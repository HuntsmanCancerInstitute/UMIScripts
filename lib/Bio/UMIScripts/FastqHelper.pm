package Bio::UMIScripts::FastqHelper;
our $VERSION = 5.01;

=head1 Bio::UMIScripts::FastqHelper - Fastq helper routines for UMIScripts

=head1 DESCRIPTION

Exported subroutines for working with Fastq files in the UMIScript package.
Not necessarily intended for general consumption.

=head2 Note on design

This package could easily be written with an object-oriented design and 
interface. However, in an interest to eke the maximum amount of performance 
with very large files and rapid throughput, it was intentionally written 
with simple exported subroutines and named indexes to array references to 
avoid as much overhead as possible. This is not necessarily the best 
approach.

=head1 METHODS

All subroutines are exported by default.

=head2 Subroutines

=over 4

=item read_fastq_filehandle($filename)

Pass filename to open and read from. Will open compressed files with external utility. 
Expects usual file extensions. Returns an L<IO::File> object;

=item write_fastq_filehandle($filename, $cpu)

Pass a filename to open a file for writing. Optionally pass the number of CPU cores
for threaded compression by external utility. The filename may include a F<.gz> 
extension, in which case it will be passed through an external compression utility.

Compressed fastq text files are written using either external F<pigz> multi-threaded 
utility with the indicated number of threads (default 1), or if not available, 
standard F<gzip> utility. 

A file name of C<stdout> opens a file handle to C<STDOUT>.

Returns L<IO::File> or L<IO::Handle> (for standard out) object.

=item write_bam_file_handle($filename, $cpu)

Pass a Bam filename to open a file for writing. The filename must have a F<.bam> 
extension. The file handle will be piped to the external F<samtools> application 
for writing. If F<samtools> is not in the current C<PATH> and found automatically, 
then the path may be explicitly set using the variable 
C<$Bio::UMIScripts::FastqHelper::SAMTOOLS_APP>. If F<samtools> is not present, 
then a compressed, text C<SAM> file will be written using the extension F<.sam.gz>.

Optionally pass an integer for the number of compression threads (default 1). 

=item get_fastq_read($fh)

Pass an opened L<IO::File> object for reading. Reads next four lines. Will die 
if four lines are not available. The header line is automatically split into C<name> 
and C<description>. All lines except for spacer are C<chomp>ed.

Returns an array reference of five items representing the read: name, description, 
sequence, spacer, and quality. These may be accessed by exported index names. See 
L</"Exported Variables"> below.

=back

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Bioinformatics Shared Resource
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  

=cut

use strict;
use IO::File;
use IO::Handle;
use List::Util qw(min);
use Bio::UMIScripts::FastqRead;
require Exporter;

our $GZIP_APP;
our $PIGZ_APP;
our $SAMTOOLS_APP;


BEGIN {
	# get paths to help applications
	$SAMTOOLS_APP = qx(which samtools);
	chomp $SAMTOOLS_APP;
	$GZIP_APP = qx(which gzip);
	chomp $GZIP_APP;
	$PIGZ_APP = qx(which pigz);
	chomp $PIGZ_APP;
}

our @ISA = qw(Exporter);
our @EXPORT_OK = qw(
	read_fastq_filehandle
	write_fastq_filehandle
	write_bam_filehandle
	get_fastq_read
	$SAMTOOLS_APP
);

sub read_fastq_filehandle {
	my $file = $_[0];
	unless ($file =~ /\.(?:fastq|fq|txt)(?:\.gz|\.bz2)?$/) {
		print STDERR "file $file does not have expected file extension!\n";
	}
	my $fh;
	if ($file =~ /\.gz$/i) {
		# decompress with gzip
		$fh = IO::File->new("$GZIP_APP -d -c $file |") or
			die "unable to open file $file! $!\n";
	}
	elsif ($file =~ /\.bz2$/) {
		# hope for the best????
		$fh = IO::File->new("bzip2 -d -c $file |") or
			die "unable to open file $file! $!\n";
	}
	else {
		# why do we have an uncompressed file???? geez! 
		$fh = IO::File->new($file) or
			die "unable to open file $file! $!\n";
	}
	return $fh;
}

sub write_fastq_filehandle {
	my ($outfile, $cpu) = @_;
	$cpu ||= 3;
	
	# open filehandle
	my $fh;
	if (lc($outfile) eq 'stdout') {
		$fh = IO::Handle->new;
		$fh->fdopen(fileno(STDOUT), 'w');
	}
	elsif ($outfile =~ /\.gz$/) {
		if ($PIGZ_APP) {
			print STDERR " Writing with $PIGZ_APP -p $cpu to $outfile\n";
			$fh = IO::File->new("| $PIGZ_APP -p $cpu -c > $outfile") or
				die "cannot write to compressed file '$outfile' $!\n";
		}
		elsif ($GZIP_APP) {
			print STDERR " Writing with $GZIP_APP to $outfile\n";
			$fh = IO::File->new("| $GZIP_APP -c > $outfile") or
				die "cannot write to compressed file '$outfile' $!\n";
		}
		else {
			print STDERR " Gzip applications not available!\n";
			$outfile =~ s/\.gz$//;
			print STDERR " Writing to $outfile\n";
			$fh = IO::File->new($outfile, 'w') or
				die "cannot write to file '$outfile' $!\n";
		}
	}
	else {
		print STDERR " Writing to $outfile\n";
		$fh = IO::File->new($outfile, 'w') or
			die "cannot write to file '$outfile' $!\n";
	}
	return $fh;
}

sub write_bam_filehandle {
	my ($outfile, $cpu) = @_;
	$cpu ||= 3;
	
	# check that we have supported external helper applications
	if ($outfile =~ /\.bam$/) {
		unless ($SAMTOOLS_APP) {
			print STDERR " Samtools not installed - writing text SAM file\n";
			if ($PIGZ_APP or $GZIP_APP) {
				$outfile =~ s/\.bam$/.sam.gz/;
			}
			else {
				# ugh, uncompressed sam file, really!!???
				$outfile =~ s/\.bam$/.sam/;
			}
			return write_fastq_filehandle($outfile, $cpu);
		}
	}
	
	# open filehandle
	if ($outfile =~ /\.bam$/) {
		print STDERR " Writing with $SAMTOOLS_APP view -b -\@ $cpu to $outfile\n";
		my $fh = IO::File->new("| $SAMTOOLS_APP view -b -@ $cpu -o $outfile - ") or
			die "unable to open output to $SAMTOOLS_APP! $!\n";
		$fh->print("\@HD\tVN:1.6\tSO:unsorted\n");
		return $fh;
	}
	else {
		# not a bam file!? use the other one
		return write_fastq_filehandle($outfile);
	}
}

sub get_fastq_read {
	my $header  = $_[0]->getline;
	return unless $header;
	unless (substr($header, 0, 1) eq '@') {
		die "Line is not a proper fastq line!\n$header\n";
	}
	my $sequence = $_[0]->getline or die "malformed file! no sequence line";
	my $spacer   = $_[0]->getline or die "malformed file! no spacer line";
	my $quality  = $_[0]->getline or die "malformed file! no quality line";
	chomp $header;
	chomp $sequence;
	chomp $spacer;
	chomp $quality;
	my ($name, $description) = split '\s+', $header, 2;
	return Bio::UMIScripts::FastqRead->new(
		$name,
		$description,
		$sequence,
		$spacer,
		$quality
	);
}

1;
