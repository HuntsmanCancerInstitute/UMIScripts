#!/usr/bin/perl

# Timothy J. Parnell, PhD
# Huntsman Cancer Institute
# University of Utah
# Salt Lake City, UT 84112
#  
# This package is free software; you can redistribute it and/or modify
# it under the terms of the Artistic License 2.0.  
# 
# Updated versions of this file may be found in the repository
# https://github.com/HuntsmanCancerInstitute/UMIScripts

use strict;
use IO::File;
use IO::Handle;
use Getopt::Long;

########################################

my $VERSION = 0.1;


### Options
my $read_file1;
my $read_file2;
my $umi_file1;
my $umi_file2;
my $outfile;
my $sam_format;
my $append_name;
my $keep_sample;
my $cpu = 4;
my $sam_app = sprintf("%s", qx(which samtools));
my $pigz_app = sprintf("%s", qx(which pigz));
my $gzip_app = sprintf("%s", qx(which gzip));
my $help;

chomp $sam_app;
chomp $pigz_app;
chomp $gzip_app;


my $description = <<END;

A script to merge a separate UMI fastq file with sequence read fastq file(s). 

Unique Molecular Indexes (UMIs) are short, random nucleotide sequences that can
help to distinguish independent DNA molecules with the same sequence
composition. By associating the UMI to a nucleotide molecule at the beginning of
library preparation, identical sequence reads can be distinguished between
unique, biologically-derived molecules and PCR-derived or sequencing artifact
duplicates.

Either paired-end or single-end reads may be merged with a UMI Fastq file.

UMIs can be merged into sequence read Fastq files in three ways:

   1 Inserted as standard SAM attribute tags RX and QX in the read header.
     This is compatible with e.g. BWA, Bowties2 aligners. DEFAULT.

   2 Exported as an unaligned SAM/BAM format with RX and QX attribute tags.
     This is compatible with e.g. Bowtie2, STAR, BWA, and Novoalign.
     Use the '--bam' option.

   3 Appended to the read name as ":UMI_sequence". This is compatible with 
     all aligners but requires special software to utilize. Non-standard. 
     Use the '--name' option.

Alignments can be de-duplicated utilizing the UMI codes with external software.
See 'bam_umi_dedup.pl' in this software package, or 
Picard 'UmiAwareMarkDuplicatesWithMateCigar' as possibilities.

External utilities may be required for execution, inlcuding pigz and/or gzip, 
and samtools. These are searched for in your environment PATH.

Read and UMI Fastq files may be empirically determined based on size and name. 

END

my $usage = <<END;

VERSION: $VERSION

USAGE: 
    merge_umi_fastq.pl *.fastq.gz
    
    merge_umi_fastq.pl -1 R1.fq.gz -2 R2.fq.gz -u UMI.fq.gz -o unaligned.bam

OPTIONS:
  Input:
    -1 --read1 <file>     First fastq read
    -2 --read2 <file>     Second fastq read, optional
    -u --umi <file>       UMI fastq read
    --umi2 <file>         Second UMI fastq file, optional
  Output:
    -o --out <file>       Path to output (input basename)
                            use 'stdout' for (interleaved) piping
    -b --bam              Write output as unaligned Bam format
                            or SAM format if stdout or samtools unavailable 
  Options:
    -n --name             Append UMI to read name instead of SAM tag RX
    -k --keepbc           Keep existing Illumina CASAVA sample barcode
                            from read description if present
  Other:
    --samtools <path>     Path to samtools ($sam_app)
    -h --help             Show full description and help

END

unless (@ARGV) {
	print $usage;
	exit;
}

# get command line options
GetOptions( 
	'1|read1=s'         => \$read_file1, # first fastq file
	'2|read2=s'         => \$read_file2, # second fastq file
	'u|umi1=s'          => \$umi_file1, # first umi fastq file
	'umi2=s'            => \$umi_file2, # second umi fastq file
	'o|out=s'           => \$outfile, # output file base name
	'b|bam!'            => \$sam_format, # output sam format
	'n|name!'           => \$append_name, # append UMI to name
	'k|keepbc!'         => \$keep_sample, # keep the sample bar code
	'samtools=s'        => \$sam_app, # samtools application
	'cpu=i'             => \$cpu, # number of CPU cores for compression
	'h|help'            => \$help, # print help
) or die "bad options!\n";

if ($help) {
	print $description;
	print $usage;
	exit 0;
}





### Check options
if (not $read_file1 and not $umi_file1 and scalar @ARGV) {
	# user provided unordered list of fastq files
	# sort them out by size
	my @list = 
		sort {$b->[1] <=> $a->[1]}      # decreasing size
		map { [$_, (stat($_))[7] ] }    # array of filename and size <- 8th element from stat array
		@ARGV;
	if (scalar @list == 2) {
		# one fastq and one umi
		$read_file1 = $list[0]->[0];
		$umi_file1  = $list[1]->[0];
	}
	elsif (scalar @list == 3) {
		# two fastq and one umi
		# sort the fastq names asciibetically to get read1 and read2
		($read_file1, $read_file2) = 
			map {$_->[0]}
			sort {$a->[0] cmp $b->[0]}
			($list[0], $list[1]);
		$umi_file1  = $list[2]->[0];
	}
	elsif (scalar @list == 4) {
		# two fastq and two umi
		# sort the fastq names asciibetically to get read1 and read2
		($read_file1, $read_file2) = 
			map {$_->[0]}
			sort {$a->[0] cmp $b->[0]}
			($list[0], $list[1]);
		($umi_file1, $umi_file2) = 
			map {$_->[0]}
			sort {$a->[0] cmp $b->[0]}
			($list[2], $list[3]);
	}
	else {
		my $i = scalar @list;
		die "ERROR! $i input files provided! Unable to identify appropriately!";
	}
	print STDERR "Empirically determined following provided files:\n";
	print STDERR " - read1: $read_file1\n";
	print STDERR " - read2: $read_file2\n" if $read_file2;
	print STDERR " - umi1:  $umi_file1\n";
	print STDERR " - umi2:  $umi_file2\n" if $umi_file2;
}

if (not $read_file1) {
	die "Must provide an input read fastq file!\n";
}
if (not $umi_file1) {
	die "Must provide a UMI fastq files!\n";
}
if ($outfile and $outfile =~ /\.[s|b]am(?:\.gz)?$/) {
	# convenience
	$sam_format = 1;
}
if ($sam_format and not $sam_app) {
	if ($outfile and $outfile =~ m/\.bam$/) {
		$outfile =~ s/\bam$/sam.gz/;
		print STDERR "samtools application is not present. Writing to $outfile\n";
	}
}





### Open input filehandles
my ($read_fh1, $read_fh2, $umi_fh1, $umi_fh2);

$read_fh1 = open_input($read_file1);
$read_fh2 = open_input($read_file2) if $read_file2;
$umi_fh1  = open_input($umi_file1);
$umi_fh2  = open_input($read_file2) if $read_file2;





### Open output filehandles
my ($out_fh1, $out_fh2);
if ($outfile =~ /^stdout$/i) {
	# writing everything to standard output
	# this is overkill, but just to make things consistent, 
	# open an IO::Handle to STDOUT
	$out_fh1 = IO::Handle->new;
	$out_fh1->fdopen(fileno(STDOUT), 'w');
	if ($read_file2) {
		# reuse the same output 
		$out_fh2 = $out_fh1;
	}
}
elsif ($outfile =~ /\.bam$/ and $sam_format) {
	# special case, write to external samtools to convert
	print STDERR " Writing with $sam_app view -b -\@ $cpu to $outfile\n";
	$out_fh1 = IO::File->new("| $sam_app view -b -@ $cpu -o $outfile - ") or 
		die "unable to open output to $sam_app! $!\n";
	if ($read_file2) {
		# reuse the same output 
		$out_fh2 = $out_fh1;
	}
}
elsif ($outfile =~ /\.sam/ and $sam_format) {
	# write to sam file
	$out_fh1 = open_output($read_file1, 1); 
	if ($read_file2) {
		# reuse the same output 
		$out_fh2 = $out_fh1;
	}
}
elsif (not $outfile and $sam_format) {
	# write to sam or bam file depending on availability of samtools
	if ($sam_app) {
		# write to bam file with samtools
		$outfile = $read_file1;
		$outfile =~ s/\.(?:fq|fastq)(?:\.gz)?$/.bam/;
		print STDERR " Writing with $sam_app view -b -\@ $cpu to $outfile\n";
		$out_fh1 = IO::File->new("| $sam_app view -b -@ $cpu -o $outfile - ") or 
			die "unable to open output to $sam_app! $!\n";
		if ($read_file2) {
			# reuse the same output 
			$out_fh2 = $out_fh1;
		}
	}
	else {
		# write to text sam file
		$out_fh1 = open_output($read_file1, 1); 
		if ($read_file2) {
			# reuse the same output 
			$out_fh2 = $out_fh1;
		}
	}
}
else {
	# writing fastq text file output
	$out_fh1 = open_output($read_file1, 1); 
	if ($read_file2 and not $sam_format) {
		# ordinary second fastq file
		$out_fh2 = open_output($read_file2, 2);
	}
}


### Write sam header as necessary
if ($sam_format) {
	$out_fh1->print("\@HD\tVN:1.6\tSO:unsorted\n");
	my $cl = "$0 --read1 $read_file1 ";
	$cl .= "--read2  $read_file2 " if $read_file2;
	$cl .= "--umi $umi_file1 ";
	$cl .= "--umi2 $umi_file2 " if $umi_file2;
	$cl .= "--out $outfile " if $outfile;
	$cl .= "--bam " if $sam_format;
	$cl .= "--name " if $append_name;
	$cl .= "--keepbc " if $keep_sample;
	$out_fh1->print("\@PG\tID:merge_umi_fastq\tVN:$VERSION\tCL:$cl\n");
}




# define flags
	# 0x4d	77	PAIRED,UNMAP,MUNMAP,READ1
	# 0x8d	141	PAIRED,UNMAP,MUNMAP,READ2
	# 0x4   4   UNMAP 
my $flag1 = $read_file2 ? 77 : 4; 
my $flag2 = 141;




#### Iterate through files
my $goodCount = 0;
while (my $header1  = $read_fh1->getline) {
	my $sequence1 = $read_fh1->getline or die "malformed file $read_file1! no sequence line";
	my $spacer1 = $read_fh1->getline or die "malformed file $read_file1! no sequence line";
	my $quality1 = $read_fh1->getline or die "malformed file $read_file1! no sequence line";
	
	# first UMI
	my $umi_head1 = $umi_fh1->getline or die "malformed file $umi_file1! no sequence line";
	my $umi_seq1 = $umi_fh1->getline or die "malformed file $umi_file1! no sequence line";
	my $umi_space1 = $umi_fh1->getline or die "malformed file $umi_file1! no sequence line";
	my $umi_qual1 = $umi_fh1->getline or die "malformed file $umi_file1! no sequence line";
	
	# start chomping as needed
	chomp $header1;
	chomp $sequence1;
	chomp $quality1;
	chomp $umi_seq1;
	chomp $umi_qual1;
	
	# split the first header
	my ($r_name1, $r_desc1) = split / +/, $header1, 2;
	my ($u_name1, $u_desc1) = split / +/, $umi_head1, 2;
	
	# basic header check
	unless ($r_name1 eq $u_name1) {
		die "malformed fastq files! sequence IDs from files don't match! compare\n" . 
			" $read_file1: $header1\n $umi_file1: $umi_head1\n";
	}
		
	# second UMI if present
	my ($umi_seq, $umi_qual);
	if ($umi_file2) {
		# read
		my $umi_head2 = $umi_fh2->getline or die "malformed file $umi_file2! no sequence line";
		my $umi_seq2 = $umi_fh2->getline or die "malformed file $umi_file2! no sequence line";
		my $umi_space2 = $umi_fh2->getline or die "malformed file $umi_file2! no sequence line";
		my $umi_qual2 = $umi_fh2->getline or die "malformed file $umi_file2! no sequence line";
		
		# check
		my ($u_name2, $u_desc2) = split / +/, $umi_head2, 2;
		unless ($u_name1 eq $u_name2) {
			die "malformed UMI fastq files! sequence IDs from files don't match! compare\n" . 
				" $umi_file1: $umi_head1\n $umi_file2: $umi_head2\n";
		}
		
		# assemble based on SAM spec
		$umi_seq = join('-', $umi_seq1, $umi_seq2);
		$umi_qual = join(' ', $umi_qual1, $umi_qual2);
	}
	else {
		# only one file, easy
		$umi_seq = $umi_seq1;
		$umi_qual = $umi_qual1;
	}
	
	# compose tags
	my $new_name = $append_name ? "$r_name1:$umi_seq" : $r_name1;
	my $tags = "RX:Z:$umi_seq\tQX:Z:$umi_qual";
	if ($keep_sample and $r_desc1 =~ m/\d+:\w:\d+:([ATGCN]+[\+\-]?[ATGCN]*)/) {
		$tags .= "\tBC:Z:$1";
	}
	
	# output first line
	if ($sam_format) {
		# compose a SAM formatted line and print
		$out_fh1->printf("%s\n", join("\t", 
			substr($r_name1,1), # QNAME minus the leading @ symbol
			$flag1,             # FLAG
			'*',                # RNAME
			0,                  # POS
			0,                  # MAPQ
			'*',                # CIGAR
			'*',                # RNEXT
			0,                  # PNEXT
			0,                  # TLEN
			$sequence1,         # SEQ
			$quality1,          # QUAL
			$tags               # TAGS
		) );
	}
	elsif ($append_name) {
		# compose a fastq line with appended read name header
		$out_fh1->printf("%s:%s %s\n%s\n%s%s\n", 
			$r_name1, $umi_seq, $r_desc1,   # HEADER
			$sequence1,                     # SEQUENCE
			$spacer1,                       # SPACER
			$quality1                       # QUALITY
		);
	}
	else {
		# compose a fastq line with SAM tags
		$out_fh1->printf("%s %s\n%s\n%s%s\n", 
			$r_name1, $tags,                # HEADER
			$sequence1,                     # SEQUENCE
			$spacer1,                       # SPACER
			$quality1                       # QUALITY
		);
	}
	
	# second fastq file stuff
	if ($read_file2) {
		# read line
		my $header2  = $read_fh2->getline or die "malformed file $read_file2! no header line";
		my $sequence2 = $read_fh2->getline or die "malformed file $read_file2! no sequence line";
		my $spacer2 = $read_fh2->getline or die "malformed file $read_file2! no spacer line";
		my $quality2 = $read_fh2->getline or die "malformed file $read_file2! no quality line";
		
		# chomp as necessary
		chomp $header2;
		chomp $sequence2;
		chomp $quality2;
		
		# split the first header
		my ($r_name2, $r_desc2) = split / +/, $header2, 2;
	
		# basic header check
		unless ($r_name1 eq $r_name2) {
			die "malformed fastq files! sequence IDs from files don't match! compare\n" . 
				" $read_file1: $header1\n $read_file2: $header2\n";
		}
		
		# output second line
		if ($sam_format) {
			# compose a SAM formatted line and print
			$out_fh2->printf("%s\n", join("\t", 
				substr($r_name2,1), # QNAME minus the leading @ symbol
				$flag2,             # FLAG
				'*',                # RNAME
				0,                  # POS
				0,                  # MAPQ
				'*',                # CIGAR
				'*',                # RNEXT
				0,                  # PNEXT
				0,                  # TLEN
				$sequence2,         # SEQ
				$quality2,          # QUAL
				$tags               # TAGS
			) );
		}
		elsif ($append_name) {
			# compose a fastq line with appended read name header
			$out_fh2->printf("%s:%s %s\n%s\n%s%s\n", 
				$r_name2, $umi_seq, $r_desc2,   # HEADER
				$sequence2,                     # SEQUENCE
				$spacer2,                       # SPACER
				$quality2                       # QUALITY
			);
		}
		else {
			# compose a fastq line with SAM tags
			$out_fh2->printf("%s %s\n%s\n%s%s\n", 
				$r_name2, $tags,                # HEADER
				$sequence2,                     # SEQUENCE
				$spacer2,                       # SPACER
				$quality2                       # QUALITY
			);
		}
	}
	
	$goodCount++;
}

# close all file handles
$read_fh1->close;
$read_fh2->close if $read_fh2;
$umi_fh1->close;
$umi_fh2->close if $umi_fh2; 
$out_fh1->close;
$out_fh2->close if $out_fh2;




#### Finish
print STDERR " $goodCount reads were processed\n";




sub open_input {
	my $file = shift;
	my $fh;
	if ($file =~ /\.gz$/i) {
		# decompress with gzip
		$fh = IO::File->new("$gzip_app -d -c $file |") or 
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

sub open_output {
	my ($file, $index) = @_;
	my $filename;
	
	# check for format and output file
	if ($sam_format) {
		# writing a SAM file, likely only one
		# user may or may not have specified an output file
		if ($outfile) {
			# great! user specified
			if ($outfile =~ /\.sam(\.gz)?$/) {
				# has proper extension
				$filename = $outfile;
			}
			else {
				# add an extension
				$outfile .= '.sam.gz';
			}
		}
		else {
			# no outfile, so make one up
			$filename = $file;
			$filename =~ s/(?:txt|fq|fastq)(?:\.gz)$/sam.gz/;
		}
	}
	elsif ($outfile) {
		if ($outfile =~ /\.(fq|fastq|fq\.gz|fastq\.gz)$/) {
			# user provided extension?
			my $ext = $1;
			$filename = $outfile;
			$filename =~ s/$ext\Z//;
			$filename .= sprintf("_%s.%s", $index, $ext);
		}
		else {
			# generate with extension, assume compression
			$filename = sprintf("%s_%s.fastq.gz", $outfile, $index);
		}
	}
	else {
		# generate from input file
		$filename = $file;
		$filename =~ s/\.(?:txt|fq|fastq)(?:\.gz)$/.umi.fastq.gz/;
	}
	
	# open filehandle
	my $fh;
	if ($filename =~ /\.gz$/) {
		if ($pigz_app) {
			print STDERR " Writing with $pigz_app -p $cpu to $filename\n";
			$fh = IO::File->new("| $pigz_app -p $cpu -c > $filename") or 
				die "cannot write to compressed file '$filename' $!\n";
		}
		else {
			print STDERR " Writing with $gzip_app to $filename\n";
			$fh = IO::File->new("| $gzip_app -c > $filename") or 
				die "cannot write to compressed file '$filename' $!\n";
		}
	}
	else {
		print STDERR " Writing to $filename\n";
		$fh = IO::File->new($filename, 'w') or 
			die "cannot write to file '$filename' $!\n";
	}
	return $fh;
}





