#!/usr/bin/perl

use strict;
use Bio::ToolBox::db_helper::bam; 

# version 1.0 - initial version
# version 1.1 - added option to write the duplicates to second bam file


unless (@ARGV) {
	print <<END;

A script to remove duplicates from ChIP-Nexus barcoded Bam alignment files.
Alignments that match the same coordinate are checked for the random barcode 
in the read name; Reads with the same barcode are sorted and selected for one 
to be retained. Selection criteria include mapping quality and the sum of read 
base qualities. Retained reads are written to the output bam file. 

See the script ChIPNexus_fastq_barcode_processer.pl for pre-processing the 
Fastq files prior to alignments, which will append the random barcode to the 
read name.

Bam files must be sorted by coordinate. Bam files will be indexed as necessary.

Optionally provide a second output file for the discarded duplicate alignments.

Usage: $0 <input.bam> <output.bam>
       $0 <input.bam> <output.bam> <duplicates.bam>

END
	exit;
}

# input bam file
my $infile = shift @ARGV;
my $sam = open_bam_db($infile) or # automatically takes care of indexing
	die "unable to read bam $infile $!";

# output bam file
my $outfile = shift @ARGV or die "no output file provided!\n";
my $outbam = write_new_bam_file($outfile) or 
	die "unable to open output bam file $outfile! $!";
	# this uses low level Bio::DB::Bam object

# duplicates bam file
my $dupfile = shift @ARGV || undef;
my $dupbam;
if ($dupfile) {
	$dupbam = write_new_bam_file($dupfile) or 
		die "unable to open output duplicate bam file $dupfile! $!";
		# this uses low level Bio::DB::Bam object
}

# write header
my $header = $sam->bam->header;
$outbam->header_write($header);
if ($dupbam) {
	$dupbam->header_write($header);
}

#initialize counters
my $totalCount = 0;
my $goodCount  = 0;
my $badCount   = 0;


# process on chromosomes
# this could be parallelized per chromosome to speed up, if this is too slow
for my $tid (0 .. $sam->n_targets - 1) {
	my $seq_length = $sam->target_len($tid);
	
	# prepare callback data structure
	# anonymous hash of 
	my $data = {
		position   => 0,
		reads      => [],
	};
	
	# walk through the reads on the chromosome
	$sam->bam_index->fetch($sam->bam, $tid, 0, $seq_length, \&callback, $data);
	
	# check to make sure we don't leave something behind
	if ($data->{reads}->[0]) {
		write_reads($data);
	}
}

# finish up
printf "
 %12s total mapped reads
 %12s unique mapped reads retained
 %12s duplicate mapped reads discarded
", $totalCount, $goodCount, $badCount;
undef $outbam;
check_bam_index($outfile);
if ($dupbam) {
	undef $dupbam;
	check_bam_index($dupfile);
}

# end


### alignment callback
sub callback {
	my ($a, $data) = @_;
	$totalCount += 1;
	
	# check position 
	my $pos = $a->pos;
	if ($pos == $data->{position}) {
		push @{ $data->{reads} }, $a;
	}
	else {
		write_reads($data);
		$data->{position} = $pos;
		push @{ $data->{reads} }, $a;
	}
}

### write passing alignments 
sub write_reads {
	my $data = shift;
	
	# put the reads into a hash based on tag identifer
	my %tag2reads;
	foreach my $a (@{ $data->{reads} }) {
		my $tag;
		if ($a->qname =~ /:([AGTC]+)$/) {
			$tag = $1;
		}
		else {
			die sprintf 
				" unable to identify bar code string at end of read '%s'! Please pre-process your reads with ChIPNexus_fastq_barcode_processer.pl to append the bar code sequence to the read name as :NNNNN", 
				$a->qname;
		}
		$tag2reads{$tag} ||= [];
		push @{ $tag2reads{$tag} }, $a;
	}
	
	# write out the unique reads
	foreach my $tag (keys %tag2reads) {
		if (scalar @{ $tag2reads{$tag} } == 1) {
			# good, only 1 read, we can write it
			$outbam->write1( $tag2reads{$tag}->[0] );
			$goodCount++; # accepted read 
		}
		else {
			# damn, more than 1 read, we gotta sort
			# sorting first on mapping quality, then sum of base qualities
			my @sorted = 
				map {$_->[0]} 
				sort { ($a->[1] <=> $b->[1]) or ($a->[2] <=> $b->[2]) } 
				map { [$_, $_->qual, sum($_->qscore)] } 
				@{ $tag2reads{$tag} };
			$outbam->write1( shift @sorted );
			$goodCount++; # accepted read 
			$badCount += scalar(@sorted); # bad reads
			if ($dupbam) {
				foreach my $a (@sorted) {
					$dupbam->write1($a);
				}
			}
		}
	}
	
	# reset
	$data->{reads} = [];
}

sub sum {
	my $sum = 0;
	foreach (@_) { $sum += $_ }
	return $sum;
}


