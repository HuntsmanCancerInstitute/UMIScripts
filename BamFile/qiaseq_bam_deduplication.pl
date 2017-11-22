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
# https://github.com/tjparnell/HCI-Scripts/

use strict;
use Getopt::Long;
use List::Util qw(sum max);
use Bio::ToolBox::db_helper 1.50 qw(
	open_db_connection 
	low_level_bam_fetch
	$BAM_ADAPTER
);
# this can import either Bio::DB::Sam or Bio::DB::HTS depending on availability
# this script is mostly bam adapter agnostic

### Get Options
my ($infile, $outfile, $dupfile);
my $min_mapq = 13;
my $max_isize = 2000;
my $bc_length = 12;

unless (@ARGV) {
print <<DOC;

A simple barcode-aware de-duplication application for bam files. 

This is intended for use with the Qiaseq application, but may work with other 
sequencing libraries with molecular barcodes. The application will identify proper 
alignment pairs with identical coordinates, and then compare barcodes of each pair 
stored in the read name. Only one unique barcode is retained for each position. 
Duplicate barcodes are considered true PCR duplicates and discarded, or optionally 
written to a second bam file. The best barcode alignment is selected by the 
highest mapping quality or the highest sum of base qualities.

Barcodes from Qiaseq libraries can be extracted from the second read using the 
companion application qiaseq_FastqBarcodeTagger.pl. Barcodes from a third sequencing 
application may be processed with the USeq application FastqBarcodeTagger. These 
insert the barcode and barcode qualities into each pair read name. 

Note that this application does not check barcode sequence quality. No attempt to 
collapse the duplicate sequence reads is performed; the "best" one is simply 
retained and the others discarded. For low-levels of duplication, this may be 
acceptable.

Supplemental and secondary alignments are silently discarded. Alignment pairs 
that do not match the expected F - R orientation on the same chromosome and within 
the indicated distance are silently discarded. Singletons are also silently discarded. 

Currently only supports paired-end alignment files.

USAGE:  qiaseq_bam_deduplication.pl --in in.bam --out out.bam
       
OPTIONS:
    --in <file>      The input bam file, should be sorted and indexed
    --out <file>     The output bam file
    --dup <file>     The duplicates bam file, optional if you want to keep 
                     the duplicate alignments. Non-standard alignments 
                     are not retained.
    --len <integer>  Length of the molecular barcode. Default $bc_length.
    --mapq <integer> The minimum mapping quality to retain. 
                     Default is $min_mapq.
    --size <integer> The minimum paired-insertion size to retain.
                     Default is $max_isize bp.

DOC
	exit 0;
}

GetOptions( 
	'in=s'       => \$infile, # the input bam file path
	'out=s'      => \$outfile, # name of output file 
	'dup=s'      => \$dupfile, # the name of the duplicates file
	'mapq=i'     => \$min_mapq, # minimum mapping quality allowed
	'size=i'     => \$max_isize, # maximum insertion size allowed
	'bam=s'      => \$BAM_ADAPTER, # specifically set the bam adapter, advanced!
) or die " unrecognized option(s)!! please refer to the help documentation\n\n";



### Open bam files
# input bam file
my $sam = open_db_connection($infile) or # automatically takes care of indexing
	die "unable to read bam $infile $!";

# read header and set adapter specific alignment writer
my ($header, $write_alignment);
if ($BAM_ADAPTER eq 'sam') {
	$header = $sam->bam->header;
	$write_alignment = \&write_sam_alignment;
}
elsif ($BAM_ADAPTER eq 'hts') {
	$header = $sam->hts_file->header_read;
	$write_alignment = \&write_hts_alignment;
}
else {
	die "unrecognized bam adapter $BAM_ADAPTER!";
}

# open bam new bam file
$outfile .= '.bam' unless $outfile =~ /\.bam$/i;
my $outbam = Bio::ToolBox::db_helper::write_new_bam_file($outfile) or 
	die "unable to open output bam file $outfile! $!";
	# this uses low level Bio::DB::Bam object
	# using an unexported subroutine imported as necessary depending on bam availability
$outbam->header_write($header);

# open duplicate bam file if requested
my $dupbam;
if ($dupfile) {
	$dupfile .= '.bam' unless $dupfile =~ /\.bam$/i;
	$dupbam = Bio::ToolBox::db_helper::write_new_bam_file($dupfile) or 
		die "unable to open duplicate bam file $dupfile! $!";
		# this uses low level Bio::DB::Bam object
		# using an unexported subroutine imported as necessary depending on bam availability
	$dupbam->header_write($header);
}

# set counters
my $totalCount     = 0;
my $uniqueCount    = 0;
my $duplicateCount = 0;
my $tossCount      = 0;
my $unpairedCount  = 0;
my %dupDepth2Count;

# walk through bam file again
for my $tid (0 .. $sam->n_targets - 1) {
	my $seq_length = $sam->target_len($tid);
	
	# prepare callback data structure
	my $data = {
		position   => -1,
		reads      => [],
		keepers    => {}, # hash of names of rev pe reads to keep
		dupkeepers => {}, # hash of names of rev pe dup reads to keep
	};
	
	# walk through the reads on the chromosome
	low_level_bam_fetch($sam, $tid, 0, $seq_length, \&pe_callback, $data);
	
	# check to make sure we don't leave something behind
	# possible with single-end, should be unlikely with paired-end
	if (defined $data->{reads}->[0]) {
		write_out_alignments($data);
	}
	if (scalar keys %{$data->{keepers}}) {
		$unpairedCount += scalar(keys %{$data->{keepers}});
	}
}




### Print results
printf "  Total mapped: %22d 
  Unique count: %22d
  Retained duplicate count: %10d 
  Removed duplicate count: %11d
  Unpaired count: %20d\n",
	$totalCount, $uniqueCount, $duplicateCount, $tossCount, $unpairedCount;
print " Duplicate Depth counts\n";
foreach (sort {$a <=> $b} keys %dupDepth2Count) {
	printf "  %d\t%d\n", $_, $dupDepth2Count{$_};
}



### Finish up
printf " Wrote %d pairs of alignments to $outfile\n", 
	($uniqueCount + $duplicateCount) * 2;
printf " Wrote %d pairs of alignments to $dupfile\n", $tossCount * 2
	if $dupfile;
exit; # bam files should automatically be closed





### paired-end alignment callback for writing
sub pe_callback {
	my ($a, $data) = @_;
	
	## check the alignment
	return unless $a->proper_pair; # consider only proper pair alignments
	return unless $a->tid == $a->mtid;
	my $isize = $a->isize; 
	if ($a->reversed) {
		# in proper FR pairs reverse alignments are negative
		return if $isize > 0; # pair is RF orientation
		$isize = abs($isize);
	}
	return if $isize > $max_isize;
	return if $a->qual < $min_mapq; # mapping quality
	my $flag = $a->flag;
	return if ($flag & 0x0100); # secondary alignment
	return if ($flag & 0x0200); # QC failed but still aligned? is this necessary?
	return if ($flag & 0x0800); # supplementary hit
	# we ignore duplicate marks, assume it's not marked anyway....
	
	# process forward reads
	$totalCount++ unless $a->reversed; # only count forward alignments
	
	# check position 
	my $pos = $a->pos;
	if ($pos == $data->{position}) {
		push @{ $data->{reads} }, $a;
	}
	else {
		# write out the reads
		write_out_alignments($data) if defined $data->{reads}->[0];
		
		# reset
		$data->{reads} = [];
		$data->{position} = $pos;
		push @{ $data->{reads} }, $a;
	}
}


### Identify and write out paired-end alignments 
sub write_out_alignments {
	my $data = shift;
	
	# split up based on reported insertion size and strand
	my %f_sizes;
	my %r_sizes;
	while (my $a = shift @{$data->{reads}}) {
		my $s = $a->isize; 
		if ($a->reversed) {
			$r_sizes{$s} ||= [];
			push @{$r_sizes{$s}}, $a;
		}
		else {
			$f_sizes{$s} ||= [];
			push @{$f_sizes{$s}}, $a;
		}
	}
	
	# write out forward alignments
	foreach my $s (sort {$a <=> $b} keys %f_sizes) {
		my $number = scalar @{ $f_sizes{$s} };
		if ($number == 1) {
			# only one, write it
			my $a = $f_sizes{$s}->[0];
			$data->{keepers}{ $a->qname } += 1;
			&$write_alignment($outbam, $a);
			$uniqueCount++;
		}
		else {
			# must sort the bams by the barcode
			$dupDepth2Count{$number} += 1;
			my %bc2a; # barcode-to-alignment hash
			for my $i (0 .. $number - 1) {
				my $a = $f_sizes{$s}->[$i];
				if ($a->qname =~ /:BMF:([ACGT]{$bc_length})/) {
					$bc2a{$1} ||= [];
					push @{ $bc2a{$1} }, $a;
				}
			}
			
			# now write one of each
			foreach my $bc (keys %bc2a) {
				if (scalar @{$bc2a{$bc}} == 1) {
					# only 1 alignment, excellent, write it
					my $a = $bc2a{$bc}->[0];
					$data->{keepers}{ $a->qname } += 1;
					&$write_alignment($outbam, $a);
					$uniqueCount++;
				}
				else {
					# more than one alignment, must rank them by sum of base qualities
					my @sorted = 
						map {$_->[0]} 
						sort { ($a->[1] <=> $b->[1]) or ($a->[2] <=> $b->[2]) } 
						map { [$_, $_->qual, sum($_->qscore)] } 
						@{ $bc2a{$bc} };
					# write the first one
					my $best = shift @sorted;
					$data->{keepers}{ $best->qname } += 1;
					&$write_alignment($outbam, $best);
					
					# increment counters
					$duplicateCount++;
					$tossCount += scalar(@sorted);
					
					# write the duplicates
					if ($dupbam) {
						while (@sorted) {
							my $a = shift @sorted;
							$data->{dupkeepers}{$a->qname};
							&$write_alignment($dupbam, $a);
						}
					}
				}
			}
		}
	}
	
	# write out reverse alignments
	foreach my $s (sort {$a <=> $b} keys %r_sizes) {
		while (my $a = shift @{ $r_sizes{$s} }) {
			my $name = $a->qname;
			if (exists $data->{keepers}{$name}) {
				# we have processed the forward alignment as a keeper
				# immediately write the reverse alignment
				&$write_alignment($outbam, $a);
				delete $data->{keepers}{$name};
			}
			elsif (exists $data->{dupkeepers}{$name}) {
				# we have processed the forward alignment as a duplicate keeper
				# immediately write the reverse alignment
				&$write_alignment($dupbam, $a);
				delete $data->{dupkeepers}{$name};
			}
		}
	}
}


sub write_sam_alignment {
	# wrapper for writing Bio::DB::Sam alignments
	return $_[0]->write1($_[1]);
}

sub write_hts_alignment {
	# wrapper for writing Bio::DB::HTS alignments
	return $_[0]->write1($header, $_[1]);
}








