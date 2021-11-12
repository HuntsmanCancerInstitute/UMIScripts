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
use File::Which;
use IO::File;
use List::Util qw(sum min max);
use List::MoreUtils qw(firstidx);
use Parallel::ForkManager;
use Text::Levenshtein::Flexible;
	# this Levenshtein module is preferred over others available specifically 
	# because we can weight insertions, deletions, and substitutions and
	# we only want substitutions to be considered for UMI codes
use Bio::ToolBox::db_helper 1.50 qw(
	open_db_connection 
	low_level_bam_fetch
	$BAM_ADAPTER
);
# this can import either Bio::DB::Sam or Bio::DB::HTS depending on availability
# this script is mostly bam adapter agnostic

my $VERSION = 2.1;
# version 1.0 - initial version
# version 1.1 - added option to write the duplicates to second bam file
# version 1.2 - add marking and parallel processing, make compatible with 
#               qiagen_umi_extractor script
# version 1.3 - remove the UMI in the UMI regex of the read name to work 
#               with Illumina UMI marked fastq files, allow N in UMI, 
#               write all reads without UMI codes instead of dying. 
# version 1.4 - add paired-end support and samtools cat for merging bams
# version 1.5 - improve pe support by checking each strand separately
#               this will find inverted pairs that might get tossed, breaking pairs
# version 1.6 - handle secondary alignments, write program line in bam header 
# version 1.7 - improve alignment counting methodologies
# version 1.8 - improve samtools cat command
# version 1.9 - improve help and option checking
# version 2.0 - deduplicate ALL alignments, support for SAM tags, tolerate mismatches
# version 2.1 - change default to SAM tag instead of read name

#### Inputs
my $infile;
my $outfile;
my $umi_location = 'RX';
my $markdups;
my $mismatch = 1;
my $indel_score = 10;
my $opt_distance = 0;
my $name_coordinates;
my $cpu = 4;
my $no_sam;
my $paired;
my $keep_secondary = 0;
my $help;
my $sam_app = sprintf("%s", which 'samtools');
my $tilepos;
my $xpos;
my $ypos;
my $start_time = time;

my $description = <<END;

A script to remove duplicates based on a Unique Molecular Index (UMI) code.
Alignments that match the same coordinate are checked for the random UMI
code in the read name; Reads with the same UMI are sorted and selected for
one to be retained.

UMI sequences may be embedded in alignments in one of two locations:

  1 Added as an alignment attribute tag, typically the RX tag per current
    SAM specifications. This is the standard method, but requires 
    compatible aligners. UMI sequence qualities (QX tag) are ignored.

  2 Appended to the alignment read name as ":NNNN", where NNNN is a
    sequence of indeterminate length comprised of [A,T,G,C]. This is a
    non-standard method but generally compatible with all aligners. 

At each chromosomal position, one representative alignment is selected
amongst all represented UMI sequences and the remainder are discarded
(default) or marked as duplicate (bit flag 0x400) and retained. UMI
sequence tags can tolerate mismatches up to the indicated number;
insertions or deletions are not tolerated. In general, longer UMI sequences
should tolerate more mismatches, but at the risk of missing optimal
matches. Increased tolerance results in decreased UMI-unique alignments.
Alignments without a detectable UMI flag are simply written out. 

Selection criteria amongst UMI-duplicates include the mapping qualities of
the alignment (and possibly mate pair using tag 'MQ' if present) or the sum
of base qualities of the read (and possibly mate pair using the samtools
fixmate tag 'ms').

Optical duplicate checking may be optionally included by specifying the
pixel distance threshold for marking as duplicates. Tile coordinates must
be present in the alignment read name. Optical de-duplication occurs before
UMI checking, so some false positives may be expected. Optical duplicates
are not distinguished from UMI duplicates when marking.

Bam files must be sorted by coordinate. Bam files may be indexed as
necessary. Unmapped (flag 0x4) alignments are silently discarded. Read
groups (tag RG) are ignored; de-duplication should be done with only one
read group per Bam file.

DISCLAIMER: Alignments are de-duplicated solely on alignment coordinates
and the UMI sequences of the alignments at the current position, as well as
properly paired mate pairs. No guarantees are made for maintaining the same
molecule between secondary, supplementary (chimeric), and mate pair
alignments on separate chromosomes. A unique Molecule Identifier (tag MI)
is not calculated. These results may be sufficient for most applications,
but not all.

END

my $usage = <<END;   

VERSION: $VERSION

USAGE:  bam_umi_dedup.pl --in in.bam --out out.bam

OPTIONS:
    -i --in <file>        The input bam file, should be sorted and indexed
    -o --out <file>       The output bam file
    -u --umi <string>     SAM tag name for UMI sequence. Default 'RX'
                            Specify 'name' when appended to read name.
    -m --mark             Mark duplicates (flag 0x400) instead of discarding
    -t --tolerance <int>  Mismatch tolerance ($mismatch)
    -d --distance <int>   Set optical duplicate distance threshold.
                            Use 100 for unpatterned flowcell (HiSeq) or 
                            2500 for patterned flowcell (NovaSeq). Default 0.
    --coord <string>      Provide the tile:X:Y integer 1-base positions in the 
                            read name for optical checking. For Illumina CASAVA 1.8 
                            7-element names, this is 5:6:7 (default)
    -c --cpu <int>        Specify the number of forks to use ($cpu) 
    --samtools <path>     Path to samtools ($sam_app)
    -h --help             Display full description and help

END


#### Parse Inputs
unless (@ARGV) {
	print $usage;
	exit 0;
}
GetOptions( 
	'i|in=s'        => \$infile, # the input bam file path
	'o|out=s'       => \$outfile, # name of output file 
	'u|umi=s'       => \$umi_location, # the location where the UMI tag is stored
	'm|mark!'       => \$markdups, # set duplicate flag instead of removing
	't|tolerance=i' => \$mismatch, # number mismatches allowed
	'indel=i'       => \$indel_score, # weight for indel scoring in calculating distance
	'd|distance=i'  => \$opt_distance, # optical pixel distance
	'coord=s'       => \$name_coordinates, # tile:X:Y name positions
	'c|cpu=i'       => \$cpu, # number of cpu cores to use
	'p|pe!'         => \$paired, # legacy flag for paired-end alignments
	's|secondary!'  => \$keep_secondary, # legacy flag to keep secondary alignments 
	'samtools=s'    => \$sam_app, # path to samtools
	'bam=s'         => \$BAM_ADAPTER, # specifically set the bam adapter, advanced!
	'nosam!'        => \$no_sam, # avoid using external sam adapter, advanced!
	'h|help!'       => \$help, # print help
) or die " Unrecognized option(s)!! please refer to the help documentation\n\n";

if ($help) {
	print $description;
	print $usage;
	exit 0;
}
unless ($infile) {
	die " Must provide an input file name!\n";
}
unless ($outfile) {
	die " Must provide an output file name!\n";
}
unless ($outfile =~ /\.bam$/) {
	$outfile .= '.bam';
}
undef $sam_app if $no_sam;
my $sam_version;
if ($sam_app) {
	my $sam_help = qx($sam_app 2>&1);
	if ($sam_help =~ /Version: 1\.(\d+) /) {
		$sam_version = $1;
	}
	else {
		print " Unrecognized samtools version! Disabling\n";
		undef $sam_app;
	}
}
if ($name_coordinates) {
	if ($name_coordinates =~ /(\d):(\d):(\d)/) {
		# coordinates must be converted to 0-based indices
		$tilepos = $1 - 1;
		$xpos = $2 - 1;
		$ypos = $3 - 1;
	}
	else {
		die " Name coordinate must be integers as tile:X:Y, such as 5:6:7\n";
	}
}
else {
	# these defaults are for Illumina CASAVA 1.8+ Fastq data in 0-based coordinates
	$tilepos = 4;
	$xpos = 5;
	$ypos = 6;
}

# which callback to use
my $get_umi;
if ($umi_location eq 'name') {
	print " Collecting UMI from name\n";
	$get_umi = \&get_umi_from_name;
}
elsif ($umi_location =~ /^[A-Z]{2}$/) {
	print " Collecting UMI from tag $umi_location\n";
	$get_umi = \&get_umi_from_tag;
}
else {
	die "unrecognized location '$umi_location' for storing UMI sequence!";
}





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

# update headers
my $htext = $header->text;
{
	my $pp;
	foreach my $line (split "\n", $htext) {
		# assume the last observed program should become the previous program ID
		if (substr($line,0,3) eq '@PG' and $line =~ /ID:([\w\.\-]+)/) {
			$pp =$1;
		}
	}
	$htext .= "\@PG\tID:bam_umi_dedup\t";
	$htext .= "PP:$pp\t" if $pp;
	$htext .= sprintf("VN:%s\tCL:%s", $VERSION, $0);
	$htext .= " --in $infile --out $outfile --tag $umi_location";
	$htext .= " --mark" if $markdups;
	$htext .= " --tolerance $mismatch --indel $indel_score";
	$htext .= " --distance $opt_distance" if $opt_distance;
	$htext .= " --coord $name_coordinates" if $name_coordinates and $opt_distance;
	$htext .= " --samtools $sam_app" if $sam_app;
	$htext .= " --bam $BAM_ADAPTER\n";
}





#### Deduplicate

#initialize counters
my $totalSingleCount    = 0;
my $opticalSingleCount  = 0;
my $uniqueSingleCount   = 0;
my $dupSingleCount      = 0;
my $totalPairedCount    = 0;
my $opticalPairedCount  = 0;
my $uniquePairedCount   = 0;
my $dupPairedCount      = 0;
my $untagCount          = 0;

# deduplicate in single or multi-thread
my $outbam; # I need to return the final output bam to main:: and let the normal 
          # perl exit close it properly, otherwise it crashes hard, and there's 
          # no proper close for the object!!!??????
deduplicate_multithread() or die " Something went wrong with de-duplication!\n";
undef $outbam if $outbam;






#### Finish up
printf "
 File %s:
 %12s total alignments
", $infile, ($totalSingleCount + $totalPairedCount);

if ($totalSingleCount) {
	printf  " %12s single-end alignments, including chimeric and paired-singletons
 %12s (%.1f%%) UMI-unique single-end retained
 %12s (%.1f%%) UMI-duplicate single-end %s
",
	$totalSingleCount, 
	$uniqueSingleCount, 
	($uniqueSingleCount/$totalSingleCount) * 100, 
	$dupSingleCount, 
	($dupSingleCount/$totalSingleCount) * 100, 
	$markdups ? 'marked' : 'discarded';
	if ($opticalSingleCount) {
		printf " %12s (%.1f%%) optical duplicate single-end alignments were %s\n", 
			$opticalSingleCount, ($opticalSingleCount/$totalSingleCount) * 100,
			$markdups ? 'marked' : 'discarded';
	}
}

if ($totalPairedCount) {
	printf " %12s paired-end alignments
 %12s (%.1f%%) UMI-unique paired-end retained
 %12s (%.1f%%) UMI-duplicate paired-end %s
", 
	$totalPairedCount, 
	$uniquePairedCount, 
	($uniquePairedCount/$totalPairedCount) * 100, 
	$dupPairedCount, 
	($dupPairedCount/$totalPairedCount) * 100, 
	$markdups ? 'marked' : 'discarded';
	if ($opticalPairedCount) {
		printf " %12s (%.1f%%) optical duplicate paired-end alignments were %s\n", 
			$opticalPairedCount, ($opticalPairedCount/$totalPairedCount) * 100,
			$markdups ? 'marked' : 'discarded';
	}
}

if ($untagCount) {
	printf " %12s (%.1f%%) alignments without UMI tags retained\n", $untagCount, 
		($untagCount/($totalSingleCount + $totalPairedCount)) * 100;
}

printf "\n Wrote $outfile in %.1f minutes \n", (time - $start_time) / 60;
#### End













######## Subroutines ###########


sub deduplicate_multithread {
	
	# set up manager
	my $pm = Parallel::ForkManager->new($cpu);
	$pm->run_on_finish( sub {
		my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data) = @_;
		
		# add child counts to global values
		$totalSingleCount   += $data->{totalSingleCount};
		$opticalSingleCount += $data->{opticalSingleCount};
		$uniqueSingleCount  += $data->{uniqueSingleCount};
		$dupSingleCount     += $data->{dupSingleCount};
		$totalPairedCount   += $data->{totalPairedCount};
		$opticalPairedCount += $data->{opticalPairedCount};
		$uniquePairedCount  += $data->{uniquePairedCount};
		$dupPairedCount     += $data->{dupPairedCount};
		$untagCount         += $data->{untagCount};
		
		# check coordinate errors
		if ($opt_distance and $data->{coordErrorCount} > 100000) {
			# 100K is entirely arbitrary, want to avoid random, occasional errors
			print " Too many alignments have read name tile coordinate errors!\n Disabling optical duplicate checking\n";
			# this won't have an effect on currently running forks, but at least we'll 
			# stop this for subsequent forks
			$opt_distance = 0;
		}
	});
	
	# prepare targets and names
	my @targets = (0 .. $sam->n_targets - 1);
	my $tempfile = $outfile;
	$tempfile =~ s/\.bam$//i;
	my @targetfiles = map {$tempfile . ".$$.$_.bam"} @targets;
	
	# walk through the file in parallel, one fork per chromosome
	for my $tid (@targets) {
		$pm->start and next;
		
		### in child
		
		# Clone Bam object
		$sam->clone;
		if ($BAM_ADAPTER eq 'sam') {
			$header = $sam->bam->header;
		}
		elsif ($BAM_ADAPTER eq 'hts') {
			$header = $sam->hts_file->header_read;
		}
		
		# prepare a temporary bam file
		my $tf = $targetfiles[$tid];
		my $tempbam = Bio::ToolBox::db_helper::write_new_bam_file($tf) or 
			die "unable to open output bam file $tf! $!";
		$tempbam->header_write($header);
		
		# prepare distance calculator
		my $Calculator;
		if ($mismatch) {
			# set maximum distance allowed to user-defined mismatch level
			# set scoring weights of insertion 10, deletion 10, substitution 1
			$Calculator = Text::Levenshtein::Flexible->new(
				$mismatch, $indel_score, $indel_score, 1) or 
				die "unable to initiate Text::Levenshtein::Flexible object! $!";
		}
		
		# prepare callback data structure
		my $data = {
			position            => -1,
			totalSingleCount    => 0,
			opticalSingleCount  => 0,
			uniqueSingleCount   => 0,
			dupSingleCount      => 0,
			totalPairedCount    => 0,
			opticalPairedCount  => 0,
			uniquePairedCount   => 0,
			dupPairedCount      => 0,
			untagCount          => 0,
			coordErrorCount     => 0,
			outbam              => $tempbam,
			single_reads        => [],
			paired_reads        => [],
			keepers             => {},
			duplicates          => {},
			opt_duplicates      => {},
			calculator          => $Calculator
		};
	
		# walk through the reads on the chromosome
		my $seq_length = $sam->target_len($tid);
		low_level_bam_fetch($sam, $tid, 0, $seq_length, \&callback, $data);
	
		# check to make sure we don't leave something behind
		if (defined $data->{single_reads}->[0] or defined $data->{paired_reads}->[0]) {
			write_reads($data);
		}
		
		# warnings
		if (scalar keys %{$data->{keepers}}) {
			printf " ! %d missing unmatched mate pair alignments on %s\n", 
				scalar keys %{$data->{duplicates}}, $sam->target_name($tid);
		}
		if (scalar keys %{$data->{duplicates}}) {
			printf " ! %d missing unmatched duplicate mate pair alignments on %s\n", 
				scalar keys %{$data->{duplicates}}, $sam->target_name($tid);
		}
		if (scalar keys %{$data->{opt_duplicates}}) {
			printf " ! %d missing unmatched optical-duplicate mate pair alignments on %s\n", 
				scalar keys %{$data->{opt_duplicates}}, $sam->target_name($tid);
		}
		
		# finish and return to parent
		delete $data->{position};
		delete $data->{outbam};
		delete $data->{single_reads};
		delete $data->{paired_reads};
		delete $data->{keepers};
		delete $data->{duplicates};
		delete $data->{opt_duplicates};
		delete $data->{calculator};
		undef $tempbam;
		undef $Calculator;
		$pm->finish(0, $data); 
	}
	$pm->wait_all_children;

	
	
	# attempt to use external samtools to merge the bam files
	# this should be faster than going through Perl and bam adapters
	if ($sam_app) {
		# write a temporary sam header with updated text
		my $samfile = $outfile;
		$samfile =~ s/\.bam$/temp.sam/;
		my $fh = IO::File->new($samfile, '>') or 
			die "unable to write temporary sam file\n";
		$fh->print($htext);
		$fh->close;
		
		# prepare a samtools concatenate command
		my $command = sprintf "%s cat -h %s -o %s ", $sam_app, $samfile, $outfile;
		if ($sam_version >= 10) {
			$command .= sprintf "--no-PG --threads %s ", $cpu;
		}
		$command .= join(' ', @targetfiles);
		print " executing '$sam_app cat' to merge children...\n";
		if (system($command)) {
			die "something went wrong with command '$command'!\n";
		}
		unlink $samfile, @targetfiles;
		return 1;
	}
	
	# otherwise we do this the old fashioned slow way
	# open final bam new bam file
	$outbam = Bio::ToolBox::db_helper::write_new_bam_file($outfile) or 
		die "unable to open output bam file $outfile! $!";
	# we can't write updated text to the header - otherwise it crashes
	# $header->text($htext);
	$outbam->header_write($header);

	# now remerge all the child files and write to main bam file
	if ($BAM_ADAPTER eq 'sam') {
		# open the low level bam file
		# avoid the usual high level bam file opening because that will force an 
		# automatic index, which we don't want or need
		# We will read the alignments as is. It should already be sorted and in order.
		# Must read the header first!
		for my $tf (@targetfiles) {
			my $inbam = Bio::DB::Bam->open($tf) or die "unable to open $tf!\n";
			my $h = $inbam->header;
			while (my $a = $inbam->read1) {
				&$write_alignment($outbam, $a);
			}
			undef $inbam;
			unlink $tf;
		} 
	}	
	elsif ($BAM_ADAPTER eq 'hts') {
		for my $tf (@targetfiles) {
			my $inbam = Bio::DB::HTSfile->open($tf) or die "unable to open $tf!\n";
			my $h = $inbam->header_read;
			while (my $a = $inbam->read1($h)) { 
				&$write_alignment($outbam, $a);
			}
			undef $inbam;
			unlink $tf;
		}
	}
	
	# finished
	return 1
}



### universal alignment callback
sub callback {
	my ($a, $data) = @_;
	
	# do we really need to check for QC fail FLAG 0x200 ????
	
	# First check position
	if ($a->pos > $data->{position}) {
		# current position is beyond last one
		# must dump reads to disk
		if (defined $data->{single_reads}->[0] or defined $data->{paired_reads}->[0]) {
			write_reads($data);
			# reset
			$data->{position} = $a->pos;
			$data->{single_reads} = [];
			$data->{paired_reads} = [];
		}
	}
	
	# Process alignment
	if ($a->paired) {
		# paired alignment
		
		# get supplementary status
		my $sup = $a->flag & 0x800 ? 1 : 0;
		
		# check alignment and process accordingly
		if ($a->proper_pair and not $sup) {
			# proper pair and not supplementary? great!
			push @{ $data->{paired_reads} }, $a;
			$data->{totalPairedCount}++;
		}
		elsif ($a->munmapped or $a->tid != $a->mtid or $sup or $a->unmapped) {
			# singletons or cross-chromosome alignments 
			# also don't trust that supplementary reads are "properly paired"
			# we really shouldn't have unmapped reads here
			# treat all of these as single-end reads
			push @{ $data->{single_reads} }, $a;
			$data->{totalSingleCount}++;
		}
		else {
			# other paired reads
			# typically same chromosome but miles apart or inverse orientation
			# we can still take them as pairs
			push @{ $data->{paired_reads} }, $a;
			$data->{totalPairedCount}++;
		}
	}
	else {
		# single-end alignment
		push @{ $data->{single_reads} }, $a;
		$data->{totalSingleCount}++;
	}
}


sub write_reads {
	my $data = shift;
	
	# arbitrarily write single end reads first
	write_se_reads($data) if defined $data->{single_reads}->[0];
	
	# then paired
	write_pe_reads($data) if defined $data->{paired_reads}->[0];
}


### write single-end alignments 
sub write_se_reads {
	my $data = shift;
	
	# split up based on end point and strand
	my %fends; # forward ends
	my %rends; # reverse ends
	while (my $a = shift @{$data->{single_reads}}) {
		my $end = $a->calend; 
		if ($a->reversed) {
			$rends{$end} ||= [];
			push @{ $rends{$end} }, $a;
		}
		else {
			$fends{$end} ||= [];
			push @{ $fends{$end} }, $a;
		}
	}
	
	# write forward reads
	foreach my $pos (sort {$a <=> $b} keys %fends){
		write_se_reads_on_strand($data, $fends{$pos});
	}

	# write reverse reads
	foreach my $pos (sort {$a <=> $b} keys %rends){
		write_se_reads_on_strand($data, $rends{$pos});
	}
}


sub write_se_reads_on_strand {
	my ($data, $reads) = @_; 
		# data reference and array reference of alignment reads

	# collect the non-optical duplicate reads
	my ($nonopt_reads, $optdup_reads, $error);
	if ($opt_distance) {
		($nonopt_reads, $optdup_reads, $error) = identify_optical_duplicates($reads);
		$data->{opticalSingleCount} += scalar @$optdup_reads;
		$data->{coordErrorCount} += $error;
	}
	else {
		# blindly take everything
		$nonopt_reads = $reads;
	}
	
	# put the reads into a hash based on tag identifer
	my %tag2reads;
	my @untagged;
	while (my $a = shift @{ $nonopt_reads }) {
	
		# collect UMI
		my $tag = &$get_umi($a);
		if (defined $tag) {
			$tag2reads{$tag} ||= [];
			push @{ $tag2reads{$tag} }, $a;
		}
		else {
			push @untagged, $a;
		}
	}

	# first collapse UMI tags based on mismatches
	if ($mismatch and scalar(keys %tag2reads) > 1) {
		collapse_umi_hash(\%tag2reads, $data->{calculator});
	}

	# then write out one of each tag
	foreach my $tag (keys %tag2reads) {
		if (scalar @{ $tag2reads{$tag} } == 1) {
			# good, only 1 read, we can write it
			&$write_alignment($data->{outbam}, $tag2reads{$tag}->[0] );
			$data->{uniqueSingleCount}++; # accepted read 
		}
		else {
			# damn, more than 1 read, we gotta sort
			# sorting first on mapping quality, then sum of base qualities
			my @sorted = 
				map {$_->[0]} 
				sort { ($a->[1] <=> $b->[1]) or ($a->[2] <=> $b->[2]) } 
				map { [$_, $_->qual, sum($_->qscore)] } 
				@{ $tag2reads{$tag} };
		
			# write the best one
			&$write_alignment($data->{outbam}, shift @sorted );
		
			# update counters
			$data->{uniqueSingleCount}++; # accepted read 
			$data->{dupSingleCount} += scalar(@sorted); # bad reads
		
			# write remaining as necessary
			if ($markdups) {
				foreach my $a (@sorted) {
					mark_alignment($a);
					&$write_alignment($data->{outbam}, $a);
				}
			}
		}
	}
	
	# write optical duplicate reads
	if ($markdups and scalar @$optdup_reads) {
		foreach my $a (@$optdup_reads) {
			mark_alignment($a);
			&$write_alignment($data->{outbam}, $a);
		}
	}
	
	# write untagged reads
	if (@untagged) {
		$data->{untagCount} += scalar(@untagged);
		foreach my $a (@untagged) {
			&$write_alignment($data->{outbam}, $a);
		}
	}
}


### Identify and write out paired-end alignments 
sub write_pe_reads {
	my $data = shift;
	
	# split up based on reported insertion size and strand
	my %f_sizes;
	my %r_sizes;
	while (my $a = shift @{$data->{paired_reads}}) {
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
		write_pe_reads_on_strand($data, $f_sizes{$s});
	}
	
	# write out reverse alignments
	foreach my $s (sort {$a <=> $b} keys %r_sizes) {
		write_pe_reads_on_strand($data, $r_sizes{$s});
	}
}

sub write_pe_reads_on_strand {
	my ($data, $reads) = @_; 
		# data reference and array reference of alignment reads
	
	
	# first identify possible mates of previously marked duplicates
	my @nondup_reads;
	while (my $a = shift @$reads) {
		my $name = $a->qname;
		if (exists $data->{keepers}{$name}) {
			&$write_alignment($data->{outbam}, $a);
			delete $data->{keepers}{$name};
			$data->{uniquePairedCount}++; 
		}
		elsif (exists $data->{duplicates}{$name}) {
			if ($markdups) {
				mark_alignment($a);
				&$write_alignment($data->{outbam}, $a);
			}
			delete $data->{duplicates}{$name};
			$data->{dupPairedCount}++; 
		}
		elsif (exists $data->{opt_duplicates}{$name}) {
			if ($markdups) {
				mark_alignment($a);
				&$write_alignment($data->{outbam}, $a);
			}
			delete $data->{opt_duplicates}{$name};
			$data->{opticalPairedCount}++; 
		}
		else {
			push @nondup_reads, $a;
		}
	}
	
	# collect the optical duplicate reads
	my ($nonopt_reads, $optdup_reads, $error);
	if ($opt_distance) {
		($nonopt_reads, $optdup_reads, $error) = identify_optical_duplicates(\@nondup_reads);
		$data->{optical} += scalar @$optdup_reads;
		$data->{coordErrorCount} += $error;
	}
	else {
		# blindly take everything
		$nonopt_reads = $reads;
	}
	
	# then collect the UMI codes from remaining alignments
	my %tag2reads; # barcode-to-alignment hash
	my @untagged;
	while (my $a = shift @{ $nonopt_reads }) {
		# get UMI
		my $tag = &$get_umi($a);
		if (defined $tag) {
			$tag2reads{$tag} ||= [];
			push @{ $tag2reads{$tag} }, $a;
		}
		else {
			push @untagged, $a;
		}
	}
	
	# first collapse UMI tags based on mismatches
	if ($mismatch and scalar(keys %tag2reads) > 1) {
		collapse_umi_hash(\%tag2reads, $data->{calculator});
	}
	
	# then write one of each tag
	foreach my $tag (keys %tag2reads) {
		if (scalar @{$tag2reads{$tag}} == 1) {
			# only 1 alignment, excellent, write it
			my $a = $tag2reads{$tag}->[0];
			$data->{keepers}{ $a->qname } += 1;
			&$write_alignment($data->{outbam}, $a);
			$data->{uniquePairedCount}++;
		}
		else {
			# more than one alignment, must rank them by sum of base qualities
			my @sorted = map {$_->[0]} 
				sort { ($a->[1] <=> $b->[1]) or ($a->[2] <=> $b->[2]) } 
				map { [ 
						$_, 
						sum($_->qual, $_->aux_get('MQ') || 0),
						sum($_->qscore, $_->aux_get('ms') || 0) 
				] } @{ $tag2reads{$tag} };
			
			# write the first one
			my $best = shift @sorted;
			$data->{keepers}{ $best->qname } += 1;
			&$write_alignment($data->{outbam}, $best);
			
			# increment counters
			$data->{uniquePairedCount}++;
			$data->{dupPairedCount} += scalar(@sorted);
			
			# record remaining duplicate names
			if ($markdups) {
				# write these as marked
				foreach my $a (@sorted) {
					$data->{duplicates}{$a->qname} += 1;
					mark_alignment($a);
					&$write_alignment($data->{outbam}, $a);
				}
			}
			else {
				foreach my $a (@sorted) {
					$data->{duplicates}{$a->qname} += 1;
				}
			}
		}
	}
	
	# write optical duplicate reads, and remember them so to write its pair
	if ($markdups and scalar @$optdup_reads) {
		foreach my $a (@$optdup_reads) {
			$data->{opt_duplicates}{$a->qname} += 1;
			mark_alignment($a);
			&$write_alignment($data->{outbam}, $a);
		}
	}
	elsif (scalar @$optdup_reads) {
		foreach my $a (@$optdup_reads) {
			$data->{opt_duplicates}{$a->qname} += 1;
		}
	}
	
	# write untagged reads, and remember them so to write its pair
	if (@untagged) {
		$data->{untagCount} += scalar(@untagged);
		foreach my $a (@untagged) {
			$data->{keepers}{$a->qname} += 1; # remember to write pair
			&$write_alignment($data->{outbam}, $a);
		}
	}
}


sub identify_optical_duplicates {
	my $alignments = shift;
	my $data = shift || undef; # this is only for reporting distances
	
	# check for only one alignment and quickly return
	if (scalar @$alignments == 1) {
		return ($alignments, []);
	}
	
	# extract tile and coordinates from each alignment given and sort 
	my %tile2alignment;
	my @dups;
	my @nondups;
	my $coord_error = 0;
	foreach my $a (@$alignments) {
		my @bits = split(':', $a->qname);
		if (
			$#bits < $ypos or (
				not defined $bits[$tilepos] and 
				not defined $bits[$xpos] and 
				not defined $bits[$ypos]
			)
		) {
			# not able to identify coordinates in the alignment name
			# record error count and place in nondups array
			$coord_error++;
			push @nondups, $a;
		}
		else {
			my $rg = $a->aux_get("RG") || '1';
			my $rgtile = join('.', $rg, $bits[$tilepos]);
			$tile2alignment{$rgtile} ||= [];
			# each item in the array is another array of [x, y, $a]
			push @{$tile2alignment{$rgtile}}, [$bits[$xpos], $bits[$ypos], $a];
		}
	}
	
	# check whether we have actionable alignments
	if (scalar keys %tile2alignment == 0) {
		# couldn't find any coordinates!
		return (\@nondups, \@dups, $coord_error);
	}
	
	
	
	# walk through each read-group tile
	while (my ($rgtile, $rgtile_alnts) = each %tile2alignment) {
		# check the number of alignments
		if (scalar @$rgtile_alnts == 1) {
			# we only have one, so it's a nondup
			push @nondups, $rgtile_alnts->[0][2];
		}
		else {
			# we have more than one so must sort through
			
			# sort by increasing X coordinate
			my @spots = sort {$a->[0] <=> $b->[0]} @$rgtile_alnts;
			# collect the deltas from one spot to the next on the X axis
			# start with second sorted spot and subtract the one before it
			my @xdiffs = map { $spots[$_]->[0] - $spots[$_-1]->[0] } (1..$#spots);
			
			
			
			## Check 
			# working on the assumption here that by definition an optical cluster is 
			# going to be below the threshold on both axes, and not just one axis
			# therefore we only need to check the x axis first
			if (min(@xdiffs) > $opt_distance) {
				# everything is ok
				foreach (@spots) {
					push @nondups, $_->[2];
				}
			}
			else {
				# two or more spots are too close on X axis, must check out
				my $first = shift @spots;
				while (@spots) {
					if ($spots[0]->[0] - $first->[0] < $opt_distance) {
						# second is too close to first
						# keep taking spots until they're no longer close
						my @closest;
						push @closest, $first, shift @spots;
						# continue comparing the next one with the previous one
						while (
							scalar @spots and 
							$spots[0]->[0] - $closest[-1]->[0] < $opt_distance
						) {
							push @closest, shift @spots;
						}
						
						## Now must process this X cluster
						# check Y coordinates by sorting and calculating deltas
						@closest = sort {$a->[1] <=> $b->[1]} @closest;
						my @ydiffs = map { $closest[$_]->[1] - $closest[$_-1]->[1] } (1..$#closest);
						
						# check the Y deltas in this X cluster
						if (min(@ydiffs) > $opt_distance) {
							# they're all good
							foreach (@closest) {
								push @nondups, $_->[2];
							}
						}
						else {
							# we definitely have some XY clusters within threshold
							my $xyfirst = shift @closest;
							while (@closest) {
								if ($closest[0]->[1] - $xyfirst->[1] < $opt_distance) {
									# these two are close
									my @clustered;
									push @clustered, $xyfirst, shift @closest;
									# continue compareing the next one with the previous one
									while (
										scalar @closest and 
										$closest[0]->[1] - $clustered[-1]->[1] < $opt_distance
									) {
										push @clustered, shift @closest;
									}
									# done, no more
									
									# take the first one as pseudo random chosen one
									push @nondups, $clustered[0]->[2];
									# remainder are optical duplicates
									for my $i (1..$#clustered) {
										push @dups, $clustered[$i]->[2];
									}
									
									## prepare for next round
									$xyfirst = shift @closest || undef;
								}
								else {
									# this one is ok
									push @nondups, $xyfirst->[2];
									$xyfirst = shift @closest;
								}
							}
							
							# check for last remaining alignment
							push @nondups, $xyfirst->[2] if defined $xyfirst;
						}
						
						## Prepare for next round
						$first = shift(@spots) || undef;
					}
					else {
						# first is ok
						push @nondups, $first->[2];
						$first = shift @spots;
					}
					# continue
				}
				
				# check for last remaining alignment
				push @nondups, $first->[2] if defined $first;
			}
		}
	}
	
	# finished
	return (\@nondups, \@dups, $coord_error);
}


sub get_umi_from_name {
	# get UMI from read name
	if ($_[0]->qname =~ /:([AGCTN]{4,25})([\+\-][AGCTN]{4,25})?/i) {
		if ($2) {
			return $1 . $2;
		}
		else {
			return $1;
		}
	}
	else {
		return undef;
	}
}

sub get_umi_from_tag {
	# get UMI from alignment tag
	return $_[0]->aux_get($umi_location) || undef;
}

sub collapse_umi_hash {
	my ($tag2a, $Calculator) = @_;
	# create first list from all the available UMI sequences
	# we always sort the list so that we can merge similar tags in a systematic 
	# fashion rather than randomly, shouldn't matter for single mismatches, but for 
	# two (or more!?) mismatches, you can imagine hierarchical trees of relatedness
	# but that's way too computationally intense for doing here, hence sorting is better 
	# than nothing
	my @list = sort {$a cmp $b} keys %$tag2a;
	while (scalar(@list) > 1) {
		# compare first tag with the remainder
		my $first = shift @list;
		my @results = $Calculator->distance_lc_all($first, @list);
		
		# check results
		if (@results) {
			# combine alignments for matching tags, then delete the hash entry
			foreach (@results) {
				my $t = $_->[0]; # result is array of [word, #mismatches]
				push @{ $tag2a->{$first} }, @{$tag2a->{$t}};
				delete $tag2a->{$t};
				# remove this from the list
				my $i = firstidx {$_ eq $t} @list;
				if (defined $i) {
					splice @list, $i, 1;
				}
			}
		}
	}
	return 1;
}

sub write_sam_alignment {
	# wrapper for writing Bio::DB::Sam alignments
	return $_[0]->write1($_[1]);
}

sub write_hts_alignment {
	# wrapper for writing Bio::DB::HTS alignments
	return $_[0]->write1($header, $_[1]);
}

sub mark_alignment {
	# set alignment flag as a duplicate
	my $f = $_[0]->flag;
	unless ($f & 0x400) {
		$f += 0x400;
		$_[0]->flag($f);
	}
}



