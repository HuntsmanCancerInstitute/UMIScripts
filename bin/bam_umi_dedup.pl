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
use List::Util qw(sum max);
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

my $VERSION = 2.0;
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

#### Inputs
my $infile;
my $outfile;
my $umi_location = 'name';
my $markdups;
my $mismatch = 1;
my $indel_score = 10;
my $paired;
my $keep_secondary = 0;
my $cpu = 4;
my $no_sam;
my $help;
my $sam_app = sprintf("%s", which 'samtools');

my $description = <<END;

A script to remove duplicates based on a Unique Molecular Index (UMI) code.
Alignments that match the same coordinate are checked for the random UMI code 
in the read name; Reads with the same UMI are sorted and selected for one 
to be retained. 

UMI sequences may be embedded in alignments in one of two locations:

  1 Added as an alignment attribute tag, typically the RX tag per 
    current SAM specifications. This is the standard method, but 
    requires compatible aligners. UMI sequence qualities (QX tag) 
    are ignored. 

  2 Appended to the alignment read name as ":NNNN", where NNNN is a 
    sequence of indeterminate length comprised of [A,T,G,C]. This 
    is a non-standard method but generally compatible with all aligners.
    Currently the default method for legacy reasons.

At each chromosomal position, one representative alignment is selected amongst
all represented UMI sequences and the remainder are discarded (default) or
marked as duplicate (bit flag 0x400) and retained. UMI sequence tags can
tolerate mismatches up to the indicated number; insertions or deletions are not
tolerated. In general, longer UMI sequences should tolerate more mismatches.
Increased tolerance results in decreased UMI-unique alignments. Alignments
without a detectable UMI flag are simply written out.

For single-end alignments, selection criteria amongst UMI-duplicates include
mapping quality and the sum of read base qualities. For paired-end alignments,
selection criteria include the mapping qualities of the forward and possibly
reverse mate alignments (using tag 'MQ' if present) and the base qualities of
the forward and possibly reverse mate alignments (using the samtools fixmate tag
'ms'). 

Bam files must be sorted by coordinate. Bam files may be indexed as necessary.
Unmapped (flag 0x4) alignments are silently discarded. Read groups (tag RG) are
ignored; de-duplication should be done with only one read group per Bam file.

DISCLAIMER: Alignments are de-duplicated solely on alignment coordinates and the
UMI sequences of the alignments at the current position, as well as properly
paired mate pairs. No guarantees are made for maintaining the same molecule
between secondary, supplementary (chimeric), and mate pair alignments on
separate chromosomes. A unique Molecule Identifier (tag MI) is not calculated.
These results may be sufficient for most applications, but not all.

END

my $usage = <<END;   

VERSION: $VERSION

USAGE:  bam_umi_dedup.pl --in in.bam --out out.bam

OPTIONS:
    -i --in <file>        The input bam file, should be sorted and indexed
    -o --out <file>       The output bam file
    -u --umi <string>     The location of the UMI sequence. See description.
                            Typically 'RX'. Default 'name'.
    -m --mark             Mark duplicates instead of discarding
    -t --tolerance <int>  Mismatch tolerance ($mismatch)
    -c --cpu <int>        Specify the number of threads to use ($cpu) 
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
	'p|pe!'         => \$paired, # legacy flag for paired-end alignments
	'm|mark!'       => \$markdups, # set duplicate flag instead of removing
	't|tolerance=i' => \$mismatch, # number mismatches allowed
	'indel=i'       => \$indel_score, # weight for indel scoring in calculating distance
	's|secondary!'  => \$keep_secondary, # legacy flag to keep secondary alignments 
	'c|cpu=i'       => \$cpu, # number of cpu cores to use
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
	$htext .= " --samtools $sam_app" if $sam_app;
	$htext .= " --bam $BAM_ADAPTER\n";
}





#### Deduplicate

#initialize counters
my $totalSingleCount  = 0;
my $uniqueSingleCount = 0;
my $dupSingleCount    = 0;
my $totalPairedCount  = 0;
my $uniquePairedCount = 0;
my $dupPairedCount    = 0;
my $untagCount        = 0;

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
}

if ($untagCount) {
	printf " %12s (%.1f%%) alignments without UMI tags retained\n", $untagCount, 
		($untagCount/($totalSingleCount + $totalPairedCount)) * 100;
}


#### End













######## Subroutines ###########


sub deduplicate_multithread {
	
	# set up manager
	my $pm = Parallel::ForkManager->new($cpu);
	$pm->run_on_finish( sub {
		my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data) = @_;
		
		# add child counts to global values
		$totalSingleCount  += $data->{totalSingleCount};
		$uniqueSingleCount += $data->{uniqueSingleCount};
		$dupSingleCount    += $data->{dupSingleCount};
		$totalPairedCount  += $data->{totalPairedCount};
		$uniquePairedCount += $data->{uniquePairedCount};
		$dupPairedCount    += $data->{dupPairedCount};
		$untagCount        += $data->{untagCount};
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
			uniqueSingleCount   => 0,
			dupSingleCount      => 0,
			totalPairedCount    => 0,
			uniquePairedCount   => 0,
			dupPairedCount      => 0,
			untagCount          => 0,
			outbam              => $tempbam,
			single_reads        => [],
			paired_reads        => [],
			keepers             => {},
			duplicates          => {},
			calculator          => $Calculator
		};
	
		# walk through the reads on the chromosome
		my $seq_length = $sam->target_len($tid);
		low_level_bam_fetch($sam, $tid, 0, $seq_length, \&callback, $data);
	
		# check to make sure we don't leave something behind
		if (defined $data->{single_reads}->[0] or defined $data->{paired_reads}->[0]) {
			write_reads($data);
		}
		
		# finish and return to parent
		delete $data->{position};
		delete $data->{outbam};
		delete $data->{single_reads};
		delete $data->{paired_reads};
		delete $data->{keepers};
		delete $data->{duplicates};
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
	
	# put the reads into a hash based on tag identifer
	my %tag2reads;
	my @untagged;
	while (my $a = shift @{ $data->{single_reads} }) {
		
		# first check for any existing if supplemental or secondary
		# otherwise collect UMI
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
	
	# must sort the bams by the barcode
	my %tag2reads; # barcode-to-alignment hash
	my @untagged;
	foreach my $a (@$reads) {
		my $name = $a->qname;
			# we will count those weird situations here where a reverse read came first
			# and was processed first, leaving a forward read name in the keeper hash
			# where it might not get counted normally - so we do so here 
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
		else {
			my $tag = &$get_umi($a);
			if (defined $tag) {
				$tag2reads{$tag} ||= [];
				push @{ $tag2reads{$tag} }, $a;
			}
			else {
				push @untagged, $a;
			}
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
	
	# write untagged reads
	if (@untagged) {
		$data->{untagCount} += scalar(@untagged);
		foreach my $a (@untagged) {
			$data->{keepers}{$a->qname} += 1; # remember to write pair
			&$write_alignment($data->{outbam}, $a);
		}
	}
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
	my @list = keys %$tag2a;
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
			}
			# make new list for next cycle - but skip the current tag
			@list = ();
			foreach (keys %tag2a) {
				push @list, $_ if $_ ne $first;
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



