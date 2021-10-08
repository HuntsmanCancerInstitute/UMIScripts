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
use Bio::ToolBox::db_helper 1.50 qw(
	open_db_connection 
	low_level_bam_fetch
	$BAM_ADAPTER
);
# this can import either Bio::DB::Sam or Bio::DB::HTS depending on availability
# this script is mostly bam adapter agnostic

my $VERSION = 1.9;
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

#### Inputs
my $parallel;
eval {
	# check for parallel support
	require Parallel::ForkManager;
	$parallel = 1;
};
my $infile;
my $outfile;
my $markdups;
my $paired;
my $keep_secondary = 0;
my $cpu = $parallel ? 4 : 1;
my $no_sam;
my $help;
my $sam_app = sprintf("%s", which 'samtools');

my $description = <<END;

A script to remove duplicates based on a Unique Molecular Index (UMI) code.
Alignments that match the same coordinate are checked for the random UMI code 
in the read name; Reads with the same UMI are sorted and selected for one 
to be retained. 

UMI codes are extracted from the end of the read name as ":NNNN", where 
NNNN is a sequence of indeterminate length comprised of [A,T,G,C,N]. 
See the other scripts in the UMIScripts project for pre-processing the Fastq 
files prior to alignments, which will append the random UMI code to the read 
name. This should also work with UMIs extracted with the new Illumina 
BCL2Fastq v2 software, which are appended to the read name.

For single-end alignments, selection criteria include mapping quality and the 
sum of read base qualities. For paired-end alignments, selection criteria include 
the mapping qualities of the forward and possibly reverse mate alignments (using tag 
'MQ' if present) and the base qualities of the forward and possibly reverse mate 
alignments (using the samtools fixmate tag 'ms'). 

Duplicate alignments by default are discarded, or they may be optionally 
marked as duplicate (bit flag 0x400) and retained.

Bam files must be sorted by coordinate. Bam files will be indexed as necessary.

Unaligned (flag 0x4) and supplementary (flag 0x800) reads are silently 
discarded. Secondary alignments may be kept and de-duplicated if specified; 
however, no guarantee is made that the same sequence read (out of multiple 
with the same UMI) will be retained across all primary and secondary alignments.

END

my $usage = <<END;   

VERSION: $VERSION

USAGE:  bam_umi_dedup.pl --in in.bam --out out.bam

OPTIONS:
    -i --in <file>      The input bam file, should be sorted and indexed
    -o --out <file>     The output bam file.
    -p --pe             The bam file consists of paired-end alignments
                          Only properly-paired alignments are retained
    -m --mark           Mark duplicates instead of discarding
    -s --secondary      Keep secondary (multi-hit) reads
    -c --cpu <int>      Specify the number of threads to use ($cpu) 
    --samtools <path>   Path to samtools >= 1.10 ($sam_app)
    -h --help           Display help

END


#### Parse Inputs
unless (@ARGV) {
	print $usage;
	exit 0;
}
GetOptions( 
	'i|in=s'        => \$infile, # the input bam file path
	'o|out=s'       => \$outfile, # name of output file 
	'p|pe!'         => \$paired, # paired-end alignments
	'm|mark!'       => \$markdups, # mark duplicates instead of remove
	's|secondary!'  => \$keep_secondary, # keep secondary alignments 
	'c|cpu=i'       => \$cpu, # number of cpu cores to use
	'samtools=s'    => \$sam_app, # path to samtools
	'bam=s'         => \$BAM_ADAPTER, # specifically set the bam adapter, advanced!
	'nosam!'        => \$no_sam, # avoid using external sam adapter, advanced!
	'h|help!'       => \$help, # print help
) or die " unrecognized option(s)!! please refer to the help documentation\n\n";

if ($help) {
	print $description;
	print $usage;
	exit 0;
}
unless ($infile) {
	die "must provide an input file name!\n";
}
unless ($outfile) {
	die "must provide an output file name!\n";
}
if ($cpu > 1 and not $parallel) {
	warn "Must install Parallel::ForkManager to run multiple forks! Disabling\n";
	$cpu = 1;
}
if (not $sam_app and $cpu > 1) {
	print "A samtools application must be available when running multiple forks! Disabling\n";
	$cpu = 1;
} 
my $sam_version;
if ($sam_app and $cpu > 1) {
	my $sam_help = qx($sam_app 2>&1);
	if ($sam_help =~ /Version: 1\.(\d+) /) {
		$sam_version = $1;
	}
	else {
		print "Unrecognized samtools version! Disabling\n";
		$sam_app = '';
		$cpu = 1;
	}
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
$htext .= sprintf("\@PG\tID:bam_umi_dedup\tVN:%s\tCL:%s", $VERSION, $0);
$htext .= " --pe" if $paired;
$htext .= " --mark" if $markdups;
$htext .= " --secondary" if $keep_secondary;
$htext .= " --samtools $sam_app" if ($cpu > 1 and $sam_app);
$htext .= " --bam $BAM_ADAPTER --in $infile --out $outfile\n";


# which callback to use
my $callback = $paired ? \&pe_callback : \&se_callback;


#initialize counters
my $totalCount = 0;
my $uniqueCount  = 0;
my $duplicateCount   = 0;
my $untagCount = 0;


# deduplicate in single or multi-thread
printf "  - %sing secondary alignments\n", $keep_secondary ? 'Keep' : 'Discard';
printf "  - Discarding single and improperly paired alignments\n" if $paired;
print "  - Discarding supplementary alignments\n  - Discarding unmapped alignments\n";
my $outbam; # I need to return the final output bam to main:: and let the normal 
          # perl exit close it properly, otherwise it crashes hard, and there's 
          # no proper close for the object!!!??????
if ($parallel and $cpu > 1) {
	deduplicate_multithread() or die " Something went wrong with de-duplication!\n";
}
else {
	deduplicate_singlethread() or die " Something went wrong with de-duplication!\n";
}
undef $outbam if $outbam;


# finish up
my $type = $paired ? 'paired-alignments' : 'alignments';
printf "
 File $infile:
 %12s total mapped $type
 %12s (%.1f%%) UMI-unique $type retained
 %12s (%.1f%%) UMI-duplicate $type %s
", 
	$totalCount, $uniqueCount, ($uniqueCount/$totalCount) * 100, $duplicateCount, 
	($duplicateCount/$totalCount) * 100, $markdups ? 'marked' : 'discarded';
if ($untagCount) {
	printf " %12s (%.1f%%) alignments without UMI tags retained\n", $untagCount, 
		($untagCount/$totalCount) * 100;
}


#### End
















# process on chromosomes
sub deduplicate_singlethread {
	
	# open bam new bam file
	$outbam = Bio::ToolBox::db_helper::write_new_bam_file($outfile) or 
		die "unable to open output bam file $outfile! $!";
		# using an unexported subroutine imported as necessary depending on bam availability
	# we can't write updated text to the header - otherwise it crashes
	# $header->text($htext);
	$outbam->header_write($header);
	
	# walk through the file
	for my $tid (0 .. $sam->n_targets - 1) {
		my $seq_length = $sam->target_len($tid);
	
		# prepare callback data structure
		# anonymous hash of 
		my $data = {
			position        => -1,
			totalCount      => 0,
			uniqueCount     => 0,
			duplicateCount  => 0,
			untagged        => 0,
			outbam          => $outbam,
			reads           => [],
			keepers         => {},
			dupkeepers      => {}
		};
	
		# walk through the reads on the chromosome
		low_level_bam_fetch($sam, $tid, 0, $seq_length, $callback, $data);
	
		# check to make sure we don't leave something behind
		if ($data->{reads}->[0]) {
			if ($paired) {
				write_pe_reads($data);
			}
			else {
				write_se_reads($data);
			}
		}
		
		# add up the local counts to global counts
		$totalCount += $data->{totalCount};
		$uniqueCount  += $data->{uniqueCount};
		$duplicateCount   += $data->{duplicateCount};
		$untagCount += $data->{untagged};
	}
	
	# finished
	return 1;
}


sub deduplicate_multithread {
	
	# set up manager
	my $pm = Parallel::ForkManager->new($cpu);
	$pm->run_on_finish( sub {
		my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data) = @_;
		# add child counts to global values
		$totalCount += $data->{totalCount};
		$uniqueCount  += $data->{uniqueCount};
		$duplicateCount   += $data->{duplicateCount};
		$untagCount += $data->{untagged};
	});
	
	# prepare targets and names
	my @targets = (0 .. $sam->n_targets - 1);
	my $tempfile = $outfile;
	$tempfile =~ s/\.bam$//i;
	my @targetfiles = map {$tempfile . ".temp.$_.bam"} @targets;
	
	# walk through the file in parallel, one fork per chromosome
	for my $tid (@targets) {
		$pm->start and next;
		
		### in child
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
		
		# prepare callback data structure
		my $data = {
			position        => -1,
			totalCount      => 0,
			uniqueCount     => 0,
			duplicateCount  => 0,
			untagged        => 0,
			outbam          => $tempbam,
			reads           => [],
			keepers         => {},
			dupkeepers      => {}
		};
	
		# walk through the reads on the chromosome
		my $seq_length = $sam->target_len($tid);
		low_level_bam_fetch($sam, $tid, 0, $seq_length, $callback, $data);
	
		# check to make sure we don't leave something behind
		if ($data->{reads}->[0]) {
			if ($paired) {
				write_pe_reads($data);
			}
			else {
				write_se_reads($data);
			}
		}
		
		# finish and return to parent
		delete $data->{outbam};
		delete $data->{keepers};
		delete $data->{reads};
		delete $data->{position};
		undef $tempbam;
		$pm->finish(0, $data); 
	}
	$pm->wait_all_children;

	# attempt to use external samtools to merge the bam files
	# this should be faster than going through Perl and bam adapters
	if ($sam_app and not $no_sam) {
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



### alignment callback
sub se_callback {
	my ($a, $data) = @_;
	my $flag = $a->flag;
	return if ($flag & 0x100 and not $keep_secondary); # secondary alignment
	return if $flag & 0x800; # supplementary hit
	# we ignore duplicate marks, assume it's not marked anyway....
	$data->{totalCount} += 1;
	
	# check position 
	my $pos = $a->pos;
	if ($pos == $data->{position}) {
		push @{ $data->{reads} }, $a;
	}
	else {
		write_se_reads($data);
		$data->{position} = $pos;
		push @{ $data->{reads} }, $a;
	}
}

### paired-end alignment callback for writing
sub pe_callback {
	my ($a, $data) = @_;
	
	## check the alignment
	return unless $a->proper_pair; # consider only proper pair alignments
	return unless $a->tid == $a->mtid; # we can't handle cross-chromosome alignments
	my $flag = $a->flag;
	return if ($flag & 0x100 and not $keep_secondary); # secondary alignment
	return if ($flag & 0x800); # supplementary hit
	# we ignore duplicate marks, assume it's not marked anyway....
	
	# check position 
	my $pos = $a->pos;
	if ($pos == $data->{position}) {
		push @{ $data->{reads} }, $a;
	}
	else {
		# write out the reads
		write_pe_reads($data) if defined $data->{reads}->[0];
		
		# reset
		$data->{reads} = [];
		$data->{position} = $pos;
		push @{ $data->{reads} }, $a;
	}
}


### write single-end alignments 
sub write_se_reads {
	my $data = shift;
	
	# put the reads into a hash based on tag identifer
	my %tag2reads;
	my @untagged;
	foreach my $a (@{ $data->{reads} }) {
		my $tag;
		if ($a->qname =~ /:([AGTCN]+)$/) {
			$tag = $1;
		}
		else {
			push @untagged, $a;
			next;
		}
		$tag2reads{$tag} ||= [];
		push @{ $tag2reads{$tag} }, $a;
	}
	
	# write out the unique reads
	foreach my $tag (keys %tag2reads) {
		if (scalar @{ $tag2reads{$tag} } == 1) {
			# good, only 1 read, we can write it
			&$write_alignment($data->{outbam}, $tag2reads{$tag}->[0] );
			$data->{uniqueCount}++; # accepted read 
		}
		else {
			# damn, more than 1 read, we gotta sort
			# sorting first on mapping quality, then sum of base qualities
			my @sorted = 
				map {$_->[0]} 
				sort { ($a->[1] <=> $b->[1]) or ($a->[2] <=> $b->[2]) } 
				map { [$_, $_->qual, sum($_->qscore)] } 
				@{ $tag2reads{$tag} };
			&$write_alignment($data->{outbam}, shift @sorted );
			$data->{uniqueCount}++; # accepted read 
			$data->{duplicateCount} += scalar(@sorted); # bad reads
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
		$data->{untagged} += scalar(@untagged);
		foreach my $a (@untagged) {
			&$write_alignment($data->{outbam}, $a);
		}
	}
	
	# reset
	$data->{reads} = [];
}

### Identify and write out paired-end alignments 
sub write_pe_reads {
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
		my ($uniq, $dup, $untag) = write_pe_reads_on_strand($data, $f_sizes{$s});
		# record up the counts
		$data->{uniqueCount} += $uniq;
		$data->{duplicateCount} += $dup;
		$data->{untagged} += $untag;
		$data->{totalCount} += ($uniq + $dup + $untag);
	}
	
	# write out reverse alignments
	foreach my $s (sort {$a <=> $b} keys %r_sizes) {
		my ($uniq, $dup, $untag) = write_pe_reads_on_strand($data, $r_sizes{$s});
		# discard these counts - we only count forward strand
		# weird situations where reverse came first still have a second chance at counts
	}
}

sub write_pe_reads_on_strand {
	my ($data, $reads) = @_; 
		# data reference and array reference of alignment reads
	
	# internal counts for alignments, only counting forward reads
	my $uniq_c   = 0;
	my $dup_c    = 0;
	my $untag_c  = 0;
	
	# must sort the bams by the barcode
	my %bc2a; # barcode-to-alignment hash
	my @untagged;
	my $number = scalar @$reads - 1;
	for my $i (0 .. $number) {
		my $a = $reads->[$i];
		my $name = $a->qname;
			# we will count those weird situations here where a reverse read came first
			# and was processed first, leaving a forward read name in the keeper hash
			# where it might not get counted normally - so we do so here 
		if (exists $data->{keepers}{$name}) {
			&$write_alignment($data->{outbam}, $a);
			delete $data->{keepers}{$name};
			$uniq_c++ if not $a->reversed; 
		}
		elsif (exists $data->{dupkeepers}{$name}) {
			mark_alignment($a);
			&$write_alignment($data->{outbam}, $a);
			delete $data->{dupkeepers}{$name};
			$dup_c++ if not $a->reversed; 
		}
		elsif ($name =~ /:([AGTCN]+)$/) {
			$bc2a{$1} ||= [];
			push @{ $bc2a{$1} }, $a;
		}
		else {
			push @untagged, $a;
		}
	}
	
	# now write one of each
	foreach my $bc (keys %bc2a) {
		if (scalar @{$bc2a{$bc}} == 1) {
			# only 1 alignment, excellent, write it
			my $a = $bc2a{$bc}->[0];
			$data->{keepers}{ $a->qname } += 1;
			&$write_alignment($data->{outbam}, $a);
			$uniq_c++ if not $a->reversed;
		}
		else {
			# more than one alignment, must rank them by sum of base qualities
			my @sorted = map {$_->[0]} 
				sort { ($a->[1] <=> $b->[1]) or ($a->[2] <=> $b->[2]) } 
				map { [ 
						$_, 
						sum($_->qual, $_->aux_get('MQ') || 0),
						sum($_->qscore, $_->aux_get('ms') || 0) 
				] } @{ $bc2a{$bc} };
			
			# write the first one
			my $best = shift @sorted;
			$data->{keepers}{ $best->qname } += 1;
			&$write_alignment($data->{outbam}, $best);
			
			# increment counters
			$uniq_c++ if not $best->reversed;
			foreach (@sorted) {
				$dup_c++ if not $_->reversed;
			}
			
			# write marked remaining duplicates
			if ($markdups) {
				foreach my $a (@sorted) {
					$data->{dupkeepers}{$a->qname} += 1;
					mark_alignment($a);
					&$write_alignment($data->{outbam}, $a);
				}
			}
		}
	}
	
	# write untagged reads
	if (@untagged) {
		foreach my $a (@untagged) {
			$data->{keepers}{$a->qname} += 1; # remember to write pair
			&$write_alignment($data->{outbam}, $a);
			$untag_c++ if not $a->reversed; # only forward reads counted
		}
	}
	
	# return the internal counts
	return ($uniq_c, $dup_c, $untag_c);
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
	# mark alignments as a duplicate
	my $f = $_[0]->flag;
	unless ($f & 0x400) {
		$f += 0x400;
		$_[0]->flag($f);
	}
}



