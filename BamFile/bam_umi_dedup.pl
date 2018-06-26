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
my $parallel;
eval {
	# check for parallel support
	require Parallel::ForkManager;
	$parallel = 1;
};

my $VERSION = 1.3;
# version 1.0 - initial version
# version 1.1 - added option to write the duplicates to second bam file
# version 1.2 - add marking and parallel processing, make compatible with 
#               qiagen_umi_extractor script
# version 1.3 - remove the UMI in the UMI regex of the read name to work 
#               with Illumina UMI marked fastq files, allow N in UMI, 
#               write all reads without UMI codes instead of dying. 

unless (@ARGV) {
	print <<END;

A script to remove duplicates based on a Unique Molecular Index (UMI) code.
Alignments that match the same coordinate are checked for the random UMI code 
in the read name; Reads with the same UMI are sorted and selected for one 
to be retained. Selection criteria include mapping quality and the sum of read 
base qualities. Retained reads are written to the output bam file. 

See the scripts ChIPNexus_fastq_barcode_processer.pl or qiagen_umi_extractor.pl 
for pre-processing the Fastq files prior to alignments, which will append the 
random UMI code to the read name. This should also work with UMIs extracted 
with the new Illumina BCL2Fastq v2 software, which are appended to the read 
name.

Bam files must be sorted by coordinate. Bam files will be indexed as necessary.
Currently only single-end alignments are accepted.

Duplicate alignments by default are discarded, or they may be optionally 
marked as duplicate (bit flag 0x400) and retained.

Unaligned reads are silently discarded.

USAGE:  bam_umi_dedup.pl --in in.bam --out out.bam
       
OPTIONS:
    -i | --in <file>    The input bam file, should be sorted and indexed
    -o | --out <file>   The output bam file.
    -m | --mark         Write all alignments to output and mark duplicates 
                        with flag bit 0x400.
    -c | --cpu <int>    Specify the number of threads to use (4) 

END
	exit;
}


#### Inputs
my $infile;
my $outfile;
my $markdups;
my $cpu;
my @program_options = @ARGV;
GetOptions( 
	'in=s'       => \$infile, # the input bam file path
	'out=s'      => \$outfile, # name of output file 
	'mark!'      => \$markdups, # mark duplicates instead of remove
	'cpu=i'      => \$cpu, # number of cpu cores to use
	'bam=s'      => \$BAM_ADAPTER, # specifically set the bam adapter, advanced!
) or die " unrecognized option(s)!! please refer to the help documentation\n\n";

if ($parallel and not defined $cpu) {
	$cpu = 4;
}
if ($cpu and not $parallel) {
	warn "must install Parallel::ForkManager to multi-thread!\n";
	$cpu = 1;
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
# my $htext = $header->text;
# $htext .= sprintf("\@PG\tID:bam_umi_dedup\tVN:%s\tCL:%s\n", 
# 	$VERSION, join(' ', $0, @program_options));
# $header->text($htext);


#initialize counters
my $totalCount = 0;
my $goodCount  = 0;
my $badCount   = 0;
my $untagCount = 0;


# deduplicate in single or multi-thread
print " Removing duplicates and writing new bam file";
my $ob; # I need to return the final output bam to main:: and let the normal 
          # perl exit close it properly, otherwise it crashes hard, and there's 
          # no proper close for the object!!!??????
if ($parallel and $cpu > 1) {
	print " in $cpu threads....\n";
	$ob = deduplicate_multithread();
}
else {
	print "....\n";
	$ob = deduplicate_singlethread();
}


# finish up
printf "
 File $infile:
 %12s total mapped alignments
 %12s (%.1f%%) UMI-unique alignments retained
 %12s (%.1f%%) UMI-duplicate alignments %s
", $totalCount, $goodCount, ($goodCount/$totalCount) * 100, $badCount, 
($badCount/$totalCount) * 100, $markdups ? 'marked' : 'discarded';
if ($untagCount) {
	printf " %12s (%.1f%%) alignments without UMI tags retained\n", $untagCount, 
		($untagCount/$totalCount) * 100;
}


# end


# process on chromosomes
sub deduplicate_singlethread {
	
	# open bam new bam file
	my $outbam = Bio::ToolBox::db_helper::write_new_bam_file($outfile) or 
		die "unable to open output bam file $outfile! $!";
		# using an unexported subroutine imported as necessary depending on bam availability
# 	$header->text($htext);
	$outbam->header_write($header);
	
	# walk through the file
	for my $tid (0 .. $sam->n_targets - 1) {
		my $seq_length = $sam->target_len($tid);
	
		# prepare callback data structure
		# anonymous hash of 
		my $data = {
			position   => 0,
			totalCount => 0,
			goodCount  => 0,
			badCount   => 0,
			untagged   => 0,
			outbam     => $outbam,
			reads      => [],
		};
	
		# walk through the reads on the chromosome
		low_level_bam_fetch($sam, $tid, 0, $seq_length, \&callback, $data);
	
		# check to make sure we don't leave something behind
		if ($data->{reads}->[0]) {
			write_reads($data);
		}
		
		# add up the local counts to global counts
		$totalCount += $data->{totalCount};
		$goodCount  += $data->{goodCount};
		$badCount   += $data->{badCount};
		$untagCount += $data->{untagged};
	}
	
	# finished
	return $outbam;
}


sub deduplicate_multithread {
	
	# set up manager
	my $pm = Parallel::ForkManager->new($cpu);
	$pm->run_on_finish( sub {
		my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data) = @_;
		# add child counts to global values
		$totalCount += $data->{totalCount};
		$goodCount  += $data->{goodCount};
		$badCount   += $data->{badCount};
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
		
		# prepare a temporary bam file
		my $tf = $targetfiles[$tid];
		my $tempbam = Bio::ToolBox::db_helper::write_new_bam_file($tf) or 
			die "unable to open output bam file $tf! $!";
		$tempbam->header_write($header);
		
		# prepare callback data structure
		my $data = {
			position   => 0,
			totalCount => 0,
			goodCount  => 0,
			badCount   => 0,
			untagged   => 0,
			outbam     => $tempbam,
			reads      => [],
		};
	
		# walk through the reads on the chromosome
		my $seq_length = $sam->target_len($tid);
		low_level_bam_fetch($sam, $tid, 0, $seq_length, \&callback, $data);
	
		# check to make sure we don't leave something behind
		if ($data->{reads}->[0]) {
			write_reads($data);
		}
		
		# finish and return to parent
		delete $data->{outbam};
		delete $data->{reads};
		delete $data->{position};
		undef $tempbam;
		$pm->finish(0, $data); 
	}
	$pm->wait_all_children;

	# open final bam new bam file
	my $outbam = Bio::ToolBox::db_helper::write_new_bam_file($outfile) or 
		die "unable to open output bam file $outfile! $!";
# 	$header->text($htext);
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
	return $outbam;
}



### alignment callback
sub callback {
	my ($a, $data) = @_;
	$data->{totalCount} += 1;
	
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
			$data->{goodCount}++; # accepted read 
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
			$data->{goodCount}++; # accepted read 
			$data->{badCount} += scalar(@sorted); # bad reads
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



