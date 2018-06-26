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
use IO::File;
use IO::Handle;
use Getopt::Long;
use List::Util qw(sum);
my $mismatch_ok = 0;
eval {
	require String::Approx; 
	String::Approx->import('aindex');
	$mismatch_ok = 1;
};

########################################


my $fixed = 'AACTGTAGGCACCATCAAT'; 
my $file1;
my $outfile;
my $failfile;
my $umi_length = 12;
my $mismatch;

my $help = <<END;

A script to strip out the Qiagen UMI from single-end reads 
prepared from a Qiagen small RNA library. Qiagen recommends doing 
a long read to capture the small RNA insert, read through the 
19 bp adapter sequence, and then through the $umi_length bp UMI sequence. 

This script will search for the adapter, optionally allowing for up to 
2 bp mismatch, remove it, and take the remainder of the read as the UMI. 
Depending on the size of the insert and the length of the read, the 
entire UMI may not be captured. A UMI of less than 3 bp is considered 
failed. UMIs with Ns are also considered failed.

The script requires that the entire fixed adapter sequence be present. 
For larger insertions where the read sequence may not reach the end of 
the adapter and to the UMI, the reads are either discarded or written 
to an optional Fastq file. These failed reads can be aligned separately 
if desired; be sure to use adapter trimming on these reads.

The UMI is appended to the name of the read (for lack of a better place 
to store it), and can be used in subsequent processing to remove PCR 
duplicates in combination with alignment information. See the script 
qiaseq_bam_deduplication.pl as the companion application to do this.

Usage: $0 -i input.fastq.gz --fail noUMI.fastq.gz | <aligner>

Options:
	-i | --input <file>  Fastq file, may be gzipped
	-o | --out <file>    Specify output filename for the output fastq files.
	                     Include a .gz extension for compression.
	                     Default is to print to STDOUT for piping
	-f | --fail <file>   Specify optional filename for failed fastq reads 
	                     where the UMI cannot be found.
	-m | --mismatch      Allow for up to 2 mismatches in adapter sequence.
	                     Requires String::Approx to be installed. 
	                     About 3X slower, but gains a few percent extra.
	-a | --adapt <ATGC>  The adapter sequence to remove. Up to two   
	                     mismatches are allowed. Default is '$fixed'.
	--len <integer>      Length of the UMI barcode. Default $umi_length.

END

unless (@ARGV) {
	print $help;
	exit;
}



### Options
GetOptions( 
	'i|input=s'         => \$file1, # input fastq file
	'o|out=s'           => \$outfile, # output fastq file 
	'f|fail=s'          => \$failfile, # failed fastq file
	'a|adapt=s'         => \$fixed, # fixed anchor sequence
	'len=i'             => \$umi_length, # length of the barcode
	'm|mismatch!'       => \$mismatch, # allow mismatches
) or die "bad options!\n";



### Preparation
# check input fastq files
if (not $file1) {
	if (scalar(@ARGV) == 1) {
		$file1 = shift @ARGV;
	}
	else {
		die "need to provide fastq input!\n";
	}
}
if (not -e $file1) {
	die "input file $file1 does not exist!\n"
}

# prepare index callback to search for position of adapter sequence
my $myindex;
if ($mismatch and $mismatch_ok) {
	# we are using mismatches
	$myindex = sub {
		return aindex($fixed, [2, 'D0', 'S2', 'I0'], $_[0]);
	};
}
elsif ($mismatch and not $mismatch_ok) {
	# we want to use mismatches but cannot
	warn "Mismatch option not allowed! Please install String::Approx module\n";
	$myindex = sub {
		return index($_[0], $fixed);
	};
}
else {
	# use complete matching only
	$myindex = sub {
		return index($_[0], $fixed);
	};
}


### open filehandles
my ($infh, $outfh, $failfh);

# input
if ($file1 =~ /gz$/i) {
	$infh = IO::File->new("gzip -dc $file1 |") or 
		die "unable to open $file1! $!\n";
}
else {
	$infh = IO::File->new("$file1") or 
		die "unable to open $file1! $!\n";
}

# output
if ($outfile) {
	if ($outfile =~ /\.gz$/) {
		$outfh = IO::File->new("| gzip >$outfile") or 
			die "unable to open output $outfile! $!\n";
	}
	else {
		$outfh = IO::File->new("$outfile", 'w') or 
			die "unable to open output $outfile! $!\n";
	}
}
else {
	# this is overkill, but just to make things consistent, 
	# open an IO::Handle to STDOUT
	$outfh = IO::Handle->new;
	$outfh->fdopen(fileno(STDOUT), 'w');
}

# fail
if ($failfile) {
	if ($failfile =~ /\.gz$/) {
		$failfh = IO::File->new("| gzip >$failfile") or 
			die "unable to open output $failfile! $!\n";
	}
	else {
		$failfh = IO::File->new("$failfile", 'w') or 
			die "unable to open output $failfile! $!\n";
	}
}


### process the files
# set variables
my $notFoundCount  = 0;
my $nCount = 0;
my $shortCount = 0;
my %umi2length;
my $fix_length = length($fixed);

# loop through input file
while (my $header1  = $infh->getline) {
	my $sequence1 = $infh->getline or die "malformed file! no sequence line";
	my $spacer1   = $infh->getline or die "malformed file! no spacer line";
	my $quality1  = $infh->getline or die "malformed file! no quality line";
	
	# find the position of the fixed adapter using the callback determined above
	my $i = &$myindex($sequence1);
	
	if ($i == -1) {
		# the adapter sequence was not found
		$notFoundCount++;
		if ($failfh) {
			$failfh->print("$header1$sequence1$spacer1$quality1");
		}
	}
	else {
		# we found the adapter sequence
		my $umi = substr($sequence1, $i + $fix_length);
		chomp($umi);
		
		# record the length
		my $len = length($umi);
		$umi2length{$len} += 1;
		
		# check length of umi
		if ($len < 3) {
			# arbitrary decision that we want a umi of at least 3 bp
			# otherwise it's not much help, I guess?
			$shortCount++;
			if ($failfh) {
				$failfh->print("$header1$sequence1$spacer1$quality1");
			}
			next;
		}
		elsif ($len > $umi_length) {
			# umi is too long!
			$umi = substr($umi, 0, $umi_length);
		}
		
		# check contents of umi
		if ($umi =~ /N/i) {
			# don't know what to do with Ns, so fail them
			$nCount++; # count how many
			if ($failfh) {
				$failfh->print("$header1$sequence1$spacer1$quality1");
			}
			next;
		}
		
		# remove barcode
		my $trimmed_seq = substr($sequence1, 0, $i);
		my $trimmed_qual = substr($quality1, 0, $i);
		
		# write output
		my ($name1, $desc1) = split /\s+/, $header1;
		chomp($desc1); # we may not have a description
		$outfh->print("$name1:UMI:$umi $desc1\n$trimmed_seq\n$spacer1$trimmed_qual\n");
	}
}


### finish
$infh->close;
$failfh->close if $failfile;
$outfh->close if $outfile;
my $goodCount = sum(values %umi2length);
my $badCount = $notFoundCount + $nCount + $shortCount;
my $totalCount = $goodCount + $badCount;
if ($totalCount == 0) {
	if ($outfile) {
		print "nothing processed!\n";
	}
	else {
		warn "nothing processed!\n";
	}
	exit;
}

# print stats
if ($outfile) {
	# safe to print to stdout
	print " Processed $file1\n";
	printf " %d reads had a discoverable UMI and were processed\n %d (%.1f%%) reads had no discoverable UMI\n", 
		$goodCount, $badCount, ($badCount / ($goodCount + $badCount)) * 100;
	printf "   %d reads did not have a complete adapter\n", $notFoundCount;
	printf "   %d reads had a UMI < 3 bp\n", $shortCount;
	printf "   %d reads had Ns in the UMI\n", $nCount;
	print " UMI length counts\n";
	foreach my $k (sort {$a <=> $b} keys %umi2length) {
		printf("%6s %d\n", $k, $umi2length{$k});
	}
	print " Wrote $outfile\n";
	print " Wrote $failfile\n" if $failfile;
}
else {
	# print to stderr instead
	warn sprintf " %d reads had a discoverable UMI and were processed\n %d (%.1f%%) reads had no discoverable UMI\n", 
		$goodCount, $badCount, ($badCount / ($goodCount + $badCount)) * 100;
	warn sprintf "   %d reads did not have a complete adapter\n", $notFoundCount;
	warn sprintf "   %d reads had a UMI < 3 bp\n", $shortCount;
	warn sprintf "   %d reads had Ns in the UMI\n", $nCount;
	warn " UMI length counts\n";
	foreach my $k (sort {$a <=> $b} keys %umi2length) {
		warn sprintf("%6s %d\n", $k, $umi2length{$k});
	}
	warn " Wrote $failfile\n" if $failfile;
}
exit;




