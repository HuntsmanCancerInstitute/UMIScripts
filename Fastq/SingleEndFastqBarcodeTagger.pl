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

########################################



my $help = <<END;

A script to append a barcode read sequence from a second fastq 
file to the read name of a single-end fastq file. The barcode 
is a unique molecular index (UMI) to properly identify PCR 
duplicates.

Compare to the USeq tool FastqBarcodeTagger for paired-end reads.

Use the bam_umi_dedup.pl in this repository for de-duplication based 
on the UMI barcode.

Usage: SingleEndFastqBarcodeTagger -f read1.fastq.gz -b bc.fastq.gz -o read1bc.fastq.gz

Options:
    -f --fastq <file>     Fastq read, may be gzipped
    -b --barcode <file>   Barcode fastq read, may be gzipped
    -o --out <name>       Specify a new basename for the output fastq files.
                            Use 'stdout' to pipe to standout out.

END

unless (@ARGV) {
	print $help;
	exit;
}

my $file1;
my $file2;
my $out;
my $print_help;

# options
GetOptions( 
	'f|fastq=s'    => \$file1, # first fastq file
	'b|barcode=s'  => \$file2, # second fastq file, with barcode
	'o|out=s'      => \$out, # output file base name
	'h|help'       => \$print_help, # output directory
) or die "bad options!\n";

if ($print_help) {
	print $help;
	exit 0;
}

# check input fastq files
if (not $file1 and not $file2) {
	die "need to provide two fastq files!\n";
}




### open filehandles
my ($infh1, $infh2, $outfh);

if ($file1 =~ /gz$/i) {
	$infh1 = IO::File->new("gzip -dc $file1 |") or 
		die "unable to open $file1! $!\n";
}
else {
	$infh1 = IO::File->new("$file1") or 
		die "unable to open $file1! $!\n";
}
if ($file2 =~ /gz$/i) {
	$infh2 = IO::File->new("gzip -dc $file2 |") or 
		die "unable to open $file2! $!\n";
}
else {
	$infh2 = IO::File->new("$file2") or 
		die "unable to open $file2! $!\n";
}
if ($out =~ /^stdout$/i) {
	# this is overkill, but just to make things consistent, 
	# open an IO::Handle to STDOUT
	$outfh = IO::Handle->new;
	$outfh->fdopen(fileno(STDOUT), 'w');
}
elsif ($out =~ /gz$/i) {
	# check for pigz to make it go faster
	my $gzipper = `which pigz`;
	chomp $gzipper;
	if ($gzipper) {
		$gzipper .= ' -p 3 '; # three threads should do it
	}
	else {
		$gzipper = 'gzip';
	}
	# open file handle
	$outfh = IO::File->new("| $gzipper >$out") or 
		die "unable to open output $out! $!\n";
}
else {
	$outfh = IO::File->new("$out", 'w') or 
		die "unable to open output $out! $!\n";
}


# process the files
my $goodCount = 0;

while (my $header1  = $infh1->getline) {
	my $sequence1 = $infh1->getline or die "malformed file! no sequence line";
	my $spacer1 = $infh1->getline or die "malformed file! no spacer line";
	my $quality1 = $infh1->getline or die "malformed file! no quality line";
	
	# second fastq file stuff
	my $header2  = $infh2->getline;
	my $sequence2 = $infh2->getline or die "malformed file! no sequence line";
	my $spacer2 = $infh2->getline or die "malformed file! no spacer line";
	my $quality2 = $infh2->getline or die "malformed file! no quality line";
	
	# split the first header, I hope that's a space in there!!!!
	my ($name1, $desc1) = split / +/, $header1;
	my ($name2, $desc2) = split / +/, $header2;
	
	# basic header check
	unless ($name1 eq $name2) {
		die "malformed fastq files! sequence IDs from files don't match! compare\n" . 
			" $file1: $header1\n $file2: $header2\n";
	}
	
	# write out
	chomp $sequence2;
	$outfh->printf("%s:%s %s%s%s%s", $name1, $sequence2, $desc2, $sequence1, 
		$spacer1, $quality1);
	$goodCount++;
}

# close file handles
$infh1->close;
$infh2->close;
$outfh->close;

# finish
if ($out =~ /^stdout$/i) {
	warn sprintf " %d reads were processed\n", $goodCount;
}
else {
	printf " %d reads were processed\n", $goodCount;
	print " wrote $out\n";
}






