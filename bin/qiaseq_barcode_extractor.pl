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
use Getopt::Long;
use File::Basename qw(fileparse);
use String::Approx qw(amatch);

########################################


my $fixed = 'ATTGGAGTCCT'; # antisequence is AGGACTCCAAT
my $help = <<END;

A script to strip out the Qiagen QiaSeq barcode from the second 
fastq file and write it out as a third fastq file file. Intended 
as a pre-processor for using the USeq tools FastqBarcodeTagger, 
MatchMates, and Consensus to properly handle PCR and potential 
duplicates with unique molecular barcode information.

Usage: $0 --out <basename> --f1 <read_1> --f2 <read_2>
       $0 --dir <directory> <read_1> <read_2>

Options:
        --f1 <file>     First read in Fastq, may be gzipped
        --f2 <file>     Second Fastq read, expected to contain the 
                        barcode. May be gzipped. Both files may also
                        be simply appended to the command line.
        --len <integer> Length of the molecular barcode. Default 12.
        --fix <ATGC>    Sequence of the fixed portion of the barcode. 
                        Up to two edits (mismatch, indel) are allowed. 
                        Default is '$fixed'.
        --out <basename> Specify a new basename for the output fastq files.
                        Default is to reuse the basename from the first 
                        fastq file.
        --dir <directory> Specify a different output directory. Default 
                        is to reuse the input files' directory.

END

unless (@ARGV) {
	print $help;
	exit;
}

my $file1;
my $file2;
my $outbase;
my $outdir;
my $bc_length = 12;

# options
GetOptions( 
	'f1=s'       => \$file1, # first fastq file
	'f2=s'       => \$file2, # second fastq file, with barcode
	'fix=s'      => \$fixed, # fixed anchor sequence
	'len=i'      => \$bc_length, # length of the barcode
	'out=s'      => \$outbase, # output file base name
	'dir=s'      => \$outdir, # output directory
) or die "bad options!\n";


# check input fastq files
if (not $file1 and not $file2) {
	if (scalar(@ARGV) == 2) {
		$file1 = shift @ARGV;
		$file2 = shift @ARGV;
	}
	else {
		die "need to provide two fastq files!\n";
	}
}

# check basename
unless ($outbase or $outdir) {
	die "need to provide either a different output file basename or an output directory!\n";
}

# set filename paths using first fastq as an example
my ($suffix, $filename, $directories);
($filename, $directories, $suffix) = fileparse($file1, 
	qw(.txt .fastq .txt.gz .fastq.gz .fq .fq.gz));
unless ($outbase) {
	$outbase = $filename;
	$outbase =~ s/_\d$//; # strip the read number suffix
}
unless ($outdir) {
	$outdir = $directories;
}


# directory
if ($outdir) {
	unless (-e $outdir) {
		mkdir $outdir or die "unable to make directory $outdir! $!\n";
	}
}

# # summary
# printf " processing:\n  %s\n  %s\n Output base name: %s\n Output directory: %s\n", 
# 	$file1, $file2, $outbase, $outdir;


### open filehandles
my ($infh1, $infh2, $outfh1, $outfh2, $outfh3, $out1, $out2, $out3);

if ($suffix =~ /gz$/i) {
	$infh1 = IO::File->new("gzip -dc $file1 |") or 
		die "unable to open $file1! $!\n";
	$out1 = $outdir . $outbase . "_1" . $suffix;
	$outfh1 = IO::File->new("| gzip >$out1") or 
		die "unable to open output $out1! $!\n";
	$infh2 = IO::File->new("gzip -dc $file2 |") or 
		die "unable to open $file2! $!\n";
	$out2 = $outdir . $outbase . "_2" . $suffix;
	$outfh2 = IO::File->new("| gzip >$out2") or 
		die "unable to open output $out2! $!\n";
	$out3 = $outdir . $outbase . "_bc" . $suffix;
	$outfh3 = IO::File->new("| gzip >$out3") or 
		die "unable to open output $out3! $!\n";
}
else {
	$infh1 = IO::File->new("$file1") or 
		die "unable to open $file1! $!\n";
	$out1 = $outdir . $outbase . "_1" . $suffix;
	$outfh1 = IO::File->new("$out1", 'w') or 
		die "unable to open output $out1! $!\n";
	$infh2 = IO::File->new("$file2") or 
		die "unable to open $file2! $!\n";
	$out2 = $outdir . $outbase . "_2" . $suffix;
	$outfh2 = IO::File->new("$out2", 'w') or 
		die "unable to open output $out2! $!\n";
	$out3 = $outdir . $outbase . "_bc" . $suffix;
	$outfh3 = IO::File->new("$out3", 'w') or 
		die "unable to open output $out3! $!\n";
}


# process the files
my $goodCount = 0;
my $badCount  = 0;
my $fix_length = length($fixed);

while (my $header1  = $infh1->getline) {
	my $sequence1 = $infh1->getline or die "malformed file! no sequence line";
	my $spacer1 = $infh1->getline or die "malformed file! no spacer line";
	my $quality1 = $infh1->getline or die "malformed file! no quality line";
	
	# second fastq file stuff
	my $header2  = $infh2->getline;
	my $sequence2 = $infh2->getline or die "malformed file! no sequence line";
	my $spacer2 = $infh2->getline or die "malformed file! no spacer line";
	my $quality2 = $infh2->getline or die "malformed file! no quality line";
	
	# basic header check
	unless (amatch($header1, [3], $header2)) {
		die "malformed fastq files! headers from files don't match! compare\n" . 
			" $file1: $header1\n $file2: $header2\n";
	}
	
	# write as appropriate
	if ( amatch($fixed, [2, 'D2', 'S2', 'I0'], substr($sequence2, $bc_length, $fix_length)) ) {
		# get barcode
		my $bc = substr($sequence2, 0, $bc_length);
		my $quality_bc = substr($quality2, 0, $bc_length);
		
		# remove barcode
		$sequence2 = substr($sequence2, $bc_length);
		$quality2  = substr($quality2, $bc_length);
		
		# write output
		# we have to write out fastq1 too to keep the files in sync
		$outfh1->print("$header1$sequence1$spacer1$quality1");
		$outfh2->print("$header2$sequence2$spacer2$quality2");
		$outfh3->print("$header2$bc\n$spacer2$quality_bc\n");
		$goodCount++;
	}
	else {
		$badCount++;
	}
}

# finish
printf " %d reads had a bar code and were processed\n %d (%.2f%%) reads had no detectable bar code\n", 
	$goodCount, $badCount, $badCount / ($goodCount + $badCount);
print " wrote $out1\n";
print " wrote $out2\n";
print " wrote $out3\n";







