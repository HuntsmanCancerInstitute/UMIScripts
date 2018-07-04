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
my $file1;
my $file2;
my $outbase;
my $outdir;
my $bc_length = 12;
my $interleave;

my $help = <<END;

A script to strip out the Qiagen QiaSeq barcode from the second 
fastq file and incorporate it into the Fastq header. Intended 
as a replacement processor for the USeq tool FastqBarcodeTagger
to properly handle PCR and potential duplicates with unique 
molecular barcode information.

IMPORTANT: Remember to trim the reverse complement of the fixed 
barcode from the first read using cutadapt or similar! Otherwise 
overlapping pairs will yield unaligned bases. Use standard 
adapter for the 2nd read trim.

Usage: $0 --out <basename> --f1 <read_1> --f2 <read_2>
       $0 --interleave <read_1> <read_2>

Options:
        --f1 <file>     First read in Fastq, may be gzipped
        --f2 <file>     Second Fastq read, expected to contain the 
                        barcode. May be gzipped. Both files may also
                        be simply appended to the command line.
        --len <integer> Length of the molecular barcode. Default $bc_length.
        --fix <ATGC>    Sequence of the fixed portion of the barcode. 
                        Up to two edits (mismatch, indel) are allowed. 
                        Default is '$fixed'.
        --out <basename> Specify a new basename for the output fastq files.
                        Default is to reuse the basename from the first 
                        fastq file.
        --dir <directory> Specify a different output directory. Default 
                        is to reuse the input files' directory.
        --interleave    Flag to write interleaved output to standard out

END

unless (@ARGV) {
	print $help;
	exit;
}

# options
GetOptions( 
	'f1=s'       => \$file1, # first fastq file
	'f2=s'       => \$file2, # second fastq file, with barcode
	'fix=s'      => \$fixed, # fixed anchor sequence
	'len=i'      => \$bc_length, # length of the barcode
	'out=s'      => \$outbase, # output file base name
	'dir=s'      => \$outdir, # output directory
	'interleave!'=> \$interleave, # interleave output to stdout
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
unless ($outbase or $outdir or $interleave) {
	die "need to provide either a output basename, output directory, or interleave option!\n";
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


### open filehandles
my ($infh1, $infh2, $outfh1, $outfh2, $out1, $out2);

if ($suffix =~ /gz$/i) {
	$infh1 = IO::File->new("gzip -dc $file1 |") or 
		die "unable to open $file1! $!\n";
	$infh2 = IO::File->new("gzip -dc $file2 |") or 
		die "unable to open $file2! $!\n";
	unless ($interleave) {
		$out1 = $outdir . $outbase . "_1" . $suffix;
		$outfh1 = IO::File->new("| gzip >$out1") or 
			die "unable to open output $out1! $!\n";
		$out2 = $outdir . $outbase . "_2" . $suffix;
		$outfh2 = IO::File->new("| gzip >$out2") or 
			die "unable to open output $out2! $!\n";
	}
}
else {
	$infh1 = IO::File->new("$file1") or 
		die "unable to open $file1! $!\n";
	$infh2 = IO::File->new("$file2") or 
		die "unable to open $file2! $!\n";
	unless ($interleave) {
		$out1 = $outdir . $outbase . "_1" . $suffix;
		$outfh1 = IO::File->new("$out1", 'w') or 
			die "unable to open output $out1! $!\n";
		$out2 = $outdir . $outbase . "_2" . $suffix;
		$outfh2 = IO::File->new("$out2", 'w') or 
			die "unable to open output $out2! $!\n";
	}
}


# process the files
my $goodCount = 0;
my $badCount  = 0;
my $fix_length = length($fixed);
my $bc_fix_length = $bc_length + $fix_length;
my $replacement = 'N' x $bc_fix_length;

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
		# This may not be exact since the random portion could be 11 bp instead of 12
		# resulting in the fixed moving a little bit, which is compensated by the errors  
		# in the match above. This means we take one bp from fix as random. I think 
		# that's ok, since it's pretty rare.
		my $bc = substr($sequence2, 0, $bc_length);
		my $quality_bc = substr($quality2, 0, $bc_length);
		
		# remove barcode
		$sequence2 = substr($sequence2, $bc_fix_length);
		$quality2  = substr($quality2, $bc_fix_length);
		
		# write output
		# we have to write out fastq1 too to keep the files in sync
		my ($name1, $desc1) = split / /, $header1;
		my ($name2, $desc2) = split / /, $header2;
		
		if ($interleave) {
			print "$name1:BMF:$bc$quality_bc/1 $desc1$sequence1$spacer1$quality1$name2:BMF:$bc$quality_bc/2 $desc2$sequence2$spacer2$quality2"
		}
		else {
			$outfh1->print("$name1:BMF:$bc$quality_bc/1 $desc1$sequence1$spacer1$quality1");
			$outfh2->print("$name2:BMF:$bc$quality_bc/2 $desc2$sequence2$spacer2$quality2");
		}
		$goodCount++;
	}
	else {
		$badCount++;
		# print $sequence2;
	}
}

# finish
if ($interleave) {
	warn sprintf " %d reads had a bar code and were processed\n %d (%.2f%%) reads had no detectable bar code\n", 
		$goodCount, $badCount, $badCount / ($goodCount + $badCount);
}
else {
	printf " %d reads had a bar code and were processed\n %d (%.2f%%) reads had no detectable bar code\n", 
		$goodCount, $badCount, $badCount / ($goodCount + $badCount);
	print " wrote $out1\n";
	print " wrote $out2\n";
}






