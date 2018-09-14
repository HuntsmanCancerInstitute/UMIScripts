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

####### Documentation
unless (@ARGV) {
	print <<END;

A script to pull out Unique Molecular Index (UMI) barcodes embedded at 
the beginning of the sequence read, as opposed to a seperate sequence read. 
The UMI barcode may or may not be followed by a short fixed sequence. 

The UMI code is appended to the read name in the format of ":UMI:NNNNN".

This script will work with single-end fastq files, paired-end fastq files 
where both have an independent UMI at the beginning, and paired-end fastq 
files where only one file has an UMI at the beginning. Fastq files are 
expected to be gzip compressed, and output files will be gzip compressed.

The sequence quality of the UMI is not checked nor retained. Mismatches in 
any fixed sequence are not tolerated.

Usage: 

  embedded_UMI_barcode_tagger.pl --input file1.fastq.gz
  
  embedded_UMI_barcode_tagger.pl --input file1.fastq.gz --umipair file2.fastq.gz \\
    --out file1.bc.fastq.gz --outpair file2.bc.fastq.gz

Options:

  -i --input   <file>     Fastq file1 with UMI
  -p --pair    <file>     Paired fastq file2 without UMI
  -u --umipair <file>     Paired fastq file2 with UMI
  -o --out     <file>     Output file1. Optional
  -s --second  <file>     Output file2. Optional
  -l --length  <integer>  Length of the randome UMI at beginning
  -f --fixed   <text>     Optional fixed sequence after UMI and before insert 
  
END
	exit;
}


#### Parameters
my $infile; 
my $inpair; 
my $umipair;  
my $outfile; 
my $outpair; 
my $fixed = ''; 
my $umi_length = 0;

GetOptions( 
	'input=s'           => \$infile, # input fastq file with umi
	'pair=s'            => \$inpair, # input paired without umi
	'umipair=s'         => \$umipair, # input paired with umi code
	'out=s'             => \$outfile, # output fastq file 
	'second=s'          => \$outpair, # output paired fastq file
	'fixed=s'           => \$fixed, # fixed portion of anchor sequence
	'length=i'          => \$umi_length, # length of the barcode
) or die "unrecognized options!\n";






#### Check parameters
unless ($infile) {
	die "must provide an input file!\n";
}
unless ($umi_length > 0) {
	die "must provide a length for the UMI barcode!\n";
}

# identify appropriate gzip program
my $zipper = `which pigz`;
chomp($zipper);
if ($zipper) {
	# three threads for parallel gzip output usually matches speed of input
	$zipper .= " -p 3 ";
}
else {
	$zipper = "gzip -c ";
}


# prepare the umi barcode search callback routine
# we're using substr and index instead of regex to speed things up just a little
my $fix_length = length($fixed);
my $barCodeLength = $umi_length + $fix_length;
my $search_callback;
if ($fixed) {
	$search_callback = sub {
		my $seq = shift;
		if (substr($seq, $umi_length, $fix_length) eq $fixed) {
			# we have the correct fixed code
			# pull the UMI sequence
			my $u = substr($seq, 0, $umi_length);
			# check for N's
			if (index($u,'N') == -1) {
				# not found
				return $u;
			}
			else {
				return;
			}
		}
		return;
	};
}
else {
	$search_callback = sub {
		my $seq = shift;
		my $u = substr($seq, 0, $umi_length);
		# check for N's
		if (index($u,'N') == -1) {
			# not found
			return $u;
		}
		else {
			return;
		}
	};
}





####### Process file based on number of inputs

# One Fastq only
if ($infile and not $inpair and not $umipair) {
	print "processing $infile....\n";
	my $goodCount = 0;
	my $badCount  = 0;
	
	# open input file
	if ($infile !~ /\.gz$/i) {
		"why is $infile not gzip compressed!!!!\n";
	}
	my $infh = IO::File->new("gzip -dc $infile |") or 
		die "unable to open $infile! $!";
	
	# check output file name
	unless ($outfile) {
		my ($basename, $path, $extension) = fileparse($infile, qw(.txt.gz .fastq.gz));
		$outfile = $path . $basename . '.bc' . $extension;
	}
	my $outfh = IO::File->new("| $zipper > $outfile") or 
		die "unable to open $outfile! $!";
	
	while (my $header = $infh->getline) {
		chomp $header; 
		# we only chomp the header, we can safely skip the others
		my $sequence = $infh->getline or die "malformed $infile file! no sequence line";
		my $spacer = $infh->getline or die "malformed $infile file! no spacer line";
		my $quality = $infh->getline or die "malformed $infile file! no quality line";
		
		# check for bar code sequence
		my $umi = &$search_callback($sequence);
		if ($umi) {
			my ($name, $desc) = split /\s+/, $header;
			$name .= ":UMI:$umi";
			$header = $name . ' ' . $desc;
			$sequence = substr($sequence, $barCodeLength);
			$quality = substr($quality, $barCodeLength);
			$outfh->print("$header\n$sequence$spacer$quality");
			$goodCount++;
		}
		else {
			$badCount++;
		}
	}
	
	# close and finish
	$infh->close;
	$outfh->close;
	printf " %12s reads had a bar code and were processed\n %12s reads had no bar code\n", 
		$goodCount, $badCount;
	print " wrote $outfile\n";
}

# Paired Fastq, first is barcoded only
elsif ($infile and $inpair and not $umipair) {
	print "processing $infile for UMIs with $inpair ....\n";
	my $goodCount = 0;
	my $badCount  = 0;
	
	# open input file
	if ($infile !~ /\.gz$/i) {
		"why is $infile not gzip compressed!!!!\n";
	}
	if ($inpair !~ /\.gz$/i) {
		"why is $inpair not gzip compressed!!!!\n";
	}
	my $infh1 = IO::File->new("gzip -dc $infile |") or 
		die "unable to open $infile! $!";
	my $infh2 = IO::File->new("gzip -dc $inpair |") or 
		die "unable to open $inpair! $!";
	
	# check output file name
	unless ($outfile) {
		my ($basename, $path, $extension) = fileparse($infile, qw(.txt.gz .fastq.gz));
		$outfile = $path . $basename . '.bc' . $extension;
	}
	unless ($outpair) {
		my ($basename, $path, $extension) = fileparse($inpair, qw(.txt.gz .fastq.gz));
		$outpair = $path . $basename . '.bc' . $extension;
	}
	
	# open output file handles
	my $outfh1 = IO::File->new("| $zipper > $outfile") or 
		die "unable to open $outfile! $!";
	my $outfh2 = IO::File->new("| $zipper > $outpair") or 
		die "unable to open $outpair! $!";
	
	while (my $header1 = $infh1->getline) {
		chomp $header1; 
		# we only chomp the header, we can safely skip the others
		my $sequence1 = $infh1->getline or die "malformed $infile file! no sequence line";
		my $spacer1 = $infh1->getline or die "malformed $infile file! no spacer line";
		my $quality1 = $infh1->getline or die "malformed $infile file! no quality line";
		
		# second file
		my $header2 = $infh2->getline or die "malformed $inpair file! no sequence line";
		my $sequence2 = $infh2->getline or die "malformed $inpair file! no sequence line";
		my $spacer2 = $infh2->getline or die "malformed $inpair file! no spacer line";
		my $quality2 = $infh2->getline or die "malformed $inpair file! no quality line";
		chomp $header2;
		
		# check for bar code sequence
		my $umi = &$search_callback($sequence1);
		if ($umi) {
			my $random = $1 or die "nothing captured!"; # grabbed from the barcode regex
			
			# first file
			my ($name, $desc) = split /\s+/, $header1;
			$name .= ":UMI:$umi";
			$header1 = $name . ' ' . $desc;
			$sequence1 = substr($sequence1, $barCodeLength);
			$quality1 = substr($quality1, $barCodeLength);
			$outfh1->print("$header1\n$sequence1$spacer1$quality1");
			
			# second file
			($name, $desc) = split /\s+/, $header2;
			$name .= ":UMI:$umi";
			$header2 = $name . ' ' . $desc;
			$sequence2 = substr($sequence2, $barCodeLength);
			$quality2 = substr($quality2, $barCodeLength);
			$outfh2->print("$header2\n$sequence2$spacer2$quality2");
			
			$goodCount++;
		}
		else {
			$badCount++;
		}
	}
	
	# close and finish
	$infh1->close;
	$infh2->close;
	$outfh1->close;
	$outfh2->close;
	printf " %12s reads had a bar code and were processed\n %12s reads had no bar code\n", 
		$goodCount, $badCount;
	print " wrote $outfile and $outpair\n";
}

# Paired Fastq, both are barcoded
elsif ($infile and not $inpair and $umipair) {
	print "processing both $infile and $umipair for UMIs....\n";
	my $goodCount = 0;
	my $badCount  = 0;
	
	# open input file
	if ($infile !~ /\.gz$/i) {
		"why is $infile not gzip compressed!!!!\n";
	}
	if ($umipair !~ /\.gz$/i) {
		"why is $umipair not gzip compressed!!!!\n";
	}
	my $infh1 = IO::File->new("gzip -dc $infile |") or 
		die "unable to open $infile! $!";
	my $infh2 = IO::File->new("gzip -dc $umipair |") or 
		die "unable to open $umipair! $!";
	
	# check output file name
	unless ($outfile) {
		my ($basename, $path, $extension) = fileparse($infile, qw(.txt.gz .fastq.gz));
		$outfile = $path . $basename . '.bc' . $extension;
	}
	unless ($outpair) {
		my ($basename, $path, $extension) = fileparse($umipair, qw(.txt.gz .fastq.gz));
		$outpair = $path . $basename . '.bc' . $extension;
	}
	
	# open output file handles
	my $outfh1 = IO::File->new("| $zipper > $outfile") or 
		die "unable to open $outfile! $!";
	my $outfh2 = IO::File->new("| $zipper > $outpair") or 
		die "unable to open $outpair! $!";
	
	while (my $header1 = $infh1->getline) {
		chomp $header1; 
		# we only chomp the header, we can safely skip the others
		my $sequence1 = $infh1->getline or die "malformed $infile file! no sequence line";
		my $spacer1 = $infh1->getline or die "malformed $infile file! no spacer line";
		my $quality1 = $infh1->getline or die "malformed $infile file! no quality line";
		
		# second file
		my $header2 = $infh2->getline or die "malformed $umipair file! no sequence line";
		my $sequence2 = $infh2->getline or die "malformed $umipair file! no sequence line";
		my $spacer2 = $infh2->getline or die "malformed $umipair file! no spacer line";
		my $quality2 = $infh2->getline or die "malformed $umipair file! no quality line";
		chomp $header2;
		
		# check for bar code sequence
		my $umi1 = &$search_callback($sequence1);
		my $umi2 = &$search_callback($sequence2);

		# entire barcode
		if ($umi1 and $umi2) {
			my $umi = $umi1 . $umi2;
		
			# first file
			my ($name, $desc) = split /\s+/, $header1;
			$name .= ":UMI:$umi";
			$header1 = $name . ' ' . $desc;
			$sequence1 = substr($sequence1, $barCodeLength);
			$quality1 = substr($quality1, $barCodeLength);
			$outfh1->print("$header1\n$sequence1$spacer1$quality1");
			
			# second file
			($name, $desc) = split /\s+/, $header2;
			$name .= ":UMI:$umi";
			$header2 = $name . ' ' . $desc;
			$sequence2 = substr($sequence2, $barCodeLength);
			$quality2 = substr($quality2, $barCodeLength);
			$outfh2->print("$header2\n$sequence2$spacer2$quality2");
			
			$goodCount++;
		}
		else {
			# what happens if only one read had a barcode? do we still use it?
			# this only matters if they had used fixed sequences
			$badCount++;
		}
	}
	
	# close and finish
	$infh1->close;
	$infh2->close;
	$outfh1->close;
	$outfh2->close;
	printf " %12s reads had a bar code and were processed\n %12s reads had no bar code\n", 
		$goodCount, $badCount;
	print " wrote $outfile and $outpair\n";
}




