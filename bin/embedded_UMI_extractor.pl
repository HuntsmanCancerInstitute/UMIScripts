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
use Getopt::Long qw(:config no_ignore_case);
use File::Basename qw(fileparse);
use String::Approx qw(aindex);
use List::Util qw(min);

my $VERSION = 2;

####### Documentation
unless (@ARGV) {
	print <<END;

A script to pull out Unique Molecular Index (UMI) barcodes embedded at 
the beginning of the sequence read, as opposed to a seperate sequence read. 
The UMI barcode may or may not be followed by a short fixed sequence. 

The UMI code is appended to the read name in the format of ":UMI:NNNNN".

This script will work with single-end fastq files, paired-end fastq files
where both have an independent UMI at the beginning, and paired-end fastq
files where only one file has an UMI at the beginning. If the UMI is only
included in read2 of a paired-end sequence, simply reverse the order and
set read2 as the input with the UMI; for purposes of extracting the UMI
only, the order is not absolutely required.

Fastq files are expected to be gzip compressed, and output files will be 
gzip compressed.

The sequence quality of the UMI is checked for a minimum acceptable base 
quality. However, it is not retained in the read name. The presence of 
any 'N' bases in the UMI are automatically considered a failure. Reads that 
do not contain the given fixed sequence(s) are considered failed. If both 
reads contain UMIs, then both both reads must pass in order to be written 
to output. Failed reads are skipped entirely and not written to output. 
Result counts are printed to standard output.

Version: $VERSION

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
  -l --length  <integer>  Length of UMI at beginning of read1
  -L --length2 <integer>  Length of UMI at beginning of read2 (if different)
  -f --fixed   <text>     Optional fixed sequence after UMI and before insert 
  -F --fixed2  <text>     Optional fixed sequence after UMI and before insert
                            for read2 (if different)
  -q --qual    <integer>  Minimum base quality for UMI sequence (default 20)
  
END
	exit;
}


#### Parameters
my $infile; 
my $inpair; 
my $umipair;  
my $outfile; 
my $outpair; 
my $fixed1 = ''; 
my $fixed2 = '';
my $umi1_length = 0;
my $umi2_length = 0;
my $min_basequal = 20;

GetOptions( 
	'input=s'           => \$infile, # input fastq file with umi
	'pair=s'            => \$inpair, # input paired without umi
	'umipair=s'         => \$umipair, # input paired with umi code
	'out=s'             => \$outfile, # output fastq file 
	'second=s'          => \$outpair, # output paired fastq file
	'f|fixed=s'         => \$fixed1, # fixed portion of anchor sequence
	'F|fixed2=s'        => \$fixed2, # fixed portion of anchor sequence
	'l|length=i'        => \$umi1_length, # length of the barcode
	'L|length2=i'       => \$umi2_length, # length of the barcode
	'q|qual=i'          => \$min_basequal, # minimum base quality for UMI seq
) or die "unrecognized options!\n";






#### Check parameters
unless ($infile) {
	die "must provide an input file!\n";
}
unless ($umi1_length > 0) {
	die "must provide a length for the UMI barcode!\n";
}





#### Prepare the search callbacks

# prepare the umi barcode search callback routine for read1
# we're using substr and index instead of regex to speed things up just a little
my $fix1_length = length($fixed1);
my $barCode1_length = $umi1_length + $fix1_length;
my $search1_callback;
if ($fixed1) {
	$search1_callback = sub {
		my ($seq, $qual) = @_;
		if (substr($seq, $umi1_length, $fix1_length) eq $fixed1) {
			# we have the correct fixed code
			# pull the UMI sequence
			my $u = substr($seq, 0, $umi1_length);
			my $q = substr($qual, 0, $umi1_length);
			return check_quality($u, $q);
		}
		
		# not found, try an imperfect match, allowing only 1 mismatch or deletion
		my $fix_i = aindex($fixed1, [1, 'D1', 'S1', 'I0'], $seq);
		if ($fix_i > 0 and abs(($umi1_length - 1) - $fix_i) <= 2) {
			# we found it and it's within the expected position
			# need to compensate that the umi1_length is essentially a 1-based position
			# but index is returned as 0-based index, plus tolerate a possible deletion
			my $u = substr($seq, 0, $fix_i);
			my $q = substr($qual, 0, $fix_i);
			return check_quality($u, $q);
		}
		
		# still not found
		return;
	};
}
else {
	$search1_callback = sub {
		my ($seq, $qual) = @_;
		my $u = substr($seq, 0, $umi1_length);
		my $q = substr($qual, 0, $umi1_length);
		return check_quality($u, $q);
	};
}

# prepare the umi callback for read2 as necessary
my ($search2_callback, $fix2_length, $barCode2_length, $umi2_length);
if ($umipair) {
	# defaults are the same as read1
	$fixed2 ||= $fixed1;
	$umi2_length ||= $umi1_length;
	
	# calculate new lengths
	$fix2_length = length($fixed2);
	$barCode2_length = $umi2_length + $fix2_length;
	if ($fixed2) {
		$search2_callback = sub {
			my ($seq, $qual) = @_;
			if (substr($seq, $umi2_length, $fix2_length) eq $fixed2) {
				# we have the correct fixed code
				# pull the UMI sequence
				my $u = substr($seq, 0, $umi2_length);
				my $q = substr($qual, 0, $umi2_length);
				return check_quality($u, $q);
			}
		
			# not found, try an imperfect match, allowing only 1 mismatch or deletion
			my $fix_i = aindex($fixed2, [1, 'D1', 'S1', 'I0'], $seq);
			if ($fix_i > 0 and abs(($umi2_length - 1) - $fix_i) <= 2) {
				# we found it and it's within the expected position
				# need to compensate that the umi1_length is essentially a 1-based position
				# but index is returned as 0-based index, plus tolerate a possible deletion
				my $u = substr($seq, 0, $fix_i);
				my $q = substr($qual, 0, $fix_i);
				return check_quality($u, $q);
			}
		
			# still not found
			return;
		};
	}
	else {
		$search2_callback = sub {
			my ($seq, $qual) = @_;
			my $u = substr($seq, 0, $umi2_length);
			my $q = substr($qual, 0, $umi2_length);
			return check_quality($u, $q);
		};
	}
}






####### Process file based on number of inputs

# One Fastq only
if ($infile and not $inpair and not $umipair) {
	print "processing $infile....\n";
	my $goodCount = 0;
	my $badCount  = 0;
	
	# open input file
	if ($infile !~ /\.gz$/i) {
		die "why is $infile not gzip compressed!!!!\n";
	}
	my $infh = IO::File->new("gzip -dc $infile |") or 
		die "unable to open $infile! $!";
	
	# check output file name
	unless ($outfile) {
		my ($basename, $path, $extension) = fileparse($infile, qw(.txt.gz .fastq.gz .fq.gz));
		$outfile = $path . $basename . '.bc' . $extension;
	}
	my $outfh = IO::File->new("| gzip -c > $outfile") or 
		die "unable to open $outfile! $!";
	
	while (my $header = $infh->getline) {
		chomp $header; 
		# we only chomp the header, we can safely skip the others
		my $sequence = $infh->getline or die "malformed $infile file! no sequence line";
		my $spacer = $infh->getline or die "malformed $infile file! no spacer line";
		my $quality = $infh->getline or die "malformed $infile file! no quality line";
		
		# check for bar code sequence
		my $umi = &$search1_callback($sequence, $quality);
		if ($umi) {
			my ($name, $desc) = split /\s+/, $header;
			$name .= ":UMI:$umi";
			$header = $name . ' ' . $desc;
			$sequence = substr($sequence, $barCode1_length);
			$quality = substr($quality, $barCode1_length);
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
	my $total = $goodCount + $badCount;
	printf " %12s reads total\n %12s reads had a bar code and were processed (%.2f%%)\n %12s reads had no bar code (%.2f%%)\n", 
		$total, $goodCount, ($goodCount/$total) * 100, $badCount, ($badCount/$total) * 100;
	print " wrote $outfile\n";
}

# Paired Fastq, first is barcoded only
elsif ($infile and $inpair and not $umipair) {
	print "processing $infile for UMIs with $inpair ....\n";
	my $goodCount = 0;
	my $badCount  = 0;
	
	# open input file
	if ($infile !~ /\.gz$/i) {
		die "why is $infile not gzip compressed!!!!\n";
	}
	if ($inpair !~ /\.gz$/i) {
		die "why is $inpair not gzip compressed!!!!\n";
	}
	my $infh1 = IO::File->new("gzip -dc $infile |") or 
		die "unable to open $infile! $!";
	my $infh2 = IO::File->new("gzip -dc $inpair |") or 
		die "unable to open $inpair! $!";
	
	# check output file name
	unless ($outfile) {
		my ($basename, $path, $extension) = fileparse($infile, qw(.txt.gz .fastq.gz .fq.gz));
		$outfile = $path . $basename . '.bc' . $extension;
	}
	unless ($outpair) {
		my ($basename, $path, $extension) = fileparse($inpair, qw(.txt.gz .fastq.gz .fq.gz));
		$outpair = $path . $basename . '.bc' . $extension;
	}
	
	# open output file handles
	my $outfh1 = IO::File->new("| gzip -c > $outfile") or 
		die "unable to open $outfile! $!";
	my $outfh2 = IO::File->new("| gzip -c > $outpair") or 
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
		my $umi = &$search1_callback($sequence1, $quality1);
		if ($umi) {
			my $random = $1 or die "nothing captured!"; # grabbed from the barcode regex
			
			# first file
			my ($name, $desc) = split /\s+/, $header1;
			$name .= ":UMI:$umi";
			$header1 = $name . ' ' . $desc;
			$sequence1 = substr($sequence1, $barCode1_length);
			$quality1 = substr($quality1, $barCode1_length);
			$outfh1->print("$header1\n$sequence1$spacer1$quality1");
			
			# second file
			($name, $desc) = split /\s+/, $header2;
			$name .= ":UMI:$umi";
			$header2 = $name . ' ' . $desc;
			$sequence2 = substr($sequence2, $barCode1_length);
			$quality2 = substr($quality2, $barCode1_length);
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
	my $total = $goodCount + $badCount;
	printf " %12s read pairs total\n %12s read pairs had a bar code and were processed (%.2f%%)\n %12s read pairs had no bar code (%.2f%%)\n", 
		$total, $goodCount, ($goodCount/$total) * 100, $badCount, ($badCount/$total) * 100;
	print " wrote $outfile and $outpair\n";
}

# Paired Fastq, both are barcoded
elsif ($infile and not $inpair and $umipair) {
	print "processing both $infile and $umipair for UMIs....\n";
	my $goodCount   = 0;
	my $badCount    = 0;
	my $singleCount = 0;
	
	# open input file
	if ($infile !~ /\.gz$/i) {
		die "why is $infile not gzip compressed!!!!\n";
	}
	if ($umipair !~ /\.gz$/i) {
		die "why is $umipair not gzip compressed!!!!\n";
	}
	my $infh1 = IO::File->new("gzip -dc $infile |") or 
		die "unable to open $infile! $!";
	my $infh2 = IO::File->new("gzip -dc $umipair |") or 
		die "unable to open $umipair! $!";
	
	# check output file name
	unless ($outfile) {
		my ($basename, $path, $extension) = fileparse($infile, qw(.txt.gz .fastq.gz .fq.gz));
		$outfile = $path . $basename . '.bc' . $extension;
	}
	unless ($outpair) {
		my ($basename, $path, $extension) = fileparse($umipair, qw(.txt.gz .fastq.gz .fq.gz));
		$outpair = $path . $basename . '.bc' . $extension;
	}
	
	# open output file handles
	my $outfh1 = IO::File->new("| gzip -c > $outfile") or 
		die "unable to open $outfile! $!";
	my $outfh2 = IO::File->new("| gzip -c > $outpair") or 
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
		my $umi1 = &$search1_callback($sequence1, $quality1);
		my $umi2 = &$search2_callback($sequence2, $quality2);

		# entire barcode
		if ($umi1 and $umi2) {
			my $umi = $umi1 . $umi2;
		
			# first file
			my ($name, $desc) = split /\s+/, $header1;
			$name .= ":UMI:$umi";
			$header1 = $name . ' ' . $desc;
			$sequence1 = substr($sequence1, $barCode1_length);
			$quality1 = substr($quality1, $barCode1_length);
			$outfh1->print("$header1\n$sequence1$spacer1$quality1");
			
			# second file
			($name, $desc) = split /\s+/, $header2;
			$name .= ":UMI:$umi";
			$header2 = $name . ' ' . $desc;
			$sequence2 = substr($sequence2, $barCode2_length);
			$quality2 = substr($quality2, $barCode2_length);
			$outfh2->print("$header2\n$sequence2$spacer2$quality2");
			
			$goodCount++;
		}
		elsif ($umi1 or $umi2) {
			# only one read had a good UMI sequence
			$singleCount++;
		}
		else {
			# neither had a good UMI sequence
			$badCount++;
		}
	}
	
	# close and finish
	$infh1->close;
	$infh2->close;
	$outfh1->close;
	$outfh2->close;
	my $total = $goodCount + $singleCount + $badCount;
	printf " %12s read pairs total\n %12s read pairs had a bar code and were processed (%.2f%%)\n %12s read pairs had one failed bar code (%.2f%%)\n %12s read pairs had no bar code in either read (%.2f%%)\n", 
		$total, $goodCount, ($goodCount/$total) * 100, $singleCount, 
		($singleCount/$total) * 100, $badCount, ($badCount/$total) * 100;
	print " wrote $outfile and $outpair\n";
}


sub check_quality {
	my ($seq, $qual) = @_;
	if (index($seq,'N') == -1) {
		# good sequence, no N's
		my @quals = map {ord($_) - 33} split q(), $qual;
		if (min(@quals) >= $min_basequal) {
			# good quality
			return $seq;
		}
	}
	return; # otherwise return nothing indicating a bad sequence
}

