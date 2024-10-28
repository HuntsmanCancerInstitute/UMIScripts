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
use English qw(-no_match_vars);
use IO::File;
use Getopt::Long;
use String::Approx qw(aindex);

########################################
our $VERSION = 0.2;

my $fixed = 'ATTGATGGTGCCTACAGTT'; # read2 adapter
my $file1;
my $file2;
my $outfile;
my $suspectfile;
my $failfile;
my $umi_length = 12;
my $print_help;

my $help = <<END;

A script to extract the Unique Molecular Index (UMI) from the second 
read of a paired-end Qiagen small RNA library. The structure of 
a typical small RNA paired-end read is as follows:

R1  tagcttatcagactgatgttga AACTGTAGGCACCATCAAT NNNNNNNNN
    <--miRNA-------------> <--adapter--------> <--UMI-->

R2  NNNNNNNNNNNN ATTGATGGTGCCTACAGTT tcaacatcagtctgataag
    <--UMI-----> <--adapter--------> <--miRNA---------->

Since the small RNA insert may be of variable length, the full UMI 
sequence ($umi_length) is not reliably present in R1, hence extracting it from 
R2 is necessary. The Qiagen adapter is identified in R2 (allowing for 
mismatches as necessary), the preceding UMI is extracted, and the UMI 
appended to the read name. The Qiagen adapter is also searched for and 
removed from R1 if present, allowing for mismatches and deletions as 
necessary. 

Only one Fastq file, R1, is written out. Reads with adapter sequence 
in both R2 and R1 is always written to output. Reads with the adapter 
identified only in R2 is written to output or optionally to a second 
file (option --suspect). Reads where the adapter sequence is not found 
in either R2 or R1 may only be written to an optional file (option --fail). 
No UMI code is appended to failed reads.

Sequencing primers are not searched for, but normally should not be 
necessary except for long insertions (suspect reads).

The UMI is appended to the read name as ":NNNNNNNNNNNN". After 
alignment, the bam file may be de-duplicated with bam_umi_dedup.pl.


Usage: 
   smallRNA_pe_umi_extractor.pl --out <basename> --f1 <read_1> --f2 <read_2>
   smallRNA_pe_umi_extractor.pl <read_1> <read_2>

Options:
    -1 |--f1 <file>      First read in Fastq, may be gzipped
    -2 |--f2 <file>      Second Fastq read, expected to contain the 
                           barcode. May be gzipped. Both files may also
                           be simply appended to the command line.
    -o |--out <file>     Specify output fastq filename. Default is STDOUT.
                           GZip compression is fully supported.
    -s |--suspect <file> Specify optional filename for suspect fastq reads
                           where the adapter is only found in Read2.
                           Default is to include in primary output.
    -f |--fail <file>    Specify optional filename for failed fastq reads 
                           where the adapter cannot be found in either read.
                           Default is to not include in primary output.
    -a |--adapt <ATGC>   The Qiagen adapter sequence to look for in Read2.   
                           Default is '$fixed'.
    -l |--len <integer>  Length of the Unique Molecular Index. Default $umi_length.
    -h |--help           This help
END

unless (@ARGV) {
	print $help;
	exit;
}

# options
GetOptions( 
	'1|f1=s'            => \$file1, # first fastq file
	'2|f2=s'            => \$file2, # second fastq file, with barcode
	'o|out=s'           => \$outfile, # output fastq file 
	's|suspect=s'       => \$suspectfile, # suspect fastq file
	'f|fail=s'          => \$failfile, # failed fastq file
	'a|adapt=s'         => \$fixed, # fixed anchor sequence
	'l|len=i'           => \$umi_length, # length of the barcode
	'h|help!'           => \$print_help, # print help
) or die "bad or unknown options! use --help for documentation\n";

if ($print_help) {
	print $help;
	exit;
}


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






### open filehandles
my ($infh1, $infh2, $outfh, $suspectfh, $failfh);

# input1
if ($file1 =~ /gz$/i) {
	$infh1 = IO::File->new("gzip -dc $file1 |") or 
		die "unable to open $file1! $OS_ERROR\n";
}
else {
	$infh1 = IO::File->new("$file1") or 
		die "unable to open $file1! $OS_ERROR\n";
}

# input2
if ($file2 =~ /gz$/i) {
	$infh2 = IO::File->new("gzip -dc $file2 |") or 
		die "unable to open $file2! $OS_ERROR\n";
}
else {
	$infh2 = IO::File->new("$file2") or 
		die "unable to open $file2! $OS_ERROR\n";
}

# output
if ($outfile) {
	if ($outfile =~ /\.gz$/) {
		$outfh = IO::File->new("| gzip >$outfile") or 
			die "unable to open output $outfile! $OS_ERROR\n";
	}
	else {
		$outfh = IO::File->new("$outfile", 'w') or 
			die "unable to open output $outfile! $OS_ERROR\n";
	}
}
else {
	# open an IO::Handle to STDOUT
	$outfh = IO::Handle->new;
	$outfh->fdopen(fileno(STDOUT), 'w');
	$outfile = 'stdout'; # just so we know
}

# suspect
if ($suspectfile) {
	if ($suspectfile =~ /\.gz$/) {
		$suspectfh = IO::File->new("| gzip >$suspectfile") or 
			die "unable to open output $suspectfile! $OS_ERROR\n";
	}
	else {
		$suspectfh = IO::File->new("$suspectfile", 'w') or 
			die "unable to open output $suspectfile! $OS_ERROR\n";
	}
}

# fail
if ($failfile) {
	if ($failfile =~ /\.gz$/) {
		$failfh = IO::File->new("| gzip >$failfile") or 
			die "unable to open output $failfile! $OS_ERROR\n";
	}
	else {
		$failfh = IO::File->new("$failfile", 'w') or 
			die "unable to open output $failfile! $OS_ERROR\n";
	}
}








### Counts and global values
my $goodCount      = 0;
my $notFoundCount  = 0;
my $missingR1Count = 0;
my $fix_length = length($fixed);
my $rfixed = reverse scalar($fixed);
$rfixed =~ tr/ATGC/TACG/;
my $rfix_short = substr($rfixed, 0, 12); # trim to the first 12 bases







### Process the files
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
# 	unless (amatch($header1, [2], $header2)) {
# 		die "malformed fastq files! headers from files don't match! compare\n" . 
# 			" $file1: $header1\n $file2: $header2\n";
# 	}
	
	### Check for adapter sequence in read2
	
	# first attempt perfect match
	my $i = index($sequence2, $fixed);
	if ($i <= 0) {
		# not found, so try an imperfect match allowing for substitutions
		# this should be entirely within the read, so there should be deletions
		$i = aindex($fixed, [2, 'D0', 'S2', 'I0'], $sequence2);
	}
	
	# process the read accordingly
	if ($i > 0) {
		# the adaptor in read2 was found
		
		# get UMI sequence
		# This may not be exact since the random portion could be 11 bp instead of 12
		# resulting in the fixed moving a little bit, which is compensated by the errors  
		# in the match above. This means we take one bp from fix as random. I think 
		# that's ok, since it's pretty rare.
		my $umi = substr($sequence2, 0, $umi_length);
		
		# now process read1 and trim its adapter
		my $i1 = index($sequence1, $rfixed);
		
		if ($i1 < 1) {
			# not found
			# try again with trimmed sequence
			$i1 = index($sequence1, $rfix_short);
			
			if ($i1 < 1) {
				# still not found? try a third time
				# with mismatches and deletions
				$i1 = aindex($rfix_short, [4, 'D4', 'S2', 'I0'], $sequence1);
			}
			# in theory I should just keep trying with shorter and shorter possibilities
			# but this seems potentially dangerous as you lose specificity for 
			# diminishing returns
			# in the end it's up to the aligner to align or soft trim
			# and it's likely not a real miRNA anyway....
		}
		
		
		# write read1 to output
		chomp($header1);
		my ($name1, $desc1) = split /\s/, $header1;
		if ($i1 > 1) {
			# found it, strip everything after the adapter
			my $trim_sequence1 = substr($sequence1, 0, $i1);
			my $trim_quality1  = substr($quality1, 0, $i1);
			# write output
			$outfh->printf("%s:%s %s\n%s\n%s%s\n", $name1, $umi, $desc1, 
				$trim_sequence1, $spacer1, $trim_quality1);
			$goodCount++;
		}
		else {
			# suspect RNA fragment
			# it's probably not a short RNA read, but does that disqualify it?
			if ($suspectfh) {
				$suspectfh->printf("%s:%s %s\n%s%s%s", $name1, $umi, $desc1, 
					$sequence1, $spacer1, $quality1);
			}
			else {
				$outfh->printf("%s:%s %s\n%s%s%s", $name1, $umi, $desc1, 
					$sequence1, $spacer1, $quality1);
			}
			$missingR1Count++;
		}
	}
	else {
		# adapter was never found in read2, so fail
		$notFoundCount++;
		if ($failfh) {
			$failfh->print("$header1$sequence1$spacer1$quality1");
		}
	}
}

# finish
if ($outfile eq 'stdout') {
	warn sprintf(" %10s spots had a detectable adapter sequence in both reads\n",
		$goodCount);
	warn sprintf(" %10s spots had a detectable adapter sequence in read2 only\n",
		$missingR1Count);
	warn sprintf(" %10s spots had no detectable adapter sequence\n",
		$notFoundCount);
	warn sprintf(" %10s spots were written to %s\n", 
		$suspectfile ? $goodCount : $goodCount + $missingR1Count, 
		$outfile);
}
else {
	printf(" %10s spots had a detectable adapter sequence in both reads\n",
		$goodCount);
	printf(" %10s spots had a detectable adapter sequence in read2 only\n",
		$missingR1Count);
	printf(" %10s spots had no detectable adapter sequence\n",
		$notFoundCount);
	printf(" %10s spots were written to %s\n", 
		$suspectfile ? $goodCount : $goodCount + $missingR1Count, 
		$outfile);
}






