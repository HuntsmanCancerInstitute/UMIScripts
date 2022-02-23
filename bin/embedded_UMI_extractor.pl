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
use Getopt::Long qw(:config no_ignore_case);
use File::Basename qw(fileparse);
use Bio::UMIScripts::FastqHelper qw(
	read_fastq_filehandle
	write_fastq_filehandle
	write_bam_filehandle
	get_fastq_read
	$SAMTOOLS_APP
);
use Bio::UMIScripts::UMIHelper qw(
	umi_sam_tags_from_fastq_read
	extract_umi_with_fixed_from_fastq
	extract_umi_from_fastq
	name_append_umi_from_fastq_read
);

my $VERSION = 3.01;

####### Documentation
my $description =  <<END;

A script to pull out an UMI embedded within a sequence read fastq file.
 
Unique Molecular Indexes (UMIs) are short, random nucleotide sequences that can
help to distinguish independent DNA molecules with the same sequence
composition. By associating the UMI to a nucleotide molecule at the beginning of
library preparation, identical sequence reads can be distinguished between
unique, biologically-derived molecules and PCR-derived or sequencing artifact
duplicates.

UMIs are expected to be at the beginning of the sequence read. The UMI barcode 
may or may not be followed by a short fixed sequence precdeding the sequence 
read. Either single-end or paired-end fastq files may be processed. Paired-end 
files may have an UMI in either one or both files. If a fixed sequence is provided, 
it must be found for the read to be written to output; paired-end UMI containing 
files must have both found.

UMIs can be merged into sequence read Fastq files in three ways:

   1 Inserted as standard SAM attribute tags RX and QX in the read header.
     This is compatible with e.g. BWA, Bowties2 aligners. DEFAULT.

   2 Exported as an unaligned SAM/BAM format with RX and QX attribute tags.
     This is compatible with e.g. Bowtie2, STAR, BWA, and Novoalign.
     Use the '--bam' option.

   3 Appended to the read name as ":UMI_sequence". This is compatible with 
     all aligners but requires special software to utilize. Non-standard. 
     Use the '--name' option.

Alignments can be de-duplicated utilizing the UMI codes with external software.
See 'bam_umi_dedup.pl' in this software package, or 
Picard 'UmiAwareMarkDuplicatesWithMateCigar' as possibilities.

External utilities may be required for execution, inlcuding pigz and/or gzip, 
and samtools. These are searched for in your environment PATH.

END

my $usage = <<END;

Version: $VERSION

Usage: 

  embedded_UMI_extractor.pl -l 8 -f GTAGTA *.fastq.gz 
  
  embedded_UMI_extractor.pl -1 file1.fastq.gz -2 file2.fastq.gz \\
    --bam --out unaligned.bam --length 8 -f GTAGTA -F ACGACG 

Options:

  Input:
    -1 --read1 <file>       First fastq read
    -2 --read2 <file>       Second fastq read, optional
    -u --umi <integer>      Read containing UMI: 1, 2, or 12 (default 1 or 12)
  
  Output:
    -o --out <file>         Path to output (input basename)
                              use 'stdout' for (interleaved) piping
    -b --bam                Write output as unaligned BAM format
                              or SAM format if stdout or samtools unavailable 
    -n --name               Append UMI to read name instead of SAM tag RX
  
  UMI Options:
    -l --length  <integer>  Length of UMI at beginning of read1
    -L --length2 <integer>  Length of UMI at beginning of read2 (if different)
    -f --fixed   <text>     Fixed sequence after UMI for read1 (optional) 
    -F --fixed2  <text>     Fixed sequence after UMI for read2 (optional, if different) 
  
  Other:
    --samtools <path>       Path to samtools ($SAMTOOLS_APP)
    -h --help               Show full description and help

END


#### Parameters
my $read_file1;
my $read_file2;
my $umi_location;
my $outfile; 
my $sam_format;
my $append_name;
my $fixed1 = q(); 
my $fixed2 = q();
my $umi1_length = 0;
my $umi2_length = 0;
my $cpu = 3;
my $help;

unless (@ARGV) {
	print $usage;
	exit 0;
}

GetOptions( 
	'1|read1=s'         => \$read_file1, # first fastq file
	'2|read2=s'         => \$read_file2, # second fastq file
	'u|umi=i'           => \$umi_location, # umi_location
	'o|out=s'           => \$outfile, # output fastq basename 
	'b|bam!'            => \$sam_format, # bam output
	'n|name'            => \$append_name, # add UMI to name
	'f|fixed=s'         => \$fixed1, # fixed portion of anchor sequence
	'F|fixed2=s'        => \$fixed2, # fixed portion of anchor sequence
	'l|length=i'        => \$umi1_length, # length of the barcode
	'L|length2=i'       => \$umi2_length, # length of the barcode
	'samtools=s'        => \$SAMTOOLS_APP, # samtools application
	'cpu=i'             => \$cpu, # number of CPU cores for compression
	'h|help'            => \$help, # print help
) or die "unrecognized options!\n";


if ($help) {
	print $description;
	print $usage;
	exit 0;
}





#### Check parameters
unless ($read_file1) {
	if (@ARGV) {
		$read_file1 = shift @ARGV;
		$read_file2 = shift @ARGV || undef;
	}
	else {
		die "must provide at least one input file!\n";
	}
}
unless ($umi1_length > 0) {
	die "must provide a length for the UMI barcode!\n";
}
unless ($umi_location) {
	$umi_location = $read_file2 ? 12 : 1;
}
unless ($umi_location == 1 or $umi_location == 2 or $umi_location == 12) {
	die "unrecognized UMI location of $umi_location! Must 1, 2, or 12\n";
}
if ($umi_location == 2 or $umi_location == 12) {
	$umi2_length ||= $umi1_length;
	$fixed2 ||= $fixed1;
}
if ($sam_format and not $SAMTOOLS_APP) {
	if ($outfile and $outfile =~ m/\.bam$/) {
		$outfile =~ s/bam$/sam.gz/;
		print STDERR "samtools application is not present. Writing to $outfile\n";
	}
}

# callback for extracting UMI
my $extract_umi;
if ($fixed1 or $fixed2) {
	# using a fixed sequence between UMI and read
	$extract_umi = \&extract_umi_with_fixed_from_fastq;
}
else {
	# no fixed sequence
	$extract_umi = \&extract_umi_from_fastq;
}




### Open input filehandles
my ($in_fh1, $in_fh2);
$in_fh1 = read_fastq_filehandle($read_file1);
if ($read_file2) {
	$in_fh2 = read_fastq_filehandle($read_file2);
}





### Open output filehandles
my ($out_fh1, $out_fh2);

if (not $outfile) {
	if ($read_file2 and not $sam_format) {
		# open two separate output handles
		$outfile = $read_file1;
		$outfile =~ s/\.(txt|fq|fastq)/.umi.$1/;
		$out_fh1 = write_fastq_filehandle($outfile, $cpu);
		my $outfile2 = $read_file2;
		$outfile2 =~ s/\.(txt|fq|fastq)/.umi.$1/;
		$out_fh2 = write_fastq_filehandle($outfile2, $cpu);
	}
	elsif (not $read_file2 and not $sam_format) {
		# just one file output
		$outfile = $read_file1;
		$outfile =~ s/\.(txt|fq|fastq)/.umi.$1/;
		$out_fh1 = write_fastq_filehandle($outfile, $cpu);
	}
	elsif ($read_file2 and $sam_format) {
		# open one bam output handle for both
		$outfile = $read_file1;
		$outfile =~ s/\.(?:txt|fq|fastq)(?:\.gz)?$/.bam/;
		$out_fh1 = write_bam_filehandle($outfile, $cpu);
		$out_fh2 = $out_fh2;
	}
	elsif (not $read_file2 and $sam_format) {
		# just one file output
		$outfile = $read_file1;
		$outfile =~ s/\.(?:fq|fastq)(?:\.gz)?$/.bam/;
		$out_fh1 = write_bam_filehandle($outfile, $cpu);
	}
	else {
		die "programming error:\n outfile $outfile\n sam format $sam_format\n first file $read_file1\n second file $read_file2\n";
	}
}
elsif ($outfile and lc($outfile) eq 'stdout') {
	$out_fh1 = write_fastq_filehandle($outfile); # this can handle stdout too
	if ($read_file2) {
		# reuse the same output 
		$out_fh2 = $out_fh1;
	}
}
elsif ($outfile and $sam_format) {
	unless ($outfile =~ /\.bam$/) {
		$outfile .= '.bam';
	}
	$out_fh1 = write_bam_filehandle($outfile, $cpu);
	if ($read_file2) {
		# reuse the same output 
		$out_fh2 = $out_fh1;
	}
}
elsif ($outfile and not $sam_format) {
	# writing fastq text file output for both files
	# use given out file as a base and make up extension
	my $outfile1 = $outfile . '.1.fastq.gz';
	$out_fh1 = write_fastq_filehandle($outfile1, $cpu);
	if ($read_file2) {
		my $outfile2 = $outfile . '.2.fastq.gz';
		$out_fh2 = write_fastq_filehandle($outfile2, $cpu);
	}
}
else {
	die "programming error:\n outfile $outfile\n sam format $sam_format\n first file $read_file1\n second file $read_file2\n";
}


### Write sam header as necessary
if ($sam_format) {
	my $cl = "$0 --read1 $read_file1 ";
	$cl .= "--read2  $read_file2 " if $read_file2;
	$cl .= "--out $outfile " if $outfile;
	$cl .= "--umi $umi_location ";
	$cl .= "--fixed $fixed1 " if $fixed1;
	$cl .= "--fixed2 $fixed2 " if $fixed2;
	$cl .= "--length $umi1_length " if $umi1_length;
	$cl .= "--length2 $umi2_length " if $umi2_length;
	$cl .= "--bam " if $sam_format;
	$cl .= "--name " if $append_name;
	$out_fh1->print("\@PG\tID:embedded_UMI_extractor\tVN:$VERSION\tCL:$cl\n");
}

# define SAM flags
	# 0x4d	77	PAIRED,UNMAP,MUNMAP,READ1
	# 0x8d	141	PAIRED,UNMAP,MUNMAP,READ2
	# 0x4   4   UNMAP 
my $flag1 = $read_file2 ? 77 : 4; 
my $flag2 = 141;






####### Process file based on inputs

# initialize global counters
my $goodCount = 0;
my $badCount  = 0;
my $singleCount = 0;

if ($read_file1 and not $read_file2) {
	# One Fastq only
	print STDERR "processing $read_file1....\n";
	if ($sam_format) {
		# write SAM file
		process_single_end_sam();
	}
	elsif ($append_name) {
		# append UMI to read name and write Fastq
		process_single_end_name_append();
	}
	else {
		# put SAM tags in Fastq comment
		process_single_end_fastq_tag();
	}
}
elsif ($read_file1 and $read_file2 and $umi_location != 12) {
	# Paired Fastq, first is barcoded only
	print STDERR "processing $read_file1 for UMIs with $read_file2 ....\n";
	if ($sam_format) {
		# write SAM file
		process_paired_end_with_one_umi_sam();
	}
	elsif ($append_name) {
		# append UMI to read name and write Fastq
		process_paired_end_with_one_umi_name_append();
	}
	else {
		# put SAM tags in Fastq comment
		process_paired_end_with_one_umi_fastq_tag();
	}
}
elsif ($read_file1 and $read_file2 and $umi_location == 12) {
	# Paired Fastq, both are barcoded
	print STDERR "processing both $read_file1 and $read_file2 for UMIs....\n";
	if ($sam_format) {
		# write SAM file
		process_paired_end_with_two_umi_sam();
	}
	elsif ($append_name) {
		# append UMI to read name and write Fastq
		process_paired_end_with_two_umi_name_append();
	}
	else {
		# put SAM tags in Fastq comment
		process_paired_end_with_two_umi_fastq_tag();
	}
}





#### Finish

# close filehandles
$in_fh1->close;
$in_fh2->close if $in_fh2;
$out_fh1->close;
$out_fh2->close if $out_fh2;

# print status
my $total = $goodCount + $badCount + $singleCount;
if ($singleCount) {
	printf STDERR " %12s read pairs total\n %12s read pairs had a bar code and were processed (%.2f%%)\n %12s read pairs had one failed bar code (%.2f%%)\n %12s read pairs had no bar code in either read (%.2f%%)\n", 
		$total, $goodCount, ($goodCount/$total) * 100, $singleCount, 
		($singleCount/$total) * 100, $badCount, ($badCount/$total) * 100;
}
else {
	printf STDERR " %12s reads total\n %12s reads had a bar code and were processed (%.2f%%)\n %12s reads had no bar code (%.2f%%)\n", 
		$total, $goodCount, ($goodCount/$total) * 100, $badCount, ($badCount/$total) * 100;
}

exit 0;








############ Subroutines

sub process_single_end_sam {
	# iterate
	while (my $read = get_fastq_read($in_fh1)) {
		
		# extract UMI sequence
		my $umi = &$extract_umi($read, $umi1_length, $fixed1);
		if ($umi) {
			my $tag = umi_sam_tags_from_fastq_read($umi);
			$out_fh1->print( $read->sam_string($flag1, $tag) );
			$goodCount++;
		}
		else {
			$badCount++;
		}
	}
}

sub process_single_end_name_append {
	# iterate
	while (my $read = get_fastq_read($in_fh1)) {
		
		# extract UMI sequence
		my $umi = &$extract_umi($read, $umi1_length, $fixed1);
		if ($umi) {
			name_append_umi_from_fastq_read($read, $umi);
			$out_fh1->print( $read->fastq_string );
			$goodCount++;
		}
		else {
			$badCount++;
		}
	}
}

sub process_single_end_fastq_tag {
	# iterate
	while (my $read = get_fastq_read($in_fh1)) {
		
		# extract UMI sequence
		my $umi = &$extract_umi($read, $umi1_length, $fixed1);
		if ($umi) {
			my $tag = umi_sam_tags_from_fastq_read($umi);
			$out_fh1->print( $read->fastq_string($tag) );
			$goodCount++;
		}
		else {
			$badCount++;
		}
	}
}

sub process_paired_end_with_one_umi_sam {
	# iterate
	while (my $read1 = get_fastq_read($in_fh1)) {
		
		# second read
		my $read2 = get_fastq_read($in_fh2) or 
			die "premature end of file for $read_file2!";
		unless ($read1->compare_names($read2)) {
			printf STDERR " Mismatching read names!\n   Compare read1 %s with read2 %s\n",
				$read1->name, $read2->name;
			die "mismatched files!";
		}
		
		# extract UMI sequence
		my $umi = &$extract_umi($umi_location == 1 ? $read1 : $read2, 
			$umi1_length, $fixed1);
		if ($umi) {
			my $tag = umi_sam_tags_from_fastq_read($umi);
			$out_fh1->print( $read1->sam_string($flag1, $tag) );
			$out_fh2->print( $read2->sam_string($flag2, $tag) );
			$goodCount++;
		}
		else {
			$badCount++;
		}
	}
}

sub process_paired_end_with_one_umi_name_append {
	# iterate
	while (my $read1 = get_fastq_read($in_fh1)) {
		
		# second read
		my $read2 = get_fastq_read($in_fh2) or 
			die "premature end of file for $read_file2!";
		unless ($read1->compare_names($read2)) {
			printf STDERR " Mismatching read names!\n   Compare read1 %s with read2 %s\n",
				$read1->name, $read2->name;
			die "mismatched files!";
		}
		
		# extract UMI sequence
		my $umi = &$extract_umi($umi_location == 1 ? $read1 : $read2, 
			$umi1_length, $fixed1);
		if ($umi) {
			name_append_umi_from_fastq_read($read1, $umi);
			name_append_umi_from_fastq_read($read2, $umi);
			$out_fh1->print( $read1->fastq_string );
			$out_fh2->print( $read2->fastq_string );
			$goodCount++;
		}
		else {
			$badCount++;
		}
	}
}

sub process_paired_end_with_one_umi_fastq_tag {
	# iterate
	while (my $read1 = get_fastq_read($in_fh1)) {
		
		# second read
		my $read2 = get_fastq_read($in_fh2) or 
			die "premature end of file for $read_file2!";
		unless ($read1->compare_names($read2)) {
			printf STDERR " Mismatching read names!\n   Compare read1 %s with read2 %s\n",
				$read1->name, $read2->name;
			die "mismatched files!";
		}
		
		# extract UMI sequence
		my $umi = &$extract_umi($umi_location == 1 ? $read1 : $read2, 
			$umi1_length, $fixed1);
		if ($umi) {
			my $tag = umi_sam_tags_from_fastq_read($umi);
			$out_fh1->print( $read1->fastq_string($tag) );
			$out_fh2->print( $read2->fastq_string($tag) );
			$goodCount++;
		}
		else {
			$badCount++;
		}
	}
}



sub process_paired_end_with_two_umi_sam {
	# iterate
	while (my $read1 = get_fastq_read($in_fh1)) {
		
		# second read
		my $read2 = get_fastq_read($in_fh2) or 
			die "premature end of file for $read_file2!";
		unless ($read1->compare_names($read2)) {
			printf STDERR " Mismatching read names!\n   Compare read1 %s with read2 %s\n",
				$read1->name, $read2->name;
			die "mismatched files!";
		}
		
		# extract UMI sequence
		my $umi1 = &$extract_umi($read1, $umi1_length, $fixed1);
		my $umi2 = &$extract_umi($read2, $umi2_length, $fixed2);
		if ($umi1 and $umi2) {
			# both found, both good
			my $umi = $umi1->concatenate_reads($umi2);
			my $tag = umi_sam_tags_from_fastq_read($umi);
			$out_fh1->print( $read1->sam_string($flag1, $tag) );
			$out_fh2->print( $read2->sam_string($flag2, $tag) );
			$goodCount++;
		}
		elsif ($umi1 or $umi2) {
			# only one found
			$singleCount++;
		}
		else {
			# none found
			$badCount++;
		}
	}
}

sub process_paired_end_with_two_umi_name_append {
	# iterate
	while (my $read1 = get_fastq_read($in_fh1)) {
		
		# second read
		my $read2 = get_fastq_read($in_fh2) or 
			die "premature end of file for $read_file2!";
		unless ($read1->compare_names($read2)) {
			printf STDERR " Mismatching read names!\n   Compare read1 %s with read2 %s\n",
				$read1->name, $read2->name;
			die "mismatched files!";
		}
		
		# extract UMI sequence
		my $umi1 = &$extract_umi($read1, $umi1_length, $fixed1);
		my $umi2 = &$extract_umi($read2, $umi2_length, $fixed2);
		if ($umi1 and $umi2) {
			# both found, both good
			my $umi = $umi1->concatenate_reads($umi2);
			name_append_umi_from_fastq_read($read1, $umi);
			name_append_umi_from_fastq_read($read2, $umi);
			$out_fh1->print( $read1->fastq_string );
			$out_fh2->print( $read2->fastq_string );
			$goodCount++;
		}
		elsif ($umi1 or $umi2) {
			# only one found
			$singleCount++;
		}
		else {
			# none found
			$badCount++;
		}
	}
}

sub process_paired_end_with_two_umi_fastq_tag {
	# iterate
	while (my $read1 = get_fastq_read($in_fh1)) {
		
		# second read
		my $read2 = get_fastq_read($in_fh2) or 
			die "premature end of file for $read_file2!";
		unless ($read1->compare_names($read2)) {
			printf STDERR " Mismatching read names!\n   Compare read1 %s with read2 %s\n",
				$read1->name, $read2->name;
			die "mismatched files!";
		}
		
		# extract UMI sequence
		my $umi1 = &$extract_umi($read1, $umi1_length, $fixed1);
		my $umi2 = &$extract_umi($read2, $umi2_length, $fixed2);
		if ($umi1 and $umi2) {
			# both found, both good
			my $umi = $umi1->concatenate_reads($umi2);
			my $tag = umi_sam_tags_from_fastq_read($umi);
			$out_fh1->print( $read1->fastq_string($tag) );
			$out_fh2->print( $read2->fastq_string($tag) );
			$goodCount++;
		}
		elsif ($umi1 or $umi2) {
			# only one found
			$singleCount++;
		}
		else {
			# none found
			$badCount++;
		}
	}
}





