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
# https://github.com/HuntsmanCancerInstitute/UMIScripts

use strict;
use IO::File;
use IO::Handle;
use Getopt::Long qw(:config no_ignore_case);
use Bio::UMIScripts::FastqHelper qw(
	read_fastq_filehandle
	write_fastq_filehandle
	write_bam_filehandle
	get_fastq_read
);
use Bio::UMIScripts::UMIHelper qw(
	umi_sam_tags_from_fastq_read
	extract_umi_with_fixed_from_fastq
	extract_umi_from_fastq
	name_append_umi_from_fastq_read
);

########################################

my $VERSION = 1;


### Options
my $read_file1;
my $read_file2;
my $umi_file1;
my $umi_file2;
my $outfile;
my $sam_format;
my $append_name;
my $keep_sample;
my $cpu = 3;
my $help;
my $start_time = time;


my $description = <<END;

A script to merge a separate UMI fastq file with sequence read fastq file(s). 

Unique Molecular Indexes (UMIs) are short, random nucleotide sequences that can
help to distinguish independent DNA molecules with the same sequence
composition. By associating the UMI to a nucleotide molecule at the beginning of
library preparation, identical sequence reads can be distinguished between
unique, biologically-derived molecules and PCR-derived or sequencing artifact
duplicates.

Either paired-end or single-end reads may be merged with a UMI Fastq file.

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

Read and UMI Fastq files may be empirically determined based on size and name. 

END

my $usage = <<END;

VERSION: $VERSION

USAGE: 
    merge_umi_fastq.pl *.fastq.gz
    
    merge_umi_fastq.pl -1 R1.fq.gz -2 R2.fq.gz -u UMI.fq.gz -o unaligned.bam

OPTIONS:
  Input:
    -1 --read1 <file>     First fastq read
    -2 --read2 <file>     Second fastq read, optional
    -u --umi <file>       UMI fastq read
    -U --umi2 <file>      Second UMI fastq file, optional
  
  Output:
    -o --out <file>       Output file (base)name (input basename)
                            use 'stdout' for (interleaved) piping
    -b --bam              Write output as unaligned BAM format
                            or SAM format if stdout or samtools unavailable 
    -n --name             Append UMI to read name instead of SAM tag RX
  
  Other:
    --samtools <path>     Path to samtools ($Bio::UMIScripts::FastqRead::SAMTOOLS_APP)
    -h --help             Show full description and help

END

unless (@ARGV) {
	print $usage;
	exit 0;
}

# get command line options
GetOptions( 
	'1|read1=s'         => \$read_file1, # first fastq file
	'2|read2=s'         => \$read_file2, # second fastq file
	'u|umi1=s'          => \$umi_file1, # first umi fastq file
	'U|umi2=s'          => \$umi_file2, # second umi fastq file
	'o|out=s'           => \$outfile, # output file base name
	'b|bam!'            => \$sam_format, # output sam format
	'n|name!'           => \$append_name, # append UMI to name
	'k|keepbc!'         => \$keep_sample, # keep the sample bar code
	'samtools=s'        => \$Bio::UMIScripts::FastqRead::SAMTOOLS_APP, # samtools application
	'cpu=i'             => \$cpu, # number of CPU cores for compression
	'h|help'            => \$help, # print help
) or die "bad options!\n";





### Check options
if (not $read_file1 and not $umi_file1 and scalar @ARGV) {
	# user provided unordered list of fastq files
	# sort them out by size
	my @list =
		sort {$b->[1] <=> $a->[1]}      # decreasing size
		map { [$_, (stat $_)[7] ] }    # array of filename and size <- 8th element from stat array
		@ARGV;
	if (scalar @list == 2) {
		# one fastq and one umi
		$read_file1 = $list[0]->[0];
		$umi_file1  = $list[1]->[0];
	}
	elsif (scalar @list == 3) {
		# two fastq and one umi
		# sort the fastq names asciibetically to get read1 and read2
		($read_file1, $read_file2) =
			map {$_->[0]}
			sort {$a->[0] cmp $b->[0]}
			($list[0], $list[1]);
		$umi_file1  = $list[2]->[0];
	}
	elsif (scalar @list == 4) {
		# two fastq and two umi
		# sort the fastq names asciibetically to get read1 and read2
		($read_file1, $read_file2) =
			map {$_->[0]}
			sort {$a->[0] cmp $b->[0]}
			($list[0], $list[1]);
		($umi_file1, $umi_file2) =
			map {$_->[0]}
			sort {$a->[0] cmp $b->[0]}
			($list[2], $list[3]);
	}
	else {
		my $i = scalar @list;
		die "ERROR! $i input files provided! Unable to identify appropriately!";
	}
	print STDERR "Empirically determined following provided files:\n";
	print STDERR " - read1: $read_file1\n";
	print STDERR " - read2: $read_file2\n" if $read_file2;
	print STDERR " - umi1:  $umi_file1\n";
	print STDERR " - umi2:  $umi_file2\n" if $umi_file2;
}

if (not $read_file1) {
	die "Must provide an input read fastq file!\n";
}
if (not $umi_file1) {
	die "Must provide a UMI fastq files!\n";
}
if ($outfile and $outfile =~ /\.[s|b]am(?:\.gz)?$/) {
	# convenience
	$sam_format = 1;
}
if ($sam_format and not $Bio::UMIScripts::FastqRead::SAMTOOLS_APP) {
	if ($outfile and $outfile =~ m/\.bam$/) {
		$outfile =~ s/\bam$/sam.gz/;
		print STDERR "samtools application is not present. Writing to $outfile\n";
	}
}



### Open input filehandles

# Read fastq files
my ($read_fh1, $read_fh2);
$read_fh1 = read_fastq_filehandle($read_file1);
$read_fh2 = read_fastq_filehandle($read_file2) if $read_file2;


# UMI fastq files
# also set the subroutine callback for working with either one or two UMI files
my ($umi_fh1, $umi_fh2, $read_umi);
$umi_fh1  = read_fastq_filehandle($umi_file1);
if ($umi_file2) {
	# working with two UMI fastq files
	$umi_fh2  = read_fastq_filehandle($umi_file2) ;
	
	# set the UMI read subroutine
	$read_umi = \&read_two_umi_fastq;
}
else {
	# working with one UMI fastq file
	$read_umi = \&read_one_umi_fastq;
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
elsif ($outfile and (lc($outfile) eq 'stdout' )) {
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
	$cl .= "--umi $umi_file1 ";
	$cl .= "--umi2 $umi_file2 " if $umi_file2;
	$cl .= "--out $outfile " if $outfile;
	$cl .= "--bam " if $sam_format;
	$cl .= "--name " if $append_name;
	$cl .= "--keepbc " if $keep_sample;
	$out_fh1->print("\@PG\tID:merge_umi_fastq\tVN:$VERSION\tCL:$cl\n");
}

# define SAM flags
	# 0x4d	77	PAIRED,UNMAP,MUNMAP,READ1
	# 0x8d	141	PAIRED,UNMAP,MUNMAP,READ2
	# 0x4   4   UNMAP 
my $flag1 = $read_file2 ? 77 : 4; 
my $flag2 = 141;







#### MAIN 

my $count;
if ($read_file2) {
	# paired-end
	if ($sam_format) {
		# write SAM file
		$count = process_paired_sam();
	}
	elsif ($append_name) {
		# append UMI to read name and write Fastq
		$count = process_paired_name_append();
	}
	else {
		# put SAM tags in Fastq comment
		$count = process_paired_fastq_tag();
	}
}
else {
	# single-end
	if ($sam_format) {
		# write SAM file
		$count = process_single_sam();
	}
	elsif ($append_name) {
		# append UMI to read name and write Fastq
		$count = process_single_name_append();
	}
	else {
		# put SAM tags in Fastq comment
		$count = process_single_fastq_tag();
	}
}

# close all file handles
$read_fh1->close;
$read_fh2->close if $read_fh2;
$umi_fh1->close;
$umi_fh2->close if $umi_fh2; 
$out_fh1->close;
$out_fh2->close if $out_fh2;


# Finish
printf STDERR " $count reads were processed in %.1f minutes\n", (time - $start_time) / 60;












################ Subroutines

sub read_two_umi_fastq {
	my $check = shift;
	my $umi1 = get_fastq_read($umi_fh1) or 
		die "premature end of file for $umi_file1!";
	my $umi2 = get_fastq_read($umi_fh2) or 
		die "premature end of file for $umi_file2!";
	if (
		$check->compare_names($umi1) and 
		$check->compare_names($umi2)
	) {
		# names match, let's merge the two
		# official SAM format says to concatenate with a +
		return $umi1->concatenate_reads($umi2);
	}
	else {
		printf STDERR " Mismatching read names!\n   Compare read %s with UMI1 %s and UMI2 %s\n",
			$check->name, $umi1->name, $umi2->name;
		return;
	}
}

sub read_one_umi_fastq {
	my $check = shift;
	my $umi1 = get_fastq_read($umi_fh1) or
		die "premature end of file for $umi_file1!";
	if ($check->compare_names($umi1)) {
		# names match
		return $umi1;
	}
	else {
		printf STDERR " Mismatching read names!\n   Compare read %s with UMI1 %s\n",
			$check->name, $umi1->name;
		return;
	}
}

sub process_paired_sam {
	# two paired-end fastq files written to sam/bam while keeping the sample barcode
	print STDERR " Processing paired-end fastq files\n";
	print STDERR " Writing output in SAM/BAM format with UMI tags\n";
	my $i = 0;
	while (my $read1 = get_fastq_read($read_fh1)) {
		# second read
		my $read2 = get_fastq_read($read_fh2) or
			die "premature end of file for $read_file2!";
		unless ($read1->compare_names($read2)) {
			printf STDERR " Mismatching read names!\n   Compare read1 %s with read2 %s\n",
				$read1->name, $read2->name;
			die "mismatched files!";
		}
		
		# UMI tags
		my $umi = &$read_umi($read1) or die "mismatched files!";
		my $tag = umi_sam_tags_from_fastq_read($umi);
		
		# first read
		$out_fh1->print( $read1->sam_string($flag1, $tag) );
		
		# second read
		$out_fh1->print( $read2->sam_string($flag2, $tag) );
		
		$i++;
	}
	return $i;
}

sub process_paired_name_append {
	# two paired-end fastq files written to fastq with UMI appended to read names 
	# sample barcode will automatically be kept, if present
	print STDERR " Processing paired-end fastq files\n";
	print STDERR " Appending UMI to read name\n";
	my $i = 0;
	while (my $read1 = get_fastq_read($read_fh1)) {
		# second read
		my $read2 = get_fastq_read($read_fh2) or
			die "premature end of file for $read_file2!";
		unless ($read1->compare_names($read2)) {
			printf STDERR " Mismatching read names!\n   Compare read1 %s with read2 %s\n",
				$read1->name, $read2->name;
			die "mismatched files!";
		}
		
		# UMI 
		my $umi = &$read_umi($read1) or die "mismatched files!";
		
		# first read
		name_append_umi_from_fastq_read($read1, $umi);
		$out_fh1->print( $read1->fastq_string );
		
		# second read
		name_append_umi_from_fastq_read($read2, $umi);
		$out_fh2->print( $read2->fastq_string );
		
		$i++;
	}
	return $i;
}

sub process_paired_fastq_tag {
	# two paired-end fastq files written to fastq with UMI as SAM tags
	print STDERR " Processing paired-end fastq files\n";
	print STDERR " Appending UMI as SAM tags to read comment\n";
	my $i = 0;
	while (my $read1 = get_fastq_read($read_fh1)) {
		# second read
		my $read2 = get_fastq_read($read_fh2) or
			die "premature end of file for $read_file2!";
		unless ($read1->compare_names($read2)) {
			printf STDERR " Mismatching read names!\n   Compare read1 %s with read2 %s\n",
				$read1->name, $read2->name;
			die "mismatched files!";
		}
		
		# UMI 
		my $umi = &$read_umi($read1) or die "mismatched files!";
		my $tag = umi_sam_tags_from_fastq_read($umi);
		
		# first read
		$out_fh1->print( $read1->fastq_string($tag) );
		
		# second read
		$out_fh2->print( $read2->fastq_string($tag) );
		
		$i++;
	}
	return $i;
}

sub process_single_name_append {
	# single-end fastq files written to fastq with UMI appended to read names 
	# sample barcode will automatically be kept, if present
	print STDERR " Processing single-end fastq file\n";
	print STDERR " Appending UMI to read name\n";
	my $i = 0;
	while (my $read1 = get_fastq_read($read_fh1)) {
		# UMI 
		my $umi = &$read_umi($read1) or die "mismatched files!";
		
		# first read
		name_append_umi_from_fastq_read($read1, $umi);
		$out_fh1->print( $read1->fastq_string );
		
		$i++;
	}
	return $i;
}

sub process_single_sam {
	# single-end fastq files written to sam/bam
	print STDERR " Processing single-end fastq file\n";
	print STDERR " Writing output in SAM/BAM format with UMI tags\n";
	my $i = 0;
	while (my $read1 = get_fastq_read($read_fh1)) {
		# UMI tags
		my $umi = &$read_umi($read1) or die "mismatched files!";
		my $tag = umi_sam_tags_from_fastq_read($umi);
		
		# first read
		$out_fh1->print( $read1->sam_string($flag1, $tag) );
		
		$i++;
	}
	return $i;
}

sub process_single_fastq_tag {
	# single-end fastq files written to fastq with UMI as SAM tags
	print STDERR " Processing single-end fastq file\n";
	print STDERR " Appending UMI as SAM tags to read comment\n";
	my $i = 0;
	while (my $read1 = get_fastq_read($read_fh1)) {
		# UMI 
		my $umi = &$read_umi($read1) or die "mismatched files!";
		my $tag = umi_sam_tags_from_fastq_read($umi);
		
		# first read
		$out_fh1->print( $read1->fastq_string($tag) );
		
		$i++;
	}
	return $i;
}




