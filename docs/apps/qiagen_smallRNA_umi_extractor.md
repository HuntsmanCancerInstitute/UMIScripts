# UMI Scripts - qiagen\_smallRNA\_umi\_extractor

|[Home](Readme.md)|[Install](Install.md)|[Usage](Usage.md)|[Applications](Applications.md)|

## qiagen\_smallRNA\_umi\_extractor

A script to strip out the Qiagen UMI from single-end reads 
prepared from a Qiagen small RNA library. Qiagen recommends doing 
a long read to capture the small RNA insert, read through the 
19 bp adapter sequence, and then through the 12 bp UMI sequence. 

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
[bam_umi_dedup.pl](bam_umi_dedup.md) as the companion application to do this.

### Usage:

	qiagen_smallRNA_umi_extractor.pl -i input.fastq.gz --fail noUMI.fastq.gz | <aligner>

### Options:

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
	                     mismatches are allowed. Default is 'AACTGTAGGCACCATCAAT'.
	--len <integer>      Length of the UMI barcode. Default 12.

