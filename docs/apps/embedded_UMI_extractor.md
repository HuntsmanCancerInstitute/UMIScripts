# UMI Scripts - embeded\_UMI\_extractor

|[Home](Readme.md)|[Install](Install.md)|[Usage](Usage.md)|[Applications](Applications.md)|

## embeded\_UMI\_extractor

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

1. SAM attribute in Fastq comment
	Inserted as standard SAM attribute tags RX and QX in the read header.
	This is compatible with e.g. BWA, Bowties2 aligners. DEFAULT.

2. SAM attribute in unaligned Bam format
	Exported as an unaligned SAM/BAM format with RX and QX attribute tags.
	This is compatible with e.g. Bowtie2, STAR, BWA, and Novoalign.
	Use the `--bam` option.

3. Append to name
	Appended to the read name as ":UMI_sequence". This is compatible with 
	all aligners but requires special software to utilize. Non-standard. 
	Use the `--name` option.

Alignments can be de-duplicated utilizing the UMI codes with external software.
See [bam_umi_dedup.pl](bam_umi_dedup.md) in this software package, or 
Picard 'UmiAwareMarkDuplicatesWithMateCigar' as possibilities.

External utilities may be required for execution, inlcuding `pigz` and/or `gzip`,
and `samtools`. These are searched for in your environment `PATH`.


Version: 3.01

### Usage

	embedded_UMI_extractor.pl -l 8 -f GTAGTA *.fastq.gz 

	embedded_UMI_extractor.pl -1 file1.fastq.gz -2 file2.fastq.gz \
	--bam --out unaligned.bam --length 8 -f GTAGTA -F ACGACG 

### Options

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

    --samtools <path>       Path to samtools (/usr/local/bin/samtools)
    -h --help               Show full description and help

