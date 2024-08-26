# UMI Scripts - merge\_UMI\_fastq

|[Home](Readme.md)|[Install](Install.md)|[Usage](Usage.md)|[Applications](Applications.md)|

## merge\_UMI\_fastq

A script to merge a separate UMI fastq file with sequence read fastq file(s). 

Unique Molecular Indexes (UMIs) are short, random nucleotide sequences that can
help to distinguish independent DNA molecules with the same sequence
composition. By associating the UMI to a nucleotide molecule at the beginning of
library preparation, identical sequence reads can be distinguished between
unique, biologically-derived molecules and PCR-derived or sequencing artifact
duplicates.

Either paired-end or single-end reads may be merged with a UMI Fastq file.

UMIs can be merged into sequence read Fastq files in three ways:

1. SAM attribute in Fastq comment
	Inserted as standard SAM attribute tags `RX` and `QX` in the read header.
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

Read and UMI Fastq files may be empirically determined based on size and name. 


VERSION: 1.01

### USAGE: 
    merge_umi_fastq.pl *.fastq.gz
    
    merge_umi_fastq.pl -1 R1.fq.gz -2 R2.fq.gz -u UMI.fq.gz -o unaligned.bam

### OPTIONS

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

    --samtools <path>     Path to samtools (/usr/local/bin/samtools)
    -h --help             Show full description and help


