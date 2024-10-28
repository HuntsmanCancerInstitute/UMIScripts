# UMI Scripts - smallRNA\_pe\_umi\_extractor

|[Home](Readme.md)|[Install](Install.md)|[Usage](Usage.md)|[Applications](Applications.md)|

## smallRNA\_pe\_umi\_extractor

A script to extract the Unique Molecular Index (UMI) from the second 
read of a short-read (50 bp) paired-end Qiagen small RNA library.
The structure of a typical small RNA paired-end read is as follows:

- R1
	`tagcttatcagactgatgttga AACTGTAGGCACCATCAAT NNNNNNNNN`
	`<--miRNA-------------> <--adapter--------> <--UMI-->`

- R2  
	`NNNNNNNNNNNN ATTGATGGTGCCTACAGTT tcaacatcagtctgataag`
    `<--UMI-----> <--adapter--------> <--miRNA---------->`

Since the small RNA insert may be of variable length, the full UMI 
sequence (12 bp by default) may not be reliably present within the full 
read length of R1, in which case extracting it from R2 becomes necessary.

The Qiagen adapter is identified in R2 (allowing for mismatches as
necessary), the preceding UMI is extracted, and the UMI appended to the
read name. The Qiagen adapter is also searched for and removed from R1 if
present, allowing for mismatches and deletions as necessary. 

Only one Fastq file, R1, is written out. Reads with adapter sequence 
in both R2 and R1 is always written to output. Reads with the adapter 
identified only in R2 is written to output or optionally to a second 
file (option `--suspect`). Reads where the adapter sequence is not found 
in either R2 or R1 may only be written to an optional file (option `--fail`). 
No UMI code is appended to failed reads.

Sequencing primers are not searched for, but normally should not be 
necessary except for long insertions (suspect reads).

The UMI is appended to the read name as ":NNNNNNNNNNNN". After 
alignment, the bam file may be de-duplicated [bam_umi_dedup.pl](bam_umi_dedup.md).


### Usage:

	smallRNA_pe_umi_extractor.pl --out <basename> --f1 <read_1> --f2 <read_2>
	
	smallRNA_pe_umi_extractor.pl <read_1> <read_2>

### Options:

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
                           Default is 'ATTGATGGTGCCTACAGTT'.
    -l |--len <integer>  Length of the Unique Molecular Index. Default 12.
    -h |--help           This help
