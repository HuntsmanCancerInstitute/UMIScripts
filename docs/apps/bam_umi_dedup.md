# UMI Scripts - bam\_umi\_dedup

|[Home](Readme.md)|[Install](Install.md)|[Usage](Usage.md)|[Applications](Applications.md)|

## bam\_umi\_dedup

A script to remove duplicates based on a Unique Molecular Index (UMI) code.
Alignments that match the same coordinate are checked for the random UMI
code in the read name; Reads with the same UMI are sorted and selected for
one to be retained.

UMI sequences may be embedded in alignments in one of two locations:

1. SAM attribute `RX`
	 Added as an alignment attribute tag, typically the `RX` tag per current
    SAM specifications. This is the standard method, but requires 
    compatible aligners. UMI sequence qualities (`QX` tag) are ignored.

2. Appended to name
	Appended to the alignment read name as `:NNNN`, where NNNN is a
	sequence of indeterminate length comprised of [A,T,G,C]. This is a
	non-standard method but generally compatible with all aligners. 

At each chromosomal position, one representative alignment is selected
amongst all represented UMI sequences and the remainder are discarded
(default) or marked as duplicate (bit flag `0x400`) and retained. 
Alignments without a detectable UMI flag are simply written out. 

UMI sequence tags can tolerate mismatches up to the indicated number.
Insertions or deletions are also tolerated, but can be ignored by 
increasing their penalty score. In general, longer UMI sequences
should tolerate more mismatches. Allowing for mismatches can increase 
runtime considerably; hence a maximum depth threshold is enforced for 
detecting mismatches to preserve execution efficiency.

Selection criteria amongst UMI-duplicates include the mapping qualities of
the alignment (and possibly mate pair using tag `MQ` if present) or the sum
of base qualities of the read (and possibly mate pair using the `samtools
fixmate` tag `ms`).

Optical duplicate checking may be optionally included by specifying the
pixel distance threshold for marking as duplicates. Tile coordinates must
be present in the alignment read name. Optical duplicates are not 
distinguished from UMI duplicates when marking. 

Alignment files must be sorted by coordinate and indexed. Unmapped
(flag `0x4`) alignments are silently discarded. Read groups (tag `RG`) are
ignored; de-duplication should be done with only one read group per
alignment file.

Both Bam and Cram file formats are supported as input and output files.
Cram files require an indexed reference file.

**DISCLAIMER**: Alignments are de-duplicated solely on alignment coordinates
and the UMI sequences of the alignments at the current position, as well as
properly paired mate pairs. No guarantees are made for maintaining the same
molecule between secondary, supplementary (chimeric), and mate pair
alignments on separate chromosomes. These results may be sufficient for most 
applications, but not all. Additionally, alignments may fail subsequent 
verification checks because of this. 

VERSION: 2.3

### USAGE:

	bam_umi_dedup.pl --in in.bam --out out.bam

### OPTIONS:

Required:

    -i --in <file>        The input bam file, should be sorted and indexed
    -o --out <file>       The output bam file

UMI options:

    -u --umi <string>     SAM tag name for UMI sequence. Default 'RX'
                            Specify 'name' when UMI appended to read name.
    -m --mark             Mark duplicates (flag 0x400) instead of discarding
    -t --tolerance <int>  UMI sequence edit distance tolerance (1)
       --indel <int>      Set insertion/deletion penalty score (1)
       --skip <int>       Skip mismatch detection if depth exceeds (5000)

Other options:

    -f --fasta <file>     Provide indexed fasta file for Cram files
    -d --distance <int>   Set optical duplicate distance threshold.
                            Use 100 for unpatterned flowcell (HiSeq) or 
                            2500 for patterned flowcell (NextSeq or NovaSeq6000)
                            or 200 for NovaseqX. Default 0.
       --coord <string>   Provide the tile:X:Y integer 1-base positions in the 
                            read name for optical checking. For Illumina CASAVA 1.8 
                            7-element names, this is 5:6:7 (default)
    -c --cpu <int>        Specify the number of forks to use (4) 
       --samtools <path>  Path to samtools (/usr/local/bin/samtools)
       --nosam            Do not use samtools for final concatenation (slower)
    -h --help             Display full description and help


