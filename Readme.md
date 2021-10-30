# UMI Scripts

Simple scripts for handling Unique Molecular Indexes in Fastq and Bam files for 
simple PCR de-duplication purposes.

## Description

Unique molecular indexes (UMIs) are short, random, nucleotide sequences incorporated 
into next generation sequence reads for uniquely identifying biological molecules of 
RNA or DNA. This aids in distinguishing biological duplicates from PCR derived duplicates.

This package of scripts aid in processing these UMIs, both from Fastq and Bam files. 
For simplicity and convenience, the UMI is appended to the read name. Other processing 
schemes may use Bam attribute tags.

These scripts are designed to be fairly simplistic and work as quickly as possible. 
As such, UMIs are tested for exact matches without regard to base qualities or mismatches, 
which is probably suitable in most cases where users are simply looking to remove PCR 
artifacts. Users needing more advanced tools or reporting may need to look elsewhere.

## Applications

Below is a description of the included programs. Execute each script without any options to 
display the internal help page and list of available options.

### Fastq

These pre-alignment scripts are intended to work on Fastq files and extract a UMI and
either append it to the read name or write it as a third file. This is usually done prior
to trimming adapters. GZip compressed files are natively handled. Processed reads in 
some cases can be piped out to downstream applications (adapter trimming or alignment).

- `embedded_UMI_extractor.pl`

    A simple script for extracting the UMI barcode from the beginning of the read. 
    A fixed sequence may or may not be present between the UMI and the sequenced 
    insert. Both single-end and paired-end (with or without a second UMI in the 
    second read) fastq files are supported. The UMI is appended to the read name.

- `merge_umi_fastq.pl`

    A script for merging a third fastq file of UMI sequences with one (single-end) 
    or two (paired-end) fastq reads. The UMI may be appended to the read name 
    (non-standard), added as SAM tags to the Fastq read header, or written out as 
    an unaligned Sam or Bam file, depending on need and available support by the aligner.

- `smallRNA_pe_umi_extractor.pl`

    A simple script for extracting the adapter and UMI sequence from a 
    paired-end implementation of the Qiagen small RNA library preparation. 
    The UMI and Qiagen adapter is identified in R2, any complete or partial 
    Qiagen adapter identified and removed from R1, and R1 is written to output 
    with the UMI appended to the read name for de-duplication later.

- `qiagen_smallRNA_umi_extractor.pl`

    A simple script for extracting the UMI from a Fastq read derived from 
    Qiagen's Small RNA library preparation kit. Typically, a 60 base or longer read 
    is needed. Since small RNAs are typically under 30 bases, the sequence goes 
    through the RNA insert, Qiagen adapter, and index code. This script searches for the 
    Qiagen adapter sequence, allowing for up to two mismatches, and extracts the remaining 
    sequence as the UMI. 
    
    For a paired-end implementation of this, see `smallRNA_pe_umi_extractor.pl`.

### Bam

These post-alignment scripts are intended to work on Bam files for deduplication based 
on coordinates and UMI sequence. They work on files only (no streaming). 

- `bam_umi_dedup.pl`

    A generic, multi-threaded, UMI-aware Bam de-duplication application. Alignments 
    may be marked or removed based on coordinate, strand, and UMI sequence. UMI 
    duplicates are selected based on highest quality scores (mapping and base quality). 
    All single-end, proper paired-end, supplementary, and secondary alignments are 
    processed, with caveats. 

### Deprecated

These are deprecated scripts, kept here for posterity, but superseded by other 
scripts. Generally don't use.

- `SingleEndFastqBarcodeTagger.pl`

    A simple script for incorporating a second Fastq read of UMIs into the first 
    Fastq file, appending the UMI to the read name.

- `qiaseq_bam_deduplication.pl`

    A paired-end, UMI-aware Bam de-duplication application designed to work with the 
    qiaseq fastq scripts above. UMI duplicates are selected for the highest quality 
    sequence. Superseded by the generic `bam_umi_dedup.pl`.

- `qiaseq_barcode_extractor.pl`

    A simple script to extract the UMI code from the beginning of the second Fastq read 
    (paired-end) derived from the Qiagen Qiaseq PCR-amplicon based library preparation 
    kit. It will remove the UMI code and write it as a third Fastq file. This is 
    designed to work as a Qiagen-specific pre-processor for using the 
    [USeq](https://github.com/HuntsmanCancerInstitute/USeq) FastqBarcodeTagger, 
    MatchMates, and Consensus applications for collapsing PCR duplicate reads into 
    a single consensus read for re-alignment. See 
    [Workflows](https://github.com/HuntsmanCancerInstitute/Workflows) for examples.

- `qiaseq_FastqBarcodeTagger.pl`

    A simple script, similar to `qiaseq_barcode_extractor.pl`, to extract the UMI 
    from Qiagen Qiaseq PCR-amplicon based library preparation and appending the 
    UMI to the read name. Intended as a Qiagen-specific replacement processor for the 
    [USeq](https://github.com/HuntsmanCancerInstitute/USeq) tool FastqBarcodeTagger 
    and analogous to `embedded_UMI_extractor.pl`. Reads may then be aligned as normal 
    and de-duplicated using the below scripts.



## Installation

The package only includes scripts in the `bin` directory, which can be used as is. 
The scripts for bam file processing require additional library requirements of 
[Bio::ToolBox](https://metacpan.org/pod/Bio::ToolBox) and 
[Bio::DB::HTS](https://metacpan.org/pod/Bio::DB::HTS). Installation instructions 
can be found at the [Bio::ToolBox repository](https://github.com/tjparnell/biotoolbox).
Adapter identification with mismatches requires 
[String::Approx](https://metacpan.org/pod/String::Approx).
A standard [Module::Build](https://metacpan.org/pod/Module::Build) script is also 
provided for automated installation.

    perl ./Build.PL
    ./Build
    ./Build install

# AUTHOR

    Timothy J. Parnell, PhD
    Bioinformatics Shared Resource
    Huntsman Cancer Institute
    University of Utah
    Salt Lake City, UT, 84112

# LICENSE

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0. For details, see the
full text of the license in the file LICENSE.

This package is distributed in the hope that it will be useful, but it
is provided "as is" and without any express or implied warranties. For
details, see the full text of the license in the file LICENSE.
