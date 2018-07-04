# UMI Scripts

Scripts for handling Unique Molecular Indexes in Fastq and Bam files

## Description

Unique molecular indexes (UMIs) are short, random, nucleotide sequences incorporated 
into next generation sequence reads for uniquely identifying biological molecules of 
RNA or DNA. This aids in distinguishing biological duplicates from PCR derived duplicates.

This package of scripts aid in processing these UMIs, both from Fastq and Bam files. 
For simplicity and convenience, the UMI is appended to the read name. Other processing 
schemes may use Bam attribute tags.

## Applications

Below is a description of the included programs. Execute each script without any options to 
display the internal help page and list of available options.

### Fastq

- `SingleEndFastqBarcodeTagger.pl`

    A simple script for incorporating a second Fastq read of UMIs into the first 
    Fastq file, appending the UMI to the read name.

- `qiagen_smallRNA_umi_extractor.pl`

    A simple script for extracting the UMI from a Fastq read derived from 
    Qiagen's Small RNA library preparation kit. Typically, a 60 base or longer read 
    is needed. Since small RNAs are typically under 30 bases, the sequence goes 
    through the RNA insert, adapter, and index code. This script searches for the 
    adapter sequence, allowing for up to two mismatches, and extracts the remaining 
    sequence as the UMI.

- `qiaseq_barcode_extractor.pl`

    A simple script to extract the UMI code from the beginning of the second Fastq read 
    (paired-end) derived from the Qiagen Qiaseq PCR-amplicon based library preparation 
    kit. It will remove the UMI code and write it as a third Fastq file. This is 
    designed to work as a pre-processor for using the 
    [USeq](https://github.com/HuntsmanCancerInstitute/USeq) FastqBarcodeTagger, 
    MatchMates, and Consensus applications.

- `qiaseq_FastqBarcodeTagger.pl`

    A simple script, similar to `qiaseq_barcode_extractor.pl`, to extract the UMI 
    from Qiagen Qiaseq PCR-amplicon based library preparation and incorporating the 
    UMI into the read name. Intended as a replacement processor for the 
    [USeq](https://github.com/HuntsmanCancerInstitute/USeq) tool FastqBarcodeTagger.

### Bam

- `bam_umi_dedup.pl`

    A multi-threaded, single-end, UMI-aware Bam de-duplication application. Alignments 
    may be marked or removed based on coordinate, strand, and UMI sequence. UMI 
    duplicates are selected based on highest quality scores.

- `qiaseq_bam_deduplication.pl`

    A paired-end, UMI-aware Bam de-duplication application designed to work with the 
    qiaseq fastq scripts above. UMI duplicates are selected for the highest quality 
    sequence.

## Installation

The package only includes scripts in the `bin` directory, which can be used as is. 
The scripts for bam file processing require additional library requirements of 
[Bio::ToolBox](https://metacpan.org/pod/Bio::ToolBox) and 
[Bio::DB::HTS](https://metacpan.org/pod/Bio::DB::HTS). Installation instructions 
can be found at the [Bio::ToolBox repository](https://github.com/tjparnell/biotoolbox).
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
