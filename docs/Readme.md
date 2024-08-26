# UMI Scripts

|[Home](Readme.md)|[Install](Install.md)|[Usage](Usage.md)|[Applications](Applications.md)|

Scripts for handling Unique Molecular Indexes in Fastq and Bam files, primarily for 
de-duplication purposes.

## Description

Unique molecular indexes (UMIs) are short, random, nucleotide sequences incorporated 
into next generation sequence reads for uniquely identifying biological molecules of 
RNA or DNA. This aids in distinguishing biological duplicates from PCR-derived
or technical (optical) duplicates.

This package of scripts aid in processing these UMIs, both from Fastq and Bam files.
Embedded UMIs can be extracted from the primary sequence read in a Fastq file, or
they may be provided as a separate Fastq file. They may be stored either appended
to the read name or assigned to SAM/Bam attribute tag.

These scripts are designed to be fairly simplistic and work as quickly as possible. 
Users needing more advanced tools or reporting may need to look elsewhere.

## Installation

See the [Installation](Install.md) documentation.

## Usage

See the [Usage](Usage.md) documentation.

## Applications

See the list of available [Applications](Applications.md).

Execute each script without options or with the `-h` or `--help` options to
print the included documentation.

## AUTHOR

    Timothy J. Parnell, PhD
    Cancer Bioinformatics Shared Resource
    Huntsman Cancer Institute
    University of Utah
    Salt Lake City, UT, 84112

## LICENSE

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0. For details, see the
full text of the license in the file LICENSE.

This package is distributed in the hope that it will be useful, but it
is provided "as is" and without any express or implied warranties. For
details, see the full text of the license in the file LICENSE.
