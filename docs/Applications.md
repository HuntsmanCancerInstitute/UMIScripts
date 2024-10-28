# UMI Scripts - Applications

|[Home](Readme.md)|[Install](Install.md)|[Usage](Usage.md)|[Applications](Applications.md)|


The following applications are included with this package free of charge.


### Fastq applications

These pre-alignment scripts are intended to work on Fastq files and extract a UMI and
either append it to the read name or write it as a third file. This is usually done prior
to trimming adapters. GZip compressed files are natively handled. Processed reads in 
some cases can be piped out to downstream applications (adapter trimming or alignment).

- [embedded_UMI_extractor.pl](apps/embedded_UMI_extractor.md)

    A simple script for extracting the UMI barcode from the beginning of the read. 
    A fixed sequence may or may not be present between the UMI and the sequenced 
    insert. Both single-end and paired-end (with or without a second UMI in the 
    second read) fastq files are supported. The UMI may added as SAM tags to the 
    Fastq read header, appended to the read name, or written out as an unaligned 
    Sam or Bam file.

- [merge_umi_fastq.pl](apps/merge_umi_fastq.md)

    A script for merging a third fastq file of UMI sequences with one (single-end) 
    or two (paired-end) fastq reads. The UMI may be appended to the read name 
    (non-standard), added as SAM tags to the Fastq read header, or written out as 
    an unaligned Sam or Bam file, depending on need and available support by the aligner.

- [qiagen_smallRNA_umi_extractor.pl](apps/qiagen_smallRNA_umi_extractor.md)

    A simple script for extracting the UMI from a Fastq read derived from
    Qiagen's Small RNA library preparation kit. Typically, a 60 base or
    longer read is needed. Since small RNAs are typically under 30 bases,
    the sequence goes through the RNA insert, Qiagen adapter, and index
    code. This script searches for the Qiagen adapter sequence, allowing
    for up to two mismatches, and extracts the remaining sequence as the
    UMI. Use the `smallRNA_pe_umi_extractor.pl` script when using using
    paired-end 50 bp reads to reliably extract the UMI from the second read.

- [smallRNA_pe_umi_extractor.pl](apps/smallRNA_pe_umi_extractor.md)

    A simple script for extracting the adapter and UMI sequence from a
    paired-end implementation of the Qiagen small RNA library preparation.
    This is intended for short (50 bp) paired-end sequence reads from the
    library, where variable insert size may preclude reading a complete UMI
    sequence in a single read. The UMI and Qiagen adapter is identified in
    R2, any complete or partial Qiagen adapter identified and removed from
    R1, and R1 is written to output with the UMI appended to the read name
    for de-duplication later.
    

### Bam applications

These post-alignment applications are intended to work on Bam files for deduplication
based on coordinates and UMI sequence. They work on files only (no streaming). 

- [bam_umi_dedup.pl](apps/bam_umi_dedup.md)

    A generic, multi-threaded, UMI-aware Bam de-duplication application. Alignments 
    may be marked or removed based on coordinate, strand, and UMI sequence. UMI 
    duplicates are selected based on highest quality scores (mapping and base quality). 
    All single-end, proper paired-end, supplementary, and secondary alignments are 
    processed, with caveats: supplementary, chimeric, secondary, and mate pairs on 
    separate chromosomes are treated independently. 

### Deprecated applications

These are deprecated scripts, kept here for posterity, but superseded by other 
scripts. Generally do not use. The scripts can be found in the `deprecated` folder.

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


