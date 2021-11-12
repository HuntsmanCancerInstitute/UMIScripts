# Introduction

Universal Molecular Indexes, or UMIs, are for uniquely identifying the molecule from
which a sequencing read is derived, and is primarily used to differentiate reads from
unique biological molecules as opposed to a duplicate molecule derived from PCR
amplification or flow-cell spot duplication (optical duplicate). In essence, they are
oligomers of any nucleotide (`NNNNN`) in the adapter oligo. UMI sequences may be
included inline with the primary sequencing read (usually at the beginning of the
read, either with or without an interposing fixed sequence), or read as a separate
sequencing read, similar to sample barcode index tags for de-multiplexing. 

UMIs are common in single-cell sequencing. Since single-cell sequencing have their
own workflows, this document deals with regular DNASeq and RNASeq alignments.

The UMI sequence typically doesn't provide enough entropy to be used solely for
de-duplication, so it must be used in conjunction with the read sequence. This is
best done after the alignment step, since the alignment coordinates can then be used
in combination with the UMI to determine uniqueness. 

# Preparing for alignment

Keeping the UMI during the alignment step is tricky. It can't be part of the primary
read, since it won't align to the genome, so it has to be associated through the use
of alignment tags. The current [SAM
specification](http://samtools.github.io/hts-specs/SAMtags.pdf) uses the `RX` tag for
unique molecular index sequences. 

Unfortunately, by default, the Fastq file format doesn't allow for the use of extra
metadata, so alternatives must be made.

1. Attach the UMI sequence to the read name. This makes it part of the read name and
will be incorporated in the alignment. This is a non-standard method, but universally
accepted by all aligners since they all keep the read name.

2. Convert the UMI to SAM formatted attribute tags and stick them in the Fastq read
comment or description field, the part after the read name, delineated by a space or
tab. Some aligners can incorporate this into the SAM output.

3. Convert the Fastq file into an unaligned Bam file, where the UMI is incorporated
as a proper attribute tag. Some aligners can read directly from unaligned Bam files.

The method chosen is mostly dependent on the aligner software to be used, and to a
lesser extent, the de-duplication software used. 


# UMIScripts tool example

On local HCI Linux machines or using our CHPC module repository, you will need to
load the module.

    module load umiscripts

### Separate UMI Fastq file

For a separate Fastq file, one can use the `merge_umi_fastq.pl` tool. It can take
multiple Fastq files and sort them appropriately as Read1, Read2, and UMI based on
name and size, or they can be explicitly set with appropriate command line options.
The output file name or base name can optionally be provided.

- By default, it will generate SAM tags in the Fastq comment. 

        merge_umi_fastq.pl *.fastq.gz

- To append the UMI sequence to the read name, use the `--name` option. 

        merge_umi_fastq.pl --name *.fastq.gz

- To write as an unaligned Bam file, use the `--bam` option. A copy of
[samtools](https://github.com/samtools/samtools) must be in the `PATH`.

        merge_umi_fastq.pl --bam *.fastq.gz

The application can process 1 million paired-end reads in about 15-20 seconds. 

### Embedded UMI Fastq file

When the UMI is embedded inline to the read, it must be extracted. It may be in the
first read, second read, or both. The `embedded_UMI_extractor.pl` application can do
this. Currently, it only appends the UMI to the read name. The
length of the UMI sequence, and the fixed sequence (if any) should be given.

    embedded_UMI_extractor.pl --input R1.fastq.gz --pair R2.fastq.gz --length 12 --fixed GATC

# Aligners

The choice of aligner dictates the method used. Below are examples for the four
common aligners used in the Core. 

Output is piped into `samtools fixmate` to append the mate score, which can be used
for ranking alignments during de-duplication later on.

### BWA

BWA will accept UMI in the read name, as a Fastq comment, and as an unaligned Bam
file with an extra step.

To align Fastq with SAM tag comments:

    bwa mem -t $CPU -v 1 -C reference.fasta *.fastq.gz | \
    samtools fixmate -m - output.bam

To align with unaligned Bam file, it must first be converted to Fastq with tags as
comments using `samtools` as a pre-step:

    samtools fastq -T RX,RQ unaligned.bam | \
    bwa mem -t $CPU -C -v 1 -p reference.fasta - | \
    samtools fixmate -m - output.bam

### STAR

STAR does not accept Fastq comments, but it will accept an unaligned Bam file, as of
version 2.6. It will also of course accept UMI in the read name. You must designate
the command for reading the input file, and designate single-end "SE" or paired-end
"PE". 

    STAR --runMode alignReads \
    --readFilesIn unaligned.bam \
    --readFilesCommand samtools view \
    --readFilesType SAM PE \
    --readFilesSAMattrKeep RX QX \
    ...

STAR will also accept a SAM text file as an alternative to a Bam file, in which 
case the `--readFilesCommand` should be changed as appropriate.

### Bowtie2

Bowtie2 will accept unaligned Bam files as of version 2.3, and Fastq comments as of
version 2.4. It will of course accept UMI in the name.

To align with unaligned Bam files:

    bowtie2 --preserve-tags --align-paired-reads \
    -x reference.fasta -b unaligned.bam | \
    samtools fixmate -m - output.bam

To align with Fastq with SAM tag comments:

    bowtie2 --sam-append-comment \
    -x reference.fasta -1 read1.fastq.gz -2 read2.fastq.gz | \
    samtools fixmate -m - output.bam

### Novoalign

Novoalign appears to have the most versatility, with the ability to handle any
format, including extracting the UMI appended to the read name, or even extracting an
embedded UMI sequences from the main sequence read. 

To extract UMI appended to the read name:

    novoalign -d index.nix -f *.fastq.gz \
    --umi2 RX --tune NOVASEQ -o SAM | \
    samtools fixmate -m - output.bam

To align with Fastq SAM comments:

    novoalign -d index.nix -f *.fastq.gz \
    -C --tune NOVASEQ -o SAM | \
    samtools fixmate -m - output.bam

To align an unaligned Bam file:

    novoalign -d index.nix -f unaligned.bam \
    -F BAM RX,QX --tune NOVASEQ -o SAM | \
    samtools fixmate -m - output.bam


# De-duplication of alignment files
Once you have alignment files, they may be de-duplicated. The Bam files must be
sorted and indexed. Duplicates at any given coordinate are checked for the UMI and
the best one is kept, while the remaining are marked or discarded. There are two
programs that can be used.

### Picard

The Picard tool
[UmiAwareMarkDuplicatesWithMateCigar](https://broadinstitute.github.io/picard/command
-line-overview.html#UmiAwareMarkDuplicatesWithMateCigar) is the most complete and
thorough application, at the great expense of speed, i.e. it is painfully slow! It
can only accommodate UMI sequences in a SAM attribute tag, RX by default.

    java -jar picard.jar UmiAwareMarkDuplicatesWithMateCigar \
    -I input.bam -O output.bam \
    -M marked_dup_metrics.txt -UMI_METRICS umi_metrics.txt \
    --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 --UMI_TAG_NAME RX

### UMIScripts

The included de-deduplication script is considerably faster, but with a notable 
caveat: no guarantee is made for retaining identical molecules at secondary, 
supplementary, and inter-chromosomal alignments (they are treated independently). 
Otherwise, for normal alignments, results are comparable to Picard. For most 
count-based applications such as ChIPSeq or RNASeq, this limitation may be acceptable.

By default, duplicates are discarded, or they can marked with the `--mark` option.
Unmapped alignments are silently discarded. 

For de-duplication with SAM attribute tag `RX` (default):

    bam_umi_dedup.pl --in input.bam --out output.bam --cpu $CPU

For de-duplication with the UMI appended to the read name:

    bam_umi_dedup.pl --in input.bam --out output.bam --umi name --cpu $CPU


For a human (hg38) WGS Bam file with 668M alignments, the Picard tool (version 2.26.3) 
completes in 7.5 hours, while the UMIScripts tool completes in 48 minutes with 16 
threads. Results may vary.



