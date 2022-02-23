#!/usr/bin/perl -w

# Test script for UMIScripts

use strict;
use Test::More;
use File::Spec;
use FindBin '$Bin';
use Bio::UMIScripts::FastqConstant;
use Bio::UMIScripts::FastqHelper qw(
	read_fastq_filehandle
	write_fastq_filehandle
	write_bam_filehandle
	get_fastq_read
);
use Bio::UMIScripts::UMIHelper qw(
	umi_sam_tags_from_fastq_read
	extract_umi_with_fixed_from_fastq
	extract_umi_from_fastq
	name_append_umi_from_fastq_read
);


### Named indexes
is(NAME, 0, 'Named index');
is(DESC, 1, 'Description index');
is(SEQ, 2, 'Sequence index');
is(SPACER, 3, 'Spacer index');
is(QUAL, 4, 'Quality index');


### Example Fastq files
my $fastq1 = File::Spec->catfile($Bin, "data", "example_R1.fastq"); # 100 bp read 1
my $fastq2 = File::Spec->catfile($Bin, "data", "example_R2.fastq"); # 100 bp read 2
my $fastq3 = File::Spec->catfile($Bin, "data", "example_R3.fastq"); # matching UMI read
my $fastq4 = File::Spec->catfile($Bin, "data", "example_R4.fastq"); # different UMI read



### Open fastq files
my $fq_fh1 = read_fastq_filehandle($fastq1);
isa_ok($fq_fh1, 'IO::File', 'open fastq file example_R1');

my $fq_fh2 = read_fastq_filehandle($fastq2);
isa_ok($fq_fh2, 'IO::File', 'open fastq file example_R2');

my $fq_fh3 = read_fastq_filehandle($fastq3);
isa_ok($fq_fh3, 'IO::File', 'open fastq file example_R3');

my $fq_fh4 = read_fastq_filehandle($fastq4);
isa_ok($fq_fh4, 'IO::File', 'open fastq file example_R4');



### Check Read
my $read1 = get_fastq_read($fq_fh1);
isa_ok($read1, 'Bio::UMIScripts::FastqRead', 'read1 object');
is($read1->name, 'A001:414:xxxxx:1:1101:2989:2566', 'read1 name');
is($read1->description, '1:N:0:CGTATCTC+GTACCTTG', 'read1 description');
is($read1->sequence, 'ATGAAGCCACAACGTACCCAAACCTATGGGACACAATGAAAGCATTTCTA', 'read1 sequence string');
is($read1->quality, 'FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF', 'read1 quality string');
is($read1->spacer, "+", 'read1 spacer string');

my $read2 = get_fastq_read($fq_fh2);
isa_ok($read2, 'Bio::UMIScripts::FastqRead', 'read2 object');
is($read2->name, 'A001:414:xxxxx:1:1101:2989:2566', 'read2 name');
is($read2->sequence, 'ACCAACTCCTCGTTTGGTTAATTCTTTGAATAGTTCTTCTTGTTTCCACT', 'read2 sequence string');


my $read3 = get_fastq_read($fq_fh3);
isa_ok($read3, 'Bio::UMIScripts::FastqRead', 'read3 object');
is($read3->name, 'A001:414:xxxxx:1:1101:2989:2566', 'read3 name');
is($read3->sequence, 'TCCGTTTCTGT', 'read3 sequence string');

my $read4 = get_fastq_read($fq_fh4);
isa_ok($read4, 'Bio::UMIScripts::FastqRead', 'read4 object');
is($read4->name, 'A001:414:xxxxx:1:1101:3070:3991', 'read4 name');
is($read4->sequence, 'GCAGGCACCTN', 'read4 sequence string');



### Compare read names
ok($read1->compare_names($read2), 'compare identical names'); # identical
ok(!$read1->compare_names($read4), 'compare different names'); # very different
$read3->[NAME] = '@A001:414:xxxxx:1:1101:2989:256'; # 1 deletion only
ok(!$read1->compare_names($read3), 'compare similar names 1D'); # 1 differences only
$read3->[NAME] = '@A001:414:xxxxx:1:1101:2989:356'; # 1 deletion 1 substitution
ok(!$read1->compare_names($read3), 'compare similar names 1D 1S'); # 2 differences only



### Check read quality
ok($read1->check_quality(20), 'check read1 minimum quality of 20 true');
ok($read1->check_quality(30), 'check read1 minimum quality of 30 true');
ok(!$read1->check_quality(40), 'check read1 minimum quality of 40 false');
ok(!$read4->check_quality(20), 'check read4 minimum quality of 20 false');



### Concatenate Read
my $cat1 = $read3->concatenate_reads($read4);
isa_ok($cat1, 'Bio::UMIScripts::FastqRead', 'concatenated read object');
is($cat1->name, 'A001:414:xxxxx:1:1101:2989:356', 'concatenated read name');
is($cat1->description, '', 'concatenated read description');
is($cat1->sequence, 'TCCGTTTCTGT+GCAGGCACCTN', 'concatenated read sequence string');
is($cat1->quality, 'FFFFF:FF:F:+,F,,,FFFFF#', 'concatenated read quality string');



### Extract Illumina barcode
is($read1->extract_illumina_sample, 'BC:Z:CGTATCTC+GTACCTTG', 'extract Illumina barcode as SAM tag');


### Convert read to Sam tags
my $sam_tags = umi_sam_tags_from_fastq_read($read3);
is($sam_tags, "RX:Z:TCCGTTTCTGT\tQX:Z:FFFFF:FF:F:", 'convert read to SAM tags');


### Duplicate a read
my $example = $read1->duplicate_read;
isa_ok($example, 'Bio::UMIScripts::FastqRead', 'duplicated read object');
is($example->name, 'A001:414:xxxxx:1:1101:2989:2566', 'duplicated read name');
is($example->sequence, 'ATGAAGCCACAACGTACCCAAACCTATGGGACACAATGAAAGCATTTCTA', 'duplicated read sequence string');
is($example->quality, 'FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF', 'duplicated read quality string');



### Extract UMI subsequence by length
my $umi = extract_umi_from_fastq($example, 10);
isa_ok($umi, 'Bio::UMIScripts::FastqRead', 'extracted UMI read object by length');
is($umi->sequence, 'ATGAAGCCAC', 'extracted UMI sequence');
is($umi->quality, 'FFFFFFFFFF', 'extracted UMI quality');



### Extract UMI subsequence with fixed sequence
# perfect match
$example = $read1->duplicate_read;
$umi = extract_umi_with_fixed_from_fastq($example, 11, 'ACGTAC');
isa_ok($umi, 'Bio::UMIScripts::FastqRead', 'extracted UMI read object by length and fixed sequence');
is($umi->sequence, 'ATGAAGCCACA', 'extracted UMI sequence');
is($umi->quality, 'FFFFFFFFFFF', 'extracted UMI quality');
is($example->sequence, 'CCAAACCTATGGGACACAATGAAAGCATTTCTA', 'trimmed sequence');
is($example->quality, 'FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF', 'trimmed quality');

# one mismatch in fixed
$example = $read1->duplicate_read;
$umi = extract_umi_with_fixed_from_fastq($example, 11, 'ACGTTC');
isa_ok($umi, 'Bio::UMIScripts::FastqRead', 'extracted UMI read object by length and fixed sequence mismatch');
is($umi->sequence, 'ATGAAGCCACA', 'extracted UMI sequence');
is($umi->quality, 'FFFFFFFFFFF', 'extracted UMI quality');
is($example->sequence, 'CCAAACCTATGGGACACAATGAAAGCATTTCTA', 'trimmed sequence');
is($example->quality, 'FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF', 'trimmed quality');

# fixed offset by extra base
$example = $read1->duplicate_read;
$umi = extract_umi_with_fixed_from_fastq($example, 10, 'ACGTAC');
isa_ok($umi, 'Bio::UMIScripts::FastqRead', 'extracted UMI read object by length and fixed sequence offset +1');
is($umi->sequence, 'ATGAAGCCAC', 'extracted UMI sequence');
is($umi->quality, 'FFFFFFFFFF', 'extracted UMI quality');
is($example->sequence, 'CCCAAACCTATGGGACACAATGAAAGCATTTCTA', 'trimmed sequence');
is($example->quality, 'FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF', 'trimmed quality');

# fixed offset by missing base
$example = $read1->duplicate_read;
$umi = extract_umi_with_fixed_from_fastq($example, 12, 'ACGTAC');
isa_ok($umi, 'Bio::UMIScripts::FastqRead', 'extracted UMI read object by length and fixed sequence offset -1');
is($umi->sequence, 'ATGAAGCCACAA', 'extracted UMI sequence');
is($umi->quality, 'FFFFFFFFFFFF', 'extracted UMI quality');
is($example->sequence, 'CAAACCTATGGGACACAATGAAAGCATTTCTA', 'trimmed sequence');
is($example->quality, 'FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF', 'trimmed quality');

# no match
$example = $read1->duplicate_read;
$umi = extract_umi_with_fixed_from_fastq($example, 11, 'GATCTA');
ok(!defined $umi, 'extract embedded UMI by length and non-matching sequence');



### Append name
$example = $read1->duplicate_read;
ok(name_append_umi_from_fastq_read($example, $read3), 'append UMI sequence to name');
is($example->name, 'A001:414:xxxxx:1:1101:2989:2566:TCCGTTTCTGT', 'appended new name');



### Fastq string
# default
my $expected = <<END;
\@A001:414:xxxxx:1:1101:2989:2566 1:N:0:CGTATCTC\+GTACCTTG
ATGAAGCCACAACGTACCCAAACCTATGGGACACAATGAAAGCATTTCTA
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
END
is($read1->fastq_string(), $expected, 'export to default fastq string');

# add sam string, generated above
$expected = <<END;
\@A001:414:xxxxx:1:1101:2989:2566 RX:Z:TCCGTTTCTGT	QX:Z:FFFFF:FF:F:
ATGAAGCCACAACGTACCCAAACCTATGGGACACAATGAAAGCATTTCTA
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
END
is($read1->fastq_string($sam_tags), $expected, 'export to fastq string with SAM tags');



### SAM string
# default
$expected = join("\t", 
	'A001:414:xxxxx:1:1101:2989:2566', # name
	4,                      # 0x4   4   UNMAP 
	'*',                    # RNAME
	0,                      # POS
	0,                      # MAPQ
	'*',                    # CIGAR
	'*',                    # RNEXT
	0,                      # PNEXT
	0,                      # TLEN
	'ATGAAGCCACAACGTACCCAAACCTATGGGACACAATGAAAGCATTTCTA', # SEQ
	'FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF', # QUAL
	"\n"                    # note extra tab at end
);
is($read1->sam_string(), $expected, 'export to default SAM string');

# read1 of two with SAM tags
$expected = join("\t", 
	'A001:414:xxxxx:1:1101:2989:2566', # name
	77,                     # 0x4d	77	PAIRED,UNMAP,MUNMAP,READ1
	'*',                    # RNAME
	0,                      # POS
	0,                      # MAPQ
	'*',                    # CIGAR
	'*',                    # RNEXT
	0,                      # PNEXT
	0,                      # TLEN
	'ATGAAGCCACAACGTACCCAAACCTATGGGACACAATGAAAGCATTTCTA', # SEQ
	'FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF', # QUAL
	'RX:Z:TCCGTTTCTGT',     # TAGS
	'QX:Z:FFFFF:FF:F:'             
) . "\n";
is($read1->sam_string(77, $sam_tags), $expected, 'export to SAM string with SAM tags');



### Finished
done_testing();







