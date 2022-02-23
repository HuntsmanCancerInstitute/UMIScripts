#!/usr/bin/perl

# Test script for UMIScripts

use strict;
use Test::More;
# use Test2::V0;
use Test::Script;
use File::Spec;
use FindBin '$Bin';

print "Using $Bin\n";

### Scripts
my @scripts = qw(
	bam_umi_dedup.pl
	embedded_UMI_extractor.pl
	merge_umi_fastq.pl
	qiagen_smallRNA_umi_extractor.pl
	smallRNA_pe_umi_extractor.pl
);



### Compilation check
foreach my $script (@scripts) {
	script_compiles("bin/$script", "$script compiles");
}


### Process fastq files
# to be implemented



### Bam deduplication
# to be implemented


done_testing();

