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
	embedded_UMI_extractor.pl
	merge_umi_fastq.pl
	qiagen_smallRNA_umi_extractor.pl
	smallRNA_pe_umi_extractor.pl
);

### Compilation check
foreach my $script (@scripts) {
	script_compiles("bin/$script", "$script compiles");
}


### check bam_umi_dedup.pl if Bio::DB::HTS is available
my $check = 0;
eval { require Bio::DB::HTS; $check = 1; };
if ($check) {
	script_compiles("bin/bam_umi_dedup.pl", "bam_umi_dedup.pl compiles");
}


### Process fastq files
# to be implemented



### Bam deduplication
# to be implemented


done_testing();

