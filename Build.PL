#!/usr/bin/perl

use strict;
use warnings;
use Module::Build;

my $build = Module::Build->new(
	dist_name         => 'UMIScripts',
	module_name       => 'UMIScripts',
	dist_version_from => 'lib/Bio/UMIScripts.pm',
	dist_abstract     => 'Scripts for handling Unique Molecular Indexes in Fastq and Bam files',
	dist_author       => 'Timothy Parnell <timothy.parnell@hci.utah.edu>',
	license           => 'artistic_2',
	meta_merge        => {
		resources     => {
			repository     => 'https://github.com/HuntsmanCancerInstitute/UMIScripts'
		}
	},
	configure_requires => {
		'Module::Build'           => 0,
	},
	requires          => {
		'List::Util'                    => 0,
		'String::Approx'                => 0,
		'Text::Levenshtein::Flexible'   => 0,
	},
	test_requires     => {
		'Test::More'                    => 0,
		'Test::Script'                  => 0
	},
	recommends        => {
		'Bio::DB::HTS'                  => 0,
		'Parallel::ForkManager'         => 0,
		'List::MoreUtils'               => 0,
	},
);

$build->create_build_script;
exit 0;
