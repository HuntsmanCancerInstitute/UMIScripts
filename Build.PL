#!/usr/bin/perl

use strict;
use warnings;
use Module::Build;

my $build = Module::Build->new(
	dist_name         => 'UMIScripts',
	dist_version      => 1,
	dist_author       => 'Timothy Parnell <parnell.tj@gmail.com>',
	license           => 'artistic_2',
	meta_merge        => {
		resources     => {
			repository     => 'https://github.com/tjparnell/'
		}
	},
	requires          => {
		'Bio::ToolBox'              => 1.60,
		'Parallel::ForkManager'     => 0,
		'Bio::DB::HTS'              => 0,
	},
	recommends        => {
		'String::Approx'            => 0,
	},
);

$build->create_build_script;
exit 0;
