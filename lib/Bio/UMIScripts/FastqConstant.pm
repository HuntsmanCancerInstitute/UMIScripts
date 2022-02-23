package Bio::UMIScripts::FastqConstant;
our $VERSION = 5;

=head1 Bio::UMIScripts::FastqConstant - Fastq Read object array indexes

=head1 DESCRIPTION

Exported constant values for working with Bio::UMIScripts::FastqRead objects.
Not intended for general consumption.

=back

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Bioinformatics Shared Resource
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  

=cut

use strict;
require Exporter;
use constant {
	NAME    => 0,
	DESC    => 1,
	SEQ     => 2,
	SPACER  => 3,
	QUAL    => 4
};

our @ISA = qw(Exporter);
our @EXPORT = qw(NAME DESC SEQ SPACER QUAL);

1;
