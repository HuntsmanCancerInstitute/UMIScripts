package Bio::UMIScripts::FastqRead;
our $VERSION = 5;

=head1 Bio::UMIScripts::FastqRead - Fastq Read object

=head1 DESCRIPTION

This module represents Fastq reads. 

B<Note> that objects have minimal error checking in the name of speed and efficiency.
It is up to the user to verify if values make sense. Do not be surprised if 
something fails if given garbage input.

=head1 METHODS

=over 4

=item new(@lines)

Create a new read object. Generally not used independently but returned as 
a read from L<Bio::UMIScripts::FastqHelper::get_fastq_read()>;

=item name()

Returns the name of the read.

=item description()

Returns the name of the read.

=item sequence()

Returns the name of the read.

=item spacer()

Returns the name of the read.

=item quality()

Returns the name of the read.

=item duplicate_read()

Returns a new read object identical to the current one.

=item compare_names($read2)

Compares the names of two fastq reads to determine if they are identical. 
Returns 1 if editing (Levenshtein) distance is 1 or less, or 0 if not.

=item check_quality($min_basequal)

Checks whether a Fastq Read contains any C<N> bases and whether all base qualities 
match or exceed the provided minimum base quality (assumes standard Sanger quality 
where Phred+33 encoding). Returns 1 if true, 0 if false.

=item concatenate_reads($read2)

Concatenates two Fastq Reads, joining sequence and quality with a C<+>. Name and 
description are empty values. This is suitable for concatenating barcode sequences 
for SAM tags per specification, but not for concatenating standard sequence reads. 
Returns a new Fastq Read array. 

=item fastq_string()

=item fastq_string($comment)

Formats and returns the Fastq Read as a proper four-line Fastq text string. 
If a comment string is provided, it is used in place of the original
comment, e.g. if a SAM formatted comment string is to be used.

=item sam_string()

=item sam_string($flag)

=item sam_string($flag, $tag)

Formats the Fastq Read as a proper SAM line and returns it. The SAM flag 
should be provided as an integer and used to indicate read and alignment 
status (default C<4> representing C<UNMAPPED>). Any additional SAM formatted 
attribute tags may be provided for appending to the line. No checks are made 
for flag or tags. 

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
use List::Util qw(min);
use String::Approx qw(amatch);
use Bio::UMIScripts::FastqConstant;

sub new {
	my $class = shift;
	my @a = @_;
	return bless \@a, ref($class) || $class;
}

sub name {
	return substr $_[0]->[NAME], 1;
}

sub description {
	return $_[0]->[DESC];
}

sub sequence {
	return $_[0]->[SEQ];
}

sub spacer {
	if (length $_[0]->[SPACER] > 1) {
		return substr $_[0]->[SPACER], 1;
	}
	else {
		return $_[0]->[SPACER];
	}
}

sub quality {
	return $_[0]->[QUAL];
}

sub compare_names {
	# read1, read2
	# do a simple compare first
	return 1 if ($_[0]->[NAME] eq $_[1]->[NAME]);
	
	# otherwise compare allowing one substitution only and case insensitive
	if ( amatch($_[0]->[NAME], ['i', 'S1', 'I0', 'D0'], $_[1]->[NAME]) ) {
		return 1;
	}
	else {
		return 0;
	}
}

sub extract_illumina_sample {
	if ($_[0]->[DESC] =~ m/\d+:\w:\d+:([ATGCN]+[\+\-]?[ATGCN]*)/) {
		return "BC:Z:$1";
	}
	else {
		return '';
	}
}


sub check_quality {
	# my ($read, $min_basequal) = @_;
	if (index($_[0]->[SEQ],'N') == -1) {
		# good sequence, no N's
		my @quals = map {ord($_) - 33} split q(), $_[0]->[QUAL];
		if (min(@quals) >= $_[1]) {
			# good quality
			return 1;
		}
	}
	return 0; # otherwise return nothing indicating a bad sequence
}

sub concatenate_reads {
	# read1, read2
	my @a = (
		$_[0]->[NAME],
		'',
		$_[0]->[SEQ] . '+' . $_[1]->[SEQ],
		'',
		$_[0]->[QUAL] . '+' . $_[1]->[QUAL]
	);
	return $_[0]->new(@a);
}

sub fastq_string {
	# my ($read, $comment) = @_;
	return sprintf "%s %s\n%s\n%s\n%s\n",
		$_[0]->[NAME],
		$_[1] || $_[0]->[DESC],
		$_[0]->[SEQ],
		$_[0]->[SPACER] || '+',
		$_[0]->[QUAL];
}

sub sam_string {
	# my ($read, $flag, $tag) = @_;
	return sprintf "%s\n", join "\t",
		$_[0]->name,            # QNAME
		$_[1] || 4,             # FLAG
		'*',                    # RNAME
		0,                      # POS
		0,                      # MAPQ
		'*',                    # CIGAR
		'*',                    # RNEXT
		0,                      # PNEXT
		0,                      # TLEN
		$_[0]->[SEQ],           # SEQ
		$_[0]->[QUAL],          # QUAL
		$_[2] || ''             # TAGS
	;
}

sub duplicate_read {
	return $_[0]->new(
		$_[0]->[NAME]   || '',
		$_[0]->[DESC]   || '',
		$_[0]->[SEQ]    || '',
		$_[0]->[SPACER] || '',
		$_[0]->[QUAL]   || ''
	);
}

1;
