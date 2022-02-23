package Bio::UMIScripts::UMIHelper;
our $VERSION = 5;

=head1 Bio::UMIScripts::UMIHelper - UMI helper routines for UMIScripts

=head1 DESCRIPTION

Exported subroutines for working with UMI sequences in the UMIScript package.

=head1 METHODS

No subroutines are exported by default.

=over 4

=item extract_umi_from_fastq($read, $umi_length)

Extracts the UMI sequence from a Fastq Read. Pass the Fastq Read and the expected integer 
length of the UMI sequence. 

The UMI sequence is extracted, and the UMI and fixed sequence removed from the original 
fastq Read. Corresponding changes are made to the Fastq Read quality as well. The UMI 
Fastq Read is generated and returned.

    printf "%s\n", $Read1->seq; # ATGAAGCCACAACGTACCCAAACCTATGGGACACAATGAAAGCA
    my $umi = extract_umi_from_fastq($Read1, 10);
    printf "%s\n", $umi->seq;   # ATGAAGCCAC
    printf "%s\n", $Read1->seq; #           AACGTACCCAAACCTATGGGACACAATGAAAGCA
    

=item extract_umi_with_fixed_from_fastq($read, $umi_length, $fixed_sequence)

Extracts the UMI sequence from a Fastq Read. Pass the Fastq Read, the expected integer 
length of the UMI sequence, and the fixed sequence between the UMI and read sequence. 
The Fixed sequence is search for in the Fastq sequence, first at the expected position. 
Failing that, it uses a fuzzy matching, tolerating one mismatch or one deletion, but 
not both or any insertions.

The UMI sequence is extracted, and the UMI and fixed sequence removed from the original 
fastq Read. Corresponding changes are made to the original Fastq Read quality as well. 
The UMI Fastq Read is generated and returned.

    printf "%s\n", $Read1->seq; # ATGAAGCCACAACGTACCCAAACCTATGGGACACAATGAAAGCA
    my $umi = extract_umi_with_fixed_from_fastq($Read1, 10, 'AACGTAC');
    printf "%s\n", $umi->seq;   # ATGAAGCCAC
    printf "%s\n", $Read1->seq; #                  CCAAACCTATGGGACACAATGAAAGCA
    

=item name_append_umi_from_fastq_read($read1, $umi)

Appends the sequence from the Read2 Fastq (UMI) to the name of the Read1 Fastq.

    printf "%s\n", $Read1->name; # A001:414:xxxxx:1:1101:2989:2566
    printf "%s\n", $umi->seq;   # ATGAAGCCAC
    name_append_umi_from_fastq_read($Read1, $umi);
    printf "%s\n", $Read1->name; # A001:414:xxxxx:1:1101:2989:2566:ATGAAGCCAC

=item umi_sam_tags_from_fastq_read($read)

Pass a fastq Read object representing a UMI sequence with quality. Returns a string 
representation of SAM tags for a UMI. Specifically, sequence is put into a C<RX> tag 
and quality is put into a C<QX> tag.

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
use String::Approx qw(amatch);
use Bio::UMIScripts::FastqConstant;
use Bio::UMIScripts::FastqRead;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw(
	umi_sam_tags_from_fastq_read
	extract_umi_with_fixed_from_fastq
	extract_umi_from_fastq
	name_append_umi_from_fastq_read
);


# the true value
1;



sub umi_sam_tags_from_fastq_read {
	return sprintf "RX:Z:%s\tQX:Z:%s", $_[0]->[SEQ], $_[0]->[QUAL];
}

sub extract_umi_with_fixed_from_fastq {
	my ($read, $umi_length, $fixed) = @_;
	
	# prepare new
	my $umi = $read->new($read->name, $read->description, '', '', '');
	
	# look for the UMI fixed sequence, if present pull the UMI sequence 
	# I have tried a number of strategies here
	
	# String::Approx aindex() doesn't seem to work well, especially with short 
	# fixed strings that we usually see, even worse if there are runs of identical 
	# bases. It always seems to match too early, even when you restrict the search 
	# space. This would've been ideal, if it worked, since the UMI may have 1 bp 
	# deletions (or insertions) that shifts the fixed code from expected position.
	
	# Text::Levenshtein::Flexible->distance() works, but fails tests when presented 
	# with indels in the primer, essentially shifts of the fixed sequence. 
	
	# Text::Levenshtein:::Flexible->l_distance() solves this problem, mostly. You 
	# can only give weights to the scoring, and then test the maximal distance. Since 
	# a deletion in the UMI shifts the fixed sequence, this requires a 1 deletion 
	# and one substitution, or a distance of 2. But this also indiscriminately 
	# allows 2 substitutions, not ideal.
	
	# String::Approx amatch() seems to be the best all around solution, since I can 
	# limit the maximum number of insertions, deletions, and substitutions. It 
	# doesn't involve scanning like aindex() so no early false positives. For 
	# UMI deletions that shift the fixed seq, this can be matched easily allowing 
	# only 1 deletion and 1 substitution, without allowing 2 substitutions. 
	
	# Searching the entire read increases run time and returns false positives
	# so limit to where we expect the fixed sequence to be
	my $seq_to_check = substr $read->[SEQ], $umi_length, length($fixed);
	
	# check fixed sequence
	if ($seq_to_check eq $fixed) {
		# we have exact match to fixed sequence
		my $len = $umi_length + length($fixed);
		my $s = substr $read->[SEQ], 0, $len, q();
		my $q = substr $read->[QUAL], 0, $len, q();
		$umi->[SEQ]  = substr $s, 0, $umi_length;
		$umi->[QUAL] = substr $q, 0, $umi_length;
		return $umi;
	}
	elsif (amatch($fixed, ['I0', 'D1', 'S1'], $seq_to_check) ) {
		# we have more-or-less the correct fixed code
		# allowing for zero insertions, one deletion, and one substitution max
		my $len = $umi_length + length $fixed;
		my $s = substr $read->[SEQ], 0, $len, q();
		my $q = substr $read->[QUAL], 0, $len, q();
		$umi->[SEQ]  = substr $s, 0, $umi_length;
		$umi->[QUAL] = substr $q, 0, $umi_length;
		return $umi;
	}
	else {
		# still not found
		return;
	}
}

sub extract_umi_from_fastq {
	my ($read, $umi_length) = @_;
	my $umi = Bio::UMIScripts::FastqRead->new($read->[NAME], $read->[DESC], q(), q(), q());
	# extract sequences while altering original sequence
	$umi->[QUAL] = substr $read->[QUAL], 0, $umi_length, q();
	$umi->[SEQ]  = substr $read->[SEQ], 0, $umi_length, q();
	return $umi;
}

sub name_append_umi_from_fastq_read {
	# read1, read2
	$_[0]->[NAME] .= sprintf ":%s", $_[1]->[SEQ];
	return 1;
}

1;
