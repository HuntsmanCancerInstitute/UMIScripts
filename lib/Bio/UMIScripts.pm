package Bio::UMIScripts;
our $VERSION = 5.01;

=head1 Bio::UMIScripts - Modules to support UMIScripts

=head1 SYNOPSIS

	use Bio::UMIScripts::FastqHelper qw(
	    read_fastq_filehandle
	    write_fastq_filehandle
	    write_bam_filehandle
	    get_fastq_read
	);
	use Bio::UMIScripts::UMIHelper qw(
	    extract_umi_from_fastq
	    extract_umi_with_fixed_from_fastq
	    umi_sam_tags_from_fastq_read
	);

	# open a Fastq IO::File handle
	my $fh = read_fastq_filehandle("file1.fastq.gz");

	# open output IO:File handles
	my $out1 = write_fastq_filehandle("output.fastq.gz");
	my $bam1 = write_bam_filehandle("unaligned.bam");

	while (my $read = get_fastq_read($fh)) {
		# each read is a Bio::UMIScripts::FastqRead object
	
		# basics
		my $name = $read->name;     # 'A001:414:xxxxx:1:1101:2989:2566'
		my $seq  = $read->sequence; # 'ATGAAGCCACAACGTACCCAAACCTATGGGACACAATGAAAGCA'
		my $qual = $read->quality;  # 'FFFFFFFF:FFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFF'
	
		# extract first 8 bp as UMI read and alter original read
		my $umi = extract_umi_from_fastq($read, 8);
	
		# or extract with fixed sequence with fuzzy matching
		$umi = extract_umi_with_fixed_from_fastq($read, 8, 'ACAACG');
	
		# export UMI as SAM tags - 'RX:Z:ATGAAGCC	RQ:Z:FFFFFFFF'
		my $tags = umi_sam_tags_from_fastq_read($umi);
	
		# write to fastq with SAM tag as comment
		$out1->print( $read->fastq_string($tags) );
	
		# or write to unaligned bam file with SAM tags
		$bam1->print( $read->sam_string(4, $tags) );
	}

=head1 DESCRIPTION

A package for working with Fastq files, and specifically Unique Molecular Index 
sequences. These are primarily intended for use with the packaged scripts, although 
they could be used as described above.


=head1 AUTHOR

 Timothy J. Parnell, PhD
 Bioinformatics Shared Resource
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  

=cut


1;
