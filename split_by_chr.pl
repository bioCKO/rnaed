#!/usr/bin/perl

use strict 'vars';

die("Usage: <file>") if (@ARGV != 1);

use Bio::SeqIO;

my $is = Bio::SeqIO->new(-file => @ARGV[0], -format => 'fasta');

while (my $obj = $is->next_seq()) {
	my $chr = $obj->display_id;

	print("$chr\n");

	my $os = Bio::SeqIO->new(-file => ">" . @ARGV[0] . "_" . $chr, -format => 'fasta');
	$os->write_seq($obj);
}
