#!/usr/bin/perl

use Bio::SeqIO;

my $file = @ARGV[0];
$file =~ /^\S+_(\S+)$/;
my $chr = $1;

my $is = Bio::SeqIO->new(-file => $file, -format => 'fasta');

while (my $obj = $is->next_seq()) {
	my $seq = $obj->seq;

	while ($seq =~ m/(A{5,})/gi) {
		my $stop = pos($seq);
		my $start = $stop - length($1) + 1;

		$start -- ;

		print("$chr\t$start\t$stop\n");
	}

	while ($seq =~ m/(C{5,})/gi) {
		my $stop = pos($seq);
		my $start = $stop - length($1) + 1;

		$start -- ;

		print("$chr\t$start\t$stop\n");
	}

	while ($seq =~ m/(G{5,})/gi) {
		my $stop = pos($seq);
		my $start = $stop - length($1) + 1;

		$start -- ;

		print("$chr\t$start\t$stop\n");
	}

	while ($seq =~ m/(T{5,})/gi) {
		my $stop = pos($seq);
		my $start = $stop - length($1) + 1;

		$start -- ;

		print("$chr\t$start\t$stop\n");
	}
}