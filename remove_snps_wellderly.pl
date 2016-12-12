#!/usr/bin/perl

use strict 'vars';

my %remove;

open(fh, "/sc/orga/projects/STARNET/oscar/rna_editing/dbSNP/scripps.wellderly.GRCh38.bed") or die("Cannot open /sc/orga/projects/STARNET/oscar/rna_editing/dbSNP/scripps.wellderly.GRCh38.bed");

while (my $line = <fh>) {
	chomp($line);

	$line =~ /^(\S+)\t\S+\t(\S+)/;
	my ($chr, $pos) = ($1,$2);

	$remove{$chr."_".$pos} = 1;
}

close(fh);

my @files = <./calls.prefinal/*.no_NHLBI>;

foreach my $file (@files) {
	print("$file\n");

	open(fh, $file);
	open(fh_out, ">" . $file . ".no_wellderly");

	while (my $line = <fh>) {
		chomp($line);

		$line =~ /^(\S+) (\S+)/;
		my ($chr, $pos) = ($1,$2);

		if ($remove{$chr."_".$pos} eq "") {
			print(fh_out "$line\n");
		}
	}

	close(fh);
	close(fh_out);
}