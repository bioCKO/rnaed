#!/usr/bin/perl

use strict 'vars';

my %exac;

open(fh, "ExAC.r0.3.1.sites.vep.vcf.GRCh38.bed") or die("Cannot open ExAC.r0.3.1.sites.vep.vcf.GRCh38.bed");

while (my $line = <fh>) {
	chomp($line);
	next if ($line =~ /^#/);

	$line =~ /^(\S+)\t\S+\t(\S+)/;
	my ($chr, $pos) = ($1,$2);

	$exac{$chr."_".$pos} = 1;
}

close(fh);

my @files = <./calls.prefinal/*.no_ucsc_141>;

foreach my $file (@files) {
	print("$file\n");

	open(fh, $file);
	open(fh_out, ">" . $file . ".no_exac");

	while (my $line = <fh>) {
		chomp($line);

		$line =~ /^(\S+) (\S+)/;
		my ($chr, $pos) = ($1,$2);

		if ($exac{$chr."_".$pos} eq "") {
			print(fh_out "$line\n");
		}
	}

	close(fh);
	close(fh_out);
}