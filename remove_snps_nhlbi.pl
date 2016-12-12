#!/usr/bin/perl

use strict 'vars';

my %exac;

open(fh, "ESP6500SI-V2-SSA137.GRCh38-liftover.vcf") or die("Cannot open ESP6500SI-V2-SSA137.GRCh38-liftover.vcf");

<fh>;

while (my $line = <fh>) {
	chomp($line);

	$line =~ /\s(\S+):(\S+)$/;
	my ($chr, $pos) = ($1,$2);

	$exac{$chr."_".$pos} = 1;
}

close(fh);

my @files = <./calls.prefinal/*.no_exac>;

foreach my $file (@files) {
	print("$file\n");

	open(fh, $file);
	open(fh_out, ">" . $file . ".no_NHLBI");

	while (my $line = <fh>) {
		chomp($line);

		$line =~ /^(\S+) (\S+)/;
		my ($chr, $pos) = ($1,$2);

		$chr =~ s/chr//;

		if ($exac{$chr."_".$pos} eq "") {
			print(fh_out "$line\n");
		}
	}

	close(fh);
	close(fh_out);
}