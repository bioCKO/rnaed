#!/usr/bin/perl

use strict 'vars';

my %remove;

open(fh, "dbSNP.cDNA_evidence_only") or die("Cannot open dbSNP.cDNA_evidence_only");

while (my $line = <fh>) {
	chomp($line);
	$line =~ /^(\S+) /;
	my $rs = $1;

	$remove{$rs} = 1;
}

close(fh);

open(fh, "All_20160407.vcf") or die("Cannot open All_20160407.vcf");

while (my $line = <fh>) {
	chomp($line);

	next if ($line =~ /^#/);

	$line =~ /^(\S+)\t(\S+)\t(\S+)\t/;
	my ($chr, $pos, $rs) = ($1,$2,$3);

	if ($remove{$rs} eq "") {
		print("$chr $pos\n");
	}
	else {
		print(STDERR "$rs removed\n");
	}
}

close(fh);