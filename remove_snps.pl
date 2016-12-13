#!/usr/bin/perl

use strict 'vars';

my %dbsnp;

my $prog = 0;

open(fh, "All_20160527.vcf.chr_and_pos") or die("Cannot open All_20160527.vcf.chr_and_pos");

while (my $line = <fh>) {
	chomp($line);
	next if ($line =~ /^#/);

	$line =~ /^(\S+) (\S+)/;
	my ($chr, $pos) = ($1,$2);

	$dbsnp{$chr."_".$pos} = 1;

	$prog ++ ;

	if ( ($prog % 200000) == 0) {
		print(STDERR "chr: $chr\n");
		print(STDERR "pos: $pos\n");
		print(STDERR "$prog / 152,300,679\n");
	}
}

close(fh);

my @files = <./mpileup.filt/*.bed.no_simple_repeats.no_homopolymers.no_SJ>;

foreach my $file (@files) {
	print("$file\n");

	open(fh, $file);
	open(fh_out, ">" . $file . ".no_dbSNP");

	while (my $line = <fh>) {
		chomp($line);

		$line =~ /^(\S+)\t\S+\t(\S+)/;
		my ($chr, $pos) = ($1,$2);
		$chr =~ s/chr//;

		if ($dbsnp{$chr."_".$pos} eq "") {
			print(fh_out "$line\n");
		}
	}

	close(fh);
	close(fh_out);
}
