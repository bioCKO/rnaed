#!/usr/bin/perl

use strict 'vars';

open(fh, @ARGV[0]) or die("Cannot open @ARGV[0]");

while (my $line = <fh>) {
	chomp($line);

	$line =~ /^\S+\t\S+\t\S+\t\S+\t(\S+)\t/;
	my $foo = $1;

	$foo =~ s/[\-\+][0-9]+[ATGC]+//gi;
	$foo =~ s/[0-9\*\;\=\-\:\\$\~\^\<\>\,\.]//g;

	if (length($foo) >= 2 && $line !~ /random/ && $line !~ /chrUn/ && $line !~ /chrM/) {
		print("$line\n");
	}
}

close(fh);
