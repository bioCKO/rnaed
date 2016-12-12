#!/usr/bin/perl

use strict 'vars';

my $prefix = @ARGV[0];

my @files = <./mpileup/*$prefix\.*>;

open(fh_out, ">./mpileup.merged/$prefix" . ".mpileup");

foreach my $file (@files) {
	print(STDERR "$file\n");

	open(fh, $file);

	while (my $line = <fh>) {
		chomp($line);
		$line =~ /\t(\S+)\t\S+$/;
 		my $depth = $1;
 		$depth =~ s/[0-9\*\;\=\-\:\$\~\^\<\>\,\.]//g;

 		if (length($depth) >= 2) {
 			print(fh_out "$line\n");
 		}
	}

	close(fh);
}

close(fh_out);