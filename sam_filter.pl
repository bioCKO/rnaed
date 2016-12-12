#!/usr/bin/perl

use strict 'vars';

open(fh, @ARGV[0]) or die("Cannot open @ARGV[0]");

while (my $line = <fh>) {
	chomp($line);

	if ($line =~ /^\@/) {
		print("$line\n");
	}
	else {

		$line =~ /^\S+\t\S+\t\S+\t\S+\t(\S+)\t/;
		my $maq = $1;

		if ($maq > 20) {
			$line =~ /NM:i:(\d+)/;
			my $edit_dist = $1;

			if ($edit_dist < 5) {
				print("$line\n");
			}
		}
	}
}

close(fh);