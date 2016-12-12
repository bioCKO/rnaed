#!/usr/bin/perl

use strict 'vars';

die("Usage: <xml file>") if (@ARGV != 1);

open(fh, @ARGV[0]) or die("Cannot open @ARGV[0]");

while (my $line = <fh>) {
	chomp($line);

	if ($line =~ /<Rs rsId=\"(\S+)\"/) {
		my $rs = $1;
		my @foo;
		my %molecule_types;
		my $validated = 0;

		while (my $line2 = <fh>) {
			chomp($line2);

			push(@foo, $line2);

			if ($line2 =~ /<\/Rs>/) {
				last;
			}
		}

		foreach my $item (@foo) {
			if ($item =~ /molType=\"(\S+)\"/) {
				my $ml = $1;
				$molecule_types{$ml} = 1;
			}

			if ($item =~ /Validation/ && $item =~ /true/) {
				$validated = 1;
			}
		}

		my $types;

		foreach my $type (keys(%molecule_types)) {
			$types .= $type . ",";
		}

		chop($types);

		if ($types eq "cDNA" && $validated == 0) {
			print("rs$rs $validated $types\n");
		}
	}
}

close(fh);