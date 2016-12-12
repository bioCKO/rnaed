#!/usr/bin/perl

use strict 'vars';

my %list;

open(fh, "./ref/gencode.v24.chr_patch_hapl_scaff.annotation.S.gtf") or die("Cannot open ./ref/gencode.v24.chr_patch_hapl_scaff.annotation.S.gtf");

while (my $line = <fh>) {
	chomp($line);
	next if ($line =~ /^#/);

	if ($line =~ /\texon\t/ && $line =~ /transcript_id/) {
		$line =~ /transcript_id \"(\S+)\"/;
		
		my $tr = $1;

		$line =~ /^(\S+)\t\S+\t\S+\t(\S+)\t(\S+)\t.\t(.)\t/;
		my ($chr, $start, $end, $strand) = ($1,$2,$3,$4);

		push(@$tr, { chr => $chr, start => $start, end => $end, strand => $strand });

		$list{$tr} = 1;
	}
}

close(fh);

foreach my $tr (keys(%list)) {
	if (@$tr > 1) {
		for (my $i=0; $i<@$tr - 1; $i++) {
			my $exon1 = @$tr[$i];
			my $exon2 = @$tr[$i+1];

			if ($exon1->{chr} eq $exon2->{chr}) {
				my $first_start = $exon1->{end} - 5;
				my $first_end = $exon1->{end} + 5;

				my $second_start = $exon2->{start} - 5;
				my $second_end = $exon2->{start} + 5;

				$first_start -- ;
				$second_start -- ;

				if ($exon1->{chr} =~ /chr[0-9XY]+/) {
					print($exon1->{chr} . "\t$first_start\t$first_end\n");
					print($exon1->{chr} . "\t$second_start\t$second_end\n");
				}
			}
		}
	}
}