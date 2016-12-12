#!/usr/bin/perl

use strict 'vars';

die("Usage: <GSNAP bed> <STAR bed>\n") if (@ARGV != 2);

my %list;

my $bed_gsnap = @ARGV[0];
my $bed_star = @ARGV[1];

open(fh, $bed_gsnap) or die("Cannot open $bed_gsnap");

while (my $line = <fh>) {
	chomp($line);

	$line =~ /^(\S+)\t(\S+)\t(\S+)\t(\S+)\t\S+\t(\S+)/;
	my ($chr, $start, $stop, $id, $strand) = ($1,$2,$3,$4,$5);

	$list{$id}{strand} = $strand;
	$list{$id}{chr} = $chr;
	$list{$id}{start} = $start;
	$list{$id}{stop} = $stop;
}

close(fh);

open(fh, $bed_star) or die("Cannot open $bed_star");

while (my $line = <fh>) {
	chomp($line);

	$line =~ /^(\S+)\t(\S+)\t(\S+)\t(\S+)\t\S+\t(\S+)/;
	my ($chr, $start, $stop, $id, $strand) = ($1,$2,$3,$4,$5);

	# Alignments have same strand and same chromosome
	if ($strand eq $list{$id}{strand} && $chr eq $list{$id}{chr}) {
		my $a_start = $start;
		my $a_stop = $stop;

		my $b_start = $list{$id}{start};
		my $b_stop = $list{$id}{stop};

		# Check if the alignments overlap
		if ($a_stop >= $b_start && $a_start <= $b_stop) {
			my @q;
			push(@q, { start => $a_start, stop => $a_stop });
			push(@q, { start => $b_start, stop => $b_stop });

			@q = sort { $a->{start} <=> $b->{start} } @q;

			my $combined_length;

			if (@q[0]->{start} == @q[1]->{start}) {
				@q = sort { $a->{stop} <=> $b->{stop} } @q;
			}

			$combined_length = @q[1]->{stop} - @q[0]->{start};
			my $overlap = @q[0]->{stop} - @q[1]->{start};

			# foreach my $item (@q) {
			# 	print($item->{start} . "<\n");
			# }

			my $percent_overlap = $overlap / $combined_length;

			if ($percent_overlap > 0.99) {
				print("$id\n");
			}
		}
	}
}

close(fh);