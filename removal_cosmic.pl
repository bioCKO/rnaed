#!/usr/bin/perl

my @files = <./COSMIC/Cosmic*C>;

my %remove;

foreach my $file (@files) {
	open(fh, $file);

	while (my $line = <fh>) {
		chomp($line);

		$line =~ /^(\S+)\t(\S+)/;
		my ($chr, $pos) = ($1,$2);

		$chr = "chr" . $chr;

		$remove{$chr."_".$pos} = 1;
	}

	close(fh);
}

my @files = <./calls.prefinal/*.no_wellderly>;

foreach my $file (@files) {
	print("$file\n");

	open(fh, $file);
	open(fh_out, ">" . $file . ".no_COSMIC");

	while (my $line = <fh>) {
		chomp($line);

		$line =~ /^(\S+) (\S+)/;
		my ($chr, $pos) = ($1,$2);

		if ($remove{$chr."_".$pos} eq "") {
			print(fh_out "$line\n");
		}
	}

	close(fh);
	close(fh_out);
}