#!/usr/bin/perl

use strict 'vars';

die("Usage: <file>") if (@ARGV != 1);

my %list;

open(fh, @ARGV[0]) or die("Cannot open @ARGV[0]");

while (my $line = <fh>) {
	chomp($line);

	$list{$line} = 1;
}

close(fh);

my $sample = @ARGV[0];

$sample =~ /\.\/read_lists\/(\S+)\.read_list/;
$sample = $1;

my $fq = "./data.merged/$sample" . ".fastq";

open(fh, $fq) or die("Cannot open $fq");

while (my $line1 = <fh>) {
	my $line2 = <fh>;
	my $line3 = <fh>;
	my $line4 = <fh>;

	chomp($line1);
	chomp($line2);
	chomp($line3);
	chomp($line4);

	$line1 =~ /^\@(\S+)/;
	my $id = $1;

	if ($list{$id} ne "") {
		print("$line1\n");
		print("$line2\n");
		print("$line3\n");
		print("$line4\n");
	}
}

close(fh);