#!/usr/bin/perl

use strict 'vars';

# Record the strand of genomic features
my %genomic_features_strand;

open(fh, "./ref/rm.S.bed") or die("Cannot open ./ref/rm.S.bed");

while (my $line = <fh>) {
	chomp($line);

	$line =~ /^(\S+)\t(\S+)\t(\S+)\t(\S+)\t\S+\t(\S+)/;
	my ($chr, $start, $stop, $name, $strand) = ($1,$2,$3,$4,$5);

	$genomic_features_strand{$chr."_".$start."_".$stop} = $strand;
}

close(fh);

open(fh, "./ref/gencode.v24.chr_patch_hapl_scaff.annotation.S.gtf") or die("Cannot open ./ref/gencode.v24.chr_patch_hapl_scaff.annotation.S.gtf");

while (my $line = <fh>) {
	chomp($line);

	next if ($line !~ /\tgene\t/);

	$line =~ /^(\S+)\t\S+\t\S+\t(\S+)\t(\S+)\t\S+\t(.)\t/;
	my ($chr, $start, $stop, $strand) = ($1,$2,$3,$4);

	$line =~ /gene_id \"(\S+?)\"\;/;
	my $gene = $1;

	$start -- ;

	$genomic_features_strand{$chr."_".$start."_".$stop} = $strand;
}

close(fh);

my @files = <./calls/*.intersected>;

foreach my $file (@files) {
	print(STDERR "$file\n");
	my %info;

	open(fh, $file);

	while (my $line = <fh>) {
		chomp($line);

		$line =~ /^(\S+)\t\S+\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)/;
		my ($variant_chr, $variant_pos, $feature_chr, $feature_start, $feature_stop, $feature_name) = ($1,$2,$3,$4,$5,$6);
		my $ref = $variant_chr . "_" . $variant_pos;

		$info{$ref}{$feature_name} = $feature_chr."_".$feature_start."_".$feature_stop;
	}

	close(fh);

	$file =~ /^.+\/(\S+?)\./;

	my $sample = $1;
	my $calls = "./calls/" . $sample . ".calls";
	
	open(fh, $calls) or die("Cannot open $calls");
	open(fh_out, ">" . $calls . ".annotated");

	while (my $line = <fh>) {
		chomp($line);

		$line =~ /^(\S+) (\S+) /;
		my ($chr, $pos) = ($1,$2);

		$line =~ /(\S+)$/;
		my $strand = $1;
		my $r = $info{$chr."_".$pos};
		my $all_overlapping_features;

		foreach my $tr (keys(%$r)) {
			$info{$chr."_".$pos}{$tr} =~ /^(\S+)_(\S+)_(\S+)/;
			my ($f_chr, $f_start, $f_stop) = ($1,$2,$3);
			my $f_strand = $genomic_features_strand{$f_chr."_".$f_start."_".$f_stop};

			$all_overlapping_features .= $tr . "_" . $f_strand . ";";
		}

		chop($all_overlapping_features);

		$all_overlapping_features = "NA" if ($all_overlapping_features eq "");

		print(fh_out "$line $all_overlapping_features\n");
	}

	close(fh);
	close(fh_out);
}