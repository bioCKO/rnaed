#!/usr/bin/perl
# OF; May 26, 2016
# Example:
# ./get_RNA_editing_variants.pl ./gsnap.aln/AOR_507.gsnap.F.S.bam.mpileup.F > AOR_507.calls
# cat AOR_317.calls | awk '{lines[toupper($3"_"$4)]++} END {for (i in lines) {print i,lines[i]}}' | sort -k 2,2nr

use strict 'vars';

#my %variants;
my %total_reads;

die("Usage: <file with samtools pileup>\n") if (@ARGV != 1);

#my $args_bed_file = @ARGV[0];
my $args_samtools_mpileup = @ARGV[0];

# Sample name
$args_samtools_mpileup =~ /^.+\/(\S+?)\./;
my $sample = uc($1);

# Information about mapping bias so that I know which samples are strand specific and what the direction is
my %seq_direction;
my %read_lengths;

open(fh, "mapping_bias.final.txt") or die("Cannot open mapping_bias.final.txt");

<fh>;

while (my $line = <fh>) {
	chomp($line);

	$line =~ /^(\S+) (\S+) (\S+) (\S+)/;
	my ($sample, $sense, $antisense, $rl) = (uc($1),$2,$3,$4);
	my $ratio = $sense / $antisense;

	$total_reads{$sample} += $sense;
	$total_reads{$sample} += $antisense;

	$read_lengths{$sample} = $rl;

	if ($ratio > 0.90 && $ratio < 1.10) {
		$seq_direction{$sample} = "non_strand_specific";
	}
	else {
		if ($sense > $antisense) {
			$seq_direction{$sample} = "sense";
		}
		else {
			$seq_direction{$sample} = "antisense";
		}
	}
}

close(fh);

my $read_length = $read_lengths{$sample};

if ($read_length eq "") {
	die("Read length not found.\n");
}

# Ignore the sample if it has less than this number of reads
if ($total_reads{$sample} < 1000000) {
	die("Fewer than 1,000,000 reads (" . $total_reads{$sample} . "). Exiting.\n");
}

# Target positions
my $targets = "./mpileup.filt/$sample" . ".mpileup.F.bed.no_simple_repeats.no_homopolymers.no_SJ.no_dbSNP";
my %target_positions;

open(fh, $targets) or die("Cannot open $targets");

while (my $line = <fh>) {
	chomp($line);
	$line =~ /^(\S+)\t\S+\t(\S+)/;
	my ($chr, $pos) = ($1,$2);

	$target_positions{$chr."_".$pos} = 1;
}

close(fh);

# Parse the mpileup
open(fh, $args_samtools_mpileup) or die("Cannot open $args_samtools_mpileup");

while (my $line = <fh>) {
	chomp($line);

	$line =~ /^(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)/;
	my ($chr, $pos, $ref_base, $depth, $pileup, $quality, $read_positions) = ($1,$2,$3,$4,$5,$6,$7);
	next if ($target_positions{$chr."_".$pos} eq "");

	my $ref = $chr."_".$pos;
	my @rp = split(/,/,$read_positions);

	# Ignore positions with indels.
	next if ($pileup =~ /[\-\+][0-9]+[ATGC]+/gi);

	$pileup =~ s/[\$]//g; # Remove end marks
	$pileup =~ s/\^.//g; # Remove start mark.

	my %mismatches;
	my %read_pos;
	my $total_depth = 0;

	# Parse the pilup
	for (my $i=0; $i<length($pileup); $i++) {
		my $char = substr($pileup, $i, 1);
		my $quality_score = substr($quality, $i, 1);
		my $phred_score = ord($quality_score) - 33;
		my $position_in_read = @rp[$i];

		if ($char =~ /([atgc])/i) {
			my $q = $1;

			if ($phred_score >= 20) {

				# Ignore this variant if it's within 10 bp of the 5' start of the read
				# This is because there is 5' positional nucleotide bias in all samples:
				# /Users/rand/Bioinformatik/Proj/rna_editing/nucleotide_bias
				# http://127.0.0.1/?q=rna_editing/node/1376

				my $accept_position = 0;

				# Read is mapped in the forward direction
				if ($char =~ /[ATGC]/) {
					if ($position_in_read > 6) {
						$accept_position = 1;
						$read_pos{$position_in_read} ++ ;
					}
				}
				else {
					my $offset = $read_length - $position_in_read + 1;

					if ($offset > 6) {
						$accept_position = 1;
						$read_pos{$offset} ++ ;
					}
				}

				if ($accept_position == 1) {
					if ($seq_direction{$sample} eq "sense" || $seq_direction{$sample} eq "antisense") {
						$mismatches{$char} ++ ;
					}
					else {
						# Unstranded sequencing
						$mismatches{uc($char)} ++ ;
					}

					$total_depth ++ ;
				}
			}
		}
		elsif ($char =~ /([\.\,])/) {
			if ($phred_score >= 20) {
				$total_depth ++ ;
			}
		}
	}

	# If more than one base was different, I ignore this position
	next if (keys(%mismatches) > 1);
	next if (keys(%read_pos) == 0);

	my $edited_base = (keys(%mismatches))[0]; # the non-genomic base
	my $edited_depth = $mismatches{$edited_base};

	if ($ref_base eq "") {
		die("Ref base is empty at $chr $pos\n");
	}

	# Ignore IUPAC nucleotides
	next if ($ref_base =~ /[RYSWKMBNDHV]/);

	# Does the signal come from the forward or reverse strand?
	# This info is only known is strand-specific libraries.
	my $strand_info = "NA";

	if ($seq_direction{$sample} eq "sense") {
		# Reads supporting transcription from the reverse strand
		if ($edited_base =~ /[atgc]/) {
			$strand_info = "-";

			if ($edited_base eq "a") {
				$edited_base = "t";
			}
			elsif ($edited_base eq "t") {
				$edited_base = "a";
			}
			elsif ($edited_base eq "g") {
				$edited_base = "c";
			}
			elsif ($edited_base eq "c") {
				$edited_base = "g";
			}

			if ($ref_base eq "A") {
				$ref_base = "T";
			}
			elsif ($ref_base eq "T") {
				$ref_base = "A";
			}
			elsif ($ref_base eq "G") {
				$ref_base = "C";
			}
			elsif ($ref_base eq "C") {
				$ref_base = "G";
			}
		}
		else {
			$strand_info = "+";
		}
	}
	elsif ($seq_direction{$sample} eq "antisense") {
		if ($edited_base =~ /[ATGC]/) {
			$strand_info = "-";

			if ($edited_base eq "A") {
				$edited_base = "T";
			}
			elsif ($edited_base eq "T") {
				$edited_base = "A";
			}
			elsif ($edited_base eq "G") {
				$edited_base = "C";
			}
			elsif ($edited_base eq "C") {
				$edited_base = "G";
			}

			if ($ref_base eq "A") {
				$ref_base = "T";
			}
			elsif ($ref_base eq "T") {
				$ref_base = "A";
			}
			elsif ($ref_base eq "G") {
				$ref_base = "C";
			}
			elsif ($ref_base eq "C") {
				$ref_base = "G";
			}
		}
		else {
			$strand_info = "+";
		}
	}
	# Not strand specific sequencing. Report the variant in relation to the forward strand.
	else {
		# if ($feature_strand eq "-") {
		# 	if ($ref_base eq "A") {
		# 		$ref_base = "T";
		# 	}
		# 	elsif ($ref_base eq "T") {
		# 		$ref_base = "A";
		# 	}
		# 	elsif ($ref_base eq "G") {
		# 		$ref_base = "C";
		# 	}
		# 	elsif ($ref_base eq "C") {
		# 		$ref_base = "G";
		# 	}

		# 	if (uc($edited_base) eq "A") {
		# 		$edited_base = "T";
		# 	}
		# 	elsif (uc($edited_base) eq "T") {
		# 		$edited_base = "A";
		# 	}
		# 	elsif (uc($edited_base) eq "G") {
		# 		$edited_base = "C";
		# 	}
		# 	elsif (uc($edited_base) eq "C") {
		# 		$edited_base = "G";
		# 	}
		# }
	}

	if ($edited_depth >= 2) {
		# If the mismatch has the same position in all reads, it's suspicious.
		if (keys(%read_pos) > 1) {
			my $positions;
			foreach my $str (keys(%read_pos)) {
				$positions .= $str.":".$read_pos{$str} . ";";
			}

			chop($positions);

			print($chr . " " . $pos . " " . $edited_depth . " " . $ref_base . " " . $edited_base . " " . $pileup . " $total_depth $positions $strand_info\n");
		}
	}
}

close(fh);