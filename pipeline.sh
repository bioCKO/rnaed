#!/bin/bash

echo RNA editing pipeline
echo Downloading the human genome GRCh38
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

echo Uncompressing.
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz