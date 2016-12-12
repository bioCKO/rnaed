#!/bin/bash

############################################################
# Assumes the following software versions are installed:   #
# GNU Wget                                                 #
# STAR v. 2.5.1b                                           #
############################################################

THREADS=3

echo RNA editing pipeline
echo Downloading the human genome GRCh38
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

echo Uncompressing.
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

echo Indexing genome
mkdir genome.index

STAR --runMode genomeGenerate --genomeDir genome.index --genomeFastaFiles GCA_000001405.15_GRCh38_no_alt_analysis_set.fna --runThreadN $THREADS --sjdbGTFfile /sc/orga/projects/STARNET/oscar/Machine_learning/ref/GENCODE_r24/gencode.v24.chr_patch_hapl_scaff.annotation.gtf