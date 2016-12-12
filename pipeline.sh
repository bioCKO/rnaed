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

gunzip ./ref/gencode.v24.chr_patch_hapl_scaff.annotation.S.gtf.gz

STAR --runMode genomeGenerate --genomeDir genome.index --genomeFastaFiles GCA_000001405.15_GRCh38_no_alt_analysis_set.fna --runThreadN $THREADS --sjdbGTFfile ./ref/gencode.v24.chr_patch_hapl_scaff.annotation.S.gtf

for file in ./fastq/*.fastq
do
	prefix=${file##*/}
	prefix=${prefix%.*}

	STAR --genomeLoad NoSharedMemory --genomeDir ./genome.index --readFilesIn $file --readFilesCommand cat --runThreadN $THREADS --outStd SAM --outSAMmode Full --outFilterMultimapNmax 1 --outFilterMismatchNoverLmax 0.1 > ${prefix}.sam
done