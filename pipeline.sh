#!/bin/bash

############################################################
# Assumes the following software versions are installed:   #
# GNU Wget                                                 #
# STAR v. 2.5.1b                                           #
# samtools                                                 #
# java v. 1.8.0                                            #
# picard v. 1.112                                          #
# bamtools                                                 #
############################################################

THREADS=3
MEM=4G

echo RNA editing pipeline
echo Downloading the human genome GRCh38
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

echo Uncompressing.
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

echo Indexing genome
mkdir genome.index

gunzip ./ref/gencode.v24.chr_patch_hapl_scaff.annotation.S.gtf.gz

STAR --runMode genomeGenerate --genomeDir genome.index --genomeFastaFiles GCA_000001405.15_GRCh38_no_alt_analysis_set.fna --runThreadN ${THREADS} --sjdbGTFfile ./ref/gencode.v24.chr_patch_hapl_scaff.annotation.S.gtf

mkdir mapped

for file in ./fastq/*.fastq
do
	prefix=${file##*/}
	prefix=${prefix%.*}

	STAR --genomeLoad NoSharedMemory --genomeDir ./genome.index --readFilesIn $file --readFilesCommand cat --runThreadN ${THREADS} --outStd SAM --outSAMmode Full --outFilterMultimapNmax 1 --outFilterMismatchNoverLmax 0.1 > ./mapped/${prefix}.sam
done

mkdir split

for file in ./mapped/*.sam
do
	samtools view -T GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -bS ${file} > ${file%.*}.bam
	samtools sort -m ${mem} ${file%.*}.bam ${file%.*}.S.bam
	samtools index ${file%.*}.S.bam
	samtools rmdup -s ${file%.*}.S.bam ${file%.*}.rmdups.bam

	java -Xms${mem} -Xmx4{mem} -XX:+UseSerialGC -jar AddOrReplaceReadGroups.jar I=${file%.*}.rmdups.bam O=${file%.*}.rg.bam SO=coordinate RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample

	samtools index ${file%.*}.rg.bam

	# split bam files by chromosome
	bamtools split -in ${file%.*}.rg.bam -reference

	# Delete sam records containing sequences mapped to contigs, etc.
	find . -type f -name "*EBV*" | xargs -I {} rm {}
	find . -type f -name "*chrUn*" | xargs -I {} rm {}
	find . -type f -name "*_random*" | xargs -I {} rm {}
	find . -type f -name "*chrM.bam" | xargs -I {} rm {}

	mv -v ./mapped/*REF* ./split
done

mkdir mpileup

for file in ./split/*.bam
do
	samtools mpileup -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna ${file} > ./mpileup/${file##*/}.mpileup.out
done

for prefix in `ls mpileup | cut -d '.' -f1 | sort -u`
do
	./filter.pl $prefix
done