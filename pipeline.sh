#!/bin/bash

############################################################
# Assumes the following software versions are installed:   #
# GNU Wget                                                 #
# STAR v. 2.5.1b                                           #
# samtools v. 1.3                                          #
# java v. 1.8.0                                            #
# picard v. 1.112                                          #
# bamtools                                                 #
# bedtools                                                 #
# GNU Awk                                                  #
# BioPerl (Bio::SeqIO)                                     #
# Perl 5                                                   #
# gmap-2016-05-25                                          #
# NCBI liftOver                                            #
# GNU parallel                                             #
############################################################

THREADS=3
MEM=4G

echo RNA editing pipeline
echo Downloading the human genome GRCh38
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

echo Uncompressing
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
	samtools sort -m ${MEM} ${file%.*}.bam -o ${file%.*}.S.bam
	samtools index ${file%.*}.S.bam
	samtools rmdup -s ${file%.*}.S.bam ${file%.*}.rmdups.bam

	# AddOrReplaceReadGroups.jar is part of Picard Tools
	java -Xms${MEM} -Xmx${MEM} -XX:+UseSerialGC -jar AddOrReplaceReadGroups.jar I=${file%.*}.rmdups.bam O=${file%.*}.rg.bam SO=coordinate RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample

	samtools index ${file%.*}.rg.bam

	# split bam files by chromosome
	bamtools split -in ${file%.*}.rg.bam -reference

	# Delete sam records containing sequences mapped to contigs, etc.
	find mapped -type f -name "*EBV*" | xargs -I {} rm {}
	find mapped -type f -name "*chrUn*" | xargs -I {} rm {}
	find mapped -type f -name "*_random*" | xargs -I {} rm {}
	find mapped -type f -name "*chrM.bam" | xargs -I {} rm {}

	mv -v ./mapped/*REF* ./split
done

mkdir mpileup

(for file in ./split/*.bam
do
	echo "samtools mpileup -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna ${file} > ./mpileup/${file##*/}.mpileup.out"
done) > cmds.txt

parallel -j ${THREADS} < cmds.txt

mkdir mpileup.merged

(for prefix in `ls mpileup | cut -d '.' -f1 | sort -u`
do
	echo "./filter.pl $prefix"
done) > cmds.txt

parallel -j ${THREADS} < cmds.txt

mkdir mpileup.filt

(for file in ./mpileup.merged/*.mpileup
do
	p=${file##*/}
	echo "./remove_indels.pl $file > ./mpileup.filt/${p}.F"
done) > cmds.txt

parallel -j ${THREADS} < cmds.txt

gunzip ./ref/simple_repeats.plus_minus_one_bp.bed.gz

(for file in ./mpileup.filt/*.mpileup.F
do
	p=${file##*/}
	echp "awk '{print $1"\t"($2-1)"\t"$2}' ${file} > ./mpileup.filt/${p}.bed"
done) > cmds.txt

parallel -j ${THREADS} < cmds.txt

(for file in ./mpileup.filt/*.bed
do
	echo "bedtools intersect -v -a ${file} -b ./ref/simple_repeats.plus_minus_one_bp.bed > ${file}.no_simple_repeats"
done) > cmds.txt

parallel -j ${THREADS} < cmds.txt

./split_by_chr.pl GCA_000001405.15_GRCh38_no_alt_analysis_set.fna

find . -type f -name "*EBV*" | xargs -I {} rm {}
find . -type f -name "*chrUn*" | xargs -I {} rm {}
find . -type f -name "*_random*" | xargs -I {} rm {}
find . -type f -name "*chrM.bam" | xargs -I {} rm {}

(for chr_no in `seq 1 22` X Y
do
	chr=chr${chr_no}
	echo "./homopolymers.pl GCA_000001405.15_GRCh38_no_alt_analysis_set.fna_${chr} > ${chr}.hp"
done) > cmds.txt

parallel -j ${THREADS} < cmds.txt

cat *.hp > homopolymers_5bp-length.grch38

rm -v *.hp
rm -v GCA_000001405.15_GRCh38_no_alt_analysis_set.fna_*

(for file in ./mpileup.filt/*.no_simple_repeats
do
	echo "bedtools intersect -v -a ${file} -b homopolymers_5bp-length.grch38 > ${file}.no_homopolymers"
done) > cmds.txt

parallel -j ${THREADS} < cmds.txt

./make_splice_junction_bed_file.pl > SJ.bed

(for file in ./mpileup.filt/*.no_homopolymers
do
	echo "bedtools intersect -v -a ${file} -b SJ.bed | grep -v chrM | grep -v chrEBV > ${file}.no_SJ"
done) > cmds.txt

parallel -j ${THREADS} < cmds.txt

wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b147_GRCh38p2/VCF/All_20160527.vcf.gz
gunzip All_20160527.vcf.gz

wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b147_GRCh38p2/XML/*

(for file in ds_*.xml.gz
do
	echo "gunzip $file"
done) > cmds.txt

parallel -j ${THREADS} < cmds.txt

(for file in *.xml
do
	echo "./xml_parser.pl ${file} > ${file}.cDNA_snps"
done) > cmds.txt

parallel -j ${THREADS} < cmds.txt

cat *.cDNA_snps > dbSNP.cDNA_evidence_only

./prepare.pl > All_20160527.vcf.chr_and_pos
./remove_snps.pl

gmap_build -D . -d hg38.gsnap -k 15 GCA_000001405.15_GRCh38_no_alt_analysis_set.fna

mkdir read_lists

for file in ./mpileup.filt/*.no_dbSNP
do
	prefix=${file##*/}
	prefix=${prefix%.*}
	prefix=${prefix%.*}
	prefix=${prefix%.*}
	prefix=${prefix%.*}
	prefix=${prefix%.*}
	prefix=${prefix%.*}
	prefix=${prefix%.*}

	samtools view -L ${file} ./mapped/${prefix}.rg.bam | awk '{print $1}' | sort -u > ./read_lists/${prefix}.read_list
done

(for file in ./read_lists/*.read_list
do
	prefix=${file##*/}
	prefix=${prefix%.*}

	echo "./extract_fastq.pl ${file} > ${file}.fq"
done) > cmds.txt

parallel -j ${THREADS} < cmds.txt

mkdir gsnap.aln

for file in ./read_lists/*.fq
do
	prefix=${file##*/}
	prefix=${prefix%.*}
	prefix=${prefix%.*}

	gsnap -t ${THREADS} --trim-mismatch-score 0 -A sam -m 5 -i 2 -N 1 -D . -d hg38.gsnap ${file} > ./gsnap.aln/${prefix}.gsnap.sam

	./sam_filter.pl ./gsnap.aln/${prefix}.gsnap.sam > ./gsnap.aln/${prefix}.gsnap.F.sam
done

for file in ./gsnap.aln/*.F.sam
do
	prefix=${file##*/}
	prefix=${prefix%.*}
	prefix=${prefix%.*}
	prefix=${prefix%.*}

	samtools view -T GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -bS ${file} > ${file%.*}.bam
done

for file in ./gsnap.aln/*.gsnap.F.bam ./mapped/*.rg.bam
do
	prefix=${file##*/}
	prefix=${prefix%.*}

	bamToBed -i ${file} > ${file}.bed
done

for file in ./gsnap.aln/*.gsnap.F.bam
do
	prefix=${file##*/}
	prefix=${prefix%.*}
	prefix=${prefix%.*}
	prefix=${prefix%.*}

	./find_concordant.pl ./gsnap.aln/${prefix}.gsnap.F.bam.bed ./mapped/${prefix}.rg.bam.bed > ./gsnap.aln/${prefix}.gsnap.bam.targets
done

for file in ./gsnap.aln/*.targets
do
	prefix=${file##*/}
	prefix=${prefix%.*}
	prefix=${prefix%.*}
	prefix=${prefix%.*}

	# picard
	java -Xms${MEM} -Xmx${MEM} -XX:+UseSerialGC -jar FilterSamReads.jar INPUT=./gsnap.aln/${prefix}.gsnap.F.bam FILTER=includeReadList READ_LIST_FILE=${file} OUTPUT=./gsnap.aln/${prefix}.final.bam
done

for file in ./gsnap.aln/*.final.bam
do
	prefix=${file##*/}
	prefix=${prefix%.*}
	prefix=${prefix%.*}

	# the --output-BP options is not in older versions of samtools

	samtools sort -m ${MEM} ${file} -o ${file%.*}.S.bam

	samtools mpileup --output-BP -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna ${file%.*}.S.bam | perl -ne ' {
     $_ =~ /\t(\S+)\t\S+\t\S+$/;
     my $q = $1;
     $q =~ s/[0-9\*\;\=\-\:\\$\~\^\<\>\,\.]//g;

     if (length($q) >= 2) {
         print("$_");
     }
 } ' > ${file%.*}.S.bam.mpileup
done

mkdir calls

for file in ./gsnap.aln/*.mpileup
do
	prefix=${file##*/}
	prefix=${prefix%.*}
	prefix=${prefix%.*}
	prefix=${prefix%.*}
	prefix=${prefix%.*}

	./get_rna_editing_variants.pl ${file} > ./calls/${prefix}.calls
done

gunzip ./ref/rm.S.bed.gz

(awk '$3 == "gene" {print $1"\t"($4-1)"\t"$5"\t"$10}' ./ref/gencode.v24.chr_patch_hapl_scaff.annotation.S.gtf | sed 's/[\";]//g'
cat ./ref/rm.S.bed | grep -v rich | grep -v -P '\([ATGC]+\)n' | cut -d$'\t' -f1-4) > genes_and_repeats.bed

for file in ./calls/*.calls
do
	awk '{print $1"\t"($2-1)"\t"$2}' $file > ${file}.bed
done

(for file in ./calls/*.bed
do
	echo "bedtools intersect -wa -wb -a ${file} -b genes_and_repeats.bed > ${file}.intersected"
done) > cmds.txt

parallel -j ${THREADS} < cmds.txt

./annotate.pl

cat ./calls/*.annotated | cut -d ' ' -f1-2,4,5,10 | awk '{$3=toupper($3); $4=toupper($4); print $0}' | sort -u > feature_overlap.txt

mv -v calls calls.prefinal
mkdir calls

for file in ./calls.prefinal/*.annotated
do
 echo $file

 awk '{
  if ($1 != "chr6") {
    print($0)
  }
  else {
    if ($2 < 28510120 || $2 > 33480577) {
        print($0)
    }
  }
 }' $file > ${file}.no_MHC
done

./remove_snps_ucsc_146.pl

./remove_snps_ucsc_141.pl

wget ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/ExAC.r0.3.1.sites.vep.vcf.gz

gunzip ExAC.r0.3.1.sites.vep.vcf.gz

awk '$7 == "PASS" {
 chr="chr"$1
 start=$2 - 1
 stop=$2

 print(chr"\t"start"\t"stop)
}' ExAC.r0.3.1.sites.vep.vcf > ExAC.r0.3.1.sites.vep.vcf.bed

liftOver ExAC.r0.3.1.sites.vep.vcf.bed hg19ToHg38.over.chain ExAC.r0.3.1.sites.vep.vcf.GRCh38.bed unmapped.txt

./remove_snps_exac.pl

wget http://evs.gs.washington.edu/evs_bulk_data/ESP6500SI-V2-SSA137.GRCh38-liftover.snps_indels.txt.tar.gz

gunzip ESP6500SI-V2-SSA137.GRCh38-liftover.snps_indels.txt.tar.gz

tar -zxvf ESP6500SI-V2-SSA137.GRCh38-liftover.snps_indels.txt.tar.gz

(awk 'NR == 8' ESP6500SI-V2-SSA137.GRCh38-liftover.chr22.snps_indels.txt
for file in ESP6500SI-V2-SSA137.GRCh38-liftover.chr*
do
 cat $file | grep -v -P '^#'
done) > ESP6500SI-V2-SSA137.GRCh38-liftover.vcf

rm -v ESP6500SI-V2-SSA137.GRCh38-liftover.chr*

./remove_snps_nhlbi.pl

wget -r https://genomics.scripps.edu/browser/files/wellderly/vcf/

for file in *.gz
do
 echo $file
 zcat ${file} > ../scripps.wellderly.uncomp/${file%.*}
done

(for file in *.vcf
do
 awk '{
  if (length($4) == 1 && length($5) == 1) {
   chr="chr"$1
   start = $2 - 1
   stop = $2

   print(chr"\t"start"\t"stop)
 }
}' $file
done) > ../scripps.wellderly.bed

liftOver scripps.wellderly.bed hg19ToHg38.over.chain scripps.wellderly.GRCh38.bed unmapped.txt

./remove_snps_wellderly.pl

./removal_cosmic.pl

for file in ./calls.prefinal/*.no_COSMIC
do
 d=${file##*/}
 d=`echo $d | sed 's/\(.*\)\.calls.*/\1/'`

 ln -sfn $file ./calls/${d}.calls
done

echo All done.