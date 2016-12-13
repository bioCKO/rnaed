# About
This is an RNA editing pipeline that calls RNA editing events from RNA-seq data. It uses the human genome as reference, but can easily be modified to use any reference genome. The aim of this repository is to provide the code used in our corresponding article.

# How to use
1. Clone the repository:

   ```
   git clone https://github.com/oscar-franzen/rnaed/
   ```
2. Install software dependencies. The tools should be in your ```$PATH```.
  * [GNU Wget](https://www.gnu.org/software/wget/)
  * [STAR v. 2.5.1b](https://github.com/alexdobin/STAR)
  * [samtools v. 1.3](http://www.htslib.org/download/)
  * [java v. 1.8.0](http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html)
  * [Picard v. 1.112](https://broadinstitute.github.io/picard/)
  * [bamtools](https://github.com/pezmaster31/bamtools)
  * [bedtools](http://bedtools.readthedocs.io/en/latest/)
  * GNU Awk
  * [BioPerl (Bio::SeqIO)](http://bioperl.org/INSTALL.html)
  * Perl 5.x
  * [gmap-2016-05-25](http://research-pub.gene.com/gmap/)
  * [NCBI liftOver](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver)
  * [GNU parallel](https://www.gnu.org/software/parallel/)

  Optional:
  * [htseq](http://www-huber.embl.de/HTSeq/doc/count.html)
3. Prepare a file called "ucsc_146.snps.GRCh38.final", goto http://genome.ucsc.edu/cgi-bin/hgTables and select group "Variation" and track "All SNPs(146)" and Assembly "Dec. 2013(GRCh38/hg38)". In output format, choose "BED". Type in the filename ucsc_146.bed and press "Get output".  

   Then execute:  
   ```bash
   awk -F '\t' '{print $1"\t"$3}' ucsc_146.bed > ucsc_146.snps.GRCh38.final
   ```
4. Repeat step 2, but instead for version 141. The final file should be called ucsc_141.snps.GRCh38.final.

5. Register and download [COSMIC](http://cancer.sanger.ac.uk/cosmic). The license of this database is not compatible with unrestricted distribution. Only the VCF files are needed. Put all files in a directory called `COSMIC`. Execute:

   ```
   for file in ./COSMIC/Cosmic*
   do
    cat ${file} | grep -v -E '^#' | awk -F '\t' '{print $1"\t"$2}' > ${file}.C
   done
   ```

6. Prepare a file named `mapping_bias.final.txt`. This file should have four columns:
  * the sample name (taken from the fastq file: `sample_name.fastq`)
  * number of reads mapping mapping in the forward direction on protein-coding exons
  * number of reads mapping in the reverse direction on protein-coding exons
  * sequence read length

   One way to retrieve the counts is to use the program htseq.
7. Put your *.fastq files in a directory called `fastq'.
8. Run:

   ```bash
   chmod +x pipeline.sh *.pl
   ./pipeline.sh
   ```

# Publication
* Manuscript in preparation.

# Contact
* Oscar Franzén, <p.oscar.franzen@gmail.com>