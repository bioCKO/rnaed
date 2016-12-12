# About
This is an RNA editing pipeline that calls RNA editing events from RNA-seq data. It uses the human genome as reference, but can easily be modified to use any reference genome. The aim of this repository is to provide the code used in our corresponding article.

# How to use
1. Clone the repository:

   ```
   git clone https://github.com/oscar-franzen/rnaed/
   ```
2. Install the necessary software components:
  * [GNU Wget](https://www.gnu.org/software/wget/)
  * [STAR v. 2.5.1b](https://github.com/alexdobin/STAR)
  * samtools v. 1.3
  * java v. 1.8.0
  * picard v. 1.112
  * bamtools
  * bedtools
  * GNU Awk
  * BioPerl (Bio::SeqIO)
  * Perl 5
  * gmap-2016-05-25
  * NCBI liftOver
  * GNU parallel
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

6. Put your *.fastq files in a directory called `fastq'.
7. Run:

   ```bash
   chmod +x pipeline.sh *.pl
   ./pipeline.sh
   ```

# Publication
* Manuscript in preparation.

# Contact
* Oscar Franzén, <p.oscar.franzen@gmail.com>