# About
This is an RNA editing pipeline that calls RNA editing events from RNA-seq data.

# How to use
1. Install the necessary software components (listed in the top of pipeline.sh).
2. Prepare a file called "ucsc_146.snps.GRCh38.final", goto http://genome.ucsc.edu/cgi-bin/hgTables and select group "Variation" and track "All SNPs(146)" and Assembly "Dec. 2013(GRCh38/hg38)". In output format, choose "BED". Type in the filename ucsc_146.bed and press "Get output".  

   Then execute:  
   ```bash
   awk -F '\t' '{print $1"\t"$3}' ucsc_146.bed > ucsc_146.snps.GRCh38.final
   ```
3. Repeat step 2, but instead for version 141. The final file should be called ucsc_141.snps.GRCh38.final.

4. Register and download [COSMIC](http://cancer.sanger.ac.uk/cosmic). The license of this database is not compatible with unrestricted distribution. Only the VCF files are needed. Execute:

   ```
   for file in ./COSMIC/Cosmic*
   do
    cat ${file} | grep -v -E '^#' | awk -F '\t' '{print $1"\t"$2}' > ${file}.C
   done
   ```

5. Put your *.fastq files in a directory called `fastq'.
6. Run:

   ```bash
   chmod +x pipeline.sh *.pl
   ./pipeline.sh
   ```

# Publication
* Manuscript in preparation.

# Contact
* Oscar Franzén, <p.oscar.franzen@gmail.com>