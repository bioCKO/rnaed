# About
This is an RNA editing pipeline that calls RNA editing events from RNA-seq data.

# How to use
1. Install the necessary software components (listed in the top of pipeline.sh).
2. Prepare a file called "ucsc_146.snps.GRCh38.final", goto http://genome.ucsc.edu/cgi-bin/hgTables and select group "Variation" and track "All SNPs(146)" and Assembly "Dec. 2013(GRCh38/hg38)". In output format, choose "BED". Type in the filename ucsc_146.bed and press "Get output". Then execute: ```bash
awk -F '\t' '{print $1"\t"$3}' ucsc_146.bed > ucsc_146.snps.GRCh38.final```
3. Put your *.fastq files in a directory called `fastq'.
4. Run ```bash
chmod +x pipeline.sh *.pl
./pipeline.sh
```