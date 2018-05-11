# Analysis pipeline for Spel *et al.* (2018)

## Summary
In this repository you will find all the command, scripts and reference files that
were used for the analysis of the two way phenotype screen in the Spel *et_al* (2018)
paper. The output of the analysis is a csv file with 

## Requirements
* A computer running some flavor of Linux with sufficient RAM, storage and 
computational power to analyze NGS data. The scripts are not tested (and likely
incompatible to due a different implementation of awk) on Mac OS X.
* Bowtie1
* R
** Bioconductor, DESeq2 and IRanges need to be installed (Bioconductor 3.5, IRanges 2.10.3 and DESeq2 1.16.1 were used)
* Python3
** Pandas and Numpy are required. Pandas 0.20.3 and Numpy 1.13.1 were used

## Walkthrough
1. Extract guides from sequence reads
2. Align extracted guides using Bowtie
3. Counting of aligned guides
4. Running DESeq2
5. Postprocessing

## 1. Extract guides from sequence reads

Run the following command for each of the 4 fastq.gz files comprising the raw sequence reads:
```bash
zcat reads_sample_x.fastq.gz | awk '{if(NR%4==2) print $0}' | sed -n -e 's/.*ACCG\([A-Z]\{19,21\}\)GTT.*/\1/p'  | awk '{print ">\n"$0}' > extracted_reads_sample_x.fasta
```
Adjust the two references 'reads_sample_x' to your actual filename.

This command does the following:
- decompresses the gz file
- extract every 4th line containing of the sequence read
- extracts any sequence of 19 - 21 bases flanked that is flanked by the ACCG and GTT sequences
- adds a ">" and new line to comply to the fasta file format

## 2. Align extracted guides using Bowtie

Again, for each of the 4 fasta you have, run the following commands:
```bash
bowtie -p 30 -v 1 -m 1 -f bowtie_1_ref/ref extracted_reads_sample_x.fasta aligned_reads_sample_x.sam --max supressed_reads_sample_x.fasta
bowtie -p 30 -v 0 -m 1 -f bowtie_1_ref/ref supressed_reads_sample_x.fasta aligned_recovered_reads_sample_x.sam
```
Decomposition of the command:
- -p: number of core
- -v: number of total mismatches
- -m: number of sequences that a reads is allowed to align against (1 for unique unambiguous alignment)
- --max: write reads that did not pass the -m threshold to a fasta file

So what the command does is that it executes Bowtie and tries to align the extracted guides to the
reference library. The first time it runs up to one mismatch is allowed. However, a read is never allowed
to align to more than one guide in the library. In case a read aligns to multiple guides in the library 
the reads is written to a new fasta file (called supressed_reads_sample_x.fasta). Bowtie is than executed
a seconds time to try to align these reads unambiguously to the library by not allowing any mismatches.

## 3. Counting of aligned guides

To determine the frequency of each of the reads, run the following command for each of the 4 samples:

```bash
cat aligned_reads_sample_x.sam aligned_recovered_reads_sample_x.sam | awk '{print $3}' | sort | uniq -c | awk '{print $1"\t"$2}'> sample_x.counts
```

## 4. Running DESeq2

First, either open the DESeq2_analysis.R script in your favorite text editor or directly in R or Rstudio. Next, alter 
the lines 19 till 22 to match the names of the 4 count files you created in the previous step. Now run the script inter-
actively using R or Rstudio or directly using Rscript.

```r
Rscript DESeq2_analysis.R
```

## 5. Postprocessing

The postprocessing basically consists of two parts. The first part deals with the NonTargetingControlGuides and assigns
them random but equally divided over 167 hypothetical NonTargetingControlGenes. The second part of the script is the hit
selection. Genes are selected as hit when at least 3 guides have an FDR corrected p-value < 0.05 and their 

```python
Python3 post_processing.py
```


