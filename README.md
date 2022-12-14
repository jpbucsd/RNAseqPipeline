# RNAseqPipeline
RNA sequencing pipeline for sun lab

# Usage
## RNAseq.sh
RNAseq.sh is the main script in this pipeline that executes all of the subscripts. This pipeline is hosted in the SDSC TSCC super computer, where a graphics user interface on the users personal computer sends tasks to the super computer that execute the scripts in this repository. Usage of the commands in this repository that are automated by the GUI are as follows:

##### bash RNAseq.sh -f <(location of fastq directory in TSCC)
##### &emsp;-o (output directory in TSCC)
#####    &emsp;-d (flag indicating differential expression analysis)
#####    &emsp;-PADJ (padj value)
#####    &emsp;-log10 (log10 value to consider a significant result)
#####    &emsp;-A (log10 value to highlight the gene in output when no gene list is available)
#####    &emsp;-s (path to slr file, a custom file format that contains the necessary information to determine grouping and comparisons fastq files in different combinations)
#####    &emsp;-G (path to txt file containing a list of genes to highlight in output charts)
#####    &emsp;-a (flag indicating alignment must be done)
#####    &emsp;-Q (flag to produce a quality report of fastq files)
#####    &emsp;-T (flag to trim fastq files of low quality)


## DifferentialExpression.R
This R script performs differential expression analysis with DESEQ2. It produces a .csv file containing the log2Fold chain for each gene between the sample sets, and it produces a .csv file containing Rlog normalized gene counts for each gene among each sample.

##### Rscript Rtest.R -1 (name of first sample set) set1file1.genes.results set1file2.genes.results ... set1fileN.genes.results 
#####    &emsp;-2 (name of second sample set) set2file1.genes.results set2file2.genes.results ... set2fileN.genes.results 
#####    &emsp;-d <path/to/directory>

## PCA.py
This python script performs principled component analysis on each sample, resulting in a chart where each point represents one sample. This can be used to confirm the efficacy of the experiment by ensuring that experimental samples cluster together and control samples cluster together.

##### python PCA.py -f <file1.genes.results> <file2.genes.results> ... <fileN.genes.results> 
#####    &emsp;--oDir <output/directory/of/PCA/chart> 
#####    &emsp;--numComps <number of comparissons for PCA, default is sample number. There must be at least 2 files!>


## RSEM.sh
A deprecated file, its functionality was incorporated directly into RNAseq.sh

## align.sh
A deprecated file, its functionality was incorporated directly into RNAseq.sh

## aligner_tools.sh
A deprecated file, its functionality was incorporated directly into RNAseq.sh

## heatmap.py
This script produces a heatmap using the Rlog values produced by differential expression analysis.
Usage of this script has not yet been implemented into the graphics user interface.
##### python heatmap.py --g genelist.txt
#####    &emsp;-f normalizedCounts.csv
#####    &emsp;-c (toggles clustering on)

## results.py
This script produces lists of up and down regulated genes, a background list, and volcano plots based on the log2Fold change calculated with differential expression analysis.
##### results.py -f file1 ... fileN
#####    &emsp;--padjusted or --padj - adjusted pvalue, a padj of 0.5 implies 50% of significant results are false positives. results with Padj above 50% are filtered out
#####    &emsp;--log10 or --psig - -log10(pvalue) to determine which results are significant
#####    &emsp;--Llog10 to determine which results are worth naming, a log 10 value as threshold for which values to name
#####    &emsp;--odir output directory


## splicingvsexpression.py Usage
This script is not a part of the pipeline, and is just an additional tool. It has not been configured for reuse and currently takes a log2Fold change .csv and a .csv containing alternative splicing counts and produces a chart comparing the two. 
##### python splicingvsexpression.py

