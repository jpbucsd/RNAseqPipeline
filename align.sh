#!/bin/bash
#arguement 1: directory with fastq files
#arguement 2: reference genome file (ends in .fna)
for filename in $1/*.fq.gz do
	tar xzvf ${filename}
done

for filename in $1/*fastq.gz do
	tar xzvf ${filename}
done

for filename in $1/*.fq do
	bwa mem $2 filename | samtools view -bS > ${filename}.bam
done

for filename in $1/*.fastq do
	bwa mem $2 filename | samtools view -bS > ${filename}.bam
done

#a reference genome index may have to be created from a reference genome fasta and GTF files. this can be done with bwa index. see lab 1 from cse 185

#the next step will be to use the tool rsem, which will output a customized txt file containing the expression of each gene
#maybe a custom tool for this part based on what dr. sun wants



