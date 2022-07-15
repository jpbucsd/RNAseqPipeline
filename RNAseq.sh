#!/bin/bash

#important variables

fqDir=""
oDir=""
indexF=0
fqFlag=0
oFlag=0
align=1

#usage of command line arguments
#-h - help //todo
#-f fastqdirectory
#-o output directory
#-index (compute index) usage of index with no other arguments will cause only indexing to occur. without this flag indexing will not occur

#input command line arguments
for var in "$@"
do
	if [[ $fqFlag == 1 ]]
	then
		fqFlag=0
		fqDir="$var"
		align=2
		#next check for other flags before -
	elif [[ $oFlag == 1 ]]
	then
		oFlag=0
		oDir="$var"/
	elif [[ "$var" == "-"* ]]
	then
		if [[ "$var" == *"-i"* ]]
		then
			#calculate index!
			indexF=1
			if [[ "$align" == 1 ]]
			then
				#it is assumed we want to align, unless index is called. if index is called then we will set align to 0 
				#if align is 2, it means there have already been parameters for alignment
				align=0
			fi
		fi
		if [[ "$var" == *"-f"* ]]
		then
			#fast q files coming up
			fqFlag=1
			align=2
		fi
		if [[ "$var" == *"-o"* ]]
		then
			oFlag=1
		fi
	else
		#default choice to assume fqdir if fqdir is empty
		if [[ "$fqDir" == "" ]]
		then
			fqDir="$var"
			align=2
		fi
	fi
done

#this script produces an index file from the referencce genome (fasta format) and a GTF annotations file.
#the produced index file will next be used in the second step written in the script align.sh
#this script will only execute if the instruction to generate an index is included.

if [[ $indexF == 1 ]]
then
	echo "indexing"
	mkdir refGen
	cd refGen
	wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
	wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz
	gzip -d GCF_000001405.40_GRCh38.p14_genomic.fna.gz
	gzip -d GCF_000001405.40_GRCh38.p14_genomic.gtf.gz
        
	#find number of chromsomes for star to work without crashing
	#this has already been calculated and for now will be set to the calculated value of 22. if this number returns a value below 18 make it 18
	#chrN=$(cat GCF_000001405.40_GRCh38.p14_genomic.fna | grep chr | wc -l)
	#baseN=$(cat GCF_000001405.40_GRCh38.p14_genomic.fna | grep -v chr | wc -m)
	#BpC=$(expr $baseN / $chrN)
	#logV=$(echo "l($BpC)/l(2)" | bc -l)
	#factor=$(printf "%.0f\n" $logV)
	#echo $factor
	factor=22

		
	mkdir genome
	chmod -R 0777 genome
	#star requires 36GB of ram and 12 threads
	STAR --runMode genomeGenerate --genomeFastaFiles GCF_000001405.40_GRCh38.p14_genomic.fna --sjdbGTFfile GCF_000001405.40_GRCh38.p14_genomic.gtf --genomeDir genome --genomeChrBinNbits $factor --runThreadN 16
	
	rsem-prepare-reference --gtf GCF_000001405.40_GRCh38.p14_genomic.gtf --runThreadN 16 GCF_000001405.40_GRCh38.p14_genomic.fna reference
		
	#leave refGen to return to normal
	cd ../
	
fi

#at this point it is assumed that there is a folder containing indexed information at refGen/genome/

#the next step in the RNA-seq pipeline is alignment of .fq to the reference genome
#because the previous step was performed in STAR, in order to get BAM files this too must be performed with STAR.

if [[ $align != 0 ]]
then
	for filename in $fqDir/*.fq.gz 
	do
		echo "unzipping $filename"
		gzip -d ${filename}
	done

	for filename in $fqDir/*.fastq.gz 
	do
		echo "unzipping $filename"
		gzip -d ${filename}
	done
	
	mkdir ${oDir%/}
	chmod -R 0777 ${oDir%/}
	mkdir "${oDir%/}/$fqDir"
	chmod -R 0777 "${oDir%/}/$fqDir"

	for filename in $fqDir/*.fq 
	do
		echo "aligning $filename"
		#the following line is the originally intended bwa command with $2 being an index genome in fna format or something.
		#bwa mem $2 filename | samtools view -bS > ${filename}.bam
		STAR --genomeDir refGen/genome --readFilesIn ${filename} --outFileNamePrefix "$oDir${filename%.*}" --runThreadN 16
	done
	
	for filename in $fqDir/*.fastq 
	do
		echo "aligning $filename"
		STAR --genomeDir refGen/genome --readFilesIn ${filename} --outFileNamePrefix "$oDir${filename%.*}" --runThreadN 16
		#bwa mem $2 filename | samtools view -bS > ${filename}.bam
	done

	#the next step will be to use the tool rsem, which will output a customized txt file containing the expression of each gene
	#maybe a custom tool for this part based on what dr. sun wants

	#quantifying gene expression
	#cd ${oDir%/}
	

	for filename in $oDir/$fqDir/*.out.sam
	do
		echo "calculating expression of ${filename}"
		sample_fname="${filename##*/}"
		sample_name="${sample_fname%.*}"
		rsem-calculate-expression --num-threads 16 "${filename}" refGen/reference "${sample_name}"
	done
	
	#cd ../
	#at this point we will have gene.results files?? and they can be used in the R language with Deseq2
fi

