#!/bin/bash

#important variables

fqDir=""
oDir=""
indexF=0
fqFlag=0
oFlag=0
align=1
readLength=50
quantify=1
rFlag=0
pFlag=0
paired=0
pair1="_1"
pair2="_2"

#usage of command line arguments
#-h - help //todo
#-f fastqdirectory
#-o output directory for alignment

#-r readLength
#-p signals paired ends (should have the same name and end in _1 and _2, otherwise type what files will end in)

#-index (compute index) usage of index with no other arguments will cause only indexing to occur. without this flag indexing will not occur

#input command line arguments
for var in "$@"
do
	if [[ $fqFlag == 1 ]]
	then
		fqFlag=0
		fqDir="$var"
		align=2
		quantify=2
		#next check for other flags before -
	elif [[ $oFlag == 1 ]]
	then
		oFlag=0
		oDir="$var"/
	elif [[ $rFlag == 1 ]]
	then
		rFlag=0
		readLength=$var
	elif [[ "$var" == "-"* ]]
	then
		#some flags may have multiple phrases following them so the -option comes first in the else ifs. 
		#unflag these first
		if [[ $pFlag != 0 ]]
		then
			pFlag=0
		fi
		if [[ "$var" == *"-i"* ]]
		then
			#calculate index!
			indexF=1
			if [[ "$align" == 1 ]]
			then
				#it is assumed we want to align, unless index is called. if index is called then we will set align to 0 
				#if align is 2, it means there have already been parameters for alignment
				align=0
				quantify=0
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
		if [[ "$var" == *"-o"* ]]
		then
			oFlag=1
		fi
		if [[ "$var" == *"-r"* ]]
		then
			rFlag=1
		fi
		if [[ "$var" == *"-p"* ]]
		then
			pFlag=1
			paired=1
		fi
	elif [[ $pFlag == 1 ]]
	then
		pFlag=2
		pair1="$var"
	elif [[ $pFlag == 2 ]]
	then
		pFlag=0
		pair2="$var"
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

#check if genome has been indexed for the proper settings
if [[ $indexF == 0 ]]
then
	tempBool=0
	if [ -d "refGen/genome$readLength" ] 
	then
		echo "genome for length $readLength has already been indexed, moving to alignment"
    		tempBool=1
	fi
	if [[ $tempBool == 0 ]]
	then
		#must index
		echo "genome for length $readLength has not been indexed, moving to indexing"
		indexF=1
	fi
fi
if [[ $indexF == 1 ]]
then
	echo "indexing"
	
	if ! [ -d refGen ]
	then
		mkdir refGen
	fi
	cd refGen
	if ! [ -f GCF_000001405.40_GRCh38.p14_genomic.fna ]
	then
		if ! [ -f GCF_000001405.40_GRCh38.p14_genomic.fna.gz ]
		then
			wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
		fi
		gzip -d GCF_000001405.40_GRCh38.p14_genomic.fna.gz
	fi
	if ! [ -f GCF_000001405.40_GRCh38.p14_genomic.gtf ]
	then
		if ! [ -f GCF_000001405.40_GRCh38.p14_genomic.gtf.gz ]
		then
			wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz
		fi
		gzip -d GCF_000001405.40_GRCh38.p14_genomic.gtf.gz
	fi
	
	
	
        
	#find number of chromsomes for star to work without crashing
	#this has already been calculated and for now will be set to the calculated value of 22. if this number returns a value below 18 make it 18
	chrN=$(cat GCF_000001405.40_GRCh38.p14_genomic.fna | grep chr | wc -l)
	baseN=$(cat GCF_000001405.40_GRCh38.p14_genomic.fna | grep -v chr | wc -m)
	BpC=$(expr $baseN / $chrN)
	logV=$(echo "l($BpC)/l(2)" | bc -l)
	factor=$(printf "%.0f\n" $logV)
	#factor=22
	
	overhang=$(expr $readLength - 1 )
	
		
	mkdir "genome${readLength}"
	chmod -R 0777 "genome${readLength}"
	#star requires 36GB of ram and 12 threads
	STAR --runMode genomeGenerate --genomeFastaFiles GCF_000001405.40_GRCh38.p14_genomic.fna --sjdbGTFfile GCF_000001405.40_GRCh38.p14_genomic.gtf --genomeDir "genome${readLength}" --genomeChrBinNbits $factor --runThreadN 16 --sjdbOverhang $overhang
	
	#we must produce an RSEM index to use RSEM for quantitative analysis later. This step runs extremely quickly compared to STAR's indexing, and finishes easily on a regular computer
	cp GCF_000001405.40_GRCh38.p14_genomic.fna GCF_000001405.40_GRCh38.p14_genomic.fa
	rsem-prepare-reference --gtf GCF_000001405.40_GRCh38.p14_genomic.gtf --num-threads 16 GCF_000001405.40_GRCh38.p14_genomic.fa GCF_000001405.40_GRCh38.p14_genomic

	#leave refGen to return to normal
	cd ../
	
fi

#at this point it is assumed that there is a folder containing indexed information at refGen/genome/

#the next step in the RNA-seq pipeline is alignment of .fq to the reference genome
#because the previous step was performed in STAR, in order to get BAM files this too must be performed with STAR.

my_array[0]=foo
my_array[1]=yobo
for i in "${!my_array[@]}"
do 
    echo "${my_array[$i]}"
done

samples[0]=""

if [[ paired == 1 ]]
then
	for filename in $fqDir/*.fastq.gz
	do
		mv $filename ${filename%.fastq.gz}.fq.gz
	done
	
	for filename in $fqDir/*_1.fq.gz
	do
		it=0
		if [ -f "${filename%_1.fq.gz}_2.fq.gz" ]
		then
			samples[$it]="${filename%_1.fq.gz}"
			it=$(expr $it + 1)
		else
			"ERROR: no matching pair for $filename ; will not be included"
		fi	
	done
	mkdir ${oDir%/}
	chmod -R 0777 ${oDir%/}
	mkdir "${oDir%/}/$fqDir"
	chmod -R 0777 "${oDir%/}/$fqDir"
		
	if [[ $align != 0 ]]
	then
		for i in "${!samples[@]}"
		do 
		   	base="${samples[$i]}"
		    	read1=$fqDir/${samples[$i]}_1.fq.gz
			read2=$fqDir/${samples[$i]}_2.fq.gz
		        echo "unzipping $read1"
			gzip -d ${read1}
			echo "unzipping $read2"
			gzip -d ${read2}

			echo "aligning $base"

			STAROPTS="--outSAMattributes NH HI AS NM MD \
				--outFilterType BySJout \
				--outFilterMultimapNmax 20 \
				--outFilterMismatchNmax 999 \
				--outFilterMismatchNoverReadLmax 0.04 \
			    	--alignIntronMin 20 \
			    	--alignIntronMax 1000000 \
			    	--alignMatesGapMax 1000000 \
			    	--alignSJoverhangMin 8 \
			    	--alignSJDBoverhangMin 1 \
			    	--sjdbScore 1 \
			    	--limitBAMsortRAM 50000000000"

			 STAR --genomeDir refGen/genome$readLength --readFilesIn ${read1%.gz} ${read2%.gz} --outFileNamePrefix "$oDir${base}" --runThreadN 16 --quantMode TranscriptomeSAM ${STAROPTS}

			 echo "calculating expression of ${base}"
			 rsem-calculate-expression --num-threads 16 --paired-end --alignments "${oDir}${base}Aligned.toTranscriptome.out.bam" refGen/GCF_000001405.40_GRCh38.p14_genomic "${base}"
		done
	fi
else
	#unpaired pipeline
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
			STAR --genomeDir refGen/genome --readFilesIn ${filename} --outFileNamePrefix "$oDir${filename%.*}" --runThreadN 16 --quantMode TranscriptomeSAM
		done

		for filename in $fqDir/*.fastq 
		do
			echo "aligning $filename"
			STAR --genomeDir refGen/genome --readFilesIn ${filename} --outFileNamePrefix "$oDir${filename%.*}" --runThreadN 16 --quantMode TranscriptomeSAM
		done

		#quantifying gene expression
		#cd ${oDir%/}

		#it is very likely that the following for loop can be removed entirely
		#for filename in $oDir$fqDir/*.out.sam
		#do
		#	echo "converting $filename to "${filename%.*}.bam""
		#	samtools view --threads 16 -S -b $filename > "${filename%.*}.bam"
		#done

		for filename in $oDir$fqDir/*.toTranscriptome.out.bam
		do
			echo "calculating expression of ${filename}"
			sample_fname="${filename##*/}"
			sample_name="${sample_fname%.*}"
			rsem-calculate-expression --num-threads 16 --alignments "${filename%.*}.bam" refGen/GCF_000001405.40_GRCh38.p14_genomic "${sample_name}"
		done

		#cd ../
		#at this point we will have gene.results files?? and they can be used in the R language with Deseq2
	fi

	if [[ $quantify != 0 ]]
	then

	fi
fi
