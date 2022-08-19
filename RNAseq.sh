#!/bin/bash

#important variables

fqDir=""
oDir=""
indexF=0
fqFlag=0
oFlag=0
align=0
readLength=50
quantify=1
rFlag=0
pFlag=0
paired=0
pair1="_1"
pair2="_2"
slr=""
sFlag=0
analysis=0
pca=0


rsd=$(pwd) #directory where RNA seq commands are stored, may not be root of fastq reads

#usage of command line arguments
#-h - help //todo
#-f fastqdirectory
#-o output directory

#-r readLength
#-p signals paired ends (should have the same name and end in _1 and _2, otherwise type what files will end in)

#-index (compute index) usage of index with no other arguments will cause only indexing to occur. without this flag indexing will not occur //depricated
#-s slr file, specifications from graphical interface
#input command line arguments
for var in "$@"
do
	if [[ $fqFlag == 1 ]]
	then
		fqFlag=0
		fqDir="$var"
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
	elif [[ $sFlag == 1 ]]
	then
		sFlag=0
		slr=$var
	elif [[ "$var" == "-"* ]]
	then
		#some flags may have multiple phrases following them so the -option comes first in the else ifs. 
		#unflag these first
		if [[ $pFlag != 0 ]]
		then
			pFlag=0
		fi
		if [[ "$var" == *"-f"* ]]
		then
			#fast q files coming up
			fqFlag=1
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
		if [[ "$var" == *"-s"* ]]
		then
			sFlag=1
		fi
		if [[ "$var" == "-"*"d"* ]]
		then
			analysis=1
		fi
		if [[ "$var" == "-"*"a"* ]]
		then
			align=2
		fi
		if [[ "$var" == "-"*"PCA"* ]]
		then
			pca=1
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
		fi
	fi
done

#this script produces an index file from the referencce genome (fasta format) and a GTF annotations file.
#the produced index file will next be used in the second step written in the script align.sh
#this script will only execute if the instruction to generate an index is included.

#check if genome has been indexed for the proper settings
if [[ $align != 0 ]]
then
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

	if [ -d "genome$readLength" ] 
	then
		echo "genome for length $readLength has already been indexed with STAR"
	else
		#must index
		echo "genome for length $readLength has not been indexed with STAR, moving to indexing"
		#indexF=1




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

	fi

	if [ -f *.n2g.idx.fa ] 
	then
		echo "genome has already been indexed with RSEM"
	else
		#we must produce an RSEM index to use RSEM for quantitative analysis later. This step runs extremely quickly compared to STAR's indexing, and finishes easily on a regular computer
		echo "genome has not yet been indexed with RSEM, indexing..."

		if ! [ -f GCF_000001405.40_GRCh38.p14_genomic.fa ]
		then
			cp GCF_000001405.40_GRCh38.p14_genomic.fna GCF_000001405.40_GRCh38.p14_genomic.fa
		fi
		rsem-prepare-reference --gtf GCF_000001405.40_GRCh38.p14_genomic.gtf --num-threads 16 GCF_000001405.40_GRCh38.p14_genomic.fa GCF_000001405.40_GRCh38.p14_genomic
	fi

	#leave refGen and return to normal
	cd ../

	#at this point it is assumed that there is a folder containing indexed information at refGen/genome/

	#the next step in the RNA-seq pipeline is alignment of .fq to the reference genome
	#because the previous step was performed in STAR, in order to get BAM files this too must be performed with STAR.

	samples[0]=""

	if [[ $paired == 1 ]]
	then
		cd $fqDir
		for filename in *.fastq.gz
		do
			mv $filename ${filename%.fastq.gz}.fq.gz
		done

		for filename in *.fq.gz
		do
			echo "unzipping $filename"
			gzip -d $filename
		done
		it=0
		for filename in *${pair1}.fq
		do
			if [ -f "${filename%${pair1}.fq}${pair2}.fq" ]
			then
				samples[$it]="${filename%${pair1}.fq}"
				it=$(expr $it + 1)
			else
				"ERROR: no matching pair for $filename ; will not be included"
			fi	
		done
		mkdir ${oDir%/}
		chmod -R 0777 ${oDir%/}
		mkdir "${oDir%/}/Alignment"
		chmod -R 0777 "${oDir%/}/Alignment"

		if [[ $align != 0 ]]
		then
			for i in "${!samples[@]}"
			do 
				base="${samples[$i]}"
				read1=${samples[$i]}${pair1}.fq
				read2=${samples[$i]}${pair2}.fq


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

				 STAR --genomeDir ${rsd}/refGen/genome$readLength --readFilesIn ${read1} ${read2} --outFileNamePrefix "${oDir}Alignment/${base}" --runThreadN 16 --quantMode TranscriptomeSAM ${STAROPTS}

				 echo "calculating expression of ${base}"
				 rsem-calculate-expression --num-threads 16 --paired-end --alignments "${oDir}Alignment/${base}Aligned.toTranscriptome.out.bam" ${rsd}/refGen/GCF_000001405.40_GRCh38.p14_genomic "${oDir}${base}"
			done
		fi
		#leave the raw read directory, where all alignments and quantifications have now been saved to
		#cd ../
		cd $rsd
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
	fi
fi
if [[ $analysis != 0 ]] || [[ $pca != 0 ]]
then
	factors[0]=""
	cat $slr | grep -n "Attributes:" > tempFile.slr
	while IFS=$'\t' read -r -a parsedArray
	do
		it=0
		for item in "${parsedArray[@]}"
		do
			if [[ it -ge 0 ]]
			then
				factors[$it]=$item
			fi
			it=$(expr $it + 1)
		done
	done < tempFile.slr
	rm tempFile.slr
	
	factor1="0"
	factor2="0"
	cat $slr | grep -n "Comparison:" > tempFile.slr
	while IFS=$'\t' read -r -a parsedArray
	do
		factor1="${parsedArray[1]}"
		factor2="${parsedArray[2]}"

		firstFactor[0]=""
		it=0
		cat $slr | grep -n "~${factor1}" > stempFile.slr
		while IFS=$'\t' read -r -a sparsedArray
		do
			firstFactor[$it]=${sparsedArray[1]}
			it=$(expr $it + 1)
		done < stempFile.slr
		rm stempFile.slr

		secondFactor[0]=""
		it=0
		cat $slr | grep -n "~${factor2}" > stempFile.slr
		while IFS=$'\t' read -r -a sparsedArray
		do
			secondFactor[$it]=${sparsedArray[1]}
			it=$(expr $it + 1)
		done < stempFile.slr
		rm stempFile.slr

		if [[ $analysis != 0 ]]
		then
			firstFactorName=${factors[$(($factor1))]}
			secondFactorName=${factors[$(($factor2))]}
			echo "Rscript DifferentialExpression.R -1 $firstFactorName ${firstFactor[@]/%/.genes.results} -2 $secondFactorName ${secondFactor[@]/%/.genes.results} -d ${fqDir}/${oDir}"

			Rscript DifferentialExpression.R -1 $firstFactorName "${firstFactor[@]/%/.genes.results}" -2 $secondFactorName "${secondFactor[@]/%/.genes.results}" -d "${fqDir}/${oDir}"

			echo "python results.py -f ${fqDir}/${oDir}/${firstFactorName}_vs_${secondFactorName}.csv --padj 0.5 --log10 5 --Llog10 30 --odir ${fqDir}/${oDir}/results"

			python results.py -f "${fqDir}/${oDir}/${firstFactorName}_vs_${secondFactorName}.csv" --padj 0.5 --log10 5 --Llog10 30 --odir "${fqDir}/${oDir}/results"
		fi
		
		if [[ $pca != 0 ]]
		then
			echo "python PCA.py -f ${firstFactor[@]/%/.genes.results} ${secondFactor[@]/%/.genes.results} --oDir ${fqDir}/${oDir} --numComps 2"
			
			python PCA.py -f "${firstFactor[@]/%/.genes.results}" "${secondFactor[@]/%/.genes.results}" --oDir "${fqDir}/${oDir}" --numComps 2
		fi
	done < tempFile.slr
	rm tempFile.slr
fi
