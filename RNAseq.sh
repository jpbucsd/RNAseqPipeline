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
padj=0.5
padFlag=0
logT=5
logTFlag=0
logA=30
logAFlag=0
threadFlag=0
threads=0
qReport=0
trim=0
gFlag=0
gList=0
geneList=""
HgeneList=""
hgList=0;
hgFlag=0;

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
	elif [[ $padFlag == 1 ]]
	then
		padFlag=0
		padj=$var
	elif [[ $logTFlag == 1 ]]
	then
		logTFlag=0
		logT=$var
	elif [[ $logAFlag == 1 ]]
	then
		logAFlag=0
		logA=$var
	elif [[ $threadFlag == 1 ]]
	then
		threadFlag=0
		threads=$var
	elif [[ $sFlag == 1 ]]
	then
		sFlag=0
		slr=$var
	elif [[ $gFlag == 1 ]]
	then
		gFlag=0
		geneList=$var
	elif [[ $hgFlag == 1 ]]
	then
		hgFlag=0
		HgeneList=$var
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
		if [[ "$var" == *"-t"* ]]
		then
			threadFlag=1
		fi
		if [[ "$var" == "-"*"PADJ" ]]
		then
			padFlag=1
		fi
		if [[ "$var" == "-"*"log10" ]]
		then
			logTFlag=1
		fi
		if [[ "$var" == "-"*"A" ]]
		then
			logAFlag=1
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
			logAFlag=0
		fi
		if [[ "$var" == "-Q"* ]]
		then
			qReport=1
		fi
		if [[ "$var" == "-T"* ]]
		then
			trim=1
		fi
		if [[ "$var" == "-G"* ]]
		then
			gFlag=1
			gList=1
		fi
		if [[ "$var" == "-H"* ]]
		then
			hgFlag=1
			hgList=1
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


if [[ $qReport != 0 ]]
then
		cd "$fqDir"
		for filename in *.fastq.gz
		do
			mv $filename ${filename%.fastq.gz}.fq.gz
		done
		
		

		for filename in *.fq.gz
		do
			if [ ! -f "${filename%.*.*}.fq" ]
			then
				echo "unzipping $filename"
				gzip -dc $filename > "${filename%.*.*}.fq"
			fi	
		done


		mkdir quality
		chmod -R 0777 quality
		
		for filename in *.fq
		do
			echo "assessing quality of $filename"
			fastqc -o quality $filename
		done
		
		cd $rsd	
fi

if [[ $trim != 0 ]]
then
	#todo
	#${pair1}
	
	cd "$fqDir"
	for filename in *.fastq.gz
	do
		mv $filename ${filename%.fastq.gz}.fq.gz
	done
	for filename in *.fq.gz
	do
		if [ ! -f "${filename%.*.*}.fq" ]
		then
			echo "unzipping $filename"
			gzip -dc $filename > "${filename%.*.*}.fq"
		fi
	done
	
	for filename in *${pair1}.fq
	do
		if [ -f "${filename%${pair1}.fq}${pair2}.fq" ]
		then
			echo "cleaning $filename and ${filename%${pair1}.fq}${pair2}.fq"
			fastp -i "${filename}" -o "${filename%${pair1}*}${pair1}_clean.fq" -I "${filename%${pair1}.fq}${pair2}.fq" -O "${filename%${pair1}*}${pair2}_clean.fq"
		else
			"ERROR: no matching pair for $filename ; will not be cleaned/included"
		fi	
	done
		
	for filename in *_clean.fq
	do	
		echo "removing ${filename%_clean.fq}.fq"
		rm "${filename%_clean.fq}.fq"
		
		echo "assessing quality of cleaned file $filename"
		fastqc -o quality $filename
		
		echo "renaming $filename to ${filename%_clean.fq}.fq"
		mv $filename "${filename%_clean.fq}.fq"
	done
fi

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
		STAR --runMode genomeGenerate --genomeFastaFiles GCF_000001405.40_GRCh38.p14_genomic.fna --sjdbGTFfile GCF_000001405.40_GRCh38.p14_genomic.gtf --genomeDir "genome${readLength}" --genomeChrBinNbits $factor --runThreadN $threads --sjdbOverhang $overhang

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
		rsem-prepare-reference --gtf GCF_000001405.40_GRCh38.p14_genomic.gtf --num-threads $threads GCF_000001405.40_GRCh38.p14_genomic.fa GCF_000001405.40_GRCh38.p14_genomic
	fi

	#leave refGen and return to normal
	cd ../

	#at this point it is assumed that there is a folder containing indexed information at refGen/genome/

	#the next step in the RNA-seq pipeline is alignment of .fq to the reference genome
	#because the previous step was performed in STAR, in order to get BAM files this too must be performed with STAR.

	samples[0]=""

	if [[ $paired == 1 ]]
	then
		cd "$fqDir"
		for filename in *.fastq.gz
		do
			mv $filename ${filename%.fastq.gz}.fq.gz
		done

		for filename in *.fq.gz
		do
			if [ ! -f "${filename%.*.*}.fq" ]
			then
				echo "unzipping $filename"
				gzip -dc $filename > "${filename%.*.*}.fq"
			fi
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
				
				 STAR --genomeDir ${rsd}/refGen/genome$readLength --readFilesIn ${read1} ${read2} --outFileNamePrefix "${oDir}Alignment/${base}" --runThreadN $threads --quantMode TranscriptomeSAM ${STAROPTS}
				 
				 echo "calculating expression of ${base}"
				 rsem-calculate-expression --num-threads $threads --paired-end --alignments "${oDir}Alignment/${base}Aligned.toTranscriptome.out.bam" ${rsd}/refGen/GCF_000001405.40_GRCh38.p14_genomic "${oDir}${base}"
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
				STAR --genomeDir refGen/genome --readFilesIn ${filename} --outFileNamePrefix "$oDir${filename%.*}" --runThreadN $threads --quantMode TranscriptomeSAM
			done

			for filename in $fqDir/*.fastq 
			do
				echo "aligning $filename"
				STAR --genomeDir refGen/genome --readFilesIn ${filename} --outFileNamePrefix "$oDir${filename%.*}" --runThreadN $threads --quantMode TranscriptomeSAM
			done

			#quantifying gene expression
			#cd ${oDir%/}

			#it is very likely that the following for loop can be removed entirely
			#for filename in $oDir$fqDir/*.out.sam
			#do
			#	echo "converting $filename to "${filename%.*}.bam""
			#	samtools view --threads $threads -S -b $filename > "${filename%.*}.bam"
			#done

			for filename in $oDir$fqDir/*.toTranscriptome.out.bam
			do
				echo "calculating expression of ${filename}"
				sample_fname="${filename##*/}"
				sample_name="${sample_fname%.*}"
				rsem-calculate-expression --num-threads $threads --alignments "${filename%.*}.bam" refGen/GCF_000001405.40_GRCh38.p14_genomic "${sample_name}"
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
			if [[ it -ge 1 ]]
			then
				factors[$(expr $it - 1)]=$item
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
		#the first element in this array will be a string, of which each index determines the type of parameter in the rest of the string
		#1 represents the value for comparison 1, 2 represents the value for comparison 2
		#A represents a parameter that is unioned with the parameter denoted by 1. therefore either 1 or A will be compared against 2
		#B represents a parameter that is unioned with the parameter denoted by 2. therefore 1  will be compared against 2 or B
		#c represents a parameter that is intersected with 1. Therefore samples within C and 1 will be compared against 2
		#D represents a paraemter that is intersected with 2, therefore samples within D and 2 will be compared against 1
		determinant="${parsedArray[1]}";
		c1union=(${parsedArray[2]});
		c2union=(${parsedArray[3]});
		c1intersect=();
		c2intersect=();
		for (( s=2; s<${#determinant}; s++ )); do
			#we want to skip the first two letters, they will always be 1 and 2
			echo "${determinant:$s:1}"
			if [ "${determinant:$s:1}" == "A" ]
			then
				c1union+=${parsedArray[$((s+2))]};
			elif [ "${determinant:$s:1}" == "B" ]
			then
				c2union+=${parsedArray[$((s+2))]};
			elif [ "${determinant:$s:1}" == "C" ]
			then
				c1intersect+=${parsedArray[$((s+2))]};
			elif [ "${determinant:$s:1}" == "D" ]
			then
				c2intersect+=${parsedArray[$((s+2))]};
			fi
  			#add to lists of union comparison one and intersect comparison 1. and also for 2. 
			#then for the union, grep from stepfile twice into stempfile13
			#for intersection grep from stepfile13 into 23
		done

		#union on 1
		for i in "${c1union[@]}"
		do
		   cat $slr | grep -n "~$i}" >> stempFile1.slr
		done
		
		#union on 2
		for i in "${c2union[@]}"
		do
		   cat $slr | grep -n "~$i}" >> stempFile2.slr
		done
		
		#intersection on 1
		for i in "${c1intersect[@]}"
		do
		   cat stempFile1 | grep -n "~$i}" > stempFile11.slr
		   cat stempFile11.slr > stempFile1.slr
		done
		
		#intersection on 2
		for i in "${c2intersect[@]}"
		do
		   cat stempFile2 | grep -n "~$i}" > stempFile22.slr
		   cat stempFile22.slr > stempFile2.slr
		done
		
		#collect the samples after unions and intersections
		
		firstFactor[0]=""
		it=0
		
		while IFS=$'\t' read -r -a sparsedArray
		do
			firstFactor[$it]=${sparsedArray[1]}
			it=$(expr $it + 1)
			echo firstFactor[$it]
		done < stempFile1.slr
		rm stempFile1.slr
		rm stempFile11.slr
		
		secondFactor[0]=""
		it=0
		
		while IFS=$'\t' read -r -a sparsedArray
		do
			secondFactor[$it]=${sparsedArray[1]}
			it=$(expr $it + 1)
			echo secondFactor[$it]
		done < stempFile2.slr
		rm stempFile2.slr
		rm stempFile22.slr

		if [[ $analysis != 0 ]]
		then
			firstFactorName=${factors[$((${parsedArray[2]}))]}
			secondFactorName=${factors[$((${parsedArray[3]}))]}
			echo "Rscript DifferentialExpression.R -1 $firstFactorName ${firstFactor[@]/%/.genes.results} -2 $secondFactorName ${secondFactor[@]/%/.genes.results} -d ${fqDir}/${oDir}"

			Rscript DifferentialExpression.R -1 $firstFactorName "${firstFactor[@]/%/.genes.results}" -2 $secondFactorName "${secondFactor[@]/%/.genes.results}" -d "${fqDir}/${oDir}"

			if [[ $gList == 0 ]]
			then
				echo "python results.py -f ${fqDir}/${oDir}/${firstFactorName}_vs_${secondFactorName}.csv --padj $padj --log10 $logT --Llog10 $logA --odir ${fqDir}/${oDir}/results"

				python results.py -f "${fqDir}/${oDir}/${firstFactorName}_vs_${secondFactorName}.csv" --padj $padj --log10 $logT --Llog10 $logA --odir "${fqDir}/${oDir}/results"
			else
				echo "python results.py -f ${fqDir}/${oDir}/${firstFactorName}_vs_${secondFactorName}.csv --padj $padj --log10 $logT --Llog10 $logA --odir ${fqDir}/${oDir}/results --g ${geneList}"

				python results.py -f "${fqDir}/${oDir}/${firstFactorName}_vs_${secondFactorName}.csv" --padj $padj --log10 $logT --Llog10 $logA --odir "${fqDir}/${oDir}/results" --g $geneList
			fi
			
		fi
		
		if [[ $pca != 0 ]]
		then
			first_PRE=("${firstFactor[@]/#/$fqDir/$oDir/}")
			second_PRE=("${secondFactor[@]/#/$fqDir/$oDir/}")
			
			echo "python PCA.py -f ${first_PRE[@]/%/.genes.results} ${second_PRE[@]/%/.genes.results} --oDir ${fqDir}/${oDir} --numComps 2"
			
			python PCA.py -f "${first_PRE[@]/%/.genes.results}" "${second_PRE[@]/%/.genes.results}" --oDir "${fqDir}/${oDir}" --numComps 2
		fi
	done < tempFile.slr
	rm tempFile.slr
	
	if [[ $analysis != 0 ]]
	then
		heat=""
		cat $slr | grep -n "Comparison:" > tempFile.slr
		IFS=$'\t' read -r -a parsedArray < tempFile.slr
		for ELEMENT in ${parsedArray[@]}; do
			heat+="${ELEMENT} "
		done
		rm tempFile.slr
		#moving into the output directory since heat maps cannot be saved into a specific directory, and just get autosaved to the cd
		cd ${fqDir}/${oDir}/
		
		if [[ $hgList == 0 ]]
		then
			echo "python ${rsd}/heatmap.py -f normalizedCounts.csv"
		
			$python {rsd}/heatmap.py -f normalizedCounts.csv
			
		else
			echo "python ${rsd}/heatmap.py -f normalizedCounts.csv --g ${HgeneList}"
		
			python ${rsd}/heatmap.py -f normalizedCounts.csv --g $HgeneList
			
		fi
		
		#the heatmap wont save to the correct place, we must move it
		echo "mv ${rsd}/heatmap.eps ${fqDir}/${oDir}/heatmap.eps"
		mv ${rsd}/heatmap.eps ${fqDir}/${oDir}/heatmap.eps
		 
		#in the future the following parameters should be added to the heatmap command through command line variables
		#--g /projects/ps-bryansunlab/labTools/RNAseq/geneListHeatTest.txt -c
		cd $rsd	
	fi
fi
