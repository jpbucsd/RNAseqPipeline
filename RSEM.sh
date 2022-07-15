#!/bin/bash
#this file is the same as RNAseq.sh however it only includes code pertinent to the usage of RSEM
#the purpose is for when the original script failed RSEM but succeeded everything else

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

if [[ $indexF == 1 ]]
then
	cd refGen
	#i think this indexing (bowtie2) is an alterative to star and it will not be used for this pipeline
	rsem-prepare-reference --gtf GCF_000001405.40_GRCh38.p14_genomic.gtf --num-threads 16 --bowtie2 GCF_000001405.40_GRCh38.p14_genomic.fna GCF_000001405.40_GRCh38.p14_genomic
	#leave refGen to return to normal
	cd ../
	
fi

if [[ $align != 0 ]]
then
	for filename in $oDir$fqDir/*.out.sam
	do
		echo "converting $filename to "${filename%.*}.bam""
		samtools view --threads 16 -S -b $filename > "${filename%.*}.bam"
		echo "calculating expression of ${filename}
		sample_fname="${filename##*/}"
		sample_name="${sample_fname%.*}"
		rsem-calculate-expression --num-threads 16 --bam "${filename%.*}.bam" refGen/genome "${sample_name}"
	done
fi
