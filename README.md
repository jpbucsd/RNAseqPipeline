# RNAseqPipeline
RNA sequencing pipeline for sun lab

# Indexing the genome
The first required step is indexing the human genome using genome annotations.
The following command will produce a folder called genome containing index information STAR will use to align reads in later steps
##### bash RNAseq.sh -index
This command should not be used a second time for future sequencing as it is repetetive and the way it is currently set up it costs 13 core hours.

