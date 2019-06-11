#!/bin/bash 
# Map the reads against reference and estimate gene and isoform expression of three Mock cell samples 

REF=/home/sucsb/Databases/genomes/mouse/mouse_ref
bowtiePath=/home/sucsb/bioinfo-tools/RSEM_tutorial-master/software/bowtie2-2.2.6

# Map reads and calculate counts for each fastq file
for item in `ls `
do
	if [ -f "$item" ] && [[ $item == *".fastq.gz" ]] ; then
		seqtk trimfq $item > ${item%%.*}.fastq
		rsem-calculate-expression -p 28 --bowtie2 --bowtie2-path $bowtiePath ${item%%.*}.fastq $REF ${item%%_*}
	fi
done

# DEG analysis: FDR < 0.05

rsem-generate-data-matrix SRR1731203.genes.results SRR1731204.genes.results SRR1731205.genes.results SRR1731206.genes.results SRR1731199.genes.results SRR1731200.genes.results SRR1731201.genes.results SRR1731202.genes.results   >  CntMat_mutNFAT1vsMock.txt 
rsem-run-ebseq CntMat_mutNFAT1vsMock.txt 4,4 geneMat.results 
rsem-control-fdr geneMat.results 0.05 mutNFAT1vsMock.deg


