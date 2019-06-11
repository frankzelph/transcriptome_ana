#!/bin/bash 
# Map the reads against reference and estimate gene and isoform expression of three Mock cell samples 

REF_dir=/home/sucsb/Databases/genomes/mouse/STAR_index
GTFfile=/home/sucsb/Databases/genomes/mouse/Mus_musculus.GRCm38.82.chr.gtf

# Map reads and calculate counts for each fastq file
for item in `ls ..`
do
	if [ -f "../$item" ] && [[ $item == *".fastq" ]] ; then
		#seqtk $item > ${item%%.*}_trimed.fastq
		# mapping reads
		STAR --runThreadN 28 --genomeDir $REF_dir \
			 --readFilesIn ../$item --genomeLoad LoadAndRemove \
			 --outFilterMultimapNmax 1 --outFilterType BySJout \
			 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
			 --alignIntronMin 20 --alignIntronMax 1000000 \
			 --alignMatesGapMax 100000 --outFilterMismatchNmax 0 \
			 --outFileNamePrefix ./${item%%.*}.
	fi
done

# Get count matrix using featureCounts
featureCounts -T 28 -g gene_name -s 0 -a $GTFfile\
			  -o CntMatrix.txt *.Aligned.out.sam

# DEG analysis





