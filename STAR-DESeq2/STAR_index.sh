#!/bin/bash

STAR --runThreadN 28 --runMode genomeGenerate \
     --genomeDir /home/sucsb/Databases/genomes/mouse/STAR_index \
     --genomeFastaFiles /home/sucsb/Databases/genomes/mouse/Mus_musculus.GRCm38.dna.toplevel.fa \
     --sjdbGTFfile /home/sucsb/Databases/genomes/mouse/Mus_musculus.GRCm38.82.chr.gtf \
     --sjdbOverhang 99 
