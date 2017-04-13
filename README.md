# transcriptome_ana
Record the process for transcriptome analysis, including codes.
RNA-seq analysis process

1. Genome dependent analysis (for human dataset). – TopHat → cufflinks → CummeRbund
1.1 Mapping reads to genome (TopHat)

Preparing the environment:
Intall tophat2 and bowtie2. (Google it)
Note: This bio-linux machine is already prepared.

To use tophat, firstly, the reference genome sequence and annotation files (.gtf) should be downloaded from Ensembl, UCSC, or NCBI database. 
I downloaded human genome reference and annotation files from Ensembl.
Homo_sapiens.GRCh38.dna.toplevel.fa.gz: ftp://ftp.ensembl.org/pub/release-85/fasta/homo_sapiens/dna/
Homo_sapiens.GRCh38.85.gtf.gz: ftp://ftp.ensembl.org/pub/release-85/gtf/homo_sapiens/

Decompress these two files. Use the sequence file to build a bowtie2 index:
	gunzip Homo_sapiens.GRCh38.dna.toplevel.fa.gz
	bowtie2-build  Homo_sapiens.GRCh38.dna.toplevel.fa genome

Use tophat2 to map reads to genome. Write a bash file (named as tophat2.sh) as below:

#!/bin/bash 

# Map calu3 mock repeat 1 paired reads 
tophat2 -p 12 -G genes.gtf -o Calu3_M_R1_thout genome WGC072920-Calu3-M-1_combined_R1.fastq.gz WGC072920-Calu3-M-1_combined_R2.fastq.gz 

# Map calu3 mock repeat 2 paired reads 
tophat2 -p 12 -G genes.gtf -o Calu3_M_R2_thout genome WGC072920-Calu3-M-2_combined_R1.fastq.gz WGC072920-Calu3-M-2_combined_R2.fastq.gz 

# Map calu3 mock repeat 3 paired reads 
tophat2 -p 12 -G genes.gtf -o Calu3_M_R3_thout genome WGC072920-Calu3-M-3_combined_R1.fastq.gz WGC072920-Calu3-M-3_combined_R2.fastq.gz 

…

tophat2 command line  explaining:
-p 12				# use 12 threads
-G genes.gtf			# annotation file (*.gtf)
-o Calu3_M_R1_thout	# output result files into folder Calu3_M_R1_thout
genome			# bowtie2 index
WGC072920-Calu3-M-1_combined_R1.fastq.gz WGC072920-Calu3-M-
1_combined_R2.fastq.gz
				# paired-end read files: read1_file, read2_file


Write all commands in the bash file to map all RNA-seq data to the reference. Then run those commands in terminal as below:


$ chmod 755 tophat2.sh	# modify the file to be capable of execution
$ ./tophat2.sh			# run the bash file


1.2 Assemble the mapped reads (cufflinks)
Write a bash file (cufflinks.sh) as below:

#!/bin/bash 
#cufflinks: assemble transcripts 

cufflinks -p 12 -o Calu3_M_R1_clout Calu3_M_R1_thout/accepted_hits.bam 
cufflinks -p 12 -o Calu3_M_R2_clout Calu3_M_R2_thout/accepted_hits.bam 
cufflinks -p 12 -o Calu3_M_R3_clout Calu3_M_R3_thout/accepted_hits.bam 
cufflinks -p 12 -o Calu3_6h_R1_clout Calu3_6h_R1_thout/accepted_hits.bam 
cufflinks -p 12 -o Calu3_6h_R2_clout Calu3_6h_R2_thout/accepted_hits.bam 
cufflinks -p 12 -o Calu3_6h_R3_clout Calu3_6h_R3_thout/accepted_hits.bam 
cufflinks -p 12 -o Calu3_12h_R1_clout Calu3_12h_R1_thout/accepted_hits.bam 
cufflinks -p 12 -o Calu3_12h_R2_clout Calu3_12h_R2_thout/accepted_hits.bam 
cufflinks -p 12 -o Calu3_12h_R3_clout Calu3_12h_R3_thout/accepted_hits.bam 
cufflinks -p 12 -o Calu3_24h_R1_clout Calu3_24h_R1_thout/accepted_hits.bam 
cufflinks -p 12 -o Calu3_24h_R2_clout Calu3_24h_R2_thout/accepted_hits.bam 
cufflinks -p 12 -o Calu3_24h_R3_clout Calu3_24h_R3_thout/accepted_hits.bam 
cufflinks -p 12 -o Calu3_36h_R1_clout Calu3_36h_R1_thout/accepted_hits.bam 
cufflinks -p 12 -o Calu3_36h_R2_clout Calu3_36h_R2_thout/accepted_hits.bam 
cufflinks -p 12 -o Calu3_36h_R3_clout Calu3_36h_R3_thout/accepted_hits.bam


Run it in terminal as below:

$ chmod 755 cufflinks.sh	# modify the file to be capable of execution
$ ./cufflinks.sh		# run the bash file


1.3 Merge all assembled reads (cuffmerge)
Write a text file to imply the files (assemlies.txt) which are going to be assembled.

./Calu3_M_R1_clout/transcripts.gtf 
./Calu3_M_R2_clout/transcripts.gtf 
./Calu3_M_R3_clout/transcripts.gtf 
./Calu3_6h_R1_clout/transcripts.gtf 
./Calu3_6h_R2_clout/transcripts.gtf 
./Calu3_6h_R3_clout/transcripts.gtf 
./Calu3_12h_R1_clout/transcripts.gtf 
./Calu3_12h_R2_clout/transcripts.gtf 
./Calu3_12h_R3_clout/transcripts.gtf 
./Calu3_24h_R1_clout/transcripts.gtf 
./Calu3_24h_R2_clout/transcripts.gtf 
./Calu3_24h_R3_clout/transcripts.gtf 
./Calu3_36h_R1_clout/transcripts.gtf 
./Calu3_36h_R2_clout/transcripts.gtf 
./Calu3_36h_R3_clout/transcripts.gtf

Run cuffmerge in terminal as below:

$ cuffmerge -g genes.gtf -s genome.fa -p 12 assemblies.txt

#******************Explaining********************#
-g	# implying the annotation file
-s	# implying the genome sequence file
-p	# the number of threads to be used


1.4 Calculate the expression level of genes and transcripts (cuffdiff)
Write a bash file (cuffdiff.sh) as below:

#!/bin/bash 
# cuffdiff: compare expression levels of genes and transcripts 

cuffdiff -o hu_6h_diff_out -b genome.fa -p 12 -u merged_asm/merged.gtf \ 
Calu3_M_R1_thout/accepted_hits.bam,Calu3_M_R2_thout/accepted_hits.bam,Calu3_M_R3_thout/accepted_hits.bam \ 
Calu3_6h_R1_thout/accepted_hits.bam,Calu3_6h_R2_thout/accepted_hits.bam,Calu3_6h_R3_thout/accepted_hits.bam \ 

cuffdiff -o hu_12h_diff_out -b genome.fa -p 12 -u merged_asm/merged.gtf \ 
Calu3_M_R1_thout/accepted_hits.bam,Calu3_M_R2_thout/accepted_hits.bam,Calu3_M_R3_thout/accepted_hits.bam \ 
Calu3_12h_R1_thout/accepted_hits.bam,Calu3_12h_R2_thout/accepted_hits.bam,Calu3_12h_R3_thout/accepted_hits.bam \ 

cuffdiff -o hu_24h_diff_out -b genome.fa -p 12 -u merged_asm/merged.gtf \ 
Calu3_M_R1_thout/accepted_hits.bam,Calu3_M_R2_thout/accepted_hits.bam,Calu3_M_R3_thout/accepted_hits.bam \ 
Calu3_24h_R1_thout/accepted_hits.bam,Calu3_24h_R2_thout/accepted_hits.bam,Calu3_24h_R3_thout/accepted_hits.bam \ 


cuffdiff -o hu_36h_diff_out -b genome.fa -p 12 -u merged_asm/merged.gtf \ 
Calu3_M_R1_thout/accepted_hits.bam,Calu3_M_R2_thout/accepted_hits.bam,Calu3_M_R3_thout/accepted_hits.bam \ 
Calu3_36h_R1_thout/accepted_hits.bam,Calu3_36h_R2_thout/accepted_hits.bam,Calu3_36h_R3_thout/accepted_hits.bam 


Run it in terminal as below:

$ chmod 755 cuffdiff.sh	# modify the file to be capable of execution
$ ./cuffdiff.sh			# run the bash file


1.5 Analysis and visualization through CummeRbund.
Open a terminal and change the directory to the parent path of cuffdiff output folder. Enter “R” terminal, and load “cummeRbund” library. Get differential expressed genes of all datasets.


# An example for getting differentially expressed genes
$ cd ~/Documents/Rs_Hu_rnaseq
$ R
> # Load cummeRbund library.
> library(cummeRbund)
> # Load differential expression data calculated by cufflinks (6 hpi vs mock).
> Mvs6h_cuff <- readCufflinks('hu_6h_diff_out')
> # Get the gene_ids of significant diff genes with FDR < 0.05.
> Mvs6h_sigGeneIds <- getSig(Mvs6h_cuff, alpha=0.05, level='genes') 
> # Get diff genes according to the gene_ids.
> Mvs6h_sigGenes <- getGenes(Mvs6h_cuff, Mvs6h_sigGeneIds) 
> # Generate a table including gene_names, log2 fold change, q-values and so on.
> Mvs6h_sigGenes_table <- data.frame(gene_name=Mvs6h_sigGenes@annotation[,4],  Mvs6h_sigGenes@diff[,5:10]) 
> # Get upregulated genes with  log2_fold_change > 1.
> upgenes_6h <- subset(Mvs6h_sigGenes_table, (log2_fold_change >1)) 
> # Get downregulated genes with  log2_fold_change < -1.
> downgenes_6h <- subset(Mvs6h_sigGenes_table, (log2_fold_change < -1)) 


Use the same method to get differentially expressed genes of the cells infected at 12, 24, 36 hpi compared to mock, repectively.

Work in R terminal. Draw a mirror bar plot to display the numbers of up- and down-regulated genes at each time point post infection.

> # Do statistics of the numbers of up- and down-regulated genes and organise them into a data.frame 
> stat <- data.frame( 
+ group=rep(c("up","down"), each=4), 
+ x=c("06h","12h","24h","36h"), 
+ y=c(nrow(upgenes_6h), nrow(upgenes_12h), nrow(upgenes_24h), nrow(upgenes_36h), nrow(downgenes_6h), nrow(downgenes_12h), nrow(downgenes_24h), nrow(downgenes_36h)) 
+ ) 
> # Set the number of down-regulated gene to negtive, to show their bars on the down side of the bar plot. 
> dat$y[dat$group=="Down"] <- -dat$y[dat$group=="Down"] 
> # Since the levels of group is in alphabet order by default, reverse it to make "up" displayed above in the figure. 
> levels(dat$group) 
> dat$group <- factor(dat$group, levels = rev(levels(dat$group))) 

> # Draw the mirror bar plot and output into a pdf file. 
> pdf(file='Hu_diff_genes/hu_updown_genes.pdf') 
> ggplot(dat, aes(x=x,y=y,fill=group, order=rev(group))) + geom_bar(stat="identity", position="identity") + scale_y_continuous(breaks=seq(-1300,1300, by =100), labels=abs(seq(-1300,1300,by=100)))+scale_fill_manual(values = alpha(c("red", "green"))) + labs(title = "Homo sapiens", x = "hpi", y = "Number DE genes") + geom_text(aes(label=abs(y)), vjust=0) 
> dev.off() 
> # Output the figure into a tiff format file. 
> tiff(filename='Hu_diff_genes/hu_updown_genes.tif') 
> ggplot(dat, aes(x=x,y=y,fill=group, order=rev(group))) + geom_bar(stat="identity", position="identity") + scale_y_continuous(breaks=seq(-1300,1300, by =100), labels=abs(seq(-1300,1300,by=100)))+scale_fill_manual(values = alpha(c("red", "green"))) + labs(title = "Homo sapiens", x = "hpi", y = "Number DE genes") + geom_text(aes(label=abs(y)), vjust=0) 
> dev.off() 


The figure is show as below.


























2. Genome independent analysis (for Rs datasets). Trinity → Annotation → RSEM

Firstly, download and install all required packages:
Trinity: https://github.com/trinityrnaseq/trinityrnaseq/wiki
blast+: Defaultly installed in bio-linux
Python: Defaultly installed in bio-linux
RSEM: https://github.com/bli25ucb/RSEM_tutorial

2.1 Trinity: assemble all Rs reads to get a reference transcriptome.
Put all sequence files of Rs into a directory. Open a terminal. Write a bash file as below and run it.




 
#!/bin/bash 

# Trinity assemble all reads 
Trinity --seqType fq --max_memory 30G --left WGC072920-RsK-M-1_combined_R1.fastq.gz,WGC072920-RsK-M-2_combined_R1.fastq.gz,WGC072920-RsK-M-3_combined_R1.fastq.gz,WGC072920-RsK-6h-1_combined_R1.fastq.gz,WGC072920-RsK-6h-2_combined_R1.fastq.gz,WGC072920-RsK-6h-3_combined_R1.fastq.gz,WGC072920-RsK-12h-1_combined_R1.fastq.gz,WGC072920-RsK-12h-2_combined_R1.fastq.gz,WGC072920-RsK-12h-3_combined_R1.fastq.gz,WGC072920-RsK-24h-1_combined_R1.fastq.gz,WGC072920-RsK-24h-2_combined_R1.fastq.gz,WGC072920-RsK-24h-3_combined_R1.fastq.gz,WGC072920-RsK-36h-1_combined_R1.fastq.gz,WGC072920-RsK-36h-2_combined_R1.fastq.gz,WGC072920-RsK-36h-3_combined_R1.fastq.gz --right WGC072920-RsK-M-1_combined_R2.fastq.gz,WGC072920-RsK-M-2_combined_R2.fastq.gz,WGC072920-RsK-M-3_combined_R2.fastq.gz,WGC072920-RsK-6h-1_combined_R2.fastq.gz,WGC072920-RsK-6h-2_combined_R2.fastq.gz,WGC072920-RsK-6h-3_combined_R2.fastq.gz,WGC072920-RsK-12h-1_combined_R2.fastq.gz,WGC072920-RsK-12h-2_combined_R2.fastq.gz,WGC072920-RsK-12h-3_combined_R2.fastq.gz,WGC072920-RsK-24h-1_combined_R2.fastq.gz,WGC072920-RsK-24h-2_combined_R2.fastq.gz,WGC072920-RsK-24h-3_combined_R2.fastq.gz,WGC072920-RsK-36h-1_combined_R2.fastq.gz,WGC072920-RsK-36h-2_combined_R2.fastq.gz,WGC072920-RsK-36h-3_combined_R2.fastq.gz  --CPU 12 --no_cleanup --normalize_reads 


Trinity will generate a directory “trinity_out_dir” by default, and put the results in it. The assembled sequences were in a fasta file “Trinity.fasta”.  Rename it to be “Rhinolophus_sinicus-Unigene.fa”.

2.2 Annotation:
2.2.1 Blastx the assembled transcriptome against the “nr” database.
Download the “nr” database from NCBI: ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz 
Build a local blast index of the “nr” database:

$ makeblastdb -in nr -dbtype prot -parse_seqids


Blastx: search homologous sequences of the assembled transcriptome against the “nr” database:

$ blastx -db nr -query Rhinolophus_sinicus-Unigene.fa -out Rhinolophus_sinicus-Unigene.blastx -evalue 1e-5 -outfmt 5 -max_target_seqs 1 -num_threads 12

Then use biopython to extract the GI numbers and accessions of the best hits of each query sequence. To learn how to use biopython . Please refer to http://biopython.org/DIST/docs/tutorial/Tutorial.html
Note: Since the sequence assembly and blastx take too many compute resources and time, we use the results generated by the NGS sequencing supplier.


2.2.2 Obtain gene names according to GI number by search “gene2accession” database.
Download “gene2accession” database: ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2accession.gz
Decompress it and convert it into a sqlite database:

$ sqlite3   gene2ac.db 	# Enter sqlite3 commandline
sqlite> .mode tabs
sqlite> .import gene2accession gene2ac
sqlite> 					# press Ctl + D to exit

Write a python code to search the sql database and output result:

#!/bin/py 
# Get the gene symbols according to the accession from gene2accession database. 

import re 
import sqlite3 
import datetime 

print "connecting gene2accession database..." 
conn = sqlite3.connect('gene2ac.db') 
c = conn.cursor() 
print "OK!" 

print datetime.datetime.now().time().isoformat() 
fano = open("All_Unigene.basic.annotation.csv", 'rU') 
genes = {} 
gis ={} 
ano_fstline = fano.readline() 
print "Begin annotation..." 
for line in fano: 
    tmp = re.split('\t', line.rstrip()) 
    if tmp[4] == '--': 
        genes[tmp[0]] = ("--", "--",tmp[0],) 
        # {unigeneID:(geneID,ac,symbol)} 
    else: 
        gi = re.split('\|', tmp[4])[1] 
        gis[tmp[0]] = gi 
GI_s = tuple(gis.values()) 

print len(GI_s), "GIs loaded." 
c.execute("""select GeneID,"protein_accession.version",Symbol,protein_gi from gene2ac where protein_gi in ({0})""".format(','.join('?' for _ in GI_s)),GI_s) 
print "Finish search the gene2ac database." 
all_symbols = c.fetchall() 
print len(all_symbols), "searched items." 

print "Construct a dictionary of searched results." 
sch_dict = {} 
for each in all_symbols: 
    sch_dict[each[3]] = each[:3] 

print "Mapping these searched items to the unigene IDs." 
for key in gis.keys(): 
    if gis[key] not in sch_dict.keys(): 
        genes[key] = ("--","--",key) 
    else: 
        if key not in genes: 
            genes[key] = sch_dict[gis[key]] 
fano.close() 
c.close() 
conn.close() 
print "Finish annotation." 

print datetime.datetime.now().time().isoformat() 
print "Writing the results..." 
fout = open("gene_symbol_map.txt", 'w') 
fout.write("UnigeneID\tgeneid\tprotein_accession.version\tgene_symbol\n") 
for each in sorted(genes.keys()): 
    fout.write(each+'\t'+'\t'.join(genes[each])+'\n') 
fout.close() 
print "Finished." 


2.3 RSEM: RNA-seq analysis. Indexing → Mapping →DEG → R calculation
2.3.1 Indexing. Build a reference:

$ rsem-prepare-reference --bowtie --bowtie-path /usr/bin --bowtie2 --bowtie2-path /usr/bin Rhinolophus_sinicus-Unigene.fa rs_unigenes

2.3.2 Mapping and calculation.Write a bash code as below and run it.

#!/bin/bash 
# Map the reads against reference and estimate gene and isoform expression of three Mock cell samples 
rsem-calculate-expression -p 10 --bowtie2 --paired-end WGC072920-RsK-M-1_combined_R1.fastq.gz WGC072920-RsK-M-1_combined_R2.fastq.gz rs_unigenes rsk_M1 
rsem-calculate-expression -p 10 --bowtie2 --paired-end WGC072920-RsK-M-2_combined_R1.fastq.gz WGC072920-RsK-M-2_combined_R2.fastq.gz rs_unigenes rsk_M2 
rsem-calculate-expression -p 10 --bowtie2 --paired-end WGC072920-RsK-M-3_combined_R1.fastq.gz WGC072920-RsK-M-3_combined_R2.fastq.gz rs_unigenes rsk_M3 
# Map the reads against reference and estimate gene and isoform expression of three 6 hpi samples 
rsem-calculate-expression -p 10 --bowtie2 --paired-end WGC072920-RsK-6h-1_combined_R1.fastq.gz WGC072920-RsK-6h-1_combined_R2.fastq.gz rs_unigenes rsk_6h1 
rsem-calculate-expression -p 10 --bowtie2 --paired-end WGC072920-RsK-6h-2_combined_R1.fastq.gz WGC072920-RsK-6h-2_combined_R2.fastq.gz rs_unigenes rsk_6h2 
rsem-calculate-expression -p 10 --bowtie2 --paired-end WGC072920-RsK-6h-3_combined_R1.fastq.gz WGC072920-RsK-6h-3_combined_R2.fastq.gz rs_unigenes rsk_6h3 
# Map the reads against reference and estimate gene and isoform expression of three 12 hpi samples 
rsem-calculate-expression -p 10 --bowtie2 --paired-end WGC072920-RsK-12h-1_combined_R1.fastq.gz WGC072920-RsK-12h-1_combined_R2.fastq.gz rs_unigenes rsk_12h1 
rsem-calculate-expression -p 10 --bowtie2 --paired-end WGC072920-RsK-12h-2_combined_R1.fastq.gz WGC072920-RsK-12h-2_combined_R2.fastq.gz rs_unigenes rsk_12h2 
rsem-calculate-expression -p 10 --bowtie2 --paired-end WGC072920-RsK-12h-3_combined_R1.fastq.gz WGC072920-RsK-12h-3_combined_R2.fastq.gz rs_unigenes rsk_12h3 
# Map the reads against reference and estimate gene and isoform expression of three 24 hpi samples 
rsem-calculate-expression -p 10 --bowtie2 --paired-end WGC072920-RsK-24h-1_combined_R1.fastq.gz WGC072920-RsK-24h-1_combined_R2.fastq.gz rs_unigenes rsk_24h1 
rsem-calculate-expression -p 10 --bowtie2 --paired-end WGC072920-RsK-24h-2_combined_R1.fastq.gz WGC072920-RsK-24h-2_combined_R2.fastq.gz rs_unigenes rsk_24h2 
rsem-calculate-expression -p 10 --bowtie2 --paired-end WGC072920-RsK-24h-3_combined_R1.fastq.gz WGC072920-RsK-24h-3_combined_R2.fastq.gz rs_unigenes rsk_24h3 
# Map the reads against reference and estimate gene and isoform expression of three 36 hpi samples 
rsem-calculate-expression -p 10 --bowtie2 --paired-end WGC072920-RsK-36h-1_combined_R1.fastq.gz WGC072920-RsK-36h-1_combined_R2.fastq.gz rs_unigenes rsk_36h1 
rsem-calculate-expression -p 10 --bowtie2 --paired-end WGC072920-RsK-36h-2_combined_R1.fastq.gz WGC072920-RsK-36h-2_combined_R2.fastq.gz rs_unigenes rsk_36h2 
rsem-calculate-expression -p 10 --bowtie2 --paired-end WGC072920-RsK-36h-3_combined_R1.fastq.gz WGC072920-RsK-36h-3_combined_R2.fastq.gz rs_unigenes rsk_36h3 


2.3.3 Get differentially expressed genes (FDR < 0.05).
Write a bash code as below and run it.

#!/bin/bash 
# Differential expression analysis using EBSeq 

# 6 hpi vs Mock 
rsem-generate-data-matrix rsk_6h1.genes.results rsk_6h2.genes.results
rsk_6h3.genes.results rsk_M1.genes.results rsk_M2.genes.results rsk_M3.genes.results > geneMat_rsk_6hvsM.txt 
rsem-run-ebseq geneMat_rsk_6hvsM.txt 3,3 geneMat_rsk_6hvsM.results 
rsem-control-fdr geneMat_rsk_6hvsM.results 0.05 geneMat_rsk_6hvsM.de.txt 

# 12 hpi vs Mock 
rsem-generate-data-matrix rsk_12h1.genes.results rsk_12h2.genes.results rsk_12h3.genes.results rsk_M1.genes.results rsk_M2.genes.results rsk_M3.genes.results > geneMat_rsk_12hvsM.txt 
rsem-run-ebseq geneMat_rsk_12hvsM.txt 3,3 geneMat_rsk_12hvsM.results 
rsem-control-fdr geneMat_rsk_12hvsM.results 0.05 geneMat_rsk_12hvsM.de.txt 

# 24 hpi vs Mock 
rsem-generate-data-matrix rsk_24h1.genes.results rsk_24h2.genes.results rsk_24h3.genes.results rsk_M1.genes.results rsk_M2.genes.results rsk_M3.genes.results > geneMat_rsk_24hvsM.txt 
rsem-run-ebseq geneMat_rsk_24hvsM.txt 3,3 geneMat_rsk_24hvsM.results 
rsem-control-fdr geneMat_rsk_24hvsM.results 0.05 geneMat_rsk_24hvsM.de.txt 

# 36 hpi vs Mock 
rsem-generate-data-matrix rsk_36h1.genes.results rsk_36h2.genes.results rsk_36h3.genes.results rsk_M1.genes.results rsk_M2.genes.results rsk_M3.genes.results > geneMat_rsk_36hvsM.txt 
rsem-run-ebseq geneMat_rsk_36hvsM.txt 3,3 geneMat_rsk_36hvsM.results 
rsem-control-fdr geneMat_rsk_36hvsM.results 0.05 geneMat_rsk_36hvsM.de.txt 


