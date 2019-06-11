#!/bin/Rscript

library(dplyr)
require(org.Mm.eg.db)

args <- commandArgs(trailingOnly=TRUE)
infiles <- args

for (infile in infiles) {
deg <- read.delim(file=infile)
id <- rownames(deg)	

deg <- data.frame(ENSEMBL=id, deg)

gene_symbol <- AnnotationDbi::select(org.Mm.eg.db, keys=id, columns=c("SYMBOL","ENSEMBL","GENENAME"), keytype="ENSEMBL")
gene_symbol <- subset(gene_symbol, !is.na(SYMBOL))
# join the gene symbol to the dif. expression gene table
final_deg <- left_join(gene_symbol, deg, by="ENSEMBL")
final_deg <- final_deg[, c(1,4:(dim(final_deg)[2]),2,3)]

outfile <- paste(infile, ".Anno",sep="")

write.table(final_deg, file=outfile, sep='\t', quote=F, row.names=F)

}
