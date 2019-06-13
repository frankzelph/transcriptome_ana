#!/usr/bin/Rscript
#
# Usage: Rscript volcano_plot.R DEG_FILE FIG_FILE

# Load packages
library(dplyr)
library(ggplot2)
library(ggrepel)

args <- commandArgs(trailingOnly=TRUE)
DATA_FILE <- args[1]
FIG_FILE <- args[2]


# create a svg figure
svg(filename=paste(FIG_FILE, ".svg", sep=""), width=8, height=7)

res <- read.table(DATA_FILE,sep='\t',header=T)
res$log2FoldChange <- -res$log2FoldChange  # use it accordingly.
res <- mutate(res, class=ifelse((log2FoldChange < -2 & padj < 0.01),"Down", ifelse((log2FoldChange > 2 & padj < 0.01),"Up","Not_sig.")))

# Select some genes to label in the figure
Targets <- c("Cblb","Cd5","Plbd1","Socs1","Pten","Rfwd2","Rasa2","Gtpbp4","Tox")

# Plot
ggplot(res, aes(x=log2FoldChange, y=-log10(padj))) +
	geom_point(aes(colour=class), size=1) +
#	scale_x_continuous(limits = c(-10, 10)) +
	scale_color_manual(values=c("Green","grey60","red")) +
	labs(title="DEG analysis", x="Log2 Fold Change", y="-Log10(Padj)") +
	theme(title = element_text(size=18), axis.text = element_text(size=12),plot.title = element_text(hjust = 0.5)) +
	geom_text_repel(data=filter(res, is.element(Geneid, Targets)), aes(label=Geneid), size=6)
	

dev.off()



