#!/usr/bin/env Rscript
#
# This script makes a pca of all samples and a heatmap of sample distances.
#

# The inputs to the script are counts file, design file and the number of samples.
#
# How to run?
# Rscript summary_plots.r <counts_file> <design_file> <number_of_samples>
# Example: Rscript summary_plots.r counts.txt design.txt 6
# Design file example:
#	Sample	Condition
#	MCVS450	MCVS
#	MCVS515	MCVS
#	MCVS520	MCVS
#	MNS456	MNS
#	MNS486	MNS
#	MNS580	MNS


read <- function(counts_file, design_file,number_of_samples){
  counts = read.table(counts_file, header=TRUE, sep="\t", row.names=1 )
  idx = ncol(counts) - number_of_samples
  
  # Cut out the valid columns.
  if (idx > 0) counts = counts[-c(1:idx)] else counts=counts

  numeric_idx = sapply(counts, mode) == 'numeric'
  counts[numeric_idx] = round(counts[numeric_idx], 0)
  
  colData = read.table(design_file, header=TRUE, sep="\t", row.names=1 )
  
  # Create DESEq2 dataset.
  dds = DESeqDataSetFromMatrix(countData=counts, colData=colData, design = ~1)
  
  # Variance Stabilizing Transformation.
  vsd = vst(dds)
  names = colnames(counts)
  groups = colnames(colData)
  
  rlist <- list("vsd" = vsd, "names"=names, "groups" =groups)
  return(rlist)
  
}

# Command line argument.
args = commandArgs(trailingOnly=TRUE)


if (length(args)!=3) {
  stop("Counts file, Design file and the number of samples must be specified at the commandline", call.=FALSE)
}


# Load the library while suppressing verbose messages.
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(ggplot2))

# Set the plot dimensions.
WIDTH = 12
HEIGHT = 8


# The first argument to the script -counts file
infile = args[1]

# The second  argument to the script - design file
coldata_file = args[2]

# The third argument to the script - total number of samples.
sno = args[3]
sno= as.numeric(sno)

res = read(infile, coldata_file, sno)
vsd= res$vsd
names = res$names
groups = res$groups

# Open the drawing device.
pdf('pca.pdf', width = WIDTH, height = HEIGHT)
par(mfrow = c(2,1))
nudge <- position_nudge(y = 0.5)

z=plotPCA(vsd, intgroup=c(groups)) 
z+ geom_text(aes(label = names), position=nudge, size = 2.5) +ggtitle(aes("PCA"))
dev.off()

#
# Plot heatmap of sample distances
#
library(pheatmap)
library("RColorBrewer")

sampleDists = dist(t(assay(vsd)))
sampleDistMatrix = as.matrix(sampleDists)

colnames(sampleDistMatrix) = NULL

colors = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)


# Open the drawing device.

pdf('heatmap.pdf', width = 8, height = HEIGHT)

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)  + geom_label(aes(label = names))

dev.off()

