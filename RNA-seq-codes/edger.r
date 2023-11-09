#
# Differential expression analysis with the edgeR package.
# 
# https://bioconductor.org/packages/release/bioc/html/edgeR.html
#

# Load the library
library(edgeR)

# The name of the file that contains the counts.
counts_file = "counts.csv"

# The sample file is in CSV format and must have the headers "sample" and "condition".
design_file = "design.csv"

# The final result file.
output_file = "results.csv"

# Read the sample file.
colData <- read.csv(design_file, stringsAsFactors=F)

# Turn conditions into factors.
colData$condition = factor(colData$condition)

# The first level should correspond to the first entry in the file!
# Required when building a model.
colData$condition = relevel(colData$condition, toString(colData$condition[1]))

# Isolate the sample names.
sample_names <- colData$sample

# Read the data from the standard input.
df = read.csv(counts_file, header=TRUE, row.names=1 )

# Created rounded integers for the count data
counts = round(df[, sample_names])

# Other columns in the dataframe that are not sample information. 
otherCols = df[!(names(df) %in% sample_names)]

# Using the same naming as in the library.
group <- colData$condition

# Creates a DGEList object from a table of counts and group.
dge <- DGEList(counts=counts, group=group)

# Maximizes the negative binomial conditional common likelihood to estimate a common dispersion value across all genes.
dis <- estimateCommonDisp(dge)

# Estimates tagwise dispersion values by an empirical Bayes method based on weighted conditional maximum likelihood.
tag <- estimateTagwiseDisp(dis)

# Compute genewise exact tests for differences in the means between the groups.
etx <- exactTest(tag)

# Extracts the most differentially expressed genes.
etp <- topTags(etx, n=nrow(counts))

# Get the scale of the data
scale = dge$samples$lib.size * dge$samples$norm.factors

# Get the normalized counts
normed = round(t(t(counts)/scale) * mean(scale))

# Turn the DESeq2 results into a data frame.
data = merge(otherCols, etp$table, by="row.names")

# Create column placeholders.
data$baseMean = 1
data$baseMeanA = 1
data$baseMeanB = 1
data$foldChange = 2 ^ data$logFC
data$falsePos = 1

# Rename the column.
names(data)[names(data)=="logFC"] <-"log2FoldChange"

# Compute the adjusted p-value
data$PAdj = p.adjust(data$PValue, method="hochberg")

# Rename the first columns for consistency with other methods.
colnames(data)[1] <- "name"

# Create a merged output that contains the normalized counts.
total <- merge(data, normed, by.x='name', by.y="row.names")

# Sort the data for the output. 
total = total[with(total, order(PValue, -foldChange)), ]

# Compute the false discovery counts on the sorted table.
data$falsePos = 1:nrow(data) * data$FDR

# Sample names for condition A
col_names_A = data.frame(split(colData, colData$condition)[1])[,1]

# Sample names for condition B
col_names_B = data.frame(split(colData, colData$condition)[2])[,1]

# Create the individual baseMean columns.
total$baseMeanA = rowMeans(total[, col_names_A])
total$baseMeanB = rowMeans(total[, col_names_B])
total$baseMean = total$baseMeanA + total$baseMeanB

# Round the numbers
total$foldChange = round(total$foldChange, 3)
total$FDR = round(total$FDR, 4)
total$PAdj = round(total$PAdj, 4)
total$logCPM = round(total$logCPM, 1)
total$log2FoldChange = round(total$log2FoldChange, 1)
total$baseMean = round(total$baseMean, 1)
total$baseMeanA = round(total$baseMeanA, 1)
total$baseMeanB = round(total$baseMeanB, 1)
total$falsePos = round(total$falsePos, 0)

# Reformat these columns as string.
total$PAdj = formatC(total$PAdj, format = "e", digits = 1)
total$PValue = formatC(total$PValue, format = "e", digits = 1)

# Reorganize columns names to make more sense.
new_cols = c("name", names(otherCols), "baseMean","baseMeanA","baseMeanB","foldChange",
         "log2FoldChange","logCPM","PValue","PAdj", "FDR","falsePos",col_names_A, col_names_B)

# Slice the dataframe with new columns.
total = total[, new_cols]

# Write the result to the standard output.
write.csv(total, file=output_file, row.names=FALSE, quote=FALSE)

# Inform the user.
print("# Tool: edgeR")
print(paste("# Design: ", design_file))
print(paste("# Input: ", counts_file))
print(paste("# Output: ", output_file))

