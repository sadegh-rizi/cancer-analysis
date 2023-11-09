#
# Combine transcripts into genes with the tximport package.
# 
# https://bioconductor.org/packages/release/bioc/html/tximport.html
#

# Load the packages
library(tximport)

# The directory where the counts for each sample are located.
data_dir <- 'salmon'

# The sample file is in CSV format and must have the headers "sample" and "condition".
design_file = "design.csv"

# What software created the mappings.
method <- "salmon"

# The ouput file name.
output_file = "counts.csv"

# The name of the file that contains transcript to gene mapping.
# See our guide on how to make a mapping file.
tx2gene_file = "tx2gene.csv"

# Inform the user.
print("# Tool: Combine genes")
print(paste("# Sample: ", design_file))
print(paste("# Data dir: ", data_dir))
print(paste("# Gene map: ", tx2gene_file))

# Read the sample file
sample_data <- read.csv(design_file, stringsAsFactors=F)

# Isolate the sample names.
sample_names <- sample_data$sample

# Generate the file names that contain the quantification data.
files <- file.path(data_dir,  sample_names, "quant.sf")

# Read the transcript to gene mapping file.
tx2gene = read.csv(tx2gene_file)

# Summarize over transcripts.
tx <- tximport(files, type=method, tx2gene=tx2gene)

# Transform into a dataframe.
df <- data.frame(tx$counts)

# Set the column names
colnames(df) <- sample_names

# Add a rowname by gene id
df$ensembl_gene_id = row.names(df)

# Create a smaller data frame that connects gene to name.
id2name = tx2gene[, c("ensembl_gene_id", "external_gene_name")]

# Need to de-duplicate so the merge won't create new entries.
id2name = id2name[!duplicated(id2name),]

# Add the gene names to the list
df = merge(df, id2name, by="ensembl_gene_id", all.x = TRUE, all.y = FALSE)

# This will be the new column order.
cols <- c("ensembl_gene_id", "external_gene_name", sample_names)

# Reorganize the columns.
df <- df[, cols]

# Save the  resulting summarized counts.
write.csv(df, file=output_file,  row.names = FALSE, quote = FALSE)

# Inform the user.
print(paste("# Results: ", output_file))
