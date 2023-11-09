#
#  Obtain gene names for transcripts with the biomaRt package.
#
# https://bioconductor.org/packages/release/bioc/html/biomaRt.html
#

# Load the biomart manager
library(biomaRt)

# The biomaRt dataset name.
dataset <- "drerio_gene_ensembl"

# The ouput file name.
output_file = "tx2gene.csv"

# Make a connection to the validated dataset.
mart <- useEnsembl(dataset = dataset, biomart = 'ensembl')

# The attributes that we want to obtain.
# The first column must match the feature id used during quantification.
attributes <- c(
  "ensembl_transcript_id_version", 
  "ensembl_gene_id",
  "ensembl_transcript_id", 
  "transcript_length",
  "external_gene_name"
)

# Perform the query.
data <- biomaRt::getBM(attributes = attributes, mart = mart)

# Save the data into a file.
write.csv(data, file=output_file, row.names=FALSE, quote=FALSE)

# Inform the user.
print("# Tool: Create tx2gene mapping")
print(paste("# Dataset: ", dataset))
print(paste("# Output: ", output_file))
