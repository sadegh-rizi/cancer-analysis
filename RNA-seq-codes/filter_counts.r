#
# Filter a count table to remove entries with low expression.
#

# The sample file is in CSV format and must have the headers "sample" and "condition".
sample_file = "samples.csv"

# The count file that contains the summarized counts.
counts_file = "combined_genes.csv"

# The result file name.
output_file = "filtered_counts.csv"

# Inform the user.
print("# Tool: Filter counts")
print(paste("# Sample: ", sample_file))
print(paste("# Counts: ", counts_file))

# Read the sample file
sd <- read.csv(sample_file, stringsAsFactors=F)

# Isolate the sample names.
sn <- sd$sample

# Read the count file.
df = read.csv(counts_file, stringsAsFactors = FALSE)

# Create a matrix from the sample names
dm = as.matrix(df[,sn])

#
# The filtering condition below selects for at least 10 reads across all conditions.
#
min_count = 10

# Apply the filter on the row sums.
keep <- (rowSums(dm) >= min_count) 

# Compute reporting counts.
count_total = length(keep)
count_good = sum(keep)
count_removed = count_total - count_good

# Notify the user on filtering.
print (paste("# Minimum:", min_count, "Total:", count_total, "Kept:", count_good, "Removed:", count_removed))

# Slice the input dataframe with the boolean vector.
df = df[keep,]

# Save the  resulting filtered counts.
write.csv(df, file=output_file,  row.names = FALSE, quote = FALSE)

# Inform the user.
print(paste("# Results:", output_file))

