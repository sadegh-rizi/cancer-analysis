#
# Compare two RNA seq data analysis result files.
#

# The results files to be compared.
file1 = "results1.csv"
file2 = "results2.csv"

# Read the data files.
data1 <- read.csv(file1)
data2 <- read.csv(file2)

# Select the rows where the FDR is under a cutoff.
sig1 <- subset(data1, FDR <= 0.5)
sig2 <- subset(data2, FDR <= 0.5)

# Extract the gene names only
name1 = sig1$name
name2 = sig2$name

# Intersect: common elements.
isect <- intersect(name1, name2)

# Elements from both.
uni <- union(name1, name2)

# Elements in file 1 only.
only1 <- setdiff(name1, name2)

# Elements in file 2 only.
only2 <- setdiff(name2, name1)

# Report the differences
print("# Tool: compare_results.r")
print("")
print(paste("# File 1:", length(name1), file1))
print(paste("# File 2:", length(name2), file2))
print("")
print(paste("# Union:", length(uni)))
print(paste("# Intersect:", length(isect)))
print(paste("# File 1 only:", length(only1)))
print(paste("# File 2 only:", length(only2)))

print("----")

print("Only 1:")
print(paste( only1))

print("----")

print("Only 2:")
print(paste( only2))

expect = subset(data1, grepl("UP|DOWN", name))$name
