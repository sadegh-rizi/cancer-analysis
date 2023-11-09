#
# This Makefile perform the Mission Impossible RNA-seq analysis.
#
#
# More info in the Biostar Handbook volume RNA-Seq by Example
#

# The reference genome.
REF = refs/genome.fa

# The file containing transcripts.
TRX = refs/transcripts.fa

# The name of the HISAT2 index.
HISAT2_INDEX = idx/genome

# The name of the salmon index
SALMON_INDEX = idx/transcripts.salmon

# The design file.
DESIGN = design.csv

# These targets are not files.
.PHONY: data align results

# Tell the user to read the source of this Makefile to understand it.
usage:
	@echo "#"
	@echo "# Use the the source, Luke!"
	@echo "#"

# Create the design file and ids.txt file.
${DESIGN}:

	# The design file
	@echo "sample,condition" > design.csv
	@echo "BORED_1,bored" >> design.csv
	@echo "BORED_2,bored" >> design.csv
	@echo "BORED_3,bored" >> design.csv
	@echo "EXCITED_1,excited" >> design.csv
	@echo "EXCITED_2,excited" >> design.csv
	@echo "EXCITED_3,excited" >> design.csv

	# Create the ids.txt file (first column only, delete first line)
	cat design.csv | cut -f1 -d , | sed 1d > ids.txt

# Download the data for the analysis.
data:
	# Download the reference genome.
	wget -nc http://data.biostarhandbook.com/books/rnaseq/data/golden.genome.tar.gz

	# Unpack the reference genome.
	tar xzvf golden.genome.tar.gz

	# Download the sequencing reads.
	wget -nc http://data.biostarhandbook.com/books/rnaseq/data/golden.reads.tar.gz

	# Unpack the sequencing reads.
	tar zxvf golden.reads.tar.gz

# Build the HISAT2 and salmon index for the reference genomes
index: ${REF}
	mkdir -p idx
	hisat2-build ${REF} ${HISAT2_INDEX}
	salmon index -t ${TRX} -i ${SALMON_INDEX}

# Runs a HISAT2 alignment.
align: ${DESIGN} ${HISAT2_INDEX_FILE}
	mkdir -p bam
	cat ids.txt | parallel --progress --verbose "hisat2 -x ${HISAT2_INDEX} -1 reads/{}_R1.fq -2 reads/{}_R2.fq | samtools sort > bam/{}.bam"
	cat ids.txt | parallel -j 1 echo "bam/{}.bam" | \
		xargs featureCounts -p -a refs/features.gff -o counts.txt
	RScript code/parse_featurecounts.r

# Run a SALMON classification.
classify: ${DESIGN} ${SALMON_INDEX}
	mkdir -p salmon
	cat ids.txt | parallel --progress --verbose "salmon quant -i ${SALMON_INDEX} -l A --validateMappings -1 reads/{}_R1.fq -2 reads/{}_R2.fq  -o salmon/{}"
	RScript code/combine_transcripts.r

# Runs an analysis on the aligned data.
results:
	mkdir -p res
	RScript code/deseq2.r
	RScript code/create_heatmap.r

# Required software
install:
	@echo ""
	@echo mamba install wget parallel samtools subread hisat2 salmon bioconductor-tximport bioconductor-edger bioconductor-biomart bioconductor-deseq2 r-gplots
	@echo ""

