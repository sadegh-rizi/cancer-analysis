# 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.17")
# 

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManSeq2")
# ager")
# 
# BiocManager::install("DE
# pkgs = c( "edgeR", "gplots", "tximport", "biomaRt")
# BiocManager::install(pkgs)


BiocManager::install("TCGAbiolinks")
browseVignettes("TCGAbiolinks")

library(TCGAbiolinks)
library(dplyr)
library(DT)
library(DESeq2)
library(ggplot2)
library("factoextra")
packageVersion("TCGAbiolinks")

getProjectSummary('TCGA-STAD')
#Query for stomach tumors RNA-seq data
query <- GDCquery(
  project = c("TCGA-STAD"),#Stomach
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  access="open",
  experimental.strategy = "RNA-Seq",
  
  
)



datatable(
  getResults(query), 
  filter = 'top',
  options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
  rownames = FALSE
)

# Choosing specific barcodes

# query <- GDCquery(
#   project = c("TCGA-STAD"),#Stomach
#   data.category = "Transcriptome Profiling",
#   data.type = "Gene Expression Quantification",
#   access="open",
#   experimental.strategy = "RNA-Seq",
#   barcode=c('TCGA-FP-7735-11A-01R-2055-13',
#             'TCGA-FP-7735-01A-11R-2055-13',
#             'TCGA-BR-6457-01A-21R-1802-13',
#             'TCGA-BR-6457-11A-01R-1802-13',
#             'TCGA-HU-8238-01A-11R-2343-13',
#             'TCGA-HU-8238-11A-01R-2343-13')


datatable(
  getResults(query), 
  filter = 'top',
  options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
  rownames = FALSE
)



getResults(query)
GDCdownload(query)

tcga_br<-GDCprepare(query,summarizedExperiment = TRUE)
tcga_br_counts <- assay(tcga_br,'unstranded')  
tcga_br_tpm <- assay(tcga_br,'tpm_unstrand') 


n<-448
sample_data<- as.data.frame(colData(tcga_br))[1:n,]
datatable(sample_data, 
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE
)

gene_data <- as.data.frame(rowRanges(tcga_br))
datatable(gene_data, 
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE
)


treatments<-as.data.frame(bind_rows(sample_data[,"treatments"]))

datatable(treatments
,  options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)


#Saving the files
write.csv(tcga_br_counts,"STAD_counts.csv")
write.csv(tcga_br_tpm,"STAD_tpm.csv")
write.csv(gene_data,"gene_data.csv")
write.csv(treatments,"treatments.csv")

  #Excluding PAR genes : 44 genes 
#Working with 20. samples for the start
#Transposing such that the samples are in rows and genes in columns
tcga_br_counts = t(tcga_br_counts[!grepl("PAR",rownames(tcga_br_counts)),1:n])
tcga_br_tpm = t(tcga_br_tpm[!grepl("PAR",rownames(tcga_br_tpm)),1:n])



# selection of normal samples "NT"
samplesNT <- TCGAquery_SampleTypes(
  barcode = rownames(tcga_br_tpm),
  typesample = c("NT")
)

# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(
  barcode = rownames(tcga_br_tpm), 
  typesample = c("TP")
)
sample_type=ifelse(rownames(tcga_br_tpm) %in% samplesNT,"Normal","Tumor")
sample_type

write.csv(sample_type,"sample_type.csv")
pca <- prcomp(tcga_br_tpm)
plot(pca)
pca.x <- as.data.frame(pca$x)
pca.var <- pca$sdev ^2
pve <- signif(100* pca.var / sum(pca.var),2)
pve

#Principal components
barplot(pve[1:10])
ggsave("pve_pc.png")



#Tumor and Normal PCA
ggplot(pca.x)+geom_point(aes(PC1,PC2,color=sample_type))
ggsave("tumor_normal_pca.png",width = 7,height = 5)


#Normal PCA
ggplot(pca.x[sample_type=="Normal",])+geom_point(aes(PC1,PC2))
ggsave("normal_pca.png",width = 7,height = 5)

#Tumor PCA
ggplot(pca.x[sample_type=="Tumor",])+geom_point(aes(PC1,PC2))
ggsave("tumor_pca.png",width = 7,height = 5)

ggplot(pca.x[sample_type=="Tumor",])+geom_point(aes(PC1,PC2,color=sample_data[sample_type=="Tumor","ajcc_pathologic_stage"]))
ggsave("tumor_stage_pca.png",width = 7,height = 5)


#Hierarchial clustering

# Creating the distance matrix
tpm.dist <- dist(tcga_br_tpm,method = "euclidean")

#Creating 
tpm.hc <- hclust(d=tpm.dist,method="ward.D2")
fviz_dend(tpm.hc,cex=0.4)

