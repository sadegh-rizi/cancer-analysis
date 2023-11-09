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
BiocManager::install("ComplexHeatmap")

BiocManager::install("ConsensusClusterPlus")


BiocManager::install("TCGAbiolinks")
browseVignettes("TCGAbiolinks")

library(TCGAbiolinks)
library(dplyr)
library(DT)
library(DESeq2)
packageVersion("TCGAbiolinks")

getProjectSummary('TCGA-STAD')

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


query <- GDCquery(
  project = c("TCGA-STAD"),#Stomach
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  access="open",
  experimental.strategy = "RNA-Seq",
  barcode=c('TCGA-FP-7735-11A-01R-2055-13',
            'TCGA-FP-7735-01A-11R-2055-13',
            'TCGA-BR-6457-01A-21R-1802-13',
            'TCGA-BR-6457-11A-01R-1802-13',
            'TCGA-HU-8238-01A-11R-2343-13',
            'TCGA-HU-8238-11A-01R-2343-13')
)

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
br_Cor <-  TCGAanalyze_Preprocessing(tcga_br)
br_Cor


# normalization of genes
dataNorm <- TCGAanalyze_Normalization(
  tabDF = br_Cor, 
  geneInfo =  geneInfoHT
)

# quantile filter of genes
# dataFilt <- TCGAanalyze_Filtering(
#   tabDF = dataNorm,
#   method = "quantile", 
#   qnt.cut =  0.25
# )
dataFilt <- t(dataNorm %>% 
  TCGAanalyze_Filtering(method = "varFilter") %>%
  TCGAanalyze_Filtering(method = "filter1") %>%  
  TCGAanalyze_Filtering(method = "filter2",foldChange = 1))
#dataFilt=dataFilt[c(sample(1:5000,2000)),]
# selection of normal samples "NT"
samplesNT <- TCGAquery_SampleTypes(
  barcode = colnames(dataFilt),
  typesample = c("NT")
)

# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(
  barcode = colnames(dataFilt), 
  typesample = c("TP")
)

# Diff.expr.analysis (DEA)
dataDEGs <- TCGAanalyze_DEA(
  mat1 = dataFilt[,samplesNT],
  mat2 = dataFilt[,samplesTP],
  metadata=TRUE,
  pipeline = "edgeR",
  
  Cond1type = "Normal",
  Cond2type = "Tumor",
  fdr.cut = 0.05 ,
  logFC.cut = 1,
  method = "glmLRT"
)


# DEGs table with expression values in normal and tumor samples
dataDEGsFiltLevel <- TCGAanalyze_LevelTab(
  FC_FDR_table_mRNA = dataDEGs,
  typeCond1 = "Tumor",
  typeCond2 = "Normal",
  TableCond1 = dataFilt[,samplesTP],
  TableCond2 = dataFilt[,samplesNT]
)


  


#BiocManager::install("EnsDb.Hsapiens.v86")
library(EnsDb.Hsapiens.v86)

Genelist <- rownames(dataDEGsFiltLevel)
geneIDs <- ensembldb::select(EnsDb.Hsapiens.v86, keys= Genelist, keytype = "GENEID", columns = c("SYMBOL","GENEID") )

ansEA <- TCGAanalyze_EAcomplete(
  TFname = "DEA genes Normal Vs Tumor",
  RegulonList = geneIDs$SYMBOL
)

# Enrichment Analysis EA (TCGAVisualize)
# Gene Ontology (GO) and Pathway enrichment barPlot

TCGAvisualize_EAbarplot(
  tf = rownames(ansEA$ResBP), 
  GOBPTab = ansEA$ResBP,
  GOCCTab = ansEA$ResCC,
  GOMFTab = ansEA$ResMF,
  PathTab = ansEA$ResPat,
  nRGTab = Genelist, 
  nBar = 10
)

#It crashes!!
data_Hc2 <- TCGAanalyze_Clustering(
  tabDF = dataFilt,method='consensus',methodHC = 'ward.D2')



colData(tcga_br)$groupsHC <- paste0("EC",data_Hc2[[4]]$consensusClass)



#PCA
rownames(dataFilt)<- geneIDs$SYMBOL

pca <- TCGAvisualize_PCA(
  dataFilt = dataFilt,
  dataDEGsFiltLevel = dataDEGsFiltLevel,
  ntopgenes = 200, 
  group1 = samplesNT,
  group2 =  samplesTP
)

