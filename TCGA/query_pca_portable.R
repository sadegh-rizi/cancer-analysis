tcga_br_counts <- read.csv("STAD_counts.csv",header=TRUE,row.names = 1)
tcga_br_tpm <- read.csv("STAD_tpm.csv",header=TRUE,row.names = 1)
n<-448
tcga_br_counts = t(tcga_br_counts[!grepl("PAR",rownames(tcga_br_counts)),1:n])
tcga_br_tpm = t(tcga_br_tpm[!grepl("PAR",rownames(tcga_br_tpm)),1:n])



sample_type <- read.csv("sample_type.csv")$x

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

#HC on TPM data
# Creating the distance matrix
tpm.dist <- dist(tcga_br_tpm,method = "euclidean")

#Creating 
tpm.hc <- hclust(d=tpm.dist,method="ward.D2")
fviz_dend(tpm.hc,cex=0.4)


#HC on 4
