library(dplyr)
library(ggplot2)
library("factoextra")
load("Last.RData")

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
pve <- data.frame(x=seq(1,10),y=pve[1:10])

#Principal components
barplot(pve[1:10])
ggplot(pve,aes(x=factor(x),y=y))+
  geom_bar(stat="identity")+
  labs(x = "PC components", y = "PVE",title="Proportion of variance explaind")
ggsave("pve_pc.png")


perform_pca <- function(data,suffix){
  pca <- prcomp(data)
  pca.x <- as.data.frame(pca$x)
  pca.var <- pca$sdev ^2
  pve <- signif(100* pca.var / sum(pca.var),2)
  pve <- data.frame(x=seq(1,10),y=pve[1:10])
  
  ggplot(pve,aes(x=factor(x),y=y))+
    geom_bar(stat="identity")+
    labs(x = "PC components", y = "PVE",title="Proportion of variance explaind")
  ggsave(paste0("pca_pve_",suffix,".png"))
  #Tumor and Normal PCA
  ggplot(pca.x)+geom_point(aes(PC1,PC2,color=sample_type))+ggtitle(paste0("Normal and Tumor -",suffix))
  ggsave(paste0("pca_tumor_normal_",suffix,".png"),width = 7,height = 5)
  
  
  #Normal PCA
  ggplot(pca.x[sample_type=="Normal",])+geom_point(aes(PC1,PC2))+ggtitle(paste0("Normal Cases -",suffix))
  ggsave(paste0("pca_normal_",suffix,".png"),width = 7,height = 5)
  
  #Tumor PCA
  ggplot(pca.x[sample_type=="Tumor",])+geom_point(aes(PC1,PC2))+ggtitle(paste0("Tumor Cases -",suffix))
  ggsave(paste0("pca_tumor_",suffix,".png"),width = 7,height = 5)
  
  ggplot(pca.x[sample_type=="Tumor",])+geom_point(aes(PC1,PC2,color=sample_data[sample_type=="Tumor","ajcc_pathologic_stage"]))+
    ggtitle(paste0("Tumor Cases- ",suffix))
  ggsave(paste0("pca_tumor_stage_",suffix,".png"),width = 7,height = 5)
  invisible(pca.x)
}

perform_pca(tcga_br_counts,"counts")
perform_pca(tcga_br_tpm,"tpm")
perform_pca(log2(1+tcga_br_tpm),"tpm_log")
 perform_hclust <- function(data,suffix){
  
}  
  
  
#Tumor and Normal PCA
ggplot(pca.x)+geom_point(aes(PC1,PC2,color=sample_type))+ggtitle("Normal and Tumor")
ggsave("tumor_normal_pca.png",width = 7,height = 5)


#Normal PCA
ggplot(pca.x[sample_type=="Normal",])+geom_point(aes(PC1,PC2))+ggtitle("Normal Cases")
ggsave("normal_pca.png",width = 7,height = 5)

#Tumor PCA
ggplot(pca.x[sample_type=="Tumor",])+geom_point(aes(PC1,PC2))+ggtitle("Tumor Cases")
ggsave("tumor_pca.png",width = 7,height = 5)

ggplot(pca.x[sample_type=="Tumor",])+geom_point(aes(PC1,PC2,color=sample_data[sample_type=="Tumor","ajcc_pathologic_stage"]))+
  ggtitle("Tumor Cases")
ggsave("tumor_stage_pca.png",width = 7,height = 5)
  





#Hierarchial clustering

#HC on TPM data
# Creating the distance matrix for tumors
tcga_br_tpm_tumor <- tcga_br_tpm[sample_type=="Tumor",]

tpm.dist <- dist(tcga_br_tpm_tumor,method = "euclidean")

#Creating 


tpm.hc <- hclust(d=tpm.dist,method="average")

fviz_dend(tpm.hc,cex=0.1)
tpm.coph <- cophenetic(tpm.hc)
cor(tpm.dist,tpm.coph)
# Choosing the linkage method -----------
#Checking different methods for the highest correlation.
methods <- c("complete", "single", "average", "mcquitty", "median", "centroid", "ward.D", "ward.D2")
hc_list <- list()
cor_list <- list()
for (method in methods){
  hc <- hclust(d=tpm.dist,method=method)
  hc.coph <- cophenetic(hc)
  hc_list[[method]]<-hc
  cor_list[[method]]<-cor(tpm.dist,hc.coph)
 
}


print(methods)
print(hc_list)
print(cor_list) 

fviz_dend(tpm.hc,cex=0.1)
fviz_dend(hc_list$average,cex=0.1)

#Average method correlation : 0.7494675
# It means that in this method distances in the height of
#dendrogram better represent the distance matrix values


# Removing outliers from clusters ------------------------
# different heights(cutoff_values should be tested)
#Here a cut_off value was choosen to remove 9 outliers 
cut_off= 95000
fviz_dend(tpm.hc,h=cut_off, # Cut in 10 groups
          cex = 0.1, # label size
          
          color_labels_by_k = TRUE, # color labels by groups
          rect = TRUE # Add rectangle around groups
)
groupindexes <- cutree(tpm.hc, h =cut_off)
tpm.dist_pruned <-dist(tcga_br_tpm_tumor[(groupindexes==1),],method = "euclidean")
tpm.hc_pruned <- hclust(d=tpm.dist_pruned,method="average")
fviz_dend(tpm.hc_pruned,cex=0.1)

# Choosing number of clusters --------------------------
# number of clusters ranging from 2 to 25 were tested
k <- 2:25
for(i in k){
  fviz_dend(tpm.hc_pruned, k= i, # Cut in i groups
            cex = 0.3, # label size
            
            color_labels_by_k = TRUE, # color labels by groups
            rect = TRUE # Add rectangle around groups
  )
  ggsave(paste0("hclust_",i,".png"),width = 20,height = 10)
}

fviz_dend(hc_list$average, k= 10, # Cut in 10 groups
          cex = 0.1, # label size
          
          color_labels_by_k = TRUE, # color labels by groups
          rect = TRUE # Add rectangle around groups
)


# ggsave("test.png",width = 20,height = 10)
# 
# ggsave("hclust_tumor.png",width = 20,height = 10)







