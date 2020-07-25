
# 1. Load data
rm(list=ls(all=TRUE))
options(stringsAsFactors = FALSE) 
library(dplyr)
library(Seurat)
library(cowplot)  
setwd("C:/Users/TPOW31714/Desktop")
da0<- readRDS("pbmc_test_20200412.rds", refhook = NULL)

# 2. Check data
Ncounts_seurat<- da0[["RNA"]]@data #Normalized counts
dim(Ncounts_seurat)
Ncounts_seurat[1:3,1:3]
table(Idents(da0)) # Check each cluster cell numbers
zz2<- Idents(da0)

# 3. Single R

# 3.1 Seurat data to singleR data
library(SingleR)
Ncounts_singleR<- matrix(Ncounts_seurat, ncol=2700)
dim(Ncounts_singleR)
Cluster_idents<- DataFrame(cluster=Idents(da0),
                     row.names=colnames(Ncounts_seurat))

da0_singleR<- SummarizedExperiment(list(counts=Ncounts_singleR), colData=Cluster_idents)
rownames(da0_singleR)
rownames(da0_singleR)<- rownames(Ncounts_seurat)
rownames(da0_singleR)
da0_singleR

# 3.2 Ref for identifying cell types
# 3.2.1.1 Human Primary Cell Atlas Data
ref_HPCA<- HumanPrimaryCellAtlasData()
ref_HPCA

pred <- SingleR(test = da0_singleR, ref = ref_HPCA, labels =ref_HPCA$label.main,  
        assay.type.test = "counts", method = c("cluster"), clusters=zz2)
pred2 <- SingleR(test = da0_singleR, ref = ref_HPCA, labels =ref_HPCA$label.fine,  
        assay.type.test = "counts", method = c("cluster"), clusters=zz2)
# 3.2.1.2 Check comparing result
pred
dim(pred)
pred[1:2,]
cbind(0:8,pred$labels ,pred2$labels )
cbind(0:8,pred2$labels ) # Exhibit each cluster cell types

pred$scores
pred2$scores # Exhibit possible cell types scores

# 3.2.2 Monaco Immune Data

ref2<- MonacoImmuneData()
ref2
pred2<- SingleR(test = da0_singleR, ref = ref2, labels =ref2$label.main,  
        assay.type.test = "counts", method = c("cluster"),clusters=zz2)
cbind(0:8,pred2$labels )

# 3.2.3 Database Immune Cell Expression Data
ref3<- DatabaseImmuneCellExpressionData()
ref3
pred3<- SingleR(test = da0_singleR, ref = ref3, labels =ref3$label.main,  
        assay.type.test = "counts", method = c("cluster"),clusters=zz2)
cbind(0:8,pred3$labels )
