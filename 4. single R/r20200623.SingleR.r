# R-3.6.2
rm(list=ls(all=TRUE))
options(stringsAsFactors = FALSE) 

library(dplyr)
library(Seurat)
library(cowplot)  

setwd("F:/zDatabase/SingleCell/Seurat/testdata")
da<- readRDS("pbmc_test_20200412.rds", refhook = NULL)

zz<- da[["RNA"]]@data
dim(zz)
zz[1:3,1:3]
table(Idents(da))
zz2<- Idents(da)




library(SingleR)

counts<- matrix(zz, ncol=2700)
dim(counts)

colData<- DataFrame(cluster=Idents(da),
                     row.names=colnames(zz))

testz<- SummarizedExperiment(list(counts=counts), colData=colData)
rownames(testz)
rownames(testz)<- rownames(zz)
testz



ref<- HumanPrimaryCellAtlasData()
ref

pred <- SingleR(test = testz, ref = ref, labels =ref$label.main,  
        assay.type.test = "counts", method = c("cluster"), clusters=zz2)
pred2 <- SingleR(test = testz, ref = ref, labels =ref$label.fine,  
        assay.type.test = "counts", method = c("cluster"), clusters=zz2)

pred
dim(pred)
pred[1:2,]

cbind(0:8,pred$labels ,pred2$labels )
cbind(0:8,pred2$labels )


pred$scores

ref2<- MonacoImmuneData()
ref2

pred2<- SingleR(test = testz, ref = ref2, labels =ref2$label.main,  
        assay.type.test = "counts", method = c("cluster"),clusters=zz2)
cbind(0:8,pred2$labels )




ref3<- DatabaseImmuneCellExpressionData()
ref3
pred3<- SingleR(test = testz, ref = ref3, labels =ref3$label.main,  
        assay.type.test = "counts", method = c("cluster"),clusters=zz2)
cbind(0:8,pred3$labels )
