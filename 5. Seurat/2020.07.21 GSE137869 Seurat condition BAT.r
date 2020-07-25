
### 0.利用seurat讀取資料並區分好組別 ###
# https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html

rm(list=ls(all=TRUE))
library(dplyr)
library(Seurat)
library(patchwork)
library(cowplot)   #有多的package

setwd("C:/Users/TPOW31714/OneDrive/!!! 2020.04.27 Discuss with Dr.Jiang/2020.04.08 GSE137869 Aging/GSE137869/box")
t0<- Sys.time()
da0<- read.table("BAT_M_comb1.txt",sep="\t",header=T,fill=T, quote = "", check.names=F )
tn<- Sys.time()
difftime(tn,t0, tz=" ",unit='mins')
dim(da0)
da0[1:5,1:7]

#Time difference of 2.875562 mins

data1 <- CollapseSpeciesExpressionMatrix(da0)
data2 <- CreateSeuratObject(counts = data1 )
data2
#An object of class Seurat 
#32884 features across 11991 samples within 1 assay 
#Active assay: RNA (32884 features)

#### mark condition cluster  ##################################
data2$sample<- rep(c("Y","O","CR"), c(3749,4601,4785))

#根據合併時的操作, 先是Y 再來是O 再來是CR, 所以可以知道data中先有O的3456 cells,
#再來是O 的4502 cells, 再來是CR的 4033 cells
#在合併檔案時的, 即上面程式(GSE137869 Aging&CR step-2)中的dim(BM1),dim(BM2),dim(BM3), 可知各條件的cell數

dim(data2)
#32884 genes, 11991 cells
# (Y,Old,CR)==(3456,4502,4033)

table(data2$sample)
#CR    O    Y 
#4033 4502 3456

str(data2)

### 2.跑分condition的seurat ###
# https://satijalab.org/seurat/v3.1/immune_alignment.html

# 2.1 分組別
data2.list <- SplitObject(data2, split.by = "sample")

# 2.2 Normalization and FindVariableFeatures
data2.list <- lapply(X = data2.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3288)
})
# 32840genes, 所以選 nfeatures = 3284

# 2.3 Perform integration
data2.anchors <- FindIntegrationAnchors(object.list = data2.list, dims = 1:50)
data2.combined <- IntegrateData(anchorset = data2.anchors, dims = 1:50)
# 2.4 Perform an integrated analysis
DefaultAssay(data2.combined) <- "integrated"
# 2.4.1 Run the standard workflow for visualization and clustering
data2.combined <- ScaleData(data2.combined, verbose = FALSE)
data2.combined <- RunPCA(data2.combined, npcs = 30, verbose = FALSE)
# 2.4.2 t-SNE and Clustering
data2.combined <- RunUMAP(data2.combined, reduction = "pca", dims = 1:30)
data2.combined <- FindNeighbors(data2.combined, reduction = "pca", dims = 1:30)
data2.combined <- FindClusters(data2.combined, resolution = 0.3)
# PCdims 1:30, resolution = 0.4, Number of communities: 16
# PCdims 1:30, resolution = 0.3, Number of communities: 15
# 2.4.3 saveRDS
saveRDS(data2.combined, file = "GSE137869_BAT_Y&O&CR_split_PC1-30_RS0.3_20200721.rds")

# 2.4.4 loadRDS
rm(list=ls(all=TRUE))
library(dplyr)
library(Seurat)
library(patchwork)
library(cowplot)
setwd("C:/Users/TPOW31714/OneDrive/!!! 2020.04.27 Discuss with Dr.Jiang/2020.04.08 GSE137869 Aging/1. File")
data2.combined<- readRDS("GSE137869_BAT_Y&O&CR_split_PC1-30_RS0.3_20200721.rds", refhook = NULL)
UMAPPlot(data2.combined)

# 2.4.5 Visualization
p1 <- DimPlot(data2.combined, reduction = "umap", group.by = "sample")
p2 <- DimPlot(data2.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

# 2.4.6 use the split.by argument to show each condition colored by cluster.
DimPlot(data2.combined, reduction = "umap", split.by = "sample", label = TRUE)

# 2.5.1 把不同組別的cluster gene叫出來
sub_8<- subset(x = data2.combined, idents =16)

# sub<- subset(x = data2.combined, idents =1)
# sub<- subset(x = data2.combined, idents =c(0,1))

# 2.5.2 Method1 for找出不同組別的cluster genes
zz0<- sub_8[["RNA"]]@data #normalizes data
dim(zz0) #32884genes * 436cells

write.table(zz0,"BAT_cluster16_24_Ncounts_rs0.3_20200720.txt", sep="\t",row.name=T, col.name=T, quote=FALSE)

# write.table(zz0,"BM_Y&O&CR_split cluster8 normalized counts for GSEA_rs0.3_20200716.txt", sep="\t",row.name=F, col.name=T, quote=FALSE)

# 2.5.2.1 Method1 extend: USE CELL ID to select specific treatment conditions
sub_8<- subset(x = data2.combined, idents =16)
dim(sub_8)
# 2000 436
table(sub_8$sample)  #condition composition CR  O   Y 
#250 56 130
cell_id<- colnames(sub_8) #ALL cell ID
xx<- sub_8$sample=="Y"
want<- cell_id[xx]     #Selected CR Only ID
head(want)
length(want)
#[1] 250
sub_8_CR<- subset(x = data2.combined, cells =want)
sub_8_CR
zz3<- sub_8_CR[["RNA"]]@data
dim(zz3)
write.table(zz3,"BAT_cluster16_Y17_Ncounts_rs0.3_20200720.txt", sep="\t",row.name=T, col.name=T, quote=FALSE)


# 2.6 Identify conserved cell type markers
DefaultAssay(data2.combined) <- "RNA"
Conserverd.markers<- FindConservedMarkers(data2.combined, ident.1 = 8, grouping.var = "sample", verbose = FALSE)
head(Conserverd.markers)

# 2.6.1 find markers for every cluster compared to all remaining cells, report only the positive ones
BM.markers <- FindAllMarkers(data2.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(BM.markers,"GSE137869_PC1-20_cluster for markers_rs0.3_2020716.txt", sep="\t",row.name=T, col.name=T, quote=FALSE)

BM.marker_1<- BM.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC)
write.table(BM.marker_1,"GSE137869_cluster for markers_top15_rs0.4_20200427.txt", sep="\t",row.name=T, col.name=T, quote=FALSE)

