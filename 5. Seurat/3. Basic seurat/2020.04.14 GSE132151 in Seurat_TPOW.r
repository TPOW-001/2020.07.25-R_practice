rm(list=ls(all=TRUE))


##########讀取GSE132151 data############################################
library(dplyr)
library(Seurat)
library(patchwork)
setwd("C:/Users/TPOW_NEW/Desktop/2020.03.31 GSE132151/Data edited")
t0<- Sys.time()

data0<- read.table("GSE132151_2847cells_20200318t.txt",sep="\t",header=T,fill=T, quote = "", check.names=F )

tn<- Sys.time()
difftime(tn,t0, tz=" ",unit='mins')
# Time difference of 3.691845 mins

names(data0)[1:25]
rownames(data0)[1:25]
dim(data0)
da0<- data0[c(19:25307),]
########################################################################

data1 <- CollapseSpeciesExpressionMatrix(da0)
data2 <- CreateSeuratObject(counts = data1 )
data2
#An object of class Seurat 
#25289 features across 2847 samples within 1 assay 
#Active assay: RNA (25289 features)

str(data2)
############################################################################
#Formal class 'Seurat' [package "Seurat"] with 12 slots
# ..@ assays      :List of 1
#  .. ..$ RNA:Formal class 'Assay' [package "Seurat"] with 8 slots
#  .. .. .. ..@ counts       :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#  .. .. .. .. .. ..@ i       : int [1:2256114] 52 58 501 527 573 624 668 744 797 898 ...
#  .. .. .. .. .. ..@ p       : int [1:2848] 0 1421 2484 3283 3809 4598 5860 6625 7396 #8593 ...
#  .. .. .. .. .. ..@ Dim     : int [1:2] 25289 2847
#  .. .. .. .. .. ..@ Dimnames:List of 2
#  .. .. .. .. .. .. ..$ : chr [1:25289] "X0610006L08Rik" "X0610007P14Rik" "X0610009B22Rik#" "X0610009E02Rik" ...
#  .. .. .. .. .. .. ..$ : chr [1:2847] "2" "3" "4" "5" ...
#  .. .. .. .. .. ..@ x       : num [1:2256114] 1 1 1 1 1 1 1 1 1 1 ...
#  .. .. .. .. .. ..@ factors : list()
#  .. .. .. ..@ data         :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#  .. .. .. .. .. ..@ i       : int [1:2256114] 52 58 501 527 573 624 668 744 797 898 ...
#  .. .. .. .. .. ..@ p       : int [1:2848] 0 1421 2484 3283 3809 4598 5860 6625 7396 #8593 ...
#  .. .. .. .. .. ..@ Dim     : int [1:2] 25289 2847
#  .. .. .. .. .. ..@ Dimnames:List of 2
#  .. .. .. .. .. .. ..$ : chr [1:25289] "X0610006L08Rik" "X0610007P14Rik" "X0610009B22Rik#" "X0610009E02Rik" ...
#  .. .. .. .. .. .. ..$ : chr [1:2847] "2" "3" "4" "5" ...
#  .. .. .. .. .. ..@ x       : num [1:2256114] 1 1 1 1 1 1 1 1 1 1 ...
#  .. .. .. .. .. ..@ factors : list()
#  .. .. .. ..@ scale.data   : num[0 , 0 ] 
#  .. .. .. ..@ key          : chr "rna_"
#  .. .. .. ..@ assay.orig   : NULL
#  .. .. .. ..@ var.features : logi(0) 
#  .. .. .. ..@ meta.features:'data.frame':      25289 obs. of  0 variables
#  .. .. .. ..@ misc         : NULL
#  ..@ meta.data   :'data.frame':        2847 obs. of  3 variables:
#  .. ..$ orig.ident  : Factor w/ 1 level "SeuratProject": 1 1 1 1 1 1 1 1 1 1 ...
#  .. ..$ nCount_RNA  : num [1:2847] 3068 1972 1350 1113 1869 ...
#  .. ..$ nFeature_RNA: int [1:2847] 1421 1063 799 526 789 1262 765 771 1197 984 ...
#  ..@ active.assay: chr "RNA"
#  ..@ active.ident: Factor w/ 1 level "SeuratProject": 1 1 1 1 1 1 1 1 1 1 ...
#  .. ..- attr(*, "names")= chr [1:2847] "2" "3" "4" "5" ...
#  ..@ graphs      : list()
#  ..@ neighbors   : list()
#  ..@ reductions  : list()
#  ..@ project.name: chr "SeuratProject"
#  ..@ misc        : list()
#  ..@ version     :Classes 'package_version', 'numeric_version'  hidden list of 1
#  .. ..$ : int [1:3] 3 1 3
#  ..@ commands    : list()
#  ..@ tools       : list()
############################################################################


#### standard log-normalization
data2 <- NormalizeData(data2)

str(data2)
data2[["RNA"]]@data[c(502,2590:2592),1]
data2$RNA@counts[c(502,2590:2592),1]
sum(data2$RNA@counts[,1])
zz<- data2$RNA@counts 
log((1/3068)*10000)
log1p ((1/3068)*10000)
#log1p(x) computes log(1+x) accurately also for |x| << 1.


####Identification of highly variable features (feature selection)
data2 <- FindVariableFeatures(data2, nfeatures = 2000)

# Identify the 20 most highly variable genes
top20 <- head(VariableFeatures(data2), 20)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(data2)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)

#plot1 + plot2
CombinePlots(plots = list(plot1, plot2))


#### Scaling the data
# Shifts the expression of each gene, so that the mean expression across cells is 0
# Scales the expression of each gene, so that the variance across cells is 1 
#   This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
#
#data[["RNA"]]@scale.data

data2 <- ScaleData(data2)

data2[["RNA"]]@data[1:3,1:3]
data2[["RNA"]]@scale.data[1:3,1:3]

#### Perform linear dimensional reduction
data2<- RunPCA(data2, features = VariableFeatures(object = data2))
print(data2[["pca"]], dims = 1:8, nfeatures = 10)

VizDimLoadings(data2, dims = 1:2, reduction = "pca")

DimPlot(data2, reduction = "pca")
DimHeatmap(data2, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(data2, dims = 1:3, cells = 500, balanced = TRUE)


#### Determine the ??dimensionality?? of the dataset
data2 <- JackStraw(data2, num.replicate = 100)
data2 <- ScoreJackStraw(data2, dims = 1:20)

JackStrawPlot(data2, dims = 1:15)
ElbowPlot(data2)


#### Cluster the cells
data2 <- FindNeighbors(data2, dims = 1:15)
data2 <- FindClusters(data2, resolution = 0.5)
#resolution=0.4, Number of communities:6
#resolution=0.5, Number of communities:7
#resolution=1, Number of communities:12

head(Idents(data2), 5)
table(Idents(data2))


#### Run non-linear dimensional reduction (UMAP/tSNE)
data2 <- RunUMAP(data2, dims = 1:15)
DimPlot(data2, reduction = "umap")
UMAPPlot(data2)


data2 <- RunTSNE(data2, dims = 1:10)
DimPlot(data2, reduction = "tsne")
TSNEPlot(data2)


setwd("C:/Users/TPOW_NEW/Desktop/2020.04.14 GSE132151")
saveRDS(data2, file = "GSE132151_seurat_20200414.rds")

###############################################################################
###############################################################################
rm(list=ls(all=TRUE))
library(dplyr)
library(Seurat)
library(patchwork)
setwd("C:/Users/TPOW_NEW/Desktop/2020.04.14 GSE132151/1. File")
data2<- readRDS("GSE132151_seurat_20200414_rs=1.rds", refhook = NULL)
UMAPPlot(data2)


#### Finding differentially expressed features (cluster biomarkers)
# find all markers of cluster 1
cluster1.markers <- FindMarkers(data2, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(data2, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
BM.markers <- FindAllMarkers(data2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

write.table(BM.markers,"GSE132151_cluster for markers_rs1_20200414.txt", sep="\t",row.name=F, col.name=T, quote=FALSE)

BM.marker_1<- BM.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC)

write.table(BM.marker_1,"GSE132151_cluster for markers_top15_rs1_20200414.txt", sep="\t",row.name=F, col.name=T, quote=FALSE)

cluster0.markers <- FindMarkers(data2, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)


#Visualizing marker: VlnPlot
VlnPlot(data2, features = c("Nt5e", "Thy1"))

#Visualizing marker: FeaturePlot

FeaturePlot(data2, features = c("Ms4a1","Fos"))
FeaturePlot(data2, features = c("Nt5e","Thy1","Eng","Itgav","Ly6a","Nes","Pdgfra"),reduction="umap")
FeaturePlot(data2, features = c("Cebpb","Pparg","Klf5","Ppargc1a","Prdm16","Prkaa1","Runx2","Sp7","Alpl"),reduction="umap")
FeaturePlot(data2, features = c("Bglap","Spp1","Fabp4","Igfbp4","Sox9","Acan","Col2A1"),reduction="umap")
FeaturePlot(data2, features = c("Jag1","Tnfaip6","Ccl2","Hmox1","Sox9","Acan","Col2a1"),reduction="umap")
FeaturePlot(data2, features = c("Tgfb1","Hgf","Fgf7","Igf1","Sox9","Acan","Col2a1"),reduction="umap")
FeaturePlot(data2, features = c("Lgals1","Lgals3","Lgals9","Entpd1","Sox9","Acan","Col2a1"),reduction="umap")
FeaturePlot(data2, features = c("Il1a","Il1b","Il1rn","Il6","Sox9","Acan","Col2a1"),reduction="umap")
FeaturePlot(data2, features = c("Ptgs1","Ptgs2","Ptges3","Il6","Sox9","Acan","Col2a1"),reduction="umap")
FeaturePlot(data2, features = c("Icam1","Vcam1","Ido1","Arg1","Sox9","Acan","Col2a1"),reduction="umap")
FeaturePlot(data2, features = c("Faslg","H2-M3","H2-T23","Arg1","Sox9","Acan","Col2a1"),reduction="umap")





#default umap, then tsne, then pca

#Visualizing marker: DoHeatmap
top20 <- BM.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
DoHeatmap(data2,features = top20$gene) + NoLegend()

all.cells<- colnames(data2)
#all cell names

length(all.cells)
#check number of cell: 2487 cell names

want_cell<- all.cells[Idents(data2)==0]
#擷取出屬第0群的cell names

want_cell<- all.cells[Idents(data2) %in% c(0,1,2)]
#擷取出屬第0,1,2群的cell names

DoHeatmap(data2, features = c("Car3","Fos"), cells=want_cell, label=FALSE) + NoLegend()
DoHeatmap(data2, features = top20$gene, cells=want_cell, label=FALSE) + NoLegend()


# label=FALSE: 預設是TRUE, 會在圖上標群別0,1,2..., 現在只畫第0群的話, 還是會標0,1,2..., 所以乾脆不要標了



#如何只選一個group???

write.table(top20,"GSE132151_cluster for markers_top20_rs1_20200414.txt", sep="\t",row.name=F, col.name=T, quote=FALSE)



#Assigning cell type identity to clusters
new.cluster.ids <- c("0", "1", "2", "3", "4", "5", 
                     "6", "7", "8", "9", "10", "11")
names(new.cluster.ids) <- levels(data2)
pbmc <- RenameIdents(data2, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


