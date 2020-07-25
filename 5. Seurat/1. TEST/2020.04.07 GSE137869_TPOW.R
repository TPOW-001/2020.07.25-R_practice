
library(dplyr)
library(Seurat)
library(patchwork)
# Load the dataset
BMMY.data <- Read10X(data.dir = "C:/Users/TPOW_NEW/Desktop/2020.04.07 GSE137869 BM/M-Y")
# Initialize the Seurat object with the raw (non-normalized data).
BMMY<- CreateSeuratObject(counts = BMMY.data, project = "GSE137869 BM", min.cells = 3, min.features = 200)
BMMY
# An object of class Seurat 
# 13099 features across 3454 samples within 1 assay 
# Active assay: RNA (13099 features)

dense.size <- object.size(as.matrix(BMMY.data))
dense.size
# 911677264 bytes

sparse.size <- object.size(BMMY.data)
sparse.size
# 60898744 bytes

dense.size/sparse.size
# 15 bytes

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
BMMY[["percent.mt"]] <- PercentageFeatureSet(BMMY, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(BMMY, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(BMMY, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(BMMY, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

BMMY <- subset(BMMY, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

BMMY <- NormalizeData(BMMY, normalization.method = "LogNormalize", scale.factor = 10000)




BMMY <- FindVariableFeatures(BMMY, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(BMMY), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(BMMY)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2




all.genes <- rownames(BMMY)
t0<- Sys.time()
BMMY <- ScaleData(BMMY, features = all.genes)
tn<- Sys.time()
difftime(tn,t0, tz=" ",unit='mins')

# Perform linear dimensional reduction
BMMY <- RunPCA(BMMY, features = VariableFeatures(object = BMMY))

# Examine and visualize PCA results a few different ways
print(BMMY[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(BMMY, dims = 1:2, reduction = "pca")
DimPlot(BMMY, reduction = "pca")


#Determine the ‘dimensionality’ of the dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
BMMY <- JackStraw(BMMY, num.replicate = 100)
BMMY <- ScoreJackStraw(BMMY, dims = 1:20)
JackStrawPlot(BMMY, dims = 1:20)
ElbowPlot(BMMY)

#Cluster the cells
BMMY <- FindNeighbors(BMMY, dims = 1:20)
BMMY <- FindClusters(BMMY, resolution = 1)
# Resolution=1, 有20個cluster
# Resolution=0.5, 有17個cluster


# Run non-linear dimensional reduction (UMAP/tSNE)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
BMMY <- RunTSNE(BMMY, dims = 1:20)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(BMMY, reduction = "tsne")


