# MCA
# https://github.com/ggjlab/scMCA

# 1. Load data
rm(list=ls(all=TRUE))
options(stringsAsFactors = FALSE) 
library(scMCA) 
setwd("C:/Users/TPOW31714/Desktop")
da0<- read.table("BAT_cluster16_24_Ncounts_rs0.3_20200720.txt", header=T, sep="\t")
data(da0)
dim(da0)
#[1] 32884   25
# 32884 genes expression value of 25 cells

# 2. Run MCA
# scMCA has two parameters , single cell expression matrix(scdata) and 
# the number of most similar cell types
mca_result <- scMCA(scdata = da0, numbers_plot = 10)

# The return of scMCA() is a list which contains 4 parts.
# cors_matrix: Pearson correlation coefficient matrix of each cell and cell type.
# top_cors: equals to numbers_plot
# scMCA: the most relevant cell type for each query cell
# scMCA_probility: the top n relevant cell types for each query cell

# 3. Visulizing MCA
# open shiny for visualize result for scMCA
scMCA_vis(mca_result)
