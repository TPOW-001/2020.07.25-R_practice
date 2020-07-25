
##
# https://www.datanovia.com/en/blog/how-to-normalize-and-standardize-data-in-r-for-great-heatmap-visualization/

# 1.Loading data
setwd("C:/Users/TPOW31714/Desktop/2020.07.16 Normalized 06.03 RMA")
da0<- read.table("5. Yen list-1 RMA.txt",sep="\t",row.name=T, col.name=T, quote=FALSE)
da0<- read.table("5. Yen list-1 RMA.txt",sep="\t",header=T)
dim(da0)
da0[1:5,1:7]
names(da0)

# Heat map
library(heatmaply)
heatmaply(da0, xlab = "Cell types", ylab = "Gene symbols", main = "Raw data")

da1<- da0[,-1]
rownames(da1)<- da0[,1]
da1<- scale(da1)
heatmaply(da1, xlab = "Cell types", ylab = "Gene symbols", main = "Data Scaling")

write.table(da1,"1. RMA standardrization.txt",sep="\t", row.names=T,col.names=T,quote=F)


da2<- da0[,-1]
da2<- normalized(da2)
heatmaply(da2, xlab = "Features", ylab = "Cars", main = "Data Normalization")

da3<- da0[,-1]
da3<- percentize(da2)
heatmaply(da3, xlab = "Features", ylab = "Cars", main = "Percentile Transformation")






# Example
# https://cran.r-project.org/web/packages/heatmaply/vignettes/heatmaply.html
library("heatmaply")
heatmaply(mtcars)
mtcars

# 1.Correlation heatmaps
heatmaply_cor(
  cor(mtcars),
  xlab = "Features",
  ylab = "Features",
  k_col = 2,
  k_row = 2
)

# 1.1 advanced correlation heatmap
r <- cor(mtcars)
## We use this function to calculate a matrix of p-values from correlation tests
## https://stackoverflow.com/a/13112337/4747043
cor.test.p <- function(x){
  FUN <- function(x, y) cor.test(x, y)[["p.value"]]
  z <- outer(
    colnames(x), 
    colnames(x), 
    Vectorize(function(i,j) FUN(x[,i], x[,j]))
  )
  dimnames(z) <- list(colnames(x), colnames(x))
  z
}
p <- cor.test.p(mtcars)

heatmaply_cor(
  r,
  node_type = "scatter",
  point_size_mat = -log10(p), 
  point_size_name = "-log10(p-value)",
  label_names = c("x", "y", "Correlation")
)

# 2. Data transformation (scaling, normalize, and percentize)
# 2.1
heatmaply(
  mtcars, 
  xlab = "Features",
  ylab = "Cars", 
  scale = "column",
  main = "Data transformation using 'scale'"
)
# 2.2
heatmaply(
  normalize(mtcars),
  xlab = "Features",
  ylab = "Cars", 
  main = "Data transformation using 'normalize'"
)
# 2.3
heatmaply(
  percentize(mtcars),
  xlab = "Features",
  ylab = "Cars", 
  main = "Data transformation using 'percentize'"
)




