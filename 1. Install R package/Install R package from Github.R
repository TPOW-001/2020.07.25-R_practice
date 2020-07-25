
#从Github上获取R包/安装
#https://blog.csdn.net/tandelin/article/details/87601729

install.packages("glue")
install.packages("usethis")
install.packages("backports")
library(backports)
install.packages("devtools")

library(glue)
library(usethis)
library(devtools)
install_github('satijalab/seurat-data')
