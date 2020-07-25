
# 1. Loading DATA
rm(list=ls(all=TRUE))
setwd("F:/J20190806_ChiaChi/myDownload/GSE132151")
da0<- read.table("GSE132151_bone_marrow_stroma.counts.tsv.gz",sep="\t",header=T,fill=T, quote = "" )

dim(da0)
table(data0[,"cluster"])
xx<- data0[,"cluster"]=="MSC"

sum(xx)
data1<- data0[xx,]
dim(data1)
yy<- colnames(data1) %in% c("Acaa2","Acaca")
sum(yy)
(1:25307)[yy]
data2<- data1[,c(1:18, (1:25307)[yy])]
write.table(data2,"output.txt", sep="\t",row.name=F, col.name=T, quote=FALSE)



#########################################################################################
#########################################################################################
rm(list=ls(all=TRUE))

setwd("F:/J20190806_ChiaChi/myDownload/GSE132151")
data0<- read.table("TestData_20x25.txt",sep="\t",header=T,fill=T, quote = "" )
dim(data0)


table(data0[,"cluster"])
xx<- data0[,"cluster"]=="MSC"
sum(xx)
data1<- data0[xx,]
dim(data1)


yy<- colnames(data1) %in% c("X0610009O20Rik","X0610010F05Rik")
sum(yy)
(1:25)[yy]
data2<- data1[,c(1:18, (1:25)[yy])]
write.table(data2,"output.txt", sep="\t",row.name=F, col.name=T, quote=FALSE)