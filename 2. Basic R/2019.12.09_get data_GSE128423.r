rm(list=ls(all=TRUE))
options(stringsAsFactors = FALSE)  

setwd("C:/Users/TPOW_NEW/Desktop/J20190806_ChiaChi/w20191120Datasets/Download/GSE128423")
da0<- read.table("stroma.tsne.txt",sep="\t",header=T)
dim(da0)
da0[1:2,]
da<- da0[-1,]
da[1:2, ]
dim(da)

plot(da[,2],da[,3], pch=20, xlab="tSNE 1", ylab="tSNE 2")


db0<- read.table("stroma.tsne.meta.txt",sep="\t",header=T)
dim(db0)
db0[1,]
db<- db0[-1,]
db[1:3,]

sum(db[,1]==da[,1])

dc0<- read.table("stroma.TP4K.txt",sep="\t",header=T)
dim(dc0)#[1] 27998 20582
dc0[1:5,1:3]

xx<- dc0[,1] %in% c("Cebpb","Pparg","Klf5","Lep","Nppa")
sum(xx)
tar<- dc0[xx,]
tar[1:2,1:4]
tar2<- t(tar)
tar2[1:2,1:4]
colnames(tar2)<- tar2[1,]
tar3<- tar2[-1,]
tar3[1:2,]


sum(db[,1]==rownames(tar3) )

box<- cbind(da, db, tar3)
box[1:3, ]
box2<- box[,-4]
box2[1:3, ]

setwd("C:/Users/TPOW_NEW/Desktop/J20190806_ChiaChi/w20191120Datasets/Download/GSE128423")
write.table(box2,"GSE128423_tSNE_cluster_19gene_20200212-1.txt", sep="\t",row.name=F, col.name=T, quote=FALSE)

