rm(list=ls(all=TRUE))
options(stringsAsFactors = FALSE)  

setwd("C:/Users/TPOW_NEW/Desktop/J20190806_ChiaChi/w20191120Datasets/Download/GSE128423")
da0<- read.table("GSE128423_tSNE_cluster_7gene_20200113-2.txt",sep="\t",header=T)
dim(da0)
da0[1:2,]

plot(da0[,"X"],da0[,"Y"], pch=20, xlab="tSNE 1", ylab="tSNE 2")




setwd("C:/Users/TPOW_NEW/Desktop/J20190806_ChiaChi/w20191120Datasets/Download/GSE128423/plot20200113-2")

tar<- 6
zz<- da0[,tar]
q<- colnames(da0)[tar]
summary(zz)
col_red <- colorRampPalette(c("white", "red"))(11)

Rplot<- paste("GSE128423_tSNE_",q,".jpg",sep="")
jpeg(Rplot, width =10, height =8, units = 'in', res = 300)

zz<- da0[,tar]
q<- colnames(da0)[tar]
zz0<- ceiling(max(zz))/9
zz3<- ceiling(zz/zz0)+1

par(fig=c(0,1,0,1),cex=1)
plot(da0[,"X"],da0[,"Y"], pch=20, xlab="tSNE 1", ylab="tSNE 2", col=col_red[zz3],
     main=paste("Gene=",q,sep=""))

par(fig=c(0,0.3,0.75,1), new=TRUE, cex=0.8)
plot(0:11,rep(1,12), xlab="", ylab="", axes=F, type="n")
for(i in 1:11){
polygon(c(i-1,i,i,i-1), c(0,0,1,1), col =col_red[i], border=NA)
}
axis(1,at=c(0.5,11),labels=c(0,ceiling(max(zz))),lwd=0.2, padj=-1)

dev.off()



tar<- 7
zz<- da0[,tar]
q<- colnames(da0)[tar]
summary(zz)
col_red <- colorRampPalette(c("white", "red"))(3)

Rplot<- paste("GSE128423_tSNE_",q,".jpg",sep="")
jpeg(Rplot, width =10, height =8, units = 'in', res = 300)


zz0<- ceiling(max(zz))/2
zz3<- ceiling(zz/zz0)+1

par(fig=c(0,1,0,1),cex=1)
plot(da0[,"X"],da0[,"Y"], pch=20, xlab="tSNE 1", ylab="tSNE 2", col=col_red[zz3],
     main=paste("Gene=",q,sep=""))

par(fig=c(0,0.3,0.75,1), new=TRUE, cex=0.8)
plot(0:3,rep(1,4), xlab="", ylab="", axes=F, type="n")
for(i in 1:3){
polygon(c(i-1,i,i,i-1), c(0,0,1,1), col =col_red[i], border=NA)
}
axis(1,at=c(0.5,3),labels=c(0,ceiling(max(zz))),lwd=0.2, padj=-1)

dev.off()




tar<- 10
zz<- da0[,tar]
q<- colnames(da0)[tar]
summary(zz)
col_red <- colorRampPalette(c("white", "red"))(3)

Rplot<- paste("GSE128423_tSNE_",q,".jpg",sep="")
jpeg(Rplot, width =10, height =8, units = 'in', res = 300)


zz0<- ceiling(max(zz))/2
zz3<- ceiling(zz/zz0)+1

par(fig=c(0,1,0,1),cex=1)
plot(da0[,"X"],da0[,"Y"], pch=20, xlab="tSNE 1", ylab="tSNE 2", col=col_red[zz3],
     main=paste("Gene=",q,sep=""))

par(fig=c(0,0.3,0.75,1), new=TRUE, cex=0.8)
plot(0:3,rep(1,4), xlab="", ylab="", axes=F, type="n")
for(i in 1:3){
polygon(c(i-1,i,i,i-1), c(0,0,1,1), col =col_red[i], border=NA)
}
axis(1,at=c(0.5,3),labels=c(0,ceiling(max(zz))),lwd=0.2, padj=-1)

dev.off()
#
