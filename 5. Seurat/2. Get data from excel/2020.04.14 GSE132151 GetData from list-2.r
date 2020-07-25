rm(list=ls(all=TRUE))
#row是讀水平, column是垂直

setwd("C:/Users/TPOW_NEW/Desktop/2020.03.18 GSE132151")
t0<- Sys.time()
data0<- read.table("bone_marrow_stroma.counts.tsv",sep="\t",header=T,fill=T, quote = "" )
tn<- Sys.time()
difftime(tn,t0, tz=" ",unit='mins')
#Time difference of 5.648314 mins
#讀取檔案以及確認檔案讀取時間

dim(data0)
# 15028 cells *25307 
data0[1:5,1:20]
colnames(data0)[1:20]
# [1] "barcode_id"             "barcode"                "library"                "sample"                
# [5] "n_genes"                "n_counts"               "mito_frac"              "pass_filter"           
# [9] "cluster"                "SPRING_x"               "SPRING_y"               "PBA_Potential"         
#[13] "PBA_Pr_Adipo"           "PBA_Pr_Osteo"           "PBA_Pr_Chondro"         "PBA_Osteo_pseudotime"  
#[17] "PBA_Chondro_pseudotime" "PBA_Adipo_pseudotime"   "X0610006L08Rik"         "X0610007P14Rik"     

table(data0[,"pass_filter"])
table(data0[,"sample"])
table(data0[,"pass_filter"],data0[,"sample"])
table(data0[,"cluster"])
# table(data0[,"cluster"])是回報column"cluster"所有row的組成比例

plot(data0[,"SPRING_x"],data0[,"SPRING_y"], pch=16)
plot(data0[,"SPRING_x"],data0[,"SPRING_y"], col= as.factor(data0[,"cluster"]), pch=16)
# 利用spring x&y座標繪製圖形.
# col= as.factor(data0[,"cluster"] 利用cluster來區分顏色
# pch=solid dots increase the readability of this data plot

xx<- data0[,"pass_filter"]=="True"
sum(xx)
# 選擇"pass_filter"標示為True的組成
da2<- data0[xx,]
dim(da2)
# 選擇"pass_filter"標示為True的組成, 並建立da2只呈現"pass_filter"標示為True的組成
da3<- t(da2)
dim(da3)
# 建立da3將da2 row以及column轉置
setwd("C:/Users/TPOW_NEW/Desktop/2020.03.18 GSE132151/Data edited")
write.table(da2,"GSE132151_2847cells_20200318.txt", sep="\t",row.name=F, col.name=T, quote=FALSE)
write.table(da3,"GSE132151_2847cells_20200318t.txt", sep="\t",row.name=T, col.name=T, quote=FALSE)
#儲存檔案

plot(da2[,"SPRING_x"],da2[,"SPRING_y"], col= as.factor(da2[,"cluster"]), pch=16)
#利用spring x&y座標, 繪製da2 plot,顏色利用cluster區分
# 如何更改顏色???????


# !!!!!利用資料塞選出想要的gene, 並利用spring x&y座標繪製圖形

xx<- colnames(da2)=="Sox9"
sum(xx)
(1:25307)[xx]
# Sox9在21406的位置

da0<- da2
tar<- 21406
zz<- da0[,tar]
# zz<- da0[,tar]是column"21406位置"所有row的組成比例
q<- colnames(da0)[tar]
summary(zz)
col_red <- colorRampPalette(c("white", "red"))(11)
# 打算將顏色分成11等份, 可以利用summary(zz)以及table(zz)決定想要區分的顏色
col_red[1]<- "gray90"
# http://iccm.cc/colors-and-palettes-in-r-language/  講解colorRampPalette
# white->red分成11個顏色
# https://blog.csdn.net/Bone_ACE/article/details/47362619 R語言色碼表

zz<- da0[,tar]
table(zz)
#    0    1    2    3    4 
# 2621  181   36    2    7
q<- colnames(da0)[tar]
summary(zz)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0000  0.0000  0.0000  0.1008  0.0000  4.0000 

# !!!數值是常態分佈的話用以下公式
zz0<- ceiling(max(zz))/9
# summary(zz)的最大值-1 來替換數字9的部分
zz3<- ceiling(zz/zz0)+1
table(zz3)
# ceiling的公式代表為天花板
# 上述兩個公式舉例:
# 數字0,1,2,3,4---100, 打算將0畫一個顏色, 1-10一個顏色----99-100一個顏色, 一共11個顏色


# !!!數值式"非"常態分佈的話用以下公式
zz3<- zz
zz3[zz>0]<- 1
zz3<- zz3+1
table(zz3)
#    1    2 
# 2621  226
## 第二組226= 181+36+2+7
col_red<- c("gray90","red")

#假設將顏色分成三群,第一群數值0, 第二群數值1, 第三群數值2-4, 可以直接寫zz3[zz>=2]<- 2,故變成
zz3<- zz
zz3[zz>=2]<- 2
zz3<- zz3+1
table(zz3)
col_red <- colorRampPalette(c("white", "red"))(3)
col_red[1]<- "gray90"

par(fig=c(0,1,0,1),cex=1)
plot(da0[,"SPRING_x"],da0[,"SPRING_y"], pch=20, type='n',
     xlab="tSNE 1", ylab="tSNE 2",main=paste("Gene=",q,sep=""))
# type='n'是為了畫出格子, 但不要標出顏色
# main=paste("Gene=",q,sep="")是什麼? 猜測是title, 但q表示?

table(zz3)

for(i in 1:3){
xx<- zz3==i
points(da0[xx,"SPRING_x"],da0[xx,"SPRING_y"], pch=20, col=col_red[i])
}
# 此公式是為了在空白格子內繪圖, 顏色由淺到深, 以這個例子來說, i=1 灰色, i=2 粉紅, i=3 紅色
# i in 1:3, 此處的3是來自於col_red <- colorRampPalette(c("white", "red"))(3)

par(fig=c(0,0.3,0.75,1), new=TRUE, cex=0.8)
plot(0:3,rep(1,4), xlab="", ylab="", axes=F, type="n")
#畫出label的格子
for(i in 1:3){
polygon(c(i-1,i,i,i-1), c(0,0,1,1), col =col_red[i], border=NA)
}
axis(1,at=c(0,1,2,3),labels=c(0,1,2,4),lwd=0.2, padj=-1)
# 因為顏色分成3群, 0,1,2-4, 故label座標為0,1,2,4





##以下公式會有顏色被覆蓋的疑慮:
par(fig=c(0,1,0,1),cex=1)
plot(da0[,"SPRING_x"],da0[,"SPRING_y"],pch=20, xlab="tSNE 1", ylab="tSNE 2", col=col_red[zz3],main=paste("Gene=",q,sep=""))






# 2020.03.22目前快速畫圖公式:
xx<- colnames(da2)=="Sox9"
sum(xx)
(1:25307)[xx]
da0<- da2

tar<- 21406

zz<- da0[,tar]
table(zz)
summary(zz)

#假設將顏色分成三群,第一群數值0, 第二群數值1, 第三群數值2-4, 可以直接寫zz3[zz>=2]<- 2,故變成
zz3<- zz
zz3[zz>=2]<- 2
zz3<- zz3+1
table(zz3)
col_red <- colorRampPalette(c("white", "red"))(3)
col_red[1]<- "gray90"

par(fig=c(0,1,0,1),cex=1)
plot(da0[,"SPRING_x"],da0[,"SPRING_y"], pch=20, type='n',
     xlab="tSNE 1", ylab="tSNE 2",main=paste("Gene=",q,sep=""))
for(i in 1:3){
  xx<- zz3==i
  points(da0[xx,"SPRING_x"],da0[xx,"SPRING_y"], pch=20, col=col_red[i])
}
# 此公式是為了在空白格子內繪圖, 顏色由淺到深, 以這個例子來說, i=1 灰色, i=2 粉紅, i=3 紅色
# i in 1:3, 此處的3是來自於col_red <- colorRampPalette(c("white", "red"))(3)

par(fig=c(0,0.3,0.75,1), new=TRUE, cex=0.8)
plot(0:3,rep(1,4), xlab="", ylab="", axes=F, type="n")
for(i in 1:3){
polygon(c(i-1,i,i,i-1), c(0,0,1,1), col =col_red[i], border=NA)
}
axis(1,at=c(0,1,2,3),labels=c(0,1,2,4),lwd=0.2, padj=-1)
# 因為顏色分成3群, 0,1,2-4, 故label座標為0,1,2,4

zz0<- ceiling(max(zz))/9
zz3<- ceiling(zz/zz0)+1
table(zz3)