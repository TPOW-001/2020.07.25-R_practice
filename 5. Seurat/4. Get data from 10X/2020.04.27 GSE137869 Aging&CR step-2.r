#Step1: 先把每一個資料轉成gene-cell matrix
#要在資料夾GSE137869建立一個名叫box的資夾
rm(list=ls(all=TRUE))
options(stringsAsFactors = FALSE) 

library(Seurat)

#資料夾裡的東西必是如此命名,.gz可不用解壓縮
#barcodes.tsv.gz, genes.tsv.gz, matrix.mtx.gz
#要把genes.tsv.gz 改成features.tsv.gz

test.data <- Read10X(data.dir="C:/Users/TPOW_NEW/Desktop/2020.04.05 GSE137869/GSE137869/GSM4331822_BM-M-Y") #資料存放位置 

da0 <- CreateSeuratObject(counts=test.data, project="test.A1", min.cells=0, min.features=0)
# min.cells,min.features 建議是0, 不然會篩data
# ex: min.cells=3,min.features=200 等於留三個以上cell表現的gene, 及有200以上gene表現的 cell

#出現以下 warning可略過
#Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

da0    #看一下有多少gene和cell

target<- da0@assays$RNA@counts  #取出count值
dim(target)
target[1:3, 1:3]
#row=gene,  col=cell
#cell 已經有加識別來自哪個sample了   ###是一開始的原始資料就加了嗎?


t0<- Sys.time()
setwd("C:/Users/TPOW_NEW/Desktop/2020.04.05 GSE137869/GSE137869/Box")   #存放位置
write.table(target,"BM-M-Y.txt",sep="\t", row.names=T,col.names=T,quote=F)  #存哪筆data與存放檔名
tn<- Sys.time()
print( difftime(tn,t0, tz=" ",unit='mins') )  #存很久是正常, 1~10分鐘不等

#如此已經將data存成 cell-gene-matrix了
#重複將想看的data轉好

# A1=GSM4331814_Aorta-F-O,  32884 features across 3849 samples within 1 assay 
# A2=GSM4331815_Aorta-F-CR, 32884 features across 2582 samples within 1 assay 


########################################################################################################
########################################################################################################
#Step2: 合併gene-cell matrix

rm(list=ls(all=TRUE))
options(stringsAsFactors = FALSE) 

setwd("C:/Users/TPOW_NEW/Desktop/2020.04.05 GSE137869/GSE137869/Box")

BM1<- read.table("BM-M-Y.txt",sep="\t",header=T) #讀很久也是正常, 要注意記憶體夠不夠, 不夠會出現ERROR
BM2<- read.table("BM-M-O.txt",sep="\t",header=T)
BM3<- read.table("BM-M-CR.txt",sep="\t",header=T)
dim(BM1)
dim(BM2)
dim(BM3)

#先確認gene的排列順序是否一樣
rownames(BM1)==rownames(BM2)
sum(rownames(BM1)==rownames(BM2))
dim(BM1)
dim(BM2)
dim(BM3)
BM1[1:5,1:20]
BM2[1:5,1:20]
BM3[1:5,1:20]
BM1[32881:32884,1:2]
BM2[32881:32884,1:2]
BM3[32881:32884,1:2]
xx<- rownames(BM1)=="Pecam1"
(1:32883)[xx]

BM1_1<- BM1[c(1:32883, 32883), ]
BM3_1<- BM3[c(1:32883, 32883), ]
rownames(BM1_1)[32884]<- rownames(BM3_1)[32884]<- "Pecam1"
BM1_1[32884,]<- BM3_1[32884,]<- 0
data<- cbind(BM1_1, BM2, BM3_1)
dim(data)
data[32884,]

#此指令的意義是確認A1的第一個gene與A2的gene是否一樣, 一樣會給TRUE=1,  不一樣會給FALSE=0, 
#全部gene數是32884, 所以會希望 此比較的總和是32884=兩個檔案gene的排列順序是一樣
#
#有多個檔案合併時, 最好都跟某一個比一下, 確定順序都一樣
#不然在順序不一樣的時候合併會出事
#順序不一樣的時候, 當然就是想辦法把它順序弄成一樣, 遇到再求救, 先略過


#####暫時未用#####
#橫向合併
data<- cbind(A1, A2)
dim(data)
#####暫時未用#####

t0<- Sys.time()
setwd("C:/Users/TPOW_NEW/Desktop/2020.04.05 GSE137869/GSE137869/box")   #存放位置
write.table(data,"BM_M_comb1.txt",sep="\t", row.names=T,col.names=T,quote=F)  #存哪筆data與存放檔名
tn<- Sys.time()
print( difftime(tn,t0, tz=" ",unit='mins') )  #存很久是正常, 1~10分鐘不等
