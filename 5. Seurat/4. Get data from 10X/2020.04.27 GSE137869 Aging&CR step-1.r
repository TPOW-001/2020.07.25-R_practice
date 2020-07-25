#Step0: 檔案分資料夾及改檔名
# 要把GSE137869_RAW (下載的data)跟GSE137869 放在同一資料夾底下
# GSE137869要拿來放分好的資料夾


rm(list=ls(all=TRUE))
options(stringsAsFactors = FALSE) 

setwd("C:/Users/TPOW_NEW/Desktop/2020.04.05 GSE137869/GSE137869_RAW")
zz0<- list.files() #下載資料所在資料夾中有那些檔案
zz0
zz<- zz0[7:168] # 前面六個先跳過
ZZ

#  [1] "GSM4321717_GL2_1_fC_all_rawcounts.p.txt.gz"
#  [2] "GSM4321718_GL2_2_fC_all_rawcounts.p.txt.gz"
#  [3] "GSM4321719_GL2_3_fC_all_rawcounts.p.txt.gz"
#  [4] "GSM4321720_YBX_1_fC_all_rawcounts.p.txt.gz"
#  [5] "GSM4321721_YBX_2_fC_all_rawcounts.p.txt.gz"
#  [6] "GSM4321722_YBX_3_fC_all_rawcounts.p.txt.gz"
# 剩下的檔名列出來時的有個規律, 同一個sample會再一起, 三個一組, 且有個規律先barcodes.tsv.gz,genes.tsv.gz,再matrix.mtx.gz  共54個samples

i<- 1
setwd("C:/Users/TPOW_NEW/Desktop/2020.04.05 GSE137869/")
for(i in 1:54){
  aa<- unlist(strsplit(zz[i*3-2],"_")) #支解檔名
  target<- paste(aa[1], aa[2], sep="_") #取出要做為資料夾名稱的名字
  dir.create( paste("GSE137869/",target, sep="") ) #建立資料夾
  file.copy(from=paste("GSE137869_RAW/",zz[i*3-2], sep='' ),to=paste("GSE137869/",target, "/barcodes.tsv.gz", sep='' ))
  file.copy(from=paste("GSE137869_RAW/",zz[i*3-1], sep='' ),to=paste("GSE137869/",target, "/features.tsv.gz", sep='' ))
  file.copy(from=paste("GSE137869_RAW/",zz[i*3], sep='' ),to=paste("GSE137869/",target, "/matrix.mtx.gz", sep='' ))
  #移動複製檔案同時改名
}
#
