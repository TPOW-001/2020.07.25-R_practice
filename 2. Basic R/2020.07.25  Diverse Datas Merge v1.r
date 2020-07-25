# 1. Load data
rm(list=ls(all=TRUE))
options(stringsAsFactors = FALSE) 

library(magrittr)
library(dplyr)
setwd("C:/Users/TPOW31714/Desktop")
da1<- read.table("1. BMMSC OvsY p under 0.05.txt", header=T, sep="\t")
da2<- read.table("2. LUMSC OvsY p under 0.05.txt", header=T, sep="\t")
da3<- read.table("3. PAMSC OvsY p under 0.05.txt", header=T, sep="\t")


# 2. Merge函數
# https://sites.google.com/site/rlearningsite/manage/merge
Merge_1<-merge(da1, da2, by="NAME")
Merge_2<-merge(Merge_1, da3, by="NAME")

# 3. Save data
write.table(Merge_2,"2020.07.25 GSEAresult_Merge_MSCs OvsY.txt", 
            sep="\t",row.name=F, col.name=T, quote=FALSE)


# 相關課程
# https://medium.com/datainpoint/r-essentials-manipulating-dataframes-advancedly-dcf0b0f76ecd
# 1. 交集
# 左邊的資料框
name <- c("蒙其·D·魯夫", "羅羅亞·索隆", "娜美", "多尼多尼·喬巴")
age <- c(19, 21, 20, 17)
left_df <- data.frame(name, age)
# 右邊的資料框
name <- c("蒙其·D·魯夫", "多尼多尼·喬巴", "妮可·羅賓")
devil_fruit <- c("橡膠果實", "人人果實", "花花果實")
right_df <- data.frame(name, devil_fruit)
# 1.1 merge
merge(left_df, right_df)
# 1.2 dplyr
library(dplyr)
inner_join(left_df, right_df)

# 2 保留所有左邊資料框的觀測值
# 2.1 Merge
merge(left_df, right_df, all.x = TRUE)
# 2.2 dplyr
left_join(left_df, right_df)

# 3 保留所有右邊資料框的觀測值
# 3.1 Merge
merge(left_df, right_df, all.y = TRUE)
# 3.2 dplyr
right_join(left_df, right_df)

# 4 保留兩邊資料框的觀測值
# 4.1 Merge
merge(left_df, right_df, all.x = TRUE, all.y = TRUE)
# 4.2 dplyr
full_join(left_df, right_df)



# TEST
# https://blog.csdn.net/woodcorpse/article/details/80494605

# 1. Load data
rm(list=ls(all=TRUE))
options(stringsAsFactors = FALSE) 
setwd("C:/Users/TPOW31714/Desktop")
da1<- read.table("1. BMMSC OvsY p under 0.05.txt", header=T, sep="\t")
da2<- read.table("2. LUMSC OvsY p under 0.05.txt", header=T, sep="\t")
da3<- read.table("3. PAMSC OvsY p under 0.05.txt", header=T, sep="\t")

# 2. 交集intersect
# intersect()

# 2.1
Intersect_1<-data.frame(intersect(da1[,1], da2[,1]))

# 2.2
library(magrittr)
library(dplyr)
Intersect_2<- intersect(da1[,1], da2[,1]) %>% 
  data.frame()
              
# 3. 聯集union
# union()
Union_1<-data.frame(union(da1[,1], da2[,1]))

# 4.求向量x与向量y中不同的元素(只取x中不同的元素)
# setdiff(x, y)
# 找da1中不同于da2的元素
Setdiff_1<-data.frame(setdiff(da1[,1], da2[,1]))
# 找da2中不同于da1的元素
Setdiff_2<-data.frame(setdiff(da2[,1], da1[,1]))

# 5. 判断相同setequal
# setequal(x, y)
x = 1:4
y = 2:6
# 判断x与y是否相同，结果为假
setequal(x, y)
# [1] FALSE
# 找y与x是否相同，结果为假
setequal(y, x)
# [1] FALSE

# 只有完全相同的才返回TRUE
y = 1:4
setequal(x, y)
# [1] TRUE


