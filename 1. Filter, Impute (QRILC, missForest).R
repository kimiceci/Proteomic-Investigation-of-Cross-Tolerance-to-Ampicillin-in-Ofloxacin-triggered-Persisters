library(tidyverse)
library(ggplot2)
library(ggrepel)
library(imputeLCMD)
library(missForest)
library(dplyr)

# 设置工作目录
setwd("D:/Final Report/6500")

filePath <- file.path(paste0("proteinGroups_6500.txt"))
# read the protein groups data and change the " " in the header into "."
PG <- read.delim(filePath, header = TRUE, stringsAsFactors = FALSE, na.strings = "") %>%
  setNames(gsub(" ", ".", names(.)))

PG_filtered<-PG


PG_filtered$Gene.names <- sub(".*?GN=([^ ]+).*","\\1", PG_filtered$Fasta.headers)

PG_filtered <- PG_filtered %>% filter(is.na(Reverse))
PG_filtered <- PG_filtered %>% filter(is.na(Potential.contaminant))
PG_filtered <- PG_filtered %>% filter(is.na(Only.identified.by.site))

PG_filtered <- PG_filtered %>%
  filter(!((LFQ.intensity.exp1.1==0|LFQ.intensity.exp1.2==0|LFQ.intensity.exp1.3==0)&
             (LFQ.intensity.exp6.1==0|LFQ.intensity.exp6.2==0|LFQ.intensity.exp6.3==0)&
             (LFQ.intensity.exp7.1==0|LFQ.intensity.exp7.2==0|LFQ.intensity.exp7.3==0)))

LFQ_all <- PG_filtered %>%
  select(starts_with("LFQ.intensity.exp"),
         Gene.names,
         Majority.protein.IDs)

#这行代码用于检查LFQ_all数据框中Majority protein IDs列是否有重复的值
LFQ_all$Majority.protein.IDs %>% duplicated() %>% any()
#这行代码用于统计LFQ_all数据框中Gene.names列中每个基因名称的出现频率，并筛选出出现频率大于1的基因名称
LFQ_all %>% group_by(Gene.names) %>% summarize(frequency = n()) %>% arrange(desc(frequency)) %>% filter(frequency > 1)
# A tibble: 0 × 2
# ℹ 2 variables: Gene.names <chr>, frequency <int>

extract_ID <- function(c1,c2,c4){
  # 定义优先级顺序
  priority <- c('P', 'Q', 'O', 'G', 'V', 'A')
  # 检查 c2 的长度和第一个元素
  if (length(c2) == 1 || c2[1] == "P") {
    c3 <- c1[1]
    c5 <- c4[1]
    }else{
      # 计算匹配的索引
      match_indices <- order(match(c2, priority))
      # 根据匹配的索引选择 c1 和 c4 的第一个元素
      c3 <- c1[match_indices][1]
      c5 <- if (length(c4) == 1) c4[1] else c4[match_indices][1]
      }
  return(c(c3,c5))
}

#创建两个字符型向量name和ID，它们的长度等于数据框LFQ_all的行数
name <- character(nrow(LFQ_all))
ID <- character(nrow(LFQ_all))

for(i in 1:length(name)){
#从LFQ_all的列Majority.protein.IDs中提取第i个元素的值，然后使用分号（;）将其分割成多个子字符串，并将这些子字符串存储到向量c1中
  c1 <- unlist(str_split(LFQ_all$Majority.protein.IDs, ";")[i])
#从LFQ_all的列Majority.protein.IDs中提取第i个元素的值，然后使用分号（;）将其分割成多个子字符串，并从每个子字符串中提取第一个字符，最后将这些第一个字符存储到向量c2中
  c2 <- unlist(lapply(str_split(LFQ_all$Majority.protein.IDs, ";"),substr,start=1,stop=1)[i])
#从LFQ_all的列Gene.names中提取第i个元素的值，然后使用分号（;）将其分割成多个子字符串，并将这些子字符串存储到向量c4中
  c4 <- unlist(str_split(LFQ_all$Gene.names, ";")[i])
  temp=extract_ID(c1,c2,c4)
  ID[i]=temp[1]
  name[i]=temp[2]
}

LFQ <- LFQ_all %>%
  mutate(
    name = name,
    ID = ID,
  )

LFQ <- LFQ[c(1:9,12:13)]

#########impute MNAR###############
# MNAR: missing in >50% in a distinct study group 
# (n=3, that means the numbers of missing value > 1)
# MCAR/MAR: missing in <50% in a distinct study group
# (n=3, that means the numbers of missing value <= 1)
LFQ_filtered1 <- LFQ %>%
  mutate(miss_num_1 = rowSums(select(., 1:3) == 0)>1,
         miss_num_6 = rowSums(select(., 4:6) == 0)>1,
         miss_num_7 = rowSums(select(., 7:9) == 0)>1,
         )
# MNAR: TRUE, fill with QRILC; MCAR/MAR: FALSE, fill with RF

LFQ_0toNAforQRILC <- LFQ_filtered1

col_S1 <- grep("^LFQ.intensity.exp1.", colnames(LFQ_0toNAforQRILC))
col_S6 <- grep("^LFQ.intensity.exp6.", colnames(LFQ_0toNAforQRILC))
col_S7 <- grep("^LFQ.intensity.exp7.", colnames(LFQ_0toNAforQRILC))

col_indices <- data.frame(S1 = col_S1, S6 = col_S6, S7 = col_S7)
miss_indices <- grep("^miss_num_", colnames(LFQ_0toNAforQRILC))

for (i in 1:3){
  # 获取当前需要处理的行和列 
  rows_to_replace <- LFQ_0toNAforQRILC[, miss_indices[i]]  # 逻辑向量，指示哪些行需要处理 
  cols_to_replace <- col_indices[[i]]       # 需要处理的列索引 
  # 提取目标子集 
  target_subset <- LFQ_0toNAforQRILC[rows_to_replace, cols_to_replace]
  # 替换 0 为 NA 
  LFQ_0toNAforQRILC[rows_to_replace, cols_to_replace] <- replace(target_subset, target_subset == 0, NA)
#  LFQ_0toNAforQRILC[LFQ_0toNAforQRILC[,miss_indices[i]],col_indices[[i]]] <- replace(LFQ_0toNAforQRILC[LFQ_0toNAforQRILC[,miss_indices[i]],col_indices[[i]]], LFQ_0toNAforQRILC[LFQ_0toNAforQRILC[,miss_indices[i]],col_indices[[i]]] == 0, NA)
}

set.seed(0)
LFQ_tobeimputed_QRILC <- log2(LFQ_0toNAforQRILC[,1:9]+1)

# 随机性来自Gibbs随机抽样
LFQ_imputed_QRILC<-impute.QRILC(LFQ_tobeimputed_QRILC)
LFQ_imputed_QRILC<-LFQ_imputed_QRILC[[1]]

# 把剩余的0值换成NA然后转置
LFQ_tobeimputed_MF <- LFQ_imputed_QRILC
LFQ_tobeimputed_MF[LFQ_tobeimputed_MF==0] <- NA
LFQ_tobeimputed_MF <- as.data.frame(t(LFQ_tobeimputed_MF))

# 添加分组信息作为missForest的特征
LFQ_tobeimputed_MF$Groups <- c(rep("S1", 3), rep("S6", 3), rep("S7", 3))
LFQ_tobeimputed_MF$Groups <- as.factor(LFQ_tobeimputed_MF$Groups)

# 随机性来自bagging
LFQ_imputed_MF<-missForest(LFQ_tobeimputed_MF)
LFQ_imputed<-t(2^LFQ_imputed_MF$ximp[1:ncol(LFQ_imputed_MF$ximp)-1]-1)

# 为消除对数转换对原始数据的影响，只把缺失值（即原为0的位置）替换成QRILC和missForest的填补值
LFQ_filtered2<-LFQ_filtered1[1:9]
LFQ_filtered2[LFQ_filtered2==0]<-LFQ_imputed[LFQ_filtered2==0]

# Exact column names as specified 
lfq_cols_exp <- grep("^LFQ.intensity.exp", colnames(LFQ), value = TRUE)

bfimputed <- PG_filtered
bfimputed$name<-name
bfimputed$ID<-ID
#保存过滤完成的output为.csv文件
bfimputedPath <- file.path("Impute", paste0("proteinGroups_6500_filtered.csv"))
write.csv(bfimputed, bfimputedPath, row.names = FALSE)

output <- PG_filtered
output$name<-name
output$ID<-ID
output[lfq_cols_exp]<-LFQ_filtered2[lfq_cols_exp]
#保存过滤且插值完成的output为.csv文件
outputPath <- file.path(paste0("proteinGroups_6500_filtered_imputed.csv"))
write.csv(output, outputPath, row.names = FALSE)