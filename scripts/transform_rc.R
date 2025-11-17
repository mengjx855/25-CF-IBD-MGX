#### Jinxin Meng, 20220529, 20250918, v0.1.2 ####

# 20230920: 修改了一些变量名称，增加rc2fpkm的函数
# 20250317: 修改函数格式
# 20250918: 修改函数格式，修改rc2fpkm的函数中1e12为1e9

library(dplyr)
library(tibble)
library(tidyr)

# TPM的计算方法：把比对到的某个基因的fragment数目，乘以1e6，除以基因的长度，再求相对丰度；
# FPKM的计算方法：把比对到的某个基因的fragment数目，乘以1e9，除以基因的长度，其比值再除以总map的reads数；
# RPM的计算方法：把比对到的某个基因的fragment数目，乘以1e6，除以此样本fragment总数；

#### rc2tpm ####
rc2tpm <- function(profile, length) {
  if (!all(c('name', 'length') %in% colnames(length) )) 
    stop('length must have columns: name, length (bp)')
  length <- tibble::column_to_rownames(length, 'name')[rownames(profile), 'length']
  tpm <- apply(profile, 2, \(x) x / length) %>% 
    apply(2, \(x) x / sum(x) * 1e6) %>% 
    data.frame(check.names = F)
  return(tpm)
}

#### rc2fpkm ####
rc2fpkm <- function(profile, length) {
  if (!all(c('name', 'length') %in% colnames(length) )) 
    stop('length must have columns: name, length (bp)')
  length <- tibble::column_to_rownames(length, 'name')[rownames(profile), 'length']
  total_reads <- colSums(profile)
  fpkm <- apply(profile, 2, \(x) x * 1e9 / length) %>% 
    apply(1, \(x) x / total_reads) %>% 
    t() %>% 
    data.frame(check.names = F)
  return(fpkm)
}

#### rc2rpkm ####
rc2rpkm <- function(profile, length) {
  if (!all(c('name', 'length') %in% colnames(length) )) 
    stop('length must have columns: name, length (bp)')
  length <- tibble::column_to_rownames(length, 'name')[rownames(profile), 'length']
  total_reads <- colSums(profile)
  rpkm <- apply(profile, 2, \(x) x * 1e9 / length) %>% 
    apply(1, \(x) x / total_reads) %>% 
    t() %>% 
    data.frame(check.names = F)
  return(rpkm)
}

#### rc2rpm ####
rc2rpm <- function(profile) {
  rpm <- apply(profile, 2, \(x) x / sum(x) * 1e6)
  return(rpm)
}
