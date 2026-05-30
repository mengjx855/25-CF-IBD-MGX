#!/data/software/miniconda3/envs/r-base-4.4/bin/Rscript
# author: Jin-Xin Meng <jinxmeng@zju.edu.cn>

args <- commandArgs(trailingOnly = T)

if (length(args) < 2) {
  message("Usage: Rscript run-randomforest.R [CF-tpm] [out_prefix]")
  message("auto read metadata from /data/mengjx/project/10.20250623_IBD_BAC_CF_Landscape/profile/")
  # status = 1 表示非正常退出（给 Linux 系统看的信号）
  # save = "no" 表示不保存工作空间
  q(save = "no", status = 1)
}

pacman::p_load(tidyverse)
source('/data/mengjx/R_proj/R_func/model_rf_Kfold.R')

file_base = '/data/mengjx/project/10.20250623_IBD_BAC_CF_Landscape/profile/'

proj_name <- c(
  'BushmanFD_2020','FranzosaEA_2018','HallAB_2017','HeQ_2017','KumbhariA_2024',
  'LloydPriceJ_2019','SchirmerM_2018','SchirmerM_2024','WengY_2019','YanQ_2023c'
)

#### read tpm ####
tpm <- data.table::fread(args[1]) %>% 
  column_to_rownames('name')

#### rf ####
result <- map_dfr(
  proj_name, ~ {
    group <- read.delim(
      paste0(file_base, .x, '.sample_group.gz')
    ) %>% 
      dplyr::select(sample, group) %>% 
      mutate(group = ifelse(group == 'HC', 'HC', 'IBD'))
    
    profile <- tpm[,group$sample]
    
    pred <- rf_Kfold(profile, group, k = 10, seed = 2025, quiet = T)
    pred <- left_join(pred, group, by = 'sample')
    roc <- pROC::roc(pred$group, pred$HC, quiet = T)
    
    data.frame(
      proj = .x, 
      auc = as.numeric(roc$auc)
    )
  }
)

write.table(
  result, paste0(args[2], '.rf.tsv'), 
  sep = '\t', quote = F, row.names = F
)
