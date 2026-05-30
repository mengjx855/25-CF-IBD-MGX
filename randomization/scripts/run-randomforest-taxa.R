#!/data/software/miniconda3/envs/r-base-4.4/bin/Rscript
# author: Jin-Xin Meng <jinxmeng@zju.edu.cn>

args <- commandArgs(trailingOnly = T)

if (length(args) < 4) {
  message("Usage: Rscript run-randomforest-taxa.R [select_taxa] [kk2_profile] [group.tsv] [out_prefix]")
  # status = 1 表示非正常退出（给 Linux 系统看的信号）
  # save = "no" 表示不保存工作空间
  q(save = "no", status = 1)
}

source('model_rf_Kfold.R')

select_taxa <- stringr::str_split_1(as.character(args[1]), ',')

profile <- data.table::fread(args[2]) |> 
  tibble::column_to_rownames('name')

profile <- profile[select_taxa,]

group <- data.table::fread(args[3]) |> 
  dplyr::select(sample, group) |> 
  dplyr::mutate(
    group = ifelse(group == 'HC', 'HC', 'IBD')
  )

pred <- rf_Kfold(profile, group, k = 10, seed = 2025, quiet = T)
pred <- dplyr::left_join(pred, group, by = 'sample')
roc <- pROC::roc(pred$group, pred$HC, quiet = T)

result <- tibble::tibble(
  name = basename(args[4]),
  auc  = round(roc$auc, 6)
)

readr::write_tsv(
  result, paste0(args[4], '.taxa.rf.tsv'), 
  col_names = F, quote = NULL
)
