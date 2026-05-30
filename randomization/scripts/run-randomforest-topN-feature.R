#!/data/software/miniconda3/envs/r-base-4.4/bin/Rscript
# author: Jin-Xin Meng <jinxmeng@zju.edu.cn>

args <- commandArgs(trailingOnly = T)

if (length(args) < 2) {
  message("Usage: run-randomforest-topN-feature.R [top_n] [out_file]")
  # status = 1 表示非正常退出（给 Linux 系统看的信号）
  # save = "no" 表示不保存工作空间
  q(save = "no", status = 1)
}

source('/data/mengjx/project/10.20250623_IBD_BAC_CF_Landscape/github/randomization/scripts/model_rf_Kfold.R')

proj_name <- c(
  'BushmanFD_2020','FranzosaEA_2018','HallAB_2017','HeQ_2017','KumbhariA_2024',
  'LloydPriceJ_2019','SchirmerM_2018','SchirmerM_2024','WengY_2019','YanQ_2023c'
)

select_features <- readr::read_rds(
  '6.feature_importance.names.rds'
)[1:as.numeric(args[1])]

file_base <- '/data/mengjx/project/10.20250623_IBD_BAC_CF_Landscape/github/data/'

tpm <- data.table::fread(
  '/data/mengjx/project/10.20250623_IBD_BAC_CF_Landscape/github/randomization/merge_tpm.fam.tpm'
  ) |> 
  tibble::column_to_rownames('name')

tpm <- tpm[grepl('CFF', rownames(tpm)),]

results <- purrr::map_dfr(
  proj_name, ~ {
    
    group <- read.delim(
      paste0(
        file_base,
        .x, 
        '.sample_group.bz2')
    ) |> 
      dplyr::select(sample, group) |> 
      dplyr::mutate(group = ifelse(group == 'HC', 'HC', 'IBD'))
    
    purrr:::map_dfr(
      2026:2030, \(x) {
        
        pred <- rf_Kfold(
          tpm[select_features, group$sample, ] |> na.omit(), 
          group, k = 10, seed = x, quiet = T)
        
        pred <- dplyr::left_join(pred, group, by = 'sample')
        
        roc <- pROC::roc(pred$group, pred$HC, quiet = T)
        
        tibble::tibble(
          proj = .x,
          seed = x,
          auc  = round(roc$auc, 6) 
        )
      }
    )
  }
) |> 
  tibble::add_column(
    top_n = paste0('top_', args[1]), 
    .before = 1
  )

readr::write_tsv(
  results, args[2], 
  col_names = F, quote = NULL
)
