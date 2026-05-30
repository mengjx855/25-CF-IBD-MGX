#!/data/software/miniconda3/envs/r-base-4.4/bin/Rscript
# author: Jin-Xin Meng <jinxmeng@zju.edu.cn>

args <- commandArgs(trailingOnly = T)

if (length(args) < 2) {
  message("Usage: Rscript mjx-adonis-by-fam.R [fam.tpm] [out_prefix]")
  message("file_base /data/mengjx/project/10.20250623_IBD_BAC_CF_Landscape/profile/")
  # status = 1 表示非正常退出（给 Linux 系统看的信号）
  # save = "no" 表示不保存工作空间
  q(save = "no", status = 1)
}

pacman::p_load(tidyverse)

get_pc_by_cumsum_var <- function(x, cumsum = .75) {
  .cumsum <- 0
  for (i in 1:length(x)) {
    .cumsum = .cumsum + x[i]
    if (.cumsum > cumsum) 
      return(i)
  }
}

calcu_adjusted_r2 <- function(adonis_object) {
  n_observations <- adonis_object$Df[3]+1
  d_freedom <- adonis_object$Df[1]
  r2 <- adonis_object$R2[1]
  vegan::RsquareAdj(r2, n_observations, d_freedom)
}

file_base <- '/data/mengjx/project/10.20250623_IBD_BAC_CF_Landscape/profile/'

proj_name <- c(
  'BushmanFD_2020','FranzosaEA_2018','HallAB_2017','HeQ_2017','KumbhariA_2024',
  'LloydPriceJ_2019','SchirmerM_2018','SchirmerM_2024','WengY_2019','YanQ_2023c'
)

tpm <- data.table::fread(args[1]) |>
  tibble::column_to_rownames('name') |>
  data.frame(check.names = FALSE)

result <- map_dfr(
  proj_name, ~ {
    group <- read.delim(
      paste0(file_base, .x, '.sample_group.gz')
    ) |> 
      select(sample, group)
  
    tpm_x <- tpm[, group$sample] |> 
      na.omit() |> 
      data.frame()
    
    data_x <- data.frame(na.omit(tpm[, group$sample])) %>%
      t() %>% 
      vegan::decostand(method = 'hellinger') %>% 
      data.frame()
  
    data_y <- data.table::fread(
      paste0(file_base, .x, '.kk2.s.gz')
    ) %>%
      column_to_rownames('name') %>%
      dplyr::select(any_of(rownames(data_x)))
  
    dist_y <- vegan::decostand(t(data_y), method = 'hellinger') %>%
      vegan::vegdist(method = 'bray')
  
    pca_result <- vegan::pca(data_x)
    pca_summary <- summary(pca_result)
    pc_axis <- get_pc_by_cumsum_var(
      pca_summary$cont$importance[2,], cumsum = .95
    )
    variables <- vegan::scores(
      pca_result, display = 'sites', choices = 1:pc_axis
    )
  
    adonis <- vegan::adonis2(
      as.formula(
        paste0(
          'dist_y ~ ', 
          paste0(colnames(variables), collapse = ' + ')
        )
      ), 
      data.frame(variables), permutations = 1, parallel = 1)
  
    tibble(
      proj = .x, 
      r2 = adonis$R2[1],
      r2adj = calcu_adjusted_r2(adonis),
      pc_axis = pc_axis,
      pc_var = 0.95,
      details = 'explain_by_CF'
    )
  } 
)

write.table(
  result, paste0(args[2], '.explain_by_CF.tsv'),
  sep = '\t', quote = F, row.names = F
)
