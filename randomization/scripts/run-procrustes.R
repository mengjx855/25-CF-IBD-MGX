#!/data/software/miniconda3/envs/r-base-4.4/bin/Rscript
# author: Jin-Xin Meng <jinxmeng@zju.edu.cn>

args <- commandArgs(trailingOnly = T)

if (length(args) < 2) {
  message("Usage: Rscript run-procrustes.R [gene-fam-tpm] [out_prefix]")
  message("auto read kk2 profile from /data/mengjx/project/10.20250623_IBD_BAC_CF_Landscape/profile/")
  # status = 1 表示非正常退出（给 Linux 系统看的信号）
  # save = "no" 表示不保存工作空间
  q(save = "no", status = 1)
}

pacman::p_load(tidyverse)

kk_base = '/data/mengjx/project/10.20250623_IBD_BAC_CF_Landscape/profile/'

proj_name <- c(
  'BushmanFD_2020','FranzosaEA_2018','HallAB_2017','HeQ_2017','KumbhariA_2024',
  'LloydPriceJ_2019','SchirmerM_2018','SchirmerM_2024','WengY_2019','YanQ_2023c'
)

#### read tpm ####
tpm <- data.table::fread(args[1]) |>
  tibble::column_to_rownames("name") |>
  data.frame(check.names = FALSE)

#### Procrustes ####
result <- purrr::map_dfr(
  proj_name, \(proj) {
     
    group_x <- read.delim(
      paste0(kk_base, proj, ".sample_group.gz"),
      check.names = FALSE
    )
  
    ## family TPM profile
    tpm_x <- tpm |>
      dplyr::select(dplyr::any_of(group_x$sample))
  
    tpm_x <- tpm_x[
      rowSums(tpm_x, na.rm = TRUE) != 0,
      colSums(tpm_x, na.rm = TRUE) != 0,
      drop = FALSE
    ]
  
    data_x <- vegan::decostand(
      t(tpm_x), 
      method = "hellinger") |>
      t() |>
      data.frame(check.names = FALSE)
  
    ## species profile
    data_y <- data.table::fread(paste0(kk_base, proj, ".kk2.s.gz")) |>
      tibble::column_to_rownames("name") |>
      data.frame(check.names = FALSE) |>
      dplyr::select(dplyr::any_of(colnames(data_x)))
  
    data_y <- data_y[
      rowSums(data_y, na.rm = TRUE) != 0,
      colSums(data_y, na.rm = TRUE) != 0,
      drop = FALSE
    ]
    
    data_y <- vegan::decostand(t(data_y), method = "hellinger") |>
      t() |>
      data.frame(check.names = FALSE)
  
    ## 保证两个矩阵样本完全一致
    common_sample <- intersect(colnames(data_x), colnames(data_y))
  
    if (length(common_sample) < 3) {
      warning(proj, ": fewer than 3 common samples.")
      return(data.frame(name = proj, n = length(common_sample), M2 = NA, r = NA, pval = NA))
    }
  
    data_x <- data_x[, common_sample, drop = FALSE]
    data_y <- data_y[, common_sample, drop = FALSE]
    
    dist_x <- vegan::vegdist(t(data_x), method = "bray")
    dist_y <- vegan::vegdist(t(data_y), method = "bray")
    
    PCoA_x <- stats::cmdscale(dist_x, k = 2)
    PCoA_y <- stats::cmdscale(dist_y, k = 2)
    
    set.seed(2026)
    proc <- vegan::procrustes(PCoA_x, PCoA_y, symmetric = TRUE)
    proc_test <- vegan::protest(PCoA_x, PCoA_y, permutations = 1)
  
    data.frame(
      name = proj,
      n = length(common_sample),
      M2 = round(proc$ss, 4),
      r = round(proc_test$t0, 4),
      pval = proc_test$signif
    )
  }
)

write.table(
  result,
  paste0(args[2], ".procrustes.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
