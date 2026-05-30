#!/data/software/miniconda3/envs/r-base-4.4/bin/Rscript
# author: Jin-Xin Meng <jinxmeng@zju.edu.cn>

args <- commandArgs(trailingOnly = T)

if (length(args) < 3) {
  message("Usage: Rscript run-DNA-RNA-correlation.R [CF-DNASeq-tpm] [CF-RNASeq-tpm] [out_prefix]")
  message("auto read metadata from /data/mengjx/project/10.20250623_IBD_BAC_CF_Landscape/profile/")
  # status = 1 表示非正常退出（给 Linux 系统看的信号）
  # save = "no" 表示不保存工作空间
  q(save = "no", status = 1)
}

pacman::p_load(tidyverse)

file_base = '/data/mengjx/project/10.20250623_IBD_BAC_CF_Landscape/profile/'

proj_name <- c('LloydPriceJ_2019','SchirmerM_2018')

#### read tpm ####
DNA_tpm <- data.table::fread(args[1]) %>% 
  column_to_rownames('name')

RNA_tpm <- data.table::fread(args[2]) %>% 
  column_to_rownames('name')

#### rf ####
result <- map_dfr(
  proj_name, ~ {
    
    group <- read.delim(
      paste0(
        file_base, 
        .x, 
        '.exp.sample_group.gz'
      )
    ) |> 
      filter(!is.na(sample_MGX))
    
    MTX_profile <- log10(RNA_tpm[group$sample] + 1)
    MGX_profile <- log10(DNA_tpm[rownames(MTX_profile), group$sample_MGX] + 1)
    
    data <- data.frame(
      name = rownames(MTX_profile),
      MTX_mean = rowMeans(MTX_profile),
      MGX_mean = rowMeans(MGX_profile),
      row.names = NULL
    )
    
    test <- cor.test(data$MTX_mean, data$MGX_mean)
    
    tibble(
      proj = .x,
      r = as.numeric(test$estimate),
      p = test$p.value
    )
  }
)

write.table(
  result, 
  paste0(args[3], '.DNA_vs_RNA.corr.tsv'), 
  sep = '\t', quote = F, row.names = F
)
