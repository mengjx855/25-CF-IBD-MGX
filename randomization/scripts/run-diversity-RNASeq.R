#!/data/software/miniconda3/envs/r-base-4.4/bin/Rscript
# author: Jin-Xin Meng <jinxmeng@zju.edu.cn>

args <- commandArgs(trailingOnly = T)

if (length(args) < 2) {
  message("Usage: Rscript run-diversity-RNASeq.R [tpm] [out_prefix]")
  message("auto read metadata from /data/mengjx/project/10.20250623_IBD_BAC_CF_Landscape/profile/")
  # status = 1 表示非正常退出（给 Linux 系统看的信号）
  # save = "no" 表示不保存工作空间
  q(save = "no", status = 1)
}

pacman::p_load(tidyverse)
source('/data/mengjx/R_proj/R_func/calcu_difference.R')

file_base = '/data/mengjx/project/10.20250623_IBD_BAC_CF_Landscape/profile/'

proj_name <- c('LloydPriceJ_2019','SchirmerM_2018')

#### read tpm ####
tpm <- data.table::fread(args[1]) |>
  tibble::column_to_rownames('name') |>
  data.frame(check.names = FALSE)

#### alpha coef ####
result <- map_dfr(
  proj_name, ~ {
    
    group <- read.delim(
      paste0(
        file_base, .x, '.exp.sample_group.gz'
        )
      ) %>% 
      dplyr::select(sample, group)
    
    comparisons <- combn(c('CD', 'UC', 'HC'), m = 2, simplify = F)
    comparisons <- set_names(
      comparisons, 
      map_vec(comparisons, ~ paste0(.x, collapse = '_vs_'))
    )
    
    profile <- select(tpm, any_of(group$sample))
    group <- group[group$sample %in% colnames(profile), , drop = FALSE]
  
    # shannon
    data <- data.frame(value = vegan::diversity(t(profile), 'shannon')) %>% 
      rownames_to_column('sample') %>% 
      left_join(group, by = 'sample') |> 
      mutate(
        group = factor(group, c('CD', 'UC', 'HC'))
      )
    
    test <- calcu_diff(data, value ~ group) %>% 
      mutate(
        proj = .x,
        plab = add_plab(pval),
        name = paste0(proj, '|', comparison)
      ) %>% 
      select(name, pval, plab)
    
    stat <- aggregate(value ~ group, data, mean) |> 
      column_to_rownames('group')
    
    stat <- map2_dfr(
      comparisons, names(comparisons), \(x, y) {
        tibble(
          name = paste0(.x, '|', y),
          log2FC = log2(stat[x[1],] / stat[x[2],])
        )
      }
    )
    
    shannon <- left_join(stat, test, by = 'name') %>% 
      add_column(type = 'shannon', .before = 1) 
  
    # richness
    data <- data.frame(value = colSums(profile > 0)) %>% 
      rownames_to_column('sample') %>% 
      left_join(group, by = 'sample') %>% 
      mutate(
        group = factor(group, c('CD', 'UC', 'HC'))
      )
    
    test <- calcu_diff(data, value ~ group) %>% 
      mutate(
        proj = .x,
        plab = add_plab(pval),
        name = paste0(proj, '|', comparison)
      ) %>% 
      select(name, pval, plab)
    
    stat <- aggregate(value ~ group, data, mean) |> 
      column_to_rownames('group')
    
    stat <- map2_dfr(
      comparisons, names(comparisons), \(x, y) {
        tibble(
          name = paste0(.x, '|', y),
          log2FC = log2(stat[x[1],] / stat[x[2],])
        )
      }
    )
    
    richness <- left_join(stat, test, by = 'name') %>% 
      add_column(type = 'richness', .before = 1) 
    
    bind_rows(shannon, richness)
  } 
)

write.table(
  result, paste0(args[2], '.alphaDiv.tsv'), 
  sep = '\t', quote = F, row.names = F
  )
