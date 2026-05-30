#!/data/software/miniconda3/envs/r-base-4.4/bin/Rscript
# author: Jin-Xin Meng <jinxmeng@zju.edu.cn>

args <- commandArgs(trailingOnly = T)

if (length(args) < 2) {
  message("Usage: Rscript run-diversity.R [tpm] [out_prefix]")
  message("auto read metadata from /data/mengjx/project/10.20250623_IBD_BAC_CF_Landscape/profile/")
  # status = 1 表示非正常退出（给 Linux 系统看的信号）
  # save = "no" 表示不保存工作空间
  q(save = "no", status = 1)
}

pacman::p_load(tidyverse)
source('/data/mengjx/R_proj/R_func/calcu_difference.R')

calcu_adjusted_r2 <- function(adonis_object) {
  n_observations <- adonis_object$Df[3]+1
  d_freedom <- adonis_object$Df[1]
  r2 <- adonis_object$R2[1]
  vegan::RsquareAdj(r2, n_observations, d_freedom)
}

file_base = '/data/mengjx/project/10.20250623_IBD_BAC_CF_Landscape/profile/'

proj_name <- c(
  'BushmanFD_2020','FranzosaEA_2018','HallAB_2017','HeQ_2017','KumbhariA_2024',
  'LloydPriceJ_2019','SchirmerM_2018','SchirmerM_2024','WengY_2019','YanQ_2023c'
)

#### read tpm ####
tpm <- data.table::fread(args[1]) |>
  tibble::column_to_rownames('name') |>
  data.frame(check.names = FALSE)

#### alpha coef ####
result <- map_dfr(
  proj_name, ~ {
    group <- read.delim(
      paste0(file_base, .x, '.sample_group.gz')
    ) %>% 
      dplyr::select(sample, group)
    
    group_level <- intersect(c('HC','IBD','CD','UC','IC'), unique(group$group))
    
    profile <- select(tpm, any_of(group$sample))
    group <- group[group$sample %in% colnames(profile), , drop = FALSE]
  
    # shannon
    data <- data.frame(value = vegan::diversity(t(profile), 'shannon')) %>% 
      rownames_to_column('sample') %>% 
      left_join(group, by = 'sample') %>% 
      mutate(group = factor(group, group_level))
    
    test <- calcu_diff(data, value ~ group) %>% 
      filter(grepl('HC', comparison)) %>% 
      mutate(
        proj = .x,
        plab = add_plab(pval),
        name = str_remove_all(comparison, 'HC_vs_|_vs_HC'),
        name = paste0(proj, '|', name)) %>% 
      select(name, pval, plab)
    
    stat <- aggregate(value ~ group, data, mean)
    stat <- filter(stat, group != 'HC') %>% 
      rename(case = value) %>% 
      add_column(
        control = stat[stat$group == 'HC', 'value'],
        proj = .x
      ) %>% 
      mutate(
        name = paste0(proj, '|', group),
        log2FC = log2(case / control)
      ) %>% 
      select(name, log2FC)
    
    shannon = left_join(stat, test, by = 'name') %>% 
      add_column(type = 'shannon', .before = 1) 
  
    # richness
    data <- data.frame(value = colSums(profile > 0)) %>% 
      rownames_to_column('sample') %>% 
      left_join(group, by = 'sample') %>% 
      mutate(group = factor(group, group_level))
    
    test <- calcu_diff(data, value ~ group) %>% 
      add_column(proj = .x) %>% 
      filter(grepl('HC', comparison)) %>% 
      mutate(
        plab = add_plab(pval),
        name = str_remove_all(comparison, 'HC_vs_|_vs_HC'),
        name = paste0(proj, '|', name)
      ) %>% 
      select(name, pval, plab)
    
    stat <- aggregate(value ~ group, data, mean)
    stat <- filter(stat, group != 'HC') %>% 
      rename(case = value) %>% 
      add_column(
        control = stat[stat$group == 'HC', 'value'],
        proj = .x
      ) %>% 
      mutate(
        name = paste0(proj, '|', group),
        log2FC = log2(case / control)
      ) %>% 
      select(name, log2FC)
    
    richness = left_join(stat, test, by = 'name') %>% 
      add_column(type = 'richness', .before = 1) 
    
    bind_rows(shannon, richness)
  } )

write.table(
  result, paste0(args[2], '.alphaDiv.tsv'), 
  sep = '\t', quote = F, row.names = F
  )

#### adonis ####
result <- map_dfr(
  proj_name, ~ {
  
    group <- read.delim(
      paste0(file_base, .x, '.sample_group.gz')
    ) %>% 
      dplyr::select(sample, group)
    
    group_level <- intersect(c('HC','IBD','CD','UC','IC'), unique(group$group))
  
    test <- map_dfr(
      setdiff(group_level, 'HC'), \(y) {
        .group <- filter(group, group %in% c(y, 'HC'))
      
        .tpm <- select(tpm, any_of(.group$sample))
        .tpm <- .tpm[
          rowSums(.tpm, na.rm = TRUE) != 0,
          colSums(.tpm, na.rm = TRUE) != 0,
          drop = FALSE
          ]
      
        .group <- .group[match(colnames(.tpm), .group$sample), , drop = FALSE]
        .group$group <- factor(.group$group, levels = c('HC', y))
        
        
        adonis <- vegan::adonis2(
          t(.tpm) ~ group, .group, permutations = 1, 
          distance = 'bray', parallel = 1)
      
        r2adj <- calcu_adjusted_r2(adonis)
        
        data.frame(
          name = paste0(.x, '|', y), 
          r2 = adonis$R2[1], 
          r2adj = r2adj, 
          pval = adonis$`Pr(>F)`[1]
        )
      }
    )
  }
)

write.table(
  result, paste0(args[2], '.betaDiv.tsv'),
  sep = '\t', quote = F, row.names = F
)
