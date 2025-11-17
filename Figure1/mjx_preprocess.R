#### Jinxin Meng, 20251028, 20251109 ####
setwd('/data/mengjx/project/10.20250623_IBD_BAC_CF_Landscape/git/Figure1/')
pacman::p_load(tidyverse, ggpubr)
source('../scripts/palette.R')
source('../scripts/transform_rc.R')

#### profile ####
proj_name <- c('BushmanFD_2020','FranzosaEA_2018','HallAB_2017','HeQ_2017','KumbhariA_2024',
               'LloydPriceJ_2019','SchirmerM_2018','SchirmerM_2024','WengY_2019','YanQ_2023c')

gene_len <- read.delim('../pipeline/uhgp.len.bz2', header = F, col.names = c('name', 'length'))
gene_info <- read.delim('../pipeline/uhgp.m8.f.bz2', header = F) %>%
  dplyr::select(name = 1, value = 2) %>%
  mutate(cFam = stringr::str_split_i(value, ',', 1),
         cGene = stringr::str_split_i(value, ',', 2))

# CF_gene tpm
tpms <- map(proj_names, ~ {
  rc <- fread(paste0('../profile/', .x, '.rc')) %>% column_to_rownames('name')
  cvg <- fread(paste0('../profile/', .x, '.cvg')) %>% column_to_rownames('name')
  tpm <- rc2tpm(rc, gene_len) * (cvg > 50)
  tpm[rowSums(tpm) != 0, colSums(tpm) != 0] } ) %>%
  set_names(proj_names)
saveRDS(tpms, 'cRaw.tpms.rds')

# cRaw rcs
rcs <- map(proj_names, ~ {
  cvg <- fread(paste0('../profile/', .x, '.cvg')) %>% column_to_rownames('name')
  rc <- fread(paste0('../profile/', .x, '.rc')) %>% column_to_rownames('name') * (cvg > 50)
  rc[rowSums(rc) != 0, colSums(rc) != 0] } ) %>%
  set_names(proj_names)
saveRDS(rcs, 'cRaw.rcs.rds')

# cGene tpms
tpms <- map(proj_names, ~ {
  rc <- fread(paste0('../profile/', .x, '.rc')) %>% column_to_rownames('name')
  cvg <- fread(paste0('../profile/', .x, '.cvg')) %>% column_to_rownames('name')
  tpm <- rc2tpm(rc, gene_len) * (cvg > 50)
  tpm[colSums(tpm) != 0] %>%
    mutate(name = gene_info$cGene[match(rownames(tpm), gene_info$name)]) %>%
    aggregate(. ~ name, ., sum) %>%
    column_to_rownames('name') %>%
    filter(rowSums(.) != 0) } ) %>%
  set_names(proj_names)
saveRDS(tpms, 'cGene.tpms.rds')

# cGene rcs
rcs <- map(proj_names, ~ {
  cvg <- fread(paste0('../profile/', .x, '.cvg')) %>% column_to_rownames('name')
  rc <- fread(paste0('../profile/', .x, '.rc')) %>% column_to_rownames('name') * (cvg > 50)
  mutate(rc, name = gene_info$cGene[match(rownames(rc), gene_info$name)]) %>%
    aggregate(. ~ name, ., sum) %>%
    column_to_rownames('name') %>%
    filter(rowSums(.) != 0) } ) %>%
  set_names(proj_names)
saveRDS(rcs, 'cGene.rcs.rds')

# cFam tpms
tpms <- map(proj_names, ~ {
  rc <- fread(paste0('../profile/', .x, '.rc')) %>% column_to_rownames('name')
  tpm <- rc2tpm(rc, gene_len)
  tpm[colSums(tpm) != 0] %>%
    mutate(name = gene_info$cFam[match(rownames(tpm), gene_info$name)]) %>%
    aggregate(. ~ name, ., sum) %>%
    column_to_rownames('name') %>%
    filter(rowSums(.) != 0) } ) %>%
  set_names(proj_names)
saveRDS(tpms, 'cFam.tpms.rds')

# cFam rcs
rcs <- map(proj_names, ~ {
  rc <- fread(paste0('../profile/', .x, '.rc')) %>% column_to_rownames('name')
  mutate(rc, name = gene_info$cFam[match(rownames(rc), gene_info$name)]) %>%
    aggregate(. ~ name, ., sum) %>%
    column_to_rownames('name') %>%
    filter(rowSums(.) != 0) } ) %>%
  set_names(proj_names)
saveRDS(rcs, 'cFam.rcs.rds')

#### contributor ####
proj_names <- c('BushmanFD_2020','FranzosaEA_2018','HallAB_2017','HeQ_2017','KumbhariA_2024',
                'LloydPriceJ_2019','SchirmerM_2018','SchirmerM_2024','WengY_2019','YanQ_2023c')

gene_len <- read.delim('../data/uhgp-90.len', header = F, col.names = c('name', 'len'))
gene_info <- read.delim('../data/uhgp-90.m8', header = F) %>%
  dplyr::select(name = 1, value = 2) %>%
  mutate(cFam = stringr::str_split_i(value, ',', 1),
         cGene = stringr::str_split_i(value, ',', 2))
taxonomy <- read.delim('../pipeline/uhgp-90.m8.genome_info.tsv', header = F) %>% 
  dplyr::select(name = V1, taxonomy = V19) %>% 
  taxa_split()

data <- map(proj_names, ~ {
  rc <- fread(paste0('../profile/', .x, '.rc')) %>% column_to_rownames('name')
  cvg <- fread(paste0('../profile/', .x, '.cvg')) %>% column_to_rownames('name')
  group <- read.delim(paste0('../profile/', .x, '.sample_group')) %>% 
    dplyr::select(sample, group)
  tpm <- profile_smp2grp(rc2tpm(rc, gene_len) * (cvg > 50), group) %>% 
    profile_transRA() %>% 
    filter(rowSums(.) != 0) %>% 
    mutate(cFam = gene_info$cFam[match(rownames(.), gene_info$name)]) %>%
    mutate(taxa = taxonomy$genus[match(rownames(.), taxonomy$name)]) %>%
    aggregate(. ~ cFam + taxa, ., sum) %>% 
    gather('group', 'value', -cFam, -taxa)
  
} ) %>%
  set_names(proj_names)
saveRDS(tpms, 'cRaw.tpms.rds')