#### Jinxin Meng, 20251028, 20251109 ####
setwd('/data/mengjx/project/10.20250623_IBD_BAC_CF_Landscape/git/Figure6/')
pacman::p_load(tidyverse, ggpubr)
source('../scripts/model_randomforest.R')
source('../scripts/plot_roc.R')

#### Fig. 6a ####
proj_name <- c('BushmanFD_2020','FranzosaEA_2018','HallAB_2017','HeQ_2017','KumbhariA_2024',
               'LloydPriceJ_2019','SchirmerM_2018','SchirmerM_2024','WengY_2019','YanQ_2023c')

tpm <- readRDS('../Figure3/CF.tpm.rds')

roc <- map(proj_name, ~ {
  group <- read.delim(paste0('../data/', .x, '.sample_group.bz2')) %>% 
    dplyr::select(sample, group) %>% 
    mutate(group = ifelse(group == 'HC', 'HC', 'IBD'))
  pred <- rf_Kfold(tpm[[.x]], group, k = 10, seed = 2025)
  pred <- left_join(pred, group, by = 'sample')
  roc(pred$group, pred$HC) }) %>% 
  set_names(proj_name)

plot_roc_multiple(roc, title = 'Random Forest Model (All CF)')
ggsave('CF.rf.roc.all_CF.pdf', width = 8, height = 4)

map2_dfr(roc, names(roc), ~ data.frame(name = .y, auc = as.numeric(.x$auc))) %>% 
  add_column(type = 'ALL_CF') %>% 
  write_tsv('CF.rf.roc.all_CF.tsv')

#### Fig. 6b ####
meta_out <- openxlsx::read.xlsx('../Figure5/meta-analysis.xlsx')
names <- filter(meta_out, padj < 0.05) %>% 
  pull(feature) %>% unique()

roc <- map(proj_name, ~ {
  group <- read.delim(paste0('../data/', .x, '.sample_group.bz2')) %>% 
    dplyr::select(sample, group) %>% 
    mutate(group = ifelse(group == 'HC', 'HC', 'IBD'))
  pred <- rf_Kfold(tpm[[.x]][names, ], group, k = 10, seed = 2025)
  pred <- left_join(pred, group, by = 'sample')
  roc(pred$group, pred$HC) }) %>% 
  set_names(proj_name)

plot_roc_multiple(roc, title = 'Random Forest Model (46 CF signatures)')
ggsave('CF.rf.roc.46_CF.pdf', width = 8, height = 4)

map2_dfr(roc, names(roc), ~ data.frame(name = .y, auc = as.numeric(.x$auc))) %>% 
  add_column(type = '46_CF') %>% 
  write_tsv('CF.rf.roc.46_CF.tsv')

#### Fig. 6c ####
roc <- map(proj_name, ~ {
  group <- read.delim(paste0('../data/', .x, '.sample_group.bz2')) %>% 
    dplyr::select(sample, group) %>% 
    mutate(group = ifelse(group == 'HC', 'HC', 'IBD'))
  profile <- data.table::fread(paste0('../data/', .x, '.kk2.s.bz2')) %>% 
    column_to_rownames('name')
  pred <- rf_Kfold(profile, group, k = 10, seed = 2025)
  pred <- left_join(pred, group, by = 'sample')
  roc(pred$group, pred$HC) }) %>% 
  set_names(proj_name)

plot_roc_multiple(roc, title = 'Random Forest Model (Species taxa)')
ggsave('CF.rf.roc.species.pdf', width = 8, height = 4)

map2_dfr(roc, names(roc), ~ data.frame(name = .y, auc = as.numeric(.x$auc))) %>% 
  add_column(type = 'species') %>% 
  write_tsv('ML/CF.rf.roc.species.tsv')

#### Fig. 6d ####
roc <- map(proj_name, ~ {
  group <- read.delim(paste0('../data/', .x, '.sample_group.bz2')) %>% 
    dplyr::select(sample, group) %>% 
    mutate(group = ifelse(group == 'HC', 'HC', 'IBD'))
  profile <- data.table::fread(paste0('../data/', .x, '.kk2.g.bz2')) %>% 
    column_to_rownames('name')
  pred <- rf_Kfold(profile, group, k = 10, seed = 2025)
  pred <- left_join(pred, group, by = 'sample')
  roc(pred$group, pred$HC) }) %>% 
  set_names(proj_name)

plot_roc_multiple(roc, title = 'Random Forest Model (Genus taxa)')
ggsave('CF.rf.roc.genus.pdf', width = 8, height = 4)

map2_dfr(roc, names(roc), ~ data.frame(name = .y, auc = as.numeric(.x$auc))) %>% 
  add_column(type = 'genus') %>% 
  write_tsv('ML/CF.rf.roc.genus.tsv')

#### Fig. 6e ####
list(
  read.delim('CF.rf.roc.all_CF.tsv') %>% 
    select(name = 1, all_CF = 2),
  read.delim('CF.rf.roc.46_CF.tsv') %>% 
    select(name = 1, markers_CF = 2),
  read.delim('CF.rf.roc.species.tsv') %>% 
    select(name = 1, species = 2),
  read.delim('CF.rf.roc.genus.tsv') %>% 
    select(name = 1, genus = 2)) %>% 
  purrr::reduce(~ left_join(.x, .y, by = 'name')) %>% 
  mutate(
    all_CF = round(all_CF, 4),
    markers_CF = round(markers_CF, 4),
    species = round(species, 4),
    genus = round(genus, 4)
    ) %>% 
  ggtexttable(theme = ttheme('light'))
ggsave('CF.rf.roc.table.pdf', width = 6, height = 5)

#### Fig. 6f ####
model <- map_dfr(proj_name, \(x) {
  group_x <- read.delim(paste0('../data/', x, '.sample_group.bz2')) %>% 
    dplyr::select(sample, group) %>% 
    mutate(group = ifelse(group == 'HC', 'HC', 'IBD'))
  profile_x <- tpm[[x]][names, ]
  map_dfr(proj_name, \(y) {
    group_y <- read.delim(paste0('../data/', y, '.sample_group.bz2')) %>% 
      dplyr::select(sample, group) %>% 
      mutate(group = ifelse(group == 'HC', 'HC', 'IBD'))
    profile_y <- tpm[[y]][names, ]
    pred <- rf_next_vaildate(profile_x, profile_y, group_x, group_y, seed = 2025)
    pred <- left_join(pred, group_y, by = 'sample')
    roc <- roc(pred$group, pred$HC) 
    data.frame(name_y = y, auc = as.numeric(roc$auc)) } ) %>% 
    add_column(name_x = x, .before = 1) })
write_tsv(model, 'CF.rf.roc.cross.tsv')

model <- rbind(
  read.delim('CF.rf.roc.cross.tsv') %>% 
    filter(name_x != name_y),
  read.delim('CF.rf.roc.46_CF.tsv',) %>% 
    select(name_x = name, auc) %>% 
    mutate(name_y = name_x, .before = 'auc')
)

ggplot(model, aes(name_x, name_y)) +
  geom_tile(aes(fill = auc)) +
  labs(x = 'Discovery dataset', y = 'Validation dataset', fill = 'AUC') +
  geom_text(aes(label = round(auc, 2)), size = 2) +
  scale_fill_gradientn(colours = c('#4F6F9D',"#94A8D3","#FFFFFF","#FEEF6B",'#F7DF5A')) +
  coord_fixed() +
  theme(axis.ticks = element_blank(),
        axis.text = element_text(color = '#000000', size = 10),
        axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = .8, color = 'black', fill = 'transparent'))
ggsave('CF.rf.roc.cross.pdf', width = 7, height = 6)
