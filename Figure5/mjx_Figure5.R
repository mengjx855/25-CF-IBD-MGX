#### Jinxin Meng, 20251028, 20251215 ####
setwd('/data/mengjx/project/10.20250623_IBD_BAC_CF_Landscape/git/Figure5/')
pacman::p_load(tidyverse, ggpubr)
source('../scripts/calcu_difference.R')
source('../scripts/calcu_metafor.R')
source('../scripts/model_randomforest.R')
source('../scripts/plot_roc.R')

#### Fig. 5b ####
proj_name <- c('BushmanFD_2020','FranzosaEA_2018','HallAB_2017','HeQ_2017','KumbhariA_2024',
               'LloydPriceJ_2019','SchirmerM_2018','SchirmerM_2024','WengY_2019','YanQ_2023c')

tpm <- readRDS('../Figure3/CF.tpm.rds')

data <- map(proj_name, ~ {
  data <- tpm[[.x]]
  group <- read.delim(paste0('../data/', .x, '.sample_group.bz2')) %>% 
    dplyr::select(sample, group)
  comparisons <- map(setdiff(unique(group$group), 'HC'), ~ c(.x, 'HC'))
  
  mean_data <- data.frame(t(data)) %>% 
    mutate(group = group$group[match(rownames(.), group$sample)]) %>% 
    group_by(group) %>% 
    summarise_all(mean) %>% 
    ungroup %>% 
    column_to_rownames('group')
  
  sd_data <- data.frame(t(data)) %>%  
    mutate(group = group$group[match(rownames(.), group$sample)]) %>% 
    group_by(group) %>% 
    summarise_all(sd) %>% 
    ungroup %>% 
    column_to_rownames('group')
  
  n_data <- data.frame(t(data)) %>% 
    mutate(group = group$group[match(rownames(.), group$sample)]) %>%
    count(group) %>% 
    column_to_rownames('group')
  
  map(comparisons, \(x) 
      data.frame(m1 = unlist(mean_data[x[1],]),
                 sd1 = unlist(sd_data[x[1],]),
                 n1 = n_data[x[1],]) %>% 
        rownames_to_column('name') %>%
        add_column(m2 = unlist(mean_data[x[2],]),
                   sd2 = unlist(sd_data[x[2],]),
                   n2 = n_data[x[2],])) %>% 
    set_names(map_vec(comparisons, \(x) paste0(x, collapse = '_vs_')))  }) %>% 
  set_names(proj_name)

names <- map(data, \(x) map(x, \(y) pull(y, name))) %>% 
  flatten() %>% 
  reduce(~ intersect(.x, .y))

meta_in <- map(names, \(x) 
               map_dfr(proj_name, \(y)
                       map2_dfr(data[[y]], names(data[[y]]), \(a, b)
                                filter(a, name == x) %>%
                                  add_column(name2 = paste0(y, '.', b)))) %>% 
                 dplyr::select(-name) %>% 
                 rename(name = name2) %>% 
                 relocate(name)) %>% 
  set_names(names)

meta_out <- map_dfr(names, ~ calcu_metafor.1(meta_in[[.x]]) %>% 
                      mutate(padj = p.adjust(pval, method = 'BH')) %>% 
                      add_column(feature = .x, .before = 1))
openxlsx::write.xlsx(meta_out, 'meta-analysis.xlsx')

data <- meta_out %>% 
  dplyr::select(feature, estimate, ci_lb, ci_ub, padj) %>% 
  distinct() %>% 
  mutate(plab = add_plab(padj)) %>% 
  arrange(estimate) %>% 
  mutate(feature = factor(feature, rev(.$feature)),
         ypos = ifelse(estimate > 0, -.3, .3),
         enriched = ifelse(estimate > 0, 'IBD', 'HC'))

CF_level <- data$feature

ggplot(data) + 
  geom_errorbar(aes(x = feature, ymin = ci_lb, ymax = ci_ub), width = .4, linewidth = .3) +
  geom_point(aes(x = feature, y = estimate, fill = enriched), size = 3, shape = 21) +
  scale_fill_manual(values = c(IBD = '#fee722', HC = '#5f7dbe')) +
  geom_hline(yintercept = 0, linetype = 'longdash', linewidth = .3) +
  geom_text(aes(x = feature, y = ypos, label = plab), color = 'red', hjust = .5, vjust = 1, angle = 90) +
  labs(x = '', y = 'estimate ± CI (95%)') +
  theme_pubr() +
  theme(aspect.ratio = 1/3.6, 
        line = element_blank(),
        axis.text = element_text(size = 10, color = '#000000'),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        legend.position = 'right',
        axis.ticks.length = unit(1, 'mm'),
        axis.ticks = element_line(linewidth = .3, color = '#000000'),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = .8, color = '#000000', fill = 'transparent') ) +
  annotate('text', x = 5, y = -.8, label = 'Enriched in IBD patient', hjust = 0) +
  annotate('text', x = 70, y = .7, label = 'Enriched in healthy control', hjust = 1) +
  annotate('rect', xmin = .5, xmax = 46.5, ymin = -1, ymax = 1, alpha = .1, fill = '#fee722') +
  annotate('rect', xmin = 46.5, xmax = 75.5, ymin = -1, ymax = 1, alpha = .1, fill = '#5f7dbe')
ggsave('meta-analysis.pdf', width = 12, height = 5)

#### Fig. 5a ####
rc <- readRDS('../Figure3/CF.rc.rds')

difference <- map(proj_name, \(x) {
  data <- rc[[x]]
  data <- data[colSums(data) > 100]
  group <- read.delim(paste0('../data/', x, '.sample_group.bz2')) %>% 
    dplyr::select(sample, group) %>% 
    dplyr::filter(sample %in% colnames(data)) %>% 
    arrange(match(sample, colnames(data))) %>% 
    column_to_rownames('sample') %>% 
    mutate(group = factor(group))
  comparisons <- map(setdiff(unique(group$group), 'HC'), ~ c(.x, 'HC'))
  
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = data, colData = group, design = ~ group)
  des <- DESeq2::DESeq(dds)
  
  map(comparisons, ~
        data.frame(DESeq2::results(des, contrast = c('group', .x))) %>%
        rownames_to_column('name') %>%
        filter(!is.na(padj)) %>%
        mutate(enriched = ifelse(log2FoldChange > 1 & padj < 0.05, .x[1], 
                                 ifelse(log2FoldChange < -1 & padj < 0.05, .x[2], 'none')) ) ) %>% 
    set_names(map_vec(comparisons, ~ .x[1]))  } ) %>% 
  set_names(proj_name) %>% 
  list_flatten(name_spec = '{outer}|{inner}')
saveRDS(difference, 'DESeq2.rds')

data <- map2_dfr(difference, names(difference), ~ add_column(.x, dataset = .y)) %>% 
  mutate(plab = add_plab(padj, format = 4)) %>% 
  mutate(log2FoldChange = ifelse(log2FoldChange > 5, 5, log2FoldChange),
         log2FoldChange = ifelse(log2FoldChange < -5, -5, log2FoldChange))

dataset_level <- data.frame(dataset = unique(data$dataset)) %>% 
  mutate(type = str_split_i(dataset, '\\|', 2), 
         type = factor(type, c('CD', 'UC', 'IC', 'IBD'))) %>% 
  arrange(type) %>% 
  pull(dataset)

select(data, dataset, name, log2FoldChange) %>% 
  spread('name', 'log2FoldChange', fill = 0) %>% 
  gather('name', 'log2FoldChange', -dataset) %>% 
  filter(name %in% CF_level) %>% 
  mutate(dataset = factor(dataset, dataset_level),
         name = factor(name, rev(CF_level)) ) %>% 
  ggplot() +
  geom_tile(aes(name, dataset, fill = log2FoldChange)) +
  geom_text(aes(name, dataset, label = plab), filter(data, name %in% CF_level), 
            color = 'red', size = 3, inherit.aes = F) +
  scale_fill_gradientn(
    colors = c('#4F6F9D',"#94A8D3","#FFFFFF","#FEEF6B",'#F7DF5A'),
    values = scales::rescale(c(-5, -1, 0, 1, 5)),
    breaks = c(-5, -1, 0, 1, 5),
    labels = c('-5', '-1', '0', '1', '5'),
    name = 'Log2FoldChange') +
  labs(x = '', y = '') +
  coord_fixed() +
  theme(line = element_blank(),
        axis.text = element_text(size = 10, color = '#000000'),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.ticks.length = unit(2, 'mm'),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = .8, color = '#000000', fill = 'transparent') )
ggsave('DESeq2.pdf', width = 15, height = 5)

#### Fig. 5c ####
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

#### Fig. 5d ####
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

#### Fig. 5e ####
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

#### Fig. 5f ####
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

#### Fig. 5g ####
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

#### Fig. 5h ####
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
