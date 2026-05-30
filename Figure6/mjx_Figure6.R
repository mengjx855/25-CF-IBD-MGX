#### Jinxin Meng, 20251028, 20260530 ####
setwd('/data/mengjx/project/10.20250623_IBD_BAC_CF_Landscape/github/Figure6/')
pacman::p_load(tidyverse, ggpubr)
source('../scripts/calcu_difference.R')
source('../scripts/calcu_metafor.R')
source('../randomization/scripts/model_randomforest.R')
source('../randomization/scripts/plot_roc.R')

#### Fig. 6b ####
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

#### Fig. 6a ####
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

#### Fig. 6c ####
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
ggsave('CF.rf.roc.all_CF_gene.pdf', width = 8, height = 4)

map2_dfr(roc, names(roc), ~ data.frame(name = .y, auc = as.numeric(.x$auc))) %>% 
  add_column(type = 'ALL_CF') %>% 
  write_tsv('CF.rf.roc.all_CF.tsv')

#### Fig. 6d ####
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

#### Fig. S6f ####
tpm <- readRDS('../Figure4/CF.tpm.rds')

importance <- map(
  proj_name, ~ {
    
    group <- read.delim(
      paste0(
        '../data/',
        .x, 
        '.sample_group.bz2')
    ) %>% 
      dplyr::select(sample, group) %>% 
      mutate(group = ifelse(group == 'HC', 'HC', 'IBD'))
    
    rf_importance(
      tpm[[.x]][, group$sample], group, seed = 2026, ntree = 1000
    ) |> 
      mutate(
        rank = row_number()
      )
    
  }, .progress = T
) %>% 
  set_names(proj_name)

rank_df <- map2(
  importance, names(importance), ~ 
    select(.x, name, !!.y := MeanDecreaseAccuracy)
) |> 
  reduce(~ full_join(.x, .y, by = 'name') ) |> 
  mutate(
    across(
      !name, ~ 
        replace_na(.x, 0)
    )
  )

write.xlsx(
  rank_df, 'feature_importance.xlsx'
)

feature_names <- column_to_rownames(rank_df, 'name') |> 
  rowMeans() |> 
  sort(decreasing = T) |> 
  names()
write_rds(feature_names, 'feature_importance.names.rds')

# 汇总
files <- list.files('feature_importance/', pattern = 'top')

results <- map_dfr(
  files, ~ 
    read.delim(
      paste0('feature_importance/', .x), 
      header = F, col.names = c('top_n', 'proj', 'seed', 'auc')
    )
)

auc_df <- group_by(results, top_n, proj) |> 
  summarise(avg = mean(auc), .groups = 'drop')

auc_df2 <- group_by(auc_df, top_n) |> 
  summarise(mean = mean(avg), .groups = 'drop') |> 
  mutate(
    top_n = sub('top_', '', top_n),
    top_n = as.numeric(top_n)
  )

list(
  auc_by_proj = auc_df, 
  auc_by_top_n = auc_df2
) |> 
  write.xlsx('feature_topN_AUC.xlsx')

ggplot(auc_df2, aes(x = top_n, y = mean)) +
  geom_vline(xintercept = 13, linetype = 'longdash', linewidth = .4) +
  geom_hline(
    yintercept = auc_df2$mean[auc_df2$top_n == 13], 
    linetype = 'longdash', linewidth = .4
  ) +
  geom_line(linewidth = .5, color = 'black') +
  geom_point(size = 3, shape = 21, fill = 'white', color = 'black') +
  scale_x_continuous(
    breaks = seq(0, 79, 10),
    labels = seq(0, 79, 10),
    expand = c(.02, .02)
  ) +
  labs(
    x = 'Top N CF-family', y = 'Average AUC', 
    title = 'Optimal number of CF-family'
  ) +
  theme_pubr() +
  theme(aspect.ratio = 1/4,
        axis.ticks.length = unit(2, 'mm'),
        plot.title = element_text(color = 'black', hjust = .5, face = 'bold'))

ggsave('feature_topN_AUC.pdf', width = 10, height = 3)

#### Fig. 6e ####
tpm <- read_rds('../Figure4/CF.tpm.rds')

feature_names <- read_rds('feature_importance.names.rds')[1:13]

roc <- map(
  proj_name, ~ {
    
    group <- read.delim(
      paste0(
        '../data/',
        .x, 
        '.sample_group.bz2')
    ) %>% 
      dplyr::select(sample, group) %>% 
      mutate(group = ifelse(group == 'HC', 'HC', 'IBD'))
    
    pred <- rf_Kfold(
      tpm[[.x]][feature_names, group$sample], 
      group, k = 10, seed = 2026)
    
    pred <- left_join(pred, group, by = 'sample')
    
    pROC::roc(pred$group, pred$HC)
    
  }
) %>% 
  set_names(proj_name)

plot_roc_multiple(roc, title = 'Random Forest Model (Top 13 CF-family)')

ggsave('roc.top_13.pdf', width = 8, height = 4)

results <- map2_dfr(
  roc, names(roc), ~ 
    data.frame(
      name = .y, 
      auc = as.numeric(.x$auc))) %>% 
  add_column(type = 'top_13')

write_tsv(results, 'CF.rf.roc.top_13.tsv')

#### Fig. 6f ####
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

#### Fig. 6g ####
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

#### Fig. 6h ####
df <- list(
  read.delim('CF.rf.roc.all_CF.tsv') %>% 
    select(name = 1, CF79 = 2),
  read.delim('CF.rf.roc.46_CF.tsv') %>% 
    select(name = 1, CF46 = 2),
  read.delim('CF.rf.roc.top_13.tsv') |> 
    select(name = 1, CF13 = 2),
  read.delim('CF.rf.roc.species.tsv') %>% 
    select(name = 1, species = 2),
  read.delim('CF.rf.roc.genus.tsv') %>% 
    select(name = 1, genus = 2)) %>% 
  purrr::reduce(~ left_join(.x, .y, by = 'name')) %>% 
  mutate(
    across(
      where(is.numeric), ~
        round(.x, 3)
    )
  )

df2 <- data.frame(
  average = colMeans(df[,-1]) |> round(3), 
  sd = apply(df[, -1], 2, \(x) sd(x)) |> round(3),
  number = c(79, 46, 13, 4324, 1059)
) |> 
  t() |> 
  data.frame() |> 
  rownames_to_column('name')

df <- bind_rows(df, df2) |> 
  rename(project = name)

ggtexttable(
  df, theme = ttheme('blank'), rows = NULL
) |> 
  tab_add_vline(
    at.column = 2:6, column.side = 'left', 
    from.row = 2, linetype = 2
  ) |> 
  tab_add_hline(
    at.row = c(1, 2), row.side = 'top',
    linewidth = 2, linetype = 1
  ) |> 
  tab_add_hline(
    at.row = c(11, 14), row.side = 'bottom',
    linewidth = 2, linetype = 1)

ggsave('CF.rf.roc.summary.pdf', width = 6, height = 5)

#### Fig. S6m ####
# top most informative 13/46/79 genus
proj_name <- c(
  'BushmanFD_2020','FranzosaEA_2018','HallAB_2017','HeQ_2017','KumbhariA_2024',
  'LloydPriceJ_2019','SchirmerM_2018','SchirmerM_2024','WengY_2019','YanQ_2023c'
)

importance <- map(
  proj_name, ~ {
    
    group <- read.delim(
      paste0(
        '../data/', .x, '.sample_group.bz2')
    ) %>% 
      dplyr::select(sample, group) %>% 
      mutate(group = ifelse(group == 'HC', 'HC', 'IBD'))
    
    kk2 <- data.table::fread(
      paste0('../data/', .x, '.kk2.g.bz2')
    ) |> 
      column_to_rownames('name')
    
    rf_importance(
      kk2[, group$sample], group, seed = 2026, ntree = 1000
    )
    
  }, .progress = T
) %>% 
  set_names(proj_name)

rank_df <- map2(
  importance, names(importance), ~ 
    select(.x, name, !!.y := MeanDecreaseAccuracy)
) |> 
  reduce(~ full_join(.x, .y, by = 'name') ) |> 
  mutate(
    across(
      !name, ~ 
        replace_na(.x, 0)
    )
  )

feature_names <- column_to_rownames(rank_df, 'name') |> 
  rowMeans() |> 
  sort(decreasing = T) |> 
  names()

roc_df <- map_dfr(
  c(13, 46, 79), \(x) {
    map_dfr(
      proj_name, \(y) {
        
        group <- read.delim(
          paste0('../data/', y, '.sample_group.bz2')
        ) %>% 
          dplyr::select(sample, group) %>% 
          mutate(group = ifelse(group == 'HC', 'HC', 'IBD'))
        
        kk2 <- data.table::fread(
          paste0('../data/', y, '.kk2.g.bz2')
        ) |> 
          column_to_rownames('name')
        
        pred <- rf_Kfold(
          kk2[feature_names[1:x], ],
          group, k = 10, seed = 2026, quiet = T
        )
        
        pred <- left_join(pred, group, by = 'sample')
        
        roc <- pROC::roc(pred$group, pred$HC, quiet = T)
        
        tibble(
          size = x, 
          proj = y, 
          auc = round(roc$auc, 6)
        )
      }, .progress = T
    )
  }
)

write.xlsx(roc_df, 'taxa_genus_top_roc.xlsx')

# top most informative 13/46/79 species
importance <- map(
  proj_name, ~ {
    
    group <- read.delim(
      paste0(
        '../data/', .x, '.sample_group.bz2')
    ) %>% 
      dplyr::select(sample, group) %>% 
      mutate(group = ifelse(group == 'HC', 'HC', 'IBD'))
    
    kk2 <- data.table::fread(
      paste0('../data/', .x, '.kk2.s.bz2')
    ) |> 
      column_to_rownames('name')
    
    rf_importance(
      kk2[, group$sample], group, seed = 2026, ntree = 1000
    )
    
  }, .progress = T
) %>% 
  set_names(proj_name)

rank_df <- map2(
  importance, names(importance), ~ 
    select(.x, name, !!.y := MeanDecreaseAccuracy)
) |> 
  reduce(~ full_join(.x, .y, by = 'name') ) |> 
  mutate(
    across(
      !name, ~ 
        replace_na(.x, 0)
    )
  )

feature_names <- column_to_rownames(rank_df, 'name') |> 
  rowMeans() |> 
  sort(decreasing = T) |> 
  names()

roc_df <- map_dfr(
  c(13, 46, 79), \(x) {
    map_dfr(
      proj_name, \(y) {
        
        group <- read.delim(
          paste0('../data/', y, '.sample_group.bz2')
        ) %>% 
          dplyr::select(sample, group) %>% 
          mutate(group = ifelse(group == 'HC', 'HC', 'IBD'))
        
        kk2 <- data.table::fread(
          paste0('../data/', y, '.kk2.s.bz2')
        ) |> 
          column_to_rownames('name')
        
        pred <- rf_Kfold(
          kk2[feature_names[1:x], ],
          group, k = 10, seed = 2026, quiet = T
        )
        
        pred <- left_join(pred, group, by = 'sample')
        
        roc <- pROC::roc(pred$group, pred$HC, quiet = T)
        
        tibble(
          size = x, 
          proj = y, 
          auc = round(roc$auc, 6)
        )
      }, .progress = T
    )
  }
)

write.xlsx(roc_df, 'taxa_species_top_roc.xlsx')

# comparison
plot_data <- bind_rows(
  
  read.xlsx('taxa_species_top_roc.xlsx') |> 
    add_column(type = 'species'),
  
  read.xlsx('taxa_genus_top_roc.xlsx') |> 
    add_column(type = 'genus'),
  
  read.delim('CF.rf.roc.top_13.tsv') |> 
    select(size = type, proj = name, auc) |> 
    mutate(
      size = 13, 
      type = 'CF'
    ), 
  
  read.delim('CF.rf.roc.46_CF.tsv') |> 
    select(size = type, proj = name, auc) |> 
    mutate(
      size = 46,
      type = 'CF'
    ),
  read.delim('CF.rf.roc.all_CF.tsv') |> 
    select(size = type, proj = name, auc) |> 
    mutate(
      size = 79,
      type = 'CF'
    )
) |> 
  group_by(type, size) |> 
  summarise(
    mean = mean(auc),
    sd = sd(auc), 
    .groups = 'drop'
  )

ggplot(plot_data, aes(x = type, y = mean)) +
  geom_errorbar(
    aes(ymin = mean - sd, ymax = mean + sd),
    width = .25, linewidth = .4, color = 'black'
  ) +
  geom_point(
    aes(fill = type),
    shape = 21, size = 5, stroke = .4, color = 'black'
  ) +
  facet_grid(cols = vars(size)) +
  labs(
    x = '', y = 'Average AUC', fill = 'feature type'
  ) +
  theme_bw() +
  theme(
    aspect.ratio = 3,
    axis.text = element_text(size = 12, color = 'black'),
    axis.text.x = element_text(size = 12, color = 'black', angle = 20),
    axis.title = element_text(size = 12, color = 'black'),
    axis.ticks = element_line(linewidth = .5, color = 'black'),
    axis.ticks.length = grid::unit(2, 'mm'),
    panel.grid.major = element_line(color = 'grey88', linewidth = .4),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = .5, color = 'black', fill = NA),
    strip.background = element_rect(fill = 'grey95', color = 'black', linewidth = .5),
    strip.text = element_text(size = 11, color = 'black'),
    legend.position = 'top',
    legend.title = element_text(size = 11, color = 'black'),
    legend.text = element_text(size = 11, color = 'black')
  )

ggsave('taxa_top_informative.pdf', width = 7, height = 4)

#### Fig. S6n ####
CF79 <- read_rds('feature_importance.names.rds')

model <- map_dfr(
  proj_name, \(x) {
    
    group_x <- read.delim(
      paste0('../data/', x, '.sample_group.bz2')
    ) %>% 
      dplyr::select(sample, group) %>% 
      mutate(group = ifelse(group == 'HC', 'HC', 'IBD'))
    
    profile_x <- tpm[[x]][CF79, ]
    
    map_dfr(
      proj_name, \(y) {
        group_y <- read.delim(
          paste0('../data/', y, '.sample_group.bz2')
        ) %>% 
          dplyr::select(sample, group) %>% 
          mutate(group = ifelse(group == 'HC', 'HC', 'IBD'))
        
        profile_y <- tpm[[y]][CF79, ]
        pred <- rf_next_vaildate(
          profile_x, profile_y, group_x, group_y, seed = 2025
        )
        
        pred <- left_join(pred, group_y, by = 'sample')
        roc <- pROC::roc(pred$group, pred$HC, quiet = T) 
        data.frame(
          name_x = x, 
          name_y = y, 
          auc = as.numeric(roc$auc)
        )
      }
    ) 
  }, .progress = T
)

write_tsv(model, 'cross.auc.CF79.tsv')

model <- rbind(
  read.delim('cross.auc.CF79.tsv') %>% 
    filter(name_x != name_y),
    read.delim('CF.rf.roc.all_CF.tsv') %>%
    select(name_x = name, auc) %>% 
    mutate(name_y = name_x, .before = 'auc')
)

ggplot(model, aes(name_x, name_y)) +
  geom_tile(aes(fill = auc)) +
  labs(
    x = 'Discovery dataset', y = 'Validation dataset', fill = 'AUC',
    title = '79 CF-family cross-cohort validation') +
  geom_text(aes(label = round(auc, 2)), size = 2) +
  scale_fill_gradientn(
    colors = c('#4F6F9D',"#94A8D3","#FFFFFF","#FEEF6B",'#F7DF5A')
  ) +
  coord_fixed() +
  theme(
    axis.ticks = element_blank(),
    axis.text = element_text(color = '#000000', size = 10),
    axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1),
    plot.title = element_text(hjust = .5, face = 'bold'),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(linewidth = .8, color = 'black', fill = NA)
  )
ggsave('cross.auc.CF79.pdf', width = 7, height = 6)

#### Fig. S6o ####
proj_name <- c(
  'BushmanFD_2020','FranzosaEA_2018','HallAB_2017','HeQ_2017','KumbhariA_2024',
  'LloydPriceJ_2019','SchirmerM_2018','SchirmerM_2024','WengY_2019','YanQ_2023c'
)

tpm <- readRDS('../Figure4/CF.tpm.rds')

CF46 <- read.xlsx('meta-analysis.xlsx') |> 
  filter(padj < 0.05) |> 
  pull(feature) |> 
  unique()

model <- map_dfr(
  proj_name, \(x) {
    
    group_x <- read.delim(
      paste0('../data/', x, '.sample_group.bz2')
    ) %>% 
      dplyr::select(sample, group) %>% 
      mutate(group = ifelse(group == 'HC', 'HC', 'IBD'))
    
    profile_x <- tpm[[x]][CF46, ]
    
    map_dfr(
      proj_name, \(y) {
        group_y <- read.delim(
          paste0('../data/', y, '.sample_group.bz2')
        ) %>% 
          dplyr::select(sample, group) %>% 
          mutate(group = ifelse(group == 'HC', 'HC', 'IBD'))
        
        profile_y <- tpm[[y]][CF46, ]
        pred <- rf_next_vaildate(
          profile_x, profile_y, group_x, group_y, seed = 2025
        )
        
        pred <- left_join(pred, group_y, by = 'sample')
        roc <- pROC::roc(pred$group, pred$HC, quiet = T) 
        data.frame(
          name_x = x, 
          name_y = y, 
          auc = as.numeric(roc$auc)
        )
      }
    ) 
  }, .progress = T
)

write_tsv(model, 'cross.auc.CF46.tsv')

model <- rbind(
  read.delim('cross.auc.CF46.tsv') %>% 
    filter(name_x != name_y),
    read.delim('CF.rf.roc.46_CF.tsv') %>%
    select(name_x = name, auc) %>% 
    mutate(name_y = name_x, .before = 'auc')
)

ggplot(model, aes(name_x, name_y)) +
  geom_tile(aes(fill = auc)) +
  labs(
    x = 'Discovery dataset', y = 'Validation dataset', fill = 'AUC',
    title = '46 CF-family cross-cohort validation') +
  geom_text(aes(label = round(auc, 2)), size = 2) +
  scale_fill_gradientn(
    colors = c('#4F6F9D',"#94A8D3","#FFFFFF","#FEEF6B",'#F7DF5A')
  ) +
  coord_fixed() +
  theme(
    axis.ticks = element_blank(),
    axis.text = element_text(color = '#000000', size = 10),
    axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1),
    plot.title = element_text(hjust = .5, face = 'bold'),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(linewidth = .8, color = 'black', fill = NA)
  )
ggsave('cross.auc.CF46.pdf', width = 7, height = 6)

#### Fig. S6p ####
CF13 <- read_rds('feature_importance.names.rds')[1:13]

model <- map_dfr(
  proj_name, \(x) {
    
    group_x <- read.delim(
      paste0('../data/', x, '.sample_group.bz2')
    ) %>% 
      dplyr::select(sample, group) %>% 
      mutate(group = ifelse(group == 'HC', 'HC', 'IBD'))
    
    profile_x <- tpm[[x]][CF13, ]
    
    map_dfr(
      proj_name, \(y) {
        group_y <- read.delim(
          paste0('../data/', y, '.sample_group.bz2')
        ) %>% 
          dplyr::select(sample, group) %>% 
          mutate(group = ifelse(group == 'HC', 'HC', 'IBD'))
        
        profile_y <- tpm[[y]][CF13, ]
        pred <- rf_next_vaildate(
          profile_x, profile_y, group_x, group_y, seed = 2025
        )
        
        pred <- left_join(pred, group_y, by = 'sample')
        roc <- pROC::roc(pred$group, pred$HC, quiet = T) 
        data.frame(
          name_x = x, 
          name_y = y, 
          auc = as.numeric(roc$auc)
        )
      }
    ) 
  }, .progress = T
)

write_tsv(model, 'cross.auc.CF13.tsv')

model <- rbind(
  read.delim('cross.auc.CF13.tsv') %>% 
    filter(name_x != name_y),
  read.delim('roc.top_13.tsv') %>%
    select(name_x = name, auc) %>% 
    mutate(name_y = name_x, .before = 'auc')
)

ggplot(model, aes(name_x, name_y)) +
  geom_tile(aes(fill = auc)) +
  labs(
    x = 'Discovery dataset', y = 'Validation dataset', fill = 'AUC',
    title = '13 CF-family cross-cohort validation') +
  geom_text(aes(label = round(auc, 2)), size = 2) +
  scale_fill_gradientn(
    colors = c('#4F6F9D',"#94A8D3","#FFFFFF","#FEEF6B",'#F7DF5A')
  ) +
  coord_fixed() +
  theme(
    axis.ticks = element_blank(),
    axis.text = element_text(color = '#000000', size = 10),
    axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1),
    plot.title = element_text(hjust = .5, face = 'bold'),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(linewidth = .8, color = 'black', fill = NA)
  )
ggsave('cross.auc.CF13.pdf', width = 7, height = 6)
