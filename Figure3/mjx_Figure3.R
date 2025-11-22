#### Jinxin Meng, 20251028, 20251122 ####
setwd('/data/mengjx/project/10.20250623_IBD_BAC_CF_Landscape/git/Figure3/')
pacman::p_load(tidyverse, ggpubr)
source('../scripts/palette.R')
source('../scripts/calcu_difference.R')
source('../scripts/transform_rc.R')

#### profile preprocess ####
proj_name <- c('BushmanFD_2020','FranzosaEA_2018','HallAB_2017','HeQ_2017','KumbhariA_2024',
               'LloydPriceJ_2019','SchirmerM_2018','SchirmerM_2024','WengY_2019','YanQ_2023c')

gene_len <- read.delim('../pipeline/uhgp.len.bz2', header = F, col.names = c('name', 'length'))
gene_info <- read.delim('../pipeline/uhgp.m8.f.drop.info.bz2', header = F) %>%
  dplyr::select(name = 1, value = 2) %>%
  mutate(CF = stringr::str_split_i(value, ',', 1),
         CF_gene = stringr::str_split_i(value, ',', 2))

# gene tpm
tpm <- map(proj_name, ~ {
  rc <- data.table::fread(paste0('../data/', .x, '.rc.bz2')) %>% column_to_rownames('name')
  cvg <- data.table::fread(paste0('../data/', .x, '.cvg.bz2')) %>% column_to_rownames('name')
  tpm <- rc2tpm(rc, gene_len) * (cvg > 50)
  tpm[rowSums(tpm) != 0, colSums(tpm) != 0] } ) %>%
  set_names(proj_name)
saveRDS(tpm, 'gene.tpm.rds')

# gene rc
rc <- map(proj_name, ~ {
  cvg <- data.table::fread(paste0('../data/', .x, '.cvg.bz2')) %>% column_to_rownames('name')
  rc <- data.table::fread(paste0('../data/', .x, '.rc.bz2')) %>% column_to_rownames('name') * (cvg > 50)
  rc[rowSums(rc) != 0, colSums(rc) != 0] } ) %>%
  set_names(proj_name)
saveRDS(rc, 'gene.rc.rds')

# CF tpm
tpm <- map(proj_name, ~ {
  rc <- data.table::fread(paste0('../data/', .x, '.rc.bz2')) %>% column_to_rownames('name')
  tpm <- rc2tpm(rc, gene_len)
  tpm[colSums(tpm) != 0, ] %>%
    mutate(name = gene_info$CF[match(rownames(.), gene_info$name)]) %>%
    aggregate(. ~ name, ., sum) %>%
    column_to_rownames('name') %>%
    filter(rowSums(.) != 0) } ) %>%
  set_names(proj_name)
saveRDS(tpm, 'CF.tpm.rds')

# CF rc
rc <- map(proj_name, ~ {
  rc <- data.table::fread(paste0('../data/', .x, '.rc.bz2')) %>% column_to_rownames('name')
  mutate(rc, name = gene_info$CF[match(rownames(rc), gene_info$name)]) %>%
    aggregate(. ~ name, ., sum) %>%
    column_to_rownames('name') %>%
    filter(rowSums(.) != 0) } ) %>%
  set_names(proj_name)
saveRDS(rc, 'CF.rc.rds')

group <- map(proj_name, ~ 
               data.table::fread(paste0('../data/', .x, '.sample_group.bz2')) %>% 
               select(sample, group)) %>% 
  set_names(proj_name)
saveRDS(group, 'group.rds')

#### Fig. 3a ####
data <- map_dfr(proj_name, ~ {
  group <- data.table::fread(paste0('../data/', .x, '.sample_group.bz2')) %>% 
    select(sample, group)
  stat <- count(group, group)
  stat <- filter(stat, group != 'HC') %>% 
    rename(case = n) %>% 
    add_column(control = unlist(stat[stat$group == 'HC', 'n'])) %>% 
    add_column(proj = .x) } ) %>% 
  mutate(group = factor(group, c('CD', 'UC', 'IC', 'IBD')),
         name = paste0(proj, '|', group)) %>% 
  arrange(group) %>% 
  mutate(name = factor(name, name)) %>% 
  select(-group, -proj) %>% 
  gather('group', 'value', -name)

name_level <- levels(data$name)

ggbarplot(data, 'name', 'value', rotate = T, fill = 'group', palette = c('#fee722', '#5f7dbe'),
          ylab = 'Number of samples', xlab = '', legend = 'right') +
  scale_y_continuous(expand = c(0, 0)) +
  theme(aspect.ratio = 3,
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.length = unit(2, 'mm'), 
        panel.grid.major = element_line(linewidth = .5, color = 'grey88'),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = .5, color = 'black', fill = 'transparent'))
ggsave('sample.stat.pdf', width = 8, height = 10)

#### Fig. 3b ####
tpm <- read_rds('CF.tpm.rds')

data <- map(proj_name, ~ {
  group <- read.delim(paste0('../data/', .x, '.sample_group.bz2')) %>% 
    dplyr::select(sample, group) %>% 
    filter(sample %in% colnames(tpm[[.x]]))
  group <- group[match(colnames(tpm[[.x]]), group$sample),]
  group_level <- intersect(c('HC','IBD','CD','UC','IC'), unique(group$group))
  
  map(setdiff(group_level, 'HC'), \(y) {
    .group <- filter(group, group %in% c(y, 'HC'))
    .tpm <- tpm[[.x]][, .group$sample]
    .name = paste0(.x, '.', y)
    data.frame(value = rowSums(.tpm)) %>% 
      rownames_to_column('name') %>%
      rename(!!.name := value) } ) %>% 
    reduce(~ inner_join(.x, .y, by = 'name'))  } ) %>% 
  reduce(~ inner_join(.x, .y, by = 'name'))

plot_data <- column_to_rownames(data, 'name') %>% 
  apply(2, \(x) x / sum(x) * 100) %>% 
  data.frame()

.names <- rowMeans(plot_data) %>% sort %>% tail(n = 11) %>% names()

rownames_to_column(plot_data, 'name') %>% 
  mutate(name = ifelse(name %in% .names, name, 'Other CF')) %>%
  aggregate(. ~ name, ., sum) %>% 
  filter(name != 'Other CF') %>% 
  gather('group', 'value', -name) %>% 
  mutate(group = sub('\\.', '|', group),
         group = factor(group, name_level)) %>% 
  ggbarplot('group', 'value', rotate = T, fill = 'name', ylab = 'Relative abundance (%)',
            xlab = '', legend = 'right') +
  scale_fill_manual(values = c('#fb8072','#80b1d3','#ffffb3','#fccde5','#ffed6f','#fdb462',
                               '#b3de69','#8dd3c7','#bebada','#bc80bd','#ccebc5','grey77')) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(aspect.ratio = 3,
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.length = unit(2, 'mm'), 
        panel.grid.major = element_line(linewidth = .5, color = 'grey88'),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = .5, color = 'black', fill = 'transparent'))
ggsave('composition.pdf', width = 8, height = 10)

#### Fig. 3c ####
tpm <- read_rds('gene.tpm.rds')

data <- map_dfr(proj_name, ~ {
  group <- read.delim(paste0('../data/', .x, '.sample_group.bz2')) %>% 
    dplyr::select(sample, group)
  group_level <- intersect(c('HC','IBD','CD','UC','IC'), unique(group$group))
  
  data <- data.frame(value = vegan::diversity(t(tpm[[.x]]), 'shannon')) %>% 
    rownames_to_column('sample') %>% 
    left_join(group, by = 'sample') %>% 
    mutate(group = factor(group, group_level))
  
  test <- calcu_diff(data, value ~ group) %>% 
    add_column(proj = .x) %>% 
    filter(grepl('HC', comparison)) %>% 
    mutate(plab = add_plab(pval),
           name = paste0(proj, '|', sub('HC_vs_', '', comparison))) %>% 
    select(name, plab)
  
  stat <- aggregate(value ~ group, data, mean)
  stat <- filter(stat, group != 'HC') %>% 
    rename(case = value) %>% 
    add_column(control = unlist(stat[stat$group == 'HC', 'value'])) %>% 
    add_column(proj = .x) %>% 
    mutate(group = factor(group, c('CD', 'UC', 'IC', 'IBD')),
           name = paste0(proj, '|', group),
           log2FC = log2(case / control)) %>% 
    select(name, log2FC)
  
  left_join(stat, test, by = 'name') } )

mutate(data,
       name = factor(name, name_level),
       group = str_split_i(name, '\\|', 2)) %>%
  ggbarplot('name', 'log2FC', rotate = T, fill = 'group',
            palette = c(UC = '#80b1d3', CD = '#b3de69', IC = '#fdb462', IBD = '#8dd3c7'),
            ylab = 'Coefficient (Shannon)', xlab = '', legend = 'none') +
  geom_text(aes(label = plab), color = 'red', vjust = .9) +
  theme(aspect.ratio = 3,
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.length = unit(2, 'mm'), 
        panel.grid.major = element_line(linewidth = .5, color = 'grey88'),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = .5, color = 'black', fill = 'transparent'))
ggsave('shannon.pdf', width = 6, height = 10)

#### Fig. 3d ####
tpm <- read_rds('gene.tpm.rds')

data <- map_dfr(proj_name, ~ {
  group <- read.delim(paste0('../data/', .x, '.sample_group.bz2')) %>% 
    dplyr::select(sample, group)
  group_level <- intersect(c('HC','IBD','CD','UC','IC'), unique(group$group))
  
  data <- data.frame(value = colSums(tpm[[.x]] > 0)) %>% 
    rownames_to_column('sample') %>% 
    left_join(group, by = 'sample') %>% 
    mutate(group = factor(group, group_level))
  
  test <- calcu_diff(data, value ~ group) %>% 
    add_column(proj = .x) %>% 
    filter(grepl('HC', comparison)) %>% 
    mutate(plab = add_plab(pval),
           name = paste0(proj, '|', sub('HC_vs_', '', comparison))) %>% 
    select(name, plab)
  
  stat <- aggregate(value ~ group, data, mean)
  stat <- filter(stat, group != 'HC') %>% 
    rename(case = value) %>% 
    add_column(control = unlist(stat[stat$group == 'HC', 'value'])) %>% 
    add_column(proj = .x) %>% 
    mutate(group = factor(group, c('CD', 'UC', 'IC', 'IBD')),
           name = paste0(proj, '|', group),
           log2FC = log2(case / control)) %>% 
    select(name, log2FC)
  
  left_join(stat, test, by = 'name') } )

mutate(data,
       name = factor(name, name_level),
       group = str_split_i(name, '\\|', 2)) %>%
  ggbarplot('name', 'log2FC', rotate = T, fill = 'group',
            palette = c(UC = '#80b1d3', CD = '#b3de69', IC = '#fdb462', IBD = '#8dd3c7'),
            ylab = 'Coefficient (# CF gene homologs)', xlab = '', legend = 'none') +
  geom_text(aes(label = plab), color = 'red', vjust = .9) +
  theme(aspect.ratio = 3,
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.length = unit(2, 'mm'), 
        panel.grid.major = element_line(linewidth = .5, color = 'grey88'),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = .5, color = 'black', fill = 'transparent'))
ggsave('richness.pdf', width = 6, height = 10)

#### Fig. 3f ####
tpm <- read_rds('gene.tpm.rds')

calcu_adjusted_r2 <- function(adonis_object) {
  n_observations <- adonis_object$Df[3]+1
  d_freedom <- adonis_object$Df[1]
  r2 <- adonis_object$R2[1]
  adjusted_r2 <- vegan::RsquareAdj(r2, n_observations, d_freedom)
  return(adjusted_r2)
}

data <- map_dfr(proj_name, ~ {
  group <- read.delim(paste0('../data/', .x, '.sample_group.bz2')) %>% 
    dplyr::select(sample, group) %>% 
    filter(sample %in% colnames(tpm[[.x]]))
  group <- group[match(colnames(tpm[[.x]]), group$sample),]
  group_level <- intersect(c('HC','IBD','CD','UC','IC'), unique(group$group))
  
  map_dfr(setdiff(group_level, 'HC'), \(y) {
    .group <- filter(group, group %in% c(y, 'HC'))
    .tpm <- tpm[[.x]][, .group$sample]
    adonis <- vegan::adonis2(t(.tpm) ~ group, .group, permutations = 999, 
                             distance = 'bray', parallel = 110)
    r2adj <- calcu_adjusted_r2(adonis)
    data.frame(name = paste0(.x, '|', y), r2 = adonis$R2[1], r2adj = r2adj, pval = adonis$`Pr(>F)`[1])
  } )
} )
openxlsx::write.xlsx(data, 'adonis.xlsx')

mutate(data,
       name = factor(name, name_level),
       plab = add_plab(pval),
       group = str_split_i(name, '\\|', 2)) %>%
  ggbarplot('name', 'r2adj', rotate = T, fill = 'group',
            palette = c(UC = '#80b1d3', CD = '#b3de69', IC = '#fdb462', IBD = '#8dd3c7'),
            ylab = 'Adjust R2 (adonis)', xlab = '', legend = 'none') +
  geom_text(aes(label = plab), color = 'red', vjust = .9) +
  theme(aspect.ratio = 3,
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.length = unit(2, 'mm'), 
        panel.grid.major = element_line(linewidth = .5, color = 'grey88'),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = .5, color = 'black', fill = 'transparent'))
ggsave('adonis.pdf', width = 6, height = 10)
