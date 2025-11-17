#### Jinxin Meng, 20251028, 20251109 ####
setwd('/data/mengjx/project/10.20250623_IBD_BAC_CF_Landscape/git/Figure7/')
pacman::p_load(tidyverse, ggpubr)
source('../scripts/palette.R')
source('../scripts/calcu_difference.R')
source('../scripts/transform_rc.R')

#### profile preprocess ####
proj_name <- c('LloydPriceJ_2019','SchirmerM_2018')

filter_samples <- list(
  LloydPriceJ_2019.exp = read.delim('../data/LloydPriceJ_2019.exp.lib_size.bz2') %>% 
    filter(reads_count < 1e6) %>% 
    pull(sample),
  SchirmerM_2018.exp = read.delim('../data/SchirmerM_2018.exp.lib_size.bz2') %>% 
    filter(reads_count < 1e6) %>% 
    pull(sample)
)

gene_len <- read.delim('../pipeline/uhgp.len.bz2', header = F, col.names = c('name', 'length'))
gene_info <- read.delim('../pipeline/uhgp.m8.info.bz2', header = F) %>%
  dplyr::select(name = 1, value = 2) %>%
  mutate(CF = stringr::str_split_i(value, ',', 1),
         CF_gene = stringr::str_split_i(value, ',', 2))

# gene tpm
tpm <- map2(proj_name, filter_samples, ~ {
  rc <- data.table::fread(paste0('../data/', .x, '.exp.rc.bz2')) %>%
    column_to_rownames('name') %>% 
    dplyr::select(!all_of(.y))
  cvg <- data.table::fread(paste0('../data/', .x, '.exp.cvg.bz2')) %>%
    column_to_rownames('name') %>% 
    dplyr::select(!all_of(.y))
  tpm <- rc2tpm(rc, gene_len) * (cvg > 5)
  tpm[rowSums(tpm) != 0, colSums(tpm) != 0] } ) %>%
  set_names(proj_name)
saveRDS(tpm, 'gene.exp.tpm.rds')

# CF tpm
drops <- c('CF40', 'CF58', 'CF62')
tpm <- map2(proj_name, filter_samples, ~ {
  rc <- data.table::fread(paste0('../data/', .x, '.exp.rc.bz2')) %>% 
    column_to_rownames('name') %>% 
    dplyr::select(!all_of(.y))
  tpm <- rc2tpm(rc, gene_len)
  tpm[colSums(tpm) != 0, ] %>%
    mutate(name = gene_info$CF[match(rownames(.), gene_info$name)]) %>%
    aggregate(. ~ name, ., sum) %>%
    column_to_rownames('name') %>%
    filter(rowSums(.) != 0) %>% 
    filter(!rownames(.) %in% drops) } ) %>%
  set_names(proj_name)
saveRDS(tpm, 'CF.exp.tpm.rds')

# CF rc
rc <- map2(proj_name, filter_samples, ~ {
  rc <- data.table::fread(paste0('../data/', .x, '.exp.rc.bz2')) %>% 
    column_to_rownames('name') %>% 
    dplyr::select(!all_of(.y))
  mutate(rc, name = gene_info$CF[match(rownames(rc), gene_info$name)]) %>%
    aggregate(. ~ name, ., sum) %>%
    column_to_rownames('name') %>%
    filter(rowSums(.) != 0) %>% 
    filter(!rownames(.) %in% drops) } ) %>%
  set_names(proj_name)
saveRDS(rc, 'CF.exp.rc.rds')

group <- map2(proj_name, filter_samples, ~ 
               data.table::fread(paste0('../data/', .x, '.exp.sample_group.bz2')) %>% 
               filter(!sample %in% .y)) %>% 
  set_names(proj_name)
saveRDS(group, 'group.exp.rds')

#### Fig. 7a-b ####
plots <- map(
  proj_name, ~ {
    .group <- group[[.x]]
    MTX_profile <- read_rds('CF.exp.tpm.rds')[[.x]][.group$sample]
    name_level <- names(sort(rowMeans(MTX_profile), decreasing = T)[1:20])
    
    data <- MTX_profile[name_level,] %>% 
      rownames_to_column('name') %>% 
      gather('sample', 'value', -name) %>% 
      mutate(name = factor(name, rev(name_level)),
             .value = log10(value)) %>% 
      filter(!is.infinite(.value))
    
    ggplot(data, aes(.value, name, fill = name)) +
      ggridges::geom_density_ridges() +
      scale_fill_manual(values = colorRampPalette(pald('Spectral'))(20)) +
      scale_y_discrete(expand = c(0.03, 0)) +
      labs(x = 'log10 TPM', y = '', title = paste0(.x, '\nTop20 expression CF')) +
      theme_bw() +
      theme(aspect.ratio = 2,
            axis.ticks.length = unit(2, 'mm'),
            axis.text = element_text(color = 'black', size = 10),
            plot.title = element_text(color = 'black', hjust = .5, face = 'bold'),
            panel.grid.major.x = element_line(linewidth = .5, color = 'grey88'),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(linewidth = .8, color = 'black')) +
      guides(fill = 'none')
  } )
cowplot::plot_grid(plotlist = plots, nrow = 1)
ggsave('CF.exp.ggridge.pdf', width = 8, height = 8)

#### Fig. 7c-d ####
group <- read_rds('group.exp.rds')

plots <- map(
  proj_name, ~ {
    .group <- filter(group[[.x]], !is.na(sample_MGX))
    MTX_profile <- log10(read_rds('CF.exp.tpm.rds')[[.x]][.group$sample] + 1)
    MGX_profile <- log10(read_rds('../Figure3/CF.tpm.rds')[[.x]][rownames(MTX_profile), .group$sample_MGX] + 1)
    
    data <- data.frame(
      name = rownames(MTX_profile),
      MTX_mean = rowMeans(MTX_profile),
      # MTX_sd = apply(MTX_profile, 1, sd),
      MGX_mean = rowMeans(MGX_profile),
      # MGX_sd = apply(MGX_profile, 1, sd),
      row.names = NULL
    )
    
    ggplot(data, aes(MTX_mean, MGX_mean)) +
      # geom_errorbar(aes(ymin = MGX_mean - MGX_sd, ymax = MGX_mean + MGX_sd), width = 0) +
      # geom_errorbar(aes(xmin = MTX_mean - MTX_sd, xmax = MTX_mean + MTX_sd), width = 0) +
      geom_point(fill = 'white', shape = 21, size = 5, show.legend = F) +
      ggpmisc::stat_poly_eq(ggpmisc::use_label('eq', 'R2', 'P'), parse = T, 
                            color = 'blue') +
      geom_text(aes(label = name)) +
      geom_smooth(method = 'lm', se = F)  + 
      geom_abline(slope = 1) +
      labs(x = 'log10 average expression (TPM)', y = 'log10 average abundance (TPM)',
           title = .x) +
      theme_gray() +
      theme(aspect.ratio = 1,
            axis.ticks.length = unit(2, 'mm'),
            plot.title = element_text(color = 'black', hjust = .5, face = 'bold'))
  } )
cowplot::plot_grid(plotlist = plots, nrow = 1)
ggsave('association.pdf', width = 10, height = 5)

#### Fig. 7e-f ####
plots <- map(proj_name, ~ {
  .group <- dplyr::select(group[[.x]], sample, group)
  group_level <- intersect(c('HC','IBD','CD','UC','IC'), unique(.group$group))
  data <- read.delim(paste0('../data/', .x, '.exp.map_rate.bz2'), header = F, 
                     col.names = c('sample', 'value')) %>% 
    filter(sample %in% .group$sample) %>% 
    mutate(group = .group$group[match(sample, .group$sample)], 
           group = factor(group, group_level),
           value = sub('%', '', value),
           value = as.numeric(value))
  ggboxplot(data, 'group', 'value', color = 'group', palette = 'npg', legend = 'none', 
            xlab = '', ylab = 'Mapping rate (%)', title = .x, width = .7, outlier.shape = NA) + 
    geom_jitter(aes(color = group), width = .25) +
    geom_signif(comparisons = list(c('CD', 'HC'), c('UC', 'HC'), c('CD', 'UC')), 
                step_increase = .08) +
    theme(axis.ticks.length = unit(2, 'mm'),
          plot.title = element_text(color = 'black', hjust = .5, face = 'bold'),
          aspect.ratio = 1) 
} )
cowplot::plot_grid(plotlist = plots, nrow = 1, align = 'v')
ggsave('mapping.boxplot.pdf', width = 6, height = 3)

#### Fig. 7g-h ####
tpm <- read_rds('gene.exp.tpm.rds')

plots <- map(proj_name, ~ {
  .tpm <- tpm[[.x]]
  .group <- group[[.x]]
  group_level <- intersect(c('HC','IBD','CD','UC','IC'), unique(.group$group))
  data <- data.frame(value = vegan::diversity(t(tpm[[.x]]), 'shannon')) %>% 
    rownames_to_column('sample') %>% 
    left_join(.group, by = 'sample') %>% 
    mutate(group = factor(group, group_level))
  
  ggboxplot(data, 'group', 'value', color = 'group', palette = 'npg', legend = 'none', 
            xlab = '', ylab = 'Shannon Index', title = .x, width = .7, outlier.shape = NA) + 
    geom_jitter(aes(color = group), width = .25) +
    geom_signif(comparisons = list(c('CD', 'HC'), c('UC', 'HC'), c('CD', 'UC')), step_increase = .08) +
    theme(axis.ticks.length = unit(2, 'mm'),
          plot.title = element_text(color = 'black', hjust = .5, face = 'bold'),
          aspect.ratio = 1) 
})
cowplot::plot_grid(plotlist = plots, nrow = 1, align = 'v')
ggsave('div.shannon.boxplot.pdf', width = 6, height = 3)

#### Fig. 7i-j ####
tpm <- read_rds('gene.exp.tpm.rds')

plots <- map(proj_name, ~ {
  .tpm <- tpm[[.x]]
  .group <- group[[.x]]
  group_level <- intersect(c('HC','IBD','CD','UC','IC'), unique(.group$group))
  data <- data.frame(value = colSums(tpm[[.x]] > 0)) %>% 
    rownames_to_column('sample') %>% 
    left_join(.group, by = 'sample') %>% 
    mutate(group = factor(group, group_level))
  
  ggboxplot(data, 'group', 'value', color = 'group', palette = 'npg', legend = 'none', 
            xlab = '', ylab = 'Richness Index', title = .x, width = .7, outlier.shape = NA) + 
    geom_jitter(aes(color = group), width = .25) +
    geom_signif(comparisons = list(c('CD', 'HC'), c('UC', 'HC'), c('CD', 'UC')), step_increase = .08) +
    theme(axis.ticks.length = unit(2, 'mm'),
          plot.title = element_text(color = 'black', hjust = .5, face = 'bold'),
          aspect.ratio = 1) 
})
cowplot::plot_grid(plotlist = plots, nrow = 1, align = 'v')
ggsave('div.richness.boxplot.pdf', width = 6, height = 3)

#### Fig. 7k ####
source('/data/mengjx/R_proj/R_func/model_randomforest.R')
source('/data/mengjx/R_proj/R_func/plot_roc.R')

tpm <- readRDS('CF.exp.tpm.rds')
group <- readRDS('group.exp.rds')

roc <- map(
  proj_name, ~ {
    .group <- group[[.x]] %>% 
      dplyr::select(sample, group) %>% 
      mutate(group = ifelse(group == 'HC', 'HC', 'IBD'))
    pred <- rf_Kfold(tpm[[.x]], .group, k = 10, seed = 2025)
    pred <- left_join(pred, .group, by = 'sample')
    roc(pred$group, pred$HC) 
  } ) %>% 
  set_names(proj_name)

plot_roc_multiple(roc, title = 'Random Forest Model')
ggsave('CF.exp.rf.roc.all_CF.pdf', width = 8, height = 4)

#### Fig. 7l ####
data <- map(
  proj_name, ~ {
    .group <- group[[.x]] %>% 
      dplyr::select(sample, group) %>% 
      mutate(group = ifelse(group == 'HC', 'HC', 'IBD'))
    rf_importance(tpm[[.x]], .group, seed = 2025)
  } ) %>% 
  set_names(proj_name)

gene_info <- read.delim('../Figure1/pfam.rename.tsv')
data$LloydPriceJ_2019 %>% 
  head(30) %>%
  mutate(name = paste0(name, ' (', gene_info$pfam[match(name, gene_info$name)], ')')) %>% 
  ggbarplot('name', 'MeanDecreaseAccuracy', fill = 'name', legend = 'none',
            x.text.angle = 90, xlab = '') +
  scale_fill_manual(values = colorRampPalette(pald('Spectral'))(30)) +
  theme(axis.ticks.length = unit(2, 'mm'),
        plot.title = element_text(color = 'black', hjust = .5, face = 'bold'),
        aspect.ratio = 1/3) 
ggsave('CF.exp.rf.importance.pdf', width = 10, height = 5)
