#### Jinxin Meng, 20251028, 20251109 ####
setwd('/data/mengjx/project/10.20250623_IBD_BAC_CF_Landscape/git/Figure4/')
pacman::p_load(tidyverse, ggpubr)
source('../scripts/palette.R')
source('../scripts/calcu_difference.R')
source('../scripts/plot_Procrustes.R')

#### Fig. 4a and Fig. S2 ####
proj_name <- c('BushmanFD_2020','FranzosaEA_2018','HallAB_2017','HeQ_2017','KumbhariA_2024',
               'LloydPriceJ_2019','SchirmerM_2018','SchirmerM_2024','WengY_2019','YanQ_2023c')

tpm <- readRDS('../Figure3/CF.tpm.rds')

plots <- map(proj_name, ~ {
  data_x <- vegan::decostand(t(apply(tpm[[.x]], 2, \(x) x / sum(x))), method = 'hellinger') %>% 
    t %>% data.frame()
  
  data_y <- data.table::fread(paste0('../data/', .x, '.kk2.s.bz2')) %>% 
    column_to_rownames('name') %>% 
    dplyr::select(any_of(colnames(data_x))) %>% 
    t() %>% vegan::decostand(method = 'hellinger') %>% 
    t %>% data.frame()
  
  plot_Procrustes(data_x, data_y, dis_method = 'bray', title = .x, show_grid = F, show_line = F,
                  colors = c(c('#fee722', '#5f7dbe'))) + 
    theme(aspect.ratio = 1,
          axis.ticks.length = unit(2, 'mm'),
          panel.grid.major = element_line(color = 'grey77', linewidth = .4, linetype = 'longdash')) } )

cowplot::plot_grid(plotlist = plots, nrow = 2, align = 'v')
ggsave('procrustes.boxplot.pdf', width = 25, height = 10)

data.frame(
  name = proj_name, 
  M2 = c(0.2524, 0.4284, 0.6289, 0.3296, 0.4890,
         0.3863, 0.1776, 0.2355, 0.4831, 0.4649),
  plab = '***') %>% 
  ggbarplot('name', 'M2', fill = 'name', palette = 'Spectral', legend = 'none', 
            x.text.angle = 60, xlab = '', ylab = 'Procrustes M2',
            label = '***', lab.col = 'red', lab.size = 5) +
  theme(aspect.ratio = 1/2,
        axis.ticks.length = unit(2, 'mm'),
        panel.grid.major = element_line(color = 'grey77', linewidth = .4, linetype = 'longdash'))
ggsave('procrustes.M2.barplot.pdf', width = 8, height = 5)

#### Fig. 4b ####
get_pc_by_cumsum_var <- function(x, cumsum = .95) {
  .cumsum <- 0
  for (i in 1:length(x)) {
    .cumsum = .cumsum + x[i]
    if (.cumsum > cumsum) 
      return(i)
  }
}

calcu_adjusted_r2 <- function(adonis_object) {
  n_observations <- adonis_object$Df[3]+1
  d_freedom <- adonis_object$Df[1]
  r2 <- adonis_object$R2[1]
  adjusted_r2 <- vegan::RsquareAdj(r2, n_observations, d_freedom)
  return(adjusted_r2)
}

test <- map(proj_name, ~ {
  data_x <- data.frame(apply(tpm[[.x]], 2, \(x) x / sum(x))) %>%
    t() %>% vegan::decostand(method = 'hellinger') %>% 
    data.frame()
  
  data_y <- data.table::fread(paste0('../data/', .x, '.kk2.s.bz2')) %>%
    column_to_rownames('name') %>%
    dplyr::select(any_of(rownames(data_x)))
  
  dist_y <- vegan::decostand(t(data_y), method = 'hellinger') %>%
    vegan::vegdist(method = 'bray')
  
  # PCA 选择一些成分
  pca_result <- pca(data_x)
  pca_summary <- summary(pca_result)
  pc_axis <- get_pc_by_cumsum_var(pca_summary$cont$importance[2,], cumsum = .95)
  variables <- scores(pca_result, display = 'sites', choices = 1:pc_axis)
  
  adonis <- vegan::adonis2(as.formula(paste0('dist_y ~ ', paste0(colnames(variables), collapse = ' + '))), 
                           data.frame(variables), permutations = 999, parallel = 80)
  adonis$r2adj <- c(calcu_adjusted_r2(adonis), NA, NA)
  adonis$pc_axis <- c(pc_axis, NA, NA)
  adonis$pc_var <- c(.95, NA, NA)
  adonis } ) %>% 
  set_names(proj_name)
write_rds(test, 'adonis.r2.BAC_by_CF.rds')

map2_dfr(test, proj_name, ~ data.frame(name = .y, r2adj = .x$r2adj[1], 
                                       pval = .x$`Pr(>F)`[1], axis = .x$pc_axis[1])) %>% 
  ggbarplot('name', 'r2adj', fill = 'name', palette = 'Spectral', legend = 'none', 
            x.text.angle = 60, xlab = '', ylab = 'Adonis r2adj',
            label = '***', lab.col = 'red', lab.size = 5) +
  theme(aspect.ratio = 1/2,
        axis.ticks.length = unit(2, 'mm'),
        panel.grid.major = element_line(color = 'grey77', linewidth = .4, linetype = 'longdash'))
ggsave('adonis.r2.BAC_by_CF.barplot.pdf', width = 8, height = 5)

#### Fig. 4c ####
get_pc_by_cumsum_var <- function(x, cumsum = .95) {
  .cumsum <- 0
  for (i in 1:length(x)) {
    .cumsum = .cumsum + x[i]
    if (.cumsum > cumsum) 
      return(i)
  }
}

calcu_adjusted_r2 <- function(adonis_object) {
  n_observations <- adonis_object$Df[3]+1
  d_freedom <- adonis_object$Df[1]
  r2 <- adonis_object$R2[1]
  adjusted_r2 <- vegan::RsquareAdj(r2, n_observations, d_freedom)
  return(adjusted_r2)
} 

test <- map(proj_name, ~ {
  data_x <- data.table::fread(paste0('../data/', .x, '.kk2.s.bz2')) %>%
    column_to_rownames('name') %>% 
    t() %>% vegan::decostand(method = 'hellinger') %>% 
    data.frame()
  
  data_y <- data.frame(apply(tpm[[.x]], 2, \(x) x / sum(x))) %>%
    dplyr::select(any_of(rownames(data_x)))
  
  dist_y <- decostand(t(data_y), method = 'hellinger') %>%
    vegan::vegdist(method = 'bray')
  
  # PCA 选择一些成分
  pca_result <- pca(data_x)
  pca_summary <- summary(pca_result)
  pc_axis <- get_pc_by_cumsum_var(pca_summary$cont$importance[2,], cumsum = .95)
  variables <- scores(pca_result, display = 'sites', choices = 1:pc_axis)
  
  adonis <- vegan::adonis2(as.formula(paste0('dist_y ~ ', paste0(colnames(variables), collapse = ' + '))), 
                           data.frame(variables), permutations = 999, parallel = 80)
  adonis$r2adj <- c(calcu_adjusted_r2(adonis), NA, NA)
  adonis$pc_axis <- c(pc_axis, NA, NA)
  adonis$pc_var <- c(.95, NA, NA)
  adonis } ) %>% 
  set_names(proj_name)
write_rds(test, 'adonis.r2.CF_by_BAC.rds')

map2_dfr(test, proj_name, ~ data.frame(name = .y, r2adj = .x$r2adj[1], 
                                       pval = .x$`Pr(>F)`[1], axis = .x$pc_axis[1])) %>% 
  ggbarplot('name', 'r2adj', fill = 'name', palette = 'Spectral', legend = 'none', 
            x.text.angle = 60, xlab = '', ylab = 'Adonis r2adj',
            label = '***', lab.col = 'red', lab.size = 5) +
  theme(aspect.ratio = 1/2,
        axis.ticks.length = unit(2, 'mm'),
        panel.grid.major = element_line(color = 'grey77', linewidth = .4, linetype = 'longdash'))
ggsave('adonis.r2.CF_by_BAC.barplot.pdf', width = 8, height = 5)

#### Fig. 4d ####
test <- map(proj_name, ~ {
  data_x <- data.frame(apply(tpm[[.x]], 2, \(x) x / sum(x))) %>%
    t() %>% vegan::decostand(method = 'hellinger') %>% 
    data.frame()
  
  data_y <- data.table::fread(paste0('../data/', .x, '.kk2.s.bz2')) %>%
    column_to_rownames('name') %>%
    dplyr::select(any_of(rownames(data_x)))
  
  dist_y <- vegan::decostand(t(data_y), method = 'hellinger') %>%
    vegan::vegdist(method = 'bray')
  
  map(colnames(data_x), \(x) {
    adonis <- vegan::adonis2(as.formula(paste0('dist_y ~ ', x )), data_x, permutations = 999, parallel = 80) 
    adonis$r2adj <- c(calcu_adjusted_r2(adonis), NA, NA)
    adonis }) %>%  
    set_names(colnames(data_x)) } ) %>% 
  set_names(proj_name)
write_rds(test, 'adonis.r2.single.rds')

data <- map2_dfr(test, names(test), \(x, y) 
                 map2_dfr(x, names(x), \(i, j)
                          data.frame(dataset = y, gene = j, r2adj = i$r2adj[1], pval = i$`Pr(>F)`[1]) ) ) %>% 
  mutate(plab = add_plab(pval, format = 4))
# openxlsx::write.xlsx(data, 'association/adonis.r2.single.xlsx')

plot_tile <- select(data, dataset, gene, r2adj) %>% 
  spread('gene', 'r2adj', fill = 0) %>% 
  column_to_rownames('dataset')

dataset_level <- rownames(plot_tile)[hclust(dist(plot_tile))$order]
gene_level <- colnames(plot_tile)[hclust(dist(t(plot_tile)))$order]

CF_info <- read.delim('../Figure1/pfam.rename.tsv')

rownames_to_column(plot_tile, 'dataset') %>% 
  gather('gene', 'value', -dataset) %>% 
  mutate(gene = factor(gene, gene_level),
         dataset = factor(dataset, dataset_level)) %>% 
  ggplot(aes(gene, dataset)) +
  geom_tile(aes(fill = value)) +
  geom_text(aes(gene, dataset, label = plab), data, inherit.aes = F, size = 3) +
  scale_x_discrete(
    breaks = gene_level, 
    labels = paste0(gene_level, ' (', CF_info$pfam[match(gene_level, CF_info$name)], ')') ) +
  scale_fill_gradientn(colors = rev(pald('Spectral')[-11])) +
  labs(x = '', y = '', fill = 'adjusted R2') +
  coord_fixed() +
  theme(axis.ticks = element_blank(),
        axis.text = element_text(color = '#000000', size = 10),
        axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = .8, color = 'black', fill = 'transparent'))
ggsave('adonis.r2.single.pdf', width = 16, height = 6)
