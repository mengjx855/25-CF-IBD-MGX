#### Jinxin Meng, 20251028, 20251122 ####
setwd('/data/mengjx/project/10.20250623_IBD_BAC_CF_Landscape/git/Figure2/')
pacman::p_load(tidyverse, ggpubr)
source('../scripts/palette.R')
source('../scripts/calcu_difference.R')

#### Fig. 2a ####
genome_info <- data.table::fread('../pipeline/genomes-all_metadata.tsv.bz2') %>% 
  mutate(genome = sub('\\.\\d+$', '', Genome),
         domain = str_split_i(Lineage, pattern = ';', 1),
         phylum = str_split_i(Lineage, pattern = ';', 2),
         class = str_split_i(Lineage, pattern = ';', 3),
         order = str_split_i(Lineage, pattern = ';', 4),
         family = str_split_i(Lineage, pattern = ';', 5),
         genus = str_split_i(Lineage, pattern = ';', 6),
         species = str_split_i(Lineage, pattern = ';', 7) )

counts <- data.table::fread('../pipeline/uhgp.m8.f.drop.spread.species.binary.bz2') %>% 
  column_to_rownames('name')

dist <- vegan::vegdist(counts, method = 'jaccard', binary = T)
pcoa <- ape::pcoa(dist)
pcoa_rela_eig <- pcoa$values$Relative_eig[1:2]
# [1] 0.1395993 0.1092038

point <- data.frame(pcoa$vectors[, 1:2]) %>% 
  dplyr::rename_with(~ c('X1', 'X2')) %>% 
  add_column(genome = rownames(counts), .before = 1) %>% 
  mutate(phylum = genome_info$phylum[match(genome, genome_info$genome)],
         phylum = fct_lump_n(phylum, n = 11, ties.method = 'first', other_level = 'Other phyla'))

ggscatter(point, 'X1', 'X2', color = 'phylum', legend = 'right', size = 1.5,
          xlab = 'PCoA1 (14.0%)', ylab = 'PCoA2 (10.9%)') +
  scale_color_manual(values = c('#fb8072','#80b1d3','#ffffb3','#fccde5','#ffed6f','#fdb462',
                                '#b3de69','#8dd3c7','#bebada','#bc80bd','#ccebc5','grey77')) +
  theme(aspect.ratio = 3/4,
        axis.line = element_blank(),
        axis.ticks.length = unit(2, 'mm'), 
        panel.grid.major = element_line(linewidth = .8, color = 'grey88'),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = .8, color = 'black', fill = 'transparent'))
ggsave('pcoa.jaccard.scatter_plot.pdf', width = 9, height = 6)

#### Fig. 2b ####
calcu_adjusted_r2 <- function(adonis_object) {
  n_observations <- adonis_object$Df[3]+1
  d_freedom <- adonis_object$Df[1]
  r2 <- adonis_object$R2[1]
  adjusted_r2 <- vegan::RsquareAdj(r2, n_observations, d_freedom)
  return(adjusted_r2)
}

metadata <- tibble(genome = rownames(counts)) %>%
  left_join(dplyr::select(genome_info, genome, phylum, class, order, family, genus, species) %>%
              distinct(), by = 'genome')

test <- map(c('phylum','class','order','family','genus'),
            ~ vegan::adonis2(as.formula(paste0('dist ~ ', .x)),
                             metadata, permutations = 999, parallel = 110) ) %>%
  set_names(c('phylum','class','order','family','genus'))

map2_vec(names(test), test,
         ~ tibble(name = str_to_title(.x), r2adj = calcu_adjusted_r2(.y), pval = .y[1, 5])) %>% 
  mutate(plab = add_plab(pval)) %>% 
  ggbarplot('name', 'r2adj', fill = 'name', xlab = '', ylab = 'Adjusted r2', 
            palette = 'Spectral', width = .6) +
  geom_text(aes(label = plab), color = 'red', hjust = .5) +
  theme_minimal() +
  theme(aspect.ratio = 3/5, 
        axis.line = element_line(linewidth = .5),
        axis.ticks = element_line(linewidth = .5),
        axis.text = element_text(color = '#000000', size = 10),
        panel.grid.major = element_line(linewidth = .5),
        panel.grid.minor = element_blank(),
        axis.ticks.length = unit(2, 'mm')) +
  guides(fill = 'none')
ggsave('adonis.test.jaccard.pdf', width = 4, height = 5)

#### Fig. 2c and Fig. S1 ####
dist <- read_rds('dist.jaccard.rds')

taxa <- data.frame(ref_species = labels(dist)) %>% 
  mutate(domain = genome_info$domain[match(ref_species, genome_info$genome)])

tr <- ape::read.tree('bac120_iqtree.nwk')
taxa_name <- filter(taxa, domain == 'd__Bacteria') %>% pull(ref_species)
taxa_name <- taxa_name[taxa_name %in% tr$tip.label]

patristic_dist <- data.frame(ape::cophenetic.phylo(tr)[taxa_name, taxa_name]) # Patristic distance
jaccard_dist <- data.frame(as.matrix(dist)[taxa_name, taxa_name])

set.seed(2026)
data.frame(x = as.vector(as.dist(patristic_dist)),
           y = as.vector(as.dist(jaccard_dist))) %>% 
  slice_sample(prop = .0001) %>% 
  ggscatter('x', 'y', color = '#3288bd', xlab = 'Patristic distance (phylogenetic)', 
            ylab = 'Jaccard distance (CF-based)', title = 'Bacterial lineage (Species=4,691)', 
            add = 'reg.line', add.params = list(color = '#000000', linewidth = .8), 
            cor.coef = T, cor.method = 'spearman') +
  theme(aspect.ratio = 1,
        axis.ticks.length = unit(2, 'mm'),
        plot.title = element_text(color = 'black', hjust = .5, face = 'bold'))
ggsave('bacterial.dist.cor.seed2026.pdf', width = 4, height = 4.5)

plots <- map2(
  2027:2046, colorRampPalette(pald('Spectral'))(20), ~ {
    set.seed(.x)
    data.frame(x = as.vector(as.dist(patristic_dist)),
               y = as.vector(as.dist(jaccard_dist))) %>% 
      slice_sample(prop = .0001) %>% 
      ggscatter('x', 'y', color = .y, xlab = 'Patristic distance (phylogenetic)', 
                ylab = 'Jaccard distance (CF-based)', 
                title = paste0('stochastic sample ', .x - 2026, ' time'),
                subtitle = paste0('(proportion: 0.0001, seed: ', .x, ')'),
                add = 'reg.line', add.params = list(color = '#000000', linewidth = .8), 
                cor.coef = T, cor.method = 'spearman') +
      theme(aspect.ratio = 1,
            axis.ticks.length = unit(2, 'mm'),
            plot.title = element_text(color = 'black', hjust = .5, face = 'bold'),
            plot.subtitle = element_text(color = 'black', hjust = .5, face = 'bold'))
  } )
cowplot::plot_grid(plotlist = plots, nrow = 4)
ggsave('bacterial.dist.cor.ramdom20.pdf', width = 20, height = 17)

c(0.42, 0.46, 0.42, 0.41, 0.34, 0.4, 0.34, 0.4, 0.44, 0.45,
  0.39, 0.41, 0.43, 0.4, 0.44, 0.46, 0.44, 0.41, 0.41, 0.4, 0.45) %>% 
  summary()

#### Fig. 2d ####
tr <- ape::read.tree('/data/database/uhgg/phylogenies/ar122_iqtree.nwk')
taxa_name <- filter(taxa, domain == 'd__Archaea') %>% pull(ref_species)
taxa_name <- taxa_name[taxa_name %in% tr$tip.label]

patristic_dist <- data.frame(ape::cophenetic.phylo(tr)[taxa_name, taxa_name])
jaccard_dist <- data.frame(as.matrix(dist)[taxa_name, taxa_name])

data.frame(x = as.vector(as.dist(patristic_dist)),
           y = as.vector(as.dist(jaccard_dist))) %>% 
  ggscatter('x', 'y', xlab = 'Patristic distance (phylogenetic)', 
            ylab = 'Jaccard distance (CF-based)', color = '#fdae61',
            title = 'Archaeal lineage (Species=25)', 
            add = 'reg.line', add.params = list(color = '#000000', linewidth = .8), 
            cor.coef = T, cor.method = 'spearman') +
  theme(aspect.ratio = 1,
        axis.ticks.length = unit(2, 'mm'),
        plot.title = element_text(color = 'black', hjust = .5, face = 'bold'))
ggsave('archaeal.dist.cor.seed2026.pdf', width = 4, height = 4.5)

#### Fig. 2e ####
PAM_clustering <- function(dist, k) {
  cluster <- as.vector(cluster::pam(as.dist(dist), k, diss = T)$clustering)
  return(cluster)
}

dist <- read_rds('dist.jaccard.rds')

ncluster <- map_vec(1:10, ~ 
                      clusterSim::index.G1(x = (counts > 0) * 1, cl = PAM_clustering(dist, .x),
                                           d = dist, centrotypes = 'medoids'), .progress = T)

data.frame(k_cluster = seq(2, 10), CH = ncluster[-1]) %>% 
  ggline('k_cluster', 'CH', xlab = 'k clusters', ylab = 'Calinski-Harabasz index',
         title = 'Optimal number of clusters', shape = 21, point.size = 3, plot_type = 'l') +
  geom_point(size = 3, shape = 21, fill = 'white') +
  scale_x_continuous(breaks = seq(1, 10)) +
  theme(aspect.ratio = 1,
        axis.ticks.length = unit(2, 'mm'),
        plot.title = element_text(color = 'black', hjust = .5, face = 'bold'))
ggsave('pam.ncluster_CH.pdf', width = 4, height = 4)

#### Fig. 2f ####
cluster <- PAM_clustering(dist, k = 3)
cluster <- data.frame(sample = rownames(counts), cluster = cluster) %>%
  left_join(distinct(genome_info, genome, domain, phylum, class, order, family, genus, species),
            by = c('sample' = 'genome'))
write.table(cluster, 'pam.3cluster.data.tsv', sep = '\t', row.names = F, quote = F)

point <- data.frame(pcoa$vectors[,1:2]) %>% 
  dplyr::rename_with(~ c('X1', 'X2')) %>% 
  add_column(genome = rownames(counts), .before = 1) %>% 
  add_column(cluster = cluster$cluster) %>% 
  mutate(cluster = factor(cluster))

ggscatter(point, 'X1', 'X2', color = 'cluster', legend = 'right', size = 1.5,
          xlab = 'PCoA1 (14.0%)', ylab = 'PCoA2 (10.9%)', palette = c('#EA8379','#7DAEE0','#B395BD')) +
  theme(aspect.ratio = 3/4,
        axis.line = element_blank(),
        axis.ticks.length = unit(2, 'mm'), 
        panel.grid.major = element_line(linewidth = .8, color = 'grey88'),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = .8, color = 'black', fill = 'transparent'))
ggsave('pam.3cluster.scatter_plot.pdf', width = 9, height = 6)

#### Fig. 2g ####
library(ggtern)

test_data <- rownames_to_column(counts, 'genome') %>% 
  mutate(cluster = cluster$cluster[match(genome, cluster$sample)],
         cluster = paste0('cluster_', cluster)) %>% 
  dplyr::select(-genome) %>% 
  aggregate(. ~ cluster, ., sum) %>% 
  column_to_rownames('cluster') %>% 
  t %>% 
  data.frame() %>% 
  rownames_to_column('CF')

plot_point <- mutate(
  test_data,
  cluster_1 = cluster_1 / sum(cluster_1) * 100,
  cluster_2 = cluster_2 / sum(cluster_2) * 100,
  cluster_3 = cluster_3 / sum(cluster_3) * 100,
  avg = (cluster_1 + cluster_2 + cluster_3) / 3)

cramer_test <- rstatix::cramer_v(column_to_rownames(test_data, 'CF'))
chi_test <- rstatix::chisq_test(column_to_rownames(test_data, 'CF'))

ggtern(plot_point, aes(cluster_1, cluster_2, cluster_3)) +
  geom_point(aes(fill = avg), shape = 21, size = 4) +
  geom_text(aes(label = CF), size = 2) +
  labs(fill = 'Avg. prevalence', 
       title = paste0("Cramér's V: ", signif(cramer_test, 3), ', ', 'P < 2.2e-16')) +
  scale_fill_gradientn(colours = rev(pald('Spectral'))) +
  labs(x = 'CF prevalence in Cluster-1',
       y = 'CF prevalence in Cluster-2',
       z = 'CF prevalence in Cluster-3') +
  ggtern::theme_bw() +
  theme(axis.ticks.length = unit(2, 'mm'))
ggsave('pam.3cluster.ternary_plot.pdf', width = 10, height = 8)

# calcu proportion difference for CF
CF_name <- filter(rowwise(test_data), sum(c(cluster_1 >=5 , cluster_2 >=5, cluster_3 >=5)) >=1 ) %>% 
  pull(CF)

chisq_test <- map_dfr(
  CF_name, ~ 
    mutate(test_data, CF = ifelse(CF != .x, 'Other_CFs', .x)) %>% 
    aggregate(. ~ CF, ., sum) %>% 
    column_to_rownames('CF') %>%
    rstatix::chisq_test() %>% 
    add_column(CF = .x,  .before = 1) ) %>% 
  mutate(padj = p.adjust(p, 'BH'), .after = 'p')

pairwise_chisq_test <- map_dfr(
  CF_name, ~ 
    mutate(test_data, CF = ifelse(CF != .x, 'Other_CFs', .x)) %>% 
    aggregate(. ~ CF, ., sum) %>% 
    mutate(CF = factor(CF, c(.x, 'Other_CFs'))) %>% 
    column_to_rownames('CF') %>%
    rstatix::pairwise_chisq_gof_test(p.adjust.method = 'BH') %>% 
    add_column(CF = .x, .before = 1) )

list(
  prop = plot_point,
  chisq_test = chisq_test, 
  pairwise_chisq_test = pairwise_chisq_test
) %>% openxlsx::write.xlsx('pam.3cluster.CF.chisq_test.xlsx')

#### Fig. 2h ####
# calcu proportion difference for phylum
test_data <- mutate(
  cluster, 
  cluster = paste0('cluster_', cluster)) %>% 
  count(phylum, cluster) %>% 
  spread('cluster', 'n', fill = 0)

cramer_test <- rstatix::cramer_v(column_to_rownames(test_data, 'phylum'))
chi_test <- rstatix::chisq_test(column_to_rownames(test_data, 'phylum'))

mutate(cluster, 
       cluster = paste0('cluster_', cluster),
       phylum = fct_lump_n(phylum, n = 11, other_level = 'Other phyla') ) %>% 
  count(phylum, cluster) %>% 
  ggbarplot('cluster', 'n', fill = 'phylum', legend = 'right', xlab = '', 
            ylab = 'Number of microbial species', width = .7,
            subtitle = paste0("Cramér's V: ", signif(cramer_test, 3), ', ', 'P < 2.2e-16')) +
  scale_fill_manual(values = c('#fb8072','#80b1d3','#ffffb3','#fccde5','#ffed6f','#fdb462',
                               '#b3de69','#8dd3c7','#bebada','#bc80bd','#ccebc5','grey77')) +
  theme(aspect.ratio = 1.7,
        axis.ticks.length = unit(2, 'mm'),
        plot.title = element_text(color = 'black', hjust = .5, face = 'bold'))
ggsave('pam.3cluster.phylum.barplot.pdf', width = 6, height = 5)

prop_data <- mutate(
  test_data,
  cluster_1 = cluster_1 / sum(cluster_1) * 100,
  cluster_2 = cluster_2 / sum(cluster_2) * 100,
  cluster_3 = cluster_3 / sum(cluster_3) * 100)

phylum_name <- filter(rowwise(test_data), sum(c(cluster_1 >=5 , cluster_2 >=5, cluster_3 >=5)) >=1 ) %>% 
  pull(phylum)

chisq_test <- map_dfr(
  phylum_name, ~ 
    mutate(test_data, phylum = ifelse(phylum != .x, 'Other_phyla', .x)) %>% 
    aggregate(. ~ phylum, ., sum) %>% 
    mutate(phylum = factor(phylum, c(.x, 'Other_phyla'))) %>% 
    arrange(phylum) %>% 
    column_to_rownames('phylum') %>%
    rstatix::chisq_test() %>% 
    add_column(phylum = .x, .before = 1) ) %>% 
  mutate(padj = p.adjust(p, 'BH'), .after = 'p')

pairwise_chisq_test <- map_dfr(
  phylum_name, ~ 
    mutate(test_data, phylum = ifelse(phylum != .x, 'Other_phyla', .x)) %>% 
    aggregate(. ~ phylum, ., sum) %>% 
    mutate(phylum = factor(phylum, c(.x, 'Other_phyla'))) %>% 
    arrange(phylum) %>% 
    column_to_rownames('phylum') %>%
    rstatix::pairwise_chisq_gof_test(p.adjust.method = 'BH') %>% 
    add_column(phylum = .x, .before = 1) ) 

list(
  prop_data  = prop_data,
  chisq_test = chisq_test, 
  pairwise_chisq_test = pairwise_chisq_test
) %>% openxlsx::write.xlsx('pam.3cluster.phylum.chisq_test.xlsx')

#### Fig. 2i ####
test_data <- mutate(
  cluster,
  cluster = paste0('cluster_', cluster)) %>% 
  count(family, cluster) %>% 
  spread('cluster', 'n', fill = 0)

cramer_test <- rstatix::cramer_v(column_to_rownames(test_data, 'family'))
chi_test <- rstatix::chisq_test(column_to_rownames(test_data, 'family'))

prop_data <- mutate(
  test_data, 
  cluster_1 = cluster_1 / sum(cluster_1) * 100,
  cluster_2 = cluster_2 / sum(cluster_2) * 100,
  cluster_3 = cluster_3 / sum(cluster_3) * 100)

family_name <- filter(rowwise(test_data), sum(c(cluster_1 >=5 , cluster_2 >=5, cluster_3 >=5)) >=1 ) %>% 
  pull(family)

chisq_test <- map_dfr(
  family_name, ~ 
    mutate(test_data, family = ifelse(family != .x, 'Other_families', .x)) %>% 
    aggregate(. ~ family, ., sum) %>% 
    mutate(family = factor(family, c(.x, 'Other_families'))) %>% 
    arrange(family) %>% 
    column_to_rownames('family') %>%
    rstatix::chisq_test() %>% 
    add_column(family = .x, .before = 1) ) %>% 
  mutate(padj = p.adjust(p, 'BH'), .after = 'p')

pairwise_chisq_test <- map_dfr(
  family_name, ~ 
    mutate(test_data, family = ifelse(family != .x, 'Other_families', .x)) %>% 
    aggregate(. ~ family, ., sum) %>% 
    mutate(family = factor(family, c(.x, 'Other_families'))) %>% 
    arrange(family) %>% 
    column_to_rownames('family') %>%
    rstatix::pairwise_chisq_gof_test(p.adjust.method = 'BH') %>% 
    add_column(family = .x, .before = 1) )

list(
  prop_data = prop_data,
  chisq_test = chisq_test, 
  pairwise_chisq_test = pairwise_chisq_test
) %>% openxlsx::write.xlsx('landscape2/pam.3cluster.family.chisq_test.xlsx')

plot_text <- filter(chisq_test, padj < 0.05) %>% 
  filter(family != 'f__') %>% 
  mutate(plab = add_plab(padj))

plot_prop <- test_data %>% 
  mutate(cluster1 = cluster_1 / sum(cluster_1),
         cluster2 = cluster_2 / sum(cluster_2),
         cluster3 = cluster_3 / sum(cluster_3)
  ) %>% 
  filter(family %in% plot_text$family) %>% 
  dplyr::select(family, cluster1, cluster2, cluster3) %>%
  rowwise() %>% 
  mutate(enriched = ifelse(cluster1 > cluster2 & cluster1 > cluster3, 'cluster1', 
                           ifelse(cluster2 > cluster1 & cluster2 > cluster3, 'cluster2', 'cluster3'))) %>% 
  group_by(enriched) %>% 
  group_modify(~ arrange(.x, desc(!!sym(.y$enriched))) %>% head(10) ) %>% 
  ungroup() %>% 
  mutate(family = factor(family, family)) %>% 
  dplyr::select(-enriched) %>% 
  gather('group', 'value', -family)

p1 <- ggbarplot(
  plot_prop, 'family', 'value', fill = 'group', position = position_dodge2(),
  palette = c('#EA8379','#7DAEE0','#B395BD'), xlab = '', ylab = 'Proportion', 
  legend = 'right', x.text.angle = 90, color = NA, width = .8,
  subtitle = paste0("Cramér's V: ", signif(cramer_test, 3), ', ', 'P < 2.2e-16')) +
  geom_text(aes(family, 0.21, label = plab), filter(plot_text, family %in% plot_prop$family), 
            inherit.aes = F, angle = 90, color = 'red', vjust = .8) +
  scale_y_continuous(expand = c(.01, .003)) +
  geom_vline(xintercept = c(10.5, 20.5), linetype = 'longdash', linewidth = .3) +
  theme(axis.line = element_blank(),
        axis.ticks.length = unit(2, 'mm'), 
        panel.grid.major = element_line(linewidth = .5, color = 'grey88'),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = .8, color = 'black', fill = 'transparent'))

gg.gap::gg.gap(plot = p1, segments = c(.25, .55), ylim = c(0, .6), 
               rel_heights = c(0.5, 0, 0.1), tick_width = .05)

ggsave('pam.3cluster.family.barplot.pdf', width = 10, height = 6.5)
