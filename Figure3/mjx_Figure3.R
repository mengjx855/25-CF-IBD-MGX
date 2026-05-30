#### Jinxin Meng, 20251028, 20260530 ####
setwd('/data/mengjx/project/10.20250623_IBD_BAC_CF_Landscape/github/Figure3/')
pacman::p_load(tidyverse, ggpubr)
source('../scripts/palette.R')
source('../scripts/calcu_difference.R')

#### Fig. S3a ####
PAM_clustering <- function(dist, k) {
  cluster <- as.vector(cluster::pam(as.dist(dist), k, diss = T)$clustering)
  return(cluster)
}

dist <- read_rds('../Figure2/dist.jaccard.rds')

ncluster <- map_vec(
  1:10, ~ 
    clusterSim::index.G1(
      x = (counts > 0) * 1, cl = PAM_clustering(dist, .x),
      d = dist, centrotypes = 'medoids'
    ),
  .progress = T
)

data.frame(
  k_cluster = seq(2, 10), 
  CH = ncluster[-1]
  ) %>% 
  ggline(
    'k_cluster', 'CH', xlab = 'k clusters', ylab = 'Calinski-Harabasz index',
    title = 'Optimal number of clusters', shape = 21, point.size = 3, 
    plot_type = 'l'
  ) +
  geom_point(size = 3, shape = 21, fill = 'white') +
  scale_x_continuous(breaks = seq(1, 10)) +
  theme(
    aspect.ratio = 1,
    axis.ticks.length = unit(2, 'mm'),
    plot.title = element_text(color = 'black', hjust = .5, face = 'bold')
  )

ggsave('pam.ncluster_CH.pdf', width = 4, height = 4)

#### Fig. 3a ####

counts <- data.table::fread(
  '../pipeline/uhgp.m8.f.drop.spread.species.binary.bz2'
) %>% 
  column_to_rownames('name')

cluster <- PAM_clustering(dist, k = 3)

cluster <- data.frame(
  sample = rownames(counts), 
  cluster = cluster
) %>%
  left_join(
    distinct(genome_info, genome, domain, phylum, class, 
             order, family, genus, species),
    by = c('sample' = 'genome')
  )

write.table(
  cluster, 'pam.3cluster.data.tsv', sep = '\t', 
  row.names = F, quote = F
)

pcoa <- readRDS('../Figure2/pcoa.jaccard.rds')

point <- data.frame(pcoa$vectors[,1:2]) %>% 
  dplyr::rename_with(~ c('X1', 'X2')) %>% 
  add_column(genome = rownames(counts), .before = 1) %>% 
  add_column(cluster = cluster$cluster) %>% 
  mutate(cluster = factor(cluster))

ggscatter(
  point, 'X1', 'X2', color = 'cluster', legend = 'right', size = 1.5,
  xlab = 'PCoA1 (14.0%)', ylab = 'PCoA2 (10.9%)', 
  palette = c('#EA8379','#7DAEE0','#B395BD')
) +
  theme(
    aspect.ratio = 3/4,
    axis.line = element_blank(),
    axis.ticks.length = unit(2, 'mm'), 
    panel.grid.major = element_line(linewidth = .8, color = 'grey88'),
    panel.background = element_blank(),
    panel.border = element_rect(linewidth = .8, color = 'black', fill = 'transparent')
  )

ggsave('pam.3cluster.scatter_plot.pdf', width = 9, height = 6)

#### Fig. 3b ####
library(ggtern)

test_data <- rownames_to_column(counts, 'genome') %>% 
  mutate(
    cluster = cluster$cluster[match(genome, cluster$sample)],
    cluster = paste0('cluster_', cluster)
  ) %>% 
  dplyr::select(-genome) %>% 
  aggregate(. ~ cluster, ., sum) %>% 
  column_to_rownames('cluster') %>% 
  t %>% 
  data.frame() %>% 
  rownames_to_column('CF')

cluster_data <- count(cluster, cluster) %>% 
  mutate(cluster = paste0('cluster_', cluster))

cluster_data <- split(cluster_data$n, cluster_data$cluster)

prop_data <- prop_data <- mutate(
  test_data,
  cluster_1 = cluster_1 / cluster_data$cluster_1 * 100,
  cluster_2 = cluster_2 / cluster_data$cluster_2 * 100,
  cluster_3 = cluster_3 / cluster_data$cluster_3 * 100,
  avg = (cluster_1 + cluster_2 + cluster_3) / 3
)

ggtern::ggtern(prop_data, aes(cluster_1, cluster_2, cluster_3)) +
  geom_point(aes(fill = avg), shape = 21, size = 4) +
  geom_text(aes(label = CF), size = 2) +
  scale_fill_gradientn(colours = rev(pald('Spectral'))) +
  labs(
    x = 'CF prevalence in Cluster-1',
    y = 'CF prevalence in Cluster-2',
    z = 'CF prevalence in Cluster-3'
  ) +
  ggtern::theme_bw() +
  theme(axis.ticks.length = unit(2, 'mm'))

ggsave('pam.3cluster.ternary_plot.pdf', width = 10, height = 8)

# difference, fisher's exact test
fisher_test <- map_dfr(
  test_data$CF, ~ 
    filter(test_data, CF == .x) %>% 
    add_row(
      CF = 'none',
      cluster_1 = cluster_data$cluster_1 - .$cluster_1, 
      cluster_2 = cluster_data$cluster_2 - .$cluster_2, 
      cluster_3 = cluster_data$cluster_3 - .$cluster_3
    ) %>% 
    column_to_rownames('CF') %>%
    rstatix::fisher_test(detailed = T, simulate.p.value = T, B = 10000) %>% 
    add_column(CF = .x,  .before = 1) ) %>% 
  mutate(padj = p.adjust(p, 'BH'), .after = 'p')

# compare with other, and enrichment fold change
pairwise_fisher_test <- map_dfr(
  test_data$CF, ~ {
    
    .data <- filter(test_data, CF == .x) %>% 
      add_row(
        CF = 'none',
        cluster_1 = cluster_data$cluster_1 - .$cluster_1, 
        cluster_2 = cluster_data$cluster_2 - .$cluster_2, 
        cluster_3 = cluster_data$cluster_3 - .$cluster_3
      ) %>% 
      column_to_rownames('CF')
    
    efc <- map_dfr(
      c('cluster_1', 'cluster_2', 'cluster_3'), \(x) {
        
        .other <- setdiff(c('cluster_1','cluster_2','cluster_3'), x)
        
        .test_data <- mutate(.data, others = .data[[.other[1]]] + .data[[.other[2]]]) %>% 
          select(all_of(x), others)
        
        .prop <- .test_data[.x,] / colSums(.test_data)
        
        rstatix::fisher_test(
          .test_data, detailed = T, 
          simulate.p.value = T, B = 10000
          ) %>% 
          add_column(
            efc = unlist(.prop[1] / .prop[2], use.names = F), 
            .after = 'estimate'
            ) %>% 
          add_column(name = .x, cluster = x, .before = 1)
      }
    ) %>%
      select(-p.signif) %>% 
      mutate(
        padj = p.adjust(p, 'BH'),
        padj.signif = add_plab(padj), .after = 'p'
      )
  }
)

list(
  prop = prop_data,
  fisher_test = fisher_test, 
  pairwise_fisher_test = pairwise_fisher_test
) %>% 
  openxlsx::write.xlsx('pam.3cluster.CF.fisher_test.xlsx')

#### Fig. 3c ####
cluster <- read_tsv('pam.3cluster.data.tsv')

profile <- data.table::fread('../pipeline/KO.profile.bz2') %>% 
  column_to_rownames('name') %>% 
  apply(2, \(x) ifelse(x > 0, 1, 0)) %>%
  data.frame() %>% 
  select(all_of(cluster$sample)) %>% 
  filter(rowSums(.) != 0)

dist_y <- vegan::vegdist(t(profile), method = 'jaccard', binary = T) # KO_profile
# write_rds(dist_y, 'KO.profile.jaccard.rds')

# Procrustes
dist_x <- read_rds('KO.profile.jaccard.rds')
dist_y <- read_rds('../Figure2/dist.jaccard.rds')

PCoA_x <- cmdscale(dist_x)
PCoA_y <- cmdscale(dist_y)

proc <- vegan::procrustes(PCoA_x, PCoA_y, symmetric = T)
proc_test <- vegan::protest(PCoA_x, PCoA_y, permutations = 999)

proc_point <- bind_cols(
  dplyr::rename(data.frame(proc$Yrot), X1_rotated = X1, X2_rotated = X2), 
  dplyr::rename(data.frame(proc$X), X1_target = Dim1, X2_target = Dim2)
) |> 
  rownames_to_column('sample') |> 
  left_join(select(cluster, sample, cluster))

proc_centroid <- proc_point |> 
  group_by(cluster) |> 
  summarise(across(where(is.numeric), ~ median(.x, na.rm = T)))

proc_coord <- data.frame(proc$rotation)

# mantel
mantel_test <- vegan::mantel(
  dist_x, dist_y, method = 'spearman', 
  permutations = 999, parallel = 110
)
write_rds(mantel_test, 'KO.profile.mantel_test.rds')

# 绘图
ggplot(proc_centroid) +
  geom_segment(
    aes(color = as.character(cluster), 
        x = X1_rotated, y = X2_rotated, 
        xend = X1_target, yend = X2_target), 
    arrow = arrow(length = unit(.3, 'cm')), linewidth = .5) + 
  geom_point(
    aes(X1, X2, shape = type),
    pivot_longer(
      proc_centroid, cols = -cluster, 
      names_to = c('.value', 'type'), names_sep = '_'
    ) |> 
      mutate(type = case_when(
        type == 'rotated' ~ 'CF_profile',
        type == 'target' ~ 'KO_profile'
      )
    ), 
    size = 3, color = 'grey66') +
  scale_shape_manual(values = c(17, 16)) +
  labs(
    x = 'Dim 1', y = 'Dim 2', 
    title = 'Concordance between CF and KO functional spaces', 
    subtitle = 'Procrustes M2 = 0.5153, p=0.001\nMantel r = 0.686, p=0.001') +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    axis.ticks = element_line(linewidth = .5, color = 'black'),
    axis.ticks.length = unit(2, 'mm'),
    axis.title = element_text(size = 12, color = 'black'),
    axis.text = element_text(size = 12, color = 'black'),
    axis.line = element_blank(),
    plot.title = element_text(hjust = .5, size = 12, face = 'bold'),
    plot.subtitle = element_text(hjust = .5, size = 12, color = 'black'),
    plot.margin = unit(c(2, 2, 2, 2), 'mm'),
    panel.border = element_rect(linewidth = .5, color = 'black', fill = NA),
    panel.background = element_blank(),
    panel.grid = element_blank(), 
    legend.background = element_blank(),
    legend.text = element_text(size = 10, color = 'black'),
    legend.title = element_text(size = 10, color = 'black')
  )

ggsave('KO.profile.mantel_test.pdf', width = 5.5, height = 5)

#### Fig. 3d ####
cluster <- read_tsv('pam.3cluster.data.tsv')

profile <- data.table::fread('../pipeline/KO.profile.bz2') %>% 
  column_to_rownames('name') %>% 
  apply(2, \(x) ifelse(x > 0, 1, 0)) %>%
  data.frame() %>% 
  select(all_of(cluster$sample)) %>% 
  filter(rowSums(.) != 0)

# difference
test_data <- data.frame(t(profile)) %>% 
  mutate(cluster = cluster$cluster[match(rownames(.), cluster$sample)]) %>% 
  aggregate(. ~ cluster, ., sum) %>% 
  mutate(cluster = paste0('cluster_', cluster)) %>% 
  column_to_rownames('cluster') %>% 
  t %>% data.frame() %>% 
  rownames_to_column('name')

prop_data <- mutate(
  test_data,
  cluster_1 = cluster_1 / cluster_data$cluster_1 * 100,
  cluster_2 = cluster_2 / cluster_data$cluster_2 * 100,
  cluster_3 = cluster_3 / cluster_data$cluster_3 * 100
)

# cluster with others 
pairwise_fisher_test <- map_dfr(
  test_data$name, \(x) {
    .data <- filter(test_data, name == x) %>% 
      add_row(name = 'none',
              cluster_1 = cluster_data$cluster_1 - .$cluster_1, 
              cluster_2 = cluster_data$cluster_2 - .$cluster_2, 
              cluster_3 = cluster_data$cluster_3 - .$cluster_3)
    
    map_dfr(c('cluster_1','cluster_2','cluster_3'), \(y) {
      
      .other <- setdiff(c('cluster_1','cluster_2','cluster_3'), y)
      
      .test_data <- column_to_rownames(.data, 'name') %>%
        mutate(others = .data[[.other[1]]] + .data[[.other[2]]]) %>% 
        select(all_of(y), others)
      
      .prop_data <- .test_data[1,] / colSums(.test_data) * 100
      
      rstatix::fisher_test(.test_data, detailed = T, simulate.p.value = T, B = 10000) %>% 
        add_column(name = x,
                   cluster = y,
                   prop_main = .prop_data[[y]],
                   prop_other = .prop_data[['others']], 
                   efc = .prop_data[[1]] / .prop_data[[2]], .before = 1) }) %>%
      select(-p.signif) %>% 
      mutate(
        padj = p.adjust(p, 'BH'),
        padj.signif = add_plab(padj), .after = 'p')
    
  }, .progress = T)

openxlsx::write.xlsx(pairwise_fisher_test, 'KO.profile.pairwise.fisher_test.xlsx')

difference <- group_by(pairwise_fisher_test, name) %>% 
  group_modify(~ filter(.x, padj < 0.05 & efc > 2) %>%
                 arrange(desc(efc)) %>% 
                 head(n = 1)) %>% 
  ungroup()

count(difference, cluster) |> 
  plot_pie(
    fill = c('#fed439','#709ae1','#d2af81'), 
    add_n = T, font_size = 4, start = 30)
ggsave('KO.profile.pairwise.pieplot.pdf', width = 4, height = 4)

#### Fig. 3e ####
# KEGG terms enrichment analysis
kegg_db <- read.delim('../pipeline/KO_level_A_B_C_D_Description.bz2') %>% 
  filter(!lvAdes %in% c('Organismal Systems', 'Human Diseases')) %>% 
  mutate(gene = lvD, term = paste0('[', lvC, '] ', lvCdes)) %>% 
  select(term, gene)

results <- map(
  c('cluster_1', 'cluster_2', 'cluster_3'), ~ 
    
    filter(difference, cluster == .x & efc > 2) %>%
    pull(name) %>% 
    clusterProfiler::enricher(
      minGSSize = 1, maxGSSize = 2000, TERM2GENE = kegg_db, 
      pvalueCutoff = .05, qvalueCutoff = .05
    ) %>% 
    data.frame
) %>% 
  set_names(c('cluster_1', 'cluster_2', 'cluster_3'))

openxlsx::write.xlsx(results, 'KO.profile.ORA.results.xlsx')

results <- list(
  cluster_1 = read.xlsx('KO.profile.ORA.results.xlsx', sheet = 1) |> 
    add_column(cluster = 'cluster_1'),
  cluster_2 = read.xlsx('KO.profile.ORA.results.xlsx', sheet = 2) |> 
    add_column(cluster = 'cluster_2'),
  cluster_3 = read.xlsx('KO.profile.ORA.results.xlsx', sheet = 3) |> 
    add_column(cluster = 'cluster_3')
)

results <- do.call(rbind, results) |> 
  data.frame(row.names = NULL) |> 
  relocate(cluster)

data <- results |> 
  mutate(
    path_id = str_extract(ID, 'ko\\d+'), 
    name = gsub('.*\\] ', '', ID),
    name = gsub(' \\[.*', '', name),
    name = paste0(path_id, ', ', name),
    .after = 1
  ) |> 
  filter(
    !is.na(path_id),
    !grepl(' - ', ID), 
    !grepl('BR', ID),  
    !grepl(c('unknown|prediction'), ID)
  )

plot_data <- slice_min(
  data, 
  n = 20,
  order_by = p.adjust, 
  by = cluster) |> 
  select(name, cluster, FoldEnrichment, p.adjust) |> 
  arrange(cluster, desc(FoldEnrichment))

# use itol vis
ggscatter(
  plot_data, 'name', 'cluster', fill = 'p.adjust', 
  shape = 21, size = 'FoldEnrichment', x.text.angle = 45,
  xlab = '', legend = 'right', title = 'enriched KEGG Terms') +
  scale_size_continuous(range = c(3, 8)) +
  scale_fill_distiller(palette = 'Spectral', labels = scales::label_scientific()) +
  theme(
    aspect.ratio = 1/4,
    axis.line = element_blank(),
    axis.ticks.length = unit(2, 'mm'),
    axis.text = element_text(color = 'black', size = 11),
    plot.title = element_text(face = 'bold', hjust = .5),
    panel.border = element_rect(fill = NA, color = 'black', linewidth = .5),
    panel.grid.major = element_line(color = 'grey88', linewidth = .5)
  )

ggsave('KO.profile.ORA.results.top20.dotchart.pdf', 
       width = 15, height = 7)

#### Fig. 3f ####
# calcu proportion difference for phylum
test_data <- mutate(
  cluster, 
  cluster = paste0('cluster_', cluster)) %>% 
  count(phylum, cluster) %>% 
  spread('cluster', 'n', fill = 0)

fisher_test <- rstatix::fisher_test(
  column_to_rownames(test_data, 'phylum'), 
  simulate.p.value = T, B = 10000, detailed = T)

cramer_test <- effectsize::cramers_v(column_to_rownames(test_data, 'phylum'))

mutate(
  cluster, 
  cluster = paste0('cluster_', cluster),
  phylum = fct_lump_n(phylum, n = 11, other_level = 'Other phyla')
) %>% 
  count(phylum, cluster) %>% 
  ggbarplot(
    'cluster', 'n', fill = 'phylum', legend = 'right', xlab = '', 
    ylab = 'Number of microbial species', width = .7,
    subtitle = paste0("Cramér's V: ", signif(cramer_test, 3), ', ', 'P < 2.2e-16')
  ) +
  scale_fill_manual(
    values = c(
      '#fb8072','#80b1d3','#ffffb3','#fccde5','#ffed6f','#fdb462',
      '#b3de69','#8dd3c7','#bebada','#bc80bd','#ccebc5','grey77'
      )
  ) +
  theme(
    aspect.ratio = 1.7,
    axis.ticks.length = unit(2, 'mm'),
    plot.title = element_text(color = 'black', hjust = .5, face = 'bold')
  )
ggsave('pam.3cluster.phylum.barplot.pdf', width = 6, height = 5)

prop_data <- mutate(
  test_data,
  cluster_1 = cluster_1 / sum(cluster_1) * 100,
  cluster_2 = cluster_2 / sum(cluster_2) * 100,
  cluster_3 = cluster_3 / sum(cluster_3) * 100)

fisher_test <- map_dfr(
  test_data$phylum, ~ 
    mutate(test_data, phylum = ifelse(phylum != .x, 'Other_phyla', .x)) %>% 
    aggregate(. ~ phylum, ., sum) %>% 
    mutate(phylum = factor(phylum, c(.x, 'Other_phyla'))) %>% 
    arrange(phylum) %>% 
    column_to_rownames('phylum') %>%
    rstatix::fisher_test(detailed = T, simulate.p.value = T, B = 10000) %>% 
    add_column(phylum = .x, .before = 1) ) %>% 
  mutate(padj = p.adjust(p, 'BH'), .after = 'p')

pairwise_fisher_test <- map_dfr(
  test_data$phylum, ~ 
    mutate(test_data, phylum = ifelse(phylum != .x, 'Other_phyla', .x)) %>% 
    aggregate(. ~ phylum, ., sum) %>% 
    mutate(phylum = factor(phylum, c(.x, 'Other_phyla'))) %>% 
    arrange(phylum) %>% 
    column_to_rownames('phylum') %>%
    rstatix::pairwise_fisher_test(p.adjust.method = 'BH', detailed = T, 
                                  simulate.p.value = T, B = 10000) %>% 
    add_column(phylum = .x, .before = 1) )

list(
  prop_data  = prop_data,
  fisher_test = fisher_test, 
  pairwise_fisher_test = pairwise_fisher_test
) %>% openxlsx::write.xlsx('landscape3/pam.3cluster.phylum.fisher_test.xlsx')

#### Fig. 3g ####
test_data <- mutate(
  cluster,
  cluster = paste0('cluster_', cluster)
) %>% 
  count(family, cluster) %>% 
  spread('cluster', 'n', fill = 0)

cramer_test <- effectsize::cramers_v(column_to_rownames(test_data, 'family'))
fisher_test <- rstatix::fisher_test(
  column_to_rownames(test_data, 'family'), 
  simulate.p.value = T, B = 10000
)

prop_data <- mutate(
  test_data, 
  cluster_1 = cluster_1 / sum(cluster_1) * 100,
  cluster_2 = cluster_2 / sum(cluster_2) * 100,
  cluster_3 = cluster_3 / sum(cluster_3) * 100)

fisher_test <- map_dfr(
  test_data$family, ~
    mutate(test_data, family = ifelse(family != .x, 'Other_families', .x)) %>%
    aggregate(. ~ family, ., sum) %>%
    mutate(family = factor(family, c(.x, 'Other_families'))) %>%
    arrange(family) %>%
    column_to_rownames('family') %>%
    rstatix::fisher_test(detailed = T, simulate.p.value = T, B = 10000) %>%
    add_column(family = .x, .before = 1) ) %>%
  mutate(padj = p.adjust(p, 'BH'), .after = 'p')

pairwise_fisher_test <- map_dfr(
  test_data$family, ~
    mutate(test_data, family = ifelse(family != .x, 'Other_families', .x)) %>%
    aggregate(. ~ family, ., sum) %>%
    mutate(family = factor(family, c(.x, 'Other_families'))) %>%
    arrange(family) %>%
    column_to_rownames('family') %>%
    rstatix::pairwise_fisher_test(p.adjust.method = 'BH', detailed = T,
                                  simulate.p.value = T, B = 10000) %>%
    add_column(family = .x, .before = 1) )

list(
  prop_data = prop_data,
  fisher_test = fisher_test,
  pairwise_fisher_test = pairwise_fisher_test
) %>% openxlsx::write.xlsx('pam.3cluster.family.fisher_test.xlsx')

plot_text <- filter(fisher_test, padj < 0.05) %>% 
  filter(family != 'f__') %>% 
  mutate(plab = add_plab(padj))

plot_prop <- prop_data %>% 
  filter(family %in% plot_text$family) %>% 
  rowwise() %>% 
  mutate(enriched = ifelse(cluster_1 > cluster_2 & cluster_1 > cluster_3, 'cluster_1', 
                           ifelse(cluster_2 > cluster_1 & cluster_2 > cluster_3, 'cluster_2', 'cluster_3'))) %>% 
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
  subtitle = paste0("Cramér's V: ", signif(cramer_test$Cramers_v_adjusted, 3), ', ', 'P < 0.0001')) +
  geom_text(
    aes(family, 21, label = plab), filter(plot_text, family %in% plot_prop$family), 
    inherit.aes = F, angle = 90, color = 'red', vjust = .8) +
  scale_y_continuous(expand = c(.01, .003)) +
  geom_vline(xintercept = c(10.5, 20.5), linetype = 'longdash', linewidth = .3) +
  theme(
    axis.line = element_blank(),
    axis.ticks.length = unit(2, 'mm'), 
    panel.grid.major = element_line(linewidth = .5, color = 'grey88'),
    panel.background = element_blank(),
    panel.border = element_rect(linewidth = .8, color = 'black', fill = 'transparent')
  )

gg.gap::gg.gap(plot = p1, segments = c(25, 55), ylim = c(0, 60), 
               rel_heights = c(0.5, 0, 0.1), tick_width = 5)
ggsave('pam.3cluster.family.barplot.pdf', width = 10, height = 6.5)
