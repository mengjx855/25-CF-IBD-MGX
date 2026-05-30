#### Jinxin Meng, 20251028, 20260530 ####
setwd('/data/mengjx/project/10.20250623_IBD_BAC_CF_Landscape/github/Figure2/')
pacman::p_load(tidyverse, ggpubr)
source('../scripts/palette.R')
source('../scripts/calcu_difference.R')

#### Fig. 2a ####
genome_info <- data.table::fread('../pipeline/genomes-all_metadata.tsv.bz2') %>% 
  mutate(
    genome = sub('\\.\\d+$', '', Genome),
    domain = str_split_i(Lineage, pattern = ';', 1),
    phylum = str_split_i(Lineage, pattern = ';', 2),
    class = str_split_i(Lineage, pattern = ';', 3),
    order = str_split_i(Lineage, pattern = ';', 4),
    family = str_split_i(Lineage, pattern = ';', 5),
    genus = str_split_i(Lineage, pattern = ';', 6),
    species = str_split_i(Lineage, pattern = ';', 7)
  )

counts <- data.table::fread(
  '../pipeline/uhgp.m8.f.drop.spread.species.binary.bz2'
  ) %>% 
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

ggscatter(
  point, 'X1', 'X2', color = 'phylum', legend = 'right', size = 1.5,
  xlab = 'PCoA1 (14.0%)', ylab = 'PCoA2 (10.9%)'
  ) +
  scale_color_manual(
    values = c(
      '#fb8072','#80b1d3','#ffffb3','#fccde5','#ffed6f','#fdb462',
      '#b3de69','#8dd3c7','#bebada','#bc80bd','#ccebc5','grey77'
      )
    ) +
  theme(
    aspect.ratio = 3/4,
    axis.line = element_blank(),
    axis.ticks.length = unit(2, 'mm'), 
    panel.grid.major = element_line(linewidth = .8, color = 'grey88'),
    panel.background = element_blank(),
    panel.border = element_rect(linewidth = .8, color = 'black', fill = 'transparent')
  )

ggsave('pcoa.jaccard.scatter_plot.pdf', width = 9, height = 6)

#### Fig. 2b ####
point <- data.frame(pcoa$vectors[, 1:2]) %>% 
  dplyr::rename_with(~ c('X1', 'X2')) %>% 
  add_column(genome = rownames(counts), .before = 1) %>% 
  mutate(family = genome_info$family[match(genome, genome_info$genome)],
         family = fct_lump_n(family, n = 19, ties.method = 'first', other_level = 'Other families'))

ggscatter(
  point, 'X1', 'X2', color = 'family', legend = 'right', size = 1.5,
  xlab = 'PCoA1 (14.0%)', ylab = 'PCoA2 (10.9%)'
  ) +
  scale_color_manual(
    values = c(
      '#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd',
      '#8c564b','#e377c2','#7f7f7f','#bcbd22','#17becf',
      '#aec7e8','#ffbb78','#98df8a','#ff9896','#c5b0d5',
      '#c49c94','#f7b6d2','#dbdb8d','#9edae5','#c7c7c7'
    )
  ) +
  theme(
    aspect.ratio = 3/4,
    axis.line = element_blank(),
    axis.ticks.length = unit(2, 'mm'), 
    panel.grid.major = element_line(linewidth = .8, color = 'grey88'),
    panel.background = element_blank(),
    panel.border = element_rect(linewidth = .8, color = 'black', fill = 'transparent')
  )

ggsave('pcoa.jaccard.scatter_plot.pdf', width = 9, height = 6)

#### Fig. 2c ####
calcu_adjusted_r2 <- function(adonis_object) {
  n_observations <- adonis_object$Df[3]+1
  d_freedom <- adonis_object$Df[1]
  r2 <- adonis_object$R2[1]
  adjusted_r2 <- vegan::RsquareAdj(r2, n_observations, d_freedom)
  return(adjusted_r2)
}

metadata <- tibble(genome = rownames(counts)) %>%
  left_join(
    dplyr::select(genome_info, genome, domain, phylum, 
                  class, order, family, genus, species) %>%
      distinct(), by = 'genome')

taxa_level <- c('domain', 'phylum','class','order','family','genus')

test <- map(
  taxa_level, ~ 
    vegan::adonis2(
      as.formula(paste0('dist ~ ', .x)),
      metadata, permutations = 999, parallel = 110)
  ) %>%
  set_names(taxa_level)

map2_vec(
  names(test), test, ~ 
    tibble(
      name = str_to_title(.x), 
      r2adj = calcu_adjusted_r2(.y), 
      pval = .y[1, 5])
  ) %>% 
  mutate(plab = add_plab(pval)) %>% 
  ggbarplot(
    'name', 'r2adj', fill = 'name', xlab = '', ylab = 'Adjusted r2', 
    palette = 'Spectral', width = .6) +
  geom_text(
    aes(label = paste0(plab, '\n', round(r2adj * 100, 2), '%')), 
    color = 'red', hjust = .5, vjust = -1) +
  theme_minimal() +
  theme(
    aspect.ratio = 4/5, 
    axis.line = element_line(linewidth = .5),
    axis.ticks = element_line(linewidth = .5),
    axis.text = element_text(color = '#000000', size = 10),
    panel.grid.major = element_line(linewidth = .5),
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(2, 'mm')
  ) +
  guides(fill = 'none')

ggsave('adonis.test.jaccard.pdf', width = 4, height = 5)

#### Fig. 2d and Fig. S2 ####
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
    data.frame(
      x = as.vector(as.dist(patristic_dist)),
      y = as.vector(as.dist(jaccard_dist))
    ) %>% 
      slice_sample(prop = .0001) %>% 
      ggscatter(
        'x', 'y', color = .y, xlab = 'Patristic distance (phylogenetic)', 
        ylab = 'Jaccard distance (CF-based)', 
        title = paste0('stochastic sample ', .x - 2026, ' time'),
        subtitle = paste0('(proportion: 0.0001, seed: ', .x, ')'),
        add = 'reg.line', add.params = list(color = '#000000', linewidth = .8), 
        cor.coef = T, cor.method = 'spearman'
      ) +
      theme(
        aspect.ratio = 1,
        axis.ticks.length = unit(2, 'mm'),
        plot.title = element_text(color = 'black', hjust = .5, face = 'bold'),
        plot.subtitle = element_text(color = 'black', hjust = .5, face = 'bold')
      )
  }
)
cowplot::plot_grid(plotlist = plots, nrow = 4)
ggsave('bacterial.dist.cor.ramdom20.pdf', width = 20, height = 17)

c(0.42, 0.46, 0.42, 0.41, 0.34, 0.4, 0.34, 0.4, 0.44, 0.45,
  0.39, 0.41, 0.43, 0.4, 0.44, 0.46, 0.44, 0.41, 0.41, 0.4, 0.45) %>% 
  summary()

#### Fig. 2e ####
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
