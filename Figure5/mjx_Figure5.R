#### Jinxin Meng, 20251028, 20251115 ####
setwd('/data/mengjx/project/10.20250623_IBD_BAC_CF_Landscape/github/Figure5/')
pacman::p_load(tidyverse, ggpubr)
source('../scripts/palette.R')
source('../scripts/calcu_difference.R')
source('../scripts/plot_Procrustes.R')
source('../scripts/calcu_metafor-0.1.1.R')

#### Fig. 5a and Fig. S5a-j ####
proj_name <- c(
  'BushmanFD_2020','FranzosaEA_2018','HallAB_2017','HeQ_2017','KumbhariA_2024',
  'LloydPriceJ_2019','SchirmerM_2018','SchirmerM_2024','WengY_2019','YanQ_2023c'
)

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

#### Fig. 5b ####
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

#### Fig. 5c ####
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

#### Fig. 5d ####
test <- map(proj_name, ~ {
  data_x <- data.frame(apply(tpm[[.x]], 2, \(x) x / sum(x))) %>%
    t() %>% vegan::decostand(method = 'hellinger') %>% 
    data.frame()
  
  data_y <- data.table::fread(paste0('../data/', .x, '.kk2.s.bz2')) %>%
    column_to_rownames('name') %>%
    dplyr::select(any_of(rownames(data_x)))
  
  dist_y <- vegan::decostand(t(data_y), method = 'hellinger') %>%
    vegan::vegdist(method = 'bray')
  
  map(
    colnames(data_x), \(x) {
      adonis <- vegan::adonis2(
        as.formula(paste0('dist_y ~ ', x )), data_x, 
        permutations = 999, parallel = 80
      ) 
      
      adonis$r2adj <- c(calcu_adjusted_r2(adonis), NA, NA)
      adonis 
    }
  ) %>%  
    set_names(colnames(data_x)) 
  }
) %>% 
  set_names(proj_name)

write_rds(test, 'adonis.r2.single.rds')

data <- map2_dfr(
  test, names(test), \(x, y) 
  map2_dfr(
    x, names(x), \(i, j) 
    data.frame(
      dataset = y, 
      gene = j, 
      r2adj = i$r2adj[1],
      pval = i$`Pr(>F)`[1]
    ) 
  )
) %>% 
  mutate(plab = add_plab(pval, format = 4))

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

#### Fig. S5n ####
proj_name <- c(
  'BushmanFD_2020','FranzosaEA_2018','HallAB_2017','HeQ_2017',
  'KumbhariA_2024','LloydPriceJ_2019','SchirmerM_2018',
  'SchirmerM_2024','WengY_2019','YanQ_2023c')

profile <- map(
  proj_name, ~ 
    data.table::fread(
      paste0(
        '../data/',
        .x, '.kk2.s.bz2')
    )
) |> 
  reduce(~ full_join(.x, .y, by = 'name')) |> 
  column_to_rownames('name') |> 
  mutate(across(everything(), ~ replace_na(.x, 0)))

genome_info <- data.table::fread('../pipeline/genomes-species_metadata.tsv.bz2') |> 
  mutate(
    genome  = sub('\\.\\d+$', '', Genome),
    genus   = stringr::str_split_i(Lineage, pattern = ';', 6),
    species = stringr::str_split_i(Lineage, pattern = ';', 7),
    genus   = sub('^g__', '', genus),
    species = sub('^s__', '', species), 
    species = case_when(
      species == '' | is.na(species) ~ genome,
      TRUE ~ species
    )
  ) |> 
  select(name = genome, genus, species)

## 3. meta 分析
stats <- map(
  proj_name, ~ {
    
    group <- data.table::fread(
      paste0(
        '../data/', 
        .x, '.sample_group.bz2'
      )
    ) |> 
      dplyr::select(sample, group) |> 
      filter(sample %in% colnames(profile))
    
    comparisons <- map(
      setdiff(unique(group$group), 'HC'), 
      ~ c(.x, 'HC')
    )
    
    ## 反正弦平方根转换
    profile_x <- profile[, group$sample, drop = FALSE]
    profile_x <- asin(sqrt(profile_x))
    
    profile_x <- profile[, group$sample, drop = FALSE]
    rownames(profile_x) <- genome_info$name[match(rownames(profile_x), genome_info$species)]
    
    mean_data <- profile_collapse(profile_x, group, method = 'mean') |> 
      t() |> 
      data.frame(check.names = FALSE)
    
    sd_data <- profile_collapse(profile_x, group, method = 'sd') |> 
      t() |> 
      data.frame(check.names = FALSE)
    
    n_data <- count(group, group) |> 
      column_to_rownames('group')
    
    map(
      comparisons, \(x) 
      data.frame(
        name = rownames(profile_x), 
        m1   = unlist(mean_data[x[1], ], use.names = FALSE),
        sd1  = unlist(sd_data[x[1], ], use.names = FALSE),
        n1   = as.numeric(n_data[x[1], 'n']),
        m2   = unlist(mean_data[x[2], ], use.names = FALSE),
        sd2  = unlist(sd_data[x[2], ], use.names = FALSE),
        n2   = as.numeric(n_data[x[2], 'n'])
      )
    ) |> 
      set_names(
        map_vec(comparisons, \(x) paste0(x, collapse = '_vs_'))
      )
  }, 
  .progress = T
) |> 
  set_names(proj_name)

feature_names <- map(
  stats, ~ 
    map(.x, \(x) pull(x, name))
) |> 
  flatten() |> 
  reduce(~ intersect(.x, .y))

meta_in <- map(
  feature_names, \(x) 
  map_dfr(
    proj_name, \(y)
    map2_dfr(
      stats[[y]], names(stats[[y]]), \(m, n) 
      filter(m, name == x) |>
        add_column(project = paste0(y, '.', n))
    )
  ) |> 
    dplyr::select(-name) |> 
    relocate(project),
  .progress = TRUE
) |> 
  set_names(feature_names)

write_rds(meta_in, 'all_species.meta_in.rds')

meta_out <- map_dfr(
  feature_names, ~ 
    calcu_metafor.1(
      meta_in[[.x]], 
      data_rename = c('name' = 'project'), 
      simplify = T, quiet = T
    ) |> 
    add_column(feature = .x, .before = 1),
  .progress = T
)

meta_out <- meta_out |> 
  mutate(
    padj = p.adjust(pval, method = 'BH'),
    enriched = case_when(
      padj < 0.05 & estimate > 0 ~ 'case',
      padj < 0.05 & estimate < 0 ~ 'control',
      TRUE ~ 'none'
    )
  )

write.xlsx(meta_out, 'all_species.meta_out.xlsx')

metadata <- read.delim('pam.3cluster.CS.tsv')

counts <- data.table::fread('../pipeline/uhgp.m8.f.drop.spread.species.binary.bz2') |> 
  column_to_rownames('name')

richness <- rowSums(counts) |> 
  data.frame(richness = _) |> 
  rownames_to_column('feature')

strategy_richness <- data.frame(t(counts)) |> 
  profile_aggregate(metadata) |> 
  t() |> 
  data.frame() |> 
  rownames_to_column('feature') |> 
  select(-Unknown)

difference <- read.xlsx('all_species.meta_out.xlsx') |> 
  select(feature, enriched)

count(difference, enriched) |> 
  plot_pie(
    fill = c('#fed439','#709ae1','#d2af81'), 
    add_n = T, font_size = 4, start = 30)
ggsave('all_species.meta_out.pie.pdf', width = 4, height = 4)

#### Fig. S5o ####
data <- difference |> 
  left_join(richness, by = 'feature') |> 
  left_join(strategy_richness, by = 'feature') |> 
  filter(enriched != 'none' & !is.na(richness))

test <- map_dfr(
  c("richness","communication","metabolism","resistance"), ~ 
    rstatix::wilcox_test(
      data, as.formula(paste0(.x, ' ~ enriched')), 
      detailed = T) |> 
    data.frame()
)

stats <- left_join(
  data |> 
    group_by(enriched) |> 
    summarise(
      across(
        where(is.numeric), ~ 
          mean(.x, na.rm = T)
      )
    ) |> 
    pivot_longer(
      cols = !enriched, 
      names_to = 'class', 
      values_to = 'mean'), 
  data |> 
    group_by(enriched) |> 
    summarise(
      across(
        where(is.numeric), ~ 
          mean(.x, na.rm = T)
      )
    ) |> 
    pivot_longer(
      cols = !enriched, 
      names_to = 'class', 
      values_to = 'sd'),
  by = c('enriched', 'class')
) |> 
  left_join(
    count(data, enriched),
    by = 'enriched'
  )

list(
  stats = stats, 
  test = test
) |> 
  write.xlsx('all_species.meta_out.species.richness.xlsx')

ggplot(data, aes(x = enriched, y = richness, fill = enriched)) +
  geom_violin(trim = F, alpha = 0.45, color = NA) +
  geom_boxplot(fill = 'white', width = 0.25, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.1, size = .8, alpha = .4, shape = 16) +
  scale_fill_manual(values = c('case' = '#fee722', 'control' = '#5f7dbe')) +
  geom_signif(comparisons = list(c('case', 'control')), step_increase = .05,
              textsize = 4, tip_length = .02, vjust = .2, size = .4) +
  theme_pubr(base_size = 14) +
  labs(
    x = NULL, y = 'CF-family richness per species'
  ) +
  theme(
    aspect.ratio = 3/2,
    axis.ticks.length = unit(2, 'mm'),
    legend.position = 'none'
  )
ggsave('all_species.meta_out.species.all_richness.pdf', width = 3, height = 4)

select(
  data, feature, enriched, communication, metabolism, resistance
) |> 
  pivot_longer(
    cols = !c(feature, enriched), 
    names_to = 'strategy', 
    values_to = 'value'
  ) |> 
  mutate(
    strategy = recode(
      strategy,
      metabolism = 'Metabolism-centric',
      resistance = 'Stress-resistance-oriented',
      communication = 'Communication-mediated'
    )
  ) |> 
  ggplot(aes(x = enriched, y = value, fill = enriched)) +
  geom_violin(trim = F, alpha = 0.45, color = NA) +
  geom_boxplot(fill = 'white', width = 0.25, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.1, size = .8, alpha = .4, shape = 16) +
  facet_wrap(~ strategy, scales = 'free_y') +
  scale_fill_manual(values = c('case' = '#fee722', 'control' = '#5f7dbe')) +
  geom_signif(comparisons = list(c('case', 'control')), step_increase = .05,
              textsize = 4, tip_length = .02, vjust = .2, size = .4) +
  theme_pubr(base_size = 14) +
  labs(
    x = NULL, y = 'Strategy-specific CF-family count',
  ) +
  theme(
    aspect.ratio = 3/2,
    axis.ticks.length = unit(2, 'mm'),
    legend.position = 'none'
  )
ggsave('all_species.meta_out.species.strategy_richness.pdf', width = 9, height = 4)

#### Fig. 5e ####
difference <- read.xlsx('all_species.meta_out.xlsx')

fgsea_in <- pull(
  difference, estimate, name = feature
) |> 
  sort(decreasing = T)

# 构建每个 CF 的 carrier species set
TERM2GENE <- counts |> 
  rownames_to_column('GENE') |> 
  pivot_longer(cols = -GENE, names_to = 'TERM', values_to = 'presence') |> 
  filter(presence == 1) |> 
  select(TERM, GENE)

gsea_result <- clusterProfiler::GSEA(
  fgsea_in, minGSSize = 1, maxGSSize = length(fgsea_in), 
  pvalueCutoff = 1, TERM2GENE = TERM2GENE, eps = 1e-100
)
write_rds(gsea_result, 'gsea_result.rds')

gsea_result <- data.frame(gsea_result, row.names = NULL)
write.xlsx(gsea_result, 'gsea_result.xlsx')

# 火山图
plot_data <- gsea_result |> 
  left_join(metadata, by = c('ID' = 'name')) |> 
  mutate(
    class = case_when(
      p.adjust < 0.05 & group == 'metabolism' ~ 'metabolism',
      p.adjust < 0.05 & group == 'communication' ~ 'communication',
      p.adjust < 0.05 & group == 'resistance' ~ 'resistance',
      p.adjust < 0.05 & is.na(group) ~ 'unclassified',
      T ~ 'other'
    ),
    .padj = -log10(p.adjust),
    .padj = case_when(.padj > 30 ~ .padj * .5, T ~ .padj)
  )

group_colors <- c(
  'metabolism' = '#58bbd0',
  'communication' = '#eb9256',
  'resistance' = '#9b99cb',
  'unclassified' = '#ad9f2a',
  'other' = '#e1ecee')

ggplot(plot_data, aes(x = NES, y = .padj)) +
  geom_point(aes(color = class), size = 3) +
  scale_color_manual(values = group_colors) +
  geom_hline(yintercept = -log10(0.05), color = '#e1ecee', linetype = 'dashed', linewidth = 1) +
  ggrepel::geom_text_repel(aes(label = ID), size = 2, min.segment.length = .1) +
  labs(x = 'Normalized enrichment score', y = '-log10 FDR') +
  theme_pubr() +  
  theme(
    aspect.ratio = 1,
    axis.ticks.length = unit(2, 'mm'),
    legend.position = 'right'
  )

ggsave('taxa_centric/gsea_result.volcano.pdf', width = 7, height = 5)

#### Fig. 5f ####
# leading-edge species
genome_info <- data.table::fread('../pipeline/genomes-species_metadata.tsv.bz2') |> 
  mutate(
    genome  = sub('\\.\\d+$', '', Genome),
    phylum  = stringr::str_split_i(Lineage, pattern = ';', 2),
    family  = stringr::str_split_i(Lineage, pattern = ';', 5),
    genus   = stringr::str_split_i(Lineage, pattern = ';', 6),
    species = stringr::str_split_i(Lineage, pattern = ';', 7)
  ) |> 
  select(genome, phylum, family, genus, species) |> 
  mutate(
    across(
      where(is.character), ~ 
        if_else(grepl('__$', .x), paste0(.x, 'Unclassified'), .x)
    )
  )

data <- gsea_result %>%
  mutate(
    direction = case_when(
      NES > 0 ~ 'case', 
      T ~ 'control'
    )
  ) |> 
  filter(p.adjust < 0.05)

df <- select(data, CF = ID, core_enrichment) |> 
  separate_rows(core_enrichment, sep = '/') |> 
  left_join(
    select(data, CF = ID, NES), 
    by = 'CF'
  ) |> 
  left_join(
    select(difference, core_enrichment = feature, estimate), 
    by = 'core_enrichment'
  ) |> 
  filter(
    NES * estimate > 0
  ) |> 
  rename(genome = core_enrichment) |> 
  left_join(genome_info, by = 'genome') |> 
  select(CF, species, taxon = family) |> 
  dplyr::distinct(CF, species, taxon)

top_taxa <- df %>%
  count(taxon, name = 'n_total', sort = T) %>%
  filter(!grepl('CAG|UBA|QA|HG|UM|DS|RU|Ana|Mar|Mon|Gas|Acid', taxon)) |> 
  slice_head(n = 20) %>%
  pull(taxon)

df2 <- df %>%
  dplyr::mutate(
    taxon = ifelse(taxon %in% top_taxa, taxon, 'f__Other')
  ) %>%
  dplyr::count(CF, taxon, name = 'n_species') |> 
  left_join(select(data, CF = ID, direction), by = 'CF')

n_taxa <- group_by(df2, CF) |> 
  summarise(n_taxa = sum(n_species), .groups = 'drop') |> 
  mutate(
    label = paste0(CF, ' (', n_taxa, ')')
  )

df2 <- left_join(df2, n_taxa, by = 'CF') |> 
  select(-CF) |> 
  rename(CF = label)

CF_levels <- df2 %>%
  group_by(CF) %>%
  summarise(total = sum(n_species), .groups = "drop") %>%
  left_join(select(data, CF = ID, direction), by = 'CF') |> 
  arrange(direction, desc(total)) %>%
  pull(CF)

df2$CF <- factor(df2$CF, levels = CF_levels)

taxa_levels <- df2 %>%
  count(taxon, wt = n_species, sort = TRUE) %>%
  pull(taxon)

df2$taxon <- factor(df2$taxon, levels = taxa_levels)

df2 <- mutate(
  df2, 
  p_species = n_species / n_taxa * 100)

pal <- c('grey92', pald('Paired2', n = 20))
names(pal) <- taxa_levels

ggbarplot(df2, 'CF', 'p_species', fill = 'taxon', width = .8, size = .3, 
          palette = pal, legend = 'right', x.text.angle = 90) +
  facet_grid(cols = vars(direction), scales = 'free', space = 'free') +
  scale_y_continuous(expand = c(.02, .02)) +
  labs(x = '', y = 'Proportion of leading-edge species', fill = 'Host genera') +
  theme(
    axis.ticks = element_line(linewidth = .5),
    axis.ticks.length = unit(1.8, 'mm'),
    axis.line = element_blank(),
    panel.grid.major = element_line(linewidth = .5),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, linewidth = .5),
    strip.background = element_rect(fill = NA, linewidth = .5),
    legend.key.spacing.y = unit(1, 'mm'),
    legend.text = element_text(face = 'italic',size = 11)
  )

ggsave('gsea_result.core.all.pdf', width = 14, height = 6)
write.xlsx(df2, 'gsea_result.core.all.xlsx')

#### Fig. 5g ####
## specific taxa
select_taxa <- c(
  'g__Streptococcus','g__Fusobacterium','g__Enterococcus',
  'g__Faecalibacterium','g__Bifidobacterium','g__Roseburia'
)

df <- select(data, CF = ID, core_enrichment) |> 
  separate_rows(core_enrichment, sep = '/') |> 
  left_join(
    select(data, CF = ID, NES), 
    by = 'CF'
  ) |> 
  left_join(
    select(difference, core_enrichment = feature, estimate), 
    by = 'core_enrichment'
  ) |> 
  filter(
    NES * estimate > 0
  ) |> 
  rename(genome = core_enrichment) |> 
  left_join(genome_info, by = 'genome') |> 
  filter(str_detect(genus, paste(select_taxa, collapse = '|'))) |> 
  mutate(
    species = if_else(
      species == 's__Unclassified', 
      paste0('s__', str_remove(genus, 'g__'), ' ', genome),
      species
    )
  )
write.xlsx(df, 'select_taxa_for_net.xlsx')
