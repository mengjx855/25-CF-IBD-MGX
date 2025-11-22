#### Jinxin Meng, 20251028, 20251122 ####
setwd('/data/mengjx/project/10.20250623_IBD_BAC_CF_Landscape/git/Figure1/')
pacman::p_load(tidyverse, ggpubr)
source('../scripts/palette.R')
source('../scripts/transform_rc.R')
source('../scripts/plot_pie.R')

#### data ####
data <- data.table::fread('../pipeline/uhgp.m8.f.drop.info.bz2') %>% 
  dplyr::select(query = V1, subject = V2, genome = V5, ref_species = V18, lineage = V19) %>% 
  mutate(CF = str_split_i(subject, pattern = ',', 1),
         CF_gene = str_split_i(subject, pattern = ',', 2),
         domain = str_split_i(lineage, pattern = ';', 1),
         phylum = str_split_i(lineage, pattern = ';', 2),
         class = str_split_i(lineage, pattern = ';', 3),
         order = str_split_i(lineage, pattern = ';', 4),
         family = str_split_i(lineage, pattern = ';', 5),
         genus = str_split_i(lineage, pattern = ';', 6),
         species = str_split_i(lineage, pattern = ';', 7) ) %>% 
  select(-lineage)

data.table::fwrite(data, 'CF_gene.metadata.tsv', sep = '\t')

#### Fig. 1a ####
select(data, CF) %>% 
  count(CF) %>% 
  rename(name = 1) %>% 
  plot_pie(top_n = 14, color = '#000000', hemi = T, start = 90,
           fill = colorRampPalette(pald('Spectral')[3:11])(15),
           title = '78,579 CF gene homologs in UHGP, across 79 CFs')
ggsave('counts.uhgp.hemi_pie.pdf', width = 4, height = 4)

#### Fig. 1b ####
data <- read.delim('../pipeline/uhgp.m8.f.drop.spread.CF_count.bz2', header = F, 
                   col.names = c('name', 'n'))
plot_pie(data, top_n = 14, color = '#000000', hemi = T, start = 90,
         fill = colorRampPalette(pald('Spectral')[3:11])(15),
         title = '7,056,385 CF gene homologs in UHGG, across 79 CFs, 
         across 4,716 species (4,716/4,745), across 289,022 genomes (289,022/289,231)')
ggsave('counts.uhgg.hemi_pie.pdf', width = 4, height = 4)

#### Fig. 1c ####
counts <- data.table::fread('../pipeline/uhgp.m8.f.drop.spread.species.binary.for_gene.bz2', check.names = F) %>% 
  column_to_rownames('name') %>% 
  t %>% data.frame()

rare_data <- vegan::specaccum(counts, permutations = 999, method = 'random')

tibble(sites = rare_data$sites,
       richness = rare_data$richness) %>% 
  ggline('sites', 'richness', xlab = 'Number of species', 
         ylab = 'Number of reference genes', plot_type = 'l', linewidth = .8) +
  theme(aspect.ratio = 1,
        axis.line = element_blank(),
        axis.ticks.length = unit(2, 'mm'), 
        panel.grid.major = element_line(linewidth = .8, color = 'grey88'),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = .8, color = 'black', fill = 'transparent'))
ggsave('rare.data.species.pdf', width = 5, height = 5)

#### Fig. 1d ####
genome_info <- data.table::fread('../pipeline/genomes-all_metadata.tsv.bz2') %>% 
  mutate(genome = sub('\\.\\d+$', '', Genome),
         phylum = str_split_i(Lineage, pattern = ';', 2),
         family = str_split_i(Lineage, pattern = ';', 5),
         genus = str_split_i(Lineage, pattern = ';', 6),
         species = str_split_i(Lineage, pattern = ';', 7))

data <- data.table::fread('../pipeline/uhgp.m8.f.drop.spread.species.binary.bz2') %>% 
  column_to_rownames('name')

# tile plot
plot_data <- data %>% 
  mutate(phylum = genome_info$phylum[match(rownames(.), genome_info$genome)]) %>% 
  aggregate(. ~ phylum, ., sum) %>% 
  column_to_rownames('phylum') %>% 
  apply(2, \(x) ifelse(x != 0, log10(x), 0)) %>% 
  data.frame()
plot_data <- plot_data[rowSums(plot_data > 0) != 0, colSums(plot_data > 0) != 0]

phylum_level <- rev(rownames(plot_data)[hclust(dist(plot_data))$order])
CF_level <- rev(colnames(plot_data)[hclust(dist(t(plot_data)))$order])

CF_info <- read_tsv('pfam.rename.tsv')

rownames_to_column(plot_data, 'phylum') %>% 
  gather('CF', 'n', -phylum) %>% 
  mutate(phylum = factor(phylum, phylum_level),
         CF = factor(CF, CF_level)) %>% 
  ggplot(aes(CF, phylum)) +
  geom_tile(aes(fill = n)) +
  scale_x_discrete(
    breaks = CF_level, 
    labels = paste0(CF_level, ' (', CF_info$pfam[match(CF_level, CF_info$name)], ')')) +
  scale_fill_gradientn(colors = rev(pald('Spectral')[-11])) +
  labs(x = '', y = '', fill = 'log10 species count') +
  coord_fixed() +
  theme(axis.ticks = element_blank(),
        axis.text = element_text(color = '#000000', size = 10),
        axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = .8, color = 'black', fill = 'transparent'))
ggsave('counts.phylum.tile.pdf', width = 16, height = 8)

# prevalence for CF
data.frame(value = colSums(plot_data > 0) / nrow(plot_data)) %>% 
  rownames_to_column('CF') %>% 
  mutate(CF = factor(CF, CF_level)) %>% 
  ggbarplot('CF', 'value', xlab = '', ylab = 'Prevalence', fill = 'grey',
            legend = 'none', width = 1, color = 'grey') +
  theme(axis.ticks.length = unit(2, 'mm'),
        axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1))
ggsave('counts.CF.upper.pdf', width = 14, height = 3)

# prevalence for phylum
data.frame(value = rowSums(plot_data > 0) / ncol(plot_data)) %>% 
  rownames_to_column('phylum') %>% 
  mutate(phylum = factor(phylum, phylum_level)) %>% 
  ggbarplot('phylum', 'value', xlab = '', ylab = 'Prevalence', fill = 'grey',
            legend = 'none', width = 1, color = 'grey',  rotate = T) +
  theme(axis.ticks.length = unit(2, 'mm'))
ggsave('counts.phylum.right.pdf', width = 5, height = 5)

#### Fig. 1e ####
plot_data <- data %>% 
  mutate(genus = genome_info$genus[match(rownames(.), genome_info$genome)]) %>% 
  aggregate(. ~ genus, ., sum) %>% 
  column_to_rownames('genus') %>% 
  data.frame()

plot_prevalence <- data.frame(value = rowSums(plot_data > 0) / ncol(plot_data)) %>% 
  rownames_to_column('genus') %>%
  filter(genus != 'g__') %>% 
  mutate(phylum = genome_info$phylum[match(genus, genome_info$genus)],
         phylum = fct_lump_n(phylum, n = 10, other_level = 'Other phyla'))

plot_text <- group_by(plot_prevalence, phylum) %>% 
  slice_max(order_by = value, n = 6)

ggboxplot(plot_prevalence, 'phylum', 'value', color = 'phylum', legend = 'none', 
          xlab = '', ylab = 'Prevalence', x.text.angle = 30, outlier.shape = NA,
          title = 'CF prevalence in genus-level taxa') +
  geom_jitter(aes(color = phylum), width = .3) +
  ggrepel::geom_text_repel(aes(phylum, value, label = genus), plot_text, size = 2) +
  stat_compare_means() +
  scale_color_manual(values = colorRampPalette(pald('Spectral')[-6])(11)) +
  theme(axis.ticks.length = unit(2, 'mm'),
        aspect.ratio = 1/2,
        plot.title = element_text(face = 'bold', hjust = .5))
ggsave('counts.genus.boxplot.pdf', width = 10, height = 6)

#### Fig. 1f ####
plot_prevalence <- data.frame(value = rowSums(data > 0) / ncol(data)) %>% 
  rownames_to_column('name') %>%
  mutate(species = genome_info$species[match(name, genome_info$genome)]) %>% 
  mutate(phylum = genome_info$phylum[match(name, genome_info$genome)],
         phylum = fct_lump_n(phylum, n = 10, other_level = 'Other phyla'))

plot_text <- group_by(plot_prevalence, phylum) %>% 
  slice_max(order_by = value, n = 6)

ggboxplot(plot_prevalence, 'phylum', 'value', color = 'phylum', legend = 'none', 
          xlab = '', ylab = 'Prevalence', x.text.angle = 30, outlier.shape = NA,
          title = 'CF prevalence in species-level taxa') +
  geom_jitter(aes(color = phylum), width = .3) +
  ggrepel::geom_text_repel(aes(phylum, value, label = species), plot_text, size = 1) +
  stat_compare_means() +
  scale_color_manual(values = colorRampPalette(pald('Spectral')[-6])(11)) +
  theme(axis.ticks.length = unit(2, 'mm'),
        aspect.ratio = 1/2,
        plot.title = element_text(face = 'bold', hjust = .5))
ggsave('counts.species.boxplot.pdf', width = 10, height = 6)
