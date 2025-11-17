#### Jinxin Meng, 20251028, 20251109 ####
setwd('/data/mengjx/project/10.20250623_IBD_BAC_CF_Landscape/git/Figure1/')
pacman::p_load(tidyverse, ggpubr)
source('../scripts/palette.R')
source('../scripts/transform_rc.R')
source('../scripts/plot_pie.R')

#### data ####
data <- data.table::fread('../pipeline/uhgp.m8.info.bz2') %>% 
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

# these CF not included
drops <- c('CF40', 'CF58', 'CF62') 
data <- filter(data, !CF %in% drops)
data.table::fwrite(data, 'CF_gene.metadata.tsv', sep = '\t')

#### Fig. 1a ####
select(data, CF) %>% 
  count(CF) %>% 
  rename(name = 1) %>% 
  plot_pie(top_n = 14, color = '#000000', hemi = T, start = 90,
           fill = colorRampPalette(pald('Spectral')[3:11])(15))
ggsave('counts.hemi_pie.pdf', width = 4, height = 4)

#### Fig. 1b ####
# for genome
counts <- select(data, genome, CF_gene) %>% 
  count(genome, CF_gene) %>% 
  spread('genome', 'n', fill = 0) %>% 
  column_to_rownames('CF_gene')

rare_data <- vegan::specaccum(t(counts), permutations = 99)

tibble(sites = rare_data$sites, richness = rare_data$richness) %>% 
  ggline('sites', 'richness', xlab = 'Number of genome', 
         ylab = 'Number of CF genes', plot_type = 'l', linewidth = .8) +
  theme(aspect.ratio = 1,
        axis.line = element_blank(),
        axis.ticks.length = unit(2, 'mm'), 
        panel.grid.major = element_line(linewidth = .8, color = 'grey88'),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = .8, color = 'black'))
ggsave('rare.data.genome.pdf', width = 5, height = 5)

# for species 
counts <- select(data, ref_species, CF_gene) %>% 
  count(ref_species, CF_gene) %>% 
  spread('ref_species', 'n', fill = 0) %>% 
  column_to_rownames('CF_gene')

rare_data <- vegan::specaccum(t(counts), permutations = 999, method = 'random')

tibble(sites = rare_data$sites, richness = rare_data$richness) %>% 
  ggline('sites', 'richness', xlab = 'Number of species', 
         ylab = 'Number of CF genes', plot_type = 'l', linewidth = .8) +
  theme(aspect.ratio = 1,
        axis.line = element_blank(),
        axis.ticks.length = unit(2, 'mm'), 
        panel.grid.major = element_line(linewidth = .8, color = 'grey88'),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = .8, color = 'black'))
ggsave('rare.data.species.pdf', width = 5, height = 5)

# for genus 
counts <- select(data, genus, CF_gene) %>% 
  count(genus, CF_gene) %>% 
  spread('genus', 'n', fill = 0) %>% 
  column_to_rownames('CF_gene')

rare_data <- vegan::specaccum(t(counts), permutations = 999, method = 'random')

tibble(sites = rare_data$sites, richness = rare_data$richness) %>% 
  ggline('sites', 'richness', xlab = 'Number of genus', 
         ylab = 'Number of CF genes', plot_type = 'l', linewidth = .8) +
  theme(aspect.ratio = 1,
        axis.line = element_blank(),
        axis.ticks.length = unit(2, 'mm'), 
        panel.grid.major = element_line(linewidth = .8, color = 'grey88'),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = .8, color = 'black'))
ggsave('rare.data.genus.pdf', width = 5, height = 5)

#### Fig. 1c ####
# tile plot
plot_tile <- select(data, phylum, CF) %>% 
  count(phylum, CF) %>% 
  mutate(.n = log10(n)) %>% 
  select(-n) %>% 
  spread('phylum', '.n', fill = 0) %>% 
  column_to_rownames('CF')

CF_level <- rownames(plot_tile)[hclust(dist(plot_tile))$order]
phylum_level <- colnames(plot_tile)[hclust(dist(t(plot_tile)))$order]

CF_info <- read_tsv('pfam.rename.tsv')

rownames_to_column(plot_tile, 'CF') %>% 
  gather('phylum', '.n', -CF) %>% 
  mutate(phylum = factor(phylum, phylum_level),
         CF = factor(CF, CF_level)) %>% 
  ggplot(aes(CF, phylum)) +
  geom_tile(aes(fill = .n)) +
  scale_x_discrete(breaks = CF_level, 
                   labels = paste0(CF_level, ' (', CF_info$pfam[match(CF_level, CF_info$name)], ')')) +
  scale_fill_gradientn(colors = rev(pald('Spectral')[-11])) +
  labs(x = '', y = '', fill = 'log10 gene count') +
  coord_fixed() +
  theme(axis.ticks = element_blank(),
        axis.text = element_text(color = '#000000', size = 10),
        axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = .8, color = 'black'))
ggsave('counts.phylum.tile.pdf', width = 16, height = 8)

# prevalence for CF
data.frame(value = rowSums(plot_tile > 0) / ncol(plot_tile)) %>% 
  rownames_to_column('CF') %>% 
  mutate(CF = factor(CF, CF_level)) %>% 
  ggbarplot('CF', 'value', xlab = '', ylab = 'Prevalence', fill = 'grey',
            legend = 'none', width = 1, color = 'grey') +
  theme(axis.ticks.length = unit(2, 'mm'),
        axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1))
ggsave('counts.phylum.upper.pdf', width = 14, height = 3)

# prevalence for phylum
data.frame(value = colSums(plot_tile > 0) / nrow(plot_tile)) %>% 
  rownames_to_column('phylum') %>% 
  mutate(phylum = factor(phylum, phylum_level)) %>% 
  ggbarplot('phylum', 'value', xlab = '', ylab = 'Prevalence', fill = 'grey',
            legend = 'none', width = 1, color = 'grey', rotate = T) +
  theme(axis.ticks.length = unit(2, 'mm'))
ggsave('counts.phylum.right.pdf', width = 5, height = 5)

#### Fig. 1d ####
plot_data <- select(data, genus, CF) %>% 
  count(genus, CF) %>% 
  mutate(.n = log10(n)) %>% 
  select(-n) %>% 
  spread('genus', '.n', fill = 0) %>% 
  column_to_rownames('CF')

plot_prevalence <- data.frame(value = colSums(plot_data > 0) / nrow(plot_data)) %>% 
  rownames_to_column('genus') %>%
  filter(genus != 'g__') %>% 
  left_join(distinct(select(data, phylum, genus)), by = 'genus') %>% 
  mutate(phylum = fct_lump_n(phylum, n = 10, other_level = 'Other phyla'))
openxlsx::write.xlsx(plot_prevalence, 'counts.genus.boxplot.xlsx')

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

#### Fig. 1e ####
plot_data <- select(data, species, CF) %>% 
  count(species, CF) %>% 
  mutate(.n = log10(n)) %>% 
  select(-n) %>% 
  spread('species', '.n', fill = 0) %>% 
  column_to_rownames('CF')

plot_prevalence <- data.frame(value = colSums(plot_data > 0) / nrow(plot_data)) %>% 
  rownames_to_column('species') %>%
  filter(species != 's__') %>% 
  left_join(distinct(select(data, phylum, species)), by = 'species') %>% 
  mutate(phylum = fct_lump_n(phylum, n = 10, other_level = 'Other phyla'))
openxlsx::write.xlsx(plot_prevalence, 'counts.species.boxplot.xlsx')

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
