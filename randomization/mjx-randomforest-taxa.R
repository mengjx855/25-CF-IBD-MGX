#### Jinxin Meng, 20250625, 20250520 ####
setwd('/data/mengjx/project/10.20250623_IBD_BAC_CF_Landscape/github/randomization/')
pacman::p_load(tidyverse, ggpubr)
source('scripts/calcu_difference.R')

#### 20260521 random select species 79/46/13 ####
proj_name <- c(
  'BushmanFD_2020','FranzosaEA_2018','HallAB_2017','HeQ_2017','KumbhariA_2024',
  'LloydPriceJ_2019','SchirmerM_2018','SchirmerM_2024','WengY_2019','YanQ_2023c'
)

file_base <- '../data/'

kk2 <- map(
  proj_name, ~
    data.table::fread(
      paste0(file_base, .x, '.kk2.s.bz2')
    )
  ) |> 
  set_names(proj_name)

feature_names <- map(
  kk2, ~ 
    pull(.x, name) 
  ) |> 
  reduce(~ intersect(.x, .y))

seed_seqs <- seq(2026, 3024, 1)

seed_names <- paste0(
  'rep', 
  sprintf('%03d', seq_along(seed_seqs)), 
  '_seed', 
  seed_seqs
)

random_select <- map(
  seed_seqs, ~ {
    set.seed(.x)
    sample(feature_names, size = 13, replace = F) 
  }
) |> 
  set_names(seed_names)

result <- tibble::enframe(
  random_select, name = 'name', value = 'taxa'
) |> 
  rowwise() |> 
  mutate(taxa = paste(taxa, collapse = ','))

write_tsv(
  result, 'results_taxa/select_taxa.species.13.tsv', 
  col_names = F, quote = NULL
)

#### 20260521 random select genus 79/46/13 ####

proj_name <- c(
  'BushmanFD_2020','FranzosaEA_2018','HallAB_2017','HeQ_2017','KumbhariA_2024',
  'LloydPriceJ_2019','SchirmerM_2018','SchirmerM_2024','WengY_2019','YanQ_2023c'
)

file_base <- '../data/'

kk2 <- map(
  proj_name, ~
    data.table::fread(
      paste0(file_base, .x, '.kk2.g.bz2')
    )
  ) |> 
  set_names(proj_name)

feature_names <- map(
  kk2, ~ 
    pull(.x, name) 
  ) |> 
  reduce(~ intersect(.x, .y))

seed_seqs <- seq(2026, 3024, 1)

seed_names <- paste0(
  'rep', 
  sprintf('%03d', seq_along(seed_seqs)), 
  '_seed', 
  seed_seqs
)

random_select <- map(
  seed_seqs, ~ {
    set.seed(.x)
    sample(feature_names, size = 13, replace = F) 
  }
) |> 
  set_names(seed_names)

result <- tibble::enframe(
  random_select, name = 'name', value = 'taxa'
) |> 
  rowwise() |> 
  mutate(taxa = paste(taxa, collapse = ','))

write_tsv(
  result, 'results_taxa/select_taxa.genus.13.tsv', 
  col_names = F, quote = NULL
)

#### 20260521 random test species 79/46/13 ####
files <- list.files('results_taxa/random_species/', pattern = '13.results')

random <- map_dfr(
  files, ~ 
    read.delim(
      paste0('results_taxa/random_species/', .x), 
      header = F, col.names = c('rep', 'auc')
    )
) |> 
  group_by(rep) |> 
  summarise(mean = mean(auc))

# obs for 46 / 79
ref <- read.delim('../Figure6/CF.rf.roc.all_CF_gene.tsv') |>
  pull(auc) |>
  mean()

# obs for 13
ref <- read.delim(
  '../Figure6/roc.top_13.tsv',
) |>
  pull(auc) |>
  mean()

test <- calcu_empirical_p(ref, random$mean, 'right')

ggdensity(
  random, x = 'mean', fill = 'grey77', 
  title = '13 random species', xlab = 'Average AUC',
  subtitle = paste0('SES=', round(test$SES, 3),'\np=', test$pval)
) +
  geom_vline(xintercept = ref, color = 'red') +
  theme(
    aspect.ratio = 1,
    axis.line = element_blank(),
    axis.ticks.length = unit(2, 'mm'), 
    plot.title = element_text(face = 'bold', colour = 'black', hjust = .5),
    plot.subtitle = element_text(colour = 'red', hjust = .5),
    plot.caption = element_text(colour = 'black', hjust = .5),
    panel.grid.major = element_line(linewidth = .5, color = 'grey88'),
    panel.background = element_blank(),
    panel.border = element_rect(linewidth = .8, color = 'black', fill = 'transparent')
  )

ggsave('results_taxa/select_taxa.species.13.pdf', width = 4, height = 4)

#### 20260521 random test genus 79/46/13 ####
files <- list.files('results_taxa/random_genus/', pattern = '13.results')

random <- map_dfr(
  files, ~ 
    read.delim(
      paste0('results_taxa/random_genus/', .x), 
      header = F, col.names = c('rep', 'auc')
    )
) |> 
  group_by(rep) |> 
  summarise(mean = mean(auc))

# obs for 46 / 79
ref <- read.delim('../Figure6/CF.rf.roc.all_CF_gene.tsv') |>
  pull(auc) |>
  mean()

# obs for 13
ref <- read.delim(
  '../Figure6/roc.top_13.tsv',
) |>
  pull(auc) |>
  mean()

test <- calcu_empirical_p(ref, random$mean, 'right')

ggdensity(
  random, x = 'mean', fill = 'grey77', 
  title = '13 random genus', xlab = 'Average AUC',
  subtitle = paste0('SES=', round(test$SES, 3),'\np=', test$pval)
) +
  geom_vline(xintercept = ref, color = 'red') +
  theme(
    aspect.ratio = 1,
    axis.line = element_blank(),
    axis.ticks.length = unit(2, 'mm'), 
    plot.title = element_text(face = 'bold', colour = 'black', hjust = .5),
    plot.subtitle = element_text(colour = 'black', hjust = .5),
    plot.caption = element_text(colour = 'black', hjust = .5),
    panel.grid.major = element_line(linewidth = .5, color = 'grey88'),
    panel.background = element_blank(),
    panel.border = element_rect(linewidth = .8, color = 'black', fill = 'transparent')
  )

ggsave('results_taxa/select_taxa.genus.13.pdf', width = 4, height = 4)
