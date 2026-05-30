#### Jinxin Meng, 20260404, 20260525 ####
setwd('/data/mengjx/project/10.20250623_IBD_BAC_CF_Landscape/github/randomization/')
pacman::p_load(tidyverse, ggpubr, rstatix, openxlsx)
source('scripts/calcu_difference.R')

#### adonis-by-fam ####
files <- list.files('random_adonis/', pattern = 'explain_by_CF')

random <- map_dfr(
  files[-1], ~ 
    read.delim(paste0('random_adonis/', .x)) %>%
    add_column(rep = str_split_i(.x, '\\.', 1), .before = 1) 
  )

ref <- read.delim('random_adonis/rep000_seed0000.fam.explain_by_CF.tsv') |> 
  pull(r2adj) |> 
  mean()

test <- calcu_empirical_p(
  ref, 
  group_by(random, rep) |> 
    summarise(mean = mean(r2adj)) %>% 
    pull(mean),
  alternative = 'left'
  )

p1 <- group_by(random, rep) %>% 
  summarise(mean = mean(r2adj)) %>% 
  ggdensity(
    x = 'mean', fill = 'grey77', xlab = 'Average adjusted R2',
    title = 'Species-level profile variations\nexplained by CF−family',
    subtitle = paste0(
      'SES=', round(test[['SES']], 3),
      '\np=', round(test[['pval']], 3))
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

#### adonis-by-species ####
files <- list.files('random_adonis/', pattern = 'explain_by_species')

random <- map_dfr(
  files[-1], ~ 
    read.delim(paste0('random_adonis/', .x)) %>%
    add_column(rep = str_split_i(.x, '\\.', 1), .before = 1) 
)

ref <- read.delim('random_adonis/rep000_seed0000.fam.explain_by_species.tsv') |> 
  pull(r2adj) |> 
  mean()

test <- calcu_empirical_p(
  ref, 
  group_by(random, rep) |> 
    summarise(mean = mean(r2adj)) %>% 
    pull(mean),
  'left'
)

p2 <- group_by(random, rep) %>% 
  summarise(mean = mean(r2adj)) %>% 
  ggdensity(
    x = 'mean', fill = 'grey77', xlab = 'Average adjusted R2',
    title = 'CF−family profile variations\nexplained by species',
    subtitle = paste0(
      'SES=', round(test[['SES']], 3),
      '\np=', round(test[['pval']], 3)
    )
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

cowplot::plot_grid(p1, p2, nrow = 1, align = 'hv')
ggsave('results/5.Figure5.adonis.pdf', width = 8, height = 4.2)
