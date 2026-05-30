#### Jinxin Meng, 20260404, 20260530 ####
setwd('/data/mengjx/project/10.20250623_IBD_BAC_CF_Landscape/github/randomization/')
pacman::p_load(tidyverse, ggpubr, rstatix, openxlsx)
source('scripts/calcu_difference.R')

#### procrustes ####
files <- list.files('random_procrustes/', pattern = 'procrustes.tsv')

random <- map_dfr(
  files[-1], ~ 
    read.delim(paste0('random_procrustes/', .x)) %>%
    add_column(rep = str_split_i(.x, '\\.', 1), .before = 1)
)

ref <- read.delim('random_procrustes/rep000_seed0000.fam.procrustes.tsv') |> 
  pull(M2) |> mean()

test <- calcu_empirical_p(
  ref, 
  random %>%
    dplyr::group_by(rep) %>%
    dplyr::summarise(mean = mean(M2), .groups = "drop") %>%
    dplyr::pull(mean),
  alternative = 'right'
)

group_by(random, rep) %>% 
  summarise(mean = mean(M2)) %>% 
  ggdensity(
    x = 'mean', fill = 'grey77', title = 'Procrustes test',
    xlab = 'Average Procrustes M2',
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
ggsave('results/5.Figure5.procrustes.pdf', width = 4, height = 4)
