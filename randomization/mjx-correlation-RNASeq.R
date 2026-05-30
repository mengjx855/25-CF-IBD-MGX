#### Jinxin Meng, 20260404, 20260523 ####
setwd('/data/mengjx/project/10.20250623_IBD_BAC_CF_Landscape/github/randomization/')
pacman::p_load(tidyverse, ggpubr, rstatix, openxlsx)
source('scripts/calcu_difference.R')

#### correlation ####
files <- list.files('random_correlation_RNASeq/', pattern = '.tsv')

random <- map_dfr(
  files[-1], ~ 
    read.delim(paste0('random_correlation_RNASeq/', .x)) %>%
    add_column(rep = str_split_i(.x, '\\.', 1), .before = 1) 
  )

ref <- read.delim('random_correlation_RNASeq/rep000_seed0000.DNA_vs_RNA.corr.tsv') |> 
  pull(r, name = proj)

test <- map(
  c('LloydPriceJ_2019', 'SchirmerM_2018'), ~
    calcu_empirical_p(
      ref[[.x]], 
      random[random[['proj']] == .x, 'r'], 
      'right'
    )
) |> 
  set_names(
    c('LloydPriceJ_2019', 'SchirmerM_2018')
  )

plots <- map(
  c('LloydPriceJ_2019', 'SchirmerM_2018'), ~ {
    ggdensity(
      filter(random, proj == .x),
      x = 'r', fill = 'grey77', 
      title = paste0(str_to_title(.x), ' spearman R'),
      subtitle = paste0(
        'SES=', round(test[[.x]][['SES']], 3), 
        '\np=', round(test[[.x]][['pval']], 3))
    ) +
      geom_vline(xintercept = ref[[.x]], color = 'red') +
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
  }
)

cowplot::plot_grid(plotlist = plots, nrow = 1, align = 'hv')
ggsave('results_RNASeq/2.Figure7.correlation.pdf', width = 8, height = 4)
