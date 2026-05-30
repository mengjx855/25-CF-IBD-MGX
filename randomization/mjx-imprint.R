#### Jin-Xin Meng, 20260404, 20260530 ####
setwd('/data/mengjx/project/10.20250623_IBD_BAC_CF_Landscape/github/randomization/')
pacman::p_load(tidyverse, ggpubr, rstatix, openxlsx)
source('scripts/calcu_difference.R')

#### Jaccard-binary-based PCoA ####
taxa_level <- c('Domain', 'Phylum','Class','Order','Family','Genus')

files <- list.files('random_imprint/', pattern = 'pcoa.tsv')

random <- map_dfr(
  files[-1], ~ 
    read.delim(paste0('random_imprint/', .x)) |> 
    add_column(rep = str_split_i(.x, '\\.', 1), .before = 1)
  )

ref <- read.delim('random_imprint/rep000_seed0000.pcoa.tsv') |> 
  pull(r2adj, name = name)

test <- map(
  taxa_level, ~ {
    .ref <- ref[[.x]]
    .random <- filter(random, name == .x) |> pull(r2adj)
    calcu_empirical_p(.ref, .random, 'left')
  }
) |> 
  set_names(taxa_level)

plots <- map(
  taxa_level, ~ 
    filter(random, name == .x) |> 
    ggdensity(
      x = 'r2adj', fill = 'grey77', title = paste0(.x, '-level'),
      xlab = 'Adonis adjusted R2',
      subtitle = paste0(
        'SES=', round(test[[.x]][['SES']], 3), 
        '\np=', round(test[[.x]][['pval']], 3)
        )
      ) +
    geom_vline(xintercept = ref[[.x]], color = 'red') +
    theme(
      aspect.ratio = 1,
      axis.line = element_blank(),
      axis.ticks.length = unit(2, 'mm'), 
      plot.title = element_text(face = 'bold', colour = 'black', hjust = .5),
      plot.subtitle = element_text(colour = 'red', hjust = .5),
      panel.grid.major = element_line(linewidth = .5, color = 'grey88'),
      panel.background = element_blank(),
      panel.border = element_rect(linewidth = .8, color = 'black', fill = 'transparent')
      )
  )

cowplot::plot_grid(plotlist = plots, nrow = 2, ncol = 4, align = 'v')
ggsave('results/4.Figure2.imprint.pcoa.pdf', width = 16, height = 8)

#### Jaccard-dist mantel #####
files <- list.files('random_imprint/', pattern = 'mantel.tsv')

random <- map_dfr(
  files[-1], ~ 
    read.delim(paste0('random_imprint/', .x)) |> 
    add_column(rep = str_split_i(.x, '\\.', 1), .before = 1) )

ref <- read.delim('random_imprint/rep000_seed0000.mantel.tsv') |> 
  pull(statistic, name = domain)

test <- map(
  c('bacteria', 'archaea'), ~ {
    .random <- filter(random, domain == .x & !is.na(statistic)) |> 
      pull(statistic)
    calcu_empirical_p(ref[[.x]], .random, 'left')
  }
) |> 
  set_names(c('bacteria', 'archaea'))

plots <- map(
  c('bacteria', 'archaea'), ~ {
    filter(random, domain == .x) |> 
      ggdensity(
        x = 'statistic', fill = 'grey77', xlab = 'mantel r',
        title = paste0(str_to_title(.x), ' lineage'), 
        subtitle = paste0(
          'SES=', round(test[[.x]][['SES']], 3),
          '\np=', round(test[[.x]][['pval']], 3)
        )
      ) +
      geom_vline(xintercept = ref[[.x]], color = 'red') +
      theme(
        aspect.ratio = 1,
        axis.line = element_blank(),
        axis.ticks.length = unit(2, 'mm'), 
        plot.title = element_text(face = 'bold', colour = 'black', hjust = .5),
        plot.subtitle = element_text(colour = 'red', hjust = .5),
        panel.grid.major = element_line(linewidth = .5, color = 'grey88'),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = .8, color = 'black', fill = 'transparent')
      )
  }
)

cowplot::plot_grid(plotlist = plots, nrow = 1, align = 'v')
ggsave('results/4.Figure2.imprint.mantel.pdf', width = 8, height = 4)
