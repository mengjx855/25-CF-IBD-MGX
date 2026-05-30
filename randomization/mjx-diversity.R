#### Jinxin Meng, 20260404, 20260530 ####
setwd('/data/mengjx/project/10.20250623_IBD_BAC_CF_Landscape/github/randomization/')
pacman::p_load(tidyverse, ggpubr, rstatix, openxlsx)
source('scripts/calcu_difference.R')

#### gene alphaDiv ####
files <- list.files('random_div/', pattern = 'gene.alphaDiv.tsv')

random <- map_dfr(
  files[-1], ~ 
    read.delim(paste0('random_div/', .x)) %>%
    add_column(rep = str_split_i(.x, '\\.', 1), .before = 1) 
  )

ref <- read.delim('random_div/rep000_seed0000.gene.alphaDiv.tsv') |> 
  group_by(type) |> 
  summarise(mean = mean(log2FC)) |> 
  pull(mean, name = type)

test <- map(
  c('shannon', 'richness'), ~ {
    .random <- filter(random, type == .x) %>% 
      group_by(rep) %>% 
      summarise(mean = mean(log2FC)) %>% 
      pull(mean)
    calcu_empirical_p(ref[[.x]], .random, 'right')
  }
) %>% 
  set_names(
    c('shannon', 'richness')
  )

plots <- map(
  c('shannon', 'richness'), ~ 
    filter(random, type == .x) %>% 
    group_by(rep) %>% 
    summarise(mean = mean(log2FC)) %>% 
    ggdensity(
      x = 'mean', fill = 'grey77', 
      xlab = paste0(
        'Average ', str_to_title(.x), ' coefficient'
      ),
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
)

#### gene betaDiv ####
files <- list.files('random_div/', pattern = 'gene.betaDiv.tsv')

random <- map_dfr(
  files[-1], ~ 
    read.delim(paste0('random_div/', .x)) %>%
    add_column(rep = str_split_i(.x, '\\.', 1), .before = 1)
)

ref <- read.delim('random_div/rep000_seed0000.gene.betaDiv.tsv') |> 
  pull(r2adj) |> 
  mean()

test <- calcu_empirical_p(
  ref, 
  random %>%
    dplyr::group_by(rep) %>%
    dplyr::summarise(mean = mean(r2adj), .groups = "drop") %>%
    dplyr::pull(mean),
  alternative = 'right'
)

plots[[3]] <- group_by(random, rep) %>% 
  summarise(mean = mean(r2adj)) %>% 
  ggdensity(
    x = 'mean', fill = 'grey77', xlab = 'Average adjusted R2',
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
      panel.grid.major = element_line(linewidth = .5, color = 'grey88'),
      panel.background = element_blank(),
      panel.border = element_rect(linewidth = .8, color = 'black', fill = 'transparent')
    )

cowplot::plot_grid(plotlist = plots, nrow = 1, align = 'hv')
ggsave('results/3.Figure4.div.gene_level.pdf', width = 12, height = 4)

#### fam alphaDiv ####
files <- list.files('random_div/', pattern = 'fam.alphaDiv.tsv')

random <- map_dfr(
  files[-1], ~ 
    read.delim(paste0('random_div/', .x)) %>%
    add_column(rep = str_split_i(.x, '\\.', 1), .before = 1) 
  )

ref <- read.delim('random_div/rep000_seed0000.fam.alphaDiv.tsv') |> 
  group_by(type) |> 
  summarise(mean = mean(log2FC)) |> 
  pull(mean, name = type)

test <- map(
  c('shannon', 'richness'), ~ {
    .random <- filter(random, type == .x) %>% 
      group_by(rep) %>% 
      summarise(mean = mean(log2FC)) %>% 
      pull(mean)
    calcu_empirical_p(ref[[.x]], .random, 'right')
  }
) %>% 
  set_names(
    c('shannon', 'richness')
  )

plots <- map(
  c('shannon', 'richness'), ~ 
    filter(random, type == .x) %>% 
    group_by(rep) %>% 
    summarise(mean = mean(log2FC)) %>% 
    ggdensity(
      x = 'mean', fill = 'grey77', 
      xlab = paste0(
        'Average ', str_to_title(.x), ' coefficient'
      ),
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
)

#### fam betaDiv ####
files <- list.files('random_div/', pattern = 'fam.betaDiv.tsv')

random <- map_dfr(
  files[-1], ~ 
    read.delim(paste0('random_div/', .x)) %>%
    add_column(rep = str_split_i(.x, '\\.', 1), .before = 1)
)

ref <- read.delim('random_div/rep000_seed0000.fam.betaDiv.tsv') |> 
  pull(r2adj) |> 
  mean()

test <- calcu_empirical_p(
  ref, 
  random %>%
    dplyr::group_by(rep) %>%
    dplyr::summarise(mean = mean(r2adj), .groups = "drop") %>%
    dplyr::pull(mean),
  alternative = 'right'
)

plots[[3]] <- group_by(random, rep) %>% 
  summarise(mean = mean(r2adj)) %>% 
  ggdensity(
    x = 'mean', fill = 'grey77', xlab = 'Average adjusted R2',
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
      panel.grid.major = element_line(linewidth = .5, color = 'grey88'),
      panel.background = element_blank(),
      panel.border = element_rect(linewidth = .8, color = 'black', fill = 'transparent')
    )

cowplot::plot_grid(plotlist = plots, nrow = 1, align = 'hv')
ggsave('results/3.Figure4.div.fam_level.pdf', width = 12, height = 4)
