#### Jinxin Meng, 20260404, 20260523 ####
setwd('/data/mengjx/project/10.20250623_IBD_BAC_CF_Landscape/github/randomization/')
pacman::p_load(tidyverse, ggpubr, rstatix, openxlsx)
source('scripts/calcu_difference.R')

#### alphaDiv ####
files <- list.files('random_div_RNASeq/', pattern = 'gene.alphaDiv.tsv')

random <- map_dfr(
  files[-1], ~ 
    read.delim(paste0('random_div_RNASeq/', .x)) %>%
    add_column(rep = str_split_i(.x, '\\.', 1), .before = 1) 
  )

ref_file <- 'random_div_RNASeq/rep000_seed0000.gene.alphaDiv.tsv'

tests <- tibble::tibble()

map(
  c('shannon', 'richness'), \(x) {
    
    map(
      c('SchirmerM_2018', 'LloydPriceJ_2019'), \(y) {
        
        map(
          c('CD_vs_UC', 'CD_vs_HC', 'UC_vs_HC'), \(z) {
            
            random_f <- filter(
              random, 
              type == x,
              grepl(y, name),
              grepl(z, name)
            )
            
            ref <- read.delim(ref_file) |> 
              filter(
                type == x,
                grepl(y, name),
                grepl(z, name)
              ) |> 
              pull(log2FC)
            
            test <- calcu_empirical_p(ref, random_f$log2FC, 'auto', simplify = T) |> 
              add_column(index = x, proj = y, comparisons = z, .before = 1)
            
            tests <<- bind_rows(tests, test)
            
            ggdensity(
              random_f, x = 'log2FC', fill = 'grey77', 
              title = paste0(c(x, y, z), collapse = ' '), 
              subtitle = paste0('SES=', round(test$SES, 3), '\np=', round(test$pval, 3))) +
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
          }
        ) |> 
          cowplot::plot_grid(plotlist = _, align = 'hv', nrow = 1)
      }
    ) |> 
      cowplot::plot_grid(plotlist = _, align = 'hv', ncol = 1)
  }
) |> 
  cowplot::plot_grid(plotlist = _, align = 'hv', ncol = 1)

write.xlsx(tests, 'results_RNASeq/6.div.gene.xlsx')

ggsave('results_RNASeq/6.div.gene.pdf', width = 12, height = 16)

tests |> 
  ggplot(aes(comparisons, index)) +
  geom_tile(aes(fill = SES)) +
  facet_grid(rows = vars(proj)) +
  scale_fill_gradientn(
    colors = c('#4F6F9D',"#94A8D3","#FFFFFF","#FEEF6B",'#F7DF5A'),
    values = scales::rescale(c(-2.4, -1.2, 0, 1.2, 2.4)),
    limits = c(-2.4, 2.4)
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = '', y = '', fill = 'SES') +
  geom_text(aes(label = scales::number_format(accuracy = 0.001)(pval)), size = 4) +
  coord_fixed() +
  theme(
    axis.ticks = element_blank(),
    axis.text = element_text(color = '#000000', size = 10),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    plot.title = element_text(hjust = .5, face = 'bold'),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(linewidth = .8, color = 'black', fill = NA),
    strip.background = element_rect(linewidth = .8, color = 'black', fill = NA)
  )
ggsave('results_RNASeq/6.Figure7.div.gene.tile.pdf', width = 5, height = 6)
