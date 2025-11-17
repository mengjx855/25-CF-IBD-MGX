#### Jinxin Meng, 20211029, 20250404, v.21 ####

# 2023-01-01: update function: calcu_alpha.
# 2023-12-04: update function: check_file_name was deprecated.
# 2024-03-04: fix some bug
# 2025-04-04: 修改函数的某些参数名称，plot_alpha()函数中 添加 add_ref_line 参数

library(dplyr)
library(tibble)
library(tidyr)
library(vegan)
library(ggplot2)
library(ggpubr)
source('/data/mengjx/R_proj/R_func/calcu_difference.R')

#### calcu_alpha ####
calcu_alpha <- function(profile, method = "richness", out_colnames = NULL, 
                        tree = NULL, base = exp(1)) {
  profile = t(profile)
  if (method == "richness") {
    result <- rowSums(profile > 0)
  } else if (method == "chao1") {
    result <- estimateR(ceiling(profile))[2,]
  } else if (method == "observed") {
    result <- estimateR(ceiling(profile))[1,]
  } else if (method == "ace") {
    result <- estimateR(ceiling(profile))[4,]
  } else if (method == "shannon") {
    result <- diversity(profile, index = "shannon", base = base)
  } else if (method == "simpson") {
    result <- diversity(profile, index = "simpson")
  } else if (method == "pielou") {
    result <- diversity(profile, index = "shannon", base = base) / log(estimateR(profile)[1, ], base)
  } else if (method == "gc") {
    result <- 1 - rowSums(profile == 1) / rowSums(profile)
  } else if (method == "pd" & !is.null(tree)) { 
    pd <- picante::pd(profile, tree, include.root = F) # 需要有根树
    result <- pd[,1]
    names(result) <- rownames(pd)
  }
  data <- data.frame(sample = names(result), value = result, row.names = NULL)
  if (!is.null(out_colnames) & is.character(out_colnames)) 
    colnames(data)[2] <- out_colnames
  return(data)
}

#### plot_alpha ####
plot_alpha <- function(data, group, group_level = NULL, group_color = NULL, 
                       data_rename = NULL, group_rename = NULL, xlab = '', 
                       ylab = '', title = '', aspect_ratio = 1, show_grid = T,
                       show_jitter = T, rotate_x_text = F, coord_flip = F, 
                       show_diff = T, method = 'wilcox', center_group = NULL,
                       sort_value = NULL, add_ref_line = NULL, ... ) {
  
  if (!all(c('sample','value') %in% colnames(data)) & is.null(data_rename)) 
    stop("data field (sample|value)")
  
  if (!all(c('sample','group') %in% colnames(group)) & is.null(group_rename)) 
    stop("group field (sample|group)")
  
  if (!is.null(data_rename)) 
    data <- data.frame(data, check.names = F) %>% 
      dplyr::rename(all_of(data_rename))
  
  if (!is.null(group_rename)) 
    group <- data.frame(group, check.names = F) %>% 
      dplyr::rename(all_of(group_rename))
  
  if (is.null(group_level)) 
    group_level <- unique(group$group)
  
  plot_data <- left_join(data, group, by = 'sample')
  
  if (!is.null(sort_value)) {
    if (!(sort_value %in% c('asc', 'desc'))) {
      message("sort_value need specify \'asc\', \'desc\' or NULL (default) ")
      sort_value <- NULL
    }
  }

  if (!is.null(sort_value)) {
    if (sort_value == 'asc') 
      group_level <- aggregate(value ~ group, plot_data, median) %>% 
        arrange(value) %>% 
        pull(group) %>% 
        as.character()
    
    if (sort_value == 'desc') 
      group_level <- aggregate(value ~ group, plot_data, median) %>% 
        arrange(desc(value)) %>% 
        pull(group) %>% 
        as.character()
  }
  
  if (!is.null(add_ref_line))
    if (!(add_ref_line %in% group_level)) {
      message('add_ref_line not existing.')
      add_ref_line <- NULL
    }
  
  if (is.null(group_color)) {
    .color <- c("#1f78b4","#33a02c","#e31a1c","#ff7f00","#6a3d9a","#ffff99",
                "#b15928","#a6cee3","#b2df8a","#fb9a99","#fdbf6f","#cab2d6")
    .len <- length(.color)
    group_color <- rep(.color, ceiling(length(group_level)/.len))[1:length(group_level)]
  }
  
  plot_data <- mutate(plot_data, group = factor(group, group_level))
  
  p <- ggplot(plot_data, aes(group, value, fill = group)) +
    geom_boxplot(width = .618, linewidth = .4, outlier.shape = NA, 
                 show.legend = F, ...) +
    scale_fill_manual(values = group_color) +
    labs(x = xlab, y = ylab, title = title) +
    theme_pubr() +
    theme(aspect.ratio = aspect_ratio,
          axis.ticks.length = unit(2, 'mm'), 
          plot.title = element_text(hjust = .5, size = 12, face = 'bold'))
  
  if (!is.null(add_ref_line)) 
    p <- p +
    geom_hline(yintercept = median(plot_data[plot_data$group == add_ref_line, 'value']),
               linetype = "dashed", linewidth = .4, color = '#000000')
  
  if (isTRUE(show_jitter)) 
    p <- p + 
    geom_jitter(aes(color = group), size = .7, width = .2, show.legend = F) +
    scale_color_manual(values = group_color)
  
  if (isTRUE(show_diff)) {
    comparisons <- calcu_diff(data, group, method = method, 
                              center_group = center_group,
                              group_level = group_level) %>% 
      filter(pval < 0.05) %>% 
      pull(comparison) %>% 
      strsplit(x = ., split = "_vs_")
    
    p <- p + 
      stat_compare_means(comparisons = comparisons, method = method, 
                         method.args = list(exact = F), label = "p.signif", 
                         tip.length = .02, step.increase = .05, 
                         vjust = .8, size = 3)
  }
  
  if (isTRUE(show_grid)) 
    p <- p + theme(panel.grid.major = element_line(color = "grey88", linewidth = .4))
  
  if (isTRUE(rotate_x_text)) 
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
  
  if (isTRUE(coord_flip)) 
    p <- p + coord_flip()
  
  message(paste0(" ggsave(file = \"alpha-div.pdf\", width = ", 
                 length(group_level)*0.6 ,", height = 4.5)"))
  return(p)
}
