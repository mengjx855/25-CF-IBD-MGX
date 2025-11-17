#### Jinxin Meng, 20230816, 20250820, v1.0.1 ####
# 20250820: update some function.

library(dplyr)
library(pROC)
library(ggplot2)

#### plot_roc ####
# 针对一个roc对象绘制图，可输入曲线的颜色
plot_roc <- function(roc, color = '#238443', plot_se = F){
  label <- paste0('AUC: ', round(roc$auc, digits = 3), 
                  '\n(95% Cl: ', paste(round(ci.auc(roc), digits = 3)[c(1,3)], collapse = ' ~ '), ')')
  p <- ggroc(roc, legacy.axes = T, lwd = .4, color = color) +
    annotate('segment', x = 0, xend = 1, y = 0, yend = 1, color = 'black', lty = 'dashed', lwd = .2) +
    labs(x = '1 - Specificity', y = 'Sensitivity') +
    annotate('text' ,x = 0.75, y = 0.125 , label = label, size = 3) +
    scale_x_continuous(expand = c(.01, .01)) +
    scale_y_continuous(expand = c(.01, .01)) +
    theme_bw() +
    theme(axis.text = element_text(size = 10, color = 'black'),
          axis.ticks = element_line(linewidth = .5, color = 'black'),
          axis.title = element_text(size = 12, color = 'black'),
          axis.line = element_blank(),
          panel.background = element_rect(linewidth = .4, color = 'black'),
          panel.grid = element_blank(),
          aspect.ratio = 1)
  if(isTRUE(plot_se)) {
    roc_se <- ci.se(roc, specificities = seq(0, 1, 0.01), conf.level = 0.95) %>% 
      data.frame(check.names = F) %>% 
      dplyr::rename(lower = all_of('2.5%'), upper = all_of('97.5%'), median = all_of('50%')) %>% 
      rownames_to_column(var = 'spec') %>% 
      mutate(spec = as.numeric(spec))
    p <- p + geom_ribbon(data = roc_se, aes(x = 1 - spec, ymin = lower, ymax = upper), fill = color, alpha = .1)
  }
  # cat('  ggsave(file = \'roc.pdf\', width = 4, height = 4)\n')
  return(p)
} 

#### plot_roc_multiple ####
# 针对一个roc列表对象绘制图，可输入曲线的颜色，用向量
plot_roc_multiple <- function(roc_list, colors = NULL, plot_se = F, title = NULL){
  if (is.null(colors)) {
    .colors <- c('#80b1d3','#b3de69','#fdb462','#8dd3c7','#bc80bd',
                 '#fb8072','#ffed6f','#fccde5','#bebada','#ccebc5')
    .color_len <- length(.colors)
    .level_len <- length(roc_list)
    colors <- rep(.colors, times = ceiling(.level_len/.color_len))[1:.level_len]
  } 
  
  labels <- map_vec(roc_list, \(x) 
                    paste0(round(x$auc, digits = 3), ' [95%-Cl ',
                           paste(round(ci.auc(x), digits = 3)[c(1,3)], 
                                 collapse = '~'), ']')) %>% 
    paste0(names(roc_list), ' ', .)
  
  p <- ggroc(roc_list, legacy.axes = T, lwd = .4) +
    annotate('segment', x = 0, y = 0, xend = 1, yend = 1, color = 'black', lty = 'longdash', lwd = .4) +
    scale_color_manual(values = colors, breaks = names(roc_list), labels = labels) +
    scale_x_continuous(expand = c(.01, .01)) +
    scale_y_continuous(expand = c(.01, .01)) +
    labs(x = '1 - Specificity', y = 'Sensitivity', color = '') +
    theme_bw() +
    theme(axis.text = element_text(size = 12, color = '#000000'),
          axis.ticks = element_line(linewidth = .4, color = '#000000'),
          axis.ticks.length = unit(2, 'mm'),
          axis.title = element_text(size = 12, color = '#000000'),
          axis.line = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_line(linewidth = .4, linetype = 'longdash', color = 'grey88'),
          panel.grid.minor = element_blank(),
          plot.title = element_text(face = 'bold', hjust = .5, size = 12, color = '#000000'),
          legend.text = element_text(size = 10, color = '#000000'),
          aspect.ratio = 1)
  
  if (!is.null(title)) p <- p + labs(title = title)
  
  if (isTRUE(plot_se)) {
    roc_se <- lapply(roc_list, \(x) ci.se(x, specificities = seq(0, 1, 0.01), conf.level = 0.95) %>% 
                       data.frame(check.names = F) %>% 
                       rename(lower = all_of('2.5%'), upper = all_of('97.5%'), median = all_of('50%')) %>% 
                       rownames_to_column(var = 'spec') %>% 
                       mutate(spec = as.numeric(spec))) %>% 
      purrr::map2_dfr(., names(.), \(x, y) add_column(x, class = y))
    p <- p + geom_ribbon(data = roc_se, aes(x = 1 - spec, ymin = lower, ymax = upper, fill = class), alpha = .1, 
                         inherit.aes = F, show.legend = F) +
      scale_fill_manual(values = colors)
  }
  cat('  ggsave(file = \'roc_list.pdf\', width = 4, height = 4)\n')
  return(p)
}
