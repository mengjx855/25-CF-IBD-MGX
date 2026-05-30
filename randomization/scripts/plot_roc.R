#### Jinxin Meng, 20230816, 20260519, v0.1.1 ####

# 20250820: update some function.
# 20260519: simplify code, add helper functions, fix CI label and ribbon order.

library(dplyr)
library(tibble)
library(ggplot2)

#### helper: ROC plot theme ####
# ROC 图主题
# style:
#   default : 黑框、无网格，适合常规论文图
#   grid    : 黑框、带浅灰虚线网格
#   classic : 只有坐标轴，无外框，比较简洁
#   minimal : 极简风格，适合展示型图
#   nature  : 类 Nature 风格，坐标轴更突出
#   lancet  : 类 Lancet 风格，黑框较明显
theme_roc <- function(base_size = 12, show_grid = NULL,
                      style = c('default', 'grid', 'classic', 
                                'minimal', 'nature', 'lancet')) {
  
  style <- match.arg(style)
  
  if (is.null(show_grid)) {
    show_grid <- style %in% c('grid', 'minimal')
  }
  
  p <- switch(
    style,
    default = theme_bw(),
    grid = theme_bw(),
    classic = theme_classic(),
    minimal = theme_minimal(),
    nature = theme_classic(),
    lancet = theme_bw()
  )
  
  p <- p +
    theme(
      axis.text = element_text(size = base_size - 1, color = 'black'),
      axis.title = element_text(size = base_size, color = 'black'),
      axis.ticks = element_line(linewidth = .4, color = 'black'),
      axis.ticks.length = unit(2, 'mm'),
      axis.line = element_blank(),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(
        face = 'bold', hjust = .5, size = base_size, color = 'black'
      ),
      plot.subtitle = element_text(
        hjust = .5, size = base_size - 1, color = 'black'
      ),
      legend.text = element_text(size = base_size - 1, color = 'black'),
      legend.title = element_text(size = base_size - 1, color = 'black'),
      aspect.ratio = 1
    )
  
  ## 不同主题的细节
  if (style %in% c('default', 'grid')) {
    p <- p +
      theme(
        panel.border = element_rect(linewidth = .4, color = 'black', fill = NA)
      )
  }
  
  if (style == 'classic') {
    p <- p +
      theme(
        axis.line = element_line(linewidth = .4, color = 'black'),
        panel.border = element_blank()
      )
  }
  
  if (style == 'minimal') {
    p <- p +
      theme(
        panel.border = element_blank(),
        axis.line = element_line(linewidth = .3, color = 'black')
      )
  }
  
  if (style == 'nature') {
    p <- p +
      theme(
        axis.line = element_line(linewidth = .5, color = 'black'),
        axis.ticks = element_line(linewidth = .5, color = 'black'),
        panel.border = element_blank(),
        legend.key = element_blank()
      )
  }
  
  if (style == 'lancet') {
    p <- p +
      theme(
        panel.border = element_rect(linewidth = .6, color = 'black', fill = NA),
        axis.ticks = element_line(linewidth = .5, color = 'black'),
        legend.key = element_blank()
      )
  }
  
  ## 是否添加网格
  if (isTRUE(show_grid)) {
    p <- p +
      theme(
        panel.grid.major = element_line(
          linewidth = .35, linetype = 'longdash', color = 'grey88'
        ),
        panel.grid.minor = element_blank()
      )
  }
  
  return(p)
}


#### helper: AUC label ####
roc_auc_label <- function(roc, digits = 3, prefix = 'AUC') {
  
  auc_value <- as.numeric(pROC::auc(roc))
  
  auc_ci <- tryCatch(
    as.numeric(pROC::ci.auc(roc)),
    error = function(e) c(NA_real_, NA_real_, NA_real_)
  )
  
  paste0(
    prefix, ': ', round(auc_value, digits),
    '\n(95% CI: ',
    paste(round(auc_ci[c(1, 3)], digits), collapse = ' ~ '),
    ')'
  )
}

#### helper: SE confidence ribbon data ####
roc_se_data <- function(roc, by = 0.01, conf.level = 0.95) {
  
  roc_se <- tryCatch(
    pROC::ci.se(
      roc,
      specificities = seq(0, 1, by),
      conf.level = conf.level
    ),
    error = function(e) NULL
  )
  
  if (is.null(roc_se)) return(NULL)
  
  roc_se <- data.frame(roc_se, check.names = FALSE) |>
    rename(
      lower = all_of('2.5%'),
      median = all_of('50%'),
      upper = all_of('97.5%')
    ) |>
    rownames_to_column('spec') |>
    mutate(spec = as.numeric(spec))
  
  return(roc_se)
}


#### plot_roc ####
# 针对一个 pROC::roc 对象绘制 ROC 曲线。
# roc: pROC::roc() 返回的对象。
# color: ROC 曲线颜色。
# plot_se: 是否添加 sensitivity 的 95% CI ribbon。
# label_pos: AUC 标签位置，格式为 c(x, y)。

plot_roc <- function(roc, color = '#238443', plot_se = FALSE,
                     label_pos = c(0.75, 0.125), title = NULL, 
                     subtitle = NULL, linewidth = .6, se_alpha = .12,
                     theme_style = 'default', base_size = 12,
                     show_grid = FALSE, quiet = FALSE) {
  
  label <- roc_auc_label(roc)
  
  ## 基础 ROC 曲线
  p <- pROC::ggroc(
    roc,
    legacy.axes = TRUE,
    linewidth = linewidth,
    color = color
  ) +
    annotate(
      'segment',
      x = 0, xend = 1, y = 0, yend = 1,
      color = 'black',
      linetype = 'dashed',
      linewidth = .25
    ) +
    annotate(
      'text',
      x = label_pos[1],
      y = label_pos[2],
      label = label,
      size = 3
    ) +
    scale_x_continuous(expand = c(.01, .01)) +
    scale_y_continuous(expand = c(.01, .01)) +
    labs(
      x = '1 - Specificity',
      y = 'Sensitivity',
      title = title,
      subtitle = subtitle
    ) +
    theme_roc(
      base_size = base_size, 
      show_grid = show_grid, 
      style = theme_style
    )
  
  ## 添加 sensitivity 置信区间
  ## 注意这里用 p$layers 调整顺序，让 ribbon 位于 ROC 曲线下方。
  if (isTRUE(plot_se)) {
    
    roc_se <- roc_se_data(roc)
    
    if (!is.null(roc_se)) {
      ribbon_layer <- geom_ribbon(
        data = roc_se,
        aes(x = 1 - spec, ymin = lower, ymax = upper),
        fill = color,
        alpha = se_alpha,
        inherit.aes = FALSE
      )
      
      p <- p + ribbon_layer
      p$layers <- c(
        p$layers[length(p$layers)],
        p$layers[-length(p$layers)]
      )
    }
  }
  
  if (isFALSE(quiet)) {
    message("  ggsave(file = 'roc.pdf', width = 4, height = 4)")
  }
  
  return(p)
  
}

#### plot_roc_multiple ####
# 针对多个 pROC::roc 对象绘制 ROC 曲线。
# roc_list: roc 对象列表，建议命名。
# colors: 曲线颜色向量；如果为空则自动生成。
# plot_se: 是否添加 sensitivity 95% CI ribbon
plot_roc_multiple <- function(roc_list, colors = NULL, plot_se = FALSE,
                              title = NULL, subtitle = NULL,
                              linewidth = .6, se_alpha = .10,
                              theme_style = 'default', base_size = 12,
                              show_grid = TRUE, quiet = FALSE) {
  
  if (is.null(names(roc_list)) || any(names(roc_list) == '')) {
    names(roc_list) <- paste0('ROC_', seq_along(roc_list))
  }
  
  if (is.null(colors)) {
    base_colors <- c(
      '#80b1d3', '#b3de69', '#fdb462', '#8dd3c7', '#bc80bd',
      '#fb8072', '#ffed6f', '#fccde5', '#bebada', '#ccebc5'
    )
    colors <- rep(base_colors, length.out = length(roc_list))
  }
  
  names(colors) <- names(roc_list)
  
  labels <- purrr::map_chr(
    roc_list, \(x) {
      auc_value <- as.numeric(pROC::auc(x))
      auc_ci <- tryCatch(
        as.numeric(pROC::ci.auc(x)),
        error = function(e) c(NA_real_, NA_real_, NA_real_)
      )
      
      paste0(
        round(auc_value, 3),
        ' [95% CI ',
        paste(round(auc_ci[c(1, 3)], 3), collapse = '~'),
        ']'
      )
    }
  )
  
  labels <- paste0(names(roc_list), ' ', labels)
  
  p <- pROC::ggroc(
    roc_list,
    legacy.axes = TRUE,
    linewidth = linewidth
  ) +
    annotate(
      'segment',
      x = 0, y = 0, xend = 1, yend = 1,
      color = 'black',
      linetype = 'longdash',
      linewidth = .35
    ) +
    scale_color_manual(
      values = colors,
      breaks = names(roc_list),
      labels = labels
    ) +
    scale_x_continuous(expand = c(.01, .01)) +
    scale_y_continuous(expand = c(.01, .01)) +
    labs(
      x = '1 - Specificity',
      y = 'Sensitivity',
      color = '',
      title = title,
      subtitle = subtitle
    ) +
    theme_roc(
      base_size = base_size, 
      show_grid = show_grid, 
      style = theme_style
    )
  
  ## 添加多个 ROC 的 sensitivity 置信区间
  if (isTRUE(plot_se)) {
    
    roc_se <- purrr::imap_dfr(
      roc_list,
      \(x, y) {
        tmp <- roc_se_data(x)
        if (is.null(tmp)) return(NULL)
        dplyr::mutate(tmp, class = y)
      }
    )
    
    if (nrow(roc_se) > 0) {
      p <- p +
        geom_ribbon(
          data = roc_se,
          aes(x = 1 - spec, ymin = lower, ymax = upper, fill = class),
          alpha = se_alpha,
          inherit.aes = FALSE,
          show.legend = FALSE
        ) +
        scale_fill_manual(values = colors)
      
      ## 让 ribbon 位于曲线下面
      p$layers <- c(
        p$layers[length(p$layers)],
        p$layers[-length(p$layers)]
      )
    }
  }
  
  if (isFALSE(quiet)) {
    message("  ggsave(file = 'roc_list.pdf', width = 4.5, height = 4)")
  }
  
  return(p)
}
