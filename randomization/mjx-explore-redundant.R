#### Jinxin Meng, 20260404, 20260530 ####
setwd('/data/mengjx/project/10.20250623_IBD_BAC_CF_Landscape/github/randomization/')
pacman::p_load(tidyverse, ggpubr, rstatix, openxlsx)
source('scripts/calcu_difference.R')
source('scripts/model_randomforest.R')

#### evaluate feature signal and redundancy ####
# 目的：
#   用于比较不同 family/profile 表的疾病信息量和特征冗余程度。
# 
# 核心思想：
#   1. feature-label signal:
#      每个 feature 单独区分 HC vs case 的 AUC 信号。
#      signal = 2 * abs(AUC - 0.5)
#      signal = 0 表示无区分能力；
#      signal = 1 表示完美区分。
#
#   2. feature-feature redundancy:
#      feature 之间丰度 profile 的相关性。
#      mean_abs_cor 越高，说明不同 feature 之间越冗余。
#
#   3. signal / redundancy ratio:
#      平均疾病信息量 / 平均冗余程度。
#      越高表示该 profile 在单位冗余下包含更多疾病相关信息。
#
# 输入要求：
#   profile 文件：行为 feature，列为 sample，第一列名为 name。
#   group 文件：至少包含 sample 和 group 两列。

## 1. 读取 profile 并按照当前 group 样本过滤
read_profile_for_group <- function(file, group, sample_col = 'sample') {
  
  ## 读取 feature × sample 矩阵
  mat <- data.table::fread(file, data.table = FALSE, check.names = FALSE) |>
    tibble::column_to_rownames('name')
  
  ## 只保留 group 表中存在、且 profile 中也存在的样本
  samples <- group[[sample_col]]
  samples <- samples[samples %in% colnames(mat)]
  
  mat <- mat[, samples, drop = FALSE]
  mat <- as.matrix(mat)
  storage.mode(mat) <- 'numeric'
  
  mat[is.na(mat)] <- 0
  
  return(mat)
}

## 2. 构建 HC vs case 的二分类标签
make_binary_label <- function(group, sample_col = 'sample', status_col = 'group',
                              control = 'HC', case_levels = NULL) {
  
  group <- group |>
    dplyr::filter(
      !is.na(.data[[sample_col]]),
      !is.na(.data[[status_col]])
    )
  
  ## 如果指定 case_levels，则只保留 control 和指定 case
  ## 如果不指定，则默认所有非 control 样本都作为 case
  if (!is.null(case_levels)) {
    group <- group |>
      dplyr::filter(.data[[status_col]] %in% c(control, case_levels))
  }
  
  label <- ifelse(group[[status_col]] == control, 0, 1)
  names(label) <- group[[sample_col]]
  
  return(label)
}

## 3. 计算每个 feature 的单特征 AUC signal
calc_auc_signal <- function(mat, label) {
  
  ## 对齐 profile 和 label 的样本
  common <- intersect(colnames(mat), names(label))
  mat <- mat[, common, drop = FALSE]
  y <- label[common]
  
  res <- lapply(rownames(mat), \(f) {
    x <- as.numeric(mat[f, ])
    ok <- is.finite(x) & !is.na(y)
    x <- x[ok]
    yy <- y[ok]
    
    ## 如果该 feature 没有变化，或者只有一个分组，则无法计算 AUC
    if (length(unique(yy)) < 2 || length(unique(x)) < 2) {
      return(data.frame(
        feature = f,
        auc = NA_real_,
        signal = NA_real_
      ))
    }
    
    ## direction = "auto" 表示自动调整方向，使 AUC 通常 >= 0.5
    roc_obj <- tryCatch(
      pROC::roc(
        response = yy,
        predictor = x,
        levels = c(0, 1),
        direction = 'auto',
        quiet = TRUE
      ),
      error = function(e) NULL
    )
    
    auc_val <- if (is.null(roc_obj)) NA_real_ else as.numeric(pROC::auc(roc_obj))
    
    data.frame(
      feature = f,
      auc = auc_val,
      signal = 2 * abs(auc_val - 0.5)
    )
  })
  
  dplyr::bind_rows(res)
}

## 4. 计算 feature-feature redundancy
calc_redundancy <- function(mat, method = 'spearman') {
  
  ## feature 少于 2 个时无法计算相关性矩阵
  if (nrow(mat) < 2) {
    return(data.frame(
      mean_abs_cor = NA_real_,
      median_abs_cor = NA_real_,
      q90_abs_cor = NA_real_,
      prop_abs_cor_gt_0.7 = NA_real_
    ))
  }
  
  ## cor(t(mat)) 表示计算 feature 与 feature 之间的相关性
  cormat <- suppressWarnings(
    cor(t(mat), method = method, use = 'pairwise.complete.obs')
  )
  
  ## 只取上三角，避免重复计算和对角线自相关
  vals <- abs(cormat[upper.tri(cormat)])
  vals <- vals[is.finite(vals)]
  
  if (length(vals) == 0) {
    return(data.frame(
      mean_abs_cor = NA_real_,
      median_abs_cor = NA_real_,
      q90_abs_cor = NA_real_,
      prop_abs_cor_gt_0.7 = NA_real_
    ))
  }
  
  data.frame(
    mean_abs_cor = mean(vals, na.rm = TRUE),
    median_abs_cor = median(vals, na.rm = TRUE),
    q90_abs_cor = as.numeric(quantile(vals, 0.9, na.rm = TRUE)),
    prop_abs_cor_gt_0.7 = mean(vals > 0.7, na.rm = TRUE)
  )
}

## 5. 对一个 profile 文件计算 signal 和 redundancy
eval_one_profile_file <- function(file, rep_name, group, label,
                                  sample_col = 'sample',
                                  log_transform = TRUE,
                                  min_sd = 0) {
  
  ## 读取并对齐样本
  mat <- read_profile_for_group(
    file = file,
    group = group,
    sample_col = sample_col
  )
  
  common <- intersect(colnames(mat), names(label))
  mat <- mat[, common, drop = FALSE]
  label <- label[common]
  
  ## 可选 log 转换
  ## 对 TPM / abundance profile，log10(x + 1) 可以降低高丰度 feature 的主导作用
  if (isTRUE(log_transform)) {
    mat <- log10(mat + 1)
  }
  
  ## 去掉零方差或近零方差 feature
  ## 这些 feature 无法提供疾病区分信息，也会影响相关性计算
  keep <- apply(mat, 1, stats::sd, na.rm = TRUE) > min_sd
  mat <- mat[keep, , drop = FALSE]
  
  ## 如果过滤后 feature 太少，直接返回空结果
  if (nrow(mat) < 2) {
    return(tibble::tibble(
      rep = rep_name,
      n_sample = length(common),
      n_feature = nrow(mat),
      mean_signal = NA_real_,
      median_signal = NA_real_,
      q90_signal = NA_real_,
      n_signal_gt_0.2 = NA_integer_,
      mean_abs_cor = NA_real_,
      median_abs_cor = NA_real_,
      q90_abs_cor = NA_real_,
      prop_abs_cor_gt_0.7 = NA_real_,
      signal_redundancy_ratio = NA_real_
    ))
  }
  
  ## 计算 feature-label signal
  auc_df <- calc_auc_signal(mat, label)
  
  ## 计算 feature-feature redundancy
  red_df <- calc_redundancy(mat, method = 'spearman')
  
  mean_signal <- mean(auc_df$signal, na.rm = TRUE)
  
  ## 汇总输出
  tibble::tibble(
    rep = rep_name,
    n_sample = length(common),
    n_feature = nrow(mat),
    
    ## 疾病区分信号
    mean_signal = mean_signal,
    median_signal = median(auc_df$signal, na.rm = TRUE),
    q90_signal = as.numeric(quantile(auc_df$signal, 0.9, na.rm = TRUE)),
    n_signal_gt_0.2 = sum(auc_df$signal > 0.2, na.rm = TRUE),
    
    ## feature 冗余
    mean_abs_cor = red_df$mean_abs_cor,
    median_abs_cor = red_df$median_abs_cor,
    q90_abs_cor = red_df$q90_abs_cor,
    prop_abs_cor_gt_0.7 = red_df$prop_abs_cor_gt_0.7,
    
    ## 单位冗余下的信息量
    signal_redundancy_ratio = mean_signal / (red_df$mean_abs_cor + 1e-6)
  )
}

run_proj <- function(proj_name, metadata_base) {
  group <- read.delim(
    paste0(
      metadata_base, proj_name, '.sample_group.gz')
    ) %>% 
    dplyr::select(sample, group) %>% 
    mutate(group = ifelse(group == 'HC', 'HC', 'IBD'))
  
  ## 读取所有 profile 文件
  all_files <- list.files("random_profile", pattern = "fam.tpm", full.names = TRUE)
  
  ref_file <- all_files[1]
  random_files <- all_files[-1]
  random_names <- map_chr(
    basename(random_files),
    ~ str_split_1(.x, "\\.")[1]
  )
  
  ## 构建标签：HC vs Other
  label <- make_binary_label(
    group = group,
    sample_col = "sample",
    status_col = "group",     # 这里根据你的 group 文件改
    control = "HC",
    case_levels = NULL        # NULL 表示所有非 HC 都是 case
  )
  
  ## 评估CF信息冗余和特征信号
  fam_result <- eval_one_profile_file(
    file = ref_file, rep_name = 'CF', group = group, 
    label = label, sample_col = "sample" ) %>%
    mutate(type = "CF")
  
  ## random 结果
  random_result <- map2_dfr(
    random_files, random_names,~ 
      eval_one_profile_file(
        file = .x, rep_name = .y, group = group, label = label, 
        sample_col = "sample" ), .progress = TRUE
  ) %>%
    mutate(type = "Random")
  
  ## 合并
  result <- bind_rows(fam_result, random_result)
  
  ## 计算经验 p 值
  cf <- filter(result, type == "CF")
  rnd <- filter(result, type == "Random")
  
  summary_result <- tibble(
    cf_mean_signal = cf$mean_signal,
    random_mean_signal = mean(rnd$mean_signal, na.rm = TRUE),
    p_signal_high = (1 + sum(rnd$mean_signal >= cf$mean_signal, na.rm = TRUE)) / # 右尾检验
      (1 + nrow(rnd)),
    
    cf_redundancy = cf$mean_abs_cor,
    random_redundancy = mean(rnd$mean_abs_cor, na.rm = TRUE),
    p_redundancy_high = (1 + sum(rnd$mean_abs_cor >= cf$mean_abs_cor, na.rm = TRUE)) / # 右尾检验
      (1 + nrow(rnd)),
    
    cf_ratio = cf$signal_redundancy_ratio,
    random_ratio = mean(rnd$signal_redundancy_ratio, na.rm = TRUE),
    p_ratio_low = (1 + sum(rnd$signal_redundancy_ratio <= cf$signal_redundancy_ratio, na.rm = TRUE)) /
      (1 + nrow(rnd))
  )
  
  ## 可视化：signal vs redundancy
  p1 <- ggplot(result, aes(x = mean_abs_cor, y = mean_signal)) +
    geom_point(
      data = filter(result, type == "Random"),
      alpha = 0.35, size = 1.8
    ) +
    geom_point(
      data = filter(result, type == "CF"),
      color = "red", size = 3, 
    ) +
    theme_pubr() +
    labs(
      x = "Feature-feature redundancy\nMean absolute Spearman correlation",
      y = "Feature-label signal\nMean 2 × |AUC - 0.5|",
      title = "Disease signal\nfeature redundancy"
    ) +
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
  
  ## 可视化：冗余分布
  p2 <- ggdensity(rnd, x = 'mean_abs_cor', fill = 'grey77') +
    geom_vline(xintercept = cf$mean_abs_cor, color = "red") +
    labs(
      x = "Mean absolute Spearman correlation",
      title = "CF feature redundancy\ncompared with random sets",
      subtitle = paste0(
        "Empirical p = ", signif(summary_result$p_redundancy_high, 3)
        )
    ) + 
    theme_pubr() +
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
  
  ## 可视化：signal / redundancy ratio
  p3 <- ggdensity(rnd, x = 'signal_redundancy_ratio', fill = 'grey77') +
    geom_vline(xintercept = cf$signal_redundancy_ratio, color = "red") +
    labs(
      x = "Signal / redundancy ratio",
      title = "Effective disease signal\nper unit redundancy",
      subtitle = paste0("Empirical p = ", signif(summary_result$p_ratio_low, 3))
    ) + 
    theme_pubr() +
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
  
  p <- cowplot::plot_grid(p1, p2, p3, nrow = 1, align = 'hv')
  
  output <- list(
    plot    = p, 
    data    = result,
    summary = summary_result
  )
  
  return(output)
}

proj_name <- c(
  'BushmanFD_2020','FranzosaEA_2018','HallAB_2017','HeQ_2017','KumbhariA_2024',
  'LloydPriceJ_2019','SchirmerM_2018','SchirmerM_2024','WengY_2019','YanQ_2023c'
)

file_base = '../data/profile/'

results <- map(
  proj_name, ~ {
    res <- run_proj(.x, file_base)
    message(
      paste0(.x, '  end ...')
    )
    res
  }
) |> 
  set_names(
    proj_name
  )

write_rds(results, 'results/2.rf.redundancy.rds')

results <- read_rds('results/2.rf.redundancy.rds')

map2(
  results, names(results), ~ 
    cowplot::plot_grid(
      ggplot() +
        geom_text(
          aes(x = 1, y = 1, label  = .y), size = 5,
        ) +
        theme_void(), 
      .x[[1]], 
      rel_widths = c(1, 3),
      nrow = 1)
) |> 
  set_names(
    names(results)
  ) |> 
  cowplot::plot_grid(
    plotlist = _, 
    align = 'hv', 
    ncol = 1
  )

ggsave('results/2.rf.redundancy.pdf', width = 16, height = 40)

data <- map2_df(
  results, names(results), ~
    add_column(.x[[3]], proj = .y, .before = 1)
  )

write.xlsx(data, 'results/2.rf.redundancy.xlsx')

# 雷达图
plot_ggradar_pair <- function(
    data, cf_col, random_col, p_col, title = NULL,  
    proj_col = "proj", cf_name = "CF", random_name = "Random",
    colors = c(CF = "#E64B35", Random = "#4DBBD5"), 
    grid_min = 0, grid_max = NULL, show_p = TRUE, p_digits = 3, 
    legend.position = "bottom", axis.label.size = 3, grid.label.size = 3
  ) {
  
  data <- data.frame(data, check.names = FALSE)
  data[[proj_col]] <- as.character(data[[proj_col]])
  
  ## 1. 生成雷达轴标签
  ## show_p = TRUE 时，把 p 值直接写到项目名下面
  if (isTRUE(show_p)) {
    data <- data |>
      dplyr::mutate(
        proj_lab = ifelse(
          is.na(.data[[p_col]]),
          .data[[proj_col]],
          paste0(.data[[proj_col]], "\n(p = ", signif(.data[[p_col]], p_digits), ")")
        )
      )
  } else {
    data <- data |>
      dplyr::mutate(proj_lab = .data[[proj_col]])
  }
  
  ## 2. 整理成 ggradar 需要的格式
  radar_data <- data |>
    dplyr::select(
      proj = proj_lab,
      CF = dplyr::all_of(cf_col),
      Random = dplyr::all_of(random_col)
    ) |>
    tidyr::pivot_longer(
      cols = c("CF", "Random"),
      names_to = "group",
      values_to = "value"
    ) |>
    tidyr::pivot_wider(
      names_from = proj,
      values_from = value
    ) |>
    dplyr::mutate(group = factor(group, levels = c("CF", "Random"))) |>
    data.frame(check.names = FALSE)
  
  ## 3. 设置刻度范围
  value_mat <- as.matrix(radar_data[, -1])
  
  if (is.null(grid_max)) {
    grid_max <- max(value_mat, na.rm = TRUE) * 1.05
    grid_max <- signif(grid_max, 2)
  }
  
  if (grid_max <= grid_min) {
    stop("grid_max should be larger than grid_min.")
  }
  
  grid_mid <- (grid_min + grid_max) / 2
  
  ## 4. ggradar 作图
  p <- ggradar::ggradar(
    radar_data,
    grid.min = grid_min,
    grid.mid = grid_mid,
    grid.max = grid_max,
    values.radar = c(
      signif(grid_min, 2),
      signif(grid_mid, 2),
      signif(grid_max, 2)
    ),
    group.colours = colors,
    group.line.width = 0.9,
    group.point.size = 2.5,
    background.circle.colour = "white",
    gridline.min.colour = "grey88",
    gridline.mid.colour = "grey82",
    gridline.max.colour = "grey70",
    axis.label.size = axis.label.size,
    grid.label.size = grid.label.size,
    legend.position = legend.position
  ) +
    labs(title = title) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
      legend.text = element_text(size = 10)
    )
  
    # coord_equal(
    #   xlim = c(-.5, .5),
    #   ylim = c(-.5, .5),
    #   clip = "off"
    # )
  
  return(p)
}

p_redundancy <- plot_ggradar_pair(
  data = data,
  cf_col = "cf_redundancy",
  random_col = "random_redundancy",
  p_col = "p_redundancy_high",
  title = "Feature-feature redundancy"
)

p_signal <- plot_ggradar_pair(
  data = data,
  cf_col = "cf_mean_signal",
  random_col = "random_mean_signal",
  p_col = "p_signal_high",
  title = "Feature-phenotype signal"
)

p_redundancy + p_signal + 
  patchwork::plot_layout(nrow = 1) &
  theme(legend.position = "top")
ggsave('results/2.Figure6.rf.radar.pdf', width = 10, height = 6)

#### POG importance ####
proj_name <- c(
  'BushmanFD_2020','FranzosaEA_2018','HallAB_2017','HeQ_2017','KumbhariA_2024',
  'LloydPriceJ_2019','SchirmerM_2018','SchirmerM_2024','WengY_2019','YanQ_2023c'
)

file_base <- '../data/'

tpm <- data.table::fread('merge_tpm.fam.tpm') |> 
  column_to_rownames('name')

tpm <- tpm[grepl('CFF', rownames(tpm)),]

importance <- map(
  proj_name, ~ {
    
    group <- read.delim(
      paste0(
        file_base, .x, 
        '.sample_group.bz2')
    ) %>% 
      dplyr::select(sample, group) %>% 
      mutate(group = ifelse(group == 'HC', 'HC', 'IBD'))
    
    profile <- tpm[, group$sample]
    
    rf_importance(
      profile, group, seed = 2026, ntree = 1000
    )
    
  }, .progress = T
) %>% 
  set_names(proj_name)

rank_df <- map2(
  importance, names(importance), ~ 
    select(.x, name, !!.y := MeanDecreaseAccuracy)
) |> 
  reduce(~ full_join(.x, .y, by = 'name') ) |> 
  mutate(
    across(
      !name, ~ 
        replace_na(.x, 0)
    )
  )

write.xlsx(
  rank_df, 'results/6.feature_importance.xlsx'
)

feature_names <- column_to_rownames(rank_df, 'name') |> 
  rowMeans() |> 
  sort(decreasing = T) |> 
  names()
write_rds(feature_names, 'results/6.feature_importance.names.rds')

# 汇总
files <- list.files('results/6.feature_importance/', pattern = 'top')

results <- map_dfr(
  files, ~ 
    read.delim(
      paste0('results/6.feature_importance/', .x), 
      header = F, col.names = c('top_n', 'proj', 'seed', 'auc')
    )
)

auc_df <- group_by(results, top_n, proj) |> 
  summarise(avg = mean(auc), .groups = 'drop')

auc_df2 <- group_by(auc_df, top_n) |> 
  summarise(mean = mean(avg), .groups = 'drop') |> 
  mutate(
    top_n = sub('top_', '', top_n),
    top_n = as.numeric(top_n)
  )

list(
  auc_by_proj = auc_df, 
  auc_by_top_n = auc_df2
) |> 
  write.xlsx('results/6.feature_topN_AUC.xlsx')

# 0.8255701
ggplot(auc_df2, aes(x = top_n, y = mean)) +
  geom_vline(xintercept = 15, linetype = 'longdash', linewidth = .4) +
  geom_hline(
    yintercept = auc_df2$mean[auc_df2$top_n == 15], 
    linetype = 'longdash', linewidth = .4
  ) +
  geom_line(linewidth = .5, color = 'black') +
  geom_point(size = 3, shape = 21, fill = 'white', color = 'black') +
  scale_x_continuous(
    breaks = seq(0, 79, 10),
    labels = seq(0, 79, 10),
    expand = c(.02, .02)
  ) +
  labs(
    x = 'Top N CF−like POG', 
    y = 'Average AUC across ten cohorts', 
    title = 'Top-ranked CF-like POG feature compression'
  ) +
  theme_pubr() +
  theme(aspect.ratio = 1/3,
        axis.ticks.length = unit(2, 'mm'),
        plot.title = element_text(color = 'black', hjust = .5, face = 'bold'))

ggsave('results/6.Figure6.feature_topN_AUC.pdf', width = 10, height = 4)

#### random test 15 ####
files <- list.files('random_rf_top15/', pattern = 'fam.rf.tsv')

random <- map_dfr(
  files[-1], ~ 
    read.delim(paste0('random_rf_top15/', .x)) |>
    add_column(rep = str_split_i(.x, '\\.', 1), .before = 1) )

ref <- read.xlsx('results/6.feature_topN_AUC.xlsx', sheet = 'auc_by_top_n') |> 
  filter(top_n == 15) |> 
  pull(mean)

.random <- random |>
  dplyr::group_by(rep) |>
  dplyr::summarise(mean = mean(auc), .groups = 'drop') |>
  dplyr::pull(mean)

test <- calcu_empirical_p(ref, .random, 'right')

group_by(random, rep) |> 
  summarise(mean = mean(auc)) |> 
  ggdensity(
    x = 'mean', fill = 'grey77', xlab = 'Average AUC', title = '15 random POG',
    subtitle = paste0('SES=', round(test$SES, 3),'\np=', test$pval)) +
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

ggsave('results/8.Figure6.top15.rf.pdf', width = 4, height = 4)
