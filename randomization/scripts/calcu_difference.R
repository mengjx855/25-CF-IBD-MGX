#### Jinxin Meng, 20220101, 20260523, v0.8.3 ####

# 20220601: 可选择'wilcox rank-sum','one-way anova','student's t test'三种方法做差异分析；
# 20230117: diff_test_profile函数对feature进行差异分析，输入的是标准otu表和group表
# 20231204: 修改diff_test_profile函数中的for循环，使用purrr::map_dfr，速度上稍微快一丢丢。
# 20231204: 修改diff_test函数中的rbind()，使用tibble::add_column()，速度上稍微快一丢丢。
# 20250109: 添加round参数，格式化pval， sig 参数，只输出显著的
# 20250404: 修改部分参数名称，sig <- filter_by_pval, xx_colnames <- xx_rename，使用 map_xx 代替 for循环，加速计算速度
# 20250417: 修改一些BUG
# 20250908: 合并difference_analysis, 改本脚本名称为calcu_difference.R
# 20250915: 更新部分函数，简化函数，使用formula代替分组文件
# 20251013: 更新difference_analysis，计算P与丰度估计使用的数据是否转换要分开；
# 20251105: 修改一些BUG
# 20260523: add new function 'calcu_empirical_p()'

library(tidyverse)

#### add_plab ####
# 根据 p 值生成显著性标签
# format:
#   1: '***', '**', '*', 'ns'
#   2: '***', '**', '*', ''
#   3: '#', '+', '*', ''
#   4: '+', '**', '*', ''
#   5: 'p<0.001', 'p<0.01', 'p<0.05', ''
#   6: round(pvalue, 3)

.add_plab <- function(pvalue, format = 2) {

  if (missing(pvalue)) {
    cat("  format-1: '***', '**', '*', 'ns'\n")
    cat("  format-2: '***', '**', '*',   ''\n")
    cat("  format-3:   '#',  '+', '*',   ''\n")
    cat("  format-4:   '+', '**', '*',   ''\n")
    cat("  format-5: 'p<0.001', 'p<0.01', 'p<0.05', ' '\n")
    cat("  format-6: round(pvalue, digits = 3)\n")
    return(invisible(NULL))
  }
  
  if (format == 6) {
    return(round(pvalue, 3))
  }
  
  labs <- switch(
    as.character(format),
    '1' = c('***', '**', '*', 'ns'),
    '2' = c('***', '**', '*', ''),
    '3' = c('#', '+', '*', ''),
    '4' = c('+', '**', '*', ''),
    '5' = c('p<0.001', 'p<0.01', 'p<0.05', ''),
    c('***', '**', '*', '')
  )
  
  label <- cut(
    pvalue,
    include.lowest = TRUE,
    breaks = c(0, 0.001, 0.01, 0.05, Inf),
    labels = labs
  ) |>
    as.character()
  
  label[is.na(label)] <- ''
  
  return(label)
}


#### calcu_empirical_p ####
# 计算经验 P 值和标准化效应量 SES
# obs: 观测值
# random: 随机 null 分布
# alternative:
#   auto: obs >= random mean 时右尾，否则左尾
#   right: P(random >= obs)
#   left:  P(random <= obs)
#   two.sided: 双尾经验 P 值

calcu_empirical_p <- function(
    obs, random, alternative = c('auto', 'right', 'left', 'two.sided'), 
    simplify = FALSE
  ) {
  
  alternative <- match.arg(alternative)
  
  random <- random[is.finite(random)]
  n <- length(random)
  
  if (n == 0 || !is.finite(obs)) {
    out <- list(
      observed = obs,
      random_mean = NA_real_,
      SES = NA_real_,
      pval = NA_real_,
      side = NA_character_
    )
    
    if (isTRUE(simplify)) {
      out <- data.frame(out, check.names = FALSE)
    }
    
    return(out)
  }
  
  random_mean <- mean(random, na.rm = TRUE)
  
  if (alternative == 'auto') {
    alternative <- ifelse(obs >= random_mean, 'right', 'left')
  }
  
  if (alternative == 'right') {
    p <- (1 + sum(random >= obs, na.rm = TRUE)) / (1 + n)
    side <- 'right.side'
  }
  
  if (alternative == 'left') {
    p <- (1 + sum(random <= obs, na.rm = TRUE)) / (1 + n)
    side <- 'left.side'
  }
  
  if (alternative == 'two.sided') {
    p1 <- (1 + sum(random >= obs, na.rm = TRUE)) / (1 + n)
    p2 <- (1 + sum(random <= obs, na.rm = TRUE)) / (1 + n)
    p <- min(1, 2 * min(p1, p2))
    side <- 'two.sided'
  }
  
  sd_random <- sd(random, na.rm = TRUE)
  ses <- ifelse(sd_random > 0, (obs - random_mean) / sd_random, NA_real_)
  
  out <- list(
    observed = obs,
    random_mean = random_mean,
    SES = ses,
    pval = p,
    side = side
  )
  
  if (isTRUE(simplify)) {
    out <- data.frame(out, check.names = FALSE)
  }
  
  return(out)
}

#### calcu_diff ####
# 根据 formula 进行两两差异检验
# data: 数据框
# formula: value ~ group
# method: wilcox / anova / t
# var_equal: t 检验是否假设方差齐性
# add_plab: 是否添加显著性标签

calcu_diff <- function(data, formula, method = c('wilcox', 'anova', 't'), 
                       var_equal = FALSE, add_plab = FALSE, 
                       plab_fmt = 2, ...) {
  
  method <- match.arg(method)
  
  terms <- terms(formula)
  response <- all.vars(terms)[1]
  group <- all.vars(terms)[2]
  
  if (!response %in% names(data)) 
    stop('variable ', response, ' not in data')
  
  if (!group %in% names(data))
    stop('variable ', group, ' not in data')
  
  data <- data |>
    dplyr::filter(!is.na(.data[[response]]), !is.na(.data[[group]]))
  
  if (is.factor(data[[group]])) {
    group_level <- levels(droplevels(data[[group]]))
  } else {
    group_level <- unique(as.character(data[[group]]))
  }
  
  if (length(group_level) < 2) {
    stop('at least two groups are required.')
  }
  
  comparison <- combn(group_level, m = 2, simplify = FALSE)
  
  difference <- purrr::map_dfr(
    comparison, \(x) {
      
      .data <- dplyr::filter(data, .data[[group]] %in% x)
    
      test <- tryCatch({
        
        if (method == 'wilcox') {
          stats::wilcox.test(formula, .data, ...) |> suppressWarnings()
          
        } else if (method == 'anova') {
          stats::oneway.test(formula, .data, ...)
          
        } else if (method == 't') {
          stats::t.test(formula, .data, var.equal = var_equal, ...) |> suppressWarnings()
        }
        
      }, error = function(e) NULL)
      
      method_name <- dplyr::case_when(
        method == 'wilcox' ~ 'Wilcoxon rank-sum test',
        method == 'anova'  ~ 'One-way ANOVA test',
        method == 't' & isTRUE(var_equal) ~ "Student's t-test",
        method == 't' & isFALSE(var_equal) ~ 'Welch t-test'
      )
      
      data.frame(
        comparison = paste0(x, collapse = '_vs_'),
        pval = ifelse(is.null(test), NA_real_, test$p.value),
        method = method_name
      )
    }
  )
  
  if (nrow(difference) >= 3) {
    difference <- difference |>
      dplyr::mutate(padj = p.adjust(pval, method = 'BH'), .after = 'pval')
  }
  
  if (isTRUE(add_plab)) {
    if ('padj' %in% colnames(difference)) {
      difference <- difference |>
        dplyr::mutate(plab = .add_plab(padj, plab_fmt), .after = 'padj')
    } else {
      difference <- difference |>
        dplyr::mutate(plab = .add_plab(pval, plab_fmt), .after = 'pval')
    }
  }
  
  return(difference)
}

#### calcu_diff_profile ####
# 对 profile 中每个 feature 进行差异检验
# profile: 行为 feature，列为 sample
# group: 样本分组信息，至少包含 sample 和 group_by
# comparison: 指定比较组，例如 c('IBD', 'HC')；NULL 时自动两两比较

calcu_diff_profile <- function(profile, group, group_by = 'group', comparison = NULL,
                               method = c('wilcox', 'anova', 't'), 
                               add_plab = FALSE, plab_fmt = 2, 
                               var_equal = FALSE, progress = TRUE, ...) {

  method <- match.arg(method)
  
  if (!all(c('sample', group_by) %in% colnames(group))) {
    stop('group field should contain: sample | ', group_by)
  }
  
  profile <- data.frame(profile, check.names = FALSE)
  group <- data.frame(group, check.names = FALSE)

  sample_use <- intersect(group$sample, colnames(profile))
  group <- group |>
    dplyr::filter(sample %in% sample_use)
  
  profile <- profile[, group$sample, drop = FALSE]
  
  if (is.null(comparison)) {
    comparison <- combn(unique(group[[group_by]]), m = 2, simplify = FALSE)
  } else if (is.vector(comparison) && !is.list(comparison)) {
    comparison <- list(comparison)
  }
  
  data <- data.frame(t(profile), check.names = FALSE) |>
    tibble::rownames_to_column('sample') |>
    dplyr::left_join(group, by = 'sample')
  
  feature_names <- rownames(profile)
  
  difference <- purrr::map_dfr(
    comparison, \(x) {
      
      purrr::map_dfr(
        feature_names, \(f) {
          
          data |>
            dplyr::select(
              group = dplyr::all_of(group_by),
              value = dplyr::all_of(f)
            ) |>
            dplyr::filter(group %in% x) |>
            calcu_diff(
              value ~ group, 
              method = method, 
              var_equal = var_equal,
              ...
            ) |>
            tibble::add_column(name = f, .before = 1)
        },
        .progress = progress
      )
    }
  )
  
  
  difference <- difference |>
    dplyr::mutate(padj = p.adjust(pval, method = 'BH'), .after = 'pval')
  
  if (isTRUE(add_plab)) {
    difference <- difference |>
      dplyr::mutate(plab = .add_plab(padj, plab_fmt), .after = 'padj')
  }
  
  return(difference)
}


#### difference_analysis ####
# 对 profile 做两组差异分析，同时输出均值、FC、log2FC、prevalence 和 p/q 值
# profile: 行为 feature，列为 sample
# group: 样本分组信息，默认包含 sample 和 group
# comparison: 长度为 2 的比较组，例如 c('IBD', 'HC')
# trans: 用于差异检验的数据转换方式，NULL / RA / LOG10 / LOG2 / SQRT
# only_test: TRUE 表示仅检验数据转换，丰度估计仍使用原始 profile
# method: wilcox / t
# min_abundance: prevalence 判断阈值
# add_enriched: 是否根据 log2FC 和 padj 添加 enriched 分组

difference_analysis <- function(profile, group, group_rename = NULL, comparison = NULL, 
                                trans = NULL, only_test = TRUE, progress = TRUE, 
                                method = c('wilcox', 't'), exact = NULL, 
                                var_equal = FALSE, min_abundance = 0, 
                                add_enriched = FALSE, add_plab = FALSE, 
                                plab_fmt = 2, log2FC = 0, padj = 0.05, ...) {

  method <- match.arg(method)
  
  if (!all(c('sample', 'group') %in% colnames(group)) && is.null(group_rename)) {
    stop('group field should contain: sample | group')
  }
  
  if (!is.null(group_rename)) {
    group <- data.frame(group, check.names = FALSE) |>
      dplyr::rename(dplyr::all_of(group_rename))
  }
  
  if (is.null(comparison) || length(comparison) != 2) {
    stop('comparison should be c(group1, group2).')
  }
  
  profile <- data.frame(profile, check.names = FALSE)
  group <- data.frame(group, check.names = FALSE)
  
  ## 样本过滤和排序
  group <- group |>
    dplyr::filter(group %in% comparison, sample %in% colnames(profile))
  
  profile <- profile[, group$sample, drop = FALSE]
  profile <- profile[rowSums(profile, na.rm = TRUE) != 0, , drop = FALSE]
  
  group$group <- factor(group$group, levels = comparison)
  
  sample_n <- table(as.character(group$group))
  
  message(
    '[', format(Sys.time()), '] Sample: ',
    paste0(paste0(names(sample_n), ' (n=', sample_n, ')'), collapse = ', ')
  )
  
  ## 内部转换函数
  trans_profile <- function(x, trans = NULL) {
    
    x <- as.matrix(x)
    storage.mode(x) <- 'numeric'
    x[is.na(x)] <- 0
    
    if (is.null(trans)) {
      out <- x
      
    } else if (trans == 'RA') {
      out <- apply(x, 2, \(v) {
        if (sum(v, na.rm = TRUE) == 0) v else v / sum(v, na.rm = TRUE)
      } )
      
    } else if (trans %in% c('LOG10', 'LOG2')) {
      min_pos <- min(x[x > 0], na.rm = TRUE)
      if (!is.finite(min_pos)) min_pos <- 1e-6
      x[x <= 0] <- min_pos * 0.001
      
      out <- if (trans == 'LOG10') log10(x) else log2(x)
      
    } else if (trans == 'SQRT') {
      out <- sqrt(x)
      
    } else {
      stop('unknown trans: ', trans)
    }
    
    data.frame(out, check.names = FALSE)
  }
  
  ## 1. 差异检验数据
  test_profile <- trans_profile(profile, trans = trans)
  
  test_data <- data.frame(t(test_profile), check.names = FALSE) |>
    dplyr::mutate(group = group$group[match(rownames(.), group$sample)])
  
  message('[', format(Sys.time()), '] Assessing P-value.')
  
  feature_names <- rownames(profile)
  
  if (method == 'wilcox') {
    
    difference <- purrr::map_dfr(
      feature_names, \(f) {
        .data <- dplyr::select(test_data, value = dplyr::all_of(f), group)
        test <- tryCatch(
          stats::wilcox.test(value ~ group, .data, exact = exact, ...) |> 
            suppressWarnings(),
          error = function(e) NULL
        )
        data.frame(
          name = f,
          pval = ifelse(is.null(test), NA_real_, test$p.value),
          method = 'Wilcoxon rank-sum test'
        )
      },
      .progress = progress
    )
  }
  
  if (method == 't') {
    
    method_name <- ifelse(isTRUE(var_equal), "Student's t-test", 'Welch t-test')
    
    difference <- purrr::map_dfr(
      feature_names, \(f) {
        .data <- dplyr::select(test_data, value = dplyr::all_of(f), group)
        test <- tryCatch(
          stats::t.test(value ~ group, .data, var.equal = var_equal, ...) |>
            suppressWarnings(),
          error = function(e) NULL
        )
        data.frame(
          name = f,
          pval = ifelse(is.null(test), NA_real_, test$p.value),
          method = method_name
        )
      },
      .progress = progress
    )
  }
  
  difference <- difference |>
    dplyr::mutate(padj = p.adjust(pval, method = 'BH'), .after = 'pval')
  
  ## 2. 丰度估计数据
  message('[', format(Sys.time()), '] Mean abundance.')
  
  if (isTRUE(only_test)) {
    ab_profile <- data.frame(profile, check.names = FALSE)
  } else {
    ab_profile <- test_profile
  }
  
  ab_data <- data.frame(t(ab_profile), check.names = FALSE) |>
    dplyr::mutate(group = group$group[match(rownames(.), group$sample)])
  
  abundance <- aggregate(. ~ group, ab_data, mean) |>
    tibble::column_to_rownames('group') |>
    t() |>
    data.frame(check.names = FALSE) |>
    dplyr::select(
      avg_ab1 = dplyr::all_of(comparison[1]),
      avg_ab2 = dplyr::all_of(comparison[2])
    ) |>
    dplyr::mutate(
      comparison = paste0(comparison, collapse = '_vs_'),
      FC = round((avg_ab1 + 1e-6) / (avg_ab2 + 1e-6), 6),
      log2FC = round(log2((avg_ab1 + 1e-6) / (avg_ab2 + 1e-6)), 6),
      avg_ab1 = round(avg_ab1, 6),
      avg_ab2 = round(avg_ab2, 6)
    ) |>
    tibble::rownames_to_column('name')
  
  ## 3. 流行率
  message('[', format(Sys.time()), '] Prevalence.')
  
  prev_data <- data.frame(t(profile), check.names = FALSE) |>
    dplyr::mutate(group = group$group[match(rownames(.), group$sample)])
  
  prevalence <- prev_data |>
    dplyr::group_by(group) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(feature_names),
        \(x) mean(x > min_abundance, na.rm = TRUE)
      ),
      .groups = 'drop'
    ) |>
    tibble::column_to_rownames('group') |>
    t() |>
    data.frame(check.names = FALSE) |>
    dplyr::select(
      prev1 = dplyr::all_of(comparison[1]),
      prev2 = dplyr::all_of(comparison[2])
    ) |>
    dplyr::mutate(
      prev1 = round(prev1, 6),
      prev2 = round(prev2, 6)
    ) |>
    tibble::rownames_to_column('name')
  
  ## 4. 合并结果
  message('[', format(Sys.time()), '] Output result.')
  
  result <- abundance |>
    dplyr::left_join(prevalence, by = 'name') |>
    dplyr::left_join(difference, by = 'name') |>
    dplyr::relocate(prev1, prev2, .after = avg_ab2)
  
  if (isTRUE(add_enriched)) {
    message(
      '[', format(Sys.time()), '] add_enriched: abs(log2FC) > ',
      log2FC, ', padj < ', padj
    )
    
    result <- result |>
      dplyr::mutate(
        enriched = dplyr::case_when(
          .data$log2FC > log2FC & .data$padj < padj ~ comparison[1],
          .data$log2FC < -log2FC & .data$padj < padj ~ comparison[2],
          TRUE ~ 'none'
        )
      )
  }
  
  if (isTRUE(add_plab)) {
    result <- result |>
      dplyr::mutate(plab = .add_plab(padj, plab_fmt), .after = 'padj')
  }
  
  message('[', format(Sys.time()), '] end ...')
  
  return(result)
}
