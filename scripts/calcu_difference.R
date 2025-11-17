#### Jinxin Meng, 20220101, 20251105, v0.8.2 ####

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

library(dplyr)
library(tibble)
library(tidyr)
source('/data/mengjx/R_proj/R_func/profile_process.R')

#### add_plab ####
add_plab <- function(pvalue, format = 2) {
  if (missing(pvalue)) {
    cat("  format-1: '***', '**', '*', 'ns'\n")
    cat("  format-2: '***', '**', '*',   ''\n")
    cat("  format-3:   '#',  '+', '*',   ''\n")
    cat("  format-4:   '+', '**', '*',   ''\n")
    cat("  format-5: 'p<0.001', 'p<0.01', 'p<0.05', ' '\n")
    cat("  format-6: round(pvalue, digits = 3)\n")
    return(invisible(NULL)) 
  } 
  
  if (format == 1) label <- cut(pvalue, include.lowest = T, breaks = c(0, 0.001, 0.01, 0.05, 1), 
                                labels = c('***', '**', '*', 'ns')) %>% as.character()
  if (format == 2) label <- cut(pvalue, include.lowest = T, breaks = c(0, 0.001, 0.01, 0.05, 1), 
                                labels = c('***', '**', '*', '')) %>% as.character()
  if (format == 3) label <- cut(pvalue, include.lowest = T, breaks = c(0, 0.001, 0.01, 0.05, 1), 
                                labels = c('#', '+', '*', '')) %>% as.character()
  if (format == 4) label <- cut(pvalue, include.lowest = T,breaks = c(0, 0.001, 0.01, 0.05, 1), 
                                labels = c('+', '**', '*', '')) %>% as.character()
  if (format == 5) label <- cut(pvalue, include.lowest = T, breaks = c(0, 0.001, 0.01, 0.05, 1), 
                                labels = c('p<0.001', 'p<0.01', 'p<0.05', '')) %>% as.character()
  if (format == 6) label <- round(pvalue, 3)
  
  return(label)
}

#### calcu_diff ####
calcu_diff <- function(data, formula, method = 'wilcox', var_equal = FALSE, 
                       add_plab = FALSE, plab_fmt = 2, ...) {
  
  terms <- terms(formula)
  response <- all.vars(terms)[1]
  group <- all.vars(terms)[2]
  
  if (!response %in% names(data)) 
    stop(paste('variable ', response, 'not in data'))
  
  if (!group %in% names(data))
    stop(paste('variable', group, 'not in data'))
  
  if ( all(is.factor(data[[group]])) ) {
    group_level <- levels(data[[group]])
  } else {
    group_level <- unique(data[[group]])
  }
  
  if (method == 'wilcox')
    difference <- purrr::map_dfr(combn(group_level, m = 2, simplify = F), ~ {
      test <- stats::wilcox.test(formula, dplyr::filter(data, .data[[group]] %in% .x)) %>% suppressWarnings()
      data.frame(comparison = paste0(.x, collapse = '_vs_'), 
                 pval = test$p.value, 
                 method = 'Wilcoxon-rank sum test') } )
  
  if (method == 'anova')
    difference <- purrr::map_dfr(combn(group_level, m = 2, simplify = F), ~ {
      test <- stats::oneway.test(formula, dplyr::filter(data, .data[[group]] %in% .x))
      data.frame(comparison = paste0(.x, collapse = '_vs_'),
                 pval = test$p.value, 
                 method = 'One-way ANOVA test') } )
  
  if (method == 't') {
    if (isTRUE(var_equal)) method_name = "student's t-test"
    if (isFALSE(var_equal)) method_name = 'Welch t-test'
    difference <- purrr::map_dfr(combn(group_level, m = 2, simplify = F), ~ {
      test <- stats::t.test(formula, dplyr::filter(data, .data[[group]] %in% .x), var.equal = var_equal)
      data.frame(comparison = paste0(.x, collapse = '_vs_'), 
                 pval = test$p.value, 
                 method = method_name) } ) 
  }
  
  if (nrow(difference) >= 3) 
    difference <- mutate(difference, padj = p.adjust(pval, method = 'BH'), .after = 'pval')
  
  if (isTRUE(add_plab)) {
    if ('padj' %in% colnames(difference)) {
      difference <- mutate(difference, plab = add_plab(padj, plab_fmt), .after = 'padj')
    } else {
      difference <- mutate(difference, plab = add_plab(pval, plab_fmt), .after = 'pval')
    }
  }
  
  return(difference)
}

#### calcu_diff_profile ####
calcu_diff_profile <- function(profile, group, group_by = 'group', comparison = NULL,
                               method = 'wilcox', add_plab = FALSE, plab_fmt = 2, 
                               var_equal = FALSE, ...) {
  
  if (!all(c('sample', group_by) %in% colnames(group))) 
    stop('group field (sample|', group_by, '|...)')
  
  if (is.null(comparison)) {
    comparison <- combn(unique(group[[group_by]]), m = 2, simplify = F)
  } else if (is.vector(comparison)) {
    comparison <- list(comparison)
  }
  
  .names <- rownames(profile)
  .length <- length(.names)
  
  data <- data.frame(t(profile), check.names = F) %>% 
    rownames_to_column('sample') %>% 
    left_join(group, by = 'sample')
  
  difference <- purrr::map_dfr(
    comparison, \(x) purrr::map_dfr(.names, ~ {
      dplyr::select(data, group, value = all_of(.x)) %>% 
      dplyr::filter(group %in% x) %>% 
      calcu_diff(value ~ group, method = method, var_equal = var_equal) %>% 
      tibble::add_column(name = .x, .before = 1)
  }, .progress = T ) )
  
  difference <- mutate(difference, padj = p.adjust(pval, method = 'BH'), .after = 'pval')
  
  if (isTRUE(add_plab))
    difference <- mutate(difference, plab = add_plab(padj, plab_fmt), .after = 'padj')
  
  return(difference)
} %>% suppressWarnings()

#### difference_analysis ####
difference_analysis <- function(profile, group, group_rename = NULL, comparison = NULL, 
                                trans = NULL, only_test = TRUE, progress = TRUE, 
                                method = 'wilcox', exact = NULL, var_equal = FALSE, 
                                min_abundance = 0, add_enriched = FALSE, add_plab = FALSE, 
                                plab_fmt = 2, log2FC = 0, padj = 0.05, ... ) {
  
  if (!all(c('sample', 'group') %in% colnames(group)) & is.null(group_rename))
    stop('group field (sample|group|...)')
  
  if (!is.null(group_rename)) 
    group <- data.frame(group, check.names = F) %>% 
      dplyr::rename(all_of(group_rename))
  
  if (is.null(comparison)) 
    stop('missing comparison')
  
  if (!all(as.character(group$group) %in% comparison)) {
    group <- dplyr::filter(group, group %in% comparison)
    profile <- dplyr::select(profile, any_of(group$sample))
  }
  
  profile <- data.frame(profile, check.names = F)
  profile <- profile[rowSums(profile) != 0,]
  
  sample_n <- table(as.character(group$group))
  
  message(paste0('[', format(Sys.time()), '] Sample: ', 
                 paste0(paste0(names(sample_n), ' (n=', sample_n, ')'), collapse = ', ')) )
  
  # 差异检测数据
  if (is.null(trans)) { 
    data <- data.frame(t(profile), check.names = F)
  } else if (trans == 'RA') { 
    data <- data.frame(t(apply(profile, 2, \(x) x/sum(x))), check.names = F)
  } else if (trans == 'LOG10') { 
    data <- apply(profile, 2, \(x) ifelse(x == 0, x + .001 * min(profile), x)) 
    data <- data.frame(t(log10(data)), check.names = F)
  } else if (trans == 'LOG2') {
    data <- apply(profile, 2, \(x) ifelse(x == 0, x + .001 * min(profile), x)) 
    data <- data.frame(t(log2(data)), check.names = F)
  } else if (trans == 'SQER') { 
    data <- data.frame(t(sqrt(profile)), check.names = F)
  }
  
  data <- dplyr::mutate(data, group = group$group[match(rownames(data), group$sample)])
  
  # 估计P值
  message(paste0('[', format(Sys.time()), '] Assessing P-value.'))
  
  .names <- rownames(profile)
  
  if (method == 'wilcox')
    difference <- map_dfr(.names, ~ {
      .data <- dplyr::select(data, value = all_of(.x), group)
      test <- stats::wilcox.test(value ~ group, .data, exact = exact) %>% suppressWarnings()
      data.frame(name = .x, pval = test$p.value, method = 'Wilcoxon rank-sum test') },
      .progress = progress)
  
  if (method == 't') {
    if (isTRUE(var_equal)) method_name = "student's t-test"
    if (isFALSE(var_equal)) method_name = 'Welch t-test'
    difference <- map_dfr(.names, ~ {
      .data <- dplyr::select(data, value = all_of(.x), group)
      test <- stats::t.test(value ~ group, .data, var.equal = var_equal) %>% suppressWarnings()
      data.frame(name = .x, pval = test$p.value, method = method_name) },
      .progress = progress)
  }
  
  difference <- dplyr::mutate(
    difference, 
    padj = p.adjust(pval, method = 'BH'), .after = 'pval',
    pval = round(pval, 6),
    padj = round(padj, 6))
  
  # 计算丰度
  message(paste0('[', format(Sys.time()), '] Mean abundance.'))
  
  if (isTRUE(only_test)) 
    data <- data.frame(t(profile), check.names = F) 
  
  data <- data.frame(t(profile), check.names = F) %>% 
    dplyr::mutate(group = group$group[match(rownames(.), group$sample)])
  
  abundance <- aggregate(. ~ group, data, mean) %>%
    tibble::column_to_rownames('group') %>% 
    t() %>% 
    data.frame(check.names = F) %>% 
    dplyr::select(avg_ab1 = all_of(comparison[1]), 
                  avg_ab2 = all_of(comparison[2])) %>% 
    dplyr::mutate(comparison = paste0(comparison, collapse = '_vs_'),
                  FC = round((avg_ab1 + 1e-6) / (avg_ab2 + 1e-6), 6),
                  log2FC = round(log2((avg_ab1 + 1e-6) / (avg_ab2 + 1e-6)), 6),
                  avg_ab1 = round(avg_ab1, 6),
                  avg_ab2 = round(avg_ab2, 6)) %>% 
    tibble::rownames_to_column('name')
  
  # 计算流行率
  message(paste0('[', format(Sys.time()), '] Prevalence.'))
  prevalence <- dplyr::group_by(data, group) %>% 
    dplyr::group_modify(~ purrr::map_df(.x, \(x) sum(x > min_abundance)/length(x))) %>% 
    tibble::column_to_rownames('group') %>% 
    t() %>% 
    data.frame(check.names = F) %>% 
    dplyr::select(prev1 = all_of(comparison[1]), 
                  prev2 = all_of(comparison[2])) %>% 
    mutate(prev1 = round(prev1, 6),
           prev2 = round(prev2, 6)) %>% 
    tibble::rownames_to_column('name')
  
  # 整理输出
  message(paste0('[', format(Sys.time()), '] Output result.'))
  result <- dplyr::left_join(abundance, prevalence, by = 'name') %>% 
    dplyr::left_join(difference, by = 'name') %>% 
    dplyr::relocate(prev1, prev2, .after = avg_ab2)
  
  if(isTRUE(add_enriched)) {
    message(paste0('[', format(Sys.time()), '] add_enriched: abs(log2FC) > ', log2FC, ', padjust < ', padj))
    result$enriched <- ifelse(result$log2FC > log2FC & result$padj < padj, comparison[1], 
                              ifelse(result$log2FC < -log2FC & result$padj < padj, comparison[2], 'none'))
  }
  
  if (isTRUE(add_plab)) result <- mutate(result, plab = add_plab(padj, plab_fmt), .after = padj)
  
  message(paste0('[', format(Sys.time()), '] end ...'))
  return(result)
}
