#### Jinxin Meng, 202209014, 20260516, v0.1.1 ####

# 20250405: update some parameter.
# 20260516: update some details.

library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)

#### calcu_metafor ####
# 多 feature 的 meta-analysis，适合比较多个项目中 case vs control 的效应方向。
#
# 参考：
#   https://www.programmingr.com/tutorial/arcsine-transformation/
#   https://bookdown.org/robcrystalornelas/meta-analysis_of_ecological_data/#recommended-citation
# 
# 输入数据格式:
#   data = sample + features + group_col + proj_col
# 
# 示例:
#   sample   feat1   feat2   ...   group   proj
#   S1       0.12    1.35          HC      Study1
#   S2       0.30    0.85          IBD     Study1
#   S3       0.08    1.10          HC      Study2
#   S4       0.42    0.60          IBD     Study2
# 
# 参数说明:
#   sample_col: 样本列名，默认 'sample'
#   group_col : 分组列名，默认 'group'
#   proj_col  : 项目/队列列名，默认 'proj'
#   comparison: c(case, control)，第一个为病例/处理组，第二个为对照组
#   measure   : 效应值类型，默认 'SMD'，即标准化均数差，metafor 中对应 Hedges' g
#   method    : 随机效应模型中 tau^2 的估计方法，默认 'REML'
# 
# 结果解释:
#   estimate > 0 表示 comparison[1] 组高于 comparison[2] 组
#   estimate < 0 表示 comparison[2] 组高于 comparison[1] 组
# 
# 输出主要字段:
#   yi         : 每个项目中的效应值
#   vi         : 每个项目中效应值的方差
#   estimate   : 随机效应模型汇总效应值
#   ci_lb/ci_ub: 汇总效应值的 95% 置信区间
#   pval       : 汇总效应值的显著性
#   tau2       : 项目间方差，认为不同研究之间有一些其他因素导致的差异，总异质性
#   I2         : 异质性比例，去判断案例间差异大小占总差异的指标之一，但是I2不可以作为选择哪种模型（固定vs.随机）的依据
#   Q/Q_pval   : 异质性检验统计量及其 p 值，效应值总体的异质性，是评价效应值的差异程度
#                Qt越大，则效应值越离散，暗示我们有些因素对效应值有强烈的影响
#                Qt的优势是可以进行显著性检验的。如果p值不显著，那么我们认为案例间的差异是随机因素造成的，这种情况下就不需要往下继续进行Meta分析
# 
# 菌群相对丰度数据的推荐转换:
#   1) 0–1 相对丰度:
#        profile_meta <- asin(sqrt(profile))
#   2) 0–100 百分比丰度:
#        profile_meta <- asin(sqrt(profile / 100))
#   3) count table:
#        profile_ra <- sweep(profile, 2, colSums(profile), '/')
#        profile_meta <- asin(sqrt(profile_ra))
#   4) 偏态明显或低丰度菌较多时，也可作为敏感性分析使用:
#        profile_meta <- log10(profile + 1e-6)
#   5) 若希望更符合 compositional data 理论，可使用 CLR 转换:
#        profile_pc <- profile + 1e-6
#        profile_meta <- sweep(log(profile_pc), 2, colMeans(log(profile_pc)), '-')
#   6) Shannon、Richness、Simpson 等 alpha diversity 指数:
#        通常直接使用原值；若 richness 偏态明显，可考虑 log 转换

calcu_metafor <- function(data, sample_col = 'sample', group_col = 'group', 
                          proj_col = 'proj', comparison = c('Case', 'Control'),
                          measure = 'SMD', method = 'REML',simplify = FALSE, 
                          quiet = FALSE, progress = TRUE, ...) {
  
  start_time <- Sys.time()
  
  if (isTRUE(quiet)) progress = FALSE
  
  data <- data.frame(data, check.names = FALSE)
  
  ## 1. 检查必要字段
  ## data 至少需要包含样本列、分组列和项目/队列列。
  need_cols <- c(sample_col, group_col, proj_col)
  if (!all(need_cols %in% colnames(data))) {
    stop(
      'data should contain columns: ',
      paste(need_cols, collapse = ' | ')
    )
  }
  
  if (length(comparison) != 2) {
    stop('comparison should be c(case, control).')
  }
  
  ## 2. 统一内部列名
  ## 后续函数内部统一使用 sample/group/proj，方便写循环和汇总。
  data <- data |>
    dplyr::rename(
      sample = dplyr::all_of(sample_col),
      group = dplyr::all_of(group_col),
      proj = dplyr::all_of(proj_col)
    )
  
  ## 3. 提取 feature 列
  ## 除 sample/group/proj 之外，其余列都作为待分析的 feature。
  feature <- setdiff(colnames(data), c('sample', 'group', 'proj'))
  
  if (length(feature) == 0) {
    stop('No feature columns found in data.')
  }
  
  case_group <- comparison[1]
  control_group <- comparison[2]
  
  meta_out <- list()
  
  if (isTRUE(progress)) {
    pb <- utils::txtProgressBar(
      min = 0, max = length(feature), style = 3,
      width = 50, char = '#')
    flag = 1
    on.exit(close(pb), add = TRUE)
  }
  
  ## 4. 逐个 feature 计算项目内效应值，再进行随机效应 meta-analysis
  for (i in feature) {
    
    ## 4.1 计算每个项目中 case/control 的均值、标准差和样本量
    ## SMD/Hedges' g 需要每组的 mean、sd 和 n
    meta_stat <- data |>
      dplyr::filter(group %in% comparison) |>
      dplyr::select(proj, group, index = dplyr::all_of(i)) |>
      dplyr::mutate(index = suppressWarnings(as.numeric(index))) |>
      dplyr::group_by(proj, group) |>
      dplyr::summarise(
        Mean = mean(index, na.rm = TRUE),
        Sd   = stats::sd(index, na.rm = TRUE),
        N    = sum(!is.na(index)),
        .groups = 'drop'
      ) |>
      tidyr::pivot_wider(
        names_from = group,
        values_from = c(Mean, Sd, N),
        names_glue = '{.value}_{group}'
      )
    
    ## 4.2 检查该 feature 是否同时包含 case 和 control 组
    need_stat <- c(
      paste0('Mean_', case_group), paste0('Sd_', case_group), paste0('N_', case_group),
      paste0('Mean_', control_group), paste0('Sd_', control_group), paste0('N_', control_group)
    )
    
    if (!all(need_stat %in% colnames(meta_stat))) {
      if (isFALSE(quiet)) message('Warning: ', i, ' lacks case/control groups.')
      next
    }
    
    ## 4.3 整理为 metafor::escalc() 需要的格式
    ## x_* 表示 case 组；y_* 表示 control 组。
    meta_in <- meta_stat |>
      dplyr::transmute(
        proj = proj,
        x_Mean = .data[[paste0('Mean_', case_group)]],
        x_Sd   = .data[[paste0('Sd_', case_group)]],
        x_N    = .data[[paste0('N_', case_group)]],
        y_Mean = .data[[paste0('Mean_', control_group)]],
        y_Sd   = .data[[paste0('Sd_', control_group)]],
        y_N    = .data[[paste0('N_', control_group)]]
      ) |>
      dplyr::filter(         # 如果某个项目每组只有 1 个样本，sd() 是 NA，SMD 不可靠
        x_N >= 2, y_N >= 2,
        !is.na(x_Sd), !is.na(y_Sd),
        x_Sd >= 0, y_Sd >= 0,
        x_Sd > 0 | y_Sd > 0  # 如果 case 和 control 两组标准差都是 0，SMD 无法稳定计算
      )
    
    if (nrow(meta_in) == 0) {
      if (isFALSE(quiet)) message('Warning: ', i, ' is not suitable.')
      next
    }
    
    ## 4.4 计算每个项目内的效应值 yi 和方差 vi
    smd_meta <- metafor::escalc(
      measure = measure, data = meta_in, append = TRUE,
      m1i = x_Mean, m2i = y_Mean, 
      sd1i = x_Sd, sd2i = y_Sd, 
      n1i = x_N, n2i = y_N
    ) |>
      dplyr::filter(!is.na(yi), !is.na(vi))
    
    if (nrow(smd_meta) == 0) {
      if (isFALSE(quiet)) message('Warning: ', i, ' has no valid effect size.')
      next
    }
    
    ## 4.5 随机效应模型汇总效应值
    ## estimate > 0 表示 case_group 高于 control_group
    smd_rma <- tryCatch(
      metafor::rma(
        yi, vi, method = method, data = smd_meta,
        control = list(stepadj = 0.5, maxiter = 10000),
        ...
      ),
      error = function(e) {
        if (isFALSE(quiet)) message('Warning: ', i, ' rma failed: ', e$message)
        return(NULL)
      }
    )
    
    if (is.null(smd_rma)) next
    
    ## 4.6 合并项目内效应值和 meta-analysis 汇总结果
    if (isTRUE(simplify)) {
      
      meta_out[[i]] <- tibble::tibble(
        feature = i,
        measure     = measure,                   ## 效应值的计算方法
        model       = 'RM',                      ## 累积效应值计算模型
        method_tau2 = method,                    ## 随机效应模型估计案例内方差（Tau^2）的计算方法
        tau2        = as.numeric(smd_rma$tau2),  ## 案例间方差的值
        I2          = round(smd_rma$I2, 2),      ## 案例间差异大小占总差异的比例
        Q           = smd_rma$QE,                ## 异质性检验
        Q_pval      = smd_rma$QEp,               ## 异质性检验P值 越显著异质性越大
        estimate    = as.numeric(smd_rma$beta),
        ci_lb       = smd_rma$ci.lb,
        ci_ub       = smd_rma$ci.ub,
        pval        = smd_rma$pval
      )
      
    } else {
      
      meta_out[[i]] <- smd_rma$data |>
        tibble::add_column(
        measure     = measure,                   ## 效应值的计算方法
        model       = 'RM',                      ## 累积效应值计算模型
        method_tau2 = method,                    ## 随机效应模型估计案例内方差（Tau^2）的计算方法
        tau2        = as.numeric(smd_rma$tau2),  ## 案例间方差的值
        I2          = round(smd_rma$I2, 2),      ## 案例间差异大小占总差异的比例
        Q           = smd_rma$QE,                ## 异质性检验
        Q_pval      = smd_rma$QEp,               ## 异质性检验P值 越显著异质性越大
        estimate    = as.numeric(smd_rma$beta),
        ci_lb       = smd_rma$ci.lb,
        ci_ub       = smd_rma$ci.ub,
        pval        = smd_rma$pval
        ) |> 
        tibble::add_column(feature = i, .before = 1)
    }
    
    if (isTRUE(progress)) {
      utils::setTxtProgressBar(pb, flag)
      flag = flag + 1
    }
    
  }

  ## 5. 合并所有 feature 的结果
  meta_out <- dplyr::bind_rows(meta_out)
  
  if (isFALSE(quiet)) {
    if (isTRUE(progress)) {
      cat('\n')
    }
    message(
      '    estimate > 0 -> ', comparison[1], '; estimate < 0 -> ', comparison[2]
    )
  }
  
  return(meta_out)
}

#### calcu_metafor.1 一种指数 ####
# 针对一个指标或一个 feature 的 meta-analysis
#
# 输入数据需要包含:
#   name, m1, sd1, n1, m2, sd2, n2
#
# 结果解释:
#   estimate > 0 表示第 1 组 m1 高于第 2 组 m2
#   estimate < 0 表示第 2 组 m2 高于第 1 组 m1
#
# 参数说明:
#   measure : 效应值类型，默认 'SMD'，即标准化均数差，metafor 中对应 Hedges' g
#   method  : 随机效应模型中 tau^2 的估计方法，默认 'REML'
#   simplify: TRUE 只返回汇总结果；FALSE 返回每个研究效应值 + 汇总结果

calcu_metafor.1 <- function(data, data_rename = NULL, measure = 'SMD', 
                            method = 'REML', simplify = FALSE, 
                            quiet = FALSE, ...) {
  
  data <- data.frame(data, check.names = FALSE)
  
  ## 1. 字段重命名
  ## data_rename 可用于把外部列名统一改为 name/m1/sd1/n1/m2/sd2/n2
  if (!is.null(data_rename)) {
    data <- data |>
      dplyr::rename(dplyr::all_of(data_rename))
  }
  
  ## 2. 检查必要字段
  need_cols <- c('name', 'm1', 'sd1', 'n1', 'm2', 'sd2', 'n2')
  
  if (!all(need_cols %in% colnames(data))) {
    stop(
      "data field should contain: name | m1 | sd1 | n1 | m2 | sd2 | n2
      
     😀 please refer to following dataset:
         name      m1     sd1   n1    m2    sd2   n2
         <chr>     <dbl>  <dbl> <int> <dbl> <dbl> <int>
      1  HB-MSP    1307.  264.  16    1051. 146.  8
      2  HL-MP     1763.  186.  20    1051. 146.  8
      3  HL2-MP    1578.  140.  22    1051. 146.  8
      4  HN-YNBP   1354.  180.  19    1051. 146.  8
      5  IM-IMBP   1485.  135.  17    1051. 146.  8
      6  ...       ...    ...   ...   ...   ...   ..."
    )
  }
  
  ## 3. 转换为数值，并过滤不适合计算 SMD 的行
  ## n1/n2 至少为 2；sd 不能缺失；两组 sd 不能同时为 0。
  data <- data |>
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(c('m1', 'sd1', 'n1', 'm2', 'sd2', 'n2')),
        \(x) suppressWarnings(as.numeric(x))
      )
    ) |>
    dplyr::filter(
      n1 >= 2, n2 >= 2,
      !is.na(m1), !is.na(m2),
      !is.na(sd1), !is.na(sd2),
      sd1 >= 0, sd2 >= 0,
      sd1 > 0 | sd2 > 0
    )
  
  if (nrow(data) == 0) {
    if (isFALSE(quiet)) warning('No valid rows for meta-analysis.')
    return(data.frame())
  }
  
  ## 4. 计算每个研究内的效应值 yi 和方差 vi
  smd <- metafor::escalc(
    measure = measure, data = data, append = TRUE, 
    m1i = m1, m2i = m2, 
    sd1i = sd1, sd2i = sd2,
    n1i = n1, n2i = n2
  ) |>
    dplyr::filter(!is.na(yi), !is.na(vi)) |>
    data.frame(check.names = FALSE)
  
  if (nrow(smd) == 0) {
    if (isFALSE(quiet)) warning('No valid effect size was calculated.')
    return(data.frame())
  }
  
  ## 5. 随机效应模型汇总效应值
  smd_rma <- tryCatch(
    metafor::rma(
      yi, vi, method = method, data = smd,
      control = list(stepadj = 0.5, maxiter = 10000),
      ...
    ),
    error = function(e) {
      if (isFALSE(quiet)) warning('rma failed: ', e$message)
      return(NULL)
    }
  )
  
  if (is.null(smd_rma)) {
    return(data.frame())
  }
  ## 6. 输出结果
  ## simplify = FALSE: 保留每个研究的 yi/vi，并附加汇总结果
  ## simplify = TRUE : 只输出一行汇总结果
  if (isFALSE(simplify)) {
    
    smd_out <- smd_rma$data |>
      tibble::add_column(
        measure     = measure,                   ## 效应值的计算方法          
        model       = 'RM',                      ## 累积效应值计算模型
        method_tau2 = method,                    ## 随机效应模型估计案例内方差（Tau^2）的计算方法
        tau2        = as.numeric(smd_rma$tau2),  ## 案例间方差的值
        I2          = round(smd_rma$I2, 2),      ## 案例间差异大小占总差异的比例
        Q           = smd_rma$QE,                ## 异质性检验
        Q_pval      = smd_rma$QEp,               ## 异质性检验P值 越显著异质性越大
        estimate    = as.numeric(smd_rma$beta),
        ci_lb       = smd_rma$ci.lb,
        ci_ub       = smd_rma$ci.ub,
        pval        = smd_rma$pval
      )
    
  } else {
    
    smd_out <- data.frame(
      measure     = measure,                    ## 效应值的计算方法
      model       = 'RM',                       ## 累积效应值计算模型
      method_tau2 = method,                     ## 随机效应模型估计案例内方差（Tau^2）的计算方法
      tau2        = as.numeric(smd_rma$tau2),   ## 案例间方差的值
      I2          = round(smd_rma$I2, 2),       ## 案例间差异大小占总差异的比例
      Q           = smd_rma$QE,                 ## 异质性检验
      Q_pval      = smd_rma$QEp,                ## 异质性检验P值 越显著异质性越大
      estimate    = as.numeric(smd_rma$beta),
      ci_lb       = smd_rma$ci.lb,
      ci_ub       = smd_rma$ci.ub,
      pval        = smd_rma$pval
    )
  }
  
  if (isFALSE(quiet)) {
    message(
      '== estimate > 0, ==> group 1 / m1 ==\n',
      '== estimate < 0, ==> group 2 / m2 =='
    )
  }
  
  return(smd_out)
}
