#### Jinxin Meng, 202205029, 20260519, v0.3.4 ####

# 20250917: update function. rename 'rf_feature_rank()' as 'rf_importance()'
# 20251106: update rf_next_vaildate()
# 20260509: update rf_Kfold() and rf_Kfold_rep()
# 20260519: simplify code

library(dplyr)
library(tibble)
library(tidyr)

#### K折随机森林建模 ####
# K折随机森林建模，返回一个数据框，输出一个预测数据框
# 输入一个profile表; group文件是样本的分组信息，分组的列名为sample和group
# k为几折检验; seed设置随机种子, ntree设置随机森林的树数量
rf_Kfold <- function(profile, group, k = 5, seed = 2026, ntree = 1000,
                     sample_col = 'sample', group_col = 'group', 
                     remove_zero_var = TRUE, na_fill = NULL, 
                     progress = TRUE, quiet = FALSE, ...) {
  start_time <- Sys.time()
  
  if (isTRUE(quiet)) progress = FALSE
  
  ## 1. 预处理数据
  X <- data.frame(t(profile), check.names = FALSE)
  group <- data.frame(group, check.names = FALSE)
  
  if (!all(c(sample_col, group_col) %in% colnames(group))) {
    stop('group should contain columns: ', sample_col, ' | ', group_col)
  }
  
  y <- group[[group_col]][match(rownames(X), group[[sample_col]])]
  
  if (any(is.na(y))) {
    stop(
      'Some samples in profile are missing in group table: ',
      paste(rownames(X)[is.na(y)], collapse = ', ')
    )
  }
  
  y <- factor(y)
  
  if (nlevels(y) < 2) {
    stop('Random forest classification requires at least two groups.')
  }
  
  if (k < 2 || k > nrow(X)) {
    stop('k should be >= 2 and <= sample number.')
  }
  
  ## 2. 转成数值矩阵
  X <- data.frame(
    lapply(X, \(x) suppressWarnings(as.numeric(x))),
    check.names = FALSE,
    row.names = rownames(X)
  )
  
  if (any(is.na(X))) {
    if (is.null(na_fill)) {
      stop('NA values found in profile. Please remove/impute them or set na_fill.')
    } else {
      X[is.na(X)] <- na_fill
    }
  }
  
  ## 3. 去除零方差 feature
  if (isTRUE(remove_zero_var)) {
    keep <- vapply(X, stats::var, numeric(1)) > 0
    X <- X[, keep, drop = FALSE]
  }
  
  if (ncol(X) == 0) {
    stop('No valid features remained for random forest.')
  }
  
  ## 4. make.names() 重命名 feature，避免复杂列名导致模型报错
  colnames(X) <- make.names(colnames(X))
  
  ## 5. 普通随机 K 折
  set.seed(seed)
  idx <- sample(seq_len(nrow(X)))
  folds <- split(idx, rep(seq_len(k), length.out = nrow(X)))
  
  pred_list <- vector('list', k)
  
  if (isTRUE(progress)) {
    pb <- utils::txtProgressBar(min = 0, max = k, style = 3, width = 50, char = '#')
    on.exit(close(pb), add = TRUE)
  }
  
  ## 6. 逐折训练和预测
  for (i in seq_len(k)) {
    
    test_idx <- folds[[i]]
    train_idx <- setdiff(seq_len(nrow(X)), test_idx)
    
    train_x <- X[train_idx, , drop = FALSE]
    test_x  <- X[test_idx, , drop = FALSE]
    train_y <- y[train_idx]
    test_y  <- y[test_idx]
    
    set.seed(seed + i)
    
    rf_model <- randomForest::randomForest(
      x = train_x, y = train_y, ntree = ntree, 
      importance = FALSE, proximity = FALSE, ...
    )
    
    pred_class <- stats::predict(rf_model, test_x, type = 'response')
    pred_prob <- data.frame(
      stats::predict(rf_model, test_x, type = 'prob'),
      check.names = FALSE
    )
    
    ## 保证概率列顺序和原始分组水平一致
    pred_list[[i]] <- pred_prob[, levels(y), drop = FALSE] |> 
      tibble::rownames_to_column('sample')
    
    if (isTRUE(progress)) {
      utils::setTxtProgressBar(pb, i)
    }
  }
  
  pred <- do.call(rbind, pred_list)
  
  if (isFALSE(quiet)) {
    if (isTRUE(progress)) cat('\n')
    message(
      'Random forest K-fold finished. Time = ', 
      round(as.numeric(difftime(Sys.time(), start_time, units = 'secs')), 2),
      ' secs.'
    )
  }
  
  return(pred)
}

#### K折rep次重复 ####
# 重复 K 折随机森林建模
# profile 输入格式同 rf_Kfold：行为 feature，列为 sample
# group 文件为样本分组信息，默认列名为 sample 和 group
# rep 为重复次数；每次使用 seed + i - 1 作为随机种子
# 返回每个样本在 rep 次重复后的平均预测概率、最终预测类别和正确性
rf_Kfold_rep <- function(profile, group, k = 5, rep = 5, seed = 2026, ntree = 1000,
                         sample_col = 'sample', group_col = 'group',
                         remove_zero_var = TRUE, na_fill = NULL,
                         progress = TRUE, quiet = FALSE, ...) {
  
  start_time <- Sys.time()
  
  if (isTRUE(quiet)) progress = FALSE
  
  pred_list <- vector('list', rep)
  
  if (isTRUE(progress)) {
    pb <- utils::txtProgressBar(min = 0, max = rep, style = 3, width = 50, char = '#')
    on.exit(close(pb), add = TRUE)
  }
  
  for (i in seq_len(rep)) {
    pred_i <- rf_Kfold(
      profile = profile, 
      group = group,
      k = k,
      seed = seed + i - 1,
      ntree = ntree,
      sample_col = sample_col,
      group_col = group_col,
      remove_zero_var = remove_zero_var,
      na_fill = na_fill,
      progress = FALSE,
      quiet = TRUE,
      ...
    )
    
    pred_i$rep <- paste0('rep_', i)
    pred_list[[i]] <- pred_i
    
    if (isTRUE(progress)) {
      utils::setTxtProgressBar(pb, i)
    }
  }
  
  pred_all <- dplyr::bind_rows(pred_list)
  
  ## 概率列：排除非概率信息列后剩下的就是各类别预测概率
  prob_cols <- setdiff(colnames(pred_all), c('sample', 'rep'))
  
  ## 对每个样本的预测概率取均值
  pred <- pred_all |>
    dplyr::group_by(sample) |>
    dplyr::summarise(
      dplyr::across(dplyr::all_of(prob_cols), \(x) mean(x, na.rm = TRUE)),
      .groups = 'drop'
    )
  
  if (isFALSE(quiet)) {
    if (isTRUE(progress)) cat('\n')
    message(
      'Repeated random forest K-fold finished. Time = ',
      round(as.numeric(difftime(Sys.time(), start_time, units = 'secs')), 2),
      ' secs.'
    )
  }
  
  return(pred)
}

#### 多数据集交叉验证 ####
# 使用一个 dataset 作为训练集，依次用其他 dataset 作为验证集。
# profile:
#   行为 feature，列为 sample 的丰度表。
# group:
#   样本分组信息表，至少包含 sample_col、group_col 和 dataset_col。
#   默认列名为 sample / group / dataset。
# dataset_col:
#   指定数据集/队列所在列名，默认 'dataset'。
# dataset_list:
#   指定需要参与交叉验证的数据集名称。
#   默认 NULL，表示使用 group[[dataset_col]] 中所有数据集。
# positive:
#   二分类时用于计算 AUC 的阳性类别。
#   默认 NULL，使用 factor 水平中的第二个类别。
# 返回:
#   train_dataset / test_dataset / seed / ntree / AUC
# 示例:
#   rf_cross_dataset_vaildate(
#     profile = profile, group = group, 
#     dataset_col = 'dataset', positive = 'IBD'
#   )

# 每次使用一个 dataset 建模，并用其他 dataset 验证
# profile: 行为 feature，列为 sample
# group: 样本分组信息，默认包含 sample / group / dataset 三列
rf_cross_dataset_vaildate <- function(profile, group, dataset_col = 'dataset',
                                      dataset_list = NULL, seed = 2026, ntree = 1000,
                                      sample_col = 'sample', group_col = 'group',
                                      positive = NULL, na_fill = 0, ...) {
  
  profile <- data.frame(profile, check.names = FALSE)
  group <- data.frame(group, check.names = FALSE)
  
  if (is.null(dataset_list)) {
    dataset_list <- unique(as.character(group[[dataset_col]]))
  }
  
  rf <- list()
  k <- 1
  
  for (i in dataset_list) {
    ## 1. 训练集
    group_i <- group[group[[dataset_col]] == i, , drop = FALSE]
    sample_i <- intersect(group_i[[sample_col]], colnames(profile))
    group_i <- group_i[match(sample_i, group_i[[sample_col]]), , drop = FALSE]
    
    train_x <- data.frame(t(profile[, sample_i, drop = FALSE]), check.names = FALSE)
    train_y <- factor(group_i[[group_col]])
    
    train_x <- data.frame(
      lapply(train_x, \(x) suppressWarnings(as.numeric(x))),
      check.names = FALSE, row.names = rownames(train_x)
    )
    
    train_x[is.na(train_x)] <- na_fill
    
    ## 去除训练集零方差 feature
    keep <- vapply(train_x, stats::var, numeric(1)) > 0
    train_x <- train_x[, keep, drop = FALSE]
    feature_keep <- colnames(train_x)
    
    if (ncol(train_x) == 0 || nlevels(train_y) < 2) next
    
    colnames(train_x) <- paste0('var_', seq_len(ncol(train_x)))
    
    set.seed(seed)
    rf_model <- randomForest::randomForest(
      x = train_x, y = train_y, ntree = ntree,
      importance = FALSE, proximity = FALSE, ...
    )
    
    ## 2. 验证集
    for (j in setdiff(dataset_list, i)) {
      
      group_j <- group[group[[dataset_col]] == j, , drop = FALSE]
      sample_j <- intersect(group_j[[sample_col]], colnames(profile))
      group_j <- group_j[match(sample_j, group_j[[sample_col]]), , drop = FALSE]
      
      test_x <- data.frame(t(profile[feature_keep, sample_j, drop = FALSE]), check.names = FALSE)
      test_y <- factor(group_j[[group_col]], levels = levels(train_y))
      
      test_x <- data.frame(
        lapply(test_x, \(x) suppressWarnings(as.numeric(x))),
        check.names = FALSE, row.names = rownames(test_x)
      )
      
      test_x[is.na(test_x)] <- na_fill
      colnames(test_x) <- colnames(train_x)
      
      pred <- data.frame(
        stats::predict(rf_model, test_x, type = 'prob'),
        check.names = FALSE
      )
      
      auc_value <- NA_real_
      
      if (nlevels(train_y) == 2 && length(unique(na.omit(test_y))) == 2) {
        
        pos <- ifelse(is.null(positive), levels(train_y)[2], positive)
        
        if (pos %in% colnames(pred)) {
          roc <- tryCatch(
            pROC::roc(test_y, pred[[pos]], levels = levels(train_y),
                      direction = '<', quiet = TRUE),
            error = function(e) NULL
          )
          
          if (!is.null(roc)) {
            auc_value <- as.numeric(pROC::auc(roc))
          }
        }
      }
      
      rf[[k]] <- data.frame(
        train_dataset = i,
        test_dataset = j,
        seed = seed,
        ntree = ntree,
        AUC = auc_value
      )
      
      k <- k + 1
    }
  }
  
  return(dplyr::bind_rows(rf))
}

#### 随机森林变量重要性排序 ####
# 使用全部样本训练随机森林模型，并输出 feature importance。
# profile:
#   行为 feature，列为 sample 的丰度表。
# group:
#   样本分组信息表，默认列名为 sample 和 group。
# seed:
#   随机种子。
# ntree:
#   随机森林树的数量。
# na_fill:
#   profile 中 NA 的填充值，默认 0。
# 输出:
#   name:
#     原始 feature 名。
#   MeanDecreaseAccuracy:
#     置换重要性，越高说明该 feature 对分类准确率贡献越大。
#   MeanDecreaseGini:
#     Gini 重要性，越高说明该 feature 在树分裂中贡献越大。
# 示例:
#   imp <- rf_importance(
#     profile = profile, group = group, ntree = 1000
#   )

rf_importance <- function(profile, group, seed = 2026, ntree = 1000,
                          sample_col = 'sample', group_col = 'group',
                          na_fill = 0, ...) {
  
  profile <- data.frame(profile, check.names = FALSE)
  group <- data.frame(group, check.names = FALSE)
  
  sample_x <- intersect(group[[sample_col]], colnames(profile))
  group <- group[match(sample_x, group[[sample_col]]), , drop = FALSE]
  
  X <- data.frame(t(profile[, sample_x, drop = FALSE]), check.names = FALSE)
  y <- factor(group[[group_col]])
  
  X <- data.frame(
    lapply(X, \(x) suppressWarnings(as.numeric(x))),
    check.names = FALSE,
    row.names = rownames(X)
  )
  
  X[is.na(X)] <- na_fill
  
  keep <- vapply(X, stats::var, numeric(1)) > 0
  X <- X[, keep, drop = FALSE]
  
  data_rename <- data.frame(
    name = colnames(X),
    rename = paste0('var_', seq_len(ncol(X)))
  )
  
  colnames(X) <- data_rename$rename
  
  set.seed(seed)
  rf_model <- randomForest::randomForest(
    x = X, y = y, ntree = ntree,
    importance = TRUE, proximity = FALSE, ...
  )
  
  importance <- randomForest::importance(rf_model) |>
    data.frame(check.names = FALSE) |>
    tibble::rownames_to_column('feature') |>
    dplyr::mutate(name = data_rename$name[match(feature, data_rename$rename)]) |>
    dplyr::select(name, dplyr::everything(), -feature)
  
  if ('MeanDecreaseAccuracy' %in% colnames(importance)) {
    importance <- importance |>
      dplyr::arrange(dplyr::desc(MeanDecreaseAccuracy))
  }
  
  return(importance)
}

#### 随机森林留一法 ####
# Leave-one-out random forest。
# 每次留出 1 个样本作为测试集，其余样本作为训练集。
# profile:
#   行为 feature，列为 sample 的丰度表。
# group:
#   样本分组信息表，默认列名为 sample 和 group。
# progress:
#   是否显示进度条，默认 TRUE。
# 返回:
#   每个样本的随机森林预测概率。
#   第一列为 sample，后续列为各类别的预测概率。
# 注意:
#   如果某一次留一后训练集中只剩一个类别，该样本预测概率会返回 NA。
# 示例:
#   pred <- rf_loom(
#     profile = profile, group = group, ntree = 1000
#   )

rf_loom <- function(profile, group, seed = 2026, ntree = 1000,
                    sample_col = 'sample', group_col = 'group',
                    na_fill = 0, progress = TRUE, ...) {
  
  profile <- data.frame(profile, check.names = FALSE)
  group <- data.frame(group, check.names = FALSE)
  
  sample_x <- intersect(group[[sample_col]], colnames(profile))
  group <- group[match(sample_x, group[[sample_col]]), , drop = FALSE]
  
  X <- data.frame(t(profile[, sample_x, drop = FALSE]), check.names = FALSE)
  y <- factor(group[[group_col]])
  
  X <- data.frame(
    lapply(X, \(x) suppressWarnings(as.numeric(x))),
    check.names = FALSE,
    row.names = rownames(X)
  )
  
  X[is.na(X)] <- na_fill
  
  keep <- vapply(X, stats::var, numeric(1)) > 0
  X <- X[, keep, drop = FALSE]
  
  colnames(X) <- paste0('var_', seq_len(ncol(X)))
  
  pred <- list()
  
  if (isTRUE(progress)) {
    pb <- txtProgressBar(style = 3, width = 50, char = '#')
    on.exit(close(pb), add = TRUE)
  }
  
  for (i in seq_len(nrow(X))) {
    
    train_x <- X[-i, , drop = FALSE]
    test_x <- X[i, , drop = FALSE]
    train_y <- y[-i]
    
    if (length(unique(train_y)) < 2) {
      pred_i <- data.frame(matrix(NA_real_, nrow = 1, ncol = nlevels(y)))
      colnames(pred_i) <- levels(y)
      rownames(pred_i) <- rownames(test_x)
    } else {
      set.seed(seed + i)
      rf_model <- randomForest::randomForest(
        x = train_x, y = train_y, ntree = ntree,
        importance = FALSE, proximity = FALSE, ...
      )
      
      pred_i <- data.frame(
        predict(rf_model, test_x, type = 'prob'),
        check.names = FALSE
      )
      
      pred_i <- pred_i[, levels(y), drop = FALSE]
    }
    
    pred[[i]] <- pred_i |>
      tibble::rownames_to_column('sample')
    
    if (isTRUE(progress)) {
      setTxtProgressBar(pb, i / nrow(X))
    }
  }
  
  return(dplyr::bind_rows(pred))
}

#### profile_x建模，profile_y验证 ####
# 使用 profile_x 对应样本训练随机森林模型，
# 然后预测 profile_y 对应样本的类别概率。
# profile_x:
#   训练集 profile，行为 feature，列为 sample。
# profile_y:
#   验证集 profile，行为 feature，列为 sample。
# group_x:
#   训练集样本分组信息，默认列名为 sample 和 group。
# group_y:
#   验证集样本分组信息，默认列名为 sample 和 group。
# 处理逻辑:
#   1. 自动取 profile_x 和 profile_y 的共有 feature；
#   2. 自动匹配 group_x/group_y 中存在于 profile 的样本；
#   3. 使用 profile_x 建模；
#   4. 输出 profile_y 样本的预测概率。
# 返回:
#   第一列为 sample，后续列为训练集中各类别的预测概率。
# 示例:
#   pred <- rf_next_vaildate(
#     profile_x = discovery_profile, profile_y = validation_profile,
#     group_x = discovery_group, group_y = validation_group
#   )

rf_next_vaildate <- function(profile_x, profile_y, group_x, group_y,
                             seed = 2026, ntree = 1000,
                             sample_col = 'sample', group_col = 'group',
                             na_fill = 0, ...) {
  
  shared_name <- intersect(rownames(profile_x), rownames(profile_y))
  
  profile_x <- profile_x[shared_name, , drop = FALSE]
  profile_y <- profile_y[shared_name, , drop = FALSE]
  
  sample_x <- intersect(group_x[[sample_col]], colnames(profile_x))
  sample_y <- intersect(group_y[[sample_col]], colnames(profile_y))
  
  group_x <- group_x[match(sample_x, group_x[[sample_col]]), , drop = FALSE]
  group_y <- group_y[match(sample_y, group_y[[sample_col]]), , drop = FALSE]
  
  train_x <- data.frame(t(profile_x[, sample_x, drop = FALSE]), check.names = FALSE)
  test_x <- data.frame(t(profile_y[, sample_y, drop = FALSE]), check.names = FALSE)
  train_y <- factor(group_x[[group_col]])
  
  train_x <- data.frame(
    lapply(train_x, \(x) suppressWarnings(as.numeric(x))),
    check.names = FALSE,
    row.names = rownames(train_x)
  )
  
  test_x <- data.frame(
    lapply(test_x, \(x) suppressWarnings(as.numeric(x))),
    check.names = FALSE,
    row.names = rownames(test_x)
  )
  
  train_x[is.na(train_x)] <- na_fill
  test_x[is.na(test_x)] <- na_fill
  
  keep <- vapply(train_x, stats::var, numeric(1)) > 0
  train_x <- train_x[, keep, drop = FALSE]
  test_x <- test_x[, keep, drop = FALSE]
  
  colnames(train_x) <- paste0('var_', seq_len(ncol(train_x)))
  colnames(test_x) <- colnames(train_x)
  
  set.seed(seed)
  rf_model <- randomForest::randomForest(
    x = train_x, y = train_y, ntree = ntree,
    importance = FALSE, proximity = FALSE, ...
  )
  
  pred <- data.frame(
    predict(rf_model, test_x, type = 'prob'),
    check.names = FALSE
  ) |>
    dplyr::select(dplyr::all_of(levels(train_y))) |>
    tibble::rownames_to_column('sample')
  
  return(pred)
}
