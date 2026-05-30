#### Jin-Xin Meng, 20220529, 20260509, v0.0.2 ####

# 20260509: updated functions.

rf_Kfold <- function(profile, group, k = 5, seed = 2025, ntree = 1000,
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
