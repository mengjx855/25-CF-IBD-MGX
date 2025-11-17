#### Jinxin Meng, 202205029, 20251106, v1.3.2 ####
# 20250917: update function. rename 'rf_feature_rank' as 'rf_importance'
# 20251106: update rf_next_vaildate

library(dplyr)
library(tibble)
library(tidyr)
library(randomForest)

#### K折随机森林建模 ####
# K折随机森林建模，返回一个数据框，输出一个预测数据框
# 输入一个profile表; group文件是样本的分组信息，分组的列名为sample和group
# k为几折检验; # seed设置随机种子, ntree设置随机森林的树数量
# make.names修改一些复杂feature的名字，主要是预防报错，二者可以降低内存。我这里就是把名字替换了，R里可以直接用make.names()函数去修改某一个表的内容
rf_Kfold <- function(profile, group, k, seed = 2024, ntree = 1000) {
  start_time <- Sys.time()
  profile <- data.frame(t(data.frame(profile, check.names = F)))
  # 随机按照K折采样
  n_sample <- nrow(profile)
  size <- round(n_sample/k, 0)
  count <- seq_len(n_sample)
  sample_result <- list()
  for (i in 1:(k-1)) {
    set.seed(seed)
    sample_i <- sample(x = count, size = size, replace = F)
    sample_result[[i]] <- sample_i
    count <- setdiff(count, sample_i)
  }
  sample_result[[i+1]] <- count
  # 训练和测试
  colnames(profile) <- paste0('V', seq_len(ncol(profile)))
  profile$group <- as.factor(group$group[match(rownames(profile), group$sample)])
  pred <- rbind()
  pb <- txtProgressBar(style = 3, width = 50, char = '#')
  for(j in 1:length(sample_result)){ # 划分测试集和训练集
    sample_j <- sample_result[[j]]
    test <- profile[sample_j,]
    train <- profile[-sample_j,]
    set.seed(seed)
    rf_model <- randomForest(group ~ ., data = train, ntree = ntree, importance = F, proximity = T)
    pred_i <- data.frame(predict(rf_model, test, type = 'prob'))
    pred <- rbind(pred, pred_i)
    setTxtProgressBar(pb, j/length(sample_result)) # 进度计算
  }
  close(pb)
  pred <- rownames_to_column(pred, var = 'sample')
  end_time <- Sys.time()
  run_time <- end_time - start_time
  cat(paste0('  Time cost: ', round(as.numeric(run_time), 4), ' secs.\n'))
  return(pred)
}

#### K折rep次重复 ####
# K折随机森林建模，返回一个数据框，输出一个预测数据框
# profile输入的profile格式的表格, 但是要求行是样本，列是特征，也就是正常profile表格的转置
# group文件是样本的分组信息，分组的列名为sample和group
# k为几折检验
# rep为几次重复，最后会计算均值
# seed设置随机种子, ntree设置随机森林的树数量
# make.names修改一些复杂feature的名字，主要是预防报错，二者可以降低内存。我这里就是把名字替换了，R里可以直接用make.names()函数去修改某一个表的内容
rf_Kfold_rep <- function(profile, group, k, rep = 5, seed = 2023, ntree = 1000) {
  # rep次重复
  pred <- rbind()
  for (i in seq_len(rep)) {
    cat(paste0('  [', Sys.time(), '] Repeat no.', i, ' time.\n'))
    pred_i <- rf_Kfold(profile, group, k = k, seed = seed + i -1, ntree = ntree)
    pred_i <- add_column(pred_i, rep = paste0('rep_', i))
    pred <- rbind(pred, pred_i)
  }
  pred <- pred %>% 
    select(-rep) %>% 
    group_by(sample) %>% 
    summarise_all(mean)
  return(pred)
}

#### 多数据集交叉验证 ####
# 多数据集交叉验证
# 原理上，就是每一个数据集进行建模，然后用其他数据集进行验证
# 输入profile表，正常的profile表格列为样本，行为feature
# group文件是样本的分组信息，分组的列名为sample/group/dataset_name,第一列必须为'sample'，第二列必须为'group', dataset_name这列是可以指定名称的
# seed设置随机种子, ntree设置随机森林的树数量
# 指定交叉验证数据集的向量
# 返回一个交叉验证结果的数据框
rf_cross_dataset_vaildate <- function(profile, group, dataset_name = 'dataset', 
                                      dataset_list = NULL, seed = 2022, ntree = 1000) {
  # 提取多个数据集的list
  if(is.null(dataset_name)) {
    dataset_list <- group %>% 
      select(all_of(dataset_name)) %>%
      unlist() %>%
      as.character() %>%
      unique()
  } else {
    dataset_list = dataset_list
  }
  rf <- rbind() # VKH使用了BD的健康人样本
  for (i in dataset_list) {
    # 训练
    group_i <- group[group[[dataset_name]]%in%i,]
    sample_i <- group_i %>% 
      select(sample) %>% 
      data.frame(check.names = F) %>% 
      unlist() %>% 
      as.character()
    profile_i <- profile %>% 
      select(all_of(sample_i)) %>% 
      t() %>% 
      data.frame(check.names = F)
    colnames(profile_i) <- paste0('var_', seq(ncol(profile_i)))
    profile_i$group <- factor(group_i$group[match(rownames(profile_i), group_i$sample)])
    set.seed(seed)
    rf_model <- randomForest(group ~ ., data = profile_i, ntree = ntree, importance = F, proximity = T)
    # 预测
    for (j in setdiff(dataset_list, i)) {
      group_j <- group[group[[dataset_name]]%in%j,]
      sample_j <- group_j %>% 
        select(sample) %>% 
        data.frame(check.names = F) %>% 
        unlist() %>% 
        as.character()
      profile_j <- profile %>% 
        select(all_of(sample_j)) %>% 
        t() %>% 
        data.frame(check.names = F)
      colnames(profile_j) <- paste0('var_', seq(ncol(profile_j)))    
      profile_j$group <- factor(group_j$group[match(rownames(profile_j), group_j$sample)])
      pred <- predict(rf_model, profile_j, type = 'prob') %>% 
        data.frame() %>% 
        rownames_to_column(var = 'sample') %>% 
        mutate(group = group_j$group[match(.$sample, group_j$sample)])
      roc <- roc(pred$group, pred[,2])
      rf <- rbind(rf, data.frame(train_dataset = i, test_dataset = j, seed = seed, ntree = ntree, AUC = as.numeric(str_replace(roc$auc, 'Area under the curve: ', ''))))
    }
  }
  return(rf)
}

#### 随机森林变量重要性排序 ####
rf_importance <- function(profile, group, seed = 2023, ntree = 1000) {
  profile <- data.frame(t(profile), check.names = F)
  data_rename <- data.frame(name = colnames(profile), rename = paste0('n_', seq_len(ncol(profile))))
  colnames(profile) <- data_rename$rename
  profile$group <- factor(group$group[match(rownames(profile), group$sample)])
  # 建模
  set.seed(seed)
  rf_model <- randomForest(group ~ ., data = profile, ntree = ntree, importance = T, proximity = T)
  # 变量重要性
  importance <- importance(rf_model) %>% 
    data.frame(check.names = F) %>% 
    tibble::rownames_to_column('feature') %>%
    mutate(name = data_rename$name[match(feature, data_rename$rename)]) %>% 
    dplyr::select(name, MeanDecreaseAccuracy, MeanDecreaseGini) %>%
    dplyr::arrange(desc(MeanDecreaseAccuracy))
  return(importance)
}

#### 随机森林留一法 ####
# leave-one-out method
# K折随机森林建模，返回一个数据框，输出一个预测数据框
# profile_t输入的profile格式的表格,t(profile)
# group文件是样本的分组信息，分组的列名为group
# k为几折检验
# seed设置随机种子
rf_loom <- function(profile, group, seed = 2022, ntree = 1000){
  profile <- select(profile, all_of(group$sample))
  set.seed(seed)
  profile <- data.frame(t(profile), check.names = F)
  colnames(profile) <- paste0('var_', seq_len(ncol(profile))) 
  profile$group <- as.factor(group$group[match(rownames(profile),group$sample)])
  # 新建进度条
  pb <- txtProgressBar(style = 3, width = 50, char = '#')
  start_time <- Sys.time()
  pred <- rbind()
  for(i in 1:length(rownames(profile))){
    fold_test <- profile[i,]   # 取1个样品作为测试集 
    fold_train <- profile[-i,]   # 剩下的样品作为训练集
    rf_model <- randomForest(group ~ ., data = fold_train, ntree = ntree, importance = T, proximity = TRUE)
    temp <- data.frame(predict(rf_model, fold_test, type = 'prob'))
    pred <- rbind(pred, temp) 
    setTxtProgressBar(pb, i/length(rownames(profile))) # 进度计算
  }
  pred <- rownames_to_column(pred, 'sample')
  end_time <- Sys.time()
  close(pb)
  run_time <- end_time - start_time
  cat(paste0('Run times: ', round(run_time, 4), '\n'))
  return(pred)
}

#### profile_x建模，profile_y验证 ####
# 两个profile表均是转置后的，第一个是发现集，第二个为验证集
# group文件是样本的分组信息，分组的列至少有sample和group,第一列必须为'sample'，第二列必须为'group'
# seed设置随机种子, ntree设置随机森林的树数量
# 返回一个预测结果的数据框
rf_next_vaildate <- function(profile_x, profile_y, group_x, group_y, seed = 2024, 
                             ntree = 1000, label = 'profile_x for modeling and profile_y for predicting'){
  shared_name <- intersect(rownames(profile_x), rownames(profile_y))
  profile_x <- profile_x[shared_name, ]
  profile_y <- profile_y[shared_name, ]
  
  # discovery
  profile_x <- data.frame(t(profile_x))
  profile_x$group <- as.factor(group_x$group[match(rownames(profile_x), group_x$sample)])
  set.seed(seed)
  rf_model <- randomForest(group ~ ., data = profile_x, ntree = ntree, importance = F, proximity = T)
  
  # validation
  profile_y <- data.frame(t(profile_y))
  profile_y$group <- as.factor(group_y$group[match(rownames(profile_y), group_y$sample)])
  pred <- predict(rf_model, profile_y, type = 'prob') %>% 
    data.frame() %>% 
    rownames_to_column(var = 'sample')
  return(pred)
}