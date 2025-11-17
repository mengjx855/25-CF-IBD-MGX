#### Jinxin Meng, 20231215, 20250623, v1.0.1 ####
# procrustes分析|普鲁克分析|两个矩阵的相关性分析
# Procrustes分析与Mantel Test的比较
# 提到群落分析，Procrustes分析与Mantel test都是用于分析物种组成和环境属性关系的常见方法，当然二者的具体关注点还是有区别的，方法各有自身的优点。
# 但如果只聚焦在评估两数据集一致性上，似乎Procrustes分析更直观一些。
# M2统计量及其显著性检验p值提供了两个数据集之间一致性的总体度量，同时数据集的图形匹配和相关残差提供了比Mantel test更丰富的信息源。
# 在对应点的坐标匹配度较好时，两个数据集表现出良好的一致性。坐标匹配度越差表明这些点与整体趋势不匹配，
# 这类似于回归分析中残差较大的点，这些点不符合样本的总体趋势。
# 此外，PROTEST的统计功效也被证明优于Mantel test的统计功效（Peres-Neto and Jackson, 2001）。
# 因此，如果两组数据之间存在潜在关系，则Procrustes分析更有能力检测到，并且鉴于结果的图形性质，它还提供了出色的解释性准则。

library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(vegan)

#### plot_Procrustes ####
# input: profile[feature, sample] | dist 
# 输入两个profile表，或者是基于两个profile表分析的聚类矩阵。 默认是 y 旋转变换到 x
plot_Procrustes <- function(profile_x, profile_y, dist_x, dist_y, dis_method = 'bray',
                            symmetric = TRUE, colors = c(x = '#9BBB59', y = '#957DB1'),
                            xlab = NULL, ylab = NULL, title = NULL, subtitle = NULL, 
                            show_grid = FALSE, show_line = TRUE, aspect_ratio = 3/4, 
                            theme = 'default', ...) {
  if (!missing(profile_x) & !missing(profile_y)) {
    dist_x <- vegdist(t(profile_x), method = dis_method)
    dist_y <- vegdist(t(profile_y), method = dis_method)
  } else if (missing(dist_x) & missing(dist_y)) {
    stop('distance matrix missing!')
  }
  
  # 降维分析
  PCoA_x <- cmdscale(dist_x)
  PCoA_y <- cmdscale(dist_y)
  
  # 选择模式：通过对该函数中各参数的了解，可知X为目标矩阵也就是降维后的环境（功能基因等）坐标，Y为降维后的物种数据的坐标，因为后续普氏分析中旋转和缩放操作是针对Y，将Y匹配给X。
  # 另外，当symmetric=FALSE时，处于'非对称'模式，X和Y的分配值调换后，普氏分析的偏差平方和（M^2）也会随之改变。而当symmetric=TRUE时，从而给出更合适比例的对称统计。“对称”模式下，X和Y的分配值调换后，普氏分析的偏差平方和（M2）不会发生改变，但注意旋转仍将是非对称的。# 以对称模式进行普氏分析（symmetric = TRUE）
  # 作者回应 https://stats.stackexchange.com/questions/563911/what-is-the-difference-between-symmetric-and-non-symmetric-in-procrustes-protest
  # 将矩阵B拟合到目标矩阵A：我想获得与A类似的旋转/缩放/平移矩阵B。如果这是您的目标，您应该使用非对称(B到A)旋转。 
  # 如果您不想获得结果（旋转），但您只是对两个矩阵之间的相似性以及该相似性的统计量感兴趣，则应该使用对称旋转。 
  proc <- procrustes(PCoA_x, PCoA_y, symmetric = symmetric)
  
  # summary(proc) # 评价一致性
  # 如果样本中物种与环境一致性（相似性）越近，则对应的残差越小，反之物种与环境的相似性越远，则残差越大（三条辅助线对应的位置分别为残差25%、50%和75%）
  # plot(proc, kind = 2)
  # residuals(proc)
  
  # 普氏分析中M2统计量的显著性检验，置换检验999次。在 Procrustes 分析中，M2 常常指代的是 Procrustes 统计量，它是度量两个形状间差异程度的一个指标。
  # 更具体地说，M2 是源数据集的点经过旋转、缩放和/或平移后与目标数据集中对应点的平方距离和。
  # 它代表了变换后源数据集的点与目标数据集点之间的不匹配程度。M2 的数值越小，表明两组数据集的形状越相似；反之，则表明它们之间的差异越大。
  set.seed(2024)
  proc_test <- protest(PCoA_x, PCoA_y, permutations = 999)
  
  if (is.null(xlab)) 
    xlab <- 'Dim 1'
  
  if (is.null(ylab)) 
    ylab <- 'Dim 2'
  
  if (is.null(title))
    title <- paste0(stringr::str_to_sentence(string = dis_method), 
                    ' distance-based Procrustes analysis')
  
  if (is.null(subtitle))
    subtitle <- substitute("coefficients: " * M^2 == a ~ ", " ~ italic(p) < b, 
                           list(a = round(proc_test$ss, 4), b = proc_test$signif))
  
  # 首先提取降维后的数据轴1和2的坐标，并且提取转换的坐标；然后进行绘制。获得x和y轴的坐标及旋转过的坐标
  proc_point <- cbind(
    # Y-矩阵旋转后的坐标，就是物种矩阵为了逼近X做了调整，调整后的坐标。
    dplyr::rename(data.frame(proc$Yrot), X1_rotated = X1, X2_rotated = X2), 
    # X-目标矩阵，就是proc分析输入的第一个矩阵，作为目标矩阵没变化。
    dplyr::rename(data.frame(proc$X), X1_target = Dim1, X2_target = Dim2))
  proc_coord <- data.frame(proc$rotation) # Y旋转后的坐标轴
  
  # 绘图
  p <- ggplot(proc_point) + # 旋转坐标到目的坐标的一半上一个颜色，另一半上不同的颜色
    geom_segment(aes(x = X1_rotated, y = X2_rotated, # 旋转矩阵
                     xend = (X1_rotated + X1_target)/2, 
                     yend = (X2_rotated + X2_target)/2), 
                 arrow = arrow(length = unit(0, 'cm')),
                 color = colors[1], linewidth = .4) + 
    geom_segment(aes(x = (X1_rotated + X1_target)/2, # 目标矩阵
                     y = (X2_rotated + X2_target)/2, 
                     xend = X1_target, yend = X2_target), 
                 arrow = arrow(length = unit(0.15, 'cm')), 
                 color = colors[2], linewidth = .4) +
    geom_point(aes(X1_rotated, X2_rotated), color = colors[1], # 旋转坐标点
               size = 1.6, shape = 16) + 
    geom_point(aes(X1_target, X2_target), color = colors[2], # 目的坐标点
               size = 1.6, shape = 16) + 
    geom_abline(intercept = 0, slope = proc_coord[1,2]/proc_coord[1,1], linewidth = .5) +
    geom_abline(intercept = 0, slope = proc_coord[2,2]/proc_coord[2,1], linewidth = .5) +
    labs(x = xlab, y = ylab, title = title, subtitle = subtitle)
  
  if (isTRUE(show_line)) {
    p <- p + 
      geom_vline(xintercept = 0, color = 'gray', linetype = 'longdash', linewidth = .5) + 
      geom_hline(yintercept = 0, color = 'gray', linetype = 'longdash', linewidth = .5)
  }
  
  if (theme == 'pubr') { # 主题为 pubr
    p <- p + 
      ggpubr::theme_pubr() + 
      theme(aspect.ratio = aspect_ratio, 
            plot.margin =  unit(c(2, 2, 2, 2), 'mm'),
            plot.title = element_text(hjust = .5, size = 12, face = 'bold'),
            legend.position = 'right')
  }
  
  if (theme != 'pubr') { # 主题为 默认
    p <- p + 
      theme_bw() +
      theme(axis.ticks = element_line(linewidth = .5, color = 'black'),
            axis.ticks.length = unit(2, 'mm'),
            axis.title = element_text(size = 12, color = 'black'),
            axis.text = element_text(size = 12, color = 'black'),
            axis.line = element_blank(),
            plot.title = element_text(hjust = .5, size = 12, face = 'bold'),
            plot.subtitle = element_text(size = 12, color = 'black'),
            plot.margin =  unit(c(2, 2, 2, 2), 'mm'),
            panel.border = element_rect(linewidth = .5, color = 'black', fill = NA),
            panel.background = element_blank(),
            panel.grid = element_blank(), 
            legend.background = element_blank(),
            legend.text = element_text(size = 10, color = 'black'),
            legend.title = element_text(size = 10, color = 'black'),
            aspect.ratio = aspect_ratio )
  }
  
  if (isTRUE(show_grid)) {
    p <- p + 
      theme(panel.grid.major = element_line(linewidth = .5, color = 'grey90'),
            panel.grid.minor = element_blank())
  }
  
  message('ggsave(file = \'procrustes_scatterplot.pdf\', width = 6, height = 4.5)') 
  # profile_x is target matrix (arrow direction) [colors-1]
  # profile_y is matrix to be rotated. [colors-2]
  return(p)
}
