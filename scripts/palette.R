#### jinxin meng, mengjx855@163.com, 20240304, 20250727, v.2.2 ####
# 20250418: update function
# 20250727: update function

library(paletteer)

#### pal_sets ####
.pal_sets <- list(
  npg =       c('#e64b35','#4dbbd5','#00a087','#3c5488','#f39b7f','#8491b4','#91d1c2','#dc0000','#7e6148','#b09c85'),
  aaas =      c('#3b4992','#ee0000','#008b45','#631879','#008280','#bb0021','#5f559b','#a20056','#808180','#1b1919'),
  nejm =      c('#bc3c29','#0072b5','#e18727','#20854e','#7876b1','#6f99ad','#ffdc91','#ee4c97'),
  lancet =    c('#00468b','#ed0000','#42b540','#0099b4','#925e9f','#fdaf91','#ad002a','#adb6b6','#1b1919'),
  jama =      c('#374e55','#df8f44','#00a1d5','#b24745','#79af97','#6a6599','#80796b'),
  jco =       c('#0073c2','#efc000','#868686','#cd534c','#7aa6dc','#003c67','#8f7700','#3b3b3b','#a73030'),
  usscgb =    c('#ff0000','#ff9900','#ffcc00','#00ff00','#6699ff','#cc33ff','#99991e','#999999','#ff00cc','#cc0000','#ffcccc','#ffff00','#ccff00','#358000','#0000cc','#99ccff','#00ffff','#ccffff','#9900cc','#cc99ff','#996600','#666600','#666666','#cccccc','#79cc3d','#cccc99'),
  d3 =        c('#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b','#e377c2','#7f7f7f','#bcbd22','#17becf','#aec7e8','#ffbb78','#98df8a','#ff9896','#c5b0d5','#c49c94','#f7b6d2','#c7c7c7','#dbdb8d','#9edae5'),
  locuszoom = c('#d43f3a','#eea236','#5cb85c','#46b8da','#357ebd','#9632b8','#b8b8b8'),
  uchicago =  c('#800000','#767676','#ffa319','#8a9045','#155f83','#c16622','#8f3931','#58593f','#350e20'),
  startrek =  c('#cc0c00','#5c88da','#84bd00','#ffcd00','#7c878e','#00b5e2','#00af66'),
  tron =      c('#ff410d','#6ee2ff','#f7c530','#95cc5e','#d0dfe6','#f79d1e','#748aa6'),
  futurama =  c('#ff6f00','#c71000','#008ea0','#8a4198','#5a9599','#ff6348','#84d7e1','#ff95a8','#3d3b25','#ade2d0','#1a5354','#3f4041'),
  simpsons =  c('#fed439','#709ae1','#8a9197','#d2af81','#fd7446','#d5e4a2','#197ec0','#f05c3b','#46732e','#71d0f5','#370335','#075149','#c80813','#91331f','#1a9993','#fd8cc1'),
  igv =       c('#5773cc','#ffb900'),
  Rainbow =   c('#FF0000','#FF4D00','#FF9900','#FFE500','#CCFF00','#80FF00','#33FF00','#00FF19','#00FF66','#00FFB2','#00FFFF','#00B3FF','#0066FF','#001AFF','#3300FF','#7F00FF','#CC00FF','#FF00E6','#FF0099','#FF004D'),
  Rainbow2 =  c('#f04e40','#f26335','#f3732d','#f68424','#f99d28','#fdb82a','#f6cc34','#e3d63c','#d4e04c','#b1d353','#8fc958','#6fc162','#60bf88','#39beab','#1db6c7','#2899ce','#3b7cbf','#3769b3','#2e56a6','#22439c'),
  Rainbow3 =  c('#f3732d','#f68424','#f99d28','#fdb82a','#f6cc34','#e3d63c','#d4e04c','#b1d353','#8fc958','#6fc162','#60bf88','#39beab','#1db6c7','#2899ce','#3b7cbf','#3769b3'),
  Rainbow4 =  c('#f99d28','#fdb82a','#f6cc34','#e3d63c','#d4e04c','#b1d353','#8fc958','#6fc162','#60bf88','#39beab','#1db6c7','#2899ce'),
  Spectral =  c('#9E0142','#D53E4F','#F46D43','#FDAE61','#FEE08B','#FFFFBF','#E6F598','#ABDDA4','#66C2A5','#3288BD','#5E4FA2'),
  Set2 =      c('#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f','#e5c494','#b3b3b3'),
  hue =       c('#f8766d','#ea8331','#d89000','#c09b00','#a3a500','#7cae00','#39b600','#00bb4e','#00bf7d','#00c1a3','#00bfc4','#00bae0','#00b0f6','#35a2ff','#9590ff','#c77cff','#e76bf3','#fa62db','#ff62bc','#ff6a98'),
  Set3 =      c('#80b1d3','#b3de69','#fdb462','#8dd3c7','#bc80bd','#fb8072','#ffed6f','#fccde5','#bebada','#ccebc5','#ffffb3','#d9d9d9'),
  Accent =    c('#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0','#f0027f','#bf5b17','#666666'),
  Paired =    c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928'),
  Paired2 =   c('#4e79a7','#a0cbe8','#f28e2b','#ffbe7d','#59a14f','#8cd17d','#b6992d','#f1ce63','#499894','#86bcb6','#e15759','#ff9d9a','#79706e','#bab0ac','#d37295','#fabfd2','#b07aa1','#d4a6c8','#9d7660','#d7b5a6'),
  Contrast2 = c('#56B4E9','#E69F00'),
  Contrast3 = c('#f46d43','#74add1'),  # 多种免疫病病毒组研究中病例和对照组用的颜色
  Contrast4 = c('#f77d4d','#1f78b4'),  # MDD, IBS, ACVD 病例对照中的颜色方案
  Contrast5 = c('#ffb900','#5773cc'),
  ReYebl =    c('#d53e4f','#fc8d59','#fee08b','#ffffbf','#e6f598','#99d594','#3288bd'), # 渐变色 红-黄-蓝
  PuWhGr =    c('#c51b7d','#e9a3c9','#fde0ef','#f7f7f7','#e6f5d0','#a1d76a','#4d9221'), # 渐变色 紫-白-绿
  ReWhBl =    c('#b2182b','#ef8a62','#fddbc7','#f7f7f7','#d1e5f0','#67a9cf','#2166ac'), # 渐变色 红-白-蓝
  BrWhCy =    c('#8c510a','#d8b365','#f6e8c3','#f5f5f5','#c7eae5','#5ab4ac','#01665e')  # 渐变色 褐-白-青
  )

#### paette_discrete ####
# name: 预设的颜色集, n: 输出颜色的数量, 不指定参数默认输出颜色集的名称
pald <- function(name = NULL, reverse = F, n = NULL, paste = FALSE, ...) {
  if (is.null(name) & is.null(n)) {
    vars <- map2_vec(.pal_sets, names(.pal_sets), ~ paste0(.y, '-', length(.x))) %>% 
      paste(collapse = ', ')
    return(vars) 
    }
  
  if (!is_null(name) & !(name %in% names(.pal_sets))) 
    stop('Palette name not found.')
  
  if (is.null(n))
    x <- .pal_sets[[name]]
    
  if (!is.null(n)) {
    x <- .pal_sets[[name]]
    .len <- length(x)
    x <- rep(x, times = ceiling(n/.len))[1:n]
  }
  
  if (isTRUE(reverse)) x <- rev(x)
  
  if (isTRUE(paste)) cat(paste0("'", paste(x, collapse = "','"), "'"))
  
  if (isFALSE(paste)) return(x)
}

#### paette_discrete_show ####
pald_show <- function(name = NULL, n = NULL, ...) {
  if (is.null(name) & is.null(n)) {
    vars <- map2_vec(.pal_sets, names(.pal_sets), ~ paste0(.y, '-', length(.x))) %>% 
      paste(collapse = ', ')
    return(vars) 
  }
  
  .color <- pald(name = name, n = n, paste = FALSE)
  if (!is.null(name)) 
    return(scales::show_col(.color))
}

#### palette_continuous ####
# palc
# name: 预设的颜色集
# n: 输出颜色的数量
# 不指定参数默认输出颜色集的名称
palc <- function(name = NULL, n = NULL, paste = FALSE, ...) {
  .pallete_names <- c('BlWhRe','BrWhCy','YeWhBl','ReWhBl','ReYeBl','ReWhGr','Rainbow','Spectral')
  
  if (is.null(name) & is.null(n)) 
    return(paste0('😊 Palette name as following: ', 
                  paste(.pallete_names, collapse = ', ')))
  
  if (!(!is.null(name) & !is.null(n))) 
    stop('Both name and n must be specified simultaneously')
  
  if (!is_null(name) & !(name %in% .pallete_names)) 
    stop('Palette name not found.')
  
  if (name == 'hue') 
    x <- scales::hue_pal()(n)
  if (name == 'YeBr') 
    x <- colorRampPalette(c('#F3DA7D','#A13043'))(n) # 渐变色 黄-褐if (name == 'BlWhRe') x <- colorRampPalette(c('#3288bd', '#ffffff', '#d53e4f'))(n) # 渐变色 蓝-白-红
  if (name == 'BrWhCy') 
    x <- colorRampPalette(c('#a6611a','#dfc27d','#f5f5f5','#80cdc1','#018571'))(n) # 渐变色 褐-白-青
  if (name == 'YeWhBl') 
    x <- colorRampPalette(c('#e66101','#fdb863','#f7f7f7','#b2abd2','#5e3c99'))(n) # 渐变色 黄-白-紫
  if (name == 'ReWhBl') 
    x <- colorRampPalette(c('#ca0020','#f4a582','#f7f7f7','#92c5de','#0571b0'))(n) # 渐变色 红-白-蓝
  if (name == 'ReYeBl') 
    x <- colorRampPalette(c('#d7191c','#fdae61','#ffffbf','#abd9e9','#2c7bb6'))(n) # 渐变色 红-黄-蓝
  if (name == 'ReWhGr') 
    x <- colorRampPalette(c('#f46d43','#fee08b','#ffffff','#d9ef8b','#66bd63'))(n) # 渐变色 红-白-绿
  if (name == 'Rainbow') 
    x <- colorRampPalette(c('#f3732d','#f68424','#f99d28','#fdb82a','#f6cc34','#e3d63c','#d4e04c','#b1d353','#8fc958','#6fc162','#60bf88','#39beab','#1db6c7','#2899ce','#3b7cbf','#3769b3'))(n)
  if (name == 'Spectral') 
    x <- colorRampPalette(c('#9E0142','#D53E4F','#F46D43','#FDAE61','#FEE08B','#FFFFBF','#E6F598','#ABDDA4','#66C2A5','#3288BD','#5E4FA2'))(n)

  return(x)
}

#### palette_continuous_show ####
palc_show <- function(name = NULL, n = 2, ...) {
  .color <- palc(name = name, n = n, paste = FALSE)
  if (!is.null(name)) 
    return(scales::show_col(.color))
}
