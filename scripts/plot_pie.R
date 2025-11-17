#### Jinxin Meng, 20231113, 20250816, v0.1.2 ####
# 20250317: update some function.
# 20250816: update some function.
# 20251030: update some function.

library(dplyr)
library(tidyr)
library(ggplot2)

#### plot_pie ####
plot_pie <- function(data, data_rename = NULL, desc = T, level = NULL,
                     top_n = NULL, add_n = F, add_perc = T, lab_cir = F, 
                     lab_cir_flip = F, title = NULL, color = 'white', 
                     fill = 'auto', font_size = 2, hemi = F, start = 0) {
  if (!all(c('name', 'n') %in% colnames(data)) & is.null(data_rename)) 
    stop('data field (name|n)')
  
  if (!is.null(data_rename)) 
    data <- dplyr::rename(data.frame(data, check.names = F), all_of(data_rename))
  
  if (!is.null(level)) 
    data <- mutate(data, name = factor(name, level)) %>% arrange(name)
  
  if (is.null(level) & isTRUE(desc)) 
    data <- arrange(data, desc(n))
  
  if (is.null(level) & isFALSE(desc)) 
    data <- arrange(data, n)
  
  if (!is.null(top_n) & is.numeric(top_n)) {
    .vec <- head(data, top_n - 1) %>% pull(name)
    
    data <- data %>% 
      mutate(name = as.character(name), 
             name = ifelse(name %in% .vec, name, 'Other')) %>% 
      group_by(name) %>% 
      summarise(n = sum(n)) %>% 
      ungroup() %>% 
      mutate(n = round(n, 2)) 
    }
  
  if (isTRUE(desc)) 
    data <- arrange(data, desc(n))
  
  if (isFALSE(desc))
    data <- arrange(data, n)

  data <- mutate(
    data, perc = n / sum(n), 
    ypos = cumsum(perc) - 0.5 * perc,
    angle = ifelse(ypos < .5, 360 * ypos + 180, 360 * ypos)
    )
  
  if (isFALSE(add_n) & isFALSE(add_perc)) 
    data <- mutate(data, lab = name)
  
  if (isTRUE(add_n) & isTRUE(add_perc)) 
    data <- mutate(data, lab = paste0(name, ', ', prettyNum(n, big.mark = ','), ', ', round(perc * 100, 1), '%'))
  
  if (isTRUE(add_n) & isFALSE(add_perc)) 
    data <- mutate(data, lab = paste0(name, ', ', prettyNum(n, big.mark = ',')))
  
  if (isFALSE(add_n) & isTRUE(add_perc)) 
    data <- mutate(data, lab = paste0(name, ', ', round(perc * 100, 1), '%'))
  
  if (length(fill) == 1) {
    if (fill == 'auto') {
      .palette <- c('#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f',
                    '#8dd3c7','#fdb462','#80b1d3','#fccde5','#d9ef8b','#fee391')
      # .palette <- c('#82ccdc','#f0ad80','#b4b3d8','#c2b75f','#87CCBA',
      #               '#F9C851','#ACDEF3','#F9C6B3','#F5EAF2','#D3EDE4')
      .fill <- rep(.palette, time = (ceiling(nrow(data)/length(.palette))))[1:nrow(data)]
    }
    
    if (fill == 'hue') .fill <- scales::hue_pal()(nrow(data))
  } else if (length(fill) != nrow(data)) {
    .fill <- rep(fill, time = (ceiling(nrow(data)/length(fill))))[1:nrow(data)]
  }
    
  p <- ggplot(data, aes(x = 3, y = perc)) +
    geom_col(width = 1 , color = color, fill = .fill, lwd = .4, show.legend = F) +
    coord_polar('y', start = start * pi/180) +
    theme_void() +
    theme(aspect.ratio = 1, 
          plot.title = element_text(color = '#000000', size = 8 + font_size, hjust = 0.5))
  
  if (!is.null(title)) 
    p <- p + labs(title = as.character(title))
  
  if (isTRUE(lab_cir) & isFALSE(lab_cir_flip)) 
    p <- p + geom_text(aes(x = 3.5, y = ypos, label = lab, angle = -angle), 
                       size = font_size, hjust = 0.5)
  
  if (isFALSE(lab_cir) & isTRUE(lab_cir_flip)) 
    p <- p + geom_text(aes(x = 3.5, y = ypos, label = lab, angle = 270 - angle), 
                       size = font_size, hjust = 0.5)
  
  if (isFALSE(lab_cir) & isFALSE(lab_cir_flip)) 
    p <- p + geom_text(aes(x = 3.5, y = ypos, label = lab), size = font_size, hjust = 0.5)
  
  if (isTRUE(hemi)) 
    p <- p + lims(x = c(0, 3.5))
  
  message('  ggsave(file = \'pie.pdf\', width = 4, height = 4)')
  
  return(p)
}