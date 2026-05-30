#!/data/software/miniconda3/envs/r-base-4.4/bin/Rscript
# author: Jin-Xin Meng <jinxmeng@zju.edu.cn>
# created date: 2026-04-06 17:27:55
# modified date: 2026-04-06 17:27:55

args <- commandArgs(trailingOnly = T)

if (length(args) < 2) {
    message("Usage: Rscript run-imprint-adonis2.R [distance] [out_prefix]")
    # status = 1 表示非正常退出（给 Linux 系统看的信号）
    # save = "no" 表示不保存工作空间
    q(save = "no", status = 1)
}

calcu_adjusted_r2 <- function(adonis) {
  n_observations <- adonis$Df[3]+1
  d_freedom <- adonis$Df[1]
  r2 <- adonis$R2[1]
  adjusted_r2 <- vegan::RsquareAdj(r2, n_observations, d_freedom)
  return(adjusted_r2)
}

########### read distance ###########
distance <- data.table::fread(args[1]) |> 
  tibble::column_to_rownames('name') |> 
  as.dist()

########### read metadata ###########
genome_info <- data.table::fread('/data/database/uhgg/genomes-all_metadata.tsv') |> 
  dplyr::mutate(
    genome = sub('\\.\\d+$', '', Genome),
    domain = stringr::str_split_i(Lineage, pattern = ';', 1),
    phylum = stringr::str_split_i(Lineage, pattern = ';', 2),
    class  = stringr::str_split_i(Lineage, pattern = ';', 3),
    order  = stringr::str_split_i(Lineage, pattern = ';', 4),
    family = stringr::str_split_i(Lineage, pattern = ';', 5),
    genus  = stringr::str_split_i(Lineage, pattern = ';', 6)
  )

metadata <- tibble::tibble(genome = labels(distance)) |>
  dplyr::left_join(
    dplyr::select(genome_info, genome, domain, phylum, 
                  class, order, family, genus) |> 
      dplyr::distinct(), 
    by = 'genome')

taxa_level <- c('domain','phylum','class','order','family','genus')

if (length(unique(metadata[['domain']])) < 2) taxa_level <- taxa_level[-1]

########### adonis test ###########
test <- purrr::map(
  taxa_level, ~
    vegan::adonis2(
      as.formula(paste0('distance ~ ', .x)), metadata,
      permutations = 1, parallel = 1) ) |>
  purrr::set_names(taxa_level)

purrr::map2_vec(
  names(test), test, ~ 
    tibble::tibble(
      name = stringr::str_to_title(.x), 
      r2 = .y$R2[1],
      r2adj = calcu_adjusted_r2(.y), 
      pval = .y[1, 5]) 
  ) |> 
  write.table(paste0(args[2], '.pcoa.tsv'), sep = '\t', quote = F, row.names = F)
