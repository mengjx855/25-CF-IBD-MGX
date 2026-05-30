#### Jinxin Meng, 20260404, 20260530 ####
setwd('/data/mengjx/project/10.20250623_IBD_BAC_CF_Landscape/github/randomization/')
pacman::p_load(tidyverse, ggpubr, rstatix, openxlsx)

#### matched families random selection 73 ####
profile <- data.table::fread('merge_tpm_RNASeq.fam.tpm') |> 
  column_to_rownames('name')

fam_size <- read.delim(
  'pogs/uhgp-90-f.clu.final.count', 
  header = F, col.names = c('fam', 'n'))

data <- data.frame(
  name = rownames(profile),
  abun = rowMeans(log10(profile + 1)),
  prev = rowSums(profile > 1e-6) / ncol(profile),
  row.names = NULL) |> 
  mutate(
    size = log10(fam_size$n[match(name, fam_size$fam)] + 1),
    size_cut = dplyr::ntile(size, 10),
    abun_cut = dplyr::ntile(abun, 10),
    prev_cut = dplyr::ntile(prev, 10),
    bin3d = paste0('s', size_cut, 'a', abun_cut, 'p', prev_cut) )

match_data <- left_join(
  filter(data, grepl('CFF', name)) |> 
    count(bin3d, name = 'num_CF'),
  filter(data, !grepl('CFF', name)) |> 
    count(bin3d, name = 'num_FAM'), 
  by = c('bin3d'))

data_nonCF <- filter(data, !grepl('CFF', name))

pools <- split(data_nonCF$name, data_nonCF$bin3d)

seed_seqs <- seq(2026, 3024, 1)

seed_names <- paste0(
  'rep', sprintf('%03d', seq_along(seed_seqs)), 
  '_seed', seed_seqs)

random_select <- map(seed_seqs, \(x) {
  set.seed(x)
  map(
    seq_len(nrow(match_data)), ~ {
      row <- match_data[.x, ]
      sample(pools[[row$bin3d]], size = row$num_CF, replace = F) 
      } )|> 
    unlist(use.names = F) 
  } 
) |> 
  set_names(seed_names)

random_select$rep000_seed0000 <- filter(data, grepl('CF', name)) |>
  pull(name)

result <- tibble::enframe(
  random_select, name = 'name', value = 'families'
) |>
  tidyr::unnest_wider(families, names_sep = '_') |> 
  arrange(name)

write_tsv(
  result, 'results_RNASeq/1.fam.random.rep999.tsv', 
  col_names = F, quote = NULL
)

write_tsv(
  select(result, name), 'results_RNASeq/1.fam.random.rep999.ids', 
  col_names = F, quote = NULL
)

#### matched families random selection 18 ####
selected <- read_rds('results_RNASeq/4.feature_importance.names.rds')[1:18]

profile <- data.table::fread('merge_tpm_RNASeq.fam.tpm') |> 
  column_to_rownames('name')

profile <- profile[
  (grepl('FAM', rownames(profile)) | rownames(profile) %in% selected)
  , ]

fam_size <- read.delim(
  'pogs/uhgp-90-f.clu.final.count', 
  header = F, col.names = c('fam', 'n'))

data <- data.frame(
  name = rownames(profile),
  abun = rowMeans(log10(profile + 1)),
  prev = rowSums(profile > 1e-6) / ncol(profile),
  row.names = NULL) |> 
  mutate(
    size = log10(fam_size$n[match(name, fam_size$fam)] + 1),
    size_cut = dplyr::ntile(size, 10),
    abun_cut = dplyr::ntile(abun, 10),
    prev_cut = dplyr::ntile(prev, 10),
    bin3d = paste0('s', size_cut, 'a', abun_cut, 'p', prev_cut)
  )

match_data <- left_join(
  filter(data, grepl('CFF', name)) |> 
    count(bin3d, name = 'num_CF'),
  filter(data, !grepl('CFF', name)) |> 
    count(bin3d, name = 'num_FAM'), 
  by = c('bin3d'))

data_nonCF <- filter(data, !grepl('CFF', name))

pools <- split(data_nonCF$name, data_nonCF$bin3d)

seed_seqs <- seq(2026, 3024, 1)

seed_names <- paste0(
  'rep', sprintf('%03d', seq_along(seed_seqs)), 
  '_seed', seed_seqs)

random_select <- map(
  seed_seqs, \(x) {
    set.seed(x)
    map(
      seq_len(nrow(match_data)), ~ {
        row <- match_data[.x, ]
        sample(pools[[row$bin3d]], size = row$num_CF, replace = F) 
      } 
    ) |> 
      unlist(use.names = F) 
  } 
) |> 
  set_names(seed_names)

random_select$rep000_seed0000 <- filter(data, grepl('CFF', name)) |> pull(name)

result <- tibble::enframe(
  random_select, name = 'name', value = 'families'
) |>
  tidyr::unnest_wider(families, names_sep = '_') |> 
  arrange(name)

write_tsv(
  result, 
  'results_RNASeq/5.fam.random.rep999.18.tsv', 
  col_names = F, quote = NULL
)

write_tsv(
  select(result, name), 
  'results_RNASeq/5.fam.random.rep999.18.ids', 
  col_names = F, quote = NULL
)

