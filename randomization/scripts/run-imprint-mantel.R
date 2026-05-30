#!/data/software/miniconda3/envs/r-base-4.4/bin/Rscript
# author: Jin-Xin Meng <jinxmeng@zju.edu.cn>
# created date: 2026-04-06 17:27:55
# modified date: 2026-04-15 08:41:19

args <- commandArgs(trailingOnly = T)

if (length(args) < 2) {
  message("Usage: Rscript run-imprint-mantel.R [jaccard distance] [out_prefix]")
  message('Bacterial phylogenetic tree: /data/database/uhgg/phylogenies/bac120_iqtree.nwk')
  message('Archaeal phylogenetic tree: /data/database/uhgg/phylogenies/ar122_iqtree.nwk')
  # status = 1 表示非正常退出（给 Linux 系统看的信号）
  # save = "no" 表示不保存工作空间
  q(save = "no", status = 1)
}

pacman::p_load(tidyverse)

########### read distance ###########
distance <- data.table::fread(args[1]) %>% 
  column_to_rownames('name') %>% 
  as.dist()

########### mantel test ###########
genome_info <- data.table::fread('/data/database/uhgg/genomes-species_metadata.tsv') %>% 
  mutate(
    genome = sub('\\.\\d+$', '', Genome),
    domain = str_split_i(Lineage, pattern = ';', 1)
  )

taxa <- data.frame(ref_species = labels(distance)) %>% 
  mutate(domain = genome_info$domain[match(ref_species, genome_info$genome)])

# bacterial
tr <- ape::read.tree('/data/database/uhgg/phylogenies/bac120_iqtree.nwk')

select_taxa <- filter(taxa, domain == 'd__Bacteria') %>% pull(ref_species)
select_taxa <- select_taxa[select_taxa %in% tr$tip.label]

patristic_dist <- data.frame(ape::cophenetic.phylo(tr)[select_taxa, select_taxa])
jaccard_dist <- data.frame(as.matrix(distance)[select_taxa, select_taxa])

bacterial_test <- vegan::mantel(
  patristic_dist, jaccard_dist, method = 'spearman', 
  permutations = 1, parallel = 1
)

# archaeal
tr <- ape::read.tree('/data/database/uhgg/phylogenies/ar122_iqtree.nwk')

select_taxa <- filter(taxa, domain == 'd__Archaea') %>% pull(ref_species)
select_taxa <- select_taxa[select_taxa %in% tr$tip.label]

if (length(select_taxa) < 3) {
  message("Warning: Number of archaeal genomes is too low (n = ", length(select_taxa), ").")
  
  data.frame(
    domain = 'bacteria',
    statistic = bacterial_test$statistic,
    pval = bacterial_test$signif,
    permutations = 1,
    method = 'mantel test (spearman)'
  ) %>% 
    write.table(paste0(args[2], '.mantel.tsv'), sep = '\t', quote = F, row.names = F)
  
  q(save = "no", status = 0)
}

patristic_dist <- data.frame(ape::cophenetic.phylo(tr)[select_taxa, select_taxa])
jaccard_dist <- data.frame(as.matrix(distance)[select_taxa, select_taxa])

archaeal_test <- vegan::mantel(
  patristic_dist, jaccard_dist, method = 'spearman', 
  permutations = 1, parallel = 1
)

data.frame(
  domain = c('bacteria', 'archaea'),
  statistic = c(bacterial_test$statistic, archaeal_test$statistic),
  pval = c(bacterial_test$signif, archaeal_test$signif),
  permutations = 1,
  method = 'mantel test (spearman)'
  ) %>% 
  write.table(paste0(args[2], '.mantel.tsv'), sep = '\t', quote = F, row.names = F)
