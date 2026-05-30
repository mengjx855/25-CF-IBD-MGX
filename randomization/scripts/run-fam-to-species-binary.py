#!/data/software/miniconda3/envs/jinxin/bin/python3
# author: Jin-Xin Meng <jinxmeng@zju.edu.cn>
# created date: 2026-04-06 11:42:51
# modified date: 2026-04-06 11:42:51

from collections import defaultdict
import pandas as pd
import re, sys

if len(sys.argv) != 5:
    sys.exit(
        f"Usage: python {sys.argv[0]} [mapping.tsv]\n"
        f"              [uhgp-90-cluster: /data/database/uhgg/protein_catalogue/uhgp-90.tsv]\n"
        f"              [uhgg-genomes-metadata: /data/database/uhgg/genomes-all_metadata.tsv]\n"
        f"              [out_file]\n"
    )

def main():
    gene_fam_file = sys.argv[1]
    cluster_file = sys.argv[2]
    metadata_file = sys.argv[3]
    out_file = sys.argv[4]
    
    # 1. 读取 gene -> fam 映射
    gene2fam = defaultdict(set)
    fam_list = set()

    with open(gene_fam_file, "r") as f:
        for row in f:
            row = row.strip()
            if not row:
                continue
            gene, fam = row.split("\t")[:2]
            gene2fam[gene].add(fam)
            fam_list.add(fam)

    # 2. 读取 UHGP-90 聚类
    fam2genome = defaultdict(set)
    with open(cluster_file, "r") as f:
        for row in f:
            row = row.strip()
            if not row:
                continue
            rep, member = row.split("\t")[:2]
            if rep not in gene2fam:
                continue
            genome = member.rpartition("_")[0]
            for fam in gene2fam[rep]:
                fam2genome[fam].add(genome)

    # 3. 读取基因组元数据
    genome2species = {}
    count_genome_in_species = defaultdict(int)

    with open(metadata_file, "r") as f:
        for row in f:
            if row.startswith("Genome"):
                continue
            row = row.rstrip("\n")
            if not row:
                continue
            splits = row.split("\t")
            genome = splits[0]
            species = splits[13]
            genome = re.sub(r"\.\d+$", "", genome)
            species = re.sub(r"\.\d+$", "", species)
            genome2species[genome] = species
            count_genome_in_species[species] += 1

    # 4. 计算每个 fam 在每个 species 中的 genome 覆盖率
    #   如果某 fam 出现在某 species >=50% 的 genome 中，则记为 1
    result = defaultdict(dict)
    for fam in fam_list:
        if fam not in fam2genome:
            continue
        count_present = defaultdict(int)
        for genome in fam2genome[fam]:
            if genome not in genome2species:
                continue
            species = genome2species[genome]
            count_present[species] += 1
        for species, present_n in count_present.items():
            total_n = count_genome_in_species[species]
            ratio = present_n / total_n
            if ratio >= 0.5:
                result[fam][species] = 1
    
    # 5. 输出矩阵
    data = pd.DataFrame(result).fillna(0).astype(int)
    data.to_csv(out_file, sep='\t', index_label='name')

if __name__ == '__main__':
    main()

