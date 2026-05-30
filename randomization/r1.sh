# UHGP-90-wide abundance mapping
bowtie2-build --large-index --quiet --threads 80 /data/database/uhgg/protein_catalogue/uhgp-90.ffn uhgp-90

## metagenome
for i in BushmanFD_2020 FranzosaEA_2018 HallAB_2017 HeQ_2017 KumbhariA_2024 LloydPriceJ_2019 SchirmerM_2018 SchirmerM_2024 WengY_2019 WengY_2019; do
    mkdir -p uhgp_cvg/$i ; 
    cat uhgp_cvg/$i.clean_fq.filelist | parallel -j 8 --colsep="\t" scripts/run_bowtie2.sh {2} uhgp-90 uhgp_cvg/$i/{1} ;
done

## RNASeq
for i in LloydPriceJ_2019 SchirmerM_2018; do
    mkdir -p uhgp_cvg_RNASeq/$i ; 
a    cat uhgp_cvg_RNASeq/$i.clean_fq.filelist | parallel -j 8 --colsep="\t" scripts/run_bowtie2.sh {2} uhgp-90 uhgp_cvg/$i/{1} ;
done



# POG randomization framework
## UHGP-90 genes detected in metagenomic samples with a sample prevalence > 1% were retained
realpath uhgp_cvg/*/*cvg > uhgp_cvg.filepath
cut -f2 uhgp_cvg.filepath | parallel -j 10 "awk '\$6>50' {} > {}.f"
ls *filelist | cut -d. -f1 | parallel -j 10 "cut -f1 {}/*.cvg.f | grep -v 'rname' | LC_ALL=C sort | uniq -c | awk '{print \$2, \$1}' > pogs/prevalence/{}.gene_count"
awk 'BEGIN{FS=" "; OFS="\t"} {if($1 in c){c[$1]+=$2}else{c[$1]=$2}} END{for(i in c){print i, c[i]}}' *count > total.gene_count
awk '$2>(2858*0.01)' total.gene_count > total.gene_count.f
cut -f1 total.gene_count.f | seqkit grep -f - /data/database/uhgg/protein_catalogue/uhgp-90.faa > uhgp-90-f.faa

## all-vs-all align
diamond makedb --in uhgp-90-f.faa --db uhgp-90-f.dmnd --threads 56 --quiet 
diamond blastp -q uhgp-90-f.faa --db uhgp-90-f.dmnd -o uhgp-90-f.m8 --outfmt 6 qseqid sseqid evalue --max-target-seqs 560000 -e 1e-5 --threads 111 --quiet

## mcl cluster
mcxload -abc uhgp-90-f.m8 --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o uhgp-90-f.mci -write-tab uhgp-90-f.tab 2>/dev/null
mcl uhgp-90-f.mci -I 2.5 -te 56 -use-tab uhgp-90-f.tab -o uhgp-90-f.mcl 2>/dev/null
awk 'BEGIN{FS=OFS="\t"} NF>0{fam=sprintf("FAM%06d", NR); for(i=1;i<=NF;i++){print $i, fam}}' uhgp-90-f.mcl > uhgp-90-f.clu

## 若一个 MCL family 至少包含 3 个基因，且其中 CF-hit genes 占该 family 总基因数的比例 ≥10%，则定义为 confident CF-associated MCL family
awk 'BEGIN{FS=OFS="\t"} NR==FNR{split($2,s,","); gene[$1]=s[1];next} {if($1 in gene){print $0, gene[$1]}else{print $0, "others"}}' ../../pipeline/uhgp.m8.f.drop uhgp-90-f.clu > uhgp-90-f.clu.lab
grep CF uhgp-90-f.clu.lab | cut -f2 | uniq | grep -f - uhgp-90-f.clu.lab | cut -f2,3 | sort | uniq -c | awk '{print $2"\t"$3"\t"$1}' > uhgp-90-f.clu.lab.stats

## manually compete CF-like POG confidence, as POG_CF_confident.tsv
awk 'BEGIN{FS=OFS="\t"} NR==FNR{h[$1]=1;next} {if($2 in h){sub(/FAM/,"CFF",$2);print $1,$2}else{print $0}}' POG_CF_confident.tsv uhgp-90-f.clu > uhgp-90-f.clu.final
awk 'BEGIN{FS=OFS="\t"} {c[$2]++} END{for(i in c){print i, c[i]}}' uhgp-90-f.clu.final | sort -rnk2 > uhgp-90-f.clu.final.count
awk '$2>2' uhgp-90-f.clu.final.count | cut -f1 | grep -f - uhgp-90-f.clu.final > uhgp-90-f.clu.final.f   ## final POG used to subsequent analyses without singleton

## tpm metagenome profile
cat uhgp_cvg.filepath | parallel -j 30 scripts/run-cvg2tpm-fam.awk pogs/uhgp-90-f.clu.final.f {} ">" merge_tpm/{/.}.fam.tpm
cat uhgp_cvg.filepath | parallel -j 30 scripts/run-cvg2tpm-gene.awk pogs/uhgp-90-f.clu.final.f {} ">" merge_tpm/{/.}.gene.tpm
scripts/run-fast-combine-matrix.py -d merge_tpm/ -o merge_tpm.fam.tpm -n 1 -v 2 -t 1 -suffix .fam.tpm
scripts/run-fast-combine-matrix.py -d merge_tpm/ -o merge_tpm.gene.tpm -n 1 -v 2 -t 1 -suffix .gene.tpm

## tpm RNASeq profile
cat uhgp_cvg_RNASeq.filepath | parallel -j 30 scripts/run-cvg2tpm-fam.awk pog/uhgp-90-f.clu.final.f {} ">" merge_tpm_RNASeq/{/.}.fam.tpm
cat uhgp_cvg_RNASeq.filepath | parallel -j 30 scripts/run-cvg2tpm-gene-RNASeq.awk pog/uhgp-90-f.clu.final.f {} ">" merge_tpm_RNASeq/{/.}.gene.tpm
scripts/run-fast-combine-matrix.py -d merge_tpm_RNASeq/ -o merge_tpm_RNASeq.fam.tpm -n 1 -v 2 -t 1 -suffix .fam.tpm
scripts/run-fast-combine-matrix.py -d merge_tpm_RNASeq/ -o merge_tpm_RNASeq.gene.tpm -n 1 -v 2 -t 1 -suffix .gene.tpm



# random metagenome profile
## mjx-random-select.R 生成随机基因家族 with 999
## abundance profile
scripts/run-split-random-tpm-fam.pl results/1.fam.random.rep999.tsv merge_tpm.fam.tpm random_profile .fam.tpm
scripts/run-split-random-tpm-gene.pl results/1.fam.random.rep999.tsv pogs/uhgp-90-f.clu.final.f merge_tpm.gene.tpm random_profile
## presence/absence profile
scripts/run-split-mapping-by-rep.pl results/1.fam.random.rep999.tsv pogs/uhgp-90-f.clu.final.f random_profile
cat results/1.fam.random.rep999.ids | parallel -j 56 scripts/run-fam-to-species-binary.py random_profile/{}.mapping \
    /data/database/uhgg/protein_catalogue/uhgp-90.tsv /data/database/uhgg/genomes-all_metadata.tsv random_profile/{}.binary

# random imprint
## calculate distance
cat results/1.fam.random.rep999.ids | parallel -j 56 calcu-distance.py random_profile/{}.binary random_profile/{}.dist --method jaccard
## random adonis
cat results/1.fam.random.rep999.ids | parallel -j 56 scripts/run-imprint-adonis2.R random_profile/{}.dist random_imprint/{}
## mantel
cat results/1.fam.random.rep999.ids | parallel -j 56 scripts/run-imprint-mantel.R random_profile/{}.dist random_imprint/{}

# random diversity
cat results/1.fam.random.rep999.ids | parallel -j 56 scripts/run-diversity.R random_profile/{}.gene.tpm random_div/{}.gene
cat results/1.fam.random.rep999.ids | parallel -j 56 scripts/run-diversity.R random_profile/{}.fam.tpm random_div/{}.fam

# random procrustes 
cat results/1.fam.random.rep999.ids | parallel -j 56 scripts/run-procrustes.R random_profile/{}.fam.tpm random_procrustes/{}.fam

# random adonis explain
cat results/1.fam.random.rep999.ids | parallel -j 56 scripts/run-adonis-by-fam.R random_profile/{}.fam.tpm random_adonis/{}.fam
cat results/1.fam.random.rep999.ids | parallel -j 56 scripts/run-adonis-by-species.R random_profile/{}.fam.tpm random_adonis/{}.fam

# random rf
## rf
cat results/1.fam.random.rep999.ids | parallel -j 56 scripts/run-randomforest.R random_profile/{}.fam.tpm random_rf/{}.fam
## feature_importance run in results/ folder
parallel -j 40 ../scripts/run-randomforest-topN-feature.R {} 6.feature_importance/top_{}.auc.tsv ::: `seq -w 1 1 73`
## 10.fam.random.rep999.15.tsv generated by mjx-random-select.R
scripts/run-split-random-tpm-fam.pl results/10.fam.random.rep999.15.tsv merge_tpm.fam.tpm random_profile_top15 .fam.tpm
cat results/10.fam.random.rep999.15.ids | parallel -j 100 scripts/run-randomforest.R random_profile_top15/{}.fam.tpm random_rf_top15/{}.fam



# random metatransciptome profile
## mjx-random-select-RNASeq.R 生成随机基因家族 999
scripts/run-split-random-tpm-fam.pl results_RNASeq/1.fam.random.rep999.tsv merge_tpm_RNASeq.fam.tpm random_profile_RNASeq .fam.tpm
scripts/run-split-random-tpm-fam.pl results_RNASeq/1.fam.random.rep999.tsv merge_tpm.fam.tpm random_profile_RNASeq .fam.tpm.DNA
scripts/run-split-random-tpm-gene.pl results_RNASeq/1.fam.random.rep999.tsv pogs/uhgp-90-f.clu.final.f merge_tpm_RNASeq.gene.tpm random_profile_RNASeq

# RNA vs DNA abundance
cat results_RNASeq/1.fam.random.rep999.ids | parallel -j 110 scripts/run-DNA-RNA-correlation.R random_profile_RNASeq/{}.fam.tpm.DNA random_profile_RNASeq/{}.fam.tpm random_correlation_RNASeq/{}

# radnom diversity
cat results_RNASeq/1.fam.random.rep999.ids | parallel -j 56 scripts/run-diversity-RNASeq.R random_profile_RNASeq/{}.gene.tpm random_div_RNASeq/{}.gene

# random rf 
## rf
cat results_RNASeq/1.fam.random.rep999.ids | parallel -j 56 scripts/run-randomforest-RNASeq.R random_profile_RNASeq/{}.fam.tpm random_rf_RNASeq/{}.fam
## feature_importance
parallel -j 73 scripts/run-randomforest-topN-feature.R {} results_RNASeq/4.feature_importance/top_{}.auc.tsv ::: `seq -w 1 1 73`
## 5.fam.random.rep999.18.tsv generated by mjx-random-select-RNASeq.R
scripts/run-split-random-tpm-fam.pl results_RNASeq/5.fam.random.rep999.18.tsv merge_tpm_RNASeq.fam.tpm random_profile_RNASeq_top18 .fam.tpm
cat results_RNASeq/5.fam.random.rep999.18.ids | parallel -j 110 scripts/run-randomforest-RNASeq.R random_profile_RNASeq_top18/{}.fam.tpm random_rf_RNASeq_top18/{}.fam



# random taxa for randomforest
## mjx-randomforest-taxa.R 生成随机物种 999 79/46/13
for i in BushmanFD_2020 FranzosaEA_2018 HallAB_2017 HeQ_2017 KumbhariA_2024 LloydPriceJ_2019 SchirmerM_2018 SchirmerM_2024 WengY_2019 YanQ_2023c;do
    parallel -j 100 --colsep="\t" \
        "scripts/run-random-taxa-rf.R {2} \
        ../data/$i.kk2.s.gz \
        ../data/$i.sample_group.gz \
        random_species/{1}" :::: select_taxa.species.79.tsv
    cat random_species/rep* > random_species/$i.roc.species.79.results
    /usr/bin/rm random_species/rep*
done
## mjx-randomforest-taxa.R 生成随机属 999 79/46/13
for i in BushmanFD_2020 FranzosaEA_2018 HallAB_2017 HeQ_2017 KumbhariA_2024 LloydPriceJ_2019 SchirmerM_2018 SchirmerM_2024 WengY_2019 YanQ_2023c;do
    parallel -j 100 --colsep="\t" \
        "scripts/run-random-taxa-rf.R {2} \
        ../data/$i.kk2.g.gz \
        ../data/$i.sample_group.gz \
        random_genus/{1}" :::: select_taxa.genus.79.tsv
    cat random_genus/rep* > random_genus/$i.roc.genus.79.results
    /usr/bin/rm random_genus/rep*
done

