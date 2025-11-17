# Pipeline for CF identification and profiling
# Jinxin Meng, mengjx855@163.com, 2025-11-09

#### CF identification ####
# https://tavazoielab.c2b2.columbia.edu/GHA/
diamond makedb --in top02prc_withCF.faa -d top02prc_withCF
diamond blastp -d top02prc_withCF.dmnd -q /data/database/uhgg/protein_catalogue/uhgp.faa -o uhgp.m8 -f 6 -p 100 --evalue 1e-5 --max-target-seqs 5 --sensitive --query-cover 60
seqkit fx2tab -n -i -l top02prc_withCF.faa > top02prc_withCF.len
# filter, scov > 70, qcov > 70, pid > 50, e-value < 1e-10
perl -e '%slen; open I, "top02prc_withCF.len"; while(<I>){chomp; @s=split/\t/; $slen{$s[0]}=$s[1]}; %qlen; open I, "/data/database/uhgg/protein_catalogue/uhgp.len"; while(<I>){chomp; @s=split/\t/; $qlen{$s[0]}=$s[1]}; while(<>){chomp; @s=split/\t/; if($s[0] ne $a){if($s[7]>$s[6]){$len=$s[7]-$s[6]+1}else{$len=$s[6]-$s[7]+1}; $qcov=$len/$qlen{$s[0]}*100; if($s[9]>$s[8]){$len=$s[9]-$s[8]+1}else{$len=$s[8]-$s[9]+1}; $scov=$len/$slen{$s[1]}*100; printf("$_\t%.2f\t%.2f\n", $qcov, $scov) if($qcov>70 && $scov>70 && $s[2]>50 && $s[10]<1e-10); $a=$s[0]}}' uhgp.m8 > uhgp.m8.f
cut -f1 uhgp.m8.f | seqkit grep -f - /data/database/uhgg/protein_catalogue/uhgp.ffn -o uhgp.fa
cat uhgp.m8.f | perl -e 'open I, "/data/database/uhgg/genomes-all_metadata.tsv"; %h; while(<I>){chomp;@s=split/\t/;$h{$s[0]}=$_}; while(<>){chomp;@s=split/\t/;$s[0]=~/(\S+?)_/;$x=$1; if(exists $h{$x}){print "$s[0]\t$s[1]\t$s[2]\t$s[12]\t$h{$x}\n"}else{$n=$x.".1";print "$s[0]\t$s[1]\t$s[2]\t$s[12]\t$h{$n}\n"}}' > uhgp.m8.info
emapper.py --cpu 0 -i top02prc_withCF.faa -o eggnog-mapper --evalue 1e-10


#### CF metagenomic profiling ####
# public dataset: BushmanFD_2020, PRJNA562600; FranzosaEA_2018, PRJNA400072; HallAB_2017, PRJNA385949; HeQ_2017, PRJEB15371; KumbhariA_2024, PRJNA993675; LloydPriceJ_2019, PRJNA398089; SchirmerM_2018, PRJNA389280; SchirmerM_2024, PRJNA436359; WengY_2019, PRJNA429990; YanQ_2023c, PRJEB67456
bowtie2-build uhgp.fa uhgp
for i in BushmanFD_2020 FranzosaEA_2018 HallAB_2017 HeQ_2017 KumbhariA_2024 LloydPriceJ_2019 SchirmerM_2018 SchirmerM_2024 WengY_2019 WengY_2019; do
    cat $i.clean_fq.filelist | parallel -j 8 --colsep="\t" run_bowtie2.sh {2} uhgp $i/{1} ;
    cut -f1 $i.clean_fq.filelist | parallel -j 8 samtools coverage -d 0 $i/{}.sort.bam -o $i/{}.cvg ;
    combine_file_zy_folder_allsample.py -D $i -suffix .cvg -t 1 -n 1 -v 6 -o $i.cvg ;
    combine_file_zy_folder_allsample.py -D $i -suffix .cvg -t 1 -n 1 -v 4 -o $i.rc ;
done

#### UHGG microbiota profiling ####
cat BushmanFD_2020.clean_fq.filelist | parallel -j 5 --colsep="\t" run_kraken2.sh {2} /data/database/uhgg/kraken2_db_uhgg_v2.0.2 BushmanFD_2020.kk2/{1} 100
cat FranzosaEA_2018.clean_fq.filelist | parallel -j 5 --colsep="\t" run_kraken2.sh {2} /data/database/uhgg/kraken2_db_uhgg_v2.0.2 FranzosaEA_2018.kk2/{1} 100
cat HallAB_2017.clean_fq.filelist | parallel -j 5 --colsep="\t" run_kraken2.sh {2} /data/database/uhgg/kraken2_db_uhgg_v2.0.2 HallAB_2017.kk2/{1} 100
cat HeQ_2017.clean_fq.filelist | parallel -j 5 --colsep="\t" run_kraken2.sh {2} /data/database/uhgg/kraken2_db_uhgg_v2.0.2 HeQ_2017.kk2/{1} 100
cat KumbhariA_2024.clean_fq.filelist | parallel -j 5 --colsep="\t" run_kraken2.sh {2} /data/database/uhgg/kraken2_db_uhgg_v2.0.2 KumbhariA_2024.kk2/{1} 100
cat LloydPriceJ_2019.clean_fq.filelist | parallel -j 5 --colsep="\t" run_kraken2.sh {2} /data/database/uhgg/kraken2_db_uhgg_v2.0.2 LloydPriceJ_2019.kk2/{1} 100
cat SchirmerM_2018.clean_fq.filelist | parallel -j 5 --colsep="\t" run_kraken2.sh {2} /data/database/uhgg/kraken2_db_uhgg_v2.0.2 SchirmerM_2018.kk2/{1} 100
cat SchirmerM_2024.clean_fq.filelist | parallel -j 5 --colsep="\t" run_kraken2.sh {2} /data/database/uhgg/kraken2_db_uhgg_v2.0.2 SchirmerM_2024.kk2/{1} 100
cat WengY_2019.clean_fq.filelist | parallel -j 5 --colsep="\t" run_kraken2.sh {2} /data/database/uhgg/kraken2_db_uhgg_v2.0.2 WengY_2019.kk2/{1} 150
cat YanQ_2023c.clean_fq.filelist | parallel -j 5 --colsep="\t" run_kraken2.sh {2} /data/database/uhgg/kraken2_db_uhgg_v2.0.2 YanQ_2023c.kk2/{1} 150

for i in BushmanFD_2020 FranzosaEA_2018 HallAB_2017 HeQ_2017 KumbhariA_2024 LloydPriceJ_2019 SchirmerM_2018 SchirmerM_2024 WengY_2019 WengY_2019; do
    combine_file_zy_folder_allsample.py -D $i.kk2/ -suffix .bk_out.s -n 1 -v 7 -t 1 -s 0 -o BushmanFD_2020.kk2.s ;
    combine_file_zy_folder_allsample.py -D $i.kk2/ -suffix .bk_out.p -n 1 -v 7 -t 1 -o BushmanFD_2020.kk2.p ;
    combine_file_zy_folder_allsample.py -D $i.kk2/ -suffix .bk_out.f -n 1 -v 7 -t 1 -o BushmanFD_2020.kk2.f ;
    combine_file_zy_folder_allsample.py -D $i.kk2/ -suffix .bk_out.g -n 1 -v 7 -t 1 -o BushmanFD_2020.kk2.g ;
done

#### MTX expression ####
for i in LloydPriceJ_2019 SchirmerM_2018; do
    cat $i.clean_fq.filelist | parallel -j 8 --colsep="\t" ../run_bowtie2.sh {2} ../uhgp $i/{1}
    cut -f1 $i.clean_fq.filelist | parallel -j 8 samtools coverage -d 0 $i/{}.sort.bam -o $i/{}.cvg
    combine_file_zy_folder_allsample.py -D $i -suffix .cvg -t 1 -n 1 -v 6 -o $i.cvg
    combine_file_zy_folder_allsample.py -D $i -suffix .cvg -t 1 -n 1 -v 4 -o $i.rc
    grep overall $i/*log | perl -ne '/\/(\S+?).log:(\S+%)/;print "$1\t$2\n"' > $i.map_rate
done

