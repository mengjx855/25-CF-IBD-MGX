#!/usr/bin/bash
shopt -s expand_aliases

usage(){
    echo "Usage: $0 [fq|fq1,fq2] [kk2_db] [out_prefix] [read_len:-150]" 
    echo "  available database:"
    echo "  1. /data/database/uhgg/kraken2_db_uhgg_v2.0.2"
    echo "  2. /data/database/kraken2/k2_standard_16_GB_20250714"
    echo "  3. /data/database/kraken2/k2_pluspf_16_GB_20250714"
    echo "  4. /data/database/kraken2/16S_Silva138_20200326"
    echo "  5. /data/database/kraken2/k2_viral_20250714"
}

( [ $# -lt 3 ]) && { usage && exit 2; }

fq=$1; db=$2; out=$3; len=${4:-150}; trds=8
alias kraken2=/data/software/miniconda3/envs/kraken2/bin/kraken2
alias bracken=/data/software/miniconda3/envs/kraken2/bin/bracken

( [ -f $out.log ] && grep -q 'Bracken complete' $out.log ) && \
    echo -e "[$(date +%Y-%m-%d\ %H:%M:%S)] Skip sample: ${out##*/}." && exit 0

if [[ $fq =~ "," ]];then
    fq1=$(echo $fq | cut -d "," -f1)
    fq2=$(echo $fq | cut -d "," -f2)
    kraken2 --db $db --output $out.kk2_out --report $out.kk2_report --threads $trds --use-names --gzip-compressed $fq1 $fq2 >> $out.log 2>&1 && rm $out.kk2_out
else
    kraken2 --db $db --output $out.kk2_out --report $out.kk2_report --threads $trds --use-names --gzip-compressed $fq >> $out.log 2>&1 && rm $out.kk2_out
fi

bracken -l P -r $len -d $db -i $out.kk2_report -o $out.bk_out.p -w $out.bk_report.p >> $out.log 2>&1 && rm $out.bk_report.p
bracken -l F -r $len -d $db -i $out.kk2_report -o $out.bk_out.f -w $out.bk_report.f >> $out.log 2>&1 && rm $out.bk_report.f
bracken -l G -r $len -d $db -i $out.kk2_report -o $out.bk_out.g -w $out.bk_report.g >> $out.log 2>&1 && rm $out.bk_report.g
bracken -l S -r $len -d $db -i $out.kk2_report -o $out.bk_out.s -w $out.bk_report.s >> $out.log 2>&1 && rm $out.bk_report.s
bracken -l S1 -r $len -d $db -i $out.kk2_report -o $out.bk_out.t -w $out.bk_report.t >> $out.log 2>&1 && rm $out.bk_report.t

