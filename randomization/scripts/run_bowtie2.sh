#!/usr/bin/bash

( [ $# -lt 3 ] ) &&  { echo -e "Usage: $0 [fq|fq1,fq2] [*.bt2] [out_prefix] [threads:16]" && exit 2; }
( [ -e $3.log ] && grep -q 'overall' $3.log ) && echo -e "[$(date "+%F %T")] Skip sample: ${3##*/} .." && exit 0

trds=${4:-20}

if [[ $1 =~ "," ]];then
    fq1=$(echo $1 | cut -d "," -f1)
    fq2=$(echo $1 | cut -d "," -f2)
    bowtie2 --end-to-end --mm --sensitive --no-unal -1 $fq1 -2 $fq2 -x $2 -p $trds 2>> $3.log |\
        samtools view -@ 10 -bS 2>/dev/null |\
        samtools sort -m 8G -@ 20 - 2>/dev/null |\
        samtools coverage -d 0 - | awk '$4>0' > $3.cvg
else
    bowtie2 --end-to-end --mm --sensitive --no-unal -U $1 -x $2 -p $trds 2>> $3.log |\
        samtools view -@ 10 -bS 2>/dev/null |\
        samtools sort -m 8G -@ 20 - 2>/dev/null |\
        samtools coverage -d 0 - | awk '$4>0' > $3.cvg
fi
