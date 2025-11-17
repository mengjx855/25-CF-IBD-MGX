#!/usr/bin/bash

( [ $# -lt 3 ] ) &&  { echo -e "Usage: $0 [fq|fq1,fq2] [*.bt2] [out_prefix] [threads:16]" && exit 2; }
( [ -e $3.sort.bam ] ) && { echo -e "Skip sample: ${3##*/} .." && exit 0; }

trds=${4:-24}

if [[ $1 =~ "," ]];then
    fq1=$(echo $1 | cut -d "," -f1)
    fq2=$(echo $1 | cut -d "," -f2)
    bowtie2 --end-to-end --sensitive --no-unal -1 $fq1 -2 $fq2 -x $2 -p $trds 2>> $3.log |\
        samtools view -@ $trds -bS 2>/dev/null | samtools sort -@ $trds -o $3.sort.bam - 2>/dev/null
else
    bowtie2 --end-to-end --sensitive --no-unal -U $1 -x $2 -p $trds 2>> $3.log |\
        samtools view -@ $trds -bS 2>/dev/null | samtools sort -@ $trds -o $3.sort.bam - 2>/dev/null
fi
