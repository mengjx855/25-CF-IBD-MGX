#!/usr/bin/awk -f
# author: Jin-Xin Meng <jinxmeng@zju.edu.cn> 
# created date: 2026-04-01 10:30:00
# modified date: 2026-04-01 10:30:50

BEGIN {
    if (ARGC < 3) {
        print "Usage: run-cvg2tpm-fam.awk [mapping.tsv] [sample.cvg]";
        print "Note: convert cvg to tpm, and then aggregate theme by higher level classification using mapping file."
        exit 0;
    }
    FS = "\t";
    OFS = "\t";
    file_num = 0;
    ARGV[3] = ARGV[2];
    ARGC = 4;
}

FNR == 1 {
    file_num++;
}

file_num == 1 {
    map[$1] = $2;
    next;
}

# RPK
file_num == 2 {
    if (FNR > 1 && $4 > 0) {
        len = $3 - $2 + 1;
        sum += $4 / len;
    }
    next;
}

file_num == 3 {
    if (FNR == 1 || $4 == 0) next;

    if ($1 in map) {
        name = $1;
        len  = $3 - $2 + 1;
        tpm  = ($4 / len / sum) * 10^6;
        fam_id = map[name];
        fam_tpm[fam_id] += tpm;
    }
}

# 输出结果
END {
    if (ARGC < 4) exit 0;
    for (f in fam_tpm) {
        printf "%s\t%.6f\n", f, fam_tpm[f];
    }
}

