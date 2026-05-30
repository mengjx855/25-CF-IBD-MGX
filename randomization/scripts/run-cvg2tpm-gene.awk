#!/usr/bin/awk -f
# author: Jin-Xin Meng <jinxmeng@zju.edu.cn>, and chat-gpt-v5.4
# modified: 2026-04-15

BEGIN {
    if (ARGC < 2) {
        print "Usage: run-cvg2tpm-gene.awk [mapping.tsv] [sample.cvg]" > "/dev/stderr";
        print "  keep gene in pog-mapping files";
        exit 1;
    }
    FS = "\t";
    OFS = "\t";
    ARGV[3] = ARGV[2];
    ARGC = 4;
}

FNR == 1 {
    file_num++;
}

file_num == 1 {
    map[$1] = $2;
    next
    }

file_num == 2 {
    if (FNR == 1) next; 
    len = $3 - $2 + 1;
    if ($4 > 0 && $6 >= 50) {
        sum_rpk += $4 / len;
    }
    next;
}

file_num == 3 {
    if (sum_rpk <= 0) {
        print "Error: No genes passed the coverage threshold in " ARGV[2] > "/dev/stderr";
        exit 1;
    }

    if (FNR == 1) next
    
    if ($1 in map) {
        len = $3 - $2 + 1;
        if ($4 > 0 && $6 >= 50) {
            tpm = ($4 / len / sum_rpk) * 10^6;
            printf "%s\t%.6f\n", $1, tpm;
        }
    }
}

