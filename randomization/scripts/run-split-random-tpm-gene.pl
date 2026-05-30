#!/usr/bin/env perl
# author: Jin-Xin Meng <jinxmeng@zju.edu.cn>

use strict;
use warnings;

# 用法检查
die "Usage: $0 <rep_family.tsv> <mapping.tsv> <gene_tpm.tsv> <outdir> [suffix:.gene.tpm]\n" unless @ARGV >= 4;

my ($repfam, $mapping, $tpm, $outdir, $suffix) = @ARGV;

$suffix //= '.gene.tpm';

# 1. family -> rep
my %fam2rep;
open my $F, "<", $repfam or die $!;
while (<$F>) {
    chomp;
    next if /^\s*$/;
    my @s = split /\t/;
    my $rep = shift @s;
    for my $fam (@s) {
        next unless $fam;
        push @{$fam2rep{$fam}}, $rep;
    }
}
close $F;


# 2. gene -> rep
my %gene2rep;
open my $M, "<", $mapping or die $!;
while (<$M>) {
    chomp;
    next if /^\s*$/;
    my ($gene, $fam) = split /\t/;
    next unless exists $fam2rep{$fam};
    for my $rep (@{$fam2rep{$fam}}) {
        $gene2rep{$gene}{$rep} = 1;
    }
}
close $M;


# 3. 拆分 TPM
my %fh;        # 文件句柄缓存
my %printed;   # 是否写过表头
open my $T, "<", $tpm or die $!;
my $head = <$T>;
while (<$T>) {
    chomp;
    next if /^\s*$/;
    my @s = split /\t/;
    my $gene = $s[0];
    next unless exists $gene2rep{$gene};

    for my $rep (keys %{$gene2rep{$gene}}) {
        # 惰性打开文件（只打开一次）
        if (!exists $fh{$rep}) {
            open my $O, ">", "$outdir/$rep$suffix" or die $!;
            $fh{$rep} = $O;
        }

        my $O = $fh{$rep};
        if (!exists $printed{$rep}) {
            print $O $head;
            $printed{$rep} = 1;
        }

        print $O "$_\n";
    }
}
close $T;


# 4. 关闭所有文件
for my $rep (keys %fh) {
    close $fh{$rep};
}

