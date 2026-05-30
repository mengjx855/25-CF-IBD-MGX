#!/usr/bin/env perl
# author: Jin-Xin Meng <jinxmeng@zju.edu.cn>

use strict;
use warnings;
use File::Path qw(make_path);

# 根据随机蛋白质家族表拆分mapping数据，
# 为下一步得到不同基因家族在物种中的0/1矩阵

die "Usage: $0 <rep_family.tsv> <mapping.tsv> <outdir>\n" unless @ARGV == 3;

my ($rep_family_file, $mapping_file, $outdir) = @ARGV;

make_path($outdir) unless -d $outdir;

# 1. 读取 rep -> family 文件，构建 family -> rep 映射
my %fam2rep;
open my $RF, "<", $rep_family_file or die $!;
while (<$RF>) {
    chomp;
    next if /^\s*$/;
    my @s = split /\t/;
    my $rep = shift @s;
    for my $fam (@s) {
        push @{$fam2rep{$fam}}, $rep;
    }
}
close $RF;

# 2. 扫描 mapping 文件并拆分
my %fh;   # rep -> file handle
open my $MAP, "<", $mapping_file or die $!;
while (<$MAP>) {
    chomp;
    next if /^\s*$/;
    my @s = split /\t/;
    my $fam = $s[1]; # mapping: fam\tgene
    for my $rep (@{$fam2rep{$fam}}) {
        if (!exists $fh{$rep}) {
            my $outfile = "$outdir/$rep.mapping";
            open my $OUT, ">", $outfile or die $!;
            $fh{$rep} = $OUT;
        }
        my $OUT = $fh{$rep};
        print $OUT "$_\n";
    }
}
close $MAP;

# 3. 关闭所有输出文件
for my $rep (keys %fh) {
    close $fh{$rep};
}
