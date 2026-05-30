#!/usr/bin/env perl
# author: Jin-Xin Meng <jinxmeng@zju.edu.cn>

use strict;
use warnings;
use File::Path qw(make_path);

die "Usage: $0 <rep_family.tsv> <fam_tpm.tsv> <outdir> [suffix:.fam.tpm]\n"  unless @ARGV >= 3;

my ($rep_family_file, $fam_tpm_file, $outdir, $suffix) = @ARGV;

$suffix //= '.fam.tpm';

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

# 2. 扫描 family TPM 总表并拆分
my %fh;
my %printed;
open my $TPM, "<", $fam_tpm_file or die $!;
my $header = <$TPM>;
while (<$TPM>) {
    chomp;
    next if /^\s*$/;
    my @s = split /\t/;
    my $fam = $s[0];
    for my $rep (@{$fam2rep{$fam}}) {
        if (!exists $fh{$rep}) {
            open my $O, ">", "$outdir/$rep$suffix" or die $!;
            $fh{$rep} = $O;
        }

        my $O = $fh{$rep};

        if (!exists $printed{$rep}) {
            print $O $header;
            $printed{$rep} = 1;
        }

        print $O "$_\n";
    }
}
close $TPM;

# 3. 关闭所有输出文件
for my $rep (keys %fh) {
    close $fh{$rep};
}

