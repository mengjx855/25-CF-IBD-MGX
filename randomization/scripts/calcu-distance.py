#!/data/software/miniconda3/envs/jinxin/bin/python
# Jin-Xin Meng, jinxmeng@zju.edu.cn
# created date: 2025-07-18, 08:20:20
# modified date: 2026-05-13 21:15:03

# 20260406: add --rotate option.
# 20260513: fix some bug. check negative value in calculating bray or jsd distance, 
#           check sample replication, check NA, check sample number (>2).

# Usage:
#   python calcu-distance.py input.tsv output.dist.tsv --method bray --normalize TSS --transform sqrt
#
# Note:
#   1. input 默认行为 sample，列为 feature；
#   2. 如果输入为常见 profile，即 feature × sample，使用 --rotate 转置；
#   3. Hellinger 转换推荐参数：--normalize TSS --transform sqrt；
#   4. Aitchison 距离可使用：--normalize none --transform clr --method euclidean。

import sys, argparse, datetime
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform

def parse_arguments():
    '''
    Calculate distance matrix from microbiome profile.

    Recommended preprocessing:
      Count data:
        Hellinger = TSS normalization + sqrt transformation.
        Example: --normalize TSS --transform sqrt --method bray

      Compositional data:
        If already relative abundance, use --normalize none --transform sqrt.

      Aitchison-like distance:
        Use CLR transformation and Euclidean distance.
        Example: --normalize none --transform clr --method euclidean

      Binary distance:
        Jaccard distance will automatically convert data into presence/absence.
    '''
    if len(sys.argv) == 1:
        sys.argv.append('--help')

    parser = argparse.ArgumentParser(
        description=parse_arguments.__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        'in_file', type=str, 
        help='input file [data, samples × variable]')
    parser.add_argument(
        'out_file', type=str, 
        help='output file [distance matrix]')
    parser.add_argument(
        '--rotate', action='store_true', 
        help = 'rotate matrix if input a profile [default: False]')
    parser.add_argument(
        '--method', choices=['bray','jsd','cosine','euclidean','jaccard'], 
        default='bray', help='the method of distance or dissimilarity [default:bray]')
    parser.add_argument(
        '--pseudocount', metavar='',type=float, default=1e-6, 
        help='specify the pseudocount [default: 1e-6]')
    parser.add_argument(
        '--normalize', type=str, choices=['none', 'TSS', 'CSS', 'RLE'], 
        default='none', help='normalization method [default: none]')
    parser.add_argument(
        '--transform', type=str, choices=['none', 'log', 'sqrt', 'clr'], 
        default='none', help='data transformation method [default: none]')
    parser.add_argument(
        '--quiet', action='store_true', 
        help='suppress message [default: False]')
    return parser.parse_args()

def message(text, quiet=False):
    if not quiet:
        print(text)


def read_table(in_file):
    '''
    Read input table and convert values to numeric.
    '''
    data = pd.read_csv(in_file, sep='\t', index_col=0)

    if data.index.has_duplicates:
        raise ValueError('输入表的行名存在重复，请先处理重复 sample/feature ID。')

    data = data.apply(pd.to_numeric, errors='coerce')

    return data


def normalize_data(data, method, pseudocount=1e-6):
    '''
    data (pd.DataFrame): samples × features
    method (str): 归一化方法，可选值包括:
    - 'TSS': 总和标准化 (Total Sum Scaling)
    - 'CSS': 累积和标准化 (Cumulative Sum Scaling)
    - 'RLE': 相对对数表达 (Relative Log Expression)
    - 'none': 不进行归一化
    pseudocount (float): 用于处理零值的伪计数
    '''
    normalized = data.copy()

    if method == 'none':
        return normalized

    elif method == 'TSS':
        row_sums = normalized.sum(axis=1)
        row_sums[row_sums == 0] = np.nan
        normalized = normalized.div(row_sums, axis=0)

    elif method == 'CSS':
        tmp = normalized.replace(0, np.nan)
        factors = []

        for _, row in tmp.iterrows():                          # 对每个样本计算分位数阈值和缩放因子
            quantile_value = row.quantile(0.5)
            factor = row[row <= quantile_value].sum()
            factors.append(factor)

        factors = pd.Series(factors, index=normalized.index)
        factors[factors == 0] = np.nan
        normalized = normalized.div(factors, axis=0)

    elif method == 'RLE':
        tmp = normalized.replace(0, pseudocount)
        log_data = np.log(tmp)
        geo_means = log_data.mean(axis=0)                      # 在样本维度上计算平均值
        factors = []
        
        for _, row in log_data.iterrows():                     # 计算该样本相对于几何平均值的中位数差异
            median_diff = np.median(row - geo_means)
            factors.append(np.exp(median_diff))

        factors = pd.Series(factors, index=normalized.index)
        factors[factors == 0] = np.nan
        normalized = normalized.div(factors, axis=0)
        
    elif method == 'TMM':
        pass

    normalized = normalized.replace([np.inf, -np.inf], np.nan).fillna(0)

    return normalized


def transform_data(data, method, pseudocount=1e-6):
    '''
    对数据应用变换, data (pd.DataFrame): samples × features
    '''
    transformed = data.copy()

    if method == 'none':
        return transformed
    
    if method == 'sqrt':
        if (transformed.values < 0).any():
            raise ValueError('sqrt transformation does not allow negative values.')
        return np.sqrt(transformed)

    if method in ['log', 'clr']:
        if (transformed.values < 0).any():
            raise ValueError(f'{method} transformation does not allow negative values.')

        transformed = transformed.mask(transformed <= 0, pseudocount)
        log_data = np.log(transformed)

        if method == 'log':
            return log_data

        if method == 'clr':
            return log_data.subtract(log_data.mean(axis=1), axis=0)

    return transformed


def calcu_distance(data, method, pseudocount=1e-6, quiet=False):
    '''
    data (pd.DataFrame): samples × features
    '''
    if data.isnull().values.any():
        message('[Warning] 输入数据包含 NaN，已用 0 替换。', quiet)
        data = data.fillna(0)

    if data.shape[0] < 2:
        raise ValueError('样本数少于 2，无法计算距离矩阵。')

    samples = data.values.astype(np.float64)

    if method in ['bray', 'jsd'] and (samples < 0).any():
        raise ValueError(f'{method} distance does not allow negative values.')

    if method == 'bray':
        dist_vec = pdist(samples, metric='braycurtis')
    
    elif method == 'cosine':
        dist_vec = pdist(samples, metric='cosine')           # 余弦距离衡量样本向量方向的差异

    elif method == 'euclidean':
        dist_vec = pdist(samples, metric='euclidean')        # 欧氏距离计算样本在特征空间的直线距离

    elif method == 'jaccard':
        binary_samples = (samples > 0).astype(float)         # Jaccard距离衡量样本集合的相似性，将数据转换为存在/不存在(1/0)
        dist_vec = pdist(binary_samples, metric='jaccard')

    elif method == 'jsd':
        samples = np.clip(samples, pseudocount, None)        # Jensen-Shannon距离需要概率分布, 添加伪计数并归一化
        row_sums = samples.sum(axis=1, keepdims=True)

        if np.any(row_sums == 0):
            raise ValueError('JSD distance found zero-sum samples.')

        samples = samples / row_sums
        dist_vec = pdist(samples, metric='jensenshannon')

    else:
        raise ValueError(f'Unsupported distance method: {method}')

    dist_matrix = pd.DataFrame(
        squareform(dist_vec),
        index=data.index,
        columns=data.index
    )

    return dist_matrix


def main():
    start_time = datetime.datetime.now()
    args = parse_arguments()
    data = read_table(args.in_file)
    
    if args.rotate:
        data = data.T
    
    # 先标准化，再转换。例如，Hellinger = TSS + sqrt。
    data = normalize_data(data, args.normalize, args.pseudocount)
    data = transform_data(data, args.transform, args.pseudocount)

    dist_matrix = calcu_distance(
        data,
        method=args.method,
        pseudocount=args.pseudocount, 
        quiet=args.quiet
    )

    dist_matrix.to_csv(
        args.out_file, sep='\t', index=True,
        index_label='name', float_format='%.6f'
    )

    end_time = datetime.datetime.now()
    elapsed = (end_time - start_time).total_seconds()

    if not args.quiet:
        print(f'程序执行完成，耗时: {elapsed:.2f} 秒')

if __name__ == '__main__':
    main()

