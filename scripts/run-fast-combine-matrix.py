#!/data/software/miniconda3/envs/jinxin/bin/python3
# author: Jin-Xin Meng, jinxmeng@zju.edu.cn
# created date: 2026-04-27 15:53:32
# modified date: 2026-04-27 15:53:32

"""
Combine two-column files into a feature x sample matrix.

Example input file:
    name    value
    gene1   10
    gene2   20

Output matrix:
    name    sample1 sample2 sample3
    gene1   10      0       5
    gene2   20      1       0

Common usage:
    python combine_matrix_fast.py -D input_dir -o merged.tsv \
        --suffix .tpm -t 1 -n 1 -v 2 -s 0

Recursive usage:
    python combine_matrix_fast.py -D input_dir -o merged.tsv \
        --suffix .tpm --recursive --name-mode relpath \
        -t 1 -n 1 -v 2 -s 0
"""

import argparse
import sys
import time
from pathlib import Path
from collections import Counter


def fraction(x):
    """Convert fraction-like value to float."""
    x = str(x)
    if "/" in x:
        a, b = x.split("/", 1)
        return float(a) / float(b)
    return float(x)


def get_args():
    if len(sys.argv) == 1:
        sys.argv.append("--help")

    parser = argparse.ArgumentParser(
        description="Combine two-column files into a feature x sample matrix.",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        "-d", "--input-dir",
        required=True,
        help="Input directory"
    )

    parser.add_argument(
        "-o", "--out-file",
        required=True,
        help="Output matrix"
    )

    parser.add_argument(
        "-suffix", "--suffix",
        required=True,
        help="File suffix, e.g. .tpm, .txt, .abund"
    )

    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Search files recursively under input directory"
    )

    parser.add_argument(
        "--name-mode",
        choices=["basename", "relpath"],
        default="basename",
        help=(
            "How to name sample columns.\n"
            "basename: use file basename without suffix [default]\n"
            "relpath : use relative path without suffix; useful for recursive mode"
        )
    )

    parser.add_argument(
        "-t",
        type=int,
        default=0,
        choices=[0, 1],
        help="Whether input file has one header line. 1=yes, 0=no [default: 0]"
    )

    parser.add_argument(
        "--skip",
        type=int,
        default=0,
        help="Number of lines to skip. If >0, overrides -t"
    )

    parser.add_argument(
        "-n",
        type=int,
        default=1,
        help="Column index for feature/name, 1-based [default: 1]"
    )

    parser.add_argument(
        "-v",
        type=int,
        default=2,
        help="Column index for value, 1-based [default: 2]"
    )

    parser.add_argument(
        "-s",
        type=int,
        default=1,
        choices=[0, 1, 2, 3],
        help=(
            "Input separator [default: 1]\n"
            "0 -> tab\n"
            "1 -> whitespace\n"
            "2 -> |\n"
            "3 -> comma"
        )
    )

    parser.add_argument(
        "-f", "--format",
        default="float",
        choices=["int", "float", "fraction", "str"],
        help="Value format, used for duplicate=sum and all-fill filtering [default: float]"
    )

    parser.add_argument(
        "-a", "--fill",
        default="0",
        help="Fill value for missing entries [default: 0]"
    )

    parser.add_argument(
        "--duplicate",
        choices=["last", "sum", "error"],
        default="last",
        help=(
            "How to handle duplicated feature names in one file [default: last]\n"
            "last  -> keep the last value\n"
            "sum   -> sum duplicated values\n"
            "error -> stop when duplicated feature is found"
        )
    )

    parser.add_argument(
        "--keep-all-fill",
        action="store_true",
        help="Keep rows whose values are all fill value"
    )

    parser.add_argument(
        "--sort-rows",
        action="store_true",
        help="Sort output rows by feature name"
    )

    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Do not print progress information"
    )

    return parser.parse_args()


def get_converter(fmt):
    return {
        "int": int,
        "float": float,
        "fraction": fraction,
        "str": str,
    }[fmt]


def split_line(line, sep_code):
    """Split one line according to separator code."""
    line = line.rstrip("\n\r")

    if sep_code == 0:
        return line.split("\t")
    if sep_code == 1:
        return line.split()
    if sep_code == 2:
        return line.split("|")
    if sep_code == 3:
        return line.split(",")

    raise ValueError("unknown separator code")


def list_files(indir, suffix, recursive):
    """List eligible input files."""
    root = Path(indir)

    if not root.is_dir():
        sys.exit(f"[ERROR] Input directory not found: {indir}")

    if recursive:
        files = [
            p for p in root.rglob("*")
            if p.is_file() and p.name.endswith(suffix)
        ]
    else:
        files = [
            p for p in root.iterdir()
            if p.is_file() and p.name.endswith(suffix)
        ]

    files = sorted(files)

    if not files:
        sys.exit(f"[ERROR] No files ending with '{suffix}' found in {indir}")

    return files


def get_sample_name(path, root, suffix, name_mode):
    """Generate sample name from file path."""
    if name_mode == "basename":
        name = path.name
    else:
        name = str(path.relative_to(root))

    if name.endswith(suffix):
        name = name[:-len(suffix)]

    # Make recursive path safe as column name
    name = name.replace("/", "__").replace("\\", "__")

    return name


def print_progress(done, total, start_time, prefix):
    """Print one-line realtime progress."""
    elapsed = time.time() - start_time
    speed = done / elapsed if elapsed > 0 else 0
    eta = (total - done) / speed if speed > 0 else 0

    sys.stderr.write(
        f"\r[INFO] {prefix}: {done}/{total} "
        f"({done / total * 100:.1f}%) | "
        f"{speed:.2f}/s | ETA {eta / 60:.1f} min"
    )
    sys.stderr.flush()

    if done == total:
        sys.stderr.write("\n")


def main():
    args = get_args()

    root = Path(args.input_dir)
    files = list_files(args.input_dir, args.suffix, args.recursive)

    samples = [
        get_sample_name(p, root, args.suffix, args.name_mode)
        for p in files
    ]

    # Check duplicated sample names
    sample_count = Counter(samples)
    duplicated_samples = [s for s, c in sample_count.items() if c > 1]

    if duplicated_samples:
        sys.exit(
            "[ERROR] Duplicated sample names detected:\n"
            + "\n".join(duplicated_samples[:20])
            + "\nUse --name-mode relpath when files in different folders share the same basename."
        )

    n_samples = len(samples)
    fill = args.fill
    init_row = [fill] * n_samples

    name_col = args.n - 1
    value_col = args.v - 1
    max_col = max(name_col, value_col) # 记录需要访问的最大列号，用来判断某一行列数够不够
    skip_n = args.skip if args.skip > 0 else args.t

    convert = get_converter(args.format)

    # feature -> list of values across samples
    result = {}

    if not args.quiet:
        sys.stderr.write(f"[INFO] Input directory: {args.input_dir}\n")
        sys.stderr.write(f"[INFO] Found files: {n_samples}\n")
        sys.stderr.write(f"[INFO] Recursive: {args.recursive}\n")
        sys.stderr.write(f"[INFO] Name mode: {args.name_mode}\n")
        sys.stderr.write(f"[INFO] Suffix: {args.suffix}\n")

    start_time = time.time()

    # Parse files sequentially
    for sample_idx, path in enumerate(files):
        seen_in_file = set()

        with open(path, "r", encoding="utf-8", errors="replace") as f:
            for _ in range(skip_n):
                next(f, None)

            for line_no, line in enumerate(f, start=skip_n + 1):
                if not line.strip():
                    continue

                parts = split_line(line, args.s)

                if len(parts) <= max_col:
                    sys.exit(
                        f"\n[ERROR] {path}:{line_no} only has {len(parts)} columns. "
                        f"Check -n {args.n}, -v {args.v}, and -s {args.s}"
                    )

                feature = parts[name_col]
                value = parts[value_col]

                if feature == "":
                    continue

                if feature not in result:
                    result[feature] = init_row.copy()

                if feature in seen_in_file:
                    if args.duplicate == "error":
                        sys.exit(
                            f"\n[ERROR] Duplicated feature in one file: {feature}\n"
                            f"File: {path}\n"
                            f"Line: {line_no}"
                        )
                    elif args.duplicate == "sum":
                        old_value = result[feature][sample_idx]
                        result[feature][sample_idx] = str(
                            convert(old_value) + convert(value)
                        )
                    else:
                        result[feature][sample_idx] = value
                else:
                    result[feature][sample_idx] = value
                    seen_in_file.add(feature)

        if not args.quiet:
            print_progress(
                sample_idx + 1,
                n_samples,
                start_time,
                "Parsed files"
            )

    features = sorted(result) if args.sort_rows else result.keys()
    total_features = len(result)

    if not args.quiet:
        sys.stderr.write(f"[INFO] Total features: {total_features}\n")
        sys.stderr.write(f"[INFO] Writing output: {args.out_file}\n")

    write_start = time.time()
    written = 0
    skipped = 0

    if args.format == "str":
        fill_value_for_check = fill
    else:
        fill_value_for_check = convert(fill)

    with open(args.out_file, "w", encoding="utf-8") as out:
        out.write("name\t" + "\t".join(samples) + "\n")

        for i, feature in enumerate(features, start=1):
            row = result[feature]

            if not args.keep_all_fill:
                try:
                    if args.format == "str":
                        is_all_fill = all(x == fill_value_for_check for x in row)
                    else:
                        is_all_fill = all(convert(x) == fill_value_for_check for x in row)

                    if is_all_fill:
                        skipped += 1
                        continue

                except Exception:
                    is_all_fill = all(x == fill for x in row)
                    if is_all_fill:
                        skipped += 1
                        continue

            out.write(feature + "\t" + "\t".join(row) + "\n")
            written += 1

            if not args.quiet:
                print_progress(
                    i,
                    total_features,
                    write_start,
                    "Written rows"
                )

    if not args.quiet:
        sys.stderr.write(f"[INFO] Written features: {written}\n")
        sys.stderr.write(f"[INFO] Skipped all-fill features: {skipped}\n")
        sys.stderr.write("[INFO] Done.\n")


if __name__ == "__main__":
    main()
