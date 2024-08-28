import argparse
import gzip
import os
from collections import defaultdict


def load_gtex_sQTL_database(gtex_sQTL_files):
    """Load GTEx sQTL database."""
    gtex_sQTL = defaultdict(set)

    for gtex_sQTL_file in gtex_sQTL_files:
        if gtex_sQTL_file.endswith('.gz'):
            open_func = gzip.open
            mode = 'rt'  # Read text mode for gzip
        else:
            open_func = open
            mode = 'r'  # Read mode for regular files

        with open_func(gtex_sQTL_file, mode) as f:
            next(f)  # Skip header
            for line in f:
                fields = line.strip().split('\t')
                gene_name = fields[1]
                variant_chr = fields[13]
                variant_pos = int(fields[14])
                gtex_sQTL[gene_name].add(f"{variant_chr}:{variant_pos}")
    return gtex_sQTL


def evaluate_candidate_sQTL(candidate_sQTL_file, gtex_sQTL_folder):
    gtex_sQTL_files = []
    for fname in os.listdir(gtex_sQTL_folder):
        if "sgenes" in fname:
            gtex_sQTL_files.append(os.path.join(gtex_sQTL_folder, fname))
    gtex_sQTL = load_gtex_sQTL_database(gtex_sQTL_files)
    """Load candidate sQTL database."""
    total_sQTL_cnt = 0
    hit_sQTL_cnt = 0
    with open(candidate_sQTL_file) as f:
        next(f)
        for line in f:
            fields = line.strip().split('\t')
            gene_name = fields[1]
            sQTL = fields[5]
            if sQTL in gtex_sQTL[gene_name]:
                hit_sQTL_cnt += 1
            total_sQTL_cnt += 1
    print(f"Total sQTL: {total_sQTL_cnt}, Hit sQTL: {hit_sQTL_cnt}")
    return total_sQTL_cnt, hit_sQTL_cnt


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare candidate sQTL with GTEx sQTL.")
    parser.add_argument("-c", "--candidate_sQTL", help="Candidate sQTL file", required=True)
    parser.add_argument("-g", "--gtex_sQTL", help="Folder to GTEx sQTL files", required=True)
    args = parser.parse_args()
    evaluate_candidate_sQTL(args.candidate_sQTL, args.gtex_sQTL)
