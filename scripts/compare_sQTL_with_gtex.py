import argparse
import gzip
import os
from collections import defaultdict


def load_gtex_sQTL_database(gtex_sQTL_files):
    """Load GTEx sQTL database."""
    gtex_sQTL = defaultdict(set)
    gtex_sQTL_pvalues = defaultdict(list)

    for gtex_sQTL_file in gtex_sQTL_files:
        if gtex_sQTL_file.endswith('.gz'):
            open_func = gzip.open
            mode = 'rt'  # Read text mode for gzip
        else:
            open_func = open
            mode = 'r'  # Read mode for regular files

        # Open the file with error handling for encoding issues
        with open_func(gtex_sQTL_file, mode, encoding='utf-8', errors='replace') as f:
            tissue_name = os.path.basename(gtex_sQTL_file).split('.')[0]
            next(f)  # Skip header
            for line in f:
                fields = line.strip().split('\t')
                variant_chr = fields[0].split('_')[0]
                variant_pos = int(fields[0].split('_')[1])
                gene_id = fields[1].split(':')[-1].split('.')[0]  # use gene_id main part
                pval_nominal = float(fields[6])
                gtex_sQTL[gene_id].add(f"{variant_chr}:{variant_pos}")
                gtex_sQTL_pvalues[f"{gene_id}:{variant_chr}:{variant_pos}"].append((tissue_name, pval_nominal))
    return gtex_sQTL, gtex_sQTL_pvalues


def evaluate_candidate_sQTL(candidate_sQTL_file, gtex_sQTL_folder):
    gtex_sQTL_files = []
    for fname in os.listdir(gtex_sQTL_folder):
        if "sqtl_signifpairs" in fname:
            gtex_sQTL_files.append(os.path.join(gtex_sQTL_folder, fname))
    gtex_sQTL, gtex_sQTL_pvalues = load_gtex_sQTL_database(gtex_sQTL_files)
    """Load candidate sQTL database."""
    all_sQTL = set()
    hit_sQTL = set()
    with open(candidate_sQTL_file) as f:
        next(f)
        for line in f:
            fields = line.strip().split('\t')
            gene_id = fields[0].split('.')[0]  # use gene_id main part
            gene_name = fields[1]
            sQTL = fields[12]
            pval = float(fields[10])
            if sQTL in gtex_sQTL[gene_id]:
                hit_sQTL.add(sQTL)
                print(f"Hit sQTL: {gene_id}, {gene_name}, {sQTL}, {pval}")
                for tissue, pval_nominal in gtex_sQTL_pvalues[f"{gene_id}:{sQTL}"]:
                    print(f"  GTEx sQTL: {tissue}, {pval_nominal}")
            # else:
            #     print(f"Missed sQTL: {gene_id}, {gene_name}, {sQTL}")
            all_sQTL.add(sQTL)
    print(f"Total sQTL: {len(all_sQTL)}, Hit sQTL: {len(hit_sQTL)}")
    return all_sQTL, hit_sQTL


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare candidate sQTL with GTEx sQTL.")
    parser.add_argument("-c", "--candidate_sQTL", help="Candidate sQTL file", required=True)
    parser.add_argument("-g", "--gtex_sQTL", help="Folder to GTEx sQTL files", required=True)
    args = parser.parse_args()
    evaluate_candidate_sQTL(args.candidate_sQTL, args.gtex_sQTL)
