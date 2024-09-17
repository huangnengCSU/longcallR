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
                intron_chr, intron_s, intron_e = fields[1].split(':')[:3]
                intron_s = int(intron_s) + 1
                intron_e = int(intron_e) - 1
                gene_id = fields[1].split(':')[-1].split('.')[0]  # use gene_id main part
                pval_nominal = float(fields[6])
                gtex_sQTL[gene_id].add(
                    f"{variant_chr}:{variant_pos}:{intron_s}:{intron_e}")  # chr:pos:intron_start:intron_end
                gtex_sQTL[gene_id].add(f"{variant_chr}:{variant_pos}:{intron_s}")  # chr:pos:intron_start
                gtex_sQTL[gene_id].add(f"{variant_chr}:{variant_pos}:{intron_e}")  # chr:pos:intron_end
                gtex_sQTL_pvalues[f"{gene_id}:{variant_chr}:{variant_pos}:{intron_s}:{intron_e}"].append(
                    (tissue_name, pval_nominal))
                gtex_sQTL_pvalues[f"{gene_id}:{variant_chr}:{variant_pos}:{intron_s}"].append(
                    (tissue_name, pval_nominal))
                gtex_sQTL_pvalues[f"{gene_id}:{variant_chr}:{variant_pos}:{intron_e}"].append(
                    (tissue_name, pval_nominal))
    return gtex_sQTL, gtex_sQTL_pvalues


def evaluate_candidate_sQTL(candidate_sQTL_file, gtex_sQTL_folder):
    gtex_sQTL_files = []
    for fname in os.listdir(gtex_sQTL_folder):
        if "sqtl_signifpairs" in fname:
            gtex_sQTL_files.append(os.path.join(gtex_sQTL_folder, fname))
    gtex_sQTL, gtex_sQTL_pvalues = load_gtex_sQTL_database(gtex_sQTL_files)
    """Load candidate sQTL database."""
    all_sQTL = set()
    double_hit_sQTL = set()
    one_side_hit_sQTL = set()
    with open(candidate_sQTL_file) as f:
        next(f)
        for line in f:
            fields = line.strip().split('\t')
            if fields[4] == "Exon" or fields[5] == "True":
                continue
            junc_s, junc_e = fields[3].split(":")[1].split("-")
            gene_id = fields[0].split('.')[0]  # use gene_id main part
            gene_name = fields[1]
            sQTL = fields[12] + ":" + junc_s + ":" + junc_e
            sQTL_left = fields[12] + ":" + junc_s
            sQTL_right = fields[12] + ":" + junc_e
            pval = float(fields[10])
            if sQTL in gtex_sQTL[gene_id]:
                double_hit_sQTL.add(sQTL)
                print(f"Hit sQTL: {gene_id}, {gene_name}, {sQTL}, {pval}")
                for tissue, pval_nominal in gtex_sQTL_pvalues[f"{gene_id}:{sQTL}"]:
                    print(f"  GTEx sQTL: {tissue}, {pval_nominal}")
            elif sQTL_left in gtex_sQTL[gene_id]:
                one_side_hit_sQTL.add(sQTL_left)
                print(f"Hit sQTL: {gene_id}, {gene_name}, {sQTL_left}, {pval}")
                for tissue, pval_nominal in gtex_sQTL_pvalues[f"{gene_id}:{sQTL_left}"]:
                    print(f"  GTEx sQTL: {tissue}, {pval_nominal}")
            elif sQTL_right in gtex_sQTL[gene_id]:
                one_side_hit_sQTL.add(sQTL_right)
                print(f"Hit sQTL: {gene_id}, {gene_name}, {sQTL_right}, {pval}")
                for tissue, pval_nominal in gtex_sQTL_pvalues[f"{gene_id}:{sQTL_right}"]:
                    print(f"  GTEx sQTL: {tissue}, {pval_nominal}")
            # else:
            #     print(f"Missed sQTL: {gene_id}, {gene_name}, {sQTL}")
            all_sQTL.add(sQTL)
    hit_sQTL = double_hit_sQTL.union(one_side_hit_sQTL)
    print(f"Total sQTL: {len(all_sQTL)}, Hit sQTL: {len(hit_sQTL)}, Double hit sQTL: {len(double_hit_sQTL)}, "
          f"One side hit sQTL: {len(one_side_hit_sQTL)}")
    return all_sQTL, hit_sQTL, double_hit_sQTL, one_side_hit_sQTL


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare candidate sQTL with GTEx sQTL.")
    parser.add_argument("-c", "--candidate_sQTL", help="Candidate sQTL file", required=True)
    parser.add_argument("-g", "--gtex_sQTL", help="Folder to GTEx sQTL files", required=True)
    args = parser.parse_args()
    evaluate_candidate_sQTL(args.candidate_sQTL, args.gtex_sQTL)
