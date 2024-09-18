import argparse
import gzip
import os
from collections import defaultdict


def load_gtex_sQTL_database(gtex_sQTL_files):
    """Load GTEx sQTL database."""
    gtex_sQTL = set()
    gtex_sQTL_pvalues = {}

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
                pval_nominal = float(fields[6])
                gtex_sQTL.add(f"{variant_chr}:{variant_pos}:{intron_s}:{intron_e}")  # chr:pos:intron_start:intron_end
                gtex_sQTL.add(f"{variant_chr}:{variant_pos}:{intron_s}")  # chr:pos:intron_start
                gtex_sQTL.add(f"{variant_chr}:{variant_pos}:{intron_e}")  # chr:pos:intron_end
                gtex_sQTL_pvalues[f"{variant_chr}:{variant_pos}:{intron_s}:{intron_e}"] = (tissue_name, pval_nominal)
                gtex_sQTL_pvalues[f"{variant_chr}:{variant_pos}:{intron_s}"] = (tissue_name, pval_nominal)
                gtex_sQTL_pvalues[f"{variant_chr}:{variant_pos}:{intron_e}"] = (tissue_name, pval_nominal)
    return gtex_sQTL, gtex_sQTL_pvalues


def evaluate_candidate_sQTL(candidate_sQTL_file, gtex_sQTL_folder):
    gtex_sQTL_files = []
    for fname in os.listdir(gtex_sQTL_folder):
        if "sqtl_signifpairs" in fname:
            gtex_sQTL_files.append(os.path.join(gtex_sQTL_folder, fname))
    gtex_sQTL, gtex_sQTL_pvalues = load_gtex_sQTL_database(gtex_sQTL_files)
    """Load candidate sQTL database."""
    all_junctions = set()
    hit_junctions = set()
    junctions_sQTLs = defaultdict(list)  # key: junction_set, value: (sQTL, junction_start, junction_end)
    with open(candidate_sQTL_file) as f:
        next(f)
        for line in f:
            fields = line.strip().split('\t')
            (junction, strand, junction_set, phase_set, hap1_absent, hap1_present, hap2_absent, hap2_present, pvalue,
             sor, novel, gene_names, sQTL, distance, ref_allele, alt_allele) = fields
            if novel == "True":
                continue
            junc_s, junc_e = junction.split(":")[1].split("-")
            junctions_sQTLs[junction_set].append((sQTL, int(junc_s), int(junc_e)))
    for junction_set in junctions_sQTLs.keys():
        found = False
        for (sQTL, junc_s, junc_e) in junctions_sQTLs[junction_set]:
            query_both_sides = f"{sQTL}:{junc_s}:{junc_e}"
            query_left = f"{sQTL}:{junc_s}"
            query_right = f"{sQTL}:{junc_e}"
            if query_both_sides in gtex_sQTL or query_left in gtex_sQTL or query_right in gtex_sQTL:
                found = True
        for (sQTL, junc_s, junc_e) in junctions_sQTLs[junction_set]:
            chr = sQTL.split(":")[0]
            if found:
                hit_junctions.add(f"{chr}:{junc_s}-{junc_e}")
            all_junctions.add(f"{chr}:{junc_s}-{junc_e}")
    print(f"Total junctions: {len(all_junctions)}, Hit junctions: {len(hit_junctions)}")
    return all_junctions, hit_junctions


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare candidate sQTL with GTEx sQTL.")
    parser.add_argument("-c", "--candidate_sQTL", help="Candidate sQTL file", required=True)
    parser.add_argument("-g", "--gtex_sQTL", help="Folder to GTEx sQTL files", required=True)
    args = parser.parse_args()
    evaluate_candidate_sQTL(args.candidate_sQTL, args.gtex_sQTL)
