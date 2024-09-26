import argparse
from collections import defaultdict

import pysam


class Variant:
    def __init__(self, chr, pos, ref_allele, alt_allele):
        self.chr = chr
        self.pos = pos
        self.ref_allele = ref_allele
        self.alt_allele = alt_allele


def load_all_variants(vcf_file):
    """Load all variants from a VCF file."""
    all_variants = defaultdict(list)
    with pysam.VariantFile(vcf_file) as vcf:
        for record in vcf:
            if record.filter.keys() != ['PASS']:
                continue
            is_indel = any(len(allele) != len(record.ref) for allele in record.alleles)
            if is_indel:
                continue
            gt = record.samples[0]['GT']
            if gt == (0, 1) or gt == (1, 0):
                chr = record.chrom
                pos = record.pos  # 1-based position
                ref_allele = record.ref
                alt_allele = record.alts[0]
                variant = Variant(chr, pos, ref_allele, alt_allele)
                if all_variants[chr]:
                    assert pos >= all_variants[chr][-1].pos  # Ensure variants are sorted by position
                all_variants[chr].append(variant)
    return all_variants


def binary_search_snp(query_position, variants_list):
    """Binary search for one closest SNP near the query position.
    query_position: 1-based position
    variants_list: list of Variant objects
    """
    left, right = 0, len(variants_list) - 1

    # Binary search to find the closest position
    while left <= right:
        mid = (left + right) // 2
        mid_pos = variants_list[mid].pos

        if mid_pos == query_position:
            return variants_list[mid]
        elif mid_pos < query_position:
            left = mid + 1
        else:
            right = mid - 1

    # Now, left is the smallest index greater than query_position
    # and right is the largest index less than query_position.
    if left >= len(variants_list):  # `query_position` is beyond the last element
        return variants_list[right]
    if right < 0:  # `query_position` is before the first element
        return variants_list[left]

    # Compare which of the two neighbors (left or right) is closer to `query_position`
    left_dist = abs(variants_list[left].pos - query_position)
    right_dist = abs(variants_list[right].pos - query_position)

    return variants_list[left] if left_dist < right_dist else variants_list[right]


def find_candidate_sQTL(vcf_file, ase_event_file, output_file, sQTL_dist_threshold):
    variants = load_all_variants(vcf_file)
    fout = open(output_file, "w")
    fout.write("#Junction\tStrand\tJunction_set\tPhase_set\tHap1_absent\tHap1_present\tHap2_absent\tHap2_present"
               "\tP_value\tSOR\tNovel\tGene_names\tsQTL\tDistance\tRef\tAlt\n")
    junctions_events = defaultdict(list)
    with open(ase_event_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            (junction, strand, junction_set, phase_set, hap1_absent, hap1_present, hap2_absent, hap2_present, pvalue,
             sor, novel, gene_names) = parts
            junctions_events[junction_set].append((junction, strand, junction_set, phase_set, hap1_absent, hap1_present,
                                                   hap2_absent, hap2_present, pvalue, sor, novel, gene_names))
    for junction_set in junctions_events.keys():
        sqtl_candidates = {}  # key: sQTL, value: (pvalue, sor)
        for event in junctions_events[junction_set]:
            (junction, strand, junction_set, phase_set, hap1_absent, hap1_present, hap2_absent, hap2_present, pvalue,
             sor, novel, gene_names) = event
            chr = junction.split(":")[0]
            start_pos, end_pos = map(int, junction.split(":")[1].split("-"))
            if len(variants[chr]) == 0:
                continue
            variant_s = binary_search_snp(start_pos, variants[chr])
            variant_e = binary_search_snp(end_pos, variants[chr])
            if abs(variant_s.pos - start_pos) < abs(variant_e.pos - end_pos):
                hit_variant = variant_s
                distance = abs(variant_s.pos - start_pos)
            else:
                hit_variant = variant_e
                distance = abs(variant_e.pos - end_pos)
            if distance <= sQTL_dist_threshold:
                sqtl_candidates[hit_variant] = (float(pvalue), float(sor))
        if len(sqtl_candidates) == 0:
            continue
        sqtl_candidates = sorted(sqtl_candidates.items(), key=lambda x: x[1][0])  # Sort by p-value
        best_variant = sqtl_candidates[0][0]
        for event in junctions_events[junction_set]:
            (junction, strand, junction_set, phase_set, hap1_absent, hap1_present, hap2_absent, hap2_present, pvalue,
             sor, novel, gene_names) = event
            start_pos, end_pos = map(int, junction.split(":")[1].split("-"))
            left_distance = start_pos - best_variant.pos  # Exon>0, Intron<=0
            right_distance = best_variant.pos - end_pos  # Exon>0, Intron<=0
            if abs(left_distance) < abs(right_distance):
                distance = left_distance
            else:
                distance = right_distance
            fout.write(f"{junction}\t{strand}\t{junction_set}\t{phase_set}\t{hap1_absent}\t{hap1_present}\t"
                       f"{hap2_absent}\t{hap2_present}\t{pvalue}\t{sor}\t{novel}\t{gene_names}\t"
                       f"{best_variant.chr}:{best_variant.pos}\t{distance}\t{best_variant.ref_allele}\t"
                       f"{best_variant.alt_allele}\n")
    fout.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Search for candidate sQTLs.")
    parser.add_argument("-v", "--vcf", type=str, required=True, help="DNA/RNA VCF file")
    parser.add_argument("-a", "--ase", type=str, required=True, help="ASE event file")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output file")
    parser.add_argument("-d", "--distance", type=int, default=4, help="sQTL distance threshold")
    args = parser.parse_args()
    find_candidate_sQTL(args.vcf, args.ase, args.output, args.distance)
