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


def find_candidate_sQTL(vcf_file, ase_event_file, output_file, sQTL_dist_threshold, sor_threshold):
    variants = load_all_variants(vcf_file)
    fout = open(output_file, "w")
    fout.write(
        "#Gene_id\tGene_name\tStrand\tEvent\tExon/Junction\tHap1_absent\tHap1_present\tHap2_absent\tHap2_present\tP_value\tSOR\tsQTL\tDistance\tRef\tAlt\n")
    with open(ase_event_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            sor = float(parts[12])
            if sor >= sor_threshold:
                gene_id, gene_name, gene_strand = parts[0], parts[1], parts[2]
                exon_junction = parts[4]
                pvalue = float(parts[11])
                hap1_absent, hap1_present, hap2_absent, hap2_present = map(int, parts[7:11])
                event = parts[3]
                chr = event.split(":")[0]
                start_pos, end_pos = map(int, event.split(":")[1].split("-"))
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
                    fout.write(f"{gene_id}\t{gene_name}\t{gene_strand}\t{event}\t{exon_junction}\t"
                               f"{hap1_absent}\t{hap1_present}\t{hap2_absent}\t{hap2_present}"
                               f"\t{pvalue}\t{sor}\t{hit_variant.chr}:{hit_variant.pos}\t{distance}\t{hit_variant.ref_allele}\t{hit_variant.alt_allele}\n")
    fout.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Search for candidate sQTLs.")
    parser.add_argument("-v", "--vcf", type=str, required=True, help="DNA/RNA VCF file")
    parser.add_argument("-a", "--ase", type=str, required=True, help="ASE event file")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output file")
    parser.add_argument("-d", "--distance", type=int, default=5, help="sQTL distance threshold")
    parser.add_argument("--sor_threshold", type=float, default=3.0,
                        help="SOR threshold for detecting sQTLs, higher values are more stringent")
    args = parser.parse_args()
    find_candidate_sQTL(args.vcf, args.ase, args.output, args.distance, args.sor_threshold)
