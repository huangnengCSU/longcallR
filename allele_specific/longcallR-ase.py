#!/usr/bin/env python

import math
import gzip
import pysam
from collections import defaultdict
import concurrent.futures
from multiprocessing import Manager
from scipy.stats import betabinom
from statsmodels.stats.multitest import multipletests
from intervaltree import Interval, IntervalTree
import argparse


def get_gene_regions(annotation_file, gene_types):
    """Parse gene, exon, and intron regions from a GFF3 or GTF file.

    Supports both GENCODE/Ensembl style (gene_id / transcript_id on every line)
    and standard GFF3 style (ID / Parent hierarchy).
    """
    assert annotation_file.endswith((".gff3", ".gtf", ".gff3.gz", ".gtf.gz")), \
        "Error: Unknown annotation file format"

    is_gff3 = ".gff3" in annotation_file
    open_func = gzip.open if annotation_file.endswith(".gz") else open

    def parse_attributes_gff3(attributes):
        attr_dict = {}
        for attr in attributes.strip().split(";"):
            attr = attr.strip()
            if not attr or "=" not in attr:
                continue
            key, value = attr.split("=", 1)
            attr_dict[key.strip()] = value.strip().strip('"')
        return attr_dict

    def parse_attributes_gtf(attributes):
        attr_dict = {}
        tags = []
        for attr in attributes.strip().split(";"):
            attr = attr.strip()
            if not attr:
                continue
            parts = attr.split(" ", 1)
            if len(parts) != 2:
                continue
            key, value = parts[0].strip(), parts[1].strip().strip('"')
            if key == "tag":
                tags.append(value)
            else:
                attr_dict[key] = value
        attr_dict["tag"] = ",".join(tags)
        return attr_dict

    parse_attrs = parse_attributes_gff3 if is_gff3 else parse_attributes_gtf

    # Pre-pass (GFF3 only): build transcript_id → gene_id from mRNA/transcript lines
    # so exon lines can be resolved when they only carry Parent= pointing to an mRNA.
    transcript_to_gene = {}
    if is_gff3:
        with open_func(annotation_file, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9 or parts[2] not in ("mRNA", "transcript"):
                    continue
                attrs = parse_attributes_gff3(parts[8])
                tx_id = attrs.get("ID", "").strip()
                if not tx_id:
                    continue
                gene_id = (attrs.get("gene_id") or attrs.get("Parent", "")).strip()
                if gene_id:
                    transcript_to_gene[tx_id] = gene_id

    # Pre-pass (all formats): build gene_id -> gene_type map for exon fallback.
    gene_type_by_gene = {}
    with open_func(annotation_file, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attrs = parse_attrs(parts[8])
            gene_id = attrs.get("gene_id", "").strip()
            if not gene_id and is_gff3:
                gene_id = attrs.get("ID", "").strip()
            if not gene_id:
                continue
            gene_type = (attrs.get("gene_type") or attrs.get("gene_biotype", "")).strip()
            gene_type_by_gene[gene_id] = gene_type

    gene_regions  = {}
    gene_names    = {}
    gene_strands  = {}
    exon_regions  = defaultdict(lambda: defaultdict(list))
    seen_gene_types = set()

    with open_func(annotation_file, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] not in ("gene", "exon"):
                continue
            feature = parts[2]
            attrs = parse_attrs(parts[8])

            # Resolve gene_id: explicit attr → (GFF3 gene) ID → (GFF3 exon) Parent→map
            gene_id = attrs.get("gene_id", "").strip()
            if not gene_id and is_gff3:
                if feature == "gene":
                    gene_id = attrs.get("ID", "").strip()
                else:
                    tx_id = (attrs.get("transcript_id") or attrs.get("Parent", "")).strip()
                    if tx_id:
                        gene_id = transcript_to_gene.get(tx_id, tx_id)
            if not gene_id:
                continue

            gene_type = (attrs.get("gene_type") or attrs.get("gene_biotype", "")).strip()
            effective_gene_type = gene_type or gene_type_by_gene.get(gene_id, "")
            seen_gene_types.add(effective_gene_type)
            tag = attrs.get("tag", "")
            if (gene_types and effective_gene_type not in gene_types) or "readthrough" in tag:
                continue

            if feature == "gene":
                try:
                    start, end = int(parts[3]), int(parts[4])
                except ValueError:
                    continue
                if start <= 0 or end < start:
                    continue
                gene_name = (
                    attrs.get("gene_name")
                    or (attrs.get("Name") if is_gff3 else None)
                    or (attrs.get("ID") if is_gff3 else None)
                    or gene_id
                )
                gene_regions[gene_id] = {"chr": parts[0], "start": start, "end": end}
                gene_names[gene_id]   = gene_name
                gene_strands[gene_id] = parts[6]
            else:
                tx_id = (
                    attrs.get("transcript_id")
                    or (attrs.get("Parent") if is_gff3 else None)
                    or gene_id
                )
                try:
                    start, end = int(parts[3]), int(parts[4])
                except ValueError:
                    continue
                if start <= 0 or end < start:
                    continue
                exon_regions[gene_id][tx_id].append((parts[0], start, end))

    if gene_types:
        missing = sorted(gene_types - seen_gene_types)
        if missing:
            if len(missing) == len(gene_types):
                print(f"Warning [annotation]: none of the requested gene types {missing} were found.")
            else:
                print(f"Warning [annotation]: gene types {missing} not found in annotation.")

    intron_regions = defaultdict(lambda: defaultdict(list))
    for gene_id, transcripts in exon_regions.items():
        for tx_id, exons in transcripts.items():
            if len(exons) <= 1:
                continue
            exons_sorted = sorted(exons, key=lambda x: x[1])
            for i in range(1, len(exons_sorted)):
                intron_start = exons_sorted[i - 1][2] + 1
                intron_end   = exons_sorted[i][1] - 1
                if intron_start < intron_end:
                    intron_regions[gene_id][tx_id].append(
                        (exons_sorted[i - 1][0], intron_start, intron_end))

    print(f"Number of genes extracted from annotation: {len(gene_regions)}")
    return gene_regions, gene_names, gene_strands, exon_regions, intron_regions


def convert_mu_rho_to_alpha_beta(mu, rho):
    """
    Convert mean (mu) and overdispersion (rho) to alpha and beta parameters.
    """
    phi = (1 - rho) / rho - 1  # Dispersion parameter
    alpha = mu * phi
    beta = (1 - mu) * phi
    return alpha, beta


def beta_binomial_p_value(k_obs, n, mu, rho, alternative="two-sided"):
    """
    Calculate the p-value for a Beta-Binomial test using scipy.stats.betabinom.

    Parameters:
        k_obs (int): Observed reference allele count.
        n (int): Total read count (k_ref + k_alt).
        mu (float): Mean parameter (expected success proportion, e.g., 0.5).
        rho (float): Overdispersion parameter.
        alternative (str): 'two-sided', 'greater', or 'less'.

    Returns:
        float: P-value for the hypothesis test.
    """
    # Convert mu and rho to alpha and beta
    alpha, beta_param = convert_mu_rho_to_alpha_beta(mu, rho)

    # Initialize the Beta-Binomial distribution
    bb = betabinom(n, alpha, beta_param)

    # Observed PMF (probability of observing k_obs)
    p_obs = bb.pmf(k_obs)

    if alternative == "two-sided":
        # Two-sided test: Sum probabilities for all outcomes with P <= P(k_obs)
        pmf_values = [bb.pmf(k) for k in range(n + 1)]
        p_value = sum(p for p in pmf_values if p <= p_obs + 1e-15)
    elif alternative == "greater":
        # One-sided (greater): Sum probabilities for k >= k_obs
        p_value = bb.sf(k_obs - 1)  # sf = 1 - cdf
    elif alternative == "less":
        # One-sided (less): Sum probabilities for k <= k_obs
        p_value = bb.cdf(k_obs)
    else:
        raise ValueError("Invalid alternative hypothesis. Choose from 'two-sided', 'greater', 'less'.")

    return p_value


def choose_best_gene(read_overlap_length, gene_starts, gene_ends):
    candidates = [gene_id for gene_id, overlap in read_overlap_length.items() if overlap > 0]
    if not candidates:
        return None
    return min(
        candidates,
        key=lambda gene_id: (
            -read_overlap_length[gene_id],
            -gene_starts.get(gene_id, 0),
            gene_ends.get(gene_id, 0) - gene_starts.get(gene_id, 0),
            gene_id,
        ),
    )


def warn_chr_name_mismatch(chrs_a, chrs_b, module_name, label_a="annotation", label_b="BAM"):
    if not chrs_a or not chrs_b:
        return
    shared = chrs_a & chrs_b
    if len(shared) == len(chrs_a):
        return

    missing_cnt = len(chrs_a) - len(shared)
    if len(shared) == 0:
        print(f"Warning [{module_name}]: no shared chromosome names between {label_a} and {label_b}")
    else:
        print(f"Warning [{module_name}]: {missing_cnt} {label_a} chromosome names are missing in {label_b}")


def run_chr_name_checks(annotation_chrs, bam_chrs, module_name):
    warn_chr_name_mismatch(annotation_chrs, bam_chrs, module_name, label_a="annotation", label_b="BAM")




def merge_gene_exon_regions(exon_regions):
    """Merge transcript exons into gene regions."""
    # merged_genes_exons_sorted_by_start = dict()  # key: chr, value: list of sorted (collapsed_exons, gene_id, gene_name)

    # merged_genes_exons, key: chr, value: dict of gene_id: [(start, end), ..., (start, end)]
    merged_genes_exons = defaultdict(lambda: defaultdict(list))
    for gene_id, transcripts in exon_regions.items():
        collapsed_exons = IntervalTree()
        chromosome = None
        chr_set = set()
        for transcript_id, exons in transcripts.items():
            for (chr, start, end) in exons:
                chr_set.add(chr)
        if len(chr_set) > 1:
            # this gene has exons on multiple chromosomes
            continue
        # Iterate over transcripts and collect intervals
        for transcript_id, exons in transcripts.items():
            for (chr, start, end) in exons:
                if chromosome is None:
                    chromosome = chr
                else:
                    assert chromosome == chr, f"Error: Inconsistent chromosome in gene {gene_id}"
                collapsed_exons.add(Interval(start, end + 1))  # Interval is left-inclusive, right-exclusive
        # Merge overlapping intervals and adjust to 1-based closed intervals
        collapsed_exons.merge_overlaps()
        collapsed_exons = sorted((interval.begin, interval.end - 1) for interval in collapsed_exons)
        merged_genes_exons[chromosome][gene_id].extend(collapsed_exons)
    return merged_genes_exons


def assign_reads_to_gene(bam_file, merged_genes_exons):
    """Assign reads to genes based on their alignment positions."""

    # read_assignment, key: read_name, value: gene_id
    read_assignment = {}

    trees = defaultdict(IntervalTree)  # key: chr, value: IntervalTree
    gene_intervals = defaultdict(lambda: defaultdict(IntervalTree))  # key: chr, value: dict of gene_id: IntervalTree
    for chrom in merged_genes_exons.keys():
        for gene_id, merged_exons in merged_genes_exons[chrom].items():
            gene_region = (merged_exons[0][0], merged_exons[-1][1])  # 1-based, start-inclusive, end-inclusive
            trees[chrom].add(Interval(gene_region[0], gene_region[1] + 1, gene_id))

            # Build IntervalTree for exon regions within the gene
            for exon_start, exon_end in merged_exons:
                gene_intervals[chrom][gene_id].add(Interval(exon_start, exon_end + 1))

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch():
            if read.is_unmapped:
                continue
            chromosome = read.reference_name
            start_pos = read.reference_start  # 0-based, inclusive
            end_pos = read.reference_end  # 0-based, exclusive
            if chromosome not in trees:
                continue
            # query should be 1-based, left-inclusive, right-exclusive
            overlapping_intervals = trees[chromosome].overlap(start_pos + 1, end_pos + 1)
            candidate_gene_ids = [interval.data for interval in overlapping_intervals]
            gene_starts = {interval.data: interval.begin for interval in overlapping_intervals}
            gene_ends = {interval.data: interval.end for interval in overlapping_intervals}

            # parse cigar string to get the splice alignment regions
            cigar = read.cigartuples
            splice_regions = []  # list of splice alignment positions of a read, 1-based and start/end inclusive
            current_pos = start_pos + 1  # 1-based
            shift = 0
            for operation, length in cigar:
                if operation in {0, 2, 7, 8}:
                    shift += length
                elif operation == 3:
                    if shift > 0:
                        splice_regions.append((current_pos, current_pos + shift - 1))
                    current_pos += (shift + length)  # Move past the skipped region
                    shift = 0
            if shift > 0:
                splice_regions.append((current_pos, current_pos + shift - 1))

            read_overlap_length = {}  # key: gene_id, value: overlap_length
            # calculate the overlap of splice regions with gene exons to determine assignment of read to gene
            for gene_id in candidate_gene_ids:
                if gene_id not in gene_intervals[chromosome]:
                    continue
                overlap_length = 0
                for splice_region in splice_regions:
                    overlap_length += sum(
                        max(0, min(splice_region[1], interval.end - 1) - max(splice_region[0], interval.begin) + 1)
                        for interval in gene_intervals[chromosome][gene_id].overlap(splice_region[0], splice_region[1] + 1)
                    )
                read_overlap_length[gene_id] = overlap_length
            if read_overlap_length:
                best_gene_id = choose_best_gene(read_overlap_length, gene_starts, gene_ends)
                if best_gene_id is not None:
                    read_assignment[read.query_name] = best_gene_id
    return read_assignment


def process_chunk(bam_file, chromosome, start, end, shared_tree, shared_gene_intervals):
    read_assignment = {}
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch(chromosome, start, end):
            if read.is_unmapped:
                continue
            start_pos = read.reference_start
            end_pos = read.reference_end

            # query should be 1-based, left-inclusive, right-exclusive
            overlapping_intervals = shared_tree.overlap(start_pos + 1, end_pos + 1)
            candidate_gene_ids = [interval.data for interval in overlapping_intervals]
            gene_starts = {interval.data: interval.begin for interval in overlapping_intervals}
            gene_ends = {interval.data: interval.end for interval in overlapping_intervals}

            # parse cigar string to get the splice alignment regions
            cigar = read.cigartuples
            splice_regions = []  # list of splice alignment positions of a read, 1-based and start/end inclusive
            current_pos = start_pos + 1  # 1-based
            shift = 0
            for operation, length in cigar:
                if operation in {0, 2, 7, 8}:
                    shift += length
                elif operation == 3:
                    if shift > 0:
                        splice_regions.append((current_pos, current_pos + shift - 1))
                    current_pos += (shift + length)  # Move past the skipped region
                    shift = 0
            if shift > 0:
                splice_regions.append((current_pos, current_pos + shift - 1))

            read_overlap_length = {}  # key: gene_id, value: overlap_length
            # calculate the overlap of splice regions with gene exons to determine assignment of read to gene
            for gene_id in candidate_gene_ids:
                if gene_id not in shared_gene_intervals:
                    continue
                overlap_length = 0
                for splice_region in splice_regions:
                    overlap_length += sum(
                        max(0, min(splice_region[1], interval.end - 1) - max(splice_region[0], interval.begin) + 1)
                        for interval in shared_gene_intervals[gene_id].overlap(splice_region[0], splice_region[1] + 1)
                    )
                read_overlap_length[gene_id] = overlap_length
            if read_overlap_length:
                best_gene_id = choose_best_gene(read_overlap_length, gene_starts, gene_ends)
                if best_gene_id is not None:
                    read_assignment[read.query_name] = best_gene_id
    return read_assignment


def assign_reads_to_gene_parallel(bam_file, merged_genes_exons, threads):
    """Assign reads to genes based on their alignment positions."""

    # read_assignment, key: read_name, value: gene_id
    read_assignment = {}

    trees_by_chr = defaultdict(IntervalTree)  # key: chr, value: IntervalTree
    gene_intervals_by_chr = defaultdict(
        lambda: defaultdict(IntervalTree))  # key: chr, value: dict of gene_id: IntervalTree
    for chrom in merged_genes_exons.keys():
        for gene_id, merged_exons in merged_genes_exons[chrom].items():
            gene_region = (merged_exons[0][0], merged_exons[-1][1])  # 1-based, start-inclusive, end-inclusive
            trees_by_chr[chrom].add(Interval(gene_region[0], gene_region[1] + 1, gene_id))

            # Build IntervalTree for exon regions within the gene
            for exon_start, exon_end in merged_exons:
                gene_intervals_by_chr[chrom][gene_id].add(Interval(exon_start, exon_end + 1))

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        chromosomes = bam.references
        chromosome_lengths = dict(zip(bam.references, bam.lengths))

    chunks = []
    for chromosome in chromosomes:
        total_length = chromosome_lengths[chromosome]  # 1-based
        chunk_size = max(1, math.ceil(total_length / threads))
        for i in range(threads):
            start = i * chunk_size  # 0-based, inclusive
            end = min((i + 1) * chunk_size, total_length)  # 0-based, exclusive
            if start >= end:
                continue
            tree = trees_by_chr[chromosome]
            gene_intervals = gene_intervals_by_chr[chromosome]
            chunks.append((bam_file, chromosome, start, end, tree, gene_intervals))

    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(process_chunk, *chunk) for chunk in chunks]
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            read_assignment.update(result)

    return read_assignment


def transform_read_assignment(read_assignment):
    """select the best gene assignment for each read"""
    gene_assigned_reads = defaultdict(list)  # key: gene_id, value: list of read_name
    for read_name, gene_id in read_assignment.items():
        gene_assigned_reads[gene_id].append(read_name)
    return gene_assigned_reads


def load_whole_genome_phased_vcf(vcf_file):
    """
    Load phased variants from a whole-genome VCF file.
    :param vcf_file:
    :return: wg_vcfs: key is chr:pos, value is directory{genotype, paternal, maternal}
    """
    wg_vcfs = {}  # key is chr:pos, value is 0|1 or 1|0
    with pysam.VariantFile(vcf_file) as vcf:
        for record in vcf.fetch():
            gt = record.samples[0]['GT']
            # Skip non-SNPs (indels and MNPs)
            if len(record.ref) != 1 or any(len(alt) != 1 for alt in record.alts):
                continue
            # Check for only phased heterozygous variants (0|1 or 1|0)
            if gt in [(0, 1), (1, 0)] and record.samples[0].phased:
                ref_allele = record.ref
                alt_allele = record.alts[0]
                if gt == (0, 1):
                    wg_vcfs[f"{record.contig}:{record.pos}"] = {"gt": gt,
                                                                "pat": alt_allele,
                                                                "mat": ref_allele}
                else:
                    wg_vcfs[f"{record.contig}:{record.pos}"] = {"gt": gt,
                                                                "pat": ref_allele,
                                                                "mat": alt_allele}
    return wg_vcfs


def load_dna_vcf(vcf_file):
    """
    Load DNA VCF file.
    :param vcf_file:
    :return: dna_vcfs: key is chr:pos, value is directory{genotype, ref, alt}
    """
    dna_vcfs = {}
    with pysam.VariantFile(vcf_file) as vcf:
        for record in vcf.fetch():
            gt = record.samples[0]['GT']
            # Skip non-SNPs (indels and MNPs)
            if len(record.ref) != 1 or any(len(alt) != 1 for alt in record.alts):
                continue
            # Check for heterozygous variants (0/1 or 1/0 or 0|1 or 1|0)
            if gt in [(0, 1), (1, 0)]:
                ref_allele = record.ref
                alt_allele = record.alts[0]
                dna_vcfs[f"{record.contig}:{record.pos}"] = {"gt": gt,
                                                             "ref": ref_allele,
                                                             "alt": alt_allele}
    return dna_vcfs


def load_longcallR_phased_vcf(vcf_file, with_dp_af = False):
    """
    get the longcallR phased vcf
    :param vcf_file:
    :param with_dp_af: whether to include depth and allele fraction
    :return: rna_vcfs: key is ps, value is variants
    """
    import pysam
    from collections import defaultdict
    rna_vcfs = defaultdict(list)  # key is ps, value is positions
    with pysam.VariantFile(vcf_file) as vcf:
        for record in vcf.fetch():
            if 'PASS' not in record.filter.keys():
                continue
            gt = record.samples[0]['GT']
            # Skip non-SNPs (indels and MNPs)
            if len(record.ref) != 1 or any(len(alt) != 1 for alt in record.alts):
                continue
            # Check for only phased heterozygous variants (0|1 or 1|0)
            if gt in [(0, 1), (1, 0)] and record.samples[0].phased:
                ps = record.samples[0].get('PS', None)
                if ps and ps!=".":
                    ps = str(ps)
                    if with_dp_af:
                        dp = record.samples[0]['DP']
                        af = record.samples[0]['AF'][0]
                        if math.isnan(af) or math.isnan(dp) or dp == 0:
                            continue
                        rna_vcfs[ps].append(f"{record.contig}:{record.pos}:{dp}:{af}") # 1-based, ctg:pos:dp:af
                    else:
                        rna_vcfs[ps].append(f"{record.contig}:{record.pos}")  # 1-based, ctg:pos
    return rna_vcfs


def get_reads_tag(bam_file, chr, start_pos, end_pos):
    reads_tag = {}
    with pysam.AlignmentFile(bam_file, "rb") as f:
        for read in f.fetch(chr, start_pos - 1, end_pos):  # start_pos is 1-based, pysam expects 0-based
            PS = str(read.get_tag("PS")) if read.has_tag("PS") else None
            HP = read.get_tag("HP") if read.has_tag("HP") else None
            reads_tag[read.query_name] = {"PS": PS, "HP": HP}
    return reads_tag


def load_main_ps_table(main_ps_file):
    main_ps = {}
    if main_ps_file is None:
        return main_ps
    with open(main_ps_file) as f:
        header = None
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if header is None:
                if fields[0].startswith("#"):
                    fields[0] = fields[0].lstrip("#")
                header = {name: idx for idx, name in enumerate(fields)}
                continue
            ps_col = header.get("Main_PS", header.get("PS"))
            if ps_col is None or ps_col >= len(fields):
                continue
            ps = fields[ps_col]
            if ps in {"", "."}:
                continue
            key = None
            for key_col in ("Gene_id", "gene_id"):
                idx = header.get(key_col)
                if idx is not None and idx < len(fields) and fields[idx]:
                    key = fields[idx]
                    break
            if key is None:
                for key_col in ("Gene_name", "gene_name"):
                    idx = header.get(key_col)
                    if idx is not None and idx < len(fields) and fields[idx]:
                        key = fields[idx]
                        break
            if key is not None:
                main_ps[key] = ps
    return main_ps


def selected_main_ps(main_ps_by_gene, gene_id, gene_name):
    if main_ps_by_gene is None:
        return None
    return main_ps_by_gene.get(gene_id) or main_ps_by_gene.get(gene_name)


def choose_phase_set(phase_set_hap_count, requested_ps=None):
    if requested_ps is not None:
        return requested_ps, phase_set_hap_count.get(requested_ps, {1: 0, 2: 0})
    ps_read_cnt = {
        ps: hap_cnt[1] + hap_cnt[2]
        for ps, hap_cnt in phase_set_hap_count.items()
    }
    if not ps_read_cnt:
        return None, None
    most_reads_ps = sorted(ps_read_cnt.items(), key=lambda x: x[1], reverse=True)[0][0]
    return most_reads_ps, phase_set_hap_count[most_reads_ps]


def calculate_ase_pvalue(bam_file, gene_id, gene_name, gene_region, min_count, overdispersion, gene_assigned_reads,
                         main_ps_by_gene=None):
    reads_tag = get_reads_tag(bam_file, gene_region["chr"], gene_region["start"], gene_region["end"])
    assigned_reads = set(gene_assigned_reads[gene_id])
    phase_set_hap_count = defaultdict(lambda: {1: 0, 2: 0})  # key: phase set, value: {haplotype: count}
    for rname in assigned_reads:
        if rname in reads_tag:
            ps = reads_tag[rname]["PS"]
            hp = reads_tag[rname]["HP"]
            if ps and hp:
                phase_set_hap_count[ps][hp] += 1
    requested_ps = selected_main_ps(main_ps_by_gene, gene_id, gene_name)
    most_reads_ps, hap_count = choose_phase_set(phase_set_hap_count, requested_ps)
    if most_reads_ps is None:
        return (gene_name, gene_region["chr"], gene_region["start"], gene_region["end"], 1.0, ".", 0, 0)
    if hap_count[1] + hap_count[2] < min_count:
        return (gene_name, gene_region["chr"], gene_region["start"], gene_region["end"], 1.0, most_reads_ps, 0, 0)
    # p_value_ase = binomtest(hap_count[1], hap_count[1] + hap_count[2], 0.5, alternative='two-sided').pvalue
    p_value_ase = beta_binomial_p_value(hap_count[1], hap_count[1] + hap_count[2], 0.5, overdispersion,
                                        alternative='two-sided')
    return (gene_name, gene_region["chr"], gene_region["start"], gene_region["end"], p_value_ase, most_reads_ps,
            hap_count[1], hap_count[2])


def calculate_ase_pvalue_pat_mat(bam_file, gene_id, gene_name, gene_region, min_count, overdispersion,
                                 gene_assigned_reads, rna_vcfs, wg_vcfs, main_ps_by_gene=None):
    reads_tag = get_reads_tag(bam_file, gene_region["chr"], gene_region["start"], gene_region["end"])
    assigned_reads = set(gene_assigned_reads[gene_id])
    phase_set_hap_count = defaultdict(lambda: {1: 0, 2: 0})  # key: phase set, value: {haplotype: count}
    for rname in assigned_reads:
        if rname in reads_tag:
            ps = reads_tag[rname]["PS"]
            hp = reads_tag[rname]["HP"]
            if ps and hp:
                phase_set_hap_count[ps][hp] += 1
    requested_ps = selected_main_ps(main_ps_by_gene, gene_id, gene_name)
    most_reads_ps, hap_count = choose_phase_set(phase_set_hap_count, requested_ps)
    if most_reads_ps is None:
        return (gene_name, gene_region["chr"], gene_region["start"], gene_region["end"], 1.0, ".", 0, 0, 0, 0, 0,
                0)
    h1_count, h2_count = hap_count[1], hap_count[2]
    if h1_count + h2_count < min_count:
        return (gene_name, gene_region["chr"], gene_region["start"], gene_region["end"], 1.0, ".", 0, 0, 0, 0, 0,
                0)
    # p_value_ase = binomtest(h1_count, h1_count + h2_count, 0.5, alternative='two-sided').pvalue
    p_value_ase = beta_binomial_p_value(hap_count[1], hap_count[1] + hap_count[2], 0.5, overdispersion,
                                        alternative='two-sided')

    # determine paternal and maternal alleles for H1 and H2
    ps_variants = rna_vcfs.get(most_reads_ps, [])
    ps_reads = [rname for rname in assigned_reads if
                rname in reads_tag and reads_tag[rname]["PS"] == most_reads_ps]
    h1_reads = [rname for rname in ps_reads if reads_tag[rname]["HP"] == 1]
    h2_reads = [rname for rname in ps_reads if reads_tag[rname]["HP"] == 2]

    ps_variant_pos = [
        int(pos.split(":")[1]) - 1
        for pos in ps_variants
        if pos.split(":")[0] == gene_region["chr"]
    ]  # 0-based
    reads_pat_mat_cnt = defaultdict(
        lambda: {"pat": 0, "mat": 0})  # key: read name, value: {paternal: count, maternal: count}
    for pileupcolumn in pysam.AlignmentFile(bam_file, "rb").pileup(gene_region["chr"], gene_region["start"] - 1,
                                                                   gene_region["end"], max_depth=100000):
        if pileupcolumn.pos not in ps_variant_pos or f"{gene_region['chr']}:{pileupcolumn.pos + 1}" not in wg_vcfs:
            continue
        for pileup_read in pileupcolumn.pileups:
            if not pileup_read.is_del and not pileup_read.is_refskip:
                read_name = pileup_read.alignment.query_name
                if read_name in ps_reads:
                    base = pileup_read.alignment.query_sequence[pileup_read.query_position].upper()
                    # wg_vcfs is 1-based, pileupcolumn.pos is 0-based
                    if base == wg_vcfs[f"{gene_region['chr']}:{pileupcolumn.pos + 1}"]["pat"]:
                        reads_pat_mat_cnt[read_name]["pat"] += 1
                    elif base == wg_vcfs[f"{gene_region['chr']}:{pileupcolumn.pos + 1}"]["mat"]:
                        reads_pat_mat_cnt[read_name]["mat"] += 1
                    else:
                        continue
    h1_pat_cnt, h1_mat_cnt = 0, 0
    for reads in h1_reads:
        if reads in reads_pat_mat_cnt:
            if reads_pat_mat_cnt[reads]["pat"] > reads_pat_mat_cnt[reads]["mat"]:
                h1_pat_cnt += 1
            elif reads_pat_mat_cnt[reads]["pat"] < reads_pat_mat_cnt[reads]["mat"]:
                h1_mat_cnt += 1
            else:
                continue
    h2_pat_cnt, h2_mat_cnt = 0, 0
    for reads in h2_reads:
        if reads in reads_pat_mat_cnt:
            if reads_pat_mat_cnt[reads]["pat"] > reads_pat_mat_cnt[reads]["mat"]:
                h2_pat_cnt += 1
            elif reads_pat_mat_cnt[reads]["pat"] < reads_pat_mat_cnt[reads]["mat"]:
                h2_mat_cnt += 1
            else:
                continue
    return (gene_name, gene_region["chr"], gene_region["start"], gene_region["end"], p_value_ase, most_reads_ps,
            h1_count, h2_count, h1_pat_cnt, h1_mat_cnt, h2_pat_cnt, h2_mat_cnt)


def calculate_ase_pvalue_filtering(bam_file, gene_id, gene_name, gene_region, min_count, overdispersion,
                                 gene_assigned_reads, rna_vcfs, dna_vcfs, main_ps_by_gene=None):
    reads_tag = get_reads_tag(bam_file, gene_region["chr"], gene_region["start"], gene_region["end"])
    assigned_reads = set(gene_assigned_reads[gene_id])
    phase_set_hap_count = defaultdict(lambda: {1: 0, 2: 0})  # key: phase set, value: {haplotype: count}
    for rname in assigned_reads:
        if rname in reads_tag:
            ps = reads_tag[rname]["PS"]
            hp = reads_tag[rname]["HP"]
            if ps and hp:
                phase_set_hap_count[ps][hp] += 1
    requested_ps = selected_main_ps(main_ps_by_gene, gene_id, gene_name)
    most_reads_ps, hap_count = choose_phase_set(phase_set_hap_count, requested_ps)
    if most_reads_ps is None:
        return (gene_name, gene_region["chr"], gene_region["start"], gene_region["end"], 1.0, ".", 0, 0)
    h1_count, h2_count = hap_count[1], hap_count[2]
    if h1_count + h2_count < min_count:
        return (gene_name, gene_region["chr"], gene_region["start"], gene_region["end"], 1.0, most_reads_ps, 0, 0)
    # p_value_ase = binomtest(h1_count, h1_count + h2_count, 0.5, alternative='two-sided').pvalue
    p_value_ase = beta_binomial_p_value(hap_count[1], hap_count[1] + hap_count[2], 0.5, overdispersion,
                                        alternative='two-sided')

    # filtering if all RNA heterozygous variants are not in DNA VCF
    overlapped_cnt = 0
    ps_variants = rna_vcfs.get(most_reads_ps, [])
    for snp in ps_variants:
        ctg_pos = snp.split(":")[0] + ":" + snp.split(":")[1]
        if ctg_pos in dna_vcfs:
            depth = int(snp.split(":")[2])
            allele_fraction = float(snp.split(":")[3])
            alt_cnt = int(round(depth * allele_fraction))
            alt_cnt = max(0, min(depth, alt_cnt))
            p_value_ase_allele = beta_binomial_p_value(alt_cnt, depth, 0.5, overdispersion, alternative='two-sided')
            if depth >= min_count and p_value_ase_allele < 0.05:
                overlapped_cnt += 1
    if overlapped_cnt == 0:
        return (gene_name, gene_region["chr"], gene_region["start"], gene_region["end"], 1.0, ".", 0, 0)
    return (gene_name, gene_region["chr"], gene_region["start"], gene_region["end"], p_value_ase, most_reads_ps,
            hap_count[1], hap_count[2])


def analyze_ase_genes(annotation_file, bam_file, out_file, threads, gene_types, min_support, overdispersion,
                      main_ps_by_gene=None,
                      ):
    gene_regions, gene_names, gene_strands, exon_regions, intron_regions = get_gene_regions(annotation_file, gene_types)
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        bam_chrs = set(bam.references)
    annotation_chrs = {r["chr"] for r in gene_regions.values()}
    run_chr_name_checks(annotation_chrs, bam_chrs, "ase")
    merged_genes_exons = merge_gene_exon_regions(exon_regions)
    read_assignment = assign_reads_to_gene_parallel(bam_file, merged_genes_exons, threads)
    gene_assigned_reads = transform_read_assignment(read_assignment)
    gene_args = [
        (bam_file, gene_id, gene_names[gene_id], gene_regions[gene_id], min_support, overdispersion)
        for gene_id in gene_regions.keys()
        if gene_id in gene_assigned_reads
        and (main_ps_by_gene is None or selected_main_ps(main_ps_by_gene, gene_id, gene_names[gene_id]) is not None)
    ]
    results = []
    with Manager() as manager:
        shared_assignments = manager.dict(gene_assigned_reads)
        shared_main_ps = manager.dict(main_ps_by_gene or {})
        with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
            futures = [executor.submit(calculate_ase_pvalue, *gene_data, shared_assignments, shared_main_ps)
                       for gene_data in gene_args]
            for future in concurrent.futures.as_completed(futures):
                results.append(future.result())
    # apply Benjamini–Hochberg correction for all genes with at least min_support reads
    pass_idx = []
    p_values = []
    print(f"total number of genes: {len(results)}")
    for idx, (gene_name, chrom, gene_start, gene_end, p_value, ps, h1, h2) in enumerate(results):
        if h1 + h2 >= min_support:
            pass_idx.append(idx)
            p_values.append(p_value)
    print(f"number of genes with at least {min_support} reads: {len(pass_idx)}")
    if p_values:
        reject, adjusted_p_values, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')
    else:
        adjusted_p_values = []
    with open(out_file, "w") as f:
        f.write("#Gene_name\tChr\tGene_start\tGene_end\tPS\tH1\tH2\tP_value\tlogFC\n")
        for pi in range(len(pass_idx)):
            idx = pass_idx[pi]
            gene_name, chrom, gene_start, gene_end, p_value, ps, h1, h2 = results[idx]
            p_value = adjusted_p_values[pi]
            logFC = math.log2((h1 + 1) / (h2 + 1))
            f.write(f"{gene_name}\t{chrom}\t{gene_start}\t{gene_end}\t{ps}\t{h1}\t{h2}\t{p_value}\t{logFC}\n")


def analyze_ase_genes_pat_mat(annotation_file, bam_file, vcf_file1, vcf_file2, out_file, threads, gene_types,
                              min_support, overdispersion, main_ps_by_gene=None):
    rna_vcfs = load_longcallR_phased_vcf(vcf_file1)
    wg_vcfs = load_whole_genome_phased_vcf(vcf_file2)
    gene_regions, gene_names, gene_strands, exon_regions, intron_regions = get_gene_regions(annotation_file, gene_types)
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        bam_chrs = set(bam.references)
    annotation_chrs = {r["chr"] for r in gene_regions.values()}
    run_chr_name_checks(annotation_chrs, bam_chrs, "ase")
    merged_genes_exons = merge_gene_exon_regions(exon_regions)
    read_assignment = assign_reads_to_gene_parallel(bam_file, merged_genes_exons, threads)
    gene_assigned_reads = transform_read_assignment(read_assignment)
    gene_args = [
        (bam_file, gene_id, gene_names[gene_id], gene_regions[gene_id], min_support, overdispersion)
        for gene_id in gene_regions.keys()
        if gene_id in gene_assigned_reads
        and (main_ps_by_gene is None or selected_main_ps(main_ps_by_gene, gene_id, gene_names[gene_id]) is not None)
    ]
    results = []
    with Manager() as manager:
        shared_assignments = manager.dict(gene_assigned_reads)
        shared_rna_vcfs = manager.dict(rna_vcfs)
        shared_wg_vcfs = manager.dict(wg_vcfs)
        shared_main_ps = manager.dict(main_ps_by_gene or {})
        with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
            futures = [executor.submit(calculate_ase_pvalue_pat_mat, *gene_data, shared_assignments, shared_rna_vcfs,
                                       shared_wg_vcfs, shared_main_ps) for gene_data in gene_args]
            for future in concurrent.futures.as_completed(futures):
                results.append(future.result())
    # apply Benjamini–Hochberg correction
    pass_idx = []
    p_values = []
    print(f"total number of genes: {len(results)}")
    for idx, (gene_name, chrom, gene_start, gene_end, p_value, ps, h1, h2, h1_pat, h1_mat, h2_pat,
              h2_mat) in enumerate(results):
        if h1 + h2 >= min_support:
            pass_idx.append(idx)
            p_values.append(p_value)
    print(f"number of genes with at least {min_support} reads: {len(pass_idx)}")
    if p_values:
        reject, adjusted_p_values, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')
    else:
        adjusted_p_values = []
    with open(out_file, "w") as f:
        f.write("#Gene_name\tChr\tGene_start\tGene_end\tPS\tH1\tH2\tP_value\tH1_Paternal\tH1_Maternal\tH2_Paternal\tH2_Maternal\tlogFC\n")
        for pi in range(len(pass_idx)):
            idx = pass_idx[pi]
            gene_name, chrom, gene_start, gene_end, p_value, ps, h1, h2, h1_pat, h1_mat, h2_pat, h2_mat = results[idx]
            p_value = adjusted_p_values[pi]
            logFC = math.log2((h1 + 1) / (h2 + 1))
            f.write(
                f"{gene_name}\t{chrom}\t{gene_start}\t{gene_end}\t{ps}\t{h1}\t{h2}\t{p_value}\t{h1_pat}\t{h1_mat}\t{h2_pat}\t{h2_mat}\t{logFC}\n")


def analyze_ase_genes_with_filtering(annotation_file, bam_file, vcf_file1, vcf_file3, out_file, threads, gene_types,
                              min_support, overdispersion, main_ps_by_gene=None):
    rna_vcfs = load_longcallR_phased_vcf(vcf_file1, with_dp_af=True)
    dna_vcfs = load_dna_vcf(vcf_file3)
    gene_regions, gene_names, gene_strands, exon_regions, intron_regions = get_gene_regions(annotation_file, gene_types)
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        bam_chrs = set(bam.references)
    annotation_chrs = {r["chr"] for r in gene_regions.values()}
    run_chr_name_checks(annotation_chrs, bam_chrs, "ase")
    merged_genes_exons = merge_gene_exon_regions(exon_regions)
    read_assignment = assign_reads_to_gene_parallel(bam_file, merged_genes_exons, threads)
    gene_assigned_reads = transform_read_assignment(read_assignment)
    gene_args = [
        (bam_file, gene_id, gene_names[gene_id], gene_regions[gene_id], min_support, overdispersion)
        for gene_id in gene_regions.keys()
        if gene_id in gene_assigned_reads
        and (main_ps_by_gene is None or selected_main_ps(main_ps_by_gene, gene_id, gene_names[gene_id]) is not None)
    ]
    results = []
    with Manager() as manager:
        shared_assignments = manager.dict(gene_assigned_reads)
        shared_rna_vcfs = manager.dict(rna_vcfs)
        shared_dna_vcfs = manager.dict(dna_vcfs)
        shared_main_ps = manager.dict(main_ps_by_gene or {})
        with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
            futures = [executor.submit(calculate_ase_pvalue_filtering, *gene_data, shared_assignments, shared_rna_vcfs,
                                       shared_dna_vcfs, shared_main_ps) for gene_data in gene_args]
            for future in concurrent.futures.as_completed(futures):
                results.append(future.result())
    # apply Benjamini–Hochberg correction for all genes with at least min_support reads
    pass_idx = []
    p_values = []
    print(f"total number of genes: {len(results)}")
    for idx, (gene_name, chrom, gene_start, gene_end, p_value, ps, h1, h2) in enumerate(results):
        if h1 + h2 >= min_support:
            pass_idx.append(idx)
            p_values.append(p_value)
    print(f"number of genes with at least {min_support} reads: {len(pass_idx)}")
    if p_values:
        reject, adjusted_p_values, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')
    else:
        adjusted_p_values = []
    with open(out_file, "w") as f:
        f.write("#Gene_name\tChr\tStart\tEnd\tPS\tH1\tH2\tP_value\tlogFC\n")
        for pi in range(len(pass_idx)):
            idx = pass_idx[pi]
            gene_name, chrom, gene_start, gene_end, p_value, ps, h1, h2 = results[idx]
            p_value = adjusted_p_values[pi]
            logFC = math.log2((h1 + 1) / (h2 + 1))
            f.write(f"{gene_name}\t{chrom}\t{gene_start}\t{gene_end}\t{ps}\t{h1}\t{h2}\t{p_value}\t{logFC}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--bam", required=True, help="phased BAM file")
    parser.add_argument("--vcf1", required=False, help="LongcallR phased vcf file", default=None)
    parser.add_argument("--vcf2", required=False, help="Whole genome haplotype phased DNA vcf file", default=None)
    parser.add_argument("--vcf3", required=False, help="DNA vcf file", default=None)
    parser.add_argument("-a", "--annotation", required=True, help="Annotation file")
    parser.add_argument("-d", "--overdispersion", type=float, default=0.001, help="Overdispersion parameter")
    parser.add_argument("-o", "--output", required=True, help="prefix of output file")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads")
    parser.add_argument("-g", "--gene_types", type=str, nargs="*", default=["protein_coding", "lncRNA"],
                        help='Gene types to be analyzed (default: protein_coding lncRNA). '
                             'Pass --gene_types with no values to include all gene types.', )
    parser.add_argument("-m", "--min_support", type=int, default=10,
                        help="Minimum support reads for counting event (default: 10)")
    parser.add_argument("--main_ps", default=None,
                        help="TSV of user-selected/main PS per gene. Recognized columns: Gene_id/Gene_name and Main_PS/PS")

    args = parser.parse_args()

    gene_types = set(args.gene_types)
    print(f"Gene types: {'all' if not gene_types else ', '.join(sorted(gene_types))}")
    main_ps_by_gene = load_main_ps_table(args.main_ps) if args.main_ps else None
    if main_ps_by_gene is not None:
        print(f"Loaded main PS genes: {len(main_ps_by_gene)}")

    if args.vcf1 and args.vcf2:
        analyze_ase_genes_pat_mat(args.annotation, args.bam, args.vcf1, args.vcf2, args.output + ".patmat_ase.tsv",
                                  args.threads, gene_types, args.min_support, args.overdispersion, main_ps_by_gene)
    elif args.vcf1 and args.vcf3:
        analyze_ase_genes_with_filtering(args.annotation, args.bam, args.vcf1, args.vcf3, args.output + ".filter_ase.tsv",
                                  args.threads, gene_types, args.min_support, args.overdispersion, main_ps_by_gene)
    else:
        analyze_ase_genes(args.annotation, args.bam, args.output + ".ase.tsv", args.threads, gene_types,
                          args.min_support, args.overdispersion, main_ps_by_gene)
