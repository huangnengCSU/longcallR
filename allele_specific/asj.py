import argparse
import concurrent.futures
from multiprocessing import Manager
import gzip
import math
from collections import defaultdict
import networkx as nx

import numpy as np
import pysam
from scipy.stats import fisher_exact, power_divergence, chi2
from statsmodels.stats.multitest import multipletests
from intervaltree import Interval, IntervalTree
import time


def get_gene_regions(annotation_file, gene_types):
    """Parse gene, exon, and intron regions from a GFF3 or GTF file.
    :param annotation_file: Path to the annotation file
    :return: Gene regions, exon regions, and intron regions
    """
    assert annotation_file.endswith((".gff3", ".gtf", ".gff3.gz", ".gtf.gz")), "Error: Unknown annotation file format"

    gene_regions = {}
    gene_names = {}
    gene_strands = {}
    exon_regions = defaultdict(lambda: defaultdict(list))
    intron_regions = defaultdict(lambda: defaultdict(list))

    def process_gene(parts, gene_id, gene_name):
        chr, start, end = parts[0], int(parts[3]), int(parts[4])
        gene_regions[gene_id] = {"chr": chr, "start": start, "end": end}  # 1-based, start-inclusive, end-inclusive
        gene_names[gene_id] = gene_name
        strand = parts[6]
        gene_strands[gene_id] = strand

    def process_exon(parts, gene_id, transcript_id):
        chr, start, end = parts[0], int(parts[3]), int(parts[4])
        exon_regions[gene_id][transcript_id].append((chr, start, end))  # 1-based, start-inclusive, end-inclusive

    def parse_attributes_gff3(attributes):
        attr_dict = {}
        for attr in attributes.strip().split(";"):
            key, value = attr.strip().split("=")
            attr_dict[key] = value.replace('"', '')
        return attr_dict

    def parse_attributes_gtf(attributes):
        attr_dict = {}
        for attr in attributes.strip().split(";"):
            if attr:
                key, value = attr.strip().split(" ")
                if key == "tag":
                    attr_dict[key] = attr_dict.get(key, []) + [value.replace('"', '')]
                else:
                    attr_dict[key] = value.replace('"', '')
        attr_dict["tag"] = ",".join(attr_dict.get("tag", []))
        return attr_dict

    def parse_file(file_handle, file_type):
        for line in file_handle:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            feature_type = parts[2]
            attributes = parts[8]

            if file_type == "gff3":
                attr_dict = parse_attributes_gff3(attributes)
            elif file_type == "gtf":
                attr_dict = parse_attributes_gtf(attributes)

            if feature_type == "gene":
                gene_id = attr_dict["gene_id"]
                gene_type = attr_dict["gene_type"]
                tag = attr_dict.get("tag", "")
                try:
                    gene_name = attr_dict["gene_name"]
                except KeyError:
                    gene_name = "."  # Use a placeholder if gene name is not available
                if gene_type in gene_types and "readthrough" not in tag:
                    process_gene(parts, gene_id, gene_name)
            elif feature_type == "exon":
                gene_type = attr_dict["gene_type"]
                transcript_id = attr_dict["transcript_id"]
                gene_id = attr_dict["gene_id"]
                tag = attr_dict.get("tag", "")
                if gene_type in gene_types and "readthrough" not in tag:
                    process_exon(parts, gene_id, transcript_id)

    open_func = gzip.open if annotation_file.endswith(".gz") else open
    file_type = "gff3" if ".gff3" in annotation_file else "gtf"

    with open_func(annotation_file, "rt") as f:
        parse_file(f, file_type)

    # Calculate intron regions based on exons
    for gene_id, transcripts in exon_regions.items():
        for transcript_id, exons in transcripts.items():
            if len(exons) == 1:
                continue
            exons_sorted = sorted(exons, key=lambda x: x[1])
            for i in range(1, len(exons_sorted)):
                intron_start = exons_sorted[i - 1][2] + 1
                intron_end = exons_sorted[i][1] - 1
                if intron_start < intron_end:
                    intron_regions[gene_id][transcript_id].append(
                        (exons_sorted[i - 1][0], intron_start, intron_end))  # 1-based, start-inclusive, end-inclusive

    return gene_regions, gene_names, gene_strands, exon_regions, intron_regions


def get_exon_intron_regions(read, ref_seq, no_gtag):
    exon_regions = []  # 1-based, start-inclusive, end-inclusive
    intron_regions = []  # 1-based, start-inclusive, end-inclusive, gt-ag tag
    reference_start = read.reference_start + 1  # 1-based
    current_position = reference_start
    for cigartuple in read.cigartuples:
        operation, length = cigartuple
        if operation in {0, 7, 8}:
            if exon_regions and exon_regions[-1][1] + 1 == current_position:
                # Extend the last exon region if it is contiguous
                exon_regions[-1] = (exon_regions[-1][0], exon_regions[-1][1] + length)
            else:
                # Start a new exon region
                exon_start = current_position
                exon_end = current_position + length - 1
                exon_regions.append((exon_start, exon_end))
            current_position += length
        elif operation == 2:  # 'D' operation represents deletions (still part of exon on the reference)
            if exon_regions and exon_regions[-1][1] + 1 == current_position:
                # Extend the last exon region if it is contiguous
                exon_regions[-1] = (exon_regions[-1][0], exon_regions[-1][1] + length)
            else:
                # Start a new exon region, case: 100M20N10D100M
                exon_start = current_position
                exon_end = current_position + length - 1
                exon_regions.append((exon_start, exon_end))
            current_position += length
        elif operation == 3:  # 'N' operation represents skipped region (intron)
            intron_start = current_position  # 1-based, start-inclusive
            intron_end = current_position + length - 1  # 1-based, end-inclusive
            intron_left_seq = ref_seq[intron_start - 1: intron_start + 1].upper()
            intron_right_seq = ref_seq[intron_end - 2: intron_end].upper()
            if no_gtag:
                intron_regions.append((intron_start, intron_end, False))
            else:
                if (intron_left_seq == "GT" and intron_right_seq == "AG") or (
                        intron_left_seq == "CT" and intron_right_seq == "AC"):
                    intron_regions.append((intron_start, intron_end, True))
                else:
                    intron_regions.append((intron_start, intron_end, False))
            current_position += length
        else:
            pass
    return exon_regions, intron_regions


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


def process_chunk(bam_file, chromosome, start, end, ref_seq, no_gtag, min_junctions, shared_tree, shared_gene_intervals):
    read_assignment = {}
    reads_positions = {}
    reads_tags = {}
    reads_exons = {}
    reads_junctions = {}

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch(chromosome, start, end):
            if read.is_unmapped:
                continue
            start_pos = read.reference_start
            end_pos = read.reference_end

            # get read positions and read tags
            # if not read.has_tag("HP"):
            #     continue
            # HP_tag = read.get_tag("HP")
            if read.has_tag("HP"):
                HP_tag = read.get_tag("HP")
            else:
                HP_tag = "."
            if read.has_tag("PS"):
                PS_tag = read.get_tag("PS")
            else:
                PS_tag = "."
            reads_tags[read.query_name] = {"PS": PS_tag, "HP": HP_tag}
            # get read start position and end position, 1-based, start-inclusive, end-inclusive
            reads_positions[read.query_name] = (read.reference_start + 1, read.reference_end)

            # get all exons and introns
            exon_regions, intron_regions = get_exon_intron_regions(read, ref_seq, no_gtag)
            # filter artifacts with too few junctions
            if len(intron_regions) <= min_junctions:
                del reads_positions[read.query_name]
                del reads_tags[read.query_name]
                continue
            reads_exons[read.query_name] = exon_regions
            reads_junctions[read.query_name] = intron_regions

            # query should be 1-based, left-inclusive, right-exclusive
            overlapping_intervals = shared_tree.overlap(start_pos + 1, end_pos + 1)
            candidate_gene_ids = [interval.data for interval in overlapping_intervals]

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
                        for interval in shared_gene_intervals[gene_id].overlap(*splice_region)
                    )
                read_overlap_length[gene_id] = overlap_length
            if read_overlap_length:
                best_gene_id = max(read_overlap_length, key=read_overlap_length.get)
                read_assignment[read.query_name] = best_gene_id
    return read_assignment, reads_positions, reads_tags, reads_exons, reads_junctions


def load_reads(bam_file, genome_dict, merged_genes_exons, threads, no_gtag, min_junctions=0):
    """Assign reads to genes based on their alignment positions."""

    # read_assignment, key: read_name, value: gene_id
    read_assignment = {}

    reads_positions = {}  # key: read_name, value: (start, end)
    reads_tags = {}  # key: read_name, value: {"PS": phase set, "HP": haplotype}

    # key: read_name, value: exons(list)/introns(list)
    reads_exons = {}
    reads_junctions = {}

    trees_by_chr = defaultdict(IntervalTree)  # key: chr, value: IntervalTree
    gene_intervals_by_chr = defaultdict(lambda: defaultdict(IntervalTree))  # key1: chr, key2: gene_id
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
        if chromosome not in genome_dict:
            continue
        total_length = chromosome_lengths[chromosome]  # 1-based
        chunk_size = max(1, math.ceil(total_length / threads))
        for i in range(threads):
            start = i * chunk_size  # 0-based, inclusive
            end = min((i + 1) * chunk_size, total_length)  # 0-based, exclusive
            if start >= end:
                continue
            tree = trees_by_chr[chromosome]
            gene_intervals = gene_intervals_by_chr[chromosome]
            chunks.append((bam_file, chromosome, start, end, genome_dict[chromosome], no_gtag, min_junctions, tree, gene_intervals))

    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(process_chunk, *chunk) for chunk in chunks]
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            read_assignment.update(result[0])
            reads_positions.update(result[1])
            reads_tags.update(result[2])
            reads_exons.update(result[3])
            reads_junctions.update(result[4])

    return read_assignment, reads_positions, reads_tags, reads_exons, reads_junctions


def transform_read_assignment(read_assignment):
    """select the best gene assignment for each read"""
    gene_assigned_reads = defaultdict(list)  # key: gene_id, value: list of read_name
    for read_name, gene_id in read_assignment.items():
        gene_assigned_reads[gene_id].append(read_name)
    return gene_assigned_reads

def cluster_junctions_connected_components(reads_junctions, min_count=10):
    junctions_clusters = []
    junctions = {}  # key: (start, end), value: count
    gt_ag_dict = {}  # key: (start, end), value: gt-ag tag
    for read_name, list_of_junctions in reads_junctions.items():
        for (start, end, tag) in list_of_junctions:
            junction_region = (start, end)
            junctions[junction_region] = junctions.get(junction_region, 0) + 1
            gt_ag_dict[junction_region] = tag
    junctions = {k: v for k, v in junctions.items() if v >= min_count}
    # Create a graph
    G = nx.Graph()
    # Add nodes for each junction
    for junction in junctions.keys():
        G.add_node((junction[0], junction[1], "junction"))  # "junction = (intron_start, intron_end)"
    # Add edges between overlapping junctions (sharing start or end positions)
    junction_list = [(junction[0], junction[1], "junction") for junction in list(junctions.keys())]
    merged_list = junction_list
    for i in range(len(merged_list)):
        for j in range(i + 1, len(merged_list)):
            start1, end1, type1 = merged_list[i]
            start2, end2, type2 = merged_list[j]
            # Check if they share a start or end position
            if type1 == type2:
                # junction-junction or exon-exon should share the donor or acceptor site
                if start1 == start2 or end1 == end2:
                    G.add_edge(merged_list[i], merged_list[j])
    connected_components = list(nx.connected_components(G))
    for component in connected_components:
        clu = []
        for node in component:
            if node[2] == "junction":
                gtag_tag = gt_ag_dict[(node[0], node[1])]
                clu.append((node[0], node[1], gtag_tag))
        if len(clu) > 0:
            junctions_clusters.append(clu)
    return junctions_clusters, junctions


def cluster_junctions_exons_connected_components(reads_junctions, reads_exons, min_count=10):
    junctions_clusters = []
    junctions = {}  # key: (start, end), value: count
    gt_ag_dict = {}  # key: (start, end), value: gt-ag tag
    for read_name, list_of_junctions in reads_junctions.items():
        for (start, end, tag) in list_of_junctions:
            junction_region = (start, end)
            junctions[junction_region] = junctions.get(junction_region, 0) + 1
            gt_ag_dict[junction_region] = tag
    junctions = {k: v for k, v in junctions.items() if v >= min_count}
    exons = {}  # key: (start, end), value: count
    for read_name, exon_regions in reads_exons.items():
        if len(exon_regions) == 0:
            pass
        # single exon
        if len(exon_regions) == 1:
            pass
        # two exons
        if len(exon_regions) == 2:
            pass
        # more than two exons
        if len(exon_regions) > 2:
            for i, exon_region in enumerate(exon_regions):
                if i == 0 or i == len(exon_regions) - 1:
                    continue
                exons[exon_region] = exons.get(exon_region, 0) + 1
    exons = {k: v for k, v in exons.items() if v >= min_count}
    # Create a graph
    G = nx.Graph()
    # Add nodes for each junction
    for junction in junctions.keys():
        G.add_node((junction[0], junction[1], "junction"))  # "junction = (intron_start, intron_end)"
    # Add nodes for each exon
    for exon in exons.keys():
        G.add_node((exon[0] - 1, exon[1] + 1, "exon"))  # "exon = (exon_start-1, exon_end+1)"
    # Add edges between overlapping junctions (sharing start or end positions)
    junction_list = [(junction[0], junction[1], "junction") for junction in list(junctions.keys())]
    # need to connected exon and junction, but exon and junction have one base difference
    exon_list = [(exon[0] - 1, exon[1] + 1, "exon") for exon in list(exons.keys())]
    merged_list = junction_list + exon_list
    for i in range(len(merged_list)):
        for j in range(i + 1, len(merged_list)):
            start1, end1, type1 = merged_list[i]
            start2, end2, type2 = merged_list[j]
            # Check if they share a start or end position
            if type1 == type2:
                # junction-junction or exon-exon should share the donor or acceptor site
                if start1 == start2 or end1 == end2:
                    G.add_edge(merged_list[i], merged_list[j])
            else:
                # junction-exon or exon-junction should connect to each other
                if start1 == end2 or end1 == start2:
                    G.add_edge(merged_list[i], merged_list[j])
    connected_components = list(nx.connected_components(G))
    for component in connected_components:
        clu = []
        for node in component:
            if node[2] == "junction":
                gtag_tag = gt_ag_dict[(node[0], node[1])]
                clu.append((node[0], node[1], gtag_tag))
        if len(clu) > 0:
            junctions_clusters.append(clu)
    return junctions_clusters, junctions


def check_absent_present(start_pos, end_pos, reads_positions, reads_junctions):
    """
    Find the reads where an exon or junction is absent or present.
    :param start_pos: 1-based, start-inclusive
    :param end_pos: 1-based, end-inclusive
    :param reads_positions:
    :param reads_junctions:
    :return:
    """
    absent_reads = []
    present_reads = []
    for read_name, (read_start, read_end) in reads_positions.items():
        # if read_start > start_pos or read_end < end_pos:
        #     continue
        ## Based on Max Marin case, we need to use the following condition
        if read_start > end_pos or read_end < start_pos:
            continue
        present = False
        for junction_start, junction_end, _ in reads_junctions[read_name]:
            if junction_start == start_pos and junction_end == end_pos:
                present_reads.append(read_name)
                present = True
                break
        if not present:
            absent_reads.append(read_name)
    return absent_reads, present_reads

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
            # Skip indels by checking if any alternate allele differs in length from the reference allele
            if any(len(record.ref) != len(alt) for alt in record.alts):
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
            # Skip indels by checking if any alternate allele differs in length from the reference allele
            if any(len(record.ref) != len(alt) for alt in record.alts):
                continue
            # Check for only phased heterozygous variants (0|1 or 1|0)
            if gt in [(0, 1), (1, 0)] and record.samples[0].phased:
                ps = record.samples[0].get('PS', None)
                if ps and ps!=".":
                    if with_dp_af:
                        dp = record.samples[0]['DP']
                        af = record.samples[0]['AF'][0]
                        if math.isnan(af) or math.isnan(dp) or dp == 0:
                            continue
                        rna_vcfs[ps].append(f"{record.contig}:{record.pos}:{dp}:{af}") # 1-based, ctg:pos:dp:af
                    else:
                        rna_vcfs[ps].append(f"{record.contig}:{record.pos}")  # 1-based, ctg:pos
    return rna_vcfs


class AseEvent:
    def __init__(self, chr, start, end, novel, gt_ag_tag, gene_name, strand, junction_set, phase_set,
                 hap1_absent, hap1_present, hap2_absent, hap2_present, p_value, sor):
        self.chr = chr
        self.start = start  # 1-based, inclusive
        self.end = end  # 1-based, inclusive
        self.strand = strand
        self.junction_set = junction_set
        self.phase_set = phase_set
        self.hap1_absent = hap1_absent
        self.hap1_present = hap1_present
        self.hap2_absent = hap2_absent
        self.hap2_present = hap2_present
        self.p_value = p_value
        self.sor = sor
        self.novel = novel
        self.gt_ag_tag = gt_ag_tag
        self.gene_name = gene_name

    @staticmethod
    def __header__():
        return ("#Junction\tStrand\tJunction_set\tPhase_set\tHap1_absent\tHap1_present\tHap2_absent\tHap2_present\t"
                "P_value\tSOR\tNovel\tGT_AG\tGene_name")

    def __str__(self):
        return (f"{self.chr}:{self.start}-{self.end}\t{self.strand}\t{self.junction_set}\t{self.phase_set}\t"
                f"{self.hap1_absent}\t{self.hap1_present}\t{self.hap2_absent}\t{self.hap2_present}\t"
                f"{self.p_value}\t{self.sor}\t{self.novel}\t{self.gt_ag_tag}\t{self.gene_name}")


def calc_sor(hap1_absent, hap1_present, hap2_absent, hap2_present):
    R = ((hap1_absent + 1) * (hap2_present + 1)) / ((hap1_present + 1) * (hap2_absent + 1))
    R_inverse = 1 / R
    sum = R + R_inverse
    SOR = math.log(sum)
    return SOR


def g_test_2x2(table, pseudocount=1e-10):
    """
    Perform a G-test on a 2x2 contingency table.
    Parameters:
        table (numpy.ndarray): A 2x2 contingency table.
        pseudocount (float): Small value to avoid log(0) errors for cells with zero counts.
    Returns:
        g_stat (float): The G-test statistic.
        p_value (float): The p-value for the test.
    """

    table = np.array(table)
    # Compute totals
    row_totals = table.sum(axis=1)
    col_totals = table.sum(axis=0)
    grand_total = table.sum()
    # Calculate expected frequencies
    expected = np.outer(row_totals, col_totals) / grand_total
    # Add pseudocount to avoid log(0)
    observed = table + pseudocount
    expected += pseudocount
    # Compute the G-statistic
    G = 2 * np.sum(observed * np.log(observed / expected))
    df = 1  # degrees of freedom
    p_value = 1 - chi2.cdf(G, df)
    return G, p_value


def haplotype_event_test(absent_reads, present_reads, reads_tags):
    """
    Perform Fisher's exact test to determine if the haplotype distribution is significantly different between absent and present reads.
    :param absent_reads:
    :param present_reads:
    :param reads_tags:
    :return:
    """
    hap_absent_counts = defaultdict(lambda: {1: 0, 2: 0})  # key: phase set, value: {hap1: count, hap2: count}
    hap_present_counts = defaultdict(lambda: {1: 0, 2: 0})  # key: phase set, value: {hap1: count, hap2: count}
    for read_name in absent_reads:
        hap = reads_tags[read_name]["HP"]
        phase_set = reads_tags[read_name]["PS"]
        hap_absent_counts[phase_set][hap] += 1
    for read_name in present_reads:
        hap = reads_tags[read_name]["HP"]
        phase_set = reads_tags[read_name]["PS"]
        hap_present_counts[phase_set][hap] += 1
    all_phase_sets = set(hap_absent_counts.keys()).union(set(hap_present_counts.keys()))
    # get the ps with the most reads
    ps_read_count = {}
    for ps in all_phase_sets:
        h1_a, h2_a = hap_absent_counts[ps][1], hap_absent_counts[ps][2]
        h1_p, h2_p = hap_present_counts[ps][1], hap_present_counts[ps][2]
        ps_read_count[ps] = h1_a + h2_a + h1_p + h2_p
    if ps_read_count:
        most_reads_ps = sorted(ps_read_count.items(), key=lambda x: x[1], reverse=True)[0][0]
    else:
        return None
    phase_set = most_reads_ps
    table = np.array([[hap_absent_counts[phase_set][1], hap_absent_counts[phase_set][2]],
                      [hap_present_counts[phase_set][1], hap_present_counts[phase_set][2]]])
    ## Fisher's exact test
    oddsratio, pvalue_fisher = fisher_exact(table)
    ## G-test
    # g_stat, pvalue_gtest = power_divergence(f_obs=table + 1e-300, lambda_="log-likelihood")
    # pvalue_gtest = np.min(pvalue_gtest)
    g_stat, pvalue_gtest = g_test_2x2(table)
    ## Use the maximum p-value from Fisher's exact test and G-test
    pvalue = max(pvalue_fisher, pvalue_gtest)
    ## Calculate SOR, refer to GATK AS_StrandOddsRatio, https://gatk.broadinstitute.org/hc/en-us/articles/360037224532-AS-StrandOddsRatio
    sor = calc_sor(hap_absent_counts[phase_set][1], hap_present_counts[phase_set][1],
                   hap_absent_counts[phase_set][2], hap_present_counts[phase_set][2])
    event = (phase_set, hap_absent_counts[phase_set][1], hap_present_counts[phase_set][1],
             hap_absent_counts[phase_set][2], hap_present_counts[phase_set][2], pvalue, sor)
    return event

    # for phase_set in all_phase_sets:
    #     if phase_set == ".":
    #         continue
    #     # take phased reads without phase set into account
    #     hap_absent_counts[phase_set][1] = hap_absent_counts[phase_set].get(1, 0) + hap_absent_counts["."].get(1, 0)
    #     hap_absent_counts[phase_set][2] = hap_absent_counts[phase_set].get(2, 0) + hap_absent_counts["."].get(2, 0)
    #     hap_present_counts[phase_set][1] = hap_present_counts[phase_set].get(1, 0) + hap_present_counts["."].get(1, 0)
    #     hap_present_counts[phase_set][2] = hap_present_counts[phase_set].get(2, 0) + hap_present_counts["."].get(2, 0)
    #     table = np.array([[hap_absent_counts[phase_set][1], hap_absent_counts[phase_set][2]],
    #                       [hap_present_counts[phase_set][1], hap_present_counts[phase_set][2]]])
    #     ## Fisher's exact test
    #     oddsratio, pvalue_fisher = fisher_exact(table)
    #     ## G-test
    #     g_stat, pvalue_gtest = power_divergence(f_obs=table + 1e-30, lambda_="log-likelihood")
    #     pvalue_gtest = np.min(pvalue_gtest)
    #
    #     ## Use the maximum p-value from Fisher's exact test and G-test
    #     pvalue = max(pvalue_fisher, pvalue_gtest)
    #
    #     ## Calculate SOR, refer to GATK AS_StrandOddsRatio, https://gatk.broadinstitute.org/hc/en-us/articles/360037224532-AS-StrandOddsRatio
    #     sor = calc_sor(hap_absent_counts[phase_set][1], hap_present_counts[phase_set][1],
    #                    hap_absent_counts[phase_set][2], hap_present_counts[phase_set][2])
    #     event = (phase_set, hap_absent_counts[phase_set][1], hap_present_counts[phase_set][1],
    #              hap_absent_counts[phase_set][2], hap_present_counts[phase_set][2], pvalue, sor)
    #     events.append(event)
    # return events


def analyze_gene(gene_name, gene_strand, annotation_exons, annotation_junctions, gene_region, gene_reads, min_count, cluster_with_exons):
    global reads_positions, reads_tags, reads_exons, reads_introns
    # Subset reads for this gene
    phased_read_names = [name for name in gene_reads if reads_tags[name]["HP"] != "."]
    sub_reads_positions = {name: reads_positions[name] for name in phased_read_names}
    sub_reads_tags = {name: reads_tags[name] for name in phased_read_names}
    sub_reads_exons = {name: reads_exons[name] for name in phased_read_names}
    sub_reads_introns = {name: reads_introns[name] for name in phased_read_names}


    chr = gene_region["chr"]
    start = gene_region["start"]
    end = gene_region["end"]
    gene_junction_set = set()
    for transcript_id, anno_junctions in annotation_junctions.items():
        for anno_junc in anno_junctions:
            gene_junction_set.add(anno_junc)
    gene_exon_set = set()
    for transcript_id, anno_exons in annotation_exons.items():
        for anno_exon in anno_exons:
            gene_exon_set.add(anno_exon)
    if not cluster_with_exons:
        junctions_clusters, read_junctions = cluster_junctions_connected_components(sub_reads_introns, min_count)
    else:
        junctions_clusters, read_junctions = cluster_junctions_exons_connected_components(sub_reads_introns, sub_reads_exons,
                                                                                      min_count)

    # filter reads which have no overlapped exons with current gene exons
    intervalt = IntervalTree()
    for anno_exon in gene_exon_set:
        exon_start, exon_end = anno_exon[1:3]
        intervalt.addi(exon_start, exon_end + 1)  # interval is half-open, left inclusive, right exclusive

    reads_to_remove = []
    for qname, read_exons in sub_reads_exons.items():
        overlapped = False
        for (exon_start, exon_end) in read_exons:
            if bool(intervalt.overlap(exon_start, exon_end + 1)):
                overlapped = True
                break
        if not overlapped:
            reads_to_remove.append(qname)

    # Remove reads after collecting them
    for qname in reads_to_remove:
        del sub_reads_positions[qname]
        del sub_reads_exons[qname]
        del sub_reads_introns[qname]
        del sub_reads_tags[qname]

    gene_ase_events = []  # each gene may have multiple allele-specific junctions

    # Analyze junction regions
    for junc_cluster in junctions_clusters:
        if len(junc_cluster) == 0:
            continue
        junction_set = f"{chr}:{junc_cluster[0][0]}-{junc_cluster[0][1]}"
        for read_junc in junc_cluster:
            junction_start = read_junc[0]
            junction_end = read_junc[1]
            gt_ag_tag = read_junc[2]
            novel = (chr, junction_start, junction_end) not in gene_junction_set
            # (extended_junction_start, extended_junction_end) = junctions_extended[(junction_start, junction_end)]
            # absences, presents = check_absent_present(extended_junction_start, extended_junction_end, reads_positions,
            #                                           reads_introns)
            absences, presents = check_absent_present(junction_start, junction_end, sub_reads_positions, sub_reads_introns)
            test_result = haplotype_event_test(absences, presents, sub_reads_tags)
            if test_result is None:
                continue
            (phase_set, h1_a, h1_p, h2_a, h2_p, pvalue, sor) = test_result
            gene_ase_events.append(AseEvent(chr, junction_start, junction_end, novel, gt_ag_tag, gene_name, gene_strand,
                                            junction_set, phase_set, h1_a, h1_p, h2_a, h2_p, pvalue, sor))
    return gene_ase_events


def analyze_gene_with_filtering(gene_name, gene_strand, annotation_exons, annotation_junctions, gene_region, gene_reads, min_count, cluster_with_exons):
    global reads_positions, reads_tags, reads_exons, reads_introns, dna_vcfs, rna_vcfs
    # Subset reads for this gene
    phased_read_names = [name for name in gene_reads if reads_tags[name]["HP"] != "."]
    sub_reads_positions = {name: reads_positions[name] for name in phased_read_names}
    sub_reads_tags = {name: reads_tags[name] for name in phased_read_names}
    sub_reads_exons = {name: reads_exons[name] for name in phased_read_names}
    sub_reads_introns = {name: reads_introns[name] for name in phased_read_names}

    chr = gene_region["chr"]
    start = gene_region["start"]
    end = gene_region["end"]
    gene_junction_set = set()
    for transcript_id, anno_junctions in annotation_junctions.items():
        for anno_junc in anno_junctions:
            gene_junction_set.add(anno_junc)
    gene_exon_set = set()
    for transcript_id, anno_exons in annotation_exons.items():
        for anno_exon in anno_exons:
            gene_exon_set.add(anno_exon)
    if not cluster_with_exons:
        junctions_clusters, read_junctions = cluster_junctions_connected_components(sub_reads_introns, min_count)
    else:
        junctions_clusters, read_junctions = cluster_junctions_exons_connected_components(sub_reads_introns, sub_reads_exons,
                                                                                      min_count)

    # filter reads which have no overlapped exons with current gene exons
    intervalt = IntervalTree()
    for anno_exon in gene_exon_set:
        exon_start, exon_end = anno_exon[1:3]
        intervalt.addi(exon_start, exon_end + 1)  # interval is half-open, left inclusive, right exclusive

    reads_to_remove = []

    # filter reads not phased by any DNA variants
    for qname in sub_reads_tags.keys():
        phase_set = sub_reads_tags[qname]["PS"]
        ps_variants = rna_vcfs.get(phase_set, [])
        overlapped_snps_cnt = 0
        for snp in ps_variants:
            ctg_pos = snp.split(":")[0] + ":" + snp.split(":")[1]
            if ctg_pos in dna_vcfs:
                overlapped_snps_cnt += 1
        if overlapped_snps_cnt == 0:
            reads_to_remove.append(qname)

    for qname, read_exons in sub_reads_exons.items():
        overlapped = False
        for (exon_start, exon_end) in read_exons:
            if bool(intervalt.overlap(exon_start, exon_end + 1)):
                overlapped = True
                break
        if not overlapped:
            reads_to_remove.append(qname)

    # Remove reads after collecting them
    for qname in set(reads_to_remove):
        del sub_reads_positions[qname]
        del sub_reads_exons[qname]
        del sub_reads_introns[qname]
        del sub_reads_tags[qname]

    gene_ase_events = []  # each gene may have multiple allele-specific junctions

    # Analyze junction regions
    for junc_cluster in junctions_clusters:
        if len(junc_cluster) == 0:
            continue
        junction_set = f"{chr}:{junc_cluster[0][0]}-{junc_cluster[0][1]}"
        for read_junc in junc_cluster:
            junction_start = read_junc[0]
            junction_end = read_junc[1]
            gt_ag_tag = read_junc[2]
            novel = (chr, junction_start, junction_end) not in gene_junction_set
            # (extended_junction_start, extended_junction_end) = junctions_extended[(junction_start, junction_end)]
            # absences, presents = check_absent_present(extended_junction_start, extended_junction_end, reads_positions,
            #                                           reads_introns)
            absences, presents = check_absent_present(junction_start, junction_end, sub_reads_positions, sub_reads_introns)
            test_result = haplotype_event_test(absences, presents, sub_reads_tags)
            if test_result is None:
                continue
            (phase_set, h1_a, h1_p, h2_a, h2_p, pvalue, sor) = test_result
            gene_ase_events.append(AseEvent(chr, junction_start, junction_end, novel, gt_ag_tag, gene_name, gene_strand,
                                            junction_set, phase_set, h1_a, h1_p, h2_a, h2_p, pvalue, sor))
    return gene_ase_events


# --- Global variables for large dicts ---
reads_positions = None
reads_tags = None
reads_exons = None
reads_introns = None
dna_vcfs = None
rna_vcfs = None

def analyze(annotation_file, bam_file, reference_file, output_prefix, min_count, gene_types, threads, no_gtag, min_junctions, cluster_with_exons):
    global reads_positions, reads_tags, reads_exons, reads_introns, dna_vcfs, rna_vcfs

    all_ase_events = {}  # key: (chr, start, end), value: {gene_name: AseEvent}
    start_time = time.time()
    anno_gene_regions, anno_gene_names, anno_gene_strands, anno_exon_regions, anno_intron_regions = get_gene_regions(
        annotation_file, gene_types)
    print(f"Annotation file parsed in {time.time() - start_time:.2f} seconds")
    merged_genes_exons = merge_gene_exon_regions(anno_exon_regions)
    genome_dict = {}
    ref_genome = pysam.FastaFile(reference_file)
    for chrom in ref_genome.references:
        genome_dict[chrom] = ref_genome.fetch(chrom)
    start_time = time.time()
    read_assignment, reads_positions_local, reads_tags_local, reads_exons_local, reads_introns_local = load_reads(bam_file,
                                                                                            genome_dict,
                                                                                            merged_genes_exons,
                                                                                            threads,
                                                                                            no_gtag,
                                                                                            min_junctions)
    print(f"Reads assigned to genes in {time.time() - start_time:.2f} seconds")
    gene_assigned_reads = transform_read_assignment(read_assignment)

    # Set globals ONCE
    reads_positions = reads_positions_local
    reads_tags = reads_tags_local
    reads_exons = reads_exons_local
    reads_introns = reads_introns_local

    with open(output_prefix + ".gene_coverage.tsv", "w") as f:
        f.write("#Gene_name\tChr\tStart\tEnd\tNum_reads\n")
        for gene_id, gene_region in anno_gene_regions.items():
            if gene_id not in gene_assigned_reads:
                gene_coverage = 0
            else:
                gene_coverage = len(gene_assigned_reads[gene_id])
            f.write(f"{anno_gene_names[gene_id]}\t{gene_region['chr']}\t{gene_region['start']}\t"
                    f"{gene_region['end']}\t{gene_coverage}\n")

    gene_data_list = []
    for gene_id, gene_region in anno_gene_regions.items():
        if (gene_region["chr"] in genome_dict) and (len(gene_assigned_reads[gene_id]) > 0):
            gene_name = anno_gene_names[gene_id]
            gene_strand = anno_gene_strands[gene_id]
            gene_anno_exons = anno_exon_regions[gene_id]
            gene_anno_introns = anno_intron_regions[gene_id]
            gene_reads = gene_assigned_reads[gene_id]
            gene_data_list.append(
                (gene_name, gene_strand, gene_anno_exons, gene_anno_introns, gene_region, gene_reads, min_count, cluster_with_exons))

    print(f"Total genes to be analyzed: {len(gene_data_list)}")
    start_time = time.time()
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(analyze_gene, *gene_data) for gene_data in gene_data_list]
        for future in concurrent.futures.as_completed(futures):
            gene_ase_events = future.result()
            for event in gene_ase_events:
                # multiple junctions in one gene_ase_events
                key = (event.chr, event.start, event.end)
                if key in all_ase_events.keys():
                    all_ase_events[key][event.gene_name] = event
                else:
                    all_ase_events[key] = {event.gene_name: event}
    print(f"All Gene processed completed in {time.time() - start_time:.2f} seconds")

    # apply Benjamini–Hochberg correction for all junctions with enough reads
    pass_idx = []  # index of junctions
    p_values = []
    junctions = []
    for key in all_ase_events.keys():
        for gname in all_ase_events[key].keys():
            junctions.append((key, gname))  # key: (chr, start, end), gname: gene name
    print(f"Total junctions: {len(junctions)}")
    for idx in range(len(junctions)):
        junc = junctions[idx][0]
        gname = junctions[idx][1]
        event = all_ase_events[junc][gname]
        if event.hap1_absent + event.hap1_present + event.hap2_absent + event.hap2_present >= min_count:
            pass_idx.append(idx)
            p_values.append(event.p_value)
    print(f"number of junctions with at least {min_count} reads: {len(pass_idx)}")
    reject, adjusted_p_values, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')
    asj_genes = {}
    with open(output_prefix + ".asj.tsv", "w") as f:
        f.write(AseEvent.__header__() + "\n")
        for pi in range(len(pass_idx)):
            junc = junctions[pass_idx[pi]][0]
            gname = junctions[pass_idx[pi]][1]
            event = all_ase_events[junc][gname]
            event.p_value = adjusted_p_values[pi]
            f.write(event.__str__() + "\n")
            if not no_gtag and not event.gt_ag_tag:
                continue
            if gname not in asj_genes:
                asj_genes[gname] = [event.chr, event.p_value, event.sor]
            else:
                if event.p_value < asj_genes[gname][1]:
                    asj_genes[gname] = [event.chr, event.p_value, event.sor]
    print(f"number of genes with allele-specific junctions: {len(asj_genes.keys())}")
    with open(output_prefix + ".asj_gene.tsv", "w") as f:
        f.write(f"#Gene_name\tChr\tP_value\tSOR\n")
        for gene_name in asj_genes:
            chr, pvalue, sor = asj_genes[gene_name]
            f.write(f"{gene_name}\t{chr}\t{pvalue}\t{sor}\n")

def analyze_with_filtering(annotation_file, bam_file, reference_file, output_prefix, min_count, gene_types, threads, no_gtag, min_junctions, cluster_with_exons, dna_vcfs_in, rna_vcfs_in):
    global reads_positions, reads_tags, reads_exons, reads_introns, dna_vcfs, rna_vcfs

    all_ase_events = {}  # key: (chr, start, end), value: {gene_name: AseEvent}
    start_time = time.time()
    anno_gene_regions, anno_gene_names, anno_gene_strands, anno_exon_regions, anno_intron_regions = get_gene_regions(
        annotation_file, gene_types)
    print(f"Annotation file parsed in {time.time() - start_time:.2f} seconds")
    merged_genes_exons = merge_gene_exon_regions(anno_exon_regions)
    genome_dict = {}
    ref_genome = pysam.FastaFile(reference_file)
    for chrom in ref_genome.references:
        genome_dict[chrom] = ref_genome.fetch(chrom)
    start_time = time.time()
    read_assignment, reads_positions_local, reads_tags_local, reads_exons_local, reads_introns_local = load_reads(bam_file,
                                                                                            genome_dict,
                                                                                            merged_genes_exons,
                                                                                            threads,
                                                                                            no_gtag,
                                                                                            min_junctions)
    print(f"Reads assigned to genes in {time.time() - start_time:.2f} seconds")
    gene_assigned_reads = transform_read_assignment(read_assignment)

    # Set globals ONCE
    reads_positions = reads_positions_local
    reads_tags = reads_tags_local
    reads_exons = reads_exons_local
    reads_introns = reads_introns_local
    dna_vcfs = dna_vcfs_in
    rna_vcfs = rna_vcfs_in

    with open(output_prefix + ".gene_coverage.tsv", "w") as f:
        f.write("#Gene_name\tChr\tStart\tEnd\tNum_reads\n")
        for gene_id, gene_region in anno_gene_regions.items():
            if gene_id not in gene_assigned_reads:
                gene_coverage = 0
            else:
                gene_coverage = len(gene_assigned_reads[gene_id])
            f.write(f"{anno_gene_names[gene_id]}\t{gene_region['chr']}\t{gene_region['start']}\t"
                    f"{gene_region['end']}\t{gene_coverage}\n")
    gene_data_list = []
    for gene_id, gene_region in anno_gene_regions.items():
        if (gene_region["chr"] in genome_dict) and (len(gene_assigned_reads[gene_id]) > 0):
            gene_name = anno_gene_names[gene_id]
            gene_strand = anno_gene_strands[gene_id]
            gene_anno_exons = anno_exon_regions[gene_id]
            gene_anno_introns = anno_intron_regions[gene_id]
            gene_reads = gene_assigned_reads[gene_id]
            gene_data_list.append((gene_name, gene_strand, gene_anno_exons, gene_anno_introns, gene_region, gene_reads, min_count, cluster_with_exons))
    print(f"Total genes to be analyzed: {len(gene_data_list)}")

    start_time = time.time()
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(analyze_gene_with_filtering, *gene_data) for gene_data in gene_data_list]
        for future in concurrent.futures.as_completed(futures):
            gene_ase_events = future.result()
            for event in gene_ase_events:
                # multiple junctions in one gene_ase_events
                key = (event.chr, event.start, event.end)
                if key in all_ase_events.keys():
                    all_ase_events[key][event.gene_name] = event
                else:
                    all_ase_events[key] = {event.gene_name: event}
    print(f"All Gene processed completed in {time.time() - start_time:.2f} seconds")

    # apply Benjamini–Hochberg correction for all junctions with enough reads
    pass_idx = []  # index of junctions
    p_values = []
    junctions = []
    for key in all_ase_events.keys():
        for gname in all_ase_events[key].keys():
            junctions.append((key, gname))  # key: (chr, start, end), gname: gene name
    print(f"Total junctions: {len(junctions)}")
    for idx in range(len(junctions)):
        junc = junctions[idx][0]
        gname = junctions[idx][1]
        event = all_ase_events[junc][gname]
        if event.hap1_absent + event.hap1_present + event.hap2_absent + event.hap2_present >= min_count:
            pass_idx.append(idx)
            p_values.append(event.p_value)
    print(f"number of junctions with at least {min_count} reads: {len(pass_idx)}")
    reject, adjusted_p_values, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')
    asj_genes = {}
    with open(output_prefix + ".asj.tsv", "w") as f:
        f.write(AseEvent.__header__() + "\n")
        for pi in range(len(pass_idx)):
            junc = junctions[pass_idx[pi]][0]
            gname = junctions[pass_idx[pi]][1]
            event = all_ase_events[junc][gname]
            event.p_value = adjusted_p_values[pi]
            f.write(event.__str__() + "\n")
            if not no_gtag and not event.gt_ag_tag:
                continue
            if gname not in asj_genes:
                asj_genes[gname] = [event.chr, event.p_value, event.sor]
            else:
                if event.p_value < asj_genes[gname][1]:
                    asj_genes[gname] = [event.chr, event.p_value, event.sor]
    print(f"number of genes with allele-specific junctions: {len(asj_genes.keys())}")
    with open(output_prefix + ".asj_gene.tsv", "w") as f:
        f.write(f"#Gene_name\tChr\tP_value\tSOR\n")
        for gene_name in asj_genes:
            chr, pvalue, sor = asj_genes[gene_name]
            f.write(f"{gene_name}\t{chr}\t{pvalue}\t{sor}\n")


if __name__ == "__main__":
    parse = argparse.ArgumentParser()
    parse.add_argument("-a", "--annotation_file", help="Annotation file in GFF3 or GTF format", required=True)
    parse.add_argument("-b", "--bam_file", help="BAM file", required=True)
    parse.add_argument("--dna_vcf", help="DNA VCF file", required=False)
    parse.add_argument("--rna_vcf", help="RNA VCF file", required=False)
    parse.add_argument("--min_junctions", help="Minimum number of junctions to be considered", default=2, type=int)
    parse.add_argument("--cluster_with_exons", action="store_true", help="Cluster junctions with exons", default=False)
    parse.add_argument("-f", "--reference", help="Reference genome file", required=True)
    parse.add_argument("-o", "--output_prefix",
                       help="prefix of output differential splicing file and allele-specific junctions file",
                       required=True)
    parse.add_argument("-t", "--threads", help="Number of threads", default=1, type=int)
    parse.add_argument("-g", "--gene_types", type=str, nargs="+", default=["protein_coding", "lncRNA"],
                       help='Gene types to be analyzed. Default is ["protein_coding", "lncRNA"]', )
    parse.add_argument("-m", "--min_sup", help="Minimum support of phased reads for exon or junction", default=10,
                       type=int)
    parse.add_argument("--no_gtag", action="store_true", help="Do not filter read junction with GT-AG signal")
    args = parse.parse_args()
    if args.dna_vcf and args.rna_vcf:
        dna_vcfs = load_dna_vcf(args.dna_vcf)
        rna_vcfs = load_longcallR_phased_vcf(args.rna_vcf, with_dp_af=False)
        analyze_with_filtering(args.annotation_file, args.bam_file, args.reference, args.output_prefix, args.min_sup, args.gene_types,
                args.threads, args.no_gtag, args.min_junctions, args.cluster_with_exons, dna_vcfs, rna_vcfs)
    else:
        analyze(args.annotation_file, args.bam_file, args.reference, args.output_prefix, args.min_sup, args.gene_types,
            args.threads, args.no_gtag, args.min_junctions, args.cluster_with_exons)

