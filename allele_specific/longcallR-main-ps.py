#!/usr/bin/env python

import argparse
import gzip
from collections import defaultdict

import pysam
from intervaltree import Interval, IntervalTree


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


def get_gene_regions(annotation_file, gene_types):
    assert annotation_file.endswith((".gff3", ".gtf", ".gff3.gz", ".gtf.gz")), \
        "Error: Unknown annotation file format"

    is_gff3 = ".gff3" in annotation_file
    open_func = gzip.open if annotation_file.endswith(".gz") else open
    parse_attrs = parse_attributes_gff3 if is_gff3 else parse_attributes_gtf

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
                gene_id = (attrs.get("gene_id") or attrs.get("Parent", "")).strip()
                if tx_id and gene_id:
                    transcript_to_gene[tx_id] = gene_id

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
            gene_type_by_gene[gene_id] = (attrs.get("gene_type") or attrs.get("gene_biotype", "")).strip()

    gene_regions = {}
    gene_names = {}
    exon_regions = defaultdict(lambda: defaultdict(list))
    with open_func(annotation_file, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] not in ("gene", "exon"):
                continue
            feature = parts[2]
            attrs = parse_attrs(parts[8])

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
            if gene_types and effective_gene_type not in gene_types:
                continue

            if feature == "gene":
                start, end = int(parts[3]), int(parts[4])
                gene_name = (
                    attrs.get("gene_name")
                    or (attrs.get("Name") if is_gff3 else None)
                    or (attrs.get("ID") if is_gff3 else None)
                    or gene_id
                )
                gene_regions[gene_id] = {"chr": parts[0], "start": start, "end": end}
                gene_names[gene_id] = gene_name
            else:
                transcript_id = (attrs.get("transcript_id") or attrs.get("Parent") or gene_id).strip()
                start, end = int(parts[3]), int(parts[4])
                exon_regions[gene_id][transcript_id].append((parts[0], start, end))

    return gene_regions, gene_names, exon_regions


def merge_gene_exon_regions(exon_regions):
    merged_genes_exons = defaultdict(lambda: defaultdict(list))
    for gene_id, transcripts in exon_regions.items():
        collapsed_exons = IntervalTree()
        chromosome = None
        chr_set = {chrom for exons in transcripts.values() for chrom, _, _ in exons}
        if len(chr_set) > 1:
            continue
        for exons in transcripts.values():
            for chrom, start, end in exons:
                if chromosome is None:
                    chromosome = chrom
                collapsed_exons.add(Interval(start, end + 1))
        collapsed_exons.merge_overlaps()
        merged_genes_exons[chromosome][gene_id].extend(
            sorted((interval.begin, interval.end - 1) for interval in collapsed_exons)
        )
    return merged_genes_exons


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


def read_splice_regions(read):
    splice_regions = []
    current_pos = read.reference_start + 1
    shift = 0
    for operation, length in read.cigartuples or []:
        if operation in {0, 2, 7, 8}:
            shift += length
        elif operation == 3:
            if shift > 0:
                splice_regions.append((current_pos, current_pos + shift - 1))
            current_pos += shift + length
            shift = 0
    if shift > 0:
        splice_regions.append((current_pos, current_pos + shift - 1))
    return splice_regions


def build_gene_trees(merged_genes_exons):
    trees = defaultdict(IntervalTree)
    gene_intervals = defaultdict(lambda: defaultdict(IntervalTree))
    for chrom, genes in merged_genes_exons.items():
        for gene_id, merged_exons in genes.items():
            if not merged_exons:
                continue
            gene_region = (merged_exons[0][0], merged_exons[-1][1])
            trees[chrom].add(Interval(gene_region[0], gene_region[1] + 1, gene_id))
            for exon_start, exon_end in merged_exons:
                gene_intervals[chrom][gene_id].add(Interval(exon_start, exon_end + 1))
    return trees, gene_intervals


def assign_read_to_gene(read, trees, gene_intervals):
    chrom = read.reference_name
    if chrom not in trees:
        return None
    start_pos = read.reference_start
    end_pos = read.reference_end
    overlapping = trees[chrom].overlap(start_pos + 1, end_pos + 1)
    gene_starts = {interval.data: interval.begin for interval in overlapping}
    gene_ends = {interval.data: interval.end for interval in overlapping}
    splice_regions = read_splice_regions(read)

    read_overlap_length = {}
    for gene_id in [interval.data for interval in overlapping]:
        overlap_length = 0
        for splice_start, splice_end in splice_regions:
            overlap_length += sum(
                max(0, min(splice_end, interval.end - 1) - max(splice_start, interval.begin) + 1)
                for interval in gene_intervals[chrom][gene_id].overlap(splice_start, splice_end + 1)
            )
        read_overlap_length[gene_id] = overlap_length
    return choose_best_gene(read_overlap_length, gene_starts, gene_ends)


def count_gene_phase_sets(bam_file, merged_genes_exons):
    trees, gene_intervals = build_gene_trees(merged_genes_exons)
    counts = defaultdict(lambda: defaultdict(lambda: {1: 0, 2: 0}))
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch():
            if read.is_unmapped or not read.has_tag("PS") or not read.has_tag("HP"):
                continue
            hp = read.get_tag("HP")
            if hp not in {1, 2}:
                continue
            gene_id = assign_read_to_gene(read, trees, gene_intervals)
            if gene_id is None:
                continue
            ps = read.get_tag("PS")
            counts[gene_id][ps][hp] += 1
    return counts


def main():
    parser = argparse.ArgumentParser(
        description="Select the main phase set (PS) per gene from a haplotagged BAM."
    )
    parser.add_argument("-b", "--bam", required=True, help="Haplotagged BAM file with HP/PS tags")
    parser.add_argument("-a", "--annotation", required=True, help="Annotation file in GTF or GFF3 format")
    parser.add_argument("-o", "--output", required=True, help="Output TSV path")
    parser.add_argument("-g", "--gene_types", type=str, nargs="*", default=["protein_coding", "lncRNA"],
                        help="Gene types to be analyzed. Pass --gene_types with no values to include all gene types.")
    parser.add_argument("-m", "--min_support", type=int, default=10,
                        help="Minimum phased reads required for reporting a selected main PS")
    args = parser.parse_args()

    gene_types = set(args.gene_types) if args.gene_types else set()
    gene_regions, gene_names, exon_regions = get_gene_regions(args.annotation, gene_types)
    merged_genes_exons = merge_gene_exon_regions(exon_regions)
    counts = count_gene_phase_sets(args.bam, merged_genes_exons)

    with open(args.output, "w") as out:
        out.write("#Gene_id\tGene_name\tChr\tStart\tEnd\tNum_PS\tPS\tPS_reads\t"
                  "H1\tH2\n")
        for gene_id, region in sorted(gene_regions.items(), key=lambda x: (x[1]["chr"], x[1]["start"], x[1]["end"], x[0])):
            ps_counts = counts.get(gene_id, {})
            if not ps_counts:
                continue
            main_ps, hap_counts = max(
                ps_counts.items(),
                key=lambda item: (item[1][1] + item[1][2], str(item[0])),
            )
            main_reads = hap_counts[1] + hap_counts[2]
            if main_reads < args.min_support:
                continue
            out.write(f"{gene_id}\t{gene_names.get(gene_id, gene_id)}\t{region['chr']}\t{region['start']}\t"
                      f"{region['end']}\t{len(ps_counts)}\t{main_ps}\t{main_reads}\t"
                      f"{hap_counts[1]}\t{hap_counts[2]}\n")


if __name__ == "__main__":
    main()
