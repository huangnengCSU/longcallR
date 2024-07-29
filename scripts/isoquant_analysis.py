import argparse
import os.path

import pysam
from scipy.stats import binomtest


def parse_isoquant_read_assignment(assignment_file, assignment_type=["unique", "unique_minor_difference"]):
    records = {}
    with open(assignment_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if parts[5] in assignment_type:
                read_name, isoform_id, gene_id = parts[0], parts[3], parts[4]
                records.setdefault(gene_id, {}).setdefault(isoform_id, []).append(read_name)
    return records


def get_gene_regions(annotation_file):
    assert ".gff3" in annotation_file or ".gtf" in annotation_file, "Error: Unsupported annotation file format"
    gene_regions = {}

    def process_line(parts, gene_id):
        chr, start, end = parts[0], int(parts[3]), int(parts[4])
        gene_regions[gene_id] = {"chr": chr, "start": start, "end": end}

    def parse_attributes(attributes, key):
        return [x for x in attributes if x.startswith(key)][0].split("=")[1]

    if annotation_file.endswith(".gz"):
        import gzip
        with gzip.open(annotation_file, "rt") as f:
            if ".gff3" in annotation_file:
                for line in f:
                    if line.startswith("#"):
                        continue
                    parts = line.strip().split("\t")
                    if parts[2] == "gene":
                        attributes = parts[8].split(";")
                        gene_id = parse_attributes(attributes, "ID=")
                        process_line(parts, gene_id)
            elif ".gtf" in annotation_file:
                for line in f:
                    if line.startswith("#"):
                        continue
                    parts = line.strip().split("\t")
                    if parts[2] == "gene":
                        attributes = parts[8].split(";")
                        tag, gene_id = attributes[0].split(" ")
                        gene_id = gene_id.replace('"', '')
                        assert tag == "gene_id"
                        process_line(parts, gene_id)
    else:
        with open(annotation_file) as f:
            if ".gff3" in annotation_file:
                for line in f:
                    if line.startswith("#"):
                        continue
                    parts = line.strip().split("\t")
                    if parts[2] == "gene":
                        attributes = parts[8].split(";")
                        gene_id = parse_attributes(attributes, "ID=")
                        process_line(parts, gene_id)
            elif ".gtf" in annotation_file:
                for line in f:
                    if line.startswith("#"):
                        continue
                    parts = line.strip().split("\t")
                    if parts[2] == "gene":
                        attributes = parts[8].split(";")
                        tag, gene_id = attributes[0].split(" ")
                        gene_id = gene_id.replace('"', '')
                        assert tag == "gene_id"
                        process_line(parts, gene_id)
    return gene_regions


def get_reads_tag(bam_file, chr, start_pos, end_pos):
    reads_tag = {}
    with pysam.AlignmentFile(bam_file, "rb") as f:
        for read in f.fetch(chr, start_pos, end_pos):
            PS = read.get_tag("PS") if read.has_tag("PS") else None
            HP = read.get_tag("HP") if read.has_tag("HP") else None
            reads_tag[read.query_name] = {"PS": PS, "HP": HP}
    return reads_tag

def write_header(out_file):
    with open(out_file, "w") as f:
        f.write("Gene\tIsoform\tChromosome\tStart\tEnd\tPS\tHP1\tHP2\tP-value\n")

def analyze_haplotype_isoform(assignment_file, annotation_file, bam_file, out_file, min_support=10, p_value=0.05):
    gene_regions = get_gene_regions(annotation_file)
    isoquant_read_assignment = parse_isoquant_read_assignment(assignment_file)
    write_header(out_file)
    with open(out_file, "a") as fwriter:
        for gene_id, isoform_records in isoquant_read_assignment.items():
            gene_region = gene_regions[gene_id]
            reads_tag = get_reads_tag(bam_file, gene_region["chr"], gene_region["start"], gene_region["end"])
            for isoform_id, read_names in isoform_records.items():
                ps_dict = {}
                for rname in read_names:
                    if rname in reads_tag:
                        read_info = reads_tag[rname]
                        if read_info['PS'] and read_info['HP']:
                            ps_dict.setdefault(read_info['PS'], [0, 0])[read_info['HP'] - 1] += 1
                for ps, hp_counts in ps_dict.items():
                    total_counts = sum(hp_counts)
                    if total_counts >= min_support:
                        p = binomtest(hp_counts[0], total_counts, 0.5, alternative='two-sided').pvalue
                        if p <= p_value:
                            fwriter.write(f"{gene_id}\t{isoform_id}\t{gene_region['chr']}\t{gene_region['start']}\t{gene_region['end']}\t{ps}\t{hp_counts[0]}\t{hp_counts[1]}\t{p}\n")


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--assignment", required=True, help="Isoquant read assignment file")
    parser.add_argument("-a", "--annotation", required=True, help="Annotation file")
    parser.add_argument("-b", "--bam", required=True, help="BAM file")
    parser.add_argument("-o", "--output", required=True, help="Output tsv file to report ASE events")
    parser.add_argument("--min_support", type=int, default=10,
                        help="Minimum support reads for counting event (default: 10)")
    parser.add_argument("--p_value", type=float, default=0.05,
                        help="P-value threshold for Binomial test (default: 0.05)")

    args = parser.parse_args()
    if os.path.exists(args.output):
        raise FileExistsError(f"Output file {args.output} already exists")
    analyze_haplotype_isoform(args.assignment, args.annotation, args.bam, args.output, args.min_support, args.p_value)
