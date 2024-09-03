import argparse

import pysam
from scipy.stats import binomtest
from collections import defaultdict


def parse_isoquant_read_assignment(assignment_file, assignment_type={"unique", "unique_minor_difference"},
                                   classification_type={"full_splice_match", "incomplete_splice_match",
                                                        "mono_exon_match"}):
    records = {}
    with open(assignment_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if parts[5] in assignment_type:
                read_name, isoform_id, gene_id, additional_info = parts[0], parts[3], parts[4], parts[8]
                additional_info_dict = {key.strip(): value.strip() for key, value in
                                        (pair.split('=') for pair in additional_info.split(';') if pair.strip())}
                try:
                    classification = additional_info_dict["Classification"]
                    if classification in classification_type:
                        records.setdefault(gene_id, {}).setdefault(isoform_id, []).append(read_name)
                except KeyError:
                    continue
    return records


def parse_isoquant_tpm(tpm_file):
    tpm_dict = {}
    with open(tpm_file, 'r') as f:
        next(f)
        for line in f:
            parts = line.strip().split("\t")
            feature_id, tpm = parts[0], float(parts[1])
            tpm_dict[feature_id] = tpm
    return tpm_dict


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


def analyze_haplotype_isoform(assignment_file, annotation_file, bam_file, out_file, assignment_type,
                              classification_type, min_support=10, p_value=0.05):
    isoquant_read_assignment = parse_isoquant_read_assignment(assignment_file, assignment_type, classification_type)
    gene_regions = get_gene_regions(annotation_file)
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
                            fwriter.write(
                                f"{gene_id}\t{isoform_id}\t{gene_region['chr']}\t{gene_region['start']}\t{gene_region['end']}\t{ps}\t{hp_counts[0]}\t{hp_counts[1]}\t{p}\n")


def analyze_tpm(assignment_file, annotation_file, bam_file, assignment_type, classification_type, isoquant_gene_tpm,
                isoquant_transcript_tpm, output_gene_tpm, output_transcript_tpm):
    isoquant_read_assignment = parse_isoquant_read_assignment(assignment_file, assignment_type, classification_type)
    gene_regions = get_gene_regions(annotation_file)
    gene_tpm = parse_isoquant_tpm(isoquant_gene_tpm)
    transcript_tpm = parse_isoquant_tpm(isoquant_transcript_tpm)
    gene_tpm_writer = open(output_gene_tpm, "w")
    gene_tpm_writer.write("Gene\tTPM\tHap1_TPM\tHap2_TPM\n")
    transcript_tpm_writer = open(output_transcript_tpm, "w")
    transcript_tpm_writer.write("Isoform\tTPM\tHap1_TPM\tHap2_TPM\n")
    for gene_id, isoform_records in isoquant_read_assignment.items():
        gene_region = gene_regions[gene_id]
        reads_tag = get_reads_tag(bam_file, gene_region["chr"], gene_region["start"], gene_region["end"])
        gene_dict = defaultdict(int)
        for isoform_id, read_names in isoform_records.items():
            ps_dict = defaultdict(int)
            for rname in read_names:
                if rname in reads_tag:
                    read_info = reads_tag[rname]
                    if read_info['HP']:
                        ps_dict[read_info['HP']] += 1
                        gene_dict[read_info['HP']] += 1
            h1_cnt = ps_dict[1]
            h2_cnt = ps_dict[2]
            tpm = transcript_tpm.get(isoform_id, 0)
            if h1_cnt + h2_cnt == 0:
                h1_tpm, h2_tpm = 0, 0
            else:
                h1_tpm = h1_cnt / (h1_cnt + h2_cnt) * tpm
                h2_tpm = h2_cnt / (h1_cnt + h2_cnt) * tpm
            transcript_tpm_writer.write(f"{isoform_id}\t{tpm}\t{h1_tpm}\t{h2_tpm}\n")
        h1_cnt = gene_dict[1]
        h2_cnt = gene_dict[2]
        tpm = gene_tpm.get(gene_id, 0)
        if h1_cnt + h2_cnt == 0:
            h1_tpm, h2_tpm = 0, 0
        else:
            h1_tpm = h1_cnt / (h1_cnt + h2_cnt) * tpm
            h2_tpm = h2_cnt / (h1_cnt + h2_cnt) * tpm
        gene_tpm_writer.write(f"{gene_id}\t{tpm}\t{h1_tpm}\t{h2_tpm}\n")
    gene_tpm_writer.close()
    transcript_tpm_writer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--assignment", required=True, help="Isoquant read assignment file")
    parser.add_argument("-g", "--gene_tpm", required=True, help="Isoquant gene TPM file")
    parser.add_argument("-t", "--transcript_tpm", required=True, help="Isoquant transcript TPM file")
    parser.add_argument("-a", "--annotation", required=True, help="Annotation file")
    parser.add_argument("-b", "--bam", required=True, help="BAM file")
    parser.add_argument("-o", "--output", required=True, help="prefix of output file")
    parser.add_argument('--assignment_type', type=str, nargs='+', default=["unique", "unique_minor_difference"],
                        help='Assignment types to include. Default is ["unique", "unique_minor_difference"].')
    parser.add_argument('--classification_type', type=str, nargs='+',
                        default=["full_splice_match", "incomplete_splice_match", "mono_exon_match"],
                        help='Classification types to include. Default is ["full_splice_match", "incomplete_splice_match", "mono_exon_match"].')
    parser.add_argument("--min_support", type=int, default=10,
                        help="Minimum support reads for counting event (default: 10)")
    parser.add_argument("--p_value", type=float, default=0.05,
                        help="P-value threshold for Binomial test (default: 0.05)")

    args = parser.parse_args()

    assignment_type = set(args.assignment_type)
    classification_type = set(args.classification_type)
    analyze_haplotype_isoform(args.assignment, args.annotation, args.bam, args.output + ".ase.tsv", assignment_type,
                              classification_type, args.min_support, args.p_value)
    analyze_tpm(args.assignment, args.annotation, args.bam, assignment_type, classification_type,
                args.gene_tpm, args.transcript_tpm,
                args.output + ".haplotype_gene_tpm.tsv",
                args.output + ".haplotype_transcript_tpm.tsv")
