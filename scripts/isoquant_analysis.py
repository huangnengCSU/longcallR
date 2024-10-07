import argparse
import gzip
import pysam
from scipy.stats import binomtest
from collections import defaultdict
from multiprocessing import Manager
import concurrent.futures


def parse_isoquant_read_assignment(assignment_file, assignment_type, classification_type):
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


# def get_gene_regions(annotation_file):
#     assert ".gff3" in annotation_file or ".gtf" in annotation_file, "Error: Unsupported annotation file format"
#     gene_regions = {}
#
#     def process_line(parts, gene_id):
#         chr, start, end = parts[0], int(parts[3]), int(parts[4])
#         gene_regions[gene_id] = {"chr": chr, "start": start, "end": end}
#
#     def parse_attributes(attributes, key):
#         return [x for x in attributes if x.startswith(key)][0].split("=")[1]
#
#     if annotation_file.endswith(".gz"):
#         import gzip
#         with gzip.open(annotation_file, "rt") as f:
#             if ".gff3" in annotation_file:
#                 for line in f:
#                     if line.startswith("#"):
#                         continue
#                     parts = line.strip().split("\t")
#                     if parts[2] == "gene":
#                         attributes = parts[8].split(";")
#                         gene_id = parse_attributes(attributes, "ID=")
#                         process_line(parts, gene_id)
#             elif ".gtf" in annotation_file:
#                 for line in f:
#                     if line.startswith("#"):
#                         continue
#                     parts = line.strip().split("\t")
#                     if parts[2] == "gene":
#                         attributes = parts[8].split(";")
#                         tag, gene_id = attributes[0].split(" ")
#                         gene_id = gene_id.replace('"', '')
#                         assert tag == "gene_id"
#                         process_line(parts, gene_id)
#     else:
#         with open(annotation_file) as f:
#             if ".gff3" in annotation_file:
#                 for line in f:
#                     if line.startswith("#"):
#                         continue
#                     parts = line.strip().split("\t")
#                     if parts[2] == "gene":
#                         attributes = parts[8].split(";")
#                         gene_id = parse_attributes(attributes, "ID=")
#                         process_line(parts, gene_id)
#             elif ".gtf" in annotation_file:
#                 for line in f:
#                     if line.startswith("#"):
#                         continue
#                     parts = line.strip().split("\t")
#                     if parts[2] == "gene":
#                         attributes = parts[8].split(";")
#                         tag, gene_id = attributes[0].split(" ")
#                         gene_id = gene_id.replace('"', '')
#                         assert tag == "gene_id"
#                         process_line(parts, gene_id)
#     return gene_regions


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
        return {key_value.split("=")[0]: key_value.split("=")[1] for key_value in attributes.split(";")}

    def parse_attributes_gtf(attributes):
        attr_dict = {}
        for attr in attributes.strip().split(";"):
            if attr:
                key, value = attr.strip().split(" ")
                attr_dict[key] = value.replace('"', '')
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
                try:
                    gene_name = attr_dict["gene_name"]
                except KeyError:
                    gene_name = "."  # Use a placeholder if gene name is not available
                if gene_type in gene_types:
                    process_gene(parts, gene_id, gene_name)
            elif feature_type == "exon":
                gene_type = attr_dict["gene_type"]
                transcript_id = attr_dict["transcript_id"]
                gene_id = attr_dict["gene_id"]
                if gene_type in gene_types:
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


def get_reads_tag(bam_file, chr, start_pos, end_pos):
    reads_tag = {}
    with pysam.AlignmentFile(bam_file, "rb") as f:
        for read in f.fetch(chr, start_pos, end_pos):
            PS = read.get_tag("PS") if read.has_tag("PS") else None
            HP = read.get_tag("HP") if read.has_tag("HP") else None
            reads_tag[read.query_name] = {"PS": PS, "HP": HP}
    return reads_tag


def calculate_ase_pvalue(bam_file, gene_id, gene_name, gene_region, min_count, isoquant_read_assignments):
    reads_tag = get_reads_tag(bam_file, gene_region["chr"], gene_region["start"], gene_region["end"])
    assigned_reads = set()
    for isoform_id, reads in isoquant_read_assignments[gene_id].items():
        assigned_reads.update(reads)
    phase_set_hap_count = defaultdict(lambda: {1: 0, 2: 0})  # key: phase set, value: {haplotype: count}
    for rname in assigned_reads:
        ps = reads_tag[rname]["PS"]
        hp = reads_tag[rname]["HP"]
        if ps and hp:
            phase_set_hap_count[ps][hp] += 1
    p_value = 1.0
    h1_count, h2_count = 0, 0
    phase_set = "."  # Placeholder for phase set
    for ps, hap_count in phase_set_hap_count.items():
        if hap_count[1] + hap_count[2] < min_count:
            continue
        # Calculate Binomial test p-value
        total_reads = hap_count[1] + hap_count[2]
        p_value_ase = binomtest(hap_count[1], total_reads, 0.5, alternative='two-sided').pvalue
        if p_value_ase < p_value:
            p_value = p_value_ase
            h1_count = hap_count[1]
            h2_count = hap_count[2]
            phase_set = ps
    return (gene_name, p_value, phase_set, h1_count, h2_count)


def analyze_ase_genes(assignment_file, annotation_file, bam_file, out_file, threads, gene_types, assignment_type,
                      classification_type, min_support):
    isoquant_read_assignments = parse_isoquant_read_assignment(assignment_file, assignment_type, classification_type)
    gene_regions, gene_names, gene_strands, exon_regions, intron_regions = get_gene_regions(annotation_file, gene_types)

    gene_args = [(bam_file, gene_id, gene_names[gene_id], gene_regions[gene_id], min_support)
                 for gene_id in gene_regions.keys()
                 if gene_id in isoquant_read_assignments]
    results = []
    # Use a Manager to share isoquant_read_assignments across processes
    with Manager() as manager:
        shared_assignments = manager.dict(isoquant_read_assignments)
        with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
            # Submit all tasks at once without chunking
            futures = [executor.submit(calculate_ase_pvalue, *gene_data, shared_assignments) for gene_data in gene_args]
            for future in concurrent.futures.as_completed(futures):
                results.append(future.result())
    with open(out_file, "w") as f:
        f.write("Gene\tPS\tH1\tH2\tP-value\n")
        for gene_name, p_value, ps, h1, h2 in results:
            f.write(f"{gene_name}\t{ps}\t{h1}\t{h2}\t{p_value}\n")


def analyze_tpm(assignment_file, annotation_file, bam_file, gene_types, assignment_type, classification_type,
                isoquant_gene_tpm,
                isoquant_transcript_tpm, output_gene_tpm, output_transcript_tpm):
    isoquant_read_assignment = parse_isoquant_read_assignment(assignment_file, assignment_type, classification_type)
    gene_regions, gene_names, gene_strands, exon_regions, intron_regions = get_gene_regions(annotation_file, gene_types)
    gene_tpm = parse_isoquant_tpm(isoquant_gene_tpm)
    transcript_tpm = parse_isoquant_tpm(isoquant_transcript_tpm)
    gene_tpm_writer = open(output_gene_tpm, "w")
    gene_tpm_writer.write("Gene\tTPM\tHap1_TPM\tHap2_TPM\n")
    transcript_tpm_writer = open(output_transcript_tpm, "w")
    transcript_tpm_writer.write("Isoform\tTPM\tHap1_TPM\tHap2_TPM\tHap1_reads\tHap2_reads\tP-value\n")
    for gene_id, isoform_records in isoquant_read_assignment.items():
        if gene_id not in gene_regions:
            continue
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
                pvalue = 1.0
            else:
                h1_tpm = h1_cnt / (h1_cnt + h2_cnt) * tpm
                h2_tpm = h2_cnt / (h1_cnt + h2_cnt) * tpm
                pvalue = binomtest(h1_cnt, h1_cnt + h2_cnt, 0.5, alternative='two-sided').pvalue
            transcript_tpm_writer.write(f"{isoform_id}\t{tpm}\t{h1_tpm}\t{h2_tpm}\t{h1_cnt}\t{h2_cnt}\t{pvalue}\n")
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
    parser.add_argument("-b", "--bam", required=True, help="BAM file")
    parser.add_argument("-a", "--annotation", required=True, help="Annotation file")
    parser.add_argument("-i", "--assignment", required=True, help="Isoquant read assignment file")
    parser.add_argument("-g", "--gene_tpm", required=True, help="Isoquant gene TPM file")
    parser.add_argument("-t", "--transcript_tpm", required=True, help="Isoquant transcript TPM file")
    parser.add_argument("-o", "--output", required=True, help="prefix of output file")
    parser.add_argument("-p", "--processes", type=int, default=1, help="Number of process to run")
    parser.add_argument("--gene_types", type=str, nargs="+", default=["protein_coding", "lncRNA"],
                        help='Gene types to be analyzed. Default is ["protein_coding", "lncRNA"]', )
    parser.add_argument('--assignment_type', type=str, nargs='+', default=["unique", "unique_minor_difference"],
                        help='Assignment types to include. Default is ["unique", "unique_minor_difference"].')
    parser.add_argument('--classification_type', type=str, nargs='+',
                        default=["full_splice_match", "incomplete_splice_match", "mono_exon_match"],
                        help='Classification types to include. Default is ["full_splice_match", "incomplete_splice_match", "mono_exon_match"].')
    parser.add_argument("--min_support", type=int, default=10,
                        help="Minimum support reads for counting event (default: 10)")

    args = parser.parse_args()

    gene_types = set(args.gene_types)
    assignment_type = set(args.assignment_type)
    classification_type = set(args.classification_type)
    analyze_ase_genes(args.assignment, args.annotation, args.bam, args.output + ".ase.tsv", args.processes, gene_types,
                      assignment_type, classification_type, args.min_support)

    analyze_tpm(args.assignment, args.annotation, args.bam, gene_types, assignment_type, classification_type,
                args.gene_tpm, args.transcript_tpm, args.output + ".haplotype_gene_tpm.tsv",
                args.output + ".haplotype_transcript_tpm.tsv")
