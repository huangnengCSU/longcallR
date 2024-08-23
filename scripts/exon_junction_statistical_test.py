import argparse
import gzip
from collections import defaultdict

import pysam
from scipy.stats import fisher_exact


def get_gene_regions(annotation_file):
    """Parse gene, exon, and intron regions from a GFF3 or GTF file.
    :param annotation_file: Path to the annotation file
    :return: Gene regions, exon regions, and intron regions
    """
    assert annotation_file.endswith(
        (".gff3", ".gtf", ".gff3.gz", ".gtf.gz")), "Error: Unsupported annotation file format"

    gene_regions = {}
    exon_regions = defaultdict(lambda: defaultdict(list))
    intron_regions = defaultdict(lambda: defaultdict(list))

    def process_gene(parts, gene_id):
        chr, start, end = parts[0], int(parts[3]), int(parts[4])
        gene_regions[gene_id] = {"chr": chr, "start": start, "end": end}  # 1-based, start-inclusive, end-inclusive

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
                process_gene(parts, gene_id)
            elif feature_type == "exon":
                transcript_id = attr_dict["transcript_id"]
                gene_id = attr_dict["gene_id"]
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

    return gene_regions, exon_regions, intron_regions


def parse_reads_from_alignment(bam_file, chr, start_pos, end_pos):
    """Parse reads from a BAM file that overlap a specific region.
    :param bam_file: Path to the BAM file
    :param chr: Chromosome
    :param start_pos: Start position, 1-based inclusive
    :param end_pos: End position, 1-based inclusive
    :return:
    """

    def get_exon_intron_regions(read):
        exon_regions = []  # 1-based, start-inclusive, end-inclusive
        intron_regions = []  # 1-based, start-inclusive, end-inclusive
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
                current_position += length
            elif operation == 3:  # 'N' operation represents skipped region (intron)
                intron_start = current_position
                intron_end = current_position + length - 1
                intron_regions.append((intron_start, intron_end))
                current_position += length
            elif operation in {1, 4, 5}:  # 'I', 'S', 'H' operations
                pass
            else:
                current_position += length
        return exon_regions, intron_regions

    reads_exons = {}
    reads_junctions = {}
    reads_positions = {}  # 1-based, start-inclusive, end-inclusive
    reads_tags = {}  # key: read name, value: {"PS": phase set, "HP": haplotype}
    samfile = pysam.AlignmentFile(bam_file, "rb")
    fetch_start = start_pos - 1  # 0-based, start-inclusive
    fetch_end = end_pos  # 0-based, end-exclusive
    for read in samfile.fetch(chr, fetch_start, fetch_end):
        # read is not phased, ignore
        if not read.has_tag("PS") or not read.get_tag("HP"):
            continue
        PS_tag = read.get_tag("PS")
        HP_tag = read.get_tag("HP")
        reads_tags[read.query_name] = {"PS": PS_tag, "HP": HP_tag}
        # get read start position and end position, 1-based, start-inclusive, end-inclusive
        read_start = read.reference_start + 1
        read_end = read.reference_end
        reads_positions[read.query_name] = (read_start, read_end)
        # get all exons and introns
        exon_regions, intron_regions = get_exon_intron_regions(read)
        reads_exons[read.query_name] = exon_regions
        reads_junctions[read.query_name] = intron_regions
    return reads_positions, reads_exons, reads_junctions, reads_tags


def cluster_exons(reads_exons, min_count=10):
    exons = {}  # key: (start, end), value: count
    for read_name, exon_regions in reads_exons.items():
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
    # delete exons with count less than 10
    exons = {k: v for k, v in exons.items() if v >= min_count}
    return exons


def cluster_junctions(reads_junctions, min_count=10):
    junctions = {}  # key: (start, end), value: count
    for read_name, junction_regions in reads_junctions.items():
        for junction_region in junction_regions:
            junctions[junction_region] = junctions.get(junction_region, 0) + 1
    # delete junctions with count less than 10
    junctions = {k: v for k, v in junctions.items() if v >= min_count}
    return junctions


def check_absent_present(start_pos, end_pos, exon_or_junction, reads_positions, reads_exons, reads_junctions):
    """
    Find the reads where an exon or junction is absent or present.
    :param chr:
    :param start_pos: 1-based, start-inclusive
    :param end_pos: 1-based, end-inclusive
    :param exon_or_junction: 0 for exon, 1 for junction
    :param reads_positions:
    :param reads_exons:
    :param reads_junctions:
    :return:
    """
    absent_reads = []
    present_reads = []
    for read_name, (read_start, read_end) in reads_positions.items():
        if read_start > start_pos or read_end < end_pos:
            continue
        if exon_or_junction == 0:
            present = False
            for exon_start, exon_end in reads_exons[read_name]:
                if exon_start == start_pos and exon_end == end_pos:
                    present_reads.append(read_name)
                    present = True
                    break
            if not present:
                absent_reads.append(read_name)
        elif exon_or_junction == 1:
            present = False
            for junction_start, junction_end in reads_junctions[read_name]:
                if junction_start == start_pos and junction_end == end_pos:
                    present_reads.append(read_name)
                    present = True
                    break
            if not present:
                absent_reads.append(read_name)
    return absent_reads, present_reads


def haplotype_event_fisher_test(gene_id, chr, start, end, exon_or_junction, absent_reads, present_reads, reads_tags,
                                p_value_threshold=0.05):
    """
    Perform Fisher's exact test to determine if the haplotype distribution is significantly different between absent and present reads.
    :param gene_id:
    :param chr:
    :param start: 1-based, start-inclusive
    :param end: 1-based, end-inclusive
    :param exon_or_junction: 0 for exon, 1 for junction
    :param absent_reads:
    :param present_reads:
    :param reads_tags:
    :param p_value_threshold:
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
    for phase_set in all_phase_sets:
        ## if not exist, set 0
        hap_absent_counts[phase_set][1] = hap_absent_counts[phase_set].get(1, 0)
        hap_absent_counts[phase_set][2] = hap_absent_counts[phase_set].get(2, 0)
        hap_present_counts[phase_set][1] = hap_present_counts[phase_set].get(1, 0)
        hap_present_counts[phase_set][2] = hap_present_counts[phase_set].get(2, 0)
        ## Fisher's exact test
        table = [[hap_absent_counts[phase_set][1], hap_absent_counts[phase_set][2]],
                 [hap_present_counts[phase_set][1], hap_present_counts[phase_set][2]]]
        oddsratio, pvalue = fisher_exact(table)
        if pvalue < p_value_threshold:
            if exon_or_junction == 0:
                print(
                    f"event:{gene_id}, exon, {chr}:{start}-{end}, PS: {phase_set}, H1: {hap_absent_counts[phase_set][1]} absent, "
                    f"{hap_present_counts[phase_set][1]} present; H2: {hap_absent_counts[phase_set][2]} absent, "
                    f"{hap_present_counts[phase_set][2]} present; p-value: {pvalue}")
            elif exon_or_junction == 1:
                print(
                    f"event:{gene_id}, junction, {chr}:{start}-{end}, PS: {phase_set}, H1: {hap_absent_counts[phase_set][1]} absent, "
                    f"{hap_present_counts[phase_set][1]} present; H2: {hap_absent_counts[phase_set][2]} absent, "
                    f"{hap_present_counts[phase_set][2]} present; p-value: {pvalue}")


def analyze(annotation_file, bam_file, min_count, p_value_threshold):
    anno_gene_regions, anno_exon_regions, anno_intron_regions = get_gene_regions(annotation_file)
    for gene_id, region in anno_gene_regions.items():
        chr = region["chr"]
        start = region["start"]
        end = region["end"]
        reads_positions, reads_exons, reads_introns, reads_tags = parse_reads_from_alignment(bam_file, chr, start, end)
        exons = cluster_exons(reads_exons, min_count)
        junctions = cluster_junctions(reads_introns, min_count)
        for exon in exons:
            exon_start = exon[0]
            exon_end = exon[1]
            absences, presents = check_absent_present(exon_start, exon_end, 0, reads_positions, reads_exons,
                                                      reads_introns)
            haplotype_event_fisher_test(gene_id, chr, exon_start, exon_end, 0, absences, presents, reads_tags,
                                        p_value_threshold)
        for junction in junctions:
            junction_start = junction[0]
            junction_end = junction[1]
            absences, presents = check_absent_present(junction_start, junction_end, 1, reads_positions, reads_exons,
                                                      reads_introns)
            haplotype_event_fisher_test(gene_id, chr, junction_start, junction_end, 1, absences, presents, reads_tags,
                                        p_value_threshold)


if __name__ == "__main__":
    parse = argparse.ArgumentParser()
    parse.add_argument("-a", "--annotation_file", help="Annotation file in GFF3 or GTF format", required=True)
    parse.add_argument("-b", "--bam_file", help="BAM file", required=True)
    parse.add_argument("-s", "--min_sup", help="Minimum support of phased reads for exon or junction", default=10)
    parse.add_argument("-p", "--p_value_threshold", help="P-value threshold for Fisher's exact test", default=0.05)
    args = parse.parse_args()
    analyze(args.annotation_file, args.bam_file, args.min_sup, args.p_value_threshold)
