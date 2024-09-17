import argparse
import concurrent.futures
import gzip
import math
from collections import defaultdict

import numpy as np
import pysam
from scipy.stats import fisher_exact, power_divergence


def get_gene_regions(annotation_file):
    """Parse gene, exon, and intron regions from a GFF3 or GTF file.
    :param annotation_file: Path to the annotation file
    :return: Gene regions, exon regions, and intron regions
    """
    assert annotation_file.endswith(
        (".gff3", ".gtf", ".gff3.gz", ".gtf.gz")), "Error: Unsupported annotation file format"

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
                if gene_type != "processed_pseudogene":
                    process_gene(parts, gene_id, gene_name)
            elif feature_type == "exon":
                gene_type = attr_dict["gene_type"]
                transcript_id = attr_dict["transcript_id"]
                gene_id = attr_dict["gene_id"]
                if gene_type != "processed_pseudogene":
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
    contigs = samfile.references
    if chr not in contigs:
        return reads_positions, reads_exons, reads_junctions, reads_tags
    for read in samfile.fetch(chr, fetch_start, fetch_end):
        # read is not phased, ignore
        # if not read.has_tag("PS") or not read.has_tag("HP"):
        #     continue
        if not read.has_tag("HP"):
            continue
        if read.has_tag("PS"):
            PS_tag = read.get_tag("PS")
        else:
            PS_tag = "."
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
        # if read_start > end_pos or read_end < start_pos:
        #     continue
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


class AseEvent:
    def __init__(self, chr, start, end, exon_or_junction, novel, gene_id, gene_name, gene_strand, phase_set,
                 hap1_absent,
                 hap1_present, hap2_absent, hap2_present, p_value, sor):
        self.chr = chr
        self.start = start  # 1-based, inclusive
        self.end = end  # 1-based, inclusive
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.gene_strand = gene_strand
        self.novel = novel
        self.exon_or_junction = exon_or_junction
        self.phase_set = phase_set
        self.hap1_absent = hap1_absent
        self.hap1_present = hap1_present
        self.hap2_absent = hap2_absent
        self.hap2_present = hap2_present
        self.p_value = p_value
        self.sor = sor

    @staticmethod
    def __header__():
        return ("#gene_id\tgene_name\tstrand\tregion\tExon/Junction\tnovel\tphase_set\t"
                "hap1_absent\thap1_present\thap2_absent\thap2_present\tp_value\tSOR")

    def __str__(self):
        return (
            f"{self.gene_id}\t{self.gene_name}\t{self.gene_strand}\t{self.chr}:{self.start}-{self.end}\t{self.exon_or_junction}\t"
            f"{self.novel}\t{self.phase_set}\t{self.hap1_absent}\t{self.hap1_present}\t{self.hap2_absent}\t"
            f"{self.hap2_present}\t{self.p_value}\t{self.sor}")


def calc_sor(hap1_absent, hap1_present, hap2_absent, hap2_present):
    R = ((hap1_absent + 1) * (hap2_present + 1)) / ((hap1_present + 1) * (hap2_absent + 1))
    R_inverse = 1 / R
    sum = R + R_inverse
    SOR = math.log(sum)
    return SOR


def haplotype_event_test(gene_id, gene_name, gene_strand, chr, start, end, exon_or_junction, novel,
                         absent_reads, present_reads, reads_tags, p_value_threshold):
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
    ase_events = []
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
        # if phase_set == ".":
        #     continue
        ## if not exist, set 0
        # take phased reads without phase set into account
        hap_absent_counts[phase_set][1] = hap_absent_counts[phase_set].get(1, 0) + hap_absent_counts["."].get(1, 0)
        hap_absent_counts[phase_set][2] = hap_absent_counts[phase_set].get(2, 0) + hap_absent_counts["."].get(2, 0)
        hap_present_counts[phase_set][1] = hap_present_counts[phase_set].get(1, 0) + hap_present_counts["."].get(1, 0)
        hap_present_counts[phase_set][2] = hap_present_counts[phase_set].get(2, 0) + hap_present_counts["."].get(2, 0)
        table = np.array([[hap_absent_counts[phase_set][1], hap_absent_counts[phase_set][2]],
                          [hap_present_counts[phase_set][1], hap_present_counts[phase_set][2]]])
        ## Fisher's exact test
        oddsratio, pvalue_fisher = fisher_exact(table)
        ## G-test
        g_stat, pvalue_gtest = power_divergence(f_obs=table + 1e-30, lambda_="log-likelihood")
        pvalue_gtest = np.min(pvalue_gtest)

        ## Use the maximum p-value from Fisher's exact test and G-test
        pvalue = max(pvalue_fisher, pvalue_gtest)

        ## Calculate SOR, refer to GATK AS_StrandOddsRatio, https://gatk.broadinstitute.org/hc/en-us/articles/360037224532-AS-StrandOddsRatio
        sor = calc_sor(hap_absent_counts[phase_set][1], hap_present_counts[phase_set][1],
                       hap_absent_counts[phase_set][2], hap_present_counts[phase_set][2])

        if pvalue < p_value_threshold:
            if exon_or_junction == 0:
                event = AseEvent(chr, start, end, 'Exon', novel, gene_id, gene_name, gene_strand, phase_set,
                                 hap_absent_counts[phase_set][1], hap_present_counts[phase_set][1],
                                 hap_absent_counts[phase_set][2], hap_present_counts[phase_set][2], pvalue, sor)
                ase_events.append(event)
            elif exon_or_junction == 1:
                event = AseEvent(chr, start, end, 'Junc', novel, gene_id, gene_name, gene_strand, phase_set,
                                 hap_absent_counts[phase_set][1], hap_present_counts[phase_set][1],
                                 hap_absent_counts[phase_set][2], hap_present_counts[phase_set][2], pvalue, sor)
                ase_events.append(event)
    return ase_events


def analyze_gene(gene_id, gene_name, gene_strand, annotation_exons, annotation_junctions, region, bam_file, min_count,
                 p_value_threshold):
    chr = region["chr"]
    start = region["start"]
    end = region["end"]

    # get all exons and junctions
    gene_exon_set = set()
    for transcript_id, anno_exons in annotation_exons.items():
        for exon in anno_exons:
            gene_exon_set.add(exon)
    gene_junction_set = set()
    for transcript_id, anno_junctions in annotation_junctions.items():
        for junction in anno_junctions:
            gene_junction_set.add(junction)

    # Extract relevant reads and regions
    reads_positions, reads_exons, reads_introns, reads_tags = parse_reads_from_alignment(bam_file, chr, start, end)
    read_clustered_exons = cluster_exons(reads_exons, min_count)
    read_clustered_junctions = cluster_junctions(reads_introns, min_count)

    gene_ase_events = []
    # Analyze exon regions
    for exon in read_clustered_exons:
        exon_start = exon[0]
        exon_end = exon[1]
        novel_exon = (chr, exon_start, exon_end) not in gene_exon_set
        absences, presents = check_absent_present(exon_start, exon_end, 0, reads_positions, reads_exons,
                                                  reads_introns)
        ase_events = haplotype_event_test(gene_id, gene_name, gene_strand, chr, exon_start, exon_end, 0, novel_exon,
                                          absences, presents, reads_tags, p_value_threshold)
        gene_ase_events.extend(ase_events)

    # Analyze junction regions
    for junction in read_clustered_junctions:
        junction_start = junction[0]
        junction_end = junction[1]
        novel_junction = (chr, junction_start, junction_end) not in gene_junction_set
        absences, presents = check_absent_present(junction_start, junction_end, 1, reads_positions,
                                                  reads_exons, reads_introns)
        ase_events = haplotype_event_test(gene_id, gene_name, gene_strand, chr, junction_start, junction_end, 1,
                                          novel_junction, absences, presents, reads_tags, p_value_threshold)
        gene_ase_events.extend(ase_events)

    return gene_ase_events


def analyze(annotation_file, bam_file, output_file, min_count, threads, p_value_threshold):
    all_ase_events = []
    anno_gene_regions, anno_gene_names, anno_gene_strands, anno_exon_regions, anno_intron_regions = get_gene_regions(
        annotation_file)

    # Prepare data for multiprocessing
    gene_data_list = [(gene_id, anno_gene_names[gene_id], anno_gene_strands[gene_id], anno_exon_regions[gene_id],
                       anno_intron_regions[gene_id], region, bam_file, min_count, p_value_threshold) for gene_id, region
                      in anno_gene_regions.items()]

    # Use ProcessPoolExecutor for multiprocessing
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        # Submit tasks and collect futures
        futures = [executor.submit(analyze_gene, *gene_data) for gene_data in gene_data_list]

        # As each future completes, collect the results
        for future in concurrent.futures.as_completed(futures):
            all_ase_events.extend(future.result())  # Waits for the future to complete and gets the result

    # write to output file
    with open(output_file, "w") as f:
        f.write(AseEvent.__header__() + "\n")
        for event in all_ase_events:
            f.write(event.__str__() + "\n")


if __name__ == "__main__":
    parse = argparse.ArgumentParser()
    parse.add_argument("-a", "--annotation_file", help="Annotation file in GFF3 or GTF format", required=True)
    parse.add_argument("-b", "--bam_file", help="BAM file", required=True)
    parse.add_argument("-o", "--output_file", help="Output file", required=True)
    parse.add_argument("-t", "--threads", help="Number of threads", default=1, type=int)
    parse.add_argument("-s", "--min_sup", help="Minimum support of phased reads for exon or junction", default=10)
    parse.add_argument("-p", "--p_value_threshold", help="P-value threshold for Fisher's exact test", default=1e-10)
    args = parse.parse_args()
    analyze(args.annotation_file, args.bam_file, args.output_file, args.min_sup, args.threads, args.p_value_threshold)
