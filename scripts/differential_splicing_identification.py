import argparse
import concurrent.futures
import gzip
import math
from collections import defaultdict
import networkx as nx

import numpy as np
import pysam
from scipy.stats import fisher_exact, power_divergence
from intervaltree import IntervalTree


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


def cluster_junctions_connected_components(reads_junctions, min_count=10):
    junctions_clusters = []
    junctions = {}  # key: (start, end), value: count

    # Count occurrences of each junction region
    for read_name, junction_regions in reads_junctions.items():
        for junction_region in junction_regions:
            junctions[junction_region] = junctions.get(junction_region, 0) + 1

    # Remove junctions with count less than min_count
    junctions = {k: v for k, v in junctions.items() if v >= min_count}

    # Create a graph
    G = nx.Graph()

    # Add nodes for each junction
    for junction in junctions.keys():
        G.add_node(junction)  # "junction = (intron_start, intron_end)"

    # Add edges between overlapping junctions (sharing start or end positions)
    junction_list = list(junctions.keys())
    for i in range(len(junction_list)):
        for j in range(i + 1, len(junction_list)):
            start1, end1 = junction_list[i]
            start2, end2 = junction_list[j]

            # Check if they share a start or end position
            if start1 == start2 or end1 == end2:
                G.add_edge(junction_list[i], junction_list[j])

    # Find connected components
    connected_components = list(nx.connected_components(G))

    # Collect the clusters of junctions
    for component in connected_components:
        junctions_clusters.append(list(component))

    return junctions_clusters, junctions


def cluster_junctions_exons_connected_components(reads_junctions, reads_exons, min_count=10):
    junctions_clusters = []
    junctions = {}  # key: (start, end), value: count

    # Count occurrences of each junction region
    for read_name, junction_regions in reads_junctions.items():
        for junction_region in junction_regions:
            junctions[junction_region] = junctions.get(junction_region, 0) + 1

    # Remove junctions with count less than min_count
    junctions = {k: v for k, v in junctions.items() if v >= min_count}

    # Count occurrences of each exon region
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

    # Remove exons with count less than min_count
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

    # Find connected components
    connected_components = list(nx.connected_components(G))

    # Collect the clusters of junctions
    for component in connected_components:
        clu = []
        for node in component:
            if node[2] == "junction":
                clu.append((node[0], node[1]))
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
        if read_start > start_pos or read_end < end_pos:
            continue
        # if read_start > end_pos or read_end < start_pos:
        #     continue
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
    def __init__(self, chr, start, end, novel, gene_name, strand, junction_set, phase_set, hap1_absent, hap1_present,
                 hap2_absent, hap2_present, p_value, sor):
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
        self.gene_names = [gene_name]

    @staticmethod
    def __header__():
        return ("#Junction\tStrand\tJunction_set\tPhase_set\tHap1_absent\tHap1_present\tHap2_absent\tHap2_present\t"
                "P_value\tSOR\tNovel\tGene_names")

    def __str__(self):
        gene_names = ",".join(set(self.gene_names))
        return (f"{self.chr}:{self.start}-{self.end}\t{self.strand}\t{self.junction_set}\t{self.phase_set}\t"
                f"{self.hap1_absent}\t{self.hap1_present}\t{self.hap2_absent}\t{self.hap2_present}\t"
                f"{self.p_value}\t{self.sor}\t{self.novel}\t{gene_names}")


def calc_sor(hap1_absent, hap1_present, hap2_absent, hap2_present):
    R = ((hap1_absent + 1) * (hap2_present + 1)) / ((hap1_present + 1) * (hap2_absent + 1))
    R_inverse = 1 / R
    sum = R + R_inverse
    SOR = math.log(sum)
    return SOR


def haplotype_event_test(absent_reads, present_reads, reads_tags):
    """
    Perform Fisher's exact test to determine if the haplotype distribution is significantly different between absent and present reads.
    :param absent_reads:
    :param present_reads:
    :param reads_tags:
    :return:
    """
    events = []
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
        if phase_set == ".":
            continue
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
        event = (phase_set, hap_absent_counts[phase_set][1], hap_present_counts[phase_set][1],
                 hap_absent_counts[phase_set][2], hap_present_counts[phase_set][2], pvalue, sor)
        events.append(event)
    return events


def analyze_gene(gene_name, gene_strand, annotation_exons, annotation_junctions, gene_region, bam_file, min_count,
                 p_value_threshold, sor_threshold):
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

    # Extract relevant reads and regions
    reads_positions, reads_exons, reads_introns, reads_tags = parse_reads_from_alignment(bam_file, chr, start, end)
    # junctions_clusters, read_junctions = cluster_junctions_connected_components(reads_introns, min_count)
    junctions_clusters, read_junctions = cluster_junctions_exons_connected_components(reads_introns, reads_exons,
                                                                                      min_count)

    # filter reads which have no overlapped exons with current gene exons
    intervalt = IntervalTree()
    for anno_exon in gene_exon_set:
        exon_start, exon_end = anno_exon[1:3]
        intervalt.addi(exon_start, exon_end + 1)  # interval is half-open, left inclusive, right exclusive

    reads_to_remove = []
    for qname, read_exons in reads_exons.items():
        overlapped = False
        for (exon_start, exon_end) in read_exons:
            if bool(intervalt.overlap(exon_start, exon_end + 1)):
                overlapped = True
                break
        if not overlapped:
            reads_to_remove.append(qname)

    # Remove reads after collecting them
    for qname in reads_to_remove:
        del reads_positions[qname]
        del reads_exons[qname]
        del reads_introns[qname]
        del reads_tags[qname]

    gene_ase_events = []

    # Analyze junction regions
    for junc_cluster in junctions_clusters:
        if len(junc_cluster) == 0:
            continue
        junction_set = f"{chr}:{junc_cluster[0][0]}-{junc_cluster[0][1]}"
        for read_junc in junc_cluster:
            junction_start = read_junc[0]
            junction_end = read_junc[1]
            novel = (chr, junction_start, junction_end) not in gene_junction_set
            absences, presents = check_absent_present(junction_start, junction_end, reads_positions, reads_introns)
            test_results = haplotype_event_test(absences, presents, reads_tags)
            for event in test_results:
                (phase_set, h1_a, h1_p, h2_a, h2_p, pvalue, sor) = event
                if pvalue == 0:
                    gene_ase_events.append(AseEvent(chr, junction_start, junction_end, novel, gene_name, gene_strand,
                                                    junction_set, phase_set, h1_a, h1_p, h2_a, h2_p, pvalue, sor))
                elif pvalue < p_value_threshold and sor >= sor_threshold:
                    gene_ase_events.append(AseEvent(chr, junction_start, junction_end, novel, gene_name, gene_strand,
                                                    junction_set, phase_set, h1_a, h1_p, h2_a, h2_p, pvalue, sor))
    return gene_ase_events


def analyze(annotation_file, bam_file, output_file, min_count, threads, p_value_threshold, sor_threshold):
    all_ase_events = {}  # key: (chr, start, end), value: AseEvent
    anno_gene_regions, anno_gene_names, anno_gene_strands, anno_exon_regions, anno_intron_regions = get_gene_regions(
        annotation_file)

    # Prepare data for multiprocessing
    gene_data_list = [(anno_gene_names[gene_id], anno_gene_strands[gene_id], anno_exon_regions[gene_id],
                       anno_intron_regions[gene_id], gene_region, bam_file, min_count, p_value_threshold, sor_threshold)
                      for gene_id, gene_region in anno_gene_regions.items()]

    # Use ProcessPoolExecutor for multiprocessing
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        # Submit tasks and collect futures
        futures = [executor.submit(analyze_gene, *gene_data) for gene_data in gene_data_list]

        # As each future completes, collect the results
        for future in concurrent.futures.as_completed(futures):
            # all_ase_events.extend(future.result().values())  # Waits for the future to complete and gets the result
            for event in future.result():
                if (event.chr, event.start, event.end) in all_ase_events.keys():
                    tmp_event = all_ase_events[(event.chr, event.start, event.end)]
                    if event.p_value < tmp_event.p_value:
                        tmp_event.p_value = event.p_value
                        tmp_event.sor = event.sor
                        tmp_event.hap1_absent = event.hap1_absent
                        tmp_event.hap1_present = event.hap1_present
                        tmp_event.hap2_absent = event.hap2_absent
                        tmp_event.hap2_present = event.hap2_present
                        tmp_event.phase_set = event.phase_set
                    tmp_event.gene_names.extend(event.gene_names)
                    if not event.novel:
                        tmp_event.novel = event.novel
                else:
                    all_ase_events[(event.chr, event.start, event.end)] = event

    # write to output file
    with open(output_file, "w") as f:
        f.write(AseEvent.__header__() + "\n")
        for key, event in all_ase_events.items():
            f.write(event.__str__() + "\n")


if __name__ == "__main__":
    parse = argparse.ArgumentParser()
    parse.add_argument("-a", "--annotation_file", help="Annotation file in GFF3 or GTF format", required=True)
    parse.add_argument("-b", "--bam_file", help="BAM file", required=True)
    parse.add_argument("-o", "--output_file", help="Output file", required=True)
    parse.add_argument("-t", "--threads", help="Number of threads", default=1, type=int)
    parse.add_argument("-m", "--min_sup", help="Minimum support of phased reads for exon or junction", default=10,
                       type=int)
    parse.add_argument("-p", "--p_value_threshold", help="P-value threshold for Fisher's exact test", default=1e-15,
                       type=float)
    parse.add_argument("-s", "--sor_threshold", help="SOR threshold for filtering, higher value is more stringent",
                       default=2.0, type=float)
    args = parse.parse_args()
    analyze(args.annotation_file, args.bam_file, args.output_file, args.min_sup, args.threads, args.p_value_threshold,
            args.sor_threshold)
