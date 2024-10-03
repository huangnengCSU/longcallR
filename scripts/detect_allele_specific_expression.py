import argparse
from collections import defaultdict
import gzip
import pysam
from scipy.stats import binomtest
import concurrent.futures


def get_gene_regions(annotation_file, gene_types):
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


def parse_reads_from_alignment(bam_file, chr, start_pos, end_pos):
    """Parse reads from a BAM file that overlap a specific region.
    :param bam_file: Path to the BAM file
    :param chr: Chromosome
    :param start_pos: Start position, 1-based inclusive
    :param end_pos: End position, 1-based inclusive
    :return:
    """

    def get_overlap_length(read, start_pos, end_pos):
        """Calculate the overlap length between a read and a region. Only consider Match and MisMatch.
            :param start_pos: Start position, 1-based inclusive
            :param end_pos: End position, 1-based inclusive
        """
        reference_start = read.reference_start + 1  # 1-based
        current_position = reference_start
        overlap_length = 0
        for cigartuple in read.cigartuples:
            operation, length = cigartuple
            if current_position > end_pos:
                break
            if operation in {0, 7, 8}:
                for i in range(length):
                    if current_position + i >= start_pos and current_position + i <= end_pos:
                        overlap_length += 1
                    if current_position + i > end_pos:
                        break
                current_position += length
            elif operation in {2, 3}:  # 'D', 'N' operations
                current_position += length
            elif operation in {1, 4, 5}:  # 'I', 'S', 'H' operations
                pass
        return overlap_length

    reads_positions = {}  # 1-based, start-inclusive, end-inclusive
    reads_tags = {}  # key: read name, value: {"PS": phase set, "HP": haplotype}
    samfile = pysam.AlignmentFile(bam_file, "rb")
    fetch_start = start_pos - 1  # 0-based, start-inclusive
    fetch_end = end_pos  # 0-based, end-exclusive
    contigs = samfile.references
    if chr not in contigs:
        return reads_positions, reads_tags
    for read in samfile.fetch(chr, fetch_start, fetch_end):
        # read is not phased, ignore
        # if not read.has_tag("PS") or not read.has_tag("HP"):
        #     continue
        if read.is_secondary or read.is_supplementary or read.is_unmapped:
            continue
        # overlap length between read and gene region should over 50% of the read length
        read_length = read.query_length
        overlap_length = get_overlap_length(read, start_pos, end_pos)
        if overlap_length < read_length * 0.5:
            continue
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
    return reads_positions, reads_tags


def calculate_ase_for_gene(gene_name, gene_region, bam_file, min_count):
    chr = gene_region["chr"]
    start = gene_region["start"]
    end = gene_region["end"]
    reads_positions, reads_tags = parse_reads_from_alignment(bam_file, chr, start, end)
    phase_set_hap_count = defaultdict(lambda: {1: 0, 2: 0})  # key: phase set, value: {haplotype: count}
    for tag in reads_tags.values():
        ps = tag["PS"]
        hp = tag["HP"]
        if ps == ".":
            continue
        phase_set_hap_count[ps][hp] += 1
    p_value = 1.0
    h1_count, h2_count = 0, 0
    phase_set = "."
    for ps in phase_set_hap_count.keys():
        if phase_set_hap_count[ps][1] + phase_set_hap_count[ps][2] < min_count:
            continue
        # Calculate Binomial test p-value
        total_reads = phase_set_hap_count[ps][1] + phase_set_hap_count[ps][2]
        p_value_ase = binomtest(phase_set_hap_count[ps][1], total_reads, 0.5, alternative='two-sided').pvalue
        if p_value_ase < p_value:
            p_value = p_value_ase
            h1_count = phase_set_hap_count[ps][1]
            h2_count = phase_set_hap_count[ps][2]
            phase_set = ps
    return (gene_name, p_value, phase_set, h1_count, h2_count)


def detect_allele_specific_expression(annotation_file, bam_file, output_file, min_count, gene_types, threads):
    gene_regions, gene_names, gene_strands, exon_regions, intron_regions = get_gene_regions(annotation_file, gene_types)
    gene_args = [(gene_names[gene_id], region, bam_file, min_count) for gene_id, region in gene_regions.items()]
    results = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(calculate_ase_for_gene, *gene_data) for gene_data in gene_args]
        for future in concurrent.futures.as_completed(futures):
            (gene_name, p_value, ps, h1, h2) = future.result()
            results.append((gene_name, p_value, ps, h1, h2))
    with open(output_file, "w") as f:
        f.write("Gene\tASE P-value\tPhase set\tHaplotype 1 count\tHaplotype 2 count\n")
        for (gene_name, p_value, ps, h1, h2) in results:
            f.write(f"{gene_name}\t{p_value}\t{ps}\t{h1}\t{h2}\n")


if __name__ == "__main__":
    parse = argparse.ArgumentParser(description="Detect allele-specific expression (ASE) from phased RNA-seq data.")
    parse.add_argument("-a", "--annotation", help="GFF3 or GTF annotation file", required=True)
    parse.add_argument("-b", "--bam", help="Phased RNA-seq BAM file", required=True)
    parse.add_argument("-o", "--output", help="Output file", required=True)
    parse.add_argument("-m", "--min_count", help="Minimum read count for ASE", type=int, default=10)
    parse.add_argument("--gene_types", help="Gene types to consider", nargs="+", default=["protein_coding", "lncRNA"])
    parse.add_argument("-t", "--threads", help="Number of threads", type=int, default=1)
    args = parse.parse_args()
    detect_allele_specific_expression(args.annotation, args.bam, args.output, args.min_count, args.gene_types,
                                      args.threads)
