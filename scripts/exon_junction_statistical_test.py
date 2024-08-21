import argparse

import gzip
from collections import defaultdict


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


if __name__ == "__main__":
    parse = argparse.ArgumentParser()
    parse.add_argument("-a", "--annotation_file", help="Annotation file in GFF3 or GTF format", required=True)
    gene_regions, exon_regions, intron_regions = get_gene_regions(parse.parse_args().annotation_file)
