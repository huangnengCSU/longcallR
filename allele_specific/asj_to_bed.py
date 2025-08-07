import sys
def convert_asj_to_bed(input_tsv, pval_threshold=1e-10):
    with open(input_tsv, "r") as infile:
        header = infile.readline().strip().split("\t")
        for line in infile:
            cols = line.strip().split("\t")
            rd = dict(zip(header, cols))
            pvalue = float(rd["P_value"])
            if pvalue >= pval_threshold:
                continue
            junction = rd["#Junction"] # Format: "chr:start-end", 1-based, both inclusive
            strand = rd["Strand"]
            gene_name = rd["Gene_name"]
            chrom, positions = junction.split(":")
            start, end = positions.split("-")
            start = str(int(start) - 1)  # Convert to 0-based BED format
            extra_info = ";".join([f"{h}={v}" for h, v in zip(header, cols)])
            print(f"{chrom}\t{start}\t{end}\t{gene_name}\t{pvalue}\t{strand}\t{extra_info}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python asj_to_bed.py <input_tsv> [p_value_threshold]")
        sys.exit(1)

    input_tsv = sys.argv[1]
    pval_threshold = float(sys.argv[2]) if len(sys.argv) > 2 else 1e-10

    convert_asj_to_bed(input_tsv, pval_threshold)