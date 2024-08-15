import argparse

import pysam


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Add PS field to T2T-Q100 VCF records.")
    parser.add_argument("-i", "--input_vcf", help="Input VCF file", required=True)
    parser.add_argument("-o", "--output_vcf", help="Output VCF file", required=True)

    return parser.parse_args()


def add_ps_tag_to_header(vcf):
    # Add the PS tag to the FORMAT field of the VCF header
    vcf.header.formats.add("PS", 1, "String", "Phase set")


def transform_records(input_vcf, output_vcf):
    vcf = pysam.VariantFile(input_vcf)
    # Add the PS tag to the header if it's not already present
    if "PS" not in vcf.header.formats:
        add_ps_tag_to_header(vcf)
    out_vcf = pysam.VariantFile(output_vcf, 'w', header=vcf.header)
    for record in vcf:
        chr = record.chrom
        chr = chr.replace("chr", "")
        record.samples[0]['PS'] = chr
        out_vcf.write(record)


def main():
    args = parse_args()
    transform_records(args.input_vcf, args.output_vcf)


if __name__ == "__main__":
    main()
