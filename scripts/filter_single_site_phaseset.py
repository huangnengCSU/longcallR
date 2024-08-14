import argparse

import pysam


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Filter and transform phased VCF records with single-site phase sets.")
    parser.add_argument("-i", "--input_vcf", help="Input VCF file", required=True)
    parser.add_argument("-o", "--output_vcf", help="Output VCF file", required=True)

    return parser.parse_args()



def collect_phase_sets(vcf):
    """Collect phase sets from the VCF file."""
    phase_sets = {}  # key: phase set id, value: list of positions
    for record in vcf:
        chr = record.chrom
        pos = record.pos
        samples = record.samples
        if "PS" not in samples[0].keys():
            continue
        phase_set = samples[0]['PS']
        phase_sets.setdefault(phase_set, []).append(f"{chr}:{pos}")
    return phase_sets


def filter_and_transform_records(input_vcf, output_vcf, phase_sets):
    """Filter specific positions and transform phased records to non-phased records."""
    vcf = pysam.VariantFile(input_vcf)
    out_vcf = pysam.VariantFile(output_vcf, 'w', header=vcf.header)

    # Identify positions with single-site phase sets
    filter_sites = [positions[0] for phase_set, positions in phase_sets.items() if len(positions) == 1]

    for record in vcf:
        chr = record.chrom
        pos = record.pos
        site_id = f"{chr}:{pos}"

        # If the site is in the list of filter sites, remove the 'PS' field
        if site_id in filter_sites:
            if "PS" in record.samples[0]:
                del record.samples[0]['PS']
                del record.format['PS']
                gt = record.samples[0]['GT']
                if gt is not None and record.samples[0].phased:
                    # Convert phased genotype to unphased
                    record.samples[0]['GT'] = gt
                    record.samples[0].phased = False

        out_vcf.write(record)

    print(f"Filtered {len(filter_sites)} sites")


def main():
    args = parse_args()
    vcf = pysam.VariantFile(args.input_vcf)

    # Collect phase sets
    phase_sets = collect_phase_sets(vcf)

    # Filter and transform records
    filter_and_transform_records(args.input_vcf, args.output_vcf, phase_sets)


if __name__ == "__main__":
    main()
