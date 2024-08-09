#!/bin/bash

# Function to display usage information
usage() {
    echo "Usage: $0 -v <vcf_file> -b <bam_file> -r <reference_file> -o <output_dir> -t <threads>"
    echo "  -v  Path to the VCF file"
    echo "  -b  Path to the BAM file"
    echo "  -r  Path to the reference genome"
    echo "  -o  Path to the output directory"
    echo "  -t  Number of threads"
    exit 1
}

# Parse command-line arguments
while getopts ":v:b:r:o:t:" opt; do
    case ${opt} in
        v )
            vcf=$OPTARG
            ;;
        b )
            bam=$OPTARG
            ;;
        r )
            ref=$OPTARG
            ;;
        o )
            out=$OPTARG
            ;;
        t )
            threads=$OPTARG
            ;;
        \? )
            usage
            ;;
    esac
done

# Check if all arguments are provided
if [ -z "$vcf" ] || [ -z "$bam" ] || [ -z "$ref" ] || [ -z "$out" ] || [ -z "$threads" ]; then
    usage
fi

# Create output directory if it doesn't exist
mkdir -p "$out"

# Chromosomes, chr1, chr2, chr3, ..., chr22
CHR=()
for i in {1..22}; do
    CHR+=("chr${i}")
done

# Split the VCF file by chromosome
parallel --joblog "${out}/splited_vcf.log" -j"$threads" \
"bcftools view -r {1} -f PASS $vcf > $out/{1}.splited.vcf" ::: "${CHR[@]}"

# Whatshap phase
parallel --joblog "${out}/whatshap_phase.log" -j"$threads" \
"whatshap phase \
--only-snvs \
--ignore-read-groups \
--reference $ref \
--chromosome {1} \
-o ${out}/phased_{1}.vcf \
$out/{1}.splited.vcf \
$bam" ::: "${CHR[@]}"

# Merge all phased VCF files
merged_vcf="${out}/merged_phased.vcf"
bcftools concat ${out}/phased_chr*.vcf > "$merged_vcf"

# Sort the merged VCF file
sorted_vcf="${out}/merged_phased_sorted.vcf"
bcftools sort "$merged_vcf" > "$sorted_vcf"

# Compress the sorted VCF file
compressed_vcf="${sorted_vcf}.gz"
bgzip -c "$sorted_vcf" > "$compressed_vcf"

# Create an index for the compressed VCF file
tabix -p vcf "$compressed_vcf"

echo "Phasing, merging, sorting, compressing, and indexing completed."