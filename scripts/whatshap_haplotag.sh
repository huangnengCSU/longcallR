#!/bin/bash

# Function to display usage information
usage() {
    echo "Usage: $0 -b <bam_file> -v <vcf_file> -r <reference_file> -o <output_dir> -t <threads>"
    exit 1
}

# Parse command-line arguments
while getopts ":b:v:r:o:t:" opt; do
    case ${opt} in
        b )
            bam=$OPTARG
            ;;
        v )
            vcf=$OPTARG
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
if [ -z "$bam" ] || [ -z "$vcf" ] || [ -z "$ref" ] || [ -z "$out" ] || [ -z "$threads" ]; then
    usage
fi

# Create output directory if it doesn't exist
mkdir -p "$out"

# Chromosomes, chr1, chr2, chr3, ..., chr22
CHR=()
for i in {1..22}; do
    CHR+=("chr${i}")
done

# Whatshap haplotag
time parallel --joblog "${out}/whatshap_haplotag.log" -j"$threads" \
"whatshap haplotag \
--output ${out}/{1}.bam \
--reference $ref \
--ignore-read-groups \
--regions {1} \
$vcf \
$bam" ::: "${CHR[@]}"

# Merge all phased BAM files
merged_bam="${out}/merged_phased.bam"
samtools merge -@ "$threads" -f "$merged_bam" ${out}/chr*.bam

# Sort the merged BAM file
sorted_bam="${out}/merged_phased_sorted.bam"
samtools sort -@ "$threads" -o "$sorted_bam" "$merged_bam"

# Create an index for the sorted BAM file
samtools index "$sorted_bam"

echo "Haplotagging, merging, sorting, and indexing of BAM files completed."