# Allele-specific junction

## Allele-Specific Junction Quantification
To quantify allele-specific junction based on the LongcallR-phased BAM, run the following command:
```
python longcallR-asj.py
-b <phased_bam>
-a <annotation>
-f <reference>
-o <output_prefix>
-t <threads>
-g <gene_types>
-m <min_sup>
```
+ -b, --bam_file: longcallR phased bam file
+ -a, --annotation_file: Annotation file GTF for GFF3 format
+ -f, --reference: Reference genome file
+ -o, --output_prefix: Prefix of output files
+ -t, --threads: Number of threads
+ -g, --gene_types: Gene types to be analyzed. Default: ["protein_coding", "lncRNA"]
+ -m, --min_sup: Minimum support of phased reads for counting event. Default: 10

## High-Confidence ASJ Using Known Genomic SNPs
To increase confidence by requiring allele-specific junction to be supported by both longcallR's RNA SNPs and known genomic SNPs, use the following command:
```
python longcallR-asj.py
-b <phased_bam>
-a <annotation>
-f <reference>
--dna_vcf <genomic_vcf>
--rna_vcf <rna_vcf>
-o <output_prefix>
-t <threads>
-g <gene_types>
-m <min_sup>
```
+ --dna_vcf: DNA vcf file
+ --rna_vcf: longcallR phased RNA vcf file

## Store Allele-Specific Junctions in BED Format for IGV Visualization
To visualize allele-specific junctions in IGV, extract regions with a p-value below a specified threshold (default: 1e-10) and save them in BED format:
```
python asj_to_bed.py output.asj.tsv [p_value_threshold] > output.asj.bed
```
- output.asj.tsv: Input file containing allele-specific junctions
- p_value_threshold (optional): Significance threshold (default: 1e-10)
- output.asj.bed: Output BED file for IGV visualization


## Output format
Running `longcallR-asj.py` will generate three output files `output.asj.tsv`, `output.asj_gene.tsv` and `output.gene_coverage.tsv`. The formats of these files are described below:

1.output.asj.tsv

| Column         | Description                                                          |
|----------------|----------------------------------------------------------------------|
| `#Junction`    | Junction coordinates in the format `chr:start-end`                   |
| `Strand`       | Strand of the junction (`+` or `-`)                                  |
| `Junction_set` | Junction set identifier (cluster of related junctions)               |
| `Phase_set`    | Phase set ID                                                         |
| `Hap1_absent`  | Number of reads from haplotype 1 where the junction is absent        |
| `Hap1_present` | Number of reads from haplotype 1 where the junction is present       |
| `Hap2_absent`  | Number of reads from haplotype 2 where the junction is absent        |
| `Hap2_present` | Number of reads from haplotype 2 where the junction is present       |
| `P_value`      | P-value for allele-specific junction                                 |
| `SOR`          | Symmetric odds ratio                                                 |
| `Novel`        | Whether the junction is novel (`TRUE` / `FALSE`)                     |
| `GT_AG`        | Whether junction has canonical GT-AG splice sites (`TRUE` / `FALSE`) |
| `Gene_name`    | Associated gene symbol                                               |

2.output.asj_gene.tsv

| Column         | Description                                                                        |
|----------------|------------------------------------------------------------------------------------|
| `#Gene_name`   | Gene symbol that containing allele-specific junctions                              |
| `Chr`          | Chromosome                                                                         |
| `P_value`      | P-value for the most significant allele-specific junction in the gene              |
| `SOR`          | Symmetric odds ratio for the most significant allele-specific junction in the gene |

3.output.gene_coverage.tsv

| Column         | Description                                |
|----------------|--------------------------------------------|
| `#Gene_name`   | Gene symbol                                |
| `Chr`          | Chromosome                                 |
| `Start`        | start position of gene, 1-based, inclusive |
| `End`          | end position of gene, 1-based, inclusive   |
| `Num_reads`    | Number of reads in the gene                |

