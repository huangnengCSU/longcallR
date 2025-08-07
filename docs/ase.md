# Allele-specific expressed gene

## Allele-Specific Expression Quantification
To quantify allele-specific gene expression based on the LongcallR-phased BAM, run the following command to get `output.ase.tsv`:
```
python longcallR-ase.py
-b <phased_bam> 
-a <annotation>
-o <output_prefix>
-t <threads>
--gene_types <gene_types>
--min_support <min_support>
```
+ -b, --bam: longcallR phased bam file
+ -a, --annotation: Annotation file GTF for GFF3 format
+ -o, --output: Prefix of output files
+ -t, --threads: Number of threads
+ --gene_types: Gene types to be analyzed. Default: ["protein_coding", "lncRNA"]
+ --min_support: Minimum support reads for counting event. Default: 10

## High-Confidence ASE Using Known Genomic SNPs
To increase confidence by requiring allele-specific expression to be supported by both LongcallR SNPs and known genomic SNPs, use the following command to get `output.filter_ase.tsv`:
```
python longcallR-ase.py
-b <phased_bam>
--vcf1 <longcallR_phased_vcf>
--vcf3 <genomic_vcf>
-a <annotation>
-o <output_prefix>
-t <threads>
--gene_types <gene_types>
--min_support <min_support>
```
+ --vcf1: LongcallR phased vcf file
+ --vcf3: DNA vcf file

## Using Phased Genomic VCF for Parental Allele Resolution
To assign allele-specific expression to paternal and maternal haplotypes, provide both the LongcallR-phased VCF and a whole-genome haplotype-phased DNA VCF, get `output.patmat_ase.tsv`:
```
python longcallR-ase.py
-b <phased_bam>
--vcf1 <longcallR_phased_vcf>
--vcf2 <phased_genomic_vcf>
-a <annotation>
-o <output_prefix>
-t <threads>
--gene_types <gene_types>
--min_support <min_support>
```
+ --vcf1: LongcallR phased vcf file
+ --vcf2: Whole genome haplotype phased DNA vcf file

## Output format
1.output.ase.tsv or output.filter_ase.tsv

| Column       | Description                               |
|--------------|-------------------------------------------|
| `#Gene_name` | Gene symbol                               |
| `Chr`        | Chromosome                                |
| `PS`         | Phase set ID                              |
| `H1`         | Read count assigned to haplotype 1        |
| `H2`         | Read count assigned to haplotype 2        |
| `P_value`    | Binomial test p-value for allelic imbalance |

2.output.patmat_ase.tsv

| Column         | Description                                          |
|----------------|------------------------------------------------------|
| `#Gene_name`   | Gene symbol                                          |
| `Chr`          | Chromosome                                           |
| `PS`           | Phase set ID                                         |
| `H1`           | Read count assigned to haplotype 1                   |
| `H2`           | Read count assigned to haplotype 2                   |
| `P_value`      | Binomial test p-value for allelic imbalance          |
| `H1_Paternal`  | Read count in H1 supporting the paternal allele      |
| `H1_Maternal`  | Read count in H1 supporting the maternal allele      |
| `H2_Paternal`  | Read count in H2 supporting the paternal allele      |
| `H2_Maternal`  | Read count in H2 supporting the maternal allele      |
