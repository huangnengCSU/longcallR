# Usage

## General Usage
To run `longcallR` for **SNP calling** and **read phasing**, use the following command:
```
longcallR
-b <input.sorted.bam>
-f <reference.fa>
-p <preset>
-t <threads>
-o <output_prefix>
```

+ -b, --bam-path: Input BAM file (must be sorted and indexed)
+ -f, --ref-path: Reference FASTA file
+ -p, --preset: Preset for sequencing platform
+ -t, --threads: Number of threads
+ -o, --output: Output file prefix


Available Presets

| **Preset**     | **Description**                                               |
|----------------|---------------------------------------------------------------|
| `hifi-masseq`  | PacBio Mas-Seq data (strand bias filtering **disabled**)      |
| `hifi-isoseq`  | PacBio Iso-Seq data (strand bias filtering **enabled**)       |
| `ont-cdna`     | ONT cDNA data (strand bias filtering **enabled**)             |
| `ont-drna`     | ONT direct RNA data (strand bias filtering **disabled**)      |

## Advanced Usage
To perform SNP calling and phasing in a **specific genomic region**, use the following command:
```
longcallR
-b <input.sorted.bam>
-f <reference.fa>
-p <preset>
-t <threads>
-r <chr:start-end>
-o <output_prefix>
```

To perform SNP calling and phasing based on **user-provided candidate SNPs** (e.g., known genomic SNPs), use the following command:
```
longcallR
-b <input.sorted.bam>
-f <reference.fa>
-v <input.vcf>
-p <preset>
-t <threads>
-o <output_prefix>
```

To customize thresholds for **read coverage**, **read length**, **mapping quality**, **base quality**, **allele fraction**, **allele fraction for low fraction candidates** (more challenge sites), **strand bias** filtering, use:
```
longcallR
-b <input.sorted.bam>
-f <reference.fa>
-p <preset>
--min-depth <min_depth>
--max-depth <max_depth>
--min-read-length <min_length>
--min-mapq <min_mapq>
--min-baseq <min_baseq>
--min-allele-freq <min_allele_freq>
--low-allele-frac-cutoff <low_allele_frac_cutoff>
--strand-bias <strand_bias>
-t <threads>
-o <output_prefix>
```



Full Arguments

+ -b, --bam-path: Input BAM file (must be sorted and indexed)
+ -f, --ref-path: Reference FASTA file
+ -p, --preset: Preset for sequencing platform
+ -t, --threads: Number of threads
+ -o, --output: Output file prefix
+ -a, --annotation: Annotation file, GFF3 or GTF format
+ -r, --region: Region to be processed. Format: chr:start-end, left-closed, right-open
+ -x, --contigs: Conitgs to be processed. Example: -x chr1 chr2 chr3
+ -v , --input-vcf: Input vcf file as Candidate SNPs
+ --exon-only: When set, only call SNPs in exons
+ --max-enum-snps: Maximum number of SNPs for enumerate haplotypes [Default: 10]
+ --min-read-length: Minimum length for reads [Default: 500]
+ --min-mapq: Minimum mapping quality for reads [Default: 20]
+ --min-baseq: Minimim base quality for allele calling [Default: 10]
+ --divergence: Max sequence divergence for valid reads [Default: 0.05]
+ --min-depth: Minimum depth for a candidate SNP [Default: 10]
+ --max-depth: Maximum depth for a candidate SNP [Default: 50000]
+ --min-allele-freq: Minimum allele frequency for high allele fraction candidate SNPs [Default: 0.20]
+ --min-allele-freq-include-intron: Minimum allele frequency for high allele fraction candidate SNPs include intron [Default: 0.0]
+ --min-qual: Minimum QUAL for candidate SNPs [Default: 2]
+ --strand-bias: Whether to use strand bias to filter SNPs [Default: false] [possible values: true, false]
+ --distance-to-read-end: Ignore bases within distance to read end [Default: 20]
+ --polya-tail-length: PolyA tail length threshold for filtering [Default: 5]
+ --dense-win-size: Window size used to identify dense regions of candidate SNPs [Default: 500]
+ --min-dense-cnt: Minimum number of candidate SNPs within the dense window to consider the region as dense [Default: 5]
+ --min-linkers: Minimum number of related candidate heterozygous SNPs required to perform phasing in a region [Default: 1]
+ --min-phase-score: Minimum phase score to filter candidate SNPs [Default: 8.0]
+ --truncation: When set, apply truncation of high coverage regions
+ --truncation-coverage: Read number threshold for region truncation [Default: 200000]
+ --downsample: When set, allow downsampling of high coverage regions
+ --downsample-depth: Downsample depth [Default: 10000]
+ --min-read-assignment-diff: Minimum absolute difference between haplotype assignment probabilities required for a read to be confidently assigned [Default: 0.15]
+ --low-allele-frac-cutoff: Minimum allele frequency for low allele fraction candidate SNPs [Default: 0.05]
+ --low-allele-cnt-cutoff: Minimum allele count for low allele fraction candidate SNPs [Default: 10]
+ --no-bam-output: When set, do not output phased bam file
+ --get-blocks: When set, show all regions to be processed.





