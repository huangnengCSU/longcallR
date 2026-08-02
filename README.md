<p align="center">
  <img src="longcallR.png" alt="longcallR logo" style="width:85%">
</p>

<p align="center">
  <a href="https://crates.io/crates/longcallR">
    <img src="https://img.shields.io/crates/v/longcallR.svg" alt="Crates.io version">
  </a>
  <a href="https://opensource.org/licenses/MIT">
    <img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="License: MIT">
  </a>
  <a href="https://crates.io/crates/longcallR">
    <img src="https://img.shields.io/crates/d/longcallR.svg" alt="Crates.io Downloads">
  </a>
  <a href="https://github.com/huangnengCSU/longcallR/releases">
    <img src="https://img.shields.io/github/downloads/huangnengCSU/longcallR/total.svg" alt="GitHub Downloads">
  </a>
  <a href="https://github.com/huangnengCSU/longcallR">
    <img src="https://img.shields.io/github/stars/huangnengCSU/longcallR.svg?style=social&label=Star" alt="GitHub Stars">
  </a>
</p>

## Getting Started
```sh
# download and build
git clone https://github.com/huangnengCSU/longcallR.git
cd longcallR
cargo build --release # The executable will be located at: target/release/longcallR

# call SNPs and phasing
longcallR -b input.bam -f ref.fa -o output -t 8 -p ont-cdna       # Nanopore cDNA reads
longcallR -b input.bam -f ref.fa -o output -t 8 -p ont-drna       # Nanopore dRNA reads
longcallR -b input.bam -f ref.fa -o output -t 8 -p hifi-isoseq    # PacBio iso-seq reads
longcallR -b input.bam -f ref.fa -o output -t 8 -p hifi-masseq    # PacBio mas-seq reads

# Note: By default, longcallR detects whether each region contains single- or double-stranded reads.
# Strand bias filtering is enabled only for confidently double-stranded regions; ambiguous regions remain disabled.
# Enabling it explicitly for single-stranded data may cause many false negatives.
longcallR -b input.bam -f ref.fa -o output -t 8 -p <preset> --strand-bias true

# Allele-specific junction analysis
longcallR asj -a annotation.gtf.gz -b phased.bam -f ref.fa -o output_prefix -t threads

# Allele-specific expression analysis
longcallR ase -a annotation.gtf.gz -b phased.bam -o output_prefix -t threads
```

## Table of Contents
- [Introduction](#introduction)
- [Installation](#installation)
- [Usage](#usage)
- [Demo](#demo)
- [Alignment](#alignment)
- [Citation](#citation)
- [Update Log](#update-log)
- [License](#license)

## Introduction
longcallR is a tool for SNP calling, haplotype phasing, and allele-specific analysis with long-read RNA-seq data. It supports ONT cDNA and direct RNA sequencing, as well as PacBio Iso-Seq and MAS-Seq.

Full documentation is available at [huangnengCSU.github.io/longcallR](https://huangnengCSU.github.io/longcallR/).

## Installation

### Install from source
LongcallR is written in [Rust](https://www.rust-lang.org) and uses [Cargo](https://doc.rust-lang.org/cargo/getting-started/installation.html) to build as follows:
```bash
git clone https://github.com/huangnengCSU/longcallR.git
cd longcallR
cargo build --release

INSTALL_DIR="$(pwd)/target/release"
export PATH="$PATH:$INSTALL_DIR"

# Check Installation
longcallR -h
```

### Install from Conda
longcallR can be installed via Conda:
```bash
mamba install longcallr

# Check Installation
longcallR -h
```

### Install from Crates.io
longcallR can be installed from crates.io:
```bash
cargo install longcallR

# Check Installation
longcallR -h
```

## Usage

### Basic Example
```bash
longcallR \
--bam-path input.bam \                  ## The alignment bam file
--ref-path ref.fa \                     ## The reference file must be indexed.
--preset ${PRESET} \                    ## option: {ont-cdna, ont-drna, hifi-isoseq, hifi-masseq}
--output ${OUTPUT_PREFIX}               ## Prefix of output files
```
For advanced options, see the [documentation](https://huangnengCSU.github.io/longcallR/)

Advanced option: if a phased VCF is already available, `--direct-haplotag` can be used with `--input-vcf` to directly haplotag reads without re-phasing SNPs.

### Allele specific analysis
```bash
# call allele-specific expression
longcallR ase -b <phased_bam> -a <annotation> -o <output_prefix> -t <threads>

# call allele-specific junction
longcallR asj -b <phased_bam> -a <annotation> -f <reference> -o <output_prefix> -t <threads>

# store allele-specific junction in BED format for IGV visualization
asj_to_bed.py output.asj.tsv [p_value_threshold] > output.asj.bed
```
For detailed information on available options and output file formats, please refer to the [documentation](https://huangnengCSU.github.io/longcallR/)

## Demo
```
longcallR -b demo/demo.bam -f demo/chr20.fa -o demo/test -t 8 -p hifi-masseq
```

## Alignment
Based on the strand orientation of reads from different PacBio and Nanopore protocols, we recommend the following alignment parameters for [**minimap2**](https://github.com/lh3/minimap2):
```sh
minimap2 -ax splice:hq -uf ref.fa query.fa > aln.sam    # PacBio Kinnex/Mas-seq
minimap2 -ax splice:hq ref.fa query.fa > aln.sam        # PacBio Iso-seq
minimap2 -ax splice -uf -k14 ref.fa reads.fa > aln.sam  # Nanopore dRNA
minimap2 -ax splice ref.fa reads.fa > aln.sam           # Nanopore cDNA
```
Note: The `-uf` option forces minimap2 to consider only a single transcript strand during alignment.
-	PacBio Kinnex/MAS-Seq and Nanopore dRNA sequencing produce single-stranded reads.
-	Nanopore cDNA sequencing produces double-stranded reads.
-	PacBio Iso-Seq datasets may be either single- or double-stranded.

We recommend using the `-uf` option for single-stranded reads and omitting it for double-stranded reads.

## Results:
The results of manuscript can be found and reproduced using the data and scripts available at: [Zenodo](https://doi.org/10.5281/zenodo.17842979).

## Citation
If you use LongcallR in your work or analysis, please cite the preprint:
> Neng Huang, Heng Li, Human Pangenome Reference Consortium, SNP calling, haplotype phasing and allele-specific analysis with long RNA-seq reads. *nature methods*, 2026. [https://doi.org/10.1038/s41592-026-03045-6](https://doi.org/10.1038/s41592-026-03045-6)

## Update Log

### Unreleased - 2026-08-02
- Added automatic per-region detection of single- and double-stranded reads using forward and reverse coverage stored in the base-frequency matrix.
- Classification requires at least 50 informative positions within a region, with a read coverage depth of at least 10 at each position.
- Strand bias filtering is now enabled automatically only for confidently double-stranded regions. It remains disabled for single-stranded or ambiguous regions to reduce false negatives.
- Explicit `--strand-bias true` or `--strand-bias false` settings override automatic detection.
- The runtime log reports whether strand bias filtering is enabled or disabled for each region and whether the decision was automatic or explicitly configured.

### 2.0.1 - 2026-04-03
- Improved robustness of `--direct-haplotag` mode: non-SNV records (indels, spanning deletions, symbolic alleles) in the input VCF are now filtered out before haplotagging to prevent mismatched interpretation at indel loci.

### 2.0.0 - 2026-03-27
- Added support for direct haplotagging from a pre-phased VCF using `--input-vcf --direct-haplotag`.
- Reimplemented the functionality of longcallR-ase.py and longcallR-asj.py in Rust and integrated them into longcallR as the `longcallR ase` and `longcallR asj` subcommands. The Python implementations will be removed in a future release.


## License
MIT License

Copyright (c) 2024 Dana-Farber Cancer Institute.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
