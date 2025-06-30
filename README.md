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
cargo build --release # The executable will be located at: longcallR/target/release/longcallR

# call SNPs
longcallR -b input.bam -f ref.fa -o output -t 8 -p ont-cdna       # Nanopore cDNA reads
longcallR -b input.bam -f ref.fa -o output -t 8 -p ont-drna       # Nanopore dRNA reads (no strand bias filtering)
longcallR -b input.bam -f ref.fa -o output -t 8 -p hifi-isoseq    # PacBio iso-seq reads
longcallR -b input.bam -f ref.fa -o output -t 8 -p hifi-masseq    # PacBio mas-seq reads (no strand bias filtering)

# Allele-specific junction analysis
python allele_specific/longcallR-asj.py -a annotation.gtf -b phased.bam -f ref.fa -o output_prefix -t threads

# Allele-specific expression analysis
python allele_specific/longcallR-ase.py -a annotation.gtf -b phased.bam -o output_prefix -t threads
```

## Table of Contents
- [Introduction](#introduction)
- [Installation](#installation)
- [Usage](#usage)
- [Demo](#demo)
- [Alignment](#alignment)
- [Citation](#citation)
- [License](#license)

## Introduction
LongcallR is a SNP caller for single molecule long-read RNA-seq data. LongcallR supports Nanopore cDNA sequecing and dRNA sequencing, PacBio Iso-Seq and MAS-Seq.

## Installation

LongcallR is written in [Rust](https://www.rust-lang.org) and uses [Cargo](https://doc.rust-lang.org/cargo/getting-started/installation.html) to build as follows:
```
git clone https://github.com/huangnengCSU/longcallR.git
cd longcallR
cargo build --release
```
Alternatively, install directly from crates.io:
```
cargo install longcallR
```

## Usage

General usage
```
./longcallR \
--bam-path input.bam \                  ## The alignment bam file
--ref-path ref.fa \                     ## The reference file must be indexed.
--preset ${PRESET} \                    ## option: {ont-cdna, ont-drna, hifi-isoseq, hifi-masseq}
--output ${OUTPUT_DIR}/${PREFIX}        ## output path and prefix of output files
```

## Demo
```
./longcallR -b demo/demo.bam -f demo/chr20.fa -o demo/test -t 8 -p hifi-masseq
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

## Citation
If you use LongcallR in your work or analysis, please cite the preprint:
> Neng Huang, Heng Li, SNP calling, haplotype phasing and allele-specific analysis with long RNA-seq reads. *bioRxiv*, 2025. [doi.org/10.1101/2025.05.26.656191](https://doi.org/10.1101/2025.05.26.656191)


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