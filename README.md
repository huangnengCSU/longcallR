## Getting Started
```sh
# download and build longcallR
git clone https://github.com/huangnengCSU/longcallR.git
cd longcallR
cargo build --release

# call small variants from Nanopore cDNA reads
./longcallR -b input.bam -f ref.fa -o output -t 8 --platform ont --preset ont-cdna

# call small variants from Nanopore dRNA reads
./longcallR -b input.bam -f ref.fa -o output -t 8 --platform ont --preset ont-drna

# call small variants from PacBio iso-seq reads
./longcallR -b input.bam -f ref.fa -o output -t 8 --platform hifi --preset hifi-isoseq

# call small variants from PacBio mas-seq reads
./longcallR -b input.bam -f ref.fa -o output -t 8 --platform hifi --preset hifi-masseq
```

## Table of Contents
- [Introduction](#introduction)
- [Compiling](#compiling)
- [Usage](#usage)
- [Results](#results)
- [Demo](#demo)
- [License](#license)

## Introduction
LongcallR is a small variant caller for single molecule long-read RNA-seq data. LongcallR supports Nanopore cDNA sequecing and dRNA sequencing, PacBio Iso-Seq and MAS-Seq.

## Compiling

LongcallR is written in [Rust](https://www.rust-lang.org) and uses [Cargo](https://doc.rust-lang.org/cargo/getting-started/installation.html) to build as follows:
```
git clone https://github.com/huangnengCSU/longcallR.git
cd longcallR
cargo build --release
```

## Usage

General usage
```
./longcallR \
--bam-path input.bam \                  ## The alignment bam file
--ref-path ref.fa \                     ## The reference file must be indexed.
--platform ${PLATFORM} \                ## options: {ont, hifi}
--preset ${PRESET} \                    ## option: {ont-cdna, ont-drna, hifi-isoseq, hifi-masseq}
--output ${OUTPUT_DIR}/${PREFIX}        ## output path and prefix of output files
```

## Results

The following table shows the result of several datasets. [WTC-11 Iso-Seq](https://zenodo.org/records/5920920), [WTC-11 ONT](https://www.encodeproject.org/experiments/ENCSR539ZXJ/), [HG002 MAS-Seq](https://downloads.pacbcloud.com/public/dataset/Kinnex-full-length-RNA/) and [HG004 MAS-Seq](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data_RNAseq/AshkenazimTrio/HG004_NA24143_mother/PacBio_Pacbio-MASseq/) are public available. The ground truth of WTC-11 were described in Mark D. Robinson et.al [paper](https://link.springer.com/article/10.1186/s13059-023-02923-y). The ground truth of HG002 and HG004 were from GIAB [NISTv4.2.1](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/). All truth variants used for benchmarking have at least 10X coverage in the input. The results in sub-table (all sites) are obtained by [hap.py](https://github.com/Illumina/hap.py.git). The sub-table (non-A-to-G sites) is achieved by removing all A-G and T-C substitutions from the results of hap.py. The sub-table (evaluation without genotype error) is for separating genotype errors from SNP discovery errors.
![alt text](image.png)

## Demo

```
./longcallR -b demo/demo.bam -f demo/chr20.fa -o demo/test -t 8 --platform hifi --preset hifi-masseq
```

## License
MIT License

Copyright (c) 2024 Heng Li, Neng Huang.

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