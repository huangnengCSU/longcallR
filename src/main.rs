use std::collections::HashMap;

use crate::thread::run;
use crate::util::*;
use clap::{ArgAction, Parser};
use rust_lapper::Interval;

mod ase;
mod asj;
mod candidate;
mod exon;
mod fragment;
mod low_frac;
mod phase;
mod snp;
mod snpfrags;
mod strandedness;
mod thread;
mod util;
mod vcf;

const ALL_GENE_TYPES_SENTINEL: &str = "__LONGCALLR_ALL_GENE_TYPES__";

fn inject_empty_multi_value_flag_sentinel(
    argv: &mut Vec<String>,
    long_flag: &str,
    short_flag: Option<&str>,
    sentinel: &str,
) {
    if argv.len() <= 1 {
        return;
    }
    let mut i: usize = 1;
    while i < argv.len() {
        let is_target_flag = argv[i] == long_flag || short_flag.map(|s| argv[i] == s).unwrap_or(false);
        if is_target_flag {
            let next_is_value = argv
                .get(i + 1)
                .map(|s| !s.starts_with('-'))
                .unwrap_or(false);
            if !next_is_value {
                argv.insert(i + 1, sentinel.to_string());
                i += 1;
            }
        }
        i += 1;
    }
}

mod constants {
    pub const MAX_BASE_QUALITY: u8 = 30;
}

pub static VALID_ALLELES: &[char] = &['A', 'C', 'G', 'T', 'a', 'c', 'g', 't'];

#[derive(clap::ValueEnum, Debug, Clone)]
pub enum Preset {
    HifiIsoseq, // PacBio HiFi, both strand
    HifiMasseq, // PacBio IsoSeq, transcript strand
    OntCdna,    // Oxford Nanopore, both strand
    OntDrna,    // Oxford Nanopore direct RNA, transcript strand
}

#[derive(clap::ValueEnum, Debug, Clone)]
pub enum Platform {
    Hifi, // PacBio long-read RNA sequencing
    Ont,  // Oxford Nanopore long-read RNA sequencing
}

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input BAM file (must be sorted and indexed)
    #[arg(short = 'b', long)]
    bam_path: String,

    /// Reference FASTA file
    #[arg(short = 'f', long)]
    ref_path: String,

    /// Annotation file, GFF3 or GTF format
    #[arg(short = 'a', long)]
    annotation: Option<String>,

    /// Output file prefix
    #[arg(short = 'o', long)]
    output: String,

    /// Region to be processed. Format: chr:start-end, left-closed, right-open
    #[arg(short = 'r', long)]
    region: Option<String>,

    /// Contigs to be processed. Example: -x chr1 chr2 chr3
    #[arg(short = 'x', long, num_args(0..))]
    contigs: Option<Vec<String>>,

    /// Input vcf file as Candidate SNPs
    #[arg(short = 'v', long)]
    input_vcf: Option<String>,

    /// Number of threads [Default: 1]
    #[arg(short = 't', long)]
    threads: Option<usize>,

    /// Preset for sequencing platform, choices: hifi-isoseq, hifi-masseq, ont-cdna, ont-drna
    #[arg(short = 'p', long)]
    preset: Preset,

    /// Minimum allele frequency for high allele fraction candidate SNPs [Default: 0.20]
    #[arg(long)]
    min_allele_freq: Option<f32>,

    /// Minimum allele frequency for high allele fraction candidate SNPs include intron [Default: 0.0]
    #[arg(long)]
    min_allele_freq_include_intron: Option<f32>,

    /// Minimum allele frequency for low allele fraction candidate SNPs [Default: 0.05]
    #[arg(long)]
    low_allele_frac_cutoff: Option<f32>,

    /// Minimum allele count for low allele fraction candidate SNPs [Default: 10]
    #[arg(long)]
    low_allele_cnt_cutoff: Option<u32>,

    /// Minimum length for reads [Default: 500]
    #[arg(long)]
    min_read_length: Option<usize>,

    /// Minimum mapping quality for reads [Default: 20]
    #[arg(long)]
    min_mapq: Option<u8>,

    /// Minimim base quality for allele calling [Default: 10]
    #[arg(long)]
    min_baseq: Option<u8>,

    /// Max sequence divergence for valid reads [Default: 0.05]
    #[arg(long)]
    divergence: Option<f32>,

    /// Minimum depth for a candidate SNP [Default: 10]
    #[arg(long)]
    min_depth: Option<u32>,

    /// Maximum depth for a candidate SNP [Default: 50000]
    #[arg(long)]
    max_depth: Option<u32>,

    /// Whether to use strand bias to filter SNPs [Default: auto-detect]
    #[arg(long)]
    strand_bias: Option<bool>,

    /// Minimum QUAL for candidate SNPs [Default: 2]
    #[arg(long)]
    min_qual: Option<u32>,

    /// Ignore bases within distance to read end [Default: 20]
    #[arg(long)]
    distance_to_read_end: Option<u32>,

    /// PolyA tail length threshold for filtering [Default: 5]
    #[arg(long)]
    polya_tail_length: Option<u32>,

    /// Window size used to identify dense regions of candidate SNPs [Default: 500]
    #[arg(long)]
    dense_win_size: Option<u32>,

    /// Minimum number of candidate SNPs within the dense window to consider the region as dense [Default: 5]
    #[arg(long)]
    min_dense_cnt: Option<u32>,

    /// Minimum number of related candidate heterozygous SNPs required to perform phasing in a region [Default: 1]
    #[arg(long)]
    min_linkers: Option<u32>,

    /// Maximum number of SNPs for enumerate haplotypes [Default: 10]
    #[arg(long)]
    max_enum_snps: Option<usize>,

    /// Minimum phase score to filter candidate SNPs [Default: 8.0]
    #[arg(long)]
    min_phase_score: Option<f32>,

    /// Minimum absolute difference between haplotype assignment probabilities required for a read to be confidently assigned [Default: 0.15]
    #[arg(long)]
    min_read_assignment_diff: Option<f64>,

    /// When set, apply truncation of high coverage regions
    #[arg(long, action = ArgAction::SetTrue, default_value = "false")]
    truncation: bool,

    /// Read number threshold for region truncation [Default: 200000]
    #[arg(long)]
    truncation_coverage: Option<u32>,

    /// When set, allow downsampling of high coverage regions
    #[clap(long, action = ArgAction::SetTrue, default_value = "false")]
    downsample: bool,

    /// Downsample depth [Default: 10000]
    #[arg(long)]
    downsample_depth: Option<u32>,

    /// When set, only call SNPs in exons
    #[arg(long, action = ArgAction::SetTrue, default_value = "false")]
    exon_only: bool,

    /// When set, do not output phased bam file
    #[arg(long, action = ArgAction::SetTrue, default_value = "false")]
    no_bam_output: bool,

    /// When set, directly haplotag reads from a pre-phased input VCF without re-phasing SNPs
    #[arg(long, action = ArgAction::SetTrue, default_value = "false")]
    direct_haplotag: bool,

    /// When set, show all regions to be processed
    #[clap(long, action = ArgAction::SetTrue, default_value = "false")]
    get_blocks: bool,
}

fn build_regions(
    bam_path: &str,
    ref_path: &str,
    threads: usize,
    input_region: Option<String>,
    input_contigs: Option<Vec<String>>,
    min_mapq: u8,
    min_read_length: usize,
    divergence: f32,
    truncation: bool,
    truncation_coverage: u32,
    anno_path: Option<String>,
    exon_only: bool,
) -> (Vec<Region>, HashMap<String, Vec<Interval<u32, u8>>>) {
    let mut regions = if let Some(region_str) = input_region {
        vec![Region::new(region_str.to_string())]
    } else {
        extract_isolated_regions_parallel(
            bam_path,
            ref_path,
            threads,
            input_contigs,
            min_mapq,
            min_read_length,
            divergence,
            truncation,
            truncation_coverage,
        )
    };
    let (anno_gene_regions, anno_exon_regions) = if let Some(anno_path_str) = anno_path.clone() {
        parse_annotation(anno_path_str)
    } else {
        (HashMap::new(), HashMap::new())
    };
    if exon_only {
        // divided region into multiple regions overlapped with gene regions
        regions = intersect_gene_regions(&regions, &anno_gene_regions, threads, true);
    }
    (regions, anno_exon_regions)
}

fn main() {
    let raw_args: Vec<String> = std::env::args().collect();
    if raw_args.get(1).map(|x| x == "ase").unwrap_or(false) {
        let mut ase_argv: Vec<String> = Vec::new();
        ase_argv.push(format!("{} ase", raw_args[0]));
        ase_argv.extend(raw_args.iter().skip(2).cloned());
        inject_empty_multi_value_flag_sentinel(
            &mut ase_argv,
            "--gene-types",
            None,
            ALL_GENE_TYPES_SENTINEL,
        );
        let ase_args = ase::AseArgs::parse_from(ase_argv);
        ase::run_ase(ase_args);
        return;
    }
    if raw_args.get(1).map(|x| x == "asj").unwrap_or(false) {
        let mut asj_argv: Vec<String> = Vec::new();
        asj_argv.push(format!("{} asj", raw_args[0]));
        asj_argv.extend(raw_args.iter().skip(2).cloned());
        inject_empty_multi_value_flag_sentinel(
            &mut asj_argv,
            "--gene-types",
            Some("-g"),
            ALL_GENE_TYPES_SENTINEL,
        );
        let asj_args = asj::AsjArgs::parse_from(asj_argv);
        asj::run_asj(asj_args);
        return;
    }

    let arg = Args::parse();
    let bam_path = arg.bam_path.as_str();
    let out_bam = (arg.output.clone() + ".phased.bam").clone();
    let out_vcf = (arg.output.clone() + ".vcf").clone();
    let ref_path = arg.ref_path.as_str();
    let anno_path = arg.annotation;
    let input_vcf = arg.input_vcf;
    let input_region = arg.region;
    let input_contigs = arg.contigs;
    let preset = arg.preset;

    let get_blocks = arg.get_blocks;
    let no_bam_output = arg.no_bam_output;
    let direct_haplotag = arg.direct_haplotag;
    let downsample = arg.downsample;
    let truncation = arg.truncation;
    let exon_only = arg.exon_only;

    let platform: Platform;
    let min_depth: Option<u32>;
    let min_phase_score: Option<f32>;
    let min_read_assignment_diff: Option<f64>;
    let min_linkers: Option<u32>;
    let min_allele_freq: Option<f32>;
    let min_allele_freq_include_intron: Option<f32>;
    let distance_to_read_end: Option<u32>;
    let dense_win_size: Option<u32>;
    let min_dense_cnt: Option<u32>;
    let strand_bias: Option<bool>;
    let threads: Option<usize>;
    let max_enum_snps: Option<usize>;
    let min_mapq: Option<u8>;
    let divergence: Option<f32>;
    let min_baseq: Option<u8>;
    let min_qual: Option<u32>;
    let polya_tail_length: Option<u32>;
    let max_depth: Option<u32>;
    let truncation_coverage: Option<u32>;
    let downsample_depth: Option<u32>;
    let min_read_length: Option<usize>;
    let low_allele_frac_cutoff: Option<f32>;
    let low_allele_cnt_cutoff: Option<u32>;

    match preset {
        Preset::OntCdna => {
            platform = Platform::Ont;
            min_depth = Option::from(arg.min_depth.unwrap_or(10));
            min_phase_score = Option::from(arg.min_phase_score.unwrap_or(13.0));
            min_read_assignment_diff = Option::from(arg.min_read_assignment_diff.unwrap_or(0.0));
            min_linkers = Option::from(arg.min_linkers.unwrap_or(1));
            min_allele_freq = Option::from(arg.min_allele_freq.unwrap_or(0.20));
            min_allele_freq_include_intron =
                Option::from(arg.min_allele_freq_include_intron.unwrap_or(0.0));
            distance_to_read_end = Option::from(arg.distance_to_read_end.unwrap_or(20));
            dense_win_size = Option::from(arg.dense_win_size.unwrap_or(100));
            min_dense_cnt = Option::from(arg.min_dense_cnt.unwrap_or(5));
            strand_bias = arg.strand_bias;

            threads = Option::from(arg.threads.unwrap_or(1));
            max_enum_snps = Option::from(arg.max_enum_snps.unwrap_or(10));
            min_mapq = Option::from(arg.min_mapq.unwrap_or(20));
            divergence = Option::from(arg.divergence.unwrap_or(0.5));
            min_baseq = Option::from(arg.min_baseq.unwrap_or(10));
            min_qual = Option::from(arg.min_qual.unwrap_or(2));
            polya_tail_length = Option::from(arg.polya_tail_length.unwrap_or(5));
            max_depth = Option::from(arg.max_depth.unwrap_or(50000));
            truncation_coverage = Option::from(arg.truncation_coverage.unwrap_or(200000));
            downsample_depth = Option::from(arg.downsample_depth.unwrap_or(10000));
            min_read_length = Option::from(arg.min_read_length.unwrap_or(500));
            low_allele_frac_cutoff = Option::from(arg.low_allele_frac_cutoff.unwrap_or(0.05));
            low_allele_cnt_cutoff = Option::from(arg.low_allele_cnt_cutoff.unwrap_or(10));
            println!("Preset: ont-cdna");
        }

        Preset::OntDrna => {
            platform = Platform::Ont;
            min_depth = Option::from(arg.min_depth.unwrap_or(10));
            min_phase_score = Option::from(arg.min_phase_score.unwrap_or(13.0));
            min_read_assignment_diff = Option::from(arg.min_read_assignment_diff.unwrap_or(0.0));
            min_linkers = Option::from(arg.min_linkers.unwrap_or(1));
            min_allele_freq = Option::from(arg.min_allele_freq.unwrap_or(0.20));
            min_allele_freq_include_intron =
                Option::from(arg.min_allele_freq_include_intron.unwrap_or(0.0));
            distance_to_read_end = Option::from(arg.distance_to_read_end.unwrap_or(20));
            dense_win_size = Option::from(arg.dense_win_size.unwrap_or(100));
            min_dense_cnt = Option::from(arg.min_dense_cnt.unwrap_or(5));
            strand_bias = arg.strand_bias;

            threads = Option::from(arg.threads.unwrap_or(1));
            max_enum_snps = Option::from(arg.max_enum_snps.unwrap_or(10));
            min_mapq = Option::from(arg.min_mapq.unwrap_or(20));
            divergence = Option::from(arg.divergence.unwrap_or(0.5));
            min_baseq = Option::from(arg.min_baseq.unwrap_or(10));
            min_qual = Option::from(arg.min_qual.unwrap_or(2));
            polya_tail_length = Option::from(arg.polya_tail_length.unwrap_or(5));
            max_depth = Option::from(arg.max_depth.unwrap_or(50000));
            truncation_coverage = Option::from(arg.truncation_coverage.unwrap_or(200000));
            downsample_depth = Option::from(arg.downsample_depth.unwrap_or(10000));
            min_read_length = Option::from(arg.min_read_length.unwrap_or(500));
            low_allele_frac_cutoff = Option::from(arg.low_allele_frac_cutoff.unwrap_or(0.05));
            low_allele_cnt_cutoff = Option::from(arg.low_allele_cnt_cutoff.unwrap_or(10));
            println!("Preset: ont-drna");
        }

        Preset::HifiIsoseq => {
            platform = Platform::Hifi;
            min_depth = Option::from(arg.min_depth.unwrap_or(6));
            min_phase_score = Option::from(arg.min_phase_score.unwrap_or(11.0));
            min_read_assignment_diff = Option::from(arg.min_read_assignment_diff.unwrap_or(0.0));
            min_linkers = Option::from(arg.min_linkers.unwrap_or(1));
            min_allele_freq = Option::from(arg.min_allele_freq.unwrap_or(0.15));
            min_allele_freq_include_intron =
                Option::from(arg.min_allele_freq_include_intron.unwrap_or(0.0));
            distance_to_read_end = Option::from(arg.distance_to_read_end.unwrap_or(40));
            dense_win_size = Option::from(arg.dense_win_size.unwrap_or(100));
            min_dense_cnt = Option::from(arg.min_dense_cnt.unwrap_or(5));
            strand_bias = arg.strand_bias;

            threads = Option::from(arg.threads.unwrap_or(1));
            max_enum_snps = Option::from(arg.max_enum_snps.unwrap_or(10));
            min_mapq = Option::from(arg.min_mapq.unwrap_or(20));
            divergence = Option::from(arg.divergence.unwrap_or(0.5));
            min_baseq = Option::from(arg.min_baseq.unwrap_or(10));
            min_qual = Option::from(arg.min_qual.unwrap_or(2));
            polya_tail_length = Option::from(arg.polya_tail_length.unwrap_or(5));
            max_depth = Option::from(arg.max_depth.unwrap_or(50000));
            truncation_coverage = Option::from(arg.truncation_coverage.unwrap_or(200000));
            downsample_depth = Option::from(arg.downsample_depth.unwrap_or(10000));
            min_read_length = Option::from(arg.min_read_length.unwrap_or(500));
            low_allele_frac_cutoff = Option::from(arg.low_allele_frac_cutoff.unwrap_or(0.05));
            low_allele_cnt_cutoff = Option::from(arg.low_allele_cnt_cutoff.unwrap_or(10));
            println!("Preset: hifi-isoseq");
        }

        Preset::HifiMasseq => {
            platform = Platform::Hifi;
            min_depth = Option::from(arg.min_depth.unwrap_or(6));
            min_phase_score = Option::from(arg.min_phase_score.unwrap_or(11.0));
            min_read_assignment_diff = Option::from(arg.min_read_assignment_diff.unwrap_or(0.0));
            min_linkers = Option::from(arg.min_linkers.unwrap_or(1));
            min_allele_freq = Option::from(arg.min_allele_freq.unwrap_or(0.15));
            min_allele_freq_include_intron =
                Option::from(arg.min_allele_freq_include_intron.unwrap_or(0.0));
            distance_to_read_end = Option::from(arg.distance_to_read_end.unwrap_or(40));
            dense_win_size = Option::from(arg.dense_win_size.unwrap_or(100));
            min_dense_cnt = Option::from(arg.min_dense_cnt.unwrap_or(5));
            strand_bias = arg.strand_bias;

            threads = Option::from(arg.threads.unwrap_or(1));
            max_enum_snps = Option::from(arg.max_enum_snps.unwrap_or(10));
            min_mapq = Option::from(arg.min_mapq.unwrap_or(20));
            divergence = Option::from(arg.divergence.unwrap_or(0.5));
            min_baseq = Option::from(arg.min_baseq.unwrap_or(10));
            min_qual = Option::from(arg.min_qual.unwrap_or(2));
            polya_tail_length = Option::from(arg.polya_tail_length.unwrap_or(5));
            max_depth = Option::from(arg.max_depth.unwrap_or(50000));
            truncation_coverage = Option::from(arg.truncation_coverage.unwrap_or(200000));
            downsample_depth = Option::from(arg.downsample_depth.unwrap_or(10000));
            min_read_length = Option::from(arg.min_read_length.unwrap_or(500));
            low_allele_frac_cutoff = Option::from(arg.low_allele_frac_cutoff.unwrap_or(0.05));
            low_allele_cnt_cutoff = Option::from(arg.low_allele_cnt_cutoff.unwrap_or(10));
            println!("Preset: hifi-masseq");
        }
    }

    // check if set exon_only, anno_path should be set
    if exon_only && anno_path.is_none() {
        panic!("--exon-only is set, but annotation file is not provided");
    }

    if get_blocks {
        let (regions, _) = build_regions(
            bam_path,
            ref_path,
            threads.unwrap(),
            input_region,
            input_contigs,
            min_mapq.unwrap(),
            min_read_length.unwrap(),
            divergence.unwrap(),
            truncation,
            truncation_coverage.unwrap(),
            anno_path,
            exon_only,
        );
        for reg in regions.iter() {
            if reg.gene_id.is_none() {
                println!(
                    "{}:{}-{} {}",
                    reg.chr,
                    reg.start,
                    reg.end,
                    reg.max_coverage.unwrap()
                );
            } else {
                println!(
                    "{}:{}-{} {} {:?}",
                    reg.chr,
                    reg.start,
                    reg.end,
                    reg.max_coverage.unwrap(),
                    reg.gene_id.clone().unwrap()
                );
            }
        }
        return;
    }

    if direct_haplotag && input_vcf.is_none() {
        panic!("direct_haplotag is set, but input_vcf is not provided");
    }

    let (regions, exon_regions) = build_regions(
        bam_path,
        ref_path,
        threads.unwrap(),
        input_region,
        input_contigs,
        min_mapq.unwrap(),
        min_read_length.unwrap(),
        divergence.unwrap(),
        truncation,
        truncation_coverage.unwrap(),
        anno_path,
        exon_only,
    );

    run(
        bam_path,
        ref_path,
        input_vcf,
        out_vcf.as_str(),
        out_bam.as_str(),
        threads.unwrap(),
        regions,
        exon_regions,
        exon_only,
        &platform,
        min_mapq.unwrap(),
        min_baseq.unwrap(),
        divergence.unwrap(),
        min_allele_freq.unwrap(),
        min_qual.unwrap(),
        min_allele_freq_include_intron.unwrap(),
        strand_bias,
        min_depth.unwrap(),
        max_depth.unwrap(),
        downsample,
        downsample_depth.unwrap(),
        min_read_length.unwrap(),
        distance_to_read_end.unwrap(),
        polya_tail_length.unwrap(),
        dense_win_size.unwrap(),
        min_dense_cnt.unwrap(),
        min_linkers.unwrap(),
        min_phase_score.unwrap(),
        max_enum_snps.unwrap(),
        min_read_assignment_diff.unwrap(),
        direct_haplotag,
        no_bam_output,
        low_allele_frac_cutoff.unwrap(),
        low_allele_cnt_cutoff.unwrap(),
    );
}
