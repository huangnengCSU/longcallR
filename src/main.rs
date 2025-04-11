use std::collections::HashMap;

use crate::thread::run;
use crate::util::*;
use clap::{ArgAction, Parser};
use rand::seq::SliceRandom;
use rust_htslib::bam::Read;
use rust_lapper::Interval;

mod candidate;
mod exon;
mod fragment;
mod phase;
mod snp;
mod snpfrags;
mod somatic;
mod thread;
mod util;
mod vcf;

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
    /// Path to input bam file
    #[arg(short = 'b', long)]
    bam_path: String,

    /// Path to reference file
    #[arg(short = 'f', long)]
    ref_path: String,

    /// annotation file
    #[arg(short = 'a', long)]
    annotation: Option<String>,

    /// Output bam file path
    #[arg(short = 'o', long)]
    output: String,

    /// Region to realign (Optional). Format: chr:start-end, left-closed, right-open.
    #[arg(short = 'r', long)]
    region: Option<String>,

    /// Contigs to be processed. Example: -x chr1 chr2 chr3
    #[arg(short = 'x', long, num_args(0..))]
    contigs: Option<Vec<String>>,

    /// Input vcf file
    #[arg(short = 'v', long)]
    input_vcf: Option<String>,

    /// Number of threads [Default: 1]
    #[arg(short = 't', long)]
    threads: Option<usize>,

    /// Preset of parameters, choices: hifi-isoseq, hifi-masseq, ont-cdna, ont-drna
    #[arg(short = 'p', long)]
    preset: Preset,

    /// Maximum number of iteration for phasing [Default: 100]
    #[arg(long)]
    max_iters: Option<i32>,

    /// Maximum number of SNPs for enumerate haplotypes [Default: 10]
    #[arg(long)]
    max_enum_snps: Option<usize>,

    /// Random flip fraction for snps and fragments [Default: 0.2]
    #[arg(long)]
    random_flip_fraction: Option<f32>,

    /// Minimum mapping quality for reads [Default: 20]
    #[arg(long)]
    min_mapq: Option<u8>,

    /// Minimim base quality for allele calling [Default: 10]
    #[arg(long)]
    min_baseq: Option<u8>,

    /// Max sequence divergence for valid reads [Default: 0.05]
    #[arg(long)]
    divergence: Option<f32>,

    /// Minimum allele frequency for candidate SNPs [Default: 0.20]
    #[arg(long)]
    min_allele_freq: Option<f32>,

    /// Minimum allele frequency for candidate SNPs include intron [Default: 0.0]
    #[arg(long)]
    min_allele_freq_include_intron: Option<f32>,

    /// Minimum QUAL for candidate SNPs [Default: 2]
    #[arg(long)]
    min_qual: Option<u32>,

    /// Whether to use strand bias to filter SNPs [Default: false]
    #[arg(long)]
    strand_bias: Option<bool>,

    /// Ignore bases within distance to read end [Default: 20]
    #[arg(long)]
    distance_to_read_end: Option<u32>,

    /// PolyA tail length threshold for filtering [Default: 5]
    #[arg(long)]
    polya_tail_length: Option<u32>,

    /// Dense window size for candidate SNPs [Default: 500]
    #[arg(long)]
    dense_win_size: Option<u32>,

    /// Minimum dense cnt for candidate SNPs [Default: 5]
    #[arg(long)]
    min_dense_cnt: Option<u32>,

    /// Minimum linked heterozygous snps for phasing [Default: 1]
    #[arg(long)]
    min_linkers: Option<u32>,

    /// Minimum phase score to filter SNPs [Default: 8.0]
    #[arg(long)]
    min_phase_score: Option<f32>,

    /// Minimum depth to filter SNPs [Default: 10]
    #[arg(long)]
    min_depth: Option<u32>,

    /// Maximum depth to filter SNPs [Default: 50000]
    #[arg(long)]
    max_depth: Option<u32>,

    /// When set, apply truncation of high coverage regions
    #[arg(long, action = ArgAction::SetTrue, default_value = "false")]
    truncation: bool,

    /// coverage for truncation [Default: 200000]
    #[arg(long)]
    truncation_coverage: Option<u32>,

    /// When set, allow downsampling of high coverage regions
    #[clap(long, action = ArgAction::SetTrue, default_value = "false")]
    downsample: bool,

    /// Downsample depth [Default: 10000]
    #[arg(long)]
    downsample_depth: Option<u32>,

    /// Minimum read length to filter reads [Default: 500]
    #[arg(long)]
    min_read_length: Option<usize>,

    /// Minimum absolute difference |p(h1) - p(h2)| to consider assignment legal [Default: 0.15]
    #[arg(long)]
    min_read_assignment_diff: Option<f64>,

    /// Somatic mutation allele fraction cutoff [Default: 0.05]
    #[arg(long)]
    somatic_allele_frac_cutoff: Option<f32>,

    /// Somatic mutation allele count cutoff [Default: 10]
    #[arg(long)]
    somatic_allele_cnt_cutoff: Option<u32>,

    /// When set, do not output phased bam file.
    #[arg(long, action = ArgAction::SetTrue, default_value = "false")]
    no_bam_output: bool,

    /// When set, show blocks information.
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
    let exon_regions = if let Some(anno_path_str) = anno_path {
        let (anno_gene_regions, anno_exon_regions) = parse_annotation(anno_path_str);
        regions = intersect_gene_regions(&regions, &anno_gene_regions, threads);
        anno_exon_regions
    } else {
        HashMap::new()
    };
    (regions, exon_regions)
}

fn main() {
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
    let downsample = arg.downsample;
    let truncation = arg.truncation;

    let mut platform = Platform::Hifi;
    let mut min_depth = arg.min_depth;
    let mut min_phase_score = arg.min_phase_score;
    let mut read_assignment_cutoff = arg.min_read_assignment_diff;
    let mut min_linkers = arg.min_linkers;
    let mut min_allele_freq = arg.min_allele_freq;
    let mut min_allele_freq_include_intron = arg.min_allele_freq_include_intron;
    let mut distance_to_read_end = arg.distance_to_read_end;
    let mut dense_win_size = arg.dense_win_size;
    let mut min_dense_cnt = arg.min_dense_cnt;
    let mut strand_bias = arg.strand_bias;



    let mut threads = arg.threads;
    let mut max_iters = arg.max_iters;
    let mut max_enum_snps = arg.max_enum_snps;
    let mut random_flip_fraction = arg.random_flip_fraction;
    let mut min_mapq = arg.min_mapq;
    let mut divergence = arg.divergence;
    let mut min_baseq = arg.min_baseq;
    let mut min_qual = arg.min_qual;
    let mut polya_tail_length = arg.polya_tail_length;
    let mut max_depth = arg.max_depth;
    let mut truncation_coverage = arg.truncation_coverage;
    let mut downsample_depth = arg.downsample_depth;
    let mut min_read_length = arg.min_read_length;
    let mut somatic_allele_frac_cutoff = arg.somatic_allele_frac_cutoff;
    let mut somatic_allele_cnt_cutoff = arg.somatic_allele_cnt_cutoff;




    match preset {
        Preset::OntCdna => {
            platform = Platform::Ont;
            min_depth  = Option::from(arg.min_depth.unwrap_or(10));
            min_phase_score = Option::from(arg.min_phase_score.unwrap_or(14.0));
            read_assignment_cutoff = Option::from(arg.min_read_assignment_diff.unwrap_or(0.15));
            min_linkers = Option::from(arg.min_linkers.unwrap_or(1));
            min_allele_freq = Option::from(arg.min_allele_freq.unwrap_or(0.20));
            min_allele_freq_include_intron = Option::from(arg.min_allele_freq_include_intron.unwrap_or(0.05));
            distance_to_read_end = Option::from(arg.distance_to_read_end.unwrap_or(20));
            dense_win_size = Option::from(arg.dense_win_size.unwrap_or(500));
            min_dense_cnt = Option::from(arg.min_dense_cnt.unwrap_or(5));
            strand_bias = Option::from(arg.strand_bias.unwrap_or(true));

            threads = Option::from(arg.threads.unwrap_or(1));
            max_iters = Option::from(arg.max_iters.unwrap_or(100));
            max_enum_snps = Option::from(arg.max_enum_snps.unwrap_or(10));
            random_flip_fraction = Option::from(arg.random_flip_fraction.unwrap_or(0.2));
            min_mapq = Option::from(arg.min_mapq.unwrap_or(20));
            divergence = Option::from(arg.divergence.unwrap_or(0.05));
            min_baseq = Option::from(arg.min_baseq.unwrap_or(10));
            min_qual = Option::from(arg.min_qual.unwrap_or(2));
            polya_tail_length = Option::from(arg.polya_tail_length.unwrap_or(5));
            max_depth = Option::from(arg.max_depth.unwrap_or(50000));
            truncation_coverage = Option::from(arg.truncation_coverage.unwrap_or(200000));
            downsample_depth = Option::from(arg.downsample_depth.unwrap_or(10000));
            min_read_length = Option::from(arg.min_read_length.unwrap_or(500));
            somatic_allele_frac_cutoff = Option::from(arg.somatic_allele_frac_cutoff.unwrap_or(0.05));
            somatic_allele_cnt_cutoff = Option::from(arg.somatic_allele_cnt_cutoff.unwrap_or(10));
            println!("Preset: ont-cdna");
        }

        Preset::OntDrna => {
            platform = Platform::Ont;
            min_depth  = Option::from(arg.min_depth.unwrap_or(10));
            min_phase_score = Option::from(arg.min_phase_score.unwrap_or(14.0));
            read_assignment_cutoff = Option::from(arg.min_read_assignment_diff.unwrap_or(0.15));
            min_linkers = Option::from(arg.min_linkers.unwrap_or(2));
            min_allele_freq = Option::from(arg.min_allele_freq.unwrap_or(0.20));
            min_allele_freq_include_intron = Option::from(arg.min_allele_freq_include_intron.unwrap_or(0.05));
            distance_to_read_end = Option::from(arg.distance_to_read_end.unwrap_or(20));
            dense_win_size = Option::from(arg.dense_win_size.unwrap_or(500));
            min_dense_cnt = Option::from(arg.min_dense_cnt.unwrap_or(5));
            strand_bias = Option::from(arg.strand_bias.unwrap_or(false));

            threads = Option::from(arg.threads.unwrap_or(1));
            max_iters = Option::from(arg.max_iters.unwrap_or(100));
            max_enum_snps = Option::from(arg.max_enum_snps.unwrap_or(10));
            random_flip_fraction = Option::from(arg.random_flip_fraction.unwrap_or(0.2));
            min_mapq = Option::from(arg.min_mapq.unwrap_or(20));
            divergence = Option::from(arg.divergence.unwrap_or(0.05));
            min_baseq = Option::from(arg.min_baseq.unwrap_or(10));
            min_qual = Option::from(arg.min_qual.unwrap_or(2));
            polya_tail_length = Option::from(arg.polya_tail_length.unwrap_or(5));
            max_depth = Option::from(arg.max_depth.unwrap_or(50000));
            truncation_coverage = Option::from(arg.truncation_coverage.unwrap_or(200000));
            downsample_depth = Option::from(arg.downsample_depth.unwrap_or(10000));
            min_read_length = Option::from(arg.min_read_length.unwrap_or(500));
            somatic_allele_frac_cutoff = Option::from(arg.somatic_allele_frac_cutoff.unwrap_or(0.05));
            somatic_allele_cnt_cutoff = Option::from(arg.somatic_allele_cnt_cutoff.unwrap_or(10));
            println!("Preset: ont-drna");
        }

        Preset::HifiIsoseq => {
            platform = Platform::Hifi;
            min_depth  = Option::from(arg.min_depth.unwrap_or(6));
            min_phase_score = Option::from(arg.min_phase_score.unwrap_or(11.0));
            read_assignment_cutoff = Option::from(arg.min_read_assignment_diff.unwrap_or(0.0));
            min_linkers = Option::from(arg.min_linkers.unwrap_or(1));
            min_allele_freq = Option::from(arg.min_allele_freq.unwrap_or(0.15));
            min_allele_freq_include_intron = Option::from(arg.min_allele_freq_include_intron.unwrap_or(0.0));
            distance_to_read_end = Option::from(arg.distance_to_read_end.unwrap_or(40));
            dense_win_size = Option::from(arg.dense_win_size.unwrap_or(100));
            min_dense_cnt = Option::from(arg.min_dense_cnt.unwrap_or(5));
            strand_bias = Option::from(arg.strand_bias.unwrap_or(true));

            threads = Option::from(arg.threads.unwrap_or(1));
            max_iters = Option::from(arg.max_iters.unwrap_or(100));
            max_enum_snps = Option::from(arg.max_enum_snps.unwrap_or(10));
            random_flip_fraction = Option::from(arg.random_flip_fraction.unwrap_or(0.2));
            min_mapq = Option::from(arg.min_mapq.unwrap_or(20));
            divergence = Option::from(arg.divergence.unwrap_or(0.05));
            min_baseq = Option::from(arg.min_baseq.unwrap_or(10));
            min_qual = Option::from(arg.min_qual.unwrap_or(2));
            polya_tail_length = Option::from(arg.polya_tail_length.unwrap_or(5));
            max_depth = Option::from(arg.max_depth.unwrap_or(50000));
            truncation_coverage = Option::from(arg.truncation_coverage.unwrap_or(200000));
            downsample_depth = Option::from(arg.downsample_depth.unwrap_or(10000));
            min_read_length = Option::from(arg.min_read_length.unwrap_or(500));
            somatic_allele_frac_cutoff = Option::from(arg.somatic_allele_frac_cutoff.unwrap_or(0.05));
            somatic_allele_cnt_cutoff = Option::from(arg.somatic_allele_cnt_cutoff.unwrap_or(10));
            println!("Preset: hifi-isoseq");
        }

        Preset::HifiMasseq => {
            platform = Platform::Hifi;
            min_depth  = Option::from(arg.min_depth.unwrap_or(6));
            min_phase_score = Option::from(arg.min_phase_score.unwrap_or(11.0));
            read_assignment_cutoff = Option::from(arg.min_read_assignment_diff.unwrap_or(0.0));
            min_linkers = Option::from(arg.min_linkers.unwrap_or(1));
            min_allele_freq = Option::from(arg.min_allele_freq.unwrap_or(0.15));
            min_allele_freq_include_intron = Option::from(arg.min_allele_freq_include_intron.unwrap_or(0.0));
            distance_to_read_end = Option::from(arg.distance_to_read_end.unwrap_or(40));
            dense_win_size = Option::from(arg.dense_win_size.unwrap_or(100));
            min_dense_cnt = Option::from(arg.min_dense_cnt.unwrap_or(5));
            strand_bias = Option::from(arg.strand_bias.unwrap_or(false));

            threads = Option::from(arg.threads.unwrap_or(1));
            max_iters = Option::from(arg.max_iters.unwrap_or(100));
            max_enum_snps = Option::from(arg.max_enum_snps.unwrap_or(10));
            random_flip_fraction = Option::from(arg.random_flip_fraction.unwrap_or(0.2));
            min_mapq = Option::from(arg.min_mapq.unwrap_or(20));
            divergence = Option::from(arg.divergence.unwrap_or(0.05));
            min_baseq = Option::from(arg.min_baseq.unwrap_or(10));
            min_qual = Option::from(arg.min_qual.unwrap_or(2));
            polya_tail_length = Option::from(arg.polya_tail_length.unwrap_or(5));
            max_depth = Option::from(arg.max_depth.unwrap_or(50000));
            truncation_coverage = Option::from(arg.truncation_coverage.unwrap_or(200000));
            downsample_depth = Option::from(arg.downsample_depth.unwrap_or(10000));
            min_read_length = Option::from(arg.min_read_length.unwrap_or(500));
            somatic_allele_frac_cutoff = Option::from(arg.somatic_allele_frac_cutoff.unwrap_or(0.05));
            somatic_allele_cnt_cutoff = Option::from(arg.somatic_allele_cnt_cutoff.unwrap_or(10));
            println!("Preset: hifi-masseq");
        }
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
        &platform,
        max_iters.unwrap(),
        min_mapq.unwrap(),
        min_baseq.unwrap(),
        divergence.unwrap(),
        min_allele_freq.unwrap(),
        min_qual.unwrap(),
        min_allele_freq_include_intron.unwrap(),
        strand_bias.unwrap(),
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
        random_flip_fraction.unwrap(),
        read_assignment_cutoff.unwrap(),
        no_bam_output,
        somatic_allele_frac_cutoff.unwrap(),
        somatic_allele_cnt_cutoff.unwrap(),
    );
}
