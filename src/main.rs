use std::collections::HashMap;

use clap::{ArgAction, Parser, ValueEnum};
use rand::seq::SliceRandom;
use rust_htslib::bam::Read;

use crate::thread::multithread_phase_haplotag;
use crate::util::*;

mod snp;
mod util;
mod snpfrags;
mod somatic;
mod exon;
mod thread;
mod vcf;
mod candidate;
mod fragment;
mod phase;

#[derive(clap::ValueEnum, Debug, Clone)]
pub enum Preset {
    hifi_isoseq,
    // PacBio HiFi, both strand
    hifi_masseq,
    // PacBio IsoSeq, transcript strand
    ont_cdna,
    // Oxford Nanopore, both strand
    ont_drna, // direct RNA, transcript strand
}

#[derive(clap::ValueEnum, Debug, Clone)]
pub enum Platform {
    hifi,
    // PacBio long-read RNA sequencing
    ont, // Oxford Nanopore long-read RNA sequencing
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

    /// Number of threads, default 1
    #[arg(short = 't', long, default_value_t = 1)]
    threads: usize,

    /// Platform for sequencing reads, choices: hifi, ont
    #[arg(short = 'p', long)]
    platform: Platform,

    /// Preset of parameters, valid choices: Preset
    #[arg(long)]
    preset: Option<Preset>,

    /// Maximum number of iteration for phasing
    #[arg(long, default_value_t = 100)]
    max_iters: i32,

    /// Maximum number of SNPs for enumerate haplotypes
    #[arg(long, default_value_t = 10)]
    max_enum_snps: usize,

    /// Random flip fraction for snps and fragments
    #[arg(long, default_value_t = 0.2)]
    random_flip_fraction: f32,

    /// Minimum mapping quality for reads
    #[arg(long, default_value_t = 20)]
    min_mapq: u8,

    /// Minimim base quality for allele
    #[arg(long, default_value_t = 10)]
    min_baseq: u8,

    // /// threshold for differental average base quality of two alleles
    // #[arg(long, default_value_t = 10)]
    // diff_baseq: u8,

    /// Minimum allele frequency for candidate SNPs
    #[arg(long, default_value_t = 0.20)]
    min_allele_freq: f32,

    /// fraction threshold of heterozygous variant used for phasing, high fraction het var
    #[arg(long, default_value_t = 0.25)]
    hetvar_high_frac_cutoff: f32,

    // /// Minimum support number for each allele
    // #[arg(long, default_value_t = 2)]
    // min_allele_cnt: u32,

    /// Minimum allele frequency for candidate SNPs include intron
    #[arg(long, default_value_t = 0.0)]
    min_allele_freq_include_intron: f32,

    // /// Minimum allele frequency for homozygous SNPs
    // #[arg(long, default_value_t = 0.75)]
    // min_homozygous_freq: f32,


    // Minimum QUAL for candidate SNPs
    #[arg(long, default_value_t = 20)]
    min_qual: u32,

    // /// Minimum variant quality for candidate SNPs
    // #[arg(long, default_value_t = 200)]
    // min_qual_for_candidate: u32,

    // /// Minimum variant quality for single snp and rna editing (higher than min_qual_for_candidate)
    // #[arg(long, default_value_t = 256)]
    // min_qual_for_singlesnp_rnaedit: u32,

    /// Whether to use strand bias to filter SNPs
    #[arg(long, action = ArgAction::SetTrue, default_value = "false")]
    use_strand_bias: bool,

    /// Variants strand bias threshold to filter SNPs, most of the variant allele appear on one strand
    #[arg(long, default_value_t = 0.9)]
    strand_bias_threshold: f32,

    /// Cover reads strand bias threshold to filter SNPs
    #[arg(long, default_value_t = 0.9)]
    cover_strand_bias_threshold: f32,

    /// Ignore with distance to splicing site
    #[arg(long, default_value_t = 20)]
    distance_to_splicing_site: u32,

    /// Window size for local error rate
    #[arg(long, default_value_t = 3)]
    window_size: u32,

    /// Ignore with distance to read end
    #[arg(long, default_value_t = 20)]
    distance_to_read_end: u32,

    // /// threshold for differental average distance to read end of two alleles
    // #[arg(long, default_value_t = 200)]
    // diff_distance_to_read_end: i64,

    /// PolyA tail length threshold
    #[arg(long, default_value_t = 5)]
    polya_tail_length: u32,

    /// Dense window size
    #[arg(long, default_value_t = 500)]
    dense_win_size: u32,

    /// Minimum dense cnt
    #[arg(long, default_value_t = 5)]
    min_dense_cnt: u32,

    // /// Average dense distance
    // #[arg(long, default_value_t = 60.0)]
    // avg_dense_dist: f32,

    /// Minimum linked heterozygous snps for phasing
    #[arg(long, default_value_t = 2)]
    min_linkers: u32,

    /// Minimum phase score to filter SNPs
    #[arg(long, default_value_t = 8.0)]
    min_phase_score: f32,

    /// Minimum depth to filter SNPs
    #[arg(long, default_value_t = 10)]
    min_depth: u32,

    /// Maximum depth to filter SNPs
    #[arg(long, default_value_t = 50000)]
    max_depth: u32,

    /// Minimum read length to filter reads
    #[arg(long, default_value_t = 500)]
    min_read_length: usize,

    /// Read assignment cutoff, the read is phased only if the probability of assignment P(hap1)-P(hap2) > cutoff or P(hap2)-P(hap1) > cutoff
    #[arg(long, default_value_t = 0.15)]
    read_assignment_cutoff: f64,

    /// Imbalance allele expression cutoff, allele1 / allele2 > cutoff or allele2 / allele1 > cutoff.
    #[arg(long, default_value_t = 2.0)]
    imbalance_allele_expression_cutoff: f32,

    // /// Allele-specific expression allele fraction cutoff
    // #[arg(long, default_value_t = 0.10)]
    // ase_allele_frac_cutoff: f32,

    // /// Allele-specific expression allele count cutoff
    // #[arg(long, default_value_t = 2)]
    // ase_allele_cnt_cutoff: u32,

    // /// Allele-specific expression phased read count cutoff
    // #[arg(long, default_value_t = 10)]
    // ase_ps_read_cutoff: u32,

    // /// Allele-specific expression phase score cutoff
    // #[arg(long, default_value_t = 20.0)]
    // ase_ps_cutoff: f32,

    /// Somatic mutation allele fraction cutoff
    #[arg(long, default_value_t = 0.05)]
    somatic_allele_frac_cutoff: f32,

    /// Somatic mutation allele count cutoff
    #[arg(long, default_value_t = 10)]
    somatic_allele_cnt_cutoff: u32,

    /// Without phasing, only using genotype probability
    #[clap(long, action = ArgAction::SetTrue, default_value = "false")]
    genotype_only: bool,

    /// When set, output vcf file does not contain phase information.
    #[clap(long, action = ArgAction::SetFalse, default_value = "true")]
    no_phase_vcf: bool,

    /// When set, do not output phased bam file.
    #[arg(long, action = ArgAction::SetTrue, default_value = "false")]
    no_bam_output: bool,

    /// When set, find haplotype-specific exons.
    #[arg(long, action = ArgAction::SetTrue, default_value = "false")]
    haplotype_specific_exon: bool,

    /// minimum number of reads to support the haplotype-specific exon
    #[arg(long, default_value_t = 8)]
    min_sup_haplotype_exon: u32,

    /// haplotype_bam_output
    #[arg(long, action = ArgAction::SetTrue, default_value = "false")]
    haplotype_bam_output: bool,

    /// output read assignment
    #[arg(long, action = ArgAction::SetTrue, default_value = "false")]
    output_read_assignment: bool,

    // /// debug SNP
    // #[clap(long, action = ArgAction::SetTrue)]
    // debug_snp: bool,

    /// get blocks
    #[clap(long, action = ArgAction::SetTrue, default_value = "false")]
    debug_block: bool,
}

fn main() {
    let arg = Args::parse();
    let bam_path = arg.bam_path.as_str();
    let out_bam = (arg.output.clone() + ".phased.bam").clone();
    let out_vcf = (arg.output.clone() + ".vcf").clone();
    let ref_path = arg.ref_path.as_str();
    let anno_path = arg.annotation;
    let input_region = arg.region;
    let input_contigs = arg.contigs;
    let threads = arg.threads;
    let debug_block = arg.debug_block;
    let preset = arg.preset;
    let platform = arg.platform;
    let max_iters = arg.max_iters;
    let max_enum_snps = arg.max_enum_snps;
    let random_flip_fraction = arg.random_flip_fraction;
    let genotype_only = arg.genotype_only;
    let phasing_output = arg.no_phase_vcf; // default=true
    let no_bam_output = arg.no_bam_output; // default=false
    let haplotype_bam_output = arg.haplotype_bam_output; // default=false
    let output_read_assignment = arg.output_read_assignment; // default=false
    let haplotype_specific_exon = arg.haplotype_specific_exon; // default=false
    let min_sup_haplotype_exon = arg.min_sup_haplotype_exon;
    let min_mapq = arg.min_mapq;
    let min_baseq = arg.min_baseq;
    let min_qual = arg.min_qual;
    let hetvar_high_frac_cutoff = arg.hetvar_high_frac_cutoff;
    let strand_bias_threshold = arg.strand_bias_threshold;
    let cover_strand_bias_threshold = arg.cover_strand_bias_threshold;
    let distance_to_splicing_site = arg.distance_to_splicing_site;
    let window_size = arg.window_size;
    let polya_tail_length = arg.polya_tail_length;
    let max_depth = arg.max_depth;
    let min_read_length = arg.min_read_length;
    let imbalance_allele_expression_cutoff = arg.imbalance_allele_expression_cutoff;
    let somatic_allele_frac_cutoff = arg.somatic_allele_frac_cutoff;
    let somatic_allele_cnt_cutoff = arg.somatic_allele_cnt_cutoff;

    let mut min_allele_freq = arg.min_allele_freq;
    let mut min_allele_freq_include_intron = arg.min_allele_freq_include_intron;
    let mut use_strand_bias = arg.use_strand_bias;
    let mut distance_to_read_end = arg.distance_to_read_end;
    let mut dense_win_size = arg.dense_win_size;
    let mut min_dense_cnt = arg.min_dense_cnt;
    let mut min_linkers = arg.min_linkers;
    let mut min_phase_score = arg.min_phase_score;
    let mut min_depth = arg.min_depth;
    let mut read_assignment_cutoff = arg.read_assignment_cutoff;


    if preset.is_some() {
        match preset.unwrap() {
            Preset::ont_cdna => {
                min_depth = 10;
                min_phase_score = 14.0;
                read_assignment_cutoff = 0.15;
                min_linkers = 1;
                min_allele_freq = 0.20;
                min_allele_freq_include_intron = 0.05;
                distance_to_read_end = 20;
                dense_win_size = 500;
                min_dense_cnt = 5;
                use_strand_bias = true;
                println!("Preset: ont-cdna");
            }

            Preset::ont_drna => {
                min_depth = 10;
                min_phase_score = 14.0;
                read_assignment_cutoff = 0.15;
                min_linkers = 2;
                min_allele_freq = 0.20;
                min_allele_freq_include_intron = 0.05;
                distance_to_read_end = 20;
                dense_win_size = 500;
                min_dense_cnt = 5;
                use_strand_bias = false;
                println!("Preset: ont-drna");
            }

            Preset::hifi_isoseq => {
                min_depth = 6;
                min_phase_score = 14.0;
                read_assignment_cutoff = 0.0;
                min_linkers = 1;
                min_allele_freq = 0.15;
                min_allele_freq_include_intron = 0.0;
                distance_to_read_end = 40;
                dense_win_size = 100;
                min_dense_cnt = 5;
                use_strand_bias = true;
                println!("Preset: hifi-isoseq");
            }

            Preset::hifi_masseq => {
                min_depth = 6;
                min_phase_score = 14.0;
                read_assignment_cutoff = 0.0;
                min_linkers = 1;
                min_allele_freq = 0.15;
                min_allele_freq_include_intron = 0.0;
                distance_to_read_end = 40;
                dense_win_size = 100;
                min_dense_cnt = 5;
                use_strand_bias = false;
                println!("Preset: hifi-masseq");
            }

            _ => {
                panic!("Preset not supported");
            }
        }
    }

    match platform {
        Platform::hifi => {
            println!("Platform: PacBio HiFi");
        }
        Platform::ont => {
            println!("Platform: Oxford Nanopore");
        }
        _ => {
            panic!("Platform not supported");
        }
    }

    if debug_block {
        let mut regions = Vec::new();
        if input_region.is_some() {
            let region = Region::new(input_region.unwrap());
            regions = vec![region];
        } else {
            regions = multithread_produce3(
                bam_path.to_string().clone(),
                ref_path.to_string().clone(),
                threads,
                input_contigs,
                min_mapq,
                min_read_length,
            );
        }

        if anno_path.is_some() {
            let (gene_regions, exon_regions) = parse_annotation(anno_path.unwrap());
            regions = intersect_gene_regions(&regions, &gene_regions, threads);
        }
        for reg in regions.iter() {
            if reg.gene_id.is_none() {
                println!("{}:{}-{}", reg.chr, reg.start, reg.end);
            } else {
                println!("{}:{}-{} {:?}", reg.chr, reg.start, reg.end, reg.gene_id.clone().unwrap());
            }
        }
        return;
    }

    let mut regions = Vec::new();
    let mut exon_regions = HashMap::new();
    if input_region.is_some() {
        let region = Region::new(input_region.unwrap());
        regions = vec![region];
    } else {
        // TODO: cut weak connected regions caused by alignment error to avoid too large regions
        regions = multithread_produce3(
            bam_path.to_string().clone(),
            ref_path.to_string().clone(),
            threads,
            input_contigs,
            min_mapq,
            min_read_length,
        );
    }

    if anno_path.is_some() {
        let (gene_regions_anno, exon_regions_anno) = parse_annotation(anno_path.unwrap());
        regions = intersect_gene_regions(&regions, &gene_regions_anno, threads);
        exon_regions = exon_regions_anno;
    }

    multithread_phase_haplotag(
        bam_path.to_string().clone(),
        ref_path.to_string().clone(),
        out_vcf.clone(),
        out_bam.clone(),
        threads,
        regions,
        exon_regions,
        genotype_only,
        &platform,
        max_iters,
        min_mapq,
        min_baseq,
        min_allele_freq,
        min_qual,
        hetvar_high_frac_cutoff,
        min_allele_freq_include_intron,
        use_strand_bias,
        strand_bias_threshold,
        cover_strand_bias_threshold,
        min_depth,
        max_depth,
        min_read_length,
        distance_to_splicing_site,
        window_size,
        distance_to_read_end,
        polya_tail_length,
        dense_win_size,
        min_dense_cnt,
        min_linkers,
        min_phase_score,
        max_enum_snps,
        random_flip_fraction,
        read_assignment_cutoff,
        imbalance_allele_expression_cutoff,
        no_bam_output,
        haplotype_bam_output,
        output_read_assignment,
        haplotype_specific_exon,
        min_sup_haplotype_exon,
        somatic_allele_frac_cutoff,
        somatic_allele_cnt_cutoff,
    );
}
