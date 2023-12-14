mod align;
mod align2;
mod bam_reader;
mod base_matrix;
mod isolated_region;
mod matrix;
mod pileup2matrix;
mod util;
mod profile;
mod runt;
mod phase;
mod vcf;

// extern crate bio;
use clap::{Parser, ArgAction};
use bam_reader::{BamReader, Region};
use rust_htslib::{bam, bam::Read, bam::Format, bam::record::Aux};
use std::time::{Duration, Instant};
// use bam_reader::{write_read_records1, write_read_records2, write_read_records3};
use crate::base_matrix::*;
use crate::matrix::ColumnBaseCount;
use crate::util::*;
use align::nw_splice_aware;
use bio::io::fasta;
use isolated_region::{find_isolated_regions};
use matrix::PileupMatrix;
use pileup2matrix::generate_pileup_matrix;
use rust_htslib::bam::record::CigarString;
use std::collections::HashMap;
use std::io;
use std::process::exit;
use std::sync::{mpsc, Arc, Mutex};
use std::thread;
use rand::seq::SliceRandom;
use crate::profile::{*};
use crate::phase::{*};
use crate::vcf::VCFRecord;


#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Path to input bam file
    #[arg(short = 'b', long)]
    bam_path: String,

    /// Path to reference file
    #[arg(short = 'f', long)]
    ref_path: String,

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

    /// Maximum number of SNPs for enumerate haplotypes
    #[arg(long, default_value_t = 10)]
    max_enum_snps: usize,

    /// Random flip fraction for snps and fragments
    #[arg(long, default_value_t = 0.2)]
    random_flip_fraction: f32,

    /// Minimum allele frequency for candidate SNPs
    #[arg(long, default_value_t = 0.25)]
    min_allele_freq: f32,

    /// Minimum allele frequency for candidate SNPs include intron
    #[arg(long, default_value_t = 0.05)]
    min_allele_freq_include_intron: f32,

    /// Minimum allele frequency for homozygous SNPs
    #[arg(long, default_value_t = 0.85)]
    min_homozygous_freq: f32,

    /// Minimum support number for each allele
    #[arg(long, default_value_t = 3)]
    min_allele_cnt: u32,

    /// Variants strand bias threshold to filter SNPs, most of the variant allele appear on one strand
    #[arg(long, default_value_t = 0.9)]
    strand_bias_threshold: f32,

    /// Cover reads strand bias threshold to filter SNPs
    #[arg(long, default_value_t = 0.9)]
    cover_strand_bias_threshold: f32,

    /// Minimum phase score to filter SNPs
    #[arg(long, default_value_t = 8.0)]
    min_phase_score: f32,

    /// Minimum depth to filter SNPs
    #[arg(long, default_value_t = 10)]
    min_depth: u32,

    /// Read assignment cutoff, the read is phased only if the probability of assignment P(hap1)-P(hap2) > cutoff or P(hap2)-P(hap1) > cutoff
    #[arg(long, default_value_t = 0.15)]
    read_assignment_cutoff: f64,

    /// Without phasing, only using genotype probability
    #[clap(long, action = ArgAction::SetTrue)]
    genotype_only: bool,

    /// When set, output vcf file does not contain phase information.
    #[clap(long, action = ArgAction::SetFalse)]
    no_phase_vcf: bool,

    /// debug SNP
    #[clap(long, action = ArgAction::SetTrue)]
    debug_snp: bool,

    /// get blocks
    #[clap(long, action = ArgAction::SetTrue)]
    debug_block: bool,
}

fn main() {
    let arg = Args::parse();
    let bam_path = arg.bam_path.as_str();
    let out_bam = (arg.output.clone() + ".phased.bam").clone();
    let out_vcf = (arg.output.clone() + ".vcf").clone();
    // let output_file = arg.output.as_str();
    let ref_path = arg.ref_path.as_str();
    let input_region = arg.region;
    let input_contigs = arg.contigs;
    let threads = arg.threads;
    let genotype_only = arg.genotype_only;
    let max_enum_snps = arg.max_enum_snps;
    let random_flip_fraction = arg.random_flip_fraction;
    let min_allele_freq = arg.min_allele_freq;
    let min_allele_freq_include_intron = arg.min_allele_freq_include_intron;
    let min_allele_cnt = arg.min_allele_cnt;
    let strand_bias_threshold = arg.strand_bias_threshold;
    let cover_strand_bias_threshold = arg.cover_strand_bias_threshold;
    let min_phase_score = arg.min_phase_score;
    let min_depth = arg.min_depth;
    let read_assignment_cutoff = arg.read_assignment_cutoff;
    let min_homozygous_freq = arg.min_homozygous_freq;
    let phasing_output = arg.no_phase_vcf;  // default=true
    let debug_snp = arg.debug_snp; // default=false
    let debug_block = arg.debug_block; // default=false

    if debug_block {
        let regions = multithread_produce3(bam_path.to_string().clone(), threads, input_contigs);
        for reg in regions.iter() {
            println!("{}:{}-{}", reg.chr, reg.start, reg.end);
        }
        return;
    }

    if debug_snp {
        let region = Region::new(input_region.unwrap());
        let mut profile = Profile::default();
        let ref_seqs = read_references(ref_path);
        profile.init_with_pileup(bam_path, &region);
        profile.append_reference(&ref_seqs);
        // for bf in profile.freq_vec.iter() {
        //     println!("bf: {:?}", bf);
        // }
        let mut snpfrag = SNPFrag::default();
        snpfrag.get_candidate_snps(&profile, min_allele_freq, min_allele_freq_include_intron, min_depth, min_homozygous_freq, strand_bias_threshold, cover_strand_bias_threshold);
        snpfrag.filter_fp_snps(strand_bias_threshold, None);
        for i in snpfrag.hete_snps.iter() {
            println!("hete snp: {:?}", snpfrag.candidate_snps[*i]);
        }

        for i in snpfrag.homo_snps.iter() {
            println!("homo snp: {:?}", snpfrag.candidate_snps[*i]);
        }
        return;
    }

    if input_region.is_some() {
        let region = Region::new(input_region.unwrap());
        let mut profile = Profile::default();
        let mut readnames: Vec<String> = Vec::new();
        let ref_seqs = read_references(ref_path);
        profile.init_with_pileup(bam_path, &region);
        profile.append_reference(&ref_seqs);
        let mut snpfrag = SNPFrag::default();
        snpfrag.get_candidate_snps(&profile, min_allele_freq, min_allele_freq_include_intron, min_depth, min_homozygous_freq, strand_bias_threshold, cover_strand_bias_threshold);
        let mut read_assignments: HashMap<String, i32> = HashMap::new();
        snpfrag.get_fragments(bam_path, &region);
        for edge in snpfrag.edges.iter() {
            println!("edge: {:?}->{:?}", snpfrag.candidate_snps[edge.0[0]], snpfrag.candidate_snps[edge.0[1]]);
        }
        snpfrag.filter_fp_snps(strand_bias_threshold, None);
        let mut vcf_records: Vec<VCFRecord> = Vec::new();
        if genotype_only {
            // without phasing
            vcf_records = snpfrag.output_vcf3(min_allele_freq, min_homozygous_freq);
        } else {
            if snpfrag.hete_snps.len() > 0 {
                // for i in snpfrag.hete_snps.iter() {
                //     println!("hete snp: {:?}", snpfrag.candidate_snps[*i]);
                // }
                let mut v: Vec<_> = snpfrag.edges.iter().collect();
                v.sort_by(|x, y| x.0[0].cmp(&y.0[0]));
                unsafe { snpfrag.init_haplotypes(); }
                unsafe { snpfrag.init_assignment(); }
                snpfrag.phase(max_enum_snps, random_flip_fraction);
                read_assignments = snpfrag.assign_reads(read_assignment_cutoff);
                snpfrag.add_phase_score(min_allele_cnt);
            }
            vcf_records = snpfrag.output_vcf2(min_phase_score, min_homozygous_freq, phasing_output);
        }

        for rd in vcf_records.iter() {
            if rd.alternative.len() == 1 {
                println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", std::str::from_utf8(&rd.chromosome).unwrap(),
                         rd.position,
                         std::str::from_utf8(&rd.id).unwrap(),
                         std::str::from_utf8(&rd.reference).unwrap(),
                         std::str::from_utf8(&rd.alternative[0]).unwrap(),
                         rd.qual,
                         std::str::from_utf8(&rd.filter).unwrap(),
                         std::str::from_utf8(&rd.info).unwrap(),
                         std::str::from_utf8(&rd.format).unwrap(),
                         rd.genotype);
            } else if rd.alternative.len() == 2 {
                println!("{}\t{}\t{}\t{}\t{},{}\t{}\t{}\t{}\t{}\t{}", std::str::from_utf8(&rd.chromosome).unwrap(),
                         rd.position,
                         std::str::from_utf8(&rd.id).unwrap(),
                         std::str::from_utf8(&rd.reference).unwrap(),
                         std::str::from_utf8(&rd.alternative[0]).unwrap(),
                         std::str::from_utf8(&rd.alternative[1]).unwrap(),
                         rd.qual,
                         std::str::from_utf8(&rd.filter).unwrap(),
                         std::str::from_utf8(&rd.info).unwrap(),
                         std::str::from_utf8(&rd.format).unwrap(),
                         rd.genotype);
            }
        }

        let mut bam_reader = bam::Reader::from_path(&bam_path).unwrap();
        let header = bam::Header::from_template(&bam_reader.header());
        let mut bam_writer = bam::Writer::from_path(out_bam, &header, Format::Bam).unwrap();
        for r in bam_reader.records() {
            let mut record = r.unwrap();
            if record.is_unmapped() || record.is_secondary() || record.is_supplementary() {
                continue;
            }
            let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
            if read_assignments.contains_key(&qname) {
                let asg = read_assignments.get(&qname).unwrap();
                if *asg != 0 {
                    let _ = record.push_aux(b"HP:i", Aux::I32(*asg));
                }
            }
            let _ = bam_writer.write(&record).unwrap();
        }
    } else {
        let regions = multithread_produce3(bam_path.to_string().clone(), threads, input_contigs);
// multithread_phase_maxcut(bam_path.to_string().clone(), ref_path.to_string().clone(), output_file.to_string().clone(), threads, regions);
        multithread_phase_haplotag(bam_path.to_string().clone(),
                                   ref_path.to_string().clone(),
                                   out_vcf.clone(),
                                   out_bam.clone(),
                                   threads,
                                   regions,
                                   genotype_only,
                                   min_allele_freq,
                                   min_allele_freq_include_intron,
                                   min_allele_cnt,
                                   strand_bias_threshold,
                                   cover_strand_bias_threshold,
                                   min_depth,
                                   min_homozygous_freq,
                                   min_phase_score,
                                   max_enum_snps,
                                   random_flip_fraction,
                                   read_assignment_cutoff,
                                   phasing_output);
    }
}
