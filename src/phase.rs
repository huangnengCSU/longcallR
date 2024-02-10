use std::collections::{HashMap, HashSet, VecDeque};
use std::{fs};
use std::fs::File;
use std::hash::Hash;
use std::io::{BufRead, BufReader, Write};
use crate::bam_reader::Region;
use crate::util::{parse_fai, Profile};
use rust_htslib::{bam, bam::Read, bam::record::Record, bam::Format, bam::record::Aux};
use rust_htslib::htslib::{drand48};
use std::sync::{Mutex};
use bio::bio_types::strand::ReqStrand::Forward;
use rayon::prelude::*;
use rand::seq::SliceRandom;
use crate::base_matrix::load_reference;
use crate::vcf::VCFRecord;
use std::cmp::{max, Ordering};
use chrono::Local;


#[derive(Debug, Clone, Default)]
pub struct CandidateSNP {
    pub chromosome: Vec<u8>,
    pub pos: i64,
    // position on the reference, 0-based
    pub alleles: [char; 2],
    // major and minor alleles
    pub allele_freqs: [f32; 2],
    // major and minor allele frequencies
    pub reference: char,
    pub depth: u32,
    pub variant_type: i32,
    // 1: heterozygous SNP, 2: homozygous SNP, 3: triallelic SNP
    pub variant_quality: f64,
    // the confidence that the variant exists at this site given the data, phred-scaled
    pub genotype_probability: [f64; 3],
    // 0th: homo var, 1st: hete var, 2nd: homo ref
    pub genotype_quality: f64,
    pub haplotype: i32,
    // delta: 1,-1: hap1, hap2 if phased, 0: unassigned
    pub phase_score: f64,
    pub snp_cover_fragments: Vec<usize>,
    // index of the fragment cover this SNP
    pub rna_editing: bool,
    // A->G, T-C, rna_editing variant does not have haplotype information, no phasing
    pub filter: bool,
    // filter out this SNP or not in phasing process
    pub allele_imbalance: bool,
    // imbalanced allele expression
}

#[derive(Debug, Clone, Default)]
pub struct Edge {
    pub snp_idxes: [usize; 2],
    // index of candidate SNPs(SNPFrag.snps), start node and end node
    pub snp_poses: [i64; 2],
    // position of candidate SNPs(SNPFrag.snps), start node and end node
    pub frag_idxes: Vec<usize>,
    // index of fragments(SNPFrag.fragments) cover this edge.
    pub w: f64,
    // weight of edge,  w_{ij}=\sum_{k}x_{ki}x_{kj}log\frac{1-\epsilon_{kij}}{\epsilon_{kij}}
}

#[derive(Debug, Clone, Default)]
pub struct FragElem {
    pub snp_idx: usize,
    // index of candidate SNPs(SNPFrag.snps)
    pub pos: i64,
    // position on the reference, 0-based
    pub base: char,
    // base pair
    pub baseq: u8,
    // base quality
    pub strand: u32,
    // read strand,  0: forward, 1: reverse
    pub p: i32,
    // base allele on alphabet  {-1, 1, 0}, 1: base==alleles[0], -1: base==alleles[1], 0: not covered (not major variant allele and reference allele, deletions or N)
    pub prob: f64,
    // error rate of observe current base
}

#[derive(Debug, Clone, Default)]
pub struct Fragment {
    pub fragment_idx: usize,
    // index of multiple fragments(SNPFrag.fragments)
    pub read_id: String,
    // read name
    pub list: Vec<FragElem>,
    // single fragment
    pub haplotag: i32,
    // sigma: 0,1,-1
    pub assignment: i32,
    // haplotype assignment of the fragment, 0,1,2. 0: unassigned, 1: hap1, 2: hap2
    pub assignment_score: f64,
    // probability of the haplotype assignment
    pub exons: Vec<[i64; 2]>,
    // exons of the read on the reference, 0-based, [start, end)
}

#[derive(Debug, Clone, Default)]
pub struct SNPFrag {
    pub region: Region,
    pub candidate_snps: Vec<CandidateSNP>,
    // candidate SNPs
    pub hete_snps: Vec<usize>,
    // index of candidate heterozygous SNPs
    pub homo_snps: Vec<usize>,
    // index of candidate homozygous SNPs
    pub fragments: Vec<Fragment>,
    // multiple fragments
    pub phased: bool,
    // haplotype is phased or not
    pub edges: HashMap<[usize; 2], Edge>,
    // edges of the graph, key is [snp_idx of start_node, snp_idx of end_node]
}

fn cmp_f64(a: &f64, b: &f64) -> Ordering {
    if a < b {
        return Ordering::Less;
    } else if a > b {
        return Ordering::Greater;
    }
    return Ordering::Equal;
}

impl SNPFrag {
    pub fn get_candidate_snps(&mut self,
                              profile: &Profile,
                              min_allele_freq: f32,
                              min_allele_freq_include_intron: f32,
                              min_qual_for_candidate: u32,
                              min_coverage: u32,
                              max_coverage: u32,
                              min_baseq: u8,
                              min_homozygous_freq: f32,
                              no_strand_bias: bool,
                              strand_bias_threshold: f32,
                              cover_strand_bias_threshold: f32,
                              distance_to_splicing_site: u32,
                              window_size: u32,
                              distance_to_read_end: u32,
                              diff_distance_to_read_end: i64,
                              diff_baseq: u8,
                              dense_win_size: u32,
                              min_dense_cnt: u32,
                              avg_dense_dist: f32) {
        // get candidate SNPs, filtering with min_coverage, deletion_freq, min_allele_freq_include_intron, cover_strand_bias_threshold
        let pileup = &profile.freq_vec;
        let mut position = profile.region.start - 1;    // 0-based
        for bfidx in 0..pileup.len() {
            let bf = &pileup[bfidx];
            // println!("{:?}", bf);
            if bf.i {
                // insertion base
                continue;
            }

            // 1.filtering with depth
            let depth = bf.get_depth_exclude_intron_deletion();
            if depth < min_coverage {
                // println!("depth: {}, {}", position, depth);
                position += 1;
                continue;
            }
            if depth > max_coverage {
                // println!("depth: {}, {}", position, depth);
                position += 1;
                continue;
            }

            let (allele1, allele1_cnt, allele2, allele2_cnt) = bf.get_two_major_alleles();

            // filtering average distance to read end is significant different for allele1 and allele2
            // filtering average base quality is significant different for allele1 and allele2
            let mut allele1_dists: Vec<i64> = Vec::new();
            let mut allele2_dists: Vec<i64> = Vec::new();
            let mut allele1_quals: Vec<u8> = Vec::new();
            let mut allele2_quals: Vec<u8> = Vec::new();
            match allele1 {
                'a' => {
                    allele1_dists = bf.distance_to_end.a.clone();
                    allele1_quals = bf.baseq.a.clone();
                }
                'A' => {
                    allele1_dists = bf.distance_to_end.a.clone();
                    allele1_quals = bf.baseq.a.clone();
                }
                'c' => {
                    allele1_dists = bf.distance_to_end.c.clone();
                    allele1_quals = bf.baseq.c.clone();
                }
                'C' => {
                    allele1_dists = bf.distance_to_end.c.clone();
                    allele1_quals = bf.baseq.c.clone();
                }
                'g' => {
                    allele1_dists = bf.distance_to_end.g.clone();
                    allele1_quals = bf.baseq.g.clone();
                }
                'G' => {
                    allele1_dists = bf.distance_to_end.g.clone();
                    allele1_quals = bf.baseq.g.clone();
                }
                't' => {
                    allele1_dists = bf.distance_to_end.t.clone();
                    allele1_quals = bf.baseq.t.clone();
                }
                'T' => {
                    allele1_dists = bf.distance_to_end.t.clone();
                    allele1_quals = bf.baseq.t.clone();
                }
                _ => {
                    println!("Error: unknown allele");
                }
            }
            match allele2 {
                'a' => {
                    allele2_dists = bf.distance_to_end.a.clone();
                    allele2_quals = bf.baseq.a.clone();
                }
                'A' => {
                    allele2_dists = bf.distance_to_end.a.clone();
                    allele2_quals = bf.baseq.a.clone();
                }
                'c' => {
                    allele2_dists = bf.distance_to_end.c.clone();
                    allele2_quals = bf.baseq.c.clone();
                }
                'C' => {
                    allele2_dists = bf.distance_to_end.c.clone();
                    allele2_quals = bf.baseq.c.clone();
                }
                'g' => {
                    allele2_dists = bf.distance_to_end.g.clone();
                    allele2_quals = bf.baseq.g.clone();
                }
                'G' => {
                    allele2_dists = bf.distance_to_end.g.clone();
                    allele2_quals = bf.baseq.g.clone();
                }
                't' => {
                    allele2_dists = bf.distance_to_end.t.clone();
                    allele2_quals = bf.baseq.t.clone();
                }
                'T' => {
                    allele2_dists = bf.distance_to_end.t.clone();
                    allele2_quals = bf.baseq.t.clone();
                }
                _ => {
                    println!("Error: unknown allele");
                }
            }


            {
                let mut avg_alleles_baseq = 0.0;
                for baseq in allele1_quals.iter() {
                    avg_alleles_baseq += *baseq as f32;
                }
                for baseq in allele2_quals.iter() {
                    avg_alleles_baseq += *baseq as f32;
                }
                avg_alleles_baseq /= (allele1_quals.len() + allele2_quals.len()) as f32;
                if avg_alleles_baseq < min_baseq as f32 {
                    println!("{}:{} low allele average baseq: {}", profile.region.chr, position, avg_alleles_baseq);
                    position += 1;
                    continue;
                }
            }

            // let mut allele1_median_dist = 0;
            // let mut allele2_median_dist = 0;
            // if allele1_dists.len() >= 5 && allele2_dists.len() >= 5 {
            //     allele1_dists.sort();
            //     allele2_dists.sort();
            //     let mid1 = allele1_dists.len() / 2;
            //     let mid2 = allele2_dists.len() / 2;
            //     allele1_median_dist = allele1_dists[mid1];
            //     allele2_median_dist = allele2_dists[mid2];
            //     if (allele1_median_dist - allele2_median_dist).abs() > diff_distance_to_read_end {
            //         println!("{} allele1_median_dist: {}, allele2_median_dist: {}", position, allele1_median_dist, allele2_median_dist);
            //         position += 1;
            //         continue;
            //     }
            // }

            // let mut allele1_median_qual = 0;
            // let mut allele2_median_qual = 0;
            // if allele1_dists.len() >= 5 && allele2_dists.len() >= 5 {
            //     allele1_quals.sort();
            //     allele2_quals.sort();
            //     let mid1 = allele1_quals.len() / 2;
            //     let mid2 = allele2_quals.len() / 2;
            //     allele1_median_qual = allele1_quals[mid1];
            //     allele2_median_qual = allele2_quals[mid2];
            //     if (allele1_median_qual as i32 - allele2_median_qual as i32).abs() as u8 > diff_baseq {
            //         println!("{} allele1_median_qual: {}, allele2_median_qual: {}", position, allele1_median_qual, allele2_median_qual);
            //         position += 1;
            //         continue;
            //     }
            // }


            // 2.filtering with depth, considering intron reads
            let depth_include_intron = bf.get_depth_include_intron();
            if (allele1_cnt as f32) / (depth_include_intron as f32) < min_allele_freq_include_intron {
                // maybe caused by erroneous intron alignment
                println!("allele freq include intron: {}:{}, {}, {}, {}", profile.region.chr, position, allele1_cnt, allele2_cnt, depth_include_intron);
                position += 1;
                continue;
            }

            // 3.filtering with deletion frequency
            if bf.d > allele1_cnt {
                println!("indels: {}:{}, {}, {}, {}", profile.region.chr, position, bf.d, allele1_cnt, allele2_cnt);
                position += 1;
                continue;
            }

            if !no_strand_bias {

                // 4.filtering snps only covered by one strand reads (may caused by intron alignment error)
                let total_cover_cnt = bf.forward_cnt + bf.backward_cnt;   // does not include intron reads
                if bf.forward_cnt as f32 / total_cover_cnt as f32 > cover_strand_bias_threshold || bf.backward_cnt as f32 / total_cover_cnt as f32 > cover_strand_bias_threshold {
                    // strand bias
                    println!("cover strand bias: {}:{}, {}, {}, {}", profile.region.chr, position, bf.forward_cnt, bf.backward_cnt, total_cover_cnt);
                    // println!("{:?}", bf);
                    position += 1;
                    continue;
                }

                // 5.filtering with strand bias
                let mut variant_allele: Vec<char> = Vec::new();
                // if allele1 != bf.ref_base {
                if allele1 != bf.ref_base && allele1_cnt >= 4 {
                    variant_allele.push(allele1);
                }
                // if allele2 != bf.ref_base {
                if allele2 != bf.ref_base && allele2_cnt >= 4 {
                    variant_allele.push(allele2);
                }
                if variant_allele.len() > 0 {
                    // variant allele
                    let mut strand_bias = false;
                    for allele_base in variant_allele.iter() {
                        let mut fcnt = 0;
                        let mut bcnt = 0;
                        match allele_base {
                            'a' => { [fcnt, bcnt] = bf.base_strands.a; }
                            'A' => { [fcnt, bcnt] = bf.base_strands.a; }
                            'c' => { [fcnt, bcnt] = bf.base_strands.c; }
                            'C' => { [fcnt, bcnt] = bf.base_strands.c; }
                            'g' => { [fcnt, bcnt] = bf.base_strands.g; }
                            'G' => { [fcnt, bcnt] = bf.base_strands.g; }
                            't' => { [fcnt, bcnt] = bf.base_strands.t; }
                            'T' => { [fcnt, bcnt] = bf.base_strands.t; }
                            _ => {
                                println!("Error: unknown allele");
                                position += 1;
                                continue;
                            }
                        }
                        let total_cnt = fcnt + bcnt;
                        if fcnt as f32 / total_cnt as f32 > strand_bias_threshold || bcnt as f32 / total_cnt as f32 > strand_bias_threshold {
                            // strand bias
                            println!("strand bias: {}:{}, {}, {}, {}", profile.region.chr, position, fcnt, bcnt, total_cnt);
                            // println!("{:?}", bf);
                            strand_bias = true;
                            continue;
                        }
                    }
                    if strand_bias {
                        position += 1;
                        continue;
                    }
                }
            }


            // // 7. filtering close to read end
            // let mut left_depths: Vec<u32> = Vec::new();
            // let mut right_depths: Vec<u32> = Vec::new();
            // let mut lext = 0;
            // let mut lidx = bfidx as i32;
            // let mut rext = 0;
            // let mut ridx = bfidx;
            // let mut num_variant_allele = 0;
            // if allele1 != bf.ref_base {
            //     num_variant_allele += allele1_cnt;
            // }
            // if allele2 != bf.ref_base {
            //     num_variant_allele += allele2_cnt;
            // }
            // while lext <= distance_to_read_end {
            //     if lidx < 0 {
            //         break;
            //     }
            //     let lbf = &pileup[lidx as usize];
            //     if lbf.i {
            //         // ignore insertion region
            //         lidx -= 1;
            //         continue;
            //     }
            //     if lbf.d as f32 / (lbf.a + lbf.c + lbf.g + lbf.t + lbf.d) as f32 >= 0.7 {
            //         // ignore deletion region
            //         lidx -= 1;
            //         continue;
            //     }
            //     left_depths.push(lbf.a + lbf.c + lbf.g + lbf.t + lbf.d + lbf.n);
            //     lidx -= 1;
            //     lext += 1;
            // }
            // println!("{} left_depths: \n{:?}", position, left_depths);
            // while rext <= distance_to_read_end {
            //     if ridx >= pileup.len() {
            //         break;
            //     }
            //     let rbf = &pileup[ridx];
            //     if rbf.i {
            //         // ignore insertion region
            //         ridx += 1;
            //         continue;
            //     }
            //     if rbf.d as f32 / (rbf.a + rbf.c + rbf.g + rbf.t + rbf.d) as f32 >= 0.7 {
            //         // ignore deletion region
            //         ridx += 1;
            //         continue;
            //     }
            //     right_depths.push(rbf.a + rbf.c + rbf.g + rbf.t + rbf.d + rbf.n);
            //     ridx += 1;
            //     rext += 1;
            // }
            // println!("{} right_depths: \n{:?}", position, right_depths);
            //
            //
            // if left_depths.len() > 0 {
            //     let max_left_depth = left_depths.iter().max().unwrap().clone() as i32;
            //     let min_left_depth = left_depths.iter().min().unwrap().clone() as i32;
            //     if ((max_left_depth - min_left_depth) > (num_variant_allele as f32 * 0.9) as i32) && (min_left_depth < left_depths[0] as i32) {
            //         println!("close to left read end: {}, {}, {}", position, max_left_depth, min_left_depth);
            //         position += 1;
            //         continue;
            //     }
            // }
            //
            //
            // if right_depths.len() > 0 {
            //     let max_right_depth = right_depths.iter().max().unwrap().clone() as i32;
            //     let min_right_depth = right_depths.iter().min().unwrap().clone() as i32;
            //     if ((max_right_depth - min_right_depth) > (num_variant_allele as f32 * 0.9) as i32) && (min_right_depth < right_depths[0] as i32) {
            //         println!("close to right read end: {}, {}, {}", position, max_right_depth, min_right_depth);
            //         position += 1;
            //         continue;
            //     }
            // }


            // filtering by local high error rate
            let mut local_misalignment_ratio: Vec<f32> = Vec::new();
            let mut lext = 1;
            let mut lidx = bfidx as i32 - 1;
            let mut rext = 1;
            let mut ridx = bfidx + 1;
            while lext <= window_size {
                if lidx < 0 {
                    break;
                }
                let lbf = &pileup[lidx as usize];
                if lbf.i {
                    // ignore insertion region
                    lidx -= 1;
                    continue;
                }
                // println!("Left bases: {}, {}, {}, {}, {}", lbf.a, lbf.c, lbf.g, lbf.t, lbf.d);
                let local_error_rate = lbf.get_none_ref_count() as f32 / (lbf.a + lbf.c + lbf.g + lbf.t + lbf.d) as f32;
                local_misalignment_ratio.push(local_error_rate);
                lext += 1;
                lidx -= 1;
            }
            while rext <= window_size {
                if ridx >= pileup.len() {
                    break;
                }
                let rbf = &pileup[ridx];
                if rbf.i {
                    // ignore insertion region
                    ridx += 1;
                    continue;
                }
                // println!("Right bases: {}, {}, {}, {}, {}", rbf.a, rbf.c, rbf.g, rbf.t, rbf.d);
                let local_error_rate = rbf.get_none_ref_count() as f32 / (rbf.a + rbf.c + rbf.g + rbf.t + rbf.d) as f32;
                local_misalignment_ratio.push(local_error_rate);
                rext += 1;
                ridx += 1;
            }
            // println!("{:?}", bf);
            // println!("{} local_misalignment_ratio: \n{:?}", position, local_misalignment_ratio);

            let mut N_cnts: Vec<u32> = Vec::new();  // number of N bases (intron)
            let mut INS_cnts: Vec<u32> = Vec::new();  // number of insertions
            let mut lext = 0;
            let mut lidx = bfidx as i32;
            let mut rext = 1;
            let mut ridx = bfidx + 1;
            let mut num_variant_allele = 0;
            if allele1 != bf.ref_base {
                num_variant_allele += allele1_cnt;
            }
            if allele2 != bf.ref_base {
                num_variant_allele += allele2_cnt;
            }
            while lext <= distance_to_splicing_site {
                if lidx < 0 {
                    break;
                }
                let lbf = &pileup[lidx as usize];
                INS_cnts.push(lbf.ni);
                if lbf.i {
                    // ignore insertion region
                    lidx -= 1;
                    continue;
                }
                N_cnts.push(lbf.n);
                lidx -= 1;
                lext += 1;
            }
            INS_cnts.reverse();
            N_cnts.reverse();
            while rext <= distance_to_splicing_site {
                if ridx >= pileup.len() {
                    break;
                }
                let rbf = &pileup[ridx];
                INS_cnts.push(rbf.ni);
                if rbf.i {
                    // ignore insertion region
                    ridx += 1;
                    continue;
                }
                N_cnts.push(rbf.n);
                ridx += 1;
                rext += 1;
            }
            // println!("{} N_cnts: \n{:?}", position, N_cnts);
            // println!("{} INS_cnts: \n{:?}", position, INS_cnts);

            let max_N_cnt = N_cnts.iter().max().unwrap().clone();
            let min_N_cnt = N_cnts.iter().min().unwrap().clone();
            if local_misalignment_ratio.len() as u32 > window_size {
                let sum_local_misalignment_ratio = local_misalignment_ratio.iter().sum::<f32>();
                let mut high_error_cnt = 0;
                for error_rate in local_misalignment_ratio.iter() {
                    // if mismatch rate > 15%, we think it is high error rate site for ONT reads
                    if *error_rate > 0.10 {
                        high_error_cnt += 1;
                    }
                }
                if high_error_cnt as f32 / local_misalignment_ratio.len() as f32 >= 0.5 || sum_local_misalignment_ratio / local_misalignment_ratio.len() as f32 > 0.20 {
                    if N_cnts.len() > 0 {
                        if max_N_cnt - min_N_cnt > (num_variant_allele as f32 * 0.8) as u32 {
                            println!("close to the splicing site: {}:{}, {}, {}", profile.region.chr, position, max_N_cnt, min_N_cnt);
                            println!("{}:{} high local error rate: {}, {}", profile.region.chr, position, high_error_cnt as f32 / local_misalignment_ratio.len() as f32, sum_local_misalignment_ratio / local_misalignment_ratio.len() as f32);
                            position += 1;
                            continue;
                        }
                    }
                }
            }

            // filtering insertion cause false variant near splicing site
            let mut insertion_concerned = false;
            let mut ins_cnt = 0;
            for INS_cnt in INS_cnts.iter() {
                if (*INS_cnt > (num_variant_allele as f32 * 0.8) as u32) && (max_N_cnt - min_N_cnt > (num_variant_allele as f32 * 0.8) as u32) {
                    insertion_concerned = true;
                    ins_cnt = *INS_cnt;
                }
            }
            if insertion_concerned {
                println!("insertion cause false variant: {}:{}, {}, {}", profile.region.chr, position, num_variant_allele, ins_cnt);
                position += 1;
                continue;
            }


            /*let mut local_misalignment_ratio: Vec<f32> = Vec::new();
            let mut lext = 1;
            let mut lidx = bfidx as i32 - 1;
            let mut rext = 1;
            let mut ridx = bfidx + 1;
            while lext <= distance_to_splicing_site {
                if lidx < 0 {
                    break;
                }
                let lbf = &pileup[lidx as usize];
                if lbf.i {
                    // ignore insertion region
                    lidx -= 1;
                    continue;
                }
                let local_error_rate = lbf.get_none_ref_count() as f32 / (lbf.a + lbf.c + lbf.g + lbf.t + lbf.d) as f32;
                local_misalignment_ratio.push(local_error_rate);
                lext += 1;
                lidx -= 1;
            }
            while rext <= distance_to_splicing_site {
                if ridx >= pileup.len() {
                    break;
                }
                let rbf = &pileup[ridx];
                if rbf.i {
                    // ignore insertion region
                    ridx += 1;
                    continue;
                }
                let local_error_rate = rbf.get_none_ref_count() as f32 / (rbf.a + rbf.c + rbf.g + rbf.t + rbf.d) as f32;
                local_misalignment_ratio.push(local_error_rate);
                rext += 1;
                ridx += 1;
            }
            println!("{} local_misalignment_ratio: \n{:?}", position, local_misalignment_ratio);
            if local_misalignment_ratio.len() as u32 > distance_to_splicing_site {
                let sum_local_misalignment_ratio = local_misalignment_ratio.iter().sum::<f32>();
                let mut high_error_cnt = 0;
                for error_rate in local_misalignment_ratio.iter() {
                    // if mismatch rate > 10%, we think it is high error rate site for ONT reads
                    if *error_rate > 0.10 {
                        high_error_cnt += 1;
                    }
                }
                if high_error_cnt as f32 / local_misalignment_ratio.len() as f32 > 0.5 || sum_local_misalignment_ratio / local_misalignment_ratio.len() as f32 > 0.10 {
                    println!("{} high local error rate: {}, {}", position, high_error_cnt as f32 / local_misalignment_ratio.len() as f32, sum_local_misalignment_ratio / local_misalignment_ratio.len() as f32);
                    position += 1;
                    continue;
                }
            }*/


            /*// 6.filtering close to the donor and acceptor site
            let mut N_cnts: Vec<u32> = Vec::new();  // number of N bases (intron)
            let mut B_cnts: Vec<u32> = Vec::new();  // number of nucleotide bases
            let mut lext = 0;
            let mut lidx = bfidx as i32;
            let mut rext = 1;
            let mut ridx = bfidx + 1;
            let mut num_variant_allele = 0;
            if allele1 != bf.ref_base {
                num_variant_allele += allele1_cnt;
            }
            if allele2 != bf.ref_base {
                num_variant_allele += allele2_cnt;
            }
            while lext <= distance_to_splicing_site {
                if lidx < 0 {
                    break;
                }
                let lbf = &pileup[lidx as usize];
                if lbf.i {
                    // ignore insertion region
                    lidx -= 1;
                    continue;
                }
                N_cnts.push(lbf.n);
                B_cnts.push(lbf.a + lbf.c + lbf.g + lbf.t + lbf.d);
                lidx -= 1;
                lext += 1;
            }
            N_cnts.reverse();
            B_cnts.reverse();
            while rext <= distance_to_splicing_site {
                if ridx >= pileup.len() {
                    break;
                }
                let rbf = &pileup[ridx];
                if rbf.i {
                    // ignore insertion region
                    ridx += 1;
                    continue;
                }
                N_cnts.push(rbf.n);
                B_cnts.push(rbf.a + rbf.c + rbf.g + rbf.t + rbf.d);
                ridx += 1;
                rext += 1;
            }
            println!("{} N_cnts: \n{:?}", position, N_cnts);
            println!("{} B_cnts: \n{:?}", position, B_cnts);
            if N_cnts.len() > 0 {
                let max_N_cnt = N_cnts.iter().max().unwrap().clone();
                let min_N_cnt = N_cnts.iter().min().unwrap().clone();
                if max_N_cnt - min_N_cnt > (num_variant_allele as f32 * 0.8) as u32 {
                    println!("close to the donor and acceptor site: {}, {}, {}", position, max_N_cnt, min_N_cnt);
                    position += 1;
                    continue;
                }
                // if min_N_cnt == 0 {
                //     if max_N_cnt - min_N_cnt > 10 {
                //         println!("close to the donor and acceptor site: {}, {}, {}", position, max_N_cnt, min_N_cnt);
                //         position += 1;
                //         continue;
                //     }
                // } else {
                //     // if (max_N_cnt - min_N_cnt) as f32 / (min_N_cnt as f32) > 0.5 {
                //     if (max_N_cnt - min_N_cnt) as f32 / (min_N_cnt as f32) > 0.3 && max_N_cnt - min_N_cnt >= 8 {
                //         println!("close to the donor and acceptor site: {}, {}, {}", position, max_N_cnt, min_N_cnt);
                //         position += 1;
                //         continue;
                //     }
                // }
            }*/


            // 8. filtering lying in long homopolymer regions


            // genotype likelihood
            let mut loglikelihood = [0.0, 0.0, 0.0];
            let mut logprob = [0.0, 0.0, 0.0];
            let theta = 0.001;  // mutation rate
            let background_prob = [theta / 2.0, theta, 1.0 - 1.5 * theta];    // background of probability of observe homo variant, hete variant and homo reference
            let mut identical_baseqs;
            let mut different_baseqs;
            if bf.ref_base == 'A' {
                identical_baseqs = &bf.baseq.a;
                different_baseqs = [&bf.baseq.c, &bf.baseq.g, &bf.baseq.t];
            } else if bf.ref_base == 'C' {
                identical_baseqs = &bf.baseq.c;
                different_baseqs = [&bf.baseq.a, &bf.baseq.g, &bf.baseq.t];
            } else if bf.ref_base == 'G' {
                identical_baseqs = &bf.baseq.g;
                different_baseqs = [&bf.baseq.a, &bf.baseq.c, &bf.baseq.t];
            } else if bf.ref_base == 'T' {
                identical_baseqs = &bf.baseq.t;
                different_baseqs = [&bf.baseq.a, &bf.baseq.c, &bf.baseq.g];
            } else if bf.ref_base == 'N' {
                println!("{}:{} ref base {}", profile.region.chr, position, bf.ref_base);
                position += 1;
                continue;
            } else {
                println!("{}:{} unknown ref base {}", profile.region.chr, position, bf.ref_base);
                position += 1;
                continue;
            }

            for bq in identical_baseqs.iter() {
                let error_rate = 0.1_f64.powf((*bq as f64) / 10.0);
                loglikelihood[0] += error_rate.log10();
                loglikelihood[2] += (1.0 - error_rate).log10();
            }

            for bq_vec in different_baseqs.iter() {
                for bq in bq_vec.iter() {
                    let error_rate = 0.1_f64.powf((*bq as f64) / 10.0);
                    loglikelihood[0] += (1.0 - error_rate).log10();
                    loglikelihood[2] += error_rate.log10();
                }
            }

            let num_reads = bf.a + bf.c + bf.g + bf.t;
            loglikelihood[1] -= (num_reads as f64) * 2.0_f64.log10();   // example: logL(0) = -1, logL(1) = -6, logL(2) = -26

            // PL: phred-scaled likelihood
            // https://gatk.broadinstitute.org/hc/en-us/articles/360035890451-Calculation-of-PL-and-GQ-by-HaplotypeCaller-and-GenotypeGVCFs
            // https://gatk.broadinstitute.org/hc/en-us/articles/360035890511
            logprob = loglikelihood.clone();
            // multiple prior probability of genotype
            logprob[0] += (background_prob[0] as f64).log10();
            logprob[1] += (background_prob[1] as f64).log10();
            logprob[2] += (background_prob[2] as f64).log10();
            let max_logprob = logprob[0].max(logprob[1]).max(logprob[2]);
            // println!("1:{}:{},{:?}", position, max_logprob, logprob);
            logprob[0] = logprob[0] - max_logprob;
            logprob[1] = logprob[1] - max_logprob;
            logprob[2] = logprob[2] - max_logprob;
            // println!("2:{}:{},{:?}", position, max_logprob, logprob);
            let mut genotype_prob = logprob.clone();
            genotype_prob[0] = 10.0_f64.powf(logprob[0]);
            genotype_prob[1] = 10.0_f64.powf(logprob[1]);
            genotype_prob[2] = 10.0_f64.powf(logprob[2]);
            let sum_genotype_prob = genotype_prob[0] + genotype_prob[1] + genotype_prob[2];
            // println!("3:{}:{},{:?}", position, sum_genotype_prob, genotype_prob);
            genotype_prob = [genotype_prob[0] / sum_genotype_prob, genotype_prob[1] / sum_genotype_prob, genotype_prob[2] / sum_genotype_prob];
            // println!("4:{}:{},{:?}", position, correction_factor, genotype_prob);
            // QUAL phred-scaled quality score for the assertion made in ALT. i.e. give -10log_10 prob(call in ALT is wrong).
            // If ALT is `.` (no variant) then this is -10log_10 p(variant), and if ALT is not `.` this is -10log_10 p(no variant).
            let variant_quality = -10.0 * ((10e-301_f64.max(genotype_prob[2])).log10());    // if variant_quality is greater than 3000, we set it to 3000

            // calculate GQ: The value of GQ is simply the difference between the second lowest PL and the lowest PL (which is always 0, normalized PL)
            let mut log10_likelihood = loglikelihood.clone();
            let max_log10_likelihood = log10_likelihood[0].max(log10_likelihood[1]).max(log10_likelihood[2]);
            log10_likelihood[0] = 10.0_f64.powf(log10_likelihood[0] - max_log10_likelihood);
            log10_likelihood[1] = 10.0_f64.powf(log10_likelihood[1] - max_log10_likelihood);
            log10_likelihood[2] = 10.0_f64.powf(log10_likelihood[2] - max_log10_likelihood);
            let sum_log10_likelihood = log10_likelihood[0] + log10_likelihood[1] + log10_likelihood[2];
            let mut sorted_pl = [log10_likelihood[0] / sum_log10_likelihood, log10_likelihood[1] / sum_log10_likelihood, log10_likelihood[2] / sum_log10_likelihood];
            sorted_pl[0] = -10.0 * sorted_pl[0].log10();    // phred scale likelihood of genotype: 1/1
            sorted_pl[1] = -10.0 * sorted_pl[1].log10();    // phred scale likelihood of genotype: 0/1
            sorted_pl[2] = -10.0 * sorted_pl[2].log10();    // phred scale likelihood of genotype: 0/0
            sorted_pl.sort_by(cmp_f64);
            let genotype_quality = sorted_pl[1] - sorted_pl[0];

            // filter low variant quality
            if variant_quality < min_qual_for_candidate as f64 {
                println!("{}:{} low variant quality", profile.region.chr, position);
                position += 1;
                continue;
            }

            if genotype_prob[0] > genotype_prob[1] && genotype_prob[0] > genotype_prob[2] {
                // candidate homozygous SNP
                let allele1_freq = (allele1_cnt as f32) / (depth as f32);
                let allele2_freq = (allele2_cnt as f32) / (depth as f32);
                let mut candidate_snp = CandidateSNP::default();
                candidate_snp.chromosome = profile.region.chr.clone().into_bytes();
                candidate_snp.pos = position as i64;
                candidate_snp.alleles = [allele1, allele2];
                candidate_snp.allele_freqs = [allele1_freq, allele2_freq];
                candidate_snp.reference = bf.ref_base;
                candidate_snp.depth = depth;

                if allele1 != bf.ref_base && allele2 != bf.ref_base && allele1_freq < min_homozygous_freq && allele2_freq > 0.0 {
                    candidate_snp.variant_type = 3; // triallelic SNP, triallelic SNP is also considered as homozygous SNP
                } else {
                    candidate_snp.variant_type = 2; // homozygous SNP
                }


                candidate_snp.variant_quality = variant_quality;
                candidate_snp.genotype_probability = genotype_prob.clone();
                candidate_snp.genotype_quality = genotype_quality;
                if candidate_snp.variant_type == 2 {
                    if bf.ref_base == 'A' && allele1 == 'G' {
                        candidate_snp.rna_editing = true;
                        candidate_snp.filter = true;
                    }
                    if bf.ref_base == 'T' && allele1 == 'C' {
                        candidate_snp.rna_editing = true;
                        candidate_snp.filter = true;
                    }
                } else if candidate_snp.variant_type == 3 {
                    if bf.ref_base == 'A' && (allele1 == 'G' || allele2 == 'G') {
                        candidate_snp.rna_editing = true;
                        candidate_snp.filter = true;
                    }
                    if bf.ref_base == 'T' && (allele1 == 'C' || allele2 == 'C') {
                        candidate_snp.rna_editing = true;
                        candidate_snp.filter = true;
                    }
                }
                // println!("homo genotype quality: {:?},{:?}", likelihood, candidate_snp.genotype_quality);
                self.candidate_snps.push(candidate_snp);
                self.homo_snps.push(self.candidate_snps.len() - 1);
            } else if genotype_prob[1] > genotype_prob[0] && genotype_prob[1] > genotype_prob[2] {
                // candidate heterozygous SNP
                let allele1_freq = (allele1_cnt as f32) / (depth as f32);
                let allele2_freq = (allele2_cnt as f32) / (depth as f32);
                if allele1 != bf.ref_base && allele1_freq < min_allele_freq {
                    println!("{}:{} allele1 {:?} low frequency {}", profile.region.chr, position, allele1, allele1_freq);
                    position += 1;
                    continue;
                } else if allele2 != bf.ref_base && allele2_freq < min_allele_freq {
                    println!("{}:{} allele2 {:?} low frequency {}", profile.region.chr, position, allele2, allele2_freq);
                    position += 1;
                    continue;
                }
                let mut candidate_snp = CandidateSNP::default();
                candidate_snp.chromosome = profile.region.chr.clone().into_bytes();
                candidate_snp.pos = position as i64;
                candidate_snp.alleles = [allele1, allele2];
                candidate_snp.allele_freqs = [allele1_freq, allele2_freq];
                candidate_snp.reference = bf.ref_base;
                candidate_snp.depth = depth;
                candidate_snp.variant_type = 1; // heterozygous SNP
                candidate_snp.variant_quality = variant_quality;
                candidate_snp.genotype_probability = genotype_prob.clone();
                candidate_snp.genotype_quality = genotype_quality;
                if bf.ref_base == 'A' {
                    if (allele1 == 'A' && allele2 == 'G') || (allele1 == 'G' && allele2 == 'A') {
                        candidate_snp.rna_editing = true;
                        candidate_snp.filter = true;
                    }
                }
                if bf.ref_base == 'T' {
                    if (allele1 == 'T' && allele2 == 'C') || (allele1 == 'C' && allele2 == 'T') {
                        candidate_snp.rna_editing = true;
                        candidate_snp.filter = true;
                    }
                }
                // println!("hete genotype quality: {:?},{:?}", likelihood, candidate_snp.genotype_quality);
                self.candidate_snps.push(candidate_snp);
                self.hete_snps.push(self.candidate_snps.len() - 1);
            }
            position += 1;
        }

        // filter dense SNPs
        for i in 0..self.candidate_snps.len() {
            for j in i..self.candidate_snps.len() {
                if self.candidate_snps[j].pos - self.candidate_snps[i].pos > dense_win_size as i64 {
                    if (j - 1 - i + 1) as u32 >= min_dense_cnt && ((self.candidate_snps[j - 1].pos - self.candidate_snps[i].pos + 1) as f32) / ((j - 1 - i + 1) as f32) <= avg_dense_dist {
                        for tk in i..j {
                            // println!("dense SNPs: {}", self.candidate_snps[tk].pos);
                            // even rna editing may be filtered by dense SNPs
                            self.candidate_snps[tk].rna_editing = false;
                            self.candidate_snps[tk].filter = true;
                        }
                    }
                    break;
                }
            }
        }

        // update hete_snps, remove filtered SNPs
        let mut tmp_hetes: Vec<usize> = Vec::new();
        for i in self.hete_snps.iter() {
            if !self.candidate_snps[*i].filter {
                tmp_hetes.push(*i);
            }
        }
        self.hete_snps = tmp_hetes;
    }

    pub unsafe fn init_haplotypes(&mut self) {
        // initialize haplotype of heterozygous snp
        for i in self.hete_snps.iter() {
            if drand48() < 0.5 {
                self.candidate_snps[*i].haplotype = -1;
            } else {
                self.candidate_snps[*i].haplotype = 1;
            }
        }
    }

    pub unsafe fn init_assignment(&mut self) {
        let mut fragment_cnt: HashMap<usize, usize> = HashMap::new();
        let mut covered_fragments: HashSet<usize> = HashSet::new();

        for i in self.hete_snps.iter() {
            let snp = &self.candidate_snps[*i];
            for k in snp.snp_cover_fragments.iter() {
                if fragment_cnt.contains_key(k) {
                    let cnt = fragment_cnt.get_mut(k).unwrap();
                    *cnt += 1;
                } else {
                    fragment_cnt.insert(*k, 1);
                }
            }
        }

        for (k, v) in fragment_cnt.iter() {
            if *v >= 2 {
                covered_fragments.insert(*k);
            }
        }

        for k in covered_fragments.iter() {
            if drand48() < 0.5 {
                self.fragments[*k].haplotag = -1;
            } else {
                self.fragments[*k].haplotag = 1;
            }
        }
    }

    pub fn get_fragments(&mut self, bam_path: &str, region: &Region) {
        let mut bam_reader: bam::IndexedReader = bam::IndexedReader::from_path(bam_path).unwrap();
        bam_reader.fetch((region.chr.as_str(), region.start, region.end)).unwrap();
        let mut record = Record::new();
        if self.hete_snps.len() == 0 {
            return;
        }
        while let Some(result) = bam_reader.read(&mut record) {
            if result.is_err() {
                panic!("BAM parsing failed...");
            }
            // TODO: filtering unmapped, secondary, supplementary reads?
            if record.is_unmapped() || record.is_secondary() || record.is_supplementary() {
                continue;
            }
            let pos = record.pos(); // 0-based
            if pos > self.candidate_snps[*self.hete_snps.last().unwrap()].pos { continue; }
            let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
            let cigar = record.cigar();
            let seq = record.seq().as_bytes();
            let strand = if record.strand() == Forward { 0 } else { 1 };
            let mut pos_on_ref = pos;   // 0-based
            let mut pos_on_query = cigar.leading_softclips();   // 0-based
            let mut idx = 0; // index in self.hete_snps
            let mut snp_pos = -1;   // pre-computed position of candidate SNPs
            let mut alleles;   // pre-computed alleles of candidate SNPs
            if pos <= self.candidate_snps[*self.hete_snps.first().unwrap()].pos {
                snp_pos = self.candidate_snps[self.hete_snps[idx]].pos;
                alleles = self.candidate_snps[self.hete_snps[idx]].alleles.clone();
            } else {
                // find the first SNP in the read
                while idx < self.hete_snps.len() {
                    if self.candidate_snps[self.hete_snps[idx]].pos >= pos {
                        break;
                    }
                    idx += 1;
                }
                assert!(idx < self.hete_snps.len(), "Error: idx < self.candidate_snps.len()");
                snp_pos = self.candidate_snps[self.hete_snps[idx]].pos;
                alleles = self.candidate_snps[self.hete_snps[idx]].alleles.clone();
            }

            let mut fragment = Fragment::default();
            fragment.read_id = qname.clone();
            fragment.fragment_idx = self.fragments.len();

            let mut exon_start = -1;
            let mut exon_end = -1;
            exon_start = pos_on_ref;
            exon_end = exon_start;

            for cg in cigar.iter() {
                // since need to parse whole cigar to obtain exons, can not break before the end of cigar
                // if pos_on_ref > self.candidate_snps[*self.hete_snps.last().unwrap()].pos {
                //     break;
                // }
                // if idx >= self.hete_snps.len() {
                //     break;
                // }
                match cg.char() as u8 {
                    b'S' | b'H' => {
                        continue;
                    }
                    b'M' | b'X' | b'=' => {
                        for _ in 0..cg.len() {
                            // assert!(pos_on_ref <= snp_pos, "Error: pos_on_ref <= snp_pos");
                            if pos_on_ref == snp_pos {
                                let mut frag_elem = FragElem::default();
                                frag_elem.snp_idx = self.hete_snps[idx];
                                frag_elem.pos = pos_on_ref;
                                frag_elem.base = seq[pos_on_query as usize] as char;
                                frag_elem.baseq = record.qual()[pos_on_query as usize];
                                frag_elem.strand = strand;
                                frag_elem.prob = 10.0_f64.powf(-(frag_elem.baseq as f64) / 10.0);
                                if frag_elem.base == alleles[0] {
                                    frag_elem.p = 1;    // reference allele
                                } else if frag_elem.base == alleles[1] {
                                    frag_elem.p = -1;   // alternate allele
                                } else {
                                    frag_elem.p = 0;    // not covered
                                }
                                // calculate the cost of links of two snps as hapcut (optioanl)
                                // if fragment.list.len() > 0 {
                                //     for prev_frag_elem in fragment.list.iter() {
                                //         if prev_frag_elem.p == 0 || frag_elem.p == 0 {
                                //             continue;
                                //         }
                                //         let q1 = 0.1_f64.powf((prev_frag_elem.baseq as f64) / 10.0); // probability of error for prev_frag_elem
                                //         let q2 = 0.1_f64.powf((frag_elem.baseq as f64) / 10.0);  // probability of error for frag_elem
                                //         let epsilon = q1 * (1.0 - q2) + (1.0 - q1) * q2; // probability of sequencing error only occurs on start node or end node.
                                //         let w = (prev_frag_elem.p as f64) * (frag_elem.p as f64) * f64::log10((1.0 - epsilon) / epsilon);
                                //         // println!("q1:{}, q2:{}, epsilon:{}, w:{}", q1, q2, epsilon, w);
                                //         if self.edges.contains_key(&[prev_frag_elem.snp_idx, frag_elem.snp_idx]) {
                                //             let edge = self.edges.get_mut(&[prev_frag_elem.snp_idx, frag_elem.snp_idx]).unwrap();
                                //             edge.frag_idxes.push(fragment.fragment_idx);
                                //             edge.w += w;
                                //         } else {
                                //             let mut edge = Edge::default();
                                //             edge.snp_idxes = [prev_frag_elem.snp_idx, frag_elem.snp_idx];
                                //             edge.snp_poses = [prev_frag_elem.pos, frag_elem.pos];
                                //             edge.w = w;
                                //             edge.frag_idxes.push(fragment.fragment_idx);
                                //             self.edges.insert(edge.snp_idxes, edge);
                                //         }
                                //     }
                                // }
                                // println!("M: {:?}", frag_elem);

                                // filtered SNP and rna editing site will not be used for haplotype phasing, not covered SNP will not be used for haplotype phasing
                                if self.candidate_snps[frag_elem.snp_idx].filter == false && self.candidate_snps[frag_elem.snp_idx].rna_editing == false && frag_elem.p != 0 {
                                    fragment.list.push(frag_elem);
                                }
                                idx += 1;
                                // since need to parse whole cigar to obtain exons, can not break before the end of cigar
                                // if idx >= self.hete_snps.len() {
                                //     pos_on_query += 1;
                                //     pos_on_ref += 1;
                                //     break;
                                // }
                                if idx < self.hete_snps.len() {
                                    snp_pos = self.candidate_snps[self.hete_snps[idx]].pos;
                                    alleles = self.candidate_snps[self.hete_snps[idx]].alleles.clone();
                                }
                            }
                            pos_on_query += 1;
                            pos_on_ref += 1;
                        }
                    }
                    b'I' => {
                        pos_on_query += cg.len() as i64;
                    }
                    b'D' => {
                        for _ in 0..cg.len() {
                            // assert!(pos_on_ref <= snp_pos, "Error: pos_on_ref <= snp_pos");
                            if pos_on_ref == snp_pos {
                                let mut frag_elem = FragElem::default();
                                frag_elem.snp_idx = self.hete_snps[idx];
                                frag_elem.pos = pos_on_ref;
                                frag_elem.base = '-';
                                frag_elem.baseq = 0;
                                frag_elem.strand = strand;
                                frag_elem.p = 0;
                                // calculate the cost of links of two snps as hapcut (optioanl)
                                // if fragment.list.len() > 0 {
                                //     for prev_frag_elem in fragment.list.iter() {
                                //         if prev_frag_elem.p == 0 || frag_elem.p == 0 {
                                //             continue;
                                //         }
                                //         let q1 = 0.1_f64.powf((prev_frag_elem.baseq as f64) / 10.0); // probability of error for prev_frag_elem
                                //         let q2 = 0.1_f64.powf((frag_elem.baseq as f64) / 10.0);  // probability of error for frag_elem
                                //         let epsilon = q1 * (1.0 - q2) + (1.0 - q1) * q2; // probability of sequencing error only occurs on start node or end node.
                                //         let w = (prev_frag_elem.p as f64) * (frag_elem.p as f64) * f64::log10((1.0 - epsilon) / epsilon);
                                //         // println!("q1:{}, q2:{}, epsilon:{}, w:{}", q1, q2, epsilon, w);
                                //         if self.edges.contains_key(&[prev_frag_elem.snp_idx, frag_elem.snp_idx]) {
                                //             let edge = self.edges.get_mut(&[prev_frag_elem.snp_idx, frag_elem.snp_idx]).unwrap();
                                //             edge.frag_idxes.push(fragment.fragment_idx);
                                //             edge.w += w;
                                //         } else {
                                //             let mut edge = Edge::default();
                                //             edge.snp_idxes = [prev_frag_elem.snp_idx, frag_elem.snp_idx];
                                //             edge.snp_poses = [prev_frag_elem.pos, frag_elem.pos];
                                //             edge.w = w;
                                //             edge.frag_idxes.push(fragment.fragment_idx);
                                //             self.edges.insert(edge.snp_idxes, edge);
                                //         }
                                //     }
                                // }
                                // println!("D: {:?}", frag_elem);

                                // Deletion will not be used for haplotype phasing
                                // filtered SNP and rna editing site will not be used for haplotype phasing
                                // if self.candidate_snps[frag_elem.snp_idx].filter == false && self.candidate_snps[frag_elem.snp_idx].rna_editing == false {
                                //     fragment.list.push(frag_elem);
                                // }

                                idx += 1;
                                // since need to parse whole cigar to obtain exons, can not break before the end of cigar
                                // if idx >= self.hete_snps.len() {
                                //     pos_on_ref += 1;
                                //     break;
                                // }
                                if idx < self.hete_snps.len() {
                                    snp_pos = self.candidate_snps[self.hete_snps[idx]].pos;
                                    alleles = self.candidate_snps[self.hete_snps[idx]].alleles.clone();
                                }
                            }
                            pos_on_ref += 1;
                        }
                    }
                    b'N' => {
                        exon_end = pos_on_ref;
                        fragment.exons.push([exon_start, exon_end]);
                        exon_start = -1;
                        exon_end = -1;
                        for _ in 0..cg.len() {
                            // assert!(pos_on_ref <= snp_pos, "Error: pos_on_ref <= snp_pos");
                            if pos_on_ref == snp_pos {
                                let mut frag_elem = FragElem::default();
                                frag_elem.snp_idx = self.hete_snps[idx];
                                frag_elem.pos = pos_on_ref;
                                frag_elem.base = b'-' as char;
                                frag_elem.baseq = 0;
                                frag_elem.strand = strand;
                                frag_elem.p = 0;
                                // calculate the cost of links of two snps as hapcut (optioanl)
                                // if fragment.list.len() > 0 {
                                //     for prev_frag_elem in fragment.list.iter() {
                                //         if prev_frag_elem.p == 0 || frag_elem.p == 0 {
                                //             continue;
                                //         }
                                //         let q1 = 0.1_f64.powf((prev_frag_elem.baseq as f64) / 10.0); // probability of error for prev_frag_elem
                                //         let q2 = 0.1_f64.powf((frag_elem.baseq as f64) / 10.0);  // probability of error for frag_elem
                                //         let epsilon = q1 * (1.0 - q2) + (1.0 - q1) * q2; // probability of sequencing error only occurs on start node or end node.
                                //         let w = (prev_frag_elem.p as f64) * (frag_elem.p as f64) * f64::log10((1.0 - epsilon) / epsilon);
                                //         // println!("q1:{}, q2:{}, epsilon:{}, w:{}", q1, q2, epsilon, w);
                                //         if self.edges.contains_key(&[prev_frag_elem.snp_idx, frag_elem.snp_idx]) {
                                //             let edge = self.edges.get_mut(&[prev_frag_elem.snp_idx, frag_elem.snp_idx]).unwrap();
                                //             edge.frag_idxes.push(fragment.fragment_idx);
                                //             edge.w += w;
                                //         } else {
                                //             let mut edge = Edge::default();
                                //             edge.snp_idxes = [prev_frag_elem.snp_idx, frag_elem.snp_idx];
                                //             edge.snp_poses = [prev_frag_elem.pos, frag_elem.pos];
                                //             edge.w = w;
                                //             edge.frag_idxes.push(fragment.fragment_idx);
                                //             self.edges.insert(edge.snp_idxes, edge);
                                //         }
                                //     }
                                // }
                                // println!("N: {:?}", frag_elem);

                                // Intron will not be used for haplotype phasing
                                // filtered SNP and rna editing site will not be used for haplotype phasing
                                // if self.candidate_snps[frag_elem.snp_idx].filter == false && self.candidate_snps[frag_elem.snp_idx].rna_editing == false {
                                //     fragment.list.push(frag_elem);
                                // }

                                idx += 1;
                                // since need to parse whole cigar to obtain exons, can not break before the end of cigar
                                // if idx >= self.hete_snps.len() {
                                //     pos_on_ref += 1;
                                //     break;
                                // }
                                if idx < self.hete_snps.len() {
                                    snp_pos = self.candidate_snps[self.hete_snps[idx]].pos;
                                    alleles = self.candidate_snps[self.hete_snps[idx]].alleles.clone();
                                }
                            }
                            pos_on_ref += 1;
                        }
                        exon_start = pos_on_ref;
                        exon_end = exon_start;
                    }
                    _ => {
                        panic!("Error: unknown cigar operation: {}", cg.char());
                    }
                }
            }
            // the whole read is a single exon, no intron.
            if exon_start != -1 && exon_end != -1 && pos_on_ref > exon_end {
                exon_end = pos_on_ref;
            }
            if exon_end != exon_start {
                fragment.exons.push([exon_start, exon_end]);
            }
            exon_start = -1;
            exon_end = -1;

            // filter fragment with no heterozygous links (deletion, intron, not reference allele and alternate allele do not count as heterozygous links)
            let mut link_hete_cnt = 0;
            for fe in fragment.list.iter() {
                if fe.p != 0 {
                    link_hete_cnt += 1;
                }
            }
            if link_hete_cnt >= 2 {
                for fe in fragment.list.iter() {
                    // record each snp cover by which fragments
                    self.candidate_snps[fe.snp_idx].snp_cover_fragments.push(fragment.fragment_idx);
                }
                // println!("Fragment: {:?}", fragment);
                self.fragments.push(fragment);
            }
        }
    }

    // pub fn filter_fp_snps(&mut self, strand_bias_threshold: f32, phase_score_cutoff: Option<f64>) {
    //     /*
    //     *? Filter false positive SNPs, which may be cause by strand bias or in a dense cluster of variants (100bp length has over 3 snps).
    //      */
    //
    //
    //     // 1. filter dense cluster of variants
    //     let mut homo_hete_snps: Vec<i64> = Vec::new();
    //     for i in self.hete_snps.iter() {
    //         homo_hete_snps.push(self.candidate_snps[*i].pos);
    //     }
    //     for i in self.homo_snps.iter() {
    //         homo_hete_snps.push(self.candidate_snps[*i].pos);
    //     }
    //     // sort homo_hete_snps in ascending order
    //     homo_hete_snps.sort();
    //
    //     let mut filter_window: HashSet<i64> = HashSet::new();   // record the SNP position in a dense cluster of variants
    //     for i in 0..homo_hete_snps.len() {
    //         for j in i..homo_hete_snps.len() {
    //             // TODO: dense cluster of snps?
    //             if homo_hete_snps[j] - homo_hete_snps[i] > 300 {
    //                 // if (j - 1) - i + 1 >= 9 {
    //                 //     for tk in i..j {
    //                 //         filter_window.insert(homo_hete_snps[tk]);
    //                 //     }
    //                 // }
    //                 if (j - 1 - i + 1) >= 5 && ((homo_hete_snps[j - 1] - homo_hete_snps[i] + 1) as f32) / ((j - 1 - i + 1) as f32) <= 66.66 {
    //                     for tk in i..j {
    //                         filter_window.insert(homo_hete_snps[tk]);
    //                     }
    //                 }
    //                 break;
    //             }
    //         }
    //     }
    //
    //     // println!("dense cluster filter: {:?}", filter_window);
    //
    //     if filter_window.len() > 0 {
    //         // filter homo snps in a dense cluster of variants
    //         let mut filtered_homo_snps: Vec<usize> = Vec::new();
    //         let mut filtered_hete_snps: Vec<usize> = Vec::new();
    //
    //         for i in self.homo_snps.iter() {
    //             if !filter_window.contains(&self.candidate_snps[*i].pos) {
    //                 filtered_homo_snps.push(*i);
    //             } else {
    //                 self.candidate_snps[*i].filter = true;
    //             }
    //         }
    //         self.homo_snps = filtered_homo_snps;
    //
    //         for i in self.hete_snps.iter() {
    //             if !filter_window.contains(&self.candidate_snps[*i].pos) {
    //                 filtered_hete_snps.push(*i);
    //             } else {
    //                 self.candidate_snps[*i].filter = true;
    //             }
    //         }
    //         self.hete_snps = filtered_hete_snps;
    //     }
    //
    //     // // 2. filtering strand bias
    //     // let mut filter_window: HashSet<usize> = HashSet::new();
    //     // for i in 0..self.candidate_snps.len() {
    //     //     let snp = &mut self.candidate_snps[i];
    //     //     if snp.filter == true { continue; }
    //     //     let mut variant_strand_cnt = [0, 0];    // variant strand [forward, reverse]
    //     //     for k in snp.snp_cover_fragments.iter() {
    //     //         for fe in self.fragments[*k].list.iter() {
    //     //             if fe.snp_idx == i {
    //     //                 if fe.p != 0 {
    //     //                     if fe.base != snp.reference {
    //     //                         if fe.strand == 0 {
    //     //                             variant_strand_cnt[0] += 1;
    //     //                         } else {
    //     //                             variant_strand_cnt[1] += 1;
    //     //                         }
    //     //                     }
    //     //                 }
    //     //             }
    //     //         }
    //     //     }
    //     //     let total_variant_cnt = variant_strand_cnt[0] + variant_strand_cnt[1];
    //     //     println!("{:?} strand bias: {}, {}, {}", snp.pos, variant_strand_cnt[0], variant_strand_cnt[1], total_variant_cnt);
    //     //     if total_variant_cnt > 0 {
    //     //         // variant strand bias
    //     //         if variant_strand_cnt[0] as f32 / total_variant_cnt as f32 >= strand_bias_threshold || variant_strand_cnt[1] as f32 / total_variant_cnt as f32 >= strand_bias_threshold {
    //     //             snp.filter = true;  // current snp is filtered because of strand bias
    //     //             println!("strand bias: {}, {}, {}, {}", snp.pos, variant_strand_cnt[0], variant_strand_cnt[1], total_variant_cnt);
    //     //             continue;
    //     //         }
    //     //     }
    //     //     filter_window.insert(i);   // index of snp
    //     // }
    //     // if filter_window.len() > 0 {
    //     //     let mut filtered_homo_snps: Vec<usize> = Vec::new();
    //     //     let mut filtered_hete_snps: Vec<usize> = Vec::new();
    //     //
    //     //     for i in self.homo_snps.iter() {
    //     //         if !filter_window.contains(i) {
    //     //             filtered_homo_snps.push(*i);
    //     //         }
    //     //     }
    //     //     self.homo_snps = filtered_homo_snps;
    //     //
    //     //     for i in self.hete_snps.iter() {
    //     //         if !filter_window.contains(i) {
    //     //             filtered_hete_snps.push(*i);
    //     //         }
    //     //     }
    //     //     self.hete_snps = filtered_hete_snps;
    //     // }
    // }

    pub fn keep_reliable_snps_in_component(&mut self) {
        let mut gqmap: HashMap<usize, f64> = HashMap::new();
        for i in self.hete_snps.iter() {
            let snp = &mut self.candidate_snps[*i];
            if snp.filter == true { continue; }
            gqmap.insert(*i, snp.genotype_quality);
        }
        // sort gqmap in descending order
        let mut gqmap_vec: Vec<(&usize, &f64)> = gqmap.iter().collect();
        gqmap_vec.sort_by(|a, b| b.1.partial_cmp(a.1).unwrap());
        if gqmap_vec.len() > 10 {
            let mut filtered_hete_snps: Vec<usize> = Vec::new();
            for i in 0..((gqmap_vec.len() as f32 * 0.7) as usize) {
                // for i in 0..gqmap_vec.len() {
                filtered_hete_snps.push(*gqmap_vec[i].0);
            }
            self.hete_snps = filtered_hete_snps;
        }
    }

    // pub fn optimization_using_maxcut(&mut self) {
    //     // optimization using maxcut
    //     let mut MAX_ITER = 3;
    //     while MAX_ITER > 0 {
    //         MAX_ITER -= 1;
    //         println!("haplotype: {:?}", self.haplotype);
    //         let mut haplotype_weights: Vec<Edge> = Vec::new();
    //         let mut nodeset: HashSet<usize> = HashSet::new(); // all SNP nodes
    //         let mut sum_hap_wts = 0.0;
    //         for edge in self.edges.iter() {
    //             let mut te = edge.1.clone();
    //             nodeset.insert(te.snp_idxes[0]);
    //             nodeset.insert(te.snp_idxes[1]);
    //             // Haplotype weight = \delta_{i}\delta_{j}W_{ij}, \delta_{i} (\delta_{j}) is haplotype of node i (j), W_{ij} is the weight of edge ij.
    //             te.w = te.w * self.haplotype[te.snp_idxes[0]] as f64 * self.haplotype[te.snp_idxes[1]] as f64;
    //             sum_hap_wts += te.w;
    //             haplotype_weights.push(te);
    //         }
    //         if haplotype_weights.len() == 0 { break; }
    //         // Sort haplotype weights in ascending order to get the lowest value
    //         haplotype_weights.sort_by(|a, b| a.w.partial_cmp(&b.w).unwrap());
    //         // TODO: select k lowest value for calculation.
    //         if haplotype_weights.first().unwrap().w >= 0.0 {
    //             println!("All edges have positive haplotype weights. Phasing is done. Haplotype = {:?}", self.haplotype);
    //             break;
    //         }
    //         let mut S1: HashSet<usize> = HashSet::new();    // whole S1
    //         let mut S2: HashSet<usize> = HashSet::new();    // whole S2
    //         let mut s1: HashSet<usize> = HashSet::new();    // s1 in current connected component
    //         let mut s2: HashSet<usize> = HashSet::new();    // s2 in current connected component
    //         // initial select the lowest weight edge
    //         s1.insert(haplotype_weights[0].snp_idxes[0]);
    //         s2.insert(haplotype_weights[0].snp_idxes[1]);
    //
    //         while s1.len() + s2.len() < nodeset.len() {
    //             let union_s1_s2 = s1.union(&s2).cloned().collect();
    //             let mut max_abs_weight: f64 = 0.0;
    //             let mut max_weight_node: i32 = -1;
    //             // Find the vertex maximizes the haplotype weights
    //             for cand_node in nodeset.difference(&union_s1_s2) {
    //                 // candidate node is from V-(s1+s2)
    //                 let mut s1_weights: f64 = 0.0;
    //                 let mut s2_weights: f64 = 0.0;
    //
    //                 // calculate the sum of weight for all edges from cand_node to set s1
    //                 for s1_node in s1.iter() {
    //                     if self.edges.contains_key(&[*cand_node, *s1_node]) {
    //                         let tn = self.edges.get(&[*cand_node, *s1_node]).unwrap();
    //                         s1_weights += tn.w * (self.haplotype[tn.snp_idxes[0]] as f64) * (self.haplotype[tn.snp_idxes[1]] as f64);
    //                     } else if self.edges.contains_key(&[*s1_node, *cand_node]) {
    //                         let tn = self.edges.get(&[*s1_node, *cand_node]).unwrap();
    //                         s1_weights += tn.w * (self.haplotype[tn.snp_idxes[0]] as f64) * (self.haplotype[tn.snp_idxes[1]] as f64);
    //                     }
    //                 }
    //
    //                 // calculate the sum of weight for all edges from cand_node to set s2
    //                 for s2_node in s2.iter() {
    //                     if self.edges.contains_key(&[*cand_node, *s2_node]) {
    //                         let tn = self.edges.get(&[*cand_node, *s2_node]).unwrap();
    //                         s2_weights += tn.w * (self.haplotype[tn.snp_idxes[0]] as f64) * (self.haplotype[tn.snp_idxes[1]] as f64);
    //                     } else if self.edges.contains_key(&[*s2_node, *cand_node]) {
    //                         let tn = self.edges.get(&[*s2_node, *cand_node]).unwrap();
    //                         s2_weights += tn.w * (self.haplotype[tn.snp_idxes[0]] as f64) * (self.haplotype[tn.snp_idxes[1]] as f64);
    //                     }
    //                 }
    //
    //                 // println!("s1_weights = {:?}", s1_weights);
    //                 // println!("s2_weights = {:?}", s2_weights);
    //
    //                 // maximum value
    //                 if (s1_weights - s2_weights).abs() > max_abs_weight.abs() {
    //                     max_abs_weight = s1_weights - s2_weights;
    //                     max_weight_node = *cand_node as i32;
    //                 }
    //             }
    //
    //             // disconnected component
    //             if max_abs_weight.abs() == 0.0 && max_weight_node == -1 {
    //                 // remove s1 and s2 nodes from nodeset and re-init s1 and s2 for next connected component
    //                 for node in s1.union(&s2) {
    //                     nodeset.remove(node);
    //                 }
    //                 // change the weight of processed edges to INF to avoid been selected in next connected component
    //                 for edge in haplotype_weights.iter_mut() {
    //                     if !nodeset.contains(&edge.snp_idxes[0]) || !nodeset.contains(&edge.snp_idxes[1]) {
    //                         edge.w = f64::INFINITY;
    //                     }
    //                 }
    //
    //                 haplotype_weights.sort_by(|a, b| a.w.partial_cmp(&b.w).unwrap());
    //                 if haplotype_weights.first().unwrap().w >= 0.0 {
    //                     break;
    //                 }
    //
    //                 // clear s1/s2 for current connected component and add the edge with the lowest haplotype weight of next connected component to s1/s2
    //                 println!("Clear s1 and s2 for next connected component");
    //                 println!("Current component s1 = {:?}", s1);
    //                 println!("Current component s2 = {:?}", s2);
    //                 for node in s1.iter() {
    //                     S1.insert(*node);
    //                 }
    //                 for node in s2.iter() {
    //                     S2.insert(*node);
    //                 }
    //                 s1.clear();
    //                 s2.clear();
    //                 s1.insert(haplotype_weights[0].snp_idxes[0]);
    //                 s2.insert(haplotype_weights[0].snp_idxes[1]);
    //                 continue;
    //             }
    //
    //             if max_abs_weight > 0.0 {
    //                 s1.insert(max_weight_node as usize);
    //                 // println!("{:?} move to s1 {:?}", max_weight_node as usize, s1);
    //             } else {
    //                 s2.insert(max_weight_node as usize);
    //                 // println!("{:?} move to s2 {:?}", max_weight_node as usize, s2);
    //             }
    //         }
    //         println!("Current component s1 = {:?}", s1);
    //         println!("Current component s2 = {:?}", s2);
    //         for node in s1.iter() {
    //             S1.insert(*node);
    //         }
    //         for node in s2.iter() {
    //             S2.insert(*node);
    //         }
    //
    //         // check sum of haplotype weight increase
    //         let mut t_haplotype = self.haplotype.clone();
    //         for i in S1.iter() {
    //             if t_haplotype[*i] == 1 {
    //                 t_haplotype[*i] = -1;
    //             } else {
    //                 t_haplotype[*i] = 1;
    //             }
    //         }
    //         let mut t_sum_hap_wts = 0.0;
    //         for edge in self.edges.iter() {
    //             let mut te = edge.1.clone();
    //             te.w = te.w * t_haplotype[te.snp_idxes[0]] as f64 * t_haplotype[te.snp_idxes[1]] as f64;
    //             t_sum_hap_wts += te.w;
    //         }
    //         println!("latest haplotype weight: {:?}, previous haplotype weight: {:?}", t_sum_hap_wts, sum_hap_wts);
    //         if t_sum_hap_wts > sum_hap_wts {
    //             self.haplotype = t_haplotype;
    //             self.phased = true;
    //         }
    //         println!("MaxCut optimization");
    //         println!("haplotype: {:?}", self.haplotype);
    //     }
    // }

    pub fn cal_sigma_delta(sigma_k: i32, delta: &Vec<i32>, ps: &Vec<i32>, probs: &Vec<f64>) -> f64 {
        // calculate P(sigma_k | delta)
        // sigma_k: the assignment of read k, 1 or -1.
        // delta: the haplotypes of the SNPs covered by read k, each haplotype is 1 or -1.
        // ps: the allele of each base, 1,-1
        // probs: the probability of observing base at each SNP for read k, equals to 10^(-Q/10).
        let mut q1: f64 = 1.0;
        let mut q2: f64 = 1.0;
        let mut q3: f64 = 1.0;

        for i in 0..delta.len() {
            if sigma_k * delta[i] == ps[i] {
                // q1 = q1 * probs[i];
                q1 = q1 * (1.0 - probs[i]);
            } else {
                // q1 = q1 * (1.0 - probs[i]);
                q1 = q1 * probs[i];
            }
        }

        for i in 0..delta.len() {
            if delta[i] == ps[i] {
                // q2 = q2 * probs[i];
                // q3 = q3 * (1.0 - probs[i]);
                q2 = q2 * (1.0 - probs[i]);
                q3 = q3 * probs[i];
            } else {
                // q2 = q2 * (1.0 - probs[i]);
                // q3 = q3 * probs[i];
                q2 = q2 * probs[i];
                q3 = q3 * (1.0 - probs[i]);
            }
        }
        return 10e-9_f64.max(q1 / (q2 + q3));
    }

    pub fn cal_log_sigma_delta(sigma_k: i32, delta: &Vec<i32>, ps: &Vec<i32>, probs: &Vec<f64>) -> f64 {
        // calculate P(sigma_k | delta)
        // sigma_k: the assignment of read k, 1 or -1.
        // delta: the haplotypes of the SNPs covered by read k, each haplotype is 1 or -1.
        // ps: the allele of each base, 1,-1
        // probs: the probability of observing base at each SNP for read k, equals to 10^(-Q/10).
        let mut q1: f64 = 0.0;
        let mut q2: f64 = 0.0;
        let mut q3: f64 = 0.0;

        for i in 0..delta.len() {
            if sigma_k * delta[i] == ps[i] {
                q1 += (1.0 - probs[i]).log10();
            } else {
                q1 += probs[i].log10();
            }
        }

        for i in 0..delta.len() {
            if delta[i] == ps[i] {
                q2 += (1.0 - probs[i]).log10();
                q3 += probs[i].log10();
            } else {
                q2 += probs[i].log10();
                q3 += (1.0 - probs[i]).log10();
            }
        }
        return 10e-9_f64.max(1.0 - q1 / (q2 + q3));
    }

    pub fn cal_delta_sigma(delta_i: i32, sigma: &Vec<i32>, ps: &Vec<i32>, probs: &Vec<f64>) -> f64 {
        // calculate P(delta_i | sigma)
        // delta_i: the haplotype of SNP i, 1 or -1.
        // sigma: the assignments of the reads cover SNP i, each haplotype is 1 or -1.
        // ps: the allele of each base, 1,-1
        // probs: the probability of observing base at SNP i for each read, equals to 10^(-Q/10).

        let mut q1: f64 = 1.0;
        let mut q2: f64 = 1.0;
        let mut q3: f64 = 1.0;

        for k in 0..sigma.len() {
            if delta_i * sigma[k] == ps[k] {
                // q1 = q1 * probs[k];
                q1 = q1 * (1.0 - probs[k]);
            } else {
                // q1 = q1 * (1.0 - probs[k]);
                q1 = q1 * probs[k];
            }
        }

        for k in 0..sigma.len() {
            if sigma[k] == ps[k] {
                // q2 = q2 * probs[k];
                // q3 = q3 * (1.0 - probs[k]);
                q2 = q2 * (1.0 - probs[k]);
                q3 = q3 * probs[k];
            } else {
                // q2 = q2 * (1.0 - probs[k]);
                // q3 = q3 * probs[k];
                q2 = q2 * probs[k];
                q3 = q3 * (1.0 - probs[k]);
            }
        }
        return 10e-9_f64.max(q1 / (q2 + q3));
    }

    pub fn cal_log_delta_sigma(delta_i: i32, sigma: &Vec<i32>, ps: &Vec<i32>, probs: &Vec<f64>) -> f64 {
        // calculate P(delta_i | sigma)
        // delta_i: the haplotype of SNP i, 1 or -1.
        // sigma: the assignments of the reads cover SNP i, each haplotype is 1 or -1.
        // ps: the allele of each base, 1,-1
        // probs: the probability of observing base at SNP i for each read, equals to 10^(-Q/10).

        let mut q1: f64 = 0.0;
        let mut q2: f64 = 0.0;
        let mut q3: f64 = 0.0;

        for k in 0..sigma.len() {
            if delta_i * sigma[k] == ps[k] {
                q1 += (1.0 - probs[k]).log10();
            } else {
                q1 += probs[k].log10();
            }
        }

        for k in 0..sigma.len() {
            if sigma[k] == ps[k] {
                q2 += (1.0 - probs[k]).log10();
                q3 += probs[k].log10();
            } else {
                q2 += probs[k].log10();
                q3 += (1.0 - probs[k]).log10();
            }
        }
        return 10e-9_f64.max(1.0 - q1 / (q2 + q3));
    }

    pub fn cal_delta_sigma_sum(delta_i: i32, sigma: &Vec<i32>, ps: &Vec<i32>, probs: &Vec<f64>) -> f64 {
        // Similar in cal_delta_sigma. Due to the numerical calculation, we replace the multiply of probabilities with the sum of the logarithms of the probabilities.

        let mut q1: f64 = 0.0;
        let mut q2: f64 = 0.0;
        let mut q3: f64 = 0.0;

        for k in 0..sigma.len() {
            if delta_i * sigma[k] == ps[k] {
                // q1 = q1 + probs[k].log10();
                q1 = q1 + (1.0 - probs[k]).log10();
            } else {
                // q1 = q1 + (1.0 - probs[k]).log10();
                q1 = q1 + probs[k].log10();
            }
        }

        for k in 0..sigma.len() {
            if sigma[k] == ps[k] {
                // q2 = q2 + probs[k].log10();
                // q3 = q3 + (1.0 - probs[k]).log10();
                q2 = q2 + (1.0 - probs[k]).log10();
                q3 = q3 + probs[k].log10();
            } else {
                // q2 = q2 + (1.0 - probs[k]).log10();
                // q3 = q3 + probs[k].log10();
                q2 = q2 + probs[k].log10();
                q3 = q3 + (1.0 - probs[k]).log10();
            }
        }
        return q1 / (q2 + q3);
    }

    pub fn cal_inconsistent_percentage(delta_i: i32, sigma: &Vec<i32>, ps: &Vec<i32>, probs: &Vec<f64>) -> f64 {
        let mut consisitent = 0;
        let mut inconsistent = 0;
        for k in 0..sigma.len() {
            if delta_i * sigma[k] == ps[k] {
                consisitent += 1;
            } else {
                inconsistent += 1;
            }
        }
        return 10e-6_f64.max((inconsistent as f64) / ((consisitent + inconsistent) as f64));
    }

    // pub fn cal_overall_probability(snpfrag: &SNPFrag, sigma: &Vec<i32>, delta: &Vec<i32>) -> f64 {
    //     let mut logp = 0.0;
    //     for k in 0..sigma.len() {
    //         for fe in snpfrag.fragments[k].list.iter() {
    //             if fe.p != 0 {
    //                 let i = fe.snp_idx;
    //                 if sigma[k] * delta[i] == fe.p {
    //                     // logp += fe.prob.log10();
    //                     logp += (1.0 - fe.prob).log10();
    //                 } else {
    //                     // logp += (1.0 - fe.prob).log10();
    //                     logp += fe.prob.log10();
    //                 }
    //             }
    //         }
    //     }
    //     return logp;
    // }

    pub fn cal_overall_probability(snpfrag: &SNPFrag, processed_snps: &HashSet<usize>, covered_fragments: &HashSet<usize>) -> f64 {
        let mut logp = 0.0;
        for k in covered_fragments.iter() {
            for fe in snpfrag.fragments[*k].list.iter() {
                if snpfrag.candidate_snps[fe.snp_idx].filter == false && processed_snps.contains(&fe.snp_idx) {
                    if fe.p != 0 {
                        if snpfrag.fragments[*k].haplotag * snpfrag.candidate_snps[fe.snp_idx].haplotype == fe.p {
                            logp += (1.0 - fe.prob).log10();
                        } else {
                            logp += fe.prob.log10();
                        }
                    }
                }
            }
        }
        return logp;
    }

    pub fn check_new_haplotag(snpfrag: &SNPFrag, processed_snps: &HashSet<usize>, updated_haplotag: &HashMap<usize, i32>) -> i32 {
        let mut logp = 0.0;
        let mut pre_logp = 0.0;
        for (k, h) in updated_haplotag.iter()
        {
            for fe in snpfrag.fragments[*k].list.iter() {
                if snpfrag.candidate_snps[fe.snp_idx].filter == false && processed_snps.contains(&fe.snp_idx) {
                    if fe.p != 0 {
                        if h * snpfrag.candidate_snps[fe.snp_idx].haplotype == fe.p {
                            logp += (1.0 - fe.prob).log10();
                        } else {
                            logp += fe.prob.log10();
                        }

                        if snpfrag.fragments[*k].haplotag * snpfrag.candidate_snps[fe.snp_idx].haplotype == fe.p {
                            pre_logp += (1.0 - fe.prob).log10();
                        } else {
                            pre_logp += fe.prob.log10();
                        }
                    }
                }
            }
        }
        // println!("check_new_haplotag: logp = {}, pre_logp = {}", logp, pre_logp);
        let mut rv = 0;
        if logp > pre_logp { rv = 1; } else if logp == pre_logp { rv = 0; } else {
            // println!("check_new_haplotag: logp = {}, pre_logp = {}", logp, pre_logp);
            rv = -1;
        }
        return rv;
    }


    pub fn check_new_haplotype(snpfrag: &SNPFrag, updated_haplotype: &HashMap<usize, i32>) -> i32 {
        let mut logp = 0.0;
        let mut pre_logp = 0.0;
        for (i, h) in updated_haplotype.iter()
        {
            for k in snpfrag.candidate_snps[*i].snp_cover_fragments.iter() {
                for fe in snpfrag.fragments[*k].list.iter() {
                    if fe.snp_idx != *i { continue; }
                    if fe.p != 0 {
                        if snpfrag.fragments[*k].haplotag * h == fe.p {
                            logp += (1.0 - fe.prob).log10();
                        } else {
                            logp += fe.prob.log10();
                        }

                        if snpfrag.fragments[*k].haplotag * snpfrag.candidate_snps[fe.snp_idx].haplotype == fe.p {
                            pre_logp += (1.0 - fe.prob).log10();
                        } else {
                            pre_logp += fe.prob.log10();
                        }
                    }
                }
            }
        }
        // println!("check_new_haplotype: logp = {}, pre_logp = {}", logp, pre_logp);
        let mut rv = 0;
        if logp > pre_logp { rv = 1; } else if logp == pre_logp { rv = 0; } else { rv = -1; }
        return rv;
    }

    pub fn cross_optimize(&mut self) -> f64 {
        // Iteration:
        //     1. evaluate the assignment of each read based on the current SNP haplotype.
        //     2. evaluate the SNP haplotype based on the read assignment.
        // If P(sigma, delta) increase, repeat Iteration;
        // Else break;

        // self.optimization_using_maxcut();   // get the initial SNP haplotype, assume most of SNP haplotype is correct.

        let mut phasing_increase: bool = true;
        let mut haplotag_increase: bool = true;
        let mut num_iters = 0;
        let mut fragments_cnt: HashMap<usize, usize> = HashMap::new();
        let mut covered_fragments: HashSet<usize> = HashSet::new();
        let mut processed_snps: HashSet<usize> = HashSet::new();
        for i in self.hete_snps.iter() {
            processed_snps.insert(*i);
            for k in self.candidate_snps[*i].snp_cover_fragments.iter() {
                // covered_fragments.insert(*k);
                if !fragments_cnt.contains_key(k) {
                    fragments_cnt.insert(*k, 1);
                } else {
                    let cnt = fragments_cnt.get_mut(k).unwrap();
                    *cnt += 1;
                }
            }
        }

        for (k, v) in fragments_cnt.iter() {
            if *v >= 2 {
                covered_fragments.insert(*k);
            }
        }

        while phasing_increase | haplotag_increase {
            println!("{} {}:{}-{} cross optimization iterations: {}", Local::now().format("%Y-%m-%d %H:%M:%S"), self.region.chr, self.region.start, self.region.end, num_iters);
            // optimize sigma
            let mut tmp_haplotag: HashMap<usize, i32> = HashMap::new();
            for k in covered_fragments.iter() {
                let sigma_k = self.fragments[*k].haplotag;
                let mut delta: Vec<i32> = Vec::new();
                let mut ps: Vec<i32> = Vec::new();
                let mut probs: Vec<f64> = Vec::new();
                for fe in self.fragments[*k].list.iter() {
                    if self.candidate_snps[fe.snp_idx].filter == false && processed_snps.contains(&fe.snp_idx) {
                        if fe.p != 0 {
                            ps.push(fe.p);
                            probs.push(fe.prob);
                            delta.push(self.candidate_snps[fe.snp_idx].haplotype);
                        }
                    }
                }

                if SNPFrag::cal_log_sigma_delta(sigma_k, &delta, &ps, &probs) < SNPFrag::cal_log_sigma_delta(sigma_k * (-1), &delta, &ps, &probs) {
                    // tmp_haplotag.push(sigma_k * (-1));  // flip sigma_k
                    tmp_haplotag.insert(*k, sigma_k * (-1));
                } else {
                    // tmp_haplotag.push(sigma_k);
                    tmp_haplotag.insert(*k, sigma_k);
                }
            }

            // assert!(SNPFrag::cal_overall_probability(&self, &processed_snps, &self.haplotype) >= SNPFrag::cal_overall_probability(&self, &self.haplotag, &self.haplotype));
            let check_val = SNPFrag::check_new_haplotag(&self, &processed_snps, &tmp_haplotag);
            // assert!(check_val >= 0, "ckeck val bug: {:?}", self.candidate_snps);

            if check_val > 0 {
                // new haplotag increases the probability P(sigma, delta)
                for (k, h) in tmp_haplotag.iter() {
                    self.fragments[*k].haplotag = *h;
                }
            } else {
                haplotag_increase = false;
            }

            // optimize delta
            let mut tmp_haplotype: HashMap<usize, i32> = HashMap::new();
            for i in self.hete_snps.iter() {
                let delta_i = self.candidate_snps[*i].haplotype;
                let mut sigma: Vec<i32> = Vec::new();
                let mut ps: Vec<i32> = Vec::new();
                let mut probs: Vec<f64> = Vec::new();
                for k in self.candidate_snps[*i].snp_cover_fragments.iter() {
                    // k is fragment index
                    for fe in self.fragments[*k].list.iter() {
                        if self.candidate_snps[fe.snp_idx].filter == false && fe.snp_idx == *i {
                            if fe.p != 0 {
                                ps.push(fe.p);
                                probs.push(fe.prob);
                                sigma.push(self.fragments[*k].haplotag);
                            }
                        }
                    }
                }

                if SNPFrag::cal_log_delta_sigma(delta_i, &sigma, &ps, &probs) < SNPFrag::cal_log_delta_sigma(delta_i * (-1), &sigma, &ps, &probs) {
                    tmp_haplotype.insert(*i, delta_i * (-1));
                    // tmp_haplotype.push(delta_i * (-1));
                } else {
                    tmp_haplotype.insert(*i, delta_i);
                    // tmp_haplotype.push(delta_i);
                }
            }

            let check_val = SNPFrag::check_new_haplotype(&self, &tmp_haplotype);
            // assert!(check_val >= 0, "ckeck val bug: {:?}", self.candidate_snps);

            if check_val > 0 {
                // new haplotype increases the probability P(sigma, delta)
                for (i, h) in tmp_haplotype.iter() {
                    self.candidate_snps[*i].haplotype = *h;
                }
            } else {
                phasing_increase = false;
            }
            num_iters += 1;
            if num_iters > 10 {
                println!("Cross optimization may not converge: {}:{}-{}", self.region.chr, self.region.start, self.region.end);
                break;
            }
        }

        let prob = SNPFrag::cal_overall_probability(&self, &processed_snps, &covered_fragments);
        // println!("cross optimization iterations: {}", num_iters);
        // println!("haplotype: {:?}", self.haplotype);
        // println!("score: {:?}", prob);
        // println!("assignment: {:?}", self.haplotag);
        return prob;
    }

    pub fn phase(&mut self, max_enum_snps: usize, random_flip_fraction: f32) {
        let mut largest_prob = f64::NEG_INFINITY;
        let mut best_haplotype: HashMap<usize, i32> = HashMap::new();
        let mut best_haplotag: HashMap<usize, i32> = HashMap::new();
        let mut covered_fragments: HashSet<usize> = HashSet::new();
        let mut fragments_cnt: HashMap<usize, usize> = HashMap::new();
        for i in self.hete_snps.iter() {
            // if self.candidate_snps[*i].filter == true { continue; }
            for k in self.candidate_snps[*i].snp_cover_fragments.iter() {
                if !fragments_cnt.contains_key(k) {
                    fragments_cnt.insert(*k, 1);
                } else {
                    let cnt = fragments_cnt.get_mut(k).unwrap();
                    *cnt += 1;
                }
            }
        }


        for (k, v) in fragments_cnt.iter() {
            // each fragment is covered by at least two SNPs
            if *v >= 2 {
                covered_fragments.insert(*k);
            }
        }

        // println!("Phasing informations:");
        // for i in self.hete_snps.iter() {
        //     println!("SNP {:?}", self.candidate_snps[*i]);
        //     for k in self.candidate_snps[*i].snp_cover_fragments.iter() {
        //         if covered_fragments.contains(k) {
        //             println!("Fragment {:?}", self.fragments[*k]);
        //         }
        //     }
        // }
        // println!();

        if self.hete_snps.len() <= max_enum_snps {
            // enumerate the haplotype, then optimize the assignment
            let mut haplotype_enum: Vec<Vec<i32>> = Vec::new();
            let init_hap: Vec<i32> = vec![1; self.hete_snps.len()];
            haplotype_enum.push(init_hap.clone());
            for ti in 0..self.hete_snps.len() {
                for tj in 0..haplotype_enum.len() {
                    let mut tmp_hap = haplotype_enum[tj].clone();
                    tmp_hap[ti] = tmp_hap[ti] * (-1);
                    assert!(tmp_hap.len() == self.hete_snps.len());
                    haplotype_enum.push(tmp_hap);
                }
            }
            for hap in haplotype_enum.iter() {
                for i in 0..self.hete_snps.len() {
                    self.candidate_snps[self.hete_snps[i]].haplotype = hap[i];
                }
                let prob = self.cross_optimize();
                if prob > largest_prob {
                    largest_prob = prob;
                    best_haplotype.clear();
                    best_haplotag.clear();
                    for i in self.hete_snps.iter() {
                        best_haplotype.insert(*i, self.candidate_snps[*i].haplotype);
                    }
                    for k in covered_fragments.iter() {
                        best_haplotag.insert(*k, self.fragments[*k].haplotag);
                    }
                }
            }
            for i in self.hete_snps.iter() {
                self.candidate_snps[*i].haplotype = best_haplotype[i];
            }
            for k in covered_fragments.iter() {
                self.fragments[*k].haplotag = best_haplotag[k];
            }
        } else {
            // optimize haplotype and read assignment alternatively
            let mut max_iter = 5;
            while max_iter >= 0 {
                unsafe { self.init_haplotypes(); }
                unsafe { self.init_assignment(); }
                let prob = self.cross_optimize();
                if prob > largest_prob {
                    largest_prob = prob;
                    best_haplotype.clear();
                    best_haplotag.clear();
                    for i in self.hete_snps.iter() {
                        best_haplotype.insert(*i, self.candidate_snps[*i].haplotype);
                    }
                    for k in covered_fragments.iter() {
                        best_haplotag.insert(*k, self.fragments[*k].haplotag);
                    }
                }
                for i in self.hete_snps.iter() {
                    self.candidate_snps[*i].haplotype = best_haplotype[i];
                }
                for k in covered_fragments.iter() {
                    self.fragments[*k].haplotag = best_haplotag[k];
                }

                // block flip: flip all the snps after a random position to jump local optimization
                // let mut unflipped_haplotype: Vec<i32> = Vec::new();
                // for i in self.hete_snps.iter() {
                //     unflipped_haplotype.push(self.candidate_snps[*i].haplotype);
                // }
                // for ti in 0..unflipped_haplotype.len() {
                //     let mut tmp_hap: Vec<i32> = Vec::new();
                //     for tj in 0..unflipped_haplotype.len() {
                //         if tj < ti {
                //             tmp_hap.push(unflipped_haplotype[tj]);
                //         } else {
                //             tmp_hap.push(unflipped_haplotype[tj] * (-1));
                //         }
                //     }
                //     assert!(tmp_hap.len() == self.hete_snps.len());
                //     for i in 0..self.hete_snps.len() {
                //         self.candidate_snps[self.hete_snps[i]].haplotype = tmp_hap[i];
                //     }
                //     let prob = self.cross_optimize();
                //     if prob > largest_prob {
                //         largest_prob = prob;
                //         best_haplotype.clear();
                //         best_haplotag.clear();
                //         for i in self.hete_snps.iter() {
                //             best_haplotype.insert(*i, self.candidate_snps[*i].haplotype);
                //         }
                //         for k in covered_fragments.iter() {
                //             best_haplotag.insert(*k, self.fragments[*k].haplotag);
                //         }
                //     }
                // }
                // for i in self.hete_snps.iter() {
                //     self.candidate_snps[*i].haplotype = best_haplotype[i];
                // }
                // for k in covered_fragments.iter() {
                //     self.fragments[*k].haplotag = best_haplotag[k];
                // }

                // flip a fraction of snps and reads
                let num_flip_haplotype = (self.hete_snps.len() as f32 * random_flip_fraction) as usize;
                let num_flip_haplotag = (covered_fragments.len() as f32 * random_flip_fraction) as usize;
                let mut rng = rand::thread_rng();
                let selected_flip_haplotype: Vec<_> = self.hete_snps.iter().collect::<Vec<_>>().choose_multiple(&mut rng, num_flip_haplotype).cloned().collect();
                for ti in selected_flip_haplotype {
                    self.candidate_snps[*ti].haplotype = self.candidate_snps[*ti].haplotype * (-1);
                }
                let selected_flip_haplotag: Vec<_> = covered_fragments.iter().collect::<Vec<_>>().choose_multiple(&mut rng, num_flip_haplotag).cloned().collect();
                for tk in selected_flip_haplotag {
                    self.fragments[*tk].haplotag = self.fragments[*tk].haplotag * (-1);
                }
                let prob = self.cross_optimize();
                if prob > largest_prob {
                    largest_prob = prob;
                    best_haplotype.clear();
                    best_haplotag.clear();
                    for i in self.hete_snps.iter() {
                        best_haplotype.insert(*i, self.candidate_snps[*i].haplotype);
                    }
                    for k in covered_fragments.iter() {
                        best_haplotag.insert(*k, self.fragments[*k].haplotag);
                    }
                }
                max_iter -= 1;
                // println!("largest_prob: {}", largest_prob);
            }
            for i in self.hete_snps.iter() {
                self.candidate_snps[*i].haplotype = best_haplotype[i];
            }
            for k in covered_fragments.iter() {
                self.fragments[*k].haplotag = best_haplotag[k];
            }
        }
    }

    pub fn add_phase_score(&mut self, min_allele_cnt: u32, imbalance_allele_expression_cutoff: f32) {
        // calculate phase score for each snp
        for ti in 0..self.candidate_snps.len() {
            let snp = &self.candidate_snps[ti];
            if snp.filter == true || snp.variant_type != 1 {
                continue;
            }
            let delta_i = snp.haplotype;
            let mut sigma: Vec<i32> = Vec::new();
            let mut ps: Vec<i32> = Vec::new();
            let mut probs: Vec<f64> = Vec::new();
            let mut num_hap1 = 0;
            let mut num_hap2 = 0;
            for k in snp.snp_cover_fragments.iter() {
                if self.fragments[*k].assignment == 0 {
                    continue;
                }
                for fe in self.fragments[*k].list.iter() {
                    if fe.snp_idx == ti {
                        if fe.base != '-' {
                            // ignore intron
                            if self.fragments[*k].assignment == 1 {
                                num_hap1 += 1;
                            } else if self.fragments[*k].assignment == 2 {
                                num_hap2 += 1;
                            }
                        }
                        if fe.p != 0 {
                            ps.push(fe.p);
                            probs.push(fe.prob);
                            sigma.push(self.fragments[*k].haplotag);
                        }
                    }
                }
            }

            let mut phase_score = 0.0;

            if num_hap1 < min_allele_cnt || num_hap2 < min_allele_cnt {
                // filter SNPs with low allele count, no confident phase
                println!("Low allele count: {}:{}, {}, {}", self.region.chr, snp.pos, num_hap1, num_hap2);
                phase_score = 0.0;
            } else {
                if sigma.len() > 0 {
                    // println!("delta_i: {:?}", delta_i);
                    // println!("sigma: {:?}", sigma);
                    // println!("ps: {:?}", ps);
                    // phase_score = -10.0_f64 * SNPFrag::cal_inconsistent_percentage(delta_i, &sigma, &ps, &probs).log10();
                    phase_score = -10.0_f64 * (1.0 - SNPFrag::cal_log_delta_sigma(delta_i, &sigma, &ps, &probs)).log10();   // calaulate assignment score
                } else {
                    phase_score = 0.0;  // all reads belong to unknown group, maybe caused by single SNP
                }
            }

            // used for calculating allele imbalance expression
            let mut hap1_allele_cnt: [u32; 2] = [0, 0];
            let mut hap2_allele_cnt: [u32; 2] = [0, 0];
            for k in 0..sigma.len() {
                if delta_i * sigma[k] == ps[k] {
                    if sigma[k] == 1 {
                        hap1_allele_cnt[0] += 1;
                    } else {
                        hap2_allele_cnt[0] += 1;
                    }
                } else {
                    if sigma[k] == 1 {
                        hap1_allele_cnt[1] += 1;
                    } else {
                        hap2_allele_cnt[1] += 1;
                    }
                }
            }
            let mut allele_imbalance: bool = false;
            if max(hap1_allele_cnt[0], hap1_allele_cnt[1]) as f32 > max(hap2_allele_cnt[0], hap2_allele_cnt[1]) as f32 * imbalance_allele_expression_cutoff {
                allele_imbalance = true;
            } else if max(hap2_allele_cnt[0], hap2_allele_cnt[1]) as f32 * imbalance_allele_expression_cutoff < max(hap1_allele_cnt[0], hap1_allele_cnt[1]) as f32 {
                allele_imbalance = true;
            }
            // if phase_score > 8.0 {
            //     println!("IMB {}: {}, {:?}, {:?}", snp.pos, phase_score, hap1_allele_cnt, hap2_allele_cnt);
            // }
            self.candidate_snps[ti].allele_imbalance = allele_imbalance;
            self.candidate_snps[ti].phase_score = phase_score;
        }
    }

    pub fn assign_reads(&mut self, read_assignment_cutoff: f64) -> HashMap<String, i32> {
        let mut read_assignments: HashMap<String, i32> = HashMap::new();
        let mut processed_snps: HashSet<usize> = HashSet::new();
        for i in self.hete_snps.iter() {
            processed_snps.insert(*i);
        }
        for k in 0..self.fragments.len() {
            let sigma_k = self.fragments[k].haplotag;
            let mut delta: Vec<i32> = Vec::new();
            let mut ps: Vec<i32> = Vec::new();
            let mut probs: Vec<f64> = Vec::new();
            for fe in self.fragments[k].list.iter() {
                if self.candidate_snps[fe.snp_idx].filter == false && processed_snps.contains(&fe.snp_idx) {
                    if fe.p != 0 {
                        ps.push(fe.p);
                        probs.push(fe.prob);
                        delta.push(self.candidate_snps[fe.snp_idx].haplotype);
                    }
                }
            }

            if sigma_k == 0 {
                // unasigned haplotag, cluster the read into unknown group
                self.fragments[k].assignment = 0;
                self.fragments[k].assignment_score = 0.0;
                read_assignments.insert(self.fragments[k].read_id.clone(), 0);
            } else {
                let q = SNPFrag::cal_log_sigma_delta(sigma_k, &delta, &ps, &probs);
                let qn = SNPFrag::cal_log_sigma_delta(sigma_k * (-1), &delta, &ps, &probs);
                if q - qn > read_assignment_cutoff {
                    if sigma_k == 1 {
                        self.fragments[k].assignment = 1;
                        self.fragments[k].assignment_score = q;
                        read_assignments.insert(self.fragments[k].read_id.clone(), 1);
                    } else {
                        self.fragments[k].assignment = 2;
                        self.fragments[k].assignment_score = q;
                        read_assignments.insert(self.fragments[k].read_id.clone(), 2);
                    }
                } else {
                    // unknown which haplotype the read belongs to, cluster the read into unknown group
                    self.fragments[k].assignment = 0;
                    self.fragments[k].assignment_score = 0.0;
                    read_assignments.insert(self.fragments[k].read_id.clone(), 0);
                }
            }
        }
        return read_assignments;
    }

    // pub fn output_vcf(&self) -> Vec<VCFRecord> {
    //     let mut records: Vec<VCFRecord> = Vec::new();
    //     assert_eq!(self.haplotype.len(), self.snps.len());
    //     for i in 0..self.haplotype.len() {
    //         let snp = &self.snps[i];
    //         let hp = self.haplotype[i];
    //         let mut phase_score = 0;
    //         // use \delta_{i}\delta_{j}W_{ij} as phased score
    //         for edge in self.edges.iter() {
    //             if edge.0[0] == i || edge.0[1] == i {
    //                 phase_score += (edge.1.w as i32) * self.haplotype[edge.0[0]] * self.haplotype[edge.0[1]];
    //             }
    //         }
    //         let mut rd: VCFRecord = VCFRecord::default();
    //         rd.chromosome = snp.chromosome.clone();
    //         rd.position = snp.pos as u64 + 1;   // position in vcf format is 1-based
    //         rd.reference = vec![snp.reference as u8];
    //         rd.id = vec!['.' as u8];
    //         if snp.alleles[0] == snp.reference {
    //             rd.alternative = vec![vec![snp.alleles[1] as u8]];
    //             rd.qual = cmp::max(1, (-10.0 * f64::log10(0.01_f64.max((0.5 - snp.allele_freqs[1] as f64).abs() / 0.5))) as i32);
    //             if hp == -1 {
    //                 rd.genotype = format!("{}:{}:{}:{:.2}:{:.2}", "0|1", rd.qual, snp.depth, snp.allele_freqs[1], phase_score);
    //             } else {
    //                 rd.genotype = format!("{}:{}:{}:{:.2}:{:.2}", "1|0", rd.qual, snp.depth, snp.allele_freqs[1], phase_score);
    //             }
    //         } else if snp.alleles[1] == snp.reference {
    //             rd.alternative = vec![vec![snp.alleles[0] as u8]];
    //             // rd.qual = -10.0 * f32::log10((0.5 - snp.allele_freqs[0]).abs() / 0.5);
    //             rd.qual = cmp::max(1, (-10.0 * f64::log10(0.01_f64.max((0.5 - snp.allele_freqs[0] as f64).abs() / 0.5))) as i32);
    //             if hp == -1 {
    //                 rd.genotype = format!("{}:{}:{}:{:.2}:{:.2}", "0|1", rd.qual, snp.depth, snp.allele_freqs[0], phase_score);
    //             } else {
    //                 rd.genotype = format!("{}:{}:{}:{:.2}:{:.2}", "1|0", rd.qual, snp.depth, snp.allele_freqs[0], phase_score);
    //             }
    //         } else {
    //             rd.alternative = vec![vec![snp.alleles[0] as u8], vec![snp.alleles[1] as u8]];
    //             let q1 = cmp::max(1, (-10.0 * f64::log10(0.01_f64.max((0.5 - snp.allele_freqs[0] as f64).abs() / 0.5))) as i32);
    //             let q2 = cmp::max(1, (-10.0 * f64::log10(0.01_f64.max((0.5 - snp.allele_freqs[1] as f64).abs() / 0.5))) as i32);
    //             // if q1 > q2 { rd.qual = q2; } else { rd.qual = q1; }
    //             rd.qual = cmp::min(q1, q2);
    //             rd.genotype = format!("{}:{}:{}:{},{:.2}:{:.2}", "1|2", rd.qual, snp.depth, snp.allele_freqs[0], snp.allele_freqs[1], phase_score);
    //         }
    //
    //         rd.filter = "PASS".to_string().into_bytes();
    //         rd.info = ".".to_string().into_bytes();
    //         rd.format = "GT:GQ:DP:AF:PQ".to_string().into_bytes();
    //         records.push(rd);
    //     }
    //     return records;
    // }

    pub fn phased_output_vcf(&mut self, min_phase_score: f32, min_homozygous_freq: f32, output_phasing: bool, min_qual_for_candidate: u32, min_qual_for_singlesnp_rnaedit: u32) -> Vec<VCFRecord> {
        let mut records: Vec<VCFRecord> = Vec::new();

        // output heterozygous SNPs
        // assert_eq!(self.haplotype.len(), self.snps.len());
        for i in 0..self.candidate_snps.len() {
            let snp = &self.candidate_snps[i];
            let hp = snp.haplotype;
            if snp.filter == true && snp.rna_editing == false {
                // dense SNP
                let mut rd: VCFRecord = VCFRecord::default();
                rd.chromosome = snp.chromosome.clone();
                rd.position = snp.pos as u64 + 1;   // position in vcf format is 1-based
                rd.reference = vec![snp.reference as u8];
                rd.id = vec!['.' as u8];
                if snp.variant_type == 1 {
                    if snp.alleles[0] != snp.reference && snp.alleles[1] == snp.reference {
                        rd.alternative = vec![vec![snp.alleles[0] as u8]];
                        rd.genotype = format!("{}:{}:{}:{:.2}", "0/1", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[0]);
                    } else if snp.alleles[1] != snp.reference && snp.alleles[0] == snp.reference {
                        rd.alternative = vec![vec![snp.alleles[1] as u8]];
                        rd.genotype = format!("{}:{}:{}:{:.2}", "0/1", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[1]);
                    } else {
                        rd.alternative = vec![vec![snp.alleles[0] as u8], vec![snp.alleles[1] as u8]];
                        rd.genotype = format!("{}:{}:{}:{:.2},{:.2}", "1/2", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[0], snp.allele_freqs[1]);
                    }
                } else if snp.variant_type == 2 {
                    rd.alternative = vec![vec![snp.alleles[0] as u8]];
                    rd.genotype = format!("{}:{}:{}:{:.2}", "1/1", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[0]);
                } else if snp.variant_type == 3 {
                    rd.alternative = vec![vec![snp.alleles[0] as u8], vec![snp.alleles[1] as u8]];
                    rd.genotype = format!("{}:{}:{}:{:.2},{:.2}", "1/2", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[0], snp.allele_freqs[1]);
                }
                rd.qual = snp.variant_quality as i32;
                rd.filter = "dn".to_string().into_bytes();
                rd.info = format!("RDS={}", "dense_snp").to_string().into_bytes();
                rd.format = "GT:GQ:DP:AF".to_string().into_bytes();
                records.push(rd);
                continue;
            }

            if snp.variant_type == 3 {
                let mut rd: VCFRecord = VCFRecord::default();
                rd.chromosome = snp.chromosome.clone();
                rd.position = snp.pos as u64 + 1;   // position in vcf format is 1-based
                rd.reference = vec![snp.reference as u8];
                rd.id = vec!['.' as u8];
                rd.alternative = vec![vec![snp.alleles[0] as u8], vec![snp.alleles[1] as u8]];
                rd.qual = snp.variant_quality as i32;
                rd.genotype = format!("{}:{}:{}:{:.2},{:.2}", "1/2", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[0], snp.allele_freqs[1]);
                if snp.rna_editing == true {
                    rd.info = format!("RDS={}", "rna_editing").to_string().into_bytes();
                    if snp.variant_quality < min_qual_for_singlesnp_rnaedit as f64 {
                        rd.filter = "RnaEdit".to_string().into_bytes();
                    } else {
                        rd.filter = "PASS".to_string().into_bytes();
                    }
                    // rd.filter = "RnaEdit".to_string().into_bytes();
                } else {
                    rd.info = "RDS=.".to_string().into_bytes();
                    if snp.variant_quality < min_qual_for_candidate as f64 {
                        rd.filter = "LowQual".to_string().into_bytes();
                    } else {
                        rd.filter = "PASS".to_string().into_bytes();
                    }
                }
                rd.format = "GT:GQ:DP:AF".to_string().into_bytes();
                records.push(rd);
            } else if snp.variant_type == 2 {
                let mut rd: VCFRecord = VCFRecord::default();
                rd.chromosome = snp.chromosome.clone();
                rd.position = snp.pos as u64 + 1;   // position in vcf format is 1-based
                rd.reference = vec![snp.reference as u8];
                rd.id = vec!['.' as u8];
                rd.alternative = vec![vec![snp.alleles[0] as u8]];
                rd.qual = snp.variant_quality as i32;
                rd.genotype = format!("{}:{}:{}:{:.2}", "1/1", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[0]);
                if snp.rna_editing == true {
                    rd.info = format!("RDS={}", "rna_editing").to_string().into_bytes();
                    if snp.variant_quality < min_qual_for_singlesnp_rnaedit as f64 {
                        rd.filter = "RnaEdit".to_string().into_bytes();
                    } else {
                        rd.filter = "PASS".to_string().into_bytes();
                    }
                    // rd.filter = "RnaEdit".to_string().into_bytes();
                } else {
                    rd.info = "RDS=.".to_string().into_bytes();
                    if snp.variant_quality < min_qual_for_candidate as f64 {
                        rd.filter = "LowQual".to_string().into_bytes();
                    } else {
                        rd.filter = "PASS".to_string().into_bytes();
                    }
                }
                rd.format = "GT:GQ:DP:AF".to_string().into_bytes();
                records.push(rd);
            } else if snp.variant_type == 1 {
                if snp.phase_score < min_phase_score as f64 && snp.allele_freqs[0] > min_homozygous_freq && snp.alleles[0] != snp.reference {
                    let mut rd: VCFRecord = VCFRecord::default();
                    rd.chromosome = snp.chromosome.clone();
                    rd.position = snp.pos as u64 + 1;   // position in vcf format is 1-based
                    rd.reference = vec![snp.reference as u8];
                    rd.id = vec!['.' as u8];
                    rd.alternative = vec![vec![snp.alleles[0] as u8]];
                    rd.qual = snp.variant_quality as i32;
                    rd.genotype = format!("{}:{}:{}:{:.2}", "1/1", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[0]);
                    if snp.rna_editing == true {
                        rd.info = format!("RDS={}", "rna_editing").to_string().into_bytes();
                        if snp.variant_quality < min_qual_for_singlesnp_rnaedit as f64 {
                            rd.filter = "RnaEdit".to_string().into_bytes();
                        } else {
                            rd.filter = "PASS".to_string().into_bytes();
                        }
                        // rd.filter = "RnaEdit".to_string().into_bytes();
                    } else {
                        rd.info = "RDS=.".to_string().into_bytes();
                        if snp.variant_quality < min_qual_for_candidate as f64 {
                            rd.filter = "LowQual".to_string().into_bytes();
                        } else {
                            rd.filter = "PASS".to_string().into_bytes();
                        }
                    }
                    rd.format = "GT:GQ:DP:AF".to_string().into_bytes();
                    records.push(rd);
                } else {
                    let mut rd: VCFRecord = VCFRecord::default();
                    rd.chromosome = snp.chromosome.clone();
                    rd.position = snp.pos as u64 + 1;   // position in vcf format is 1-based
                    rd.reference = vec![snp.reference as u8];
                    rd.id = vec!['.' as u8];

                    // output information containing phasing information
                    if output_phasing {
                        if snp.rna_editing == true || snp.phase_score == 0.0 {
                            // rna editing no phase information or phase_score == 0: single SNP, no phasing information
                            if snp.alleles[0] == snp.reference {
                                rd.alternative = vec![vec![snp.alleles[1] as u8]];
                                rd.genotype = format!("{}:{}:{}:{:.2}", "0/1", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[1]);
                            } else if snp.alleles[1] == snp.reference {
                                rd.alternative = vec![vec![snp.alleles[0] as u8]];
                                rd.genotype = format!("{}:{}:{}:{:.2}", "0/1", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[0]);
                            } else {
                                rd.alternative = vec![vec![snp.alleles[0] as u8], vec![snp.alleles[1] as u8]];
                                rd.genotype = format!("{}:{}:{}:{:.2},{:.2}", "1/2", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[0], snp.allele_freqs[1]);
                            }
                            if snp.rna_editing == true {
                                rd.info = format!("RDS={}", "rna_editing").to_string().into_bytes();
                                if snp.variant_quality < min_qual_for_singlesnp_rnaedit as f64 {
                                    rd.filter = "RnaEdit".to_string().into_bytes();
                                } else {
                                    rd.filter = "PASS".to_string().into_bytes();
                                }
                                // rd.filter = "RnaEdit".to_string().into_bytes();
                            } else if snp.phase_score == 0.0 {
                                rd.info = format!("RDS={}", "single_snp").to_string().into_bytes();
                                if snp.variant_quality < min_qual_for_singlesnp_rnaedit as f64 {
                                    rd.filter = "LowQual".to_string().into_bytes();
                                } else {
                                    rd.filter = "PASS".to_string().into_bytes();
                                }
                            } else {
                                rd.info = "RDS=.".to_string().into_bytes();
                                if snp.variant_quality < min_qual_for_candidate as f64 {
                                    rd.filter = "LowQual".to_string().into_bytes();
                                } else {
                                    rd.filter = "PASS".to_string().into_bytes();
                                }
                            }
                            rd.qual = snp.variant_quality as i32;
                            rd.format = "GT:GQ:DP:AF".to_string().into_bytes();
                        } else {
                            if snp.alleles[0] == snp.reference {
                                rd.alternative = vec![vec![snp.alleles[1] as u8]];
                                rd.qual = snp.variant_quality as i32;
                                if hp == -1 {
                                    rd.genotype = format!("{}:{}:{}:{:.2}:{:.2}", "0|1", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[1], snp.phase_score);
                                } else {
                                    rd.genotype = format!("{}:{}:{}:{:.2}:{:.2}", "1|0", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[1], snp.phase_score);
                                }
                            } else if snp.alleles[1] == snp.reference {
                                rd.alternative = vec![vec![snp.alleles[0] as u8]];
                                rd.qual = snp.variant_quality as i32;
                                if hp == -1 {
                                    rd.genotype = format!("{}:{}:{}:{:.2}:{:.2}", "0|1", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[0], snp.phase_score);
                                } else {
                                    rd.genotype = format!("{}:{}:{}:{:.2}:{:.2}", "1|0", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[0], snp.phase_score);
                                }
                            } else {
                                rd.alternative = vec![vec![snp.alleles[0] as u8], vec![snp.alleles[1] as u8]];
                                rd.qual = snp.variant_quality as i32;
                                rd.genotype = format!("{}:{}:{}:{:.2},{:.2}:{:.2}", "1|2", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[0], snp.allele_freqs[1], snp.phase_score);
                            }
                            if snp.phase_score < min_phase_score as f64 || snp.variant_quality < min_qual_for_candidate as f64 {
                                // not confident phase or low variant quality
                                rd.filter = "LowQual".to_string().into_bytes();
                            } else {
                                rd.filter = "PASS".to_string().into_bytes();
                            }
                            if snp.allele_imbalance == true {
                                rd.info = "RDS=.;AlleleImbalance".to_string().into_bytes();
                            } else {
                                rd.info = "RDS=.".to_string().into_bytes();
                            }
                            rd.format = "GT:GQ:DP:AF:PQ".to_string().into_bytes();
                        }
                    }

                    // output information does not contain phasing information, but still filtered by phase score
                    if !output_phasing {
                        if snp.rna_editing == true || snp.phase_score == 0.0 {
                            // rna editing no phase information or phase_score == 0: single SNP, no phasing information
                            if snp.alleles[0] == snp.reference {
                                rd.alternative = vec![vec![snp.alleles[1] as u8]];
                                rd.genotype = format!("{}:{}:{}:{:.2}", "0/1", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[1]);
                            } else if snp.alleles[1] == snp.reference {
                                rd.alternative = vec![vec![snp.alleles[0] as u8]];
                                rd.genotype = format!("{}:{}:{}:{:.2}", "0/1", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[0]);
                            } else {
                                rd.alternative = vec![vec![snp.alleles[0] as u8], vec![snp.alleles[1] as u8]];
                                rd.genotype = format!("{}:{}:{}:{:.2},{:.2}", "1/2", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[0], snp.allele_freqs[1]);
                            }
                            if snp.rna_editing == true {
                                rd.info = format!("RDS={}", "rna_editing").to_string().into_bytes();
                                if snp.variant_quality < min_qual_for_singlesnp_rnaedit as f64 {
                                    rd.filter = "RnaEdit".to_string().into_bytes();
                                } else {
                                    rd.filter = "PASS".to_string().into_bytes();
                                }
                                // rd.filter = "RnaEdit".to_string().into_bytes();
                            } else if snp.phase_score == 0.0 {
                                rd.info = format!("RDS={}", "single_snp").to_string().into_bytes();
                                if snp.variant_quality < min_qual_for_singlesnp_rnaedit as f64 {
                                    rd.filter = "LowQual".to_string().into_bytes();
                                } else {
                                    rd.filter = "PASS".to_string().into_bytes();
                                }
                            } else {
                                rd.info = "RDS=.".to_string().into_bytes();
                                if snp.variant_quality < min_qual_for_candidate as f64 {
                                    rd.filter = "LowQual".to_string().into_bytes();
                                } else {
                                    rd.filter = "PASS".to_string().into_bytes();
                                }
                            }
                            rd.qual = snp.variant_quality as i32;
                            rd.format = "GT:GQ:DP:AF".to_string().into_bytes();
                        } else {
                            if snp.alleles[0] == snp.reference {
                                rd.alternative = vec![vec![snp.alleles[1] as u8]];
                                rd.qual = snp.variant_quality as i32;
                                rd.genotype = format!("{}:{}:{}:{:.2}", "0/1", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[1]);
                            } else if snp.alleles[1] == snp.reference {
                                rd.alternative = vec![vec![snp.alleles[0] as u8]];
                                rd.qual = snp.variant_quality as i32;
                                rd.genotype = format!("{}:{}:{}:{:.2}", "0/1", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[0]);
                            } else {
                                rd.alternative = vec![vec![snp.alleles[0] as u8], vec![snp.alleles[1] as u8]];
                                rd.qual = snp.variant_quality as i32;
                                rd.genotype = format!("{}:{}:{}:{:.2},{:.2}", "1/2", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[0], snp.allele_freqs[1]);
                            }
                            if snp.phase_score < min_phase_score as f64 || snp.variant_quality < min_qual_for_candidate as f64 {
                                // not confident phase or low variant quality
                                rd.filter = "LowQual".to_string().into_bytes();
                            } else {
                                rd.filter = "PASS".to_string().into_bytes();
                            }
                            rd.info = "RDS=.".to_string().into_bytes();
                            rd.format = "GT:GQ:DP:AF:PQ".to_string().into_bytes();
                        }
                    }
                    records.push(rd);
                }
            } else {
                println!("Unknown variant type: {:?}:{:?}", snp.chromosome, snp.pos);
                continue;
            }
        }
        return records;
    }

    pub fn output_vcf(&mut self, min_qual: u32) -> Vec<VCFRecord> {
        let mut records: Vec<VCFRecord> = Vec::new();

        // output heterozygous SNPs
        // assert_eq!(self.haplotype.len(), self.snps.len());
        for i in 0..self.candidate_snps.len() {
            let snp = &self.candidate_snps[i];
            if snp.filter == true && snp.rna_editing == false {
                // dense SNP
                let mut rd: VCFRecord = VCFRecord::default();
                rd.chromosome = snp.chromosome.clone();
                rd.position = snp.pos as u64 + 1;   // position in vcf format is 1-based
                rd.reference = vec![snp.reference as u8];
                rd.id = vec!['.' as u8];
                if snp.variant_type == 1 {
                    if snp.alleles[0] != snp.reference && snp.alleles[1] == snp.reference {
                        rd.alternative = vec![vec![snp.alleles[0] as u8]];
                        rd.genotype = format!("{}:{}:{}:{:.2}", "0/1", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[0]);
                    } else if snp.alleles[1] != snp.reference && snp.alleles[0] == snp.reference {
                        rd.alternative = vec![vec![snp.alleles[1] as u8]];
                        rd.genotype = format!("{}:{}:{}:{:.2}", "0/1", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[1]);
                    } else {
                        rd.alternative = vec![vec![snp.alleles[0] as u8], vec![snp.alleles[1] as u8]];
                        rd.genotype = format!("{}:{}:{}:{:.2},{:.2}", "1/2", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[0], snp.allele_freqs[1]);
                    }
                } else if snp.variant_type == 2 {
                    rd.alternative = vec![vec![snp.alleles[0] as u8]];
                    rd.genotype = format!("{}:{}:{}:{:.2}", "1/1", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[0]);
                }
                rd.qual = snp.variant_quality as i32;
                rd.filter = "dn".to_string().into_bytes();
                rd.info = format!("RDS={}", "dense_snp").to_string().into_bytes();
                rd.format = "GT:GQ:DP:AF".to_string().into_bytes();
                records.push(rd);
                continue;
            }

            if snp.variant_type == 2 {
                let mut rd: VCFRecord = VCFRecord::default();
                rd.chromosome = snp.chromosome.clone();
                rd.position = snp.pos as u64 + 1;   // position in vcf format is 1-based
                rd.reference = vec![snp.reference as u8];
                rd.id = vec!['.' as u8];


                rd.alternative = vec![vec![snp.alleles[0] as u8]];
                rd.qual = snp.variant_quality as i32;
                rd.genotype = format!("{}:{}:{}:{:.2}", "1/1", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[0]);
                if snp.variant_quality < min_qual as f64 {
                    rd.filter = "LowQual".to_string().into_bytes();
                } else {
                    rd.filter = "PASS".to_string().into_bytes();
                }
                rd.info = "RDS=.".to_string().into_bytes();
                rd.format = "GT:GQ:DP:AF".to_string().into_bytes();
                records.push(rd);
            } else if snp.variant_type == 1 {
                let mut rd: VCFRecord = VCFRecord::default();
                rd.chromosome = snp.chromosome.clone();
                rd.position = snp.pos as u64 + 1;   // position in vcf format is 1-based
                rd.reference = vec![snp.reference as u8];
                rd.id = vec!['.' as u8];
                if snp.alleles[0] == snp.reference {
                    rd.alternative = vec![vec![snp.alleles[1] as u8]];
                    rd.qual = snp.variant_quality as i32;
                    rd.genotype = format!("{}:{}:{}:{:.2}", "0/1", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[1]);
                } else if snp.alleles[1] == snp.reference {
                    rd.alternative = vec![vec![snp.alleles[0] as u8]];
                    rd.qual = snp.variant_quality as i32;
                    rd.genotype = format!("{}:{}:{}:{:.2}", "0/1", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[0]);
                } else {
                    rd.alternative = vec![vec![snp.alleles[0] as u8], vec![snp.alleles[1] as u8]];
                    rd.qual = snp.variant_quality as i32;
                    rd.genotype = format!("{}:{}:{}:{:.2},{:.2}", "1/2", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[0], snp.allele_freqs[1]);
                }
                if snp.variant_quality < min_qual as f64 {
                    rd.filter = "LowQual".to_string().into_bytes();
                } else {
                    rd.filter = "PASS".to_string().into_bytes();
                }
                rd.info = "RDS=.".to_string().into_bytes();
                rd.format = "GT:GQ:DP:AF".to_string().into_bytes();
                records.push(rd);
            } else {
                println!("Unknown variant type: {:?}:{:?}", snp.chromosome, snp.pos);
                continue;
            }
        }
        return records;
    }
}

// fn parse_fai(fai_path: &str) -> Vec<(String, u32)> {
//     let mut contig_lengths: Vec<(String, u32)> = Vec::new();
//     let file = File::open(fai_path).unwrap();
//     let reader = BufReader::new(file);
//     for r in reader.lines() {
//         let line = r.unwrap().clone();
//         let parts: Vec<&str> = line.split('\t').collect();
//         contig_lengths.push((parts[0].to_string(), parts[1].parse().unwrap()));
//     }
//     return contig_lengths;
// }

// pub fn multithread_phase_maxcut(bam_file: String, ref_file: String, vcf_file: String, thread_size: usize, isolated_regions: Vec<Region>) {
//     let pool = rayon::ThreadPoolBuilder::new().num_threads(thread_size - 1).build().unwrap();
//     let vcf_records_queue = Mutex::new(VecDeque::new());
//     let bam: bam::IndexedReader = bam::IndexedReader::from_path(&bam_file).unwrap();
//     // let header = bam::Header::from_template(bam.header());
//     // let mut bam_writer = Mutex::new(bam::Writer::from_path(&out_bam, &header, Format::Bam).unwrap());
//     let ref_seqs = load_reference(ref_file.clone());
//     let fai_path = ref_file + ".fai";
//     if fs::metadata(&fai_path).is_err() {
//         panic!("Reference index file .fai does not exist.");
//     }
//     let contig_lengths = parse_fai(fai_path.as_str());
//
//     // multiplethreads for low coverage regions
//     pool.install(|| {
//         isolated_regions.par_iter().for_each(|reg| {
//             println!("Start {:?}", reg);
//             let mut profile = Profile::default();
//             let mut readnames: Vec<String> = Vec::new();
//             profile.init_with_pileup(&bam_file.as_str(), &reg);
//             profile.append_reference(&ref_seqs);
//             let mut snpfrag = SNPFrag::default();
//             snpfrag.get_candidate_snps(&profile, 0.3, 0.01, 10, 0.8, 0.9);
//             if snpfrag.snps.len() == 0 { return; }
//             for snp in snpfrag.snps.iter() {
//                 println!("snp: {:?}", snp);
//             }
//             snpfrag.get_fragments(&bam_file, &reg);
//             unsafe { snpfrag.init_haplotypes(); }
//             snpfrag.optimization_using_maxcut();
//             let vcf_records = snpfrag.output_vcf();
//             {
//                 let mut queue = vcf_records_queue.lock().unwrap();
//                 for rd in vcf_records.iter() {
//                     queue.push_back(rd.clone());
//                 }
//             }
//         });
//     });
//
//     let mut vf = File::create(vcf_file).unwrap();
//     vf.write("##fileformat=VCFv4.3\n".as_bytes()).unwrap();
//     vf.write("##FILTER=<ID=PASS,Description=\"All filters passed\">\n".as_bytes()).unwrap();
//     for ctglen in contig_lengths.iter() {
//         let chromosome = ctglen.0.clone();
//         let chromosome_len = ctglen.1.clone();
//         vf.write(format!("##contig=<ID={},length={}>\n", chromosome, chromosome_len).as_bytes()).unwrap();
//     }
//     vf.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n".as_bytes()).unwrap();
//     vf.write("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n".as_bytes()).unwrap();
//     vf.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n".as_bytes()).unwrap();
//     vf.write("##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n".as_bytes()).unwrap();
//     vf.write("##FORMAT=<ID=PQ,Number=1,Type=Integer,Description=\"Phasing Quality\">\n".as_bytes()).unwrap();
//     vf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample\n".as_bytes()).unwrap();
//     for rd in vcf_records_queue.lock().unwrap().iter() {
//         if rd.alternative.len() == 1 {
//             vf.write(format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n", std::str::from_utf8(&rd.chromosome).unwrap(),
//                              rd.position,
//                              std::str::from_utf8(&rd.id).unwrap(),
//                              std::str::from_utf8(&rd.reference).unwrap(),
//                              std::str::from_utf8(&rd.alternative[0]).unwrap(),
//                              rd.qual,
//                              std::str::from_utf8(&rd.filter).unwrap(),
//                              std::str::from_utf8(&rd.info).unwrap(),
//                              std::str::from_utf8(&rd.format).unwrap(),
//                              rd.genotype).as_bytes()).unwrap();
//         } else if rd.alternative.len() == 2 {
//             vf.write(format!("{}\t{}\t{}\t{}\t{},{}\t{}\t{}\t{}\t{}\t{}\n", std::str::from_utf8(&rd.chromosome).unwrap(),
//                              rd.position,
//                              std::str::from_utf8(&rd.id).unwrap(),
//                              std::str::from_utf8(&rd.reference).unwrap(),
//                              std::str::from_utf8(&rd.alternative[0]).unwrap(),
//                              std::str::from_utf8(&rd.alternative[1]).unwrap(),
//                              rd.qual,
//                              std::str::from_utf8(&rd.filter).unwrap(),
//                              std::str::from_utf8(&rd.info).unwrap(),
//                              std::str::from_utf8(&rd.format).unwrap(),
//                              rd.genotype).as_bytes()).unwrap();
//         }
//     }
// }

pub fn multithread_phase_haplotag(bam_file: String,
                                  ref_file: String,
                                  vcf_file: String,
                                  phased_bam_file: String,
                                  thread_size: usize,
                                  isolated_regions: Vec<Region>,
                                  genotype_only: bool,
                                  platform: &String,
                                  min_mapq: u8,
                                  min_baseq: u8,
                                  diff_baseq: u8,
                                  min_allele_freq: f32,
                                  min_allele_freq_include_intron: f32,
                                  min_qual_for_candidate: u32,
                                  min_qual_for_singlesnp_rnaedit: u32,
                                  min_allele_cnt: u32,
                                  no_strand_bias: bool,
                                  strand_bias_threshold: f32,
                                  cover_strand_bias_threshold: f32,
                                  min_depth: u32,
                                  max_depth: u32,
                                  min_read_length: usize,
                                  distance_to_splicing_site: u32,
                                  window_size: u32,
                                  distance_to_read_end: u32,
                                  diff_distance_to_read_end: i64,
                                  polya_tail_len: u32,
                                  dense_win_size: u32,
                                  min_dense_cnt: u32,
                                  avg_dense_dist: f32,
                                  min_homozygous_freq: f32,
                                  min_phase_score: f32,
                                  max_enum_snps: usize,
                                  random_flip_fraction: f32,
                                  read_assignment_cutoff: f64,
                                  imbalance_allele_expression_cutoff: f32,
                                  output_phasing: bool,
                                  no_bam_output: bool) {
    let pool = rayon::ThreadPoolBuilder::new().num_threads(thread_size).build().unwrap();
    let vcf_records_queue = Mutex::new(VecDeque::new());
    let read_haplotag_queue = Mutex::new(VecDeque::new());
    // let bam: bam::IndexedReader = bam::IndexedReader::from_path(&bam_file).unwrap();
    // let header = bam::Header::from_template(bam.header());
    // let mut bam_writer = Mutex::new(bam::Writer::from_path(&out_bam, &header, Format::Bam).unwrap());
    let ref_seqs = load_reference(ref_file.clone());
    let fai_path = ref_file + ".fai";
    if fs::metadata(&fai_path).is_err() {
        panic!("Reference index file .fai does not exist.");
    }
    let contig_lengths = parse_fai(fai_path.as_str());

    // multiplethreads for low coverage regions
    pool.install(|| {
        isolated_regions.par_iter().for_each(|reg| {
            // println!("Start {:?}", reg);
            let mut profile = Profile::default();
            let ref_seq = ref_seqs.get(&reg.chr).unwrap();
            profile.init_with_pileup(&bam_file.as_str(), &reg, ref_seq, platform, min_mapq, min_baseq, min_read_length, min_depth, max_depth, distance_to_read_end, polya_tail_len);
            // profile.append_reference(&ref_seqs);
            let mut snpfrag = SNPFrag::default();
            snpfrag.region = reg.clone();
            snpfrag.get_candidate_snps(&profile, min_allele_freq, min_allele_freq_include_intron, min_qual_for_candidate, min_depth, max_depth, min_baseq, min_homozygous_freq, no_strand_bias, strand_bias_threshold, cover_strand_bias_threshold, distance_to_splicing_site, window_size, distance_to_read_end, diff_distance_to_read_end, diff_baseq, dense_win_size, min_dense_cnt, avg_dense_dist);
            // for snp in snpfrag.candidate_snps.iter() {
            //     println!("snp: {:?}", snp);
            // }
            snpfrag.get_fragments(&bam_file, &reg);
            // snpfrag.filter_fp_snps(strand_bias_threshold, None);
            if genotype_only {
                // without phasing
                let vcf_records = snpfrag.output_vcf(min_qual_for_singlesnp_rnaedit);
                {
                    let mut queue = vcf_records_queue.lock().unwrap();
                    for rd in vcf_records.iter() {
                        queue.push_back(rd.clone());
                    }
                }
            } else {
                if snpfrag.hete_snps.len() > 0 {
                    unsafe { snpfrag.init_haplotypes(); }
                    unsafe { snpfrag.init_assignment(); }
                    // snpfrag.keep_reliable_snps_in_component();
                    snpfrag.phase(max_enum_snps, random_flip_fraction);
                    let read_assignments = snpfrag.assign_reads(read_assignment_cutoff);
                    // println!("{}:{}-{} hap1:", snpfrag.region.chr, snpfrag.region.start, snpfrag.region.end);
                    // for frag in snpfrag.fragments.iter() {
                    //     if frag.haplotag == -1 {
                    //         println!("{:?}", frag.exons);
                    //     }
                    // }
                    // println!("{}:{}-{} hap2:", snpfrag.region.chr, snpfrag.region.start, snpfrag.region.end);
                    // for frag in snpfrag.fragments.iter() {
                    //     if frag.haplotag == 1 {
                    //         println!("{:?}", frag.exons);
                    //     }
                    // }
                    snpfrag.add_phase_score(min_allele_cnt, imbalance_allele_expression_cutoff);
                    {
                        let mut queue = read_haplotag_queue.lock().unwrap();
                        for a in read_assignments.iter() {
                            queue.push_back((a.0.clone(), a.1.clone()));
                        }
                    }
                }
                let vcf_records = snpfrag.phased_output_vcf(min_phase_score, min_homozygous_freq, output_phasing, min_qual_for_candidate, min_qual_for_singlesnp_rnaedit);
                {
                    let mut queue = vcf_records_queue.lock().unwrap();
                    for rd in vcf_records.iter() {
                        queue.push_back(rd.clone());
                    }
                }
            }
        });
    });

    println!("{} Start writing VCF file.", Local::now().format("%Y-%m-%d %H:%M:%S"));
    let mut vf = File::create(vcf_file).unwrap();
    vf.write("##fileformat=VCFv4.3\n".as_bytes()).unwrap();
    for ctglen in contig_lengths.iter() {
        let chromosome = ctglen.0.clone();
        let chromosome_len = ctglen.1.clone();
        vf.write(format!("##contig=<ID={},length={}>\n", chromosome, chromosome_len).as_bytes()).unwrap();
    }
    vf.write("##FILTER=<ID=PASS,Description=\"All filters passed\">\n".as_bytes()).unwrap();
    vf.write("##FILTER=<ID=LowQual,Description=\"Low phasing quality\">\n".as_bytes()).unwrap();
    vf.write("##FILTER=<ID=RnaEdit,Description=\"RNA editing\">\n".as_bytes()).unwrap();
    vf.write("##FILTER=<ID=dn,Description=\"Dense cluster of variants\">\n".as_bytes()).unwrap();
    vf.write("##INFO=<ID=RDS,Number=1,Type=String,Description=\"RNA editing or Dense SNP or Single SNP.\">\n".as_bytes()).unwrap();
    vf.write("##INFO=<ID=AlleleImbalance,Number=0,Type=Flag,Description=\"Imbalanced allele expression.\">\n".as_bytes()).unwrap();
    vf.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n".as_bytes()).unwrap();
    vf.write("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n".as_bytes()).unwrap();
    vf.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n".as_bytes()).unwrap();
    vf.write("##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n".as_bytes()).unwrap();
    vf.write("##FORMAT=<ID=PQ,Number=1,Type=Float,Description=\"Phasing Quality\">\n".as_bytes()).unwrap();
    vf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample\n".as_bytes()).unwrap();

    for rd in vcf_records_queue.lock().unwrap().iter() {
        if rd.alternative.len() == 1 {
            vf.write(format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n", std::str::from_utf8(&rd.chromosome).unwrap(),
                             rd.position,
                             std::str::from_utf8(&rd.id).unwrap(),
                             std::str::from_utf8(&rd.reference).unwrap(),
                             std::str::from_utf8(&rd.alternative[0]).unwrap(),
                             rd.qual,
                             std::str::from_utf8(&rd.filter).unwrap(),
                             std::str::from_utf8(&rd.info).unwrap(),
                             std::str::from_utf8(&rd.format).unwrap(),
                             rd.genotype).as_bytes()).unwrap();
        } else if rd.alternative.len() == 2 {
            vf.write(format!("{}\t{}\t{}\t{}\t{},{}\t{}\t{}\t{}\t{}\t{}\n", std::str::from_utf8(&rd.chromosome).unwrap(),
                             rd.position,
                             std::str::from_utf8(&rd.id).unwrap(),
                             std::str::from_utf8(&rd.reference).unwrap(),
                             std::str::from_utf8(&rd.alternative[0]).unwrap(),
                             std::str::from_utf8(&rd.alternative[1]).unwrap(),
                             rd.qual,
                             std::str::from_utf8(&rd.filter).unwrap(),
                             std::str::from_utf8(&rd.info).unwrap(),
                             std::str::from_utf8(&rd.format).unwrap(),
                             rd.genotype).as_bytes()).unwrap();
        }
    }
    drop(vf);
    println!("{} Finish writing VCF file.", Local::now().format("%Y-%m-%d %H:%M:%S"));

    if !no_bam_output {
        println!("{} Start writing phased BAM file.", Local::now().format("%Y-%m-%d %H:%M:%S"));
        let mut read_assignments: HashMap<String, i32> = HashMap::new();
        for rd in read_haplotag_queue.lock().unwrap().iter() {
            read_assignments.insert(rd.0.clone(), rd.1.clone());
        }
        let mut bam_reader = bam::IndexedReader::from_path(&bam_file).unwrap();
        let header = bam::Header::from_template(&bam_reader.header());
        let mut bam_writer = bam::Writer::from_path(phased_bam_file, &header, Format::Bam).unwrap();
        bam_writer.set_threads(thread_size).unwrap();
        for region in isolated_regions.iter() {
            bam_reader.fetch((region.chr.as_str(), region.start, region.end)).unwrap(); // set region
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
        }
        drop(bam_writer);
        println!("{} Finish writing phased BAM file.", Local::now().format("%Y-%m-%d %H:%M:%S"));
    }
}

