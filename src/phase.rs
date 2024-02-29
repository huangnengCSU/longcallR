use std::collections::{HashMap, HashSet, VecDeque};
use std::{fs};
use std::fs::File;
use std::hash::Hash;
use std::io::{Write};
use crate::util::{Region, independent_test, parse_fai, Profile, load_reference, VCFRecord};
use rust_htslib::{bam, bam::Read, bam::record::Record, bam::Format, bam::record::Aux};
use rust_htslib::htslib::{drand48};
use std::sync::{Mutex};
use bio::bio_types::strand::ReqStrand::Forward;
use rayon::prelude::*;
use rand::seq::SliceRandom;
use std::cmp::{Ordering};
use crate::Platform;
use probability::distribution::{Binomial, Distribution};


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
    pub haplotype_expression: [u32; 2],
    // [allele1, allele2], the number of reads supporting two alleles expreessing on two haplotype
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
    // base allele on alphabet  {-1, 1, 0}, 1: base==alleles[0], -1: base==alleles[1], 0: not covered (bases except allele1 and allele2, deletions or N)
    pub prob: f64,
    // error rate of observe current base
}

#[derive(Debug, Clone, Default, PartialEq, Eq, Hash)]
pub struct Exon {
    pub chr: String,
    pub start: i64,
    // start position on the reference, 0-based, inclusive
    pub end: i64,
    // end position on the reference, 0-based, exclusive
    pub state: u8, // 0: start exon, 1: internal exon, 2: end exon, 3: whole read is a single exon
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
    pub exons: Vec<Exon>,
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
    pub min_linkers: i32,
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
            if bf.i {
                continue;
            }

            // 1.filtering with depth
            let depth = bf.get_depth_exclude_intron_deletion();
            if depth < min_coverage {
                position += 1;
                continue;
            }
            if depth > max_coverage {
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
                    position += 1;
                    continue;
                }
            }

            // filtering with depth, considering intron reads
            let depth_include_intron = bf.get_depth_include_intron();
            if (allele1_cnt as f32) / (depth_include_intron as f32) < min_allele_freq_include_intron {
                // maybe caused by erroneous intron alignment
                position += 1;
                continue;
            }

            // filtering with deletion frequency
            if bf.d > allele1_cnt {
                position += 1;
                continue;
            }

            if !no_strand_bias {

                // filtering snps only covered by one strand reads (may caused by intron alignment error)
                let total_cover_cnt = bf.forward_cnt + bf.backward_cnt;   // does not include intron reads
                if bf.forward_cnt as f32 / total_cover_cnt as f32 > cover_strand_bias_threshold || bf.backward_cnt as f32 / total_cover_cnt as f32 > cover_strand_bias_threshold {
                    position += 1;
                    continue;
                }

                // filtering with strand bias
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
                            strand_bias = true;
                            continue;
                        }
                    }
                    if strand_bias {
                        position += 1;
                        continue;
                    }
                }

                // // use fisher's test or chi-square test to test strand bias
                // let mut allele_freq_mat = [0; 4];
                // match allele1 {
                //     'a' => {
                //         allele_freq_mat[0] = bf.base_strands.a[0];
                //         allele_freq_mat[1] = bf.base_strands.a[1];
                //     }
                //     'A' => {
                //         allele_freq_mat[0] = bf.base_strands.a[0];
                //         allele_freq_mat[1] = bf.base_strands.a[1];
                //     }
                //     'c' => {
                //         allele_freq_mat[0] = bf.base_strands.c[0];
                //         allele_freq_mat[1] = bf.base_strands.c[1];
                //     }
                //     'C' => {
                //         allele_freq_mat[0] = bf.base_strands.c[0];
                //         allele_freq_mat[1] = bf.base_strands.c[1];
                //     }
                //     'g' => {
                //         allele_freq_mat[0] = bf.base_strands.g[0];
                //         allele_freq_mat[1] = bf.base_strands.g[1];
                //     }
                //     'G' => {
                //         allele_freq_mat[0] = bf.base_strands.g[0];
                //         allele_freq_mat[1] = bf.base_strands.g[1];
                //     }
                //     't' => {
                //         allele_freq_mat[0] = bf.base_strands.t[0];
                //         allele_freq_mat[1] = bf.base_strands.t[1];
                //     }
                //     'T' => {
                //         allele_freq_mat[0] = bf.base_strands.t[0];
                //         allele_freq_mat[1] = bf.base_strands.t[1];
                //     }
                //     _ => {
                //         println!("Error: unknown allele");
                //         position += 1;
                //         continue;
                //     }
                // }
                // match allele2 {
                //     'a' => {
                //         allele_freq_mat[2] = bf.base_strands.a[0];
                //         allele_freq_mat[3] = bf.base_strands.a[1];
                //     }
                //     'A' => {
                //         allele_freq_mat[2] = bf.base_strands.a[0];
                //         allele_freq_mat[3] = bf.base_strands.a[1];
                //     }
                //     'c' => {
                //         allele_freq_mat[2] = bf.base_strands.c[0];
                //         allele_freq_mat[3] = bf.base_strands.c[1];
                //     }
                //     'C' => {
                //         allele_freq_mat[2] = bf.base_strands.c[0];
                //         allele_freq_mat[3] = bf.base_strands.c[1];
                //     }
                //     'g' => {
                //         allele_freq_mat[2] = bf.base_strands.g[0];
                //         allele_freq_mat[3] = bf.base_strands.g[1];
                //     }
                //     'G' => {
                //         allele_freq_mat[2] = bf.base_strands.g[0];
                //         allele_freq_mat[3] = bf.base_strands.g[1];
                //     }
                //     't' => {
                //         allele_freq_mat[2] = bf.base_strands.t[0];
                //         allele_freq_mat[3] = bf.base_strands.t[1];
                //     }
                //     'T' => {
                //         allele_freq_mat[2] = bf.base_strands.t[0];
                //         allele_freq_mat[3] = bf.base_strands.t[1];
                //     }
                //     _ => {
                //         println!("Error: unknown allele");
                //         position += 1;
                //         continue;
                //     }
                // }
                //
                // if (allele_freq_mat[0] + allele_freq_mat[2]) > 0 && (allele_freq_mat[1] + allele_freq_mat[3]) > 0 {
                //     let phred_pvalue = independent_test([allele_freq_mat[0] as u32, allele_freq_mat[1] as u32, allele_freq_mat[2] as u32, allele_freq_mat[3] as u32]);
                //     if phred_pvalue > 100.0 {
                //         position += 1;
                //         continue;
                //     }
                // } else {
                //     // filtering with strand bias
                //     let mut variant_allele: Vec<char> = Vec::new();
                //     // if allele1 != bf.ref_base {
                //     if allele1 != bf.ref_base && allele1_cnt >= 4 {
                //         variant_allele.push(allele1);
                //     }
                //     // if allele2 != bf.ref_base {
                //     if allele2 != bf.ref_base && allele2_cnt >= 4 {
                //         variant_allele.push(allele2);
                //     }
                //     if variant_allele.len() > 0 {
                //         // variant allele
                //         let mut strand_bias = false;
                //         for allele_base in variant_allele.iter() {
                //             let mut fcnt = 0;
                //             let mut bcnt = 0;
                //             match allele_base {
                //                 'a' => { [fcnt, bcnt] = bf.base_strands.a; }
                //                 'A' => { [fcnt, bcnt] = bf.base_strands.a; }
                //                 'c' => { [fcnt, bcnt] = bf.base_strands.c; }
                //                 'C' => { [fcnt, bcnt] = bf.base_strands.c; }
                //                 'g' => { [fcnt, bcnt] = bf.base_strands.g; }
                //                 'G' => { [fcnt, bcnt] = bf.base_strands.g; }
                //                 't' => { [fcnt, bcnt] = bf.base_strands.t; }
                //                 'T' => { [fcnt, bcnt] = bf.base_strands.t; }
                //                 _ => {
                //                     println!("Error: unknown allele");
                //                     position += 1;
                //                     continue;
                //                 }
                //             }
                //             let total_cnt = fcnt + bcnt;
                //             if fcnt as f32 / total_cnt as f32 > strand_bias_threshold || bcnt as f32 / total_cnt as f32 > strand_bias_threshold {
                //                 if allele1 != bf.ref_base {
                //                     println!("strand bias: {}:{}", profile.region.chr, position);
                //                 }
                //                 // strand bias
                //                 strand_bias = true;
                //                 continue;
                //             }
                //         }
                //         if strand_bias {
                //             position += 1;
                //             continue;
                //         }
                //     }
                // }
            }

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
                let local_error_rate = rbf.get_none_ref_count() as f32 / (rbf.a + rbf.c + rbf.g + rbf.t + rbf.d) as f32;
                local_misalignment_ratio.push(local_error_rate);
                rext += 1;
                ridx += 1;
            }

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
                position += 1;
                continue;
            }

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
            let mut variant_prob = logprob.clone();
            variant_prob[0] = 10.0_f64.powf(logprob[0]);
            variant_prob[1] = 10.0_f64.powf(logprob[1]);
            variant_prob[2] = 10.0_f64.powf(logprob[2]);
            let sum_variant_prob = variant_prob[0] + variant_prob[1] + variant_prob[2];
            // println!("3:{}:{},{:?}", position, sum_variant_prob, variant_prob);
            variant_prob = [variant_prob[0] / sum_variant_prob, variant_prob[1] / sum_variant_prob, variant_prob[2] / sum_variant_prob];
            // println!("4:{}:{},{:?}", position, correction_factor, variant_prob);
            // QUAL phred-scaled quality score for the assertion made in ALT. i.e. give -10log_10 prob(call in ALT is wrong).
            // If ALT is `.` (no variant) then this is -10log_10 p(variant), and if ALT is not `.` this is -10log_10 p(no variant).
            let variant_quality = -10.0 * ((10e-301_f64.max(variant_prob[2])).log10());    // if variant_quality is greater than 3000, we set it to 3000

            // calculate GQ: The value of GQ is simply the difference between the second lowest PL and the lowest PL (which is always 0, normalized PL)
            let mut log10_likelihood = loglikelihood.clone();
            let max_log10_likelihood = log10_likelihood[0].max(log10_likelihood[1]).max(log10_likelihood[2]);
            log10_likelihood[0] = 10.0_f64.powf(log10_likelihood[0] - max_log10_likelihood);
            log10_likelihood[1] = 10.0_f64.powf(log10_likelihood[1] - max_log10_likelihood);
            log10_likelihood[2] = 10.0_f64.powf(log10_likelihood[2] - max_log10_likelihood);
            let sum_log10_likelihood = log10_likelihood[0] + log10_likelihood[1] + log10_likelihood[2];
            let mut genotype_prob = [log10_likelihood[0] / sum_log10_likelihood, log10_likelihood[1] / sum_log10_likelihood, log10_likelihood[2] / sum_log10_likelihood];
            // println!("{}:{},{},{},{}", profile.region.chr, position, genotype_prob[0], genotype_prob[1], genotype_prob[2]);
            let mut phred_genotype_prob = [0.0, 0.0, 0.0];
            phred_genotype_prob[0] = -10.0 * genotype_prob[0].log10();    // phred scale likelihood of genotype: 1/1
            phred_genotype_prob[1] = -10.0 * genotype_prob[1].log10();    // phred scale likelihood of genotype: 0/1
            phred_genotype_prob[2] = -10.0 * genotype_prob[2].log10();    // phred scale likelihood of genotype: 0/0
            phred_genotype_prob.sort_by(cmp_f64);
            let genotype_quality = phred_genotype_prob[1] - phred_genotype_prob[0];

            // filter low variant quality
            if variant_quality < min_qual_for_candidate as f64 {
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
                    position += 1;
                    continue;
                } else if allele2 != bf.ref_base && allele2_freq < min_allele_freq {
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
            if *v as i32 >= self.min_linkers {
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
        assert!(self.min_linkers > 0, "Error: min_linkers <= 0");
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
                match cg.char() as u8 {
                    b'S' | b'H' => {
                        continue;
                    }
                    b'M' | b'X' | b'=' => {
                        for _ in 0..cg.len() {
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
                                // filtered SNP and rna editing site will not be used for haplotype phasing, not covered SNP will not be used for haplotype phasing
                                if self.candidate_snps[frag_elem.snp_idx].filter == false && self.candidate_snps[frag_elem.snp_idx].rna_editing == false && frag_elem.p != 0 {
                                    fragment.list.push(frag_elem);
                                }
                                idx += 1;
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
                            if pos_on_ref == snp_pos {
                                let mut frag_elem = FragElem::default();
                                frag_elem.snp_idx = self.hete_snps[idx];
                                frag_elem.pos = pos_on_ref;
                                frag_elem.base = '-';
                                frag_elem.baseq = 0;
                                frag_elem.strand = strand;
                                frag_elem.p = 0;
                                idx += 1;
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
                        if fragment.exons.len() == 0 {
                            fragment.exons.push(Exon { chr: region.chr.clone(), start: exon_start, end: exon_end, state: 0 });  // start exon
                        } else {
                            fragment.exons.push(Exon { chr: region.chr.clone(), start: exon_start, end: exon_end, state: 1 });  // internal exon
                        }
                        exon_start = -1;
                        exon_end = -1;
                        for _ in 0..cg.len() {
                            if pos_on_ref == snp_pos {
                                let mut frag_elem = FragElem::default();
                                frag_elem.snp_idx = self.hete_snps[idx];
                                frag_elem.pos = pos_on_ref;
                                frag_elem.base = b'-' as char;
                                frag_elem.baseq = 0;
                                frag_elem.strand = strand;
                                frag_elem.p = 0;
                                idx += 1;
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
                if fragment.exons.len() > 0 {
                    fragment.exons.push(Exon { chr: region.chr.clone(), start: exon_start, end: exon_end, state: 2 });  // end exon
                } else {
                    fragment.exons.push(Exon { chr: region.chr.clone(), start: exon_start, end: exon_end, state: 3 });  // single exon
                }
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
            if link_hete_cnt >= self.min_linkers {
                for fe in fragment.list.iter() {
                    // record each snp cover by which fragments
                    self.candidate_snps[fe.snp_idx].snp_cover_fragments.push(fragment.fragment_idx);
                }
                self.fragments.push(fragment);
            }
        }
    }

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
                q1 = q1 * (1.0 - probs[i]);
            } else {
                q1 = q1 * probs[i];
            }
        }

        for i in 0..delta.len() {
            if delta[i] == ps[i] {
                q2 = q2 * (1.0 - probs[i]);
                q3 = q3 * probs[i];
            } else {
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
                q1 = q1 * (1.0 - probs[k]);
            } else {
                q1 = q1 * probs[k];
            }
        }

        for k in 0..sigma.len() {
            if sigma[k] == ps[k] {
                q2 = q2 * (1.0 - probs[k]);
                q3 = q3 * probs[k];
            } else {
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
                q1 = q1 + (1.0 - probs[k]).log10();
            } else {
                q1 = q1 + probs[k].log10();
            }
        }

        for k in 0..sigma.len() {
            if sigma[k] == ps[k] {
                q2 = q2 + (1.0 - probs[k]).log10();
                q3 = q3 + probs[k].log10();
            } else {
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
        let mut rv = 0;
        if logp > pre_logp { rv = 1; } else if logp == pre_logp { rv = 0; } else {
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

        let mut phasing_increase: bool = true;
        let mut haplotag_increase: bool = true;
        let mut num_iters = 0;
        let mut fragments_cnt: HashMap<usize, usize> = HashMap::new();
        let mut covered_fragments: HashSet<usize> = HashSet::new();
        let mut processed_snps: HashSet<usize> = HashSet::new();
        for i in self.hete_snps.iter() {
            processed_snps.insert(*i);
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
            if *v as i32 >= self.min_linkers {
                covered_fragments.insert(*k);
            }
        }

        while phasing_increase | haplotag_increase {
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
                    tmp_haplotag.insert(*k, sigma_k * (-1));
                } else {
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
                } else {
                    tmp_haplotype.insert(*i, delta_i);
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
                break;
            }
        }

        let prob = SNPFrag::cal_overall_probability(&self, &processed_snps, &covered_fragments);
        return prob;
    }

    pub fn phase(&mut self, max_enum_snps: usize, random_flip_fraction: f32, max_iters: i32) {
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
            if *v as i32 >= self.min_linkers {
                covered_fragments.insert(*k);
            }
        }


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
            let mut max_iter: i32 = max_iters;
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
                phase_score = 0.0;
            } else {
                if sigma.len() > 0 {
                    phase_score = -10.0_f64 * (1.0 - SNPFrag::cal_log_delta_sigma(delta_i, &sigma, &ps, &probs)).log10();   // calaulate assignment score
                } else {
                    phase_score = 0.0;  // all reads belong to unknown group, maybe caused by single SNP
                }
            }


            // use binomial test to check allele imbalance expression
            let mut hap1_allele_cnt: [u32; 2] = [0, 0];
            let mut hap2_allele_cnt: [u32; 2] = [0, 0];

            for k in 0..sigma.len() {
                if sigma[k] == 1 {
                    // hap1
                    if ps[k] == 1 {
                        hap1_allele_cnt[0] += 1;    // hap1 allele1
                    } else if ps[k] == -1 {
                        hap1_allele_cnt[1] += 1;    // hap1 allele2
                    }
                } else if sigma[k] == -1 {
                    // hap2
                    if ps[k] == 1 {
                        hap2_allele_cnt[0] += 1;    // hap2 allele1
                    } else if ps[k] == -1 {
                        hap2_allele_cnt[1] += 1;    // hap2 allele2
                    }
                }
            }

            // let mut observed_cnt = 0;    // the number of observation
            // let mut trails_cnt = 0;  // the number of trials
            // if (hap1_allele_cnt[0] + hap2_allele_cnt[0]) > (hap1_allele_cnt[1] + hap2_allele_cnt[1]) {
            //     if hap1_allele_cnt[0] > hap2_allele_cnt[0] {
            //         observed_cnt = hap2_allele_cnt[1];
            //         trails_cnt += hap2_allele_cnt[1] + hap1_allele_cnt[0];
            //     } else {
            //         observed_cnt = hap1_allele_cnt[1];
            //         trails_cnt += hap1_allele_cnt[1] + hap2_allele_cnt[0];
            //     }
            // } else {
            //     if hap1_allele_cnt[1] > hap2_allele_cnt[1] {
            //         observed_cnt = hap2_allele_cnt[0];
            //         trails_cnt += hap2_allele_cnt[0] + hap1_allele_cnt[1];
            //     } else {
            //         observed_cnt = hap1_allele_cnt[0];
            //         trails_cnt += hap1_allele_cnt[0] + hap2_allele_cnt[1];
            //     }
            // }
            // let mut allele_imbalance: bool = false;
            // let binom = Binomial::new(trails_cnt as usize, 0.5);
            // let less_p = binom.distribution(observed_cnt as f64);
            // let phred_pvalue = -10.0_f64 * less_p.log10();
            // if phred_pvalue > 50.0 { allele_imbalance = true; }
            // self.candidate_snps[ti].allele_imbalance = allele_imbalance;

            let mut expression_cnt = [0; 2];
            if (hap1_allele_cnt[0] + hap2_allele_cnt[0]) > (hap1_allele_cnt[1] + hap2_allele_cnt[1]) {
                if hap1_allele_cnt[0] > hap2_allele_cnt[0] {
                    expression_cnt[0] = hap1_allele_cnt[0];
                    expression_cnt[1] = hap2_allele_cnt[1];
                } else {
                    expression_cnt[0] = hap2_allele_cnt[0];
                    expression_cnt[1] = hap1_allele_cnt[1];
                }
            } else {
                if hap1_allele_cnt[1] > hap2_allele_cnt[1] {
                    expression_cnt[0] = hap2_allele_cnt[0];
                    expression_cnt[1] = hap1_allele_cnt[1];
                } else {
                    expression_cnt[0] = hap1_allele_cnt[0];
                    expression_cnt[1] = hap2_allele_cnt[1];
                }
            }
            self.candidate_snps[ti].haplotype_expression = expression_cnt;
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
                if snp.phase_score < min_phase_score as f64 && snp.allele_freqs[0] > min_homozygous_freq && snp.alleles[0] != snp.reference && snp.filter == false {
                    let mut rd: VCFRecord = VCFRecord::default();
                    rd.chromosome = snp.chromosome.clone();
                    rd.position = snp.pos as u64 + 1;   // position in vcf format is 1-based
                    rd.reference = vec![snp.reference as u8];
                    rd.id = vec!['.' as u8];
                    rd.alternative = vec![vec![snp.alleles[0] as u8]];
                    rd.qual = snp.variant_quality as i32;
                    rd.genotype = format!("{}:{}:{}:{:.2}:{:.2}", "1/1", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[0], snp.phase_score);
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
                    rd.format = "GT:GQ:DP:AF:PQ".to_string().into_bytes();
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
                                    rd.genotype = format!("{}:{}:{}:{:.2}:{:.2}:{},{}", "0|1", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[1], snp.phase_score, snp.haplotype_expression[0], snp.haplotype_expression[1]);
                                } else {
                                    rd.genotype = format!("{}:{}:{}:{:.2}:{:.2}:{},{}", "1|0", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[1], snp.phase_score, snp.haplotype_expression[0], snp.haplotype_expression[1]);
                                }
                            } else if snp.alleles[1] == snp.reference {
                                rd.alternative = vec![vec![snp.alleles[0] as u8]];
                                rd.qual = snp.variant_quality as i32;
                                if hp == -1 {
                                    rd.genotype = format!("{}:{}:{}:{:.2}:{:.2}:{},{}", "0|1", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[0], snp.phase_score, snp.haplotype_expression[0], snp.haplotype_expression[1]);
                                } else {
                                    rd.genotype = format!("{}:{}:{}:{:.2}:{:.2}:{},{}", "1|0", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[0], snp.phase_score, snp.haplotype_expression[0], snp.haplotype_expression[1]);
                                }
                            } else {
                                rd.alternative = vec![vec![snp.alleles[0] as u8], vec![snp.alleles[1] as u8]];
                                rd.qual = snp.variant_quality as i32;
                                rd.genotype = format!("{}:{}:{}:{:.2},{:.2}:{:.2}:{},{}", "1|2", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[0], snp.allele_freqs[1], snp.phase_score, snp.haplotype_expression[0], snp.haplotype_expression[1]);
                            }
                            if snp.phase_score < min_phase_score as f64 || snp.variant_quality < min_qual_for_candidate as f64 {
                                // not confident phase or low variant quality
                                rd.filter = "LowQual".to_string().into_bytes();
                            } else {
                                rd.filter = "PASS".to_string().into_bytes();
                            }
                            rd.info = "RDS=.".to_string().into_bytes();
                            rd.format = "GT:GQ:DP:AF:PQ:AE".to_string().into_bytes();
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

fn exon_cluster(mut exons: Vec<Exon>, smallest_start: i64, largest_end: i64, min_sup: i32) -> HashMap<Exon, Vec<Exon>> {
    let mut cover_vec = vec![0; (largest_end - smallest_start) as usize];
    let mut cover_exon_idx: Vec<Vec<usize>> = vec![Vec::new(); (largest_end - smallest_start) as usize];
    for idx in 0..exons.len() {
        let e = &exons[idx];
        for j in e.start..e.end {
            cover_vec[j as usize - smallest_start as usize] += 1;
        }
        cover_exon_idx[e.start as usize - smallest_start as usize].push(idx);
    }

    // for overlapped exons, merge them into one cluster
    let mut clustersI: Vec<Vec<Exon>> = Vec::new();
    let mut ts = -1;    // index in cover_vec
    let mut te = -1;    // index in cover_vec
    for i in 0..cover_vec.len() {
        if cover_vec[i] == 0 {
            if ts >= 0 && te >= 0 {
                let mut t_cluster: Vec<Exon> = Vec::new();
                for ti in (ts as usize)..=(te as usize) {
                    for idx in cover_exon_idx[ti].iter() {
                        t_cluster.push(exons[*idx].clone());
                    }
                }
                clustersI.push(t_cluster);
                ts = -1;
                te = -1;
            }
            continue;
        } else {
            if ts == -1 && te == -1 {
                ts = i as i64;
                te = i as i64;
            } else {
                te = i as i64;
            }
        }
    }
    if ts >= 0 && te >= 0 {
        let mut t_cluster: Vec<Exon> = Vec::new();
        for ti in (ts as usize)..=(te as usize) {
            for idx in cover_exon_idx[ti].iter() {
                t_cluster.push(exons[*idx].clone());
            }
        }
        clustersI.push(t_cluster);
    }

    // for each level I cluster, divide them into level II clusters (as hierarchical clustering)
    let mut clusters: HashMap<Exon, Vec<Exon>> = HashMap::new();
    for c in clustersI.iter() {
        let mut exon_hashmap: HashMap<Exon, Vec<Exon>> = HashMap::new();
        for e in c.iter() {
            let mut texon = e.clone();
            if e.state == 0 {
                // start exon ignores the start position when matching exons, because the start position of start exon is not consistent
                texon.start = -1;
            } else if e.state == 2 {
                // end exon ignores the end position when matching exons, because the end position of end exon is not consistent
                texon.end = -1;
            }
            if exon_hashmap.contains_key(&texon) {
                let exon_vec = exon_hashmap.get_mut(&texon).unwrap();
                exon_vec.push(e.clone());
            } else {
                exon_hashmap.insert(texon.clone(), vec![e.clone()]);
            }
        }

        // ignore exons without enough support
        let mut filter_exons: Vec<Exon> = Vec::new();
        for (e, v) in exon_hashmap.iter() {
            if (v.len() as i32) < min_sup {
                filter_exons.push(e.clone());
            }
        }
        for e in filter_exons.iter() {
            exon_hashmap.remove(e);
        }
        for (e, v) in exon_hashmap.iter() {
            clusters.insert(e.clone(), v.clone());
        }
    }
    return clusters;
}

pub fn multithread_phase_haplotag(bam_file: String,
                                  ref_file: String,
                                  vcf_file: String,
                                  phased_bam_file: String,
                                  thread_size: usize,
                                  isolated_regions: Vec<Region>,
                                  genotype_only: bool,
                                  platform: &Platform,
                                  max_iters: i32,
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
                                  min_linkers: i32,
                                  min_phase_score: f32,
                                  max_enum_snps: usize,
                                  random_flip_fraction: f32,
                                  read_assignment_cutoff: f64,
                                  imbalance_allele_expression_cutoff: f32,
                                  output_phasing: bool,
                                  no_bam_output: bool,
                                  haplotype_bam_output: bool,
                                  output_read_assignment: bool,
                                  haplotype_specific_exon: bool,
                                  min_sup_haplotype_exon: u32) {
    let pool = rayon::ThreadPoolBuilder::new().num_threads(thread_size).build().unwrap();
    let vcf_records_queue = Mutex::new(VecDeque::new());
    let read_haplotag1_queue = Mutex::new(VecDeque::new());
    let read_haplotag2_queue = Mutex::new(VecDeque::new());
    let read_haplotag_queue = Mutex::new(VecDeque::new());
    let haplotype_exon_queue = Mutex::new(VecDeque::new());
    let ref_seqs = load_reference(ref_file.clone());
    let fai_path = ref_file + ".fai";
    if fs::metadata(&fai_path).is_err() {
        panic!("Reference index file .fai does not exist.");
    }
    let contig_lengths = parse_fai(fai_path.as_str());
    let mut contig_order = Vec::new();
    for (k, _) in contig_lengths.iter() {
        contig_order.push(k.clone());
    }

    pool.install(|| {
        isolated_regions.par_iter().for_each(|reg| {
            let mut profile = Profile::default();
            let ref_seq = ref_seqs.get(&reg.chr).unwrap();
            profile.init_with_pileup(&bam_file.as_str(), &reg, ref_seq, platform, min_mapq, min_baseq, min_read_length, min_depth, max_depth, distance_to_read_end, polya_tail_len);
            let mut snpfrag = SNPFrag::default();
            snpfrag.region = reg.clone();
            snpfrag.min_linkers = min_linkers;
            snpfrag.get_candidate_snps(&profile, min_allele_freq, min_allele_freq_include_intron, min_qual_for_candidate, min_depth, max_depth, min_baseq, min_homozygous_freq, no_strand_bias, strand_bias_threshold, cover_strand_bias_threshold, distance_to_splicing_site, window_size, distance_to_read_end, diff_distance_to_read_end, diff_baseq, dense_win_size, min_dense_cnt, avg_dense_dist);
            snpfrag.get_fragments(&bam_file, &reg);
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
                    snpfrag.phase(max_enum_snps, random_flip_fraction, max_iters);
                    let read_assignments = snpfrag.assign_reads(read_assignment_cutoff);
                    snpfrag.add_phase_score(min_allele_cnt, imbalance_allele_expression_cutoff);
                    let mut haplotype_exons: Vec<(Exon, i32, i32)> = Vec::new();
                    {
                        if haplotype_bam_output || haplotype_specific_exon {
                            let mut queue1 = read_haplotag1_queue.lock().unwrap();
                            let mut queue2 = read_haplotag2_queue.lock().unwrap();
                            let mut hap1_read_count = 0;
                            let mut hap2_read_count = 0;
                            let mut haplotype_read_count_pass = false;
                            for a in read_assignments.iter() {
                                if *a.1 == 1 {
                                    hap1_read_count += 1;
                                } else if *a.1 == 2 {
                                    hap2_read_count += 1;
                                }
                                if hap1_read_count >= 10 && hap2_read_count >= 10 {
                                    haplotype_read_count_pass = true;
                                    break;
                                }
                            }
                            if haplotype_bam_output && haplotype_read_count_pass {
                                for a in read_assignments.iter() {
                                    if *a.1 == 1 {
                                        queue1.push_back(a.0.clone());
                                    } else if *a.1 == 2 {
                                        queue2.push_back(a.0.clone());
                                    }
                                }
                            }
                            if haplotype_specific_exon && haplotype_read_count_pass {
                                let mut hap1_exons: Vec<Exon> = Vec::new();
                                let mut hap2_exons: Vec<Exon> = Vec::new();
                                let mut hap1_smallest_start = 0;
                                let mut hap1_largest_end = 0;
                                let mut hap2_smallest_start = 0;
                                let mut hap2_largest_end = 0;
                                // collect exons
                                for frag in snpfrag.fragments.iter() {
                                    if frag.assignment == 1 {
                                        for e in frag.exons.iter() {
                                            hap1_exons.push(e.clone());
                                            if hap1_smallest_start == 0 || e.start < hap1_smallest_start {
                                                hap1_smallest_start = e.start;
                                            }
                                            if e.end > hap1_largest_end {
                                                hap1_largest_end = e.end;
                                            }
                                        }
                                    } else if frag.assignment == 2 {
                                        for e in frag.exons.iter() {
                                            hap2_exons.push(e.clone());
                                            if hap2_smallest_start == 0 || e.start < hap2_smallest_start {
                                                hap2_smallest_start = e.start;
                                            }
                                            if e.end > hap2_largest_end {
                                                hap2_largest_end = e.end;
                                            }
                                        }
                                    }
                                }
                                // consensus exons
                                let hap1_consensus_exons = exon_cluster(hap1_exons.clone(), hap1_smallest_start, hap1_largest_end, 0);
                                let hap2_consensus_exons = exon_cluster(hap2_exons.clone(), hap2_smallest_start, hap2_largest_end, 0);
                                let mut combined_consensus_exons: HashMap<Exon, (i32, i32)> = HashMap::new();
                                for (e, v) in hap1_consensus_exons.iter() {
                                    if combined_consensus_exons.contains_key(e) {
                                        let (c1, c2) = combined_consensus_exons.get_mut(e).unwrap();
                                        *c1 += v.len() as i32;
                                    } else {
                                        combined_consensus_exons.insert(e.clone(), (v.len() as i32, 0));
                                    }
                                }
                                for (e, v) in hap2_consensus_exons.iter() {
                                    if combined_consensus_exons.contains_key(e) {
                                        let (c1, c2) = combined_consensus_exons.get_mut(e).unwrap();
                                        *c2 += v.len() as i32;
                                    } else {
                                        combined_consensus_exons.insert(e.clone(), (0, v.len() as i32));
                                    }
                                }
                                for (e, counts) in combined_consensus_exons.iter() {
                                    if counts.0 * counts.1 == 0 && counts.0 + counts.1 >= min_sup_haplotype_exon as i32 {
                                        if e.state == 1 {
                                            // println!("exon1: {}:{}-{}, hap1:{:?}, hap2:{:?}", e.chr, e.start + 1, e.end + 1, counts.0, counts.1);
                                            haplotype_exons.push((e.clone(), counts.0, counts.1));
                                        }
                                        if e.state == 0 {
                                            let mut start_sum = 0;
                                            let mut start_mean = 0;
                                            let mut exon_cnt = 0;
                                            if counts.0 > 0 {
                                                for ex in hap1_consensus_exons.get(e).unwrap().iter() {
                                                    start_sum += ex.start;
                                                    exon_cnt += 1;
                                                }
                                                start_mean = start_sum / exon_cnt;
                                            } else {
                                                for ex in hap2_consensus_exons.get(e).unwrap().iter() {
                                                    start_sum += ex.start;
                                                    exon_cnt += 1;
                                                }
                                                start_mean = start_sum / exon_cnt;
                                            }
                                            // println!("exon0: {}:{}-{}, hap1:{:?}, hap2:{:?}", e.chr, start_mean + 1, e.end + 1, counts.0, counts.1);
                                            let mut ec = e.clone();
                                            ec.start = start_mean;
                                            haplotype_exons.push((ec, counts.0, counts.1));
                                        }
                                        if e.state == 2 {
                                            let mut end_sum = 0;
                                            let mut end_mean = 0;
                                            let mut exon_cnt = 0;
                                            if counts.0 > 0 {
                                                for ex in hap1_consensus_exons.get(e).unwrap().iter() {
                                                    end_sum += ex.end;
                                                    exon_cnt += 1;
                                                }
                                                end_mean = end_sum / exon_cnt;
                                            } else {
                                                for ex in hap2_consensus_exons.get(e).unwrap().iter() {
                                                    end_sum += ex.end;
                                                    exon_cnt += 1;
                                                }
                                                end_mean = end_sum / exon_cnt;
                                            }
                                            // println!("exon2: {}:{}-{}, hap1:{:?}, hap2:{:?}", e.chr, e.start + 1, end_mean + 1, counts.0, counts.1);
                                            let mut ec = e.clone();
                                            ec.end = end_mean;
                                            haplotype_exons.push((ec, counts.0, counts.1));
                                        }
                                        if e.state == 3 {
                                            let relaxed_start = [e.start - 20, e.start + 20];
                                            let relaxed_end = [e.end - 20, e.end + 20];
                                            let mut unique_flag = true;
                                            if counts.0 > 0 {
                                                for ex in hap2_consensus_exons.iter() {
                                                    if ex.0.state != 3 {
                                                        continue;
                                                    }
                                                    if ex.0.start >= relaxed_start[0] && ex.0.start <= relaxed_start[1] && ex.0.end >= relaxed_end[0] && ex.0.end <= relaxed_end[1] {
                                                        // highly overlapped, not unique
                                                        unique_flag = false;
                                                        break;
                                                    }
                                                }
                                            } else {
                                                for ex in hap1_consensus_exons.iter() {
                                                    if ex.0.state != 3 {
                                                        continue;
                                                    }
                                                    if ex.0.start >= relaxed_start[0] && ex.0.start <= relaxed_start[1] && ex.0.end >= relaxed_end[0] && ex.0.end <= relaxed_end[1] {
                                                        // highly overlapped, not unique
                                                        unique_flag = false;
                                                        break;
                                                    }
                                                }
                                            }
                                            if unique_flag {
                                                // println!("exon3: {}:{}-{}, hap1:{:?}, hap2:{:?}", e.chr, e.start + 1, e.end + 1, counts.0, counts.1);
                                                haplotype_exons.push((e.clone(), counts.0, counts.1));
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        if !no_bam_output {
                            let mut queue = read_haplotag_queue.lock().unwrap();
                            for a in read_assignments.iter() {
                                queue.push_back((a.0.clone(), a.1.clone()));
                            }
                        }
                    }
                    if haplotype_specific_exon {
                        let mut queue = haplotype_exon_queue.lock().unwrap();
                        for e in haplotype_exons.iter() {
                            queue.push_back(e.clone());
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
    vf.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n".as_bytes()).unwrap();
    vf.write("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n".as_bytes()).unwrap();
    vf.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n".as_bytes()).unwrap();
    vf.write("##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n".as_bytes()).unwrap();
    vf.write("##FORMAT=<ID=PQ,Number=1,Type=Float,Description=\"Phasing Quality\">\n".as_bytes()).unwrap();
    vf.write("##FORMAT=<ID=AE,Number=A,Type=Integer,Description=\"Haplotype expression of two alleles\">\n".as_bytes()).unwrap();
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

    if output_read_assignment {
        let mut assignment_writer = File::create(phased_bam_file.replace(".phased.bam", ".assignment.tsv")).unwrap();
        for rd in read_haplotag_queue.lock().unwrap().iter() {
            // println!("{}:{}", rd.0, rd.1);
            assignment_writer.write(format!("{}\t{}\n", rd.0, rd.1).as_bytes()).unwrap();
        }
        drop(assignment_writer);
    }

    if haplotype_specific_exon {
        let mut exon_hashmap: HashMap<String, Vec<(Exon, i32, i32)>> = HashMap::new();
        let mut exon_writer = File::create(phased_bam_file.replace(".phased.bam", ".haplotype_exon.tsv")).unwrap();
        exon_writer.write("#Chromosome\tExon start\tExon end\tExon state\tHap1 expression\tHap2 expression\n".as_bytes()).unwrap();
        for rd in haplotype_exon_queue.lock().unwrap().iter() {
            let e = &rd.0;
            if exon_hashmap.contains_key(&e.chr) {
                let exon_vec = exon_hashmap.get_mut(&e.chr).unwrap();
                exon_vec.push(rd.clone());
            } else {
                exon_hashmap.insert(e.chr.clone(), vec![rd.clone()]);
            }
        }
        for chr in contig_order.iter() {
            if !exon_hashmap.contains_key(chr) {
                continue;
            }
            let mut exons_sorted = exon_hashmap.get(chr).unwrap().clone();
            exons_sorted.sort_by(|a, b| a.0.start.cmp(&b.0.start));
            for rd in exons_sorted.iter() {
                let e = &rd.0;
                exon_writer.write(format!("{}\t{}\t{}\t{}\t{}\t{}\n", e.chr, e.start + 1, e.end, e.state, rd.1, rd.2).as_bytes()).unwrap();   // 1-based, start inclusive, end inclusive
            }
        }
        drop(exon_writer);
    }

    if !no_bam_output {
        if !haplotype_bam_output {
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
        } else {
            let mut hap1_read_assignments: HashSet<String> = HashSet::new();
            let mut hap2_read_assignments: HashSet<String> = HashSet::new();
            for rname in read_haplotag1_queue.lock().unwrap().iter() {
                hap1_read_assignments.insert(rname.clone());
            }
            for rname in read_haplotag2_queue.lock().unwrap().iter() {
                hap2_read_assignments.insert(rname.clone());
            }
            let mut bam_reader = bam::IndexedReader::from_path(&bam_file).unwrap();
            let header = bam::Header::from_template(&bam_reader.header());
            let mut hap1_bam_writer = bam::Writer::from_path(phased_bam_file.replace("phased", "hap1"), &header, Format::Bam).unwrap();
            let mut hap2_bam_writer = bam::Writer::from_path(phased_bam_file.replace("phased", "hap2"), &header, Format::Bam).unwrap();
            hap1_bam_writer.set_threads(thread_size).unwrap();
            hap2_bam_writer.set_threads(thread_size).unwrap();
            for region in isolated_regions.iter() {
                bam_reader.fetch((region.chr.as_str(), region.start, region.end)).unwrap(); // set region
                for r in bam_reader.records() {
                    let mut record = r.unwrap();
                    if record.is_unmapped() || record.is_secondary() || record.is_supplementary() {
                        continue;
                    }
                    let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
                    if hap1_read_assignments.contains(&qname) {
                        let _ = hap1_bam_writer.write(&record).unwrap();
                    } else if hap2_read_assignments.contains(&qname) {
                        let _ = hap2_bam_writer.write(&record).unwrap();
                    }
                }
            }
            drop(hap1_bam_writer);
            drop(hap2_bam_writer);
        }
    }
}

