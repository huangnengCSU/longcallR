use std::cmp::Ordering;
use std::collections::{HashMap, HashSet, VecDeque};
use std::fs;
use std::fs::File;
use std::hash::Hash;
use std::io::Write;
use std::sync::Mutex;

use bio::bio_types::strand::ReqStrand::Forward;
use petgraph::algo::kosaraju_scc;
use petgraph::prelude::*;
use probability::distribution::Distribution;
use rand::Rng;
use rayon::prelude::*;
use rust_htslib::{bam, bam::ext::BamRecordExtensions, bam::Format, bam::Read, bam::record::Aux, bam::record::Record};
use rust_lapper::{Interval, Lapper};

use crate::Platform;
use crate::util::{load_reference, parse_fai, Profile, Region, VCFRecord};

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
    // 0: homo ref, 1: heterozygous SNP, 2: homozygous SNP, 3: triallelic SNP
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
    pub ase: bool,
    // allele specific expressed
    pub single: bool,
    // current snp has surrounding haplotype links or not, only works for heterozygous snps
    pub somatic: bool,
    // somatic mutation
    pub phase_set: u32,
    // phase set id is the position of the first snp in the phase set
    pub haplotype_expression: [u32; 4],
    // hap1_ref, hap1_alt, hap2_ref, hap2_alt
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
    // base allele on alphabet  {-1, 1, 0}, 1: base==ref, -1: base==alt, 0: not covered (bases except ref allele and alt allele, deletions or N)
    pub prob: f64,
    // error rate of observe current base
    pub ase_snp: bool,
    // this is a potential allele specific expressed SNP
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
    pub num_hete_links: u32,
    // number of linked heterozygous snps in the fragment
    pub num_ase_links: u32,
    // number of linked allele specific expressed snps in the fragment
}

#[derive(Debug, Clone, Default)]
pub struct SNPFrag {
    pub region: Region,
    pub candidate_snps: Vec<CandidateSNP>,
    // candidate SNPs
    pub hete_snps: Vec<usize>,
    // index of candidate heterozygous SNPs, used for phasing of non-ase snps, non-filter snps, non-rna-edit snps
    pub homo_snps: Vec<usize>,
    // index of candidate homozygous SNPs
    pub ase_snps: Vec<usize>,
    // index of potential allele specific expressed SNPs
    pub hete_homo_snps: Vec<usize>,
    // index of candidate heterozygous SNPs and homozygous SNPs, used for dense snps filter
    pub ase_hete_snps: Vec<usize>,
    // index of potential allele specific expressed SNPs and heterozygous SNPs, used for construct fragment of hete snps and ase snps
    pub somatic_snps: Vec<usize>,
    // index of candidate somatic mutation
    pub fragments: Vec<Fragment>,
    // multiple fragments
    pub phased: bool,
    // haplotype is phased or not
    pub edges: HashMap<[usize; 2], Edge>,
    // edges of the graph, key is [snp_idx of start_node, snp_idx of end_node]
    pub min_linkers: u32,
    // the number of links for snps can be phased
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
    pub fn get_candidate_snps(
        &mut self,
        profile: &Profile,
        exon_region_vec: Vec<Interval<usize, u8>>,
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
        avg_dense_dist: f32,
        ase_allele_frac_cutoff: f32,
        ase_allele_cnt_cutoff: u32,
    ) {
        // get candidate SNPs, filtering with min_coverage, deletion_freq, min_allele_freq_include_intron, cover_strand_bias_threshold
        let pileup = &profile.freq_vec;
        let mut use_annotation: bool = false;
        if exon_region_vec.len() > 0 { use_annotation = true; }
        let mut exon_intervaltree = Lapper::new(exon_region_vec);
        let mut position = profile.region.start - 1; // 0-based
        for bfidx in 0..pileup.len() {
            let bf = &pileup[bfidx];
            if bf.i {
                continue;
            }

            if use_annotation && exon_intervaltree.find((position + 1) as usize, (position + 2) as usize).count() == 0 {
                // filter, not covered by exon
                position += 1;
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
                let total_cover_cnt = bf.forward_cnt + bf.backward_cnt; // does not include intron reads
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
                            'a' => {
                                [fcnt, bcnt] = bf.base_strands.a;
                            }
                            'A' => {
                                [fcnt, bcnt] = bf.base_strands.a;
                            }
                            'c' => {
                                [fcnt, bcnt] = bf.base_strands.c;
                            }
                            'C' => {
                                [fcnt, bcnt] = bf.base_strands.c;
                            }
                            'g' => {
                                [fcnt, bcnt] = bf.base_strands.g;
                            }
                            'G' => {
                                [fcnt, bcnt] = bf.base_strands.g;
                            }
                            't' => {
                                [fcnt, bcnt] = bf.base_strands.t;
                            }
                            'T' => {
                                [fcnt, bcnt] = bf.base_strands.t;
                            }
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

            let mut N_cnts: Vec<u32> = Vec::new(); // number of N bases (intron)
            let mut INS_cnts: Vec<u32> = Vec::new(); // number of insertions
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
            let theta = 0.001; // mutation rate
            let background_prob = [theta / 2.0, theta, 1.0 - 1.5 * theta]; // background of probability of observe homo variant, hete variant and homo reference
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
                println!(
                    "{}:{} unknown ref base {}",
                    profile.region.chr, position, bf.ref_base
                );
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
            loglikelihood[1] -= (num_reads as f64) * 2.0_f64.log10(); // example: logL(0) = -1, logL(1) = -6, logL(2) = -26

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
            variant_prob = [
                variant_prob[0] / sum_variant_prob,
                variant_prob[1] / sum_variant_prob,
                variant_prob[2] / sum_variant_prob,
            ];
            // println!("4:{}:{},{:?}", position, correction_factor, variant_prob);
            // QUAL phred-scaled quality score for the assertion made in ALT. i.e. give -10log_10 prob(call in ALT is wrong).
            // If ALT is `.` (no variant) then this is -10log_10 p(variant), and if ALT is not `.` this is -10log_10 p(no variant).
            let variant_quality = -10.0 * ((10e-301_f64.max(variant_prob[2])).log10()); // if variant_quality is greater than 3000, we set it to 3000

            // calculate GQ: The value of GQ is simply the difference between the second lowest PL and the lowest PL (which is always 0, normalized PL)
            let mut log10_likelihood = loglikelihood.clone();
            let max_log10_likelihood = log10_likelihood[0].max(log10_likelihood[1]).max(log10_likelihood[2]);
            log10_likelihood[0] = 10.0_f64.powf(log10_likelihood[0] - max_log10_likelihood);
            log10_likelihood[1] = 10.0_f64.powf(log10_likelihood[1] - max_log10_likelihood);
            log10_likelihood[2] = 10.0_f64.powf(log10_likelihood[2] - max_log10_likelihood);
            let sum_log10_likelihood = log10_likelihood[0] + log10_likelihood[1] + log10_likelihood[2];
            let mut genotype_prob = [
                log10_likelihood[0] / sum_log10_likelihood,
                log10_likelihood[1] / sum_log10_likelihood,
                log10_likelihood[2] / sum_log10_likelihood,
            ];
            // println!("{}:{},{},{},{}", profile.region.chr, position, genotype_prob[0], genotype_prob[1], genotype_prob[2]);
            let mut phred_genotype_prob = [0.0, 0.0, 0.0];
            phred_genotype_prob[0] = -10.0 * genotype_prob[0].log10(); // phred scale likelihood of genotype: 1/1
            phred_genotype_prob[1] = -10.0 * genotype_prob[1].log10(); // phred scale likelihood of genotype: 0/1
            phred_genotype_prob[2] = -10.0 * genotype_prob[2].log10(); // phred scale likelihood of genotype: 0/0
            phred_genotype_prob.sort_by(cmp_f64);
            let genotype_quality = phred_genotype_prob[1] - phred_genotype_prob[0];

            if genotype_prob[0] > genotype_prob[1] && genotype_prob[0] > genotype_prob[2] {
                if variant_quality < min_qual_for_candidate as f64 {
                    // TODO: keep the site to find somatic mutation
                    position += 1;
                    continue;
                }
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
                candidate_snp.variant_type = 2;

                if allele1 != bf.ref_base && allele2 != bf.ref_base {
                    if allele1_freq < min_homozygous_freq && allele2_freq > 0.0 {
                        candidate_snp.variant_type = 3; // triallelic SNP, triallelic SNP is also considered as homozygous SNP, e.g. ref: A, alt: C, G
                    }
                }

                candidate_snp.variant_quality = variant_quality;
                candidate_snp.genotype_probability = genotype_prob.clone();
                candidate_snp.genotype_quality = genotype_quality;
                if candidate_snp.variant_type == 2 {
                    if bf.ref_base == 'A' && allele1 == 'G' {
                        candidate_snp.rna_editing = true;
                    }
                    if bf.ref_base == 'T' && allele1 == 'C' {
                        candidate_snp.rna_editing = true;
                    }
                } else if candidate_snp.variant_type == 3 {
                    if bf.ref_base == 'A' && (allele1 == 'G' || allele2 == 'G') {
                        candidate_snp.rna_editing = true;
                    }
                    if bf.ref_base == 'T' && (allele1 == 'C' || allele2 == 'C') {
                        candidate_snp.rna_editing = true;
                    }
                }
                // println!("homo genotype quality: {:?},{:?}", likelihood, candidate_snp.genotype_quality);
                self.candidate_snps.push(candidate_snp);
                self.homo_snps.push(self.candidate_snps.len() - 1);
                self.hete_homo_snps.push(self.candidate_snps.len() - 1);
            } else if genotype_prob[1] > genotype_prob[0] && genotype_prob[1] > genotype_prob[2] {
                // candidate heterozygous SNP
                let allele1_freq = (allele1_cnt as f32) / (depth as f32);
                let allele2_freq = (allele2_cnt as f32) / (depth as f32);
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
                    }
                }
                if bf.ref_base == 'T' {
                    if (allele1 == 'T' && allele2 == 'C') || (allele1 == 'C' && allele2 == 'T') {
                        candidate_snp.rna_editing = true;
                    }
                }

                if allele1 != bf.ref_base && allele2 != bf.ref_base {
                    candidate_snp.variant_type = 3; // triallelic SNP
                    self.candidate_snps.push(candidate_snp);
                    self.homo_snps.push(self.candidate_snps.len() - 1);
                    self.hete_homo_snps.push(self.candidate_snps.len() - 1);
                    position += 1;
                    continue;
                }

                if allele1 != bf.ref_base && allele1_freq < min_allele_freq {
                    // alternative allele frequency is smaller than min allele freq
                    if allele1_freq >= ase_allele_frac_cutoff && allele1_cnt >= ase_allele_cnt_cutoff {
                        candidate_snp.ase = true;
                        candidate_snp.variant_type = 0;
                        self.candidate_snps.push(candidate_snp);
                        self.ase_snps.push(self.candidate_snps.len() - 1);
                        self.ase_hete_snps.push(self.candidate_snps.len() - 1);
                    }
                    // TODO: keep the site to find somatic mutation
                    position += 1;
                    continue;
                } else if allele2 != bf.ref_base && allele2_freq < min_allele_freq {
                    // alternative allele frequency is smaller than min allele freq
                    if allele2_freq >= ase_allele_frac_cutoff && allele2_cnt >= ase_allele_cnt_cutoff {
                        candidate_snp.ase = true;
                        candidate_snp.variant_type = 0;
                        self.candidate_snps.push(candidate_snp);
                        self.ase_snps.push(self.candidate_snps.len() - 1);
                        self.ase_hete_snps.push(self.candidate_snps.len() - 1);
                    }
                    // TODO: keep the site to find somatic mutation
                    position += 1;
                    continue;
                }
                if variant_quality < min_qual_for_candidate as f64 {
                    // heterozygous SNP may have low variant quality
                    if allele1 != bf.ref_base && allele1_freq >= ase_allele_frac_cutoff && allele1_cnt >= ase_allele_cnt_cutoff {
                        candidate_snp.ase = true;
                        candidate_snp.variant_type = 0;
                        self.candidate_snps.push(candidate_snp);
                        self.ase_snps.push(self.candidate_snps.len() - 1);
                        self.ase_hete_snps.push(self.candidate_snps.len() - 1);
                    } else if allele2 != bf.ref_base && allele2_freq >= ase_allele_frac_cutoff && allele2_cnt >= ase_allele_cnt_cutoff {
                        candidate_snp.ase = true;
                        candidate_snp.variant_type = 0;
                        self.candidate_snps.push(candidate_snp);
                        self.ase_snps.push(self.candidate_snps.len() - 1);
                        self.ase_hete_snps.push(self.candidate_snps.len() - 1);
                    }
                    // TODO: keep the site to find somatic mutation
                    position += 1;
                    continue;
                }
                // println!("hete genotype quality: {:?},{:?}", likelihood, candidate_snp.genotype_quality);
                self.candidate_snps.push(candidate_snp);
                self.hete_snps.push(self.candidate_snps.len() - 1);
                self.hete_homo_snps.push(self.candidate_snps.len() - 1);
                self.ase_hete_snps.push(self.candidate_snps.len() - 1);
            } else if genotype_prob[2] > genotype_prob[0] && genotype_prob[2] > genotype_prob[1] {
                let allele1_freq = (allele1_cnt as f32) / (depth as f32);
                let allele2_freq = (allele2_cnt as f32) / (depth as f32);
                let mut candidate_snp = CandidateSNP::default();
                candidate_snp.chromosome = profile.region.chr.clone().into_bytes();
                candidate_snp.pos = position as i64;
                candidate_snp.alleles = [allele1, allele2];
                candidate_snp.allele_freqs = [allele1_freq, allele2_freq];
                candidate_snp.reference = bf.ref_base;
                candidate_snp.depth = depth;
                candidate_snp.variant_type = 0;
                candidate_snp.variant_quality = variant_quality;
                candidate_snp.genotype_probability = genotype_prob.clone();
                candidate_snp.genotype_quality = genotype_quality;
                if allele1 != bf.ref_base && allele2 != bf.ref_base {
                    position += 1;
                    continue;
                }
                if allele1 != bf.ref_base && allele1_freq >= ase_allele_frac_cutoff && allele1_cnt >= ase_allele_cnt_cutoff {
                    candidate_snp.ase = true;
                    self.candidate_snps.push(candidate_snp);
                    self.ase_snps.push(self.candidate_snps.len() - 1);
                    self.ase_hete_snps.push(self.candidate_snps.len() - 1);
                } else if allele2 != bf.ref_base && allele2_freq >= ase_allele_frac_cutoff && allele2_cnt >= ase_allele_cnt_cutoff {
                    candidate_snp.ase = true;
                    self.candidate_snps.push(candidate_snp);
                    self.ase_snps.push(self.candidate_snps.len() - 1);
                    self.ase_hete_snps.push(self.candidate_snps.len() - 1);
                }
                // TODO: keep the site to find somatic mutation
            }
            position += 1;
        }

        // filter dense SNPs
        // for i in 0..self.candidate_snps.len() {
        //     for j in i..self.candidate_snps.len() {
        //         if self.candidate_snps[j].pos - self.candidate_snps[i].pos > dense_win_size as i64 {
        //             if (j - 1 - i + 1) as u32 >= min_dense_cnt && ((self.candidate_snps[j - 1].pos - self.candidate_snps[i].pos + 1) as f32) / ((j - 1 - i + 1) as f32) <= avg_dense_dist {
        //                 for tk in i..j {
        //                     // println!("dense SNPs: {}", self.candidate_snps[tk].pos);
        //                     // even rna editing may be filtered by dense SNPs
        //                     self.candidate_snps[tk].rna_editing = false;
        //                     self.candidate_snps[tk].filter = true;
        //                 }
        //             }
        //             break;
        //         }
        //     }
        // }

        // filter dense homozygous and heterozygous SNPs
        for i in 0..self.hete_homo_snps.len() {
            for j in i..self.hete_homo_snps.len() {
                if j == self.hete_homo_snps.len() - 1 {
                    // distance from snp i to end snp is smaller than dense_win_size
                    if self.candidate_snps[self.hete_homo_snps[j]].pos - self.candidate_snps[self.hete_homo_snps[i]].pos <= dense_win_size as i64 && (j - i + 1) as u32 >= min_dense_cnt {
                        for tk in i..j {
                            self.candidate_snps[self.hete_homo_snps[tk]].filter = true;
                        }
                    }
                }
                if self.candidate_snps[self.hete_homo_snps[j]].pos - self.candidate_snps[self.hete_homo_snps[i]].pos > dense_win_size as i64 {
                    if (j - 1 - i + 1) as u32 >= min_dense_cnt && ((self.candidate_snps[self.hete_homo_snps[j - 1]].pos - self.candidate_snps[self.hete_homo_snps[i]].pos + 1) as f32) / ((j - 1 - i + 1) as f32) <= avg_dense_dist {
                        for tk in i..j {
                            self.candidate_snps[self.hete_homo_snps[tk]].filter = true;
                        }
                    }
                    break;
                }
            }
        }

        // filter dense allele-specific expressed SNPs
        for i in 0..self.ase_snps.len() {
            for j in i..self.ase_snps.len() {
                if j == self.ase_snps.len() - 1 {
                    // distance from snp i to end snp is smaller than dense_win_size
                    if self.candidate_snps[self.ase_snps[j]].pos - self.candidate_snps[self.ase_snps[i]].pos <= dense_win_size as i64 && (j - i + 1) as u32 >= min_dense_cnt {
                        for tk in i..j {
                            self.candidate_snps[self.ase_snps[tk]].filter = true;
                        }
                    }
                }
                if self.candidate_snps[self.ase_snps[j]].pos - self.candidate_snps[self.ase_snps[i]].pos > dense_win_size as i64 {
                    if (j - 1 - i + 1) as u32 >= min_dense_cnt {
                        for tk in i..j {
                            self.candidate_snps[self.ase_snps[tk]].filter = true;
                        }
                    }
                    break;
                }
            }
        }

        // update homo_snps, remove filtered homo SNPs
        let mut tmp_homo: Vec<usize> = Vec::new();
        for i in self.homo_snps.iter() {
            if !self.candidate_snps[*i].filter {
                tmp_homo.push(*i);
            }
        }
        self.homo_snps = tmp_homo;

        // update hete_snps, remove filtered hete SNPs
        let mut tmp_hetes: Vec<usize> = Vec::new();
        for i in self.hete_snps.iter() {
            if !self.candidate_snps[*i].filter {
                tmp_hetes.push(*i);
            }
        }
        self.hete_snps = tmp_hetes;

        // update hete_homo_snps, remove filtered hete SNPs and homo SNPs
        let mut tmp_hete_homo: Vec<usize> = Vec::new();
        for i in self.hete_homo_snps.iter() {
            if !self.candidate_snps[*i].filter {
                tmp_hete_homo.push(*i);
            }
        }
        self.hete_homo_snps = tmp_hete_homo;

        // update ase_hete_snps, remove filtered hete SNPs
        let mut tmp_ase_hete: Vec<usize> = Vec::new();
        for i in self.ase_hete_snps.iter() {
            if !self.candidate_snps[*i].filter {
                tmp_ase_hete.push(*i);
            }
        }
        self.ase_hete_snps = tmp_ase_hete;
    }

    pub unsafe fn init_haplotypes(&mut self) {
        // initialize haplotype of heterozygous snp
        let mut rng = rand::thread_rng();
        for i in self.hete_snps.iter() {
            let rg: f64 = rng.gen();
            if rg < 0.5 {
                self.candidate_snps[*i].haplotype = -1;
            } else {
                self.candidate_snps[*i].haplotype = 1;
            }
        }
    }

    pub unsafe fn init_haplotypes_ase(&mut self) {
        // initialize haplotype of heterozygous snp
        let mut rng = rand::thread_rng();
        for i in self.ase_snps.iter() {
            let rg: f64 = rng.gen();
            if rg < 0.5 {
                self.candidate_snps[*i].haplotype = -1;
            } else {
                self.candidate_snps[*i].haplotype = 1;
            }
        }
    }

    pub unsafe fn init_assignment(&mut self) {
        let mut rng = rand::thread_rng();
        for k in 0..self.fragments.len() {
            if self.fragments[k].num_hete_links < self.min_linkers {
                continue;
            }
            let rg: f64 = rng.gen();
            if rg < 0.5 {
                self.fragments[k].haplotag = -1;
            } else {
                self.fragments[k].haplotag = 1;
            }
        }
    }

    pub fn get_fragments_with_ase(&mut self, bam_path: &str, region: &Region) {
        let mut bam_reader: bam::IndexedReader = bam::IndexedReader::from_path(bam_path).unwrap();
        bam_reader.fetch((region.chr.as_str(), region.start, region.end)).unwrap();
        let mut record = Record::new();
        if self.ase_hete_snps.len() == 0 {
            return;
        }
        // assert!(self.min_linkers >= 0, "Error: min_linkers <= 0");
        while let Some(result) = bam_reader.read(&mut record) {
            if result.is_err() {
                panic!("BAM parsing failed...");
            }
            // TODO: filtering unmapped, secondary, supplementary reads?
            if record.is_unmapped() || record.is_secondary() || record.is_supplementary() {
                continue;
            }
            let pos = record.pos(); // 0-based
            if pos > self.candidate_snps[*self.ase_hete_snps.last().unwrap()].pos {
                continue;
            }
            let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
            let cigar = record.cigar();
            let seq = record.seq().as_bytes();
            let strand = if record.strand() == Forward { 0 } else { 1 };
            let mut pos_on_ref = pos; // 0-based
            let mut pos_on_query = cigar.leading_softclips(); // 0-based
            let mut idx = 0; // index in self.ase_hete_snps
            let mut snp_pos = -1; // pre-computed position of candidate SNPs
            let mut alleles; // pre-computed alleles of candidate SNPs
            if pos <= self.candidate_snps[*self.ase_hete_snps.first().unwrap()].pos {
                snp_pos = self.candidate_snps[self.ase_hete_snps[idx]].pos;
                alleles = self.candidate_snps[self.ase_hete_snps[idx]].alleles.clone();
            } else {
                // find the first SNP in the read
                while idx < self.ase_hete_snps.len() {
                    if self.candidate_snps[self.ase_hete_snps[idx]].pos >= pos {
                        break;
                    }
                    idx += 1;
                }
                assert!(
                    idx < self.ase_hete_snps.len(),
                    "Error: idx < self.candidate_snps.len()"
                );
                snp_pos = self.candidate_snps[self.ase_hete_snps[idx]].pos;
                alleles = self.candidate_snps[self.ase_hete_snps[idx]].alleles.clone();
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
                                frag_elem.snp_idx = self.ase_hete_snps[idx];
                                frag_elem.pos = pos_on_ref;
                                frag_elem.base = seq[pos_on_query as usize] as char;
                                frag_elem.baseq = record.qual()[pos_on_query as usize];
                                frag_elem.strand = strand;
                                frag_elem.prob = 10.0_f64.powf(-(frag_elem.baseq as f64) / 10.0);
                                // if frag_elem.base == alleles[0] {
                                //     frag_elem.p = 1;    // reference allele
                                // } else if frag_elem.base == alleles[1] {
                                //     frag_elem.p = -1;   // alternate allele
                                // } else {
                                //     frag_elem.p = 0;    // not covered
                                // }
                                if frag_elem.base == self.candidate_snps[frag_elem.snp_idx].reference {
                                    frag_elem.p = 1; // reference allele
                                } else if (frag_elem.base == alleles[0] || frag_elem.base == alleles[1]) && frag_elem.base != self.candidate_snps[frag_elem.snp_idx].reference {
                                    frag_elem.p = -1; // alternate allele
                                } else {
                                    frag_elem.p = 0; // not covered
                                }
                                if self.candidate_snps[frag_elem.snp_idx].ase == true {
                                    frag_elem.ase_snp = true;
                                }
                                // filtered SNP will not be used for haplotype phasing, ase snp will still be used for construct fragment.
                                if self.candidate_snps[frag_elem.snp_idx].filter == false && frag_elem.p != 0 {
                                    fragment.list.push(frag_elem);
                                }
                                idx += 1;
                                if idx < self.ase_hete_snps.len() {
                                    snp_pos = self.candidate_snps[self.ase_hete_snps[idx]].pos;
                                    alleles = self.candidate_snps[self.ase_hete_snps[idx]].alleles.clone();
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
                                frag_elem.snp_idx = self.ase_hete_snps[idx];
                                frag_elem.pos = pos_on_ref;
                                frag_elem.base = '-';
                                frag_elem.baseq = 0;
                                frag_elem.strand = strand;
                                frag_elem.p = 0;
                                if self.candidate_snps[frag_elem.snp_idx].ase == true {
                                    frag_elem.ase_snp = true;
                                }
                                idx += 1;
                                if idx < self.ase_hete_snps.len() {
                                    snp_pos = self.candidate_snps[self.ase_hete_snps[idx]].pos;
                                    alleles = self.candidate_snps[self.ase_hete_snps[idx]].alleles.clone();
                                }
                            }
                            pos_on_ref += 1;
                        }
                    }
                    b'N' => {
                        exon_end = pos_on_ref;
                        if fragment.exons.len() == 0 {
                            fragment.exons.push(Exon {
                                chr: region.chr.clone(),
                                start: exon_start,
                                end: exon_end,
                                state: 0,
                            }); // start exon
                        } else {
                            fragment.exons.push(Exon {
                                chr: region.chr.clone(),
                                start: exon_start,
                                end: exon_end,
                                state: 1,
                            }); // internal exon
                        }
                        exon_start = -1;
                        exon_end = -1;
                        for _ in 0..cg.len() {
                            if pos_on_ref == snp_pos {
                                let mut frag_elem = FragElem::default();
                                frag_elem.snp_idx = self.ase_hete_snps[idx];
                                frag_elem.pos = pos_on_ref;
                                frag_elem.base = b'-' as char;
                                frag_elem.baseq = 0;
                                frag_elem.strand = strand;
                                frag_elem.p = 0;
                                if self.candidate_snps[frag_elem.snp_idx].ase == true {
                                    frag_elem.ase_snp = true;
                                }
                                idx += 1;
                                if idx < self.ase_hete_snps.len() {
                                    snp_pos = self.candidate_snps[self.ase_hete_snps[idx]].pos;
                                    alleles = self.candidate_snps[self.ase_hete_snps[idx]].alleles.clone();
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
                    fragment.exons.push(Exon {
                        chr: region.chr.clone(),
                        start: exon_start,
                        end: exon_end,
                        state: 2,
                    }); // end exon
                } else {
                    fragment.exons.push(Exon {
                        chr: region.chr.clone(),
                        start: exon_start,
                        end: exon_end,
                        state: 3,
                    }); // single exon
                }
            }
            exon_start = -1;
            exon_end = -1;

            // hete snps >= 1 || ase snps >= 2
            let mut hete_links = 0;
            let mut ase_links = 0;
            for fe in fragment.list.iter() {
                if fe.p != 0 {
                    if fe.ase_snp == true {
                        ase_links += 1;
                    } else {
                        hete_links += 1;
                    }
                }
            }
            fragment.num_hete_links = hete_links;
            fragment.num_ase_links = ase_links;
            // For hifi data, min_linkers is 1, for nanopore data, min_linkers is 2 (preset). For phasing, at least min_linkers hete snps or at least 2 ase snps.
            assert!(self.min_linkers > 0, "Error: min_linkers <= 0");
            if hete_links >= self.min_linkers || ase_links >= 2 {
                for fe in fragment.list.iter() {
                    // record each snp cover by which fragments
                    self.candidate_snps[fe.snp_idx].snp_cover_fragments.push(fragment.fragment_idx);
                }
                self.fragments.push(fragment);
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
        // assert!(self.min_linkers > 0, "Error: min_linkers <= 0");
        while let Some(result) = bam_reader.read(&mut record) {
            if result.is_err() {
                panic!("BAM parsing failed...");
            }
            // TODO: filtering unmapped, secondary, supplementary reads?
            if record.is_unmapped() || record.is_secondary() || record.is_supplementary() {
                continue;
            }
            let pos = record.pos(); // 0-based
            if pos > self.candidate_snps[*self.hete_snps.last().unwrap()].pos {
                continue;
            }
            let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
            let cigar = record.cigar();
            let seq = record.seq().as_bytes();
            let strand = if record.strand() == Forward { 0 } else { 1 };
            let mut pos_on_ref = pos; // 0-based
            let mut pos_on_query = cigar.leading_softclips(); // 0-based
            let mut idx = 0; // index in self.hete_snps
            let mut snp_pos = -1; // pre-computed position of candidate SNPs
            let mut alleles; // pre-computed alleles of candidate SNPs
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
                assert!(
                    idx < self.hete_snps.len(),
                    "Error: idx < self.candidate_snps.len()"
                );
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
                                // if frag_elem.base == alleles[0] {
                                //     frag_elem.p = 1;    // reference allele
                                // } else if frag_elem.base == alleles[1] {
                                //     frag_elem.p = -1;   // alternate allele
                                // } else {
                                //     frag_elem.p = 0;    // not covered
                                // }
                                if frag_elem.base == self.candidate_snps[frag_elem.snp_idx].reference {
                                    frag_elem.p = 1; // reference allele
                                } else if (frag_elem.base == alleles[0] || frag_elem.base == alleles[1]) && frag_elem.base != self.candidate_snps[frag_elem.snp_idx].reference {
                                    frag_elem.p = -1; // alternate allele
                                } else {
                                    frag_elem.p = 0; // not covered
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
                            fragment.exons.push(Exon {
                                chr: region.chr.clone(),
                                start: exon_start,
                                end: exon_end,
                                state: 0,
                            }); // start exon
                        } else {
                            fragment.exons.push(Exon {
                                chr: region.chr.clone(),
                                start: exon_start,
                                end: exon_end,
                                state: 1,
                            }); // internal exon
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
                    fragment.exons.push(Exon {
                        chr: region.chr.clone(),
                        start: exon_start,
                        end: exon_end,
                        state: 2,
                    }); // end exon
                } else {
                    fragment.exons.push(Exon {
                        chr: region.chr.clone(),
                        start: exon_start,
                        end: exon_end,
                        state: 3,
                    }); // single exon
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
            if snp.filter == true {
                continue;
            }
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
        let mut log_q1: f64 = 0.0;
        let mut log_q2: f64 = 0.0;
        let mut log_q3: f64 = 0.0;

        for i in 0..delta.len() {
            if sigma_k * delta[i] == ps[i] {
                log_q1 += (1.0 - probs[i]).log10();
            } else {
                log_q1 += probs[i].log10();
            }
        }

        for i in 0..delta.len() {
            if delta[i] == ps[i] {
                log_q2 += (1.0 - probs[i]).log10();
                log_q3 += probs[i].log10();
            } else {
                log_q2 += probs[i].log10();
                log_q3 += (1.0 - probs[i]).log10();
            }
        }
        let max_logq = log_q1.max(log_q2.max(log_q3));
        let q1 = 10.0_f64.powf(log_q1 - max_logq);
        let q2 = 10.0_f64.powf(log_q2 - max_logq);
        let q3 = 10.0_f64.powf(log_q3 - max_logq);
        // println!("sigma delta q1:{:?}, q2+q3:{:?}", q1, q2 + q3);
        return q1 / (q2 + q3);
    }

    pub fn cal_sigma_delta_log(
        sigma_k: i32,
        delta: &Vec<i32>,
        ps: &Vec<i32>,
        probs: &Vec<f64>,
    ) -> f64 {
        // same as call_sigma_delta, but return log10 value to avoid underflow
        let mut log_q1: f64 = 0.0;
        let mut log_q2: f64 = 0.0;
        let mut log_q3: f64 = 0.0;

        for i in 0..delta.len() {
            if sigma_k * delta[i] == ps[i] {
                log_q1 += (1.0 - probs[i]).log10();
            } else {
                log_q1 += probs[i].log10();
            }
        }

        for i in 0..delta.len() {
            if delta[i] == ps[i] {
                log_q2 += (1.0 - probs[i]).log10();
                log_q3 += probs[i].log10();
            } else {
                log_q2 += probs[i].log10();
                log_q3 += (1.0 - probs[i]).log10();
            }
        }

        // The exact formula: logP = log(\frac{A}{A+B})=logA-log(A+B)=logA-log(10^{logA}+10^{logB})
        // let log_p = log_q1 - f64::log10(10.0_f64.powf(log_q2) + 10.0_f64.powf(log_q3));
        // to avoid underflow, use approximate 1.0-log(A)/(log(A)+log(B)) as A/(A+B).
        // 0.99/(0.99+0.01) = 0.99, log(0.99)/(log(0.99)+log(0.01)) = 0.00217765
        // 0.01/(0.99+0.01) = 0.01, log(0.01)/(log(0.99)+log(0.01)) = 0.99782235
        let log_p = 1.0 - log_q1 / (log_q2 + log_q3);
        return log_p;
    }

    pub fn cal_delta_sigma(delta_i: i32, sigma: &Vec<i32>, ps: &Vec<i32>, probs: &Vec<f64>) -> f64 {
        // calculate P(delta_i | sigma)
        // delta_i: the haplotype of SNP i, 1 or -1.
        // sigma: the assignments of the reads cover SNP i, each haplotype is 1 or -1.
        // ps: the allele of each base, 1,-1
        // probs: the probability of observing base at SNP i for each read, equals to 10^(-Q/10).

        let mut log_q1: f64 = 0.0;
        let mut log_q2: f64 = 0.0;
        let mut log_q3: f64 = 0.0;

        for k in 0..sigma.len() {
            if delta_i * sigma[k] == ps[k] {
                log_q1 += (1.0 - probs[k]).log10();
            } else {
                log_q1 += probs[k].log10();
            }
        }

        for k in 0..sigma.len() {
            if sigma[k] == ps[k] {
                log_q2 += (1.0 - probs[k]).log10();
                log_q3 += probs[k].log10();
            } else {
                log_q2 += probs[k].log10();
                log_q3 += (1.0 - probs[k]).log10();
            }
        }
        let max_logq = log_q1.max(log_q2.max(log_q3));
        let q1 = 10.0_f64.powf(log_q1 - max_logq);
        let q2 = 10.0_f64.powf(log_q2 - max_logq);
        let q3 = 10.0_f64.powf(log_q3 - max_logq);
        // println!("delta sigma q1:{:?}, q2+q3:{:?}", q1, q2 + q3);
        return q1 / (q2 + q3);
    }

    pub fn cal_delta_sigma_log(
        delta_i: i32,
        sigma: &Vec<i32>,
        ps: &Vec<i32>,
        probs: &Vec<f64>,
    ) -> f64 {
        // same as call_delta_sigma, but return log10 value to avoid underflow
        let mut log_q1: f64 = 0.0;
        let mut log_q2: f64 = 0.0;
        let mut log_q3: f64 = 0.0;

        for k in 0..sigma.len() {
            if delta_i * sigma[k] == ps[k] {
                log_q1 += (1.0 - probs[k]).log10();
            } else {
                log_q1 += probs[k].log10();
            }
        }

        for k in 0..sigma.len() {
            if sigma[k] == ps[k] {
                log_q2 += (1.0 - probs[k]).log10();
                log_q3 += probs[k].log10();
            } else {
                log_q2 += probs[k].log10();
                log_q3 += (1.0 - probs[k]).log10();
            }
        }

        // logP = log(\frac{A}{A+B})=logA-log(A+B)=logA-log(10^{logA}+10^{logB})
        // let log_p = log_q1 - f64::log10(10.0_f64.powf(log_q2) + 10.0_f64.powf(log_q3));
        // to avoid underflow, use approximate 1.0-log(A)/(log(A)+log(B)) as A/(A+B).
        // 0.99/(0.99+0.01) = 0.99, log(0.99)/(log(0.99)+log(0.01)) = 0.00217765
        // 0.01/(0.99+0.01) = 0.01, log(0.01)/(log(0.99)+log(0.01)) = 0.99782235
        let log_p = 1.0 - log_q1 / (log_q2 + log_q3);
        return log_p;
    }

    pub fn cal_inconsistent_percentage(
        delta_i: i32,
        sigma: &Vec<i32>,
        ps: &Vec<i32>,
        probs: &Vec<f64>,
    ) -> f64 {
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

    pub fn cal_overall_probability(snpfrag: &SNPFrag) -> f64 {
        // calculate the log10 probability of the current configuration of sigma and delta
        let mut logp = 0.0;
        for k in 0..snpfrag.fragments.len() {
            if snpfrag.fragments[k].haplotag == 0 {
                continue;
            }
            for fe in snpfrag.fragments[k].list.iter() {
                if fe.ase_snp == true {
                    continue;
                }
                assert_ne!(fe.p, 0, "Error: phasing with unexpected hete SNP.");
                if snpfrag.fragments[k].haplotag * snpfrag.candidate_snps[fe.snp_idx].haplotype == fe.p {
                    logp += (1.0 - fe.prob).log10();
                } else {
                    logp += fe.prob.log10();
                }
            }
        }
        return logp;
    }

    pub fn cal_overall_probability_ase(snpfrag: &SNPFrag) -> f64 {
        // calculate the log10 probability of the current configuration of sigma and delta
        let mut logp = 0.0;
        for k in 0..snpfrag.fragments.len() {
            for fe in snpfrag.fragments[k].list.iter() {
                assert_ne!(fe.p, 0, "Error: phasing with unexpected hete SNP.");
                if snpfrag.fragments[k].haplotag * snpfrag.candidate_snps[fe.snp_idx].haplotype == fe.p {
                    logp += (1.0 - fe.prob).log10();
                } else {
                    logp += fe.prob.log10();
                }
            }
        }
        return logp;
    }

    pub fn check_new_haplotag(snpfrag: &SNPFrag, updated_haplotag: &HashMap<usize, i32>) -> i32 {
        // updated_haplotag: the index of the fragments will be updated
        let mut logp = 0.0;
        let mut pre_logp = 0.0;
        for (k, h) in updated_haplotag.iter() {
            let mut delta: Vec<i32> = Vec::new();
            let mut ps: Vec<i32> = Vec::new();
            let mut probs: Vec<f64> = Vec::new();
            if snpfrag.fragments[*k].haplotag == 0 {
                continue;
            }
            for fe in snpfrag.fragments[*k].list.iter() {
                if fe.ase_snp == true {
                    continue;
                }
                assert_ne!(fe.p, 0, "Error: phasing with unexpected hete SNP.");
                ps.push(fe.p);
                probs.push(fe.prob);
                delta.push(snpfrag.candidate_snps[fe.snp_idx].haplotype);
            }
            if delta.len() == 0 {
                continue;
            }
            logp += SNPFrag::cal_sigma_delta_log(*h, &delta, &ps, &probs);
            pre_logp += SNPFrag::cal_sigma_delta_log(snpfrag.fragments[*k].haplotag, &delta, &ps, &probs);
        }

        let p = logp;
        let pre_p = pre_logp;
        // println!("haplotag p:{}, pre_p:{}", p, pre_p);
        let mut rv = 0;
        if p > pre_p {
            rv = 1;
        } else if p == pre_p {
            rv = 0;
        } else {
            rv = -1;
        }
        assert!(
            rv >= 0,
            "Error: update haplotag should not decrease the probability. {} -> {}",
            pre_logp,
            logp
        );
        return rv;
    }

    pub fn check_new_haplotag_ase(
        snpfrag: &SNPFrag,
        updated_haplotag: &HashMap<usize, i32>,
    ) -> i32 {
        // updated_haplotag: the index of the fragments will be updated
        let mut logp = 0.0;
        let mut pre_logp = 0.0;
        for (k, h) in updated_haplotag.iter() {
            let mut delta: Vec<i32> = Vec::new();
            let mut ps: Vec<i32> = Vec::new();
            let mut probs: Vec<f64> = Vec::new();
            for fe in snpfrag.fragments[*k].list.iter() {
                ps.push(fe.p);
                probs.push(fe.prob);
                delta.push(snpfrag.candidate_snps[fe.snp_idx].haplotype);
            }
            if delta.len() == 0 {
                continue;
            }
            logp += SNPFrag::cal_sigma_delta_log(*h, &delta, &ps, &probs);
            pre_logp += SNPFrag::cal_sigma_delta_log(snpfrag.fragments[*k].haplotag, &delta, &ps, &probs);
        }

        let p = logp;
        let pre_p = pre_logp;
        // println!("haplotag p:{}, pre_p:{}", p, pre_p);
        let mut rv = 0;
        if p > pre_p {
            rv = 1;
        } else if p == pre_p {
            rv = 0;
        } else {
            rv = -1;
        }
        assert!(
            rv >= 0,
            "Error: update haplotag should not decrease the probability. {} -> {}",
            pre_logp,
            logp
        );
        return rv;
    }

    pub fn check_new_haplotype(snpfrag: &SNPFrag, updated_haplotype: &HashMap<usize, i32>) -> i32 {
        let mut logp = 0.0;
        let mut pre_logp = 0.0;
        for (i, h) in updated_haplotype.iter() {
            let mut sigma: Vec<i32> = Vec::new();
            let mut ps: Vec<i32> = Vec::new();
            let mut probs: Vec<f64> = Vec::new();
            for k in snpfrag.candidate_snps[*i].snp_cover_fragments.iter() {
                if snpfrag.fragments[*k].haplotag == 0 {
                    continue;
                }
                for fe in snpfrag.fragments[*k].list.iter() {
                    if fe.snp_idx != *i {
                        continue;
                    }
                    if fe.ase_snp == true {
                        continue;
                    }
                    assert_ne!(fe.p, 0, "Error: phasing with unexpected hete SNP.");
                    ps.push(fe.p);
                    probs.push(fe.prob);
                    sigma.push(snpfrag.fragments[*k].haplotag);
                }
            }
            if sigma.len() == 0 {
                // println!("SNP: {:?}", snpfrag.candidate_snps[*i]);
                // println!("SNP {} is not covered by any fragment.", snpfrag.candidate_snps[*i].pos);
                continue;
            }
            logp += SNPFrag::cal_delta_sigma_log(*h, &sigma, &ps, &probs);
            pre_logp += SNPFrag::cal_delta_sigma_log(
                snpfrag.candidate_snps[*i].haplotype,
                &sigma,
                &ps,
                &probs,
            );
        }
        let p = logp;
        let pre_p = pre_logp;
        // println!("haplotype p:{}, pre_p:{}", p, pre_p);
        let mut rv = 0;
        if p > pre_p {
            rv = 1;
        } else if p == pre_p {
            rv = 0;
        } else {
            rv = -1;
        }
        assert!(
            rv >= 0,
            "Error: update haplotype should not decrease the probability. {} -> {}",
            pre_logp,
            logp
        );
        return rv;
    }

    pub fn check_new_haplotype_ase(
        snpfrag: &SNPFrag,
        updated_haplotype: &HashMap<usize, i32>,
    ) -> i32 {
        let mut logp = 0.0;
        let mut pre_logp = 0.0;
        for (i, h) in updated_haplotype.iter() {
            let mut sigma: Vec<i32> = Vec::new();
            let mut ps: Vec<i32> = Vec::new();
            let mut probs: Vec<f64> = Vec::new();
            for k in snpfrag.candidate_snps[*i].snp_cover_fragments.iter() {
                for fe in snpfrag.fragments[*k].list.iter() {
                    if fe.snp_idx != *i {
                        continue;
                    }
                    assert_ne!(fe.p, 0, "Error: phasing with unexpected hete SNP.");
                    ps.push(fe.p);
                    probs.push(fe.prob);
                    sigma.push(snpfrag.fragments[*k].haplotag);
                }
            }
            if sigma.len() == 0 {
                // println!("SNP: {:?}", snpfrag.candidate_snps[*i]);
                // println!("SNP {} is not covered by any fragment.", snpfrag.candidate_snps[*i].pos);
                continue;
            }
            logp += SNPFrag::cal_delta_sigma_log(*h, &sigma, &ps, &probs);
            pre_logp += SNPFrag::cal_delta_sigma_log(
                snpfrag.candidate_snps[*i].haplotype,
                &sigma,
                &ps,
                &probs,
            );
        }
        let p = logp;
        let pre_p = pre_logp;
        // println!("haplotype p:{}, pre_p:{}", p, pre_p);
        let mut rv = 0;
        if p > pre_p {
            rv = 1;
        } else if p == pre_p {
            rv = 0;
        } else {
            rv = -1;
        }
        assert!(
            rv >= 0,
            "Error: update haplotype should not decrease the probability. {} -> {}",
            pre_logp,
            logp
        );
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

        while phasing_increase | haplotag_increase {
            // optimize sigma
            let mut tmp_haplotag: HashMap<usize, i32> = HashMap::new();
            let mut processed_snps = HashSet::new(); // some snps in self.hete_snps may be filtered by previous steps, record the snps that covered by the fragments
            for k in 0..self.fragments.len() {
                let sigma_k = self.fragments[k].haplotag;
                let mut delta: Vec<i32> = Vec::new();
                let mut ps: Vec<i32> = Vec::new();
                let mut probs: Vec<f64> = Vec::new();
                if sigma_k == 0 {
                    continue;
                }
                for fe in self.fragments[k].list.iter() {
                    if fe.ase_snp == true {
                        continue;
                    }
                    assert_ne!(fe.p, 0, "Error: phase for unexpected allele.");
                    ps.push(fe.p);
                    probs.push(fe.prob);
                    delta.push(self.candidate_snps[fe.snp_idx].haplotype);
                    processed_snps.insert(fe.snp_idx);
                }

                let q = SNPFrag::cal_sigma_delta_log(sigma_k, &delta, &ps, &probs);
                let qn = SNPFrag::cal_sigma_delta_log(sigma_k * (-1), &delta, &ps, &probs);
                // println!("optimize sigma {} q:{}, qn:{}, sigma: {}", k, q, qn, sigma_k);

                if q < qn {
                    tmp_haplotag.insert(k, sigma_k * (-1));
                } else {
                    tmp_haplotag.insert(k, sigma_k);
                }
            }

            // assert!(SNPFrag::cal_overall_probability(&self, &processed_snps, &self.haplotype) >= SNPFrag::cal_overall_probability(&self, &self.haplotag, &self.haplotype));
            let check_val = SNPFrag::check_new_haplotag(&self, &tmp_haplotag);
            assert!(check_val >= 0, "ckeck val bug: {:?}", self.candidate_snps);
            for (k, h) in tmp_haplotag.iter() {
                // when prob is equal, we still perform the flip to avoid bug of underflow
                self.fragments[*k].haplotag = *h;
            }
            if check_val == 0 {
                haplotag_increase = false;
            } else {
                haplotag_increase = true;
                phasing_increase = true;
            }
            self.check_local_optimal_configuration(false, true);

            // optimize delta
            let mut tmp_haplotype: HashMap<usize, i32> = HashMap::new();
            for i in self.hete_snps.iter() {
                let delta_i = self.candidate_snps[*i].haplotype;
                let mut sigma: Vec<i32> = Vec::new();
                let mut ps: Vec<i32> = Vec::new();
                let mut probs: Vec<f64> = Vec::new();
                for k in self.candidate_snps[*i].snp_cover_fragments.iter() {
                    if self.fragments[*k].haplotag == 0 {
                        continue;
                    }
                    // k is fragment index
                    for fe in self.fragments[*k].list.iter() {
                        if fe.snp_idx == *i {
                            assert_ne!(fe.ase_snp, true, "Error: phase for unexpected ase SNP.");
                            assert_ne!(fe.p, 0, "Error: phase for unexpected allele.");
                            ps.push(fe.p);
                            probs.push(fe.prob);
                            sigma.push(self.fragments[*k].haplotag);
                        }
                    }
                }

                let q = SNPFrag::cal_delta_sigma_log(delta_i, &sigma, &ps, &probs);
                let qn = SNPFrag::cal_delta_sigma_log(delta_i * (-1), &sigma, &ps, &probs);
                // println!("optimize delta {} q:{:?}, qn:{:?}, delta: {}", i, q, qn, delta_i);
                if q < qn {
                    tmp_haplotype.insert(*i, delta_i * (-1));
                } else {
                    tmp_haplotype.insert(*i, delta_i);
                }
            }
            let check_val = SNPFrag::check_new_haplotype(&self, &tmp_haplotype);
            assert!(check_val >= 0, "ckeck val bug: {:?}", self.candidate_snps);
            for (i, h) in tmp_haplotype.iter() {
                // when prob is equal, we still perform the flip to avoid bug of underflow
                self.candidate_snps[*i].haplotype = *h;
            }
            if check_val == 0 {
                phasing_increase = false;
            } else {
                phasing_increase = true;
                haplotag_increase = true;
            }
            self.check_local_optimal_configuration(true, false);
            num_iters += 1;
            if num_iters > 20 {
                break;
            }
        }
        // sigma reaches the optimal solution first and then delta reaches the optimal solution. After this, equal probability flip of delta may destroy the optimum of sigma again.
        // self.check_local_optimal_configuration(true, true);
        let prob = SNPFrag::cal_overall_probability(&self);
        return prob;
    }

    pub fn cross_optimize_ase_snps(&mut self) -> f64 {
        // Iteration:
        //     1. evaluate the assignment of each read based on the current SNP haplotype.
        //     2. evaluate the SNP haplotype based on the read assignment.
        // If P(sigma, delta) increase, repeat Iteration;
        // Else break;

        let mut phasing_increase: bool = true;
        let mut haplotag_increase: bool = true;
        let mut num_iters = 0;

        while phasing_increase | haplotag_increase {
            // optimize sigma
            let mut tmp_haplotag: HashMap<usize, i32> = HashMap::new();
            let mut processed_snps = HashSet::new(); // some snps in self.hete_snps may be filtered by previous steps, record the snps that covered by the fragments
            for k in 0..self.fragments.len() {
                let mut sigma_k = self.fragments[k].haplotag;
                if sigma_k == 0 {
                    // the fragment may be unphased in the first step phasing with heterozygous snps, random assign sigma
                    let mut rng = rand::thread_rng();
                    let rg: f64 = rng.gen();
                    if rg < 0.5 {
                        self.fragments[k].haplotag = -1;
                        sigma_k = -1;
                    } else {
                        self.fragments[k].haplotag = 1;
                        sigma_k = 1;
                    }
                }
                let mut delta: Vec<i32> = Vec::new();
                let mut ps: Vec<i32> = Vec::new();
                let mut probs: Vec<f64> = Vec::new();
                for fe in self.fragments[k].list.iter() {
                    ps.push(fe.p);
                    probs.push(fe.prob);
                    delta.push(self.candidate_snps[fe.snp_idx].haplotype);
                    processed_snps.insert(fe.snp_idx);
                }

                let q = SNPFrag::cal_sigma_delta_log(sigma_k, &delta, &ps, &probs);
                let qn = SNPFrag::cal_sigma_delta_log(sigma_k * (-1), &delta, &ps, &probs);
                // println!("optimize sigma {} q:{}, qn:{}, sigma: {}", k, q, qn, sigma_k);

                if q < qn {
                    tmp_haplotag.insert(k, sigma_k * (-1));
                } else {
                    tmp_haplotag.insert(k, sigma_k);
                }
            }

            // assert!(SNPFrag::cal_overall_probability(&self, &processed_snps, &self.haplotype) >= SNPFrag::cal_overall_probability(&self, &self.haplotag, &self.haplotype));
            let check_val = SNPFrag::check_new_haplotag_ase(&self, &tmp_haplotag);
            assert!(check_val >= 0, "ckeck val bug: {:?}", self.candidate_snps);
            for (k, h) in tmp_haplotag.iter() {
                // when prob is equal, we still perform the flip to avoid bug of underflow
                self.fragments[*k].haplotag = *h;
            }
            if check_val == 0 {
                haplotag_increase = false;
            } else {
                haplotag_increase = true;
                phasing_increase = true;
            }
            self.check_local_optimal_configuration_ase(false, true);

            // optimize delta
            let mut tmp_haplotype: HashMap<usize, i32> = HashMap::new();
            for i in self.ase_snps.iter() {
                // flip ase snps and keep heterozygous snps unchanged
                let delta_i = self.candidate_snps[*i].haplotype;
                let mut sigma: Vec<i32> = Vec::new();
                let mut ps: Vec<i32> = Vec::new();
                let mut probs: Vec<f64> = Vec::new();
                for k in self.candidate_snps[*i].snp_cover_fragments.iter() {
                    // k is fragment index
                    for fe in self.fragments[*k].list.iter() {
                        if fe.snp_idx == *i {
                            assert_ne!(fe.p, 0, "Error: phase for unexpected allele.");
                            ps.push(fe.p);
                            probs.push(fe.prob);
                            sigma.push(self.fragments[*k].haplotag);
                        }
                    }
                }

                let q = SNPFrag::cal_delta_sigma_log(delta_i, &sigma, &ps, &probs);
                let qn = SNPFrag::cal_delta_sigma_log(delta_i * (-1), &sigma, &ps, &probs);
                // println!("optimize delta {} q:{:?}, qn:{:?}, delta: {}", i, q, qn, delta_i);
                if q < qn {
                    tmp_haplotype.insert(*i, delta_i * (-1));
                } else {
                    tmp_haplotype.insert(*i, delta_i);
                }
            }
            let check_val = SNPFrag::check_new_haplotype_ase(&self, &tmp_haplotype);
            assert!(check_val >= 0, "ckeck val bug: {:?}", self.candidate_snps);
            for (i, h) in tmp_haplotype.iter() {
                // when prob is equal, we still perform the flip to avoid bug of underflow
                self.candidate_snps[*i].haplotype = *h;
            }
            if check_val == 0 {
                phasing_increase = false;
            } else {
                phasing_increase = true;
                haplotag_increase = true;
            }
            self.check_local_optimal_configuration_ase(true, false); // maybe failed
            num_iters += 1;
            if num_iters > 20 {
                break;
            }
        }
        // sigma reaches the optimal solution first and then delta reaches the optimal solution. After this, equal probability flip of delta may destroy the optimum of sigma again.
        // self.check_local_optimal_configuration(true, true);
        let prob = SNPFrag::cal_overall_probability_ase(&self);
        return prob;
    }

    fn save_best_configuration(
        &self,
        best_haplotype: &mut HashMap<usize, i32>,
        best_haplotag: &mut HashMap<usize, i32>,
    ) {
        best_haplotype.clear();
        best_haplotag.clear();
        for i in self.hete_snps.iter() {
            best_haplotype.insert(*i, self.candidate_snps[*i].haplotype);
        }
        for k in 0..self.fragments.len() {
            best_haplotag.insert(k, self.fragments[k].haplotag);
        }
    }

    fn load_best_configuration(
        &mut self,
        best_haplotype: &HashMap<usize, i32>,
        best_haplotag: &HashMap<usize, i32>,
    ) {
        for i in self.hete_snps.iter() {
            self.candidate_snps[*i].haplotype = best_haplotype.get(i).unwrap().clone();
        }
        for k in 0..self.fragments.len() {
            self.fragments[k].haplotag = best_haplotag.get(&k).unwrap().clone();
        }
    }

    pub fn phase(&mut self, max_enum_snps: usize, random_flip_fraction: f32, max_iters: i32) {
        let mut largest_prob = f64::NEG_INFINITY;
        let mut best_haplotype: HashMap<usize, i32> = HashMap::new();
        let mut best_haplotag: HashMap<usize, i32> = HashMap::new();

        if self.hete_snps.len() <= max_enum_snps {
            // enumerate the haplotype, then optimize the assignment
            let mut haplotype_enum: Vec<Vec<i32>> = Vec::new();
            let init_hap: Vec<i32> = vec![1; self.hete_snps.len()];
            haplotype_enum.push(init_hap.clone());
            for ti in 0..self.hete_snps.len() {
                for tj in 0..haplotype_enum.len() {
                    let mut tmp_hap = haplotype_enum[tj].clone();
                    tmp_hap[ti] = tmp_hap[ti] * (-1);
                    haplotype_enum.push(tmp_hap);
                }
            }
            assert!(
                haplotype_enum.len() == 2_usize.pow(self.hete_snps.len() as u32),
                "Error: Not all combinations included"
            );
            for hap in haplotype_enum.iter() {
                for i in 0..self.hete_snps.len() {
                    self.candidate_snps[self.hete_snps[i]].haplotype = hap[i];
                }
                unsafe {
                    self.init_assignment();
                }
                let prob = self.cross_optimize();
                if prob > largest_prob {
                    largest_prob = prob;
                    self.save_best_configuration(&mut best_haplotype, &mut best_haplotag);
                }
            }
            self.load_best_configuration(&best_haplotype, &best_haplotag);
        } else {
            // optimize haplotype and read assignment alternatively
            let mut max_iter: i32 = max_iters;
            while max_iter >= 0 {
                // random initialization of haplotype and haplotag at each iteration
                unsafe {
                    self.init_haplotypes();
                }
                unsafe {
                    self.init_assignment();
                }
                let prob = self.cross_optimize();
                if prob > largest_prob {
                    largest_prob = prob;
                    self.save_best_configuration(&mut best_haplotype, &mut best_haplotag);
                }
                self.load_best_configuration(&best_haplotype, &best_haplotag);

                // when initial setting has reached to local optimal, flip all the snps after a specific position to jump out local optimization
                let mut unflipped_haplotype: Vec<i32> = Vec::new();
                for i in self.hete_snps.iter() {
                    unflipped_haplotype.push(self.candidate_snps[*i].haplotype);
                }
                for ti in 0..unflipped_haplotype.len() {
                    let mut tmp_hap: Vec<i32> = Vec::new();
                    for tj in 0..unflipped_haplotype.len() {
                        if tj < ti {
                            tmp_hap.push(unflipped_haplotype[tj]);
                        } else {
                            tmp_hap.push(unflipped_haplotype[tj] * (-1));
                        }
                    }
                    // block flip
                    {
                        assert_eq!(tmp_hap.len(), self.hete_snps.len());
                        for i in 0..self.hete_snps.len() {
                            self.candidate_snps[self.hete_snps[i]].haplotype = tmp_hap[i];
                        }
                    }
                    let prob = self.cross_optimize();
                    if prob > largest_prob {
                        largest_prob = prob;
                        self.save_best_configuration(&mut best_haplotype, &mut best_haplotag);
                    }
                    self.load_best_configuration(&best_haplotype, &best_haplotag);

                    // when current block flip has reached to local optimal, flip a fraction of snps and reads to jump out local optimization
                    {
                        let mut rng = rand::thread_rng();
                        for ti in 0..self.hete_snps.len() {
                            let rg: f64 = rng.gen();
                            if rg < random_flip_fraction as f64 {
                                self.candidate_snps[self.hete_snps[ti]].haplotype = self.candidate_snps[self.hete_snps[ti]].haplotype * (-1);
                            }
                        }
                        for tk in 0..self.fragments.len() {
                            if self.fragments[tk].haplotag == 0 {
                                continue;
                            }
                            let rg: f64 = rng.gen();
                            if rg < random_flip_fraction as f64 {
                                self.fragments[tk].haplotag = self.fragments[tk].haplotag * (-1);
                            }
                        }
                    }
                    let prob = self.cross_optimize();
                    if prob > largest_prob {
                        largest_prob = prob;
                        self.save_best_configuration(&mut best_haplotype, &mut best_haplotag);
                    }
                    self.load_best_configuration(&best_haplotype, &best_haplotag);
                }
                self.load_best_configuration(&best_haplotype, &best_haplotag);
                max_iter -= 1;
            }
            self.load_best_configuration(&best_haplotype, &best_haplotag);
        }
    }

    pub fn phase_ase_hete_snps(
        &mut self,
        max_enum_snps: usize,
        random_flip_fraction: f32,
        max_iters: i32,
    ) {
        let mut largest_prob = f64::NEG_INFINITY;
        let mut best_haplotype: HashMap<usize, i32> = HashMap::new();
        let mut best_haplotag: HashMap<usize, i32> = HashMap::new();
        if self.ase_snps.len() <= max_enum_snps {
            // enumerate the haplotype, then optimize the assignment
            let mut haplotype_enum: Vec<Vec<i32>> = Vec::new();
            let init_hap: Vec<i32> = vec![1; self.ase_snps.len()];
            haplotype_enum.push(init_hap.clone());
            for ti in 0..self.ase_snps.len() {
                for tj in 0..haplotype_enum.len() {
                    let mut tmp_hap = haplotype_enum[tj].clone();
                    tmp_hap[ti] = tmp_hap[ti] * (-1);
                    haplotype_enum.push(tmp_hap);
                }
            }
            assert!(
                haplotype_enum.len() == 2_usize.pow(self.ase_snps.len() as u32),
                "Error: Not all combinations included"
            );
            for hap in haplotype_enum.iter() {
                for i in 0..self.ase_snps.len() {
                    assert!(
                        self.candidate_snps[self.ase_snps[i]].ase == true,
                        "Error: only ase snps can be flipped."
                    );
                    self.candidate_snps[self.ase_snps[i]].haplotype = hap[i];
                }
                let prob = self.cross_optimize_ase_snps();
                if prob > largest_prob {
                    largest_prob = prob;
                    self.save_best_configuration(&mut best_haplotype, &mut best_haplotag);
                }
            }
            self.load_best_configuration(&best_haplotype, &best_haplotag);
        } else {
            // optimize haplotype and read assignment alternatively
            let mut max_iter: i32 = max_iters;
            while max_iter >= 0 {
                // random initialization of haplotype for ase snps
                unsafe {
                    self.init_haplotypes_ase();
                }
                let prob = self.cross_optimize_ase_snps();
                if prob > largest_prob {
                    largest_prob = prob;
                    self.save_best_configuration(&mut best_haplotype, &mut best_haplotag);
                }
                self.load_best_configuration(&best_haplotype, &best_haplotag);

                // when initial setting has reached to local optimal, flip all the snps after a specific position to jump out local optimization
                let mut unflipped_haplotype: Vec<i32> = Vec::new();
                for i in self.ase_snps.iter() {
                    unflipped_haplotype.push(self.candidate_snps[*i].haplotype);
                }
                for ti in 0..unflipped_haplotype.len() {
                    let mut tmp_hap: Vec<i32> = Vec::new();
                    for tj in 0..unflipped_haplotype.len() {
                        if tj < ti {
                            tmp_hap.push(unflipped_haplotype[tj]);
                        } else {
                            tmp_hap.push(unflipped_haplotype[tj] * (-1));
                        }
                    }
                    // block flip
                    {
                        assert_eq!(tmp_hap.len(), self.ase_snps.len());
                        for i in 0..self.ase_snps.len() {
                            self.candidate_snps[self.ase_snps[i]].haplotype = tmp_hap[i];
                        }
                    }
                    let prob = self.cross_optimize_ase_snps();
                    if prob > largest_prob {
                        largest_prob = prob;
                        self.save_best_configuration(&mut best_haplotype, &mut best_haplotag);
                    }
                    self.load_best_configuration(&best_haplotype, &best_haplotag);

                    // when current block flip has reached to local optimal, flip a fraction of snps and reads to jump out local optimization
                    {
                        let mut rng = rand::thread_rng();
                        for ti in 0..self.ase_snps.len() {
                            let rg: f64 = rng.gen();
                            if rg < random_flip_fraction as f64 {
                                self.candidate_snps[self.ase_snps[ti]].haplotype = self.candidate_snps[self.ase_snps[ti]].haplotype * (-1);
                            }
                        }
                        for tk in 0..self.fragments.len() {
                            let rg: f64 = rng.gen();
                            if rg < random_flip_fraction as f64 {
                                self.fragments[tk].haplotag = self.fragments[tk].haplotag * (-1);
                            }
                        }
                    }
                    let prob = self.cross_optimize_ase_snps();
                    if prob > largest_prob {
                        largest_prob = prob;
                        self.save_best_configuration(&mut best_haplotype, &mut best_haplotag);
                    }
                    self.load_best_configuration(&best_haplotype, &best_haplotag);
                }
                self.load_best_configuration(&best_haplotype, &best_haplotag);
                max_iter -= 1;
            }
            self.load_best_configuration(&best_haplotype, &best_haplotag);
        }
    }

    pub fn add_phase_score(
        &mut self,
        min_allele_cnt: u32,
        min_homozygous_freq: f32,
        min_phase_score: f32,
    ) {
        // calculate phase score for each snp
        for ti in 0..self.candidate_snps.len() {
            let snp = &mut self.candidate_snps[ti];
            if snp.filter == true || snp.variant_type != 1 || snp.ase == true {
                continue;
            }
            if snp.snp_cover_fragments.len() == 0 {
                // no surranding haplotype links
                snp.single = true;
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
                if self.fragments[*k].num_hete_links < self.min_linkers {
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
                        assert_ne!(fe.ase_snp, true, "Error: phase for unexpected ase SNP.");
                        assert_ne!(fe.p, 0, "Error: phase for unexpected allele.");
                        ps.push(fe.p);
                        probs.push(fe.prob);
                        sigma.push(self.fragments[*k].haplotag);
                    }
                }
            }

            let mut phase_score = 0.0;

            if num_hap1 < min_allele_cnt || num_hap2 < min_allele_cnt {
                // filter SNPs with low allele count, unconfident phase
                phase_score = 0.0;
            } else {
                if sigma.len() > 0 {
                    phase_score = -10.0_f64 * (1.0 - SNPFrag::cal_delta_sigma_log(delta_i, &sigma, &ps, &probs)).log10(); // calaulate assignment score
                } else {
                    phase_score = 0.0; // all reads belong to unknown group, unconfident phase
                }
            }

            if phase_score < min_phase_score as f64 && snp.allele_freqs[0] > min_homozygous_freq && snp.alleles[0] != snp.reference && snp.filter == false {
                // transfer from heterozygous to homozygous
                snp.variant_type = 2;
                snp.ase = false;
            }

            let mut haplotype_allele_expression: [u32; 4] = [0, 0, 0, 0];   // hap1_ref, hap1_alt, hap2_ref, hap2_alt
            for k in 0..sigma.len() {
                if sigma[k] == 1 {
                    // hap1
                    if ps[k] == 1 {
                        haplotype_allele_expression[0] += 1; // hap1 allele1, reference allele
                    } else if ps[k] == -1 {
                        haplotype_allele_expression[1] += 1; // hap1 allele2, alternative allele
                    }
                } else if sigma[k] == -1 {
                    // hap2
                    if ps[k] == 1 {
                        haplotype_allele_expression[2] += 1; // hap2 allele1, reference allele
                    } else if ps[k] == -1 {
                        haplotype_allele_expression[3] += 1; // hap2 allele2, alternative allele
                    }
                }
            }
            // TODO: detect whether candidate with low phase score is a potential somatic mutation
            self.candidate_snps[ti].haplotype_expression = haplotype_allele_expression;
            self.candidate_snps[ti].phase_score = phase_score;
        }
    }

    pub fn assign_reads(&mut self, read_assignment_cutoff: f64) -> HashMap<String, i32> {
        let mut read_assignments: HashMap<String, i32> = HashMap::new();
        for k in 0..self.fragments.len() {
            let sigma_k = self.fragments[k].haplotag;
            let mut delta: Vec<i32> = Vec::new();
            let mut ps: Vec<i32> = Vec::new();
            let mut probs: Vec<f64> = Vec::new();
            for fe in self.fragments[k].list.iter() {
                if fe.ase_snp == true {
                    continue;
                }
                assert_ne!(fe.p, 0, "Error: phase for unexpected allele.");
                ps.push(fe.p);
                probs.push(fe.prob);
                delta.push(self.candidate_snps[fe.snp_idx].haplotype);
            }
            if sigma_k == 0 {
                // unasigned haplotag, cluster the read into unknown group
                self.fragments[k].assignment = 0;
                self.fragments[k].assignment_score = 0.0;
                read_assignments.insert(self.fragments[k].read_id.clone(), 0);
            } else {
                let mut q = 0.0;
                let mut qn = 0.0;
                if delta.len() > 0 {
                    q = SNPFrag::cal_sigma_delta_log(sigma_k, &delta, &ps, &probs);
                    qn = SNPFrag::cal_sigma_delta_log(sigma_k * (-1), &delta, &ps, &probs);
                }

                if (q - qn).abs() > read_assignment_cutoff {
                    if q > qn {
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
                        if sigma_k == 1 {
                            self.fragments[k].assignment = 2;
                            self.fragments[k].assignment_score = qn;
                            self.fragments[k].haplotag = -1;
                            read_assignments.insert(self.fragments[k].read_id.clone(), 2);
                        } else {
                            self.fragments[k].assignment = 1;
                            self.fragments[k].assignment_score = qn;
                            self.fragments[k].haplotag = 1;
                            read_assignments.insert(self.fragments[k].read_id.clone(), 1);
                        }
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

    pub fn assign_reads_ase(&mut self, read_assignment_cutoff: f64) -> HashMap<String, i32> {
        let mut read_assignments: HashMap<String, i32> = HashMap::new();
        for k in 0..self.fragments.len() {
            let sigma_k = self.fragments[k].haplotag;
            let mut delta_hete: Vec<i32> = Vec::new();
            let mut ps_hete: Vec<i32> = Vec::new();
            let mut probs_hete: Vec<f64> = Vec::new();
            let mut delta_ase: Vec<i32> = Vec::new();
            let mut ps_ase: Vec<i32> = Vec::new();
            let mut probs_ase: Vec<f64> = Vec::new();

            for fe in self.fragments[k].list.iter() {
                assert_ne!(fe.p, 0, "Error: phase for unexpected allele.");
                if fe.ase_snp == true {
                    ps_ase.push(fe.p);
                    probs_ase.push(fe.prob);
                    delta_ase.push(self.candidate_snps[fe.snp_idx].haplotype);
                } else {
                    ps_hete.push(fe.p);
                    probs_hete.push(fe.prob);
                    delta_hete.push(self.candidate_snps[fe.snp_idx].haplotype);
                }
            }
            if sigma_k == 0 {
                // unasigned haplotag, cluster the read into unknown group
                self.fragments[k].assignment = 0;
                self.fragments[k].assignment_score = 0.0;
                read_assignments.insert(self.fragments[k].read_id.clone(), 0);
            } else {
                let mut q_hete = 0.0;
                let mut qn_hete = 0.0;
                let mut q_ase = 0.0;
                let mut qn_ase = 0.0;
                if delta_hete.len() > 0 {
                    q_hete = SNPFrag::cal_sigma_delta_log(sigma_k, &delta_hete, &ps_hete, &probs_hete);
                    qn_hete = SNPFrag::cal_sigma_delta_log(
                        sigma_k * (-1),
                        &delta_hete,
                        &ps_hete,
                        &probs_hete,
                    );
                }
                if delta_ase.len() > 0 {
                    q_ase = SNPFrag::cal_sigma_delta_log(sigma_k, &delta_ase, &ps_ase, &probs_ase);
                    qn_ase = SNPFrag::cal_sigma_delta_log(
                        sigma_k * (-1),
                        &delta_ase,
                        &ps_ase,
                        &probs_ase,
                    );
                }

                if (q_hete - qn_hete).abs() > read_assignment_cutoff && (q_ase - qn_ase).abs() > read_assignment_cutoff {
                    // consider both ase and hete snps
                    if q_hete >= qn_hete && q_ase >= qn_ase {
                        // both ase and hete snps support the same haplotype
                        if sigma_k == 1 {
                            self.fragments[k].assignment = 1;
                            self.fragments[k].assignment_score = q_hete.max(q_ase);
                            read_assignments.insert(self.fragments[k].read_id.clone(), 1);
                        } else {
                            self.fragments[k].assignment = 2;
                            self.fragments[k].assignment_score = q_hete.max(q_ase);
                            read_assignments.insert(self.fragments[k].read_id.clone(), 2);
                        }
                    } else if q_hete < qn_hete && q_ase < qn_ase {
                        // both ase and hete snps support the opposite haplotype
                        if sigma_k == 1 {
                            self.fragments[k].assignment = 2;
                            self.fragments[k].assignment_score = qn_hete.max(qn_ase);
                            self.fragments[k].haplotag = -1;
                            read_assignments.insert(self.fragments[k].read_id.clone(), 2);
                        } else {
                            self.fragments[k].assignment = 1;
                            self.fragments[k].assignment_score = qn_hete.max(qn_ase);
                            self.fragments[k].haplotag = 1;
                            read_assignments.insert(self.fragments[k].read_id.clone(), 1);
                        }
                    } else if q_hete >= qn_hete && q_ase < qn_ase {
                        // hete snps and ase snps have conflict read assignment
                        if delta_hete.len() >= 2 {
                            // if more than 2 hete snps, use hete snps to assign the read
                            if sigma_k == 1 {
                                self.fragments[k].assignment = 1;
                                self.fragments[k].assignment_score = q_hete;
                                read_assignments.insert(self.fragments[k].read_id.clone(), 1);
                            } else {
                                self.fragments[k].assignment = 2;
                                self.fragments[k].assignment_score = q_hete;
                                read_assignments.insert(self.fragments[k].read_id.clone(), 2);
                            }
                        } else {
                            // if only 1 hete snp, use the higher probability to assign the read
                            if q_hete >= qn_ase {
                                if sigma_k == 1 {
                                    self.fragments[k].assignment = 1;
                                    self.fragments[k].assignment_score = q_hete;
                                    read_assignments.insert(self.fragments[k].read_id.clone(), 1);
                                } else {
                                    self.fragments[k].assignment = 2;
                                    self.fragments[k].assignment_score = q_hete;
                                    read_assignments.insert(self.fragments[k].read_id.clone(), 2);
                                }
                            } else {
                                if sigma_k == 1 {
                                    self.fragments[k].assignment = 2;
                                    self.fragments[k].assignment_score = qn_ase;
                                    self.fragments[k].haplotag = -1;
                                    read_assignments.insert(self.fragments[k].read_id.clone(), 2);
                                } else {
                                    self.fragments[k].assignment = 1;
                                    self.fragments[k].assignment_score = qn_ase;
                                    self.fragments[k].haplotag = 1;
                                    read_assignments.insert(self.fragments[k].read_id.clone(), 1);
                                }
                            }
                        }
                    } else if q_hete < qn_hete && q_ase >= qn_ase {
                        // hete snps and ase snps have conflict read assignment
                        if delta_hete.len() >= 2 {
                            // if more than 2 hete snps, use hete snps to assign the read
                            if sigma_k == 1 {
                                self.fragments[k].assignment = 2;
                                self.fragments[k].assignment_score = qn_hete;
                                self.fragments[k].haplotag = -1;
                                read_assignments.insert(self.fragments[k].read_id.clone(), 2);
                            } else {
                                self.fragments[k].assignment = 1;
                                self.fragments[k].assignment_score = qn_hete;
                                self.fragments[k].haplotag = 1;
                                read_assignments.insert(self.fragments[k].read_id.clone(), 1);
                            }
                        } else {
                            // if only 1 hete snp, use the higher probability to assign the read
                            if qn_hete >= q_ase {
                                if sigma_k == 1 {
                                    self.fragments[k].assignment = 2;
                                    self.fragments[k].assignment_score = qn_hete;
                                    self.fragments[k].haplotag = -1;
                                    read_assignments.insert(self.fragments[k].read_id.clone(), 2);
                                } else {
                                    self.fragments[k].assignment = 1;
                                    self.fragments[k].assignment_score = qn_hete;
                                    self.fragments[k].haplotag = 1;
                                    read_assignments.insert(self.fragments[k].read_id.clone(), 1);
                                }
                            } else {
                                if sigma_k == 1 {
                                    self.fragments[k].assignment = 1;
                                    self.fragments[k].assignment_score = q_ase;
                                    read_assignments.insert(self.fragments[k].read_id.clone(), 1);
                                } else {
                                    self.fragments[k].assignment = 2;
                                    self.fragments[k].assignment_score = q_ase;
                                    read_assignments.insert(self.fragments[k].read_id.clone(), 2);
                                }
                            }
                        }
                    } else {
                        // unknown which haplotype the read belongs to, cluster the read into unknown group
                        self.fragments[k].assignment = 0;
                        self.fragments[k].assignment_score = 0.0;
                        read_assignments.insert(self.fragments[k].read_id.clone(), 0);
                    }
                } else if (q_hete - qn_hete).abs() > read_assignment_cutoff {
                    // only consider hete snps
                    if q_hete > qn_hete {
                        if sigma_k == 1 {
                            self.fragments[k].assignment = 1;
                            self.fragments[k].assignment_score = q_hete;
                            read_assignments.insert(self.fragments[k].read_id.clone(), 1);
                        } else {
                            self.fragments[k].assignment = 2;
                            self.fragments[k].assignment_score = q_hete;
                            read_assignments.insert(self.fragments[k].read_id.clone(), 2);
                        }
                    } else {
                        if sigma_k == 1 {
                            self.fragments[k].assignment = 2;
                            self.fragments[k].assignment_score = qn_hete;
                            self.fragments[k].haplotag = -1;
                            read_assignments.insert(self.fragments[k].read_id.clone(), 2);
                        } else {
                            self.fragments[k].assignment = 1;
                            self.fragments[k].assignment_score = qn_hete;
                            self.fragments[k].haplotag = 1;
                            read_assignments.insert(self.fragments[k].read_id.clone(), 1);
                        }
                    }
                } else if (q_ase - qn_ase).abs() > read_assignment_cutoff {
                    // only consider ase snps
                    if q_ase > qn_ase {
                        if sigma_k == 1 {
                            self.fragments[k].assignment = 1;
                            self.fragments[k].assignment_score = q_ase;
                            read_assignments.insert(self.fragments[k].read_id.clone(), 1);
                        } else {
                            self.fragments[k].assignment = 2;
                            self.fragments[k].assignment_score = q_ase;
                            read_assignments.insert(self.fragments[k].read_id.clone(), 2);
                        }
                    } else {
                        if sigma_k == 1 {
                            self.fragments[k].assignment = 2;
                            self.fragments[k].assignment_score = qn_ase;
                            self.fragments[k].haplotag = -1;
                            read_assignments.insert(self.fragments[k].read_id.clone(), 2);
                        } else {
                            self.fragments[k].assignment = 1;
                            self.fragments[k].assignment_score = qn_ase;
                            self.fragments[k].haplotag = 1;
                            read_assignments.insert(self.fragments[k].read_id.clone(), 1);
                        }
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


    pub fn assign_phase_set(&mut self, min_qual_for_candidate: u32, min_phase_score: f32, ase_ps_cutoff: f32) -> HashMap<String, u32> {
        let mut phase_set: HashMap<String, u32> = HashMap::new();
        let mut graph: GraphMap<usize, Vec<usize>, Undirected> = GraphMap::new();  // node is index in candidate snp, edge is index in fragments
        // construct graph for hete snps
        for i in 0..self.candidate_snps.len() {
            let snp = &self.candidate_snps[i];
            if snp.variant_type == 0 || snp.variant_type == 2 || snp.variant_type == 3 {
                continue;
            }
            if snp.variant_type == 1 {
                // hete snps
                let mut is_low_qual = false;
                let mut is_dense = false;
                let mut is_rna_edit = false;
                let mut is_single_snp = false;
                let mut is_unconfident_phased_snp = false;
                let mut is_ase_snp = false;
                if snp.filter == true {
                    is_dense = true;
                }
                if snp.single == true {
                    is_single_snp = true;
                }
                if snp.ase == true {
                    is_ase_snp = true;
                }
                if !is_dense && !is_single_snp && snp.phase_score == 0.0 {
                    is_unconfident_phased_snp = true;
                }
                if is_dense || is_single_snp || is_unconfident_phased_snp || snp.haplotype == 0 {
                    continue;
                }
                // ase snp
                if is_ase_snp && snp.phase_score < ase_ps_cutoff as f64 {
                    continue;
                }
                // hete snp
                if !is_ase_snp {
                    if snp.variant_quality < min_qual_for_candidate as f64 {
                        continue;
                    }
                    if snp.phase_score < min_phase_score as f64 {
                        continue;
                    }
                }
                // ase snps > ase_ps_cutoff or hete snps > min_phase_score, construct graph
                graph.add_node(i);
            }
        }
        for k in 0..self.fragments.len() {
            let frag = &self.fragments[k];
            if frag.assignment == 0 { continue; }
            let mut node_snps = Vec::new();
            for fe in frag.list.iter() {
                if graph.contains_node(fe.snp_idx) {
                    node_snps.push(fe.snp_idx);
                }
            }
            if node_snps.len() >= 2 {
                for j in 0..node_snps.len() - 1 {
                    if !graph.contains_edge(node_snps[j], node_snps[j + 1]) {
                        graph.add_edge(node_snps[j], node_snps[j + 1], vec![k]);    // weight is a vector of fragment index, which is covered by the edge
                    } else {
                        graph.edge_weight_mut(node_snps[j], node_snps[j + 1]).unwrap().push(k);
                    }
                }
            }
        }
        let scc = kosaraju_scc(&graph);
        let region = self.region.clone().to_string();
        for component_nodes in scc.iter() {
            if component_nodes.len() <= 1 {
                continue;
            }
            let mut phase_id = 0;
            for node in component_nodes.iter() {
                if phase_id == 0 {
                    phase_id = (self.candidate_snps[*node].pos + 1) as u32;  // 1-based;
                }
                self.candidate_snps[*node].phase_set = phase_id;
                for edge in graph.edges(*node) {
                    let frag_idxes = edge.2;
                    for k in frag_idxes.iter() {
                        let fragment = &self.fragments[*k];
                        let read_id = fragment.read_id.clone();
                        if phase_set.contains_key(&read_id) {
                            continue;
                        }
                        phase_set.insert(read_id, phase_id);
                    }
                }
            }
        }
        return phase_set;
    }

    fn check_local_optimal_configuration(&self, used_for_haplotype: bool, used_for_haplotag: bool) {
        // check sigma
        if used_for_haplotag {
            for k in 0..self.fragments.len() {
                let sigma_k = self.fragments[k].haplotag;
                let mut delta: Vec<i32> = Vec::new();
                let mut ps: Vec<i32> = Vec::new();
                let mut probs: Vec<f64> = Vec::new();
                if sigma_k == 0 {
                    continue;
                }
                for fe in self.fragments[k].list.iter() {
                    if fe.ase_snp == true {
                        continue;
                    }
                    assert_ne!(fe.p, 0, "Error: phase for unexpected allele.");
                    ps.push(fe.p);
                    probs.push(fe.prob);
                    delta.push(self.candidate_snps[fe.snp_idx].haplotype);
                }
                if delta.len() == 0 {
                    continue;
                }
                let q = SNPFrag::cal_sigma_delta_log(sigma_k, &delta, &ps, &probs);
                let qn = SNPFrag::cal_sigma_delta_log(sigma_k * (-1), &delta, &ps, &probs);
                // println!("q:{}, qn:{}", q, qn);
                assert!(q >= qn, "{} Error: read assignment is not local optimal. {}->{}\n{:?}\ndelta:{:?}\nps:{:?}\nprobs:{:?}\nsigma:{}\n{:?}\n{:?}", k, q, qn, self.region, delta, ps, probs, sigma_k, used_for_haplotype, used_for_haplotag);
            }
        }

        // check delta
        if used_for_haplotype {
            for i in self.hete_snps.iter() {
                let delta_i = self.candidate_snps[*i].haplotype;
                let mut sigma: Vec<i32> = Vec::new();
                let mut ps: Vec<i32> = Vec::new();
                let mut probs: Vec<f64> = Vec::new();
                for k in self.candidate_snps[*i].snp_cover_fragments.iter() {
                    if self.fragments[*k].haplotag == 0 {
                        continue;
                    }
                    // k is fragment index
                    for fe in self.fragments[*k].list.iter() {
                        if fe.snp_idx == *i {
                            assert_ne!(fe.ase_snp, true, "Error: phase for unexpected ase SNP.");
                            assert_ne!(fe.p, 0, "Error: phase for unexpected allele.");
                            ps.push(fe.p);
                            probs.push(fe.prob);
                            sigma.push(self.fragments[*k].haplotag);
                        }
                    }
                }
                if sigma.len() == 0 {
                    continue;
                }
                let q = SNPFrag::cal_delta_sigma_log(delta_i, &sigma, &ps, &probs);
                let qn = SNPFrag::cal_delta_sigma_log(delta_i * (-1), &sigma, &ps, &probs);
                assert!(q >= qn, "{} Error: phase is not local optimal. {}->{}\n{:?}\nsigma:{:?}\nps:{:?}\nprobs:{:?}\ndelta:{}\n{:?}\n{:?}", i, q, qn, self.region, sigma, ps, probs, delta_i, used_for_haplotype, used_for_haplotag);
            }
        }
    }

    fn check_local_optimal_configuration_ase(
        &self,
        used_for_haplotype: bool,
        used_for_haplotag: bool,
    ) {
        // check sigma
        if used_for_haplotag {
            for k in 0..self.fragments.len() {
                let sigma_k = self.fragments[k].haplotag;
                let mut delta: Vec<i32> = Vec::new();
                let mut ps: Vec<i32> = Vec::new();
                let mut probs: Vec<f64> = Vec::new();
                for fe in self.fragments[k].list.iter() {
                    ps.push(fe.p);
                    probs.push(fe.prob);
                    delta.push(self.candidate_snps[fe.snp_idx].haplotype);
                }
                if delta.len() == 0 {
                    continue;
                }
                let q = SNPFrag::cal_sigma_delta_log(sigma_k, &delta, &ps, &probs);
                let qn = SNPFrag::cal_sigma_delta_log(sigma_k * (-1), &delta, &ps, &probs);
                // println!("q:{}, qn:{}", q, qn);
                assert!(q >= qn, "{} Error: read assignment is not local optimal. {}->{}\n{:?}\ndelta:{:?}\nps:{:?}\nprobs:{:?}\nsigma:{}\n{:?}\n{:?}", k, q, qn, self.region, delta, ps, probs, sigma_k, used_for_haplotype, used_for_haplotag);
            }
        }

        // check delta
        if used_for_haplotype {
            for i in self.ase_snps.iter() {
                let delta_i = self.candidate_snps[*i].haplotype;
                let mut sigma: Vec<i32> = Vec::new();
                let mut ps: Vec<i32> = Vec::new();
                let mut probs: Vec<f64> = Vec::new();
                for k in self.candidate_snps[*i].snp_cover_fragments.iter() {
                    // k is fragment index
                    for fe in self.fragments[*k].list.iter() {
                        if fe.snp_idx == *i {
                            assert_ne!(fe.p, 0, "Error: phase for unexpected allele.");
                            ps.push(fe.p);
                            probs.push(fe.prob);
                            sigma.push(self.fragments[*k].haplotag);
                        }
                    }
                }
                if sigma.len() == 0 {
                    continue;
                }
                let q = SNPFrag::cal_delta_sigma_log(delta_i, &sigma, &ps, &probs);
                let qn = SNPFrag::cal_delta_sigma_log(delta_i * (-1), &sigma, &ps, &probs);
                assert!(q >= qn, "{} Error: phase is not local optimal. {}->{}\n{:?}\nsigma:{:?}\nps:{:?}\nprobs:{:?}\ndelta:{}\n{:?}\n{:?}", i, q, qn, self.region, sigma, ps, probs, delta_i, used_for_haplotype, used_for_haplotag);
            }
        }
    }

    pub fn rescue_ase_snps(&mut self) {
        // use surrounding heterozygous snps to evaluate ase snps
        // TODO: add phase score still has problem
        for i in self.ase_hete_snps.iter() {
            if self.candidate_snps[*i].ase == true {
                let mut sigma: Vec<i32> = Vec::new();
                let mut ps: Vec<i32> = Vec::new();
                let mut probs: Vec<f64> = Vec::new();
                let mut num_hap1 = 0;
                let mut num_hap2 = 0;
                for k in self.candidate_snps[*i].snp_cover_fragments.iter() {
                    // k is fragment index
                    if self.fragments[*k].assignment == 0 {
                        continue;
                    }
                    for fe in self.fragments[*k].list.iter() {
                        if fe.snp_idx == *i {
                            if self.fragments[*k].assignment == 1 {
                                num_hap1 += 1;
                            } else if self.fragments[*k].assignment == 2 {
                                num_hap2 += 1;
                            }
                            assert_ne!(fe.p, 0, "Error: phase for unexpected allele.");
                            ps.push(fe.p);
                            probs.push(fe.prob);
                            sigma.push(self.fragments[*k].haplotag);
                        }
                    }
                }
                if num_hap1 < 3 || num_hap2 < 3 {
                    continue;
                }
                if sigma.len() > 10 {
                    let phase_score1 = -10.0_f64 * (1.0 - SNPFrag::cal_delta_sigma_log(1, &sigma, &ps, &probs)).log10();
                    let phase_score2 = -10.0_f64 * (1.0 - SNPFrag::cal_delta_sigma_log(-1, &sigma, &ps, &probs)).log10();
                    if f64::max(phase_score1, phase_score2) > 40.0 {
                        if phase_score1 > phase_score2 {
                            self.candidate_snps[*i].haplotype = 1;
                        } else {
                            self.candidate_snps[*i].haplotype = -1;
                        }
                        println!(
                            "Rescue ASE SNP: {} {} {} {}",
                            String::from_utf8(self.candidate_snps[*i].chromosome.clone()).unwrap(),
                            self.candidate_snps[*i].pos,
                            phase_score1,
                            phase_score2
                        );
                        // self.candidate_snps[*i].ase = false;
                    }
                }
            }
        }
    }

    pub fn rescue_ase_snps_v2(
        &mut self,
        ase_allele_cnt_cutoff: u32,
        ase_ps_count_cutoff: u32,
        ase_ps_cutoff: f32,
    ) {
        // keep the previously calculated phasing for heterozygous snps and only flip the ase snps to get the phasing of ase snps
        for i in self.ase_snps.iter() {
            assert!(self.candidate_snps[*i].ase == true, "rescue not ase snp.");
            let mut sigma: Vec<i32> = Vec::new();
            let mut ps: Vec<i32> = Vec::new();
            let mut probs: Vec<f64> = Vec::new();
            let mut num_hap1 = 0;
            let mut num_hap2 = 0;
            for k in self.candidate_snps[*i].snp_cover_fragments.iter() {
                // k is fragment index
                if self.fragments[*k].assignment == 0 {
                    continue;
                }
                for fe in self.fragments[*k].list.iter() {
                    if fe.snp_idx == *i {
                        if self.fragments[*k].assignment == 1 {
                            num_hap1 += 1;
                        } else if self.fragments[*k].assignment == 2 {
                            num_hap2 += 1;
                        }
                        assert_ne!(fe.p, 0, "Error: phase for unexpected allele.");
                        ps.push(fe.p);
                        probs.push(fe.prob);
                        sigma.push(self.fragments[*k].haplotag);
                    }
                }
            }
            if num_hap1 < ase_allele_cnt_cutoff || num_hap2 < ase_allele_cnt_cutoff {
                // each haplotype should have at least 3 reads
                continue;
            }
            if sigma.len() > ase_ps_count_cutoff as usize {
                // each site should have at least 10 reads
                let phase_score1 = -10.0_f64 * (1.0 - SNPFrag::cal_delta_sigma_log(1, &sigma, &ps, &probs)).log10();
                let phase_score2 = -10.0_f64 * (1.0 - SNPFrag::cal_delta_sigma_log(-1, &sigma, &ps, &probs)).log10();
                if f64::max(phase_score1, phase_score2) > ase_ps_cutoff as f64 {
                    if phase_score1 > phase_score2 {
                        self.candidate_snps[*i].haplotype = 1;
                        self.candidate_snps[*i].phase_score = phase_score1;
                        self.candidate_snps[*i].variant_type = 1;
                    } else {
                        self.candidate_snps[*i].haplotype = -1;
                        self.candidate_snps[*i].phase_score = phase_score2;
                        self.candidate_snps[*i].variant_type = 1;
                    }
                }
                let mut haplotype_allele_expression: [u32; 4] = [0, 0, 0, 0];   // hap1_ref, hap1_alt, hap2_ref, hap2_alt
                for k in 0..sigma.len() {
                    if sigma[k] == 1 {
                        // hap1
                        if ps[k] == 1 {
                            haplotype_allele_expression[0] += 1;
                        } else if ps[k] == -1 {
                            haplotype_allele_expression[1] += 1;
                        }
                    } else if sigma[k] == -1 {
                        // hap2
                        if ps[k] == 1 {
                            haplotype_allele_expression[2] += 1;
                        } else if ps[k] == -1 {
                            haplotype_allele_expression[3] += 1;
                        }
                    }
                }
                self.candidate_snps[*i].haplotype_expression = haplotype_allele_expression;
                // TODO: detect whether candidate with low phase score is a potential somatic mutation
            }
        }
    }

    pub fn phased_output_vcf(
        &mut self,
        min_phase_score: f32,
        min_homozygous_freq: f32,
        output_phasing: bool,
        min_qual_for_candidate: u32,
        min_qual_for_singlesnp_rnaedit: u32,
    ) -> Vec<VCFRecord> {
        let mut records: Vec<VCFRecord> = Vec::new();

        // output heterozygous SNPs
        // assert_eq!(self.haplotype.len(), self.snps.len());
        for i in 0..self.candidate_snps.len() {
            let snp = &self.candidate_snps[i];
            let hp = snp.haplotype;
            // if snp.ase == true {
            //     continue;
            // }
            if snp.filter == true && snp.rna_editing == false {
                // Dense SNP. filter==true && rna_editing==false: dense SNP. filter==true && rna_editing==true: rna editing site
                let mut rd: VCFRecord = VCFRecord::default();
                rd.chromosome = snp.chromosome.clone();
                rd.position = snp.pos as u64 + 1; // position in vcf format is 1-based
                rd.reference = vec![snp.reference as u8];
                rd.id = vec!['.' as u8];
                if snp.variant_type == 1 {
                    if snp.alleles[0] != snp.reference && snp.alleles[1] == snp.reference {
                        rd.alternative = vec![vec![snp.alleles[0] as u8]];
                        rd.genotype = format!(
                            "{}:{}:{}:{:.2}",
                            "0/1", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[0]
                        );
                    } else if snp.alleles[1] != snp.reference && snp.alleles[0] == snp.reference {
                        rd.alternative = vec![vec![snp.alleles[1] as u8]];
                        rd.genotype = format!(
                            "{}:{}:{}:{:.2}",
                            "0/1", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[1]
                        );
                    } else {
                        rd.alternative = vec![vec![snp.alleles[0] as u8], vec![snp.alleles[1] as u8]];
                        rd.genotype = format!(
                            "{}:{}:{}:{:.2},{:.2}",
                            "1/2",
                            snp.genotype_quality as i32,
                            snp.depth,
                            snp.allele_freqs[0],
                            snp.allele_freqs[1]
                        );
                    }
                } else if snp.variant_type == 2 {
                    rd.alternative = vec![vec![snp.alleles[0] as u8]];
                    rd.genotype = format!(
                        "{}:{}:{}:{:.2}",
                        "1/1", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[0]
                    );
                } else if snp.variant_type == 3 {
                    rd.alternative = vec![vec![snp.alleles[0] as u8], vec![snp.alleles[1] as u8]];
                    rd.genotype = format!(
                        "{}:{}:{}:{:.2},{:.2}",
                        "1/2",
                        snp.genotype_quality as i32,
                        snp.depth,
                        snp.allele_freqs[0],
                        snp.allele_freqs[1]
                    );
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
                rd.position = snp.pos as u64 + 1; // position in vcf format is 1-based
                rd.reference = vec![snp.reference as u8];
                rd.id = vec!['.' as u8];
                rd.alternative = vec![vec![snp.alleles[0] as u8], vec![snp.alleles[1] as u8]];
                rd.qual = snp.variant_quality as i32;
                rd.genotype = format!(
                    "{}:{}:{}:{:.2},{:.2}",
                    "1/2",
                    snp.genotype_quality as i32,
                    snp.depth,
                    snp.allele_freqs[0],
                    snp.allele_freqs[1]
                );
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
                rd.position = snp.pos as u64 + 1; // position in vcf format is 1-based
                rd.reference = vec![snp.reference as u8];
                rd.id = vec!['.' as u8];
                rd.alternative = vec![vec![snp.alleles[0] as u8]];
                rd.qual = snp.variant_quality as i32;
                rd.genotype = format!(
                    "{}:{}:{}:{:.2}",
                    "1/1", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[0]
                );
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
                    rd.position = snp.pos as u64 + 1; // position in vcf format is 1-based
                    rd.reference = vec![snp.reference as u8];
                    rd.id = vec!['.' as u8];
                    rd.alternative = vec![vec![snp.alleles[0] as u8]];
                    rd.qual = snp.variant_quality as i32;
                    rd.genotype = format!(
                        "{}:{}:{}:{:.2}:{:.2}",
                        "1/1",
                        snp.genotype_quality as i32,
                        snp.depth,
                        snp.allele_freqs[0],
                        snp.phase_score
                    );
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
                    rd.position = snp.pos as u64 + 1; // position in vcf format is 1-based
                    rd.reference = vec![snp.reference as u8];
                    rd.id = vec!['.' as u8];

                    // output information containing phasing information
                    if output_phasing {
                        if snp.rna_editing == true || snp.phase_score == 0.0 {
                            // rna editing no phase information or phase_score == 0: single SNP, no phasing information
                            if snp.alleles[0] == snp.reference {
                                rd.alternative = vec![vec![snp.alleles[1] as u8]];
                                rd.genotype = format!(
                                    "{}:{}:{}:{:.2}",
                                    "0/1",
                                    snp.genotype_quality as i32,
                                    snp.depth,
                                    snp.allele_freqs[1]
                                );
                            } else if snp.alleles[1] == snp.reference {
                                rd.alternative = vec![vec![snp.alleles[0] as u8]];
                                rd.genotype = format!(
                                    "{}:{}:{}:{:.2}",
                                    "0/1",
                                    snp.genotype_quality as i32,
                                    snp.depth,
                                    snp.allele_freqs[0]
                                );
                            } else {
                                rd.alternative = vec![vec![snp.alleles[0] as u8], vec![snp.alleles[1] as u8]];
                                rd.genotype = format!(
                                    "{}:{}:{}:{:.2},{:.2}",
                                    "1/2",
                                    snp.genotype_quality as i32,
                                    snp.depth,
                                    snp.allele_freqs[0],
                                    snp.allele_freqs[1]
                                );
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
                                    rd.genotype = format!(
                                        "{}:{}:{}:{:.2}:{:.2}:{},{},{},{}",
                                        "0|1",
                                        snp.genotype_quality as i32,
                                        snp.depth,
                                        snp.allele_freqs[1],
                                        snp.phase_score,
                                        snp.haplotype_expression[0],
                                        snp.haplotype_expression[1],
                                        snp.haplotype_expression[2],
                                        snp.haplotype_expression[3]
                                    );
                                } else if hp == 1 {
                                    rd.genotype = format!(
                                        "{}:{}:{}:{:.2}:{:.2}:{},{},{},{}",
                                        "1|0",
                                        snp.genotype_quality as i32,
                                        snp.depth,
                                        snp.allele_freqs[1],
                                        snp.phase_score,
                                        snp.haplotype_expression[0],
                                        snp.haplotype_expression[1],
                                        snp.haplotype_expression[2],
                                        snp.haplotype_expression[3]
                                    );
                                }
                            } else if snp.alleles[1] == snp.reference {
                                rd.alternative = vec![vec![snp.alleles[0] as u8]];
                                rd.qual = snp.variant_quality as i32;
                                if hp == -1 {
                                    rd.genotype = format!(
                                        "{}:{}:{}:{:.2}:{:.2}:{},{},{},{}",
                                        "0|1",
                                        snp.genotype_quality as i32,
                                        snp.depth,
                                        snp.allele_freqs[0],
                                        snp.phase_score,
                                        snp.haplotype_expression[0],
                                        snp.haplotype_expression[1],
                                        snp.haplotype_expression[2],
                                        snp.haplotype_expression[3]
                                    );
                                } else if hp == 1 {
                                    rd.genotype = format!(
                                        "{}:{}:{}:{:.2}:{:.2}:{},{},{},{}",
                                        "1|0",
                                        snp.genotype_quality as i32,
                                        snp.depth,
                                        snp.allele_freqs[0],
                                        snp.phase_score,
                                        snp.haplotype_expression[0],
                                        snp.haplotype_expression[1],
                                        snp.haplotype_expression[2],
                                        snp.haplotype_expression[3]
                                    );
                                }
                            } else {
                                rd.alternative = vec![vec![snp.alleles[0] as u8], vec![snp.alleles[1] as u8]];
                                rd.qual = snp.variant_quality as i32;
                                rd.genotype = format!(
                                    "{}:{}:{}:{:.2},{:.2}:{:.2}:{},{},{},{}",
                                    "1|2",
                                    snp.genotype_quality as i32,
                                    snp.depth,
                                    snp.allele_freqs[0],
                                    snp.allele_freqs[1],
                                    snp.phase_score,
                                    snp.haplotype_expression[0],
                                    snp.haplotype_expression[1],
                                    snp.haplotype_expression[2],
                                    snp.haplotype_expression[3]
                                );
                            }
                            if snp.phase_score < min_phase_score as f64 || snp.variant_quality < min_qual_for_candidate as f64 {
                                // not confident phase or low variant quality
                                rd.filter = "LowQual".to_string().into_bytes();
                            } else {
                                rd.filter = "PASS".to_string().into_bytes();
                            }
                            if snp.ase == true {
                                rd.info = format!("RDS={}", "ase_snp").to_string().into_bytes();
                            } else {
                                rd.info = "RDS=.".to_string().into_bytes();
                            }
                            rd.format = "GT:GQ:DP:AF:PQ:AE".to_string().into_bytes();
                        }
                    }

                    // output information does not contain phasing information, but still filtered by phase score
                    if !output_phasing {
                        if snp.rna_editing == true || snp.phase_score == 0.0 {
                            // rna editing no phase information or phase_score == 0: single SNP, no phasing information
                            if snp.alleles[0] == snp.reference {
                                rd.alternative = vec![vec![snp.alleles[1] as u8]];
                                rd.genotype = format!(
                                    "{}:{}:{}:{:.2}",
                                    "0/1",
                                    snp.genotype_quality as i32,
                                    snp.depth,
                                    snp.allele_freqs[1]
                                );
                            } else if snp.alleles[1] == snp.reference {
                                rd.alternative = vec![vec![snp.alleles[0] as u8]];
                                rd.genotype = format!(
                                    "{}:{}:{}:{:.2}",
                                    "0/1",
                                    snp.genotype_quality as i32,
                                    snp.depth,
                                    snp.allele_freqs[0]
                                );
                            } else {
                                rd.alternative = vec![vec![snp.alleles[0] as u8], vec![snp.alleles[1] as u8]];
                                rd.genotype = format!(
                                    "{}:{}:{}:{:.2},{:.2}",
                                    "1/2",
                                    snp.genotype_quality as i32,
                                    snp.depth,
                                    snp.allele_freqs[0],
                                    snp.allele_freqs[1]
                                );
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
                                if snp.ase == true {
                                    rd.info = format!("RDS={}", "ase_snp").to_string().into_bytes();
                                } else {
                                    rd.info = "RDS=.".to_string().into_bytes();
                                }
                                // rd.info = "RDS=.".to_string().into_bytes();
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
                                rd.genotype = format!(
                                    "{}:{}:{}:{:.2}",
                                    "0/1",
                                    snp.genotype_quality as i32,
                                    snp.depth,
                                    snp.allele_freqs[1]
                                );
                            } else if snp.alleles[1] == snp.reference {
                                rd.alternative = vec![vec![snp.alleles[0] as u8]];
                                rd.qual = snp.variant_quality as i32;
                                rd.genotype = format!(
                                    "{}:{}:{}:{:.2}",
                                    "0/1",
                                    snp.genotype_quality as i32,
                                    snp.depth,
                                    snp.allele_freqs[0]
                                );
                            } else {
                                rd.alternative = vec![vec![snp.alleles[0] as u8], vec![snp.alleles[1] as u8]];
                                rd.qual = snp.variant_quality as i32;
                                rd.genotype = format!(
                                    "{}:{}:{}:{:.2},{:.2}",
                                    "1/2",
                                    snp.genotype_quality as i32,
                                    snp.depth,
                                    snp.allele_freqs[0],
                                    snp.allele_freqs[1]
                                );
                            }
                            if snp.phase_score < min_phase_score as f64 || snp.variant_quality < min_qual_for_candidate as f64 {
                                // not confident phase or low variant quality
                                rd.filter = "LowQual".to_string().into_bytes();
                            } else {
                                rd.filter = "PASS".to_string().into_bytes();
                            }
                            if snp.ase == true {
                                rd.info = format!("RDS={}", "ase_snp").to_string().into_bytes();
                            } else {
                                rd.info = "RDS=.".to_string().into_bytes();
                            }
                            // rd.info = "RDS=.".to_string().into_bytes();
                            rd.format = "GT:GQ:DP:AF:PQ".to_string().into_bytes();
                        }
                    }
                    records.push(rd);
                }
            } else if snp.variant_type == 0 {
                // homo ref, ase snp
                continue;
            } else {
                println!("Unknown variant type: {:?}:{:?}", snp.chromosome, snp.pos);
                continue;
            }
        }
        return records;
    }

    pub fn phased_output_vcf2(&mut self,
                              min_phase_score: f32,
                              ase_ps_cutoff: f32,
                              min_qual_for_candidate: u32,
                              min_qual_for_singlesnp_rnaedit: u32) -> Vec<VCFRecord> {
        let mut records: Vec<VCFRecord> = Vec::new();
        for i in 0..self.candidate_snps.len() {
            let snp = &self.candidate_snps[i];
            if snp.variant_type == 0 {
                continue;
            }
            if snp.variant_type == 1 {
                let mut is_dense = false;
                let mut is_rna_edit = false;
                let mut is_single_snp = false;
                let mut is_ase_snp = false;

                // dense snp?
                if snp.filter == true && snp.rna_editing == false {
                    is_dense = true;
                }

                // rna edit?
                if snp.rna_editing == true {
                    is_rna_edit = true;
                }

                // single snp?
                if snp.filter == false && snp.rna_editing == false && snp.phase_score == 0.0 {
                    is_single_snp = true;
                }

                // ase snp?
                if snp.ase == true {
                    is_ase_snp = true;
                }

                let mut rd: VCFRecord = VCFRecord::default();
                rd.chromosome = snp.chromosome.clone();
                rd.position = snp.pos as u64 + 1; // position in vcf format is 1-based
                rd.id = vec!['.' as u8];
                rd.reference = vec![snp.reference as u8];
                if snp.alleles[0] != snp.reference {
                    rd.alternative = vec![vec![snp.alleles[0] as u8]];
                } else if snp.alleles[1] != snp.reference {
                    rd.alternative = vec![vec![snp.alleles[1] as u8]];
                }
                rd.qual = snp.variant_quality as i32;
                if is_dense {
                    rd.filter = "dn".to_string().into_bytes();
                    rd.info = format!("RDS={}", "dense_snp").to_string().into_bytes();
                }
                if is_rna_edit {
                    if snp.variant_quality < min_qual_for_singlesnp_rnaedit as f64 {
                        rd.filter = "RnaEdit".to_string().into_bytes();
                        rd.info = format!("RDS={}", "rna_editing").to_string().into_bytes();
                    } else {
                        rd.filter = "PASS".to_string().into_bytes();
                        rd.info = format!("RDS={}", "rna_editing").to_string().into_bytes();
                    }
                }
                if is_single_snp {
                    if snp.variant_quality < min_qual_for_singlesnp_rnaedit as f64 {
                        rd.filter = "LowQual".to_string().into_bytes();
                        rd.info = format!("RDS={}", "single_snp").to_string().into_bytes();
                    } else {
                        rd.filter = "PASS".to_string().into_bytes();
                        rd.info = format!("RDS={}", "single_snp").to_string().into_bytes();
                    }
                }
                if is_ase_snp {
                    if snp.phase_score < ase_ps_cutoff as f64 {
                        rd.filter = "LowQual".to_string().into_bytes();
                        rd.info = format!("RDS={}", "ase_snp").to_string().into_bytes();
                    } else {
                        rd.filter = "PASS".to_string().into_bytes();
                        rd.info = format!("RDS={}", "ase_snp").to_string().into_bytes();
                    }
                }
                if !is_dense && !is_rna_edit && !is_single_snp && !is_ase_snp {
                    if snp.variant_quality < min_qual_for_candidate as f64 || snp.phase_score < min_phase_score as f64 {
                        rd.filter = "LowQual".to_string().into_bytes();
                        rd.info = "RDS=.".to_string().into_bytes();
                    } else {
                        rd.filter = "PASS".to_string().into_bytes();
                        rd.info = "RDS=.".to_string().into_bytes();
                    }
                }

                if is_dense || is_rna_edit || is_single_snp || snp.phase_score == 0.0 || snp.haplotype == 0 {
                    rd.genotype = format!(
                        "{}:{}:{}:{:.2}",
                        "0/1",
                        snp.genotype_quality as i32,
                        snp.depth,
                        snp.allele_freqs[1]
                    );
                    rd.format = "GT:GQ:DP:AF".to_string().into_bytes();
                } else {
                    let mut af = 0.0;
                    if snp.alleles[0] == snp.reference {
                        af = snp.allele_freqs[1];
                    } else if snp.alleles[1] == snp.reference {
                        af = snp.allele_freqs[0];
                    } else {
                        panic!("Error: unexpected allele. ref: {}, alt1: {}, alt2: {}, {}:{}", snp.reference, snp.alleles[0], snp.alleles[1], String::from_utf8_lossy(&snp.chromosome), snp.pos);
                    }
                    if snp.haplotype == -1 {
                        rd.genotype = format!(
                            "{}:{}:{}:{:.2}:{:.2}:{},{},{},{}",
                            "0|1",
                            snp.genotype_quality as i32,
                            snp.depth,
                            af,
                            snp.phase_score,
                            snp.haplotype_expression[0],
                            snp.haplotype_expression[1],
                            snp.haplotype_expression[2],
                            snp.haplotype_expression[3]
                        );
                    } else if snp.haplotype == 1 {
                        rd.genotype = format!(
                            "{}:{}:{}:{:.2}:{:.2}:{},{},{},{}",
                            "1|0",
                            snp.genotype_quality as i32,
                            snp.depth,
                            af,
                            snp.phase_score,
                            snp.haplotype_expression[0],
                            snp.haplotype_expression[1],
                            snp.haplotype_expression[2],
                            snp.haplotype_expression[3]
                        );
                    }
                    rd.format = "GT:GQ:DP:AF:PQ:AE".to_string().into_bytes();
                }
                records.push(rd);
            }
            if snp.variant_type == 2 {
                let mut is_dense = false;
                let mut is_rna_edit = false;

                if snp.filter == true && snp.rna_editing == false {
                    is_dense = true;
                }
                if snp.rna_editing == true {
                    is_rna_edit = true;
                }

                let mut rd: VCFRecord = VCFRecord::default();
                rd.chromosome = snp.chromosome.clone();
                rd.position = snp.pos as u64 + 1; // position in vcf format is 1-based
                rd.id = vec!['.' as u8];
                rd.reference = vec![snp.reference as u8];
                rd.alternative = vec![vec![snp.alleles[0] as u8]];
                rd.qual = snp.variant_quality as i32;

                if is_dense {
                    rd.filter = "dn".to_string().into_bytes();
                    rd.info = format!("RDS={}", "dense_snp").to_string().into_bytes();
                }
                if is_rna_edit {
                    if snp.variant_quality < min_qual_for_singlesnp_rnaedit as f64 {
                        rd.filter = "RnaEdit".to_string().into_bytes();
                        rd.info = format!("RDS={}", "rna_editing").to_string().into_bytes();
                    } else {
                        rd.filter = "PASS".to_string().into_bytes();
                        rd.info = format!("RDS={}", "rna_editing").to_string().into_bytes();
                    }
                }
                if !is_dense && !is_rna_edit {
                    if snp.variant_quality < min_qual_for_candidate as f64 {
                        rd.filter = "LowQual".to_string().into_bytes();
                        rd.info = "RDS=.".to_string().into_bytes();
                    } else {
                        rd.filter = "PASS".to_string().into_bytes();
                        rd.info = "RDS=.".to_string().into_bytes();
                    }
                }
                rd.genotype = format!(
                    "{}:{}:{}:{:.2}",
                    "1/1",
                    snp.genotype_quality as i32,
                    snp.depth,
                    snp.allele_freqs[0]
                );
                rd.format = "GT:GQ:DP:AF".to_string().into_bytes();
                records.push(rd);
            }
            if snp.variant_type == 3 {
                let mut is_dense = false;
                let mut is_rna_edit = false;

                if snp.filter == true && snp.rna_editing == false {
                    is_dense = true;
                }
                if snp.rna_editing == true {
                    is_rna_edit = true;
                }

                let mut rd: VCFRecord = VCFRecord::default();
                rd.chromosome = snp.chromosome.clone();
                rd.position = snp.pos as u64 + 1; // position in vcf format is 1-based
                rd.id = vec!['.' as u8];
                rd.reference = vec![snp.reference as u8];
                rd.alternative = vec![vec![snp.alleles[0] as u8], vec![snp.alleles[1] as u8]];
                rd.qual = snp.variant_quality as i32;

                if is_dense {
                    rd.filter = "dn".to_string().into_bytes();
                    rd.info = format!("RDS={}", "dense_snp").to_string().into_bytes();
                }
                if is_rna_edit {
                    if snp.variant_quality < min_qual_for_singlesnp_rnaedit as f64 {
                        rd.filter = "RnaEdit".to_string().into_bytes();
                        rd.info = format!("RDS={}", "rna_editing").to_string().into_bytes();
                    } else {
                        rd.filter = "PASS".to_string().into_bytes();
                        rd.info = format!("RDS={}", "rna_editing").to_string().into_bytes();
                    }
                }
                if !is_dense && !is_rna_edit {
                    if snp.variant_quality < min_qual_for_candidate as f64 {
                        rd.filter = "LowQual".to_string().into_bytes();
                        rd.info = "RDS=.".to_string().into_bytes();
                    } else {
                        rd.filter = "PASS".to_string().into_bytes();
                        rd.info = "RDS=.".to_string().into_bytes();
                    }
                }
                rd.genotype = format!(
                    "{}:{}:{}:{:.2},{:.2}",
                    "1/2",
                    snp.genotype_quality as i32,
                    snp.depth,
                    snp.allele_freqs[0],
                    snp.allele_freqs[1]
                );
                rd.format = "GT:GQ:DP:AF".to_string().into_bytes();
                records.push(rd);
            }
            if snp.variant_type != 0 && snp.variant_type != 1 && snp.variant_type != 2 && snp.variant_type != 3 {
                panic!("Unknown variant type: {} {}:{}", snp.variant_type, String::from_utf8_lossy(&snp.chromosome), snp.pos);
            }
        }
        return records;
    }

    pub fn phased_output_vcf3(&mut self,
                              min_phase_score: f32,
                              ase_ps_cutoff: f32,
                              min_qual_for_candidate: u32,
                              min_qual_for_singlesnp_rnaedit: u32) -> Vec<VCFRecord> {
        let mut records: Vec<VCFRecord> = Vec::new();
        for i in 0..self.candidate_snps.len() {
            let snp = &self.candidate_snps[i];
            if snp.variant_type == 0 {
                continue;
            }
            if snp.variant_type == 1 {
                let mut is_dense = false;
                let mut is_rna_edit = false;
                let mut is_single_snp = false;
                let mut is_unconfident_phased_snp = false;
                let mut is_ase_snp = false;

                // dense snp?
                if snp.filter == true {
                    is_dense = true;
                }

                // rna edit?
                if snp.rna_editing == true {
                    is_rna_edit = true;
                }

                // single snp?
                if snp.single == true {
                    is_single_snp = true;
                }

                // unconfident phase?
                if !is_dense && !is_single_snp && snp.phase_score == 0.0 {
                    is_unconfident_phased_snp = true;
                }

                // ase snp?
                if snp.ase == true {
                    is_ase_snp = true;
                }

                let mut rd: VCFRecord = VCFRecord::default();
                rd.chromosome = snp.chromosome.clone();
                rd.position = snp.pos as u64 + 1; // position in vcf format is 1-based
                rd.id = vec!['.' as u8];
                rd.reference = vec![snp.reference as u8];
                if snp.alleles[0] != snp.reference {
                    rd.alternative = vec![vec![snp.alleles[0] as u8]];
                } else if snp.alleles[1] != snp.reference {
                    rd.alternative = vec![vec![snp.alleles[1] as u8]];
                }
                rd.qual = snp.variant_quality as i32;
                if is_dense {
                    rd.filter = "dn".to_string().into_bytes();
                    rd.info = format!("RDS={}", "dense_snp").to_string().into_bytes();
                }
                if is_single_snp {
                    if snp.variant_quality < min_qual_for_singlesnp_rnaedit as f64 {
                        rd.filter = "LowQual".to_string().into_bytes();
                        rd.info = format!("RDS={}", "single_snp").to_string().into_bytes();
                    } else {
                        rd.filter = "PASS".to_string().into_bytes();
                        rd.info = format!("RDS={}", "single_snp").to_string().into_bytes();
                    }
                }
                if is_ase_snp {
                    if snp.phase_score < ase_ps_cutoff as f64 {
                        rd.filter = "LowQual".to_string().into_bytes();
                        if is_rna_edit {
                            rd.info = format!("RDS={}", "ase_snp,rna_editing").to_string().into_bytes();    // rna editing snp has bad phase
                        } else {
                            rd.info = format!("RDS={}", "ase_snp").to_string().into_bytes();
                        }
                    } else {
                        rd.filter = "PASS".to_string().into_bytes();
                        rd.info = format!("RDS={}", "ase_snp").to_string().into_bytes();
                    }
                }
                if !is_dense && !is_single_snp && !is_ase_snp {
                    if snp.variant_quality < min_qual_for_candidate as f64 {
                        rd.filter = "LowQual".to_string().into_bytes();
                        rd.info = "RDS=.".to_string().into_bytes();
                    } else if snp.phase_score < min_phase_score as f64 {
                        // bad phase, unconfindent phase
                        rd.filter = "LowQual".to_string().into_bytes();
                        if is_rna_edit {
                            rd.info = format!("RDS={}", "rna_editing").to_string().into_bytes();    // rna editing snp has bad phase
                        } else {
                            rd.info = "RDS=.".to_string().into_bytes();
                        }
                    } else {
                        // good phase, germline snp
                        rd.filter = "PASS".to_string().into_bytes();
                        rd.info = "RDS=.".to_string().into_bytes();
                    }
                }

                if is_dense || is_single_snp || is_unconfident_phased_snp || snp.phase_score == 0.0 || snp.haplotype == 0 {
                    rd.genotype = format!(
                        "{}:{}:{}:{:.2}",
                        "0/1",
                        snp.genotype_quality as i32,
                        snp.depth,
                        snp.allele_freqs[1]
                    );
                    rd.format = "GT:GQ:DP:AF".to_string().into_bytes();
                } else {
                    let mut af = 0.0;
                    if snp.alleles[0] == snp.reference {
                        af = snp.allele_freqs[1];
                    } else if snp.alleles[1] == snp.reference {
                        af = snp.allele_freqs[0];
                    } else {
                        panic!("Error: unexpected allele. ref: {}, alt1: {}, alt2: {}, {}:{}", snp.reference, snp.alleles[0], snp.alleles[1], String::from_utf8_lossy(&snp.chromosome), snp.pos);
                    }
                    if snp.phase_set != 0 {
                        if snp.haplotype == -1 {
                            rd.genotype = format!(
                                "{}:{}:{}:{}:{:.2}:{:.2}:{},{},{},{}",
                                "0|1",
                                snp.phase_set,
                                snp.genotype_quality as i32,
                                snp.depth,
                                af,
                                snp.phase_score,
                                snp.haplotype_expression[0],
                                snp.haplotype_expression[1],
                                snp.haplotype_expression[2],
                                snp.haplotype_expression[3]
                            );
                        } else if snp.haplotype == 1 {
                            rd.genotype = format!(
                                "{}:{}:{}:{}:{:.2}:{:.2}:{},{},{},{}",
                                "1|0",
                                snp.phase_set,
                                snp.genotype_quality as i32,
                                snp.depth,
                                af,
                                snp.phase_score,
                                snp.haplotype_expression[0],
                                snp.haplotype_expression[1],
                                snp.haplotype_expression[2],
                                snp.haplotype_expression[3]
                            );
                        }
                        rd.format = "GT:PS:GQ:DP:AF:PQ:AE".to_string().into_bytes();
                    } else {
                        if snp.haplotype == -1 {
                            rd.genotype = format!(
                                "{}:{}:{}:{:.2}:{:.2}:{},{},{},{}",
                                "0|1",
                                snp.genotype_quality as i32,
                                snp.depth,
                                af,
                                snp.phase_score,
                                snp.haplotype_expression[0],
                                snp.haplotype_expression[1],
                                snp.haplotype_expression[2],
                                snp.haplotype_expression[3]
                            );
                        } else if snp.haplotype == 1 {
                            rd.genotype = format!(
                                "{}:{}:{}:{:.2}:{:.2}:{},{},{},{}",
                                "1|0",
                                snp.genotype_quality as i32,
                                snp.depth,
                                af,
                                snp.phase_score,
                                snp.haplotype_expression[0],
                                snp.haplotype_expression[1],
                                snp.haplotype_expression[2],
                                snp.haplotype_expression[3]
                            );
                        }
                        rd.format = "GT:GQ:DP:AF:PQ:AE".to_string().into_bytes();
                    }
                }
                records.push(rd);
            }
            if snp.variant_type == 2 {
                let mut is_dense = false;

                if snp.filter == true {
                    is_dense = true;
                }

                let mut rd: VCFRecord = VCFRecord::default();
                rd.chromosome = snp.chromosome.clone();
                rd.position = snp.pos as u64 + 1; // position in vcf format is 1-based
                rd.id = vec!['.' as u8];
                rd.reference = vec![snp.reference as u8];
                rd.alternative = vec![vec![snp.alleles[0] as u8]];
                rd.qual = snp.variant_quality as i32;

                if is_dense {
                    rd.filter = "dn".to_string().into_bytes();
                    rd.info = format!("RDS={}", "dense_snp").to_string().into_bytes();
                }
                if !is_dense {
                    if snp.variant_quality < min_qual_for_candidate as f64 {
                        rd.filter = "LowQual".to_string().into_bytes();
                        rd.info = "RDS=.".to_string().into_bytes();
                    } else {
                        rd.filter = "PASS".to_string().into_bytes();
                        rd.info = "RDS=.".to_string().into_bytes();
                    }
                }
                rd.genotype = format!(
                    "{}:{}:{}:{:.2}",
                    "1/1",
                    snp.genotype_quality as i32,
                    snp.depth,
                    snp.allele_freqs[0]
                );
                rd.format = "GT:GQ:DP:AF".to_string().into_bytes();
                records.push(rd);
            }
            if snp.variant_type == 3 {
                let mut is_dense = false;

                if snp.filter == true {
                    is_dense = true;
                }

                let mut rd: VCFRecord = VCFRecord::default();
                rd.chromosome = snp.chromosome.clone();
                rd.position = snp.pos as u64 + 1; // position in vcf format is 1-based
                rd.id = vec!['.' as u8];
                rd.reference = vec![snp.reference as u8];
                rd.alternative = vec![vec![snp.alleles[0] as u8], vec![snp.alleles[1] as u8]];
                rd.qual = snp.variant_quality as i32;

                if is_dense {
                    rd.filter = "dn".to_string().into_bytes();
                    rd.info = format!("RDS={}", "dense_snp").to_string().into_bytes();
                }

                if !is_dense {
                    if snp.variant_quality < min_qual_for_candidate as f64 {
                        rd.filter = "LowQual".to_string().into_bytes();
                        rd.info = "RDS=.".to_string().into_bytes();
                    } else {
                        rd.filter = "PASS".to_string().into_bytes();
                        rd.info = "RDS=.".to_string().into_bytes();
                    }
                }
                rd.genotype = format!(
                    "{}:{}:{}:{:.2},{:.2}",
                    "1/2",
                    snp.genotype_quality as i32,
                    snp.depth,
                    snp.allele_freqs[0],
                    snp.allele_freqs[1]
                );
                rd.format = "GT:GQ:DP:AF".to_string().into_bytes();
                records.push(rd);
            }
            if snp.variant_type != 0 && snp.variant_type != 1 && snp.variant_type != 2 && snp.variant_type != 3 {
                panic!("Unknown variant type: {} {}:{}", snp.variant_type, String::from_utf8_lossy(&snp.chromosome), snp.pos);
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
            if snp.ase == true {
                continue;
            }
            if snp.filter == true && snp.rna_editing == false {
                // dense SNP
                let mut rd: VCFRecord = VCFRecord::default();
                rd.chromosome = snp.chromosome.clone();
                rd.position = snp.pos as u64 + 1; // position in vcf format is 1-based
                rd.reference = vec![snp.reference as u8];
                rd.id = vec!['.' as u8];
                if snp.variant_type == 1 {
                    if snp.alleles[0] != snp.reference && snp.alleles[1] == snp.reference {
                        rd.alternative = vec![vec![snp.alleles[0] as u8]];
                        rd.genotype = format!(
                            "{}:{}:{}:{:.2}",
                            "0/1", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[0]
                        );
                    } else if snp.alleles[1] != snp.reference && snp.alleles[0] == snp.reference {
                        rd.alternative = vec![vec![snp.alleles[1] as u8]];
                        rd.genotype = format!(
                            "{}:{}:{}:{:.2}",
                            "0/1", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[1]
                        );
                    } else {
                        rd.alternative = vec![vec![snp.alleles[0] as u8], vec![snp.alleles[1] as u8]];
                        rd.genotype = format!(
                            "{}:{}:{}:{:.2},{:.2}",
                            "1/2",
                            snp.genotype_quality as i32,
                            snp.depth,
                            snp.allele_freqs[0],
                            snp.allele_freqs[1]
                        );
                    }
                } else if snp.variant_type == 2 {
                    rd.alternative = vec![vec![snp.alleles[0] as u8]];
                    rd.genotype = format!(
                        "{}:{}:{}:{:.2}",
                        "1/1", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[0]
                    );
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
                rd.position = snp.pos as u64 + 1; // position in vcf format is 1-based
                rd.reference = vec![snp.reference as u8];
                rd.id = vec!['.' as u8];

                rd.alternative = vec![vec![snp.alleles[0] as u8]];
                rd.qual = snp.variant_quality as i32;
                rd.genotype = format!(
                    "{}:{}:{}:{:.2}",
                    "1/1", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[0]
                );
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
                rd.position = snp.pos as u64 + 1; // position in vcf format is 1-based
                rd.reference = vec![snp.reference as u8];
                rd.id = vec!['.' as u8];
                if snp.alleles[0] == snp.reference {
                    rd.alternative = vec![vec![snp.alleles[1] as u8]];
                    rd.qual = snp.variant_quality as i32;
                    rd.genotype = format!(
                        "{}:{}:{}:{:.2}",
                        "0/1", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[1]
                    );
                } else if snp.alleles[1] == snp.reference {
                    rd.alternative = vec![vec![snp.alleles[0] as u8]];
                    rd.qual = snp.variant_quality as i32;
                    rd.genotype = format!(
                        "{}:{}:{}:{:.2}",
                        "0/1", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[0]
                    );
                } else {
                    rd.alternative = vec![vec![snp.alleles[0] as u8], vec![snp.alleles[1] as u8]];
                    rd.qual = snp.variant_quality as i32;
                    rd.genotype = format!(
                        "{}:{}:{}:{:.2},{:.2}",
                        "1/2",
                        snp.genotype_quality as i32,
                        snp.depth,
                        snp.allele_freqs[0],
                        snp.allele_freqs[1]
                    );
                }
                if snp.variant_quality < min_qual as f64 {
                    rd.filter = "LowQual".to_string().into_bytes();
                } else {
                    rd.filter = "PASS".to_string().into_bytes();
                }
                rd.info = "RDS=.".to_string().into_bytes();
                rd.format = "GT:GQ:DP:AF".to_string().into_bytes();
                records.push(rd);
            } else if snp.variant_type == 0 {
                // homo ref. ase snp
                continue;
            } else {
                println!("Unknown variant type: {:?}:{:?}", snp.chromosome, snp.pos);
                continue;
            }
        }
        return records;
    }

    pub fn output_vcf2(&mut self, min_qual: u32) -> Vec<VCFRecord> {
        let mut records: Vec<VCFRecord> = Vec::new();

        // output heterozygous SNPs
        // assert_eq!(self.haplotype.len(), self.snps.len());
        for i in 0..self.candidate_snps.len() {
            let snp = &self.candidate_snps[i];
            if snp.ase == true {
                continue;
            }
            if snp.filter == true {
                // dense SNP
                let mut rd: VCFRecord = VCFRecord::default();
                rd.chromosome = snp.chromosome.clone();
                rd.position = snp.pos as u64 + 1; // position in vcf format is 1-based
                rd.reference = vec![snp.reference as u8];
                rd.id = vec!['.' as u8];
                if snp.variant_type == 1 {
                    if snp.alleles[0] != snp.reference && snp.alleles[1] == snp.reference {
                        rd.alternative = vec![vec![snp.alleles[0] as u8]];
                        rd.genotype = format!(
                            "{}:{}:{}:{:.2}",
                            "0/1", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[0]
                        );
                    } else if snp.alleles[1] != snp.reference && snp.alleles[0] == snp.reference {
                        rd.alternative = vec![vec![snp.alleles[1] as u8]];
                        rd.genotype = format!(
                            "{}:{}:{}:{:.2}",
                            "0/1", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[1]
                        );
                    } else {
                        rd.alternative = vec![vec![snp.alleles[0] as u8], vec![snp.alleles[1] as u8]];
                        rd.genotype = format!(
                            "{}:{}:{}:{:.2},{:.2}",
                            "1/2",
                            snp.genotype_quality as i32,
                            snp.depth,
                            snp.allele_freqs[0],
                            snp.allele_freqs[1]
                        );
                    }
                } else if snp.variant_type == 2 {
                    rd.alternative = vec![vec![snp.alleles[0] as u8]];
                    rd.genotype = format!(
                        "{}:{}:{}:{:.2}",
                        "1/1", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[0]
                    );
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
                rd.position = snp.pos as u64 + 1; // position in vcf format is 1-based
                rd.reference = vec![snp.reference as u8];
                rd.id = vec!['.' as u8];

                rd.alternative = vec![vec![snp.alleles[0] as u8]];
                rd.qual = snp.variant_quality as i32;
                rd.genotype = format!(
                    "{}:{}:{}:{:.2}",
                    "1/1", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[0]
                );
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
                rd.position = snp.pos as u64 + 1; // position in vcf format is 1-based
                rd.reference = vec![snp.reference as u8];
                rd.id = vec!['.' as u8];
                if snp.alleles[0] == snp.reference {
                    rd.alternative = vec![vec![snp.alleles[1] as u8]];
                    rd.qual = snp.variant_quality as i32;
                    rd.genotype = format!(
                        "{}:{}:{}:{:.2}",
                        "0/1", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[1]
                    );
                } else if snp.alleles[1] == snp.reference {
                    rd.alternative = vec![vec![snp.alleles[0] as u8]];
                    rd.qual = snp.variant_quality as i32;
                    rd.genotype = format!(
                        "{}:{}:{}:{:.2}",
                        "0/1", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[0]
                    );
                } else {
                    rd.alternative = vec![vec![snp.alleles[0] as u8], vec![snp.alleles[1] as u8]];
                    rd.qual = snp.variant_quality as i32;
                    rd.genotype = format!(
                        "{}:{}:{}:{:.2},{:.2}",
                        "1/2",
                        snp.genotype_quality as i32,
                        snp.depth,
                        snp.allele_freqs[0],
                        snp.allele_freqs[1]
                    );
                }
                if snp.variant_quality < min_qual as f64 {
                    rd.filter = "LowQual".to_string().into_bytes();
                } else {
                    rd.filter = "PASS".to_string().into_bytes();
                }
                rd.info = "RDS=.".to_string().into_bytes();
                rd.format = "GT:GQ:DP:AF".to_string().into_bytes();
                records.push(rd);
            } else if snp.variant_type == 0 {
                // homo ref. ase snp
                continue;
            } else {
                println!("Unknown variant type: {:?}:{:?}", snp.chromosome, snp.pos);
                continue;
            }
        }
        return records;
    }
}

fn exon_cluster(
    mut exons: Vec<Exon>,
    smallest_start: i64,
    largest_end: i64,
    min_sup: i32,
) -> HashMap<Exon, Vec<Exon>> {
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
    let mut ts = -1; // index in cover_vec
    let mut te = -1; // index in cover_vec
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

pub fn multithread_phase_haplotag(
    bam_file: String,
    ref_file: String,
    vcf_file: String,
    phased_bam_file: String,
    thread_size: usize,
    isolated_regions: Vec<Region>,
    exon_regions: HashMap<String, Vec<Interval<usize, u8>>>,
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
    min_linkers: u32,
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
    min_sup_haplotype_exon: u32,
    ase_allele_frac_cutoff: f32,
    ase_allele_cnt_cutoff: u32,
    ase_ps_count_cutoff: u32,
    ase_ps_cutoff: f32,
) {
    let pool = rayon::ThreadPoolBuilder::new().num_threads(thread_size).build().unwrap();
    let vcf_records_queue = Mutex::new(VecDeque::new());
    let read_haplotag1_queue = Mutex::new(VecDeque::new());
    let read_haplotag2_queue = Mutex::new(VecDeque::new());
    let read_haplotag_queue = Mutex::new(VecDeque::new());
    let read_phaseset_queue = Mutex::new(VecDeque::new());
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
            let mut exon_region_vec = Vec::new();
            if !reg.gene_id.is_none() {
                let gene_id_field = reg.gene_id.clone().unwrap();
                for gene_id in gene_id_field.split(",").collect::<Vec<&str>>() {
                    if exon_regions.contains_key(gene_id) {
                        exon_region_vec.extend(exon_regions.get(gene_id).unwrap().clone());
                    }
                }
                if exon_region_vec.len() == 0 {
                    // this region is done, no exon region covered
                    return;
                }
            }
            profile.init_with_pileup(
                &bam_file.as_str(),
                &reg,
                ref_seq,
                platform,
                min_mapq,
                min_baseq,
                min_read_length,
                min_depth,
                max_depth,
                distance_to_read_end,
                polya_tail_len,
            );
            let mut snpfrag = SNPFrag::default();
            snpfrag.region = reg.clone();
            snpfrag.min_linkers = min_linkers;
            snpfrag.get_candidate_snps(
                &profile,
                exon_region_vec,
                min_allele_freq,
                min_allele_freq_include_intron,
                min_qual_for_candidate,
                min_depth,
                max_depth,
                min_baseq,
                min_homozygous_freq,
                no_strand_bias,
                strand_bias_threshold,
                cover_strand_bias_threshold,
                distance_to_splicing_site,
                window_size,
                distance_to_read_end,
                diff_distance_to_read_end,
                diff_baseq,
                dense_win_size,
                min_dense_cnt,
                avg_dense_dist,
                ase_allele_frac_cutoff,
                ase_allele_cnt_cutoff,
            );
            // TODO: for very high depth region, down-sampling the reads
            // snpfrag.get_fragments(&bam_file, &reg);
            snpfrag.get_fragments_with_ase(&bam_file, &reg);
            if genotype_only {
                // without phasing
                let vcf_records = snpfrag.output_vcf2(min_qual_for_singlesnp_rnaedit);
                {
                    let mut queue = vcf_records_queue.lock().unwrap();
                    for rd in vcf_records.iter() {
                        queue.push_back(rd.clone());
                    }
                }
            } else {
                if snpfrag.hete_snps.len() > 0 {
                    unsafe {
                        snpfrag.init_haplotypes();
                    }
                    unsafe {
                        snpfrag.init_assignment();
                    }
                    snpfrag.phase(max_enum_snps, random_flip_fraction, max_iters);
                    let read_assignments = snpfrag.assign_reads(read_assignment_cutoff);
                    snpfrag.add_phase_score(min_allele_cnt, min_homozygous_freq, min_phase_score);
                    // TODO: assign phased fragments to somatic mutations
                    snpfrag.phase_ase_hete_snps(max_enum_snps, random_flip_fraction, max_iters);
                    // assign reads to haplotypes, filter reads having conflicted ase snps and heterozygous snps
                    let read_assignments_ase = snpfrag.assign_reads_ase(read_assignment_cutoff);
                    // snpfrag.rescue_ase_snps();
                    snpfrag.rescue_ase_snps_v2(ase_allele_cnt_cutoff, ase_ps_count_cutoff, ase_ps_cutoff);

                    // merge read_assignments and read_assignments_ase, read_assignments_ase first, then read_assignments
                    let mut merge_reads_assignments = read_assignments_ase.clone();
                    for (k, v) in read_assignments.iter() {
                        if !read_assignments_ase.contains_key(k) {
                            merge_reads_assignments.insert(k.clone(), v.clone());
                        }
                    }
                    let phase_sets = snpfrag.assign_phase_set(min_qual_for_candidate, min_phase_score, ase_ps_cutoff);

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
                                let hap1_consensus_exons = exon_cluster(
                                    hap1_exons.clone(),
                                    hap1_smallest_start,
                                    hap1_largest_end,
                                    0,
                                );
                                let hap2_consensus_exons = exon_cluster(
                                    hap2_exons.clone(),
                                    hap2_smallest_start,
                                    hap2_largest_end,
                                    0,
                                );
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
                                                haplotype_exons.push((
                                                    e.clone(),
                                                    counts.0,
                                                    counts.1,
                                                ));
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        // if !no_bam_output {
                        //     let mut queue = read_haplotag_queue.lock().unwrap();
                        //     for a in read_assignments.iter() {
                        //         queue.push_back((a.0.clone(), a.1.clone()));
                        //     }
                        // }

                        // output assignment both for ase snps and heterozygous snps
                        if !no_bam_output {
                            let mut queue = read_haplotag_queue.lock().unwrap();
                            for a in merge_reads_assignments.iter() {
                                queue.push_back((a.0.clone(), a.1.clone()));
                            }
                            let mut queue = read_phaseset_queue.lock().unwrap();
                            for a in phase_sets.iter() {
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
                // let vcf_records = snpfrag.phased_output_vcf(
                //     min_phase_score,
                //     min_homozygous_freq,
                //     output_phasing,
                //     min_qual_for_candidate,
                //     min_qual_for_singlesnp_rnaedit,
                // );
                let vcf_records = snpfrag.phased_output_vcf3(min_phase_score, ase_ps_cutoff, min_qual_for_candidate, min_qual_for_singlesnp_rnaedit);
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
    vf.write("##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase Set\">\n".as_bytes()).unwrap();
    vf.write("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n".as_bytes()).unwrap();
    vf.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n".as_bytes()).unwrap();
    vf.write("##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n".as_bytes()).unwrap();
    vf.write("##FORMAT=<ID=PQ,Number=1,Type=Float,Description=\"Phasing Quality\">\n".as_bytes()).unwrap();
    vf.write("##FORMAT=<ID=AE,Number=A,Type=Integer,Description=\"Haplotype expression of two alleles\">\n".as_bytes()).unwrap();
    vf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample\n".as_bytes()).unwrap();

    for rd in vcf_records_queue.lock().unwrap().iter() {
        if rd.alternative.len() == 1 {
            vf.write(
                format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                    std::str::from_utf8(&rd.chromosome).unwrap(),
                    rd.position,
                    std::str::from_utf8(&rd.id).unwrap(),
                    std::str::from_utf8(&rd.reference).unwrap(),
                    std::str::from_utf8(&rd.alternative[0]).unwrap(),
                    rd.qual,
                    std::str::from_utf8(&rd.filter).unwrap(),
                    std::str::from_utf8(&rd.info).unwrap(),
                    std::str::from_utf8(&rd.format).unwrap(),
                    rd.genotype
                ).as_bytes(),
            ).unwrap();
        } else if rd.alternative.len() == 2 {
            vf.write(
                format!(
                    "{}\t{}\t{}\t{}\t{},{}\t{}\t{}\t{}\t{}\t{}\n",
                    std::str::from_utf8(&rd.chromosome).unwrap(),
                    rd.position,
                    std::str::from_utf8(&rd.id).unwrap(),
                    std::str::from_utf8(&rd.reference).unwrap(),
                    std::str::from_utf8(&rd.alternative[0]).unwrap(),
                    std::str::from_utf8(&rd.alternative[1]).unwrap(),
                    rd.qual,
                    std::str::from_utf8(&rd.filter).unwrap(),
                    std::str::from_utf8(&rd.info).unwrap(),
                    std::str::from_utf8(&rd.format).unwrap(),
                    rd.genotype
                ).as_bytes(),
            ).unwrap();
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
        exon_writer.write(
            "#Chromosome\tExon start\tExon end\tExon state\tHap1 expression\tHap2 expression\n".as_bytes(),
        ).unwrap();
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
                exon_writer.write(
                    format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\n",
                        e.chr,
                        e.start + 1,
                        e.end,
                        e.state,
                        rd.1,
                        rd.2
                    ).as_bytes(),
                ).unwrap(); // 1-based, start inclusive, end inclusive
            }
        }
        drop(exon_writer);
    }

    if !no_bam_output {
        if !haplotype_bam_output {
            let mut read_assignments: HashMap<String, i32> = HashMap::new();
            for rd in read_haplotag_queue.lock().unwrap().iter() {
                if read_assignments.contains_key(&rd.0) {
                    read_assignments.remove(&rd.0);   // one read belongs to at least two regions
                } else {
                    read_assignments.insert(rd.0.clone(), rd.1);
                }
            }
            let mut read_phasesets: HashMap<String, u32> = HashMap::new();
            for rd in read_phaseset_queue.lock().unwrap().iter() {
                if read_phasesets.contains_key(&rd.0) {
                    read_phasesets.remove(&rd.0);   // one read belongs to at least two regions or two phase sets
                } else {
                    read_phasesets.insert(rd.0.clone(), rd.1);
                }
            }
            let mut bam_reader = bam::IndexedReader::from_path(&bam_file).unwrap();
            let header = bam::Header::from_template(&bam_reader.header());
            let mut bam_writer = bam::Writer::from_path(phased_bam_file, &header, Format::Bam).unwrap();
            bam_writer.set_threads(thread_size).unwrap();
            for region in isolated_regions.iter() {
                // TODO: duplicate reads in different regions
                bam_reader.fetch((region.chr.as_str(), region.start, region.end)).unwrap(); // set region
                for r in bam_reader.records() {
                    let mut record = r.unwrap();
                    if record.is_unmapped() || record.is_secondary() || record.is_supplementary() {
                        continue;
                    }
                    if record.reference_start() + 1 < region.start as i64 || record.reference_end() + 1 > region.end as i64 {
                        // reads beyond the region boundary will be ignored to provent duplicated reads
                        continue;
                    }
                    let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
                    if read_assignments.contains_key(&qname) {
                        let asg = read_assignments.get(&qname).unwrap();
                        if *asg != 0 {
                            let _ = record.push_aux(b"HP:i", Aux::I32(*asg));
                        }
                    }
                    if read_phasesets.contains_key(&qname) {
                        let ps = read_phasesets.get(&qname).unwrap();
                        let _ = record.push_aux(b"PS:i", Aux::U32(*ps));
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
            let mut hap1_bam_writer = bam::Writer::from_path(
                phased_bam_file.replace("phased", "hap1"),
                &header,
                Format::Bam,
            ).unwrap();
            let mut hap2_bam_writer = bam::Writer::from_path(
                phased_bam_file.replace("phased", "hap2"),
                &header,
                Format::Bam,
            ).unwrap();
            hap1_bam_writer.set_threads(thread_size).unwrap();
            hap2_bam_writer.set_threads(thread_size).unwrap();
            for region in isolated_regions.iter() {
                bam_reader.fetch((region.chr.as_str(), region.start, region.end)).unwrap(); // set region
                for r in bam_reader.records() {
                    let mut record = r.unwrap();
                    if record.is_unmapped() || record.is_secondary() || record.is_supplementary() {
                        continue;
                    }
                    if record.reference_start() + 1 < region.start as i64 || record.reference_end() + 1 > region.end as i64 {
                        // reads beyond the region boundary will be ignored to provent duplicated reads
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
