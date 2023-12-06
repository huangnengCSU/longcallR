use std::collections::{HashMap, HashSet, VecDeque};
use std::{cmp, fs};
use std::fs::File;
use std::hash::Hash;
use std::io::{BufRead, BufReader, Write};
use crate::bam_reader::Region;
use crate::profile::Profile;
use rust_htslib::{bam, bam::Read, bam::record::Record, bam::Format, bam::record::Aux};
use rust_htslib::htslib::{drand48, fai_load};
use std::sync::{mpsc, Arc, Mutex, Condvar};
use bio::bio_types::strand::ReqStrand::Forward;
use rayon::prelude::*;
use rand::seq::SliceRandom;
use log::error;
use crate::base_matrix::load_reference;
use crate::vcf::VCFRecord;


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
    // 1: heterozygous SNP, 2: homozygous SNP
    pub genotype_likelihood: [f64; 3],
    // 0th: homo var, 1st: hete var, 2nd: homo ref
    pub genotype_quality: f64,
    pub haplotype: i32,
    // delta: 1,-1: hap1, hap2 if phased, 0: unassigned
    pub phase_score: f64,
    pub snp_cover_fragments: Vec<usize>,
    // index of the fragment cover this SNP
    pub filter: bool,
    // filter out this SNP or not
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
    // haplotype of base on the alphabet {-1, 1, 0}, 1: base==alleles[0], -1: base==alleles[1], 0: not covered
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

impl SNPFrag {
    pub fn get_candidate_snps(&mut self, profile: &Profile, min_allele_freq: f32, min_allele_freq_include_intron: f32, min_coverage: u32, min_homozygous_freq: f32, cover_strand_bias_threshold: f32) {
        // get candidate SNPs, filtering with min_coverage, deletion_freq, min_allele_freq_include_intron, cover_strand_bias_threshold
        let pileup = &profile.freq_vec;
        let mut position = profile.region.start - 1;    // 0-based
        for bf in pileup.iter() {
            // println!("{:?}", bf);
            if bf.i {
                // insertion base
                continue;
            }

            // 1.filtering with depth
            let depth = bf.get_depth_exclude_intron_deletion();
            if depth < min_coverage {
                position += 1;
                continue;
            }

            let (allele1, allele1_cnt, allele2, allele2_cnt) = bf.get_two_major_alleles();

            // 2.filtering with depth, considering intron reads
            let depth_include_intron = bf.get_depth_include_intron();
            if (allele1_cnt as f32) / (depth_include_intron as f32) < min_allele_freq_include_intron {
                // maybe caused by erroneous intron alignment
                println!("allele freq include intron: {}, {}, {}, {}", position, allele1_cnt, allele2_cnt, depth_include_intron);
                position += 1;
                continue;
            }

            // 3.filtering with deletion frequency
            if bf.d > allele1_cnt {
                println!("indels: {}, {}, {}, {}", position, bf.d, allele1_cnt, allele2_cnt);
                position += 1;
                continue;
            }

            // 4.filtering snps only covered by one strand reads (may caused by intron alignment error)
            let total_cover_cnt = bf.forward_cnt + bf.backward_cnt;   // does not include intron reads
            if bf.forward_cnt as f32 / total_cover_cnt as f32 > cover_strand_bias_threshold || bf.backward_cnt as f32 / total_cover_cnt as f32 > cover_strand_bias_threshold {
                // strand bias
                println!("cover strand bias: {}, {}, {}, {}", position, bf.forward_cnt, bf.backward_cnt, total_cover_cnt);
                position += 1;
                continue;
            }


            // genotype likelihood
            let ploidy = 2;
            let mut likelihood = [1.0, 1.0, 1.0];
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
            } else {
                panic!("Error: unknown reference base");
            }

            for bq in identical_baseqs.iter() {
                let error_rate = 0.1_f64.powf((*bq as f64) / 10.0);
                likelihood[0] *= (ploidy) as f64 * error_rate;
                likelihood[1] *= (ploidy - 1) as f64 * error_rate + (1.0 - error_rate);
                likelihood[2] *= (ploidy - 2) as f64 * error_rate + 2.0 * (1.0 - error_rate);
            }

            for bq_vec in different_baseqs.iter() {
                for bq in bq_vec.iter() {
                    let error_rate = 0.1_f64.powf((*bq as f64) / 10.0);
                    likelihood[0] *= ploidy as f64 * (1.0 - error_rate);
                    likelihood[1] *= (ploidy - 1) as f64 * (1.0 - error_rate) + error_rate;
                    likelihood[2] *= (ploidy - 2) as f64 * (1.0 - error_rate) + 2.0 * error_rate;
                }
            }

            let num_reads = bf.a + bf.c + bf.g + bf.t;
            likelihood[0] = likelihood[0] / ((ploidy as f64).powf(num_reads as f64));
            likelihood[1] = likelihood[1] / ((ploidy as f64).powf(num_reads as f64));
            likelihood[2] = likelihood[2] / ((ploidy as f64).powf(num_reads as f64));

            // prior probability of genotype
            likelihood[0] *= background_prob[0];
            likelihood[1] *= background_prob[1];
            likelihood[2] *= background_prob[2];
            let likelihood_sum = likelihood[0] + likelihood[1] + likelihood[2];
            likelihood[0] = 10e-9_f64.max(likelihood[0] / likelihood_sum);
            likelihood[1] = 10e-9_f64.max(likelihood[1] / likelihood_sum);
            likelihood[2] = 10e-9_f64.max(likelihood[2] / likelihood_sum);

            // println!("likelihood: {}, {}, {}", likelihood[0], likelihood[1], likelihood[2]);

            if likelihood[0] > likelihood[1] && likelihood[0] > likelihood[2] {
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
                candidate_snp.variant_type = 2; // homozygous SNP
                candidate_snp.genotype_likelihood = likelihood.clone();
                candidate_snp.genotype_quality = -10.0_f64 * (1.0 - likelihood[0]).log10();
                self.candidate_snps.push(candidate_snp);
                self.homo_snps.push(self.candidate_snps.len() - 1);
            } else if likelihood[1] > likelihood[0] && likelihood[1] > likelihood[2] {
                // candidate heterozygous SNP
                let allele1_freq = (allele1_cnt as f32) / (depth as f32);
                let allele2_freq = (allele2_cnt as f32) / (depth as f32);
                if allele2_freq < min_allele_freq {
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
                candidate_snp.genotype_likelihood = likelihood.clone();
                candidate_snp.genotype_quality = -10.0_f64 * (1.0 - likelihood[1]).log10();
                self.candidate_snps.push(candidate_snp);
                self.hete_snps.push(self.candidate_snps.len() - 1);
            }


            // let allele1_freq = (allele1_cnt as f32) / (depth as f32);
            // let allele2_freq = (allele2_cnt as f32) / (depth as f32);
            // if allele1_freq >= min_homozygous_freq && allele1 != bf.ref_base {
            //     // candidate homozygous SNP
            //     let mut candidate_snp = CandidateSNP::default();
            //     candidate_snp.chromosome = profile.region.chr.clone().into_bytes();
            //     candidate_snp.pos = position as i64;
            //     candidate_snp.alleles = [allele1, allele2];
            //     candidate_snp.allele_freqs = [allele1_freq, allele2_freq];
            //     candidate_snp.reference = bf.ref_base;
            //     candidate_snp.depth = depth;
            //     candidate_snp.variant_type = 2; // homozygous SNP
            //     self.candidate_snps.push(candidate_snp);
            //     self.homo_snps.push(self.candidate_snps.len() - 1);
            // } else if allele2_freq > min_allele_freq {
            //     // candidate heterozgous SNP
            //     let mut candidate_snp = CandidateSNP::default();
            //     candidate_snp.chromosome = profile.region.chr.clone().into_bytes();
            //     candidate_snp.pos = position as i64;
            //     candidate_snp.alleles = [allele1, allele2];
            //     candidate_snp.allele_freqs = [allele1_freq, allele2_freq];
            //     candidate_snp.reference = bf.ref_base;
            //     candidate_snp.depth = depth;
            //     candidate_snp.variant_type = 1; // heterozygous SNP
            //     self.candidate_snps.push(candidate_snp);
            //     self.hete_snps.push(self.candidate_snps.len() - 1);
            // }
            position += 1;
        }
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
        for i in 0..self.fragments.len() {
            if drand48() < 0.5 {
                self.fragments[i].haplotag = -1;
            } else {
                self.fragments[i].haplotag = 1;
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
            for cg in cigar.iter() {
                if pos_on_ref > self.candidate_snps[*self.hete_snps.last().unwrap()].pos {
                    break;
                }
                if idx >= self.hete_snps.len() {
                    break;
                }
                match cg.char() as u8 {
                    b'S' | b'H' => {
                        continue;
                    }
                    b'M' => {
                        for _ in 0..cg.len() {
                            assert!(pos_on_ref <= snp_pos, "Error: pos_on_ref <= snp_pos");
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
                                if fragment.list.len() > 0 {
                                    for prev_frag_elem in fragment.list.iter() {
                                        if prev_frag_elem.p == 0 || frag_elem.p == 0 {
                                            continue;
                                        }
                                        let q1 = 0.1_f64.powf((prev_frag_elem.baseq as f64) / 10.0); // probability of error for prev_frag_elem
                                        let q2 = 0.1_f64.powf((frag_elem.baseq as f64) / 10.0);  // probability of error for frag_elem
                                        let epsilon = q1 * (1.0 - q2) + (1.0 - q1) * q2; // probability of sequencing error only occurs on start node or end node.
                                        let w = (prev_frag_elem.p as f64) * (frag_elem.p as f64) * f64::log10((1.0 - epsilon) / epsilon);
                                        // println!("q1:{}, q2:{}, epsilon:{}, w:{}", q1, q2, epsilon, w);
                                        if self.edges.contains_key(&[prev_frag_elem.snp_idx, frag_elem.snp_idx]) {
                                            let edge = self.edges.get_mut(&[prev_frag_elem.snp_idx, frag_elem.snp_idx]).unwrap();
                                            edge.frag_idxes.push(fragment.fragment_idx);
                                            edge.w += w;
                                        } else {
                                            let mut edge = Edge::default();
                                            edge.snp_idxes = [prev_frag_elem.snp_idx, frag_elem.snp_idx];
                                            edge.snp_poses = [prev_frag_elem.pos, frag_elem.pos];
                                            edge.w = w;
                                            edge.frag_idxes.push(fragment.fragment_idx);
                                            self.edges.insert(edge.snp_idxes, edge);
                                        }
                                    }
                                }
                                // println!("M: {:?}", frag_elem);
                                fragment.list.push(frag_elem);
                                idx += 1;
                                if idx >= self.hete_snps.len() {
                                    pos_on_query += 1;
                                    pos_on_ref += 1;
                                    break;
                                }
                                snp_pos = self.candidate_snps[self.hete_snps[idx]].pos;
                                alleles = self.candidate_snps[self.hete_snps[idx]].alleles.clone();
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
                            assert!(pos_on_ref <= snp_pos, "Error: pos_on_ref <= snp_pos");
                            if pos_on_ref == snp_pos {
                                let mut frag_elem = FragElem::default();
                                frag_elem.snp_idx = self.hete_snps[idx];
                                frag_elem.pos = pos_on_ref;
                                frag_elem.base = '-';
                                frag_elem.baseq = 0;
                                frag_elem.strand = strand;
                                frag_elem.p = 0;
                                // calculate the cost of links of two snps as hapcut (optioanl)
                                if fragment.list.len() > 0 {
                                    for prev_frag_elem in fragment.list.iter() {
                                        if prev_frag_elem.p == 0 || frag_elem.p == 0 {
                                            continue;
                                        }
                                        let q1 = 0.1_f64.powf((prev_frag_elem.baseq as f64) / 10.0); // probability of error for prev_frag_elem
                                        let q2 = 0.1_f64.powf((frag_elem.baseq as f64) / 10.0);  // probability of error for frag_elem
                                        let epsilon = q1 * (1.0 - q2) + (1.0 - q1) * q2; // probability of sequencing error only occurs on start node or end node.
                                        let w = (prev_frag_elem.p as f64) * (frag_elem.p as f64) * f64::log10((1.0 - epsilon) / epsilon);
                                        // println!("q1:{}, q2:{}, epsilon:{}, w:{}", q1, q2, epsilon, w);
                                        if self.edges.contains_key(&[prev_frag_elem.snp_idx, frag_elem.snp_idx]) {
                                            let edge = self.edges.get_mut(&[prev_frag_elem.snp_idx, frag_elem.snp_idx]).unwrap();
                                            edge.frag_idxes.push(fragment.fragment_idx);
                                            edge.w += w;
                                        } else {
                                            let mut edge = Edge::default();
                                            edge.snp_idxes = [prev_frag_elem.snp_idx, frag_elem.snp_idx];
                                            edge.snp_poses = [prev_frag_elem.pos, frag_elem.pos];
                                            edge.w = w;
                                            edge.frag_idxes.push(fragment.fragment_idx);
                                            self.edges.insert(edge.snp_idxes, edge);
                                        }
                                    }
                                }
                                // println!("D: {:?}", frag_elem);
                                fragment.list.push(frag_elem);
                                idx += 1;
                                if idx >= self.hete_snps.len() {
                                    pos_on_ref += 1;
                                    break;
                                }
                                snp_pos = self.candidate_snps[self.hete_snps[idx]].pos;
                                alleles = self.candidate_snps[self.hete_snps[idx]].alleles.clone();
                            }
                            pos_on_ref += 1;
                        }
                    }
                    b'N' => {
                        for _ in 0..cg.len() {
                            assert!(pos_on_ref <= snp_pos, "Error: pos_on_ref <= snp_pos");
                            if pos_on_ref == snp_pos {
                                let mut frag_elem = FragElem::default();
                                frag_elem.snp_idx = self.hete_snps[idx];
                                frag_elem.pos = pos_on_ref;
                                frag_elem.base = b'-' as char;
                                frag_elem.baseq = 0;
                                frag_elem.strand = strand;
                                frag_elem.p = 0;
                                // calculate the cost of links of two snps as hapcut (optioanl)
                                if fragment.list.len() > 0 {
                                    for prev_frag_elem in fragment.list.iter() {
                                        if prev_frag_elem.p == 0 || frag_elem.p == 0 {
                                            continue;
                                        }
                                        let q1 = 0.1_f64.powf((prev_frag_elem.baseq as f64) / 10.0); // probability of error for prev_frag_elem
                                        let q2 = 0.1_f64.powf((frag_elem.baseq as f64) / 10.0);  // probability of error for frag_elem
                                        let epsilon = q1 * (1.0 - q2) + (1.0 - q1) * q2; // probability of sequencing error only occurs on start node or end node.
                                        let w = (prev_frag_elem.p as f64) * (frag_elem.p as f64) * f64::log10((1.0 - epsilon) / epsilon);
                                        // println!("q1:{}, q2:{}, epsilon:{}, w:{}", q1, q2, epsilon, w);
                                        if self.edges.contains_key(&[prev_frag_elem.snp_idx, frag_elem.snp_idx]) {
                                            let edge = self.edges.get_mut(&[prev_frag_elem.snp_idx, frag_elem.snp_idx]).unwrap();
                                            edge.frag_idxes.push(fragment.fragment_idx);
                                            edge.w += w;
                                        } else {
                                            let mut edge = Edge::default();
                                            edge.snp_idxes = [prev_frag_elem.snp_idx, frag_elem.snp_idx];
                                            edge.snp_poses = [prev_frag_elem.pos, frag_elem.pos];
                                            edge.w = w;
                                            edge.frag_idxes.push(fragment.fragment_idx);
                                            self.edges.insert(edge.snp_idxes, edge);
                                        }
                                    }
                                }
                                // println!("N: {:?}", frag_elem);
                                fragment.list.push(frag_elem);
                                idx += 1;
                                if idx >= self.hete_snps.len() {
                                    pos_on_ref += 1;
                                    break;
                                }
                                snp_pos = self.candidate_snps[self.hete_snps[idx]].pos;
                                alleles = self.candidate_snps[self.hete_snps[idx]].alleles.clone();
                            }
                            pos_on_ref += 1;
                        }
                    }
                    _ => {
                        panic!("Error: unknown cigar operation");
                    }
                }
            }
            for fe in fragment.list.iter() {
                // record each snp cover by which fragments
                self.candidate_snps[fe.snp_idx].snp_cover_fragments.push(fragment.fragment_idx);
            }
            if fragment.list.len() > 0 {
                self.fragments.push(fragment);
            }
        }
    }

    pub fn filter_fp_snps(&mut self, strand_bias_threshold: f32, phase_score_cutoff: Option<f64>) {
        /*
        *? Filter false positive SNPs, which may be cause by strand bias or in a dense cluster of variants (100bp length has over 3 snps).
         */


        // 1. filter dense cluster of variants
        let mut homo_hete_snps: Vec<i64> = Vec::new();
        for i in self.hete_snps.iter() {
            homo_hete_snps.push(self.candidate_snps[*i].pos);
        }
        for i in self.homo_snps.iter() {
            homo_hete_snps.push(self.candidate_snps[*i].pos);
        }
        // sort homo_hete_snps in ascending order
        homo_hete_snps.sort();

        let mut filter_window: HashSet<i64> = HashSet::new();   // record the SNP position in a dense cluster of variants
        for i in 0..homo_hete_snps.len() {
            for j in i..homo_hete_snps.len() {
                // TODO: dense cluster of snps?
                if homo_hete_snps[j] - homo_hete_snps[i] > 300 {
                    if (j - 1) - i + 1 >= 9 {
                        for tk in i..j {
                            filter_window.insert(homo_hete_snps[tk]);
                        }
                    }
                    break;
                }
            }
        }

        println!("dense cluster filter: {:?}", filter_window);

        if filter_window.len() > 0 {
            // filter homo snps in a dense cluster of variants
            let mut filtered_homo_snps: Vec<usize> = Vec::new();
            let mut filtered_hete_snps: Vec<usize> = Vec::new();

            for i in self.homo_snps.iter() {
                if !filter_window.contains(&self.candidate_snps[*i].pos) {
                    filtered_homo_snps.push(*i);
                } else {
                    self.candidate_snps[*i].filter = true;
                }
            }
            self.homo_snps = filtered_homo_snps;

            for i in self.hete_snps.iter() {
                if !filter_window.contains(&self.candidate_snps[*i].pos) {
                    filtered_hete_snps.push(*i);
                } else {
                    self.candidate_snps[*i].filter = true;
                }
            }
            self.hete_snps = filtered_hete_snps;
        }

        // 2. filtering strand bias
        let mut filtered_hete_snps: Vec<usize> = Vec::new();
        for i in self.hete_snps.iter() {
            let snp = &mut self.candidate_snps[*i];
            if snp.filter == true { continue; }
            let mut variant_strand_cnt = [0, 0];    // variant strand [forward, reverse]
            for k in snp.snp_cover_fragments.iter() {
                for fe in self.fragments[*k].list.iter() {
                    if fe.snp_idx == *i {
                        if fe.p != 0 {
                            if fe.base != snp.reference {
                                if fe.strand == 0 {
                                    variant_strand_cnt[0] += 1;
                                } else {
                                    variant_strand_cnt[1] += 1;
                                }
                            }
                        }
                    }
                }
            }
            let total_variant_cnt = variant_strand_cnt[0] + variant_strand_cnt[1];
            if total_variant_cnt > 0 {
                // variant strand bias
                if variant_strand_cnt[0] as f32 / total_variant_cnt as f32 >= strand_bias_threshold || variant_strand_cnt[1] as f32 / total_variant_cnt as f32 >= strand_bias_threshold {
                    snp.filter = true;  // current snp is filtered because of strand bias
                    println!("strand bias: {}, {}, {}, {}", snp.pos, variant_strand_cnt[0], variant_strand_cnt[1], total_variant_cnt);
                    continue;
                }
            }
            filtered_hete_snps.push(*i);
        }
        if filtered_hete_snps.len() != self.hete_snps.len() {
            self.hete_snps = filtered_hete_snps;
        }
    }

    pub fn keep_reliable_snps_in_component(&mut self) {}

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
        if logp > pre_logp { rv = 1; } else if logp == pre_logp { rv = 0; } else { rv = -1; }
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
        let mut covered_fragments: HashSet<usize> = HashSet::new();
        let mut processed_snps: HashSet<usize> = HashSet::new();
        for i in self.hete_snps.iter() {
            processed_snps.insert(*i);
            for k in self.candidate_snps[*i].snp_cover_fragments.iter() {
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
                    // tmp_haplotag.push(sigma_k * (-1));  // flip sigma_k
                    tmp_haplotag.insert(*k, sigma_k * (-1));
                } else {
                    // tmp_haplotag.push(sigma_k);
                    tmp_haplotag.insert(*k, sigma_k);
                }
            }

            // assert!(SNPFrag::cal_overall_probability(&self, &processed_snps, &self.haplotype) >= SNPFrag::cal_overall_probability(&self, &self.haplotag, &self.haplotype));
            let check_val = SNPFrag::check_new_haplotag(&self, &processed_snps, &tmp_haplotag);
            assert!(check_val >= 0);

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
            assert!(check_val >= 0);

            if check_val > 0 {
                // new haplotype increases the probability P(sigma, delta)
                for (i, h) in tmp_haplotype.iter() {
                    self.candidate_snps[*i].haplotype = *h;
                }
            } else {
                phasing_increase = false;
            }
            num_iters += 1;
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
        for i in self.hete_snps.iter() {
            for k in self.candidate_snps[*i].snp_cover_fragments.iter() {
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
                assert!(tmp_hap.len() == self.hete_snps.len());
                for i in 0..self.hete_snps.len() {
                    self.candidate_snps[self.hete_snps[i]].haplotype = tmp_hap[i];
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
            for i in self.hete_snps.iter() {
                self.candidate_snps[*i].haplotype = best_haplotype[i];
            }
            for k in covered_fragments.iter() {
                self.fragments[*k].haplotag = best_haplotag[k];
            }
        }
    }

    pub fn add_phase_score(&mut self, min_allele_cnt: u32) {
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
                // filter SNPs with low allele count
                println!("Low allele count: {}, {}, {}", snp.pos, num_hap1, num_hap2);
                phase_score = 0.0;
            } else {
                if sigma.len() > 0 {
                    // println!("delta_i: {:?}", delta_i);
                    // println!("sigma: {:?}", sigma);
                    // println!("ps: {:?}", ps);
                    // phase_score = -10.0_f64 * SNPFrag::cal_inconsistent_percentage(delta_i, &sigma, &ps, &probs).log10();
                    phase_score = -10.0_f64 * (1.0 - SNPFrag::cal_log_delta_sigma(delta_i, &sigma, &ps, &probs)).log10();
                } else {
                    phase_score = 0.0;
                }
            }
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

    pub fn output_vcf2(&mut self, min_phase_score: f32, min_homozygous_freq: f32, output_phasing: bool) -> Vec<VCFRecord> {
        let mut records: Vec<VCFRecord> = Vec::new();

        // output heterozygous SNPs
        // assert_eq!(self.haplotype.len(), self.snps.len());
        for i in 0..self.candidate_snps.len() {
            let snp = &self.candidate_snps[i];
            let hp = snp.haplotype;
            if snp.filter == true {
                continue;
            }

            if snp.variant_type == 2 {
                let mut rd: VCFRecord = VCFRecord::default();
                rd.chromosome = snp.chromosome.clone();
                rd.position = snp.pos as u64 + 1;   // position in vcf format is 1-based
                rd.reference = vec![snp.reference as u8];
                rd.id = vec!['.' as u8];


                rd.alternative = vec![vec![snp.alleles[0] as u8]];
                rd.qual = cmp::max(1, (-10.0 * f64::log10(0.01_f64.max(1.0 - snp.allele_freqs[0] as f64))) as i32);
                rd.genotype = format!("{}:{}:{}:{:.2}", "1/1", rd.qual, snp.depth, snp.allele_freqs[1]);

                rd.filter = "PASS".to_string().into_bytes();
                rd.info = ".".to_string().into_bytes();
                rd.format = "GT:GQ:DP:AF".to_string().into_bytes();
                records.push(rd);
            } else if snp.variant_type == 1 {
                if snp.phase_score < min_phase_score as f64 {
                    if snp.allele_freqs[0] > min_homozygous_freq && snp.alleles[0] != snp.reference {
                        let mut rd: VCFRecord = VCFRecord::default();
                        rd.chromosome = snp.chromosome.clone();
                        rd.position = snp.pos as u64 + 1;   // position in vcf format is 1-based
                        rd.reference = vec![snp.reference as u8];
                        rd.id = vec!['.' as u8];


                        rd.alternative = vec![vec![snp.alleles[0] as u8]];
                        rd.qual = cmp::max(1, (-10.0 * f64::log10(0.01_f64.max(1.0 - snp.allele_freqs[0] as f64))) as i32);
                        rd.genotype = format!("{}:{}:{}:{:.2}", "1/1", rd.qual, snp.depth, snp.allele_freqs[1]);

                        rd.filter = "PASS".to_string().into_bytes();
                        rd.info = ".".to_string().into_bytes();
                        rd.format = "GT:GQ:DP:AF".to_string().into_bytes();
                        records.push(rd);
                    } else {
                        continue;
                    }
                } else {
                    let mut rd: VCFRecord = VCFRecord::default();
                    rd.chromosome = snp.chromosome.clone();
                    rd.position = snp.pos as u64 + 1;   // position in vcf format is 1-based
                    rd.reference = vec![snp.reference as u8];
                    rd.id = vec!['.' as u8];
                    if output_phasing {
                        if snp.alleles[0] == snp.reference {
                            rd.alternative = vec![vec![snp.alleles[1] as u8]];
                            rd.qual = cmp::max(1, (-10.0 * f64::log10(0.01_f64.max((0.5 - snp.allele_freqs[1] as f64).abs() / 0.5))) as i32);
                            if hp == -1 {
                                rd.genotype = format!("{}:{}:{}:{:.2}:{:.2}", "0|1", rd.qual, snp.depth, snp.allele_freqs[1], snp.phase_score);
                            } else {
                                rd.genotype = format!("{}:{}:{}:{:.2}:{:.2}", "1|0", rd.qual, snp.depth, snp.allele_freqs[1], snp.phase_score);
                            }
                        } else if snp.alleles[1] == snp.reference {
                            rd.alternative = vec![vec![snp.alleles[0] as u8]];
                            // rd.qual = -10.0 * f32::log10((0.5 - snp.allele_freqs[0]).abs() / 0.5);
                            rd.qual = cmp::max(1, (-10.0 * f64::log10(0.01_f64.max((0.5 - snp.allele_freqs[0] as f64).abs() / 0.5))) as i32);
                            if hp == -1 {
                                rd.genotype = format!("{}:{}:{}:{:.2}:{:.2}", "0|1", rd.qual, snp.depth, snp.allele_freqs[0], snp.phase_score);
                            } else {
                                rd.genotype = format!("{}:{}:{}:{:.2}:{:.2}", "1|0", rd.qual, snp.depth, snp.allele_freqs[0], snp.phase_score);
                            }
                        } else {
                            rd.alternative = vec![vec![snp.alleles[0] as u8], vec![snp.alleles[1] as u8]];
                            let q1 = cmp::max(1, (-10.0 * f64::log10(0.01_f64.max((0.5 - snp.allele_freqs[0] as f64).abs() / 0.5))) as i32);
                            let q2 = cmp::max(1, (-10.0 * f64::log10(0.01_f64.max((0.5 - snp.allele_freqs[1] as f64).abs() / 0.5))) as i32);
                            // if q1 > q2 { rd.qual = q2; } else { rd.qual = q1; }
                            rd.qual = cmp::min(q1, q2);
                            rd.genotype = format!("{}:{}:{}:{:.2},{:.2}:{:.2}", "1|2", rd.qual, snp.depth, snp.allele_freqs[0], snp.allele_freqs[1], snp.phase_score);
                        }

                        rd.filter = "PASS".to_string().into_bytes();
                        rd.info = ".".to_string().into_bytes();
                        rd.format = "GT:GQ:DP:AF:PQ".to_string().into_bytes();
                    } else {
                        if snp.alleles[0] == snp.reference {
                            rd.alternative = vec![vec![snp.alleles[1] as u8]];
                            rd.qual = cmp::max(1, (-10.0 * f64::log10(0.01_f64.max((0.5 - snp.allele_freqs[1] as f64).abs() / 0.5))) as i32);
                            rd.genotype = format!("{}:{}:{}:{:.2}", "0/1", rd.qual, snp.depth, snp.allele_freqs[1]);
                        } else if snp.alleles[1] == snp.reference {
                            rd.alternative = vec![vec![snp.alleles[0] as u8]];
                            rd.qual = cmp::max(1, (-10.0 * f64::log10(0.01_f64.max((0.5 - snp.allele_freqs[0] as f64).abs() / 0.5))) as i32);
                            rd.genotype = format!("{}:{}:{}:{:.2}", "0/1", rd.qual, snp.depth, snp.allele_freqs[0]);
                        } else {
                            rd.alternative = vec![vec![snp.alleles[0] as u8], vec![snp.alleles[1] as u8]];
                            let q1 = cmp::max(1, (-10.0 * f64::log10(0.01_f64.max((0.5 - snp.allele_freqs[0] as f64).abs() / 0.5))) as i32);
                            let q2 = cmp::max(1, (-10.0 * f64::log10(0.01_f64.max((0.5 - snp.allele_freqs[1] as f64).abs() / 0.5))) as i32);
                            // if q1 > q2 { rd.qual = q2; } else { rd.qual = q1; }
                            rd.qual = cmp::min(q1, q2);
                            rd.genotype = format!("{}:{}:{}:{:.2},{:.2}", "1/2", rd.qual, snp.depth, snp.allele_freqs[0], snp.allele_freqs[1]);
                        }

                        rd.filter = "PASS".to_string().into_bytes();
                        rd.info = ".".to_string().into_bytes();
                        rd.format = "GT:GQ:DP:AF".to_string().into_bytes();
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
}

fn parse_fai(fai_path: &str) -> Vec<(String, u32)> {
    let mut contig_lengths: Vec<(String, u32)> = Vec::new();
    let file = File::open(fai_path).unwrap();
    let reader = BufReader::new(file);
    for r in reader.lines() {
        let line = r.unwrap().clone();
        let parts: Vec<&str> = line.split('\t').collect();
        contig_lengths.push((parts[0].clone().to_string(), parts[1].clone().parse().unwrap()));
    }
    return contig_lengths;
}

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
                                  min_allele_freq: f32,
                                  min_allele_freq_include_intron: f32,
                                  min_allele_cnt: u32,
                                  strand_bias_threshold: f32,
                                  cover_strand_bias_threshold: f32,
                                  min_depth: u32,
                                  min_homozygous_freq: f32,
                                  min_phase_score: f32,
                                  max_enum_snps: usize,
                                  random_flip_fraction: f32,
                                  read_assignment_cutoff: f64,
                                  output_phasing: bool) {
    let pool = rayon::ThreadPoolBuilder::new().num_threads(thread_size - 1).build().unwrap();
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
            println!("Start {:?}", reg);
            let mut profile = Profile::default();
            profile.init_with_pileup(&bam_file.as_str(), &reg);
            profile.append_reference(&ref_seqs);
            let mut snpfrag = SNPFrag::default();
            snpfrag.get_candidate_snps(&profile, min_allele_freq, min_allele_freq_include_intron, min_depth, min_homozygous_freq, cover_strand_bias_threshold);
            for snp in snpfrag.candidate_snps.iter() {
                println!("snp: {:?}", snp);
            }
            snpfrag.get_fragments(&bam_file, &reg);
            snpfrag.filter_fp_snps(strand_bias_threshold, None);
            if snpfrag.hete_snps.len() > 0 {
                unsafe { snpfrag.init_haplotypes(); }
                unsafe { snpfrag.init_assignment(); }
                snpfrag.phase(max_enum_snps, random_flip_fraction);
                let read_assignments = snpfrag.assign_reads(read_assignment_cutoff);
                snpfrag.add_phase_score(min_allele_cnt);


                // // second round phase
                // {
                //     let mut second_round_hete_snps: Vec<usize> = Vec::new();
                //     for ti in 0..snpfrag.candidate_snps.len() {
                //         let snp = &snpfrag.candidate_snps[ti];
                //         if snp.filter == true || snp.variant_type != 1 {
                //             continue;
                //         }
                //
                //         if snp.phase_score >= min_phase_score as f64 {
                //             second_round_hete_snps.push(ti);
                //         }
                //     }
                //     snpfrag.hete_snps = second_round_hete_snps;
                //     if snpfrag.hete_snps.len() > 0 {
                //         unsafe { snpfrag.init_haplotypes(); }
                //         // remove the haplotag of each fragment
                //         for tk in 0..snpfrag.fragments.len() {
                //             snpfrag.fragments[tk].haplotag = 0;
                //         }
                //         unsafe { snpfrag.init_assignment(); }
                //         snpfrag.phase(max_enum_snps, random_flip_fraction);
                //         let read_assignments = snpfrag.assign_reads(read_assignment_cutoff);
                //         snpfrag.add_phase_score(min_allele_cnt);
                //     }
                // }


                {
                    let mut queue = read_haplotag_queue.lock().unwrap();
                    for a in read_assignments.iter() {
                        queue.push_back((a.0.clone(), a.1.clone()));
                    }
                }
            }
            // }
            let vcf_records = snpfrag.output_vcf2(min_phase_score, min_homozygous_freq, output_phasing);
            {
                let mut queue = vcf_records_queue.lock().unwrap();
                for rd in vcf_records.iter() {
                    queue.push_back(rd.clone());
                }
            }
        });
    });

    let mut vf = File::create(vcf_file).unwrap();
    vf.write("##fileformat=VCFv4.3\n".as_bytes()).unwrap();
    vf.write("##FILTER=<ID=PASS,Description=\"All filters passed\">\n".as_bytes()).unwrap();
    for ctglen in contig_lengths.iter() {
        let chromosome = ctglen.0.clone();
        let chromosome_len = ctglen.1.clone();
        vf.write(format!("##contig=<ID={},length={}>\n", chromosome, chromosome_len).as_bytes()).unwrap();
    }
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

    let mut read_assignments: HashMap<String, i32> = HashMap::new();
    for rd in read_haplotag_queue.lock().unwrap().iter() {
        read_assignments.insert(rd.0.clone(), rd.1.clone());
    }

    let mut bam_reader = bam::Reader::from_path(&bam_file).unwrap();
    let header = bam::Header::from_template(&bam_reader.header());
    let mut bam_writer = bam::Writer::from_path(phased_bam_file, &header, Format::Bam).unwrap();
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

