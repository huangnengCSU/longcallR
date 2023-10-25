use std::collections::HashMap;
use crate::bam_reader::Region;
use crate::profile::Profile;
use rust_htslib::{bam, bam::Read, bam::record::Record};
use rust_htslib::htslib::drand48;

#[derive(Debug, Clone, Default)]
pub struct CandidateSNP {
    pos: i64,
    // position on the reference, 0-based
    alleles: [char; 2],
    // major and minor alleles
    allele_freqs: [f32; 2],
    // major and minor allele frequencies
}

#[derive(Debug, Clone, Default)]
pub struct Edge {
    pub snp_idxes: [usize; 2],
    // index of candidate SNPs(SNPFrag.snps), start node and end node
    pub frag_idxes: Vec<usize>,
    // index of fragments(SNPFrag.fragments) cover this edge.
    pub w: f64,
    // weight of edge,  w_{ij}=\sum_{k}x_{ki}x_{kj}log\frac{1-\epsilon_{kij}}{\epsilon_{kij}}
}

#[derive(Debug, Clone, Default)]
struct FragElem {
    snp_idx: usize,
    // index of candidate SNPs(SNPFrag.snps)
    pos: i64,
    // position on the reference, 0-based
    base: char,
    // base pair
    baseq: u8,
    // base quality
    p: i32,
    // haplotype of base on the alphabet {-1, 1, 0}, 1: base==alleles[0], -1: base==alleles[1], 0: not covered
}

#[derive(Debug, Clone, Default)]
pub struct Fragment {
    fragment_idx: usize,
    // index of multiple fragments(SNPFrag.fragments)
    read_id: String,
    // read name
    list: Vec<FragElem>,
    // single fragment
}

#[derive(Debug, Clone, Default)]
pub struct SNPFrag {
    pub region: Region,
    pub snps: Vec<CandidateSNP>,
    // candidate SNPs
    pub fragments: Vec<Fragment>,
    // multiple fragments
    pub haplotype: Vec<i32>,
    // hap1. hap2 is bitwised complement of hap1
    pub edges: HashMap<[usize; 2], Edge>,
    // edges of the graph, key is [snp_idx of start_node, snp_idx of end_node]
}

impl SNPFrag {
    pub fn get_candidate_snps(&mut self, profile: &Profile, min_allele_freq: f32, min_coverage: u32) {
        // get candidate SNPs
        let pileup = &profile.freq_vec;
        let mut position = profile.region.start - 1;    // 0-based
        for bf in pileup.iter() {
            if bf.i {
                // insertion base
                continue;
            }
            // check current position is a candidate SNP
            let depth = bf.get_depth_exclude_intron();
            if depth < min_coverage {
                position += 1;
                continue;
            }
            let (allele1, allele1_cnt, allele2, allele2_cnt) = bf.get_two_major_alleles();
            let allele1_freq = (allele1_cnt as f32) / (depth as f32);
            let allele2_freq = (allele2_cnt as f32) / (depth as f32);
            if allele2_freq > min_allele_freq {
                // candidate SNP
                let mut candidate_snp = CandidateSNP::default();
                candidate_snp.pos = position as i64;
                // candidate_snp.ref_base = allele1 as u8;
                candidate_snp.alleles = [allele1, allele2];
                candidate_snp.allele_freqs = [allele1_freq, allele2_freq];
                self.snps.push(candidate_snp);
            }
            position += 1;
        }
    }

    unsafe fn init_haplotypes(&mut self) {
        // initialize haplotypes
        for i in 0..self.snps.len() {
            if drand48() < 0.5 {
                self.haplotype.push(-1);
            } else {
                self.haplotype.push(1);
            }
        }
    }

    pub fn get_fragments(&mut self, bam_path: &str, region: &Region) {
        let mut bam_reader: bam::IndexedReader = bam::IndexedReader::from_path(bam_path).unwrap();
        bam_reader.fetch((region.chr.as_str(), region.start, region.end)).unwrap();
        let mut record = Record::new();
        while let Some(result) = bam_reader.read(&mut record) {
            if result.is_err() {
                panic!("BAM parsing failed...");
            }
            if record.is_unmapped() || record.is_secondary() || record.is_supplementary() {
                continue;
            }
            let pos = record.pos(); // 0-based
            if pos > self.snps.last().unwrap().pos { continue; }
            let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
            // println!("{} {} {}", qname, pos, record.cigar());
            let cigar = record.cigar();
            let seq = record.seq().as_bytes();
            let mut pos_on_ref = pos;   // 0-based
            let mut pos_on_query = cigar.leading_softclips();   // 0-based
            let mut snp_offset = 0; // index of candidate SNPs
            let mut snp_pos = -1;   // pre-computed position of candidate SNPs
            let mut alleles;   // pre-computed alleles of candidate SNPs
            if pos <= self.snps.first().unwrap().pos {
                snp_pos = self.snps[snp_offset].pos;
                alleles = self.snps[snp_offset].alleles.clone();
            } else {
                // find the first SNP in the read
                while snp_offset < self.snps.len() {
                    if self.snps[snp_offset].pos >= pos {
                        break;
                    }
                    snp_offset += 1;
                }
                assert!(snp_offset < self.snps.len(), "Error: snp_offset < self.snps.len()");
                snp_pos = self.snps[snp_offset].pos;
                alleles = self.snps[snp_offset].alleles.clone();
            }

            let mut fragment = Fragment::default();
            fragment.read_id = qname.clone();
            fragment.fragment_idx = self.fragments.len();
            for cg in cigar.iter() {
                if pos_on_ref > self.snps.last().unwrap().pos {
                    break;
                }
                if snp_offset >= self.snps.len() {
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
                                frag_elem.snp_idx = snp_offset;
                                frag_elem.pos = pos_on_ref;
                                frag_elem.base = seq[pos_on_query as usize] as char;
                                frag_elem.baseq = record.qual()[pos_on_query as usize];
                                if frag_elem.base == alleles[0] {
                                    frag_elem.p = 1;    // reference allele
                                } else if frag_elem.base == alleles[1] {
                                    frag_elem.p = -1;   // alternate allele
                                } else {
                                    frag_elem.p = 0;    // not covered
                                }
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
                                            edge.w = w;
                                            edge.frag_idxes.push(fragment.fragment_idx);
                                            self.edges.insert(edge.snp_idxes, edge);
                                        }
                                    }
                                }
                                // println!("M: {:?}", frag_elem);
                                fragment.list.push(frag_elem);
                                snp_offset += 1;
                                if snp_offset >= self.snps.len() {
                                    pos_on_query += 1;
                                    pos_on_ref += 1;
                                    break;
                                }
                                snp_pos = self.snps[snp_offset].pos;
                                alleles = self.snps[snp_offset].alleles.clone();
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
                                frag_elem.snp_idx = snp_offset;
                                frag_elem.pos = pos_on_ref;
                                frag_elem.base = '-';
                                frag_elem.baseq = 0;
                                frag_elem.p = 0;
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
                                            edge.w = w;
                                            edge.frag_idxes.push(fragment.fragment_idx);
                                            self.edges.insert(edge.snp_idxes, edge);
                                        }
                                    }
                                }
                                // println!("D: {:?}", frag_elem);
                                fragment.list.push(frag_elem);
                                snp_offset += 1;
                                if snp_offset >= self.snps.len() {
                                    pos_on_ref += 1;
                                    break;
                                }
                                snp_pos = self.snps[snp_offset].pos;
                                alleles = self.snps[snp_offset].alleles.clone();
                            }
                            pos_on_ref += 1;
                        }
                    }
                    b'N' => {
                        for _ in 0..cg.len() {
                            assert!(pos_on_ref <= snp_pos, "Error: pos_on_ref <= snp_pos");
                            if pos_on_ref == snp_pos {
                                let mut frag_elem = FragElem::default();
                                frag_elem.snp_idx = snp_offset;
                                frag_elem.pos = pos_on_ref;
                                frag_elem.base = b'-' as char;
                                frag_elem.baseq = 0;
                                frag_elem.p = 0;
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
                                            edge.w = w;
                                            edge.frag_idxes.push(fragment.fragment_idx);
                                            self.edges.insert(edge.snp_idxes, edge);
                                        }
                                    }
                                }
                                // println!("N: {:?}", frag_elem);
                                fragment.list.push(frag_elem);
                                snp_offset += 1;
                                if snp_offset >= self.snps.len() {
                                    pos_on_ref += 1;
                                    break;
                                }
                                snp_pos = self.snps[snp_offset].pos;
                                alleles = self.snps[snp_offset].alleles.clone();
                            }
                            pos_on_ref += 1;
                        }
                    }
                    _ => {
                        panic!("Error: unknown cigar operation");
                    }
                }
            }
            if fragment.list.len() > 0 {
                self.fragments.push(fragment);
            }
        }
    }

    fn optimization_using_maxcut(&mut self) {
        // optimization using maxcut
    }
}

