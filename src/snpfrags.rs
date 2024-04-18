use std::collections::{HashMap, HashSet};

use petgraph::algo::kosaraju_scc;
use petgraph::graphmap::GraphMap;
use petgraph::Undirected;
use rand::Rng;
use rust_htslib::{bam, bam::Read, bam::record::Record};

use crate::snp::{CandidateSNP, Edge, Fragment};
use crate::somatic::calculate_prob_somatic;
use crate::util::Region;

#[derive(Debug, Clone, Default)]
pub struct SNPFrag {
    pub region: Region,
    pub candidate_snps: Vec<CandidateSNP>,
    // candidate SNPs
    pub homo_snps: Vec<usize>,
    // index of candidate homozygous SNPs
    pub edit_snps: Vec<usize>,
    // index of candidate rna editing SNPs
    pub low_frac_het_snps: Vec<usize>,
    // index of candidate low fraction het SNPs
    pub high_frac_het_snps: Vec<usize>,
    // index of candidate high fraction het SNPs
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

impl SNPFrag {
    pub unsafe fn init_haplotypes(&mut self) {
        // initialize haplotype of heterozygous snp
        let mut rng = rand::thread_rng();
        for i in self.high_frac_het_snps.iter() {
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

    pub fn get_somatic_haplotype_baseqs(&mut self, bam_path: &str, region: &Region, phased_fragments: &HashMap<String, i32>) {
        let mut bam_reader: bam::IndexedReader = bam::IndexedReader::from_path(bam_path).unwrap();
        bam_reader.fetch((region.chr.as_str(), region.start, region.end)).unwrap();
        let mut record = Record::new();
        // assert!(self.min_linkers >= 0, "Error: min_linkers <= 0");
        while let Some(result) = bam_reader.read(&mut record) {
            if self.somatic_snps.len() == 0 {
                continue;
            }
            if result.is_err() {
                panic!("BAM parsing failed...");
            }
            // TODO: filtering unmapped, secondary, supplementary reads?
            if record.is_unmapped() || record.is_secondary() || record.is_supplementary() {
                continue;
            }
            let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
            if !phased_fragments.contains_key(&qname) {
                continue;
            }
            let assignment = phased_fragments[&qname];
            let pos = record.pos(); // 0-based
            if pos > self.candidate_snps[*self.somatic_snps.last().unwrap()].pos {
                continue;
            }
            let cigar = record.cigar();
            let seq = record.seq().as_bytes();
            let mut pos_on_ref = pos; // 0-based
            let mut pos_on_query = cigar.leading_softclips(); // 0-based
            let mut idx = 0; // index in self.somatic_snps
            let mut snp_pos = -1; // pre-computed position of candidate SNPs
            let mut alleles; // pre-computed alleles of candidate SNPs
            if pos <= self.candidate_snps[*self.somatic_snps.first().unwrap()].pos {
                snp_pos = self.candidate_snps[self.somatic_snps[idx]].pos;
                alleles = self.candidate_snps[self.somatic_snps[idx]].alleles.clone();
            } else {
                // find the first SNP in the read
                while idx < self.somatic_snps.len() {
                    if self.candidate_snps[self.somatic_snps[idx]].pos >= pos {
                        break;
                    }
                    idx += 1;
                }
                assert!(
                    idx < self.somatic_snps.len(),
                    "Error: idx < self.candidate_snps.len()"
                );
                snp_pos = self.candidate_snps[self.somatic_snps[idx]].pos;
                alleles = self.candidate_snps[self.somatic_snps[idx]].alleles.clone();
            }

            for cg in cigar.iter() {
                match cg.char() as u8 {
                    b'S' | b'H' => {
                        continue;
                    }
                    b'M' | b'X' | b'=' => {
                        for _ in 0..cg.len() {
                            if pos_on_ref == snp_pos {
                                let somatic_cand = &mut self.candidate_snps[self.somatic_snps[idx]];
                                let base = seq[pos_on_query as usize] as char;
                                let baseq = record.qual()[pos_on_query as usize];
                                let [allele1, allele2] = somatic_cand.alleles.clone();
                                let ref_allele = somatic_cand.reference;
                                if allele1 == ref_allele || allele2 == ref_allele {
                                    if base == allele1 || base == allele2 {
                                        if base == ref_allele {
                                            if assignment == 1 {
                                                somatic_cand.hap_quals.hap1_ref_baseqs.push(baseq);
                                            } else {
                                                somatic_cand.hap_quals.hap2_ref_baseqs.push(baseq);
                                            }
                                        } else {
                                            if assignment == 1 {
                                                somatic_cand.hap_quals.hap1_alt_baseqs.push(baseq);
                                            } else {
                                                somatic_cand.hap_quals.hap2_alt_baseqs.push(baseq);
                                            }
                                        }
                                    }
                                }
                                idx += 1;
                                if idx < self.somatic_snps.len() {
                                    snp_pos = self.candidate_snps[self.somatic_snps[idx]].pos;
                                    alleles = self.candidate_snps[self.somatic_snps[idx]].alleles.clone();
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
                                idx += 1;
                                if idx < self.somatic_snps.len() {
                                    snp_pos = self.candidate_snps[self.somatic_snps[idx]].pos;
                                    alleles = self.candidate_snps[self.somatic_snps[idx]].alleles.clone();
                                }
                            }
                            pos_on_ref += 1;
                        }
                    }
                    b'N' => {
                        for _ in 0..cg.len() {
                            if pos_on_ref == snp_pos {
                                idx += 1;
                                if idx < self.somatic_snps.len() {
                                    snp_pos = self.candidate_snps[self.somatic_snps[idx]].pos;
                                    alleles = self.candidate_snps[self.somatic_snps[idx]].alleles.clone();
                                }
                            }
                            pos_on_ref += 1;
                        }
                    }
                    _ => {
                        panic!("Error: unknown cigar operation: {}", cg.char());
                    }
                }
            }
        }
        // for idx in self.somatic_snps.iter() {
        //     println!("somatic snp:{:?}\n, hap1_ref_baseqs:{:?}\n, hap1_alt_baseqs:{:?}\n, hap2_ref_baseqs:{:?}\n, hap2_alt_baseqs:{:?}\n",
        //              self.candidate_snps[*idx].pos,
        //              self.candidate_snps[*idx].hap_quals.hap1_ref_baseqs,
        //              self.candidate_snps[*idx].hap_quals.hap1_alt_baseqs,
        //              self.candidate_snps[*idx].hap_quals.hap2_ref_baseqs,
        //              self.candidate_snps[*idx].hap_quals.hap2_alt_baseqs);
        // }
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
                if fe.phase_site == false {
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
                if fe.phase_site == false {
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
                    if fe.phase_site == false {
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
                    if fe.phase_site == false {
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
            for i in self.high_frac_het_snps.iter() {
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
                            assert_ne!(fe.phase_site, false, "Error: phase for non-phase site.");
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

    // pub fn cross_optimize_ase_snps(&mut self) -> f64 {
    //     // Iteration:
    //     //     1. evaluate the assignment of each read based on the current SNP haplotype.
    //     //     2. evaluate the SNP haplotype based on the read assignment.
    //     // If P(sigma, delta) increase, repeat Iteration;
    //     // Else break;
    //
    //     let mut phasing_increase: bool = true;
    //     let mut haplotag_increase: bool = true;
    //     let mut num_iters = 0;
    //
    //     while phasing_increase | haplotag_increase {
    //         // optimize delta first with the pre-phased fragments in step of phasing hete snps
    //         // optimize delta
    //         let mut tmp_haplotype: HashMap<usize, i32> = HashMap::new();
    //         for i in self.ase_snps.iter() {
    //             // flip ase snps and keep heterozygous snps unchanged
    //             let delta_i = self.candidate_snps[*i].haplotype;
    //             let mut sigma: Vec<i32> = Vec::new();
    //             let mut ps: Vec<i32> = Vec::new();
    //             let mut probs: Vec<f64> = Vec::new();
    //             for k in self.candidate_snps[*i].snp_cover_fragments.iter() {
    //                 // k is fragment index
    //                 for fe in self.fragments[*k].list.iter() {
    //                     if fe.snp_idx == *i {
    //                         assert_ne!(fe.p, 0, "Error: phase for unexpected allele.");
    //                         ps.push(fe.p);
    //                         probs.push(fe.prob);
    //                         sigma.push(self.fragments[*k].haplotag);
    //                     }
    //                 }
    //             }
    //
    //             let q = SNPFrag::cal_delta_sigma_log(delta_i, &sigma, &ps, &probs);
    //             let qn = SNPFrag::cal_delta_sigma_log(delta_i * (-1), &sigma, &ps, &probs);
    //             // println!("optimize delta {} q:{:?}, qn:{:?}, delta: {}", i, q, qn, delta_i);
    //             if q < qn {
    //                 tmp_haplotype.insert(*i, delta_i * (-1));
    //             } else {
    //                 tmp_haplotype.insert(*i, delta_i);
    //             }
    //         }
    //         let check_val = SNPFrag::check_new_haplotype_ase(&self, &tmp_haplotype);
    //         assert!(check_val >= 0, "ckeck val bug: {:?}", self.candidate_snps);
    //         for (i, h) in tmp_haplotype.iter() {
    //             // when prob is equal, we still perform the flip to avoid bug of underflow
    //             self.candidate_snps[*i].haplotype = *h;
    //         }
    //         if check_val == 0 {
    //             phasing_increase = false;
    //         } else {
    //             phasing_increase = true;
    //             haplotag_increase = true;
    //         }
    //         self.check_local_optimal_configuration_ase(true, false); // maybe failed
    //
    //         // optimize sigma
    //         let mut tmp_haplotag: HashMap<usize, i32> = HashMap::new();
    //         let mut processed_snps = HashSet::new(); // some snps in self.hete_snps may be filtered by previous steps, record the snps that covered by the fragments
    //         for k in 0..self.fragments.len() {
    //             let mut sigma_k = self.fragments[k].haplotag;
    //             if sigma_k == 0 {
    //                 // the fragment may be unphased in the first step phasing with heterozygous snps, random assign sigma
    //                 let mut rng = rand::thread_rng();
    //                 let rg: f64 = rng.gen();
    //                 if rg < 0.5 {
    //                     self.fragments[k].haplotag = -1;
    //                     sigma_k = -1;
    //                 } else {
    //                     self.fragments[k].haplotag = 1;
    //                     sigma_k = 1;
    //                 }
    //             }
    //             let mut delta: Vec<i32> = Vec::new();
    //             let mut ps: Vec<i32> = Vec::new();
    //             let mut probs: Vec<f64> = Vec::new();
    //             for fe in self.fragments[k].list.iter() {
    //                 ps.push(fe.p);
    //                 probs.push(fe.prob);
    //                 delta.push(self.candidate_snps[fe.snp_idx].haplotype);
    //                 processed_snps.insert(fe.snp_idx);
    //             }
    //
    //             let q = SNPFrag::cal_sigma_delta_log(sigma_k, &delta, &ps, &probs);
    //             let qn = SNPFrag::cal_sigma_delta_log(sigma_k * (-1), &delta, &ps, &probs);
    //             // println!("optimize sigma {} q:{}, qn:{}, sigma: {}", k, q, qn, sigma_k);
    //
    //             if q < qn {
    //                 tmp_haplotag.insert(k, sigma_k * (-1));
    //             } else {
    //                 tmp_haplotag.insert(k, sigma_k);
    //             }
    //         }
    //
    //         // assert!(SNPFrag::cal_overall_probability(&self, &processed_snps, &self.haplotype) >= SNPFrag::cal_overall_probability(&self, &self.haplotag, &self.haplotype));
    //         let check_val = SNPFrag::check_new_haplotag_ase(&self, &tmp_haplotag);
    //         assert!(check_val >= 0, "ckeck val bug: {:?}", self.candidate_snps);
    //         for (k, h) in tmp_haplotag.iter() {
    //             // when prob is equal, we still perform the flip to avoid bug of underflow
    //             self.fragments[*k].haplotag = *h;
    //         }
    //         if check_val == 0 {
    //             haplotag_increase = false;
    //         } else {
    //             haplotag_increase = true;
    //             phasing_increase = true;
    //         }
    //         self.check_local_optimal_configuration_ase(false, true);
    //
    //
    //         num_iters += 1;
    //         if num_iters > 20 {
    //             break;
    //         }
    //     }
    //     // sigma reaches the optimal solution first and then delta reaches the optimal solution. After this, equal probability flip of delta may destroy the optimum of sigma again.
    //     // self.check_local_optimal_configuration(true, true);
    //     let prob = SNPFrag::cal_overall_probability_ase(&self);
    //     return prob;
    // }

    fn save_best_configuration(
        &self,
        best_haplotype: &mut HashMap<usize, i32>,
        best_haplotag: &mut HashMap<usize, i32>,
    ) {
        best_haplotype.clear();
        best_haplotag.clear();
        for i in self.high_frac_het_snps.iter() {
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
        for i in self.high_frac_het_snps.iter() {
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

        if self.high_frac_het_snps.len() <= max_enum_snps {
            // enumerate the haplotype, then optimize the assignment
            let mut haplotype_enum: Vec<Vec<i32>> = Vec::new();
            let init_hap: Vec<i32> = vec![1; self.high_frac_het_snps.len()];
            haplotype_enum.push(init_hap.clone());
            for ti in 0..self.high_frac_het_snps.len() {
                for tj in 0..haplotype_enum.len() {
                    let mut tmp_hap = haplotype_enum[tj].clone();
                    tmp_hap[ti] = tmp_hap[ti] * (-1);
                    haplotype_enum.push(tmp_hap);
                }
            }
            assert!(haplotype_enum.len() == 2_usize.pow(self.high_frac_het_snps.len() as u32), "Error: Not all combinations included");
            for hap in haplotype_enum.iter() {
                for i in 0..self.high_frac_het_snps.len() {
                    self.candidate_snps[self.high_frac_het_snps[i]].haplotype = hap[i];
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
                for i in self.high_frac_het_snps.iter() {
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
                        assert_eq!(tmp_hap.len(), self.high_frac_het_snps.len());
                        for i in 0..self.high_frac_het_snps.len() {
                            self.candidate_snps[self.high_frac_het_snps[i]].haplotype = tmp_hap[i];
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
                        for ti in 0..self.high_frac_het_snps.len() {
                            let rg: f64 = rng.gen();
                            if rg < random_flip_fraction as f64 {
                                self.candidate_snps[self.high_frac_het_snps[ti]].haplotype = self.candidate_snps[self.high_frac_het_snps[ti]].haplotype * (-1);
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

    // pub fn phase_ase_hete_snps(
    //     &mut self,
    //     max_enum_snps: usize,
    //     random_flip_fraction: f32,
    //     max_iters: i32,
    // ) {
    //     let mut largest_prob = f64::NEG_INFINITY;
    //     let mut best_haplotype: HashMap<usize, i32> = HashMap::new();
    //     let mut best_haplotag: HashMap<usize, i32> = HashMap::new();
    //     if self.ase_snps.len() <= max_enum_snps {
    //         // enumerate the haplotype, then optimize the assignment
    //         let mut haplotype_enum: Vec<Vec<i32>> = Vec::new();
    //         let init_hap: Vec<i32> = vec![1; self.ase_snps.len()];
    //         haplotype_enum.push(init_hap.clone());
    //         for ti in 0..self.ase_snps.len() {
    //             for tj in 0..haplotype_enum.len() {
    //                 let mut tmp_hap = haplotype_enum[tj].clone();
    //                 tmp_hap[ti] = tmp_hap[ti] * (-1);
    //                 haplotype_enum.push(tmp_hap);
    //             }
    //         }
    //         assert!(haplotype_enum.len() == 2_usize.pow(self.ase_snps.len() as u32), "Error: Not all combinations included");
    //         for hap in haplotype_enum.iter() {
    //             for i in 0..self.ase_snps.len() {
    //                 assert!(self.candidate_snps[self.ase_snps[i]].ase == true, "Error: only ase snps can be flipped.");
    //                 self.candidate_snps[self.ase_snps[i]].haplotype = hap[i];
    //             }
    //             let prob = self.cross_optimize_ase_snps();
    //             if prob > largest_prob {
    //                 largest_prob = prob;
    //                 self.save_best_configuration(&mut best_haplotype, &mut best_haplotag);
    //             }
    //         }
    //         self.load_best_configuration(&best_haplotype, &best_haplotag);
    //     } else {
    //         // optimize haplotype and read assignment alternatively
    //         let mut max_iter: i32 = max_iters;
    //         while max_iter >= 0 {
    //             // random initialization of haplotype for ase snps
    //             // if all fragments are not phased, randomly assign the haplotag, otherwise calculate haplotype based on pre-phased fragments and then calculate the haplotag of unphased fragments with newly calculated haplotype
    //             unsafe {
    //                 self.init_haplotypes_ase();
    //             }
    //             let mut phased_fragment_cnt = 0;
    //             for frag in self.fragments.iter() {
    //                 if frag.haplotag != 0 {
    //                     phased_fragment_cnt += 1;
    //                 }
    //             }
    //             if phased_fragment_cnt == 0 {
    //                 unsafe {
    //                     self.init_assignment();
    //                 }
    //             }
    //             let prob = self.cross_optimize_ase_snps();
    //             if prob > largest_prob {
    //                 largest_prob = prob;
    //                 self.save_best_configuration(&mut best_haplotype, &mut best_haplotag);
    //             }
    //             self.load_best_configuration(&best_haplotype, &best_haplotag);
    //
    //             // when initial setting has reached to local optimal, flip all the snps after a specific position to jump out local optimization
    //             let mut unflipped_haplotype: Vec<i32> = Vec::new();
    //             for i in self.ase_snps.iter() {
    //                 unflipped_haplotype.push(self.candidate_snps[*i].haplotype);
    //             }
    //             for ti in 0..unflipped_haplotype.len() {
    //                 let mut tmp_hap: Vec<i32> = Vec::new();
    //                 for tj in 0..unflipped_haplotype.len() {
    //                     if tj < ti {
    //                         tmp_hap.push(unflipped_haplotype[tj]);
    //                     } else {
    //                         tmp_hap.push(unflipped_haplotype[tj] * (-1));
    //                     }
    //                 }
    //                 // block flip
    //                 {
    //                     assert_eq!(tmp_hap.len(), self.ase_snps.len());
    //                     for i in 0..self.ase_snps.len() {
    //                         self.candidate_snps[self.ase_snps[i]].haplotype = tmp_hap[i];
    //                     }
    //                 }
    //                 let prob = self.cross_optimize_ase_snps();
    //                 if prob > largest_prob {
    //                     largest_prob = prob;
    //                     self.save_best_configuration(&mut best_haplotype, &mut best_haplotag);
    //                 }
    //                 self.load_best_configuration(&best_haplotype, &best_haplotag);
    //
    //                 // when current block flip has reached to local optimal, flip a fraction of snps and reads to jump out local optimization
    //                 {
    //                     let mut rng = rand::thread_rng();
    //                     for ti in 0..self.ase_snps.len() {
    //                         let rg: f64 = rng.gen();
    //                         if rg < random_flip_fraction as f64 {
    //                             self.candidate_snps[self.ase_snps[ti]].haplotype = self.candidate_snps[self.ase_snps[ti]].haplotype * (-1);
    //                         }
    //                     }
    //                     for tk in 0..self.fragments.len() {
    //                         let rg: f64 = rng.gen();
    //                         if rg < random_flip_fraction as f64 {
    //                             self.fragments[tk].haplotag = self.fragments[tk].haplotag * (-1);
    //                         }
    //                     }
    //                 }
    //                 let prob = self.cross_optimize_ase_snps();
    //                 if prob > largest_prob {
    //                     largest_prob = prob;
    //                     self.save_best_configuration(&mut best_haplotype, &mut best_haplotag);
    //                 }
    //                 self.load_best_configuration(&best_haplotype, &best_haplotag);
    //             }
    //             self.load_best_configuration(&best_haplotype, &best_haplotag);
    //             max_iter -= 1;
    //         }
    //         self.load_best_configuration(&best_haplotype, &best_haplotag);
    //     }
    // }

    pub fn assign_het_var_haplotype(
        &mut self,
        min_homozygous_freq: f32,
        min_phase_score: f32,
        somatic_allele_frac_cutoff: f32,
        somatic_allele_cnt_cutoff: u32,
    ) {
        // calculate phase score for each snp
        for ti in self.high_frac_het_snps.iter() {
            let snp = &mut self.candidate_snps[*ti];
            if !snp.for_phasing { continue; }
            if snp.snp_cover_fragments.len() == 0 {
                // no surranding haplotype links
                snp.single = true;
                continue;
            }
            let delta_i = snp.haplotype;
            let mut sigma: Vec<i32> = Vec::new();
            let mut ps: Vec<i32> = Vec::new();
            let mut probs: Vec<f64> = Vec::new();
            let mut assigns: Vec<i32> = Vec::new();
            let mut baseqs = Vec::new();
            let mut hap1_reads_num = 0;
            let mut hap2_reads_num = 0;
            for k in snp.snp_cover_fragments.iter() {
                if self.fragments[*k].assignment == 0 { continue; }
                if self.fragments[*k].num_hete_links < self.min_linkers { continue; }
                for fe in self.fragments[*k].list.iter() {
                    if fe.snp_idx == *ti {
                        if fe.base != '-' {
                            // ignore intron
                            if self.fragments[*k].assignment == 1 {
                                hap1_reads_num += 1;
                            } else if self.fragments[*k].assignment == 2 {
                                hap2_reads_num += 1;
                            }
                        }
                        assert_ne!(fe.phase_site, false, "Error: phase for non-phase site.");
                        assert_ne!(fe.p, 0, "Error: phase for unexpected allele.");
                        ps.push(fe.p);
                        probs.push(fe.prob);
                        sigma.push(self.fragments[*k].haplotag);
                        assigns.push(self.fragments[*k].assignment);
                        baseqs.push(fe.baseq);
                    }
                }
            }

            let debug_pos = 2341049;
            if snp.pos == debug_pos - 1 {
                println!("add phase score:");
                println!("ps:{:?}\nprobs:{:?}\nsigma:{:?}\nassigns:{:?}\n", ps, probs, sigma, assigns);
            }

            let mut phase_score = 0.0;
            if sigma.len() > 0 || hap1_reads_num >= 2 || hap2_reads_num >= 2 {
                // each haplotype should have at least 2 reads
                phase_score = -10.0_f64 * (1.0 - SNPFrag::cal_delta_sigma_log(delta_i, &sigma, &ps, &probs)).log10(); // calaulate assignment score
                if snp.pos == debug_pos - 1 {
                    println!("phase score:{}", phase_score);
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
                self.candidate_snps[*ti].germline = true;
                self.candidate_snps[*ti].haplotype_expression = haplotype_allele_expression;
                self.candidate_snps[*ti].phase_score = phase_score;
            }

            // TODO: het var with low phase score transfer to hom var
            // if phase_score < min_phase_score as f64 && snp.allele_freqs[0] > min_homozygous_freq && snp.alleles[0] != snp.reference && snp.filter == false {
            //     // transfer from heterozygous to homozygous
            //     snp.variant_type = 2;
            //     snp.ase = false;
            // }

            // TODO: het var with low phase score transfer to som var
            // if phase_score < min_phase_score as f64 {
            //     // record HapQuals for somatic mutation detection
            //     for k in 0..assigns.len() {
            //         if assigns[k] == 1 {
            //             if ps[k] == 1 {
            //                 snp.hap_quals.hap1_ref_baseqs.push(baseqs[k]);
            //             } else if ps[k] == -1 {
            //                 snp.hap_quals.hap1_alt_baseqs.push(baseqs[k]);
            //             }
            //         } else if assigns[k] == 2 {
            //             if ps[k] == 1 {
            //                 snp.hap_quals.hap2_ref_baseqs.push(baseqs[k]);
            //             } else if ps[k] == -1 {
            //                 snp.hap_quals.hap2_alt_baseqs.push(baseqs[k]);
            //             }
            //         }
            //     }
            //     let ref_allele_cnt = snp.hap_quals.hap1_ref_baseqs.len() + snp.hap_quals.hap2_ref_baseqs.len();
            //     let alt_allele_cnt = snp.hap_quals.hap1_alt_baseqs.len() + snp.hap_quals.hap2_alt_baseqs.len();
            //     if ref_allele_cnt + alt_allele_cnt > 0 && alt_allele_cnt as u32 >= somatic_allele_cnt_cutoff && alt_allele_cnt as f32 / (ref_allele_cnt + alt_allele_cnt) as f32 >= somatic_allele_frac_cutoff {
            //         // calculate somatic mutation probability
            //         let (hap1_allele_class, hap2_allele_class) = calculate_prob_somatic(&snp.hap_quals.hap1_ref_baseqs, &snp.hap_quals.hap1_alt_baseqs, &snp.hap_quals.hap2_ref_baseqs, &snp.hap_quals.hap2_alt_baseqs, 0.3);
            //         if hap1_allele_class.allcls == 0 && hap2_allele_class.allcls == 2 {
            //             let somatic_score = -10.0_f64 * (1.0 - hap2_allele_class.prob).log10();
            //             snp.cand_somatic = true;
            //             snp.somatic = true;
            //             snp.variant_type = 1;
            //             snp.somatic_score = somatic_score;
            //             snp.phase_score = 0.0;
            //             // println!("somatic snp:{}, score: {}", snp.pos, somatic_score);
            //             // println!("{:?},{:?}", hap1_allele_class, hap2_allele_class);
            //             // println!("hap1_ref_baseqs:{:?}\nhap1_alt_baseqs:{:?}\nhap2_ref_baseqs:{:?}\nhap2_alt_baseqs:{:?}", snp.hap_quals.hap1_ref_baseqs, snp.hap_quals.hap1_alt_baseqs, snp.hap_quals.hap2_ref_baseqs, snp.hap_quals.hap2_alt_baseqs);
            //         } else if hap1_allele_class.allcls == 2 && hap2_allele_class.allcls == 0 {
            //             let somatic_score = -10.0_f64 * (1.0 - hap1_allele_class.prob).log10();
            //             snp.cand_somatic = true;
            //             snp.somatic = true;
            //             snp.variant_type = 1;
            //             snp.somatic_score = somatic_score;
            //             snp.phase_score = 0.0;
            //             // println!("somatic snp:{}, score: {}", snp.pos, somatic_score);
            //             // println!("{:?},{:?}", hap1_allele_class, hap2_allele_class);
            //             // println!("hap1_ref_baseqs:{:?}\nhap1_alt_baseqs:{:?}\nhap2_ref_baseqs:{:?}\nhap2_alt_baseqs:{:?}", snp.hap_quals.hap1_ref_baseqs, snp.hap_quals.hap1_alt_baseqs, snp.hap_quals.hap2_ref_baseqs, snp.hap_quals.hap2_alt_baseqs);
            //         }
            //     }
            // }
        }
    }

    pub fn eval_low_frac_het_var_phase(&mut self, min_phase_score: f32,
                                       somatic_allele_frac_cutoff: f32,
                                       somatic_allele_cnt_cutoff: u32) {
        for ti in self.low_frac_het_snps.iter() {
            let snp = &mut self.candidate_snps[*ti];
            if snp.snp_cover_fragments.len() == 0 {
                // no surranding haplotype links
                snp.single = true;
                continue;
            }
            let mut sigma: Vec<i32> = Vec::new();
            let mut ps: Vec<i32> = Vec::new();
            let mut probs: Vec<f64> = Vec::new();
            let mut assigns: Vec<i32> = Vec::new();
            let mut baseqs = Vec::new();
            let mut hap1_reads_num = 0;
            let mut hap2_reads_num = 0;
            for k in snp.snp_cover_fragments.iter() {
                if self.fragments[*k].assignment == 0 { continue; }
                if self.fragments[*k].num_hete_links < self.min_linkers { continue; }
                for fe in self.fragments[*k].list.iter() {
                    if fe.snp_idx == *ti {
                        if fe.base != '-' {
                            // ignore intron
                            if self.fragments[*k].assignment == 1 {
                                hap1_reads_num += 1;
                            } else if self.fragments[*k].assignment == 2 {
                                hap2_reads_num += 1;
                            }
                        }
                        assert_ne!(fe.p, 0, "Error: phase for unexpected allele.");
                        ps.push(fe.p);
                        probs.push(fe.prob);
                        sigma.push(self.fragments[*k].haplotag);
                        assigns.push(self.fragments[*k].assignment);
                        baseqs.push(fe.baseq);
                    }
                }
            }

            let debug_pos = 2341049;
            if snp.pos == debug_pos - 1 {
                println!("eval_low_frac_het_var_phase:");
                println!("ps:{:?}\nprobs:{:?}\nsigma:{:?}\nassigns:{:?}\n", ps, probs, sigma, assigns);
            }

            let mut phase_score1 = 0.0;
            let mut phase_score2 = 0.0;
            let mut phase_score = 0.0;
            if sigma.len() == 0 || hap1_reads_num < 2 || hap2_reads_num < 2 {
                // each haplotype should have at least 2 reads
                continue;
            } else {
                phase_score1 = -10.0_f64 * (1.0 - SNPFrag::cal_delta_sigma_log(1, &sigma, &ps, &probs)).log10(); // calaulate assignment score
                phase_score2 = -10.0_f64 * (1.0 - SNPFrag::cal_delta_sigma_log(-1, &sigma, &ps, &probs)).log10(); // calaulate assignment score
            }

            if phase_score1.max(phase_score2) >= min_phase_score as f64 {
                phase_score = phase_score1.max(phase_score2);
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
                self.candidate_snps[*ti].germline = true;
                self.candidate_snps[*ti].haplotype = if phase_score1 >= phase_score2 { 1 } else { -1 };
                self.candidate_snps[*ti].haplotype_expression = haplotype_allele_expression;
                self.candidate_snps[*ti].phase_score = phase_score;
            }

            // TODO: het var with low phase score transfer to som var
            // if phase_score < min_phase_score as f64 {}
        }
    }
    pub fn eval_rna_edit_var_phase(&mut self, min_phase_score: f32) {
        for ti in self.edit_snps.iter() {
            let snp = &mut self.candidate_snps[*ti];
            if snp.snp_cover_fragments.len() == 0 {
                // no surranding haplotype links
                snp.single = true;
                continue;
            }
            let debug_pos = 2341049;
            if snp.pos == debug_pos - 1 {
                println!("eval_rna_edit_var_phase:");
            }
            let mut sigma: Vec<i32> = Vec::new();
            let mut ps: Vec<i32> = Vec::new();
            let mut probs: Vec<f64> = Vec::new();
            let mut assigns: Vec<i32> = Vec::new();
            let mut baseqs = Vec::new();
            let mut hap1_reads_num = 0;
            let mut hap2_reads_num = 0;
            for k in snp.snp_cover_fragments.iter() {
                if self.fragments[*k].assignment == 0 { continue; }
                if self.fragments[*k].num_hete_links < self.min_linkers { continue; }
                for fe in self.fragments[*k].list.iter() {
                    if fe.snp_idx == *ti {
                        if fe.base != '-' {
                            // ignore intron
                            if self.fragments[*k].assignment == 1 {
                                hap1_reads_num += 1;
                            } else if self.fragments[*k].assignment == 2 {
                                hap2_reads_num += 1;
                            }
                        }
                        assert_ne!(fe.p, 0, "Error: phase for unexpected allele.");
                        ps.push(fe.p);
                        probs.push(fe.prob);
                        sigma.push(self.fragments[*k].haplotag);
                        assigns.push(self.fragments[*k].assignment);
                        baseqs.push(fe.baseq);
                    }
                }
            }

            let debug_pos = 2341049;
            if snp.pos == debug_pos - 1 {
                println!("eval_rna_edit_var_phase:");
                println!("ps:{:?}\nprobs:{:?}\nsigma:{:?}\nassigns:{:?}\n", ps, probs, sigma, assigns);
            }

            let mut phase_score1 = 0.0;
            let mut phase_score2 = 0.0;
            let mut phase_score = 0.0;
            if sigma.len() == 0 || hap1_reads_num < 2 || hap2_reads_num < 2 {
                // each haplotype should have at least 2 reads
                continue;
            } else {
                phase_score1 = -10.0_f64 * (1.0 - SNPFrag::cal_delta_sigma_log(1, &sigma, &ps, &probs)).log10(); // calaulate assignment score
                phase_score2 = -10.0_f64 * (1.0 - SNPFrag::cal_delta_sigma_log(-1, &sigma, &ps, &probs)).log10(); // calaulate assignment score
            }

            if phase_score1.max(phase_score2) >= min_phase_score as f64 {
                phase_score = phase_score1.max(phase_score2);
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
                self.candidate_snps[*ti].germline = true;
                self.candidate_snps[*ti].rna_editing = false;
                self.candidate_snps[*ti].variant_type = 1;
                self.candidate_snps[*ti].haplotype = if phase_score1 >= phase_score2 { 1 } else { -1 };
                self.candidate_snps[*ti].haplotype_expression = haplotype_allele_expression;
                self.candidate_snps[*ti].phase_score = phase_score;
            }
        }
    }
    pub fn eval_som_var_phase(&mut self) {}
    pub fn eval_hom_var_phase(&mut self, min_phase_score: f32) {
        for ti in self.homo_snps.iter() {
            let snp = &mut self.candidate_snps[*ti];
            if snp.snp_cover_fragments.len() == 0 {
                // no surranding haplotype links
                snp.single = true;
                continue;
            }
            let mut sigma: Vec<i32> = Vec::new();
            let mut ps: Vec<i32> = Vec::new();
            let mut probs: Vec<f64> = Vec::new();
            let mut assigns: Vec<i32> = Vec::new();
            let mut baseqs = Vec::new();
            let mut hap1_reads_num = 0;
            let mut hap2_reads_num = 0;
            for k in snp.snp_cover_fragments.iter() {
                if self.fragments[*k].assignment == 0 { continue; }
                if self.fragments[*k].num_hete_links < self.min_linkers { continue; }
                for fe in self.fragments[*k].list.iter() {
                    if fe.snp_idx == *ti {
                        if fe.base != '-' {
                            // ignore intron
                            if self.fragments[*k].assignment == 1 {
                                hap1_reads_num += 1;
                            } else if self.fragments[*k].assignment == 2 {
                                hap2_reads_num += 1;
                            }
                        }
                        assert_ne!(fe.p, 0, "Error: phase for unexpected allele.");
                        ps.push(fe.p);
                        probs.push(fe.prob);
                        sigma.push(self.fragments[*k].haplotag);
                        assigns.push(self.fragments[*k].assignment);
                        baseqs.push(fe.baseq);
                    }
                }
            }

            let debug_pos = 2341049;
            if snp.pos == debug_pos - 1 {
                println!("eval_hom_var_phase:");
                println!("ps:{:?}\nprobs:{:?}\nsigma:{:?}\nassigns:{:?}\n", ps, probs, sigma, assigns);
            }

            let mut phase_score1 = 0.0;
            let mut phase_score2 = 0.0;
            let mut phase_score = 0.0;
            if sigma.len() == 0 || hap1_reads_num < 2 || hap2_reads_num < 2 {
                // each haplotype should have at least 2 reads
                continue;
            } else {
                phase_score1 = -10.0_f64 * (1.0 - SNPFrag::cal_delta_sigma_log(1, &sigma, &ps, &probs)).log10(); // calaulate assignment score
                phase_score2 = -10.0_f64 * (1.0 - SNPFrag::cal_delta_sigma_log(-1, &sigma, &ps, &probs)).log10(); // calaulate assignment score
            }

            if phase_score1.max(phase_score2) >= 3.0 * min_phase_score as f64 {
                // correct the genotype of this hom_var site, this site should be het_var
                phase_score = phase_score1.max(phase_score2);
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
                self.candidate_snps[*ti].germline = true;
                self.candidate_snps[*ti].hom_var = false;
                self.candidate_snps[*ti].variant_type = 1;
                self.candidate_snps[*ti].haplotype = if phase_score1 >= phase_score2 { 1 } else { -1 };
                self.candidate_snps[*ti].haplotype_expression = haplotype_allele_expression;
                self.candidate_snps[*ti].phase_score = phase_score;
            }
        }
    }

    pub fn assign_reads_haplotype(&mut self, read_assignment_cutoff: f64) -> HashMap<String, i32> {
        let mut read_assignments: HashMap<String, i32> = HashMap::new();
        for k in 0..self.fragments.len() {
            let sigma_k = self.fragments[k].haplotag;
            let mut delta: Vec<i32> = Vec::new();
            let mut ps: Vec<i32> = Vec::new();
            let mut probs: Vec<f64> = Vec::new();
            if self.fragments[k].read_id == "SRR18130587.2988112".to_string() {
                println!("{:?}", self.fragments[k]);
            }
            for fe in self.fragments[k].list.iter() {
                if fe.phase_site == false { continue; }
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

                if (q - qn).abs() >= read_assignment_cutoff {
                    if q >= qn {
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
                    // panic!("Error: unexpected condition.");
                    self.fragments[k].assignment = 0;
                    self.fragments[k].assignment_score = 0.0;
                    read_assignments.insert(self.fragments[k].read_id.clone(), 0);
                }
            }
        }
        return read_assignments;
    }

    // pub fn assign_reads_ase(&mut self, read_assignment_cutoff: f64) -> HashMap<String, i32> {
    //     let mut read_assignments: HashMap<String, i32> = HashMap::new();
    //     for k in 0..self.fragments.len() {
    //         let sigma_k = self.fragments[k].haplotag;
    //         let mut delta_hete: Vec<i32> = Vec::new();
    //         let mut ps_hete: Vec<i32> = Vec::new();
    //         let mut probs_hete: Vec<f64> = Vec::new();
    //         let mut delta_ase: Vec<i32> = Vec::new();
    //         let mut ps_ase: Vec<i32> = Vec::new();
    //         let mut probs_ase: Vec<f64> = Vec::new();
    //
    //         for fe in self.fragments[k].list.iter() {
    //             assert_ne!(fe.p, 0, "Error: phase for unexpected allele.");
    //             if fe.ase_snp == true {
    //                 ps_ase.push(fe.p);
    //                 probs_ase.push(fe.prob);
    //                 delta_ase.push(self.candidate_snps[fe.snp_idx].haplotype);
    //             } else {
    //                 ps_hete.push(fe.p);
    //                 probs_hete.push(fe.prob);
    //                 delta_hete.push(self.candidate_snps[fe.snp_idx].haplotype);
    //             }
    //         }
    //         if self.fragments[k].read_id == "m84036_230422_223801_s1/42010068/ccs/6216_8325".to_string() {
    //             println!("S1");
    //             println!("sigma_k:{:?}", self.fragments[k].haplotag);
    //             println!("delta_hete:{:?}\nps_hete:{:?}\nprobs_hete:{:?}", delta_hete, ps_hete, probs_hete);
    //             println!("delta_ase:{:?}\nps_ase:{:?}\nprobs_ase:{:?}", delta_ase, ps_ase, probs_ase);
    //             println!("ase snps: {:?}", self.ase_snps);
    //         }
    //         if sigma_k == 0 {
    //             // unasigned haplotag, cluster the read into unknown group
    //             self.fragments[k].assignment = 0;
    //             self.fragments[k].assignment_score = 0.0;
    //             read_assignments.insert(self.fragments[k].read_id.clone(), 0);
    //         } else {
    //             let mut q_hete = 0.0;
    //             let mut qn_hete = 0.0;
    //             let mut q_ase = 0.0;
    //             let mut qn_ase = 0.0;
    //             if delta_hete.len() > 0 {
    //                 q_hete = SNPFrag::cal_sigma_delta_log(sigma_k, &delta_hete, &ps_hete, &probs_hete);
    //                 qn_hete = SNPFrag::cal_sigma_delta_log(
    //                     sigma_k * (-1),
    //                     &delta_hete,
    //                     &ps_hete,
    //                     &probs_hete,
    //                 );
    //             }
    //             if delta_ase.len() > 0 {
    //                 q_ase = SNPFrag::cal_sigma_delta_log(sigma_k, &delta_ase, &ps_ase, &probs_ase);
    //                 qn_ase = SNPFrag::cal_sigma_delta_log(
    //                     sigma_k * (-1),
    //                     &delta_ase,
    //                     &ps_ase,
    //                     &probs_ase,
    //                 );
    //             }
    //
    //             if (q_hete - qn_hete).abs() >= read_assignment_cutoff && (q_ase - qn_ase).abs() >= read_assignment_cutoff {
    //                 // consider both ase and hete snps
    //                 if self.fragments[k].read_id == "m84036_230422_223801_s1/42010068/ccs/6216_8325".to_string() {
    //                     println!("S2");
    //                     println!("q_hete:{}, qn_hete:{}, q_ase:{}, qn_ase:{}", q_hete, qn_hete, q_ase, qn_ase);
    //                 }
    //                 if q_hete >= qn_hete && q_ase >= qn_ase {
    //                     // both ase and hete snps support the same haplotype
    //                     if sigma_k == 1 {
    //                         self.fragments[k].assignment = 1;
    //                         self.fragments[k].assignment_score = q_hete.max(q_ase);
    //                         read_assignments.insert(self.fragments[k].read_id.clone(), 1);
    //                         if self.fragments[k].read_id == "m84036_230422_223801_s1/42010068/ccs/6216_8325".to_string() {
    //                             println!("S2.1");
    //                         }
    //                     } else {
    //                         self.fragments[k].assignment = 2;
    //                         self.fragments[k].assignment_score = q_hete.max(q_ase);
    //                         read_assignments.insert(self.fragments[k].read_id.clone(), 2);
    //                         if self.fragments[k].read_id == "m84036_230422_223801_s1/42010068/ccs/6216_8325".to_string() {
    //                             println!("S2.2");
    //                         }
    //                     }
    //                 } else if q_hete < qn_hete && q_ase < qn_ase {
    //                     // both ase and hete snps support the opposite haplotype
    //                     if sigma_k == 1 {
    //                         self.fragments[k].assignment = 2;
    //                         self.fragments[k].assignment_score = qn_hete.max(qn_ase);
    //                         self.fragments[k].haplotag = -1;
    //                         read_assignments.insert(self.fragments[k].read_id.clone(), 2);
    //                     } else {
    //                         self.fragments[k].assignment = 1;
    //                         self.fragments[k].assignment_score = qn_hete.max(qn_ase);
    //                         self.fragments[k].haplotag = 1;
    //                         read_assignments.insert(self.fragments[k].read_id.clone(), 1);
    //                     }
    //                 } else if q_hete >= qn_hete && q_ase < qn_ase {
    //                     // hete snps and ase snps have conflict read assignment
    //                     if delta_hete.len() >= 2 {
    //                         // if more than 2 hete snps, use hete snps to assign the read
    //                         if sigma_k == 1 {
    //                             self.fragments[k].assignment = 1;
    //                             self.fragments[k].assignment_score = q_hete;
    //                             read_assignments.insert(self.fragments[k].read_id.clone(), 1);
    //                         } else {
    //                             self.fragments[k].assignment = 2;
    //                             self.fragments[k].assignment_score = q_hete;
    //                             read_assignments.insert(self.fragments[k].read_id.clone(), 2);
    //                         }
    //                     } else {
    //                         // if only 1 hete snp, use the higher probability to assign the read
    //                         if q_hete >= qn_ase {
    //                             if sigma_k == 1 {
    //                                 self.fragments[k].assignment = 1;
    //                                 self.fragments[k].assignment_score = q_hete;
    //                                 read_assignments.insert(self.fragments[k].read_id.clone(), 1);
    //                             } else {
    //                                 self.fragments[k].assignment = 2;
    //                                 self.fragments[k].assignment_score = q_hete;
    //                                 read_assignments.insert(self.fragments[k].read_id.clone(), 2);
    //                             }
    //                         } else {
    //                             if sigma_k == 1 {
    //                                 self.fragments[k].assignment = 2;
    //                                 self.fragments[k].assignment_score = qn_ase;
    //                                 self.fragments[k].haplotag = -1;
    //                                 read_assignments.insert(self.fragments[k].read_id.clone(), 2);
    //                             } else {
    //                                 self.fragments[k].assignment = 1;
    //                                 self.fragments[k].assignment_score = qn_ase;
    //                                 self.fragments[k].haplotag = 1;
    //                                 read_assignments.insert(self.fragments[k].read_id.clone(), 1);
    //                             }
    //                         }
    //                     }
    //                 } else if q_hete < qn_hete && q_ase >= qn_ase {
    //                     // hete snps and ase snps have conflict read assignment
    //                     if delta_hete.len() >= 2 {
    //                         // if more than 2 hete snps, use hete snps to assign the read
    //                         if sigma_k == 1 {
    //                             self.fragments[k].assignment = 2;
    //                             self.fragments[k].assignment_score = qn_hete;
    //                             self.fragments[k].haplotag = -1;
    //                             read_assignments.insert(self.fragments[k].read_id.clone(), 2);
    //                         } else {
    //                             self.fragments[k].assignment = 1;
    //                             self.fragments[k].assignment_score = qn_hete;
    //                             self.fragments[k].haplotag = 1;
    //                             read_assignments.insert(self.fragments[k].read_id.clone(), 1);
    //                         }
    //                     } else {
    //                         // if only 1 hete snp, use the higher probability to assign the read
    //                         if qn_hete >= q_ase {
    //                             if sigma_k == 1 {
    //                                 self.fragments[k].assignment = 2;
    //                                 self.fragments[k].assignment_score = qn_hete;
    //                                 self.fragments[k].haplotag = -1;
    //                                 read_assignments.insert(self.fragments[k].read_id.clone(), 2);
    //                             } else {
    //                                 self.fragments[k].assignment = 1;
    //                                 self.fragments[k].assignment_score = qn_hete;
    //                                 self.fragments[k].haplotag = 1;
    //                                 read_assignments.insert(self.fragments[k].read_id.clone(), 1);
    //                             }
    //                         } else {
    //                             if sigma_k == 1 {
    //                                 self.fragments[k].assignment = 1;
    //                                 self.fragments[k].assignment_score = q_ase;
    //                                 read_assignments.insert(self.fragments[k].read_id.clone(), 1);
    //                             } else {
    //                                 self.fragments[k].assignment = 2;
    //                                 self.fragments[k].assignment_score = q_ase;
    //                                 read_assignments.insert(self.fragments[k].read_id.clone(), 2);
    //                             }
    //                         }
    //                     }
    //                 } else {
    //                     // panic!("Error: unexpected condition.");
    //                     // unknown which haplotype the read belongs to, cluster the read into unknown group
    //                     self.fragments[k].assignment = 0;
    //                     self.fragments[k].assignment_score = 0.0;
    //                     read_assignments.insert(self.fragments[k].read_id.clone(), 0);
    //                 }
    //             } else if (q_hete - qn_hete).abs() >= read_assignment_cutoff {
    //                 // only consider hete snps
    //                 if self.fragments[k].read_id == "m84036_230422_223801_s1/42010068/ccs/6216_8325".to_string() {
    //                     println!("S3");
    //                     println!("q_hete:{}, qn_hete:{}, q_ase:{}, qn_ase:{}", q_hete, qn_hete, q_ase, qn_ase);
    //                 }
    //                 if q_hete > qn_hete {
    //                     if sigma_k == 1 {
    //                         self.fragments[k].assignment = 1;
    //                         self.fragments[k].assignment_score = q_hete;
    //                         read_assignments.insert(self.fragments[k].read_id.clone(), 1);
    //                     } else {
    //                         self.fragments[k].assignment = 2;
    //                         self.fragments[k].assignment_score = q_hete;
    //                         read_assignments.insert(self.fragments[k].read_id.clone(), 2);
    //                     }
    //                 } else {
    //                     if sigma_k == 1 {
    //                         self.fragments[k].assignment = 2;
    //                         self.fragments[k].assignment_score = qn_hete;
    //                         self.fragments[k].haplotag = -1;
    //                         read_assignments.insert(self.fragments[k].read_id.clone(), 2);
    //                     } else {
    //                         self.fragments[k].assignment = 1;
    //                         self.fragments[k].assignment_score = qn_hete;
    //                         self.fragments[k].haplotag = 1;
    //                         read_assignments.insert(self.fragments[k].read_id.clone(), 1);
    //                     }
    //                 }
    //             } else if (q_ase - qn_ase).abs() >= read_assignment_cutoff {
    //                 // only consider ase snps
    //                 if self.fragments[k].read_id == "m84036_230422_223801_s1/42010068/ccs/6216_8325".to_string() {
    //                     println!("S4");
    //                     println!("q_hete:{}, qn_hete:{}, q_ase:{}, qn_ase:{}", q_hete, qn_hete, q_ase, qn_ase);
    //                 }
    //                 if q_ase > qn_ase {
    //                     if sigma_k == 1 {
    //                         self.fragments[k].assignment = 1;
    //                         self.fragments[k].assignment_score = q_ase;
    //                         read_assignments.insert(self.fragments[k].read_id.clone(), 1);
    //                     } else {
    //                         self.fragments[k].assignment = 2;
    //                         self.fragments[k].assignment_score = q_ase;
    //                         read_assignments.insert(self.fragments[k].read_id.clone(), 2);
    //                     }
    //                 } else {
    //                     if sigma_k == 1 {
    //                         self.fragments[k].assignment = 2;
    //                         self.fragments[k].assignment_score = qn_ase;
    //                         self.fragments[k].haplotag = -1;
    //                         read_assignments.insert(self.fragments[k].read_id.clone(), 2);
    //                     } else {
    //                         self.fragments[k].assignment = 1;
    //                         self.fragments[k].assignment_score = qn_ase;
    //                         self.fragments[k].haplotag = 1;
    //                         read_assignments.insert(self.fragments[k].read_id.clone(), 1);
    //                     }
    //                 }
    //             } else {
    //                 // unknown which haplotype the read belongs to, cluster the read into unknown group
    //                 if self.fragments[k].read_id == "m84036_230422_223801_s1/42010068/ccs/6216_8325".to_string() {
    //                     println!("S5");
    //                     println!("q_hete:{}, qn_hete:{}, q_ase:{}, qn_ase:{}", q_hete, qn_hete, q_ase, qn_ase);
    //                 }
    //                 // panic!("Error: unexpected condition.");
    //                 self.fragments[k].assignment = 0;
    //                 self.fragments[k].assignment_score = 0.0;
    //                 read_assignments.insert(self.fragments[k].read_id.clone(), 0);
    //             }
    //         }
    //     }
    //     return read_assignments;
    // }


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
                // let mut is_low_qual = false;
                // let mut is_dense = false;
                // let mut is_rna_edit = false;
                // let mut is_single_snp = false;
                // let mut is_unconfident_phased_snp = false;
                // let mut is_ase_snp = false;
                if snp.dense || snp.single || snp.rna_editing || snp.somatic || snp.phase_score == 0.0 {
                    continue;
                }
                // if snp.single == true {
                //     is_single_snp = true;
                // }
                // if snp.ase == true {
                //     is_ase_snp = true;
                // }
                // if !is_dense && !is_single_snp && snp.phase_score == 0.0 {
                //     is_unconfident_phased_snp = true;
                // }
                // if is_dense || is_single_snp || is_unconfident_phased_snp || snp.haplotype == 0 {
                //     continue;
                // }
                // // ase snp
                // if is_ase_snp && snp.phase_score < ase_ps_cutoff as f64 {
                //     continue;
                // }
                // // hete snp
                // if !is_ase_snp {
                //     if snp.variant_quality < min_qual_for_candidate as f64 {
                //         continue;
                //     }
                //     if snp.phase_score < min_phase_score as f64 {
                //         continue;
                //     }
                // }
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

    pub fn detect_somatic_by_het(&mut self, bam_path: &str, region: &Region) {
        if self.somatic_snps.len() == 0 {
            return;
        }
        // detect confident somatic mutation with phasing result of high allele fraction het snps.
        // 1. assign phased result for candidate somatic sites.
        let mut phased_fragments: HashMap<String, i32> = HashMap::new(); // read id, assignment
        for k in 0..self.fragments.len() {
            let frag = &self.fragments[k];
            if frag.assignment == 1 || frag.assignment == 2 {
                phased_fragments.insert(frag.read_id.clone(), frag.assignment);
            }
        }
        self.get_somatic_haplotype_baseqs(bam_path, region, &phased_fragments);
        // 2. find candidates meet the criteria of somatic mutation. haplotype-specific
        for i in 0..self.somatic_snps.len() {
            let som_cand = &mut self.candidate_snps[self.somatic_snps[i]];
            let (hap1_allele_class, hap2_allele_class) = calculate_prob_somatic(&som_cand.hap_quals.hap1_ref_baseqs, &som_cand.hap_quals.hap1_alt_baseqs, &som_cand.hap_quals.hap2_ref_baseqs, &som_cand.hap_quals.hap2_alt_baseqs, 0.3);
            if hap1_allele_class.allcls == 0 && hap2_allele_class.allcls == 2 {
                let somatic_score = -10.0_f64 * (1.0 - hap2_allele_class.prob).log10();
                som_cand.somatic = true;
                som_cand.variant_type = 1;
                som_cand.somatic_score = somatic_score;
                // println!("somatic snp:{}, score: {}", som_cand.pos, somatic_score);
                // println!("{:?},{:?}", hap1_allele_class, hap2_allele_class);
                // println!("hap1_ref_baseqs:{:?}\nhap1_alt_baseqs:{:?}\nhap2_ref_baseqs:{:?}\nhap2_alt_baseqs:{:?}", som_cand.hap_quals.hap1_ref_baseqs, som_cand.hap_quals.hap1_alt_baseqs, som_cand.hap_quals.hap2_ref_baseqs, som_cand.hap_quals.hap2_alt_baseqs);
            } else if hap1_allele_class.allcls == 2 && hap2_allele_class.allcls == 0 {
                let somatic_score = -10.0_f64 * (1.0 - hap1_allele_class.prob).log10();
                som_cand.somatic = true;
                som_cand.variant_type = 1;
                som_cand.somatic_score = somatic_score;
                // println!("somatic snp:{}, score: {}", som_cand.pos, somatic_score);
                // println!("{:?},{:?}", hap1_allele_class, hap2_allele_class);
                // println!("hap1_ref_baseqs:{:?}\nhap1_alt_baseqs:{:?}\nhap2_ref_baseqs:{:?}\nhap2_alt_baseqs:{:?}", som_cand.hap_quals.hap1_ref_baseqs, som_cand.hap_quals.hap1_alt_baseqs, som_cand.hap_quals.hap2_ref_baseqs, som_cand.hap_quals.hap2_alt_baseqs);
            }
        }
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
                    if fe.phase_site == false { continue; }
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
            for i in self.high_frac_het_snps.iter() {
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
                            assert_ne!(fe.phase_site, false, "Error: phase for non-phase site.");
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
            for i in self.high_frac_het_snps.iter() {
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

    // pub fn rescue_ase_snps(&mut self) {
    //     // use surrounding heterozygous snps to evaluate ase snps
    //     // TODO: add phase score still has problem
    //     for i in self.ase_hete_snps.iter() {
    //         if self.candidate_snps[*i].ase == true {
    //             let mut sigma: Vec<i32> = Vec::new();
    //             let mut ps: Vec<i32> = Vec::new();
    //             let mut probs: Vec<f64> = Vec::new();
    //             let mut num_hap1 = 0;
    //             let mut num_hap2 = 0;
    //             for k in self.candidate_snps[*i].snp_cover_fragments.iter() {
    //                 // k is fragment index
    //                 if self.fragments[*k].assignment == 0 {
    //                     continue;
    //                 }
    //                 for fe in self.fragments[*k].list.iter() {
    //                     if fe.snp_idx == *i {
    //                         if self.fragments[*k].assignment == 1 {
    //                             num_hap1 += 1;
    //                         } else if self.fragments[*k].assignment == 2 {
    //                             num_hap2 += 1;
    //                         }
    //                         assert_ne!(fe.p, 0, "Error: phase for unexpected allele.");
    //                         ps.push(fe.p);
    //                         probs.push(fe.prob);
    //                         sigma.push(self.fragments[*k].haplotag);
    //                     }
    //                 }
    //             }
    //             if num_hap1 < 3 || num_hap2 < 3 {
    //                 continue;
    //             }
    //             if sigma.len() > 10 {
    //                 let phase_score1 = -10.0_f64 * (1.0 - SNPFrag::cal_delta_sigma_log(1, &sigma, &ps, &probs)).log10();
    //                 let phase_score2 = -10.0_f64 * (1.0 - SNPFrag::cal_delta_sigma_log(-1, &sigma, &ps, &probs)).log10();
    //                 if f64::max(phase_score1, phase_score2) > 40.0 {
    //                     if phase_score1 > phase_score2 {
    //                         self.candidate_snps[*i].haplotype = 1;
    //                     } else {
    //                         self.candidate_snps[*i].haplotype = -1;
    //                     }
    //                     println!(
    //                         "Rescue ASE SNP: {} {} {} {}",
    //                         String::from_utf8(self.candidate_snps[*i].chromosome.clone()).unwrap(),
    //                         self.candidate_snps[*i].pos,
    //                         phase_score1,
    //                         phase_score2
    //                     );
    //                     // self.candidate_snps[*i].ase = false;
    //                 }
    //             }
    //         }
    //     }
    // }

    // pub fn rescue_ase_snps_v2(
    //     &mut self,
    //     ase_allele_cnt_cutoff: u32,
    //     ase_ps_count_cutoff: u32,
    //     ase_ps_cutoff: f32,
    // ) {
    //     // keep the previously calculated phasing for heterozygous snps and only flip the ase snps to get the phasing of ase snps
    //     for i in self.ase_snps.iter() {
    //         assert!(self.candidate_snps[*i].ase == true, "rescue not ase snp.");
    //         let mut sigma: Vec<i32> = Vec::new();
    //         let mut ps: Vec<i32> = Vec::new();
    //         let mut probs: Vec<f64> = Vec::new();
    //         let mut num_hap1 = 0;
    //         let mut num_hap2 = 0;
    //         let debug_pos = 2341049;
    //         if self.candidate_snps[*i].pos == debug_pos - 1 {
    //             println!("rescue_ase_snps_v2");
    //         }
    //         if self.candidate_snps[*i].pos == debug_pos - 1 {
    //             for k in self.candidate_snps[*i].snp_cover_fragments.iter() {
    //                 println!("cover: {:?}", self.fragments[*k]);
    //             }
    //         }
    //         for k in self.candidate_snps[*i].snp_cover_fragments.iter() {
    //             // k is fragment index
    //             if self.fragments[*k].assignment == 0 {
    //                 continue;
    //             }
    //             for fe in self.fragments[*k].list.iter() {
    //                 if fe.snp_idx == *i {
    //                     if self.fragments[*k].assignment == 1 {
    //                         num_hap1 += 1;
    //                     } else if self.fragments[*k].assignment == 2 {
    //                         num_hap2 += 1;
    //                     }
    //                     assert_ne!(fe.p, 0, "Error: phase for unexpected allele.");
    //                     ps.push(fe.p);
    //                     probs.push(fe.prob);
    //                     sigma.push(self.fragments[*k].haplotag);
    //                 }
    //             }
    //         }
    //         if self.candidate_snps[*i].pos == debug_pos - 1 {
    //             println!("rescue_ase_snps_v2:");
    //             println!("ps:{:?}\nprobs:{:?}\nsigma:{:?}\n", ps, probs, sigma);
    //             println!("num_hap1:{}, num_hap2:{}", num_hap1, num_hap2);
    //         }
    //         if num_hap1 < ase_allele_cnt_cutoff || num_hap2 < ase_allele_cnt_cutoff {
    //             // each haplotype should have at least 3 reads
    //             continue;
    //         }
    //         if sigma.len() >= ase_ps_count_cutoff as usize {
    //             // each site should have at least 10 reads
    //             let phase_score1 = -10.0_f64 * (1.0 - SNPFrag::cal_delta_sigma_log(1, &sigma, &ps, &probs)).log10();
    //             let phase_score2 = -10.0_f64 * (1.0 - SNPFrag::cal_delta_sigma_log(-1, &sigma, &ps, &probs)).log10();
    //             if self.candidate_snps[*i].pos == debug_pos - 1 {
    //                 println!("phase score: {},{}", phase_score1, phase_score2);
    //             }
    //             if f64::max(phase_score1, phase_score2) > ase_ps_cutoff as f64 {
    //                 if phase_score1 > phase_score2 {
    //                     self.candidate_snps[*i].haplotype = 1;
    //                     self.candidate_snps[*i].phase_score = phase_score1;
    //                     self.candidate_snps[*i].variant_type = 1;
    //                 } else {
    //                     self.candidate_snps[*i].haplotype = -1;
    //                     self.candidate_snps[*i].phase_score = phase_score2;
    //                     self.candidate_snps[*i].variant_type = 1;
    //                 }
    //             }
    //             let mut haplotype_allele_expression: [u32; 4] = [0, 0, 0, 0];   // hap1_ref, hap1_alt, hap2_ref, hap2_alt
    //             for k in 0..sigma.len() {
    //                 if sigma[k] == 1 {
    //                     // hap1
    //                     if ps[k] == 1 {
    //                         haplotype_allele_expression[0] += 1;
    //                     } else if ps[k] == -1 {
    //                         haplotype_allele_expression[1] += 1;
    //                     }
    //                 } else if sigma[k] == -1 {
    //                     // hap2
    //                     if ps[k] == 1 {
    //                         haplotype_allele_expression[2] += 1;
    //                     } else if ps[k] == -1 {
    //                         haplotype_allele_expression[3] += 1;
    //                     }
    //                 }
    //             }
    //             self.candidate_snps[*i].haplotype_expression = haplotype_allele_expression;
    //             // TODO: If phase score does not meet ase_ps_cutoff, detect whether candidate with low phase score is a potential somatic mutation
    //         }
    //     }
    // }
}