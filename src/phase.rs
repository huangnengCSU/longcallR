use std::collections::{HashMap, HashSet};

use rand::Rng;

use crate::snpfrags::SNPFrag;

pub fn qki(sigma: i32, delta: i32, alt_fraction: f32, base_allele: i32, error_rate: f64) -> f64 {
    let mut x: i32;
    if delta == 0 {
        if alt_fraction < 0.5 {
            x = 1;
        } else {
            x = -1;
        }
    } else {
        x = sigma * delta;
    }
    if base_allele == x {
        return 1.0 - error_rate;
    } else {
        return error_rate;
    }
}

pub fn aki(sigma: i32, delta: i32, eta: i32, base_allele: i32, error_rate: f64) -> f64 {
    // sigma: the assignment of read k, 1 or -1.
    // delta: the haplotype of SNP i, 1 or -1.
    // eta: the genotype of SNP i, 0 (hetVar), 1 (homRef), -1 (homVar).
    // base_allele: the allele of the base, 1 (ref_allele), -1 (alt_allele).
    // error_rate: the probability of observing base at SNP i for read k, equals to 10^(-Q/10).
    let mut x: i32;
    if eta == 0 {
        x = sigma * delta;
    } else {
        x = eta;
    }
    if base_allele == x {
        return 1.0 - error_rate;
    } else {
        return error_rate;
    }
}

pub fn cal_sigma_delta(sigma_k: i32, delta: &Vec<i32>, altfrac: &Vec<f32>, ps: &Vec<i32>, probs: &Vec<f64>) -> f64 {
    // calculate P(sigma_k | delta)
    // sigma_k: the assignment of read k, 1 or -1.
    // delta: the haplotypes of the SNPs covered by read k, each haplotype is 1 or -1 or 0.
    // ps: the allele of each base, 1,-1
    // probs: the probability of observing base at each SNP for read k, equals to 10^(-Q/10).
    let mut log_q1: f64 = 0.0;
    let mut log_q2: f64 = 0.0;
    let mut log_q3: f64 = 0.0;

    for i in 0..delta.len() {
        log_q1 += qki(sigma_k, delta[i], altfrac[i], ps[i], probs[i]).log10();
    }

    for i in 0..delta.len() {
        log_q2 += qki(1, delta[i], altfrac[i], ps[i], probs[i]).log10();
        log_q3 += qki(-1, delta[i], altfrac[i], ps[i], probs[i]).log10();
    }
    let max_logq = log_q1.max(log_q2.max(log_q3));
    let q1 = 10.0_f64.powf(log_q1 - max_logq);
    let q2 = 10.0_f64.powf(log_q2 - max_logq);
    let q3 = 10.0_f64.powf(log_q3 - max_logq);
    // println!("sigma delta q1:{:?}, q2+q3:{:?}", q1, q2 + q3);
    return q1 / (q2 + q3);
}

pub fn cal_sigma_delta_eta_log(sigma_k: i32, delta: &Vec<i32>, eta: &Vec<i32>, ps: &Vec<i32>, probs: &Vec<f64>) -> f64 {
    // same as call_sigma_delta, but use log to avoid underflow
    let mut log_q1: f64 = 0.0;
    let mut log_q2: f64 = 0.0;
    let mut log_q3: f64 = 0.0;

    for i in 0..delta.len() {
        log_q1 += aki(sigma_k, delta[i], eta[i], ps[i], probs[i]).log10();
    }

    for i in 0..delta.len() {
        log_q2 += aki(1, delta[i], eta[i], ps[i], probs[i]).log10();
        log_q3 += aki(-1, delta[i], eta[i], ps[i], probs[i]).log10();
    }
    /*
    Given 0<=A<=1, 0<=B<=1 and A/(A+B) > B/(A+B), we have logA/(logA+logB) < logB/(logA+logB).
    When comparing A/(A+B) and B/(A+B), we can approximately use 1-log(A)/(log(A)+log(B)) and 1-log(B)/(log(A)+log(B)) to avoid underflow.
    */
    return 1.0 - log_q1 / (log_q2 + log_q3);
}

pub fn cal_delta_sigma(delta_i: i32, alt_fraction_i: f32, sigma: &Vec<i32>, ps: &Vec<i32>, probs: &Vec<f64>) -> f64 {
    // calculate P(delta_i | sigma)
    // delta_i: the haplotype of SNP i, 1 or -1.
    // sigma: the assignments of the reads cover SNP i, each haplotype is 1 or -1.
    // ps: the allele of each base, 1,-1
    // probs: the probability of observing base at SNP i for each read, equals to 10^(-Q/10).

    let mut log_q1: f64 = 0.0;
    let mut log_q2: f64 = 0.0;
    let mut log_q3: f64 = 0.0;
    let mut log_q4: f64 = 0.0;

    for k in 0..sigma.len() {
        log_q1 += qki(sigma[k], delta_i, alt_fraction_i, ps[k], probs[k]).log10();
    }

    for k in 0..sigma.len() {
        log_q2 += qki(sigma[k], 1, alt_fraction_i, ps[k], probs[k]).log10();
        log_q3 += qki(sigma[k], -1, alt_fraction_i, ps[k], probs[k]).log10();
        log_q4 += qki(sigma[k], 0, alt_fraction_i, ps[k], probs[k]).log10();
    }
    let max_logq = log_q1.max(log_q2.max(log_q3.max(log_q4)));
    let q1 = 10.0_f64.powf(log_q1 - max_logq);
    let q2 = 10.0_f64.powf(log_q2 - max_logq);
    let q3 = 10.0_f64.powf(log_q3 - max_logq);
    let q4 = 10.0_f64.powf(log_q4 - max_logq);
    // println!("delta sigma q1:{:?}, q2+q3:{:?}", q1, q2 + q3);
    return q1 / (q2 + q3 + q4);
}

pub fn cal_delta_eta_sigma_log(delta_i: i32, eta_i: i32, sigma: &Vec<i32>, ps: &Vec<i32>, probs: &Vec<f64>) -> f64 {
    // same as call_delta_sigma, but use log to avoid underflow
    let mut log_q1: f64 = 0.0;
    let mut log_q2: f64 = 0.0;
    let mut log_q3: f64 = 0.0;


    for k in 0..sigma.len() {
        log_q1 += aki(sigma[k], delta_i, eta_i, ps[k], probs[k]).log10();
    }

    for k in 0..sigma.len() {
        log_q2 += aki(sigma[k], 1, eta_i, ps[k], probs[k]).log10();
        log_q3 += aki(sigma[k], -1, eta_i, ps[k], probs[k]).log10();
    }

    /*
    Given 0<=A<=1, 0<=B<=1 and A/(A+B) > B/(A+B), we have logA/(logA+logB) < logB/(logA+logB).
    When comparing A/(A+B) and B/(A+B), we can approximately use 1-log(A)/(log(A)+log(B)) and 1-log(B)/(log(A)+log(B)) to avoid underflow.
    */
    return 1.0 - log_q1 / (log_q2 + log_q3);
}

pub fn cal_eta_delta_sigma_log(eta_i: i32, coverage_i: u32, delta_i: i32, sigma: &Vec<i32>, ps: &Vec<i32>, probs: &Vec<f64>) -> f64 {
    let mut log_q1: f64 = 0.0;
    let mut log_q2: f64 = 0.0;
    let mut log_q3: f64 = 0.0;
    let mut log_q4: f64 = 0.0;

    let prior_homref_log: f64 = (1.0 - 1.5 * 0.001 as f64).log10();
    let prior_homvar_log: f64 = (0.5 * 0.001 as f64).log10();
    let prior_hetvar_log: f64;
    if coverage_i == 0 {
        prior_hetvar_log = 0.001_f64.log10();
    } else {
        prior_hetvar_log = 0.001_f64.log10() - (coverage_i as f64) * 2.0_f64.log10();
    }

    for k in 0..sigma.len() {
        log_q1 += aki(sigma[k], delta_i, eta_i, ps[k], probs[k]).log10();
    }
    if eta_i == 0 {
        log_q1 += prior_hetvar_log;
    } else if eta_i == 1 {
        log_q1 += prior_homref_log;
    } else {
        log_q1 += prior_homvar_log;
    }

    for k in 0..sigma.len() {
        log_q2 += aki(sigma[k], delta_i, -1, ps[k], probs[k]).log10();
        log_q3 += aki(sigma[k], delta_i, 0, ps[k], probs[k]).log10();
        log_q4 += aki(sigma[k], delta_i, 1, ps[k], probs[k]).log10();
    }
    log_q2 += prior_homvar_log;
    log_q3 += prior_hetvar_log;
    log_q4 += prior_homref_log;

    /*
    Given 0<=A<=1, 0<=B<=1 and A/(A+B) > B/(A+B), we have logA/(logA+logB) < logB/(logA+logB).
    When comparing A/(A+B) and B/(A+B), we can approximately use 1-log(A)/(log(A)+log(B)) and 1-log(B)/(log(A)+log(B)) to avoid underflow.
    */
    return 1.0 - log_q1 / (log_q2 + log_q3 + log_q4);
}

pub fn cal_delta_sigma_prior_log(delta_i: i32, alt_fraction_i: f32, coverage_i: u32, sigma: &Vec<i32>, ps: &Vec<i32>, probs: &Vec<f64>) -> f64 {
    // same as call_delta_sigma, but use log to avoid underflow
    let mut log_q1: f64 = 0.0;
    let mut log_q2: f64 = 0.0;
    let mut log_q3: f64 = 0.0;
    let mut log_q4: f64 = 0.0;


    let prior_homref_log: f64 = (1.0 - 1.5 * 0.001 as f64).log10();
    let prior_homvar_log: f64 = (0.5 * 0.001 as f64).log10();
    let prior_hetvar_log: f64;
    if coverage_i == 0 {
        prior_hetvar_log = 0.001_f64.log10();
    } else {
        prior_hetvar_log = 0.001_f64.log10() - (coverage_i as f64) * 2.0_f64.log10();
    }
    // let prior_homref_log: f64 = (1.0 - 1.5 * 0.001 as f64).log10();
    // let prior_homvar_log: f64 = (0.5 * 0.001 as f64).log10();
    // let prior_hetvar_log: f64 = 0.001_f64.log10();


    for k in 0..sigma.len() {
        log_q1 += qki(sigma[k], delta_i, alt_fraction_i, ps[k], probs[k]).log10();
    }

    if delta_i == 0 {
        if alt_fraction_i < 0.5 {
            log_q1 += prior_homref_log;
        } else {
            log_q1 += prior_homvar_log;
        }
    } else {
        log_q1 += prior_hetvar_log;
    }

    for k in 0..sigma.len() {
        log_q2 += qki(sigma[k], 1, alt_fraction_i, ps[k], probs[k]).log10();
        log_q3 += qki(sigma[k], -1, alt_fraction_i, ps[k], probs[k]).log10();
        log_q4 += qki(sigma[k], 0, alt_fraction_i, ps[k], probs[k]).log10();
    }

    log_q2 += prior_hetvar_log;
    log_q3 += prior_hetvar_log;
    if alt_fraction_i < 0.5 {
        log_q4 += prior_homref_log;
    } else {
        log_q4 += prior_homvar_log;
    }

    /*
    Given 0<=A<=1, 0<=B<=1 and A/(A+B) > B/(A+B), we have logA/(logA+logB) < logB/(logA+logB).
    When comparing A/(A+B) and B/(A+B), we can approximately use 1-log(A)/(log(A)+log(B)) and 1-log(B)/(log(A)+log(B)) to avoid underflow.
    */
    return 1.0 - log_q1 / (log_q2 + log_q3 + log_q4);
}

pub fn cal_phase_score_log(delta_i: i32, eta_i: i32, sigma: &Vec<i32>, ps: &Vec<i32>, probs: &Vec<f64>) -> f64 {
    // same as call_delta_sigma, but use log to avoid underflow
    let mut log_q1: f64 = 0.0;
    let mut log_q2: f64 = 0.0;
    let mut log_q3: f64 = 0.0;

    assert!(delta_i != 0, "Error: phase for unexpected allele.");
    /*// theoretical calculation
    for k in 0..sigma.len() {
        log_q1 += qki(sigma[k], delta_i, alt_fraction_i, ps[k], probs[k]).log10();
    }
    for k in 0..sigma.len() {
        log_q2 += qki(sigma[k], 1, alt_fraction_i, ps[k], probs[k]).log10();
        log_q3 += qki(sigma[k], -1, alt_fraction_i, ps[k], probs[k]).log10();
    }
    let max_logq = log_q1.max(log_q2.max(log_q3));
    log_q1 = log_q1 - max_logq;
    log_q2 = log_q2 - max_logq;
    log_q3 = log_q3 - max_logq;
    let q1 = 10.0_f64.powf(log_q1);
    let q2 = 10.0_f64.powf(log_q2);
    let q3 = 10.0_f64.powf(log_q3);
    return q1 / (q2 + q3);*/
    // calculation in practice
    for k in 0..sigma.len() {
        log_q1 += aki(sigma[k], delta_i, eta_i, ps[k], probs[k]).log10();
    }
    for k in 0..sigma.len() {
        log_q2 += aki(sigma[k], 1, eta_i, ps[k], probs[k]).log10();
        log_q3 += aki(sigma[k], -1, eta_i, ps[k], probs[k]).log10();
    }
    return 1.0 - log_q1 / (log_q2 + log_q3);
}

pub fn cal_overall_probability(snpfrag: &SNPFrag) -> f64 {
    // calculate the log10 probability of the current configuration of sigma and delta
    let mut logp = 0.0;
    for k in 0..snpfrag.fragments.len() {
        if snpfrag.fragments[k].discarded == true {
            continue;
        }
        if snpfrag.fragments[k].haplotag == 0 {
            // unassigned fragment
            continue;
        }
        for fe in snpfrag.fragments[k].list.iter() {
            if fe.phase_site == false {
                continue;
            }
            assert_ne!(fe.p, 0, "Error: phasing with unexpected hete SNP.");
            logp += aki(snpfrag.fragments[k].haplotag, snpfrag.candidate_snps[fe.snp_idx].haplotype, snpfrag.candidate_snps[fe.snp_idx].genotype, fe.p, fe.prob).log10();
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
        let mut eta: Vec<i32> = Vec::new();
        let mut ps: Vec<i32> = Vec::new();
        let mut probs: Vec<f64> = Vec::new();
        if snpfrag.fragments[*k].discarded == true {
            continue;
        }
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
            eta.push(snpfrag.candidate_snps[fe.snp_idx].genotype);
        }
        if delta.len() == 0 {
            continue;
        }
        logp += cal_sigma_delta_eta_log(*h, &delta, &eta, &ps, &probs);
        pre_logp += cal_sigma_delta_eta_log(snpfrag.fragments[*k].haplotag, &delta, &eta, &ps, &probs);
    }

    let mut flag = 0;
    if logp > pre_logp {
        flag = 1;
    } else if logp == pre_logp {
        flag = 0;
    } else {
        flag = -1;
    }
    assert!(flag >= 0, "Error: new haplotag decrease the probability. {} -> {}", pre_logp, logp);
    return flag;
}

pub fn check_new_haplotype(snpfrag: &SNPFrag, updated_haplotype: &HashMap<usize, i32>) -> i32 {
    let mut logp = 0.0;
    let mut pre_logp = 0.0;
    for (i, h) in updated_haplotype.iter() {
        let eta_i = snpfrag.candidate_snps[*i].genotype;
        let coverage_i = snpfrag.candidate_snps[*i].depth;
        let mut sigma: Vec<i32> = Vec::new();
        let mut ps: Vec<i32> = Vec::new();
        let mut probs: Vec<f64> = Vec::new();
        for k in snpfrag.candidate_snps[*i].snp_cover_fragments.iter() {
            if snpfrag.fragments[*k].discarded == true {
                continue;
            }
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
            continue;
        }
        logp += cal_delta_eta_sigma_log(*h, eta_i, &sigma, &ps, &probs);
        pre_logp += cal_delta_eta_sigma_log(snpfrag.candidate_snps[*i].haplotype, eta_i, &sigma, &ps, &probs);
    }
    let mut flag = 0;
    if logp > pre_logp {
        flag = 1;
    } else if logp == pre_logp {
        flag = 0;
    } else {
        flag = -1;
    }
    assert!(flag >= 0, "Error: new haplotype decrease the probability. {} -> {}", pre_logp, logp);
    return flag;
}

pub fn check_new_genotype(snpfrag: &SNPFrag, updated_genotype: &HashMap<usize, i32>) -> i32 {
    let mut logp = 0.0;
    let mut pre_logp = 0.0;
    for (i, g) in updated_genotype.iter() {
        let delta_i = snpfrag.candidate_snps[*i].haplotype;
        let coverage_i = snpfrag.candidate_snps[*i].depth;
        let mut sigma: Vec<i32> = Vec::new();
        let mut ps: Vec<i32> = Vec::new();
        let mut probs: Vec<f64> = Vec::new();
        for k in snpfrag.candidate_snps[*i].snp_cover_fragments.iter() {
            if snpfrag.fragments[*k].discarded == true {
                continue;
            }
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
            continue;
        }
        logp += cal_eta_delta_sigma_log(*g, coverage_i, delta_i, &sigma, &ps, &probs);
        pre_logp += cal_eta_delta_sigma_log(snpfrag.candidate_snps[*i].genotype, coverage_i, delta_i, &sigma, &ps, &probs);
    }
    let mut flag = 0;
    if logp > pre_logp {
        flag = 1;
    } else if logp == pre_logp {
        flag = 0;
    } else {
        flag = -1;
    }
    assert!(flag >= 0, "Error: new genotype decrease the probability. {} -> {}", pre_logp, logp);
    return flag;
}

pub fn cal_inconsistent_percentage(delta_i: i32, sigma: &Vec<i32>, ps: &Vec<i32>) -> f64 {
    let mut consisitent = 0;
    let mut inconsistent = 0;
    for k in 0..sigma.len() {
        if delta_i * sigma[k] == ps[k] { consisitent += 1; } else { inconsistent += 1; }
    }
    if consisitent + inconsistent == 0 { return 0.0; }
    return (inconsistent as f64) / ((consisitent + inconsistent) as f64);
}

impl SNPFrag {
    pub unsafe fn init_haplotypes(&mut self) {
        // initialize haplotype of heterozygous snp
        for i in 0..self.candidate_snps.len() {
            let mut rng = rand::thread_rng();
            let rg: f64 = rng.gen();
            if rg < 0.5 {
                self.candidate_snps[i].haplotype = 1; // HetVar hap1
            } else {
                self.candidate_snps[i].haplotype = -1;  // HetVar hap2
            }
        }
    }

    pub unsafe fn init_assignment(&mut self) {
        for k in 0..self.fragments.len() {
            let mut rng = rand::thread_rng();
            if self.fragments[k].discarded == true { 
                continue; 
            }
            // if self.fragments[k].num_hete_links < self.min_linkers {
            //     continue;
            // }
            let rg: f64 = rng.gen();
            if rg < 0.5 {
                self.fragments[k].haplotag = -1;
            } else {
                self.fragments[k].haplotag = 1;
            }
        }
    }

    pub fn init_genotype(&mut self) {
        for i in 0..self.candidate_snps.len() {
            if self.candidate_snps[i].variant_type == 0 {
                self.candidate_snps[i].genotype = 1; // HomRef
            } else if self.candidate_snps[i].variant_type == 1 {
                self.candidate_snps[i].genotype = 0; // HetVar
            } else if self.candidate_snps[i].variant_type == 2 {
                self.candidate_snps[i].genotype = -1; // HomVar
            } else if self.candidate_snps[i].variant_type == 3 {
                self.candidate_snps[i].genotype = -1; // triallelic SNP
            }
        }
    }


    pub fn cross_optimize_sigma(&mut self) -> f64 {
        // Iteration:
        //     1. evaluate the assignment of each read based on the current SNP haplotype.
        //     2. evaluate the SNP haplotype based on the read assignment.
        // If P(sigma, delta) increase, repeat Iteration;
        // Else break;

        let mut phasing_increase: bool = true;
        let mut haplotag_increase: bool = true;
        let mut genotype_increase: bool = true;
        let mut num_iters = 0;

        while phasing_increase | haplotag_increase {
            // optimize sigma
            let mut tmp_haplotag: HashMap<usize, i32> = HashMap::new();
            let mut processed_snps = HashSet::new(); // some snps in self.hete_snps may be filtered by previous steps, record the snps that covered by the fragments
            for k in 0..self.fragments.len() {
                if self.fragments[k].discarded == true { 
                    continue; 
                }
                let sigma_k = self.fragments[k].haplotag;
                let mut delta: Vec<i32> = Vec::new();
                let mut eta: Vec<i32> = Vec::new();
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
                    eta.push(self.candidate_snps[fe.snp_idx].genotype);
                    processed_snps.insert(fe.snp_idx);
                }

                if delta.len() == 0 {
                    continue;
                }
                let q = cal_sigma_delta_eta_log(sigma_k, &delta, &eta, &ps, &probs);
                let qn = cal_sigma_delta_eta_log(sigma_k * (-1), &delta, &eta, &ps, &probs);
                // println!("optimize sigma {} q:{}, qn:{}, sigma: {}", k, q, qn, sigma_k);

                if q < qn {
                    tmp_haplotag.insert(k, sigma_k * (-1));
                } else {
                    tmp_haplotag.insert(k, sigma_k);
                }
            }

            // assert!(SNPFrag::cal_overall_probability(&self, &processed_snps, &self.haplotype) >= SNPFrag::cal_overall_probability(&self, &self.haplotag, &self.haplotype));
            let check_val = check_new_haplotag(&self, &tmp_haplotag);
            assert!(check_val >= 0, "Error: check new haplotag: {:?}", self.candidate_snps);
            for (k, h) in tmp_haplotag.iter() {
                // when prob is equal, we still perform the flip to avoid bug of underflow
                self.fragments[*k].haplotag = *h;
            }
            if check_val == 0 {
                haplotag_increase = false;
            } else {
                haplotag_increase = true;
                phasing_increase = true;
                genotype_increase = true;
            }
            self.check_local_optimal_configuration(false, true, false);

            // optimize delta
            let mut tmp_haplotype: HashMap<usize, i32> = HashMap::new();
            for i in 0..self.candidate_snps.len() {
                if self.candidate_snps[i].for_phasing == false { continue; }
                let delta_i = self.candidate_snps[i].haplotype;
                let eta_i = self.candidate_snps[i].genotype;
                let mut sigma: Vec<i32> = Vec::new();
                let mut ps: Vec<i32> = Vec::new();
                let mut probs: Vec<f64> = Vec::new();
                for k in self.candidate_snps[i].snp_cover_fragments.iter() {
                    if self.fragments[*k].discarded == true {
                        continue;
                    }
                    if self.fragments[*k].haplotag == 0 {
                        continue;
                    }
                    // k is fragment index
                    for fe in self.fragments[*k].list.iter() {
                        if fe.snp_idx == i {
                            if fe.phase_site == false { continue; }
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
                let q = cal_delta_eta_sigma_log(delta_i, eta_i, &sigma, &ps, &probs);
                let qn = cal_delta_eta_sigma_log(delta_i * (-1), eta_i, &sigma, &ps, &probs);

                if q < qn {
                    tmp_haplotype.insert(i, delta_i * (-1));
                } else {
                    tmp_haplotype.insert(i, delta_i);
                }
            }
            let check_val = check_new_haplotype(&self, &tmp_haplotype);
            assert!(check_val >= 0, "Error: check new haplotype: {:?}", self.candidate_snps);
            for (idx, hap) in tmp_haplotype.iter() {
                // when prob is equal, we still perform the flip to avoid bug of underflow
                self.candidate_snps[*idx].haplotype = *hap;
            }
            if check_val == 0 {
                phasing_increase = false;
            } else {
                phasing_increase = true;
                haplotag_increase = true;
                genotype_increase = true;
            }
            self.check_local_optimal_configuration(true, false, false);

            // optimize genotype
            let mut tmp_genotype: HashMap<usize, i32> = HashMap::new();
            for i in 0..self.candidate_snps.len() {
                if self.candidate_snps[i].for_phasing == false { continue; }
                let delta_i = self.candidate_snps[i].haplotype;
                let eta_i = self.candidate_snps[i].genotype;
                let coverage_i = self.candidate_snps[i].depth;
                let mut sigma: Vec<i32> = Vec::new();
                let mut ps: Vec<i32> = Vec::new();
                let mut probs: Vec<f64> = Vec::new();
                for k in self.candidate_snps[i].snp_cover_fragments.iter() {
                    if self.fragments[*k].discarded == true {
                        continue;
                    }
                    if self.fragments[*k].haplotag == 0 {
                        continue;
                    }
                    // k is fragment index
                    for fe in self.fragments[*k].list.iter() {
                        if fe.snp_idx == i {
                            if fe.phase_site == false { continue; }
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
                let q1 = cal_eta_delta_sigma_log(-1, coverage_i, delta_i, &sigma, &ps, &probs);
                let q2 = cal_eta_delta_sigma_log(0, coverage_i, delta_i, &sigma, &ps, &probs);
                let q3 = cal_eta_delta_sigma_log(1, coverage_i, delta_i, &sigma, &ps, &probs);

                if q1 >= q2 && q1 >= q3 {
                    tmp_genotype.insert(i, -1);
                } else if q2 >= q1 && q2 >= q3 {
                    tmp_genotype.insert(i, 0);
                } else if q3 >= q1 && q3 >= q2 {
                    tmp_genotype.insert(i, 1);
                } else {
                    panic!("Error: genotype optimization failed. {}->{}->{}", q1, q2, q3);
                }
            }
            let check_val = check_new_genotype(&self, &tmp_genotype);
            assert!(check_val >= 0, "Error: check new genotype: {:?}", self.candidate_snps);
            for (idx, gt) in tmp_genotype.iter() {
                // when prob is equal, we still perform the flip to avoid bug of underflow
                self.candidate_snps[*idx].genotype = *gt;
            }
            if check_val == 0 {
                genotype_increase = false;
            } else {
                phasing_increase = true;
                haplotag_increase = true;
                genotype_increase = true;
            }
            self.check_local_optimal_configuration(false, false, true);

            num_iters += 1;
            if num_iters > 20 {
                println!("Warning: iteration inside cross_optimize exceeds 20 times.");
                break;
            }
        }
        let prob = cal_overall_probability(&self);
        return prob;
    }

    fn check_local_optimal_configuration(&self, used_for_haplotype: bool, used_for_haplotag: bool, used_for_genotype: bool) {
        // check sigma
        if used_for_haplotag {
            for k in 0..self.fragments.len() {
                if self.fragments[k].discarded == true { 
                    continue; 
                }
                let sigma_k = self.fragments[k].haplotag;
                let mut delta: Vec<i32> = Vec::new();
                let mut eta: Vec<i32> = Vec::new();
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
                    eta.push(self.candidate_snps[fe.snp_idx].genotype);
                }
                if delta.len() == 0 {
                    continue;
                }
                let q = cal_sigma_delta_eta_log(sigma_k, &delta, &eta, &ps, &probs);
                let qn = cal_sigma_delta_eta_log(sigma_k * (-1), &delta, &eta, &ps, &probs);
                // println!("q:{}, qn:{}", q, qn);
                assert!(q >= qn, "{} Error: haplotag flipping is not local optimal. {}->{}\n{:?}\ndelta:{:?}\nps:{:?}\nprobs:{:?}\nsigma:{}\n{:?}\n{:?}\n{:?}", k, q, qn, self.region, delta, ps, probs, sigma_k, used_for_haplotype, used_for_haplotag, used_for_genotype);
            }
        }

        // check delta
        if used_for_haplotype {
            for i in 0..self.candidate_snps.len() {
                if self.candidate_snps[i].for_phasing == false { continue; }
                let delta_i = self.candidate_snps[i].haplotype;
                let eta_i = self.candidate_snps[i].genotype;
                let coverage_i = self.candidate_snps[i].depth;
                let mut sigma: Vec<i32> = Vec::new();
                let mut ps: Vec<i32> = Vec::new();
                let mut probs: Vec<f64> = Vec::new();
                for k in self.candidate_snps[i].snp_cover_fragments.iter() {
                    if self.fragments[*k].discarded == true {
                        continue;
                    }
                    if self.fragments[*k].haplotag == 0 {
                        continue;
                    }
                    // k is fragment index
                    for fe in self.fragments[*k].list.iter() {
                        if fe.snp_idx == i {
                            if fe.phase_site == false { continue; }
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
                let q = cal_delta_eta_sigma_log(delta_i, eta_i, &sigma, &ps, &probs);
                let qn = cal_delta_eta_sigma_log(delta_i * (-1), eta_i, &sigma, &ps, &probs);
                assert!(q >= qn, "{} Error: haplotype flipping is not local optimal. {}->{}\n{:?}\nsigma:{:?}\nps:{:?}\nprobs:{:?}\n{:?}\n{:?}\n{:?}", i, q, qn, self.region, sigma, ps, probs, used_for_haplotype, used_for_haplotag, used_for_genotype);
            }
        }

        // check eta
        if used_for_genotype {
            for i in 0..self.candidate_snps.len() {
                if self.candidate_snps[i].for_phasing == false { continue; }
                let delta_i = self.candidate_snps[i].haplotype;
                let eta_i = self.candidate_snps[i].genotype;
                let coverage_i = self.candidate_snps[i].depth;
                let mut sigma: Vec<i32> = Vec::new();
                let mut ps: Vec<i32> = Vec::new();
                let mut probs: Vec<f64> = Vec::new();
                for k in self.candidate_snps[i].snp_cover_fragments.iter() {
                    if self.fragments[*k].discarded == true {
                        continue;
                    }
                    if self.fragments[*k].haplotag == 0 {
                        continue;
                    }
                    // k is fragment index
                    for fe in self.fragments[*k].list.iter() {
                        if fe.snp_idx == i {
                            if fe.phase_site == false { continue; }
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
                if eta_i == 0 {
                    let q1 = cal_eta_delta_sigma_log(0, coverage_i, delta_i, &sigma, &ps, &probs);
                    let q2 = cal_eta_delta_sigma_log(1, coverage_i, delta_i, &sigma, &ps, &probs);
                    let q3 = cal_eta_delta_sigma_log(-1, coverage_i, delta_i, &sigma, &ps, &probs);
                    assert!(q1 >= q2 && q1 >= q3, "{} Error: genotype flipping is not local optimal. {}->{},{}\n{:?}\nsigma:{:?}\nps:{:?}\nprobs:{:?}\n{:?}\n{:?}\n{:?}", i, q1, q2, q3, self.region, sigma, ps, probs, used_for_haplotype, used_for_haplotag, used_for_genotype)
                } else if eta_i == 1 {
                    let q1 = cal_eta_delta_sigma_log(1, coverage_i, delta_i, &sigma, &ps, &probs);
                    let q2 = cal_eta_delta_sigma_log(0, coverage_i, delta_i, &sigma, &ps, &probs);
                    let q3 = cal_eta_delta_sigma_log(-1, coverage_i, delta_i, &sigma, &ps, &probs);
                    assert!(q1 >= q2 && q1 >= q3, "{} Error: genotype flipping is not local optimal. {}->{},{}\n{:?}\nsigma:{:?}\nps:{:?}\nprobs:{:?}\n{:?}\n{:?}\n{:?}", i, q1, q2, q3, self.region, sigma, ps, probs, used_for_haplotype, used_for_haplotag, used_for_genotype)
                } else if eta_i == -1 {
                    let q1 = cal_eta_delta_sigma_log(-1, coverage_i, delta_i, &sigma, &ps, &probs);
                    let q2 = cal_eta_delta_sigma_log(0, coverage_i, delta_i, &sigma, &ps, &probs);
                    let q3 = cal_eta_delta_sigma_log(1, coverage_i, delta_i, &sigma, &ps, &probs);
                    assert!(q1 >= q2 && q1 >= q3, "{} Error: genotype flipping is not local optimal. {}->{},{}\n{:?}\nsigma:{:?}\nps:{:?}\nprobs:{:?}\n{:?}\n{:?}\n{:?}", i, q1, q2, q3, self.region, sigma, ps, probs, used_for_haplotype, used_for_haplotag, used_for_genotype)
                }
            }
        }
    }

    fn save_best_configuration(&self, best_haplotype: &mut HashMap<usize, i32>, best_haplotag: &mut HashMap<usize, i32>, best_genotype: &mut HashMap<usize, i32>) {
        best_haplotype.clear();
        best_haplotag.clear();
        best_genotype.clear();
        for i in 0..self.candidate_snps.len() {
            best_haplotype.insert(i, self.candidate_snps[i].haplotype);
            best_genotype.insert(i, self.candidate_snps[i].genotype);
        }
        for k in 0..self.fragments.len() {
            best_haplotag.insert(k, self.fragments[k].haplotag);
        }
    }

    fn load_best_configuration(&mut self, best_haplotype: &HashMap<usize, i32>, best_haplotag: &HashMap<usize, i32>, best_genotype: &HashMap<usize, i32>) {
        for i in 0..self.candidate_snps.len() {
            self.candidate_snps[i].haplotype = best_haplotype.get(&i).unwrap().clone();
            self.candidate_snps[i].genotype = best_genotype.get(&i).unwrap().clone();
        }
        for k in 0..self.fragments.len() {
            self.fragments[k].haplotag = best_haplotag.get(&k).unwrap().clone();
        }
    }


    pub fn chain_phase(&mut self, chain_size: usize) {
        let mut chain_snps: Vec<usize> = Vec::new();
        let mut variant_scores: Vec<(usize, f64)> = Vec::new();
        for i in 0..self.candidate_snps.len() {
            if self.candidate_snps[i].for_phasing == false {
                continue;
            }
            variant_scores.push((i, self.candidate_snps[i].variant_quality));
        }
        // sort by variant quality
        variant_scores.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
        if variant_scores.len() < chain_size {
            for (i, _) in variant_scores.iter() {
                chain_snps.push(*i);
            }
        } else {
            for (i, _) in variant_scores.iter().take(chain_size) {
                chain_snps.push(*i);
            }
        }

        let mut largest_prob = f64::NEG_INFINITY;
        let mut best_haplotype: HashMap<usize, i32> = HashMap::new();
        let mut best_haplotag: HashMap<usize, i32> = HashMap::new();
        let mut best_genotype: HashMap<usize, i32> = HashMap::new();
        let mut haplotype_enum: Vec<Vec<i32>> = Vec::new();
        let init_hap: Vec<i32> = vec![1; chain_snps.len()];
        haplotype_enum.push(init_hap.clone());
        for ti in 0..chain_snps.len() {
            for tj in 0..haplotype_enum.len() {
                let mut tmp_hap = haplotype_enum[tj].clone();
                tmp_hap[ti] = tmp_hap[ti] * (-1);
                haplotype_enum.push(tmp_hap);
            }
        }
        assert!(haplotype_enum.len() == 2_usize.pow(chain_snps.len() as u32), "Error: Not all combinations included");
        for hap in haplotype_enum.iter() {
            for i in 0..chain_snps.len() {
                self.candidate_snps[i].haplotype = hap[i];
            }
            unsafe {
                self.init_assignment();
            }
            let prob = self.cross_optimize_sigma();
            if prob > largest_prob {
                largest_prob = prob;
                self.save_best_configuration(&mut best_haplotype, &mut best_haplotag, &mut best_genotype);
            }
        }
        self.load_best_configuration(&best_haplotype, &best_haplotag, &best_genotype);
    }

    pub fn phase(&mut self, max_enum_snps: usize, random_flip_fraction: f32, max_iters: i32) {
        let mut largest_prob = f64::NEG_INFINITY;
        let mut best_haplotype: HashMap<usize, i32> = HashMap::new();
        let mut best_haplotag: HashMap<usize, i32> = HashMap::new();
        let mut best_genotype: HashMap<usize, i32> = HashMap::new();

        if self.candidate_snps.len() <= max_enum_snps {
            // enumerate the haplotype, then optimize the assignment
            let mut haplotype_enum: Vec<Vec<i32>> = Vec::new();
            let init_hap: Vec<i32> = vec![1; self.candidate_snps.len()];
            haplotype_enum.push(init_hap.clone());
            for ti in 0..self.candidate_snps.len() {
                for tj in 0..haplotype_enum.len() {
                    let mut tmp_hap = haplotype_enum[tj].clone();
                    tmp_hap[ti] = tmp_hap[ti] * (-1);
                    haplotype_enum.push(tmp_hap);
                }
            }
            assert!(haplotype_enum.len() == 2_usize.pow(self.candidate_snps.len() as u32), "Error: Not all combinations included");
            for hap in haplotype_enum.iter() {
                for i in 0..self.candidate_snps.len() {
                    self.candidate_snps[i].haplotype = hap[i];
                }
                unsafe {
                    self.init_assignment();
                    self.init_genotype();
                }
                let prob = self.cross_optimize_sigma();
                if prob > largest_prob {
                    largest_prob = prob;
                    self.save_best_configuration(&mut best_haplotype, &mut best_haplotag, &mut best_genotype);
                }
            }
            self.load_best_configuration(&best_haplotype, &best_haplotag, &best_genotype);
        } else {
            // optimize haplotype and read assignment alternatively
            let mut max_iter: i32 = max_iters;
            // random initialization of haplotype and haplotag at each iteration
            unsafe {
                self.init_haplotypes();
                self.init_genotype();
            }
            unsafe {
                self.init_assignment();
            }
            let prob = self.cross_optimize_sigma();
            if prob > largest_prob {
                largest_prob = prob;
                self.save_best_configuration(&mut best_haplotype, &mut best_haplotag, &mut best_genotype);
            }
            self.load_best_configuration(&best_haplotype, &best_haplotag, &best_genotype);
            // for _ in 0..1 {
            //     {
            //         let mut rng = rand::thread_rng();
            //         for ti in 0..self.candidate_snps.len() {
            //             let rg: f64 = rng.gen();
            //             if rg < 0.33 {
            //                 self.candidate_snps[ti].haplotype = -1;
            //             } else if rg < 0.66 {
            //                 self.candidate_snps[ti].haplotype = 1;
            //             } else if rg < 1.0 {
            //                 self.candidate_snps[ti].haplotype = 0;
            //             }
            //         }
            //         for tk in 0..self.fragments.len() {
            //             if self.fragments[tk].discarded == true {
            //                 continue;
            //             }
            //             if self.fragments[tk].haplotag == 0 {
            //                 continue;
            //             }
            //             let rg: f64 = rng.gen();
            //             if rg < 0.5 as f64 {
            //                 self.fragments[tk].haplotag = self.fragments[tk].haplotag * (-1);
            //             }
            //         }
            //     }
            //     let prob = self.cross_optimize();
            //     if prob > largest_prob {
            //         largest_prob = prob;
            //         self.save_best_configuration(&mut best_haplotype, &mut best_haplotag);
            //     }
            //     self.load_best_configuration(&best_haplotype, &best_haplotag);
            // }


            // for _ in 0..1 {
            //     {
            //         let mut rng = rand::thread_rng();
            //         for tk in 0..self.fragments.len() {
            //             if self.fragments[tk].discarded == true {
            //                 continue;
            //             }
            //             if self.fragments[tk].haplotag == 0 {
            //                 continue;
            //             }
            //             let rg: f64 = rng.gen();
            //             if rg < 0.2 as f64 {
            //                 self.fragments[tk].haplotag = self.fragments[tk].haplotag * (-1);
            //             }
            //         }
            //     }
            //     let prob = self.cross_optimize_delta();
            //     if prob > largest_prob {
            //         largest_prob = prob;
            //         self.save_best_configuration(&mut best_haplotype, &mut best_haplotag);
            //     }
            //     self.load_best_configuration(&best_haplotype, &best_haplotag);
            // }

            for tidx in 0..self.candidate_snps.len() {
                {
                    for ti in 0..self.candidate_snps.len() {
                        let mut rng = rand::thread_rng();
                        let rg: f64 = rng.gen();
                        if tidx % 2 == 0 {
                            if rg < 0.1 {
                                self.candidate_snps[ti].haplotype = -1;
                            } else if rg >= 0.9 {
                                self.candidate_snps[ti].haplotype = 1;
                            }
                        } else if tidx % 2 == 1 {
                            if rg < 0.1 {
                                self.candidate_snps[ti].haplotype = 1;
                            } else if rg >= 0.9 {
                                self.candidate_snps[ti].haplotype = -1;
                            }
                        }
                    }
                }
                let prob = self.cross_optimize_sigma();
                if prob > largest_prob {
                    largest_prob = prob;
                    self.save_best_configuration(&mut best_haplotype, &mut best_haplotag, &mut best_genotype);
                }
                self.load_best_configuration(&best_haplotype, &best_haplotag, &best_genotype);

                {
                    let mut rng = rand::thread_rng();
                    for tk in 0..self.fragments.len() {
                        if self.fragments[tk].discarded == true {
                            continue;
                        }
                        if self.fragments[tk].haplotag == 0 {
                            continue;
                        }
                        let rg: f64 = rng.gen();
                        if rg < 0.1 as f64 {
                            self.fragments[tk].haplotag = self.fragments[tk].haplotag * (-1);
                        }
                    }
                }
                let prob = self.cross_optimize_sigma();
                if prob > largest_prob {
                    largest_prob = prob;
                    self.save_best_configuration(&mut best_haplotype, &mut best_haplotag, &mut best_genotype);
                }
                self.load_best_configuration(&best_haplotype, &best_haplotag, &best_genotype);
            }


            // // when initial setting has reached to local optimal, flip all the snps after a specific position to jump out local optimization
            // let mut unflipped_haplotype: Vec<i32> = Vec::new();
            // for i in 0..self.candidate_snps.len() {
            //     unflipped_haplotype.push(self.candidate_snps[i].haplotype);
            // }
            // for ti in 0..unflipped_haplotype.len() {
            //     let mut tmp_hap: Vec<i32> = Vec::new();
            //     for tj in 0..unflipped_haplotype.len() {
            //         if tj < ti {
            //             tmp_hap.push(unflipped_haplotype[tj]);
            //         } else {
            //             // TODO: only flip between 1 and -1, no jump out 0 state or jump to 0 state
            //             tmp_hap.push(unflipped_haplotype[tj] * (-1));
            //         }
            //     }
            //     // block flip
            //     {
            //         assert_eq!(tmp_hap.len(), self.candidate_snps.len());
            //         for i in 0..self.candidate_snps.len() {
            //             self.candidate_snps[i].haplotype = tmp_hap[i];
            //         }
            //     }
            //     let prob = self.cross_optimize();
            //     if prob > largest_prob {
            //         largest_prob = prob;
            //         self.save_best_configuration(&mut best_haplotype, &mut best_haplotag);
            //     }
            //     self.load_best_configuration(&best_haplotype, &best_haplotag);
            //
            //     // when current block flip has reached to local optimal, flip a fraction of snps and reads to jump out local optimization
            //     {
            //         let mut rng = rand::thread_rng();
            //         for ti in 0..self.candidate_snps.len() {
            //             let rg: f64 = rng.gen();
            //             if rg < (random_flip_fraction as f64) / 2.0 {
            //                 self.candidate_snps[ti].haplotype = self.candidate_snps[ti].haplotype * (-1);
            //             } else if rg < (random_flip_fraction as f64) {
            //                 self.candidate_snps[ti].haplotype = self.candidate_snps[ti].haplotype * 0;
            //             }
            //         }
            //         for tk in 0..self.fragments.len() {
            //             if self.fragments[tk].discarded == true {
            //                 continue;
            //             }
            //             if self.fragments[tk].haplotag == 0 {
            //                 continue;
            //             }
            //             let rg: f64 = rng.gen();
            //             if rg < random_flip_fraction as f64 {
            //                 self.fragments[tk].haplotag = self.fragments[tk].haplotag * (-1);
            //             }
            //         }
            //     }
            //     let prob = self.cross_optimize();
            //     if prob > largest_prob {
            //         largest_prob = prob;
            //         self.save_best_configuration(&mut best_haplotype, &mut best_haplotag);
            //     }
            //     self.load_best_configuration(&best_haplotype, &best_haplotag);
            // }
            self.load_best_configuration(&best_haplotype, &best_haplotag, &best_genotype);
        }
        // println!("largest_prob: {}", largest_prob);
    }
}