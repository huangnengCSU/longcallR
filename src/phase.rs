use std::collections::{HashMap, HashSet};

use itertools::Itertools;
use petgraph::algo::kosaraju_scc;
use petgraph::graphmap::GraphMap;
use petgraph::Undirected;
use petgraph::visit::Bfs;
use rand::Rng;
use rand::{SeedableRng, seq::SliceRandom};
use rand::rngs::StdRng;

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
    let mut log_q4: f64 = 0.0;
    let mut log_q5: f64 = 0.0;

    let prior_homref_log: f64 = (1.0 - 1.5 * 0.001 as f64).log10();
    let prior_homvar_log: f64 = (0.5 * 0.001 as f64).log10();
    let prior_hetvar_log: f64;
    let coverage_i = sigma.len() as u32;
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
        log_q5 += aki(sigma[k], delta_i * (-1), 0, ps[k], probs[k]).log10();
    }

    log_q2 += prior_homvar_log;
    log_q3 += prior_hetvar_log;
    log_q4 += prior_homref_log;
    log_q5 += prior_hetvar_log;

    /*
    Given 0<=A<=1, 0<=B<=1 and A/(A+B) > B/(A+B), we have logA/(logA+logB) < logB/(logA+logB).
    When comparing A/(A+B) and B/(A+B), we can approximately use 1-log(A)/(log(A)+log(B)) and 1-log(B)/(log(A)+log(B)) to avoid underflow.
    */
    return 1.0 - log_q1 / (log_q2 + log_q3 + log_q4 + log_q5);
}

pub fn cal_block_delta_eta_sigma_log(block_delta: &Vec<i32>, block_eta: &Vec<i32>, block_sigma: &Vec<Vec<i32>>, block_ps: &Vec<Vec<i32>>, block_probs: &Vec<Vec<f64>>) -> f64 {
    // same as call_delta_sigma, but use log to avoid underflow
    let mut sum_prob: f64 = 0.0;
    for i in 0..block_delta.len() {
        let delta_i = block_delta[i];
        let eta_i = block_eta[i];
        let sigma = &block_sigma[i];
        let ps = &block_ps[i];
        let probs = &block_probs[i];

        let mut log_q1: f64 = 0.0;
        let mut log_q2: f64 = 0.0;
        let mut log_q3: f64 = 0.0;
        let mut log_q4: f64 = 0.0;
        let mut log_q5: f64 = 0.0;

        let prior_homref_log: f64 = (1.0 - 1.5 * 0.001 as f64).log10();
        let prior_homvar_log: f64 = (0.5 * 0.001 as f64).log10();
        let prior_hetvar_log: f64;
        let coverage_i = sigma.len() as u32;
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
            log_q5 += aki(sigma[k], delta_i * (-1), 0, ps[k], probs[k]).log10();
        }

        log_q2 += prior_homvar_log;
        log_q3 += prior_hetvar_log;
        log_q4 += prior_homref_log;
        log_q5 += prior_hetvar_log;

        /*
        Given 0<=A<=1, 0<=B<=1 and A/(A+B) > B/(A+B), we have logA/(logA+logB) < logB/(logA+logB).
        When comparing A/(A+B) and B/(A+B), we can approximately use 1-log(A)/(log(A)+log(B)) and 1-log(B)/(log(A)+log(B)) to avoid underflow.
        */
        sum_prob += 1.0 - log_q1 / (log_q2 + log_q3 + log_q4 + log_q5);
    }
    return sum_prob;
}

pub fn cal_phase_score_log(delta_i: i32, eta_i: i32, sigma: &Vec<i32>, ps: &Vec<i32>, probs: &Vec<f64>) -> f64 {
    // same as call_delta_sigma, but use log to avoid underflow
    let mut log_q1: f64 = 0.0;
    let mut log_q2: f64 = 0.0;
    let mut log_q3: f64 = 0.0;
    assert!(delta_i != 0, "Error: phase for unexpected allele.");
    assert!(eta_i == 0, "Error: phase for unexpected genotype.");

    for k in 0..sigma.len() {
        log_q1 += aki(sigma[k], delta_i, eta_i, ps[k], probs[k]).log10();
    }

    for k in 0..sigma.len() {
        log_q2 += aki(sigma[k], 1, eta_i, ps[k], probs[k]).log10();
        log_q3 += aki(sigma[k], -1, eta_i, ps[k], probs[k]).log10();
    }
    return 1.0 - log_q1 / (log_q2 + log_q3);
}

pub fn cal_overall_probability(snpfrag: &SNPFrag, apply_downsampling: bool) -> f64 {
    // calculate the log10 probability of the current configuration of sigma and delta
    let mut logp = 0.0;
    for k in 0..snpfrag.fragments.len() {
        if !snpfrag.fragments[k].for_phasing 
            || (apply_downsampling && !snpfrag.fragments[k].downsampled)
            || snpfrag.fragments[k].haplotag == 0 {
            // unassigned fragment
            continue;
        }
        for fe in snpfrag.fragments[k].list.iter() {
            if !fe.phase_site {
                continue;
            }
            assert_ne!(fe.p, 0, "Error: phasing with unexpected hete SNP.");
            logp += aki(snpfrag.fragments[k].haplotag, snpfrag.candidate_snps[fe.snp_idx].haplotype, snpfrag.candidate_snps[fe.snp_idx].genotype, fe.p, fe.prob).log10();
        }
    }
    return logp;
}

pub fn check_new_haplotag(snpfrag: &SNPFrag, updated_haplotag: &HashMap<usize, i32>, apply_downsampling: bool) -> i32 {
    // updated_haplotag: the index of the fragments will be updated
    let mut logp = 0.0;
    let mut pre_logp = 0.0;
    for (k, h) in updated_haplotag.iter() {
        let mut delta: Vec<i32> = Vec::new();
        let mut eta: Vec<i32> = Vec::new();
        let mut ps: Vec<i32> = Vec::new();
        let mut probs: Vec<f64> = Vec::new();
        if snpfrag.fragments[*k].haplotag == 0 
            || (apply_downsampling && !snpfrag.fragments[*k].downsampled) {
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
    let flag = match logp.partial_cmp(&pre_logp).unwrap() {
        std::cmp::Ordering::Greater => 1,
        std::cmp::Ordering::Equal => 0,
        std::cmp::Ordering::Less => -1,
    };
    assert!(flag >= 0, "Error: new haplotag decrease the probability. {} -> {}", pre_logp, logp);
    flag
}

pub fn check_new_haplotype_genotype(snpfrag: &SNPFrag, updated_haplotype_genotype: &HashMap<usize, (i32, i32)>, apply_downsampling: bool) -> i32 {
    let mut logp = 0.0;
    let mut pre_logp = 0.0;
    for (i, h) in updated_haplotype_genotype.iter() {
        let delta_i = h.0;
        let eta_i = h.1;
        let mut sigma: Vec<i32> = Vec::new();
        let mut ps: Vec<i32> = Vec::new();
        let mut probs: Vec<f64> = Vec::new();
        for k in snpfrag.candidate_snps[*i].snp_cover_fragments.iter() {
            if !snpfrag.fragments[*k].for_phasing
                || (apply_downsampling && !snpfrag.fragments[*k].downsampled)
                || snpfrag.fragments[*k].haplotag == 0 {
                continue;
            }
            for fe in snpfrag.fragments[*k].list.iter() {
                if fe.snp_idx != *i || !fe.phase_site {
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
        logp += cal_delta_eta_sigma_log(delta_i, eta_i, &sigma, &ps, &probs);
        pre_logp += cal_delta_eta_sigma_log(snpfrag.candidate_snps[*i].haplotype, snpfrag.candidate_snps[*i].genotype, &sigma, &ps, &probs);
    }
    let flag = match logp.partial_cmp(&pre_logp).unwrap() {
        std::cmp::Ordering::Greater => 1,
        std::cmp::Ordering::Equal => 0,
        std::cmp::Ordering::Less => -1,
    };

    assert!(flag >= 0, "Error: new haplotype decrease the probability. {} -> {}", pre_logp, logp);
    return flag;
}

pub fn check_block_new_haplotype_haplotag(snpfrag: &SNPFrag, updated_haplotype: &HashMap<usize, i32>, updated_haplotag: &HashMap<usize, i32>, apply_downsampling: bool) -> i32 {
    let mut logp = 0.0;
    let mut pre_logp = 0.0;
    for (i, h) in updated_haplotype.iter() {
        let delta_i = *h;
        let eta_i = snpfrag.candidate_snps[*i].genotype;
        let mut sigma: Vec<i32> = Vec::new();
        let mut tmp_sigma: Vec<i32> = Vec::new();
        let mut ps: Vec<i32> = Vec::new();
        let mut probs: Vec<f64> = Vec::new();
        for k in snpfrag.candidate_snps[*i].snp_cover_fragments.iter() {
            if snpfrag.fragments[*k].haplotag == 0 
                || (apply_downsampling && !snpfrag.fragments[*k].downsampled){
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
                tmp_sigma.push(updated_haplotag[k]);
            }
        }
        if sigma.len() == 0 || tmp_sigma.len() == 0 {
            continue;
        }
        logp += cal_delta_eta_sigma_log(delta_i, eta_i, &tmp_sigma, &ps, &probs);
        pre_logp += cal_delta_eta_sigma_log(snpfrag.candidate_snps[*i].haplotype, eta_i, &sigma, &ps, &probs);
    }
    let flag = match logp.partial_cmp(&pre_logp).unwrap() {
        std::cmp::Ordering::Greater => 1,
        std::cmp::Ordering::Equal => 0,
        std::cmp::Ordering::Less => -1,
    };
    assert!(flag >= 0, "Error: new haplotype decrease the probability. {} -> {}", pre_logp, logp);
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

pub fn get_blocks(block_links: &HashMap<[usize; 2], i32>, weight_threshold: u32) -> (Vec<Vec<usize>>, GraphMap<usize, i32, Undirected>) {
    let mut graph: GraphMap<usize, i32, Undirected> = GraphMap::new();  // node is index in candidate snp, edge is index in fragments
    for (edge, w) in block_links.iter() {
        if graph.contains_edge(edge[0], edge[1]) {
            let weight = graph.edge_weight_mut(edge[0], edge[1]).unwrap();
            *weight += *w;
        } else {
            graph.add_edge(edge[0], edge[1], *w);
        }
    }
    // delete edges with weight < weight_threshold
    let mut lw_edges: Vec<(usize, usize)> = Vec::new();
    for edge in graph.all_edges() {
        if (edge.2.abs() as u32) < weight_threshold {
            lw_edges.push((edge.0, edge.1));
        }
    }
    for edge in lw_edges.iter() {
        graph.remove_edge(edge.0, edge.1);
    }
    let scc = kosaraju_scc(&graph);
    // for block in scc.iter() {
    //     let mut bfs = Bfs::new(&graph, block[0]);
    //     print!("{}->", block[0]);
    //     while let Some(nx) = bfs.next(&graph) {
    //         print!("{},", nx);
    //     }
    // }
    return (scc, graph);
}

impl SNPFrag {
    pub unsafe fn init_haplotypes(&mut self) {
        let mut rng = rand::thread_rng();
        for snp in &mut self.candidate_snps {
            snp.haplotype = if rng.gen::<f64>() < 0.5 { 1 } else { -1 };
        }
    }

    // pub unsafe fn init_haplotypes_LD(&mut self, ld_weight_threshold: u32) -> (HashSet<usize>, Vec<HashSet<usize>>) {
    //     for ti in 0..self.candidate_snps.len() {
    //         let mut rng = rand::thread_rng();
    //         let rg: f64 = rng.gen();
    //         if rg < 0.5 {
    //             self.candidate_snps[ti].haplotype = 1; // HetVar hap1
    //         } else {
    //             self.candidate_snps[ti].haplotype = -1;  // HetVar hap2
    //         }
    //     }
    //
    //     let mut ld_idxes: Vec<usize> = Vec::new();
    //     for ti in 0..self.candidate_snps.len() {
    //         ld_idxes.push(ti);
    //     }
    //
    //     // let mut r2_map = HashMap::new();
    //     let mut block_links: HashMap<[usize; 2], i32> = HashMap::new();
    //     for i in 0..ld_idxes.len() {
    //         for j in i + 1..ld_idxes.len() {
    //             if i == j { continue; }
    //             let idx1 = ld_idxes[i];
    //             let idx2 = ld_idxes[j];
    //             let snp1 = &self.candidate_snps[idx1];
    //             let snp2 = &self.candidate_snps[idx2];
    //             let (mut snp1_ref, mut snp1_alt, mut snp2_ref, mut snp2_alt);
    //             let (mut snp1_ref_frac, mut snp1_alt_frac, mut snp2_ref_frac, mut snp2_alt_frac) = (0.0, 0.0, 0.0, 0.0);
    //             if snp1.alleles[0] == snp1.reference && snp1.alleles[1] != snp1.reference {
    //                 snp1_ref = snp1.alleles[0] as u8;
    //                 snp1_ref_frac = snp1.allele_freqs[0];
    //                 snp1_alt = snp1.alleles[1] as u8;
    //                 snp1_alt_frac = snp1.allele_freqs[1];
    //             } else if snp1.alleles[0] != snp1.reference && snp1.alleles[1] == snp1.reference {
    //                 snp1_ref = snp1.alleles[1] as u8;
    //                 snp1_ref_frac = snp1.allele_freqs[1];
    //                 snp1_alt = snp1.alleles[0] as u8;
    //                 snp1_alt_frac = snp1.allele_freqs[0];
    //             } else {
    //                 continue;
    //             }
    //             if snp2.alleles[0] == snp2.reference && snp2.alleles[1] != snp2.reference {
    //                 snp2_ref = snp2.alleles[0] as u8;
    //                 snp2_ref_frac = snp2.allele_freqs[0];
    //                 snp2_alt = snp2.alleles[1] as u8;
    //                 snp2_alt_frac = snp2.allele_freqs[1];
    //             } else if snp2.alleles[0] != snp2.reference && snp2.alleles[1] == snp2.reference {
    //                 snp2_ref = snp2.alleles[1] as u8;
    //                 snp2_ref_frac = snp2.allele_freqs[1];
    //                 snp2_alt = snp2.alleles[0] as u8;
    //                 snp2_alt_frac = snp2.allele_freqs[0];
    //             } else {
    //                 continue;
    //             }
    //             if idx1 < idx2 {
    //                 if !self.allele_pairs.contains_key(&[idx1, idx2]) { continue; }
    //                 // let r2 = self.allele_pairs.get(&[idx1, idx2]).unwrap().calculate_LD_R2(snp1.alleles[0] as u8, snp1.alleles[1] as u8, snp2.alleles[0] as u8, snp2.alleles[1] as u8);
    //                 if snp1_ref_frac == 0.0 || snp1_alt_frac == 0.0 || snp2_ref_frac == 0.0 || snp2_alt_frac == 0.0 { continue; }
    //                 let (is_block, score, weight) = self.allele_pairs.get(&[idx1, idx2]).unwrap().is_same_block(snp1_ref, snp1_alt, snp2_ref, snp2_alt);
    //                 if is_block {
    //                     block_links.insert([idx1, idx2], weight);
    //                 }
    //                 // r2_map.insert([idx1, idx2], r2);
    //             } else if idx2 < idx1 {
    //                 if !self.allele_pairs.contains_key(&[idx2, idx1]) { continue; }
    //                 // let r2 = self.allele_pairs.get(&[idx2, idx1]).unwrap().calculate_LD_R2(snp2.alleles[0] as u8, snp2.alleles[1] as u8, snp1.alleles[0] as u8, snp1.alleles[1] as u8);
    //                 if snp1_ref_frac == 0.0 || snp1_alt_frac == 0.0 || snp2_ref_frac == 0.0 || snp2_alt_frac == 0.0 { continue; }
    //                 let (is_block, score, weight) = self.allele_pairs.get(&[idx2, idx1]).unwrap().is_same_block(snp2_ref, snp2_alt, snp1_ref, snp1_alt);
    //                 if is_block {
    //                     block_links.insert([idx2, idx1], weight);
    //                 }
    //                 // r2_map.insert([idx2, idx1], r2);
    //             }
    //         }
    //     }
    //
    //     // let r2_map_high: HashMap<[usize; 2], f32> = r2_map.iter().filter(|(_, &v)| v > 0.9).map(|(&k, &v)| (k, v)).collect();
    //     // let mut ld_sites: HashSet<usize> = HashSet::new();
    //     // for idxes in r2_map_high.keys() {
    //     //     ld_sites.insert(idxes[0]);
    //     //     ld_sites.insert(idxes[1]);
    //     // }
    //     // for ti in 0..self.candidate_snps.len() {
    //     //     if ld_sites.contains(&ti) {
    //     //         self.candidate_snps[ti].haplotype = 1;
    //     //     } else {
    //     //         let mut rng = rand::thread_rng();
    //     //         let rg: f64 = rng.gen();
    //     //         if rg < 0.5 {
    //     //             self.candidate_snps[ti].haplotype = 1; // HetVar hap1
    //     //         } else {
    //     //             self.candidate_snps[ti].haplotype = -1;  // HetVar hap2
    //     //         }
    //     //     }
    //     // }
    //     let (blocks_vec, block_graph) = get_blocks(&block_links, ld_weight_threshold);
    //     let mut block_snps: HashSet<usize> = HashSet::new();
    //     let mut blocks_set: Vec<HashSet<usize>> = Vec::new();
    //     for block in blocks_vec.iter() {
    //         if block.len() < 2 { continue; }
    //         // use BFS to assign snps in the phased block
    //         let mut bfs = Bfs::new(&block_graph, block[0]);
    //         let mut visited_nodes: Vec<usize> = Vec::new();
    //         self.candidate_snps[block[0]].haplotype = 1;    // set the first snp as hap1
    //         visited_nodes.push(block[0]);
    //         while let Some(nx) = bfs.next(&block_graph) {
    //             for idx in visited_nodes.iter() {
    //                 if *idx < nx {
    //                     let Some(weight) = block_links.get(&[*idx, nx]) else { continue; };
    //                     if *weight >= ld_weight_threshold as i32 {
    //                         self.candidate_snps[nx].haplotype = 1 * self.candidate_snps[*idx].haplotype;  // positive value, same haplotype
    //                         break;
    //                     } else if *weight <= -1 * ld_weight_threshold as i32 {
    //                         self.candidate_snps[nx].haplotype = -1 * self.candidate_snps[*idx].haplotype; // negative value, different haplotype
    //                         break;
    //                     }
    //                     // if block_links.contains_key(&[*idx, nx]) {
    //                     //     if *block_links.get(&[*idx, nx]).unwrap() >= ld_weight_threshold as i32 {
    //                     //         self.candidate_snps[nx].haplotype = 1 * self.candidate_snps[*idx].haplotype;
    //                     //         break;
    //                     //     } else if *block_links.get(&[*idx, nx]).unwrap() <= -1 * ld_weight_threshold as i32 {
    //                     //         self.candidate_snps[nx].haplotype = -1 * self.candidate_snps[*idx].haplotype;
    //                     //         break;
    //                     //     }
    //                     // }
    //                 } else if *idx > nx {
    //                     let Some(weight) = block_links.get(&[nx, *idx]) else { continue; };
    //                     if *weight >= ld_weight_threshold as i32 {
    //                         self.candidate_snps[nx].haplotype = 1 * self.candidate_snps[*idx].haplotype;  // positive value, same haplotype
    //                         break;
    //                     } else if *weight <= -1 * ld_weight_threshold as i32 {
    //                         self.candidate_snps[nx].haplotype = -1 * self.candidate_snps[*idx].haplotype; // negative value, different haplotype
    //                         break;
    //                     }
    //                     // if block_links.contains_key(&[nx, *idx]) {
    //                     //     if *block_links.get(&[nx, *idx]).unwrap() >= ld_weight_threshold as i32 {
    //                     //         self.candidate_snps[nx].haplotype = 1 * self.candidate_snps[*idx].haplotype;
    //                     //         break;
    //                     //     } else if *block_links.get(&[nx, *idx]).unwrap() <= -1 * ld_weight_threshold as i32 {
    //                     //         self.candidate_snps[nx].haplotype = -1 * self.candidate_snps[*idx].haplotype;
    //                     //         break;
    //                     //     }
    //                     // }
    //                 }
    //             }
    //             visited_nodes.push(nx);
    //         }
    //         let mut rng = rand::thread_rng();
    //         let rg: f64 = rng.gen();
    //         // print!("block:");
    //         blocks_set.push(HashSet::from_iter(block.iter().cloned()));
    //         for idx in block.iter().sorted() {
    //             // print!("({},{})", self.candidate_snps[*idx].pos, self.candidate_snps[*idx].haplotype);
    //             block_snps.insert(*idx);
    //         }
    //         // println!();
    //     }
    //     return (block_snps, blocks_set);
    // }

    pub unsafe fn init_haplotypes_LD2(&mut self, ld_graph: &GraphMap<usize, i32, Undirected>, ld_weight_threshold: u32) -> HashSet<usize> {
        // randomly assign haplotypes for candidate snps
        let mut rng = rand::thread_rng();
        for candidate_snp in &mut self.candidate_snps {
            candidate_snp.haplotype = if rng.gen::<f64>() < 0.5 { 1 } else { -1 };
        }

        // use BFS to assign haplotype of snps in each block
        let mut conserved_snps: HashSet<usize> = HashSet::new();    // snps in the block with more than 2 snps
        // println!("{}", format!("{}", Dot::new(&ld_graph)));
        for block in self.ld_blocks.iter() {
            if block.snp_idxes.len() < 2 { continue; }
            // use BFS to assign snps in the phased block
            let mut bfs = Bfs::new(&ld_graph, block.snp_idxes[0]);
            let mut visited_nodes: Vec<usize> = Vec::new();
            self.candidate_snps[block.snp_idxes[0]].haplotype = 1;    // set the first snp as hap1
            visited_nodes.push(block.snp_idxes[0]);
            while let Some(nx) = bfs.next(&ld_graph) {
                for visited_idx in visited_nodes.iter() {
                    let (from_idx, to_idx);
                    if *visited_idx < nx {
                        from_idx = *visited_idx;
                        to_idx = nx;
                    } else if *visited_idx > nx {
                        from_idx = nx;
                        to_idx = *visited_idx;
                    } else {
                        continue;
                    }
                    let mut weight = 0;
                    if self.allele_pairs.contains_key(&[from_idx, to_idx]) {
                        if !self.allele_pairs.get(&[from_idx, to_idx]).unwrap().valid {
                            continue;   // Score and weight are merely initial values when declared
                        }
                        if self.allele_pairs.get(&[from_idx, to_idx]).unwrap().score != 0.0 {
                            continue;   // not perfect LD
                        }
                        weight = self.allele_pairs.get(&[from_idx, to_idx]).unwrap().weight;
                    } else {
                        continue;
                    }
                    if weight >= ld_weight_threshold as i32 {
                        self.candidate_snps[nx].haplotype = 1 * self.candidate_snps[*visited_idx].haplotype;  // positive value, same haplotype
                        break;
                    } else if weight <= -1 * ld_weight_threshold as i32 {
                        self.candidate_snps[nx].haplotype = -1 * self.candidate_snps[*visited_idx].haplotype; // negative value, different haplotype
                        break;
                    }
                }
                visited_nodes.push(nx);
            }
            // print!("block:");
            if block.snp_idxes.len() >= 2 {
                // blocks with only one snp are not conserved
                for idx in block.snp_idxes.iter() {
                    // print!("({},{})", self.candidate_snps[*idx].pos, self.candidate_snps[*idx].haplotype);
                    conserved_snps.insert(*idx);
                }
            }
            // println!();
        }
        return conserved_snps;
    }

    pub unsafe fn init_assignment(&mut self) {
        let mut rng = rand::thread_rng();
        self.fragments.iter_mut()
            .filter(|fragment| fragment.for_phasing)
            .for_each(|fragment| {
                fragment.haplotag = if rng.gen::<f64>() < 0.5 { -1 } else { 1 };
            });
    }

    pub fn init_genotype(&mut self) {
        for snp in &mut self.candidate_snps {
            snp.genotype = match snp.variant_type {
                0 => 1,  // HomRef
                1 => 0,  // HetVar
                2 | 3 => -1,  // HomVar or triallelic SNP
                _ => snp.genotype,  // Default case, if needed
            };
        }
    }

    pub unsafe fn downsample_fragments(&mut self, downsample_depth: u32, seed: u64) {
        let mut rng = StdRng::seed_from_u64(seed);
        let mut shuffled_indices = (0..self.fragments.len()).collect::<Vec<usize>>();
        shuffled_indices.shuffle(&mut rng);
        let selected_indices = shuffled_indices.into_iter().take(downsample_depth as usize).collect::<Vec<usize>>();
        for i in selected_indices.iter() {
            self.fragments[*i].downsampled = true;
        }
    }

    // pub fn get_ld_blocks(&mut self, ld_weight_threshold: u32) -> (HashMap<[usize; 2], (i32, f32)>, GraphMap<usize, i32, Undirected>) {
    //     let mut ld_idxes: Vec<usize> = Vec::new();
    //     for ti in 0..self.candidate_snps.len() {
    //         if self.candidate_snps[ti].for_phasing {
    //             ld_idxes.push(ti);
    //         }
    //     }
    //
    //     let mut ld_links: HashMap<[usize; 2], (i32, f32)> = HashMap::new();
    //     for i in 0..ld_idxes.len() {
    //         for j in i + 1..ld_idxes.len() {
    //             if i == j { continue; }
    //             let idx1 = ld_idxes[i];
    //             let idx2 = ld_idxes[j];
    //             let snp1 = &self.candidate_snps[idx1];
    //             let snp2 = &self.candidate_snps[idx2];
    //             let (mut snp1_ref, mut snp1_alt, mut snp2_ref, mut snp2_alt);
    //             let (mut snp1_ref_frac, mut snp1_alt_frac, mut snp2_ref_frac, mut snp2_alt_frac) = (0.0, 0.0, 0.0, 0.0);
    //             if snp1.alleles[0] == snp1.reference && snp1.alleles[1] != snp1.reference {
    //                 snp1_ref = snp1.alleles[0] as u8;
    //                 snp1_ref_frac = snp1.allele_freqs[0];
    //                 snp1_alt = snp1.alleles[1] as u8;
    //                 snp1_alt_frac = snp1.allele_freqs[1];
    //             } else if snp1.alleles[0] != snp1.reference && snp1.alleles[1] == snp1.reference {
    //                 snp1_ref = snp1.alleles[1] as u8;
    //                 snp1_ref_frac = snp1.allele_freqs[1];
    //                 snp1_alt = snp1.alleles[0] as u8;
    //                 snp1_alt_frac = snp1.allele_freqs[0];
    //             } else {
    //                 continue;
    //             }
    //             if snp2.alleles[0] == snp2.reference && snp2.alleles[1] != snp2.reference {
    //                 snp2_ref = snp2.alleles[0] as u8;
    //                 snp2_ref_frac = snp2.allele_freqs[0];
    //                 snp2_alt = snp2.alleles[1] as u8;
    //                 snp2_alt_frac = snp2.allele_freqs[1];
    //             } else if snp2.alleles[0] != snp2.reference && snp2.alleles[1] == snp2.reference {
    //                 snp2_ref = snp2.alleles[1] as u8;
    //                 snp2_ref_frac = snp2.allele_freqs[1];
    //                 snp2_alt = snp2.alleles[0] as u8;
    //                 snp2_alt_frac = snp2.allele_freqs[0];
    //             } else {
    //                 continue;
    //             }
    //             assert!(idx1 < idx2, "Error: unexpected index order.");
    //             if !self.allele_pairs.contains_key(&[idx1, idx2]) { continue; }
    //             if snp1_ref_frac == 0.0 || snp1_alt_frac == 0.0 || snp2_ref_frac == 0.0 || snp2_alt_frac == 0.0 { continue; }
    //             // same_block is whether two snps belongs to the same block.
    //             // score is like the ratio of reads conflict the main LD pair.
    //             // weight is like the number of reads supporting the LD.
    //             let (same_block, score, weight) = self.allele_pairs.get(&[idx1, idx2]).unwrap().is_same_block(snp1_ref, snp1_alt, snp2_ref, snp2_alt);
    //             if same_block {
    //                 ld_links.insert([idx1, idx2], (weight, score));
    //             }
    //         }
    //     }
    //
    //     // construct the graph to find the connected components, which are the LD blocks
    //     let mut ld_graph: GraphMap<usize, i32, Undirected> = GraphMap::new();  // node is index in candidate snp, edge is index in fragments
    //     for (edge, (weight, score)) in ld_links.iter() {
    //         if ld_graph.contains_edge(edge[0], edge[1]) {
    //             let w = ld_graph.edge_weight_mut(edge[0], edge[1]).unwrap();
    //             *w += weight;
    //         } else {
    //             ld_graph.add_edge(edge[0], edge[1], *weight);
    //         }
    //     }
    //
    //     // delete edges with weight < weight_threshold
    //     let mut low_w_edges: Vec<(usize, usize)> = Vec::new();
    //     for edge in ld_graph.all_edges() {
    //         if (edge.2.abs() as u32) < ld_weight_threshold {
    //             low_w_edges.push((edge.0, edge.1));
    //         }
    //     }
    //     for edge in low_w_edges.iter() {
    //         ld_graph.remove_edge(edge.0, edge.1);
    //     }
    //
    //     // visualize the LD graph, node is position of snp, edge is the weight of LD
    //     // let mut vis_graph: GraphMap<i64, i32, Undirected> = GraphMap::new();
    //     // for (edge, (weight, score)) in ld_links.iter() {
    //     //     if *score != 0.0 { continue; }  // not perfect LD, make sure each edge is perfect LD
    //     //     if vis_graph.contains_edge(self.candidate_snps[edge[0]].pos, self.candidate_snps[edge[1]].pos) {
    //     //         let w = vis_graph.edge_weight_mut(self.candidate_snps[edge[0]].pos, self.candidate_snps[edge[1]].pos).unwrap();
    //     //         *w += weight;
    //     //     } else {
    //     //         vis_graph.add_edge(self.candidate_snps[edge[0]].pos, self.candidate_snps[edge[1]].pos, *weight);
    //     //     }
    //     // }
    //     // let mut low_w_edges2: Vec<(i64, i64)> = Vec::new();
    //     // for edge in vis_graph.all_edges() {
    //     //     if (edge.2.abs() as u32) < ld_weight_threshold {
    //     //         low_w_edges2.push((edge.0, edge.1));
    //     //     }
    //     // }
    //     // for edge in low_w_edges2.iter() {
    //     //     vis_graph.remove_edge(edge.0, edge.1);
    //     // }
    //     // println!("{}", format!("{}", Dot::new(&vis_graph)));
    //
    //
    //     self.ld_blocks = kosaraju_scc(&ld_graph);
    //     return (ld_links, ld_graph);
    // }


    pub fn cross_optimize(&mut self, conserved_snps: &HashSet<usize>, keep_conserved: bool, with_genotype: bool, apply_downsampling: bool) -> f64 {
        // Iteration:
        //     1. evaluate the assignment of each read based on the current SNP haplotype.
        //     2. evaluate the SNP haplotype based on the read assignment.
        // If P(sigma, delta) increase, repeat Iteration;
        // Else break;

        let mut haplotype_genotype_increase: bool = true;
        let mut haplotag_increase: bool = true;
        let mut num_iters = 0;

        while haplotype_genotype_increase | haplotag_increase {
            // optimize sigma
            let mut tmp_haplotag: HashMap<usize, i32> = HashMap::new();
            for k in 0..self.fragments.len() {
                if !self.fragments[k].for_phasing
                    || (apply_downsampling && !self.fragments[k].downsampled)
                    || self.fragments[k].haplotag == 0 {
                    continue;
                }
                let sigma_k = self.fragments[k].haplotag;
                let mut delta: Vec<i32> = Vec::new();
                let mut eta: Vec<i32> = Vec::new();
                let mut ps: Vec<i32> = Vec::new();
                let mut probs: Vec<f64> = Vec::new();
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

                if q < qn {
                    tmp_haplotag.insert(k, sigma_k * (-1));
                } else {
                    tmp_haplotag.insert(k, sigma_k);
                }
            }

            let check_val = check_new_haplotag(&self, &tmp_haplotag, apply_downsampling);
            assert!(check_val >= 0, "Error: check new haplotag: {:?}", self.candidate_snps);
            for (k, h) in tmp_haplotag.iter() {
                // when prob is equal, we still perform the flip to avoid bug of underflow
                self.fragments[*k].haplotag = *h;
            }
            if check_val == 0 {
                haplotag_increase = false;
            } else {
                haplotag_increase = true;
                haplotype_genotype_increase = true;
            }
            self.check_local_optimal_configuration(false, true, apply_downsampling);

            // optimize delta
            let mut tmp_haplotype_genotype: HashMap<usize, (i32, i32)> = HashMap::new();
            for i in 0..self.candidate_snps.len() {
                if self.candidate_snps[i].for_phasing == false { continue; }
                if keep_conserved {
                    if conserved_snps.contains(&i) { continue; }
                }
                let delta_i = self.candidate_snps[i].haplotype;
                let eta_i = self.candidate_snps[i].genotype;
                let mut sigma: Vec<i32> = Vec::new();
                let mut ps: Vec<i32> = Vec::new();
                let mut probs: Vec<f64> = Vec::new();
                for k in self.candidate_snps[i].snp_cover_fragments.iter() {
                    if !self.fragments[*k].for_phasing 
                        || (apply_downsampling && !self.fragments[*k].downsampled)
                        || self.fragments[*k].haplotag == 0 {
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
                let q1 = cal_delta_eta_sigma_log(delta_i, 0, &sigma, &ps, &probs);  // hetvar
                let q2 = cal_delta_eta_sigma_log(delta_i * (-1), 0, &sigma, &ps, &probs);   // hetvar
                let q3 = cal_delta_eta_sigma_log(delta_i, 1, &sigma, &ps, &probs);  // homref
                let q4 = cal_delta_eta_sigma_log(delta_i, -1, &sigma, &ps, &probs); // homvar

                if with_genotype {
                    let max_q = q1.max(q2.max(q3.max(q4)));
                    if q1 == max_q {
                        tmp_haplotype_genotype.insert(i, (delta_i, 0));
                    } else if q2 == max_q {
                        tmp_haplotype_genotype.insert(i, (delta_i * (-1), 0));
                    } else if q3 == max_q {
                        tmp_haplotype_genotype.insert(i, (delta_i, 1));
                    } else if q4 == max_q {
                        tmp_haplotype_genotype.insert(i, (delta_i, -1));
                    } else {
                        panic!("Error: haplotype genotype optimization failed. {}->{}->{}->{}", q1, q2, q3, q4);
                    }
                } else {
                    let mut max_q = 0.0;
                    if eta_i == 0 {
                        max_q = q1.max(q2);
                        if q1 == max_q {
                            tmp_haplotype_genotype.insert(i, (delta_i, 0));
                        } else if q2 == max_q {
                            tmp_haplotype_genotype.insert(i, (delta_i * (-1), 0));
                        }
                    } else {
                        max_q = q3.max(q4);
                        if q3 == max_q {
                            tmp_haplotype_genotype.insert(i, (delta_i, 1));
                        } else if q4 == max_q {
                            tmp_haplotype_genotype.insert(i, (delta_i, -1));
                        }
                    }
                }
            }
            let check_val = check_new_haplotype_genotype(&self, &tmp_haplotype_genotype, apply_downsampling);
            assert!(check_val >= 0, "Error: check new haplotype: {:?}", self.candidate_snps);
            for (idx, hap) in tmp_haplotype_genotype.iter() {
                // when prob is equal, we still perform the flip to avoid bug of underflow
                self.candidate_snps[*idx].haplotype = hap.0;
                self.candidate_snps[*idx].genotype = hap.1;
                // // update the haplotype of the snps in the same block
                // if block_snps.contains(&idx) {
                //     for block in blocks.iter() {
                //         if block.contains(&idx) {
                //             for ti in block.iter() {
                //                 print!("{}:{}->{},", self.candidate_snps[*ti].pos, self.candidate_snps[*ti].haplotype, hap.0);
                //                 self.candidate_snps[*ti].haplotype = hap.0;
                //             }
                //             println!();
                //         }
                //     }
                // }
            }
            if check_val == 0 {
                haplotype_genotype_increase = false;
            } else {
                haplotype_genotype_increase = true;
                haplotag_increase = true;
            }
            // self.check_local_optimal_configuration(true, false);

            num_iters += 1;
            if num_iters > 20 {
                println!("Warning: iteration inside cross_optimize exceeds 20 times.");
                break;
            }
        }
        let prob = cal_overall_probability(&self, apply_downsampling);
        prob
    }

    fn check_local_optimal_configuration(&self, used_for_haplotype_genotype: bool, used_for_haplotag: bool, apply_downsampling: bool) {
        // check sigma
        if used_for_haplotag {
            for k in 0..self.fragments.len() {
                if !self.fragments[k].for_phasing
                    || (apply_downsampling && !self.fragments[k].downsampled)
                    || self.fragments[k].haplotag == 0 {
                    continue;
                }
                let sigma_k = self.fragments[k].haplotag;
                let mut delta: Vec<i32> = Vec::new();
                let mut eta: Vec<i32> = Vec::new();
                let mut ps: Vec<i32> = Vec::new();
                let mut probs: Vec<f64> = Vec::new();
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
                assert!(q >= qn, "Error: haplotag flipping is not local optimal. {}->{}\n{:?}", q, qn, self.region);
            }
        }

        // check delta
        if used_for_haplotype_genotype {
            for i in 0..self.candidate_snps.len() {
                if self.candidate_snps[i].for_phasing == false { continue; }
                let delta_i = self.candidate_snps[i].haplotype;
                let eta_i = self.candidate_snps[i].genotype;
                let mut sigma: Vec<i32> = Vec::new();
                let mut ps: Vec<i32> = Vec::new();
                let mut probs: Vec<f64> = Vec::new();
                for k in self.candidate_snps[i].snp_cover_fragments.iter() {
                    if !self.fragments[*k].for_phasing
                        || (apply_downsampling && !self.fragments[*k].downsampled)
                        || self.fragments[*k].haplotag == 0 {
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
                let (mut q1, mut q2, mut q3, mut q4) = (0.0, 0.0, 0.0, 0.0);
                if eta_i == 0 {
                    q1 = cal_delta_eta_sigma_log(delta_i, eta_i, &sigma, &ps, &probs);
                    q2 = cal_delta_eta_sigma_log(delta_i * (-1), eta_i, &sigma, &ps, &probs);
                    // q3 = cal_delta_eta_sigma_log(delta_i, -1, &sigma, &ps, &probs);
                    // q4 = cal_delta_eta_sigma_log(delta_i, 1, &sigma, &ps, &probs);
                    assert!(q1 >= q2, "Error: haplotype and genotype flipping is not local optimal. {}->{}, {:?}, eta:{}", q1, q2, self.region, eta_i)
                } else if eta_i == 1 {
                    q1 = cal_delta_eta_sigma_log(delta_i, eta_i, &sigma, &ps, &probs);
                    q2 = cal_delta_eta_sigma_log(delta_i, eta_i * (-1), &sigma, &ps, &probs);
                    // q3 = cal_delta_eta_sigma_log(delta_i, 0, &sigma, &ps, &probs);
                    // q4 = cal_delta_eta_sigma_log(delta_i * (-1), 0, &sigma, &ps, &probs);
                    assert!(q1 >= q2, "Error: haplotype and genotype flipping is not local optimal. {}->{}, {:?}, eta:{}", q1, q2, self.region, eta_i)
                } else if eta_i == -1 {
                    q1 = cal_delta_eta_sigma_log(delta_i, eta_i, &sigma, &ps, &probs);
                    q2 = cal_delta_eta_sigma_log(delta_i, eta_i * (-1), &sigma, &ps, &probs);
                    // q3 = cal_delta_eta_sigma_log(delta_i, 0, &sigma, &ps, &probs);
                    // q4 = cal_delta_eta_sigma_log(delta_i * (-1), 0, &sigma, &ps, &probs);
                    assert!(q1 >= q2, "Error: haplotype and genotype flipping is not local optimal. {}->{}, {:?}, eta:{}", q1, q2, self.region, eta_i)
                }
                // assert!(q1 >= q2 && q1 >= q3 && q1 >= q4, "Error: haplotype and genotype flipping is not local optimal. {}->{}, {}, {}\n{:?}, eta:{}", q1, q2, q3, q4, self.region, eta_i);
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

    pub fn phase(&mut self, ld_weight_threshold: u32, max_enum_snps: usize, apply_downsampling: bool) {
        let mut largest_prob = f64::NEG_INFINITY;
        let mut best_haplotype: HashMap<usize, i32> = HashMap::new();
        let mut best_haplotag: HashMap<usize, i32> = HashMap::new();
        let mut best_genotype: HashMap<usize, i32> = HashMap::new();
        let mut conserved_snps: HashSet<usize> = HashSet::new();    //SNPs will not be flipped

        let ld_graph = self.divide_snps_into_blocks(ld_weight_threshold);
        // println!("{}", format!("{}", Dot::new(&ld_graph)));

        if self.candidate_snps.len() <= max_enum_snps {
            // enumerate the haplotype, then optimize the assignment
            let mut haplotype_enum = vec![vec![1; self.candidate_snps.len()]];
            for ti in 0..self.candidate_snps.len() {
                for tj in 0..haplotype_enum.len() {
                    let mut tmp_hap = haplotype_enum[tj].clone();
                    tmp_hap[ti] = tmp_hap[ti] * (-1);
                    haplotype_enum.push(tmp_hap);
                }
            }
            assert!(haplotype_enum.len() == 2_usize.pow(self.candidate_snps.len() as u32), "Error: Not all combinations included");
            for hap in haplotype_enum.iter() {
                self.candidate_snps.iter_mut()
                    .zip(hap)
                    .for_each(|(snp, &h)| snp.haplotype = h);
                unsafe {
                    self.init_assignment();
                    self.init_genotype();
                }
                let prob = self.cross_optimize(&conserved_snps, false, true, apply_downsampling);
                if prob > largest_prob {
                    largest_prob = prob;
                    self.save_best_configuration(&mut best_haplotype, &mut best_haplotag, &mut best_genotype);
                }
            }
            self.load_best_configuration(&best_haplotype, &best_haplotag, &best_genotype);
        } else {
            // optimize haplotype and read assignment alternatively
            // random initialization of haplotype and haplotag at each iteration
            unsafe {
                // self.init_haplotypes();
                conserved_snps = self.init_haplotypes_LD2(&ld_graph, ld_weight_threshold);
                self.init_genotype();
                self.init_assignment();
            }
            let prob = self.cross_optimize(&conserved_snps, true, false, apply_downsampling);
            if prob > largest_prob {
                largest_prob = prob;
                self.save_best_configuration(&mut best_haplotype, &mut best_haplotag, &mut best_genotype);
            }
            self.load_best_configuration(&best_haplotype, &best_haplotag, &best_genotype);

            let prob = self.cross_optimize_by_block(apply_downsampling);
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

            let mut rng = rand::thread_rng();
            for tidx in 0..=self.candidate_snps.len() / 4 {
                {
                    let flip = tidx % 2 == 1;
                    for candidate_snp in &mut self.candidate_snps {
                        let rg: f64 = rng.gen();
                        if rg < 0.1 {
                            candidate_snp.haplotype = if flip { 1 } else { -1 };
                        } else if rg >= 0.9 {
                            candidate_snp.haplotype = if flip { -1 } else { 1 };
                        }
                    }
                }
                let prob = self.cross_optimize(&conserved_snps, false, false, apply_downsampling);
                if prob > largest_prob {
                    largest_prob = prob;
                    self.save_best_configuration(&mut best_haplotype, &mut best_haplotag, &mut best_genotype);
                }
                self.load_best_configuration(&best_haplotype, &best_haplotag, &best_genotype);

                for fragment in &mut self.fragments {
                    if !fragment.for_phasing || fragment.haplotag == 0 {
                        continue;
                    }
                    if rng.gen::<f64>() < 0.1 {
                        fragment.haplotag *= -1;
                    }
                }

                let prob = self.cross_optimize(&conserved_snps, false, false, apply_downsampling);
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

    pub fn cross_optimize_by_block(&mut self, apply_downsampling: bool) -> f64 {
        // optimize delta
        let mut tmp_haplotype: HashMap<usize, i32> = HashMap::new();
        let mut tmp_haplotag: HashMap<usize, i32> = HashMap::new();
        for block in self.ld_blocks.iter() {
            let mut delta_block = Vec::new();
            let mut delta_block_flip = Vec::new();
            let mut eta_block = Vec::new();
            let mut sigma_block = Vec::new();
            let mut sigma_block_flip_vec = Vec::new();
            let mut sigma_block_flip_hashmap = HashMap::new();
            let mut ps_block = Vec::new();
            let mut probs_block = Vec::new();
            let mut block_snp_idx_set = HashSet::new();
            let mut flip_read;
            for idx in block.snp_idxes.iter() {
                block_snp_idx_set.insert(*idx);
            }
            for idx in block.snp_idxes.iter() {
                delta_block.push(self.candidate_snps[*idx].haplotype);
                delta_block_flip.push(self.candidate_snps[*idx].haplotype * (-1));
                eta_block.push(self.candidate_snps[*idx].genotype);
                let mut sigma: Vec<i32> = Vec::new();
                let mut sigma_flip: Vec<i32> = Vec::new();
                let mut ps: Vec<i32> = Vec::new();
                let mut probs: Vec<f64> = Vec::new();
                for k in self.candidate_snps[*idx].snp_cover_fragments.iter() {
                    if !self.fragments[*k].for_phasing 
                        || (apply_downsampling && !self.fragments[*k].downsampled)
                        || self.fragments[*k].haplotag == 0 {
                        continue;
                    }
                    // k is fragment index
                    flip_read = true;
                    for fe in self.fragments[*k].list.iter() {
                        if !block_snp_idx_set.contains(&fe.snp_idx) {
                            // if the read is entirely covered by the block, we will flip the read together with flipping the snp
                            flip_read = false;
                        }
                        if fe.snp_idx == *idx {
                            if fe.phase_site == false { continue; }
                            assert_ne!(fe.p, 0, "Error: phase for unexpected allele.");
                            ps.push(fe.p);
                            probs.push(fe.prob);
                            if flip_read {
                                sigma_flip.push(self.fragments[*k].haplotag * (-1));
                                sigma_block_flip_hashmap.insert(*k, self.fragments[*k].haplotag * (-1));
                            } else {
                                sigma_flip.push(self.fragments[*k].haplotag);
                                sigma_block_flip_hashmap.insert(*k, self.fragments[*k].haplotag);
                            }
                            sigma.push(self.fragments[*k].haplotag);
                        }
                    }
                }
                sigma_block.push(sigma);
                sigma_block_flip_vec.push(sigma_flip);
                ps_block.push(ps);
                probs_block.push(probs);
            }
            let q = cal_block_delta_eta_sigma_log(&delta_block, &eta_block, &sigma_block, &ps_block, &probs_block);  // hetvar
            let q_flip = cal_block_delta_eta_sigma_log(&delta_block_flip, &eta_block, &sigma_block_flip_vec, &ps_block, &probs_block);  // hetvar
            // println!("q:{}, q_flip:{}", q, q_flip);
            if q < q_flip {
                for i in 0..block.snp_idxes.len() {
                    tmp_haplotype.insert(block.snp_idxes[i], delta_block_flip[i]);
                }
                for k in 0..self.fragments.len() {
                    if sigma_block_flip_hashmap.contains_key(&k) {
                        tmp_haplotag.insert(k, sigma_block_flip_hashmap.get(&k).unwrap().clone());
                    } else {
                        tmp_haplotag.insert(k, self.fragments[k].haplotag);
                    }
                }
            } else {
                for i in 0..block.snp_idxes.len() {
                    tmp_haplotype.insert(block.snp_idxes[i], delta_block[i]);
                }
                for k in 0..self.fragments.len() {
                    tmp_haplotag.insert(k, self.fragments[k].haplotag);
                }
            }
        }
        // cal_block_delta_eta_sigma_log is quite different from cal_delta_eta_sigma_log, so the probability is not comparable
        // let check_val = check_block_new_haplotype_haplotag(&self, &tmp_haplotype, &tmp_haplotag);
        // assert!(check_val >= 0, "Error: check new haplotype: {:?}", self.candidate_snps);
        for (i, h) in tmp_haplotype.iter() {
            // when prob is equal, we still perform the flip to avoid bug of underflow
            self.candidate_snps[*i].haplotype = *h;
        }
        for (k, h) in tmp_haplotag.iter() {
            // when prob is equal, we still perform the flip to avoid bug of underflow
            self.fragments[*k].haplotag = *h;
        }
        let prob = cal_overall_probability(&self, apply_downsampling);
        prob
    }
}