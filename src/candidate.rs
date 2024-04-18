use std::cmp::Ordering;

use rust_lapper::{Interval, Lapper};

use crate::snp::CandidateSNP;
use crate::Platform;
use crate::snpfrags::SNPFrag;
use crate::util::Profile;

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
        platform: &Platform,
        exon_region_vec: Vec<Interval<usize, u8>>,
        min_allele_freq: f32,
        hetvar_high_frac_cutoff: f32,
        min_allele_freq_include_intron: f32,
        min_coverage: u32,
        max_coverage: u32,
        min_baseq: u8,
        use_strand_bias: bool,
        strand_bias_threshold: f32,
        cover_strand_bias_threshold: f32,
        distance_to_splicing_site: u32,
        window_size: u32,
        dense_win_size: u32,
        min_dense_cnt: u32,
        somatic_allele_frac_cutoff: f32,
        somatic_allele_cnt_cutoff: u32,
    ) {
        // get candidate SNPs, filtering with min_coverage, deletion_freq, min_allele_freq_include_intron, cover_strand_bias_threshold
        let pileup = &profile.freq_vec;
        let mut use_annotation: bool = false;
        if exon_region_vec.len() > 0 { use_annotation = true; }
        let mut exon_intervaltree = Lapper::new(exon_region_vec);
        let mut position = profile.region.start - 1; // 0-based
        for bfidx in 0..pileup.len() {
            let debug_pos = 1182052;
            let bf = &pileup[bfidx];
            if bf.i {
                continue;
            }
            if use_annotation && exon_intervaltree.find((position + 1) as usize, (position + 2) as usize).count() == 0 {
                // filter, not covered by exon
                position += 1;
                continue;
            }
            let depth = bf.get_depth_exclude_intron_deletion();
            if depth < min_coverage || depth > max_coverage {
                position += 1;
                continue;
            }
            let (allele1, allele1_cnt, allele2, allele2_cnt) = bf.get_two_major_alleles(bf.ref_base);
            if allele1 != bf.ref_base {
                if bf.d >= allele1_cnt {
                    position += 1;
                    continue;
                }
            } else if allele2 != bf.ref_base {
                if bf.d >= allele2_cnt {
                    position += 1;
                    continue;
                }
            }
            let depth_include_intron = bf.get_depth_include_intron();
            if (allele1_cnt as f32 + allele2_cnt as f32) / (depth_include_intron as f32) < min_allele_freq_include_intron {
                // only ont reads have this filter, hifi reads don't have this filter
                // maybe caused by erroneous intron alignment
                position += 1;
                continue;
            }
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
                    println!("Error: unknown allele: {}", allele1);
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
                    println!("Error: unknown allele: {}", allele2);
                }
            }
            // filter if all alt alleles have low base quality
            if allele1 != bf.ref_base {
                let mut allele1_bq_pass_cnt = 0;
                for bq in allele1_quals.iter() {
                    if *bq >= min_baseq {
                        allele1_bq_pass_cnt += 1;
                    }
                }
                if allele1_cnt > 0 && allele1_bq_pass_cnt < 2 {
                    position += 1;
                    continue;
                }
            } else if allele2 != bf.ref_base {
                let mut allele2_bq_pass_cnt = 0;
                for bq in allele2_quals.iter() {
                    if *bq >= min_baseq {
                        allele2_bq_pass_cnt += 1;
                    }
                }
                if allele2_cnt > 0 && allele2_bq_pass_cnt < 2 {
                    position += 1;
                    continue;
                }
            }

            if use_strand_bias {
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
                    let freads = bf.forward_cnt as i32;
                    let breads = bf.backward_cnt as i32;
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
                        if fcnt as f32 / total_cnt as f32 > strand_bias_threshold && 0.5_f32.powf((breads - bcnt) as f32) < 0.002 {
                            println!("strand bias, fread:{}, fcnt:{}, bread:{}, bcnt:{}, total_cnt:{}", freads, fcnt, breads, bcnt, total_cnt);
                            strand_bias = true;
                        } else if bcnt as f32 / total_cnt as f32 > strand_bias_threshold && 0.5_f32.powf((freads - fcnt) as f32) < 0.002 {
                            println!("strand bias, fread:{}, fcnt:{}, bread:{}, bcnt:{}, total_cnt:{}", freads, fcnt, breads, bcnt, total_cnt);
                            strand_bias = true;
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

            // if platform == Platform::ont {
            match platform {
                Platform::ont => {
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
                }
                _ => {}
            }

            if position == debug_pos - 1 {
                println!("pass filter");
            }

            // genotype likelihood
            let mut loglikelihood = [0.0, 0.0, 0.0];
            let mut logprob = [0.0, 0.0, 0.0];
            let theta = 0.001; // mutation rate
            let background_prob: [f64; 3] = [theta / 2.0, theta, 1.0 - 1.5 * theta]; // background of probability of observe homo variant, hete variant and homo reference
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
                let bq = if *bq < 33 { *bq } else { 33 };
                let error_rate = 0.1_f64.powf((bq as f64) / 10.0);
                loglikelihood[0] += error_rate.log10();
                loglikelihood[2] += (1.0 - error_rate).log10();
            }

            for bq_vec in different_baseqs.iter() {
                for bq in bq_vec.iter() {
                    let bq = if *bq < 33 { *bq } else { 33 };
                    let error_rate = 0.1_f64.powf((bq as f64) / 10.0);
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
            logprob[0] += background_prob[0].log10();
            logprob[1] += background_prob[1].log10();
            logprob[2] += background_prob[2].log10();
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
            let variant_quality = -10.0 * (10e-301_f64.max(variant_prob[2]).log10()); // if variant_quality is greater than 3000, we set it to 3000

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

            if position == debug_pos - 1 {
                println!("variant quality: {:?}, genotype quality: {:?}", variant_quality, genotype_prob);
            }

            let allele1_freq = (allele1_cnt as f32) / (depth as f32);
            let allele2_freq = (allele2_cnt as f32) / (depth as f32);
            let mut candidate_snp = CandidateSNP::default();
            candidate_snp.chromosome = profile.region.chr.clone().into_bytes();
            candidate_snp.pos = position as i64;
            candidate_snp.alleles = [allele1, allele2];
            candidate_snp.allele_freqs = [allele1_freq, allele2_freq];
            candidate_snp.reference = bf.ref_base;
            candidate_snp.depth = depth;
            candidate_snp.variant_quality = variant_quality;
            candidate_snp.genotype_probability = genotype_prob.clone();
            candidate_snp.genotype_quality = genotype_quality;
            if genotype_prob[0] > genotype_prob[1] && genotype_prob[0] > genotype_prob[2] {
                candidate_snp.variant_type = 2;
            } else if genotype_prob[1] > genotype_prob[0] && genotype_prob[1] > genotype_prob[2] {
                candidate_snp.variant_type = 1;
            } else {
                candidate_snp.variant_type = 0;
            }

            if allele1 == bf.ref_base && allele2 != bf.ref_base {
                if allele2_cnt < somatic_allele_cnt_cutoff && allele2_freq < somatic_allele_frac_cutoff {
                    position += 1;
                    continue;
                }
            } else if allele2 == bf.ref_base && allele1 != bf.ref_base {
                if allele1_cnt < somatic_allele_cnt_cutoff && allele1_freq < somatic_allele_frac_cutoff {
                    position += 1;
                    continue;
                }
            } else if allele1 != bf.ref_base && allele2 != bf.ref_base {
                if allele1_cnt < somatic_allele_cnt_cutoff &&
                    allele1_freq < somatic_allele_frac_cutoff &&
                    allele2_cnt < somatic_allele_cnt_cutoff &&
                    allele2_freq < somatic_allele_frac_cutoff
                {
                    position += 1;
                    continue;
                }
            }

            // candidate rna editing site
            if (bf.ref_base == 'A' && ((allele1 == 'G' && allele1_cnt > 0) || (allele2 == 'G' && allele2_cnt > 0))) || (bf.ref_base == 'T' && ((allele1 == 'C' && allele1_cnt > 0) || (allele2 == 'C' && allele2_cnt > 0))) {
                if candidate_snp.variant_type != 2 {
                    candidate_snp.rna_editing = true;
                    self.candidate_snps.push(candidate_snp);
                    self.edit_snps.push(self.candidate_snps.len() - 1);
                    position += 1;
                    continue;
                }
            }

            // candidate somatic mutation
            if allele1 == bf.ref_base && allele2 != bf.ref_base {
                if allele2_cnt >= somatic_allele_cnt_cutoff && allele2_freq >= somatic_allele_frac_cutoff && allele2_freq < min_allele_freq {
                    candidate_snp.cand_somatic = true;
                    self.candidate_snps.push(candidate_snp);
                    self.somatic_snps.push(self.candidate_snps.len() - 1);
                    position += 1;
                    continue;
                }
            } else if allele2 == bf.ref_base && allele1 != bf.ref_base {
                if allele1_cnt >= somatic_allele_cnt_cutoff && allele1_freq >= somatic_allele_frac_cutoff && allele1_freq < min_allele_freq {
                    candidate_snp.cand_somatic = true;
                    self.candidate_snps.push(candidate_snp);
                    self.somatic_snps.push(self.candidate_snps.len() - 1);
                    position += 1;
                    continue;
                }
            }

            if candidate_snp.variant_type == 2 {
                // hom_var
                if allele1 != bf.ref_base && allele2 != bf.ref_base && allele1_freq >= min_allele_freq && allele2_freq >= min_allele_freq {
                    candidate_snp.variant_type = 3; // triallelic SNP
                }
                candidate_snp.hom_var = true;
                candidate_snp.germline = true;
                if candidate_snp.pos == (debug_pos - 1) as i64 {
                    println!("candidate_snp1: {:?}", candidate_snp);
                }
                self.candidate_snps.push(candidate_snp);
                self.homo_snps.push(self.candidate_snps.len() - 1);
                position += 1;
                continue;
            }

            if candidate_snp.variant_type == 1 {
                // het_var
                if allele1 != bf.ref_base && allele2 != bf.ref_base {
                    candidate_snp.variant_type = 3; // triallelic SNP
                    candidate_snp.hom_var = true;
                    candidate_snp.germline = true;
                    self.candidate_snps.push(candidate_snp);
                    self.homo_snps.push(self.candidate_snps.len() - 1);
                    position += 1;
                    continue;
                    // if allele1_freq >= min_allele_freq && allele2_freq >= min_allele_freq {
                    //     candidate_snp.variant_type = 3; // triallelic SNP
                    //     candidate_snp.hom_var = true;
                    //     self.candidate_snps.push(candidate_snp);
                    //     self.homo_snps.push(self.candidate_snps.len() - 1);
                    //     position += 1;
                    //     continue;
                    // } else if allele1_freq >= min_allele_freq {
                    //     if allele1_freq >= hetvar_high_frac_cutoff {
                    //         candidate_snp.high_frac_het = true;
                    //         candidate_snp.for_phasing = true;
                    //         self.candidate_snps.push(candidate_snp);
                    //         self.high_frac_het_snps.push(self.candidate_snps.len() - 1);
                    //         position += 1;
                    //         continue;
                    //     } else {
                    //         candidate_snp.low_frac_het = true;
                    //         self.candidate_snps.push(candidate_snp);
                    //         self.low_frac_het_snps.push(self.candidate_snps.len() - 1);
                    //         position += 1;
                    //         continue;
                    //     }
                    // } else if allele2_freq >= min_allele_freq {
                    //     if allele2_freq >= hetvar_high_frac_cutoff {
                    //         candidate_snp.high_frac_het = true;
                    //         candidate_snp.for_phasing = true;
                    //         self.candidate_snps.push(candidate_snp);
                    //         self.high_frac_het_snps.push(self.candidate_snps.len() - 1);
                    //         position += 1;
                    //         continue;
                    //     } else {
                    //         candidate_snp.low_frac_het = true;
                    //         self.candidate_snps.push(candidate_snp);
                    //         self.low_frac_het_snps.push(self.candidate_snps.len() - 1);
                    //         position += 1;
                    //         continue;
                    //     }
                    // } else { panic!("Error: unexpected condition"); }
                } else if allele1 != bf.ref_base && allele2 == bf.ref_base {
                    if allele1_freq >= hetvar_high_frac_cutoff {
                        candidate_snp.high_frac_het = true;
                        candidate_snp.for_phasing = true;
                        self.candidate_snps.push(candidate_snp);
                        self.high_frac_het_snps.push(self.candidate_snps.len() - 1);
                        position += 1;
                        continue;
                    } else {
                        candidate_snp.low_frac_het = true;
                        self.candidate_snps.push(candidate_snp);
                        self.low_frac_het_snps.push(self.candidate_snps.len() - 1);
                        position += 1;
                        continue;
                    }
                } else if allele2 != bf.ref_base && allele1 == bf.ref_base {
                    if allele2_freq >= hetvar_high_frac_cutoff {
                        candidate_snp.high_frac_het = true;
                        candidate_snp.for_phasing = true;
                        self.candidate_snps.push(candidate_snp);
                        self.high_frac_het_snps.push(self.candidate_snps.len() - 1);
                        position += 1;
                        continue;
                    } else {
                        candidate_snp.low_frac_het = true;
                        self.candidate_snps.push(candidate_snp);
                        self.low_frac_het_snps.push(self.candidate_snps.len() - 1);
                        position += 1;
                        continue;
                    }
                } else { panic!("Error: unexpected condition"); }
            }

            if candidate_snp.variant_type == 0 {
                if allele1 != bf.ref_base {
                    assert!(allele1_freq < somatic_allele_frac_cutoff || allele1_cnt < somatic_allele_cnt_cutoff, "candidate: {:?}", candidate_snp);
                }
                if allele2 != bf.ref_base {
                    assert!(allele2_freq < somatic_allele_frac_cutoff || allele2_cnt < somatic_allele_cnt_cutoff, "candidate: {:?}", candidate_snp);
                }
                position += 1;
                continue;
            }

            position += 1;
            println!("Strange candidate: {:?}", candidate_snp);
        }

        // filter dense region, hom_var + high_frac_het + low_frac_het
        let mut concat_idxes = Vec::new();
        concat_idxes.extend(self.homo_snps.clone());
        concat_idxes.extend(self.high_frac_het_snps.clone());
        concat_idxes.extend(self.low_frac_het_snps.clone());
        concat_idxes.sort();
        for i in 0..concat_idxes.len() {
            for j in i..concat_idxes.len() {
                if j == concat_idxes.len() - 1 {
                    // distance from snp i to end snp is smaller than dense_win_size
                    if self.candidate_snps[concat_idxes[j]].pos - self.candidate_snps[concat_idxes[i]].pos <= dense_win_size as i64 && (j - i + 1) as u32 >= min_dense_cnt {
                        for tk in i..j {
                            self.candidate_snps[concat_idxes[tk]].dense = true;
                            self.candidate_snps[concat_idxes[tk]].for_phasing = false;
                        }
                    }
                }
                if self.candidate_snps[concat_idxes[j]].pos - self.candidate_snps[concat_idxes[i]].pos > dense_win_size as i64 {
                    if (j - 1 - i + 1) as u32 >= min_dense_cnt {
                        for tk in i..j {
                            self.candidate_snps[concat_idxes[tk]].dense = true;
                            self.candidate_snps[concat_idxes[tk]].for_phasing = false;
                        }
                    }
                    break;
                }
            }
        }

        let mut tmp_idxes = Vec::new();
        for i in self.homo_snps.iter() {
            if self.candidate_snps[*i].dense {
                continue;
            }
            tmp_idxes.push(*i);
        }
        self.homo_snps = tmp_idxes;

        let mut tmp_idxes = Vec::new();
        for i in self.low_frac_het_snps.iter() {
            if self.candidate_snps[*i].dense {
                continue;
            }
            tmp_idxes.push(*i);
        }
        self.low_frac_het_snps = tmp_idxes;

        let mut tmp_idxes = Vec::new();
        for i in self.high_frac_het_snps.iter() {
            if self.candidate_snps[*i].dense {
                continue;
            }
            tmp_idxes.push(*i);
        }
        self.high_frac_het_snps = tmp_idxes;

        let debug_pos = 1182052;
        for s in self.candidate_snps.iter() {
            if s.pos == debug_pos - 1 {
                println!("candidate_snp2: {:?}\n", s);
            }
        }
    }
}