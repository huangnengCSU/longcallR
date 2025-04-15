use lazy_static::lazy_static;
use std::cmp::{min, Ordering};
use std::collections::{HashMap, HashSet};

use petgraph::{algo::kosaraju_scc, graphmap::GraphMap, Undirected};
use rust_lapper::{Interval, Lapper};
use statrs::distribution::{Binomial, DiscreteCDF};

use crate::snp::{AlternateAllele, Block, CandidateSNP, ReferenceAllele};
use crate::snpfrags::SNPFrag;
use crate::util::Profile;
use crate::vcf::GenotypeAndQuality;
use crate::{Platform, VALID_ALLELES};

fn cmp_f64(a: &f64, b: &f64) -> Ordering {
    if a < b {
        return Ordering::Less;
    } else if a > b {
        return Ordering::Greater;
    }
    return Ordering::Equal;
}

fn cal_strand_odds_ratio(ref_fw: i32, ref_rv: i32, alt_fw: i32, alt_rv: i32) -> f32 {
    // Only test strand bias for alternate allele, not for reference allele
    let x00 = (ref_fw + 1) as f32;
    let x01 = (ref_rv + 1) as f32;
    let x10 = (alt_fw + 1) as f32;
    let x11 = (alt_rv + 1) as f32;
    let symmetrical_ratio = (x00 * x11) / (x01 * x10) + (x01 * x10) / (x00 * x11);
    let ref_ratio = f32::min(x00, x01) / f32::max(x00, x01);
    let alt_ratio = f32::min(x10, x11) / f32::max(x10, x11);
    let strand_odds_ratio = f32::ln(symmetrical_ratio) + f32::ln(ref_ratio) - f32::ln(alt_ratio);
    strand_odds_ratio
}

fn binomial_two_tailed(successes: u64, trials: u64, p: f64) -> f64 {
    let binom = Binomial::new(p, trials).unwrap();

    if successes == 0 {
        2.0 * binom.cdf(0)
    } else if successes == trials {
        2.0 * (1.0 - binom.cdf(trials - 1))
    } else {
        2.0 * f64::min(binom.cdf(successes), 1.0 - binom.cdf(successes - 1))
    }
}

lazy_static! {
    static ref SOR_THRESHOLD: f32 = cal_strand_odds_ratio(5, 5, 9, 1);
}

impl SNPFrag {
    pub fn get_candidate_snps(
        &mut self,
        profile: &Profile,
        exon_region_vec: Vec<Interval<u32, u8>>,
        exon_only: bool,
        min_allele_freq: f32,
        min_variant_qual: u32,
        min_allele_freq_include_intron: f32,
        min_coverage: u32,
        max_coverage: u32,
        min_baseq: u8,
        use_strand_bias: bool,
        dense_win_size: u32,
        min_dense_cnt: u32,
        somatic_allele_frac_cutoff: f32,
        somatic_allele_cnt_cutoff: u32,
    ) {
        // get candidate SNPs, filtering with min_coverage, deletion_freq, min_allele_freq_include_intron, cover_strand_bias_threshold
        let pileup = &profile.freq_vec;
        let exon_intervaltree = Lapper::new(exon_region_vec);
        let mut position = profile.region.start - 1; // 0-based
        for bfidx in 0..pileup.len() {
            let bf = &pileup[bfidx];
            if bf.i {
                continue;
            }
            if exon_only
                && exon_intervaltree
                    .find(position + 1, position + 2)
                    .next()
                    .is_none()
            {
                // filter, not covered by exon
                position += 1;
                continue;
            }
            let total_allele_count = bf.get_allele_counts();
            if total_allele_count < min_coverage || total_allele_count > max_coverage {
                position += 1;
                continue;
            }
            let (allele1, allele1_cnt, allele2, allele2_cnt) =
                bf.get_two_major_alleles(bf.ref_base);
            let allele1_freq = (allele1_cnt as f32) / (total_allele_count as f32);
            let allele2_freq = (allele2_cnt as f32) / (total_allele_count as f32);

            let mut reference_allele = ReferenceAllele::default();
            let mut alternate_alleles = AlternateAllele::default();

            if allele1 == bf.ref_base {
                reference_allele.base = allele1;
                reference_allele.frequency = allele1_freq;
                reference_allele.count = allele1_cnt;
                alternate_alleles.num = 1;
                alternate_alleles.base[0] = allele2;
                alternate_alleles.frequency[0] = allele2_freq;
                alternate_alleles.count[0] = allele2_cnt;
            } else if allele2 == bf.ref_base {
                reference_allele.base = allele2;
                reference_allele.frequency = allele2_freq;
                reference_allele.count = allele2_cnt;
                alternate_alleles.num = 1;
                alternate_alleles.base[0] = allele1;
                alternate_alleles.frequency[0] = allele1_freq;
                alternate_alleles.count[0] = allele1_cnt;
            } else {
                reference_allele.base = bf.ref_base;
                reference_allele.frequency = 0.0;
                reference_allele.count = 0;
                alternate_alleles.num = 2;
                alternate_alleles.base[0] = allele1;
                alternate_alleles.frequency[0] = allele1_freq;
                alternate_alleles.count[0] = allele1_cnt;
                alternate_alleles.base[1] = allele2;
                alternate_alleles.frequency[1] = allele2_freq;
                alternate_alleles.count[1] = allele2_cnt;
            }

            if !VALID_ALLELES.contains(&reference_allele.base) {
                position += 1;
                continue;
            }

            if alternate_alleles.num == 0 {
                position += 1;
                continue;
            }

            if alternate_alleles.num == 1 {
                if total_allele_count < 200
                    && alternate_alleles.frequency[0] < somatic_allele_frac_cutoff
                {
                    position += 1;
                    continue;
                }
                if total_allele_count >= 200
                    && alternate_alleles.count[0] < somatic_allele_cnt_cutoff
                {
                    position += 1;
                    continue;
                }
            }

            // if alternate_alleles.num == 1 && bf.d >= alternate_alleles.count[0] {
            //     position += 1;
            //     continue;
            // } else if alternate_alleles.num == 2 && bf.d >= alternate_alleles.count[0] {
            //     position += 1;
            //     continue;
            // }

            if bf.d >= alternate_alleles.count[0] {
                position += 1;
                continue;
            }

            if (allele1_cnt + allele2_cnt) as f32 / (bf.get_depth_include_intron() as f32)
                < min_allele_freq_include_intron
            {
                position += 1;
                continue;
            }

            let allele1_quals: Vec<u8> = bf.get_allele_baseq(allele1);
            let allele2_quals: Vec<u8> = bf.get_allele_baseq(allele2);
            // filter if all alt alleles have low base quality
            if allele1 != bf.ref_base {
                let allele1_bq_pass_cnt =
                    allele1_quals.iter().filter(|&&bq| bq >= min_baseq).count();
                if allele1_cnt > 0 && allele1_bq_pass_cnt < 2 {
                    position += 1;
                    continue;
                }
            } else if allele2 != bf.ref_base {
                let allele2_bq_pass_cnt =
                    allele2_quals.iter().filter(|&&bq| bq >= min_baseq).count();
                if allele2_cnt > 0 && allele2_bq_pass_cnt < 2 {
                    position += 1;
                    continue;
                }
            }

            // strand bias
            // https://gatk.broadinstitute.org/hc/en-us/articles/360036361772-StrandOddsRatio
            // https://gatk.broadinstitute.org/hc/en-us/articles/360036728631-AS-StrandOddsRatio
            if use_strand_bias {
                let sor = if alternate_alleles.num == 1 {
                    let [ref_fw, ref_rv] = bf.get_allele_base_strands(reference_allele.base);
                    let [alt_fw, alt_rv] = bf.get_allele_base_strands(alternate_alleles.base[0]);
                    cal_strand_odds_ratio(ref_fw, ref_rv, alt_fw, alt_rv)
                } else if alternate_alleles.num == 2 {
                    let [ref_fw, ref_rv] = bf.get_allele_base_strands(reference_allele.base);
                    let [alt1_fw, alt1_rv] = bf.get_allele_base_strands(alternate_alleles.base[0]);
                    let [alt2_fw, alt2_rv] = bf.get_allele_base_strands(alternate_alleles.base[1]);
                    let sor1 = cal_strand_odds_ratio(ref_fw, ref_rv, alt1_fw, alt1_rv);
                    let sor2 = cal_strand_odds_ratio(ref_fw, ref_rv, alt2_fw, alt2_rv);
                    sor1.max(sor2)
                } else {
                    0.0
                };
                if sor > *SOR_THRESHOLD {
                    position += 1;
                    continue;
                }
                if alternate_alleles.num == 1 {
                    let [alt_fw, alt_rv] = bf.get_allele_base_strands(alternate_alleles.base[0]);
                    // less than 30x, use binomial test
                    if alt_fw + alt_rv <= 30 {
                        let p = binomial_two_tailed(alt_fw as u64, (alt_fw + alt_rv) as u64, 0.5);
                        if p < 0.05 {
                            position += 1;
                            continue;
                        }
                    }
                    // alternative allele only occurs on one strand
                    if alt_fw * alt_rv == 0 {
                        position += 1;
                        continue;
                    }
                }
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
            let max_log10_likelihood = log10_likelihood[0]
                .max(log10_likelihood[1])
                .max(log10_likelihood[2]);
            log10_likelihood[0] = 10.0_f64.powf(log10_likelihood[0] - max_log10_likelihood);
            log10_likelihood[1] = 10.0_f64.powf(log10_likelihood[1] - max_log10_likelihood);
            log10_likelihood[2] = 10.0_f64.powf(log10_likelihood[2] - max_log10_likelihood);
            let sum_log10_likelihood =
                log10_likelihood[0] + log10_likelihood[1] + log10_likelihood[2];
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

            let mut candidate_snp = CandidateSNP::default();
            candidate_snp.chromosome = profile.region.chr.clone();
            candidate_snp.pos = position as i64;
            candidate_snp.alleles = [allele1, allele2];
            candidate_snp.allele_freqs = [allele1_freq, allele2_freq];
            candidate_snp.reference = bf.ref_base;
            candidate_snp.depth = total_allele_count;
            candidate_snp.variant_quality = variant_quality;
            candidate_snp.genotype_probability = genotype_prob.clone();
            candidate_snp.genotype_quality = genotype_quality;
            if alternate_alleles.num == 1 {
                candidate_snp.alt_allele_fraction[0] = alternate_alleles.frequency[0];
            } else if alternate_alleles.num == 2 {
                candidate_snp.alt_allele_fraction[0] = alternate_alleles.frequency[0];
                candidate_snp.alt_allele_fraction[1] = alternate_alleles.frequency[1];
            } else {
                println!(
                    "Warning, abnormal number of alternate alleles {} at {}:{}",
                    alternate_alleles.num, profile.region.chr, position
                );
            }

            if genotype_prob[0] > genotype_prob[1] && genotype_prob[0] > genotype_prob[2] {
                // homozygous variant
                candidate_snp.variant_type = 2;
                candidate_snp.genotype = -1;
            } else if genotype_prob[1] > genotype_prob[0] && genotype_prob[1] > genotype_prob[2] {
                // heterozygous variant
                candidate_snp.variant_type = 1;
                candidate_snp.genotype = 0;
            } else {
                // homozygous reference
                candidate_snp.variant_type = 0;
                candidate_snp.genotype = 1;
            }

            // Low QUAL sites are filtered out
            if variant_quality < min_variant_qual as f64 {
                position += 1;
                continue;
            }

            // candidate rna editing site
            let forward_transcript_cnt = bf.transcript_strands[0];
            let reverse_transcript_cnt = bf.transcript_strands[1];
            if reference_allele.base == 'A'
                && alternate_alleles.base[0] == 'G'
                && (forward_transcript_cnt > reverse_transcript_cnt * 2
                    || (forward_transcript_cnt == 0 && reverse_transcript_cnt == 0))  // sometimes for missing `ts` tag
                && candidate_snp.variant_type != 2
            {
                candidate_snp.rna_editing = true;
                candidate_snp.for_phasing = false;
                self.candidate_snps.push(candidate_snp);
                self.edit_snps.push(self.candidate_snps.len() - 1);
                position += 1;
                continue;
            }
            if reference_allele.base == 'T'
                && alternate_alleles.base[0] == 'C'
                && (reverse_transcript_cnt > forward_transcript_cnt * 2
                    || (forward_transcript_cnt == 0 && reverse_transcript_cnt == 0))  // sometimes for missing `ts` tag
                && candidate_snp.variant_type != 2
            {
                candidate_snp.rna_editing = true;
                candidate_snp.for_phasing = false;
                self.candidate_snps.push(candidate_snp);
                self.edit_snps.push(self.candidate_snps.len() - 1);
                position += 1;
                continue;
            }

            // candidate somatic mutation
            if alternate_alleles.num == 1 && alternate_alleles.frequency[0] < min_allele_freq {
                candidate_snp.cand_somatic = true;
                candidate_snp.for_phasing = false;
                self.candidate_snps.push(candidate_snp);
                self.somatic_snps.push(self.candidate_snps.len() - 1);
                position += 1;
                continue;
            }

            if candidate_snp.variant_type == 2 {
                // hom_var
                if alternate_alleles.num == 2
                    && alternate_alleles.frequency[0] >= min_allele_freq
                    && alternate_alleles.frequency[1] >= min_allele_freq
                {
                    candidate_snp.variant_type = 3; // triallelic SNP
                    candidate_snp.genotype = -1;
                }
                candidate_snp.hom_var = true;
                candidate_snp.for_phasing = true;
                self.candidate_snps.push(candidate_snp);
                self.homo_snps.push(self.candidate_snps.len() - 1);
                position += 1;
                continue;
            }

            if candidate_snp.variant_type == 1 {
                if alternate_alleles.num == 2 {
                    candidate_snp.variant_type = 3; // triallelic SNP
                    candidate_snp.genotype = -1;
                    candidate_snp.hom_var = true;
                    candidate_snp.for_phasing = true;
                    self.candidate_snps.push(candidate_snp);
                    self.homo_snps.push(self.candidate_snps.len() - 1);
                    position += 1;
                    continue;
                }
                if alternate_alleles.num == 1 {
                    candidate_snp.het_var = true;
                    candidate_snp.for_phasing = true;
                    self.candidate_snps.push(candidate_snp);
                    self.het_snps.push(self.candidate_snps.len() - 1);
                    position += 1;
                    continue;
                }
            }

            if candidate_snp.variant_type == 0 {
                position += 1;
                continue;
            }

            position += 1;
        }

        let mut concat_idxes = Vec::new();
        concat_idxes.extend(self.homo_snps.clone());
        concat_idxes.extend(self.het_snps.clone());
        // concat_idxes.extend(self.edit_snps.clone());
        concat_idxes.sort();

        // filter dense, in dense_win_size, over min_dense_cnt
        for i in 0..concat_idxes.len() {
            let start_pos = self.candidate_snps[concat_idxes[i]].pos;
            for j in i..concat_idxes.len() {
                let current_pos = self.candidate_snps[concat_idxes[j]].pos;
                let diff = current_pos - start_pos;

                // If the window exceeds the allowed size, check if the previous window had enough SNPs
                if diff > dense_win_size as i64 {
                    if (j - i) as u32 >= min_dense_cnt {
                        for tk in i..j {
                            self.candidate_snps[concat_idxes[tk]].dense = true;
                            self.candidate_snps[concat_idxes[tk]].for_phasing = false;
                        }
                    }
                    break;
                }

                // If reach the last element and still within the window, check the count.
                if j == concat_idxes.len() - 1 && (j - i + 1) as u32 >= min_dense_cnt {
                    for tk in i..j {
                        self.candidate_snps[concat_idxes[tk]].dense = true;
                        self.candidate_snps[concat_idxes[tk]].for_phasing = false;
                    }
                }
            }
        }

        // filter dense, 5 bp window, over 3 snps
        for i in 0..concat_idxes.len() {
            let start_pos = self.candidate_snps[concat_idxes[i]].pos;
            for j in i..concat_idxes.len() {
                let current_pos = self.candidate_snps[concat_idxes[j]].pos;
                let diff = current_pos - start_pos;
                if diff >= 5 {
                    if (j - i) as u32 >= 3 {
                        for tk in i..j {
                            self.candidate_snps[concat_idxes[tk]].dense = true;
                            self.candidate_snps[concat_idxes[tk]].for_phasing = false;
                        }
                    }
                    break;
                }

                // If reach the last element and still within the window, check the count.
                if j == concat_idxes.len() - 1 && (j - i + 1) as u32 >= 3 {
                    for tk in i..j {
                        self.candidate_snps[concat_idxes[tk]].dense = true;
                        self.candidate_snps[concat_idxes[tk]].for_phasing = false;
                    }
                }
            }
        }

        self.homo_snps.retain(|&i| !self.candidate_snps[i].dense);
        self.het_snps.retain(|&i| !self.candidate_snps[i].dense);
        // self.edit_snps.retain(|&i| !self.candidate_snps[i].dense);
    }

    pub fn import_external_candidates(
        &mut self,
        profile: &Profile,
        ref_seq: &Vec<u8>,
        chr_candidates_genotype_quality: &HashMap<usize, GenotypeAndQuality>,
        min_variant_qual: f32,
    ) {
        let pileup = &profile.freq_vec;
        let mut position = profile.region.start - 1; // 0-based
        for bfidx in 0..pileup.len() {
            let bf = &pileup[bfidx];
            if bf.i {
                continue;
            }
            let (allele1, allele1_cnt, allele2, allele2_cnt) =
                bf.get_two_major_alleles(bf.ref_base);
            if chr_candidates_genotype_quality.contains_key(&(position as usize)) {
                let ref_pos = position as usize; // 0-based
                let ref_base = ref_seq[ref_pos] as char;
                let genotype_quality = chr_candidates_genotype_quality.get(&(ref_pos)).unwrap();
                if genotype_quality.quality < min_variant_qual {
                    position += 1;
                    continue;
                }
                let total_allele_count = bf.get_allele_counts();
                let allele1_freq = (allele1_cnt as f32) / (total_allele_count as f32);
                let allele2_freq = (allele2_cnt as f32) / (total_allele_count as f32);
                let mut candidate_snp = CandidateSNP::default();
                candidate_snp.chromosome = profile.region.chr.clone();
                candidate_snp.pos = ref_pos as i64;
                candidate_snp.reference = ref_base;
                candidate_snp.alleles = [allele1, allele2];
                candidate_snp.allele_freqs = [allele1_freq, allele2_freq];
                candidate_snp.depth = total_allele_count;
                candidate_snp.variant_quality = genotype_quality.quality as f64;
                candidate_snp.genotype_quality = genotype_quality.quality as f64;

                // genotype: 0: 0|0, 1: 0|1, 2: 1|1, 3: 1|2
                match genotype_quality.genotype {
                    0 => {
                        // 0|0
                        candidate_snp.variant_type = 0;
                        candidate_snp.genotype = 1;
                    }
                    1 => {
                        // 0|1
                        candidate_snp.variant_type = 1;
                        candidate_snp.genotype = 0;
                        candidate_snp.for_phasing = true;
                        candidate_snp.het_var = true;
                        self.candidate_snps.push(candidate_snp);
                        self.het_snps.push(self.candidate_snps.len() - 1);
                    }
                    2 => {
                        // 1|1
                        candidate_snp.variant_type = 2;
                        candidate_snp.genotype = -1;
                        candidate_snp.for_phasing = true;
                        candidate_snp.hom_var = true;
                        self.candidate_snps.push(candidate_snp);
                        self.homo_snps.push(self.candidate_snps.len() - 1);
                    }
                    3 => {
                        candidate_snp.variant_type = 3;
                        candidate_snp.genotype = -1;
                        candidate_snp.hom_var = true;
                        self.candidate_snps.push(candidate_snp);
                        self.het_snps.push(self.candidate_snps.len() - 1);
                    }
                    _ => {
                        println!(
                            "Error: unknown genotype {:?} at {:?}:{:?}",
                            genotype_quality.genotype,
                            candidate_snp.chromosome,
                            candidate_snp.pos + 1 // vcf position 1-based
                        );
                        position += 1;
                        continue;
                    }
                }
            }
            position += 1;
        }
    }

    pub fn divide_snps_into_blocks(
        &mut self,
        ld_weight_threshold: u32,
    ) -> GraphMap<usize, i32, Undirected> {
        let ld_idxes: Vec<usize> = self
            .candidate_snps
            .iter()
            .enumerate()
            .filter(|(_, snp)| snp.for_phasing)
            .map(|(index, _)| index)
            .collect();

        let mut pass_ld_pair: Vec<(usize, usize)> = Vec::new();
        for i in 0..ld_idxes.len() {
            for j in i + 1..ld_idxes.len() {
                if i == j {
                    continue;
                }
                let idx1 = ld_idxes[i];
                let idx2 = ld_idxes[j];
                let snp1 = &self.candidate_snps[idx1];
                let snp2 = &self.candidate_snps[idx2];
                let (mut snp1_ref, mut snp1_alt, mut snp2_ref, mut snp2_alt);
                let (mut snp1_ref_frac, mut snp1_alt_frac, mut snp2_ref_frac, mut snp2_alt_frac) =
                    (0.0, 0.0, 0.0, 0.0);
                if snp1.alleles[0] == snp1.reference && snp1.alleles[1] != snp1.reference {
                    snp1_ref = snp1.alleles[0] as u8;
                    snp1_ref_frac = snp1.allele_freqs[0];
                    snp1_alt = snp1.alleles[1] as u8;
                    snp1_alt_frac = snp1.allele_freqs[1];
                } else if snp1.alleles[0] != snp1.reference && snp1.alleles[1] == snp1.reference {
                    snp1_ref = snp1.alleles[1] as u8;
                    snp1_ref_frac = snp1.allele_freqs[1];
                    snp1_alt = snp1.alleles[0] as u8;
                    snp1_alt_frac = snp1.allele_freqs[0];
                } else {
                    continue;
                }
                if snp2.alleles[0] == snp2.reference && snp2.alleles[1] != snp2.reference {
                    snp2_ref = snp2.alleles[0] as u8;
                    snp2_ref_frac = snp2.allele_freqs[0];
                    snp2_alt = snp2.alleles[1] as u8;
                    snp2_alt_frac = snp2.allele_freqs[1];
                } else if snp2.alleles[0] != snp2.reference && snp2.alleles[1] == snp2.reference {
                    snp2_ref = snp2.alleles[1] as u8;
                    snp2_ref_frac = snp2.allele_freqs[1];
                    snp2_alt = snp2.alleles[0] as u8;
                    snp2_alt_frac = snp2.allele_freqs[0];
                } else {
                    continue;
                }
                assert!(idx1 < idx2, "Error: unexpected index order.");
                if !self.allele_pairs.contains_key(&[idx1, idx2]) {
                    continue;
                }
                if snp1_ref_frac == 0.0
                    || snp1_alt_frac == 0.0
                    || snp2_ref_frac == 0.0
                    || snp2_alt_frac == 0.0
                {
                    continue;
                }
                // score is like the ratio of reads conflict the main LD pair.
                // weight is like the number of reads supporting the LD.
                let (score, weight) = self
                    .allele_pairs
                    .get(&[idx1, idx2])
                    .unwrap()
                    .calculate_ld(snp1_ref, snp1_alt, snp2_ref, snp2_alt);
                self.allele_pairs
                    .get_mut(&[idx1, idx2])
                    .unwrap()
                    .set_ld(score, weight);
                if score == 0.0 {
                    pass_ld_pair.push((idx1, idx2));
                }
            }
        }

        // construct the graph to find the connected components, which are the LD blocks
        let mut ld_graph: GraphMap<usize, i32, Undirected> = GraphMap::new(); // node is index in candidate snp, edge is index in fragments
        for (node1, node2) in pass_ld_pair.iter() {
            let ld_weight = self.allele_pairs.get(&[*node1, *node2]).unwrap().weight;
            if ld_graph.contains_edge(*node1, *node2) {
                *ld_graph.edge_weight_mut(*node1, *node2).unwrap() += ld_weight;
            } else {
                ld_graph.add_edge(*node1, *node2, ld_weight);
            }
        }

        // delete edges with weight < weight_threshold
        let low_w_edges: Vec<(usize, usize)> = ld_graph
            .all_edges()
            .filter(|ld_edge| (ld_edge.2.abs() as u32) < ld_weight_threshold)
            .map(|ld_edge| (ld_edge.0, ld_edge.1))
            .collect();
        for edge in low_w_edges.iter() {
            ld_graph.remove_edge(edge.0, edge.1);
        }

        // // visualize the LD graph, node is position of snp, edge is the weight of LD
        // let mut vis_graph: GraphMap<i64, i32, Undirected> = GraphMap::new();
        // for (edge, pair) in self.allele_pairs.iter() {
        //     if pair.score != 0.0 { continue; }  // not perfect LD, make sure each edge is perfect LD
        //     if vis_graph.contains_edge(self.candidate_snps[edge[0]].pos, self.candidate_snps[edge[1]].pos) {
        //         let w = vis_graph.edge_weight_mut(self.candidate_snps[edge[0]].pos, self.candidate_snps[edge[1]].pos).unwrap();
        //         *w += pair.weight;
        //     } else {
        //         vis_graph.add_edge(self.candidate_snps[edge[0]].pos, self.candidate_snps[edge[1]].pos, pair.weight);
        //     }
        // }
        // let mut low_w_edges2: Vec<(i64, i64)> = Vec::new();
        // for edge in vis_graph.all_edges() {
        //     if (edge.2.abs() as u32) < ld_weight_threshold {
        //         low_w_edges2.push((edge.0, edge.1));
        //     }
        // }
        // for edge in low_w_edges2.iter() {
        //     vis_graph.remove_edge(edge.0, edge.1);
        // }
        // println!("{}", format!("{}", Dot::new(&vis_graph)));

        let component = kosaraju_scc(&ld_graph); // each component is an LD block
        component.iter().enumerate().for_each(|(cid, indices)| {
            let block = Block {
                block_id: cid as u32,
                snp_idxes: indices.clone(),
                ..Default::default()
            };
            self.ld_blocks.push(block);
        });
        return ld_graph;
    }
}
