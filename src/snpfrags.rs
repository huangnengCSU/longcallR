use std::collections::HashMap;

use petgraph::algo::kosaraju_scc;
use petgraph::graphmap::GraphMap;
use petgraph::Undirected;
use rust_htslib::{bam, bam::Read, bam::record::Record};

use crate::phase::{cal_delta_eta_sigma_log, cal_phase_score_log, cal_sigma_delta_eta_log};
use crate::snp::{CandidateSNP, Fragment, LD_Pair};
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
    pub allele_pairs: HashMap<[usize; 2], LD_Pair>,
    // allele pair at two snp sites, key is [snp_idx of start_node, snp_idx of end_node], start < end
    pub ld_blocks: Vec<Vec<usize>>,
    // index of candidate SNPs in each LD block
    pub min_linkers: u32,
    // the number of links for snps can be phased
}

#[derive(Debug, Clone, Default)]
pub struct Edge {
    pub snp_idxes: [usize; 2],
    // index of candidate SNPs(SNPFrag.snps), start node and end node
    pub snp_poses: [i64; 2],
    // position of candidate SNPs(SNPFrag.snps), start node and end node
    pub frag_idxes: Vec<usize>,
    // index of fragments(SNPFrag.fragments) cover this edge.
    // pub w: f64,
    // // weight of edge,  w_{ij}=\sum_{k}x_{ki}x_{kj}log\frac{1-\epsilon_{kij}}{\epsilon_{kij}}
    pub w: i32,
    // support number of the edge, each read support the edge will add 1
}

impl SNPFrag {
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
                                let mut baseq = record.qual()[pos_on_query as usize];
                                baseq = if baseq < 30 { baseq } else { 30 };
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

    /*pub fn assign_het_var_haplotype(
        &mut self,
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

            let mut phase_score = 0.0;
            if sigma.len() > 0 || hap1_reads_num >= 2 || hap2_reads_num >= 2 {
                // each haplotype should have at least 2 reads
                phase_score = -10.0_f64 * (1.0 - cal_delta_sigma_log(delta_i, &sigma, &ps, &probs)).log10(); // calaulate assignment score
                if phase_score >= min_phase_score as f64 {
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
                } else {
                    self.candidate_snps[*ti].phase_score = phase_score;
                }
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
    }*/

    /*pub fn eval_low_frac_het_var_phase(&mut self, min_phase_score: f32,
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


            let mut phase_score1 = 0.0;
            let mut phase_score2 = 0.0;
            let mut phase_score = 0.0;
            if sigma.len() == 0 || hap1_reads_num < 2 || hap2_reads_num < 2 {
                // each haplotype should have at least 2 reads
                snp.single = true; // no surranding high_frac_het_snps
                continue;
            } else {
                phase_score1 = -10.0_f64 * (1.0 - cal_delta_sigma_log(1, &sigma, &ps, &probs)).log10(); // calaulate assignment score
                phase_score2 = -10.0_f64 * (1.0 - cal_delta_sigma_log(-1, &sigma, &ps, &probs)).log10(); // calaulate assignment score
                snp.single = false;
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
    }*/
    pub fn eval_rna_edit_var_phase(&mut self, min_phase_score: f32) {
        for ti in self.edit_snps.iter() {
            let snp = &mut self.candidate_snps[*ti];
            if snp.snp_cover_fragments.len() == 0 {
                // no surranding haplotype links
                snp.single = true;
                continue;
            }
            if snp.variant_type != 1 {
                snp.non_selected = true;
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

            let mut phase_score1 = 0.0;
            let mut phase_score2 = 0.0;
            let mut phase_score = 0.0;
            if sigma.len() == 0 || hap1_reads_num < 2 || hap2_reads_num < 2 {
                // each haplotype should have at least 2 reads
                snp.single = true; // no surranding high_frac_het_snps
                continue;
            } else {
                phase_score1 = -10.0_f64 * (1.0 - cal_phase_score_log(1, 0, &sigma, &ps, &probs)).log10(); // calaulate assignment score
                phase_score2 = -10.0_f64 * (1.0 - cal_phase_score_log(-1, 0, &sigma, &ps, &probs)).log10(); // calaulate assignment score
                snp.single = false;
            }

            if phase_score1.max(phase_score2) >= min_phase_score as f64 {
                // in-phase
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
                snp.non_selected = false;
                snp.germline = true;
                snp.rna_editing = false;
                // if self.candidate_snps[*ti].alleles[0] != self.candidate_snps[*ti].reference && self.candidate_snps[*ti].alleles[1] != self.candidate_snps[*ti].reference {
                //     // tri-allelic site
                //     self.candidate_snps[*ti].variant_type = 3;
                //     self.candidate_snps[*ti].hom_var = true;
                // } else {
                //     self.candidate_snps[*ti].variant_type = 1;
                // }
                snp.haplotype = if phase_score1 >= phase_score2 { 1 } else { -1 };
                snp.genotype = 0;
                snp.variant_type = 1;
                snp.haplotype_expression = haplotype_allele_expression;
                snp.phase_score = phase_score;
            } else {
                // out-of-phase
                snp.non_selected = true;
                snp.germline = false;
                snp.rna_editing = true;
            }
        }
        // resolve single snps
        // let mut ld_idxes: Vec<usize> = Vec::new();
        // for ti in self.edit_snps.iter() {
        //     let snp = &self.candidate_snps[*ti];
        //     if snp.single {
        //         ld_idxes.push(*ti);
        //     }
        // }
        // for i in 0..ld_idxes.len() {
        //     let mut r2_map = HashMap::new();
        //     for j in 0..ld_idxes.len() {
        //         if i == j { continue; }
        //         let idx1 = ld_idxes[i];
        //         let idx2 = ld_idxes[j];
        //         let snp1 = &self.candidate_snps[idx1];
        //         let snp2 = &self.candidate_snps[idx2];
        //         if idx1 < idx2 {
        //             if !self.allele_pairs.contains_key(&[idx1, idx2]) { continue; }
        //             let r2 = self.allele_pairs.get(&[idx1, idx2]).unwrap().calculate_LD_R2(snp1.alleles[0] as u8, snp1.alleles[1] as u8, snp2.alleles[0] as u8, snp2.alleles[1] as u8);
        //             r2_map.insert([idx1, idx2], r2);
        //         } else if idx2 < idx1 {
        //             if !self.allele_pairs.contains_key(&[idx2, idx1]) { continue; }
        //             let r2 = self.allele_pairs.get(&[idx2, idx1]).unwrap().calculate_LD_R2(snp2.alleles[0] as u8, snp2.alleles[1] as u8, snp1.alleles[0] as u8, snp1.alleles[1] as u8);
        //             r2_map.insert([idx2, idx1], r2);
        //         }
        //     }
        //     let r2_map_high: HashMap<[usize; 2], f32> = r2_map.iter().filter(|(_, &v)| v > 0.9).map(|(&k, &v)| (k, v)).collect();
        //     if r2_map_high.len() >= 1 {
        //         // high LD
        //         // println!("single snp:{}", self.candidate_snps[ld_idxes[i]].pos);
        //         // for (idxes, r2) in r2_map.iter() {
        //         // println!("\t{}, {}, r2={}", self.candidate_snps[idxes[0]].pos, self.candidate_snps[idxes[1]].pos, r2);
        //         // }
        //         for (idxes, r2) in r2_map_high.iter() {
        //             let idx1 = idxes[0];
        //             let idx2 = idxes[1];
        //             self.candidate_snps[idx1].germline = true;
        //             self.candidate_snps[idx2].germline = true;
        //             self.candidate_snps[idx1].rna_editing = false;
        //             self.candidate_snps[idx2].rna_editing = false;
        //             self.candidate_snps[idx1].phase_score = min_phase_score as f64; // high LD, set phase score to min_phase_score
        //             self.candidate_snps[idx2].phase_score = min_phase_score as f64; // high LD, set phase score to min_phase_score
        //             if self.candidate_snps[idx1].alleles[0] != self.candidate_snps[idx1].reference && self.candidate_snps[idx1].alleles[1] != self.candidate_snps[idx1].reference {
        //                 // tri-allelic site
        //                 self.candidate_snps[idx1].variant_type = 3;
        //                 self.candidate_snps[idx1].hom_var = true;
        //             } else {
        //                 self.candidate_snps[idx1].variant_type = 1;
        //             }
        //             if self.candidate_snps[idx2].alleles[0] != self.candidate_snps[idx2].reference && self.candidate_snps[idx2].alleles[1] != self.candidate_snps[idx2].reference {
        //                 // tri-allelic site
        //                 self.candidate_snps[idx2].variant_type = 3;
        //                 self.candidate_snps[idx2].hom_var = true;
        //             } else {
        //                 self.candidate_snps[idx2].variant_type = 1;
        //             }
        //         }
        //     }
        // }
    }
    /*pub fn eval_som_var_phase(&mut self) {}*/
    /*pub fn eval_hom_var_phase(&mut self, min_phase_score: f32) {
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

            let mut phase_score1 = 0.0;
            let mut phase_score2 = 0.0;
            let mut phase_score = 0.0;
            if sigma.len() == 0 || hap1_reads_num < 2 || hap2_reads_num < 2 {
                // each haplotype should have at least 2 reads
                continue;
            } else {
                phase_score1 = -10.0_f64 * (1.0 - cal_delta_sigma_log(1, &sigma, &ps, &probs)).log10(); // calaulate assignment score
                phase_score2 = -10.0_f64 * (1.0 - cal_delta_sigma_log(-1, &sigma, &ps, &probs)).log10(); // calaulate assignment score
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
                if self.candidate_snps[*ti].variant_type != 3 {
                    // tri-allelic site will not be changed to het_var
                    self.candidate_snps[*ti].hom_var = false;
                    self.candidate_snps[*ti].variant_type = 1;
                }
                self.candidate_snps[*ti].haplotype = if phase_score1 >= phase_score2 { 1 } else { -1 };
                self.candidate_snps[*ti].haplotype_expression = haplotype_allele_expression;
                self.candidate_snps[*ti].phase_score = phase_score;
            }
        }
    }*/

    pub fn assign_snp_haplotype_genotype(&mut self, min_phase_score: f32) {
        // calculate phase score for each snp
        for ti in 0..self.candidate_snps.len() {
            let snp = &mut self.candidate_snps[ti];
            if !snp.for_phasing {
                snp.non_selected = true;
                continue;
            }
            if snp.snp_cover_fragments.len() == 0 {
                // no surranding haplotype links
                snp.single = true;
                continue;
            }
            // if snp.variant_type != 1 {
            //     snp.non_selected = true;
            //     continue;
            // }
            let delta_i = snp.haplotype;
            let eta_i = snp.genotype;
            let mut sigma: Vec<i32> = Vec::new();
            let mut ps: Vec<i32> = Vec::new();
            let mut probs: Vec<f64> = Vec::new();
            let mut assigns: Vec<i32> = Vec::new();
            let mut baseqs = Vec::new();
            let mut hap1_reads_num = 0;
            let mut hap2_reads_num = 0;
            // if eta_i != 0 {
            //     snp.non_selected = true;
            //     continue;
            // }
            for k in snp.snp_cover_fragments.iter() {
                // heterozygous use phased reads, otherwise use all reads
                if snp.variant_type == 1 {
                    if self.fragments[*k].assignment == 0 {
                        continue;
                    }
                }
                if self.fragments[*k].num_hete_links < self.min_linkers { continue; }
                for fe in self.fragments[*k].list.iter() {
                    if fe.snp_idx == ti {
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

            if sigma.len() == 0 {
                snp.non_selected = true;
                continue;
            }

            let q1 = cal_delta_eta_sigma_log(delta_i, 0, &sigma, &ps, &probs);
            let q2 = cal_delta_eta_sigma_log(delta_i * (-1), 0, &sigma, &ps, &probs);
            let q3 = cal_delta_eta_sigma_log(delta_i, 1, &sigma, &ps, &probs);
            let q4 = cal_delta_eta_sigma_log(delta_i, -1, &sigma, &ps, &probs);

            let max_q = q1.max(q2.max(q3.max(q4)));
            if q1 == max_q {
                snp.haplotype = delta_i;
                snp.genotype = 0;
                snp.variant_type = 1;
            } else if q2 == max_q {
                snp.haplotype = delta_i * (-1);
                snp.genotype = 0;
                snp.variant_type = 1;
            } else if q3 == max_q {
                snp.haplotype = delta_i;
                snp.genotype = 1;
                snp.variant_type = 0;
            } else if q4 == max_q {
                snp.haplotype = delta_i;
                snp.genotype = -1;
                if snp.variant_type != 2 && snp.variant_type != 3 {
                    snp.variant_type = 2;
                }
            } else {
                panic!("Error: genotype optimization failed. {},{},{},{}", q1, q2, q3, q4);
            }

            if snp.genotype != 0 {
                snp.non_selected = true;
                continue;
            }

            let mut phase_score = 0.0;
            if sigma.len() > 0 && hap1_reads_num >= 1 && hap2_reads_num >= 1 {
                // each haplotype should have at least 2 reads
                phase_score = -10.0_f64 * (1.0 - cal_phase_score_log(snp.haplotype, snp.genotype, &sigma, &ps, &probs)).log10(); // calaulate assignment score
                if phase_score >= min_phase_score as f64 {
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
                    snp.germline = true;
                    snp.haplotype_expression = haplotype_allele_expression;
                    snp.phase_score = phase_score;
                } else {
                    snp.phase_score = phase_score;
                }
            } else {
                snp.phase_score = 0.19940219;
            }

            // TODO: het var with low phase score transfer to hom var
            // if phase_score < min_phase_score as f64 && snp.allele_freqs[0] >= 0.9 && snp.alleles[0] != snp.reference && snp.dense == false {
            //     // transfer from heterozygous to homozygous
            //     snp.variant_type = 2;
            //     snp.hom_var = true;
            //     snp.germline = true;
            //     snp.for_phasing = true;
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

    pub fn assign_reads_haplotype(&mut self, read_assignment_cutoff: f64) -> HashMap<String, i32> {
        let mut read_assignments: HashMap<String, i32> = HashMap::new();
        for k in 0..self.fragments.len() {
            let sigma_k = self.fragments[k].haplotag;
            let mut delta: Vec<i32> = Vec::new();
            let mut eta: Vec<i32> = Vec::new();
            let mut ps: Vec<i32> = Vec::new();
            let mut probs: Vec<f64> = Vec::new();
            for fe in self.fragments[k].list.iter() {
                if fe.phase_site == false { continue; }
                if self.candidate_snps[fe.snp_idx].haplotype == 0 { continue; }
                if self.candidate_snps[fe.snp_idx].genotype != 0 { continue; }
                assert_ne!(fe.p, 0, "Error: phase for unexpected allele.");
                ps.push(fe.p);
                probs.push(fe.prob);
                delta.push(self.candidate_snps[fe.snp_idx].haplotype);
                eta.push(self.candidate_snps[fe.snp_idx].genotype);
            }
            if sigma_k == 0 {
                // unasigned haplotag, cluster the read into unknown group
                self.fragments[k].assignment = 0;
                self.fragments[k].haplotag = 0;
                self.fragments[k].assignment_score = 0.0;
                read_assignments.insert(self.fragments[k].read_id.clone(), 0);
            } else {
                let mut q = 0.0;
                let mut qn = 0.0;
                if delta.len() > 0 {
                    q = cal_sigma_delta_eta_log(sigma_k, &delta, &eta, &ps, &probs);
                    qn = cal_sigma_delta_eta_log(sigma_k * (-1), &delta, &eta, &ps, &probs);
                } else {
                    self.fragments[k].assignment = 0;
                    self.fragments[k].haplotag = 0;
                    self.fragments[k].assignment_score = 0.0;
                    read_assignments.insert(self.fragments[k].read_id.clone(), 0);
                    continue;
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
                    self.fragments[k].haplotag = 0;
                    self.fragments[k].assignment_score = 0.0;
                    read_assignments.insert(self.fragments[k].read_id.clone(), 0);
                }
            }
        }
        return read_assignments;
    }


    pub fn assign_phase_set(&mut self, min_phase_score: f32) -> HashMap<String, u32> {
        let mut phase_set: HashMap<String, u32> = HashMap::new();
        let mut graph: GraphMap<usize, Vec<usize>, Undirected> = GraphMap::new();  // node is index in candidate snp, edge is index in fragments
        // construct graph for hete snps
        for i in 0..self.candidate_snps.len() {
            let snp = &self.candidate_snps[i];
            if snp.genotype != 0 || snp.variant_type != 1 { continue; }
            if snp.dense || snp.rna_editing { continue; }
            if snp.phase_score < min_phase_score as f64 { continue; }
            graph.add_node(i);
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
                    let haplotype_pair = [self.candidate_snps[node_snps[j]].haplotype, self.candidate_snps[node_snps[j + 1]].haplotype];
                    let mut allele_pair = [0, 0];
                    for fe in frag.list.iter() {
                        if fe.snp_idx == node_snps[j] {
                            allele_pair[0] = fe.p;
                        } else if fe.snp_idx == node_snps[j + 1] {
                            allele_pair[1] = fe.p;
                        }
                    }
                    if haplotype_pair[0] * haplotype_pair[1] != allele_pair[0] * allele_pair[1] { continue; }
                    if !graph.contains_edge(node_snps[j], node_snps[j + 1]) {
                        graph.add_edge(node_snps[j], node_snps[j + 1], vec![k]);    // weight is a vector of fragment index, which is covered by the edge
                    } else {
                        graph.edge_weight_mut(node_snps[j], node_snps[j + 1]).unwrap().push(k);
                    }
                }
            }
        }
        let mut low_w_edges: Vec<(usize, usize)> = Vec::new();
        for edge in graph.all_edges() {
            if edge.2.len() < 2 {
                low_w_edges.push((edge.0, edge.1));
            }
        }

        for edge in low_w_edges.iter() {
            graph.remove_edge(edge.0, edge.1);
        }

        let scc = kosaraju_scc(&graph);
        // println!("{:?}", scc);
        // println!("{:?}", Dot::with_config(&graph, &[Config::EdgeNoLabel]));
        let region = self.region.clone().to_string();
        for component_nodes in scc.iter() {
            // if component_nodes.len() <= 1 {
            //     continue;
            // }
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
}