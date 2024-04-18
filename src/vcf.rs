use crate::snpfrags::SNPFrag;

#[derive(Debug, Default, Clone)]
pub struct VCFRecord {
    pub chromosome: Vec<u8>,
    pub position: u64,
    pub id: Vec<u8>,
    pub reference: Vec<u8>,
    pub alternative: Vec<Vec<u8>>,
    pub qual: i32,
    pub filter: Vec<u8>,
    pub info: Vec<u8>,
    pub format: Vec<u8>,
    pub genotype: String,
    /* private fields */
}


impl SNPFrag {
    pub fn output_phased_vcf(&mut self, min_qual_for_candidate: u32) -> Vec<VCFRecord> {
        let mut records: Vec<VCFRecord> = Vec::new();
        for i in 0..self.candidate_snps.len() {
            let snp = &self.candidate_snps[i];
            let debug_pos = 1182052;
            if snp.pos == debug_pos - 1 {
                println!("phased output: {:?}", snp);
            }

            if snp.rna_editing {
                let mut rd: VCFRecord = VCFRecord::default();
                rd.chromosome = snp.chromosome.clone();
                rd.position = snp.pos as u64 + 1; // position in vcf format is 1-based
                rd.id = vec!['.' as u8];
                rd.reference = vec![snp.reference as u8];
                if snp.variant_type == 1 || snp.variant_type == 2 {
                    if snp.alleles[0] != snp.reference {
                        rd.alternative = vec![vec![snp.alleles[0] as u8]];
                    } else if snp.alleles[1] != snp.reference {
                        rd.alternative = vec![vec![snp.alleles[1] as u8]];
                    }
                } else if snp.variant_type == 3 {
                    rd.alternative = vec![vec![snp.alleles[0] as u8], vec![snp.alleles[1] as u8]];
                }
                rd.qual = snp.variant_quality as i32;
                rd.filter = "RnaEdit".to_string().into_bytes();
                rd.info = "RDS=.".to_string().into_bytes();
                let gt;
                if snp.variant_type == 1 {
                    gt = "0/1"
                } else if snp.variant_type == 2 {
                    gt = "1/1"
                } else if snp.variant_type == 3 {
                    gt = "1/2"
                } else {
                    continue;
                }
                rd.genotype = format!(
                    "{}:{}:{}:{:.2}",
                    gt,
                    snp.genotype_quality as i32,
                    snp.depth,
                    snp.allele_freqs[1]
                );
                rd.format = "GT:GQ:DP:AF".to_string().into_bytes();
                records.push(rd);
                continue;
            }

            if snp.dense {
                let mut rd: VCFRecord = VCFRecord::default();
                rd.chromosome = snp.chromosome.clone();
                rd.position = snp.pos as u64 + 1; // position in vcf format is 1-based
                rd.id = vec!['.' as u8];
                rd.reference = vec![snp.reference as u8];
                if snp.variant_type == 1 || snp.variant_type == 2 {
                    if snp.alleles[0] != snp.reference {
                        rd.alternative = vec![vec![snp.alleles[0] as u8]];
                    } else if snp.alleles[1] != snp.reference {
                        rd.alternative = vec![vec![snp.alleles[1] as u8]];
                    }
                } else if snp.variant_type == 3 {
                    rd.alternative = vec![vec![snp.alleles[0] as u8], vec![snp.alleles[1] as u8]];
                }
                rd.qual = snp.variant_quality as i32;
                rd.filter = "dn".to_string().into_bytes();
                rd.info = format!("RDS={}", "dense_snp").to_string().into_bytes();
                let gt;
                if snp.variant_type == 1 {
                    gt = "0/1"
                } else if snp.variant_type == 2 {
                    gt = "1/1"
                } else if snp.variant_type == 3 {
                    gt = "1/2"
                } else {
                    continue;
                }
                rd.genotype = format!(
                    "{}:{}:{}:{:.2}",
                    gt,
                    snp.genotype_quality as i32,
                    snp.depth,
                    snp.allele_freqs[1]
                );
                rd.format = "GT:GQ:DP:AF".to_string().into_bytes();
                records.push(rd);
                continue;
            }

            if snp.variant_type == 0 {
                continue;
            }

            if snp.variant_type == 1 {
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
                if snp.variant_quality < min_qual_for_candidate as f64 {
                    rd.filter = "LowQual".to_string().into_bytes();
                } else {
                    rd.filter = "PASS".to_string().into_bytes();
                }
                if snp.single {
                    rd.info = format!("RDS={}", "single_snp").to_string().into_bytes();
                } else {
                    rd.info = "RDS=.".to_string().into_bytes();
                }
                let gt;
                if snp.germline {
                    if snp.haplotype == -1 {
                        gt = "0|1"
                    } else if snp.haplotype == 1 {
                        gt = "1|0"
                    } else {
                        panic!("Error: unknown haplotype: {:?}", snp);
                    }
                } else {
                    // TODO: som var?
                    continue;
                }
                let mut af = 0.0;
                if snp.alleles[0] == snp.reference {
                    af = snp.allele_freqs[1];
                } else if snp.alleles[1] == snp.reference {
                    af = snp.allele_freqs[0];
                } else {
                    panic!("Error: unexpected allele. ref: {}, alt1: {}, alt2: {}, {}:{}", snp.reference, snp.alleles[0], snp.alleles[1], String::from_utf8_lossy(&snp.chromosome), snp.pos);
                }
                if snp.phase_set != 0 {
                    rd.genotype = format!(
                        "{}:{}:{}:{}:{:.2}:{:.2}:{},{},{},{}",
                        gt,
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
                    rd.format = "GT:PS:GQ:DP:AF:PQ:AE".to_string().into_bytes();
                } else {
                    rd.genotype = format!(
                        "{}:{}:{}:{:.2}:{:.2}:{},{},{},{}",
                        gt,
                        snp.genotype_quality as i32,
                        snp.depth,
                        af,
                        snp.phase_score,
                        snp.haplotype_expression[0],
                        snp.haplotype_expression[1],
                        snp.haplotype_expression[2],
                        snp.haplotype_expression[3]
                    );
                    rd.format = "GT:GQ:DP:AF:PQ:AE".to_string().into_bytes();
                }
                records.push(rd);
                continue;
            }

            if snp.variant_type == 2 || snp.variant_type == 3 {
                let mut rd: VCFRecord = VCFRecord::default();
                rd.chromosome = snp.chromosome.clone();
                rd.position = snp.pos as u64 + 1; // position in vcf format is 1-based
                rd.id = vec!['.' as u8];
                rd.reference = vec![snp.reference as u8];
                if snp.variant_type == 2 {
                    rd.alternative = vec![vec![snp.alleles[0] as u8]];
                } else if snp.variant_type == 3 {
                    rd.alternative = vec![vec![snp.alleles[0] as u8], vec![snp.alleles[1] as u8]];
                }
                // if snp.alleles[0] != snp.reference {
                //     rd.alternative = vec![vec![snp.alleles[0] as u8]];
                // } else if snp.alleles[1] != snp.reference {
                //     rd.alternative = vec![vec![snp.alleles[1] as u8]];
                // }
                rd.qual = snp.variant_quality as i32;
                if snp.variant_quality < min_qual_for_candidate as f64 {
                    rd.filter = "LowQual".to_string().into_bytes();
                } else {
                    rd.filter = "PASS".to_string().into_bytes();
                }
                rd.info = "RDS=.".to_string().into_bytes();
                let gt;
                if snp.germline {
                    if snp.variant_type == 2 {
                        gt = "1/1"
                    } else if snp.variant_type == 3 {
                        gt = "1/2"
                    } else {
                        panic!("Error: unknown haplotype: {:?}", snp);
                    }
                } else {
                    // TODO: som var?
                    continue;
                }
                rd.genotype = format!(
                    "{}:{}:{}:{:.2},{:.2}",
                    gt,
                    snp.genotype_quality as i32,
                    snp.depth,
                    snp.allele_freqs[0],
                    snp.allele_freqs[1]
                );
                rd.format = "GT:GQ:DP:AF".to_string().into_bytes();
                records.push(rd);
                continue;
            }

            // if snp.variant_type == 1 {
            //     let mut is_dense = false;
            //     let mut is_rna_edit = false;
            //     let mut is_single_snp = false;
            //     let mut is_unconfident_phased_snp = false;
            //     let mut is_ase_snp = false;
            //
            //     // dense snp?
            //     if snp.filter == true {
            //         is_dense = true;
            //     }
            //
            //     // rna edit?
            //     if snp.rna_editing == true {
            //         is_rna_edit = true;
            //     }
            //
            //     // single snp?
            //     if snp.single == true {
            //         is_single_snp = true;
            //     }
            //
            //     // unconfident phase?
            //     if !is_dense && !is_single_snp && snp.phase_score == 0.0 {
            //         is_unconfident_phased_snp = true;
            //     }
            //
            //     // ase snp?
            //     if snp.ase == true {
            //         is_ase_snp = true;
            //     }
            //
            //     let mut rd: VCFRecord = VCFRecord::default();
            //     rd.chromosome = snp.chromosome.clone();
            //     rd.position = snp.pos as u64 + 1; // position in vcf format is 1-based
            //     rd.id = vec!['.' as u8];
            //     rd.reference = vec![snp.reference as u8];
            //     if snp.alleles[0] != snp.reference {
            //         rd.alternative = vec![vec![snp.alleles[0] as u8]];
            //     } else if snp.alleles[1] != snp.reference {
            //         rd.alternative = vec![vec![snp.alleles[1] as u8]];
            //     }
            //     rd.qual = snp.variant_quality as i32;
            //     if is_dense {
            //         rd.filter = "dn".to_string().into_bytes();
            //         rd.info = format!("RDS={}", "dense_snp").to_string().into_bytes();
            //     }
            //     if is_single_snp {
            //         if snp.variant_quality < min_qual_for_singlesnp_rnaedit as f64 {
            //             rd.filter = "LowQual".to_string().into_bytes();
            //             rd.info = format!("RDS={}", "single_snp").to_string().into_bytes();
            //         } else {
            //             if (rd.reference[0] == b'A' && rd.alternative[0][0] == b'G') || (rd.reference[0] == b'T' && rd.alternative[0][0] == b'C') {
            //                 continue;
            //             }
            //             rd.filter = "PASS".to_string().into_bytes();
            //             rd.info = format!("RDS={}", "single_snp").to_string().into_bytes();
            //         }
            //     }
            //     if is_ase_snp {
            //         if snp.phase_score < ase_ps_cutoff as f64 {
            //             rd.filter = "LowQual".to_string().into_bytes();
            //             if is_rna_edit {
            //                 rd.info = format!("RDS={}", "ase_snp,rna_editing").to_string().into_bytes();    // rna editing snp has bad phase
            //             } else {
            //                 rd.info = format!("RDS={}", "ase_snp").to_string().into_bytes();
            //             }
            //         } else {
            //             rd.filter = "PASS".to_string().into_bytes();
            //             rd.info = format!("RDS={}", "ase_snp").to_string().into_bytes();
            //         }
            //     }
            //     if !is_dense && !is_single_snp && !is_ase_snp {
            //         if snp.variant_quality < min_qual_for_candidate as f64 {
            //             rd.filter = "LowQual".to_string().into_bytes();
            //             rd.info = "RDS=.".to_string().into_bytes();
            //         } else if snp.phase_score < min_phase_score as f64 {
            //             // bad phase, unconfindent phase
            //             rd.filter = "LowQual".to_string().into_bytes();
            //             if is_rna_edit {
            //                 rd.info = format!("RDS={}", "rna_editing").to_string().into_bytes();    // rna editing snp has bad phase
            //             } else {
            //                 rd.info = "RDS=.".to_string().into_bytes();
            //             }
            //         } else {
            //             // good phase, germline snp
            //             rd.filter = "PASS".to_string().into_bytes();
            //             rd.info = "RDS=.".to_string().into_bytes();
            //         }
            //     }
            //
            //     if is_dense || is_single_snp || is_unconfident_phased_snp || snp.phase_score == 0.0 || snp.haplotype == 0 {
            //         rd.genotype = format!(
            //             "{}:{}:{}:{:.2}",
            //             "0/1",
            //             snp.genotype_quality as i32,
            //             snp.depth,
            //             snp.allele_freqs[1]
            //         );
            //         rd.format = "GT:GQ:DP:AF".to_string().into_bytes();
            //     } else {
            //         let mut af = 0.0;
            //         if snp.alleles[0] == snp.reference {
            //             af = snp.allele_freqs[1];
            //         } else if snp.alleles[1] == snp.reference {
            //             af = snp.allele_freqs[0];
            //         } else {
            //             panic!("Error: unexpected allele. ref: {}, alt1: {}, alt2: {}, {}:{}", snp.reference, snp.alleles[0], snp.alleles[1], String::from_utf8_lossy(&snp.chromosome), snp.pos);
            //         }
            //         if snp.phase_set != 0 {
            //             if snp.haplotype == -1 {
            //                 rd.genotype = format!(
            //                     "{}:{}:{}:{}:{:.2}:{:.2}:{},{},{},{}",
            //                     "0|1",
            //                     snp.phase_set,
            //                     snp.genotype_quality as i32,
            //                     snp.depth,
            //                     af,
            //                     snp.phase_score,
            //                     snp.haplotype_expression[0],
            //                     snp.haplotype_expression[1],
            //                     snp.haplotype_expression[2],
            //                     snp.haplotype_expression[3]
            //                 );
            //             } else if snp.haplotype == 1 {
            //                 rd.genotype = format!(
            //                     "{}:{}:{}:{}:{:.2}:{:.2}:{},{},{},{}",
            //                     "1|0",
            //                     snp.phase_set,
            //                     snp.genotype_quality as i32,
            //                     snp.depth,
            //                     af,
            //                     snp.phase_score,
            //                     snp.haplotype_expression[0],
            //                     snp.haplotype_expression[1],
            //                     snp.haplotype_expression[2],
            //                     snp.haplotype_expression[3]
            //                 );
            //             }
            //             rd.format = "GT:PS:GQ:DP:AF:PQ:AE".to_string().into_bytes();
            //         } else {
            //             if snp.haplotype == -1 {
            //                 rd.genotype = format!(
            //                     "{}:{}:{}:{:.2}:{:.2}:{},{},{},{}",
            //                     "0|1",
            //                     snp.genotype_quality as i32,
            //                     snp.depth,
            //                     af,
            //                     snp.phase_score,
            //                     snp.haplotype_expression[0],
            //                     snp.haplotype_expression[1],
            //                     snp.haplotype_expression[2],
            //                     snp.haplotype_expression[3]
            //                 );
            //             } else if snp.haplotype == 1 {
            //                 rd.genotype = format!(
            //                     "{}:{}:{}:{:.2}:{:.2}:{},{},{},{}",
            //                     "1|0",
            //                     snp.genotype_quality as i32,
            //                     snp.depth,
            //                     af,
            //                     snp.phase_score,
            //                     snp.haplotype_expression[0],
            //                     snp.haplotype_expression[1],
            //                     snp.haplotype_expression[2],
            //                     snp.haplotype_expression[3]
            //                 );
            //             }
            //             rd.format = "GT:GQ:DP:AF:PQ:AE".to_string().into_bytes();
            //         }
            //     }
            //     records.push(rd);
            // }
            // if snp.variant_type == 2 {
            //     let mut is_dense = false;
            //
            //     if snp.filter == true {
            //         is_dense = true;
            //     }
            //
            //     let mut rd: VCFRecord = VCFRecord::default();
            //     rd.chromosome = snp.chromosome.clone();
            //     rd.position = snp.pos as u64 + 1; // position in vcf format is 1-based
            //     rd.id = vec!['.' as u8];
            //     rd.reference = vec![snp.reference as u8];
            //     rd.alternative = vec![vec![snp.alleles[0] as u8]];
            //     rd.qual = snp.variant_quality as i32;
            //
            //     if is_dense {
            //         rd.filter = "dn".to_string().into_bytes();
            //         rd.info = format!("RDS={}", "dense_snp").to_string().into_bytes();
            //     }
            //     if !is_dense {
            //         if snp.variant_quality < min_qual_for_candidate as f64 {
            //             rd.filter = "LowQual".to_string().into_bytes();
            //             rd.info = "RDS=.".to_string().into_bytes();
            //         } else {
            //             rd.filter = "PASS".to_string().into_bytes();
            //             rd.info = "RDS=.".to_string().into_bytes();
            //         }
            //     }
            //     rd.genotype = format!(
            //         "{}:{}:{}:{:.2}",
            //         "1/1",
            //         snp.genotype_quality as i32,
            //         snp.depth,
            //         snp.allele_freqs[0]
            //     );
            //     rd.format = "GT:GQ:DP:AF".to_string().into_bytes();
            //     records.push(rd);
            // }
            // if snp.variant_type == 3 {
            //     let mut is_dense = false;
            //
            //     if snp.filter == true {
            //         is_dense = true;
            //     }
            //
            //     let mut rd: VCFRecord = VCFRecord::default();
            //     rd.chromosome = snp.chromosome.clone();
            //     rd.position = snp.pos as u64 + 1; // position in vcf format is 1-based
            //     rd.id = vec!['.' as u8];
            //     rd.reference = vec![snp.reference as u8];
            //     rd.alternative = vec![vec![snp.alleles[0] as u8], vec![snp.alleles[1] as u8]];
            //     rd.qual = snp.variant_quality as i32;
            //
            //     if is_dense {
            //         rd.filter = "dn".to_string().into_bytes();
            //         rd.info = format!("RDS={}", "dense_snp").to_string().into_bytes();
            //     }
            //
            //     if !is_dense {
            //         if snp.variant_quality < min_qual_for_candidate as f64 {
            //             rd.filter = "LowQual".to_string().into_bytes();
            //             rd.info = "RDS=.".to_string().into_bytes();
            //         } else {
            //             rd.filter = "PASS".to_string().into_bytes();
            //             rd.info = "RDS=.".to_string().into_bytes();
            //         }
            //     }
            //     rd.genotype = format!(
            //         "{}:{}:{}:{:.2},{:.2}",
            //         "1/2",
            //         snp.genotype_quality as i32,
            //         snp.depth,
            //         snp.allele_freqs[0],
            //         snp.allele_freqs[1]
            //     );
            //     rd.format = "GT:GQ:DP:AF".to_string().into_bytes();
            //     records.push(rd);
            // }
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
            if snp.somatic || snp.rna_editing {
                continue;
            }
            if snp.dense == true {
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

