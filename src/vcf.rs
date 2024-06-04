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
}


impl SNPFrag {
    pub fn output_phased_vcf(&mut self, min_phase_score: f32) -> Vec<VCFRecord> {
        let mut records: Vec<VCFRecord> = Vec::new();
        for i in 0..self.candidate_snps.len() {
            let snp = &self.candidate_snps[i];
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
                let mut gt = "0/0";
                let mut af = 0.0;
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
                    af
                );
                rd.format = "GT:GQ:DP:AF".to_string().into_bytes();
                records.push(rd);
                continue;
            }

            if snp.non_selected {
                let mut rd: VCFRecord = VCFRecord::default();
                rd.chromosome = snp.chromosome.clone();
                rd.position = snp.pos as u64 + 1; // position in vcf format is 1-based
                rd.id = vec!['.' as u8];
                rd.reference = vec![snp.reference as u8];
                let mut gt = "0/0";
                let mut af = 0.0;
                if snp.rna_editing {
                    if snp.variant_type == 0 {
                        continue;
                    } else if snp.variant_type == 1 {
                        if snp.alleles[0] != snp.reference {
                            rd.alternative = vec![vec![snp.alleles[0] as u8]];
                            af = snp.allele_freqs[0];
                        } else if snp.alleles[1] != snp.reference {
                            rd.alternative = vec![vec![snp.alleles[1] as u8]];
                            af = snp.allele_freqs[1];
                        }
                        gt = "0/1";
                        rd.filter = "RnaEdit".to_string().into_bytes();
                    } else if snp.variant_type == 2 {
                        if snp.alleles[0] != snp.reference {
                            rd.alternative = vec![vec![snp.alleles[0] as u8]];
                            af = snp.allele_freqs[0];
                        } else if snp.alleles[1] != snp.reference {
                            rd.alternative = vec![vec![snp.alleles[1] as u8]];
                            af = snp.allele_freqs[1];
                        }
                        gt = "1/1";
                        rd.filter = "RnaEdit".to_string().into_bytes();
                    } else {
                        continue;
                    }
                    rd.qual = snp.variant_quality as i32;
                    rd.info = "RDS=noselect".to_string().into_bytes();
                    rd.genotype = format!("{}:{}:{}:{:.2}", gt, snp.genotype_quality as i32, snp.depth, af);
                    rd.format = "GT:GQ:DP:AF".to_string().into_bytes();
                    records.push(rd);
                    continue;
                }
                if snp.variant_type == 0 {
                    if snp.alleles[0] != snp.reference {
                        rd.alternative = vec![vec![snp.alleles[0] as u8]];
                        af = snp.allele_freqs[0];
                    } else if snp.alleles[1] != snp.reference {
                        rd.alternative = vec![vec![snp.alleles[1] as u8]];
                        af = snp.allele_freqs[1];
                    } else {
                        continue;
                    }
                    gt = "0/0";
                    rd.filter = "HomRef".to_string().into_bytes();
                } else if snp.variant_type == 1 {
                    // continue;
                    if snp.alleles[0] != snp.reference {
                        rd.alternative = vec![vec![snp.alleles[0] as u8]];
                        af = snp.allele_freqs[0];
                    } else if snp.alleles[1] != snp.reference {
                        rd.alternative = vec![vec![snp.alleles[1] as u8]];
                        af = snp.allele_freqs[1];
                    }
                    gt = "0/1";
                    rd.filter = "LowQual".to_string().into_bytes();
                } else if snp.variant_type == 2 {
                    if snp.alleles[0] != snp.reference {
                        rd.alternative = vec![vec![snp.alleles[0] as u8]];
                        af = snp.allele_freqs[0];
                    } else if snp.alleles[1] != snp.reference {
                        rd.alternative = vec![vec![snp.alleles[1] as u8]];
                        af = snp.allele_freqs[1];
                    }
                    gt = "1/1";
                    rd.filter = "PASS".to_string().into_bytes();
                } else {
                    continue;
                }
                rd.qual = snp.variant_quality as i32;
                rd.info = "RDS=noselect".to_string().into_bytes();
                rd.genotype = format!("{}:{}:{}:{:.2}", gt, snp.genotype_quality as i32, snp.depth, af);
                rd.format = "GT:GQ:DP:AF".to_string().into_bytes();
                records.push(rd);
            } else {
                let mut rd: VCFRecord = VCFRecord::default();
                rd.chromosome = snp.chromosome.clone();
                rd.position = snp.pos as u64 + 1; // position in vcf format is 1-based
                rd.id = vec!['.' as u8];
                rd.reference = vec![snp.reference as u8];
                let mut gt = "0/0";
                let mut af = 0.0;
                if snp.phase_score >= min_phase_score as f64 {
                    if snp.variant_type == 1 {
                        if snp.alleles[0] != snp.reference {
                            rd.alternative = vec![vec![snp.alleles[0] as u8]];
                            af = snp.allele_freqs[0];
                        } else if snp.alleles[1] != snp.reference {
                            rd.alternative = vec![vec![snp.alleles[1] as u8]];
                            af = snp.allele_freqs[1];
                        }
                        if snp.haplotype == 1 {
                            gt = "0|1"
                        } else {
                            gt = "1|0"
                        }
                        rd.filter = "PASS".to_string().into_bytes();
                    }
                } else {
                    if snp.variant_type == 0 {
                        if snp.alleles[0] != snp.reference {
                            rd.alternative = vec![vec![snp.alleles[0] as u8]];
                            af = snp.allele_freqs[0];
                        } else if snp.alleles[1] != snp.reference {
                            rd.alternative = vec![vec![snp.alleles[1] as u8]];
                            af = snp.allele_freqs[1];
                        } else {
                            continue;
                        }
                        gt = "0/0";
                        rd.filter = "HomRef".to_string().into_bytes();
                    } else if snp.variant_type == 1 {
                        if snp.alleles[0] != snp.reference {
                            rd.alternative = vec![vec![snp.alleles[0] as u8]];
                            af = snp.allele_freqs[0];
                        } else if snp.alleles[1] != snp.reference {
                            rd.alternative = vec![vec![snp.alleles[1] as u8]];
                            af = snp.allele_freqs[1];
                        }
                        gt = "0/1";
                        rd.filter = "LowQual".to_string().into_bytes();
                    } else if snp.variant_type == 2 {
                        if snp.alleles[0] != snp.reference {
                            rd.alternative = vec![vec![snp.alleles[0] as u8]];
                            af = snp.allele_freqs[0];
                        } else if snp.alleles[1] != snp.reference {
                            rd.alternative = vec![vec![snp.alleles[1] as u8]];
                            af = snp.allele_freqs[1];
                        }
                        gt = "1/1";
                        rd.filter = "PASS".to_string().into_bytes();
                    } else {
                        continue;
                    }
                }
                rd.qual = snp.variant_quality as i32;
                rd.info = "RDS=select".to_string().into_bytes();
                if snp.phase_set != 0 {
                    rd.genotype = format!("{}:{}:{}:{}:{:.2}:{:.2}", gt, snp.genotype_quality as i32, snp.phase_set, snp.depth, af, snp.phase_score);
                } else {
                    rd.genotype = format!("{}:{}:{}:{}:{:.2}:{:.2}", gt, snp.genotype_quality as i32, ".", snp.depth, af, snp.phase_score);
                }
                rd.format = "GT:GQ:PS:DP:AF:PQ".to_string().into_bytes();
                records.push(rd);
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
            // if snp.somatic || snp.rna_editing {
            //     continue;
            // }
            if snp.cand_somatic {
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
                            "1/2", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[0], snp.allele_freqs[1]
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
                        "1/2", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[0], snp.allele_freqs[1]
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
                    "1/2", snp.genotype_quality as i32, snp.depth, snp.allele_freqs[0], snp.allele_freqs[1]
                );
                if snp.variant_quality < min_qual as f64 {
                    rd.filter = "LowQual".to_string().into_bytes();
                } else {
                    rd.filter = "PASS".to_string().into_bytes();
                }
                rd.info = "RDS=.".to_string().into_bytes();
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

