use noodles_vcf::variant::record::samples::keys::key;
use noodles_vcf::variant::record::samples::series::Value;
use noodles_vcf::variant::record::samples::Sample;
use noodles_vcf::variant::record::{AlternateBases, Filters};
use rust_htslib::bcf;
use rust_htslib::bcf::record::GenotypeAllele;
use rust_htslib::bcf::Read;
use std::collections::{HashMap, HashSet};

use crate::snpfrags::SNPFrag;

#[derive(Debug, Default, Clone)]
pub struct VCFRecord {
    pub chromosome: String,
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
                let mut af = [0.0, 0.0];
                if snp.variant_type == 1 || snp.variant_type == 2 {
                    if snp.alleles[0] != snp.reference {
                        rd.alternative = vec![vec![snp.alleles[0] as u8]];
                        af[0] = snp.allele_freqs[0];
                    } else if snp.alleles[1] != snp.reference {
                        rd.alternative = vec![vec![snp.alleles[1] as u8]];
                        af[0] = snp.allele_freqs[1];
                    }
                } else if snp.variant_type == 3 {
                    rd.alternative = vec![vec![snp.alleles[0] as u8], vec![snp.alleles[1] as u8]];
                    af[0] = snp.allele_freqs[0];
                    af[1] = snp.allele_freqs[1];
                }
                rd.qual = snp.variant_quality as i32;
                rd.filter = "dn".to_string().into_bytes();
                rd.info = format!("RDS={}", "dense_snp").to_string().into_bytes();
                let mut gt = "0/0";
                if snp.variant_type == 1 {
                    gt = "0/1"
                } else if snp.variant_type == 2 {
                    gt = "1/1"
                } else if snp.variant_type == 3 {
                    gt = "1/2"
                } else {
                    continue;
                }
                if snp.variant_type == 3 {
                    rd.genotype = format!(
                        "{}:{}:{}:{:.2},{:.2}",
                        gt, snp.genotype_quality as i32, snp.depth, af[0], af[1]
                    );
                } else {
                    rd.genotype = format!(
                        "{}:{}:{}:{:.2}",
                        gt, snp.genotype_quality as i32, snp.depth, af[0]
                    );
                }
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
                let mut af = [0.0, 0.0];
                if snp.rna_editing {
                    if snp.variant_type == 1 || snp.variant_type == 2 {
                        if snp.alleles[0] != snp.reference {
                            rd.alternative = vec![vec![snp.alleles[0] as u8]];
                            af[0] = snp.allele_freqs[0];
                        } else if snp.alleles[1] != snp.reference {
                            rd.alternative = vec![vec![snp.alleles[1] as u8]];
                            af[0] = snp.allele_freqs[1];
                        }
                    } else {
                        continue;
                    }
                    rd.filter = "RnaEdit".to_string().into_bytes();
                    rd.qual = snp.variant_quality as i32;
                    rd.info = "RDS=noselect".to_string().into_bytes();
                    if snp.variant_type == 1 {
                        gt = "0/1"
                    } else if snp.variant_type == 2 {
                        gt = "1/1"
                    }
                    rd.genotype = format!(
                        "{}:{}:{}:{:.2}",
                        gt, snp.genotype_quality as i32, snp.depth, af[0]
                    );
                    rd.format = "GT:GQ:DP:AF".to_string().into_bytes();
                    records.push(rd);
                    continue;
                }
                if snp.variant_type == 0 || snp.variant_type == 1 || snp.variant_type == 2 {
                    if snp.alleles[0] != snp.reference {
                        rd.alternative = vec![vec![snp.alleles[0] as u8]];
                        af[0] = snp.allele_freqs[0];
                    } else if snp.alleles[1] != snp.reference {
                        rd.alternative = vec![vec![snp.alleles[1] as u8]];
                        af[0] = snp.allele_freqs[1];
                    }
                    if snp.variant_type == 0 {
                        gt = "0/0";
                        rd.filter = "HomRef".to_string().into_bytes();
                    } else if snp.variant_type == 1 {
                        gt = "0/1";
                        rd.filter = "LowQual".to_string().into_bytes();
                    } else {
                        gt = "1/1";
                        rd.filter = "PASS".to_string().into_bytes();
                    }
                } else {
                    // continue;
                    if snp.genotype == -1 || snp.genotype == 1 {
                        if snp.alleles[0] != snp.reference {
                            rd.alternative = vec![vec![snp.alleles[0] as u8]];
                            af[0] = snp.allele_freqs[0];
                        } else if snp.alleles[1] != snp.reference {
                            rd.alternative = vec![vec![snp.alleles[1] as u8]];
                            af[0] = snp.allele_freqs[1];
                        }
                        if snp.genotype == -1 {
                            gt = "1/1";
                            rd.filter = "PASS".to_string().into_bytes();
                        } else {
                            gt = "0/0";
                            rd.filter = "HomRef".to_string().into_bytes();
                        }
                    } else if snp.genotype == 0 {
                        rd.alternative =
                            vec![vec![snp.alleles[0] as u8], vec![snp.alleles[1] as u8]];
                        af[0] = snp.allele_freqs[0];
                        af[1] = snp.allele_freqs[1];
                        gt = "1/2";
                        rd.filter = "Multiallelic".to_string().into_bytes();
                    }
                }
                rd.qual = snp.variant_quality as i32;
                rd.info = "RDS=noselect".to_string().into_bytes();
                if gt == "0/0" || gt == "0/1" || gt == "1/1" {
                    rd.genotype = format!(
                        "{}:{}:{}:{:.2}",
                        gt, snp.genotype_quality as i32, snp.depth, af[0]
                    );
                } else {
                    rd.genotype = format!(
                        "{}:{}:{}:{:.2},{:.2}",
                        gt, snp.genotype_quality as i32, snp.depth, af[0], af[1]
                    );
                }
                rd.format = "GT:GQ:DP:AF".to_string().into_bytes();
                records.push(rd);
            } else {
                let mut rd: VCFRecord = VCFRecord::default();
                rd.chromosome = snp.chromosome.clone();
                rd.position = snp.pos as u64 + 1; // position in vcf format is 1-based
                rd.id = vec!['.' as u8];
                rd.reference = vec![snp.reference as u8];
                let mut gt = "0/0";
                let mut af = [0.0, 0.0];
                if snp.phase_score >= min_phase_score as f64 {
                    if snp.variant_type == 1 {
                        if snp.alleles[0] != snp.reference {
                            rd.alternative = vec![vec![snp.alleles[0] as u8]];
                            af[0] = snp.allele_freqs[0];
                        } else if snp.alleles[1] != snp.reference {
                            rd.alternative = vec![vec![snp.alleles[1] as u8]];
                            af[0] = snp.allele_freqs[1];
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
                            af[0] = snp.allele_freqs[0];
                        } else if snp.alleles[1] != snp.reference {
                            rd.alternative = vec![vec![snp.alleles[1] as u8]];
                            af[0] = snp.allele_freqs[1];
                        }
                        gt = "0/0";
                        rd.filter = "HomRef".to_string().into_bytes();
                    } else if snp.variant_type == 1 {
                        if snp.alleles[0] != snp.reference {
                            rd.alternative = vec![vec![snp.alleles[0] as u8]];
                            af[0] = snp.allele_freqs[0];
                        } else if snp.alleles[1] != snp.reference {
                            rd.alternative = vec![vec![snp.alleles[1] as u8]];
                            af[0] = snp.allele_freqs[1];
                        }
                        gt = "0/1";
                        rd.filter = "LowQual".to_string().into_bytes();
                    } else if snp.variant_type == 2 {
                        if snp.alleles[0] != snp.reference {
                            rd.alternative = vec![vec![snp.alleles[0] as u8]];
                            af[0] = snp.allele_freqs[0];
                        } else if snp.alleles[1] != snp.reference {
                            rd.alternative = vec![vec![snp.alleles[1] as u8]];
                            af[0] = snp.allele_freqs[1];
                        }
                        gt = "1/1";
                        rd.filter = "PASS".to_string().into_bytes();
                    } else {
                        // continue;
                        if snp.genotype == -1 || snp.genotype == 1 {
                            if snp.alleles[0] != snp.reference {
                                rd.alternative = vec![vec![snp.alleles[0] as u8]];
                                af[0] = snp.allele_freqs[0];
                            } else if snp.alleles[1] != snp.reference {
                                rd.alternative = vec![vec![snp.alleles[1] as u8]];
                                af[0] = snp.allele_freqs[1];
                            }
                            if snp.genotype == -1 {
                                gt = "1/1";
                                rd.filter = "PASS".to_string().into_bytes();
                            } else {
                                gt = "0/0";
                                rd.filter = "HomRef".to_string().into_bytes();
                            }
                        } else if snp.genotype == 0 {
                            rd.alternative =
                                vec![vec![snp.alleles[0] as u8], vec![snp.alleles[1] as u8]];
                            af[0] = snp.allele_freqs[0];
                            af[1] = snp.allele_freqs[1];
                            gt = "1/2";
                            rd.filter = "Multiallelic".to_string().into_bytes();
                        }
                    }
                }
                rd.qual = snp.variant_quality as i32;
                rd.info = "RDS=select".to_string().into_bytes();
                if gt == "0/0" || gt == "0/1" || gt == "1/1" || gt == "0|1" || gt == "1|0" {
                    if snp.phase_set != 0 {
                        rd.genotype = format!(
                            "{}:{}:{}:{}:{:.2}:{:.2}",
                            gt,
                            snp.genotype_quality as i32,
                            snp.phase_set,
                            snp.depth,
                            af[0],
                            snp.phase_score
                        );
                    } else {
                        rd.genotype = format!(
                            "{}:{}:{}:{}:{:.2}:{:.2}",
                            gt, snp.genotype_quality as i32, ".", snp.depth, af[0], snp.phase_score
                        );
                    }
                } else {
                    if snp.phase_set != 0 {
                        rd.genotype = format!(
                            "{}:{}:{}:{}:{:.2},{:.2}:{:.2}",
                            gt,
                            snp.genotype_quality as i32,
                            snp.phase_set,
                            snp.depth,
                            af[0],
                            af[1],
                            snp.phase_score
                        );
                    } else {
                        rd.genotype = format!(
                            "{}:{}:{}:{}:{:.2},{:.2}:{:.2}",
                            gt,
                            snp.genotype_quality as i32,
                            ".",
                            snp.depth,
                            af[0],
                            af[1],
                            snp.phase_score
                        );
                    }
                }
                rd.format = "GT:GQ:PS:DP:AF:PQ".to_string().into_bytes();
                records.push(rd);
            }
        }
        records
    }
}

pub fn load_vcf(
    vcf_path: &str,
) -> (
    HashMap<String, HashSet<usize>>,
    HashMap<String, HashMap<usize, (u8, f32)>>,
) {
    let mut input_candidates: HashMap<String, HashSet<usize>> = HashMap::new();
    let mut input_candidates_genotype_qual: HashMap<String, HashMap<usize, (u8, f32)>> =
        HashMap::new();
    let mut reader = noodles_vcf::io::reader::Builder::default()
        .build_from_path(vcf_path)
        .unwrap();
    let header = reader.read_header().unwrap();
    for result in reader.records() {
        let record = result.unwrap();
        if record.filters().is_pass(&header).unwrap() && record.reference_bases().len() == 1 {
            for alt in record.alternate_bases().iter() {
                if alt.unwrap().len() != 1 {
                    continue;
                }
            }
            let chr = record.reference_sequence_name();
            let pos = record.variant_start().unwrap().unwrap().get();
            let qual = record.quality_score().unwrap().unwrap();
            // store position by chromosome
            let samples = record.samples();
            for (_, sample) in samples.iter().enumerate() {
                let gt = sample
                    .get(&header, key::GENOTYPE)
                    .unwrap()
                    .unwrap()
                    .unwrap();
                match gt {
                    Value::Genotype(genotype) => {
                        let mut gt_vec = Vec::new();
                        for vi in genotype.iter() {
                            let (allele, _) = vi.unwrap();
                            let allele = allele.unwrap();
                            gt_vec.push(allele);
                        }
                        if gt_vec.len() != 2 {
                            continue;
                        }
                        if gt_vec[0] + gt_vec[1] == 1 {
                            // gt: 0/1
                            input_candidates_genotype_qual
                                .entry(chr.to_string())
                                .or_insert(HashMap::new())
                                .insert(pos, (1, qual));
                            input_candidates
                                .entry(chr.to_string())
                                .or_insert(HashSet::new())
                                .insert(pos);
                        } else if gt_vec[0] + gt_vec[1] == 2 {
                            // gt: 1/1
                            input_candidates_genotype_qual
                                .entry(chr.to_string())
                                .or_insert(HashMap::new())
                                .insert(pos, (2, qual));
                            input_candidates
                                .entry(chr.to_string())
                                .or_insert(HashSet::new())
                                .insert(pos);
                        } else if gt_vec[0] + gt_vec[1] == 3 {
                            // gt: 1/2
                            input_candidates_genotype_qual
                                .entry(chr.to_string())
                                .or_insert(HashMap::new())
                                .insert(pos, (3, qual));
                            input_candidates
                                .entry(chr.to_string())
                                .or_insert(HashSet::new())
                                .insert(pos);
                        }
                    }
                    _ => {}
                }
            }
        }
    }

    (input_candidates, input_candidates_genotype_qual)
}

#[derive(Debug)]
pub struct GenotypeAndQuality {
    pub genotype: u8, // 1 for 0/1, 2 for 1/1, 3 for 1/2
    pub(crate) quality: f32,
    pub phased: bool,
}

pub fn get_genotype_quality_phase_from_vcf(
    vcf_path: &str,
) -> HashMap<String, HashMap<usize, GenotypeAndQuality>> {
    let mut genotype_quality: HashMap<String, HashMap<usize, GenotypeAndQuality>> = HashMap::new();
    let mut reader = bcf::Reader::from_path(vcf_path).unwrap();
    for result in reader.records() {
        let record = result.unwrap();
        let rid = record.rid().unwrap();
        let chr = std::str::from_utf8(record.header().rid2name(rid).unwrap()).unwrap();
        let pos = record.pos(); // 0-based position
        let qual = record.qual();
        let mut phased = false;
        let sample_count = usize::try_from(record.sample_count()).unwrap();
        let gts = record.genotypes().expect("Error reading genotypes");
        for sample_index in 0..sample_count {
            let gt = gts.get(sample_index);
            if gt.len() != 2 {
                // Skip records that don't have two alleles.
                continue;
            }
            let gt0 = match gt[0] {
                GenotypeAllele::Unphased(n) | GenotypeAllele::Phased(n) => n,
                GenotypeAllele::UnphasedMissing | GenotypeAllele::PhasedMissing => 3,
            };
            // For gt[1], update `phased` if the allele is phased.
            let gt1 = match gt[1] {
                GenotypeAllele::Unphased(n) => n,
                GenotypeAllele::Phased(n) => {
                    phased = true;
                    n
                }
                GenotypeAllele::UnphasedMissing => 3,
                GenotypeAllele::PhasedMissing => {
                    phased = true;
                    3
                }
            };

            let genotype = match (gt0, gt1) {
                (0, 0) => 0,          // 0/0
                (0, 1) | (1, 0) => 1, // 0/1 or 1/0
                (1, 1) => 2,          // 1/1
                (1, 2) | (2, 1) => 3, // 1/2 or 2/1
                _ => 4,               // Other cases
            };

            genotype_quality.entry(chr.to_string()).or_default().insert(
                pos as usize,
                GenotypeAndQuality {
                    genotype,
                    quality: qual,
                    phased,
                },
            );

            // println!(
            //     "{}:{} Qual {} Sample {} phase: {} genotype: {:?}, gt: {:?}",
            //     chr, pos, qual, sample_index, phased, genotype, gt
            // );
        }
    }
    genotype_quality
}
