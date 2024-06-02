use std::collections::HashMap;
use std::hash::Hash;

use probability::distribution::Distribution;

use crate::exon::Exon;

#[derive(Debug, Clone, Default)]
pub struct HapQuals {
    pub hap1_ref_baseqs: Vec<u8>,
    pub hap1_alt_baseqs: Vec<u8>,
    pub hap2_ref_baseqs: Vec<u8>,
    pub hap2_alt_baseqs: Vec<u8>,
}

#[derive(Debug, Clone, Default)]
pub struct AlleleClass {
    pub allcls: u8,
    // 0: ref, 1: hete, 2: somatic
    pub prob: f64,
    // probability of allele class
}

#[derive(Debug, Clone, Default)]
pub struct CandidateSNP {
    pub chromosome: Vec<u8>,
    pub pos: i64,
    // position on the reference, 0-based
    pub alleles: [char; 2],
    // major and minor alleles
    pub allele_freqs: [f32; 2],
    // major and minor allele frequencies
    pub reference: char,
    pub alt_allele_fraction: f32,
    pub depth: u32,
    pub variant_type: i32,
    // 0: homo ref, 1: heterozygous SNP, 2: homozygous SNP, 3: triallelic SNP
    pub variant_quality: f64,
    // the confidence that the variant exists at this site given the data, phred-scaled
    pub genotype_probability: [f64; 3],
    // 0th: homo var, 1st: hete var, 2nd: homo ref
    pub genotype_quality: f64,
    pub genotype: i32,
    // eta: -1 (homVar), 0 (hetVar), 1 (homRef)
    pub haplotype: i32,
    // delta: 1,-1: hap1, hap2 if phased, 0: unassigned
    pub phase_score: f64,
    pub snp_cover_fragments: Vec<usize>,
    // index of the fragment cover this SNP
    pub rna_editing: bool,
    // A->G, T-C, rna_editing variant does not have haplotype information, no phasing
    pub germline: bool,
    // germline snp or not after evaluation of phased high frac het snps
    pub dense: bool,
    // lie in the dense snp region
    pub low_frac_het: bool,
    // low fraction het_var
    pub high_frac_het: bool,
    // high fraction het_var
    pub for_phasing: bool,
    // whether used in phasing
    pub hom_var: bool,
    // homo var
    pub single: bool,
    // current snp has surrounding haplotype links or not, only works for heterozygous snps
    pub non_selected: bool,
    // current site is not selected by phasing, it maybe HomRef, Somatic, Rna_edit or HomVar.
    pub cand_somatic: bool,
    // candidate site for somatic mutation detection
    pub somatic: bool,
    // detected somatic mutation by model
    pub somatic_score: f64,
    // phred somatic score
    pub hap_quals: HapQuals,
    // base qualities for identifying somatic mutation
    pub phase_set: u32,
    // phase set id is the position of the first snp in the phase set
    pub haplotype_expression: [u32; 4],
    // hap1_ref, hap1_alt, hap2_ref, hap2_alt
}

#[derive(Debug, Clone, Default)]
pub struct LD_Pair {
    // pub snp_idxes: [usize; 2],
    // // index of candidate SNPs (SNPFrag.snps), start node and end node
    pub ld_pairs: HashMap<[u8; 2], u32>,
    // support number of pair of alleles at two snp sites
}

impl LD_Pair {
    pub fn calculate_LD_R2(&self, A: u8, a: u8, B: u8, b: u8) -> f32 {
        // calculate r2 for two snps
        // A1, B1: alleles of snp1, A2, B2: alleles of snp2, https://en.wikipedia.org/wiki/Linkage_disequilibrium
        let (mut AB, mut Ab, mut aB, mut ab, mut sum) = (0.0, 0.0, 0.0, 0.0, 0.0);
        let (mut pA, mut pa, mut pB, mut pb) = (0.0, 0.0, 0.0, 0.0);
        if self.ld_pairs.contains_key(&[A, B]) {
            AB = self.ld_pairs[&[A, B]] as f32;
            sum += AB;
        }
        if self.ld_pairs.contains_key(&[A, b]) {
            Ab = self.ld_pairs[&[A, b]] as f32;
            sum += Ab;
        }
        if self.ld_pairs.contains_key(&[a, B]) {
            aB = self.ld_pairs[&[a, B]] as f32;
            sum += aB;
        }
        if self.ld_pairs.contains_key(&[a, b]) {
            ab = self.ld_pairs[&[a, b]] as f32;
            sum += ab;
        }

        // set minimum covered reads
        // if ((x11 < 2.0 || x22 < 2.0) && (x21 < 2.0 || x12 < 2.0)) || (sum < 4.0) { return 0.0; }

        if sum == 0.0 { return 0.0; }
        AB = AB / sum;    // P(AB)
        Ab = Ab / sum;    // P(Ab)
        aB = aB / sum;    // P(aB)
        ab = ab / sum;    // P(ab)

        pA = AB + Ab;    // P(A)
        pa = aB + ab;    // P(a)
        pB = AB + aB;    // P(B)
        pb = Ab + ab;    // P(b)

        let mut d = AB - pA * pB;  // D= P(AB) - P(A)P(B)
        let p = pA * (1.0 - pA) * pB * (1.0 - pB);  // P(A)*(1-P(A))*P(B)*(1-P(B))
        if pA == 1.0 && pB == 1.0 { return 1.0; }
        if p == 0.0 {
            // println!("AB: {}, Ab: {}, aB: {}, ab: {}", AB, Ab, aB, ab);
            return 0.0;
        }
        let r2 = d * d / p;   // r2 == 0, no correlation, r2 == 1, perfect positive correlation, r2 == -1. perfect negative correlation
        return r2;
    }
}

#[derive(Debug, Clone, Default)]
pub struct FragElem {
    pub snp_idx: usize,
    // index of candidate SNPs(SNPFrag.snps)
    pub pos: i64,
    // position on the reference, 0-based
    pub base: char,
    // base pair
    pub baseq: u8,
    // base quality
    pub strand: u32,
    // read strand,  0: forward, 1: reverse
    pub p: i32,
    // base allele on alphabet  {-1, 1, 0}, 1: base==ref, -1: base==alt, 0: not covered (bases except ref allele and alt allele, deletions or N)
    pub prob: f64,
    // error rate of observe current base
    pub phase_site: bool,
    // whether the snp used in phasing
}

#[derive(Debug, Clone, Default)]
pub struct Fragment {
    pub fragment_idx: usize,
    // index of multiple fragments(SNPFrag.fragments)
    pub read_id: String,
    // read name
    pub list: Vec<FragElem>,
    // single fragment
    pub haplotag: i32,
    // sigma: 0,1,-1
    pub assignment: i32,
    // haplotype assignment of the fragment, 0,1,2. 0: unassigned, 1: hap1, 2: hap2
    pub assignment_score: f64,
    // probability of the haplotype assignment
    pub exons: Vec<Exon>,
    // exons of the read on the reference, 0-based, [start, end)
    pub num_hete_links: u32,
    // number of linked heterozygous snps in the fragment
}

