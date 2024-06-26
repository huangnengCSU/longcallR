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
    pub depth: u32,
    pub variant_type: i32,
    // 0: homo ref, 1: heterozygous SNP, 2: homozygous SNP, 3: triallelic SNP
    pub variant_quality: f64,
    // the confidence that the variant exists at this site given the data, phred-scaled
    pub genotype_probability: [f64; 3],
    // 0th: homo var, 1st: hete var, 2nd: homo ref
    pub genotype_quality: f64,
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
pub struct Edge {
    pub snp_idxes: [usize; 2],
    // index of candidate SNPs(SNPFrag.snps), start node and end node
    pub snp_poses: [i64; 2],
    // position of candidate SNPs(SNPFrag.snps), start node and end node
    pub frag_idxes: Vec<usize>,
    // index of fragments(SNPFrag.fragments) cover this edge.
    pub w: f64,
    // weight of edge,  w_{ij}=\sum_{k}x_{ki}x_{kj}log\frac{1-\epsilon_{kij}}{\epsilon_{kij}}
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

