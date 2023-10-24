use crate::bam_reader::Region;
use rust_htslib::{bam, bam::Read, bam::record::Record};
use rust_htslib::htslib::drand48;

#[derive(Debug, Clone, Default)]
struct CandidateSNP {
    pos: i64,
    // position on the reference, 0-based
    ref_base: u8,
    // reference base
    alleles: [u8; 2],
    allele_freqs: [f32; 2],
}

#[derive(Debug, Clone, Default)]
struct Edge {
    snp_idx: [usize; 2],
    // index of candidate SNPs(SNPFrag.snps), start node and end node
    fragment_idx: usize,
    // index of multiple fragments(SNPFrag.fragments)
    p: [u8; 2],
    // haplotype of base of start node and end node from the ternary matrix on the alphabet {0, 1, -}(SNPFrag.reads_matrix)
    w: [f64; 2],
    // weight of edge, w[0] is the weight consisitent with haplotype, w[1] is the weight inconsistent with haplotype.
}

#[derive(Debug, Clone, Default)]
struct FragElem {
    snp_idx: usize,
    // index of candidate SNPs(SNPFrag.snps)
    pos: i64,
    // position on the reference, 0-based
    base: u8,
    // base pair
    baseq: u8,
    // base quality
    p: u8,
    // haplotype of base on the alphabet {0, 1, -}, 0: base==alleles[0], 1: base==alleles[1], -(2): not covered
}

#[derive(Debug, Clone, Default)]
struct Fragment {
    fragment_idx: usize,
    // index of multiple fragments(SNPFrag.fragments)
    read_id: String,
    // read name
    list: Vec<FragElem>,
    // single fragment
}

#[derive(Debug, Clone, Default)]
struct SNPFrag {
    region: Region,
    snps: Vec<CandidateSNP>,
    // candidate SNPs
    fragments: Vec<Fragment>,
    // multiple fragments
    haplotype: Vec<u8>,
    // hap1. hap2 is bitwised complement of hap1
    edges: Vec<Edge>,
    // edges of the graph
}

impl SNPFrag {
    fn get_candidate_snps(&mut self) {
        // get candidate SNPs
    }

    unsafe fn init_haplotypes(&mut self) {
        // initialize haplotypes
        for i in 0..self.snps.len() {
            if drand48() < 0.5 {
                self.haplotype.push(0);
            } else {
                self.haplotype.push(1);
            }
        }
    }

    fn get_fragments(&mut self, bam_path: &str, region: &Region) {
        let mut bam_reader: bam::IndexedReader = bam::IndexedReader::from_path(bam_path).unwrap();
        bam_reader.fetch((region.chr.as_str(), region.start, region.end)).unwrap();
        let mut record = Record::new();
        while let Some(result) = bam_reader.read(&mut record) {
            if result.is_err() {
                panic!("BAM parsing failed...");
            }
            if record.is_unmapped() || record.is_secondary() || record.is_supplementary() {
                continue;
            }
            let pos = record.pos(); // 0-based
            if pos > self.snps.last().unwrap().pos { continue; }
            let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
            let cigar = record.cigar();
            let seq = record.seq().as_bytes();
            let mut pos_on_ref = pos;
            let mut pos_on_query = cigar.leading_softclips();
            let mut snp_offset = 0; // index of candidate SNPs
            let mut snp_pos = -1;   // pre-computed position of candidate SNPs
            let mut alleles = [2, 2];   // pre-computed alleles of candidate SNPs
            if pos <= self.snps.first().unwrap().pos {
                snp_pos = self.snps[snp_offset].pos;
                alleles = self.snps[snp_offset].alleles.clone();
            } else {
                // find the first SNP in the read
                while snp_offset < self.snps.len() {
                    if self.snps[snp_offset].pos >= pos {
                        break;
                    }
                    snp_offset += 1;
                }
                assert!(snp_offset < self.snps.len(), "Error: snp_offset < self.snps.len()");
                snp_pos = self.snps[snp_offset].pos;
                alleles = self.snps[snp_offset].alleles.clone();
            }

            let mut fragment = Fragment::default();
            fragment.read_id = qname.clone();
            fragment.fragment_idx = self.fragments.len();
            for cg in cigar.iter() {
                if pos_on_ref > self.snps.last().unwrap().pos {
                    break;
                }
                match cg.char() as u8 {
                    b'S' | b'H' => {
                        continue;
                    }
                    b'M' => {
                        for _ in 0..cg.len() {
                            assert!(pos_on_ref <= snp_pos, "Error: pos_on_ref <= snp_pos");
                            if pos_on_ref == snp_pos {
                                let mut frag_elem = FragElem::default();
                                frag_elem.snp_idx = snp_offset;
                                frag_elem.pos = pos_on_ref;
                                frag_elem.base = seq[pos_on_query as usize];
                                frag_elem.baseq = record.qual()[pos_on_query as usize];
                                if frag_elem.base == alleles[0] {
                                    frag_elem.p = 0;
                                } else if frag_elem.base == alleles[1] {
                                    frag_elem.p = 1;
                                } else {
                                    frag_elem.p = 2;
                                }
                                if fragment.list.len() > 0 {
                                    let mut edge = Edge::default();
                                    edge.snp_idx = [fragment.list.last().unwrap().snp_idx, frag_elem.snp_idx];
                                    edge.fragment_idx = fragment.fragment_idx;
                                    edge.p = [fragment.list.last().unwrap().p, frag_elem.p];
                                    let q1 = 0.1_f64.powf((fragment.list.last().unwrap().baseq as f64 - 33.0) / 10.0);
                                    let q2 = 0.1_f64.powf((frag_elem.baseq as f64 - 33.0) / 10.0);
                                    let w1 = q1 * q2 + (1.0 - q1) * (1.0 - q2);
                                    let w2 = q1 * (1.0 - q2) + (1.0 - q1) * q2;
                                    edge.w = [(w1 / w2).log10(), (w2 / w1).log10()];
                                    self.edges.push(edge);
                                }
                                fragment.list.push(frag_elem);
                                snp_offset += 1;
                                snp_pos = self.snps[snp_offset].pos;
                                alleles = self.snps[snp_offset].alleles.clone();
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
                            assert!(pos_on_ref <= snp_pos, "Error: pos_on_ref <= snp_pos");
                            if pos_on_ref == snp_pos {
                                let mut frag_elem = FragElem::default();
                                frag_elem.snp_idx = snp_offset;
                                frag_elem.pos = pos_on_ref;
                                frag_elem.base = b'-';
                                frag_elem.baseq = 0;
                                frag_elem.p = 2;
                                if fragment.list.len() > 0 {
                                    let mut edge = Edge::default();
                                    edge.snp_idx = [fragment.list.last().unwrap().snp_idx, frag_elem.snp_idx];
                                    edge.fragment_idx = fragment.fragment_idx;
                                    edge.p = [fragment.list.last().unwrap().p, frag_elem.p];
                                    let q1 = 0.1_f64.powf((fragment.list.last().unwrap().baseq as f64 - 33.0) / 10.0);
                                    let q2 = 0.1_f64.powf((frag_elem.baseq as f64 - 33.0) / 10.0);
                                    let w1 = q1 * q2 + (1.0 - q1) * (1.0 - q2);
                                    let w2 = q1 * (1.0 - q2) + (1.0 - q1) * q2;
                                    edge.w = [(w1 / w2).log10(), (w2 / w1).log10()];
                                    self.edges.push(edge);
                                }
                                fragment.list.push(frag_elem);
                                snp_offset += 1;
                                snp_pos = self.snps[snp_offset].pos;
                                alleles = self.snps[snp_offset].alleles.clone();
                            }
                            pos_on_ref += 1;
                        }
                    }
                    b'N' => {
                        for _ in 0..cg.len() {
                            assert!(pos_on_ref <= snp_pos, "Error: pos_on_ref <= snp_pos");
                            if pos_on_ref == snp_pos {
                                let mut frag_elem = FragElem::default();
                                frag_elem.snp_idx = snp_offset;
                                frag_elem.pos = pos_on_ref;
                                frag_elem.base = b'-';
                                frag_elem.baseq = 0;
                                frag_elem.p = 2;
                                if fragment.list.len() > 0 {
                                    let mut edge = Edge::default();
                                    edge.snp_idx = [fragment.list.last().unwrap().snp_idx, frag_elem.snp_idx];
                                    edge.fragment_idx = fragment.fragment_idx;
                                    edge.p = [fragment.list.last().unwrap().p, frag_elem.p];
                                    let q1 = 0.1_f64.powf((fragment.list.last().unwrap().baseq as f64 - 33.0) / 10.0);
                                    let q2 = 0.1_f64.powf((frag_elem.baseq as f64 - 33.0) / 10.0);
                                    let w1 = q1 * q2 + (1.0 - q1) * (1.0 - q2);
                                    let w2 = q1 * (1.0 - q2) + (1.0 - q1) * q2;
                                    edge.w = [(w1 / w2).log10(), (w2 / w1).log10()];
                                    self.edges.push(edge);
                                }
                                fragment.list.push(frag_elem);
                                snp_offset += 1;
                                snp_pos = self.snps[snp_offset].pos;
                                alleles = self.snps[snp_offset].alleles.clone();
                            }
                            pos_on_ref += 1;
                        }
                    }
                    _ => {
                        panic!("Error: unknown cigar operation");
                    }
                }
            }
            if fragment.list.len() > 0 {
                self.fragments.push(fragment);
            }
        }
    }

    fn optimization_using_maxcut(&mut self) {
        // optimization using maxcut

    }
}

