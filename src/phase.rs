use crate::bam_reader::Region;
use rust_htslib::{bam, bam::Read, bam::record::Record};

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
    snp_idx: usize,
    // index of candidate SNPs(SNPFrag.snps)
    fragment_idx: usize,
    // index of multiple fragments(SNPFrag.fragments)
    p: [u8; 2],
    // haplotype of base of start node and end node from the ternary matrix on the alphabet {0, 1, -}(SNPFrag.reads_matrix)
    w: f32,
    // weight of edge
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

    fn init_haplotypes(&mut self) {
        // initialize haplotypes
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
}

