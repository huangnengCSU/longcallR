use bio::bio_types::strand::ReqStrand::Forward;
use rust_htslib::{bam, bam::Read, bam::record::Record};

use crate::exon::Exon;
use crate::phase::{FragElem, Fragment};
use crate::snpfrags::SNPFrag;
use crate::util::Region;

impl SNPFrag {
    pub fn get_fragments(&mut self, bam_path: &str, region: &Region) {
        let mut bam_reader: bam::IndexedReader = bam::IndexedReader::from_path(bam_path).unwrap();
        bam_reader.fetch((region.chr.as_str(), region.start, region.end)).unwrap();
        let mut record = Record::new();
        if self.candidate_snps.len() == 0 {
            return;
        }
        // assert!(self.min_linkers >= 0, "Error: min_linkers <= 0");
        while let Some(result) = bam_reader.read(&mut record) {
            if result.is_err() {
                panic!("BAM parsing failed...");
            }
            // TODO: filtering unmapped, secondary, supplementary reads?
            if record.is_unmapped() || record.is_secondary() || record.is_supplementary() {
                continue;
            }
            let pos = record.pos(); // 0-based
            if pos > self.candidate_snps.last().unwrap().pos {
                continue;
            }
            let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
            let cigar = record.cigar();
            let seq = record.seq().as_bytes();
            let strand = if record.strand() == Forward { 0 } else { 1 };
            let mut pos_on_ref = pos; // 0-based
            let mut pos_on_query = cigar.leading_softclips(); // 0-based
            let mut idx = 0; // index in self.ase_hete_snps
            let mut snp_pos = -1; // pre-computed position of candidate SNPs
            let mut alleles; // pre-computed alleles of candidate SNPs
            if pos <= self.candidate_snps.first().unwrap().pos {
                snp_pos = self.candidate_snps[idx].pos;
                alleles = self.candidate_snps[idx].alleles.clone();
            } else {
                // find the first SNP in the read
                while idx < self.candidate_snps.len() {
                    if self.candidate_snps[idx].pos >= pos {
                        break;
                    }
                    idx += 1;
                }
                assert!(idx < self.candidate_snps.len(), "Error: idx < self.candidate_snps.len()");
                snp_pos = self.candidate_snps[idx].pos;
                alleles = self.candidate_snps[idx].alleles.clone();
            }

            let mut fragment = Fragment::default();
            fragment.read_id = qname.clone();
            fragment.fragment_idx = self.fragments.len();

            let mut exon_start = -1;
            let mut exon_end = -1;
            exon_start = pos_on_ref;
            exon_end = exon_start;

            for cg in cigar.iter() {
                match cg.char() as u8 {
                    b'S' | b'H' => {
                        continue;
                    }
                    b'M' | b'X' | b'=' => {
                        for _ in 0..cg.len() {
                            if pos_on_ref == snp_pos {
                                let mut frag_elem = FragElem::default();
                                frag_elem.snp_idx = idx;
                                frag_elem.pos = pos_on_ref;
                                frag_elem.base = seq[pos_on_query as usize] as char;
                                frag_elem.baseq = record.qual()[pos_on_query as usize];
                                frag_elem.strand = strand;
                                frag_elem.prob = 10.0_f64.powf(-(frag_elem.baseq as f64) / 10.0);
                                if frag_elem.base == self.candidate_snps[idx].reference {
                                    frag_elem.p = 1; // reference allele
                                } else if (frag_elem.base == alleles[0] || frag_elem.base == alleles[1]) && frag_elem.base != self.candidate_snps[idx].reference {
                                    frag_elem.p = -1; // alternate allele
                                } else {
                                    frag_elem.p = 0; // not covered
                                }
                                if self.candidate_snps[frag_elem.snp_idx].for_phasing {
                                    frag_elem.phase_site = true;
                                }
                                // filtered SNP will not be used for haplotype phasing, ase snp will still be used for construct fragment.
                                if self.candidate_snps[frag_elem.snp_idx].dense == false && frag_elem.p != 0 {
                                    fragment.list.push(frag_elem);
                                }
                                idx += 1;
                                if idx < self.candidate_snps.len() {
                                    snp_pos = self.candidate_snps[idx].pos;
                                    alleles = self.candidate_snps[idx].alleles.clone();
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
                                if idx < self.candidate_snps.len() {
                                    snp_pos = self.candidate_snps[idx].pos;
                                    alleles = self.candidate_snps[idx].alleles.clone();
                                }
                            }
                            pos_on_ref += 1;
                        }
                    }
                    b'N' => {
                        exon_end = pos_on_ref;
                        if fragment.exons.len() == 0 {
                            fragment.exons.push(Exon {
                                chr: region.chr.clone(),
                                start: exon_start,
                                end: exon_end,
                                state: 0,
                            }); // start exon
                        } else {
                            fragment.exons.push(Exon {
                                chr: region.chr.clone(),
                                start: exon_start,
                                end: exon_end,
                                state: 1,
                            }); // internal exon
                        }
                        exon_start = -1;
                        exon_end = -1;
                        for _ in 0..cg.len() {
                            if pos_on_ref == snp_pos {
                                idx += 1;
                                if idx < self.candidate_snps.len() {
                                    snp_pos = self.candidate_snps[idx].pos;
                                    alleles = self.candidate_snps[idx].alleles.clone();
                                }
                            }
                            pos_on_ref += 1;
                        }
                        exon_start = pos_on_ref;
                        exon_end = exon_start;
                    }
                    _ => {
                        panic!("Error: unknown cigar operation: {}", cg.char());
                    }
                }
            }
            // the whole read is a single exon, no intron.
            if exon_start != -1 && exon_end != -1 && pos_on_ref > exon_end {
                exon_end = pos_on_ref;
            }
            if exon_end != exon_start {
                if fragment.exons.len() > 0 {
                    fragment.exons.push(Exon {
                        chr: region.chr.clone(),
                        start: exon_start,
                        end: exon_end,
                        state: 2,
                    }); // end exon
                } else {
                    fragment.exons.push(Exon {
                        chr: region.chr.clone(),
                        start: exon_start,
                        end: exon_end,
                        state: 3,
                    }); // single exon
                }
            }
            exon_start = -1;
            exon_end = -1;

            // hete snps >= 1
            let mut hete_links = 0;
            for fe in fragment.list.iter() {
                assert_ne!(fe.p, 0, "Error: fe.p==0");
                if fe.phase_site {
                    hete_links += 1;
                }
            }
            fragment.num_hete_links = hete_links;

            if fragment.read_id == "SRR18130587.4057475".to_string() {
                println!("hete_links: {}", hete_links);
                println!("fragment: {:?}", fragment);
            }

            // For hifi data, min_linkers is 1, for nanopore data, min_linkers is 2 (preset). For phasing, at least min_linkers hete snps or at least 2 ase snps.
            assert!(self.min_linkers > 0, "Error: min_linkers <= 0");
            if hete_links >= self.min_linkers {
                for fe in fragment.list.iter() {
                    // record each snp cover by which fragments
                    self.candidate_snps[fe.snp_idx].snp_cover_fragments.push(fragment.fragment_idx);
                }
                self.fragments.push(fragment);
            }
        }
    }
}