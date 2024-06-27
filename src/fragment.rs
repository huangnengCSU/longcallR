use bio::bio_types::strand::ReqStrand::Forward;
use rust_htslib::{bam, bam::Read, bam::record::Record};

use crate::exon::Exon;
use crate::snp::{FragElem, Fragment, LD_Pair};
use crate::snpfrags::{Edge, SNPFrag};
use crate::util::Region;

impl SNPFrag {
    pub fn get_fragments(&mut self, bam_path: &str, region: &Region, ref_seq: &Vec<u8>) {
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
            let mut truncated_regions: Vec<(usize, usize)> = Vec::new();    // left-closed right-closed
            let truncated_threshold = 200.0;
            let tr_match = 1.0;
            let (mut tr_start, mut tr_end) = (0, 0);
            let mut tr_score = 0.0;
            let mut tr_flag = false;
            for cg in cigar.iter() {
                match cg.char() as u8 {
                    b'S' | b'H' => {
                        continue;
                    }
                    b'M' | b'X' | b'=' => {
                        for _ in 0..cg.len() {
                            // mask the dense mismatch region in a read
                            let ref_base = ref_seq[pos_on_ref as usize] as char;
                            let base = seq[pos_on_query as usize] as char;
                            let mut baseq = record.qual()[pos_on_query as usize];
                            baseq = if baseq < 30 { baseq } else { 30 };
                            if base == ref_base {
                                tr_score -= tr_match;
                                if tr_score < 0.0 {
                                    if tr_flag {
                                        truncated_regions.push((tr_start, tr_end));
                                        tr_flag = false;
                                    }
                                    tr_score = 0.0;
                                    tr_start = pos_on_ref as usize;
                                    tr_end = pos_on_ref as usize;
                                }
                            } else {
                                tr_score += baseq as f64;
                                if tr_score > truncated_threshold {
                                    tr_flag = true;
                                    tr_end = pos_on_ref as usize;
                                }
                            }
                            if pos_on_ref == snp_pos {
                                let mut frag_elem = FragElem::default();
                                frag_elem.snp_idx = idx;
                                frag_elem.pos = pos_on_ref;
                                frag_elem.base = seq[pos_on_query as usize] as char;
                                let bq = record.qual()[pos_on_query as usize];
                                frag_elem.baseq = if bq < 30 { bq } else { 30 };
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
            if tr_flag && tr_end - tr_start > 0 {
                truncated_regions.push((tr_start, tr_end));
            }
            // remove frag_elem which is in truncated regions from fragment.list
            for frag_elem in fragment.list.clone().iter() {
                for tr in truncated_regions.iter() {
                    if frag_elem.pos >= tr.0 as i64 && frag_elem.pos <= tr.1 as i64 {
                        fragment.list.remove(fragment.list.iter().position(|x| x.pos == frag_elem.pos).unwrap());
                        break;
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

            // accomplish the pair wise LD_pair
            for i in 0..fragment.list.len() {
                for j in i + 1..fragment.list.len() {
                    let mut pair_start_idx: usize;
                    let mut pair_end_idx: usize;
                    let mut pair_start_base: u8;
                    let mut pair_end_base: u8;
                    if fragment.list[i].snp_idx < fragment.list[j].snp_idx {
                        pair_start_idx = fragment.list[i].snp_idx;
                        pair_end_idx = fragment.list[j].snp_idx;
                        pair_start_base = fragment.list[i].base as u8;
                        pair_end_base = fragment.list[j].base as u8;
                    } else {
                        pair_start_idx = fragment.list[j].snp_idx;
                        pair_end_idx = fragment.list[i].snp_idx;
                        pair_start_base = fragment.list[j].base as u8;
                        pair_end_base = fragment.list[i].base as u8;
                    }
                    let idx_key = &[pair_start_idx, pair_end_idx];  // snp index of start node, snp index of end node
                    if self.allele_pairs.contains_key(idx_key) {
                        let base_key = &[pair_start_base, pair_end_base];   // allele of start node, allele of end node
                        self.allele_pairs.get_mut(idx_key).unwrap().ld_pairs.entry(base_key.clone()).and_modify(|e| *e += 1).or_insert(1);
                    } else {
                        let mut ld_pair = LD_Pair::default();
                        ld_pair.ld_pairs.insert([pair_start_base, pair_end_base], 1);
                        self.allele_pairs.insert(idx_key.clone(), ld_pair);
                    }
                }
            }
            // hete snps >= 1
            let mut hete_links = 0;
            for fe in fragment.list.iter() {
                assert_ne!(fe.p, 0, "Error: fe.p==0");
                if fe.phase_site {
                    hete_links += 1;
                }
            }
            fragment.num_hete_links = hete_links;

            // For hifi data, min_linkers is 1, for nanopore data, min_linkers is 2 (preset). For phasing, at least min_linkers hete snps or at least 2 ase snps.
            assert!(self.min_linkers > 0, "Error: min_linkers <= 0");
            if hete_links >= self.min_linkers {
                // record edge count
                let mut phased_sites = Vec::new();
                for fe in fragment.list.iter() {
                    if self.candidate_snps[fe.snp_idx].for_phasing {
                        phased_sites.push(fe.clone());
                    }
                }
                for preidx in 0..phased_sites.len() - 1 {
                    let pre_elem = &phased_sites[preidx];
                    let next_elem = &phased_sites[preidx + 1];
                    if self.candidate_snps[pre_elem.snp_idx].for_phasing == false || self.candidate_snps[next_elem.snp_idx].for_phasing == false {
                        continue;
                    }
                    if self.edges.contains_key(&[pre_elem.snp_idx, next_elem.snp_idx]) {
                        let edge = self.edges.get_mut(&[pre_elem.snp_idx, next_elem.snp_idx]).unwrap();
                        edge.frag_idxes.push(fragment.fragment_idx);
                        edge.w += 1;
                    } else {
                        let mut edge = Edge::default();
                        edge.snp_idxes = [pre_elem.snp_idx, next_elem.snp_idx];
                        edge.snp_poses = [pre_elem.pos, next_elem.pos];
                        edge.frag_idxes.push(fragment.fragment_idx);
                        edge.w += 1;
                        self.edges.insert([pre_elem.snp_idx, next_elem.snp_idx], edge);
                    }
                }
                for fe in fragment.list.iter() {
                    // record each snp cover by which fragments
                    self.candidate_snps[fe.snp_idx].snp_cover_fragments.push(fragment.fragment_idx);
                }
                self.fragments.push(fragment);
            }
        }
    }

    pub fn clean_fragments(&mut self) {
        // trim fragment with edge support < 2
        for (k, e) in self.edges.iter() {
            if e.w < 2 {
                // if e.w == 0 { continue; }
                for idx in 0..e.frag_idxes.len() {
                    let frag_idx = e.frag_idxes[idx];
                    let mut frag = &mut self.fragments[frag_idx];
                    let pre_snp_idx = e.snp_idxes[0];
                    let next_snp_idx = e.snp_idxes[1];
                    self.candidate_snps[pre_snp_idx].snp_cover_fragments.retain(|&x| x != frag_idx);
                    self.candidate_snps[next_snp_idx].snp_cover_fragments.retain(|&x| x != frag_idx);
                    frag.list.retain(|x| x.snp_idx != pre_snp_idx && x.snp_idx != next_snp_idx);
                }
            }
        }
    }
}