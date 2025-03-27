use crate::exon::Exon;
use crate::snp::{FragElem, Fragment, LD_Pair};
use crate::snpfrags::{Edge, SNPFrag};
use crate::util::Region;
use bio::bio_types::strand::ReqStrand::Forward;
use rust_htslib::bam::record::Aux;
use rust_htslib::{bam, bam::record::Record, bam::Read};

impl SNPFrag {
    pub fn get_fragments(
        &mut self,
        bam_path: &str,
        region: &Region,
        ref_seq: &Vec<u8>,
        min_mapq: u8,
        min_read_length: usize,
        divergence: f32,
    ) {
        let mut bam_reader: bam::IndexedReader = bam::IndexedReader::from_path(bam_path).unwrap();
        bam_reader
            .fetch((region.chr.as_str(), region.start, region.end))
            .unwrap();
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
            if record.mapq() < min_mapq
                || record.seq_len() < min_read_length
                || record.is_unmapped()
                || record.is_secondary()
                || record.is_supplementary()
            {
                continue;
            }

            match record.aux(b"de") {
                Ok(Aux::Float(f)) => {
                    if f >= divergence {
                        continue;
                    }
                }
                _ => {}
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
                assert!(
                    idx < self.candidate_snps.len(),
                    "Error: idx < self.candidate_snps.len()"
                );
                snp_pos = self.candidate_snps[idx].pos;
                alleles = self.candidate_snps[idx].alleles.clone();
            }

            let mut fragment = Fragment::default();
            fragment.read_id = qname.clone();
            fragment.fragment_idx = self.fragments.len();

            // let mut truncation_regions: Vec<(usize, usize)> = Vec::new(); // left-closed right-closed
            // let truncated_threshold = 200.0;
            // let match_socre = 1.0;
            // let (mut truncation_start, mut truncation_end) = (0, 0);
            // let mut truncation_score = 0.0;
            // let mut truncation_flag = false;

            for cg in cigar.iter() {
                match cg.char() as u8 {
                    b'S' | b'H' => {
                        continue;
                    }
                    b'M' | b'X' | b'=' => {
                        for _ in 0..cg.len() {
                            // mask the dense mismatch region in a read
                            let ref_base = ref_seq[pos_on_ref as usize] as char;
                            let read_base = seq[pos_on_query as usize] as char;
                            let read_baseq = record.qual()[pos_on_query as usize];
                            // if read_base == ref_base {
                            //     truncation_score -= match_socre;
                            //     if truncation_score < 0.0 {
                            //         if truncation_flag {
                            //             truncation_regions.push((truncation_start, truncation_end));
                            //             truncation_flag = false;
                            //         }
                            //         truncation_score = 0.0;
                            //         truncation_start = pos_on_ref as usize;
                            //         truncation_end = pos_on_ref as usize;
                            //     }
                            // } else {
                            //     truncation_score += read_baseq as f64;
                            //     if truncation_score > truncated_threshold {
                            //         truncation_flag = true;
                            //         truncation_end = pos_on_ref as usize;
                            //     }
                            // }
                            if pos_on_ref == snp_pos {
                                let mut frag_elem = FragElem::default();
                                frag_elem.snp_idx = idx;
                                frag_elem.pos = pos_on_ref;
                                frag_elem.base = seq[pos_on_query as usize] as char;
                                frag_elem.baseq = if record.qual()[pos_on_query as usize] < 30 {
                                    record.qual()[pos_on_query as usize]
                                } else {
                                    30
                                };
                                frag_elem.strand = strand;
                                frag_elem.prob = 10.0_f64.powf(-(frag_elem.baseq as f64) / 10.0);
                                if frag_elem.base == self.candidate_snps[idx].reference {
                                    frag_elem.p = 1; // reference allele
                                } else if (frag_elem.base == alleles[0]
                                    || frag_elem.base == alleles[1])
                                    && frag_elem.base != self.candidate_snps[idx].reference
                                {
                                    frag_elem.p = -1; // alternate allele
                                } else {
                                    frag_elem.p = 0; // not covered
                                }
                                if self.candidate_snps[frag_elem.snp_idx].for_phasing {
                                    frag_elem.phase_site = true;
                                }
                                // filtered SNP will not be used for haplotype phasing, ase snp will still be used for construct fragment.
                                if self.candidate_snps[frag_elem.snp_idx].dense == false
                                    && frag_elem.p != 0
                                {
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
                    _ => {
                        panic!("Error: unknown cigar operation: {}", cg.char());
                    }
                }
            }

            if fragment.read_id == "m84036_230523_222603_s1/132912932/ccs/1625_5860" {
                println!("\nInit fragment: {:?}", fragment);
            }

            // if truncation_flag && truncation_end - truncation_start > 0 {
            //     truncation_regions.push((truncation_start, truncation_end));
            // }

            if fragment.read_id == "m84036_230523_222603_s1/132912932/ccs/1625_5860" {
                println!("\nModified fragment: {:?}", fragment);
            }

            // remove frag_elem which is in truncated regions from fragment.list
            // fragment.list.retain(|frag_elem| {
            //     !truncation_regions
            //         .iter()
            //         .any(|tr| frag_elem.pos >= tr.0 as i64 && frag_elem.pos <= tr.1 as i64)
            // });

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
                    self.allele_pairs
                        .entry([pair_start_idx, pair_end_idx])
                        .and_modify(|ld_pair| {
                            ld_pair
                                .ld_pairs
                                .entry([pair_start_base, pair_end_base])
                                .and_modify(|count| *count += 1)
                                .or_insert(1);
                        })
                        .or_insert_with(|| {
                            let mut ld_pair = LD_Pair::default();
                            ld_pair.ld_pairs.insert([pair_start_base, pair_end_base], 1);
                            ld_pair
                        });
                }
            }

            let hete_links = fragment
                .list
                .iter()
                .inspect(|fe| assert_ne!(fe.p, 0, "Error: fe.p == 0"))
                .filter(|fe| fe.phase_site)
                .count() as u32;
            // let hete_links = fragment.list.iter()
            //     .inspect(|fe| assert_ne!(fe.p, 0, "Error: fe.p == 0"))
            //     .count() as u32;
            fragment.num_hete_links = hete_links;

            assert!(self.min_linkers > 0, "Error: min_linkers <= 0");
            if hete_links >= self.min_linkers {
                fragment.for_phasing = true;
                let phased_sites: Vec<_> = fragment
                    .list
                    .iter()
                    .filter(|fe| self.candidate_snps[fe.snp_idx].for_phasing)
                    .cloned()
                    .collect();
                // let phased_sites: Vec<_> = fragment.list.iter()
                //     .cloned()
                //     .collect();

                for preidx in 0..phased_sites.len() - 1 {
                    let (pre_elem, next_elem) =
                        if phased_sites[preidx].pos < phased_sites[preidx + 1].pos {
                            (&phased_sites[preidx], &phased_sites[preidx + 1])
                        } else {
                            (&phased_sites[preidx + 1], &phased_sites[preidx])
                        };
                    if self.candidate_snps[pre_elem.snp_idx].for_phasing == false
                        || self.candidate_snps[next_elem.snp_idx].for_phasing == false
                    {
                        continue;
                    }
                    self.edges
                        .entry([pre_elem.snp_idx, next_elem.snp_idx])
                        .and_modify(|edge| {
                            edge.frag_idxes.push(fragment.fragment_idx);
                            edge.w += 1;
                        })
                        .or_insert_with(|| {
                            let mut edge = Edge::default();
                            edge.snp_idxes = [pre_elem.snp_idx, next_elem.snp_idx];
                            edge.snp_poses = [pre_elem.pos, next_elem.pos];
                            edge.frag_idxes.push(fragment.fragment_idx);
                            edge.w += 1;
                            edge
                        });
                }
                for fe in fragment.list.iter() {
                    self.candidate_snps[fe.snp_idx]
                        .snp_cover_fragments
                        .push(fragment.fragment_idx);
                }
                self.fragments.push(fragment);
            } else {
                fragment.for_phasing = false;
                for fe in fragment.list.iter() {
                    self.candidate_snps[fe.snp_idx]
                        .snp_cover_fragments
                        .push(fragment.fragment_idx);
                }
                self.fragments.push(fragment);
            }
        }
    }

    pub fn clean_fragments(&mut self) {
        // trim fragment with edge support < 2
        for (_, e) in self.edges.iter() {
            if e.w < 2 {
                // if e.w == 0 { continue; }
                for idx in 0..e.frag_idxes.len() {
                    let frag_idx = e.frag_idxes[idx];
                    let mut frag = &mut self.fragments[frag_idx];
                    let pre_snp_idx = e.snp_idxes[0];
                    let next_snp_idx = e.snp_idxes[1];
                    self.candidate_snps[pre_snp_idx]
                        .snp_cover_fragments
                        .retain(|&x| x != frag_idx);
                    self.candidate_snps[next_snp_idx]
                        .snp_cover_fragments
                        .retain(|&x| x != frag_idx);
                    frag.list
                        .retain(|x| x.snp_idx != pre_snp_idx && x.snp_idx != next_snp_idx);
                }
            }
        }
    }
}
