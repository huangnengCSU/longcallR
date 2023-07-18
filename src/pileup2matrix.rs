use std::collections::HashMap;
use crate::matrix::PileupMatrix;
pub use crate::bam_reader::Region;
use rust_htslib::{bam, bam::Read};
use rust_htslib::faidx::Reader;
use bio::io::fasta;
use std::io;

pub fn generate_pileup_matrix(bam_file: &String, ref_path: &String, region_string: &String, matrices_vec: &mut Vec<PileupMatrix>) {
    let mut ref_reader = fasta::Reader::from_file(ref_path).unwrap();
    let mut ref_seqs = HashMap::new();
    for r in ref_reader.records() {
        let ref_record = r.unwrap();
        ref_seqs.insert(ref_record.id().to_string(), std::str::from_utf8(ref_record.seq()).unwrap().to_string());
    }

    let mut bam_reader: bam::IndexedReader = bam::IndexedReader::from_path(bam_file).unwrap();
    let region = Region::new(region_string.clone());
    if region.start == 0 && region.end == 0 {
        bam_reader.fetch(region.chr.as_str()).unwrap(); // set region
    } else {
        bam_reader.fetch((region.chr.as_str(), region.start, region.end)).unwrap(); // set region
    }

    let mut pre_pos = -1;
    let mut pre_tid = -1;
    let mut seq_lengths: HashMap<String, u32> = HashMap::new();
    let mut max_insertion_size = 0;
    let mut pm = PileupMatrix::new();
    // store the prevois accessed qpos of each read. key: qname, value: qpos.
    // For ref_skip and del, qpos is none. So if ref_skip or del is followed by insertion, we cannot get the qpos of inserted segment.
    let mut qpos_record_map = HashMap::new();
    for p in bam_reader.pileup() {
        let pileup = p.unwrap();
        let pos = pileup.pos(); // 0-based
        let tid = pileup.tid();
        if pos < region.start - 1 || pos > region.end - 1 { // region is 1-based and pos is 0-based
            continue;
        }
        let depth = pileup.depth();

        if pos as i32 != pre_pos + 1 || tid as i32 != pre_tid {
            if pm.base_matrix.len() > 1 {
                pm.padding();
                matrices_vec.push(pm.clone());
                pm.clear();
                qpos_record_map.clear();
            }
        }

        pre_pos = pos as i32;
        pre_tid = tid as i32;

        if pm.region.start == 0 && pm.region.end == 0 {
            pm.region.chr = region.chr.clone();
            pm.region.start = pos;
        }
        if pos > pm.region.end {
            pm.region.end = pos;
        }
        max_insertion_size = 0;
        seq_lengths.clear();
        let ref_allele_seq = ref_seqs.get(&pm.region.chr).unwrap()[pos as usize..pos as usize + 1].as_bytes();
        pm.insert(&"ref".to_string(), ref_allele_seq, &pos);
        seq_lengths.insert("ref".to_string(), ref_allele_seq.len() as u32);
        max_insertion_size = if ref_allele_seq.len() as i32 > max_insertion_size {
            ref_allele_seq.len() as i32
        } else {
            max_insertion_size
        };

        for alignment in pileup.alignments() {
            let q_pos = alignment.qpos();
            let record = alignment.record();
            let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
            let seq = record.seq();
            let mut allele_base = Vec::new();
            if !alignment.is_del() {
                allele_base.push(seq[q_pos.unwrap()]);
                qpos_record_map.insert(qname.clone(), q_pos.unwrap());
            } else if alignment.is_refskip() {
                allele_base.push(b'N');
            } else if alignment.is_del() {
                allele_base.push(b'-');
            }
            match alignment.indel() {
                bam::pileup::Indel::Ins(len) => {
                    if alignment.is_del() || alignment.is_refskip() {
                        if !qpos_record_map.get(&qname).is_none() {
                            // If the interval starts from the middle of intron, because there is no non-del or non-ref_skip position,
                            // the read may have no record in qpos_record_map, and it may return None.
                            // For this case, we just ignore the insertion segment of this read (rare).
                            let tmp_qpos = qpos_record_map.get(&qname).unwrap().clone();
                            for tmpi in 1..=len {
                                allele_base.push(seq[tmp_qpos + tmpi as usize]);
                                qpos_record_map.insert(qname.clone(), tmp_qpos + tmpi as usize);
                            }
                        }
                    } else {
                        for tmpi in 1..=len {
                            allele_base.push(seq[q_pos.unwrap() + tmpi as usize]);
                            qpos_record_map.insert(qname.clone(), q_pos.unwrap() + tmpi as usize);
                        }
                    }
                }
                bam::pileup::Indel::Del(len) => {}
                bam::pileup::Indel::None => {}
            }

            pm.insert(&qname.clone(), allele_base.as_slice(), &pos);
            max_insertion_size = if allele_base.len() as i32 > max_insertion_size {
                allele_base.len() as i32
            } else {
                max_insertion_size
            };
            seq_lengths.insert(qname.clone(), allele_base.len() as u32);
        }
        if max_insertion_size > 0 {
            pm.expand(&seq_lengths, &pos, max_insertion_size);
        }
    }
    if pm.base_matrix.len() > 1 {
        pm.padding();
        matrices_vec.push(pm.clone());
        pm.clear();
    }
}
