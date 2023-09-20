// use noodles::{bam, core::Region, sam, sam::alignment::record};
use crate::bam_reader::Region;
use rust_htslib::{bam, bam::Read};

#[derive(Debug, Clone)]
pub struct ParsedRead {
    pub bam_record: bam::Record,
    // read from bam file
    pub pased_seq: Vec<u8>,
    // parsed sequence, expand insertion
    pub ref_pos: u32, // start position on reference based on the expanded insertion coordinate
}

impl ParsedRead {}

#[derive(Default, Debug, Clone)]
pub struct BaseFreq {
    pub a: u32,
    pub c: u32,
    pub g: u32,
    pub t: u32,
    pub n: u32,
    // number of introns
    pub d: u32,
    // number of deletions
    pub i: bool, // whether this position falls in an insertion
}

impl BaseFreq {
    pub fn get_score1(&self, query: u8) {}
    pub fn get_score2(&self, query: u8) {}
    pub fn get_score(&self, query: u8) {}
}

#[derive(Default, Debug, Clone)]
pub struct Profile {
    pub freq_vec: Vec<BaseFreq>,
}

impl Profile {
    pub fn init_with_pileup(&mut self, bam_path: String, ref_path: String, region: Region) {
        /*
        Generate the initial profile from the pileup of input bam file
        bam_path: bam file path
        ref_path: reference file path
        region: region to extract, 1-based, left closed, right open
         */
        let mut bam: bam::IndexedReader = bam::IndexedReader::from_path(bam_path).unwrap();
        bam.fetch((region.chr.as_str(), region.start, region.end)).unwrap(); // set region
        for p in bam.pileup() {
            let pileup = p.unwrap();
            let pos = pileup.pos(); // 0-based
            if pos + 1 < region.start || pos + 1 >= region.end {
                continue;
            }
            let mut bf = BaseFreq::default();
            let mut insert_bf: Vec<BaseFreq> = Vec::new();
            for alignment in pileup.alignments() {
                if alignment.is_refskip() {
                    bf.n += 1;
                } else if alignment.is_del() {
                    bf.d += 1;
                } else {
                    let q_pos = alignment.qpos().unwrap();
                    let record = alignment.record();
                    let seq = record.seq();
                    let base = seq[q_pos] as char;
                    match base {
                        'A' => bf.a += 1,
                        'a' => bf.a += 1,
                        'C' => bf.c += 1,
                        'c' => bf.c += 1,
                        'G' => bf.g += 1,
                        'g' => bf.g += 1,
                        'T' => bf.t += 1,
                        't' => bf.t += 1,
                        _ => {
                            panic!("Invalid nucleotide base: {}", base);
                        }
                    }

                    match alignment.indel() {
                        bam::pileup::Indel::Ins(len) => {
                            if len > insert_bf.len() as u32 {
                                for _ in 0..(len - insert_bf.len() as u32) {
                                    insert_bf.push(BaseFreq { a: 0, c: 0, g: 0, t: 0, n: 0, d: 0, i: true });   // fall in insertion
                                }
                            }
                            for tmpi in 1..=len {
                                // insert_segment.push(seq[q_pos.unwrap() + tmpi as usize] as char);
                                let base = seq[q_pos + tmpi as usize] as char;
                                let bf = &mut insert_bf[tmpi as usize - 1];
                                match base {
                                    'A' => bf.a += 1,
                                    'a' => bf.a += 1,
                                    'C' => bf.c += 1,
                                    'c' => bf.c += 1,
                                    'G' => bf.g += 1,
                                    'g' => bf.g += 1,
                                    'T' => bf.t += 1,
                                    't' => bf.t += 1,
                                    _ => {
                                        panic!("Invalid nucleotide base: {}", base);
                                    }
                                }
                            }
                        }
                        bam::pileup::Indel::Del(_len) => {}
                        bam::pileup::Indel::None => {}
                    }
                }
            }
            self.freq_vec.push(bf);
            if insert_bf.len() > 0 {
                for tmpi in 0..insert_bf.len() {
                    self.freq_vec.push(insert_bf[tmpi].clone());
                }
            }
        }
        // TODO: add reference sequence
    }
    pub fn subtract(&mut self, start_pos: u32, parsed_seq: &Vec<u8>) {}
    pub fn add(&mut self, start_pos: u32, parsed_seq: &Vec<u8>) {}
}
