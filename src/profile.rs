use std::collections::HashMap;
use crate::bam_reader::Region;
use rust_htslib::{bam, bam::Read};
use seq_io::fasta::{Reader, Record};

#[derive(Debug, Clone)]
pub struct ParsedRead {
    pub bam_record: bam::record::Record,
    // read from bam file
    pub pased_seq: Vec<u8>,
    // parsed sequence, expand insertion
    pub ref_pos: i64, // start position on reference based on the expanded insertion coordinate
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
    pub i: bool,
    // whether this position falls in an insertion
    pub ref_base: char,   // reference base, A,C,G,T,N,-
}

impl BaseFreq {
    pub fn get_score1(&self, query: u8) {}
    pub fn get_score2(&self, query: u8) {}
    pub fn get_score(&self, query: u8) {}
}

#[derive(Default, Debug, Clone)]
pub struct Profile {
    pub freq_vec: Vec<BaseFreq>,
    pub region: Region, // 1-based, left closed, right open
}

impl Profile {
    pub fn init_with_pileup(&mut self, bam_path: &str, region: &Region) {
        /*
        Generate the initial profile from the pileup of input bam file
        bam_path: bam file path
        ref_path: reference file path
        region: region to extract, 1-based, left closed, right open
         */
        self.region = region.clone();
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
                                    insert_bf.push(BaseFreq { a: 0, c: 0, g: 0, t: 0, n: 0, d: 0, i: true, ref_base: '\x00' });   // fall in insertion
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
    }
    pub fn append_reference(&mut self, references: &HashMap<String, Vec<u8>>) {
        /*
        Fill the ``ref_base`` field in each BaseFreq.
        Optional. If not called, the ``ref_base`` field in each BaseFreq will be '\0'
         */
        let chr = &self.region.chr;
        let s = self.region.start - 1;    // 0-based, inclusive
        let mut p = s as usize;
        for i in 0..self.freq_vec.len() {
            if self.freq_vec[i].i {
                self.freq_vec[i].ref_base = '-'; // insertion
            } else {
                self.freq_vec[i].ref_base = references.get(chr).unwrap()[p] as char;
                p += 1;
            }
        }
    }
    pub fn subtract(&mut self, start_pos: u32, parsed_seq: &Vec<u8>) {}
    pub fn add(&mut self, start_pos: u32, parsed_seq: &Vec<u8>) {}
}

pub fn read_references(ref_path: &str) -> HashMap<String, Vec<u8>> {
    let mut references: HashMap<String, Vec<u8>> = HashMap::new();
    let mut reader = Reader::from_path(ref_path).unwrap();
    while let Some(record) = reader.next() {
        let record = record.expect("Error reading record");
        references.insert(record.id().unwrap().to_string(), record.full_seq().to_vec());
    }
    return references;
}

pub fn read_bam(bam_path: &str, region: &Region) -> HashMap<String, ParsedRead> {
    /*
    Read each record in the bam file and store the record in PasedRead.bam_record.
    */
    let mut parsed_reads: HashMap<String, ParsedRead> = HashMap::new();
    let mut bam: bam::IndexedReader = bam::IndexedReader::from_path(bam_path).unwrap();
    bam.fetch((region.chr.as_str(), region.start, region.end)).unwrap();
    for r in bam.records() {
        let record = r.unwrap();
        let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
        let parsed_read = ParsedRead { bam_record: record.clone(), pased_seq: Vec::new(), ref_pos: -1 };
        parsed_reads.insert(qname, parsed_read);
    }
    return parsed_reads;
}
