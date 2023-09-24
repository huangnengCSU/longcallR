use std::collections::HashMap;
use crate::bam_reader::Region;
use rust_htslib::{bam, bam::Read};
use seq_io::fasta::{Reader, Record};

#[derive(Debug, Clone)]
pub struct ParsedRead {
    pub bam_record: bam::record::Record,
    // read from bam file
    pub parsed_seq: Vec<u8>,
    // parsed sequence, expand insertion
    pub pos_on_profile: i64, // start position of the first aligned base on the profile, 0-based
}

impl ParsedRead {
    pub fn init_parsed_seq(&mut self, profile: &Profile) {
        /*
        Parse cigar sequence and get the parsed sequence aligned to the profile (expanded insertion).
        */
        let ref_start = self.bam_record.pos();  // 0-based
        let mut i = 0;  // the index of the first base in the profile (with expanded insertion), offset on profile
        let mut j = profile.region.start - 1;   // 0-based, offset on reference
        if j as i64 > ref_start {
            panic!("read start position exceeds the profile region!");
        }

        while i < profile.freq_vec.len() {
            if !profile.freq_vec[i].i {
                if j as i64 == ref_start {
                    break;
                }
                j += 1;
            }
            i += 1;
        }

        println!("readname:{}, i = {}, j = {}, refbase: {}", std::str::from_utf8(self.bam_record.qname()).unwrap(), i, j, profile.freq_vec[i].ref_base);

        // TODO: parse cigar and get parsed sequence which is also aligned to the profile frequency vector
        self.pos_on_profile = i as i64;
        let mut pos_on_read = 0;
        let seq = self.bam_record.seq();
        for cg in self.bam_record.cigar().iter() {
            match cg.char() as u8 {
                b'S' => {
                    pos_on_read = cg.len();
                    continue;
                }
                b'H' => {
                    continue;
                }
                b'M' => {
                    let mut k = cg.len() as i32;
                    while k > 0 {
                        if !profile.freq_vec[i].i {
                            self.parsed_seq.push(seq[pos_on_read as usize]);
                            k -= 1;
                            pos_on_read += 1;
                        } else {
                            self.parsed_seq.push('-' as u8);
                        }
                        i += 1;
                    }
                }
                b'I' => {
                    let mut k = cg.len() as i32;
                    while k > 0 {
                        self.parsed_seq.push(seq[pos_on_read as usize]);
                        assert_eq!(profile.freq_vec[i].i, true);
                        k -= 1;
                        i += 1;
                        pos_on_read += 1;
                    }
                }
                b'D' => {
                    let mut k = cg.len() as i32;
                    while k > 0 {
                        if !profile.freq_vec[i].i {
                            self.parsed_seq.push('-' as u8);
                            k -= 1;
                        } else {
                            self.parsed_seq.push('-' as u8);
                        }
                        i += 1;
                    }
                }
                b'N' => {
                    let mut k = cg.len() as i32;
                    while k > 0 {
                        if !profile.freq_vec[i].i {
                            self.parsed_seq.push('N' as u8);
                            k -= 1;
                        } else {
                            self.parsed_seq.push('N' as u8);
                        }
                        i += 1;
                    }
                }
                _ => {
                    panic!("Invalid cigar character: {}", cg.char());
                }
            }
        }
    }
}

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
    pub region: Region, // position of aligned base (instead of input reigon), 1-based, left closed, right open
}

impl Profile {
    pub fn init_with_pileup(&mut self, bam_path: &str, region: &Region) {
        /*
        Generate the initial profile from the pileup of input bam file
        bam_path: bam file path
        ref_path: reference file path
        region: region to extract, 1-based, left closed, right open
         */
        // self.region = region.clone();
        let mut freq_vec_start = u32::MAX;  // 1-based
        let mut freq_vec_end = u32::MIN;    // 1-based
        let mut bam: bam::IndexedReader = bam::IndexedReader::from_path(bam_path).unwrap();
        bam.fetch((region.chr.as_str(), region.start, region.end)).unwrap(); // set region
        let mut read_positions: HashMap<String, u32> = HashMap::new();
        for p in bam.pileup() {
            let pileup = p.unwrap();
            let pos = pileup.pos(); // 0-based
            if pos + 1 < region.start || pos + 1 >= region.end {
                continue;
            }
            if pos + 1 < freq_vec_start {
                freq_vec_start = pos + 1;
            }
            if pos + 1 > freq_vec_end {
                freq_vec_end = pos + 1;
            }
            let mut bf = BaseFreq::default();
            let mut insert_bf: Vec<BaseFreq> = Vec::new();
            for alignment in pileup.alignments() {
                if alignment.is_refskip() {
                    bf.n += 1;
                    match alignment.indel() {
                        bam::pileup::Indel::Ins(len) => {
                            let record = alignment.record();
                            let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
                            let seq = record.seq();
                            if len > insert_bf.len() as u32 {
                                for _ in 0..(len - insert_bf.len() as u32) {
                                    insert_bf.push(BaseFreq { a: 0, c: 0, g: 0, t: 0, n: 0, d: 0, i: true, ref_base: '\x00' });   // fall in insertion
                                }
                            }
                            let q_pos = read_positions.get(&qname).unwrap();
                            for tmpi in 1..=len {
                                // insert_segment.push(seq[q_pos.unwrap() + tmpi as usize] as char);
                                let base = seq[(q_pos + tmpi) as usize] as char;
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
                            read_positions.insert(qname.clone(), *q_pos + len);
                        }
                        bam::pileup::Indel::Del(_len) => {}
                        bam::pileup::Indel::None => {}
                    }
                } else if alignment.is_del() {
                    bf.d += 1;
                } else {
                    let q_pos = alignment.qpos().unwrap();
                    let record = alignment.record();
                    let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
                    read_positions.insert(qname.clone(), q_pos as u32);
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
                                let ibf = &mut insert_bf[tmpi as usize - 1];
                                match base {
                                    'A' => ibf.a += 1,
                                    'a' => ibf.a += 1,
                                    'C' => ibf.c += 1,
                                    'c' => ibf.c += 1,
                                    'G' => ibf.g += 1,
                                    'g' => ibf.g += 1,
                                    'T' => ibf.t += 1,
                                    't' => ibf.t += 1,
                                    _ => {
                                        panic!("Invalid nucleotide base: {}", base);
                                    }
                                }
                            }
                            read_positions.insert(qname.clone(), q_pos as u32 + len);
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
        self.region = Region { chr: region.chr.clone(), start: freq_vec_start, end: freq_vec_end + 1 };
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
    Read each record in the bam file and store the record in ParsedRead.bam_record.
    */
    let mut parsed_reads: HashMap<String, ParsedRead> = HashMap::new();
    let mut bam: bam::IndexedReader = bam::IndexedReader::from_path(bam_path).unwrap();
    bam.fetch((region.chr.as_str(), region.start, region.end)).unwrap();
    for r in bam.records() {
        let record = r.unwrap();
        let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
        let parsed_read = ParsedRead { bam_record: record.clone(), parsed_seq: Vec::new(), pos_on_profile: -1 };
        parsed_reads.insert(qname, parsed_read);
    }
    return parsed_reads;
}
