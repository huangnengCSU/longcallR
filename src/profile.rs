use std::collections::HashMap;
use crate::bam_reader::Region;
use rust_htslib::{bam, bam::Read};
use seq_io::fasta::{Reader, Record};
use rust_lapper::{Interval, Lapper};

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
                            self.parsed_seq.push('*' as u8);
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
                            self.parsed_seq.push('*' as u8);
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
                            self.parsed_seq.push('*' as u8);
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
    pub fn get_score1(&self, query: u8) -> i32 {
        let max_cnt = self.a.max(self.c.max(self.g.max(self.t.max(self.n.max(self.d)))));
        match query {
            b'A' => if self.a == max_cnt { 0 } else { 1 },
            b'C' => if self.c == max_cnt { 0 } else { 1 },
            b'G' => if self.g == max_cnt { 0 } else { 1 },
            b'T' => if self.t == max_cnt { 0 } else { 1 },
            b'N' => if self.n == max_cnt { 0 } else { 1 },
            b'-' => if self.d == max_cnt { 0 } else { 1 },
            _ => {
                panic!("Invalid query base: {}", query as char);
            }
        }
    }
    pub fn get_score2(&self, query: u8) -> f64 {
        let s = (self.a + self.c + self.g + self.t + self.d + self.n) as f64;
        if s == 0.0 {
            return 0.0;
        }
        match query {
            b'A' => (s - self.a as f64) / s,
            b'C' => (s - self.c as f64) / s,
            b'G' => (s - self.g as f64) / s,
            b'T' => (s - self.t as f64) / s,
            b'N' => (s - self.n as f64) / s,
            b'-' => (s - self.d as f64) / s,
            _ => {
                panic!("Invalid query base: {}", query as char);
            }
        }
    }
    pub fn get_score(&self, query: u8) -> f64 {
        (self.get_score1(query) as f64 + self.get_score2(query)) / 2.0
    }

    pub fn get_depth_exclude_intron(&self) -> u32 {
        self.a + self.c + self.g + self.t + self.d
    }

    pub fn get_depth_include_intron(&self) -> u32 {
        self.a + self.c + self.g + self.t + self.d + self.n
    }

    pub fn get_intron_ratio(&self) -> (u32, f64) {
        (self.n, (self.n as f64) / (self.get_depth_include_intron() as f64))
    }
}

#[derive(Default, Debug, Clone)]
pub struct Profile {
    pub freq_vec: Vec<BaseFreq>,
    pub region: Region,
    // position of aligned base (instead of input reigon), 1-based, left closed, right open
    pub forward_donor_penalty: Vec<f64>,
    pub forward_acceptor_penalty: Vec<f64>,
    pub reverse_donor_penalty: Vec<f64>,
    pub reverse_acceptor_penalty: Vec<f64>,
    pub intron_intervals: Vec<Interval<usize, bool>>,    // intron regions, the base of each read in this region is 'N'.
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
        let mut read_positions: HashMap<String, u32> = HashMap::new();  // store offset on profile of each read
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
    pub fn subtract(&mut self, offset_on_profile: u32, parsed_seq: &Vec<u8>) {
        /*
        Subtract the parsed sequence from the profile.
        offset_on_profile: offset on profile, 0-based
        parsed_seq: parsed sequence, expand insertion
        */
        for i in 0..parsed_seq.len() {
            match parsed_seq[i] {
                b'A' => {
                    self.freq_vec[offset_on_profile as usize + i].a -= 1;
                }
                b'C' => {
                    self.freq_vec[offset_on_profile as usize + i].c -= 1;
                }
                b'G' => {
                    self.freq_vec[offset_on_profile as usize + i].g -= 1;
                }
                b'T' => {
                    self.freq_vec[offset_on_profile as usize + i].t -= 1;
                }
                b'N' => {
                    if !self.freq_vec[offset_on_profile as usize + i].i {
                        self.freq_vec[offset_on_profile as usize + i].n -= 1;
                    }
                }
                b'-' => {
                    if !self.freq_vec[offset_on_profile as usize + i].i {
                        self.freq_vec[offset_on_profile as usize + i].d -= 1;
                    }
                }
                b'*' => {
                    continue;
                }
                _ => {
                    panic!("Invalid base: {}", parsed_seq[i] as char);
                }
            }
        }
    }
    pub fn add(&mut self, offset_on_profile: u32, parsed_seq: &Vec<u8>) {
        /*
        Add the parsed sequence to the profile.
        offset_on_profile: offset on profile, 0-based
        parsed_seq: parsed sequence, expand insertion
        */
        for i in 0..parsed_seq.len() {
            match parsed_seq[i] {
                b'A' => {
                    self.freq_vec[offset_on_profile as usize + i].a += 1;
                }
                b'C' => {
                    self.freq_vec[offset_on_profile as usize + i].c += 1;
                }
                b'G' => {
                    self.freq_vec[offset_on_profile as usize + i].g += 1;
                }
                b'T' => {
                    self.freq_vec[offset_on_profile as usize + i].t += 1;
                }
                b'N' => {
                    if !self.freq_vec[offset_on_profile as usize + i].i {
                        self.freq_vec[offset_on_profile as usize + i].n += 1;
                    }
                }
                b'-' => {
                    if !self.freq_vec[offset_on_profile as usize + i].i {
                        self.freq_vec[offset_on_profile as usize + i].d += 1;
                    }
                }
                b'*' => {
                    continue;
                }
                _ => {
                    panic!("Invalid base: {}", parsed_seq[i] as char);
                }
            }
        }
    }

    pub fn cal_intron_penalty(&mut self) {
        let mut ref_sliding_window: Vec<u8> = Vec::new();
        for i in 0..self.freq_vec.len() + 1 {
            if i as i32 - 3 < 0 {
                self.forward_acceptor_penalty.push(30.0);
                self.reverse_acceptor_penalty.push(30.0);
            } else {
                let mut j = i as i32 - 1;
                ref_sliding_window.clear();
                while ref_sliding_window.len() < 3 && j >= 0 {
                    if !self.freq_vec[j as usize].i {
                        ref_sliding_window.push(self.freq_vec[j as usize].ref_base as u8);
                    }
                    j -= 1;
                }
                let mut tstr = String::new();
                for c in ref_sliding_window.iter().rev() {
                    tstr.push(*c as char);
                }
                if tstr.len() < 3 {
                    self.forward_acceptor_penalty.push(30.0);
                    self.reverse_acceptor_penalty.push(30.0);
                } else {
                    if tstr[1..3] == "AC".to_string() {
                        self.forward_acceptor_penalty.push(21.0);
                    } else if tstr == "AAG".to_string() || tstr == "GAG".to_string() {
                        self.forward_acceptor_penalty.push(8.0);
                    } else if tstr == "CAG".to_string() || tstr == "TAG".to_string() {
                        self.forward_acceptor_penalty.push(0.0);
                    } else {
                        self.forward_acceptor_penalty.push(30.0);
                    }

                    if tstr[1..3] == "AT".to_string() {
                        self.reverse_acceptor_penalty.push(21.0);
                    } else if tstr[1..3] == "GC".to_string() {
                        self.reverse_acceptor_penalty.push(15.0);
                    } else if tstr == "GAC".to_string() || tstr == "AAC".to_string() {
                        self.reverse_acceptor_penalty.push(8.0);
                    } else if tstr == "TAC".to_string() || tstr == "CAC".to_string() {
                        self.reverse_acceptor_penalty.push(0.0);
                    } else {
                        self.reverse_acceptor_penalty.push(30.0);
                    }
                }
            }
        }

        for i in 0..self.freq_vec.len() + 1 {
            if i < self.freq_vec.len() && self.freq_vec[i].i {
                self.forward_donor_penalty.push(30.0);
                // forward_acceptor_penalty.push(standed_penalty);
                self.reverse_donor_penalty.push(30.0);
                // reverse_acceptor_penalty.push(standed_penalty);
                continue;
            }
            if i + 2 >= self.freq_vec.len() {
                self.forward_donor_penalty.push(30.0);
                self.reverse_donor_penalty.push(30.0);
            } else {
                let mut j = i;
                ref_sliding_window.clear();
                while ref_sliding_window.len() < 3 && j < self.freq_vec.len() {
                    if !self.freq_vec[j].i {
                        ref_sliding_window.push(self.freq_vec[j].ref_base as u8);
                    }
                    j += 1;
                }
                let mut tstr = String::new();
                for c in ref_sliding_window.iter() {
                    tstr.push(*c as char);
                }
                if tstr.len() < 3 {
                    self.forward_donor_penalty.push(30.0);
                    self.reverse_donor_penalty.push(30.0);
                } else {
                    if tstr[0..2] == "AT".to_string() {
                        self.forward_donor_penalty.push(21.0);
                    } else if tstr[0..2] == "GC".to_string() {
                        self.forward_donor_penalty.push(15.0);
                    } else if tstr == "GTC".to_string() || tstr == "GTT".to_string() {
                        self.forward_donor_penalty.push(8.0);
                    } else if tstr == "GTA".to_string() || tstr == "GTG".to_string() {
                        self.forward_donor_penalty.push(0.0);
                    } else {
                        self.forward_donor_penalty.push(30.0);
                    }

                    if tstr[0..2] == "GT".to_string() {
                        self.reverse_donor_penalty.push(21.0);
                    } else if tstr == "CTT".to_string() || tstr == "CTC".to_string() {
                        self.reverse_donor_penalty.push(8.0);
                    } else if tstr == "CTG".to_string() || tstr == "CTA".to_string() {
                        self.reverse_donor_penalty.push(0.0);
                    } else {
                        self.reverse_donor_penalty.push(30.0);
                    }
                }
            }
        }
    }

    pub fn cal_intron_intervals(&mut self) {
        let mut extend_size = 96; // L>(d(i)+a(a)+telda(q)-q)/e
        let mut s: usize = 0;
        let mut e: usize = 0;
        for i in 0..self.freq_vec.len() {
            if self.freq_vec[i].get_depth_exclude_intron() == 0 && self.freq_vec[i].n > 0 {
                if extend_size > 0 {
                    extend_size -= 1;
                } else if s == 0 && e == 0 {
                    s = i;
                    e = i;
                } else {
                    e = i;
                }
            } else {
                if e > s {
                    self.intron_intervals.push(Interval { start: s, stop: e + 1, val: true });  // 0-based, left closed, right open
                }
                s = 0;
                e = 0;
                extend_size = 96;
            }
        }
    }
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