use std::process;
use std::collections::HashMap;
use crate::bam_reader::Region;
use rust_htslib::{bam, bam::Read, bam::record::Cigar, bam::record::CigarString, bam::Format};
use seq_io::fasta::{Reader, Record};
use rust_lapper::{Interval, Lapper};
use crate::align2::SpliceMatrixElement;

#[derive(Debug, Clone)]
pub struct ParsedRead {
    pub bam_record: bam::record::Record,
    // read from bam file
    pub parsed_seq: Vec<u8>,
    // parsed sequence, expand insertion
    pub pos_on_profile: i64,
    // start position of the first aligned base on the profile, 0-based
    pub alignment_score: f64,
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
            panic!("read start position exceeds the profile region!\n Readname: {}, ref_start: {}, region: {}-{}", std::str::from_utf8(self.bam_record.qname()).unwrap(), ref_start, profile.region.start, profile.region.end);
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

        // println!("readname:{}, i = {}, j = {}, refbase: {}", std::str::from_utf8(self.bam_record.qname()).unwrap(), i, j, profile.freq_vec[i].ref_base);

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
                        assert_eq!(profile.freq_vec[i].i, true, "{},{},{}", std::str::from_utf8(self.bam_record.qname()).unwrap(), self.bam_record.pos(), i);
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

    pub fn update_bam_record(&mut self, profile: &Profile) {
        /*
        Update the bam record with the parsed sequence after realignment
        profile: profile without any removal
        */
        let mut qb: char;
        let mut tb: char;
        let p = self.pos_on_profile;
        let mut shift = 0;
        let mut head_flag = true;
        let mut pre_state: u8 = 0;  // 1: M, 2:I, 3:D, 4:N
        let mut op_len = 0;

        let mut left_soft_clip = 0;
        let mut right_soft_clip = 0;
        let mut left_hard_clip = 0;
        let mut right_hard_clip = 0;
        let mut new_cigar: Vec<Cigar> = Vec::new();

        let cg = self.bam_record.cigar().iter().next().unwrap().clone();
        if cg.char() == 'S' {
            left_soft_clip = cg.len();
            new_cigar.push(Cigar::SoftClip(left_soft_clip));
        } else if cg.char() == 'H' {
            left_hard_clip = cg.len();
            new_cigar.push(Cigar::HardClip(left_hard_clip));
        }

        let cg = self.bam_record.cigar().iter().last().unwrap().clone();
        if cg.char() == 'S' {
            right_soft_clip = cg.len();
        } else if cg.char() == 'H' {
            right_hard_clip = cg.len();
        }

        for i in 0..self.parsed_seq.len() {
            qb = self.parsed_seq[i] as char;
            tb = profile.freq_vec[p as usize + i].ref_base;
            if qb != '-' && qb != '*' {
                head_flag = false;
            }
            if head_flag {
                if tb != '-' {
                    // non insertion
                    shift += 1;
                }
            } else {
                // M
                let mut cur_state: u8;  // 1: M, 2:I, 3:D, 4:N

                if qb == '*' {
                    continue;
                }

                // 1: M
                if qb != '-' && qb != 'N' && tb != '-' && tb != 'N' {
                    cur_state = 1;
                    if pre_state == 0 {
                        pre_state = cur_state;
                        op_len += 1;
                    } else {
                        if cur_state == pre_state {
                            op_len += 1;
                        } else {
                            match pre_state {
                                1 => {
                                    new_cigar.push(Cigar::Match(op_len));
                                }
                                2 => {
                                    new_cigar.push(Cigar::Ins(op_len));
                                }
                                3 => {
                                    new_cigar.push(Cigar::Del(op_len));
                                }
                                4 => {
                                    new_cigar.push(Cigar::RefSkip(op_len));
                                }
                                _ => {
                                    panic!("Invalid pre_state: {}", pre_state);
                                }
                            }
                            pre_state = cur_state;
                            op_len = 1;
                        }
                    }
                }
                // 2: I
                if tb == '-' {
                    if qb != 'N' && qb != '-' {
                        cur_state = 2;
                        if cur_state == pre_state {
                            op_len += 1;
                        } else {
                            match pre_state {
                                1 => {
                                    new_cigar.push(Cigar::Match(op_len));
                                }
                                2 => {
                                    new_cigar.push(Cigar::Ins(op_len));
                                }
                                3 => {
                                    new_cigar.push(Cigar::Del(op_len));
                                }
                                4 => {
                                    new_cigar.push(Cigar::RefSkip(op_len));
                                }
                                _ => {
                                    panic!("Invalid pre_state: {}", pre_state);
                                }
                            }
                            pre_state = cur_state;
                            op_len = 1;
                        }
                    } else if qb == 'N' {
                        // insertion in intron
                        continue;
                    } else {
                        panic!("Unknow state: qb = {}, tb = {}", qb, tb);
                    }
                }

                // 3: D
                if qb == '-' && tb != '-' {
                    cur_state = 3;
                    if cur_state == pre_state {
                        op_len += 1;
                    } else {
                        match pre_state {
                            1 => {
                                new_cigar.push(Cigar::Match(op_len));
                            }
                            2 => {
                                new_cigar.push(Cigar::Ins(op_len));
                            }
                            3 => {
                                new_cigar.push(Cigar::Del(op_len));
                            }
                            4 => {
                                new_cigar.push(Cigar::RefSkip(op_len));
                            }
                            _ => {
                                panic!("Invalid pre_state: {}", pre_state);
                            }
                        }
                        pre_state = cur_state;
                        op_len = 1;
                    }
                }

                // 4: N
                if qb == 'N' {
                    if tb == '-' {
                        continue;
                    } else {
                        cur_state = 4;
                        if cur_state == pre_state {
                            op_len += 1;
                        } else {
                            match pre_state {
                                1 => {
                                    new_cigar.push(Cigar::Match(op_len));
                                }
                                2 => {
                                    new_cigar.push(Cigar::Ins(op_len));
                                }
                                3 => {
                                    new_cigar.push(Cigar::Del(op_len));
                                }
                                4 => {
                                    new_cigar.push(Cigar::RefSkip(op_len));
                                }
                                _ => {
                                    panic!("Invalid pre_state: {}", pre_state);
                                }
                            }
                            pre_state = cur_state;
                            op_len = 1;
                        }
                    }
                }
            }
        }

        match pre_state {
            1 => {
                new_cigar.push(Cigar::Match(op_len));
            }
            2 => {
                new_cigar.push(Cigar::Ins(op_len));
            }
            3 => {
                new_cigar.push(Cigar::Del(op_len));
            }
            4 => {
                new_cigar.push(Cigar::RefSkip(op_len));
            }
            _ => {
                panic!("Invalid pre_state: {}", pre_state);
            }
        }

        if right_hard_clip > 0 {
            new_cigar.push(Cigar::HardClip(right_hard_clip));
        } else if right_soft_clip > 0 {
            new_cigar.push(Cigar::SoftClip(right_soft_clip));
        }

        let new_cigar_string = CigarString(new_cigar);
        let mut cglen = 0;
        for cg in new_cigar_string.iter() {
            if cg.char() == 'M' || cg.char() == 'I' {
                cglen += cg.len();
            }
        }
        let num_base = std::str::from_utf8(self.parsed_seq.as_slice()).unwrap()
            .replace("*", "")
            .replace("N", "")
            .replace("-", "").len() as u32; // does not contain left and right soft clip, intron.

        // println!();
        // println!("readname: {}, cigar len: {} read len: {}",
        //          std::str::from_utf8(self.bam_record.qname()).unwrap(),
        //          cglen,
        //          num_base
        // );
        if cglen != num_base {
            println!("cglen error: {}, cigar len: {}, num base: {}", std::str::from_utf8(self.bam_record.qname()).unwrap(), cglen, num_base);
            println!("{}", new_cigar_string.to_string());
        }
        assert!(cglen == num_base);

        let record = bam::Record::from(self.bam_record.clone());
        self.bam_record.set_pos(record.pos() + shift as i64);
        self.bam_record.set(record.qname(),
                            Some(&new_cigar_string),
                            record.seq().as_bytes().as_slice(),
                            record.qual());
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
    pub ref_base: char,
    // reference base, A,C,G,T,N,-
    pub intron: bool, // whether this position is an intron
}

impl BaseFreq {
    pub fn subtract(&mut self, base: u8) {
        match base {
            b'A' => {
                assert!(self.a > 0);
                self.a -= 1;
            }
            b'C' => {
                assert!(self.c > 0);
                self.c -= 1;
            }
            b'G' => {
                assert!(self.g > 0);
                self.g -= 1;
            }
            b'T' => {
                assert!(self.t > 0);
                self.t -= 1;
            }
            b'N' => {
                assert!(self.n > 0);
                self.n -= 1;
            }
            b'-' => {
                assert!(self.d > 0);
                self.d -= 1;
            }
            b'*' => {
                return;
            }
            _ => {
                panic!("Invalid base: {}", base as char);
            }
        }
    }

    pub fn add(&mut self, base: u8) {
        match base {
            b'A' => {
                self.a += 1;
            }
            b'C' => {
                self.c += 1;
            }
            b'G' => {
                self.g += 1;
            }
            b'T' => {
                self.t += 1;
            }
            b'N' => {
                self.n += 1;
            }
            b'-' => {
                self.d += 1;
            }
            b'*' => {
                return;
            }
            _ => {
                panic!("Invalid base: {}", base as char);
            }
        }
    }

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
        let mut pileups = bam.pileup();
        pileups.set_max_depth(1000000);
        for p in pileups {
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
                                    insert_bf.push(BaseFreq { a: 0, c: 0, g: 0, t: 0, n: 0, d: 0, i: true, ref_base: '\x00', intron: false });   // fall in insertion
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
                                    insert_bf.push(BaseFreq { a: 0, c: 0, g: 0, t: 0, n: 0, d: 0, i: true, ref_base: '\x00', intron: false });   // fall in insertion
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
            // update the number of dash in insertion position
            let cur_depth_with_intron = bf.get_depth_include_intron();
            self.freq_vec.push(bf);
            if insert_bf.len() > 0 {
                for tmpi in 0..insert_bf.len() {
                    insert_bf[tmpi].d = cur_depth_with_intron - insert_bf[tmpi].get_depth_include_intron();
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
        for iv in self.intron_intervals.iter() {
            for i in iv.start..iv.stop {
                self.freq_vec[i].intron = true;
            }
        }
    }

    pub fn generate_read_profile_for_realign(&mut self, parsed_read: &ParsedRead) -> (i64, i64, Vec<u8>, Profile) {
        /*
        1. generate parsed seq without intron and profile without intron as well
        2. subtract parsed seq from profile
        returns: (s, e, parsed_seq_without_intron, profile_without_intron)
        s: start position of parsed read without introns on the profile, 0-based, inclusive
        e: end position of parsed read without introns on the profile, 0-based, inclusive
        parsed_seq_without_intron: parsed sequence without introns
        profile_without_intron: profile without introns
         */
        let p = parsed_read.pos_on_profile as usize;
        let mut s: i64 = -1;   // 0-based, include
        let mut e: i64 = -1;   // 0-based, include
        let mut parsed_seq_without_intron: Vec<u8> = Vec::new();
        let mut profile_without_intron = Profile::default();
        for i in 0..parsed_read.parsed_seq.len() {
            let mut j = p + i;
            if !self.freq_vec[j].intron {
                if s == -1 && e == -1 {
                    s = j as i64;
                }
                parsed_seq_without_intron.push(parsed_read.parsed_seq[i]);
                self.freq_vec[j].subtract(parsed_read.parsed_seq[i]);   // minus parsed base from freq_vec
                profile_without_intron.freq_vec.push(self.freq_vec[j].clone());
                profile_without_intron.forward_acceptor_penalty.push(self.forward_acceptor_penalty[j]);
                profile_without_intron.forward_donor_penalty.push(self.forward_donor_penalty[j]);
                profile_without_intron.reverse_acceptor_penalty.push(self.reverse_acceptor_penalty[j]);
                profile_without_intron.reverse_donor_penalty.push(self.reverse_donor_penalty[j]);
                e = j as i64;
            }
        }
        profile_without_intron.forward_acceptor_penalty.push(self.forward_acceptor_penalty[e as usize + 1]); // acceptor penalty stored in the next position
        profile_without_intron.reverse_acceptor_penalty.push(self.reverse_acceptor_penalty[e as usize + 1]); // acceptor penalty stored in the next position
        (s, e, parsed_seq_without_intron, profile_without_intron)
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
        let parsed_read = ParsedRead { bam_record: record.clone(), parsed_seq: Vec::new(), pos_on_profile: -1, alignment_score: 0.0 };
        parsed_reads.insert(qname, parsed_read);
    }
    return parsed_reads;
}

pub fn get_bam_header(input_bam: &str) -> bam::header::Header {
    let bam: bam::IndexedReader = bam::IndexedReader::from_path(input_bam).unwrap();
    bam::Header::from_template(bam.header())
}

pub fn write_bam(out_bam_path: &str, parsed_reads: &HashMap<String, ParsedRead>, bam_header: &bam::header::Header) {
    let mut bam_writer = bam::Writer::from_path(out_bam_path, bam_header, Format::Bam).unwrap();
    for (_, pr) in parsed_reads.iter() {
        let re = bam_writer.write(&pr.bam_record);
        if re != Ok(()) {
            println!("write failed");
            process::exit(1);
        }
    }
}

pub fn banded_nw(query: &Vec<u8>, profile: &Profile, width: usize, reverse_strand: bool) -> (f64, Vec<u8>, Vec<u8>) {
    /*
    Needleman-Wunsch alignment with banding.
    query: query sequence without introns
    profile: profile without introns
    width: band width
    reverse_strand: use forward donor/acceptor penalty or reverse donor/acceptor penalty
    */
    let h = 2.0; // short gap open, q in minimap2
    let g = 1.0; // short gap extend, e in minimap2
    let h2 = 32.0; // long gap open, q hat in minimap2
    // let p = 9.0; // splice junction site penalty
    let b: f64 = 34.0; // splice junction site penalty
    let match_score = 2.0;
    let mismatch_score = -2.0;
    let standard_donor_penalty: &Vec<f64>;
    let standard_acceptor_penalty: &Vec<f64>;

    if reverse_strand {
        standard_donor_penalty = &profile.reverse_donor_penalty;
        standard_acceptor_penalty = &profile.reverse_acceptor_penalty;
    } else {
        standard_donor_penalty = &profile.forward_donor_penalty;
        standard_acceptor_penalty = &profile.forward_acceptor_penalty;
    }

    // let query_without_gap = query.replace("-", "").replace("N", "");
    let query_t = std::str::from_utf8(query.clone().as_slice()).unwrap().replace("-", "").replace("N", "").replace("*", "");
    let query_without_gap = query_t.as_bytes();

    let q_len = query_without_gap.len();
    let t_len = profile.freq_vec.len();

    // for d in reduced_donor_penalty.iter() {
    //     print!("{}\t", d);
    // }
    // println!();
    //
    // for d in reduced_acceptor_penalty.iter() {
    //     print!("{}\t", d);
    // }
    // println!();
    //
    // for d in 0..t_len {
    //     print!("{}\t", profile[d].get_ref_base() as char);
    // }
    // println!();
    //
    // for i in 0..t_len {
    //     print!("{}\t", profile[i].n_n);
    // }
    // println!();
    //
    // for i in 0..t_len {
    //     print!("{}\t", profile[i].get_depth() + profile[i].n_n);
    // }
    // println!();


    let mut mat: Vec<Vec<SpliceMatrixElement>> = vec![vec![SpliceMatrixElement { ..Default::default() }; (2 * width) + 3]; t_len + 1];
    // println!("mat size: {} x {}", mat.len(), mat[0].len());

    let mut k_vec: Vec<usize> = vec![0; t_len + 1];

    // Initialize first row
    mat[0][width + 1].ix = -h - g;
    mat[0][width + 1].ix2 = -f64::INFINITY;
    mat[0][width + 1].m = mat[0][width + 1].ix;

    let mut i = 1; // index on target
    let mut j = 1; // index on query
    let mut k = 0; // index on query_without_gap

    let mut left_bound: usize; // include this index
    let mut right_bound: usize; // exclude this index
    while i < t_len + 1 && k < q_len + 1 {
        let mut offset: usize;
        if query[j - 1] != b'*' && query[j - 1] != b'-' && query[j - 1] != b'N' {
            k += 1;
            offset = 0;
        } else {
            offset = 1;
        }
        k_vec[i] = k; // store k index of each row

        if k as i32 - width as i32 > 0 {
            left_bound = (k as i32 - width as i32) as usize;
        } else {
            left_bound = 1;
        }
        right_bound = if k + width + 1 <= q_len + 1 {
            k + width + 1
        } else {
            q_len + 1
        };
        for u in left_bound..right_bound {
            // index on alignment matrix: (i,v), v = (w+1)+(u-k), u in [left_bound, right_bound)
            let v = (width as i32 + 1 + (u as i32 - k as i32)) as usize;
            let col = &profile.freq_vec[i - 1];
            let qbase = query_without_gap[u - 1];
            let sij: f64;
            if col.get_depth_exclude_intron() <= 5 {
                let tbase = col.ref_base as u8;
                if tbase == b'-' {
                    sij = match_score;
                } else if qbase == tbase {
                    sij = match_score;
                } else {
                    sij = mismatch_score;
                }
            } else {
                sij = 2.0 - 4.0 * col.get_score(qbase);
            }

            // position specific donor penalty (mat[i], taget[i-1], profile[i-1])
            let mut dp_f = 1.0;
            let mut ap_f = 1.0;
            if i as i32 - 1 >= 0 && i < t_len {
                //if profile[i - 1].get_ref_base() == b'-' {
                if profile.freq_vec[i - 1].i {
                    dp_f = 1.0;
                } else {
                    let (curr_intron_cnt, curr_intron_percentage) = profile.freq_vec[i - 1].get_intron_ratio();
                    let mut ai = i as i32 - 2;
                    while ai >= 0 {
                        if profile.freq_vec[ai as usize].i {
                            ai -= 1;
                        } else {
                            break;
                        }
                    }
                    if ai >= 0 {
                        let (prev_intron_cnt, prev_intron_percentage) = profile.freq_vec[ai as usize].get_intron_ratio();
                        // if (curr_intron_cnt as i32 - prev_intron_cnt as i32).abs() > 10 {
                        //     dp_f = 0.0;
                        // } else {
                        //     dp_f = 1.0 - (curr_intron_percentage - prev_intron_percentage).abs();
                        // }
                        dp_f = 1.0 - (curr_intron_percentage - prev_intron_percentage).abs();
                    }
                }
            }

            // position specific acceptor penalty (mat[i], taget[i-1], profile[i-1]), previous position of current target position i-1, and not dash on reference
            if i as i32 - 1 >= 0 && i < t_len {
                if profile.freq_vec[i - 1].i {
                    ap_f = 1.0;
                } else {
                    let (curr_intron_cnt, curr_intron_percentage) = profile.freq_vec[i - 1].get_intron_ratio();
                    let mut ai = i;
                    while ai < t_len {
                        if profile.freq_vec[ai].i {
                            ai += 1;
                        } else {
                            break;
                        }
                    }
                    if ai < t_len {
                        let (later_intron_cnt, later_intron_percentage) = profile.freq_vec[ai].get_intron_ratio();
                        // if (later_intron_cnt as i32 - curr_intron_cnt as i32).abs() > 10 {
                        //     ap_f = 0.0;
                        // } else {
                        //     ap_f = 1.0 - (later_intron_percentage - curr_intron_percentage).abs();
                        // }
                        ap_f = 1.0 - (later_intron_percentage - curr_intron_percentage).abs();
                    }
                }
            }

            // println!("dp_f: {}, ap_f: {}", dp_f, ap_f);
            if col.i || col.ref_base == 'N' {
                if mat[i - 1][v + 1 - offset].m >= mat[i - 1][v + 1 - offset].ix {
                    mat[i][v].ix = mat[i - 1][v + 1 - offset].m;
                    mat[i][v].ix_s = 1; //mat[i][v].ix_prev_m = true;
                } else {
                    mat[i][v].ix = mat[i - 1][v + 1 - offset].ix;
                    mat[i][v].ix_s = 2; //mat[i][v].ix_prev_ix = true;
                }
            } else {
                if !col.i && i as i32 - 2 >= 0 && profile.freq_vec[i - 2].i {
                    if mat[i - 1][v + 1 - offset].m - h - g >= mat[i - 1][v + 1 - offset].ix - h - g {
                        mat[i][v].ix = mat[i - 1][v + 1 - offset].m - h - g;
                        mat[i][v].ix_s = 1; //mat[i][v].ix_prev_m = true;
                    } else {
                        mat[i][v].ix = mat[i - 1][v + 1 - offset].ix - h - g;
                        mat[i][v].ix_s = 2; //mat[i][v].ix_prev_ix = true;
                    }
                } else {
                    if mat[i - 1][v + 1 - offset].m - h - g >= mat[i - 1][v + 1 - offset].ix - g {
                        mat[i][v].ix = mat[i - 1][v + 1 - offset].m - h - g;
                        mat[i][v].ix_s = 1; //mat[i][v].ix_prev_m = true;
                    } else {
                        mat[i][v].ix = mat[i - 1][v + 1 - offset].ix - g;
                        mat[i][v].ix_s = 2; //mat[i][v].ix_prev_ix = true;
                    }
                }
            }

            if mat[i - 1][v + 1 - offset].m - standard_donor_penalty[i - 1] - h2 - dp_f * 28.0 >= mat[i - 1][v + 1 - offset].ix2 {
                mat[i][v].ix2 = mat[i - 1][v + 1 - offset].m - standard_donor_penalty[i - 1] - h2 - dp_f * 28.0;
                mat[i][v].ix2_s = 1; //mat[i][v].ix2_prev_m = true;
            } else {
                mat[i][v].ix2 = mat[i - 1][v + 1 - offset].ix2;
                mat[i][v].ix2_s = 3; //mat[i][v].ix2_prev_ix2 = true;
            }


            mat[i][v].m = (mat[i - 1][v - offset].m + sij)
                .max(mat[i][v].ix)
                .max(mat[i][v].ix2 - standard_acceptor_penalty[i] - ap_f * 28.0); // index i in matrix is corresponding to the index i-1 in reference (donor penalty and acceptor penalty)
            if mat[i][v].m == mat[i - 1][v - offset].m + sij {
                mat[i][v].m_s = 1; //mat[i][v].m_prev_m = true;
            } else if mat[i][v].m == mat[i][v].ix {
                mat[i][v].m_s = 2; //mat[i][v].m_prev_ix = true;
            } else if mat[i][v].m == mat[i][v].ix2 - standard_acceptor_penalty[i] - ap_f * 28.0 {
                mat[i][v].m_s = 3; //mat[i][v].m_prev_ix2 = true;
            }
        }
        i += 1;
        j += 1;
    }

    // trace back
    let mut aligned_query: Vec<u8> = Vec::new();
    let mut ref_target: Vec<u8> = Vec::new();
    let mut alignment_score = 0.0;

    let mut u: usize;
    let mut v: usize;

    i = t_len;

    // let mut score_vec = Vec::new();
    // for vv in 0..2 * width + 3 {
    //     score_vec.push(mat[i][vv].m);
    // }
    // find max value and index in score_vec
    // let (max_score, max_index) = find_max_value_and_index(&score_vec);
    let vv = (width as i32 + 1 + (q_len as i32 - k_vec[i] as i32)) as usize;
    let max_score = mat[i][vv].m;

    alignment_score = max_score;
    v = vv;
    k = k_vec[i];
    u = (v as i32 - (width as i32 + 1) + k as i32) as usize; // index on query_without_gap

    let mut trace_back_stat: u8 = 0;
    if mat[i][v].m_s == 2 {
        //mat[i][v].m_prev_ix
        trace_back_stat = 2;
        // trace_back_stat = TraceBack::IX;
    } else if mat[i][v].m_s == 3 {
        //mat[i][v].m_prev_ix2
        trace_back_stat = 3;
        // trace_back_stat = TraceBack::IX2;
    } else if mat[i][v].m_s == 1 {
        // mat[i][v].m_prev_m
        trace_back_stat = 1;
        // trace_back_stat = TraceBack::M;
    } else {
        panic!("Error: no traceback");
    }

    while i > 0 && u > 0 {
        // println!("i: {}, k: {}, m:{}, ix: {}, iy:{}, ix2:{}", i, k, mat[i][k].m, mat[i][k].ix, mat[i][k].iy, mat[i][k].ix2);
        k = k_vec[i];
        v = (width as i32 + 1 + (u as i32 - k as i32)) as usize;
        let qbase = query_without_gap[u - 1];
        let ref_base = profile.freq_vec[i - 1].ref_base;
        if trace_back_stat == 1 {
            // trace_back_stat == TraceBack::M
            if mat[i][v].m_s == 2 {
                //mat[i][v].m_prev_ix
                trace_back_stat = 2;
                // trace_back_stat = TraceBack::IX;
            } else if mat[i][v].m_s == 3 {
                // mat[i][v].m_prev_ix2
                trace_back_stat = 3;
                // trace_back_stat = TraceBack::IX2;
            } else if mat[i][v].m_s == 1 {
                // mat[i][v].m_prev_m
                aligned_query.push(qbase);
                ref_target.push(ref_base as u8);
                i -= 1;
                u -= 1;
                trace_back_stat = 1;
                // trace_back_stat = TraceBack::M;
            } else {
                panic!("Error: no traceback");
            }
        } else if trace_back_stat == 2 {
            // trace_back_stat == TraceBack::IX
            if mat[i][v].ix_s == 2 {
                // mat[i][v].ix_prev_ix
                aligned_query.push(b'-');
                ref_target.push(ref_base as u8);
                i -= 1;
                trace_back_stat = 2;
                // trace_back_stat = TraceBack::IX;
            } else if mat[i][v].ix_s == 1 {
                //mat[i][v].ix_prev_m
                if ref_base == '-' {
                    // current position is in the insertion region
                    aligned_query.push(b'*');
                    ref_target.push(ref_base as u8);
                } else {
                    aligned_query.push(b'-');
                    ref_target.push(ref_base as u8);
                }
                i -= 1;
                trace_back_stat = 1;
                // trace_back_stat = TraceBack::M;
            } else {
                panic!("Error: no traceback");
            }
        } else if trace_back_stat == 3 {
            // trace_back_stat == TraceBack::IX2
            if mat[i][v].ix2_s == 3 {
                // mat[i][v].ix2_prev_ix2
                aligned_query.push(b'N');
                ref_target.push(ref_base as u8);
                // major_target.push(major_base);
                i -= 1;
                trace_back_stat = 3;
                // trace_back_stat = TraceBack::IX2;
            } else if mat[i][v].ix2_s == 1 {
                // mat[i][v].ix2_prev_m
                aligned_query.push(b'N');
                ref_target.push(ref_base as u8);
                i -= 1;
                trace_back_stat = 1;
                // trace_back_stat = TraceBack::M;
            } else {
                panic!("Error: no traceback");
            }
        } else {
            panic!("Error: no traceback");
        }
    }

    while i > 0 {
        let ref_base = profile.freq_vec[i - 1].ref_base;
        // let major_base = profile[i - 1].get_major_base();
        if ref_base == '-' {
            aligned_query.push(b'*');
            ref_target.push(ref_base as u8);
        } else {
            aligned_query.push(b'-');
            ref_target.push(ref_base as u8);
        }
        // major_target.push(major_base);
        i -= 1;
    }
    while u > 0 {
        let qbase = query_without_gap[u - 1];
        aligned_query.push(qbase);
        ref_target.push(b'-');
        // major_target.push(b'-');
        u -= 1;
    }

    aligned_query.reverse();
    ref_target.reverse();
    // major_target.reverse();
    // println!("original:\n {:?}", std::str::from_utf8(query).unwrap());
    // println!("ref:\n {:?}", std::str::from_utf8(ref_target.as_slice()).unwrap());
    // println!("aligned:\n {:?}", std::str::from_utf8(aligned_query.as_slice()).unwrap());
    (alignment_score, aligned_query, ref_target)
}

pub fn realign(profile: &mut Profile, parsed_reads: &mut HashMap<String, ParsedRead>, readnames: &Vec<String>) {
    /*
    Realign reads based on the profile.
    profile: complete profile without any removal
    parsed_reads: parsed reads without any removal
    readnames: read names to be realigned
    */
    let mut prev_score: f64 = f64::NEG_INFINITY;
    let mut max_iter = 10;
    while max_iter > 0 {
        let mut total_score = 0.0;
        for rname in readnames.iter() {
            let pr = parsed_reads.get(rname).unwrap();
            let (s, e, read_without_intron, profile_without_intron) = profile.generate_read_profile_for_realign(pr);
            let (f_score, f_query, f_target) = banded_nw(&read_without_intron, &profile_without_intron, 10, false);
            let (r_score, r_query, r_target) = banded_nw(&read_without_intron, &profile_without_intron, 10, true);
            // println!("s = {}, e = {}, rname = {}, f_score = {}, r_score = {}", s, e, rname, f_score, r_score);
            if r_score > f_score {
                // println!("reverse strand");
                if r_score > pr.alignment_score || pr.alignment_score == 0.0 {
                    update_profile_parsed_reads(profile, parsed_reads, &rname, &r_query, &r_target);
                    parsed_reads.get_mut(rname).unwrap().alignment_score = r_score;
                } else {
                    // println!("update profile with original parsed seq");
                    update_profile_parsed_reads(profile, parsed_reads, &rname, &read_without_intron, &r_target);
                }
                total_score += r_score;
            } else {
                // println!("forward strand");
                if f_score > pr.alignment_score || pr.alignment_score == 0.0 {
                    update_profile_parsed_reads(profile, parsed_reads, &rname, &f_query, &f_target);
                    parsed_reads.get_mut(rname).unwrap().alignment_score = f_score;
                } else {
                    // println!("update profile with original parsed seq");
                    update_profile_parsed_reads(profile, parsed_reads, &rname, &read_without_intron, &r_target);
                }
                total_score += f_score;
            }
        }
        max_iter -= 1;
        // println!("iteration = {}, total_score = {}, prev_score = {}", 10 - max_iter, total_score, prev_score);
        if total_score <= prev_score {
            break;
        } else {
            prev_score = total_score;
        }
    }

    for rname in readnames.iter() {
        let mpr = parsed_reads.get_mut(rname).unwrap();
        mpr.update_bam_record(&profile);
    }
}

pub fn update_profile_parsed_reads(profile: &mut Profile, parsed_reads: &mut HashMap<String, ParsedRead>, qname: &String, aligned_query: &Vec<u8>, aligned_target: &Vec<u8>) {
    /*
    Update profile and parsed_reads based on the alignment result.
    profile: complete profile without any removal
    parsed_reads: parsed reads without any removal
    qname: query name
    aligned_query: aligned query after removing introns
    aligned_target: aligned target after removing introns
    */
    let mut new_parsed_seq: Vec<u8> = Vec::new();
    let mut front_gap_num = 0;
    let mut non_gap_flag = false;
    let p = parsed_reads.get(qname).unwrap().pos_on_profile;
    let mut i = 0;  // index on aligned_query
    let mut j = p;  // index on profile
    while i < aligned_query.len() {
        if !non_gap_flag {
            if aligned_query[i] == b'-' || aligned_query[i] == b'*' {
                front_gap_num += 1;
                assert!(profile.freq_vec[j as usize].ref_base == aligned_target[i] as char);
                j += 1;
            } else {
                non_gap_flag = true;
                if !profile.freq_vec[j as usize].intron {
                    new_parsed_seq.push(aligned_query[i]);
                    profile.freq_vec[j as usize].add(aligned_query[i]); // add base to freq_vec
                    assert!(profile.freq_vec[j as usize].ref_base == aligned_target[i] as char);
                    i += 1;
                    j += 1;
                } else {
                    while j < profile.freq_vec.len() as i64 {
                        if profile.freq_vec[j as usize].intron {
                            new_parsed_seq.push(b'N');
                            j += 1;
                        } else {
                            break;
                        }
                    }
                    new_parsed_seq.push(aligned_query[i]);
                    profile.freq_vec[j as usize].add(aligned_query[i]); // add base to freq_vec
                    assert!(profile.freq_vec[j as usize].ref_base == aligned_target[i] as char);
                    i += 1;
                    j += 1;
                }
            }
        } else {
            if !profile.freq_vec[j as usize].intron {
                new_parsed_seq.push(aligned_query[i]);
                profile.freq_vec[j as usize].add(aligned_query[i]); // add base to freq_vec
                assert!(profile.freq_vec[j as usize].ref_base == aligned_target[i] as char);
                i += 1;
                j += 1;
            } else {
                while j < profile.freq_vec.len() as i64 {
                    if profile.freq_vec[j as usize].intron {
                        new_parsed_seq.push(b'N');
                        j += 1;
                    } else {
                        break;
                    }
                }
                new_parsed_seq.push(aligned_query[i]);
                profile.freq_vec[j as usize].add(aligned_query[i]); // add base to freq_vec
                assert!(profile.freq_vec[j as usize].ref_base == aligned_target[i] as char);
                i += 1;
                j += 1;
            }
        }
    }
    // if front_gap_num > 0 {
    //     parsed_reads.get_mut(qname).unwrap().pos_on_profile = p + front_gap_num;
    //     println!("read = {} front_gap_num = {}", qname, front_gap_num);
    // }
    parsed_reads.get_mut(qname).unwrap().parsed_seq = new_parsed_seq;
}