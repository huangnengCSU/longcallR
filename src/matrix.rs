use std::process;
use std::cmp::max;
use std::collections::HashMap;
use crate::bam_reader::Region;

#[derive(Clone)]
pub struct PileupMatrix {
    pub positions: HashMap<u32, u32>,
    // key: position on reference 0-based, value: index in matrix 0-based
    pub base_matrix: HashMap<String, Vec<u8>>,
    // key: read name, value: base sequence
    pub current_pos: i64,
    // current index of matrix, 0-based
    pub max_idx: i32,
    // max index of matrix, 0-based
    pub region: Region,
}

impl PileupMatrix {
    pub fn new() -> PileupMatrix {
        PileupMatrix {
            positions: HashMap::new(),
            base_matrix: HashMap::new(),
            current_pos: -1,
            max_idx: -1,
            region: Region::new("c:0-0".to_string()),
        }
    }

    pub fn insert(&mut self, readname: &String, base_seq: &[u8], pos: &u32) {
        if pos.clone() != self.current_pos as u32 {
            self.max_idx += 1;
            self.positions.insert(pos.clone(), self.max_idx as u32);
            self.current_pos = pos.clone() as i64;
        }

        if self.base_matrix.get(readname).is_none() {
            let mut base_vec = Vec::new();
            if self.max_idx == 0 {
                for i in 0..base_seq.len() {
                    base_vec.push(base_seq[i]);
                }
                self.base_matrix.insert(readname.clone(), base_vec);
            } else {
                for i in 0..self.positions.get(pos).unwrap().clone() as i32 {
                    base_vec.push(b' ');
                }
                for i in 0..base_seq.len() {
                    base_vec.push(base_seq[i]);
                }
                self.base_matrix.insert(readname.clone(), base_vec);
            }
        } else {
            for i in 0..base_seq.len() {
                self.base_matrix.get_mut(readname).unwrap().push(base_seq[i]);
            }
        }
    }

    pub fn expand(&mut self, seq_lengths: &HashMap<String, u32>, pos: &u32, max_insertion_size: i32) {
        for (readname, seq_len) in seq_lengths.iter() {
            let expand_size = max_insertion_size - seq_len.clone() as i32;
            if expand_size > 0 {
                for i in 0..expand_size {
                    self.base_matrix.get_mut(readname).unwrap().push(b' ');
                }
            }
        }
        self.max_idx += max_insertion_size - 1;
        self.positions.insert(pos.clone(), self.max_idx as u32);
    }

    pub fn padding(&mut self) {
        for (readname, base_vec) in self.base_matrix.iter_mut() {
            if base_vec.len() < self.max_idx as usize + 1 {
                for i in base_vec.len()..self.max_idx as usize + 1 {
                    base_vec.push(b' ');
                }
            }
        }
    }

    pub fn clear(&mut self) {
        self.positions.clear();
        self.base_matrix.clear();
        self.current_pos = -1;
        self.max_idx = -1;
        self.region = Region::new("c:0-0".to_string());
    }

    pub fn generate_column_profile(base_matrix: &HashMap<String, Vec<u8>>, column_base_counts: &mut Vec<ColumnBaseCount>) {
        assert!(base_matrix.len() > 0);
        let ncols = base_matrix.iter().next().unwrap().1.len();
        let mut ref_base: u8 = 0;
        let mut column_bases: Vec<u8> = Vec::new();
        for i in 0..ncols {
            for (readname, base_vec) in base_matrix.iter() {
                if *readname == "ref".to_string() {
                    ref_base = base_vec[i];
                    continue;
                }
                column_bases.push(base_vec[i]);
            }
            let cbc = ColumnBaseCount::new_from_column(&column_bases, ref_base);
            column_base_counts.push(cbc);
            column_bases.clear();
        }
    }

    pub fn generate_reduced_profile(base_matrix: &HashMap<String, Vec<u8>>,
                                    column_base_counts: &mut Vec<ColumnBaseCount>,
                                    column_indexes: &mut Vec<usize>,
                                    reduced_base_matrix: &mut HashMap<String, Vec<u8>>) {
        assert!(base_matrix.len() > 0);
        let ncols = base_matrix.iter().next().unwrap().1.len();
        let mut ref_base: u8 = 0;
        let mut column_bases: Vec<u8> = Vec::new();
        for i in 0..ncols {
            for (readname, base_vec) in base_matrix.iter() {
                if *readname == "ref".to_string() {
                    ref_base = base_vec[i];
                    continue;
                }
                column_bases.push(base_vec[i]);
            }
            let cbc = ColumnBaseCount::new_from_column(&column_bases, ref_base);
            if (cbc.n_a + cbc.n_c + cbc.n_g + cbc.n_t + cbc.n_dash) == 0 && cbc.n_n == 0 {
                column_bases.clear();
                continue;
            }
            column_base_counts.push(cbc);
            column_indexes.push(i);
            column_bases.clear();
        }
        reduced_base_matrix.clear();
        for i in column_indexes.iter() {
            for (readname, base_vec) in base_matrix.iter() {
                if reduced_base_matrix.get(readname).is_none() {
                    reduced_base_matrix.insert(readname.clone(), vec![base_vec[*i]]);
                } else {
                    reduced_base_matrix.get_mut(readname).unwrap().push(base_vec[*i]);
                }
            }
        }
    }
}


pub struct ColumnBaseCount {
    pub n_a: u16,
    pub n_c: u16,
    pub n_g: u16,
    pub n_t: u16,
    pub n_n: u16,
    pub n_dash: u16,
    pub n_blank: u16,
    pub max_count: u16,
    pub ref_base: u8,
}

impl ColumnBaseCount {
    pub fn new() -> ColumnBaseCount {
        ColumnBaseCount {
            n_a: 0,
            n_c: 0,
            n_g: 0,
            n_t: 0,
            n_n: 0,
            n_dash: 0,
            n_blank: 0,
            max_count: 0,
            ref_base: 0,
        }
    }

    pub fn new_from_column(matrix_column: &Vec<u8>, ref_base: u8) -> ColumnBaseCount {
        let mut cbc = ColumnBaseCount {
            n_a: 0,
            n_c: 0,
            n_g: 0,
            n_t: 0,
            n_n: 0,
            n_dash: 0,
            n_blank: 0,
            max_count: 0,
            ref_base: ref_base,
        };

        for b in matrix_column.iter() {
            match b {
                b'A' => cbc.n_a += 1,
                b'a' => cbc.n_a += 1,
                b'C' => cbc.n_c += 1,
                b'c' => cbc.n_c += 1,
                b'G' => cbc.n_g += 1,
                b'g' => cbc.n_g += 1,
                b'T' => cbc.n_t += 1,
                b't' => cbc.n_t += 1,
                b'N' => cbc.n_n += 1,
                b'n' => cbc.n_n += 1,
                b'-' => cbc.n_dash += 1,
                b' ' => cbc.n_blank += 1,
                _ => (),
            }
        }

        // TODO: reduce the effect of N when calculate score.
        cbc.n_n = 1;

        // max_count does not consider padding (blank).
        cbc.max_count = max(max(max(max(max(cbc.n_a, cbc.n_c), cbc.n_g), cbc.n_t), cbc.n_n), cbc.n_dash);
        return cbc;
    }

    pub fn get_major_base(&self) -> u8 {
        if self.n_a == self.max_count {
            return b'A';
        } else if self.n_c == self.max_count {
            return b'C';
        } else if self.n_g == self.max_count {
            return b'G';
        } else if self.n_t == self.max_count {
            return b'T';
        } else if self.n_n == self.max_count {
            return b'N';
        } else if self.n_dash == self.max_count {
            return b'-';
        } else if self.n_blank == self.max_count {
            return b' ';
        } else {
            println!("Error: get_major_base() failed, exit.");
            process::exit(1);
            return 0;
        }
    }

    pub fn get_ref_base(&self) -> u8 {
        self.ref_base
    }

    pub fn get_score1(&self, x: &u8) -> i32 {
        match x {
            b'A' => 1 - (self.n_a == self.max_count) as i32,
            b'a' => 1 - (self.n_a == self.max_count) as i32,
            b'C' => 1 - (self.n_c == self.max_count) as i32,
            b'c' => 1 - (self.n_c == self.max_count) as i32,
            b'G' => 1 - (self.n_g == self.max_count) as i32,
            b'g' => 1 - (self.n_g == self.max_count) as i32,
            b'T' => 1 - (self.n_t == self.max_count) as i32,
            b't' => 1 - (self.n_t == self.max_count) as i32,
            b'N' => 1 - (self.n_n == self.max_count) as i32,
            b'n' => 1 - (self.n_n == self.max_count) as i32,
            b'-' => 1 - (self.n_dash == self.max_count) as i32,
            _ => {
                println!("Error: get_score1() failed, exit.");
                process::exit(1);
                return 0;
            }
        }
    }

    pub fn get_score2(&self, x: &u8) -> f64 {
        let s: f64 = (self.n_a + self.n_c + self.n_g + self.n_t + self.n_dash) as f64;
        if s == 0.0 {
            return 0.0;
        } else {
            match x {
                b'A' => (s - self.n_a as f64) / s,
                b'a' => (s - self.n_a as f64) / s,
                b'C' => (s - self.n_c as f64) / s,
                b'c' => (s - self.n_c as f64) / s,
                b'G' => (s - self.n_g as f64) / s,
                b'g' => (s - self.n_g as f64) / s,
                b'T' => (s - self.n_t as f64) / s,
                b't' => (s - self.n_t as f64) / s,
                b'N' => (s - self.n_n as f64) / s,
                b'n' => (s - self.n_n as f64) / s,
                b'-' => (s - self.n_dash as f64) / s,
                _ => {
                    println!("Error: get_score2() failed, exit.");
                    process::exit(1);
                    return 0.0;
                }
            }
        }
    }

    pub fn get_score(&self, x: &u8) -> f64 {
        let s1 = self.get_score1(x) as f64;
        let s2 = self.get_score2(x);
        (s1 + s2) / 2.0
    }
}