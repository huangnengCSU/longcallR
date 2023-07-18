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
}


pub struct ColumnBaseCount {
    n_a: u16,
    n_c: u16,
    n_g: u16,
    n_t: u16,
    n_n: u16,
    n_dash: u16,
    n_blank: u16,
    max_count: u16,
    ref_base: u8,
}