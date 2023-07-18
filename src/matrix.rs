use std::collections::HashMap;
use bam_reader::Region;

pub struct PileupMatrix {
    positions: HashMap<u32, u32>,
    // key: position on reference, value: index in matrix
    base_matrix: HashMap<String, vec<u8>>,
    // key: read name, value: base sequence
    current_pos: i64,
    // current index of matrix
    max_idx: i64,
    region: Region,
}

impl PileupMatrix {
    pub fn new() -> PileupMatrix {
        PileupMatrix {
            positions: HashMap::new(),
            base_matrix: HashMap::new(),
            current_pos: -1,
            max_idx: -1,
            region: Region::new("a:b-c".to_string()),
        }
    }

    pub fn insert(&self, readname: &String, base_seq: &[u8], pos: &u32) {
        if pos as i64 != self.current_pos {
            self.max_idx += 1;
            self.positions.insert(&pos, self.max_idx);
            self.current_pos = pos as i64;
        }

        if self.base_matrix.get(readname).is_none() {
            let base_vec = Vec::new();
            if self.max_idx == 0 {
                for i in 0..base_seq.len() {
                    base_vec.push(base_seq[i]);
                }
                self.base_matrix.insert(readname, base_vec);
            } else {
                for i in 0..self.positions.get(&pos) {
                    base_vec.push(' ');
                }
                for i in 0..base_seq.len() {
                    base_vec.push(base_seq[i]);
                }
            }
        } else {
            for i in 0..base_seq.len() {
                self.base_matrix.get_mut(readname).push(base_seq[i]);
            }
        }
    }

    pub fn expand(&self, seq_lengths: &HashMap<String, u32>, pos: &u32, max_insertion_size: u32) {
        for (readname, seq_len) in seq_lengths.iter_mut() {
            let expand_size = max_insertion_size - seq_len;
            if expand_size > 0 {
                for i in 0..expand_size {
                    self.base_matrix.get_mut(readname).push(' ');
                }
            }
        }
        self.max_idx += max_insertion_size - 1;
        self.positions[pos] = self.max_idx;
    }

    pub fn padding(&self) {
        for (readname, base_vec) in self.base_matrix.iter_mut() {
            if base_vec.len() < self.max_idx + 1 {
                for i in base_vec.len()..self.max_idx + 1 {
                    base_vec.push(' ');
                }
            }
        }
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