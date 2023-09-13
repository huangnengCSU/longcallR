use rust_htslib::{bam, bam::Read, bam::Record, bam::Format};
use rust_htslib::bam::record::{Cigar, CigarString};
use rust_htslib::bam::ext::BamRecordExtensions;
use bio::io::fasta;
use std::collections::HashMap;
use crate::bam_reader::Region;
use crate::matrix::ColumnBaseCount;
// use crate::align::{banded_nw_splice_aware3};
use crate::align2::{banded_nw_splice_aware3};
use std::time::{Duration, Instant};
use std::process;
use std::fs::File;
use std::io::Write;


pub struct MatElement {
    pub curr: u8,
    pub insert: Vec<u8>,
}


pub struct BaseMatrix {
    pub base_matrix: HashMap<String, Vec<u8>>,
    pub insertion_table: HashMap<usize, HashMap<String, String>>,
    pub expanded_matrix: HashMap<String, Vec<u8>>,
    pub bam_records: HashMap<String, bam::Record>,
    pub contig_name: String,
    pub start_position: u32,
    // 1-based
    pub end_position: u32,      // 1-based
}

impl BaseMatrix {
    pub fn new() -> BaseMatrix {
        BaseMatrix {
            base_matrix: HashMap::new(),
            insertion_table: HashMap::new(),
            expanded_matrix: HashMap::new(),
            bam_records: HashMap::new(),
            contig_name: String::new(),
            start_position: 0,
            end_position: 0,
        }
    }

    pub fn load_data(&mut self, bam_file: String, region: Region) {
        /*
        load reads in the given region and extend each read to the left and right
        */
        let mut bam_reader = bam::IndexedReader::from_path(bam_file).unwrap();
        let result = bam_reader.fetch((region.chr.as_str(), region.start, region.end + 1));   // fetch: left include, right exclude
        self.contig_name = region.chr.clone();
        if result.is_err() {
            panic!("fetch region failed...");
        }
        let mut record = Record::new();
        let mut max_length = 0;
        // let mut left_most_position: i64 = region.start as i64 -1;   // 0-based
        let mut left_most_position: i64 = i64::MAX;   // 0-based
        // get the left most position for all reads in this region
        while let Some(result) = bam_reader.read(&mut record) {
            if result.is_err() {
                panic!("BAM parsing failed...");
            }
            if record.is_unmapped() || record.is_secondary() || record.is_supplementary() {
                continue;
            }
            let pos = record.pos(); // 0-based
            if pos < left_most_position {
                left_most_position = pos;
            }
        }
        self.start_position = (left_most_position + 1) as u32;
        let result = bam_reader.fetch((region.chr.as_str(), region.start, region.end + 1));   // fetch: left include, right exclude
        if result.is_err() {
            panic!("fetch region failed...");
        }
        while let Some(result) = bam_reader.read(&mut record) {
            if result.is_err() {
                panic!("BAM parsing failed...");
            }
            if record.is_unmapped() || record.is_secondary() || record.is_supplementary() {
                continue;
            }
            let cigar = record.cigar();
            let seq = record.seq().as_bytes();
            let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
            let pos = record.pos(); // 0-based
            if self.bam_records.get(qname.as_str()).is_none() {
                self.bam_records.insert(qname.clone(), record.clone());
            }
            if self.base_matrix.get_mut(&qname).is_none() {
                self.base_matrix.insert(qname.clone(), Vec::new());
            }
            let mut pos_on_ref = pos;
            let mut pos_on_query = cigar.leading_softclips();
            let row = self.base_matrix.get_mut(&qname).unwrap();
            // left padding
            for _ in 0..pos_on_ref - left_most_position {
                row.push(b' ');
                // row.push(MatElement {
                //     curr: b' ',
                //     insert: Vec::new(),
                // });
            }
            for cg in cigar.iter() {
                match cg.char() as u8 {
                    b'S' | b'H' => {
                        continue;
                    }
                    b'M' => {
                        for _ in 0..cg.len() {
                            let qbase = seq[pos_on_query as usize];
                            // row.push(MatElement {
                            //     curr: qbase,
                            //     insert: Vec::new(),
                            // });
                            row.push(qbase);
                            pos_on_query += 1;
                            pos_on_ref += 1;
                        }
                    }
                    b'I' => {
                        let insert_seq = std::str::from_utf8(seq[pos_on_query as usize..(pos_on_query + cg.len() as i64) as usize].to_vec().as_slice()).unwrap().to_string();
                        if self.insertion_table.get(&(row.len() - 1)).is_none() {
                            self.insertion_table.insert(row.len().clone() - 1, HashMap::new());
                        }
                        self.insertion_table.get_mut(&(row.len() - 1)).unwrap().insert(qname.clone(), insert_seq);
                        // row.last_mut().unwrap().insert.extend(insert_seq);
                        pos_on_query += cg.len() as i64;
                    }
                    b'D' => {
                        for _ in 0..cg.len() {
                            // row.push(MatElement {
                            //     curr: b'-',
                            //     insert: Vec::new(),
                            // });
                            row.push(b'-');
                            pos_on_ref += 1;
                        }
                    }
                    b'N' => {
                        for _ in 0..cg.len() {
                            // row.push(MatElement {
                            //     curr: b'N',
                            //     insert: Vec::new(),
                            // });
                            row.push(b'N');
                            pos_on_ref += 1;
                        }
                    }
                    _ => {
                        panic!("Invalid cigar char: {}", cg.char());
                    }
                }
            }
            if pos_on_ref - left_most_position as i64 > max_length {
                max_length = pos_on_ref - left_most_position as i64;
            }
        }
        self.end_position = self.start_position + max_length as u32 - 1;
        // right padding
        for (_, row) in self.base_matrix.iter_mut() {
            while row.len() < max_length as usize {
                // row.push(MatElement {
                //     curr: b' ',
                //     insert: Vec::new(),
                // });
                row.push(b' ');
            }
        }
    }

    pub fn load_ref_data(&mut self, ref_seqs: HashMap<String, Vec<u8>>) {
        let contig_seq = ref_seqs.get(&self.contig_name).unwrap();
        let start_pos = self.start_position - 1;    // 0-based
        let end_pos = self.end_position - 1;    // 0-based, include
        let region_seq = &contig_seq[start_pos as usize..=end_pos as usize];
        let mut row: Vec<u8> = Vec::new();
        for nucle in region_seq.iter() {
            row.push(*nucle);
        }
        assert_eq!(row.len() as u32, self.end_position - self.start_position + 1);
        self.base_matrix.insert("ref".to_string(), row);
    }

    pub fn expand_insertion(&mut self) {
        let orig_matrix_len = self.end_position - self.start_position + 1;
        let mut read_names: Vec<String> = Vec::new();
        for qname in self.base_matrix.keys() {
            read_names.push(qname.clone());
            self.expanded_matrix.insert(qname.clone(), Vec::new());
        }
        for i in 0..orig_matrix_len {
            let mut max_insertion_len = 0;
            // fill current base to the expanded matrix and get the max insertion length
            for qname in read_names.iter() {
                let row = self.base_matrix.get(qname).unwrap();
                // self.expanded_matrix.get_mut(qname).unwrap().push(row[i as usize].curr);
                self.expanded_matrix.get_mut(qname).unwrap().push(row[i as usize]);
                if self.insertion_table.get(&(i as usize)).is_some() && self.insertion_table.get(&(i as usize)).unwrap().get(qname).is_some() {
                    let insert_seq = self.insertion_table.get(&(i as usize)).unwrap().get(qname).unwrap();
                    if insert_seq.len() > max_insertion_len {
                        max_insertion_len = insert_seq.len();
                    }
                }
                // if row[i as usize].insert.len() > max_insertion_len {
                //     max_insertion_len = row[i as usize].insert.len();
                // }
            }
            // fill insertion to the expanded matrix
            if max_insertion_len == 0 {
                continue;
            }
            for qname in read_names.iter() {
                let row = self.base_matrix.get(qname).unwrap();
                let mut insert_size = 0;
                if self.insertion_table.get(&(i as usize)).unwrap().get(qname).is_some() {
                    let insert_seq = self.insertion_table.get(&(i as usize)).unwrap().get(qname).unwrap();
                    for nucle in insert_seq.as_bytes() {
                        self.expanded_matrix.get_mut(qname).unwrap().push(*nucle);
                        insert_size += 1;
                    }
                }
                // for nucle in row[i as usize].insert.iter() {
                //     self.expanded_matrix.get_mut(qname).unwrap().push(*nucle);
                // }

                if insert_size == max_insertion_len {
                    continue;
                } else {
                    if i as i32 - 1 >= 0 && row[i as usize - 1] == b' ' {
                        for _ in 0..max_insertion_len - insert_size {
                            self.expanded_matrix.get_mut(qname).unwrap().push(b' ');
                        }
                        // continue;
                    } else if i + 1 < orig_matrix_len && row[i as usize + 1] == b' ' {
                        for _ in 0..max_insertion_len - insert_size {
                            self.expanded_matrix.get_mut(qname).unwrap().push(b' ');
                        }
                        // continue;
                    } else {
                        for _ in 0..max_insertion_len - insert_size {
                            self.expanded_matrix.get_mut(qname).unwrap().push(b'-');
                        }
                    }
                }
                // if row[i as usize].insert.len() == max_insertion_len {
                //     continue;
                // } else {
                //     if i as i32 - 1 >= 0 && row[i as usize - 1].curr == b' ' {
                //         for _ in 0..max_insertion_len - row[i as usize].insert.len() {
                //             self.expanded_matrix.get_mut(qname).unwrap().push(b' ');
                //         }
                //         // continue;
                //     } else if i + 1 < orig_matrix_len && row[i as usize + 1].curr == b' ' {
                //         for _ in 0..max_insertion_len - row[i as usize].insert.len() {
                //             self.expanded_matrix.get_mut(qname).unwrap().push(b' ');
                //         }
                //         // continue;
                //     } else {
                //         for _ in 0..max_insertion_len - row[i as usize].insert.len() {
                //             self.expanded_matrix.get_mut(qname).unwrap().push(b'-');
                //         }
                //     }
                // }
            }
        }
    }

    pub fn load_data_without_extension(&mut self, bam_file: String, region: Region) {
        /*
        load reads which are entirely in the given region (without extension)
        */
        let mut bam_reader = bam::IndexedReader::from_path(bam_file).unwrap();
        let result = bam_reader.fetch((region.chr.as_str(), region.start, region.end + 1)); // left include, right exclude
        self.contig_name = region.chr.clone();
        if result.is_err() {
            panic!("fetch region failed...");
        }
        let mut record = Record::new();
        self.start_position = region.start;
        self.end_position = region.end;
        let max_length = region.end - region.start + 1;
        while let Some(result) = bam_reader.read(&mut record) {
            if result.is_err() {
                panic!("BAM parsing failed...");
            }
            if record.is_unmapped() || record.is_secondary() || record.is_supplementary() {
                continue;
            }
            let cigar = record.cigar();
            let seq = record.seq().as_bytes();
            let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
            let pos = record.pos(); // 0-based
            if pos < region.start as i64 - 1 || pos > region.end as i64 - 1 {
                // left side of the read is not in the region
                continue;
            }
            if record.reference_end() > region.end as i64 - 1 {
                // right side of the read is not in the region
                continue;
            }

            if self.bam_records.get(qname.as_str()).is_none() {
                self.bam_records.insert(qname.clone(), record.clone());
            }
            if self.base_matrix.get_mut(&qname).is_none() {
                self.base_matrix.insert(qname.clone(), Vec::new());
            }
            let mut pos_on_ref = pos;
            let mut pos_on_query = cigar.leading_softclips();
            let row = self.base_matrix.get_mut(&qname).unwrap();
            // left padding
            for _ in 0..pos_on_ref - region.start as i64 + 1 {
                row.push(b' ');
                // row.push(MatElement {
                //     curr: b' ',
                //     insert: Vec::new(),
                // });
            }
            for cg in cigar.iter() {
                match cg.char() as u8 {
                    b'S' | b'H' => {
                        continue;
                    }
                    b'M' => {
                        for _ in 0..cg.len() {
                            let qbase = seq[pos_on_query as usize];
                            row.push(qbase);
                            // row.push(MatElement {
                            //     curr: qbase,
                            //     insert: Vec::new(),
                            // });
                            pos_on_query += 1;
                            pos_on_ref += 1;
                        }
                    }
                    b'I' => {
                        // let insert_seq = seq[pos_on_query as usize..(pos_on_query + cg.len() as i64) as usize].to_vec();
                        // row.last_mut().unwrap().insert.extend(insert_seq);
                        // pos_on_query += cg.len() as i64;
                        let insert_seq = std::str::from_utf8(seq[pos_on_query as usize..(pos_on_query + cg.len() as i64) as usize].to_vec().as_slice()).unwrap().to_string();
                        if self.insertion_table.get(&row.len()).is_none() {
                            self.insertion_table.insert(row.len().clone(), HashMap::new());
                        }
                        self.insertion_table.get_mut(&row.len()).unwrap().insert(qname.clone(), insert_seq);
                        pos_on_query += cg.len() as i64;
                    }
                    b'D' => {
                        for _ in 0..cg.len() {
                            // row.push(MatElement {
                            //     curr: b'-',
                            //     insert: Vec::new(),
                            // });
                            row.push(b'-');
                            pos_on_ref += 1;
                        }
                    }
                    b'N' => {
                        for _ in 0..cg.len() {
                            // row.push(MatElement {
                            //     curr: b'N',
                            //     insert: Vec::new(),
                            // });
                            row.push(b'N');
                            pos_on_ref += 1;
                        }
                    }
                    _ => {
                        panic!("Invalid cigar char: {}", cg.char());
                    }
                }
            }
        }
        // right padding
        for (_, row) in self.base_matrix.iter_mut() {
            while row.len() < max_length as usize {
                // row.push(MatElement {
                //     curr: b' ',
                //     insert: Vec::new(),
                // });
                row.push(b' ');
            }
        }
    }
}

pub fn load_reference(ref_path: String) -> HashMap<String, Vec<u8>> {
    let mut ref_seqs: HashMap<String, Vec<u8>> = HashMap::new();
    let reader = fasta::Reader::from_file(ref_path).unwrap();
    for r in reader.records() {
        let ref_record = r.unwrap();
        ref_seqs.insert(ref_record.id().to_string(), ref_record.seq().to_vec());
    }
    return ref_seqs;
}


pub fn generate_column_profile(expanded_matrix: &HashMap<String, Vec<u8>>, column_base_counts: &mut Vec<ColumnBaseCount>) {
    assert!(expanded_matrix.len() > 0);
    let ncols = expanded_matrix.iter().next().unwrap().1.len();
    let mut ref_base: u8 = 0;
    let mut column_bases: Vec<u8> = Vec::new();
    for i in 0..ncols {
        for (readname, base_vec) in expanded_matrix.iter() {
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

pub fn generate_reduced_profile(
    expanded_matrix: &HashMap<String, Vec<u8>>,
    forward_donor_penalty: &Vec<f64>,
    forward_acceptor_penalty: &Vec<f64>,
    reverse_donor_penalty: &Vec<f64>,
    reverse_acceptor_penalty: &Vec<f64>,
    column_base_counts: &mut Vec<ColumnBaseCount>,
    column_indexes: &mut Vec<usize>,
    reduced_expanded_matrix: &mut HashMap<String, Vec<u8>>,
    forward_reduced_donor_penalty: &mut Vec<f64>,
    forward_reduced_acceptor_penalty: &mut Vec<f64>,
    reverse_reduced_donor_penalty: &mut Vec<f64>,
    reverse_reduced_acceptor_penalty: &mut Vec<f64>,
    splice_boundary: &mut Vec<bool>,
) {
    assert!(expanded_matrix.len() > 0);
    let ncols = expanded_matrix.iter().next().unwrap().1.len();
    assert!(forward_donor_penalty.len() == ncols + 1);
    assert!(forward_acceptor_penalty.len() == ncols + 1);
    assert!(reverse_donor_penalty.len() == ncols + 1);
    assert!(reverse_acceptor_penalty.len() == ncols + 1);
    let mut ref_base: u8 = 0;
    let mut column_bases: Vec<u8> = Vec::new();
    // TODO: if the splicing cut is short than extend size?
    let mut extend_size = 96; // L>(d(i)+a(a)+telda(q)-q)/e
    for i in 0..ncols {
        for (readname, base_vec) in expanded_matrix.iter() {
            if *readname == "ref".to_string() {
                ref_base = base_vec[i];
                continue;
            }
            column_bases.push(base_vec[i]);
        }
        let cbc = ColumnBaseCount::new_from_column(&column_bases, ref_base);
        if (cbc.n_a + cbc.n_c + cbc.n_g + cbc.n_t + cbc.n_dash) == 0 && cbc.n_n > 0 {
            if extend_size > 0 {
                // extend 32 base at the beginning of the reduced region
                column_base_counts.push(cbc);
                column_indexes.push(i);
                column_bases.clear();
                forward_reduced_donor_penalty.push(forward_donor_penalty[i]);
                forward_reduced_acceptor_penalty.push(forward_acceptor_penalty[i]);
                reverse_reduced_donor_penalty.push(reverse_donor_penalty[i]);
                reverse_reduced_acceptor_penalty.push(reverse_acceptor_penalty[i]);
                extend_size -= 1;
            }
            column_bases.clear();
            continue;
        }
        column_base_counts.push(cbc);
        column_indexes.push(i);
        column_bases.clear();
        forward_reduced_donor_penalty.push(forward_donor_penalty[i]);
        forward_reduced_acceptor_penalty.push(forward_acceptor_penalty[i]);
        reverse_reduced_donor_penalty.push(reverse_donor_penalty[i]);
        reverse_reduced_acceptor_penalty.push(reverse_acceptor_penalty[i]);
        extend_size = 96;
    }
    // trick: donor[i] store the penalty of ref[i], acceptor[i] store the penalty of ref[i-1].
    // The size of donor and acceptor is ref_base_vec.len() + 1.
    // The last element of donor is meaningless, but the last element of acceptor is meaningful.
    forward_reduced_acceptor_penalty.push(forward_acceptor_penalty[ncols]);
    reverse_reduced_acceptor_penalty.push(reverse_acceptor_penalty[ncols]);
    reduced_expanded_matrix.clear();
    for i in column_indexes.iter() {
        for (readname, base_vec) in expanded_matrix.iter() {
            if reduced_expanded_matrix.get(readname).is_none() {
                reduced_expanded_matrix.insert(readname.clone(), vec![base_vec[*i]]);
            } else {
                reduced_expanded_matrix.get_mut(readname).unwrap().push(base_vec[*i]);
            }
        }
    }

    for i in 0..column_indexes.len() - 1 {
        if column_indexes[i + 1] - column_indexes[i] != 1 {
            splice_boundary.push(true);
        } else {
            splice_boundary.push(false);
        }
    }
    splice_boundary.push(false);
}

pub fn get_donor_acceptor_penalty(
    expanded_matrix: &HashMap<String, Vec<u8>>,
    standed_penalty: f64,
) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let mut forward_donor_penalty: Vec<f64> = Vec::new();
    let mut forward_acceptor_penalty: Vec<f64> = Vec::new();
    let mut reverse_donor_penalty: Vec<f64> = Vec::new();
    let mut reverse_acceptor_penalty: Vec<f64> = Vec::new();
    let mut ref_sliding_window: Vec<u8> = Vec::new();
    let ref_base_vec = expanded_matrix.get("ref").unwrap();
    for i in 0..ref_base_vec.len() + 1 {
        // trick: donor[i] store the penalty of ref[i], acceptor[i] store the penalty of ref[i-1].
        // The size of donor and acceptor is ref_base_vec.len() + 1.
        // The last element of donor is meaningless, but the last element of acceptor is meaningful.
        if i < ref_base_vec.len() && ref_base_vec[i] == b'-' {
            forward_donor_penalty.push(standed_penalty);
            forward_acceptor_penalty.push(standed_penalty);
            reverse_donor_penalty.push(standed_penalty);
            reverse_acceptor_penalty.push(standed_penalty);
            continue;
        }

        // trick: donor[i] store the penalty of ref[i], acceptor[i] store the penalty of ref[i-1].
        if i as i32 - 3 < 0 {
            forward_acceptor_penalty.push(standed_penalty);
            reverse_acceptor_penalty.push(standed_penalty);
        } else {
            let mut j = i as i32 - 1;
            ref_sliding_window.clear();
            while ref_sliding_window.len() < 3 && j >= 0 {
                if ref_base_vec[j as usize] != b'-' {
                    ref_sliding_window.push(ref_base_vec[j as usize]);
                }
                j -= 1;
            }
            let mut tstr = String::new();
            for c in ref_sliding_window.iter().rev() {
                tstr.push(*c as char);
            }

            if tstr.len() < 3 {
                forward_acceptor_penalty.push(standed_penalty);
                reverse_acceptor_penalty.push(standed_penalty);
            } else {
                if tstr[1..3] == "AC".to_string() {
                    forward_acceptor_penalty.push(standed_penalty - 9.0);
                } else if tstr == "AAG".to_string() || tstr == "GAG".to_string() {
                    forward_acceptor_penalty.push(standed_penalty - 22.0);
                } else if tstr == "CAG".to_string() || tstr == "TAG".to_string() {
                    forward_acceptor_penalty.push(0.0);
                } else {
                    forward_acceptor_penalty.push(standed_penalty);
                }

                if tstr[1..3] == "AT".to_string() {
                    reverse_acceptor_penalty.push(standed_penalty - 9.0);
                } else if tstr[1..3] == "GC".to_string() {
                    reverse_acceptor_penalty.push(standed_penalty - 15.0);
                } else if tstr == "GAC".to_string() || tstr == "AAC".to_string() {
                    reverse_acceptor_penalty.push(standed_penalty - 22.0);
                } else if tstr == "TAC".to_string() || tstr == "CAC".to_string() {
                    reverse_acceptor_penalty.push(0.0);
                } else {
                    reverse_acceptor_penalty.push(standed_penalty);
                }
            }
        }

        if i + 2 >= ref_base_vec.len() {
            forward_donor_penalty.push(standed_penalty);
            reverse_donor_penalty.push(standed_penalty);
        } else {
            let mut j = i;
            ref_sliding_window.clear();
            while ref_sliding_window.len() < 3 && j < ref_base_vec.len() {
                if ref_base_vec[j] != b'-' {
                    ref_sliding_window.push(ref_base_vec[j]);
                }
                j += 1;
            }
            let mut tstr = String::new();
            for c in ref_sliding_window.iter() {
                tstr.push(*c as char);
            }
            if tstr.len() < 3 {
                forward_donor_penalty.push(standed_penalty);
                reverse_donor_penalty.push(standed_penalty);
            } else {
                if tstr[0..2] == "AT".to_string() {
                    forward_donor_penalty.push(standed_penalty - 9.0);
                } else if tstr[0..2] == "GC".to_string() {
                    forward_donor_penalty.push(standed_penalty - 15.0);
                } else if tstr == "GTC".to_string() || tstr == "GTT".to_string() {
                    forward_donor_penalty.push(standed_penalty - 22.0);
                } else if tstr == "GTA".to_string() || tstr == "GTG".to_string() {
                    forward_donor_penalty.push(0.0);
                } else {
                    forward_donor_penalty.push(standed_penalty);
                }

                if tstr[0..2] == "GT".to_string() {
                    reverse_donor_penalty.push(standed_penalty - 9.0);
                } else if tstr == "CTT".to_string() || tstr == "CTC".to_string() {
                    reverse_donor_penalty.push(standed_penalty - 22.0);
                } else if tstr == "CTG".to_string() || tstr == "CTA".to_string() {
                    reverse_donor_penalty.push(0.0);
                } else {
                    reverse_donor_penalty.push(standed_penalty);
                }
            }
        }
    }
    (forward_donor_penalty, forward_acceptor_penalty, reverse_donor_penalty, reverse_acceptor_penalty)
}

pub fn profile_realign(
    expanded_matrix: &HashMap<String, Vec<u8>>,
    forward_donor_penalty: &Vec<f64>,
    forward_acceptor_penalty: &Vec<f64>,
    reverse_donor_penalty: &Vec<f64>,
    reverse_acceptor_penalty: &Vec<f64>,
    best_reduced_expanded_matrix: &mut HashMap<String, Vec<u8>>,
    best_column_indexes: &mut Vec<usize>,
) {
    let mut old_score = -f64::INFINITY;
    let mut new_score = 0.0;
    let mut profile: Vec<ColumnBaseCount> = Vec::new();
    let mut column_indexes: Vec<usize> = Vec::new();
    let mut reduced_expanded_matrix: HashMap<String, Vec<u8>> = HashMap::new();
    let mut forward_reduced_donor_penalty: Vec<f64> = Vec::new();
    let mut forward_reduced_acceptor_penalty: Vec<f64> = Vec::new();
    let mut reverse_reduced_donor_penalty: Vec<f64> = Vec::new();
    let mut reverse_reduced_acceptor_penalty: Vec<f64> = Vec::new();
    // let mut forward_hidden_splice_penalty: Vec<f64> = Vec::new();
    // let mut reverse_hidden_splice_penalty: Vec<f64> = Vec::new();
    let mut splice_boundary: Vec<bool> = Vec::new();
    let mut prev_aligned_seq: Vec<u8> = Vec::new();
    generate_reduced_profile(
        expanded_matrix,
        forward_donor_penalty,
        forward_acceptor_penalty,
        reverse_donor_penalty,
        reverse_acceptor_penalty,
        &mut profile,
        &mut column_indexes,
        &mut reduced_expanded_matrix,
        &mut forward_reduced_donor_penalty,
        &mut forward_reduced_acceptor_penalty,
        &mut reverse_reduced_donor_penalty,
        &mut reverse_reduced_acceptor_penalty,
        &mut splice_boundary,
    );
    best_column_indexes.clear();
    *best_column_indexes = column_indexes.clone();
    for i in 0..profile.len() {
        prev_aligned_seq.push(profile[i].get_major_base());
    }
    // println!("major sequence:");
    // for i in 0..profile.len() {
    //     print!("{}", profile[i].get_major_base() as char);
    // }
    // println!();
    // println!("ref sequence:");
    // for i in 0..profile.len() {
    //     print!("{}", profile[i].get_ref_base() as char);
    // }
    // println!();
    for i in 0..profile.len() {
        let cbc = &profile[i];
        new_score += cbc.get_score(&cbc.get_major_base());
    }
    let mut iteration = 0;
    let mut cnt = 0;
    while new_score > old_score {
        let realign_s = Instant::now();
        let mut double_strand_align_runtime: u128 = 0;
        // let mut insert_runtime: u128 = 0;
        let mut generate_runtime: u128 = 0;
        iteration += 1;
        println!("Iteration: {}, old_score: {}, new_score: {}", iteration, old_score, new_score);
        let mut sorted_expanded_matrix: Vec<(&String, &Vec<u8>)> = expanded_matrix.iter().collect();
        sorted_expanded_matrix.sort_by(|a, b| a.0.cmp(b.0));
        // for (readname, base_vec) in expanded_matrix.iter() {
        for (readname, base_vec) in sorted_expanded_matrix {
            if *readname == "ref".to_string() {
                continue;
            }

            cnt += 1;
            // let query = String::from_utf8(base_vec.clone()).unwrap().replace(" ", "").replace("-", "").replace("N", "");
            let query =
                std::str::from_utf8(reduced_expanded_matrix.get(readname).unwrap()).unwrap().as_bytes();

            let mut left_most_query_pos = 0;    // not include the padding
            let mut right_most_query_pos = 0;   // not include the padding
            for i in 0..query.len() {
                if query[i] != b' ' {
                    left_most_query_pos = i;
                    break;
                }
            }

            if left_most_query_pos as i32 - 11 >= 0 {
                left_most_query_pos -= 11;
            } else {
                left_most_query_pos = 0;
            }

            for i in (0..query.len()).rev() {
                if query[i] != b' ' {
                    right_most_query_pos = i;
                    break;
                }
            }

            if right_most_query_pos + 12 <= query.len() - 1 {
                right_most_query_pos += 12;
            } else {
                right_most_query_pos = query.len() - 1;
            }

            // TODO: 1. get the banded start position and banded end position on the reference coordinate and record the related column index.
            //       2. choose the sub_matrix from base_matrix between the banded start position and banded end position
            //       3. calculate the alignment (reduce the size of target by ignoring the front padding and back padding)
            //       4. update the base_matrix according the realignment and recorded related column index.

            // TODO: update reduced_expanded_matrix after realigning a group of reads instead of realigning each read, reduce the times of generate new profile.
            // println!("query: \n{}", query);
            // println!("query len:{}, profile len:{}", query.len(), profile.len());
            // println!("align begin, profile length: {}", &profile.len());
            // let (alignment_score, aligned_query, ref_target, major_target) = nw_splice_aware(&query.as_bytes().to_vec(), &profile);
            // let (alignment_score, aligned_query, ref_target, major_target) = semi_nw_splice_aware(&query.as_bytes().to_vec(), &profile);
            // let (alignment_score, aligned_query, ref_target, major_target) = banded_nw_splice_aware(&query.as_bytes().to_vec(), &profile, 20);
            // let (alignment_score, aligned_query, ref_target, major_target) = banded_nw_splice_aware2(&query.as_bytes().to_vec(), &profile, 20);
            let a_s = Instant::now();
            println!("qname: {}", readname);
            println!("query size: {}, left_most_query_pos: {}, right_most_query_pos: {}", query.len(), left_most_query_pos, right_most_query_pos);
            let (
                reverse_alignment_score,
                reverse_aligned_query,
                reverse_ref_target,
                reverse_major_target,
            ) = banded_nw_splice_aware3(
                &query[left_most_query_pos..=right_most_query_pos].to_vec(),
                &profile[left_most_query_pos..=right_most_query_pos].to_vec(),
                &reverse_reduced_donor_penalty[left_most_query_pos..=right_most_query_pos].to_vec(),
                &reverse_reduced_acceptor_penalty[left_most_query_pos..=right_most_query_pos + 1].to_vec(),
                &splice_boundary[left_most_query_pos..=right_most_query_pos].to_vec(),
                10,
            );
            let (
                forward_alignment_score,
                forward_aligned_query,
                forward_ref_target,
                forward_major_target,
            ) = banded_nw_splice_aware3(
                &query[left_most_query_pos..=right_most_query_pos].to_vec(),
                &profile[left_most_query_pos..=right_most_query_pos].to_vec(),
                &forward_reduced_donor_penalty[left_most_query_pos..=right_most_query_pos].to_vec(),
                &forward_reduced_acceptor_penalty[left_most_query_pos..=right_most_query_pos + 1].to_vec(),
                &splice_boundary[left_most_query_pos..=right_most_query_pos].to_vec(),
                10,
            );
            let a_d = a_s.elapsed().as_millis();
            double_strand_align_runtime += a_d;

            let alignment_score: f64;
            let aligned_query: Vec<u8>;
            let ref_target: Vec<u8>;
            let major_target: Vec<u8>;
            if forward_alignment_score >= reverse_alignment_score {
                alignment_score = forward_alignment_score;
                aligned_query = forward_aligned_query;
                ref_target = forward_ref_target;
                major_target = forward_major_target;
            } else {
                alignment_score = reverse_alignment_score;
                aligned_query = reverse_aligned_query;
                ref_target = reverse_ref_target;
                major_target = reverse_major_target;
            }

            // println!("align end");
            // println!("iter: {}, qname: {}", iteration, readname);
            // println!("original query: \n{}", query);
            // println!("ref target: \n{}", std::str::from_utf8(&ref_target).unwrap());
            // println!("major target: \n{}", std::str::from_utf8(&major_target).unwrap());
            // println!("aligned query: \n{}", std::str::from_utf8(&aligned_query).unwrap());
            // println!("alignment score: {}", alignment_score);
            assert!(aligned_query.len() + left_most_query_pos + query.len() - 1 - right_most_query_pos == reduced_expanded_matrix.get(readname).unwrap().len());
            // let insert_s = Instant::now();
            // reduced_expanded_matrix.insert(readname.clone(), aligned_query);
            // insert_runtime += insert_s.elapsed().as_millis();
            // profile.clear();
            let generate_s = Instant::now();
            // PileupMatrix::generate_column_profile(&reduced_expanded_matrix, &mut profile);
            update_expanded_matrix_profile(&mut reduced_expanded_matrix, &mut profile, readname.clone(), &aligned_query, left_most_query_pos, right_most_query_pos);
            generate_runtime += generate_s.elapsed().as_millis();
        }
        old_score = new_score;
        new_score = 0.0;
        for i in 0..profile.len() {
            let cbc = &profile[i];
            new_score += cbc.get_score(&cbc.get_major_base());
        }
        if new_score > old_score {
            prev_aligned_seq.clear();
            for i in 0..profile.len() {
                prev_aligned_seq.push(profile[i].get_major_base());
            }
            best_reduced_expanded_matrix.clear();
            *best_reduced_expanded_matrix = reduced_expanded_matrix.clone();
        }
        let realign_duration = realign_s.elapsed().as_millis();
        println!("realign duration: {}", realign_duration);
        println!("double strand align runtime: {}", double_strand_align_runtime);
        // println!("insert runtime: {}", insert_runtime);
        println!("generate runtime: {}", generate_runtime);
    }
    println!("total cnt: {}", cnt);
    // println!("new major sequence:");
    // println!("{}", String::from_utf8(prev_aligned_seq).unwrap());
}

pub fn update_expanded_matrix_profile(
    expanded_matrix: &mut HashMap<String, Vec<u8>>,
    column_base_counts: &mut Vec<ColumnBaseCount>,
    readname: String, readseq: &Vec<u8>,
    left_most_pos: usize,
    right_most_pos: usize,
) {
    let prev_readseq = expanded_matrix.get_mut(&readname).unwrap();
    let prev_readseq_len = prev_readseq.len();
    assert!(prev_readseq_len == readseq.len() + left_most_pos + prev_readseq_len - 1 - right_most_pos);
    assert!(prev_readseq_len == column_base_counts.len());
    for i in 0..readseq.len() {
        let prev_base = prev_readseq[left_most_pos + i];
        let new_base = readseq[i];
        if prev_base == new_base {
            continue;
        } else {
            if prev_base != b'N' {
                match prev_base {
                    b'A' => {
                        column_base_counts[left_most_pos + i].n_a -= 1;
                    }
                    b'C' => {
                        column_base_counts[left_most_pos + i].n_c -= 1;
                    }
                    b'G' => {
                        column_base_counts[left_most_pos + i].n_g -= 1;
                    }
                    b'T' => {
                        column_base_counts[left_most_pos + i].n_t -= 1;
                    }
                    b' ' => {
                        column_base_counts[left_most_pos + i].n_blank -= 1;
                    }
                    b'-' => {
                        column_base_counts[left_most_pos + i].n_dash -= 1;
                    }
                    _ => {
                        panic!("Invalid base: {}", prev_base);
                    }
                }
            }
            match new_base {
                b'A' => {
                    column_base_counts[left_most_pos + i].n_a += 1;
                    prev_readseq[left_most_pos + i] = b'A';
                }
                b'C' => {
                    column_base_counts[left_most_pos + i].n_c += 1;
                    prev_readseq[left_most_pos + i] = b'C';
                }
                b'G' => {
                    column_base_counts[left_most_pos + i].n_g += 1;
                    prev_readseq[left_most_pos + i] = b'G';
                }
                b'T' => {
                    column_base_counts[left_most_pos + i].n_t += 1;
                    prev_readseq[left_most_pos + i] = b'T';
                }
                b' ' => {
                    column_base_counts[left_most_pos + i].n_blank += 1;
                    prev_readseq[left_most_pos + i] = b' ';
                }
                b'-' => {
                    column_base_counts[left_most_pos + i].n_dash += 1;
                    prev_readseq[left_most_pos + i] = b'-';
                }
                b'N' => {
                    prev_readseq[left_most_pos + i] = b'N';
                }
                _ => {
                    panic!("Invalid base: {}", new_base);
                }
            }
            column_base_counts[left_most_pos + i].update_max_count();
        }
    }
}

pub fn update_expanded_matrix_from_realign(
    expanded_matrix: &mut HashMap<String, Vec<u8>>,
    realign_base_matrix: &HashMap<String, Vec<u8>>,
    column_indexes: &Vec<usize>,
) {
    for (readname, base_vec) in realign_base_matrix.iter() {
        if *readname == "ref".to_string() {
            continue;
        }
        // println!("readname: {}, base_vec length: {}, column_indexes length: {}", readname, base_vec.len(), column_indexes.len());
        assert!(base_vec.len() == column_indexes.len());
        for i in 0..column_indexes.len() {
            if expanded_matrix.get(readname).unwrap()[column_indexes[i]] != base_vec[i] {
                if expanded_matrix.get(readname).unwrap()[column_indexes[i]] as char != ' '
                    && base_vec[i] as char != 'N'
                {
                    println!("Modified: readname: {}, column index: {}, old base: {}, new base: {}",
                             readname,
                             column_indexes[i],
                             expanded_matrix.get(readname).unwrap()[column_indexes[i]] as char,
                             base_vec[i] as char
                    );
                }
            }
            expanded_matrix.get_mut(readname).unwrap()[column_indexes[i]] = base_vec[i];
        }
    }
}

pub fn update_bam_records_from_realign(
    expanded_matrix: &mut HashMap<String, Vec<u8>>,
    bam_records: &mut HashMap<String, bam::Record>,
    start_position: u32,
    end_position: u32,
) {
    let ref_seq = expanded_matrix.get("ref").unwrap().clone();
    for (readname, base_vec) in expanded_matrix.iter_mut() {
        let mut new_cigar: Vec<Cigar> = Vec::new();
        if *readname == "ref".to_string() {
            continue;
        }

        // get the left clip and right clip, the updated cigar will have the same left clip and right clip
        let mut left_soft_clip = 0;
        let mut right_soft_clip = 0;
        let mut left_hard_clip = 0;
        let mut right_hard_clip = 0;

        let cg = bam_records
            .get(readname)
            .unwrap()
            .cigar()
            .iter()
            .next()
            .unwrap()
            .clone();
        if cg.char() == 'S' {
            left_soft_clip = cg.len() as u32;
            new_cigar.push(Cigar::SoftClip(cg.len() as u32));
        } else if cg.char() == 'H' {
            left_hard_clip = cg.len() as u32;
            new_cigar.push(Cigar::HardClip(cg.len() as u32));
        }

        let cg = bam_records
            .get(readname)
            .unwrap()
            .cigar()
            .iter()
            .last()
            .unwrap()
            .clone();
        if cg.char() == 'S' {
            right_soft_clip = cg.len() as u32;
        } else if cg.char() == 'H' {
            right_hard_clip = cg.len() as u32;
        }

        // update the cigar
        let mut first_base_pair_index = 0;
        let mut last_base_pair_index = 0;
        for i in 0..base_vec.len() {
            if base_vec[i] != b' ' && base_vec[i] != b'N' && base_vec[i] != b'-' {
                first_base_pair_index = i;
                break;
            }
        }
        for i in (0..base_vec.len()).rev() {
            if base_vec[i] != b' ' && base_vec[i] != b'N' && base_vec[i] != b'-' {
                last_base_pair_index = i;
                break;
            }
        }
        // begin:       NNNNNNNACGT
        // begin:--NNNNNNNNNNNNACGT
        for i in 0..first_base_pair_index {
            if base_vec[i] == b'N' || base_vec[i] == b'-' {
                base_vec[i] = b' ';
            }
        }
        // ACGTNNNNNNNNN          :end
        // ACGTNNNNNNNNN----------:end
        for i in last_base_pair_index..base_vec.len() {
            if base_vec[i] == b'N' || base_vec[i] == b'-' {
                base_vec[i] = b' ';
            }
        }

        //  ACTGNNNNN         GTACNN         NNNNNNNN
        for i in first_base_pair_index..=last_base_pair_index {
            if base_vec[i] == b' ' {
                base_vec[i] = b'N';
            }
        }

        let mut prev_op: char = ' ';
        let mut prev_len = 0;
        let mut op: char = ' ';
        for i in first_base_pair_index..=last_base_pair_index {
            if prev_op == ' ' {
                if base_vec[i] != b'-' && ref_seq[i] != b'-' {
                    prev_op = 'M';
                    prev_len += 1;
                    continue;
                } else if base_vec[i] == b'-' && ref_seq[i] != b'-' {
                    prev_op = 'D';
                    prev_len += 1;
                    continue;
                } else if base_vec[i] != b'-' && ref_seq[i] == b'-' {
                    prev_op = 'I';
                    prev_len += 1;
                    continue;
                } else {
                    panic!("should not happen1");
                    process::exit(1);
                }
            }

            if base_vec[i] != b'-'
                && base_vec[i] != b'N'
                && ref_seq[i] != b'-'
                && ref_seq[i] != b'N'
            {
                op = 'M';
            } else if base_vec[i] == b'-' && ref_seq[i] != b'-' && ref_seq[i] != b'N' {
                op = 'D';
            } else if base_vec[i] != b'-' && base_vec[i] != b'N' && ref_seq[i] == b'-' {
                op = 'I';
            } else if base_vec[i] == b'N' && ref_seq[i] != b'N' && ref_seq[i] != b'-' {
                op = 'N';
            } else if base_vec[i] == b'N' && ref_seq[i] == b'-' {
                continue;
            } else if base_vec[i] == b'-' && ref_seq[i] == b'-' {
                continue;
            } else if ref_seq[i] == b'N' {
                // TODO: why reference can be N and query is dash or nucleotide at the same time?
                if base_vec[i] == b'N' {
                    op = 'N';
                } else if base_vec[i] == b'-' {
                    op = 'D';
                } else if base_vec[i] != b'N' && base_vec[i] != b'-' {
                    op = 'M';
                }
            } else {
                println!("readname: {}, base_vec[i]: {}, ref_seq[i]: {}", readname, base_vec[i] as char, ref_seq[i] as char);
                panic!("should not happen2");
                process::exit(1);
            }

            if op == prev_op {
                prev_len += 1;
            } else {
                if prev_op == 'M' {
                    new_cigar.push(Cigar::Match(prev_len));
                    prev_op = op;
                    prev_len = 1;
                } else if prev_op == 'I' {
                    new_cigar.push(Cigar::Ins(prev_len));
                    prev_op = op;
                    prev_len = 1;
                } else if prev_op == 'D' {
                    new_cigar.push(Cigar::Del(prev_len));
                    prev_op = op;
                    prev_len = 1;
                } else if prev_op == 'N' {
                    new_cigar.push(Cigar::RefSkip(prev_len));
                    prev_op = op;
                    prev_len = 1;
                } else {
                    panic!("should not happen3");
                    process::exit(1);
                }
            }
        }

        if prev_op == 'M' {
            new_cigar.push(Cigar::Match(prev_len));
        } else if prev_op == 'I' {
            new_cigar.push(Cigar::Ins(prev_len));
        } else if prev_op == 'D' {
            new_cigar.push(Cigar::Del(prev_len));
        } else if prev_op == 'N' {
            new_cigar.push(Cigar::RefSkip(prev_len));
        } else {
            panic!("should not happen4");
            process::exit(1);
        }

        if right_hard_clip > 0 {
            new_cigar.push(Cigar::HardClip(right_hard_clip));
        } else if right_soft_clip > 0 {
            new_cigar.push(Cigar::SoftClip(right_soft_clip));
        }

        let new_cigar_string = CigarString(new_cigar);
        let mut cglen = 0;
        for cg in new_cigar_string.iter() {
            // print!("{}{}", cg.len(), cg.char());
            if cg.char() == 'M' || cg.char() == 'I' || cg.char() == 'S' {
                cglen += cg.len();
            }
        }
        println!();
        println!("readname: {}, cigar len: {} read len: {}",
                 readname,
                 cglen,
                 bam_records.get(readname).unwrap().seq_len()
        );
        if cglen != bam_records.get(readname).unwrap().seq_len() as u32 {
            println!("cglen error: {}, cigar len: {}, read len: {}", readname, cglen, bam_records.get(readname).unwrap().seq_len());
            println!("{}", new_cigar_string.to_string());
            println!("{}", std::str::from_utf8(base_vec).unwrap());
        }
        assert!(cglen == bam_records.get(readname).unwrap().seq_len() as u32);

        // println!(
        //     "ref: \n{}",
        //     String::from_utf8(ref_seq.to_vec()).unwrap().clone()
        // );
        // println!("read: \n{}", String::from_utf8(base_vec.to_vec()).unwrap());

        // count how many blank bases in the beginning, this will add to the ref_start
        let mut blank_count = 0;
        let mut read_ref_start = start_position as i64 - 1; // start_position is 1-based, pos in bam is 0-based
        // print!("{},read ref start: {}, ", readname, read_ref_start);
        for i in 0..base_vec.len() {
            if base_vec[i] == b' ' {
                if ref_seq[i] != b' ' && ref_seq[i] != b'-' {
                    blank_count += 1;
                }
            } else {
                break;
            }
        }
        if blank_count > 0 {
            read_ref_start += blank_count as i64;
        }
        // print!("blank count: {}, read ref start: {}", blank_count, read_ref_start);
        // println!();
        let record = bam_records.get(readname).unwrap();
        let mut out_record = bam::Record::from(record.clone());
        out_record.set_pos(read_ref_start as i64);
        out_record.set(
            record.qname(),
            Some(&new_cigar_string),
            &record.seq().as_bytes(),
            record.qual(),
        );
        bam_records.insert(readname.clone(), out_record);
    }
}

pub fn write_bam_records(
    bam_records: &mut HashMap<String, bam::Record>,
    bam_file: &str,
    bam_header: &bam::header::Header,
) {
    let mut bam_writer = bam::Writer::from_path(bam_file, bam_header, Format::Bam).unwrap();
    for (_, record) in bam_records.iter() {
        let re = bam_writer.write(&record);
        if re != Ok(()) {
            println!("write failed");
            process::exit(1);
        }
    }
}

pub fn get_coverage_intervals(bam_path: String, depth_threshold: u32) -> (Vec<(String, i64, i64)>, Vec<(String, i64, i64)>) {
    let mut normal_depth_regions: Vec<(String, i64, i64)> = Vec::new();
    let mut high_depth_regions: Vec<(String, i64, i64)> = Vec::new();
    let mut ctg = String::new();
    let mut pre_ctg = String::new();
    let mut s: i64 = -1;
    let mut e: i64 = -1;
    let mut state = 0;  // 0 is normal depth, 1 is high depth
    let mut pre_state = 0;   // 0 is normal depth, 1 is high depth
    let mut pre_depth = 0;
    let mut pre_pos: u32 = 0;


    let mut bam = bam::IndexedReader::from_path(bam_path).unwrap();
    let header = bam.header().clone();
    for p in bam.pileup() {
        let pileup = p.unwrap();
        // let depth = pileup.alignments().filter(|alignment| {
        //     let record = alignment.record();
        //     !record.is_duplicate() && !record.is_secondary() && !record.is_supplementary() && !record.is_unmapped()
        // }).count() as u32;
        let depth = pileup.depth();
        let pos = pileup.pos(); //0-based
        // println!("pos: {}, depth: {}", pos+1, depth);

        if s == -1 && e == -1 {
            s = pos as i64;
            e = pos as i64;
            ctg = std::str::from_utf8(&header.tid2name(pileup.tid())).unwrap().to_string();
            pre_ctg = ctg.clone();
            if depth >= depth_threshold {
                state = 1;
            } else {
                state = 0;
            }
            pre_state = state;
            pre_depth = depth;
            pre_pos = pos;
            continue;
        }

        if depth >= depth_threshold {
            state = 1;
        } else {
            state = 0;
        }

        if state != pre_state || ctg != pre_ctg || pos != pre_pos + 1 {
            if e > s {
                // println!("{}, depth: {}, pre_depth: {}", pos+1, depth, pre_depth);
                if pre_state == 1 {
                    high_depth_regions.push((ctg.clone(), s + 1, e + 1));
                } else {
                    normal_depth_regions.push((ctg.clone(), s + 1, e + 1));
                }
            }
            ctg = std::str::from_utf8(&header.tid2name(pileup.tid())).unwrap().to_string();
            pre_ctg = ctg.clone();
            s = pos as i64;
            e = pos as i64;
            pre_state = state;
            pre_depth = depth;
            pre_pos = pos;
        } else {
            e = pos as i64;
            pre_depth = depth;
            pre_pos = pos;
        }
    }
    if pre_state == 1 && e > s {
        high_depth_regions.push((ctg.clone(), s + 1, e + 1));
    } else if pre_state == 0 && e > s {
        normal_depth_regions.push((ctg.clone(), s + 1, e + 1));
    }
    return (normal_depth_regions, high_depth_regions);
}


pub fn get_chrom_coverage_intervals(bam_path: String, ctgname: &str, depth_threshold: u32) -> (Vec<Region>, Vec<Region>) {
    let mut normal_depth_regions: Vec<Region> = Vec::new();
    let mut high_depth_regions: Vec<Region> = Vec::new();
    let mut ctg = String::new();
    let mut pre_ctg = String::new();
    let mut s: i64 = -1;
    let mut e: i64 = -1;
    let mut state = 0;  // 0 is normal depth, 1 is high depth
    let mut pre_state = 0;   // 0 is normal depth, 1 is high depth
    let mut pre_depth = 0;
    let mut pre_pos: u32 = 0;


    let mut bam = bam::IndexedReader::from_path(bam_path).unwrap();
    let header = bam.header().clone();
    bam.fetch(ctgname).unwrap();
    for p in bam.pileup() {
        let pileup = p.unwrap();
        // let depth = pileup.alignments().filter(|alignment| {
        //     let record = alignment.record();
        //     !record.is_duplicate() && !record.is_secondary() && !record.is_supplementary() && !record.is_unmapped()
        // }).count() as u32;
        let depth = pileup.depth();
        let pos = pileup.pos(); //0-based
        // println!("pos: {}, depth: {}", pos+1, depth);

        if s == -1 && e == -1 {
            s = pos as i64;
            e = pos as i64;
            ctg = std::str::from_utf8(&header.tid2name(pileup.tid())).unwrap().to_string();
            pre_ctg = ctg.clone();
            if depth >= depth_threshold {
                state = 1;
            } else {
                state = 0;
            }
            pre_state = state;
            pre_depth = depth;
            pre_pos = pos;
            continue;
        }

        if depth >= depth_threshold {
            state = 1;
        } else {
            state = 0;
        }

        if state != pre_state || ctg != pre_ctg || pos != pre_pos + 1 {
            if e > s {
                // println!("{}, depth: {}, pre_depth: {}", pos+1, depth, pre_depth);
                if pre_state == 1 {
                    high_depth_regions.push(Region { chr: ctg.clone(), start: s as u32 + 1, end: e as u32 + 1 });
                    // high_depth_regions.push((ctg.clone(), s + 1, e + 1));
                } else {
                    normal_depth_regions.push(Region { chr: ctg.clone(), start: s as u32 + 1, end: e as u32 + 1 });
                    // normal_depth_regions.push((ctg.clone(), s + 1, e + 1));
                }
            }
            ctg = std::str::from_utf8(&header.tid2name(pileup.tid())).unwrap().to_string();
            pre_ctg = ctg.clone();
            s = pos as i64;
            e = pos as i64;
            pre_state = state;
            pre_depth = depth;
            pre_pos = pos;
        } else {
            e = pos as i64;
            pre_depth = depth;
            pre_pos = pos;
        }
    }
    if pre_state == 1 && e > s {
        high_depth_regions.push(Region { chr: ctg.clone(), start: s as u32 + 1, end: e as u32 + 1 });
    } else if pre_state == 0 && e > s {
        normal_depth_regions.push(Region { chr: ctg.clone(), start: s as u32 + 1, end: e as u32 + 1 });
    }
    return (normal_depth_regions, high_depth_regions);
}

pub fn get_regions_by_coverage(bam_path: String, bed_path: String, depth_threshold: u32) {
    let mut pre_ctg = String::new();
    let mut s: u32 = u32::MAX;
    let mut e: u32 = u32::MAX;
    let mut pre_pos: u32 = u32::MAX;

    let mut bam = bam::IndexedReader::from_path(bam_path).unwrap();
    let header = bam.header().clone();
    let mut writer = File::create(bed_path).unwrap();
    for p in bam.pileup() {
        let pileup = p.unwrap();
        let depth = pileup.alignments().filter(|alignment| {
            let record = alignment.record();
            !alignment.is_refskip() && !record.is_duplicate() && !record.is_secondary() && !record.is_supplementary() && !record.is_unmapped()
        }).count() as u32;
        let pos = pileup.pos(); //0-based
        let ctg = std::str::from_utf8(&header.tid2name(pileup.tid())).unwrap().to_string();

        // first hit
        if depth >= depth_threshold && s == u32::MAX && e == u32::MAX {
            s = pos;
            e = pos;
            pre_pos = pos;
            pre_ctg = ctg;
            continue;
        }

        if depth >= depth_threshold {
            if pos != pre_pos + 1 || ctg != pre_ctg {
                if s != u32::MAX && e != u32::MAX {
                    writeln!(writer, "{}\t{}\t{}", pre_ctg, s, e).unwrap();    // bed is 0-based
                    s = u32::MAX;
                    e = u32::MAX;
                }
            }
            if s == u32::MAX && e == u32::MAX {
                s = pos;
                e = pos;
                pre_pos = pos;
                pre_ctg = ctg;
            } else {
                e = pos;
                pre_pos = pos;
                pre_ctg = ctg;
            }
        } else {
            if s != u32::MAX && e != u32::MAX {
                writeln!(writer, "{}\t{}\t{}", pre_ctg, s, e).unwrap(); // bed is 0-based
            }
            s = u32::MAX;
            e = u32::MAX;
            pre_pos = pos;
            pre_ctg = ctg;
        }
    }
    if s != u32::MAX && e != u32::MAX {
        writeln!(writer, "{}\t{}\t{}", pre_ctg, s, e).unwrap();    // bed is 0-based
    }
}