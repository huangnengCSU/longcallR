use std::cmp::min;
use rust_htslib::{bam, bam::Format, bam::{Read, ext::BamRecordExtensions}};
use std::sync::{mpsc, Arc, Mutex, Condvar};
use rayon::prelude::*;
use std::{fs, process};
use std::collections::{VecDeque, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::time::Instant;
use bio::bio_types::strand::ReqStrand::Forward;
use chrono::Local;
use clap::builder::Str;
use seq_io::fasta::{Reader, Record};
use threadpool::ThreadPool;
use crate::bam_reader::Region;
use crate::base_matrix::BaseMatrix;
use crate::base_matrix::{*};
use crate::isolated_region::find_isolated_regions;
use crate::profile::{BaseFreq, BaseQual, BaseStrands, DistanceToEnd};


pub fn parse_fai(fai_path: &str) -> Vec<(String, u32)> {
    let mut contig_lengths: Vec<(String, u32)> = Vec::new();
    let file = File::open(fai_path).unwrap();
    let reader = BufReader::new(file);
    for r in reader.lines() {
        let line = r.unwrap().clone();
        let parts: Vec<&str> = line.split('\t').collect();
        contig_lengths.push((parts[0].to_string(), parts[1].parse().unwrap()));
    }
    return contig_lengths;
}

pub fn find_isolated_regions_with_depth(bam_path: &str, chr: &str, ref_len: u32) -> Vec<Region> {
    let mut isolated_regions: Vec<Region> = Vec::new();
    let mut depth_vec: Vec<u32> = vec![0; ref_len as usize];
    let mut bam: bam::IndexedReader = bam::IndexedReader::from_path(bam_path).unwrap();
    let header = bam.header().clone();
    bam.fetch(chr).unwrap();
    for r in bam.records() {
        let record = r.unwrap();
        let ref_start = record.reference_start();   // 0-based, left-closed
        let ref_end = record.reference_end();   // 0-based, right-open
        for i in ref_start..ref_end {
            depth_vec[i as usize] += 1;
        }
    }
    let mut region_start = -1;
    let mut region_end = -1;
    for i in 0..ref_len {
        if depth_vec[i as usize] == 0 {
            if region_end > region_start {
                assert!(region_start >= 0);
                assert!(region_end >= 0);
                isolated_regions.push(Region { chr: chr.to_string(), start: (region_start + 1) as u32, end: (region_end + 2) as u32 });
                region_start = -1;
                region_end = -1;
            }
        } else {
            if region_start == -1 {
                region_start = i as i32;
                region_end = i as i32;
            } else {
                region_end = i as i32;
            }
        }
    }
    if region_end > region_start {
        isolated_regions.push(Region { chr: chr.to_string(), start: (region_start + 1) as u32, end: (region_end + 2) as u32 });
        region_start = -1;
        region_end = -1;
    }
    return isolated_regions;
}


pub fn multithread_produce(bam_file: String, thread_size: usize, tx_low: mpsc::Sender<Region>, tx_high: mpsc::Sender<Region>) {
    let pool = ThreadPool::new(thread_size);
    let bam = bam::IndexedReader::from_path(bam_file.clone()).unwrap();
    let bam_header = bam.header().clone();
    let mut contig_names: VecDeque<String> = VecDeque::new();
    for ctg in bam_header.target_names() {
        contig_names.push_back(std::str::from_utf8(ctg).unwrap().to_string().clone());
    }
    while !contig_names.is_empty() {
        let ctg = contig_names.pop_front().unwrap();
        let tx_l = tx_low.clone();
        let tx_h = tx_high.clone();
        let bam_file_clone = bam_file.clone();
        pool.execute(move || {
            let (normal_depth_regions, high_depth_regions) = get_chrom_coverage_intervals(bam_file_clone.clone(), ctg.as_str(), 1000);
            for region in normal_depth_regions {
                tx_l.send(region).unwrap();
            }
            for region in high_depth_regions {
                tx_h.send(region).unwrap();
            }
        });
    }
    pool.join();
}

pub fn multithread_work(bam_file: String, ref_file: String, out_bam: String, thread_size: usize, rx_low: mpsc::Receiver<Region>, rx_high: mpsc::Receiver<Region>) {
    let pool = ThreadPool::new(thread_size);
    let bam_records_queue: Arc<Mutex<VecDeque<bam::Record>>> = Arc::new(Mutex::new(VecDeque::new()));
    let bam: bam::IndexedReader = bam::IndexedReader::from_path(bam_file.clone()).unwrap();
    let header = bam::Header::from_template(bam.header());
    let bam_writer = bam::Writer::from_path(out_bam.clone(), &header, Format::Bam).unwrap();
    let mut bam_writer = Arc::new(Mutex::new(bam_writer));
    let ref_seqs = load_reference(ref_file.to_string());

    // multiplethreads for low coverage regions
    for normal_depth_region in rx_low {
        let bam_file_clone = bam_file.clone();
        let ref_seqs_clone = ref_seqs.clone();
        let bam_records_queue = Arc::clone(&bam_records_queue);
        let bam_writer_clone = Arc::clone(&bam_writer);
        pool.execute(move || {
            let mut base_matrix = BaseMatrix::new();
            base_matrix.load_data(bam_file_clone.clone().to_string(), normal_depth_region);
            base_matrix.load_ref_data(ref_seqs_clone.clone());
            base_matrix.expand_insertion();
            let (forward_donor_penalty,
                forward_acceptor_penalty,
                reverse_donor_penalty,
                reverse_acceptor_penalty) = get_donor_acceptor_penalty(&base_matrix.expanded_matrix, 30.0);
            let mut best_reduced_expanded_matrix: HashMap<String, Vec<u8>> = HashMap::new();
            let mut best_column_indexes: Vec<usize> = Vec::new();
            profile_realign(&base_matrix.expanded_matrix,
                            &forward_donor_penalty,
                            &forward_acceptor_penalty,
                            &reverse_donor_penalty,
                            &reverse_acceptor_penalty,
                            &mut best_reduced_expanded_matrix,
                            &mut best_column_indexes);
            update_expanded_matrix_from_realign(&mut base_matrix.expanded_matrix, &best_reduced_expanded_matrix, &best_column_indexes);
            update_bam_records_from_realign(&mut base_matrix.expanded_matrix, &mut base_matrix.bam_records, base_matrix.start_position, base_matrix.end_position);
            let mut queue = bam_records_queue.lock().unwrap();
            for (_, record) in base_matrix.bam_records {
                queue.push_back(record);
            }
            if queue.len() > 100000 {
                for record in queue.iter() {
                    let re = bam_writer_clone.lock().unwrap().write(&record);
                    if re != Ok(()) {
                        println!("write failed");
                        process::exit(1);
                    }
                }
                queue.clear();
            }
        });
    }
    pool.join();
    if bam_records_queue.lock().unwrap().len() > 0 {
        for record in bam_records_queue.lock().unwrap().iter() {
            let re = bam_writer.lock().unwrap().write(&record);
            if re != Ok(()) {
                println!("write failed");
                process::exit(1);
            }
        }
        bam_records_queue.lock().unwrap().clear();
    }

    let pool = ThreadPool::new(thread_size);
    for high_depth_region in rx_high {
        let bam_file_clone = bam_file.clone();
        let ref_seqs_clone = ref_seqs.clone();
        let bam_records_queue = Arc::clone(&bam_records_queue);
        let bam_writer_clone = Arc::clone(&bam_writer);
        pool.execute(move || {
            let mut base_matrix = BaseMatrix::new();
            base_matrix.load_data_without_extension(bam_file_clone.clone().to_string(), high_depth_region);
            base_matrix.load_ref_data(ref_seqs_clone.clone());
            base_matrix.expand_insertion();
            let (forward_donor_penalty,
                forward_acceptor_penalty,
                reverse_donor_penalty,
                reverse_acceptor_penalty) = get_donor_acceptor_penalty(&base_matrix.expanded_matrix, 30.0);
            let mut best_reduced_expanded_matrix: HashMap<String, Vec<u8>> = HashMap::new();
            let mut best_column_indexes: Vec<usize> = Vec::new();
            profile_realign(&base_matrix.expanded_matrix,
                            &forward_donor_penalty,
                            &forward_acceptor_penalty,
                            &reverse_donor_penalty,
                            &reverse_acceptor_penalty,
                            &mut best_reduced_expanded_matrix,
                            &mut best_column_indexes);
            update_expanded_matrix_from_realign(&mut base_matrix.expanded_matrix, &best_reduced_expanded_matrix, &best_column_indexes);
            update_bam_records_from_realign(&mut base_matrix.expanded_matrix, &mut base_matrix.bam_records, base_matrix.start_position, base_matrix.end_position);
            let mut queue = bam_records_queue.lock().unwrap();
            for (_, record) in base_matrix.bam_records {
                queue.push_back(record);
            }
            if queue.len() > 10000 {
                for record in queue.iter() {
                    let re = bam_writer_clone.lock().unwrap().write(&record);
                    if re != Ok(()) {
                        println!("write failed");
                        process::exit(1);
                    }
                }
                queue.clear();
            }
        });
    }
    pool.join();
    if bam_records_queue.lock().unwrap().len() > 0 {
        for record in bam_records_queue.lock().unwrap().iter() {
            let re = bam_writer.lock().unwrap().write(&record);
            if re != Ok(()) {
                println!("write failed");
                process::exit(1);
            }
        }
        bam_records_queue.lock().unwrap().clear();
    }
}


pub fn multithread_produce2(bam_file: String, thread_size: usize, tx_isolated_regions: mpsc::Sender<Region>) {
    let pool = ThreadPool::new(&thread_size - 1);
    let cond = Arc::new((Mutex::new(()), Condvar::new()));
    let bam = bam::IndexedReader::from_path(bam_file.clone()).unwrap();
    let bam_header = bam.header().clone();
    let mut contig_names: VecDeque<String> = VecDeque::new();
    for ctg in bam_header.target_names() {
        contig_names.push_back(std::str::from_utf8(ctg).unwrap().to_string().clone());
    }
    while !contig_names.is_empty() {
        let c = cond.clone();
        let ctg = contig_names.pop_front().unwrap();
        let tx_ir = tx_isolated_regions.clone();
        let bam_file_clone = bam_file.clone();
        pool.execute(move || {
            {
                let (lock, cvar) = &*c;
                let _guard = lock.lock().unwrap();
                cvar.notify_one();
            }
            let isolated_regions = find_isolated_regions(bam_file_clone.clone().as_str(), 1, Some(ctg.as_str()));
            for region in isolated_regions {
                tx_ir.send(region).unwrap();
            }
            std::thread::sleep(std::time::Duration::from_millis(1));
        });
        let (lock, cvar) = &*cond;
        let guard = lock.lock().unwrap();
        if pool.queued_count() >= &thread_size - 1 {
            let tt = cvar.wait(guard).unwrap();
        }
    }
    pool.join();
}

pub fn multithread_produce3(bam_file: String, ref_file: String, thread_size: usize, contigs: Option<Vec<String>>) -> Vec<Region> {
    let results: Mutex<Vec<Region>> = Mutex::new(Vec::new());
    let pool = rayon::ThreadPoolBuilder::new().num_threads(thread_size - 1).build().unwrap();
    let bam = bam::IndexedReader::from_path(bam_file.clone()).unwrap();
    let bam_header = bam.header().clone();
    let fai_path = ref_file + ".fai";
    if fs::metadata(&fai_path).is_err() {
        panic!("Reference index file .fai does not exist.");
    }
    let contig_lengths = parse_fai(fai_path.as_str());
    let mut contig_names: VecDeque<String> = VecDeque::new();
    if contigs.is_some() {
        for ctg in contigs.unwrap().iter() {
            contig_names.push_back(ctg.clone());
        }
    } else {
        // for ctg in bam_header.target_names() {
        //     contig_names.push_back(std::str::from_utf8(ctg).unwrap().to_string().clone());
        // }
        for (ctg, _) in contig_lengths.iter() {
            contig_names.push_back(ctg.clone());
        }
    }
    pool.install(|| {
        contig_names.par_iter().for_each(|ctg| {
            let isolated_regions = find_isolated_regions_with_depth(bam_file.as_str(), ctg, contig_lengths.iter().find(|(chr, _)| chr == ctg).unwrap().1);
            for region in isolated_regions {
                results.lock().unwrap().push(region);
            }
            println!("{:?} finished getting isolated regions", ctg);
        });
    });
    return results.into_inner().unwrap().clone();
}

// pub fn multithread_work2(bam_file: String, ref_file: String, out_bam: String, thread_size: usize, rx_isolated_regions: mpsc::Receiver<Region>, min_mapq: u8) {
//     let pool = ThreadPool::new(&thread_size - 1);
//     let cond = Arc::new((Mutex::new(()), Condvar::new()));
//     let bam_records_queue: Arc<Mutex<VecDeque<bam::Record>>> = Arc::new(Mutex::new(VecDeque::new()));
//     let bam: bam::IndexedReader = bam::IndexedReader::from_path(bam_file.clone()).unwrap();
//     let header = bam::Header::from_template(bam.header());
//     let bam_writer = bam::Writer::from_path(out_bam.clone(), &header, Format::Bam).unwrap();
//     let mut bam_writer = Arc::new(Mutex::new(bam_writer));
//     let ref_seqs = load_reference(ref_file.to_string());
//
//     // multiplethreads for low coverage regions
//     for reg in rx_isolated_regions {
//         let c = cond.clone();
//         let bam_file_clone = bam_file.clone();
//         let ref_seqs_clone = ref_seqs.clone();
//         let bam_records_queue = Arc::clone(&bam_records_queue);
//         let bam_writer_clone = Arc::clone(&bam_writer);
//         pool.execute(move || {
//             {
//                 let (lock, cvar) = &*c;
//                 let _guard = lock.lock().unwrap();
//                 cvar.notify_one();
//             }
//             println!("Start {:?}", reg);
//             let mut profile = Profile::default();
//             let mut readnames: Vec<String> = Vec::new();
//             profile.init_with_pileup(bam_file_clone.clone().as_str(), &reg, ref_seqs_clone.get(&reg.chr).unwrap(), min_mapq, 7, 10, 10000, 2000, 20, 5);
//             profile.append_reference(&ref_seqs_clone);
//             let mut parsed_reads = read_bam(bam_file_clone.clone().as_str(), &reg);
//             for (rname, pr) in parsed_reads.iter_mut() {
//                 pr.init_parsed_seq(&profile);
//                 readnames.push(rname.clone());
//             }
//             // println!("Total number of reads: {}", &parsed_reads.len());
//             profile.cal_intron_penalty();
//             profile.cal_intron_intervals();
//             // for bf in profile.freq_vec.iter() {
//             //     println!("{:?}", bf);
//             // }
//             realign(&mut profile, &mut parsed_reads, &readnames);
//             // let header = get_bam_header(bam_file_clone.clone().as_str());
//             // write_bam(out_bam, &parsed_reads, &header);
//             {
//                 let mut queue = bam_records_queue.lock().unwrap();
//                 for (_, pr) in parsed_reads.iter() {
//                     queue.push_back(pr.bam_record.clone());
//                 }
//                 if queue.len() > 100000 {
//                     for record in queue.iter() {
//                         let re = bam_writer_clone.lock().unwrap().write(&record);
//                         if re != Ok(()) {
//                             println!("write failed");
//                             process::exit(1);
//                         }
//                     }
//                     queue.clear();
//                 }
//             }
//             println!("End {:?}", reg);
//             std::thread::sleep(std::time::Duration::from_millis(1));
//         });
//         let (lock, cvar) = &*cond;
//         let guard = lock.lock().unwrap();
//         if pool.queued_count() >= &thread_size - 1 {
//             let tt = cvar.wait(guard).unwrap();
//         }
//     }
//     pool.join();
//     if bam_records_queue.lock().unwrap().len() > 0 {
//         for record in bam_records_queue.lock().unwrap().iter() {
//             let re = bam_writer.lock().unwrap().write(&record);
//             if re != Ok(()) {
//                 println!("write failed");
//                 process::exit(1);
//             }
//         }
//         bam_records_queue.lock().unwrap().clear();
//     }
// }


// pub fn multithread_work3(bam_file: String, ref_file: String, out_bam: String, thread_size: usize, isolated_regions: Vec<Region>, min_mapq: u8) {
//     let pool = rayon::ThreadPoolBuilder::new().num_threads(thread_size - 1).build().unwrap();
//     let bam_records_queue = Mutex::new(VecDeque::new());
//     let bam: bam::IndexedReader = bam::IndexedReader::from_path(&bam_file).unwrap();
//     let header = bam::Header::from_template(bam.header());
//     let mut bam_writer = Mutex::new(bam::Writer::from_path(&out_bam, &header, Format::Bam).unwrap());
//     let ref_seqs = load_reference(ref_file);
//
//     // multiplethreads for low coverage regions
//     pool.install(|| {
//         isolated_regions.par_iter().for_each(|reg| {
//             println!("Start {:?}", reg);
//             let mut profile = Profile::default();
//             let mut readnames: Vec<String> = Vec::new();
//             profile.init_with_pileup(&bam_file.as_str(), &reg, ref_seqs.get(&reg.chr).unwrap(), min_mapq, 10, 10000, 2000, 20, 5);
//             profile.append_reference(&ref_seqs);
//             let mut parsed_reads = read_bam(&bam_file.as_str(), &reg);
//             for (rname, pr) in parsed_reads.iter_mut() {
//                 pr.init_parsed_seq(&profile);
//                 readnames.push(rname.clone());
//             }
//             profile.cal_intron_penalty();
//             profile.cal_intron_intervals();
//             realign(&mut profile, &mut parsed_reads, &readnames);
//             {
//                 let mut queue = bam_records_queue.lock().unwrap();
//                 for (_, pr) in parsed_reads.iter() {
//                     queue.push_back(pr.bam_record.clone());
//                 }
//                 if queue.len() > 100000 {
//                     for record in queue.iter() {
//                         let re = bam_writer.lock().unwrap().write(&record);
//                         if re != Ok(()) {
//                             println!("write failed");
//                             process::exit(1);
//                         }
//                     }
//                     queue.clear();
//                 }
//             }
//             println!("End {:?}", reg);
//         });
//     });
//
//     if bam_records_queue.lock().unwrap().len() > 0 {
//         for record in bam_records_queue.lock().unwrap().iter() {
//             let re = bam_writer.lock().unwrap().write(&record);
//             if re != Ok(()) {
//                 println!("write failed");
//                 process::exit(1);
//             }
//         }
//         bam_records_queue.lock().unwrap().clear();
//     }
// }

pub fn read_references(ref_path: &str) -> HashMap<String, Vec<u8>> {
    let mut references: HashMap<String, Vec<u8>> = HashMap::new();
    let mut reader = Reader::from_path(ref_path).unwrap();
    while let Some(record) = reader.next() {
        let record = record.expect("Error reading record");
        references.insert(record.id().unwrap().to_string(), record.full_seq().to_vec());
    }
    return references;
}

#[derive(Default, Debug, Clone)]
pub struct Profile {
    pub freq_vec: Vec<BaseFreq>,
    pub region: Region,
}

impl Profile {
    pub fn init_with_pileup(&mut self, bam_path: &str, region: &Region, ref_seq: &Vec<u8>, platform: &String, min_mapq: u8, min_baseq: u8, min_read_length: usize, min_depth: u32, max_depth: u32, distance_to_read_end: u32, polya_tail_length: u32) {
        // When region is large and the number of reads is large, the runtime of init_profile_with_pileup is time-consuming.
        // This function is used to fill the profile by parsing each read in the bam file instead of using pileup.

        let start_time = Instant::now();
        println!("{} Start to construct FreqVec for region: {}:{}-{}", Local::now().format("%Y-%m-%d %H:%M:%S"), region.chr, region.start, region.end);
        let mut bam: bam::IndexedReader = bam::IndexedReader::from_path(bam_path).unwrap();
        bam.fetch((region.chr.as_str(), region.start, region.end)).unwrap();
        let vec_size = (region.end - region.start) as usize;    // end is exclusive
        self.freq_vec = vec![BaseFreq::default(); vec_size];
        self.region = region.clone();
        let freq_vec_pos = region.start as usize - 1;    // the first position on reference, 0-based, inclusive
        let polyA_win = polya_tail_length as i64;

        // fill the ref_base field in each BaseFreq
        for i in 0..vec_size {
            self.freq_vec[i].ref_base = ref_seq[freq_vec_pos + i] as char;
        }

        for r in bam.records() {
            let record = r.unwrap();
            if record.mapq() < min_mapq || record.seq_len() < min_read_length || record.is_unmapped() || record.is_secondary() || record.is_supplementary() {
                continue;
            }
            let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
            let seq = record.seq();
            let base_qual = record.qual();
            let strand = if record.strand() == Forward { 0 } else { 1 };
            let start_pos = record.pos() as usize;  // 0-based
            let cigar = record.cigar();
            let leading_softclips = cigar.leading_softclips();
            let trailing_softclips = cigar.trailing_softclips();

            let mut pos_in_freq_vec = start_pos - freq_vec_pos;
            let mut pos_in_read = if leading_softclips > 0 { leading_softclips as usize } else { 0 };
            for cg in cigar.iter() {
                match cg.char() as u8 {
                    b'S' | b'H' => {
                        continue;
                    }
                    b'M' | b'X' | b'=' => {
                        for cgi in 0..cg.len() {
                            let base = seq[pos_in_read] as char;
                            let baseq = base_qual[pos_in_read];

                            // close to left read end or right read end, check whether current position is in polyA tail
                            // let ref_base = ref_seq[pos_in_freq_vec] as char;
                            let ref_base = self.freq_vec[pos_in_freq_vec].ref_base;
                            let mut polyA_flag = false;
                            let mut homopolymer_flag = false;
                            let mut trimed_flag = false;    // trime the end of ont reads
                            // the end of ont reads are trimed since the accuracy of the ends is low
                            if platform == "ont" && ((pos_in_read as i64 - leading_softclips).abs() < distance_to_read_end as i64 || (pos_in_read as i64 - (seq.len() as i64 - trailing_softclips)).abs() < distance_to_read_end as i64) {
                                trimed_flag = true;
                            }
                            if !trimed_flag && ((pos_in_read as i64 - leading_softclips).abs() < distance_to_read_end as i64 || (pos_in_read as i64 - (seq.len() as i64 - trailing_softclips)).abs() < distance_to_read_end as i64) {
                                for tmpi in (pos_in_read as i64 - polyA_win)..=(pos_in_read as i64 + 1) {
                                    // pos_in_read is the current position, and the position 1-base to the left of polyA tail is often false positive variant allele. So the end for loop is pos_in_read+1 instead of pos_in_read.
                                    // same reason for pos_in_read - polyA_win instead fo pos_in_read - polyA_win + 1
                                    if tmpi < 0 || tmpi + polyA_win - 1 >= seq.len() as i64 {
                                        continue;
                                    }
                                    let mut polyA_cnt = 0;
                                    let mut polyT_cnt = 0;
                                    let mut polyC_cnt = 0;
                                    let mut polyG_cnt = 0;
                                    for tmpj in 0..polyA_win {
                                        if seq[(tmpi + tmpj) as usize] == b'A' && ref_base != 'A' {
                                            polyA_cnt += 1;
                                        } else if seq[(tmpi + tmpj) as usize] == b'T' && ref_base != 'T' {
                                            polyT_cnt += 1;
                                        } else if seq[(tmpi + tmpj) as usize] == b'C' && ref_base != 'C' {
                                            polyC_cnt += 1;
                                        } else if seq[(tmpi + tmpj) as usize] == b'G' && ref_base != 'G' {
                                            polyG_cnt += 1;
                                        }
                                    }
                                    if polyA_cnt >= polyA_win || polyT_cnt >= polyA_win {
                                        polyA_flag = true;
                                    }
                                    if polyC_cnt >= polyA_win || polyG_cnt >= polyA_win {
                                        homopolymer_flag = true;
                                    }
                                }
                            }

                            if !trimed_flag && !polyA_flag && !homopolymer_flag {

                                // calculate distance to read end of each allele, for filtering variants that the average distance of each allele is significantly different
                                let mut dist = 0;
                                if (pos_in_read as i64 - leading_softclips).abs() < (pos_in_read as i64 - (seq.len() as i64 - trailing_softclips)).abs() {
                                    dist = pos_in_read as i64 - leading_softclips;    // positive value
                                } else {
                                    dist = (pos_in_read as i64 - (seq.len() as i64 - trailing_softclips));    // negative value
                                }

                                match base {
                                    'A' | 'a' => {
                                        self.freq_vec[pos_in_freq_vec].a += 1;
                                        self.freq_vec[pos_in_freq_vec].baseq.a.push(baseq);
                                        if strand == 0 {
                                            self.freq_vec[pos_in_freq_vec].base_strands.a[0] += 1;
                                        } else {
                                            self.freq_vec[pos_in_freq_vec].base_strands.a[1] += 1;
                                        }
                                        self.freq_vec[pos_in_freq_vec].distance_to_end.a.push(dist);
                                    }
                                    'C' | 'c' => {
                                        self.freq_vec[pos_in_freq_vec].c += 1;
                                        self.freq_vec[pos_in_freq_vec].baseq.c.push(baseq);
                                        if strand == 0 {
                                            self.freq_vec[pos_in_freq_vec].base_strands.c[0] += 1;
                                        } else {
                                            self.freq_vec[pos_in_freq_vec].base_strands.c[1] += 1;
                                        }
                                        self.freq_vec[pos_in_freq_vec].distance_to_end.c.push(dist);
                                    }
                                    'G' | 'g' => {
                                        self.freq_vec[pos_in_freq_vec].g += 1;
                                        self.freq_vec[pos_in_freq_vec].baseq.g.push(baseq);
                                        if strand == 0 {
                                            self.freq_vec[pos_in_freq_vec].base_strands.g[0] += 1;
                                        } else {
                                            self.freq_vec[pos_in_freq_vec].base_strands.g[1] += 1;
                                        }
                                        self.freq_vec[pos_in_freq_vec].distance_to_end.g.push(dist);
                                    }
                                    'T' | 't' => {
                                        self.freq_vec[pos_in_freq_vec].t += 1;
                                        self.freq_vec[pos_in_freq_vec].baseq.t.push(baseq);
                                        if strand == 0 {
                                            self.freq_vec[pos_in_freq_vec].base_strands.t[0] += 1;
                                        } else {
                                            self.freq_vec[pos_in_freq_vec].base_strands.t[1] += 1;
                                        }
                                        self.freq_vec[pos_in_freq_vec].distance_to_end.t.push(dist);
                                    }
                                    _ => {
                                        panic!("Invalid nucleotide base: {}", base);
                                    }
                                }
                                if strand == 0 {
                                    self.freq_vec[pos_in_freq_vec].forward_cnt += 1;
                                } else {
                                    self.freq_vec[pos_in_freq_vec].backward_cnt += 1;
                                }
                            }

                            pos_in_freq_vec += 1;
                            pos_in_read += 1;
                        }
                    }
                    b'D' => {
                        for _ in 0..cg.len() {
                            self.freq_vec[pos_in_freq_vec].d += 1;
                            pos_in_freq_vec += 1;
                        }
                    }
                    b'I' => {
                        self.freq_vec[pos_in_freq_vec - 1].ni += 1; // insertion is counted as the previous position
                        pos_in_read += cg.len() as usize;
                    }
                    b'N' => {
                        for _ in 0..cg.len() {
                            self.freq_vec[pos_in_freq_vec].n += 1;
                            pos_in_freq_vec += 1;
                        }
                    }
                    _ => {
                        panic!("Error: unknown cigar operation: {}", cg.char());
                    }
                }
            }
        }
        let end_time = Instant::now();
        println!("Init profile with pileup: {}:{}-{} {} seconds", self.region.chr, self.region.start, self.region.end, end_time.duration_since(start_time).as_secs());
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
}