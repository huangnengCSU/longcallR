mod bam_reader;
mod matrix;
mod pileup2matrix;
mod align;
mod isolated_region;
mod base_matrix;
mod util;

// extern crate bio;
use std::time::{Duration, Instant};
use rust_htslib::{bam, bam::Read};
use bam_reader::{BamReader, Region};
use bam_reader::{write_read_records1, write_read_records2, write_read_records3};
use std::collections::HashMap;
use std::process::exit;
use rust_htslib::bam::record::CigarString;
use bio::io::fasta;
use std::io;
use matrix::PileupMatrix;
use pileup2matrix::generate_pileup_matrix;
use crate::matrix::ColumnBaseCount;
use align::nw_splice_aware;
use isolated_region::find_isolated_regions;
use crate::base_matrix::{*};
use crate::util::{*};
use std::sync::{mpsc, Arc, Mutex};
use std::thread;

fn main3() {
    let main_s = Instant::now();
    let bam_path = "wtc11_ont_grch38.chr22.bam";
    // let region = "chr22:30425877-30425912"; // 1-based
    // let region = "chr22:30425831-30425912";
    // let region = "chr22:37063258-37063422";
    // let region = "chr22:37051423-37063584";
    // let region = "chr22:36651200-36651536";
    // let region = "chr22:34628310-34634051";
    // let region = "chr22:37051212-37073437";
    // let region = "chr22:36961654-37042120";
    // let region = "chr22:26483883-26512499";
    // let region = "chr22:20302014-20324303";
    // let region = "chr22:26483883-26512499";
    let region = "chr22:37009116-37030993";
    // let region = "chr22:20302014-20324303";
    // let region = "chr22:26483883-26512499";
    // let region = "chr22:37830851-38181804";
    let ref_path = "GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.chr22.fna";
    // let out_path = "new.bam";
    // let out_path2 = "new2.bam";
    // let out_path3 = "new3.bam";
    // let bam_reader = BamReader::new(bam_path.to_string(), region.to_string());
    // let mut reader = bam::IndexedReader::from_path(bam_path).unwrap();
    // let header = bam::Header::from_template(reader.header());
    // reader.fetch(("chr22", 18994296, 18994359)).unwrap();
    // let mut writer = bam::Writer::from_path("new.bam", &header, bam::Format::Sam).unwrap();
    // for r in reader.records() {
    //     let record = r.unwrap();
    //     // let mut out_record = bam::Record::from(record.clone());
    //     let mut out_record = bam::Record::default();
    //     let mut new_cigar = CigarString(vec![Cigar::Match(1), Cigar::Ins(1), Cigar::Match(1)]);
    //     out_record.set(record.qname(), Some(&new_cigar), &record.seq().as_bytes(), record.qual());
    //     writer.write(&out_record);
    // }

    let reader_s = Instant::now();
    let bam_reader = BamReader::new(bam_path.to_string(), region.to_string());
    // let mut read_records: HashMap<String, Vec<bam::Record>> = HashMap::new();
    //
    let header = bam_reader.get_bam_header();
    // bam_reader.get_read_records(&mut read_records);

    // let mut writer = bam::Writer::from_path(out_path, &header, bam::Format::Bam).unwrap();
    //
    // for (qname, record_vec) in read_records.iter_mut() {
    //     for record in record_vec.iter_mut() {
    //         let mut out_record = bam::Record::from(record.clone());
    //         out_record.set_pos(123456);
    //         // let mut new_cigar = CigarString(vec![Cigar::Match(1), Cigar::Ins(1), Cigar::Match(1)]);
    //         // out_record.set(record.qname(),
    //         //                Some(&new_cigar),
    //         //                &record.seq().as_bytes(),
    //         //                record.qual());
    //         let new_cigar = record.cigar().take();
    //         out_record.set(record.qname(),
    //                        Some(&new_cigar),
    //                        &record.seq().as_bytes(),
    //                        record.qual());
    //         println!("qname: {:?}, pos: {:?}", qname, record.pos());
    //         let re = writer.write(&out_record);
    //         if re.is_err() {
    //             println!("error: {:?}", re);
    //             exit(1);
    //         }
    //     }
    // }
    // write_read_records1(&read_records, &header, String::from(out_path));
    // write_read_records2(&mut read_records, &header, String::from(out_path2));
    // write_read_records3(&mut read_records, &header, String::from(out_path3));

    let mut reader = fasta::Reader::from_file(ref_path).unwrap();
    for record in reader.records() {
        let record = record.unwrap();
        // println!("{}", record.id());
        // println!("{}", std::str::from_utf8(record.seq()).unwrap());
        // println!("{}", record.seq().len());
    }

    let mut matrices_vec: Vec<PileupMatrix> = Vec::new();
    generate_pileup_matrix(&bam_path.to_string(), &ref_path.to_string(), &region.to_string(), &mut matrices_vec);
    println!("matrices_vec.len(): {:?}", matrices_vec.len());
    // for (readname, seq) in matrices_vec[0].base_matrix.iter() {
    // println!("readname: {:?}", readname);
    // println!("seq: {:?}", String::from_utf8(seq.clone()).unwrap());
    // }
    // let mut v: Vec<_> = map.into_iter().collect();
    // v.sort_by(|x,y| x.0.cmp(&y.0));
    let mut sorted: Vec<_> = matrices_vec[0].positions.clone().into_iter().collect();
    sorted.sort_by(|a, b| a.1.cmp(&b.1));
    println!("max_idx: {:?}", matrices_vec[0].max_idx);
    // println!("current pos: {:?}", matrices_vec[0].current_pos);
    // for (ref_pos, idx) in sorted.iter() {
    //     print!("idx: {:?}\t", idx);
    //     println!("ref_pos: {:?}", ref_pos);
    // }

    let mut column_base_counts: Vec<ColumnBaseCount> = Vec::new();
    let mut column_indexes: Vec<usize> = Vec::new();
    let mut reduced_base_matrix: HashMap<String, Vec<u8>> = HashMap::new();
    let reader_duration = reader_s.elapsed();

    // PileupMatrix::generate_reduced_profile(&matrices_vec[0].base_matrix, &mut column_base_counts, &mut column_indexes, &mut reduced_base_matrix);
    // for cbc in column_base_counts.iter() {
    //     println!("A:{}, C:{}, G:{}, T:{}, N:{}, D:{}, B:{}", cbc.n_a, cbc.n_c, cbc.n_g, cbc.n_t, cbc.n_n, cbc.n_dash, cbc.n_blank);
    // }
    // for i in 0..column_base_counts.len() {
    //     let cbc = &column_base_counts[i];
    //     println!("idx:{}\tA:{}, C:{}, G:{}, T:{}, N:{}, D:{}, B:{}", &column_indexes[i], cbc.n_a, cbc.n_c, cbc.n_g, cbc.n_t, cbc.n_n, cbc.n_dash, cbc.n_blank);
    // }

    // let query = String::from_utf8(matrices_vec[0].base_matrix.iter().next().unwrap().1.to_vec()).unwrap().replace("N", "").replace(" ", "");
    // let (alignment_score, aligned_query, ref_target, major_target) = nw_splice_aware(&query, &column_base_counts);
    // println!("alignment_score: {:?}", alignment_score);
    // println!("aligned_query:\n {:?}", aligned_query);
    // println!("ref_target:\n {:?}", ref_target);
    // println!("major_target:\n {:?}", major_target);
    let mut total_align_runtime: u64 = 0;
    let mut total_update_runtime: u64 = 0;
    for i in 0..matrices_vec.len() {
        let (forward_donor_penalty,
            forward_acceptor_penalty,
            reverse_donor_penalty,
            reverse_acceptor_penalty) = matrices_vec[i].get_donor_acceptor_penalty(30.0);
        let mut best_reduced_base_matrix: HashMap<String, Vec<u8>> = HashMap::new();
        let mut best_column_indexes: Vec<usize> = Vec::new();
        let align_s = Instant::now();
        PileupMatrix::profile_realign(&matrices_vec[i].base_matrix,
                                      &forward_donor_penalty,
                                      &forward_acceptor_penalty,
                                      &reverse_donor_penalty,
                                      &reverse_acceptor_penalty,
                                      &mut best_reduced_base_matrix,
                                      &mut best_column_indexes);
        let align_duration = align_s.elapsed().as_secs();
        total_align_runtime += align_duration;
        // let seq = std::str::from_utf8(matrices_vec[i].base_matrix.iter().next().unwrap().1).unwrap().to_string();
        // println!("seq: \n{:?}", seq);
        let update_s = Instant::now();
        matrices_vec[i].update_base_matrix_from_realign(&best_reduced_base_matrix, &best_column_indexes);
        let it = matrices_vec[i].base_matrix.iter().next().unwrap();
        let qname = it.0;
        let seq = std::str::from_utf8(it.1).unwrap().to_string();
        // println!("new seq: {} \n{:?}", qname, seq);
        // println!("start pos: {}, end pos: {}", matrices_vec[i].region.start, matrices_vec[i].region.end);

        matrices_vec[i].update_bam_records_from_realign();
        let update_duration = update_s.elapsed().as_secs();
        total_update_runtime += update_duration;
        matrices_vec[i].write_bam_records("new4.bam", &header);
    }
    let main_duration = main_s.elapsed();
    println!("main: {:?}", main_duration);
    println!("reader: {:?}", reader_duration);
    println!("align: {:?}", total_align_runtime);
    println!("update: {:?}", total_update_runtime);
}

fn main2() {
    let bam_path = "wtc11_ont_grch38.chr22.bam";
    let region = "chr22:37009116-37030993";
    // let region = "chr22:26483883-26512499";
    let ref_path = "GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.chr22.fna";
    let mut matrices_vec: Vec<PileupMatrix> = Vec::new();
    generate_pileup_matrix(&bam_path.to_string(), &ref_path.to_string(), &region.to_string(), &mut matrices_vec);

    for i in 0..matrices_vec.len() {
        let mut column_base_counts: Vec<ColumnBaseCount> = Vec::new();
        let mut column_indexes: Vec<usize> = Vec::new();
        let mut reduced_base_matrix: HashMap<String, Vec<u8>> = HashMap::new();
        let mut forward_reduced_donor_penalty: Vec<f64> = Vec::new();
        let mut forward_reduced_acceptor_penalty: Vec<f64> = Vec::new();
        let mut reverse_reduced_donor_penalty: Vec<f64> = Vec::new();
        let mut reverse_reduced_acceptor_penalty: Vec<f64> = Vec::new();
        let mut splice_boundary: Vec<bool> = Vec::new();
        let (forward_donor_penalty, forward_acceptor_penalty, reverse_donor_penalty, reverse_acceptor_penalty) = matrices_vec[i].get_donor_acceptor_penalty(30.0);
        PileupMatrix::generate_reduced_profile(&matrices_vec[i].base_matrix,
                                               &forward_donor_penalty,
                                               &forward_acceptor_penalty,
                                               &reverse_donor_penalty,
                                               &reverse_acceptor_penalty,
                                               &mut column_base_counts,
                                               &mut column_indexes,
                                               &mut reduced_base_matrix,
                                               &mut forward_reduced_donor_penalty,
                                               &mut forward_reduced_acceptor_penalty,
                                               &mut reverse_reduced_donor_penalty,
                                               &mut reverse_reduced_acceptor_penalty,
                                               &mut splice_boundary);
        println!("reference:");
        for d in matrices_vec[i].base_matrix.get("ref").unwrap().iter() {
            print!("{}\t", *d as char);
        }
        println!();
        println!("forward donor_penalty:");
        for d in forward_donor_penalty.iter() {
            print!("{}\t", d);
        }
        println!();
        println!("forward acceptor_penalty:");
        for d in forward_acceptor_penalty.iter() {
            print!("{}\t", d);
        }
        println!();
        println!("reverse donor_penalty:");
        for d in reverse_donor_penalty.iter() {
            print!("{}\t", d);
        }
        println!();
        println!("reverse acceptor_penalty:");
        for d in reverse_acceptor_penalty.iter() {
            print!("{}\t", d);
        }
        println!();
        println!("reduced reference:");
        for d in reduced_base_matrix.get("ref").unwrap().iter() {
            print!("{}\t", *d as char);
        }
        println!();
        println!("forward reduced donor_penalty:");
        for d in forward_reduced_donor_penalty.iter() {
            print!("{}\t", d);
        }
        println!();
        println!("forward reduced acceptor_penalty:");
        for d in forward_reduced_acceptor_penalty.iter() {
            print!("{}\t", d);
        }
        println!();
        println!("reverse reduced donor_penalty:");
        for d in reverse_reduced_donor_penalty.iter() {
            print!("{}\t", d);
        }
        println!();
        println!("reverse reduced acceptor_penalty:");
        for d in reverse_reduced_acceptor_penalty.iter() {
            print!("{}\t", d);
        }
        println!();
        for d in splice_boundary.iter() {
            print!("{}\t", *d as u8);
        }
        println!();
    }
}

fn main5() {
    let bam_path = "wtc11_ont_grch38.chr22.bam";
    let regions = find_isolated_regions(bam_path, 2000);
    for region in regions.iter() {
        println!("{}:{}-{}", region.chr, region.start, region.end);
    }
}


fn main4() {
    let main_s = Instant::now();
    let bam_path = "wtc11_ont_grch38.chr22.bam";
    // let input_region = "chr22:37009116-37030993";
    let input_region = "chr22:42619102-43955456";
    let ref_path = "GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.chr22.fna";
    let reader_s = Instant::now();

    let mut base_matrix = BaseMatrix::new();
    let region = Region::new(input_region.to_string());
    // base_matrix.load_data(bam_path.to_string(), region);
    base_matrix.load_data_without_extension(bam_path.to_string(), region);
    let ref_seqs = load_reference(ref_path.to_string());
    base_matrix.load_ref_data(ref_seqs);
    base_matrix.expand_insertion();
    let reader_duration = reader_s.elapsed();
    let mut total_align_runtime: u64 = 0;
    let mut total_update_runtime: u64 = 0;

    let (forward_donor_penalty,
        forward_acceptor_penalty,
        reverse_donor_penalty,
        reverse_acceptor_penalty) = get_donor_acceptor_penalty(&base_matrix.expanded_matrix, 30.0);
    let mut best_reduced_expanded_matrix: HashMap<String, Vec<u8>> = HashMap::new();
    let mut best_column_indexes: Vec<usize> = Vec::new();
    let align_s = Instant::now();
    profile_realign(&base_matrix.expanded_matrix,
                    &forward_donor_penalty,
                    &forward_acceptor_penalty,
                    &reverse_donor_penalty,
                    &reverse_acceptor_penalty,
                    &mut best_reduced_expanded_matrix,
                    &mut best_column_indexes);
    let align_duration = align_s.elapsed().as_secs();
    total_align_runtime += align_duration;
    let update_s = Instant::now();
    update_expanded_matrix_from_realign(&mut base_matrix.expanded_matrix, &best_reduced_expanded_matrix, &best_column_indexes);
    update_bam_records_from_realign(&mut base_matrix.expanded_matrix, &mut base_matrix.bam_records, base_matrix.start_position, base_matrix.end_position);
    let update_duration = update_s.elapsed().as_secs();
    total_update_runtime += update_duration;
    let bam: bam::IndexedReader = bam::IndexedReader::from_path(bam_path).unwrap();
    let header = bam::Header::from_template(bam.header());
    write_bam_records(&mut base_matrix.bam_records, "new4.bam", &header);
    let main_duration = main_s.elapsed();
    println!("main: {:?}", main_duration);
    println!("reader: {:?}", reader_duration);
    println!("align: {:?}", total_align_runtime);
    println!("update: {:?}", total_update_runtime);
    println!("start pos: {}, end pos: {}", base_matrix.start_position, base_matrix.end_position);
}


fn main() {
    let bam_path = "wtc11_ont_grch38.chr22.bam";
    let (tx_l, rx_l) = mpsc::channel();
    let (tx_h, rx_h) = mpsc::channel();
    let prod_v = multithread_produce(bam_path.to_string(), 4, tx_l, tx_h);

    for handle in prod_v {
        handle.join().unwrap();
    }

    for region in rx_l {
        println!("low: {:?}", region);
    }
    for region in rx_h {
        println!("high: {:?}", region);
    }
}
