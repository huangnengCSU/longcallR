mod bam_reader;
mod matrix;
mod pileup2matrix;

// extern crate bio;

use rust_htslib::{bam, bam::Read};
use bam_reader::{BamReader};
use bam_reader::{write_read_records1, write_read_records2, write_read_records3};
use std::collections::HashMap;
use std::process::exit;
use rust_htslib::bam::record::CigarString;
use bio::io::fasta;
use std::io;
use matrix::PileupMatrix;
use pileup2matrix::generate_pileup_matrix;

fn main() {
    let bam_path = "wtc11_ont_grch38.chr22.bam";
    // let region = "chr22:30425877-30425912"; // 1-based
    // let region = "chr22:30425831-30425912";
    let region = "chr22:37063258-37063422";
    let ref_path = "GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.chr22.fna";
    let out_path = "new.bam";
    let out_path2 = "new2.bam";
    let out_path3 = "new3.bam";
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


    let bam_reader = BamReader::new(bam_path.to_string(), region.to_string());
    let mut read_records: HashMap<String, Vec<bam::Record>> = HashMap::new();

    let header = bam_reader.get_bam_header();
    bam_reader.get_read_records(&mut read_records);

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
    write_read_records1(&read_records, &header, String::from(out_path));
    write_read_records2(&mut read_records, &header, String::from(out_path2));
    write_read_records3(&mut read_records, &header, String::from(out_path3));

    let mut reader = fasta::Reader::from_file(ref_path).unwrap();
    for record in reader.records() {
        let record = record.unwrap();
        println!("{}", record.id());
        // println!("{}", std::str::from_utf8(record.seq()).unwrap());
        println!("{}", record.seq().len());
    }

    let mut matrices_vec: Vec<PileupMatrix> = Vec::new();
    generate_pileup_matrix(&bam_path.to_string(), &ref_path.to_string(), &region.to_string(), &mut matrices_vec);
    println!("matrices_vec.len(): {:?}", matrices_vec.len());
    for (readname, seq) in matrices_vec[0].base_matrix.iter() {
        println!("readname: {:?}", readname);
        println!("seq: {:?}", String::from_utf8(seq.clone()).unwrap());
    }
    // let mut v: Vec<_> = map.into_iter().collect();
    // v.sort_by(|x,y| x.0.cmp(&y.0));
    let mut sorted: Vec<_> = matrices_vec[0].positions.clone().into_iter().collect();
    sorted.sort_by(|a, b| a.1.cmp(&b.1));
    println!("max_idx: {:?}", matrices_vec[0].max_idx);
    println!("current pos: {:?}", matrices_vec[0].current_pos);
    for (ref_pos, idx) in sorted.iter() {
        print!("idx: {:?}\t", idx);
        println!("ref_pos: {:?}", ref_pos);
    }


}
