mod bam_reader;


use rust_htslib::{bam, bam::Read};
use bam_reader::{BamReader};
use bam_reader::{write_read_records1, write_read_records2, write_read_records3};
use std::collections::HashMap;
use std::process::exit;
use rust_htslib::bam::record::CigarString;

fn main() {
    let bam_path = "/mnt/f/postdoc/LR-RNA-seq/wtc11_ont_grch38.chr22.bam";
    let region = "chr22:18994296-18994359";
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
}
