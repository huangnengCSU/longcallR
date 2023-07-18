use std::collections::HashMap;
use std::process::exit;
use rust_htslib::{bam, bam::Read};
use rust_htslib::bam::Format;
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::record::CigarString;

#[derive(Clone)]
pub struct Region {
    pub(crate) chr: String,
    pub(crate) start: u32,
    // 1-based
    pub(crate) end: u32,   // 1-based
}

impl Region {
    pub fn new(region: String) -> Region {
        // region format: chr:start-end
        if !region.contains(":") {
            let chr = region;
            return Region {
                chr,
                start: 0,
                end: 0,
            };
        } else if region.contains(":") && region.contains("-") {
            let region_vec: Vec<&str> = region.split(":").collect();
            let chr = region_vec[0].to_string();
            let pos_vec: Vec<&str> = region_vec[1].split("-").collect();
            let start = pos_vec[0].parse::<u32>().unwrap();
            let end = pos_vec[1].parse::<u32>().unwrap();
            assert!(start <= end);
            return Region {
                chr,
                start,
                end,
            };
        } else {
            panic!("region format error!");
        }
    }
}


pub struct BamReader {
    bam_path: String,
    region: Region,
}

impl BamReader {
    pub fn new(bam_path: String, region: String) -> BamReader {
        let region = Region::new(region);
        BamReader {
            bam_path,
            region,
        }
    }

    pub fn read_pileup(&self) {
        let mut bam: bam::IndexedReader = bam::IndexedReader::from_path(self.bam_path.clone()).unwrap();
        if self.region.start == 0 && self.region.end == 0 {
            bam.fetch(self.region.chr.as_str()).unwrap(); // set region
        } else {
            bam.fetch((self.region.chr.as_str(), self.region.start, self.region.end)).unwrap(); // set region
        }
        for p in bam.pileup() {
            let pileup = p.unwrap();
            let pos = pileup.pos();
            if pos < self.region.start || pos > self.region.end {
                continue;
            }
            let depth = pileup.depth();
            println!("depth of pos {} = {}", pos, depth);

            for alignment in pileup.alignments() {
                let q_pos = alignment.qpos();
                let record = alignment.record();
                let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
                let seq = record.seq();
                if !alignment.is_del() {
                    println!(
                        "readname: {:?}, pos: {:?}, base: {:?}",
                        qname,
                        pos,
                        seq[q_pos.unwrap()] as char
                    );
                } else if alignment.is_refskip() {
                    println!("readname: {:?}, pos: {:?}, base: {:?}", qname, pos, 'N');
                } else if alignment.is_del() {
                    println!("readname: {:?}, pos: {:?}, base: {:?}", qname, pos, 'D');
                }

                match alignment.indel() {
                    bam::pileup::Indel::Ins(len) => {
                        let mut insert_segment = String::from("");
                        for tmpi in 1..=len {
                            insert_segment.push(seq[q_pos.unwrap() + tmpi as usize] as char);
                        }
                        println!("Insertion of length {} ({}) between this and next position.", len, insert_segment);
                    }
                    bam::pileup::Indel::Del(len) => {
                        println!("readname: {:?}, pos: {:?}, base: {:?}", qname, pos, 'D');
                    }
                    bam::pileup::Indel::None => {}
                }
            }
        }
    }

    pub fn get_read_records(&self, read_records: &mut HashMap<String, Vec<bam::Record>>) {
        let mut bam: bam::IndexedReader = bam::IndexedReader::from_path(self.bam_path.clone()).unwrap();
        if self.region.start == 0 && self.region.end == 0 {
            bam.fetch(self.region.chr.as_str()).unwrap(); // set region
        } else {
            bam.fetch((self.region.chr.as_str(), self.region.start, self.region.end)).unwrap(); // set region
        }
        for r in bam.records() {
            let record = r.unwrap();
            let qname = std::str::from_utf8(record.qname()).unwrap();
            if read_records.get(qname) != None {
                read_records.get_mut(qname).unwrap().push(record);
            } else {
                read_records.insert(qname.to_string(), vec![record]);
            }
        }
        for (rname, record_vec) in read_records.iter_mut() {
            println!("{}:", rname);
            for record in record_vec.iter() {
                for cg in record.cigar().iter() {
                    print!("{}{}", cg.len(), cg.char());
                }
                println!("");
            }
        }
    }

    pub fn get_bam_header(&self) -> bam::header::Header {
        let bam: bam::IndexedReader = bam::IndexedReader::from_path(self.bam_path.clone()).unwrap();
        bam::Header::from_template(bam.header())
    }
}

// write read records to a new bam file
pub fn write_read_records1(read_records: &HashMap<String, Vec<bam::Record>>, bam_header: &bam::header::Header, bam_path: String) {
    let mut bam_writer = bam::Writer::from_path(&bam_path, &bam_header, Format::Bam).unwrap();
    for (rname, record_vec) in read_records.iter() {
        for r in record_vec.iter() {
            let re = bam_writer.write(r);
            if re == Ok(()) {
                println!("write success");
            } else {
                println!("write failed");
            }
        }
    }
}

pub fn write_read_records2(read_records: &mut HashMap<String, Vec<bam::Record>>, bam_header: &bam::header::Header, bam_path: String) {
    let mut bam_writer = bam::Writer::from_path(&bam_path, &bam_header, Format::Bam).unwrap();
    for (rname, record_vec) in read_records.iter_mut() {
        for record in record_vec.iter_mut() {
            let mut out_record = bam::Record::from(record.clone());
            out_record.set_pos(123456);
            // let mut new_cigar = CigarString(vec![Cigar::Match(1), Cigar::Ins(1), Cigar::Match(1)]);
            // Warning: cigar length should be the same as the read length
            // out_record.set(record.qname(),
            //                Some(&new_cigar),
            //                &record.seq().as_bytes(),
            //                record.qual());
            let origin_cigar = record.cigar().take();
            out_record.set(record.qname(),
                           Some(&origin_cigar),
                           &record.seq().as_bytes(),
                           record.qual());
            println!("rname: {:?}, pos: {:?}", rname, record.pos());
            let re = bam_writer.write(&out_record);
            if re == Ok(()) {
                println!("write success");
            } else {
                println!("write failed");
                exit(1);
            }
        }
    }
}

pub fn write_read_records3(read_records: &mut HashMap<String, Vec<bam::Record>>, bam_header: &bam::header::Header, bam_path: String) {
    let mut bam_writer = bam::Writer::from_path(&bam_path, &bam_header, Format::Bam).unwrap();
    for (rname, record_vec) in read_records.iter_mut() {
        for record in record_vec.iter_mut() {
            let mut out_record = bam::Record::from(record.clone());
            out_record.set_pos(123456);
            let mut new_cigar = CigarString(vec![Cigar::Match(1), Cigar::Ins(1), Cigar::Match(1)]);
            // Warning: cigar length should be the same as the read length
            out_record.set(record.qname(),
                           Some(&new_cigar),
                           &record.seq().as_bytes(),
                           record.qual());
            println!("rname: {:?}, pos: {:?}", rname, record.pos());
            let re = bam_writer.write(&out_record);
            if re == Ok(()) {
                println!("write success");
            } else {
                println!("write failed");
                exit(1);
            }
        }
    }
}