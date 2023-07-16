use rust_htslib::{bam, bam::Read};

fn main() {
    let bam_path = "/mnt/f/postdoc/LR-RNA-seq/wtc11_ont_grch38.chr22.bam";
    let mut bam: bam::IndexedReader = bam::IndexedReader::from_path(bam_path).unwrap();
    // let header = bam::Header::from_template(bam.header());
    let chr = "chr22";
    let start_pos = 18994296;
    let end_pos = 18994359;
    bam.fetch((chr, start_pos, end_pos)).unwrap(); // set region
    for p in bam.pileup() {
        let pileup = p.unwrap();
        let pos = pileup.pos();
        if pos < start_pos || pos > end_pos {
            continue;
        }
        let depth = pileup.depth();
        println!("depth of pos {} = {}", pos, depth);

        for alignment in pileup.alignments() {
            let q_pos = alignment.qpos();
            let record = alignment.record();
            let qname = String::from_utf8_lossy(record.qname());
            let seq = record.seq();
            if !alignment.is_del() {
                println!(
                    "read: {:?}, pos: {:?}, base: {:?}",
                    qname,
                    pos,
                    seq[q_pos.unwrap()] as char
                );
            } else if alignment.is_refskip() {
                println!("read: {:?}, pos: {:?}, base: {:?}", qname, pos, 'N');
            } else if alignment.is_del() {
                println!("read: {:?}, pos: {:?}, base: {:?}", qname, pos, 'D');
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
                    println!("read: {:?}, pos: {:?}, base: {:?}", qname, pos, 'D');
                }
                bam::pileup::Indel::None => {}
            }
        }
    }
    println!("Hello, world!");
}
