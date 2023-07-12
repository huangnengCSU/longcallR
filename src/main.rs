use rust_htslib::{bam, bam::Read};

fn main() {
    let mut bam: bam::IndexedReader =
        bam::IndexedReader::from_path(&"/mnt/f/postdoc/LR-RNA-seq/wtc11_ont_grch38.chr22.bam")
            .unwrap();
    // let header = bam::Header::from_template(bam.header());
    bam.fetch(("chr22", 37051506, 37052135)).unwrap(); // set region
    for p in bam.pileup() {
        let pileup = p.unwrap();
        let pos = pileup.pos();
        let depth = pileup.depth();
        println!("{} {}", pos, depth);

        for alignment in pileup.alignments() {
            let q_pos = alignment.qpos();
            let record = alignment.record();
            let qname = String::from_utf8_lossy(record.qname());
            let seq = record.seq();
            if !alignment.is_del() {
                println!(
                    "reads: {:?}, pos: {:?}, base: {:?}",
                    qname,
                    pos,
                    seq[q_pos.unwrap()] as char
                );
            } else if (alignment.is_refskip()) {
                println!("reads: {:?}, pos: {:?}, base: {:?}", qname, pos, 'N');
            }
            match alignment.indel() {
                bam::pileup::Indel::Ins(len) => println!(
                    "Insertion of length {} between this and next position.",
                    len
                ),
                bam::pileup::Indel::Del(len) => {
                    println!("Deletion of length {} between this and next position.", len)
                }
                bam::pileup::Indel::None => (),
            }
        }
    }
    println!("Hello, world!");
}
