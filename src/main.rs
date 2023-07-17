mod bam_reader;


use bam_reader::BamReader;

fn main() {
    let bam_path = "/mnt/f/postdoc/LR-RNA-seq/wtc11_ont_grch38.chr22.bam";
    let region = "chr22:18994296-18994359";
    let bam_reader = BamReader::new(bam_path.to_string(), region.to_string());
    // bam_reader.read_pileup();
    let records = bam_reader.get_read_records();
}
