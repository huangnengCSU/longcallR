mod bam_pileup;

use bam_pileup::read_pileup;

fn main() {
    let bam_path = "/mnt/f/postdoc/LR-RNA-seq/wtc11_ont_grch38.chr22.bam";
    let chr = "chr22";
    let start_pos = 18994296;
    let end_pos = 18994359;
    read_pileup(bam_path, chr, start_pos, end_pos);
    println!("Hello, world!");
}
