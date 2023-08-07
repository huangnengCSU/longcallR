use rust_htslib::{bam, bam::Read};
use crate::bam_reader::Region;

pub fn find_isolated_regions(bam_path: &str, min_depth: u32) -> Vec<Region> {
    let mut bam: bam::IndexedReader = bam::IndexedReader::from_path(bam_path).unwrap();
    let header = bam.header().clone();
    let mut isolated_regions: Vec<Region> = Vec::new();
    let mut pre_tid: u32 = 0;
    let mut start_pos: u32 = 0;
    let mut end_pos: u32 = 0;
    let mut max_depth: u32 = 0;

    for p in bam.pileup() {
        let pileup = p.unwrap();
        let pos = pileup.pos(); // 0-based
        let tid = pileup.tid();
        let depth = pileup.depth();
        // println!("{} {} {}", pos, tid, depth);
        if pre_tid != 0 && tid != pre_tid {
            if start_pos != 0 && end_pos != 0 {
                let region = Region {
                    chr: std::str::from_utf8(&header.tid2name(pre_tid as u32)).unwrap().to_string(),
                    start: start_pos,
                    end: end_pos + 1,
                };
                if max_depth >= min_depth {
                    isolated_regions.push(region);
                }
                start_pos = 0;
                end_pos = 0;
                pre_tid = 0;
            }
        }

        if pre_tid == 0 && start_pos == 0 && end_pos == 0 {
            pre_tid = tid;
            start_pos = pos;
            end_pos = pos;
            max_depth = depth;
        } else if pos != end_pos + 1 && pre_tid != 0 && start_pos != 0 && end_pos != 0 {
            let region = Region {
                chr: std::str::from_utf8(&header.tid2name(pre_tid as u32)).unwrap().to_string(),
                start: start_pos,
                end: end_pos + 1,
            };
            if max_depth >= min_depth {
                isolated_regions.push(region);
            }
            start_pos = 0;
            end_pos = 0;
            pre_tid = 0;
        } else {
            end_pos = pos;
            if depth > max_depth {
                max_depth = depth;
            }
        }
    }
    return isolated_regions;
}