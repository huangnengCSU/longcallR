use crate::bam_reader::Region;
use rust_htslib::{bam, bam::Read};

pub fn find_isolated_regions(bam_path: &str, min_depth: u32, chr: Option<&str>) -> Vec<Region> {
    // Output the region of each isolated block: 1-based, left-closed, right-open
    let mut bam: bam::IndexedReader = bam::IndexedReader::from_path(bam_path).unwrap();
    let header = bam.header().clone();
    if chr.is_some() {
        bam.fetch(chr.unwrap()).unwrap();
    }
    let mut isolated_regions: Vec<Region> = Vec::new();
    let mut pre_tid: i32 = -1;
    let mut start_pos: u32 = 0;
    let mut end_pos: u32 = 0;
    let mut max_depth: u32 = 0;
    let mut pileups = bam.pileup();
    pileups.set_max_depth(1000000);

    for p in pileups {
        let pileup = p.unwrap();
        let pos = pileup.pos(); // 0-based
        let tid = pileup.tid();
        let depth = pileup.depth();
        // println!("{} {} {}", pos, tid, depth);

        if pre_tid == -1 {
            pre_tid = tid as i32;
            start_pos = pos;
            end_pos = pos;
            max_depth = depth;
            continue;
        }

        if tid as i32 != pre_tid && pre_tid != -1 {
            if end_pos > start_pos && max_depth >= min_depth {
                let region = Region {
                    chr: std::str::from_utf8(&header.tid2name(pre_tid as u32))
                        .unwrap()
                        .to_string(),
                    start: start_pos + 1,
                    end: end_pos + 2,
                };
                isolated_regions.push(region);
            }
            start_pos = pos;
            end_pos = pos;
            pre_tid = tid as i32;
            max_depth = depth;
            continue;
        }

        if pos != end_pos + 1 {
            if end_pos > start_pos && max_depth >= min_depth {
                let region = Region {
                    chr: std::str::from_utf8(&header.tid2name(pre_tid as u32))
                        .unwrap()
                        .to_string(),
                    start: start_pos + 1,
                    end: end_pos + 2,
                };
                isolated_regions.push(region);
            }
            start_pos = pos;
            end_pos = pos;
            pre_tid = tid as i32;
            max_depth = depth;
        } else {
            end_pos = pos;
            if depth > max_depth {
                max_depth = depth;
            }
        }
    }

    if end_pos > start_pos && max_depth >= min_depth {
        let region = Region {
            chr: std::str::from_utf8(&header.tid2name(pre_tid as u32))
                .unwrap()
                .to_string(),
            start: start_pos + 1,
            end: end_pos + 2,
        };
        isolated_regions.push(region);
    }
    return isolated_regions;
}
