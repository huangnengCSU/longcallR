use bio::bio_types::strand::ReqStrand::Forward;
use bio::io::fasta;
use std::collections::{HashMap, VecDeque};
use std::io::{BufRead, BufReader};
use std::sync::Mutex;
use std::{fs, fs::File};
// use fishers_exact::fishers_exact;
// use mathru::statistics::test::{ChiSquare, Test};
use rayon::prelude::*;
use rust_htslib::bam::record::Aux;
use rust_htslib::{
    bam,
    bam::{ext::BamRecordExtensions, Read},
};
use rust_lapper::{Interval, Lapper};
// use seq_io::fasta::Record;

use crate::constants::MAX_BASE_QUALITY;
use crate::Platform;

#[derive(Default, Clone, Debug)]
pub struct Region {
    pub(crate) chr: String,
    pub(crate) start: u32,
    // 1-based, inclusive
    pub(crate) end: u32,
    // 1-based, exclusive
    pub(crate) gene_id: Option<String>,
    // if load annotation, this field will tell which gene this region covers. Multiple gene separated by comma
}

impl Region {
    pub fn new(region: String) -> Region {
        // region format: chr:start-end
        if !region.contains(":") {
            let chr = region;
            Region {
                chr,
                start: 0,
                end: 0,
                gene_id: None,
            }
        } else if region.contains(":") && region.contains("-") {
            let region_vec: Vec<&str> = region.split(":").collect();
            let chr = region_vec[0].to_string();
            let pos_vec: Vec<&str> = region_vec[1].split("-").collect();
            let start = pos_vec[0].parse::<u32>().unwrap();
            let end = pos_vec[1].parse::<u32>().unwrap();
            let gene_id = None;
            assert!(start <= end);
            return Region {
                chr,
                start,
                end,
                gene_id,
            };
        } else {
            panic!("region format error!");
        }
    }
    pub fn to_string(&self) -> String {
        format!("{}:{}-{}", self.chr, self.start, self.end)
    }
}

#[derive(Default, Debug, Clone)]
pub struct BaseQual {
    pub a: Vec<u8>,
    pub c: Vec<u8>,
    pub g: Vec<u8>,
    pub t: Vec<u8>,
}

#[derive(Default, Debug, Clone)]
pub struct BaseStrands {
    pub a: [i32; 2],
    // [forward, backward]
    pub c: [i32; 2],
    pub g: [i32; 2],
    pub t: [i32; 2],
}

#[derive(Default, Debug, Clone)]
pub struct DistanceToEnd {
    pub a: Vec<i64>,
    // allele A to the end of read for every read
    pub c: Vec<i64>,
    // allele C to the end of read for every read
    pub g: Vec<i64>,
    // allele G to the end of read for every read
    pub t: Vec<i64>,
    // allele T to the end of read for every read
}

#[derive(Default, Debug, Clone)]
pub struct BaseFreq {
    pub a: u32,
    pub c: u32,
    pub g: u32,
    pub t: u32,
    pub n: u32,
    // number of introns
    pub d: u32,
    // number of deletions
    pub ni: u32,
    // number of insertions
    pub i: bool,
    // whether this position falls in an insertion
    pub ref_base: char,
    // reference base, A,C,G,T,N,-
    // pub intron: bool,
    // // whether this position is an intron
    pub forward_cnt: u32,
    // number of forward reads covering this position, excluding intron
    pub backward_cnt: u32,
    // number of backward reads covering this position, excluding intron
    pub baseq: BaseQual,
    pub base_strands: BaseStrands,
    pub transcript_strands: [i32; 2],
    // [forward, backward]
    pub distance_to_end: DistanceToEnd,
}

impl BaseFreq {
    pub fn get_allele_baseq(&self, allele: char) -> Vec<u8> {
        match allele {
            'A' | 'a' => self.baseq.a.clone(),
            'C' | 'c' => self.baseq.c.clone(),
            'G' | 'g' => self.baseq.g.clone(),
            'T' | 't' => self.baseq.t.clone(),
            _ => {
                panic!("Invalid allele: {}", allele);
            }
        }
    }

    pub fn get_allele_base_strands(&self, allele: char) -> [i32; 2] {
        match allele {
            'A' | 'a' => self.base_strands.a.clone(),
            'C' | 'c' => self.base_strands.c.clone(),
            'G' | 'g' => self.base_strands.g.clone(),
            'T' | 't' => self.base_strands.t.clone(),
            _ => {
                panic!("Invalid allele: {}", allele);
            }
        }
    }

    pub fn get_depth_include_intron(&self) -> u32 {
        self.a + self.c + self.g + self.t + self.d + self.n
    }

    pub fn get_allele_counts(&self) -> u32 {
        self.a + self.c + self.g + self.t
    }

    pub fn get_two_major_alleles(&self, ref_base: char) -> (char, u32, char, u32) {
        let mut x = vec![('A', self.a), ('C', self.c), ('G', self.g), ('T', self.t)];
        x.sort_by(|a, b| b.1.cmp(&a.1));
        if x[0].0 != ref_base && x[1].0 != ref_base {
            if x[2].1 == x[1].1 && x[2].0 == ref_base {
                (x[0].0, x[0].1, x[2].0, x[2].1)
            } else if x[3].1 == x[1].1 && x[3].0 == ref_base {
                (x[0].0, x[0].1, x[3].0, x[3].1)
            } else {
                (x[0].0, x[0].1, x[1].0, x[1].1)
            }
        } else {
            (x[0].0, x[0].1, x[1].0, x[1].1)
        }
    }

    pub fn get_main_alleles(&self) -> Vec<(char, u32)> {
        let mut allele_counts = vec![('A', self.a), ('C', self.c), ('G', self.g), ('T', self.t)];
        allele_counts.sort_by(|a, b| b.1.cmp(&a.1));
        let mut result = Vec::new();
        if allele_counts[0].1 == 0 {
            return result; // No valid alleles if all counts are zero
        }
        result.push(allele_counts[0].clone());
        for item in allele_counts.iter().skip(1) {
            if item.1 == 0 {
                break;
            }
            if item.1 == result.last().unwrap().1 {
                result.push(item.clone());
            } else {
                if result.len() == 1 {
                    result.push(item.clone());
                } else {
                    break;
                }
            }
        }
        result
    }

    pub fn get_none_ref_count(&self) -> u32 {
        match self.ref_base {
            'A' => self.c + self.g + self.t + self.d,
            'C' => self.a + self.g + self.t + self.d,
            'G' => self.a + self.c + self.t + self.d,
            'T' => self.a + self.c + self.g + self.d,
            _ => 0,
        }
    }
}

pub fn load_reference(ref_path: &str) -> HashMap<String, Vec<u8>> {
    let mut ref_seqs: HashMap<String, Vec<u8>> = HashMap::new();
    let reader = fasta::Reader::from_file(ref_path).unwrap();
    for r in reader.records() {
        let ref_record = r.unwrap();
        ref_seqs.insert(ref_record.id().to_string(), ref_record.seq().to_vec());
    }
    ref_seqs
}

pub fn parse_fai(fai_path: &str) -> Vec<(String, u32)> {
    let mut contig_lengths: Vec<(String, u32)> = Vec::new();
    let file = File::open(fai_path).unwrap();
    let reader = BufReader::new(file);
    for r in reader.lines() {
        let line = r.unwrap().clone();
        let parts: Vec<&str> = line.split('\t').collect();
        contig_lengths.push((parts[0].to_string(), parts[1].parse().unwrap()));
    }
    contig_lengths
}

pub fn find_isolated_regions_with_depth(
    bam_path: &str,
    chr: &str,
    ref_len: u32,
    min_mapq: u8,
    min_read_length: usize,
    divergence: f32,
) -> Vec<Region> {
    let mut isolated_regions: Vec<Region> = Vec::new();
    let mut depth_vec: Vec<u32> = vec![0; ref_len as usize];
    let mut bam: bam::IndexedReader = bam::IndexedReader::from_path(bam_path).unwrap();
    // let header = bam.header().clone();
    match bam.fetch(chr) {
        Ok(_) => {}
        Err(e) => {
            eprintln!(
                "Error fetching chromosome {}: {}. Skipping this chromosome.",
                chr, e
            );
            return isolated_regions;
        }
    }
    // bam.fetch(chr).unwrap();
    for r in bam.records() {
        let record = r.unwrap();
        if record.mapq() < min_mapq
            || record.seq_len() < min_read_length
            || record.is_unmapped()
            || record.is_secondary()
            || record.is_supplementary()
        {
            continue;
        }

        match record.aux(b"de") {
            Ok(Aux::Float(f)) => {
                if f >= divergence {
                    continue;
                }
            }
            _ => {}
        }

        let ref_start = record.reference_start(); // 0-based, left-closed
        let ref_end = record.reference_end(); // 0-based, right-open
        for i in ref_start..ref_end {
            depth_vec[i as usize] += 1;
        }
    }
    let mut region_start = -1;
    let mut region_end = -1;
    for i in 0..ref_len {
        if depth_vec[i as usize] == 0 {
            if region_end > region_start {
                assert!(region_start >= 0);
                assert!(region_end >= 0);
                isolated_regions.push(Region {
                    chr: chr.to_string(),
                    start: (region_start + 1) as u32,
                    end: (region_end + 2) as u32,
                    gene_id: None,
                });
                region_start = -1;
                region_end = -1;
            }
        } else {
            if region_start == -1 {
                region_start = i as i32;
                region_end = i as i32;
            } else {
                region_end = i as i32;
            }
        }
    }
    if region_end > region_start {
        isolated_regions.push(Region {
            chr: chr.to_string(),
            start: (region_start + 1) as u32,
            end: (region_end + 2) as u32,
            gene_id: None,
        });
        // region_start = -1;
        // region_end = -1;
    }
    isolated_regions
}

pub fn parse_annotation(
    anno_path: String,
) -> (
    HashMap<String, VecDeque<Region>>,
    HashMap<String, Vec<Interval<u32, u8>>>,
) {
    let mut gene_regions: HashMap<String, VecDeque<Region>> = HashMap::new(); // key is chr, value is a stack of gene regions
    let mut exon_regions: HashMap<String, Vec<Interval<u32, u8>>> = HashMap::new(); // key is gene id, value is exon regions
    let file = File::open(anno_path).unwrap();
    let reader = BufReader::new(file);
    let mut invs: Vec<Interval<u32, u8>> = Vec::new(); // for merging gene exons
    let mut gene_id: String = String::new();
    for line in reader.lines() {
        let line = line.unwrap();
        if line.starts_with("#") {
            continue;
        }
        let parts = line.split('\t').collect::<Vec<&str>>();
        let seqname = parts[0].to_string();
        let feature = parts[2];
        let start = parts[3].parse::<u32>().unwrap(); // 1-based, inclusive
        let end = parts[4].parse::<u32>().unwrap(); // 1-based, inclusive
        if feature == "gene" {
            if invs.len() > 0 {
                exon_regions.insert(gene_id.clone(), invs.clone());
                invs.clear();
            }
            if !gene_regions.contains_key(&seqname) {
                gene_regions.insert(seqname.clone(), VecDeque::new());
            }
            for subpart in parts[8].trim_end().split(";").collect::<Vec<&str>>() {
                if subpart.starts_with("gene_id") {
                    gene_id = subpart.replace("gene_id=", ""); // gff3 format
                    break;
                }
                if subpart.starts_with("gene_id") {
                    gene_id = subpart.replace("gene_id ", "").replace("\"", ""); // gtf format
                    break;
                }
            }
            let regions = gene_regions.get_mut(&seqname).unwrap();
            if let Some(mut top) = regions.pop_back() {
                assert!(
                    start >= top.start,
                    "Error: annotation file is not sorted. {}:{}-{}",
                    seqname,
                    start,
                    end
                );
                // if overlap, merge the overlapped regions
                if top.end <= start {
                    // end of top region is exclusive, so top.end == start is not overlap
                    regions.push_back(top);
                    regions.push_back(Region {
                        chr: seqname.clone(),
                        start,
                        end: end + 1,
                        gene_id: Some(gene_id.clone()),
                    });
                } else if top.end < end + 1 {
                    // top.end is exclusive, end is inclusive
                    // merge two overlapped regions
                    top.end = end + 1;
                    top.gene_id = Some(top.gene_id.unwrap() + "," + &gene_id);
                    regions.push_back(top);
                } else {
                    // equal end or contained
                    top.gene_id = Some(top.gene_id.unwrap() + "," + &gene_id);
                    regions.push_back(top);
                }
            } else {
                // first gene region in stack
                regions.push_back(Region {
                    chr: seqname.clone(),
                    start,
                    end: end + 1,
                    gene_id: Option::from(gene_id.clone()),
                });
            }
        } else if feature == "CDS" {
            let mut exon_gene_id = String::new();
            for subpart in parts[8].trim_end().split(";").collect::<Vec<&str>>() {
                if subpart.starts_with("gene_id") {
                    exon_gene_id = subpart.replace("gene_id=", ""); // gff3 format
                    break;
                }
                if subpart.starts_with("gene_id") {
                    exon_gene_id = subpart.replace("gene_id ", "").replace("\"", ""); // gtf format
                    break;
                }
            }
            assert_eq!(
                exon_gene_id, gene_id,
                "Error: gene_id in gene and exon are different: gene_id:{}, exon_gene_id:{}",
                gene_id, exon_gene_id
            );
            invs.push(Interval {
                start,
                stop: end + 1,
                val: 0,
            }); // 1-based,
        } else {
            continue;
        }
    }
    if invs.len() > 0 {
        exon_regions.insert(gene_id.clone(), invs.clone());
        invs.clear();
    }
    (gene_regions, exon_regions)
}

pub fn lapper_intervals(
    query_regions: &Vec<Region>,
    target_regions: &VecDeque<Region>,
) -> Vec<Region> {
    let mut result_regions = Vec::new();
    let mut invs: Vec<Interval<usize, String>> = Vec::new();
    for region in target_regions {
        invs.push(Interval {
            start: region.start as usize,
            stop: region.end as usize,
            val: region.gene_id.clone().unwrap(),
        });
    }
    let interval_tree = Lapper::new(invs);
    for q in query_regions {
        let q_inv = Interval {
            start: q.start as usize,
            stop: q.end as usize,
            val: 0,
        };
        for h_inv in interval_tree.find(q.start as usize, q.end as usize) {
            let intersected_start = q_inv.start.max(h_inv.start);
            let intersected_end = q_inv.stop.min(h_inv.stop);
            let h_gene = h_inv.val.clone();
            assert!(
                intersected_start < intersected_end,
                "Error: intersected_start >= intersected_end, query:{:?}",
                q_inv
            );
            result_regions.push(Region {
                chr: q.chr.clone(),
                start: intersected_start as u32,
                end: intersected_end as u32,
                gene_id: Option::from(h_gene),
            });
        }
    }
    result_regions
}

pub fn intersect_gene_regions(
    alignment_regions: &Vec<Region>,
    gene_regions: &HashMap<String, VecDeque<Region>>,
    thread_size: usize,
) -> Vec<Region> {
    let intersected_regions: Mutex<Vec<Region>> = Mutex::new(Vec::new());
    let mut region_map = HashMap::new();
    for region in alignment_regions {
        if !region_map.contains_key(&region.chr) {
            region_map.insert(region.chr.clone(), Vec::new());
        }
        region_map
            .get_mut(&region.chr)
            .unwrap()
            .push(region.clone());
    }
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(thread_size - 1)
        .build()
        .unwrap();
    let contig_names = region_map.keys().collect::<Vec<&String>>();
    pool.install(|| {
        contig_names.par_iter().for_each(|ctg| {
            if gene_regions.contains_key(*ctg) {
                let chr_intersected_region = lapper_intervals(
                    region_map.get(*ctg).unwrap(),
                    gene_regions.get(*ctg).unwrap(),
                );
                for region in chr_intersected_region {
                    intersected_regions.lock().unwrap().push(region);
                }
            }
        });
    });

    intersected_regions.into_inner().unwrap()
}

pub fn extract_isolated_regions_parallel(
    bam_file: &str,
    ref_file: &str,
    thread_size: usize,
    contigs: Option<Vec<String>>,
    min_mapq: u8,
    min_read_length: usize,
    divergence: f32,
) -> Vec<Region> {
    let results: Mutex<Vec<Region>> = Mutex::new(Vec::new());
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(thread_size - 1)
        .build()
        .unwrap();
    let fai_path = String::from(ref_file) + ".fai";
    if fs::metadata(&fai_path).is_err() {
        panic!("Reference index file .fai does not exist.");
    }
    let contig_lengths = parse_fai(fai_path.as_str());
    let contig_names: VecDeque<String> = match contigs {
        Some(names) => names.iter().cloned().collect(),
        None => contig_lengths.iter().map(|(chr, _)| chr.clone()).collect(),
    };
    pool.install(|| {
        contig_names.par_iter().for_each(|ctg| {
            let length = contig_lengths.iter().find(|(chr, _)| chr == ctg).unwrap().1;
            results
                .lock()
                .unwrap()
                .extend(find_isolated_regions_with_depth(
                    bam_file,
                    ctg,
                    length,
                    min_mapq,
                    min_read_length,
                    divergence,
                ));
        });
    });
    results.into_inner().unwrap()
}

// pub fn read_references(ref_path: &str) -> HashMap<String, Vec<u8>> {
//     let mut references: HashMap<String, Vec<u8>> = HashMap::new();
//     let mut reader = Reader::from_path(ref_path).unwrap();
//     while let Some(record) = reader.next() {
//         let record = record.expect("Error reading record");
//         references.insert(record.id().unwrap().to_string(), record.full_seq().to_vec());
//     }
//     references
// }

#[derive(Default, Debug, Clone)]
pub struct Profile {
    pub freq_vec: Vec<BaseFreq>,
    pub region: Region,
}

impl Profile {
    pub fn fill_data_into_freq_vec(
        &mut self,
        bam_path: &str,
        region: &Region,
        ref_seq: &Vec<u8>,
        platform: &Platform,
        min_mapq: u8,
        min_read_length: usize,
        divergence: f32,
        distance_to_read_end: u32,
        polya_tail_length: u32,
    ) {
        // When region is large and the number of reads is large, the runtime of init_profile_with_pileup is time-consuming.
        // This function is used to fill the profile by parsing each read in the bam file instead of using pileup.

        let mut bam: bam::IndexedReader = bam::IndexedReader::from_path(bam_path).unwrap();
        bam.fetch((region.chr.as_str(), region.start, region.end))
            .unwrap();
        let vec_size = (region.end - region.start) as usize; // end is exclusive
        self.freq_vec = vec![BaseFreq::default(); vec_size];
        self.region = region.clone();
        let freq_vec_start_pos = region.start as usize - 1; // the first position on reference, 0-based, inclusive
        let polya_tail_length = polya_tail_length as i64;

        // fill the ref_base field in each BaseFreq
        for i in 0..vec_size {
            self.freq_vec[i].ref_base = ref_seq[freq_vec_start_pos + i] as char;
        }

        for r in bam.records() {
            let record = r.unwrap();
            if record.mapq() < min_mapq
                || record.seq_len() < min_read_length
                || record.is_unmapped()
                || record.is_secondary()
                || record.is_supplementary()
            {
                continue;
            }

            match record.aux(b"de") {
                Ok(Aux::Float(f)) => {
                    if f >= divergence {
                        continue;
                    }
                }
                _ => {}
            }

            // let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
            let seq = record.seq();
            let base_qual = record.qual();
            let strand = if record.strand() == Forward { 0 } else { 1 };
            let mut ts = Aux::Char(b'*');
            match record.aux(b"ts") {
                Ok(value) => {
                    ts = value;
                }
                Err(_) => {}
            }
            let start_pos = record.pos() as usize; // 0-based
            let cigar = record.cigar();
            let leading_softclips = cigar.leading_softclips();
            let trailing_softclips = cigar.trailing_softclips();

            let mut pos_in_freq_vec: i32 = start_pos as i32 - freq_vec_start_pos as i32;
            let mut pos_in_read = if leading_softclips > 0 {
                leading_softclips as usize
            } else {
                0
            };
            let cigars = cigar.to_vec();
            for cg_idx in 0..cigars.len() {
                let cg = cigars[cg_idx];
                match cg.char() as u8 {
                    b'S' | b'H' => {
                        continue;
                    }
                    b'M' | b'X' | b'=' => {
                        for cgi in 0..cg.len() {
                            if pos_in_freq_vec < 0 {
                                // head of sequence which not in window region
                                pos_in_freq_vec += 1;
                                pos_in_read += 1;
                                continue;
                            }
                            if pos_in_freq_vec >= vec_size as i32 {
                                // tail of sequence which not in window region
                                break;
                            }
                            let base = seq[pos_in_read] as char;
                            let baseq = if base_qual[pos_in_read] < MAX_BASE_QUALITY {
                                base_qual[pos_in_read]
                            } else {
                                MAX_BASE_QUALITY
                            };

                            let ref_base = self.freq_vec[pos_in_freq_vec as usize].ref_base;
                            // filter alternate allele near the splicing site
                            // if base != ref_base {
                            //     if cgi <= 3
                            //         && cg_idx >= 1
                            //         && cigars[cg_idx - 1].char() as u8 == b'N'
                            //     {
                            //         pos_in_freq_vec += 1;
                            //         pos_in_read += 1;
                            //         continue;
                            //     } else if cgi >= cg.len() - 4
                            //         && cg_idx + 1 < cigars.len()
                            //         && cigars[cg_idx + 1].char() as u8 == b'N'
                            //     {
                            //         pos_in_freq_vec += 1;
                            //         pos_in_read += 1;
                            //         continue;
                            //     }
                            // }

                            let mut poly_a_flag = false;
                            let mut homopolymer_flag = false;
                            let mut trim_flag = false;
                            let distance_to_read_end = distance_to_read_end as i64;
                            let curr_pos = pos_in_read as i64;
                            let read_end_boundary = seq.len() as i64 - trailing_softclips;

                            // trim the end of ont reads due to accuracy
                            if matches!(platform, Platform::Ont) {
                                if (curr_pos - leading_softclips).abs() < distance_to_read_end
                                    || (curr_pos - read_end_boundary).abs() < distance_to_read_end
                                {
                                    trim_flag = true;
                                }
                            }

                            // skip poly_a, poly_t, poly_c, poly_g with length `polya_tail_length` within `distance_to_read_end` from read end
                            if !trim_flag {
                                if (curr_pos - leading_softclips).abs() < distance_to_read_end
                                    || (curr_pos - read_end_boundary).abs() < distance_to_read_end
                                {
                                    for tmpi in curr_pos - polya_tail_length..=(curr_pos + 1) {
                                        if tmpi < 0
                                            || tmpi + polya_tail_length - 1 >= seq.len() as i64
                                        {
                                            continue;
                                        }

                                        let mut poly_counts = [0; 4]; // A, T, C, G
                                        for tmpj in 0..polya_tail_length {
                                            let base = seq[(tmpi + tmpj) as usize];
                                            match base {
                                                b'A' if ref_base != 'A' => poly_counts[0] += 1,
                                                b'T' if ref_base != 'T' => poly_counts[1] += 1,
                                                b'C' if ref_base != 'C' => poly_counts[2] += 1,
                                                b'G' if ref_base != 'G' => poly_counts[3] += 1,
                                                _ => {}
                                            }
                                        }

                                        if poly_counts[0] >= polya_tail_length
                                            || poly_counts[1] >= polya_tail_length
                                        {
                                            poly_a_flag = true;
                                        }
                                        if poly_counts[2] >= polya_tail_length
                                            || poly_counts[3] >= polya_tail_length
                                        {
                                            homopolymer_flag = true;
                                        }
                                    }
                                }
                            }

                            if !trim_flag && !poly_a_flag && !homopolymer_flag {
                                // calculate distance to read end of each allele, for filtering variants that the average distance of each allele is significantly different
                                let pos = pos_in_read as i64;
                                let seq_len = seq.len() as i64;
                                let candidate1 = pos - leading_softclips;
                                let candidate2 = pos - (seq_len - trailing_softclips);
                                let dist = if candidate1.abs() < candidate2.abs() {
                                    candidate1
                                } else {
                                    candidate2
                                };

                                if strand == 0 {
                                    if ts == Aux::Char(b'+') {
                                        self.freq_vec[pos_in_freq_vec as usize]
                                            .transcript_strands[0] += 1; // read +, ts +, transcript +
                                    } else if ts == Aux::Char(b'-') {
                                        self.freq_vec[pos_in_freq_vec as usize]
                                            .transcript_strands[1] += 1; // read +, ts -, transcript -
                                    }
                                } else if strand == 1 {
                                    if ts == Aux::Char(b'+') {
                                        self.freq_vec[pos_in_freq_vec as usize]
                                            .transcript_strands[1] += 1; // read -, ts +, transcript -
                                    } else if ts == Aux::Char(b'-') {
                                        self.freq_vec[pos_in_freq_vec as usize]
                                            .transcript_strands[0] += 1; // read -, ts -, transcript +
                                    }
                                }

                                match base {
                                    'A' | 'a' => {
                                        self.freq_vec[pos_in_freq_vec as usize].a += 1;
                                        self.freq_vec[pos_in_freq_vec as usize].baseq.a.push(baseq);
                                        if strand == 0 {
                                            self.freq_vec[pos_in_freq_vec as usize]
                                                .base_strands
                                                .a[0] += 1;
                                        } else {
                                            self.freq_vec[pos_in_freq_vec as usize]
                                                .base_strands
                                                .a[1] += 1;
                                        }
                                        self.freq_vec[pos_in_freq_vec as usize]
                                            .distance_to_end
                                            .a
                                            .push(dist);
                                    }
                                    'C' | 'c' => {
                                        self.freq_vec[pos_in_freq_vec as usize].c += 1;
                                        self.freq_vec[pos_in_freq_vec as usize].baseq.c.push(baseq);
                                        if strand == 0 {
                                            self.freq_vec[pos_in_freq_vec as usize]
                                                .base_strands
                                                .c[0] += 1;
                                        } else {
                                            self.freq_vec[pos_in_freq_vec as usize]
                                                .base_strands
                                                .c[1] += 1;
                                        }
                                        self.freq_vec[pos_in_freq_vec as usize]
                                            .distance_to_end
                                            .c
                                            .push(dist);
                                    }
                                    'G' | 'g' => {
                                        self.freq_vec[pos_in_freq_vec as usize].g += 1;
                                        self.freq_vec[pos_in_freq_vec as usize].baseq.g.push(baseq);
                                        if strand == 0 {
                                            self.freq_vec[pos_in_freq_vec as usize]
                                                .base_strands
                                                .g[0] += 1;
                                        } else {
                                            self.freq_vec[pos_in_freq_vec as usize]
                                                .base_strands
                                                .g[1] += 1;
                                        }
                                        self.freq_vec[pos_in_freq_vec as usize]
                                            .distance_to_end
                                            .g
                                            .push(dist);
                                    }
                                    'T' | 't' => {
                                        self.freq_vec[pos_in_freq_vec as usize].t += 1;
                                        self.freq_vec[pos_in_freq_vec as usize].baseq.t.push(baseq);
                                        if strand == 0 {
                                            self.freq_vec[pos_in_freq_vec as usize]
                                                .base_strands
                                                .t[0] += 1;
                                        } else {
                                            self.freq_vec[pos_in_freq_vec as usize]
                                                .base_strands
                                                .t[1] += 1;
                                        }
                                        self.freq_vec[pos_in_freq_vec as usize]
                                            .distance_to_end
                                            .t
                                            .push(dist);
                                    }
                                    _ => {
                                        println!("Invalid nucleotide base: {}", base);
                                    }
                                }
                                if strand == 0 {
                                    self.freq_vec[pos_in_freq_vec as usize].forward_cnt += 1;
                                } else {
                                    self.freq_vec[pos_in_freq_vec as usize].backward_cnt += 1;
                                }
                            }

                            pos_in_freq_vec += 1;
                            pos_in_read += 1;
                        }
                    }
                    b'D' => {
                        for _ in 0..cg.len() {
                            if pos_in_freq_vec < 0 {
                                pos_in_freq_vec += 1;
                                continue;
                            }
                            if pos_in_freq_vec >= vec_size as i32 {
                                break;
                            }
                            self.freq_vec[pos_in_freq_vec as usize].d += 1;
                            pos_in_freq_vec += 1;
                        }
                    }
                    b'I' => {
                        if pos_in_freq_vec < 1 {
                            // smaller than 1 instead of 0, because insertion is counted as the previous position
                            pos_in_read += cg.len() as usize;
                            continue;
                        }
                        if pos_in_freq_vec >= vec_size as i32 {
                            break;
                        }
                        self.freq_vec[(pos_in_freq_vec - 1) as usize].ni += 1; // insertion is counted as the previous position
                        pos_in_read += cg.len() as usize;
                    }
                    b'N' => {
                        for _ in 0..cg.len() {
                            if pos_in_freq_vec < 0 {
                                pos_in_freq_vec += 1;
                                continue;
                            }
                            if pos_in_freq_vec >= vec_size as i32 {
                                break;
                            }
                            self.freq_vec[pos_in_freq_vec as usize].n += 1;
                            pos_in_freq_vec += 1;
                        }
                    }
                    _ => {
                        panic!("Error: unknown cigar operation: {}", cg.char());
                    }
                }
            }
        }
    }
}

// pub fn independent_test(data: [u32; 4]) -> f64 {
//     // use chi-square test or fisher's exact test to test whether two variables are independent
//     // N<40 or T<1, use fisher's exact test
//     // N>=40, T>=5 in >80% of cells, use chi-square test
//     // N>=40, T>=5 in <=80% of cells, use fisher's exact test
//     // return true if independent, false if not independent
//     let mut phred_pvalue = 0.0;
//     let total = data[0] + data[1] + data[2] + data[3];
//     let mut tl = [0.0; 4];
//     tl[0] = (((data[0] + data[1]) * (data[0] + data[2])) as f32) / (total as f32);
//     tl[1] = (((data[0] + data[1]) * (data[1] + data[3])) as f32) / (total as f32);
//     tl[2] = (((data[2] + data[3]) * (data[0] + data[2])) as f32) / (total as f32);
//     tl[3] = (((data[2] + data[3]) * (data[1] + data[3])) as f32) / (total as f32);
//     let mut t = 0.0;
//     for i in 0..data.len() {
//         if tl[i] >= 5.0 {
//             t += 1.0;
//         }
//     }
//     t = t / data.len() as f32;
//
//     if total > 40 && t > 0.8 {
//         // use chi-square test
//         let x = vec![data[0] as f64, data[1] as f64];
//         let y = vec![data[2] as f64, data[3] as f64];
//         let p = ChiSquare::test_vector(&x, &y).p_value();
//         phred_pvalue = -10.0 * p.log10();
//     }
//
//     if total <= 40 {
//         // use fisher's exact test
//         let p = fishers_exact(&data).unwrap().two_tail_pvalue;
//         phred_pvalue = -10.0 * p.log10();
//     }
//
//     phred_pvalue
// }
