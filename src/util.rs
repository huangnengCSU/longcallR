use std::{fs, fs::File};
use std::collections::{HashMap, VecDeque};
use std::io::{BufRead, BufReader};
use std::sync::Mutex;
use std::time::Instant;

use bio::bio_types::strand::ReqStrand::Forward;
use bio::io::fasta;
use fishers_exact::fishers_exact;
use mathru::statistics::test::{ChiSquare, Test};
use rayon::prelude::*;
use rust_htslib::{bam, bam::{ext::BamRecordExtensions, Read}};
use rust_lapper::{Interval, Lapper};
use seq_io::fasta::{Reader, Record};

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
            return Region {
                chr,
                start: 0,
                end: 0,
                gene_id: None,
            };
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
        return format!("{}:{}-{}", self.chr, self.start, self.end);
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
    pub intron: bool,
    // whether this position is an intron
    pub forward_cnt: u32,
    // number of forward reads covering this position, excluding intron
    pub backward_cnt: u32,
    // number of backward reads covering this position, excluding intron
    pub baseq: BaseQual,
    pub base_strands: BaseStrands,
    pub distance_to_end: DistanceToEnd,
}

impl BaseFreq {
    pub fn subtract(&mut self, base: u8) {
        match base {
            b'A' => {
                assert!(self.a > 0);
                self.a -= 1;
            }
            b'C' => {
                assert!(self.c > 0);
                self.c -= 1;
            }
            b'G' => {
                assert!(self.g > 0);
                self.g -= 1;
            }
            b'T' => {
                assert!(self.t > 0);
                self.t -= 1;
            }
            b'N' => {
                assert!(self.n > 0);
                self.n -= 1;
            }
            b'-' => {
                assert!(self.d > 0);
                self.d -= 1;
            }
            b'*' => {
                return;
            }
            _ => {
                panic!("Invalid base: {}", base as char);
            }
        }
    }

    pub fn add(&mut self, base: u8) {
        match base {
            b'A' => {
                self.a += 1;
            }
            b'C' => {
                self.c += 1;
            }
            b'G' => {
                self.g += 1;
            }
            b'T' => {
                self.t += 1;
            }
            b'N' => {
                self.n += 1;
            }
            b'-' => {
                self.d += 1;
            }
            b'*' => {
                return;
            }
            _ => {
                panic!("Invalid base: {}", base as char);
            }
        }
    }

    pub fn get_depth_include_intron(&self) -> u32 {
        self.a + self.c + self.g + self.t + self.d + self.n
    }

    pub fn get_depth_exclude_intron_deletion(&self) -> u32 {
        self.a + self.c + self.g + self.t
    }

    pub fn get_two_major_alleles(&self, ref_base: char) -> (char, u32, char, u32) {
        let mut x: Vec<(char, u32)> = [('A', self.a), ('C', self.c), ('G', self.g), ('T', self.t)].iter().cloned().collect();
        // sort by count: u32
        x.sort_by(|a, b| b.1.cmp(&a.1));
        if x[0].0 != ref_base && x[1].0 != ref_base {
            if x[2].0 == ref_base && x[1].1 == x[2].1 {
                return (x[0].0, x[0].1, x[2].0, x[2].1);
            } else if x[3].0 == ref_base && x[1].1 == x[3].1 {
                return (x[0].0, x[0].1, x[3].0, x[3].1);
            } else {
                return (x[0].0, x[0].1, x[1].0, x[1].1);
            }
        } else {
            return (x[0].0, x[0].1, x[1].0, x[1].1);
        }
    }

    pub fn get_none_ref_count(&self) -> u32 {
        match self.ref_base {
            'A' => self.c + self.g + self.t + self.d,
            'C' => self.a + self.g + self.t + self.d,
            'G' => self.a + self.c + self.t + self.d,
            'T' => self.a + self.c + self.g + self.d,
            _ => {
                0
            }
        }
    }
}


pub fn load_reference(ref_path: String) -> HashMap<String, Vec<u8>> {
    let mut ref_seqs: HashMap<String, Vec<u8>> = HashMap::new();
    let reader = fasta::Reader::from_file(ref_path).unwrap();
    for r in reader.records() {
        let ref_record = r.unwrap();
        ref_seqs.insert(ref_record.id().to_string(), ref_record.seq().to_vec());
    }
    return ref_seqs;
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
    return contig_lengths;
}

pub fn find_isolated_regions_with_depth(bam_path: &str, chr: &str, ref_len: u32, min_mapq: u8, min_read_length: usize) -> Vec<Region> {
    let mut isolated_regions: Vec<Region> = Vec::new();
    let mut depth_vec: Vec<u32> = vec![0; ref_len as usize];
    let mut bam: bam::IndexedReader = bam::IndexedReader::from_path(bam_path).unwrap();
    let header = bam.header().clone();
    bam.fetch(chr).unwrap();
    for r in bam.records() {
        let record = r.unwrap();
        if record.mapq() < min_mapq || record.seq_len() < min_read_length || record.is_unmapped() || record.is_secondary() || record.is_supplementary() {
            continue;
        }
        let ref_start = record.reference_start();   // 0-based, left-closed
        let ref_end = record.reference_end();   // 0-based, right-open
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
                isolated_regions.push(Region { chr: chr.to_string(), start: (region_start + 1) as u32, end: (region_end + 2) as u32, gene_id: None });
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
        isolated_regions.push(Region { chr: chr.to_string(), start: (region_start + 1) as u32, end: (region_end + 2) as u32, gene_id: None });
        region_start = -1;
        region_end = -1;
    }
    return isolated_regions;
}

pub fn parse_annotation(anno_path: String) -> (HashMap<String, VecDeque<Region>>, HashMap<String, Vec<Interval<usize, u8>>>) {
    let mut gene_regions: HashMap<String, VecDeque<Region>> = HashMap::new();    // key is chr, value is a stack of gene regions
    let mut exon_regions: HashMap<String, Vec<Interval<usize, u8>>> = HashMap::new();  // key is gene id, value is exon regions
    let file = File::open(anno_path).unwrap();
    let reader = BufReader::new(file);
    let mut invs: Vec<Interval<usize, u8>> = Vec::new(); // for merging gene exons
    let mut gene_id: String = String::new();
    for line in reader.lines() {
        let line = line.unwrap();
        if line.starts_with("#") { continue; }
        let parts = line.split('\t').collect::<Vec<&str>>();
        let seqname = parts[0].to_string();
        let feature = parts[2];
        let start = parts[3].parse::<u32>().unwrap();   // 1-based, inclusive
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
                    gene_id = subpart.replace("gene_id=", "");  // gff3 format
                    break;
                }
                if subpart.starts_with("gene_id") {
                    gene_id = subpart.replace("gene_id ", "").replace("\"", "");  // gtf format
                    break;
                }
            }
            let mut top = gene_regions.get_mut(&seqname).unwrap().pop_back();
            if top.is_some() {
                let mut top = top.unwrap();
                assert!(start >= top.start, "Error: annotation file is not sorted. {}:{}-{}", seqname, start, end);
                // if overlap, merge the overlapped regions
                if top.end <= start {
                    // end of top region is exclusive, so top.end == start is not overlap
                    gene_regions.get_mut(&seqname).unwrap().push_back(top);
                    gene_regions.get_mut(&seqname).unwrap().push_back(Region { chr: seqname.clone(), start: start, end: end + 1, gene_id: Option::from(gene_id.clone()) });
                } else if top.end < end + 1 {
                    // top.end is exclusive, end is inclusive
                    // merge two overlapped regions
                    top.end = end + 1;
                    top.gene_id = Option::from(top.gene_id.unwrap() + "," + &gene_id);
                    gene_regions.get_mut(&seqname).unwrap().push_back(top);
                } else {
                    // equal end or contained
                    top.gene_id = Option::from(top.gene_id.unwrap() + "," + &gene_id);
                    gene_regions.get_mut(&seqname).unwrap().push_back(top);
                }
            } else {
                // first gene region in stack
                gene_regions.get_mut(&seqname).unwrap().push_back(Region { chr: seqname.clone(), start: start, end: end + 1, gene_id: Option::from(gene_id.clone()) });
            }
        } else if feature == "exon" {
            let mut exon_gene_id = String::new();
            for subpart in parts[8].trim_end().split(";").collect::<Vec<&str>>() {
                if subpart.starts_with("gene_id") {
                    exon_gene_id = subpart.replace("gene_id=", "");  // gff3 format
                    break;
                }
                if subpart.starts_with("gene_id") {
                    exon_gene_id = subpart.replace("gene_id ", "").replace("\"", "");  // gtf format
                    break;
                }
            }
            assert!(exon_gene_id == gene_id, "Error: gene_id in gene and exon are different: gene_id:{}, exon_gene_id:{}", gene_id, exon_gene_id);
            invs.push(Interval { start: start as usize, stop: (end + 1) as usize, val: 0 });    // 1-based,
        } else {
            continue;
        }
    }
    if invs.len() > 0 {
        exon_regions.insert(gene_id.clone(), invs.clone());
        invs.clear();
    }
    return (gene_regions, exon_regions);
}

pub fn lapper_intervals(query_regions: &Vec<Region>, target_regions: &VecDeque<Region>) -> Vec<Region> {
    let mut result_regions = Vec::new();
    let mut invs: Vec<Interval<usize, String>> = Vec::new();
    for region in target_regions {
        invs.push(Interval { start: region.start as usize, stop: region.end as usize, val: region.gene_id.clone().unwrap() });
    }
    let mut interval_tree = Lapper::new(invs);
    for q in query_regions {
        let q_inv = Interval { start: q.start as usize, stop: q.end as usize, val: 0 };
        for h_inv in interval_tree.find(q.start as usize, q.end as usize) {
            let intersected_start = q_inv.start.max(h_inv.start);
            let intersected_end = q_inv.stop.min(h_inv.stop);
            let h_gene = h_inv.val.clone();
            assert!(intersected_start < intersected_end, "Error: intersected_start >= intersected_end, query:{:?}", q_inv);
            result_regions.push(Region { chr: q.chr.clone(), start: intersected_start as u32, end: intersected_end as u32, gene_id: Option::from(h_gene) });
        }
    }
    return result_regions;
}

pub fn intersect_gene_regions(alignment_regions: &Vec<Region>, gene_regions: &HashMap<String, VecDeque<Region>>, thread_size: usize) -> Vec<Region> {
    let intersected_regions: Mutex<Vec<Region>> = Mutex::new(Vec::new());
    let mut region_map = HashMap::new();
    for region in alignment_regions {
        if !region_map.contains_key(&region.chr) {
            region_map.insert(region.chr.clone(), Vec::new());
        }
        region_map.get_mut(&region.chr).unwrap().push(region.clone());
    }
    let pool = rayon::ThreadPoolBuilder::new().num_threads(thread_size - 1).build().unwrap();
    let contig_names = region_map.keys().collect::<Vec<&String>>();
    pool.install(|| {
        contig_names.par_iter().for_each(|ctg| {
            if gene_regions.contains_key(*ctg) {
                let chr_intersected_region = lapper_intervals(region_map.get(*ctg).unwrap(), gene_regions.get(*ctg).unwrap());
                for region in chr_intersected_region {
                    intersected_regions.lock().unwrap().push(region);
                }
            }
        });
    });

    return intersected_regions.into_inner().unwrap();
}

pub fn multithread_produce3(bam_file: String, ref_file: String, thread_size: usize, contigs: Option<Vec<String>>, min_mapq: u8, min_read_length: usize) -> Vec<Region> {
    let results: Mutex<Vec<Region>> = Mutex::new(Vec::new());
    let pool = rayon::ThreadPoolBuilder::new().num_threads(thread_size - 1).build().unwrap();
    let bam = bam::IndexedReader::from_path(bam_file.clone()).unwrap();
    let bam_header = bam.header().clone();
    let fai_path = ref_file + ".fai";
    if fs::metadata(&fai_path).is_err() {
        panic!("Reference index file .fai does not exist.");
    }
    let contig_lengths = parse_fai(fai_path.as_str());
    let mut contig_names: VecDeque<String> = VecDeque::new();
    if contigs.is_some() {
        for ctg in contigs.unwrap().iter() {
            contig_names.push_back(ctg.clone());
        }
    } else {
        // for ctg in bam_header.target_names() {
        //     contig_names.push_back(std::str::from_utf8(ctg).unwrap().to_string().clone());
        // }
        for (ctg, _) in contig_lengths.iter() {
            contig_names.push_back(ctg.clone());
        }
    }
    pool.install(|| {
        contig_names.par_iter().for_each(|ctg| {
            let isolated_regions = find_isolated_regions_with_depth(bam_file.as_str(), ctg, contig_lengths.iter().find(|(chr, _)| chr == ctg).unwrap().1, min_mapq, min_read_length);
            for region in isolated_regions {
                results.lock().unwrap().push(region);
            }
        });
    });
    return results.into_inner().unwrap().clone();
}

pub fn read_references(ref_path: &str) -> HashMap<String, Vec<u8>> {
    let mut references: HashMap<String, Vec<u8>> = HashMap::new();
    let mut reader = Reader::from_path(ref_path).unwrap();
    while let Some(record) = reader.next() {
        let record = record.expect("Error reading record");
        references.insert(record.id().unwrap().to_string(), record.full_seq().to_vec());
    }
    return references;
}

#[derive(Default, Debug, Clone)]
pub struct Profile {
    pub freq_vec: Vec<BaseFreq>,
    pub region: Region,
}

impl Profile {
    pub fn init_with_pileup(&mut self, bam_path: &str, region: &Region, ref_seq: &Vec<u8>, platform: &Platform, min_mapq: u8, min_baseq: u8, min_read_length: usize, min_depth: u32, max_depth: u32, distance_to_read_end: u32, polya_tail_length: u32) {
        // When region is large and the number of reads is large, the runtime of init_profile_with_pileup is time-consuming.
        // This function is used to fill the profile by parsing each read in the bam file instead of using pileup.

        let start_time = Instant::now();
        let mut bam: bam::IndexedReader = bam::IndexedReader::from_path(bam_path).unwrap();
        bam.fetch((region.chr.as_str(), region.start, region.end)).unwrap();
        let vec_size = (region.end - region.start) as usize;    // end is exclusive
        self.freq_vec = vec![BaseFreq::default(); vec_size];
        self.region = region.clone();
        let freq_vec_pos = region.start as usize - 1;    // the first position on reference, 0-based, inclusive
        let polyA_win = polya_tail_length as i64;

        // fill the ref_base field in each BaseFreq
        for i in 0..vec_size {
            self.freq_vec[i].ref_base = ref_seq[freq_vec_pos + i] as char;
        }

        for r in bam.records() {
            let record = r.unwrap();
            if record.mapq() < min_mapq || record.seq_len() < min_read_length || record.is_unmapped() || record.is_secondary() || record.is_supplementary() {
                continue;
            }
            let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
            let seq = record.seq();
            let base_qual = record.qual();
            let strand = if record.strand() == Forward { 0 } else { 1 };
            let start_pos = record.pos() as usize;  // 0-based
            let cigar = record.cigar();
            let leading_softclips = cigar.leading_softclips();
            let trailing_softclips = cigar.trailing_softclips();

            let mut pos_in_freq_vec: i32 = start_pos as i32 - freq_vec_pos as i32;
            let mut pos_in_read = if leading_softclips > 0 { leading_softclips as usize } else { 0 };
            for cg in cigar.iter() {
                match cg.char() as u8 {
                    b'S' | b'H' => {
                        continue;
                    }
                    b'M' | b'X' | b'=' => {
                        for cgi in 0..cg.len() {
                            if pos_in_freq_vec < 0 {
                                pos_in_freq_vec += 1;
                                pos_in_read += 1;
                                continue;
                            }
                            if pos_in_freq_vec >= vec_size as i32 {
                                break;
                            }
                            let base = seq[pos_in_read] as char;
                            let baseq = base_qual[pos_in_read];

                            // close to left read end or right read end, check whether current position is in polyA tail
                            let ref_base = self.freq_vec[pos_in_freq_vec as usize].ref_base;
                            let mut polyA_flag = false;
                            let mut homopolymer_flag = false;
                            let mut trimed_flag = false;    // trime the end of ont reads
                            // the end of ont reads are trimed since the accuracy of the ends is low
                            match platform {
                                Platform::ont => {
                                    if (pos_in_read as i64 - leading_softclips).abs() < distance_to_read_end as i64 || (pos_in_read as i64 - (seq.len() as i64 - trailing_softclips)).abs() < distance_to_read_end as i64 {
                                        trimed_flag = true;
                                    }
                                }
                                _ => {}
                            }
                            if !trimed_flag && ((pos_in_read as i64 - leading_softclips).abs() < distance_to_read_end as i64 || (pos_in_read as i64 - (seq.len() as i64 - trailing_softclips)).abs() < distance_to_read_end as i64) {
                                for tmpi in (pos_in_read as i64 - polyA_win)..=(pos_in_read as i64 + 1) {
                                    // pos_in_read is the current position, and the position 1-base to the left of polyA tail is often false positive variant allele. So the end for loop is pos_in_read+1 instead of pos_in_read.
                                    // same reason for pos_in_read - polyA_win instead fo pos_in_read - polyA_win + 1
                                    if tmpi < 0 || tmpi + polyA_win - 1 >= seq.len() as i64 {
                                        continue;
                                    }
                                    let mut polyA_cnt = 0;
                                    let mut polyT_cnt = 0;
                                    let mut polyC_cnt = 0;
                                    let mut polyG_cnt = 0;
                                    for tmpj in 0..polyA_win {
                                        if seq[(tmpi + tmpj) as usize] == b'A' && ref_base != 'A' {
                                            polyA_cnt += 1;
                                        } else if seq[(tmpi + tmpj) as usize] == b'T' && ref_base != 'T' {
                                            polyT_cnt += 1;
                                        } else if seq[(tmpi + tmpj) as usize] == b'C' && ref_base != 'C' {
                                            polyC_cnt += 1;
                                        } else if seq[(tmpi + tmpj) as usize] == b'G' && ref_base != 'G' {
                                            polyG_cnt += 1;
                                        }
                                    }
                                    if polyA_cnt >= polyA_win || polyT_cnt >= polyA_win {
                                        polyA_flag = true;
                                    }
                                    if polyC_cnt >= polyA_win || polyG_cnt >= polyA_win {
                                        homopolymer_flag = true;
                                    }
                                }
                            }

                            if !trimed_flag && !polyA_flag && !homopolymer_flag {

                                // calculate distance to read end of each allele, for filtering variants that the average distance of each allele is significantly different
                                let mut dist = 0;
                                if (pos_in_read as i64 - leading_softclips).abs() < (pos_in_read as i64 - (seq.len() as i64 - trailing_softclips)).abs() {
                                    dist = pos_in_read as i64 - leading_softclips;    // positive value
                                } else {
                                    dist = pos_in_read as i64 - (seq.len() as i64 - trailing_softclips);    // negative value
                                }

                                match base {
                                    'A' | 'a' => {
                                        self.freq_vec[pos_in_freq_vec as usize].a += 1;
                                        self.freq_vec[pos_in_freq_vec as usize].baseq.a.push(baseq);
                                        if strand == 0 {
                                            self.freq_vec[pos_in_freq_vec as usize].base_strands.a[0] += 1;
                                        } else {
                                            self.freq_vec[pos_in_freq_vec as usize].base_strands.a[1] += 1;
                                        }
                                        self.freq_vec[pos_in_freq_vec as usize].distance_to_end.a.push(dist);
                                    }
                                    'C' | 'c' => {
                                        self.freq_vec[pos_in_freq_vec as usize].c += 1;
                                        self.freq_vec[pos_in_freq_vec as usize].baseq.c.push(baseq);
                                        if strand == 0 {
                                            self.freq_vec[pos_in_freq_vec as usize].base_strands.c[0] += 1;
                                        } else {
                                            self.freq_vec[pos_in_freq_vec as usize].base_strands.c[1] += 1;
                                        }
                                        self.freq_vec[pos_in_freq_vec as usize].distance_to_end.c.push(dist);
                                    }
                                    'G' | 'g' => {
                                        self.freq_vec[pos_in_freq_vec as usize].g += 1;
                                        self.freq_vec[pos_in_freq_vec as usize].baseq.g.push(baseq);
                                        if strand == 0 {
                                            self.freq_vec[pos_in_freq_vec as usize].base_strands.g[0] += 1;
                                        } else {
                                            self.freq_vec[pos_in_freq_vec as usize].base_strands.g[1] += 1;
                                        }
                                        self.freq_vec[pos_in_freq_vec as usize].distance_to_end.g.push(dist);
                                    }
                                    'T' | 't' => {
                                        self.freq_vec[pos_in_freq_vec as usize].t += 1;
                                        self.freq_vec[pos_in_freq_vec as usize].baseq.t.push(baseq);
                                        if strand == 0 {
                                            self.freq_vec[pos_in_freq_vec as usize].base_strands.t[0] += 1;
                                        } else {
                                            self.freq_vec[pos_in_freq_vec as usize].base_strands.t[1] += 1;
                                        }
                                        self.freq_vec[pos_in_freq_vec as usize].distance_to_end.t.push(dist);
                                    }
                                    _ => {
                                        panic!("Invalid nucleotide base: {}", base);
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
        let end_time = Instant::now();
    }
    pub fn append_reference(&mut self, references: &HashMap<String, Vec<u8>>) {
        /*
        Fill the ``ref_base`` field in each BaseFreq.
        Optional. If not called, the ``ref_base`` field in each BaseFreq will be '\0'
         */
        let chr = &self.region.chr;
        let s = self.region.start - 1;    // 0-based, inclusive
        let mut p = s as usize;
        for i in 0..self.freq_vec.len() {
            if self.freq_vec[i].i {
                self.freq_vec[i].ref_base = '-'; // insertion
            } else {
                self.freq_vec[i].ref_base = references.get(chr).unwrap()[p] as char;
                p += 1;
            }
        }
    }
}

pub fn independent_test(data: [u32; 4]) -> f64 {
    // use chi-square test or fisher's exact test to test whether two variables are independent
    // N<40 or T<1, use fisher's exact test
    // N>=40, T>=5 in >80% of cells, use chi-square test
    // N>=40, T>=5 in <=80% of cells, use fisher's exact test
    // return true if independent, false if not independent
    let mut phred_pvalue = 0.0;
    let N = data[0] + data[1] + data[2] + data[3];
    let mut TL = [0.0; 4];
    TL[0] = (((data[0] + data[1]) * (data[0] + data[2])) as f32) / (N as f32);
    TL[1] = (((data[0] + data[1]) * (data[1] + data[3])) as f32) / (N as f32);
    TL[2] = (((data[2] + data[3]) * (data[0] + data[2])) as f32) / (N as f32);
    TL[3] = (((data[2] + data[3]) * (data[1] + data[3])) as f32) / (N as f32);
    let mut T = 0.0;
    for i in 0..data.len() {
        if TL[i] >= 5.0 {
            T += 1.0;
        }
    }
    T = T / data.len() as f32;

    if N > 40 && T > 0.8 {
        // use chi-square test
        let x = vec![data[0] as f64, data[1] as f64];
        let y = vec![data[2] as f64, data[3] as f64];
        let p = ChiSquare::test_vector(&x, &y).p_value();
        phred_pvalue = -10.0 * p.log10();
    }

    if N <= 40 {
        // use fisher's exact test
        let p = fishers_exact(&data).unwrap().two_tail_pvalue;
        phred_pvalue = -10.0 * p.log10();
    }

    return phred_pvalue;
}