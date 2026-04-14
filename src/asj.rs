use flate2::read::MultiGzDecoder;
use rayon::prelude::*;
use rust_htslib::bam;
use rust_htslib::bam::record::{Aux, Cigar};
use rust_htslib::bam::Read;
use rust_lapper::{Interval, Lapper};
use statrs::distribution::{ChiSquared, ContinuousCDF};
use statrs::function::gamma::ln_gamma;
use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

use crate::util::load_reference;

#[derive(clap::Parser, Debug, Clone)]
pub struct AsjArgs {
    /// Annotation file in GFF3 or GTF format
    #[arg(short = 'a', long = "annotation-file")]
    annotation_file: String,

    /// BAM file
    #[arg(short = 'b', long = "bam-file")]
    bam_file: String,

    /// DNA VCF for optional filtering (must be used together with --rna-vcf)
    #[arg(long = "dna-vcf")]
    dna_vcf: Option<String>,

    /// LongcallR RNA phased VCF for optional filtering (must be used together with --dna-vcf)
    #[arg(long = "rna-vcf")]
    rna_vcf: Option<String>,

    /// Remove reads with <= this many introns before ASJ analysis
    #[arg(long = "min-junctions", default_value_t = 2)]
    min_junctions: usize,

    /// Build Junction_set using exon-bridged connectivity (broader local splice clusters)
    ///
    /// Default (off): junctions are clustered only when they share donor/acceptor.
    /// Enabled: internal exons are added as bridge nodes, which can merge more nearby junctions.
    #[arg(long = "cluster-with-exons", default_value_t = false)]
    cluster_with_exons: bool,

    /// Reference genome FASTA file
    #[arg(short = 'f', long = "reference")]
    reference: String,

    /// Prefix of output files
    #[arg(short = 'o', long = "output-prefix")]
    output_prefix: String,

    /// Number of threads
    #[arg(short = 't', long = "threads", default_value_t = 1)]
    threads: usize,

    /// Gene types to be analyzed; pass --gene-types with no values to include all gene types [default: protein_coding lncRNA]
    #[arg(short = 'g', long = "gene-types", num_args(0..))]
    gene_types: Option<Vec<String>>,

    /// Minimum total phased reads required to test an ASJ event
    #[arg(short = 'm', long = "min-sup", default_value_t = 10)]
    min_sup: usize,

    /// Disable canonical splice-site check (GT-AG) for read introns
    ///
    /// Default (off): non-canonical introns are marked and can be excluded from gene summary.
    /// Enabled: skip GT-AG checking and treat all introns as allowed.
    #[arg(long = "no-gtag", default_value_t = false)]
    no_gtag: bool,
}

#[derive(Debug, Clone)]
struct GeneRegion {
    chr: String,
    start: u32,
    end: u32,
}

type TranscriptRegions = HashMap<String, Vec<(String, u32, u32)>>;
type GeneRegionsByTranscript = HashMap<String, TranscriptRegions>;

#[derive(Debug, Clone)]
struct ReadTag {
    ps: Option<i32>,
    hp: Option<i32>,
}

#[derive(Debug, Clone)]
struct AseEvent {
    chr: String,
    start: u32,
    end: u32,
    strand: String,
    junction_set: String,
    phase_set: String,
    hap1_absent: usize,
    hap1_present: usize,
    hap2_absent: usize,
    hap2_present: usize,
    p_value: f64,
    sor: f64,
    novel: bool,
    gt_ag_tag: bool,
    gene_name: String,
}

impl AseEvent {
    fn header() -> &'static str {
        "#Junction\tStrand\tJunction_set\tPhase_set\tHap1_absent\tHap1_present\tHap2_absent\tHap2_present\tP_value\tSOR\tNovel\tGT_AG\tGene_name"
    }

    fn line(&self) -> String {
        format!(
            "{}:{}-{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.chr,
            self.start,
            self.end,
            self.strand,
            self.junction_set,
            self.phase_set,
            self.hap1_absent,
            self.hap1_present,
            self.hap2_absent,
            self.hap2_present,
            self.p_value,
            self.sor,
            if self.novel { "True" } else { "False" },
            if self.gt_ag_tag { "True" } else { "False" },
            self.gene_name
        )
    }
}

fn event_total_coverage(event: &AseEvent) -> usize {
    event.hap1_absent + event.hap1_present + event.hap2_absent + event.hap2_present
}

fn event_stable_key(event: &AseEvent) -> (&str, u32, u32, &str, &str) {
    (
        event.chr.as_str(),
        event.start,
        event.end,
        event.phase_set.as_str(),
        event.gene_name.as_str(),
    )
}

fn is_better_gene_event(candidate: &AseEvent, current: &AseEvent) -> bool {
    let eps = 1e-12;
    if candidate.p_value + eps < current.p_value {
        return true;
    }
    if current.p_value + eps < candidate.p_value {
        return false;
    }

    let cand_cov = event_total_coverage(candidate);
    let curr_cov = event_total_coverage(current);
    if cand_cov > curr_cov {
        return true;
    }
    if cand_cov < curr_cov {
        return false;
    }

    if candidate.sor > current.sor + eps {
        return true;
    }
    if current.sor > candidate.sor + eps {
        return false;
    }

    event_stable_key(candidate) < event_stable_key(current)
}

#[derive(Debug, Clone)]
struct ReadDataBundle {
    read_assignment: HashMap<String, String>,       // junction-filtered, for ASJ analysis
    all_read_assignment: HashMap<String, String>,   // all reads, for coverage
    reads_positions: HashMap<String, (u32, u32)>,
    reads_tags: HashMap<String, ReadTag>,
    reads_exons: HashMap<String, Vec<(u32, u32)>>,
    reads_introns: HashMap<String, Vec<(u32, u32, bool)>>,
}

fn open_text_reader(path: &str) -> BufReader<Box<dyn std::io::Read>> {
    let file = File::open(path).unwrap_or_else(|e| panic!("failed to open {}: {}", path, e));
    if path.ends_with(".gz") {
        BufReader::new(Box::new(MultiGzDecoder::new(file)))
    } else {
        BufReader::new(Box::new(file))
    }
}

fn parse_attributes_gtf(attrs: &str) -> HashMap<String, String> {
    let mut attr_map: HashMap<String, String> = HashMap::new();
    let mut tags: Vec<String> = Vec::new();
    for item in attrs.split(';') {
        let trimmed = item.trim();
        if trimmed.is_empty() {
            continue;
        }
        let mut parts = trimmed.splitn(2, ' ');
        let key = parts.next().unwrap_or("").trim();
        let value = parts
            .next()
            .unwrap_or("")
            .trim()
            .trim_matches('"')
            .to_string();
        if key == "tag" {
            if !value.is_empty() {
                tags.push(value);
            }
        } else if !key.is_empty() {
            attr_map.insert(key.to_string(), value);
        }
    }
    attr_map.insert("tag".to_string(), tags.join(","));
    attr_map
}

fn parse_attributes_gff3(attrs: &str) -> HashMap<String, String> {
    let mut attr_map: HashMap<String, String> = HashMap::new();
    for item in attrs.split(';') {
        let trimmed = item.trim();
        if trimmed.is_empty() {
            continue;
        }
        let mut parts = trimmed.splitn(2, '=');
        let key = parts.next().unwrap_or("").trim();
        let value = parts
            .next()
            .unwrap_or("")
            .trim()
            .trim_matches('"')
            .to_string();
        if !key.is_empty() {
            attr_map.insert(key.to_string(), value);
        }
    }
    attr_map
}

#[allow(clippy::type_complexity)]
fn get_gene_regions(
    annotation_file: &str,
    gene_types: &HashSet<String>,
) -> (
    HashMap<String, GeneRegion>,
    HashMap<String, String>,
    HashMap<String, String>,
    GeneRegionsByTranscript,
    GeneRegionsByTranscript,
) {
    assert!(
        annotation_file.ends_with(".gff3")
            || annotation_file.ends_with(".gtf")
            || annotation_file.ends_with(".gff3.gz")
            || annotation_file.ends_with(".gtf.gz"),
        "Error: Unknown annotation file format"
    );

    let file_type = if annotation_file.contains(".gff3") {
        "gff3"
    } else {
        "gtf"
    };

    let mut gene_regions: HashMap<String, GeneRegion> = HashMap::new();
    let mut gene_names: HashMap<String, String> = HashMap::new();
    let mut gene_strands: HashMap<String, String> = HashMap::new();
    let mut exon_regions: GeneRegionsByTranscript = HashMap::new();

    let reader = open_text_reader(annotation_file);
    for line in reader.lines() {
        let line = line.unwrap();
        if line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 9 {
            continue;
        }
        let feature = parts[2];
        let attrs = if file_type == "gff3" {
            parse_attributes_gff3(parts[8])
        } else {
            parse_attributes_gtf(parts[8])
        };

        let gene_id = match attrs.get("gene_id") {
            Some(v) => v.clone(),
            None => continue,
        };
        let gene_type = attrs
            .get("gene_type")
            .or_else(|| attrs.get("gene_biotype"))
            .cloned()
            .unwrap_or_else(|| "".to_string());
        let tag = attrs.get("tag").cloned().unwrap_or_else(|| "".to_string());
        if (!gene_types.is_empty() && !gene_types.contains(&gene_type)) || tag.contains("readthrough") {
            continue;
        }

        if feature == "gene" {
            let chr = parts[0].to_string();
            let start = parts[3].parse::<u32>().unwrap_or(0);
            let end = parts[4].parse::<u32>().unwrap_or(0);
            if start == 0 || end == 0 || end < start {
                continue;
            }
            let strand = parts[6].to_string();
            let gene_name = attrs
                .get("gene_name")
                .cloned()
                .unwrap_or_else(|| ".".to_string());
            gene_regions.insert(gene_id.clone(), GeneRegion { chr, start, end });
            gene_names.insert(gene_id.clone(), gene_name);
            gene_strands.insert(gene_id, strand);
        } else if feature == "exon" {
            let transcript_id = match attrs.get("transcript_id") {
                Some(v) => v.clone(),
                None => continue,
            };
            let chr = parts[0].to_string();
            let start = parts[3].parse::<u32>().unwrap_or(0);
            let end = parts[4].parse::<u32>().unwrap_or(0);
            if start == 0 || end == 0 || end < start {
                continue;
            }
            exon_regions
                .entry(gene_id)
                .or_default()
                .entry(transcript_id)
                .or_default()
                .push((chr, start, end));
        }
    }

    let mut intron_regions: GeneRegionsByTranscript = HashMap::new();
    for (gene_id, tx_exons) in &exon_regions {
        for (tx_id, exons) in tx_exons {
            if exons.len() <= 1 {
                continue;
            }
            let mut sorted_exons = exons.clone();
            sorted_exons.sort_by_key(|x| x.1);
            for i in 1..sorted_exons.len() {
                let intron_start = sorted_exons[i - 1].2 + 1;
                let intron_end = sorted_exons[i].1.saturating_sub(1);
                if intron_start < intron_end {
                    intron_regions
                        .entry(gene_id.clone())
                        .or_default()
                        .entry(tx_id.clone())
                        .or_default()
                        .push((sorted_exons[i - 1].0.clone(), intron_start, intron_end));
                }
            }
        }
    }

    (
        gene_regions,
        gene_names,
        gene_strands,
        exon_regions,
        intron_regions,
    )
}

fn merge_gene_exon_regions(
    exon_regions: &GeneRegionsByTranscript,
) -> HashMap<String, HashMap<String, Vec<(u32, u32)>>> {
    let mut merged: HashMap<String, HashMap<String, Vec<(u32, u32)>>> = HashMap::new();
    for (gene_id, transcripts) in exon_regions {
        let mut all_exons: Vec<(String, u32, u32)> = Vec::new();
        let mut chr_set: HashSet<String> = HashSet::new();
        for exons in transcripts.values() {
            for (chr, s, e) in exons {
                chr_set.insert(chr.clone());
                all_exons.push((chr.clone(), *s, *e));
            }
        }
        if all_exons.is_empty() || chr_set.len() != 1 {
            continue;
        }
        let chr = chr_set.into_iter().next().unwrap();
        let mut spans: Vec<(u32, u32)> = all_exons.into_iter().map(|(_, s, e)| (s, e)).collect();
        spans.sort_by_key(|x| x.0);
        let mut merged_spans: Vec<(u32, u32)> = Vec::new();
        for (s, e) in spans {
            if let Some(last) = merged_spans.last_mut() {
                if s <= last.1 + 1 {
                    if e > last.1 {
                        last.1 = e;
                    }
                } else {
                    merged_spans.push((s, e));
                }
            } else {
                merged_spans.push((s, e));
            }
        }
        merged
            .entry(chr)
            .or_default()
            .insert(gene_id.clone(), merged_spans);
    }
    merged
}

type GeneSpanTrees = HashMap<String, Lapper<u32, String>>;
type GeneExonTrees = HashMap<String, HashMap<String, Lapper<u32, ()>>>;

fn build_gene_trees(
    merged_genes_exons: &HashMap<String, HashMap<String, Vec<(u32, u32)>>>,
) -> (GeneSpanTrees, GeneExonTrees) {
    let mut span_trees: GeneSpanTrees = HashMap::new();
    let mut exon_trees: GeneExonTrees = HashMap::new();

    for (chr, gene_map) in merged_genes_exons {
        let mut span_intervals: Vec<Interval<u32, String>> = Vec::new();
        let mut chr_exon_trees: HashMap<String, Lapper<u32, ()>> = HashMap::new();
        for (gene_id, exons) in gene_map {
            if exons.is_empty() {
                continue;
            }
            span_intervals.push(Interval {
                start: exons.first().unwrap().0,
                stop: exons.last().unwrap().1 + 1,
                val: gene_id.clone(),
            });

            let mut ex_intervals: Vec<Interval<u32, ()>> = exons
                .iter()
                .map(|(s, e)| Interval {
                    start: *s,
                    stop: *e + 1,
                    val: (),
                })
                .collect();
            ex_intervals.sort_by_key(|x| x.start);
            chr_exon_trees.insert(gene_id.clone(), Lapper::new(ex_intervals));
        }
        span_intervals.sort_by_key(|x| x.start);
        span_trees.insert(chr.clone(), Lapper::new(span_intervals));
        exon_trees.insert(chr.clone(), chr_exon_trees);
    }
    (span_trees, exon_trees)
}

fn aux_to_i32(aux: Aux<'_>) -> Option<i32> {
    match aux {
        Aux::I8(v) => Some(v as i32),
        Aux::U8(v) => Some(v as i32),
        Aux::I16(v) => Some(v as i32),
        Aux::U16(v) => Some(v as i32),
        Aux::I32(v) => Some(v),
        Aux::U32(v) => Some(v as i32),
        _ => None,
    }
}

fn parse_splice_regions(record: &bam::Record) -> Vec<(u32, u32)> {
    let mut regions: Vec<(u32, u32)> = Vec::new();
    let mut current_pos = record.pos() as u32 + 1;
    let mut shift: u32 = 0;
    for op in &record.cigar() {
        match op {
            Cigar::Match(l) | Cigar::Del(l) | Cigar::Equal(l) | Cigar::Diff(l) => {
                shift += *l;
            }
            Cigar::RefSkip(l) => {
                if shift > 0 {
                    regions.push((current_pos, current_pos + shift - 1));
                }
                current_pos += shift + *l;
                shift = 0;
            }
            _ => {}
        }
    }
    if shift > 0 {
        regions.push((current_pos, current_pos + shift - 1));
    }
    regions
}

fn intron_gtag_tag(seq: &[u8], intron_start: u32, intron_end: u32) -> bool {
    if intron_start < 1 || intron_end < 2 {
        return false;
    }
    let left_start = (intron_start - 1) as usize;
    let left_end = (intron_start + 1) as usize;
    let right_start = (intron_end.saturating_sub(2)) as usize;
    let right_end = intron_end as usize;
    if left_end > seq.len() || right_end > seq.len() || right_start >= right_end {
        return false;
    }
    let left = &seq[left_start..left_end];
    let right = &seq[right_start..right_end];
    let l0 = (left[0] as char).to_ascii_uppercase();
    let l1 = (left[1] as char).to_ascii_uppercase();
    let r0 = (right[0] as char).to_ascii_uppercase();
    let r1 = (right[1] as char).to_ascii_uppercase();
    (l0 == 'G' && l1 == 'T' && r0 == 'A' && r1 == 'G')
        || (l0 == 'C' && l1 == 'T' && r0 == 'A' && r1 == 'C')
}

fn get_exon_intron_regions(
    record: &bam::Record,
    ref_seq: &[u8],
    no_gtag: bool,
) -> (Vec<(u32, u32)>, Vec<(u32, u32, bool)>) {
    let mut exon_regions: Vec<(u32, u32)> = Vec::new();
    let mut intron_regions: Vec<(u32, u32, bool)> = Vec::new();
    let mut current_position = record.pos() as u32 + 1;

    for cg in &record.cigar() {
        match cg {
            Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) => {
                if let Some(last) = exon_regions.last_mut() {
                    if last.1 + 1 == current_position {
                        last.1 += *l;
                    } else {
                        exon_regions.push((current_position, current_position + *l - 1));
                    }
                } else {
                    exon_regions.push((current_position, current_position + *l - 1));
                }
                current_position += *l;
            }
            Cigar::Del(l) => {
                if let Some(last) = exon_regions.last_mut() {
                    if last.1 + 1 == current_position {
                        last.1 += *l;
                    } else {
                        exon_regions.push((current_position, current_position + *l - 1));
                    }
                } else {
                    exon_regions.push((current_position, current_position + *l - 1));
                }
                current_position += *l;
            }
            Cigar::RefSkip(l) => {
                let intron_start = current_position;
                let intron_end = current_position + *l - 1;
                let tag = if no_gtag {
                    false
                } else {
                    intron_gtag_tag(ref_seq, intron_start, intron_end)
                };
                intron_regions.push((intron_start, intron_end, tag));
                current_position += *l;
            }
            _ => {}
        }
    }
    (exon_regions, intron_regions)
}

fn load_reads(
    bam_file: &str,
    genome_dict: &HashMap<String, Vec<u8>>,
    span_trees: &GeneSpanTrees,
    exon_trees: &GeneExonTrees,
    no_gtag: bool,
    min_junctions: usize,
) -> ReadDataBundle {
    let mut read_assignment: HashMap<String, String> = HashMap::new();
    let mut all_read_assignment: HashMap<String, String> = HashMap::new();
    let mut reads_positions: HashMap<String, (u32, u32)> = HashMap::new();
    let mut reads_tags: HashMap<String, ReadTag> = HashMap::new();
    let mut reads_exons: HashMap<String, Vec<(u32, u32)>> = HashMap::new();
    let mut reads_introns: HashMap<String, Vec<(u32, u32, bool)>> = HashMap::new();

    let mut bam = bam::Reader::from_path(bam_file)
        .unwrap_or_else(|e| panic!("failed to open BAM {}: {}", bam_file, e));
    let header = bam.header().to_owned();
    for rec in bam.records() {
        let record = rec.unwrap();
        if record.is_unmapped() || record.tid() < 0 {
            continue;
        }
        let chr = String::from_utf8_lossy(header.tid2name(record.tid() as u32)).to_string();
        let ref_seq = match genome_dict.get(&chr) {
            Some(s) => s,
            None => continue,
        };

        let qname = String::from_utf8_lossy(record.qname()).to_string();
        let ps = record.aux(b"PS").ok().and_then(aux_to_i32);
        let hp = record.aux(b"HP").ok().and_then(aux_to_i32);
        let read_start_1 = record.pos() as u32 + 1;
        let read_end_1 = record.cigar().end_pos() as u32;

        reads_positions.insert(qname.clone(), (read_start_1, read_end_1));
        reads_tags.insert(qname.clone(), ReadTag { ps, hp });

        let (exon_regions, intron_regions) = get_exon_intron_regions(&record, ref_seq, no_gtag);

        // compute gene assignment for all reads (used for coverage)
        let span_tree = match span_trees.get(&chr) {
            Some(v) => v,
            None => {
                if intron_regions.len() <= min_junctions {
                    reads_positions.remove(&qname);
                    reads_tags.remove(&qname);
                }
                continue;
            }
        };
        let chr_exons = match exon_trees.get(&chr) {
            Some(v) => v,
            None => {
                if intron_regions.len() <= min_junctions {
                    reads_positions.remove(&qname);
                    reads_tags.remove(&qname);
                }
                continue;
            }
        };
        let start_pos0 = record.pos() as u32;
        let end_pos0 = record.cigar().end_pos() as u32;
        let candidates: Vec<(String, u32, u32)> = span_tree
            .find(start_pos0 + 1, end_pos0 + 1)
            .map(|x| (x.val.clone(), x.start, x.stop))
            .collect();

        let splice_regions = parse_splice_regions(&record);
        let mut best_gene: Option<String> = None;
        let mut best_overlap = 0u32;
        let mut best_start = 0u32;
        let mut best_span = u32::MAX;
        for (gene_id, gene_start, gene_stop) in candidates {
            let exon_tree = match chr_exons.get(&gene_id) {
                Some(v) => v,
                None => continue,
            };
            let mut overlap_sum = 0u32;
            for (sp_s, sp_e) in &splice_regions {
                for exon in exon_tree.find(*sp_s, *sp_e + 1) {
                    let left = (*sp_s).max(exon.start);
                    let right = (*sp_e + 1).min(exon.stop);
                    if right > left {
                        overlap_sum += right - left;
                    }
                }
            }

            if overlap_sum == 0 {
                continue;
            }

            let gene_span = gene_stop.saturating_sub(gene_start);
            let is_better = if overlap_sum > best_overlap {
                true
            } else if overlap_sum < best_overlap {
                false
            } else if gene_start > best_start {
                true
            } else if gene_start < best_start {
                false
            } else if gene_span < best_span {
                true
            } else if gene_span > best_span {
                false
            } else {
                match &best_gene {
                    Some(current_best) => gene_id < *current_best,
                    None => true,
                }
            };

            if is_better {
                best_overlap = overlap_sum;
                best_start = gene_start;
                best_span = gene_span;
                best_gene = Some(gene_id);
            }
        }
        if let Some(gid) = best_gene {
            all_read_assignment.insert(qname.clone(), gid.clone());
            // only include in the junction-filtered assignment if the read passes the filter
            if intron_regions.len() > min_junctions {
                read_assignment.insert(qname.clone(), gid);
            }
        }

        // filter artifacts with too few junctions for ASJ analysis
        if intron_regions.len() <= min_junctions {
            reads_positions.remove(&qname);
            reads_tags.remove(&qname);
            continue;
        }
        reads_exons.insert(qname.clone(), exon_regions);
        reads_introns.insert(qname.clone(), intron_regions);
    }

    ReadDataBundle {
        read_assignment,
        all_read_assignment,
        reads_positions,
        reads_tags,
        reads_exons,
        reads_introns,
    }
}

fn transform_read_assignment(
    read_assignment: &HashMap<String, String>,
) -> HashMap<String, Vec<String>> {
    let mut gene_assigned_reads: HashMap<String, Vec<String>> = HashMap::new();
    for (read_name, gene_id) in read_assignment {
        gene_assigned_reads
            .entry(gene_id.clone())
            .or_default()
            .push(read_name.clone());
    }
    gene_assigned_reads
}

fn cluster_junctions_connected_components(
    reads_junctions: &HashMap<String, Vec<(u32, u32, bool)>>,
    min_count: usize,
) -> Vec<Vec<(u32, u32, bool)>> {
    let mut counts: HashMap<(u32, u32), usize> = HashMap::new();
    let mut gt_ag_dict: HashMap<(u32, u32), bool> = HashMap::new();
    for juncs in reads_junctions.values() {
        for (s, e, tag) in juncs {
            *counts.entry((*s, *e)).or_insert(0) += 1;
            gt_ag_dict.insert((*s, *e), *tag);
        }
    }
    let nodes: Vec<(u32, u32)> = counts
        .iter()
        .filter_map(|(k, c)| if *c >= min_count { Some(*k) } else { None })
        .collect();
    if nodes.is_empty() {
        return Vec::new();
    }

    let mut parent: Vec<usize> = (0..nodes.len()).collect();
    let find = |x: usize, parent: &mut Vec<usize>| -> usize {
        let mut node = x;
        while parent[node] != node {
            parent[node] = parent[parent[node]];
            node = parent[node];
        }
        node
    };
    let union = |a: usize, b: usize, parent: &mut Vec<usize>| {
        let ra = find(a, parent);
        let rb = find(b, parent);
        if ra != rb {
            parent[rb] = ra;
        }
    };

    for i in 0..nodes.len() {
        for j in (i + 1)..nodes.len() {
            if nodes[i].0 == nodes[j].0 || nodes[i].1 == nodes[j].1 {
                union(i, j, &mut parent);
            }
        }
    }

    let mut comp: HashMap<usize, Vec<(u32, u32, bool)>> = HashMap::new();
    for (idx, (s, e)) in nodes.iter().enumerate() {
        let root = find(idx, &mut parent);
        let tag = *gt_ag_dict.get(&(*s, *e)).unwrap_or(&false);
        comp.entry(root).or_default().push((*s, *e, tag));
    }
    comp.into_values().collect()
}

fn cluster_junctions_exons_connected_components(
    reads_junctions: &HashMap<String, Vec<(u32, u32, bool)>>,
    reads_exons: &HashMap<String, Vec<(u32, u32)>>,
    min_count: usize,
) -> Vec<Vec<(u32, u32, bool)>> {
    let mut junction_counts: HashMap<(u32, u32), usize> = HashMap::new();
    let mut gt_ag_dict: HashMap<(u32, u32), bool> = HashMap::new();
    for juncs in reads_junctions.values() {
        for (s, e, tag) in juncs {
            *junction_counts.entry((*s, *e)).or_insert(0) += 1;
            gt_ag_dict.insert((*s, *e), *tag);
        }
    }
    let junctions: Vec<(u32, u32)> = junction_counts
        .iter()
        .filter_map(|(k, v)| if *v >= min_count { Some(*k) } else { None })
        .collect();

    let mut exon_counts: HashMap<(u32, u32), usize> = HashMap::new();
    for exon_regions in reads_exons.values() {
        if exon_regions.len() > 2 {
            for (idx, exon) in exon_regions.iter().enumerate() {
                if idx == 0 || idx == exon_regions.len() - 1 {
                    continue;
                }
                *exon_counts.entry(*exon).or_insert(0) += 1;
            }
        }
    }
    let exons: Vec<(u32, u32)> = exon_counts
        .iter()
        .filter_map(|(k, v)| if *v >= min_count { Some(*k) } else { None })
        .collect();

    #[derive(Clone)]
    struct Node {
        start: u32,
        end: u32,
        node_type: u8,
    }

    let mut nodes: Vec<Node> = Vec::new();
    for (s, e) in &junctions {
        nodes.push(Node {
            start: *s,
            end: *e,
            node_type: 0,
        });
    }
    for (s, e) in &exons {
        nodes.push(Node {
            start: s.saturating_sub(1),
            end: *e + 1,
            node_type: 1,
        });
    }
    if nodes.is_empty() {
        return Vec::new();
    }

    let mut parent: Vec<usize> = (0..nodes.len()).collect();
    let find = |x: usize, parent: &mut Vec<usize>| -> usize {
        let mut node = x;
        while parent[node] != node {
            parent[node] = parent[parent[node]];
            node = parent[node];
        }
        node
    };
    let union = |a: usize, b: usize, parent: &mut Vec<usize>| {
        let ra = find(a, parent);
        let rb = find(b, parent);
        if ra != rb {
            parent[rb] = ra;
        }
    };

    for i in 0..nodes.len() {
        for j in (i + 1)..nodes.len() {
            let n1 = &nodes[i];
            let n2 = &nodes[j];
            if n1.node_type == n2.node_type {
                if n1.start == n2.start || n1.end == n2.end {
                    union(i, j, &mut parent);
                }
            } else if n1.start == n2.end || n1.end == n2.start {
                union(i, j, &mut parent);
            }
        }
    }

    let mut comp: HashMap<usize, Vec<(u32, u32, bool)>> = HashMap::new();
    for (idx, node) in nodes.iter().enumerate() {
        if node.node_type != 0 {
            continue;
        }
        let root = find(idx, &mut parent);
        let tag = *gt_ag_dict.get(&(node.start, node.end)).unwrap_or(&false);
        comp.entry(root)
            .or_default()
            .push((node.start, node.end, tag));
    }
    comp.into_values().collect()
}

fn check_absent_present(
    start_pos: u32,
    end_pos: u32,
    reads_positions: &HashMap<String, (u32, u32)>,
    reads_junctions: &HashMap<String, Vec<(u32, u32, bool)>>,
) -> (Vec<String>, Vec<String>) {
    let mut absent_reads: Vec<String> = Vec::new();
    let mut present_reads: Vec<String> = Vec::new();
    for (read_name, (read_start, read_end)) in reads_positions {
        if *read_start > end_pos || *read_end < start_pos {
            continue;
        }
        let mut present = false;
        if let Some(juncs) = reads_junctions.get(read_name) {
            for (js, je, _) in juncs {
                if *js == start_pos && *je == end_pos {
                    present_reads.push(read_name.clone());
                    present = true;
                    break;
                }
            }
        }
        if !present {
            absent_reads.push(read_name.clone());
        }
    }
    (absent_reads, present_reads)
}

fn calc_sor(
    hap1_absent: usize,
    hap1_present: usize,
    hap2_absent: usize,
    hap2_present: usize,
) -> f64 {
    let r = ((hap1_absent as f64 + 1.0) * (hap2_present as f64 + 1.0))
        / ((hap1_present as f64 + 1.0) * (hap2_absent as f64 + 1.0));
    (r + 1.0 / r).ln()
}

fn ln_comb(n: usize, k: usize) -> f64 {
    if k > n {
        return f64::NEG_INFINITY;
    }
    ln_gamma((n + 1) as f64) - ln_gamma((k + 1) as f64) - ln_gamma((n - k + 1) as f64)
}

fn hypergeom_pmf(x: usize, row1: usize, col1: usize, n: usize) -> f64 {
    if x > row1 || x > col1 || row1 > n || col1 > n {
        return 0.0;
    }
    let col2 = n - col1;
    if row1 < x || row1 - x > col2 {
        return 0.0;
    }
    let ln_p = ln_comb(col1, x) + ln_comb(col2, row1 - x) - ln_comb(n, row1);
    ln_p.exp()
}

fn fisher_exact_two_sided(a: usize, b: usize, c: usize, d: usize) -> f64 {
    let row1 = a + b;
    let row2 = c + d;
    let col1 = a + c;
    let n = row1 + row2;
    let low = row1.saturating_sub(n - col1);
    let high = row1.min(col1);
    let p_obs = hypergeom_pmf(a, row1, col1, n);
    let tie_tol = (p_obs.abs() * 1e-12).max(1e-15);
    let mut p = 0.0;
    for x in low..=high {
        let px = hypergeom_pmf(x, row1, col1, n);
        if px <= p_obs + tie_tol {
            p += px;
        }
    }
    p.min(1.0)
}

fn g_test_2x2(a: usize, b: usize, c: usize, d: usize, pseudocount: f64) -> f64 {
    let row1 = a as f64 + b as f64;
    let row2 = c as f64 + d as f64;
    let col1 = a as f64 + c as f64;
    let col2 = b as f64 + d as f64;
    let grand = row1 + row2;
    if grand <= 0.0 {
        return 1.0;
    }
    let e11 = row1 * col1 / grand + pseudocount;
    let e12 = row1 * col2 / grand + pseudocount;
    let e21 = row2 * col1 / grand + pseudocount;
    let e22 = row2 * col2 / grand + pseudocount;
    let o11 = a as f64 + pseudocount;
    let o12 = b as f64 + pseudocount;
    let o21 = c as f64 + pseudocount;
    let o22 = d as f64 + pseudocount;
    let g = 2.0
        * (o11 * (o11 / e11).ln()
            + o12 * (o12 / e12).ln()
            + o21 * (o21 / e21).ln()
            + o22 * (o22 / e22).ln());
    let chi = ChiSquared::new(1.0).unwrap();
    1.0 - chi.cdf(g)
}

fn haplotype_event_test(
    absent_reads: &[String],
    present_reads: &[String],
    reads_tags: &HashMap<String, ReadTag>,
) -> Option<(i32, usize, usize, usize, usize, f64, f64)> {
    let mut absent_counts: HashMap<i32, (usize, usize)> = HashMap::new();
    let mut present_counts: HashMap<i32, (usize, usize)> = HashMap::new();

    for read_name in absent_reads {
        if let Some(tag) = reads_tags.get(read_name) {
            if let (Some(ps), Some(hp)) = (tag.ps, tag.hp) {
                let entry = absent_counts.entry(ps).or_insert((0, 0));
                if hp == 1 {
                    entry.0 += 1;
                } else if hp == 2 {
                    entry.1 += 1;
                }
            }
        }
    }
    for read_name in present_reads {
        if let Some(tag) = reads_tags.get(read_name) {
            if let (Some(ps), Some(hp)) = (tag.ps, tag.hp) {
                let entry = present_counts.entry(ps).or_insert((0, 0));
                if hp == 1 {
                    entry.0 += 1;
                } else if hp == 2 {
                    entry.1 += 1;
                }
            }
        }
    }

    let mut all_ps: HashSet<i32> = HashSet::new();
    all_ps.extend(absent_counts.keys().copied());
    all_ps.extend(present_counts.keys().copied());
    if all_ps.is_empty() {
        return None;
    }

    let mut ps_read_count: Vec<(i32, usize)> = all_ps
        .iter()
        .map(|ps| {
            let (h1a, h2a) = absent_counts.get(ps).copied().unwrap_or((0, 0));
            let (h1p, h2p) = present_counts.get(ps).copied().unwrap_or((0, 0));
            (*ps, h1a + h2a + h1p + h2p)
        })
        .collect();
    ps_read_count.sort_by(|a, b| b.1.cmp(&a.1));
    let phase_set = ps_read_count[0].0;

    let (h1a, h2a) = absent_counts.get(&phase_set).copied().unwrap_or((0, 0));
    let (h1p, h2p) = present_counts.get(&phase_set).copied().unwrap_or((0, 0));

    let p_fisher = fisher_exact_two_sided(h1a, h2a, h1p, h2p);
    let p_gtest = g_test_2x2(h1a, h2a, h1p, h2p, 1e-10);
    let p_value = p_fisher.max(p_gtest);
    let sor = calc_sor(h1a, h1p, h2a, h2p);

    Some((phase_set, h1a, h1p, h2a, h2p, p_value, sor))
}

fn benjamini_hochberg(p_values: &[f64]) -> Vec<f64> {
    let m = p_values.len();
    if m == 0 {
        return Vec::new();
    }
    let mut indexed: Vec<(usize, f64)> = p_values.iter().copied().enumerate().collect();
    indexed.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(Ordering::Equal));
    let mut adjusted_sorted = vec![1.0; m];
    let mut min_coeff = 1.0;
    for (rank_from_end, (_, p)) in indexed.iter().enumerate().rev() {
        let rank = rank_from_end + 1;
        let adj = (*p * m as f64 / rank as f64).min(1.0);
        if adj < min_coeff {
            min_coeff = adj;
        }
        adjusted_sorted[rank_from_end] = min_coeff;
    }
    let mut adjusted = vec![1.0; m];
    for (sorted_idx, (orig_idx, _)) in indexed.iter().enumerate() {
        adjusted[*orig_idx] = adjusted_sorted[sorted_idx];
    }
    adjusted
}

fn load_dna_vcf(vcf_file: &str) -> HashSet<String> {
    let mut dna_vcfs: HashSet<String> = HashSet::new();
    let reader = open_text_reader(vcf_file);
    for line in reader.lines() {
        let line = line.unwrap();
        if line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 10 {
            continue;
        }
        let chr = fields[0];
        let pos = match fields[1].parse::<u32>() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let ref_allele = fields[3];
        let alt = fields[4].split(',').next().unwrap_or("");
        if ref_allele.len() != 1 || alt.len() != 1 {
            continue;
        }
        let keys: Vec<&str> = fields[8].split(':').collect();
        let vals: Vec<&str> = fields[9].split(':').collect();
        let mut fmt: HashMap<&str, &str> = HashMap::new();
        for (k, v) in keys.iter().zip(vals.iter()) {
            fmt.insert(*k, *v);
        }
        let gt = fmt.get("GT").copied().unwrap_or("");
        if gt == "0/1" || gt == "1/0" || gt == "0|1" || gt == "1|0" {
            dna_vcfs.insert(format!("{}:{}", chr, pos));
        }
    }
    dna_vcfs
}

fn load_longcallr_phased_vcf(vcf_file: &str) -> HashMap<i32, Vec<(String, u32)>> {
    let mut rna_vcfs: HashMap<i32, Vec<(String, u32)>> = HashMap::new();
    let reader = open_text_reader(vcf_file);
    for line in reader.lines() {
        let line = line.unwrap();
        if line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 10 {
            continue;
        }
        if !fields[6].split(';').any(|x| x == "PASS") {
            continue;
        }
        let chr = fields[0].to_string();
        let pos = match fields[1].parse::<u32>() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let ref_allele = fields[3];
        let alt = fields[4].split(',').next().unwrap_or("");
        if ref_allele.len() != 1 || alt.len() != 1 {
            continue;
        }
        let keys: Vec<&str> = fields[8].split(':').collect();
        let vals: Vec<&str> = fields[9].split(':').collect();
        let mut fmt: HashMap<&str, &str> = HashMap::new();
        for (k, v) in keys.iter().zip(vals.iter()) {
            fmt.insert(*k, *v);
        }
        let gt = fmt.get("GT").copied().unwrap_or("");
        if gt != "0|1" && gt != "1|0" {
            continue;
        }
        let ps = match fmt.get("PS").copied() {
            Some(v) if v != "." => match v.parse::<i32>() {
                Ok(x) => x,
                Err(_) => continue,
            },
            _ => continue,
        };
        rna_vcfs.entry(ps).or_default().push((chr, pos));
    }
    rna_vcfs
}

struct AnalyzeContext<'a> {
    reads_positions: &'a HashMap<String, (u32, u32)>,
    reads_tags: &'a HashMap<String, ReadTag>,
    reads_exons: &'a HashMap<String, Vec<(u32, u32)>>,
    reads_introns: &'a HashMap<String, Vec<(u32, u32, bool)>>,
    dna_vcfs: Option<&'a HashSet<String>>,
    rna_vcfs: Option<&'a HashMap<i32, Vec<(String, u32)>>>,
}

#[allow(clippy::too_many_arguments)]
fn analyze_gene(
    gene_name: &str,
    gene_strand: &str,
    annotation_exons: &TranscriptRegions,
    annotation_junctions: &TranscriptRegions,
    gene_region: &GeneRegion,
    gene_reads: &[String],
    min_count: usize,
    cluster_with_exons: bool,
    ctx: &AnalyzeContext,
) -> Vec<AseEvent> {
    let valid_read_names: HashSet<String> = gene_reads
        .iter()
        .filter(|r| ctx.reads_tags.contains_key(*r))
        .cloned()
        .collect();
    let phased_read_names: Vec<String> = valid_read_names
        .into_iter()
        .filter(|name| {
            ctx.reads_tags
                .get(name)
                .map(|t| t.hp.is_some())
                .unwrap_or(false)
        })
        .collect();

    let mut sub_reads_positions: HashMap<String, (u32, u32)> = HashMap::new();
    let mut sub_reads_tags: HashMap<String, ReadTag> = HashMap::new();
    let mut sub_reads_exons: HashMap<String, Vec<(u32, u32)>> = HashMap::new();
    let mut sub_reads_introns: HashMap<String, Vec<(u32, u32, bool)>> = HashMap::new();

    for name in &phased_read_names {
        if let (Some(p), Some(t), Some(e), Some(i)) = (
            ctx.reads_positions.get(name),
            ctx.reads_tags.get(name),
            ctx.reads_exons.get(name),
            ctx.reads_introns.get(name),
        ) {
            sub_reads_positions.insert(name.clone(), *p);
            sub_reads_tags.insert(name.clone(), t.clone());
            sub_reads_exons.insert(name.clone(), e.clone());
            sub_reads_introns.insert(name.clone(), i.clone());
        }
    }

    let chr = gene_region.chr.clone();

    let mut gene_junction_set: HashSet<(String, u32, u32)> = HashSet::new();
    for anno_junctions in annotation_junctions.values() {
        for anno_junc in anno_junctions {
            gene_junction_set.insert(anno_junc.clone());
        }
    }
    let mut gene_exon_set: HashSet<(String, u32, u32)> = HashSet::new();
    for anno_exons in annotation_exons.values() {
        for anno_exon in anno_exons {
            gene_exon_set.insert(anno_exon.clone());
        }
    }

    let junction_clusters = if !cluster_with_exons {
        cluster_junctions_connected_components(&sub_reads_introns, min_count)
    } else {
        cluster_junctions_exons_connected_components(
            &sub_reads_introns,
            &sub_reads_exons,
            min_count,
        )
    };

    let mut anno_exon_intervals: Vec<(u32, u32)> = Vec::new();
    for anno_exon in &gene_exon_set {
        anno_exon_intervals.push((anno_exon.1, anno_exon.2));
    }

    let mut reads_to_remove: HashSet<String> = HashSet::new();

    if let (Some(dna_vcfs), Some(rna_vcfs)) = (ctx.dna_vcfs, ctx.rna_vcfs) {
        for (qname, tag) in &sub_reads_tags {
            let mut overlapped_snps_cnt = 0;
            if let Some(ps) = tag.ps {
                if let Some(ps_vars) = rna_vcfs.get(&ps) {
                    for (vchr, vpos) in ps_vars {
                        if dna_vcfs.contains(&format!("{}:{}", vchr, vpos)) {
                            overlapped_snps_cnt += 1;
                        }
                    }
                }
            }
            if overlapped_snps_cnt == 0 {
                reads_to_remove.insert(qname.clone());
            }
        }
    }

    for (qname, read_exons) in &sub_reads_exons {
        let mut overlapped = false;
        for (es, ee) in read_exons {
            for (as_, ae) in &anno_exon_intervals {
                let left = (*es).max(*as_);
                let right = (*ee).min(*ae);
                if right >= left {
                    overlapped = true;
                    break;
                }
            }
            if overlapped {
                break;
            }
        }
        if !overlapped {
            reads_to_remove.insert(qname.clone());
        }
    }

    for q in reads_to_remove {
        sub_reads_positions.remove(&q);
        sub_reads_exons.remove(&q);
        sub_reads_introns.remove(&q);
        sub_reads_tags.remove(&q);
    }

    let mut gene_ase_events: Vec<AseEvent> = Vec::new();

    for junc_cluster in &junction_clusters {
        if junc_cluster.is_empty() {
            continue;
        }
        let mut junc_cluster_sorted = junc_cluster.clone();
        junc_cluster_sorted.sort_by_key(|(start, end, _)| (*start, *end));
        let junction_set = format!(
            "{}:{}-{}",
            chr, junc_cluster_sorted[0].0, junc_cluster_sorted[0].1
        );
        for (junction_start, junction_end, gt_ag_tag) in &junc_cluster_sorted {
            let novel = !gene_junction_set.contains(&(chr.clone(), *junction_start, *junction_end));
            let (absent_reads, present_reads) = check_absent_present(
                *junction_start,
                *junction_end,
                &sub_reads_positions,
                &sub_reads_introns,
            );
            let test = haplotype_event_test(&absent_reads, &present_reads, &sub_reads_tags);
            if let Some((phase_set, h1a, h1p, h2a, h2p, pvalue, sor)) = test {
                gene_ase_events.push(AseEvent {
                    chr: chr.clone(),
                    start: *junction_start,
                    end: *junction_end,
                    strand: gene_strand.to_string(),
                    junction_set: junction_set.clone(),
                    phase_set: phase_set.to_string(),
                    hap1_absent: h1a,
                    hap1_present: h1p,
                    hap2_absent: h2a,
                    hap2_present: h2p,
                    p_value: pvalue,
                    sor,
                    novel,
                    gt_ag_tag: *gt_ag_tag,
                    gene_name: gene_name.to_string(),
                });
            }
        }
    }

    gene_ase_events
}

#[allow(clippy::too_many_arguments)]
fn analyze(
    annotation_file: &str,
    bam_file: &str,
    reference_file: &str,
    output_prefix: &str,
    min_count: usize,
    gene_types: HashSet<String>,
    threads: usize,
    no_gtag: bool,
    min_junctions: usize,
    cluster_with_exons: bool,
    dna_vcfs: Option<HashSet<String>>,
    rna_vcfs: Option<HashMap<i32, Vec<(String, u32)>>>,
) {
    let (
        anno_gene_regions,
        anno_gene_names,
        anno_gene_strands,
        anno_exon_regions,
        anno_intron_regions,
    ) = get_gene_regions(annotation_file, &gene_types);
    let merged_genes_exons = merge_gene_exon_regions(&anno_exon_regions);
    let (span_trees, exon_trees) = build_gene_trees(&merged_genes_exons);
    let genome_dict = load_reference(reference_file)
        .unwrap_or_else(|e| panic!("failed loading reference {}: {}", reference_file, e));
    let read_bundle = load_reads(
        bam_file,
        &genome_dict,
        &span_trees,
        &exon_trees,
        no_gtag,
        min_junctions,
    );
    let gene_assigned_reads = transform_read_assignment(&read_bundle.read_assignment);
    let all_gene_assigned_reads = transform_read_assignment(&read_bundle.all_read_assignment);

    let mut cov_writer = BufWriter::new(
        File::create(format!("{}.gene_coverage.tsv", output_prefix))
            .unwrap_or_else(|e| panic!("failed to create gene coverage file: {}", e)),
    );
    writeln!(cov_writer, "#Gene_name\tChr\tStart\tEnd\tNum_reads").unwrap();
    for (gene_id, region) in &anno_gene_regions {
        let cov = all_gene_assigned_reads
            .get(gene_id)
            .map(|x| x.len())
            .unwrap_or(0);
        let gname = anno_gene_names
            .get(gene_id)
            .map(|x| x.as_str())
            .unwrap_or(".");
        writeln!(
            cov_writer,
            "{}\t{}\t{}\t{}\t{}",
            gname, region.chr, region.start, region.end, cov
        )
        .unwrap();
    }

    let mut gene_data_ids: Vec<String> = anno_gene_regions
        .iter()
        .filter_map(|(gene_id, region)| {
            if genome_dict.contains_key(&region.chr)
                && gene_assigned_reads
                    .get(gene_id)
                    .map(|x| !x.is_empty())
                    .unwrap_or(false)
            {
                Some(gene_id.clone())
            } else {
                None
            }
        })
        .collect();
    gene_data_ids.sort();
    println!("Total genes to be analyzed: {}", gene_data_ids.len());

    let ctx = AnalyzeContext {
        reads_positions: &read_bundle.reads_positions,
        reads_tags: &read_bundle.reads_tags,
        reads_exons: &read_bundle.reads_exons,
        reads_introns: &read_bundle.reads_introns,
        dna_vcfs: dna_vcfs.as_ref(),
        rna_vcfs: rna_vcfs.as_ref(),
    };

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads.max(1))
        .build()
        .unwrap();

    let all_gene_events: Vec<Vec<AseEvent>> = pool.install(|| {
        gene_data_ids
            .par_iter()
            .map(|gene_id| {
                let gname = anno_gene_names.get(gene_id).unwrap();
                let gstrand = anno_gene_strands.get(gene_id).unwrap();
                let gexons = anno_exon_regions.get(gene_id).unwrap();
                let empty_introns: TranscriptRegions = HashMap::new();
                let gintrons = anno_intron_regions.get(gene_id).unwrap_or(&empty_introns);
                let gregion = anno_gene_regions.get(gene_id).unwrap();
                let greads = gene_assigned_reads.get(gene_id).unwrap();
                analyze_gene(
                    gname,
                    gstrand,
                    gexons,
                    gintrons,
                    gregion,
                    greads,
                    min_count,
                    cluster_with_exons,
                    &ctx,
                )
            })
            .collect()
    });

    let mut all_ase_events: HashMap<(String, u32, u32), HashMap<String, AseEvent>> = HashMap::new();
    for events in all_gene_events {
        for event in events {
            let key = (event.chr.clone(), event.start, event.end);
            all_ase_events
                .entry(key)
                .or_default()
                .insert(event.gene_name.clone(), event);
        }
    }

    let mut junctions: Vec<((String, u32, u32), String)> = Vec::new();
    for (k, gmap) in &all_ase_events {
        for gname in gmap.keys() {
            junctions.push((k.clone(), gname.clone()));
        }
    }
    println!("Total junctions: {}", junctions.len());

    let mut pass_idx: Vec<usize> = Vec::new();
    let mut p_values: Vec<f64> = Vec::new();
    for (idx, (junc, gname)) in junctions.iter().enumerate() {
        let event = all_ase_events.get(junc).unwrap().get(gname).unwrap();
        let total = event.hap1_absent + event.hap1_present + event.hap2_absent + event.hap2_present;
        if total >= min_count {
            pass_idx.push(idx);
            p_values.push(event.p_value);
        }
    }
    println!(
        "number of junctions with at least {} reads: {}",
        min_count,
        pass_idx.len()
    );
    let adjusted = benjamini_hochberg(&p_values);

    let mut asj_writer = BufWriter::new(
        File::create(format!("{}.asj.tsv", output_prefix))
            .unwrap_or_else(|e| panic!("failed to create asj file: {}", e)),
    );
    writeln!(asj_writer, "{}", AseEvent::header()).unwrap();

    let mut asj_genes: HashMap<String, AseEvent> = HashMap::new();
    for (pi, idx) in pass_idx.iter().enumerate() {
        let (junc, gname) = &junctions[*idx];
        let event = all_ase_events
            .get_mut(junc)
            .unwrap()
            .get_mut(gname)
            .unwrap();
        event.p_value = adjusted[pi];
        writeln!(asj_writer, "{}", event.line()).unwrap();
        if !no_gtag && !event.gt_ag_tag {
            continue;
        }
        match asj_genes.get(gname) {
            Some(current) if !is_better_gene_event(event, current) => {}
            _ => {
                asj_genes.insert(gname.clone(), event.clone());
            }
        }
    }
    println!(
        "number of genes with allele-specific junctions: {}",
        asj_genes.len()
    );

    let mut gene_writer = BufWriter::new(
        File::create(format!("{}.asj_gene.tsv", output_prefix))
            .unwrap_or_else(|e| panic!("failed to create asj_gene file: {}", e)),
    );
    writeln!(gene_writer, "#Gene_name\tChr\tP_value\tSOR").unwrap();
    let mut gene_names: Vec<String> = asj_genes.keys().cloned().collect();
    gene_names.sort();
    for gname in gene_names {
        if let Some(event) = asj_genes.get(&gname) {
            writeln!(gene_writer, "{}\t{}\t{}\t{}", gname, event.chr, event.p_value, event.sor)
                .unwrap();
        }
    }
}

pub fn run_asj(args: AsjArgs) {
    let has_filter = args.dna_vcf.is_some() && args.rna_vcf.is_some();
    if args.dna_vcf.is_some() ^ args.rna_vcf.is_some() {
        panic!("--dna-vcf and --rna-vcf must be provided together");
    }

    let gene_types_vec = args
        .gene_types
        .as_ref()
        .cloned()
        .unwrap_or_else(|| vec!["protein_coding".to_string(), "lncRNA".to_string()]);
    let gene_types: HashSet<String> = if gene_types_vec
        .iter()
        .any(|x| x == "__LONGCALLR_ALL_GENE_TYPES__")
    {
        HashSet::new()
    } else {
        gene_types_vec.into_iter().collect()
    };
    if gene_types.is_empty() {
        eprintln!("Gene types: all");
    } else {
        let mut sorted: Vec<&String> = gene_types.iter().collect();
        sorted.sort();
        eprintln!("Gene types: {}", sorted.iter().map(|s| s.as_str()).collect::<Vec<_>>().join(", "));
    }
    let dna_vcfs = if has_filter {
        Some(load_dna_vcf(args.dna_vcf.as_ref().unwrap()))
    } else {
        None
    };
    let rna_vcfs = if has_filter {
        Some(load_longcallr_phased_vcf(args.rna_vcf.as_ref().unwrap()))
    } else {
        None
    };

    analyze(
        &args.annotation_file,
        &args.bam_file,
        &args.reference,
        &args.output_prefix,
        args.min_sup,
        gene_types,
        args.threads,
        args.no_gtag,
        args.min_junctions,
        args.cluster_with_exons,
        dna_vcfs,
        rna_vcfs,
    );
}
