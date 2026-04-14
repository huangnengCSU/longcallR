use flate2::read::MultiGzDecoder;
use rayon::prelude::*;
use rust_htslib::bam;
use rust_htslib::bam::record::{Aux, Cigar};
use rust_htslib::bam::Read;
use rust_lapper::{Interval, Lapper};
use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

#[derive(clap::Parser, Debug, Clone)]
pub struct AseArgs {
    /// Phased BAM file (with PS/HP tags)
    #[arg(short = 'b', long)]
    bam: String,

    /// LongcallR phased RNA VCF file
    #[arg(long)]
    vcf1: Option<String>,

    /// Whole-genome haplotype-phased DNA VCF file
    #[arg(long)]
    vcf2: Option<String>,

    /// DNA VCF file for ASE filtering
    #[arg(long)]
    vcf3: Option<String>,

    /// Annotation file (GTF/GFF3, optionally gzipped)
    #[arg(short = 'a', long)]
    annotation: String,

    /// Overdispersion parameter for Beta-Binomial test
    #[arg(short = 'd', long, default_value_t = 0.001)]
    overdispersion: f64,

    /// Output prefix
    #[arg(short = 'o', long)]
    output: String,

    /// Number of threads
    #[arg(short = 't', long, default_value_t = 1)]
    threads: usize,

    /// Gene types to analyze; pass --gene-types with no values to include all gene types [default: protein_coding lncRNA]
    #[arg(long, num_args(0..))]
    gene_types: Option<Vec<String>>,

    /// Minimum support reads for counting ASE
    #[arg(long, default_value_t = 10)]
    min_support: u32,
}

#[derive(Debug, Clone)]
struct GeneRegion {
    chr: String,
    start: u32,
    end: u32,
}

#[derive(Debug, Clone)]
struct GeneInfo {
    name: String,
    region: GeneRegion,
}

#[derive(Debug, Clone)]
struct WgVariant {
    pat: char,
    mat: char,
}

#[derive(Debug, Clone)]
struct RnaVariant {
    chr: String,
    pos: u32,
    dp: Option<u32>,
    af: Option<f64>,
}

#[derive(Debug, Clone)]
struct BasicResult {
    gene_name: String,
    chr: String,
    p_value: f64,
    ps: String,
    h1: u32,
    h2: u32,
}

#[derive(Debug, Clone)]
struct PatMatResult {
    gene_name: String,
    chr: String,
    p_value: f64,
    ps: String,
    h1: u32,
    h2: u32,
    h1_pat: u32,
    h1_mat: u32,
    h2_pat: u32,
    h2_mat: u32,
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

type GeneExons = HashMap<String, HashMap<String, Vec<(String, u32, u32)>>>;

fn parse_gene_regions(
    annotation_file: &str,
    gene_types: &HashSet<String>,
) -> (HashMap<String, GeneInfo>, GeneExons) {
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

    let mut gene_infos: HashMap<String, GeneInfo> = HashMap::new();
    let mut exon_regions: GeneExons = HashMap::new();

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
        let feature_type = parts[2];
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

        if feature_type == "gene" {
            let chr = parts[0].to_string();
            let start = parts[3].parse::<u32>().unwrap_or(0);
            let end = parts[4].parse::<u32>().unwrap_or(0);
            if start == 0 || end == 0 || end < start {
                continue;
            }
            let gene_name = attrs
                .get("gene_name")
                .cloned()
                .unwrap_or_else(|| ".".to_string());
            gene_infos.insert(
                gene_id,
                GeneInfo {
                    name: gene_name,
                    region: GeneRegion { chr, start, end },
                },
            );
        } else if feature_type == "exon" {
            let gene_type = attrs
                .get("gene_type")
                .or_else(|| attrs.get("gene_biotype"))
                .cloned()
                .unwrap_or_else(|| "".to_string());
            let tag = attrs.get("tag").cloned().unwrap_or_else(|| "".to_string());
            if (!gene_types.is_empty() && !gene_types.contains(&gene_type)) || tag.contains("readthrough") {
                continue;
            }
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

    (gene_infos, exon_regions)
}

fn merge_gene_exons(exon_regions: &GeneExons) -> HashMap<String, HashMap<String, Vec<(u32, u32)>>> {
    let mut merged: HashMap<String, HashMap<String, Vec<(u32, u32)>>> = HashMap::new();
    for (gene_id, transcripts) in exon_regions {
        let mut chr_set: HashSet<String> = HashSet::new();
        let mut all_exons: Vec<(String, u32, u32)> = Vec::new();
        for exons in transcripts.values() {
            for (chr, start, end) in exons {
                chr_set.insert(chr.clone());
                all_exons.push((chr.clone(), *start, *end));
            }
        }
        if chr_set.len() != 1 || all_exons.is_empty() {
            continue;
        }
        let chr = chr_set.into_iter().next().unwrap();
        let mut intervals: Vec<(u32, u32)> =
            all_exons.into_iter().map(|(_, s, e)| (s, e)).collect();
        intervals.sort_by_key(|x| x.0);
        let mut merged_intervals: Vec<(u32, u32)> = Vec::new();
        for (start, end) in intervals {
            if let Some(last) = merged_intervals.last_mut() {
                if start <= last.1 {
                    if end > last.1 {
                        last.1 = end;
                    }
                } else {
                    merged_intervals.push((start, end));
                }
            } else {
                merged_intervals.push((start, end));
            }
        }
        if !merged_intervals.is_empty() {
            merged
                .entry(chr)
                .or_default()
                .insert(gene_id.clone(), merged_intervals);
        }
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
        let mut chr_exons: HashMap<String, Lapper<u32, ()>> = HashMap::new();
        for (gene_id, exons) in gene_map {
            if exons.is_empty() {
                continue;
            }
            let gene_start = exons.first().unwrap().0;
            let gene_end = exons.last().unwrap().1;
            span_intervals.push(Interval {
                start: gene_start,
                stop: gene_end + 1,
                val: gene_id.clone(),
            });

            let mut exon_intervals: Vec<Interval<u32, ()>> = exons
                .iter()
                .map(|(s, e)| Interval {
                    start: *s,
                    stop: e + 1,
                    val: (),
                })
                .collect();
            exon_intervals.sort_by_key(|i| i.start);
            chr_exons.insert(gene_id.clone(), Lapper::new(exon_intervals));
        }
        span_intervals.sort_by_key(|i| i.start);
        span_trees.insert(chr.clone(), Lapper::new(span_intervals));
        exon_trees.insert(chr.clone(), chr_exons);
    }

    (span_trees, exon_trees)
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

fn overlap_len(splice_start: u32, splice_end: u32, exon: &Interval<u32, ()>) -> u32 {
    let left = splice_start.max(exon.start);
    let right = (splice_end + 1).min(exon.stop);
    right.saturating_sub(left)
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

fn assign_reads_to_gene(
    bam_file: &str,
    span_trees: &GeneSpanTrees,
    exon_trees: &GeneExonTrees,
) -> (
    HashMap<String, String>,
    HashMap<String, (Option<i32>, Option<i32>)>,
) {
    let mut read_assignment: HashMap<String, String> = HashMap::new();
    let mut read_tags: HashMap<String, (Option<i32>, Option<i32>)> = HashMap::new();

    let mut bam = bam::Reader::from_path(bam_file)
        .unwrap_or_else(|e| panic!("failed to open BAM {}: {}", bam_file, e));
    let header = bam.header().to_owned();
    for rec in bam.records() {
        let record = rec.unwrap();
        if record.is_unmapped() || record.tid() < 0 {
            continue;
        }
        let chr = String::from_utf8_lossy(header.tid2name(record.tid() as u32)).to_string();
        let span_tree = match span_trees.get(&chr) {
            Some(t) => t,
            None => continue,
        };
        let chr_exons = match exon_trees.get(&chr) {
            Some(v) => v,
            None => continue,
        };

        let start_pos = record.pos() as u32;
        let end_pos = record.cigar().end_pos() as u32;
        if end_pos <= start_pos {
            continue;
        }

        let candidate_genes: Vec<(String, u32, u32)> = span_tree
            .find(start_pos + 1, end_pos + 1)
            .map(|x| (x.val.clone(), x.start, x.stop))
            .collect();
        if candidate_genes.is_empty() {
            continue;
        }

        let splice_regions = parse_splice_regions(&record);

        let mut best_gene: Option<String> = None;
        let mut best_overlap: u32 = 0;
        let mut best_start: u32 = 0;
        let mut best_span: u32 = u32::MAX;
        for (gene_id, gene_start, gene_stop) in candidate_genes {
            let exon_tree = match chr_exons.get(&gene_id) {
                Some(t) => t,
                None => continue,
            };
            let mut total_overlap: u32 = 0;
            for (sp_start, sp_end) in &splice_regions {
                let mut this_overlap = 0;
                for exon in exon_tree.find(*sp_start, sp_end + 1) {
                    this_overlap += overlap_len(*sp_start, *sp_end, exon);
                }
                total_overlap += this_overlap;
            }
            if total_overlap == 0 {
                continue;
            }

            let gene_span = gene_stop.saturating_sub(gene_start);
            let is_better = if total_overlap > best_overlap {
                true
            } else if total_overlap < best_overlap {
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
                best_overlap = total_overlap;
                best_start = gene_start;
                best_span = gene_span;
                best_gene = Some(gene_id);
            }
        }

        if let Some(gene_id) = best_gene {
            let read_name = String::from_utf8_lossy(record.qname()).to_string();
            read_assignment.insert(read_name.clone(), gene_id);
            let ps = record.aux(b"PS").ok().and_then(aux_to_i32);
            let hp = record.aux(b"HP").ok().and_then(aux_to_i32);
            read_tags.insert(read_name, (ps, hp));
        }
    }

    (read_assignment, read_tags)
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

fn convert_mu_rho_to_alpha_beta(mu: f64, rho: f64) -> (f64, f64) {
    let safe_rho = if rho <= 0.0 { 1e-9 } else { rho };
    let phi = (1.0 - safe_rho) / safe_rho - 1.0;
    (mu * phi, (1.0 - mu) * phi)
}

fn beta_binomial_log_pmf_series(n: u32, alpha: f64, beta: f64) -> Vec<f64> {
    let n_usize = n as usize;
    let n_f64 = n as f64;
    let mut log_pmf = vec![0.0; n_usize + 1];
    log_pmf[0] = 0.0;

    for k in 0..n_usize {
        let kf = k as f64;
        let delta = (n_f64 - kf).ln() - (kf + 1.0).ln() + (kf + alpha).ln()
            - (n_f64 - kf - 1.0 + beta).ln();
        log_pmf[k + 1] = log_pmf[k] + delta;
    }

    log_pmf
}

fn beta_binomial_p_value(k_obs: u32, n: u32, mu: f64, rho: f64) -> f64 {
    if k_obs > n {
        return 1.0;
    }
    let (alpha, beta) = convert_mu_rho_to_alpha_beta(mu, rho);
    let log_pmf = beta_binomial_log_pmf_series(n, alpha, beta);

    let max_log = log_pmf
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let pmf_scaled: Vec<f64> = log_pmf.iter().map(|lp| (lp - max_log).exp()).collect();
    let sum_scaled: f64 = pmf_scaled.iter().sum();
    if !sum_scaled.is_finite() || sum_scaled <= 0.0 {
        return 1.0;
    }

    let p_obs = pmf_scaled[k_obs as usize] / sum_scaled;
    let mut p_value = 0.0;
    for p_scaled in pmf_scaled {
        let p = p_scaled / sum_scaled;
        if p <= p_obs + 1e-15 {
            p_value += p;
        }
    }
    p_value.clamp(0.0, 1.0)
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

fn select_most_supported_ps(
    assigned_reads: &[String],
    read_tags: &HashMap<String, (Option<i32>, Option<i32>)>,
) -> Option<(i32, u32, u32)> {
    let mut counts: HashMap<i32, (u32, u32)> = HashMap::new();
    for read_name in assigned_reads {
        if let Some((Some(ps), Some(hp))) = read_tags.get(read_name) {
            if *hp == 1 {
                counts.entry(*ps).or_insert((0, 0)).0 += 1;
            } else if *hp == 2 {
                counts.entry(*ps).or_insert((0, 0)).1 += 1;
            }
        }
    }

    counts
        .into_iter()
        .max_by_key(|(_, (h1, h2))| h1 + h2)
        .map(|(ps, (h1, h2))| (ps, h1, h2))
}

fn calculate_basic_result(
    gene: &GeneInfo,
    assigned_reads: &[String],
    read_tags: &HashMap<String, (Option<i32>, Option<i32>)>,
    min_support: u32,
    overdispersion: f64,
) -> BasicResult {
    if let Some((ps, h1, h2)) = select_most_supported_ps(assigned_reads, read_tags) {
        if h1 + h2 < min_support {
            return BasicResult {
                gene_name: gene.name.clone(),
                chr: gene.region.chr.clone(),
                p_value: 1.0,
                ps: ps.to_string(),
                h1: 0,
                h2: 0,
            };
        }
        let p = beta_binomial_p_value(h1, h1 + h2, 0.5, overdispersion);
        BasicResult {
            gene_name: gene.name.clone(),
            chr: gene.region.chr.clone(),
            p_value: p,
            ps: ps.to_string(),
            h1,
            h2,
        }
    } else {
        BasicResult {
            gene_name: gene.name.clone(),
            chr: gene.region.chr.clone(),
            p_value: 1.0,
            ps: ".".to_string(),
            h1: 0,
            h2: 0,
        }
    }
}

fn load_longcallr_vcf(vcf_file: &str, with_dp_af: bool) -> HashMap<i32, Vec<RnaVariant>> {
    let mut rna_vcfs: HashMap<i32, Vec<RnaVariant>> = HashMap::new();
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
        let filter = fields[6];
        if !filter.split(';').any(|x| x == "PASS") {
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

        let (dp, af) = if with_dp_af {
            let dp = fmt.get("DP").and_then(|x| x.parse::<u32>().ok());
            let af = fmt
                .get("AF")
                .and_then(|x| x.split(',').next())
                .and_then(|x| x.parse::<f64>().ok());
            if dp.is_none() || af.is_none() {
                continue;
            }
            if dp.unwrap() == 0 || af.unwrap().is_nan() {
                continue;
            }
            (dp, af)
        } else {
            (None, None)
        };

        rna_vcfs
            .entry(ps)
            .or_default()
            .push(RnaVariant { chr, pos, dp, af });
    }
    rna_vcfs
}

fn load_wg_phased_vcf(vcf_file: &str) -> HashMap<String, WgVariant> {
    let mut wg_vcfs: HashMap<String, WgVariant> = HashMap::new();
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
        let (pat, mat) = if gt == "0|1" {
            (
                alt.chars().next().unwrap().to_ascii_uppercase(),
                ref_allele.chars().next().unwrap().to_ascii_uppercase(),
            )
        } else if gt == "1|0" {
            (
                ref_allele.chars().next().unwrap().to_ascii_uppercase(),
                alt.chars().next().unwrap().to_ascii_uppercase(),
            )
        } else {
            continue;
        };
        wg_vcfs.insert(format!("{}:{}", chr, pos), WgVariant { pat, mat });
    }
    wg_vcfs
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

fn calc_pat_mat_counts(
    bam_file: &str,
    chr: &str,
    start: u32,
    end: u32,
    ps_reads: &HashSet<String>,
    ps_variant_pos0: &HashSet<u32>,
    wg_vcfs: &HashMap<String, WgVariant>,
) -> HashMap<String, (u32, u32)> {
    const MIN_BASE_QUAL: u8 = 13;
    const MAX_PILEUP_DEPTH: u32 = 100_000;
    let mut counts: HashMap<String, (u32, u32)> = HashMap::new();
    let mut bam = bam::IndexedReader::from_path(bam_file)
        .unwrap_or_else(|e| panic!("failed to open indexed BAM {}: {}", bam_file, e));
    let tid = bam
        .header()
        .tid(chr.as_bytes())
        .unwrap_or_else(|| panic!("chromosome {} not found in BAM header", chr));
    bam.fetch((tid, (start - 1) as i64, end as i64))
        .unwrap_or_else(|e| panic!("failed to fetch region {}:{}-{}: {}", chr, start, end, e));

    let mut pileups = bam.pileup();
    pileups.set_max_depth(MAX_PILEUP_DEPTH);
    for p in pileups {
        let pileup = p.unwrap();
        let pos0 = pileup.pos() as u32;
        if !ps_variant_pos0.contains(&pos0) {
            continue;
        }
        let key = format!("{}:{}", chr, pos0 + 1);
        let variant = match wg_vcfs.get(&key) {
            Some(v) => v,
            None => continue,
        };

        for aln in pileup.alignments() {
            if aln.is_del() || aln.is_refskip() {
                continue;
            }
            let qpos = match aln.qpos() {
                Some(v) => v,
                None => continue,
            };
            let rec = aln.record();
            let read_name = String::from_utf8_lossy(rec.qname()).to_string();
            if !ps_reads.contains(&read_name) {
                continue;
            }
            let seq = rec.seq().as_bytes();
            if qpos >= seq.len() {
                continue;
            }
            let quals = rec.qual();
            if qpos >= quals.len() || quals[qpos] < MIN_BASE_QUAL {
                continue;
            }
            let base = (seq[qpos] as char).to_ascii_uppercase();
            let entry = counts.entry(read_name).or_insert((0, 0));
            if base == variant.pat {
                entry.0 += 1;
            } else if base == variant.mat {
                entry.1 += 1;
            }
        }
    }
    counts
}

fn calculate_patmat_result(
    bam_file: &str,
    gene: &GeneInfo,
    assigned_reads: &[String],
    read_tags: &HashMap<String, (Option<i32>, Option<i32>)>,
    rna_vcfs: &HashMap<i32, Vec<RnaVariant>>,
    wg_vcfs: &HashMap<String, WgVariant>,
    min_support: u32,
    overdispersion: f64,
) -> PatMatResult {
    let selected = select_most_supported_ps(assigned_reads, read_tags);
    if selected.is_none() {
        return PatMatResult {
            gene_name: gene.name.clone(),
            chr: gene.region.chr.clone(),
            p_value: 1.0,
            ps: ".".to_string(),
            h1: 0,
            h2: 0,
            h1_pat: 0,
            h1_mat: 0,
            h2_pat: 0,
            h2_mat: 0,
        };
    }

    let (ps, h1, h2) = selected.unwrap();
    if h1 + h2 < min_support {
        return PatMatResult {
            gene_name: gene.name.clone(),
            chr: gene.region.chr.clone(),
            p_value: 1.0,
            ps: ".".to_string(),
            h1: 0,
            h2: 0,
            h1_pat: 0,
            h1_mat: 0,
            h2_pat: 0,
            h2_mat: 0,
        };
    }

    let p_value = beta_binomial_p_value(h1, h1 + h2, 0.5, overdispersion);
    let ps_reads: Vec<String> = assigned_reads
        .iter()
        .filter_map(|r| {
            read_tags.get(r).and_then(|(rps, _)| {
                if rps == &Some(ps) {
                    Some(r.clone())
                } else {
                    None
                }
            })
        })
        .collect();
    let h1_reads: Vec<String> = ps_reads
        .iter()
        .filter_map(|r| {
            read_tags.get(r).and_then(|(_, hp)| {
                if hp == &Some(1) {
                    Some(r.clone())
                } else {
                    None
                }
            })
        })
        .collect();
    let h2_reads: Vec<String> = ps_reads
        .iter()
        .filter_map(|r| {
            read_tags.get(r).and_then(|(_, hp)| {
                if hp == &Some(2) {
                    Some(r.clone())
                } else {
                    None
                }
            })
        })
        .collect();
    let ps_reads_set: HashSet<String> = ps_reads.into_iter().collect();

    let ps_variant_pos0: HashSet<u32> = rna_vcfs
        .get(&ps)
        .unwrap_or(&Vec::new())
        .iter()
        .filter(|v| v.chr == gene.region.chr)
        .map(|v| v.pos - 1)
        .collect();

    let read_pat_mat_cnt = calc_pat_mat_counts(
        bam_file,
        &gene.region.chr,
        gene.region.start,
        gene.region.end,
        &ps_reads_set,
        &ps_variant_pos0,
        wg_vcfs,
    );

    let mut h1_pat = 0;
    let mut h1_mat = 0;
    for r in &h1_reads {
        if let Some((pat, mat)) = read_pat_mat_cnt.get(r) {
            if pat > mat {
                h1_pat += 1;
            } else if mat > pat {
                h1_mat += 1;
            }
        }
    }

    let mut h2_pat = 0;
    let mut h2_mat = 0;
    for r in &h2_reads {
        if let Some((pat, mat)) = read_pat_mat_cnt.get(r) {
            if pat > mat {
                h2_pat += 1;
            } else if mat > pat {
                h2_mat += 1;
            }
        }
    }

    PatMatResult {
        gene_name: gene.name.clone(),
        chr: gene.region.chr.clone(),
        p_value,
        ps: ps.to_string(),
        h1,
        h2,
        h1_pat,
        h1_mat,
        h2_pat,
        h2_mat,
    }
}

fn calculate_filtered_result(
    gene: &GeneInfo,
    assigned_reads: &[String],
    read_tags: &HashMap<String, (Option<i32>, Option<i32>)>,
    rna_vcfs: &HashMap<i32, Vec<RnaVariant>>,
    dna_vcfs: &HashSet<String>,
    min_support: u32,
    overdispersion: f64,
) -> BasicResult {
    if let Some((ps, h1, h2)) = select_most_supported_ps(assigned_reads, read_tags) {
        if h1 + h2 < min_support {
            return BasicResult {
                gene_name: gene.name.clone(),
                chr: gene.region.chr.clone(),
                p_value: 1.0,
                ps: ps.to_string(),
                h1: 0,
                h2: 0,
            };
        }
        let p = beta_binomial_p_value(h1, h1 + h2, 0.5, overdispersion);

        let mut overlapped_cnt = 0;
        if let Some(variants) = rna_vcfs.get(&ps) {
            for v in variants {
                let key = format!("{}:{}", v.chr, v.pos);
                if !dna_vcfs.contains(&key) {
                    continue;
                }
                let (dp, af) = match (v.dp, v.af) {
                    (Some(dp), Some(af)) => (dp, af),
                    _ => continue,
                };
                let mut alt_cnt = (dp as f64 * af).round() as i64;
                if alt_cnt < 0 {
                    alt_cnt = 0;
                } else if alt_cnt > dp as i64 {
                    alt_cnt = dp as i64;
                }
                let alt_cnt = alt_cnt as u32;
                let p_ase_allele = beta_binomial_p_value(alt_cnt, dp, 0.5, overdispersion);
                if dp >= min_support && p_ase_allele < 0.05 {
                    overlapped_cnt += 1;
                }
            }
        }
        if overlapped_cnt == 0 {
            return BasicResult {
                gene_name: gene.name.clone(),
                chr: gene.region.chr.clone(),
                p_value: 1.0,
                ps: ".".to_string(),
                h1: 0,
                h2: 0,
            };
        }

        BasicResult {
            gene_name: gene.name.clone(),
            chr: gene.region.chr.clone(),
            p_value: p,
            ps: ps.to_string(),
            h1,
            h2,
        }
    } else {
        BasicResult {
            gene_name: gene.name.clone(),
            chr: gene.region.chr.clone(),
            p_value: 1.0,
            ps: ".".to_string(),
            h1: 0,
            h2: 0,
        }
    }
}

fn write_basic_results(output_path: &str, results: &[BasicResult], min_support: u32) {
    let mut pass_idx: Vec<usize> = Vec::new();
    let mut p_values: Vec<f64> = Vec::new();
    for (idx, r) in results.iter().enumerate() {
        if r.h1 + r.h2 >= min_support {
            pass_idx.push(idx);
            p_values.push(r.p_value);
        }
    }
    println!("total number of genes: {}", results.len());
    println!(
        "number of genes with at least {} reads: {}",
        min_support,
        pass_idx.len()
    );
    let adjusted = benjamini_hochberg(&p_values);

    let f = File::create(output_path)
        .unwrap_or_else(|e| panic!("failed to create {}: {}", output_path, e));
    let mut writer = BufWriter::new(f);
    writeln!(writer, "#Gene_name\tChr\tPS\tH1\tH2\tP_value\tlogFC").unwrap();
    for (pi, idx) in pass_idx.iter().enumerate() {
        let r = &results[*idx];
        let logfc = ((r.h1 as f64 + 1.0) / (r.h2 as f64 + 1.0)).log2();
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.gene_name, r.chr, r.ps, r.h1, r.h2, adjusted[pi], logfc
        )
        .unwrap();
    }
}

fn write_patmat_results(output_path: &str, results: &[PatMatResult], min_support: u32) {
    let mut pass_idx: Vec<usize> = Vec::new();
    let mut p_values: Vec<f64> = Vec::new();
    for (idx, r) in results.iter().enumerate() {
        if r.h1 + r.h2 >= min_support {
            pass_idx.push(idx);
            p_values.push(r.p_value);
        }
    }
    println!("total number of genes: {}", results.len());
    println!(
        "number of genes with at least {} reads: {}",
        min_support,
        pass_idx.len()
    );
    let adjusted = benjamini_hochberg(&p_values);

    let f = File::create(output_path)
        .unwrap_or_else(|e| panic!("failed to create {}: {}", output_path, e));
    let mut writer = BufWriter::new(f);
    writeln!(writer, "#Gene_name\tChr\tPS\tH1\tH2\tP_value\tH1_Paternal\tH1_Maternal\tH2_Paternal\tH2_Maternal\tlogFC").unwrap();
    for (pi, idx) in pass_idx.iter().enumerate() {
        let r = &results[*idx];
        let logfc = ((r.h1 as f64 + 1.0) / (r.h2 as f64 + 1.0)).log2();
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.gene_name,
            r.chr,
            r.ps,
            r.h1,
            r.h2,
            adjusted[pi],
            r.h1_pat,
            r.h1_mat,
            r.h2_pat,
            r.h2_mat,
            logfc
        )
        .unwrap();
    }
}

fn run_basic(
    args: &AseArgs,
    gene_infos: &HashMap<String, GeneInfo>,
    gene_assigned_reads: &HashMap<String, Vec<String>>,
    read_tags: &HashMap<String, (Option<i32>, Option<i32>)>,
) {
    let mut gene_ids: Vec<String> = gene_infos
        .keys()
        .filter(|g| gene_assigned_reads.contains_key(*g))
        .cloned()
        .collect();
    gene_ids.sort();

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()
        .unwrap();
    let results: Vec<BasicResult> = pool.install(|| {
        gene_ids
            .par_iter()
            .map(|gene_id| {
                let gene = gene_infos.get(gene_id).unwrap();
                let reads = gene_assigned_reads.get(gene_id).unwrap();
                calculate_basic_result(
                    gene,
                    reads,
                    read_tags,
                    args.min_support,
                    args.overdispersion,
                )
            })
            .collect()
    });
    write_basic_results(
        &(args.output.clone() + ".ase.tsv"),
        &results,
        args.min_support,
    );
}

fn run_filter(
    args: &AseArgs,
    gene_infos: &HashMap<String, GeneInfo>,
    gene_assigned_reads: &HashMap<String, Vec<String>>,
    read_tags: &HashMap<String, (Option<i32>, Option<i32>)>,
    vcf1: &str,
    vcf3: &str,
) {
    let rna_vcfs = load_longcallr_vcf(vcf1, true);
    let dna_vcfs = load_dna_vcf(vcf3);

    let mut gene_ids: Vec<String> = gene_infos
        .keys()
        .filter(|g| gene_assigned_reads.contains_key(*g))
        .cloned()
        .collect();
    gene_ids.sort();

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()
        .unwrap();
    let results: Vec<BasicResult> = pool.install(|| {
        gene_ids
            .par_iter()
            .map(|gene_id| {
                let gene = gene_infos.get(gene_id).unwrap();
                let reads = gene_assigned_reads.get(gene_id).unwrap();
                calculate_filtered_result(
                    gene,
                    reads,
                    read_tags,
                    &rna_vcfs,
                    &dna_vcfs,
                    args.min_support,
                    args.overdispersion,
                )
            })
            .collect()
    });
    write_basic_results(
        &(args.output.clone() + ".filter_ase.tsv"),
        &results,
        args.min_support,
    );
}

fn run_patmat(
    args: &AseArgs,
    gene_infos: &HashMap<String, GeneInfo>,
    gene_assigned_reads: &HashMap<String, Vec<String>>,
    read_tags: &HashMap<String, (Option<i32>, Option<i32>)>,
    vcf1: &str,
    vcf2: &str,
) {
    let rna_vcfs = load_longcallr_vcf(vcf1, false);
    let wg_vcfs = load_wg_phased_vcf(vcf2);

    let mut gene_ids: Vec<String> = gene_infos
        .keys()
        .filter(|g| gene_assigned_reads.contains_key(*g))
        .cloned()
        .collect();
    gene_ids.sort();

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()
        .unwrap();
    let results: Vec<PatMatResult> = pool.install(|| {
        gene_ids
            .par_iter()
            .map(|gene_id| {
                let gene = gene_infos.get(gene_id).unwrap();
                let reads = gene_assigned_reads.get(gene_id).unwrap();
                calculate_patmat_result(
                    &args.bam,
                    gene,
                    reads,
                    read_tags,
                    &rna_vcfs,
                    &wg_vcfs,
                    args.min_support,
                    args.overdispersion,
                )
            })
            .collect()
    });
    write_patmat_results(
        &(args.output.clone() + ".patmat_ase.tsv"),
        &results,
        args.min_support,
    );
}

pub fn run_ase(args: AseArgs) {
    if args.vcf1.is_none() && (args.vcf2.is_some() || args.vcf3.is_some()) {
        panic!("--vcf2/--vcf3 requires --vcf1");
    }
    if args.vcf2.is_some() && args.vcf3.is_some() {
        panic!("--vcf2 and --vcf3 are mutually exclusive");
    }

    let gene_types_vec = args
        .gene_types
        .as_ref()
        .cloned()
        .unwrap_or_else(|| vec!["protein_coding".to_string(), "lncRNA".to_string()]);
    let gene_type_set: HashSet<String> = if gene_types_vec
        .iter()
        .any(|x| x == "__LONGCALLR_ALL_GENE_TYPES__")
    {
        HashSet::new()
    } else {
        gene_types_vec.into_iter().collect()
    };
    if gene_type_set.is_empty() {
        eprintln!("Gene types: all");
    } else {
        let mut sorted: Vec<&String> = gene_type_set.iter().collect();
        sorted.sort();
        eprintln!("Gene types: {}", sorted.iter().map(|s| s.as_str()).collect::<Vec<_>>().join(", "));
    }
    let (gene_infos, exon_regions) = parse_gene_regions(&args.annotation, &gene_type_set);
    let merged_exons = merge_gene_exons(&exon_regions);
    let (span_trees, exon_trees) = build_gene_trees(&merged_exons);
    let (read_assignment, read_tags) = assign_reads_to_gene(&args.bam, &span_trees, &exon_trees);
    let gene_assigned_reads = transform_read_assignment(&read_assignment);

    if let (Some(vcf1), Some(vcf2)) = (args.vcf1.as_deref(), args.vcf2.as_deref()) {
        run_patmat(
            &args,
            &gene_infos,
            &gene_assigned_reads,
            &read_tags,
            vcf1,
            vcf2,
        );
    } else if let (Some(vcf1), Some(vcf3)) = (args.vcf1.as_deref(), args.vcf3.as_deref()) {
        run_filter(
            &args,
            &gene_infos,
            &gene_assigned_reads,
            &read_tags,
            vcf1,
            vcf3,
        );
    } else {
        run_basic(&args, &gene_infos, &gene_assigned_reads, &read_tags);
    }
}
