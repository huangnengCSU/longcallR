use std::collections::{HashMap, HashSet, VecDeque};
use std::fs;
use std::fs::File;
use std::io::Write;
use std::sync::Mutex;

use rayon::iter::IntoParallelRefIterator;
use rayon::prelude::*;
use rust_htslib::{bam, bam::ext::BamRecordExtensions, bam::Format, bam::Read, bam::record::Aux};
use rust_lapper::Interval;

use crate::exon::{Exon, exon_cluster};
use crate::Platform;
use crate::snpfrags::SNPFrag;
use crate::util::{load_reference, parse_fai, Profile, Region};

pub fn multithread_phase_haplotag(
    bam_file: String,
    ref_file: String,
    vcf_file: String,
    phased_bam_file: String,
    thread_size: usize,
    isolated_regions: Vec<Region>,
    exon_regions: HashMap<String, Vec<Interval<usize, u8>>>,
    genotype_only: bool,
    platform: &Platform,
    max_iters: i32,
    min_mapq: u8,
    min_baseq: u8,
    min_allele_freq: f32,
    hetvar_high_frac_cutoff: f32,
    min_allele_freq_include_intron: f32,
    min_qual_for_candidate: u32,
    use_strand_bias: bool,
    strand_bias_threshold: f32,
    cover_strand_bias_threshold: f32,
    min_depth: u32,
    max_depth: u32,
    min_read_length: usize,
    distance_to_splicing_site: u32,
    window_size: u32,
    distance_to_read_end: u32,
    polya_tail_len: u32,
    dense_win_size: u32,
    min_dense_cnt: u32,
    min_linkers: u32,
    min_phase_score: f32,
    max_enum_snps: usize,
    random_flip_fraction: f32,
    read_assignment_cutoff: f64,
    imbalance_allele_expression_cutoff: f32,
    no_bam_output: bool,
    haplotype_bam_output: bool,
    output_read_assignment: bool,
    haplotype_specific_exon: bool,
    min_sup_haplotype_exon: u32,
    somatic_allele_frac_cutoff: f32,
    somatic_allele_cnt_cutoff: u32,
) {
    let pool = rayon::ThreadPoolBuilder::new().num_threads(thread_size).build().unwrap();
    let vcf_records_queue = Mutex::new(VecDeque::new());
    let read_haplotag1_queue = Mutex::new(VecDeque::new());
    let read_haplotag2_queue = Mutex::new(VecDeque::new());
    let read_haplotag_queue = Mutex::new(VecDeque::new());
    let read_phaseset_queue = Mutex::new(VecDeque::new());
    let haplotype_exon_queue = Mutex::new(VecDeque::new());
    let ref_seqs = load_reference(ref_file.clone());
    let fai_path = ref_file + ".fai";
    if fs::metadata(&fai_path).is_err() {
        panic!("Reference index file .fai does not exist.");
    }
    let contig_lengths = parse_fai(fai_path.as_str());
    let mut contig_order = Vec::new();
    for (k, _) in contig_lengths.iter() {
        contig_order.push(k.clone());
    }

    pool.install(|| {
        isolated_regions.par_iter().for_each(|reg| {
            let mut profile = Profile::default();
            let ref_seq = ref_seqs.get(&reg.chr).unwrap();
            let mut exon_region_vec = Vec::new();
            if !reg.gene_id.is_none() {
                let gene_id_field = reg.gene_id.clone().unwrap();
                for gene_id in gene_id_field.split(",").collect::<Vec<&str>>() {
                    if exon_regions.contains_key(gene_id) {
                        exon_region_vec.extend(exon_regions.get(gene_id).unwrap().clone());
                    }
                }
                if exon_region_vec.len() == 0 {
                    // this region is done, no exon region covered
                    return;
                }
            }
            profile.init_with_pileup(
                &bam_file.as_str(),
                &reg,
                ref_seq,
                platform,
                min_mapq,
                min_baseq,
                min_read_length,
                min_depth,
                max_depth,
                distance_to_read_end,
                polya_tail_len,
            );
            let mut snpfrag = SNPFrag::default();
            snpfrag.region = reg.clone();
            snpfrag.min_linkers = min_linkers;
            snpfrag.get_candidate_snps(
                &profile,
                &platform,
                exon_region_vec,
                min_allele_freq,
                hetvar_high_frac_cutoff,
                min_allele_freq_include_intron,
                min_depth,
                max_depth,
                min_baseq,
                use_strand_bias,
                strand_bias_threshold,
                cover_strand_bias_threshold,
                distance_to_splicing_site,
                window_size,
                dense_win_size,
                min_dense_cnt,
                somatic_allele_frac_cutoff,
                somatic_allele_cnt_cutoff,
                genotype_only,
            );
            // TODO: for very high depth region, down-sampling the reads
            snpfrag.get_fragments(&bam_file, &reg);
            if genotype_only {
                // without phasing
                let vcf_records = snpfrag.output_vcf(min_qual_for_candidate);
                {
                    let mut queue = vcf_records_queue.lock().unwrap();
                    for rd in vcf_records.iter() {
                        queue.push_back(rd.clone());
                    }
                }
            } else {
                if snpfrag.candidate_snps.len() >= 0 {
                    unsafe {
                        snpfrag.init_haplotypes();
                    }
                    unsafe {
                        snpfrag.init_assignment();
                    }
                    snpfrag.phase(max_enum_snps, random_flip_fraction, max_iters);
                    let read_assignments = snpfrag.assign_reads_haplotype(read_assignment_cutoff);
                    snpfrag.assign_snp_haplotype(min_phase_score);
                    // snpfrag.assign_het_var_haplotype(min_phase_score, somatic_allele_frac_cutoff, somatic_allele_cnt_cutoff);
                    // snpfrag.eval_low_frac_het_var_phase(min_phase_score, somatic_allele_frac_cutoff, somatic_allele_cnt_cutoff);
                    snpfrag.eval_rna_edit_var_phase(min_phase_score);
                    // snpfrag.eval_hom_var_phase(min_phase_score);
                    // assign phased fragments to somatic mutations and detect condifent somatic mutations
                    // println!("somatic: {}", snpfrag.somatic_snps.len());
                    // snpfrag.detect_somatic_by_het(&bam_file.as_str(), &reg);
                    // snpfrag.phase_ase_hete_snps(max_enum_snps, random_flip_fraction, max_iters);
                    // assign reads to haplotypes, filter reads having conflicted ase snps and heterozygous snps
                    // let read_assignments_ase = snpfrag.assign_reads_ase(read_assignment_cutoff);
                    // snpfrag.rescue_ase_snps();
                    // snpfrag.rescue_ase_snps_v2(ase_allele_cnt_cutoff, ase_ps_count_cutoff, ase_ps_cutoff);

                    // merge read_assignments and read_assignments_ase, read_assignments_ase first, then read_assignments
                    // let mut merge_reads_assignments = read_assignments.clone();
                    // for (k, v) in read_assignments.iter() {
                    //     if !read_assignments_ase.contains_key(k) {
                    //         merge_reads_assignments.insert(k.clone(), v.clone());
                    //     }
                    // }
                    let phase_sets = snpfrag.assign_phase_set();

                    let mut haplotype_exons: Vec<(Exon, i32, i32)> = Vec::new();
                    {
                        if haplotype_bam_output || haplotype_specific_exon {
                            let mut queue1 = read_haplotag1_queue.lock().unwrap();
                            let mut queue2 = read_haplotag2_queue.lock().unwrap();
                            let mut hap1_read_count = 0;
                            let mut hap2_read_count = 0;
                            let mut haplotype_read_count_pass = false;
                            for a in read_assignments.iter() {
                                if *a.1 == 1 {
                                    hap1_read_count += 1;
                                } else if *a.1 == 2 {
                                    hap2_read_count += 1;
                                }
                                if hap1_read_count >= 10 && hap2_read_count >= 10 {
                                    haplotype_read_count_pass = true;
                                    break;
                                }
                            }
                            if haplotype_bam_output && haplotype_read_count_pass {
                                for a in read_assignments.iter() {
                                    if *a.1 == 1 {
                                        queue1.push_back(a.0.clone());
                                    } else if *a.1 == 2 {
                                        queue2.push_back(a.0.clone());
                                    }
                                }
                            }
                            if haplotype_specific_exon && haplotype_read_count_pass {
                                let mut hap1_exons: Vec<Exon> = Vec::new();
                                let mut hap2_exons: Vec<Exon> = Vec::new();
                                let mut hap1_smallest_start = 0;
                                let mut hap1_largest_end = 0;
                                let mut hap2_smallest_start = 0;
                                let mut hap2_largest_end = 0;
                                // collect exons
                                for frag in snpfrag.fragments.iter() {
                                    if frag.assignment == 1 {
                                        for e in frag.exons.iter() {
                                            hap1_exons.push(e.clone());
                                            if hap1_smallest_start == 0 || e.start < hap1_smallest_start {
                                                hap1_smallest_start = e.start;
                                            }
                                            if e.end > hap1_largest_end {
                                                hap1_largest_end = e.end;
                                            }
                                        }
                                    } else if frag.assignment == 2 {
                                        for e in frag.exons.iter() {
                                            hap2_exons.push(e.clone());
                                            if hap2_smallest_start == 0 || e.start < hap2_smallest_start {
                                                hap2_smallest_start = e.start;
                                            }
                                            if e.end > hap2_largest_end {
                                                hap2_largest_end = e.end;
                                            }
                                        }
                                    }
                                }
                                // consensus exons
                                let hap1_consensus_exons = exon_cluster(
                                    hap1_exons.clone(),
                                    hap1_smallest_start,
                                    hap1_largest_end,
                                    0,
                                );
                                let hap2_consensus_exons = exon_cluster(
                                    hap2_exons.clone(),
                                    hap2_smallest_start,
                                    hap2_largest_end,
                                    0,
                                );
                                let mut combined_consensus_exons: HashMap<Exon, (i32, i32)> = HashMap::new();
                                for (e, v) in hap1_consensus_exons.iter() {
                                    if combined_consensus_exons.contains_key(e) {
                                        let (c1, c2) = combined_consensus_exons.get_mut(e).unwrap();
                                        *c1 += v.len() as i32;
                                    } else {
                                        combined_consensus_exons.insert(e.clone(), (v.len() as i32, 0));
                                    }
                                }
                                for (e, v) in hap2_consensus_exons.iter() {
                                    if combined_consensus_exons.contains_key(e) {
                                        let (c1, c2) = combined_consensus_exons.get_mut(e).unwrap();
                                        *c2 += v.len() as i32;
                                    } else {
                                        combined_consensus_exons.insert(e.clone(), (0, v.len() as i32));
                                    }
                                }
                                for (e, counts) in combined_consensus_exons.iter() {
                                    if counts.0 * counts.1 == 0 && counts.0 + counts.1 >= min_sup_haplotype_exon as i32 {
                                        if e.state == 1 {
                                            // println!("exon1: {}:{}-{}, hap1:{:?}, hap2:{:?}", e.chr, e.start + 1, e.end + 1, counts.0, counts.1);
                                            haplotype_exons.push((e.clone(), counts.0, counts.1));
                                        }
                                        if e.state == 0 {
                                            let mut start_sum = 0;
                                            let mut start_mean = 0;
                                            let mut exon_cnt = 0;
                                            if counts.0 > 0 {
                                                for ex in hap1_consensus_exons.get(e).unwrap().iter() {
                                                    start_sum += ex.start;
                                                    exon_cnt += 1;
                                                }
                                                start_mean = start_sum / exon_cnt;
                                            } else {
                                                for ex in hap2_consensus_exons.get(e).unwrap().iter() {
                                                    start_sum += ex.start;
                                                    exon_cnt += 1;
                                                }
                                                start_mean = start_sum / exon_cnt;
                                            }
                                            // println!("exon0: {}:{}-{}, hap1:{:?}, hap2:{:?}", e.chr, start_mean + 1, e.end + 1, counts.0, counts.1);
                                            let mut ec = e.clone();
                                            ec.start = start_mean;
                                            haplotype_exons.push((ec, counts.0, counts.1));
                                        }
                                        if e.state == 2 {
                                            let mut end_sum = 0;
                                            let mut end_mean = 0;
                                            let mut exon_cnt = 0;
                                            if counts.0 > 0 {
                                                for ex in hap1_consensus_exons.get(e).unwrap().iter() {
                                                    end_sum += ex.end;
                                                    exon_cnt += 1;
                                                }
                                                end_mean = end_sum / exon_cnt;
                                            } else {
                                                for ex in hap2_consensus_exons.get(e).unwrap().iter() {
                                                    end_sum += ex.end;
                                                    exon_cnt += 1;
                                                }
                                                end_mean = end_sum / exon_cnt;
                                            }
                                            // println!("exon2: {}:{}-{}, hap1:{:?}, hap2:{:?}", e.chr, e.start + 1, end_mean + 1, counts.0, counts.1);
                                            let mut ec = e.clone();
                                            ec.end = end_mean;
                                            haplotype_exons.push((ec, counts.0, counts.1));
                                        }
                                        if e.state == 3 {
                                            let relaxed_start = [e.start - 20, e.start + 20];
                                            let relaxed_end = [e.end - 20, e.end + 20];
                                            let mut unique_flag = true;
                                            if counts.0 > 0 {
                                                for ex in hap2_consensus_exons.iter() {
                                                    if ex.0.state != 3 {
                                                        continue;
                                                    }
                                                    if ex.0.start >= relaxed_start[0] && ex.0.start <= relaxed_start[1] && ex.0.end >= relaxed_end[0] && ex.0.end <= relaxed_end[1] {
                                                        // highly overlapped, not unique
                                                        unique_flag = false;
                                                        break;
                                                    }
                                                }
                                            } else {
                                                for ex in hap1_consensus_exons.iter() {
                                                    if ex.0.state != 3 {
                                                        continue;
                                                    }
                                                    if ex.0.start >= relaxed_start[0] && ex.0.start <= relaxed_start[1] && ex.0.end >= relaxed_end[0] && ex.0.end <= relaxed_end[1] {
                                                        // highly overlapped, not unique
                                                        unique_flag = false;
                                                        break;
                                                    }
                                                }
                                            }
                                            if unique_flag {
                                                // println!("exon3: {}:{}-{}, hap1:{:?}, hap2:{:?}", e.chr, e.start + 1, e.end + 1, counts.0, counts.1);
                                                haplotype_exons.push((
                                                    e.clone(),
                                                    counts.0,
                                                    counts.1,
                                                ));
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        // if !no_bam_output {
                        //     let mut queue = read_haplotag_queue.lock().unwrap();
                        //     for a in read_assignments.iter() {
                        //         queue.push_back((a.0.clone(), a.1.clone()));
                        //     }
                        // }

                        // output assignment both for ase snps and heterozygous snps
                        if !no_bam_output {
                            let mut queue = read_haplotag_queue.lock().unwrap();
                            for a in read_assignments.iter() {
                                queue.push_back((a.0.clone(), a.1.clone()));
                            }
                            let mut queue = read_phaseset_queue.lock().unwrap();
                            for a in phase_sets.iter() {
                                queue.push_back((a.0.clone(), a.1.clone()));
                            }
                        }
                    }
                    if haplotype_specific_exon {
                        let mut queue = haplotype_exon_queue.lock().unwrap();
                        for e in haplotype_exons.iter() {
                            queue.push_back(e.clone());
                        }
                    }
                }

                let vcf_records = snpfrag.output_phased_vcf(min_phase_score, min_qual_for_candidate);
                {
                    let mut queue = vcf_records_queue.lock().unwrap();
                    for rd in vcf_records.iter() {
                        queue.push_back(rd.clone());
                    }
                }
            }
        });
    });
    let mut vf = File::create(vcf_file).unwrap();
    vf.write("##fileformat=VCFv4.3\n".as_bytes()).unwrap();
    for ctglen in contig_lengths.iter() {
        let chromosome = ctglen.0.clone();
        let chromosome_len = ctglen.1.clone();
        vf.write(format!("##contig=<ID={},length={}>\n", chromosome, chromosome_len).as_bytes()).unwrap();
    }
    vf.write("##FILTER=<ID=PASS,Description=\"All filters passed\">\n".as_bytes()).unwrap();
    vf.write("##FILTER=<ID=LowQual,Description=\"Low phasing quality\">\n".as_bytes()).unwrap();
    vf.write("##FILTER=<ID=RnaEdit,Description=\"RNA editing\">\n".as_bytes()).unwrap();
    vf.write("##FILTER=<ID=dn,Description=\"Dense cluster of variants\">\n".as_bytes()).unwrap();
    vf.write("##INFO=<ID=RDS,Number=1,Type=String,Description=\"RNA editing or Dense SNP or Single SNP.\">\n".as_bytes()).unwrap();
    vf.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n".as_bytes()).unwrap();
    vf.write("##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase Set\">\n".as_bytes()).unwrap();
    vf.write("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n".as_bytes()).unwrap();
    vf.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n".as_bytes()).unwrap();
    vf.write("##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n".as_bytes()).unwrap();
    vf.write("##FORMAT=<ID=PQ,Number=1,Type=Float,Description=\"Phasing Quality\">\n".as_bytes()).unwrap();
    vf.write("##FORMAT=<ID=AE,Number=A,Type=Integer,Description=\"Haplotype expression of two alleles\">\n".as_bytes()).unwrap();
    vf.write("##FORMAT=<ID=SQ,Number=1,Type=Float,Description=\"Somatic Score\">\n".as_bytes()).unwrap();
    vf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample\n".as_bytes()).unwrap();

    for rd in vcf_records_queue.lock().unwrap().iter() {
        if rd.alternative.len() == 1 {
            vf.write(
                format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                    std::str::from_utf8(&rd.chromosome).unwrap(),
                    rd.position,
                    std::str::from_utf8(&rd.id).unwrap(),
                    std::str::from_utf8(&rd.reference).unwrap(),
                    std::str::from_utf8(&rd.alternative[0]).unwrap(),
                    rd.qual,
                    std::str::from_utf8(&rd.filter).unwrap(),
                    std::str::from_utf8(&rd.info).unwrap(),
                    std::str::from_utf8(&rd.format).unwrap(),
                    rd.genotype
                ).as_bytes(),
            ).unwrap();
        } else if rd.alternative.len() == 2 {
            vf.write(
                format!(
                    "{}\t{}\t{}\t{}\t{},{}\t{}\t{}\t{}\t{}\t{}\n",
                    std::str::from_utf8(&rd.chromosome).unwrap(),
                    rd.position,
                    std::str::from_utf8(&rd.id).unwrap(),
                    std::str::from_utf8(&rd.reference).unwrap(),
                    std::str::from_utf8(&rd.alternative[0]).unwrap(),
                    std::str::from_utf8(&rd.alternative[1]).unwrap(),
                    rd.qual,
                    std::str::from_utf8(&rd.filter).unwrap(),
                    std::str::from_utf8(&rd.info).unwrap(),
                    std::str::from_utf8(&rd.format).unwrap(),
                    rd.genotype
                ).as_bytes(),
            ).unwrap();
        }
    }
    drop(vf);

    if output_read_assignment {
        let mut assignment_writer = File::create(phased_bam_file.replace(".phased.bam", ".assignment.tsv")).unwrap();
        for rd in read_haplotag_queue.lock().unwrap().iter() {
            // println!("{}:{}", rd.0, rd.1);
            assignment_writer.write(format!("{}\t{}\n", rd.0, rd.1).as_bytes()).unwrap();
        }
        drop(assignment_writer);
    }

    if haplotype_specific_exon {
        let mut exon_hashmap: HashMap<String, Vec<(Exon, i32, i32)>> = HashMap::new();
        let mut exon_writer = File::create(phased_bam_file.replace(".phased.bam", ".haplotype_exon.tsv")).unwrap();
        exon_writer.write(
            "#Chromosome\tExon start\tExon end\tExon state\tHap1 expression\tHap2 expression\n".as_bytes(),
        ).unwrap();
        for rd in haplotype_exon_queue.lock().unwrap().iter() {
            let e = &rd.0;
            if exon_hashmap.contains_key(&e.chr) {
                let exon_vec = exon_hashmap.get_mut(&e.chr).unwrap();
                exon_vec.push(rd.clone());
            } else {
                exon_hashmap.insert(e.chr.clone(), vec![rd.clone()]);
            }
        }
        for chr in contig_order.iter() {
            if !exon_hashmap.contains_key(chr) {
                continue;
            }
            let mut exons_sorted = exon_hashmap.get(chr).unwrap().clone();
            exons_sorted.sort_by(|a, b| a.0.start.cmp(&b.0.start));
            for rd in exons_sorted.iter() {
                let e = &rd.0;
                exon_writer.write(
                    format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\n",
                        e.chr,
                        e.start + 1,
                        e.end,
                        e.state,
                        rd.1,
                        rd.2
                    ).as_bytes(),
                ).unwrap(); // 1-based, start inclusive, end inclusive
            }
        }
        drop(exon_writer);
    }

    if !no_bam_output {
        if !haplotype_bam_output {
            let mut read_assignments: HashMap<String, i32> = HashMap::new();
            for rd in read_haplotag_queue.lock().unwrap().iter() {
                if read_assignments.contains_key(&rd.0) {
                    read_assignments.remove(&rd.0);   // one read belongs to at least two regions
                } else {
                    read_assignments.insert(rd.0.clone(), rd.1);
                }
            }
            let mut read_phasesets: HashMap<String, u32> = HashMap::new();
            for rd in read_phaseset_queue.lock().unwrap().iter() {
                if read_phasesets.contains_key(&rd.0) {
                    read_phasesets.remove(&rd.0);   // one read belongs to at least two regions or two phase sets
                } else {
                    read_phasesets.insert(rd.0.clone(), rd.1);
                }
            }
            let mut bam_reader = bam::IndexedReader::from_path(&bam_file).unwrap();
            let header = bam::Header::from_template(&bam_reader.header());
            let mut bam_writer = bam::Writer::from_path(phased_bam_file, &header, Format::Bam).unwrap();
            bam_writer.set_threads(thread_size).unwrap();
            for region in isolated_regions.iter() {
                // TODO: duplicate reads in different regions
                bam_reader.fetch((region.chr.as_str(), region.start, region.end)).unwrap(); // set region
                for r in bam_reader.records() {
                    let mut record = r.unwrap();
                    if record.is_unmapped() || record.is_secondary() || record.is_supplementary() {
                        continue;
                    }
                    if record.reference_start() + 1 < region.start as i64 || record.reference_end() + 1 > region.end as i64 {
                        // reads beyond the region boundary will be ignored to provent duplicated reads
                        continue;
                    }
                    let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
                    if read_assignments.contains_key(&qname) {
                        let asg = read_assignments.get(&qname).unwrap();
                        if *asg != 0 {
                            let _ = record.push_aux(b"HP:i", Aux::I32(*asg));
                        }
                    }
                    if read_phasesets.contains_key(&qname) {
                        let ps = read_phasesets.get(&qname).unwrap();
                        let _ = record.push_aux(b"PS:i", Aux::U32(*ps));
                    }
                    let _ = bam_writer.write(&record).unwrap();
                }
            }
            drop(bam_writer);
        } else {
            let mut hap1_read_assignments: HashSet<String> = HashSet::new();
            let mut hap2_read_assignments: HashSet<String> = HashSet::new();
            for rname in read_haplotag1_queue.lock().unwrap().iter() {
                hap1_read_assignments.insert(rname.clone());
            }
            for rname in read_haplotag2_queue.lock().unwrap().iter() {
                hap2_read_assignments.insert(rname.clone());
            }
            let mut bam_reader = bam::IndexedReader::from_path(&bam_file).unwrap();
            let header = bam::Header::from_template(&bam_reader.header());
            let mut hap1_bam_writer = bam::Writer::from_path(
                phased_bam_file.replace("phased", "hap1"),
                &header,
                Format::Bam,
            ).unwrap();
            let mut hap2_bam_writer = bam::Writer::from_path(
                phased_bam_file.replace("phased", "hap2"),
                &header,
                Format::Bam,
            ).unwrap();
            hap1_bam_writer.set_threads(thread_size).unwrap();
            hap2_bam_writer.set_threads(thread_size).unwrap();
            for region in isolated_regions.iter() {
                bam_reader.fetch((region.chr.as_str(), region.start, region.end)).unwrap(); // set region
                for r in bam_reader.records() {
                    let mut record = r.unwrap();
                    if record.is_unmapped() || record.is_secondary() || record.is_supplementary() {
                        continue;
                    }
                    if record.reference_start() + 1 < region.start as i64 || record.reference_end() + 1 > region.end as i64 {
                        // reads beyond the region boundary will be ignored to provent duplicated reads
                        continue;
                    }
                    let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
                    if hap1_read_assignments.contains(&qname) {
                        let _ = hap1_bam_writer.write(&record).unwrap();
                    } else if hap2_read_assignments.contains(&qname) {
                        let _ = hap2_bam_writer.write(&record).unwrap();
                    }
                }
            }
            drop(hap1_bam_writer);
            drop(hap2_bam_writer);
        }
    }
}
