use std::collections::{HashMap, VecDeque};
use std::fs;
use std::fs::File;
use std::io::Write;
use std::sync::Mutex;

use rayon::iter::IntoParallelRefIterator;
use rayon::prelude::*;
use rust_htslib::{bam, bam::ext::BamRecordExtensions, bam::record::Aux, bam::Format, bam::Read};
use rust_lapper::Interval;

use crate::snpfrags::SNPFrag;
use crate::strandedness::{
    infer_dataset_strandedness, infer_read_strandedness, ReadStrandedness,
};
use crate::util::{load_reference, parse_fai, Profile, Region};
use crate::vcf::get_genotype_quality_phase_from_vcf;
use crate::Platform;

pub fn run(
    bam_file: &str,
    ref_file: &str,
    input_vcf_file: Option<String>,
    output_vcf_file: &str,
    phased_bam_file: &str,
    thread_size: usize,
    isolated_regions: Vec<Region>,
    exon_regions: HashMap<String, Vec<Interval<u32, u8>>>,
    exon_only: bool,
    platform: &Platform,
    min_mapq: u8,
    min_baseq: u8,
    divergence: f32,
    min_allele_freq: f32,
    min_qual: u32,
    min_allele_freq_include_intron: f32,
    use_strand_bias: Option<bool>,
    min_depth: u32,
    max_depth: u32,
    downsample: bool,
    downsample_depth: u32,
    min_read_length: usize,
    distance_to_read_end: u32,
    polya_tail_len: u32,
    dense_win_size: u32,
    min_dense_cnt: u32,
    min_linkers: u32,
    min_phase_score: f32,
    max_enum_snps: usize,
    read_assignment_cutoff: f64,
    direct_haplotag_from_vcf: bool,
    no_bam_output: bool,
    low_allele_frac_cutoff: f32,
    low_allele_cnt_cutoff: u32,
) {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(thread_size)
        .build()
        .unwrap();
    let vcf_records_queue = Mutex::new(VecDeque::new());
    let read_haplotag_queue = Mutex::new(VecDeque::new());
    let read_phaseset_queue = Mutex::new(VecDeque::new());
    let ref_seqs = load_reference(ref_file)
        .unwrap_or_else(|e| panic!("{}", e));
    let fai_path = String::from(ref_file) + ".fai";
    if fs::metadata(&fai_path).is_err() {
        panic!("Reference index file .fai does not exist.");
    }
    let contig_lengths = parse_fai(fai_path.as_str());
    let chr_order: HashMap<String, usize> = contig_lengths.iter()
        .enumerate()
        .map(|(i, (chr, _))| (chr.clone(), i))
        .collect();
    let candidates_genotype_quality = if let Some(vcf_file) = input_vcf_file.clone() {
        get_genotype_quality_phase_from_vcf(vcf_file.as_str())
    } else {
        HashMap::new()
    };

    let (use_strand_bias_for_dataset, strand_bias_source, library_strand) = match use_strand_bias {
        Some(value) => (
            value,
            "explicit --strand-bias setting".to_string(),
            "not-detected",
        ),
        None => {
            const MAX_STRAND_DETECTION_REGIONS: usize = 20;
            let region_count = isolated_regions.len();
            let sample_count = region_count.min(MAX_STRAND_DETECTION_REGIONS);
            let sampled_regions: Vec<&Region> = if sample_count == 0 {
                Vec::new()
            } else {
                (0..sample_count)
                    .map(|sample_index| {
                        &isolated_regions[sample_index * region_count / sample_count]
                    })
                    .collect()
            };
            let region_results: Vec<ReadStrandedness> = pool.install(|| {
                sampled_regions
                    .par_iter()
                    .map(|reg| {
                        let mut profile = Profile::default();
                        let ref_seq = ref_seqs.get(&reg.chr).unwrap_or_else(|| {
                            panic!(
                                "Chromosome '{}' from target region is missing in reference FASTA",
                                reg.chr
                            )
                        });
                        profile.fill_data_into_freq_vec(
                            bam_file,
                            reg,
                            ref_seq,
                            platform,
                            min_mapq,
                            min_read_length,
                            divergence,
                            distance_to_read_end,
                            polya_tail_len,
                        );
                        infer_read_strandedness(&profile)
                    })
                    .collect()
            });
            let (strandedness, informative_regions) =
                infer_dataset_strandedness(&region_results);
            (
                matches!(strandedness, ReadStrandedness::DoubleStrand),
                format!(
                    "auto-detected {} library from {} informative sampled regions",
                    strandedness.as_str(),
                    informative_regions
                ),
                strandedness.short_name(),
            )
        }
    };
    eprintln!(
        "Strand bias filtering: {} for all regions ({})",
        if use_strand_bias_for_dataset {
            "enabled"
        } else {
            "disabled"
        },
        strand_bias_source
    );

    pool.install(|| {
        isolated_regions.par_iter().for_each(|reg| {
            let mut profile = Profile::default();
            let ref_seq = ref_seqs.get(&reg.chr).unwrap_or_else(|| {
                panic!(
                    "Chromosome '{}' from target region is missing in reference FASTA",
                    reg.chr
                )
            });
            let mut exon_region_vec = Vec::new();
            if exon_only {
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
            profile.fill_data_into_freq_vec(
                bam_file,
                &reg,
                ref_seq,
                platform,
                min_mapq,
                min_read_length,
                divergence,
                distance_to_read_end,
                polya_tail_len,
            );
            let mut snpfrag = SNPFrag::default();
            snpfrag.region = reg.clone();
            snpfrag.min_linkers = min_linkers;
            if input_vcf_file.is_some() {
                let chr_candidates_genotype_quality = candidates_genotype_quality.get(&reg.chr);    // 0-based
                if let Some(chr_candidates_genotype_quality) = chr_candidates_genotype_quality {
                    snpfrag.import_external_candidates(
                        &profile,
                        ref_seq,
                        chr_candidates_genotype_quality,
                        0.0,
                        direct_haplotag_from_vcf,
                        exon_region_vec.clone(),
                        exon_only,
                        min_depth,
                        min_allele_freq,
                    );
                }
            } else {
                snpfrag.get_candidate_snps(
                    &profile,
                    exon_region_vec,
                    exon_only,
                    min_allele_freq,
                    min_qual,
                    min_allele_freq_include_intron,
                    min_depth,
                    max_depth,
                    min_baseq,
                    use_strand_bias_for_dataset,
                    dense_win_size,
                    min_dense_cnt,
                    low_allele_frac_cutoff,
                    low_allele_cnt_cutoff,
                );
            }
            // TODO: for very high depth region, down-sampling the reads
            snpfrag.get_fragments(
                &bam_file,
                &reg,
                ref_seq,
                min_mapq,
                min_read_length,
                divergence,
            );
            let apply_downsampling = downsample 
                && downsample_depth > 0 
                && snpfrag.fragments.len() >= downsample_depth as usize;
            if apply_downsampling {
                unsafe {
                    snpfrag.downsample_fragments(downsample_depth, 2025);
                }
            }
            if snpfrag.fragments.len() > 0 {
                println!(
                    "number of fragments: {:?} in {:?}, apply downsampling: {:?}, strand: {}",
                    snpfrag.fragments.len(),
                    reg,
                    apply_downsampling,
                    library_strand
                );
            }
            // snpfrag.clean_fragments();

            if direct_haplotag_from_vcf {
                snpfrag.seed_haplotag_from_phased_snps();
                let read_assignments = snpfrag.assign_reads_haplotype(read_assignment_cutoff, apply_downsampling);
                let phase_sets = snpfrag.assign_phase_set(0.0);

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

                return;
            }

            snpfrag.init_haplotypes();
            snpfrag.init_assignment();
            // snpfrag.chain_phase(max_enum_snps);
            snpfrag.phase(1, max_enum_snps, apply_downsampling);
            // let read_assignments = snpfrag.assign_reads_haplotype(read_assignment_cutoff);
            snpfrag.assign_reads_haplotype(read_assignment_cutoff, apply_downsampling); // run1
            snpfrag.assign_snp_haplotype_genotype(apply_downsampling);  // run1
            // let read_assignments = snpfrag.assign_reads_haplotype(read_assignment_cutoff);
            snpfrag.assign_reads_haplotype(read_assignment_cutoff, apply_downsampling); // run2
            snpfrag.assign_snp_haplotype_genotype(apply_downsampling);  // run2
            // snpfrag.assign_het_var_haplotype(min_phase_score, low_allele_frac_cutoff, low_allele_cnt_cutoff);
            // snpfrag.eval_low_frac_het_var_phase(min_phase_score, low_allele_frac_cutoff, low_allele_cnt_cutoff);

            // Evaluate RNA editing sites and low allele fraction variants using a relaxed phase score threshold.
            // This step reinserts these candidates into the phasing process to ensure that phase sites are properly updated.
            // Without this, the subsequent `assign_reads_haplotype` call would be ineffective, as it depends on an up-to-date set of phase sites.
            snpfrag.eval_rna_edit_var_phase(min_phase_score - 3.0, apply_downsampling);
            snpfrag.eval_low_frac_var_phase(min_phase_score - 3.0, apply_downsampling);
            let read_assignments = snpfrag.assign_reads_haplotype(read_assignment_cutoff, false);   // run3
            snpfrag.assign_snp_haplotype_genotype(false);   // run3

            // snpfrag.eval_hom_var_phase(min_phase_score);
            // assign phased fragments to low fraction variants and detect confident low fraction variants
            // println!("low_frac: {}", snpfrag.low_frac_snps.len());
            // snpfrag.detect_low_frac_by_het(&bam_file.as_str(), &reg);
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
            let phase_sets = snpfrag.assign_phase_set(min_phase_score);

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

            let vcf_records = snpfrag.output_phased_vcf(min_phase_score);
            {
                let mut queue = vcf_records_queue.lock().unwrap();
                for rd in vcf_records.iter() {
                    queue.push_back(rd.clone());
                }
            }
        });
    });
    if !direct_haplotag_from_vcf {
        let mut vf = File::create(output_vcf_file).unwrap();
        vf.write("##fileformat=VCFv4.3\n".as_bytes()).unwrap();
        for ctglen in contig_lengths.iter() {
            let chromosome = ctglen.0.clone();
            let chromosome_len = ctglen.1.clone();
            vf.write(format!("##contig=<ID={},length={}>\n", chromosome, chromosome_len).as_bytes())
                .unwrap();
        }
        vf.write("##FILTER=<ID=PASS,Description=\"All filters passed\">\n".as_bytes())
            .unwrap();
        vf.write("##FILTER=<ID=LowQual,Description=\"Low phasing quality\">\n".as_bytes())
            .unwrap();
        vf.write("##FILTER=<ID=HomRef,Description=\"Homo reference\">\n".as_bytes())
            .unwrap();
        vf.write("##FILTER=<ID=RnaEdit,Description=\"RNA editing\">\n".as_bytes())
            .unwrap();
        vf.write("##FILTER=<ID=Multiallelic,Description=\"Multiallelic SNP\">\n".as_bytes())
            .unwrap();
        vf.write("##FILTER=<ID=Dense,Description=\"Dense cluster of variants\">\n".as_bytes())
            .unwrap();
        vf.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n".as_bytes())
            .unwrap();
        vf.write("##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase Set\">\n".as_bytes())
            .unwrap();
        vf.write(
            "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n".as_bytes(),
        )
        .unwrap();
        vf.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n".as_bytes())
            .unwrap();
        vf.write("##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n".as_bytes())
            .unwrap();
        vf.write("##FORMAT=<ID=PQ,Number=1,Type=Integer,Description=\"Phasing Quality\">\n".as_bytes())
            .unwrap();
        vf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample\n".as_bytes())
            .unwrap();

        let mut vcf_records: Vec<_> = vcf_records_queue.lock().unwrap().iter().cloned().collect();
        vcf_records.sort_by_key(|rd| (*chr_order.get(&rd.chromosome).unwrap_or(&usize::MAX), rd.position));
        for rd in vcf_records.iter() {
            if rd.alternative.len() == 1 {
                vf.write(
                    format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                        &rd.chromosome,
                        rd.position,
                        std::str::from_utf8(&rd.id).unwrap(),
                        std::str::from_utf8(&rd.reference).unwrap(),
                        std::str::from_utf8(&rd.alternative[0]).unwrap(),
                        rd.qual,
                        std::str::from_utf8(&rd.filter).unwrap(),
                        ".",
                        std::str::from_utf8(&rd.format).unwrap(),
                        rd.genotype
                    )
                    .as_bytes(),
                )
                .unwrap();
            } else if rd.alternative.len() == 2 {
                vf.write(
                    format!(
                        "{}\t{}\t{}\t{}\t{},{}\t{}\t{}\t{}\t{}\t{}\n",
                        &rd.chromosome,
                        rd.position,
                        std::str::from_utf8(&rd.id).unwrap(),
                        std::str::from_utf8(&rd.reference).unwrap(),
                        std::str::from_utf8(&rd.alternative[0]).unwrap(),
                        std::str::from_utf8(&rd.alternative[1]).unwrap(),
                        rd.qual,
                        std::str::from_utf8(&rd.filter).unwrap(),
                        ".",
                        std::str::from_utf8(&rd.format).unwrap(),
                        rd.genotype
                    )
                    .as_bytes(),
                )
                .unwrap();
            }
        }
        drop(vf);
    }

    if !no_bam_output {
        let mut read_assignments: HashMap<String, i32> = HashMap::new();
        for rd in read_haplotag_queue.lock().unwrap().iter() {
            if read_assignments.contains_key(&rd.0) {
                // read_assignments.remove(&rd.0);   // one read belongs to at least two regions
                continue;
            } else {
                read_assignments.insert(rd.0.clone(), rd.1);
            }
        }
        let mut read_phasesets: HashMap<String, u32> = HashMap::new();
        for rd in read_phaseset_queue.lock().unwrap().iter() {
            if read_phasesets.contains_key(&rd.0) {
                // read_phasesets.remove(&rd.0);   // one read belongs to at least two regions or two phase sets
                continue;
            } else {
                read_phasesets.insert(rd.0.clone(), rd.1);
            }
        }
        let mut sorted_regions = isolated_regions.clone();
        sorted_regions.sort_by(|a, b| {
            chr_order.get(&a.chr).cmp(&chr_order.get(&b.chr))
                .then(a.start.cmp(&b.start))
        });
        let mut bam_reader = bam::IndexedReader::from_path(&bam_file).unwrap();
        let header = bam::Header::from_template(&bam_reader.header());
        let mut bam_writer = bam::Writer::from_path(phased_bam_file, &header, Format::Bam).unwrap();
        bam_writer.set_threads(thread_size).unwrap();
        for region in sorted_regions.iter() {
            // TODO: duplicate reads in different regions
            bam_reader
                .fetch((region.chr.as_str(), region.start, region.end))
                .unwrap(); // set region
            for r in bam_reader.records() {
                let mut record = r.unwrap();
                if record.is_unmapped() || record.is_secondary() || record.is_supplementary() {
                    continue;
                }
                if record.reference_start() + 1 < region.start as i64
                    || record.reference_end() + 1 > region.end as i64
                {
                    // reads beyond the region boundary will be ignored to provent duplicated reads
                    continue;
                }
                let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
                // if read_assignments.contains_key(&qname) {
                //     let asg = read_assignments.get(&qname).unwrap();
                //     if *asg != 0 {
                //         let _ = record.push_aux(b"HP:i", Aux::I32(*asg));
                //     }
                // }
                // if read_phasesets.contains_key(&qname) {
                //     let ps = read_phasesets.get(&qname).unwrap();
                //     let _ = record.push_aux(b"PS:i", Aux::U32(*ps));
                // }
                if let (Some(asg), Some(ps)) = (read_assignments.get(&qname), read_phasesets.get(&qname)) {
                    if *asg != 0 {
                        let _ = record.push_aux(b"HP:i", Aux::I32(*asg));
                        let _ = record.push_aux(b"PS:i", Aux::U32(*ps));
                    }
                }
                let _ = bam_writer.write(&record).unwrap();
            }
        }
        drop(bam_writer);
    }
}
