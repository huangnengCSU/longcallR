use mathru::algebra::abstr::cast::ToPrimitive;

use crate::phase::AlleleClass;

pub fn calculate_prob_somatic(hap1_ref_baseqs: &Vec<u8>, hap1_alt_baseqs: &Vec<u8>, hap2_ref_baseqs: &Vec<u8>, hap2_alt_baseqs: &Vec<u8>, purity: f64) -> (AlleleClass, AlleleClass) {
    let mut hap1_allele_class;
    let mut hap2_allele_class;
    let som_rate = 5.0 / 1000000.0;    // each haplotype
    let het_rate = 1.0 / 2000.0;    // each haplotype
    let ref_rate = 1.0 - het_rate - som_rate;

    // for Hap1
    let mut prob_read_ref = 1.0;
    let mut prob_read_het = 1.0;
    let mut prob_read_som = 1.0;

    // P(read|ref), P(read|het), P(read|som)
    for q in hap1_ref_baseqs.iter() {
        let epsilon = 10.0_f64.powf(-(q.to_f64() / 10.0));  // error rate
        prob_read_ref *= 1.0 - epsilon; // ref->ref
        prob_read_het *= epsilon;   // alt->ref
        prob_read_som *= purity * epsilon + (1.0 - purity) * (1.0 - epsilon);   // purity: alt->ref, 1.0-purity: ref->ref
    }
    for q in hap1_alt_baseqs.iter() {
        let epsilon = 10.0_f64.powf(-(q.to_f64() / 10.0));
        prob_read_ref *= epsilon;   // ref->alt
        prob_read_het *= 1.0 - epsilon; // alt->alt
        prob_read_som *= purity * (1.0 - epsilon) + (1.0 - purity) * epsilon; // purity: alt->alt, 1.0-purity: ref->alt
    }

    // prob * prior
    let prob_read_ref_with_prior = prob_read_ref * ref_rate;
    let prob_read_het_with_prior = prob_read_het * het_rate;
    let prob_read_som_with_prior = prob_read_som * som_rate;
    let hap1_prob_ref = prob_read_ref_with_prior / (prob_read_ref_with_prior + prob_read_het_with_prior + prob_read_som_with_prior);
    let hap1_prob_het = prob_read_het_with_prior / (prob_read_ref_with_prior + prob_read_het_with_prior + prob_read_som_with_prior);
    let hap1_prob_som = prob_read_som_with_prior / (prob_read_ref_with_prior + prob_read_het_with_prior + prob_read_som_with_prior);
    if hap1_prob_som > hap1_prob_ref && hap1_prob_som > hap1_prob_het {
        hap1_allele_class = AlleleClass { allcls: 2, prob: hap1_prob_som };
    } else if hap1_prob_het > hap1_prob_ref && hap1_prob_het > hap1_prob_som {
        hap1_allele_class = AlleleClass { allcls: 1, prob: hap1_prob_het };
    } else {
        hap1_allele_class = AlleleClass { allcls: 0, prob: hap1_prob_ref };
    }

    // for Hap2
    let mut prob_read_ref = 1.0;
    let mut prob_read_het = 1.0;
    let mut prob_read_som = 1.0;

    // P(read|ref), P(read|het), P(read|som)
    for q in hap2_ref_baseqs.iter() {
        let epsilon = 10.0_f64.powf(-(q.to_f64() / 10.0));  // error rate
        prob_read_ref *= 1.0 - epsilon; // ref->ref
        prob_read_het *= epsilon;   // alt->ref
        prob_read_som *= purity * epsilon + (1.0 - purity) * (1.0 - epsilon);   // purity: alt->ref, 1.0-purity: ref->ref
    }
    for q in hap2_alt_baseqs.iter() {
        let epsilon = 10.0_f64.powf(-(q.to_f64() / 10.0));
        prob_read_ref *= epsilon;   // ref->alt
        prob_read_het *= 1.0 - epsilon; // alt->alt
        prob_read_som *= purity * (1.0 - epsilon) + (1.0 - purity) * epsilon; // purity: alt->alt, 1.0-purity: ref->alt
    }

    // prob * prior
    let prob_read_ref_with_prior = prob_read_ref * ref_rate;
    let prob_read_het_with_prior = prob_read_het * het_rate;
    let prob_read_som_with_prior = prob_read_som * som_rate;
    let hap2_prob_ref = prob_read_ref_with_prior / (prob_read_ref_with_prior + prob_read_het_with_prior + prob_read_som_with_prior);
    let hap2_prob_het = prob_read_het_with_prior / (prob_read_ref_with_prior + prob_read_het_with_prior + prob_read_som_with_prior);
    let hap2_prob_som = prob_read_som_with_prior / (prob_read_ref_with_prior + prob_read_het_with_prior + prob_read_som_with_prior);
    if hap2_prob_som > hap2_prob_ref && hap2_prob_som > hap2_prob_het {
        hap2_allele_class = AlleleClass { allcls: 2, prob: hap2_prob_som };
    } else if hap2_prob_het > hap2_prob_ref && hap2_prob_het > hap2_prob_som {
        hap2_allele_class = AlleleClass { allcls: 1, prob: hap2_prob_het };
    } else {
        hap2_allele_class = AlleleClass { allcls: 0, prob: hap2_prob_ref };
    }
    return (hap1_allele_class, hap2_allele_class);
}