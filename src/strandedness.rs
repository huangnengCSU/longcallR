use crate::util::Profile;

const MIN_POSITION_DEPTH: u32 = 10;
const MIN_INFORMATIVE_POSITIONS: usize = 50;
const SINGLE_STRAND_MINOR_FRACTION: f32 = 0.05;
const DOUBLE_STRAND_MINOR_FRACTION: f32 = 0.20;
const REQUIRED_POSITION_FRACTION: f32 = 0.80;
const REQUIRED_REGION_FRACTION: f32 = 0.80;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReadStrandedness {
    SingleStrand,
    DoubleStrand,
    Unknown,
}

impl ReadStrandedness {
    pub fn as_str(self) -> &'static str {
        match self {
            ReadStrandedness::SingleStrand => "single-strand",
            ReadStrandedness::DoubleStrand => "double-strand",
            ReadStrandedness::Unknown => "unknown",
        }
    }
}

/// Infer library strandedness from the forward/reverse read counts at covered
/// positions in one genomic region. Ambiguous or sparse profiles are reported
/// as `Unknown` so callers can conservatively avoid strand-bias filtering.
pub fn infer_read_strandedness(profile: &Profile) -> ReadStrandedness {
    let minor_fractions: Vec<f32> = profile
        .freq_vec
        .iter()
        .filter_map(|base_freq| {
            let forward = base_freq.forward_cnt;
            let reverse = base_freq.backward_cnt;
            let depth = forward + reverse;
            (depth >= MIN_POSITION_DEPTH).then_some(forward.min(reverse) as f32 / depth as f32)
        })
        .collect();

    if minor_fractions.len() < MIN_INFORMATIVE_POSITIONS {
        return ReadStrandedness::Unknown;
    }

    let informative_count = minor_fractions.len() as f32;
    let single_fraction = minor_fractions
        .iter()
        .filter(|&&fraction| fraction <= SINGLE_STRAND_MINOR_FRACTION)
        .count() as f32
        / informative_count;
    let double_fraction = minor_fractions
        .iter()
        .filter(|&&fraction| fraction >= DOUBLE_STRAND_MINOR_FRACTION)
        .count() as f32
        / informative_count;

    if double_fraction >= REQUIRED_POSITION_FRACTION {
        ReadStrandedness::DoubleStrand
    } else if single_fraction >= REQUIRED_POSITION_FRACTION {
        ReadStrandedness::SingleStrand
    } else {
        ReadStrandedness::Unknown
    }
}

/// Combine local observations into one library-level classification. Local
/// disagreements are treated as artifacts unless at least 80% of informative
/// regions support the same result.
pub fn infer_dataset_strandedness(
    region_results: &[ReadStrandedness],
) -> (ReadStrandedness, usize) {
    let single_count = region_results
        .iter()
        .filter(|&&result| result == ReadStrandedness::SingleStrand)
        .count();
    let double_count = region_results
        .iter()
        .filter(|&&result| result == ReadStrandedness::DoubleStrand)
        .count();
    let informative_count = single_count + double_count;

    if informative_count == 0 {
        return (ReadStrandedness::Unknown, 0);
    }

    let result = if double_count as f32 / informative_count as f32 >= REQUIRED_REGION_FRACTION {
        ReadStrandedness::DoubleStrand
    } else if single_count as f32 / informative_count as f32 >= REQUIRED_REGION_FRACTION {
        ReadStrandedness::SingleStrand
    } else {
        ReadStrandedness::Unknown
    };
    (result, informative_count)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::util::BaseFreq;

    fn profile_with_counts(counts: &[(u32, u32)]) -> Profile {
        let mut profile = Profile::default();
        profile.freq_vec = counts
            .iter()
            .map(|&(forward, reverse)| BaseFreq {
                forward_cnt: forward,
                backward_cnt: reverse,
                ..BaseFreq::default()
            })
            .collect();
        profile
    }

    #[test]
    fn identifies_single_strand_profile() {
        let profile = profile_with_counts(&vec![(20, 0); 50]);
        assert_eq!(
            infer_read_strandedness(&profile),
            ReadStrandedness::SingleStrand
        );
    }

    #[test]
    fn identifies_double_strand_profile() {
        let profile = profile_with_counts(&vec![(12, 8); 50]);
        assert_eq!(
            infer_read_strandedness(&profile),
            ReadStrandedness::DoubleStrand
        );
    }

    #[test]
    fn sparse_profile_is_unknown() {
        let profile = profile_with_counts(&vec![(12, 8); 49]);
        assert_eq!(infer_read_strandedness(&profile), ReadStrandedness::Unknown);
    }

    #[test]
    fn mixed_profile_is_unknown() {
        let mut counts = vec![(20, 0); 25];
        counts.extend(vec![(12, 8); 25]);
        let profile = profile_with_counts(&counts);
        assert_eq!(infer_read_strandedness(&profile), ReadStrandedness::Unknown);
    }

    #[test]
    fn local_artifact_does_not_change_dataset_classification() {
        let mut results = vec![ReadStrandedness::SingleStrand; 9];
        results.push(ReadStrandedness::DoubleStrand);
        assert_eq!(
            infer_dataset_strandedness(&results),
            (ReadStrandedness::SingleStrand, 10)
        );
    }

    #[test]
    fn conflicting_regions_make_dataset_unknown() {
        let results = vec![
            ReadStrandedness::SingleStrand,
            ReadStrandedness::DoubleStrand,
        ];
        assert_eq!(
            infer_dataset_strandedness(&results),
            (ReadStrandedness::Unknown, 2)
        );
    }
}
