use crate::matrix::ColumnBaseCount;

#[derive(Clone)]
struct SpliceMatrixElement {
    m: f64,
    // not matching
    ix: f64,
    // gap in query: deletion
    ix2: f64,
    // gap in query: introns
    m_s: u8,    // previous state of m
    ix_s: u8,   // previous state of ix
    ix2_s: u8,  // previous state of ix2

}

impl Default for SpliceMatrixElement {
    fn default() -> SpliceMatrixElement {
        SpliceMatrixElement {
            m: -f64::INFINITY,
            ix: -f64::INFINITY,
            ix2: -f64::INFINITY,
            m_s: 0,
            ix_s: 0,
            ix2_s: 0,
        }
    }
}


pub fn banded_nw_splice_aware3(
    query: &Vec<u8>,
    profile: &Vec<ColumnBaseCount>,
    reduced_donor_penalty: &Vec<f64>,
    reduced_acceptor_penalty: &Vec<f64>,
    splice_boundary: &Vec<bool>,
    width: usize,
) -> (f64, Vec<u8>, Vec<u8>, Vec<u8>) {
    // TODO: calculate penalty of passing a intron region, where the region is cut by the process of cutting all reads introns.
    // TODO: useless of hidden_splice_penalty, remove later.
    // Case: wtc11_ont_grch38: chr22:37024802-37025006
    let h = 2.0; // short gap open, q in minimap2
    let g = 1.0; // short gap extend, e in minimap2
    let h2 = 32.0; // long gap open, q hat in minimap2
    // let p = 9.0; // splice junction site penalty
    let b: f64 = 34.0; // splice junction site penalty
    let match_score = 2.0;
    let mismatch_score = -2.0;

    // let query_without_gap = query.replace("-", "").replace("N", "");
    let query_t = std::str::from_utf8(query.clone().as_slice()).unwrap().replace("-", "").replace("N", "");
    let query_without_gap = query_t.as_bytes();

    let q_len = query_without_gap.len();
    let t_len = profile.len();

    let mut mat: Vec<Vec<SpliceMatrixElement>> = vec![vec![SpliceMatrixElement { ..Default::default() }; (2 * width) + 3]; t_len + 1];
    // println!("mat size: {} x {}", mat.len(), mat[0].len());

    let mut k_vec: Vec<usize> = vec![0; t_len + 1];

    // Initialize first row
    mat[0][width + 1].ix = -h - g;
    mat[0][width + 1].ix2 = -f64::INFINITY;
    mat[0][width + 1].m = mat[0][width + 1].ix;

    let mut i = 1; // index on target
    let mut j = 1; // index on query
    let mut k = 0; // index on query_without_gap

    let mut left_bound: usize; // include this index
    let mut right_bound: usize; // exclude this index
    while i < t_len + 1 && k < q_len + 1 {
        let mut offset: usize;
        if query[j - 1] != b'-' && query[j - 1] != b'N' {
            k += 1;
            offset = 0;
        } else {
            offset = 1;
        }
        k_vec[i] = k; // store k index of each row

        if k as i32 - width as i32 > 0 {
            left_bound = (k as i32 - width as i32) as usize;
        } else {
            left_bound = 1;
        }
        right_bound = if k + width + 1 <= q_len + 1 {
            k + width + 1
        } else {
            q_len + 1
        };
        for u in left_bound..right_bound {
            // index on alignment matrix: (i,v), v = (w+1)+(u-k), u in [left_bound, right_bound)
            let v = (width as i32 + 1 + (u as i32 - k as i32)) as usize;
            let col = &profile[i - 1];
            let qbase = query_without_gap[u - 1];
            let sij: f64;
            if qbase == b' ' {
                sij = 0.0;
            } else {
                if col.get_depth() <= 5 {
                    let tbase = col.get_ref_base();
                    if tbase == b'-' {
                        sij = match_score;
                    } else if qbase == tbase {
                        sij = match_score;
                    } else {
                        sij = mismatch_score;
                    }
                } else {
                    sij = 2.0 - 4.0 * col.get_score(&qbase);
                }
            }

            if col.get_major_base() == b'-' || col.get_major_base() == b'N' {
                if mat[i - 1][v + 1 - offset].m >= mat[i - 1][v + 1 - offset].ix {
                    mat[i][v].ix = mat[i - 1][v + 1 - offset].m;
                    mat[i][v].ix_s = 1; //mat[i][v].ix_prev_m = true;
                } else {
                    mat[i][v].ix = mat[i - 1][v + 1 - offset].ix;
                    mat[i][v].ix_s = 2; //mat[i][v].ix_prev_ix = true;
                }
            } else {
                if mat[i - 1][v + 1 - offset].m - h - g >= mat[i - 1][v + 1 - offset].ix - g {
                    mat[i][v].ix = mat[i - 1][v + 1 - offset].m - h - g;
                    mat[i][v].ix_s = 1; //mat[i][v].ix_prev_m = true;
                } else {
                    mat[i][v].ix = mat[i - 1][v + 1 - offset].ix - g;
                    mat[i][v].ix_s = 2; //mat[i][v].ix_prev_ix = true;
                }
            }

            if mat[i - 1][v + 1 - offset].m - reduced_donor_penalty[i - 1] - h2 >= mat[i - 1][v + 1 - offset].ix2 {
                mat[i][v].ix2 = mat[i - 1][v + 1 - offset].m - reduced_donor_penalty[i - 1] - h2;
                mat[i][v].ix2_s = 1; //mat[i][v].ix2_prev_m = true;
            } else {
                mat[i][v].ix2 = mat[i - 1][v + 1 - offset].ix2;
                mat[i][v].ix2_s = 3; //mat[i][v].ix2_prev_ix2 = true;
            }


            mat[i][v].m = (mat[i - 1][v - offset].m + sij)
                .max(mat[i][v].ix)
                .max(mat[i][v].ix2 - reduced_acceptor_penalty[i]); // index i in matrix is corresponding to the index i-1 in reference (donor penalty and acceptor penalty)
            if mat[i][v].m == mat[i - 1][v - offset].m + sij {
                mat[i][v].m_s = 1; //mat[i][v].m_prev_m = true;
            } else if mat[i][v].m == mat[i][v].ix {
                mat[i][v].m_s = 2; //mat[i][v].m_prev_ix = true;
            } else if mat[i][v].m == mat[i][v].ix2 - reduced_acceptor_penalty[i] {
                mat[i][v].m_s = 3; //mat[i][v].m_prev_ix2 = true;
            }
        }
        i += 1;
        j += 1;
    }

    // trace back
    let mut aligned_query: Vec<u8> = Vec::new();
    let mut ref_target: Vec<u8> = Vec::new();
    let mut major_target: Vec<u8> = Vec::new();
    let mut alignment_score = 0.0;

    let mut u: usize;
    let mut v: usize;

    i = t_len;

    // let mut score_vec = Vec::new();
    // for vv in 0..2 * width + 3 {
    //     score_vec.push(mat[i][vv].m);
    // }
    // find max value and index in score_vec
    // let (max_score, max_index) = find_max_value_and_index(&score_vec);
    let vv = (width as i32 + 1 + (q_len as i32 - k_vec[i] as i32)) as usize;
    let max_score = mat[i][vv].m;

    alignment_score = max_score;
    v = vv;
    k = k_vec[i];
    u = (v as i32 - (width as i32 + 1) + k as i32) as usize; // index on query_without_gap

    let mut trace_back_stat: u8 = 0;
    if mat[i][v].m_s == 2 {
        //mat[i][v].m_prev_ix
        trace_back_stat = 2;
        // trace_back_stat = TraceBack::IX;
    } else if mat[i][v].m_s == 3 {
        //mat[i][v].m_prev_ix2
        trace_back_stat = 3;
        // trace_back_stat = TraceBack::IX2;
    } else if mat[i][v].m_s == 1 {
        // mat[i][v].m_prev_m
        trace_back_stat = 1;
        // trace_back_stat = TraceBack::M;
    } else {
        panic!("Error: no traceback");
    }

    while i > 0 && u > 0 {
        // println!("i: {}, k: {}, m:{}, ix: {}, iy:{}, ix2:{}", i, k, mat[i][k].m, mat[i][k].ix, mat[i][k].iy, mat[i][k].ix2);
        k = k_vec[i];
        v = (width as i32 + 1 + (u as i32 - k as i32)) as usize;
        let qbase = query_without_gap[u - 1];
        let ref_base = profile[i - 1].get_ref_base();
        let major_base = profile[i - 1].get_major_base();
        if trace_back_stat == 1 {
            // trace_back_stat == TraceBack::M
            if mat[i][v].m_s == 2 {
                //mat[i][v].m_prev_ix
                trace_back_stat = 2;
                // trace_back_stat = TraceBack::IX;
            } else if mat[i][v].m_s == 3 {
                // mat[i][v].m_prev_ix2
                trace_back_stat = 3;
                // trace_back_stat = TraceBack::IX2;
            } else if mat[i][v].m_s == 1 {
                // mat[i][v].m_prev_m
                aligned_query.push(qbase);
                ref_target.push(ref_base);
                major_target.push(major_base);
                i -= 1;
                u -= 1;
                trace_back_stat = 1;
                // trace_back_stat = TraceBack::M;
            } else {
                panic!("Error: no traceback");
            }
        } else if trace_back_stat == 2 {
            // trace_back_stat == TraceBack::IX
            if mat[i][v].ix_s == 2 {
                // mat[i][v].ix_prev_ix
                aligned_query.push(b'-');
                ref_target.push(ref_base);
                major_target.push(major_base);
                i -= 1;
                trace_back_stat = 2;
                // trace_back_stat = TraceBack::IX;
            } else if mat[i][v].ix_s == 1 {
                //mat[i][v].ix_prev_m
                aligned_query.push(b'-');
                ref_target.push(ref_base);
                major_target.push(major_base);
                i -= 1;
                trace_back_stat = 1;
                // trace_back_stat = TraceBack::M;
            } else {
                panic!("Error: no traceback");
            }
        } else if trace_back_stat == 3 {
            // trace_back_stat == TraceBack::IX2
            if mat[i][v].ix2_s == 3 {
                // mat[i][v].ix2_prev_ix2
                aligned_query.push(b'N');
                ref_target.push(ref_base);
                major_target.push(major_base);
                i -= 1;
                trace_back_stat = 3;
                // trace_back_stat = TraceBack::IX2;
            } else if mat[i][v].ix2_s == 1 {
                // mat[i][v].ix2_prev_m
                aligned_query.push(b'N');
                ref_target.push(ref_base);
                major_target.push(major_base);
                i -= 1;
                trace_back_stat = 1;
                // trace_back_stat = TraceBack::M;
            } else {
                panic!("Error: no traceback");
            }
        } else {
            panic!("Error: no traceback");
        }
    }

    while i > 0 {
        let ref_base = profile[i - 1].get_ref_base();
        let major_base = profile[i - 1].get_major_base();
        aligned_query.push(b'-');
        ref_target.push(ref_base);
        major_target.push(major_base);
        i -= 1;
    }
    while u > 0 {
        let qbase = query_without_gap[u - 1];
        aligned_query.push(qbase);
        ref_target.push(b'-');
        major_target.push(b'-');
        u -= 1;
    }

    aligned_query.reverse();
    ref_target.reverse();
    major_target.reverse();
    // println!("original:\n {:?}", std::str::from_utf8(query).unwrap());
    // println!("ref:\n {:?}", std::str::from_utf8(ref_target.as_slice()).unwrap());
    // println!("aligned:\n {:?}", std::str::from_utf8(aligned_query.as_slice()).unwrap());
    (alignment_score, aligned_query, ref_target, major_target)
}