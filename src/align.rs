use std::process;
use crate::matrix::ColumnBaseCount;

struct SpliceMatrixElement {
    m: f64,
    // not matching
    ix: f64,
    // gap in query: deletion
    iy: f64,
    // gap in target: insertion
    ix2: f64,
    // gap in query: introns
    m_prev_m: bool,
    m_prev_ix: bool,
    m_prev_iy: bool,
    m_prev_ix2: bool,
    ix_prev_m: bool,
    ix_prev_ix: bool,
    iy_prev_m: bool,
    iy_prev_iy: bool,
    ix2_prev_m: bool,
    ix2_prev_ix2: bool,
}

#[derive(Eq, PartialEq)]
enum TraceBack {
    M,
    IX,
    IY,
    IX2,
}


pub fn nw_splice_aware(query: &String, profile: &Vec<ColumnBaseCount>) -> (f64, String, String, String) {
    let h = 2.0;
    let g = 1.0;
    let h2 = 12.0;
    let p = 4.0;

    let q_len = query.len();
    let t_len = profile.len();


    let mut mat: Vec<Vec<SpliceMatrixElement>> = Vec::new();
    for _ in 0..t_len + 1 {
        let mut row: Vec<SpliceMatrixElement> = Vec::new();
        for _ in 0..q_len + 1 {
            row.push(SpliceMatrixElement {
                m: 0.0,
                ix: 0.0,
                iy: 0.0,
                ix2: 0.0,
                m_prev_m: false,
                m_prev_ix: false,
                m_prev_iy: false,
                m_prev_ix2: false,
                ix_prev_m: false,
                ix_prev_ix: false,
                iy_prev_m: false,
                iy_prev_iy: false,
                ix2_prev_m: false,
                ix2_prev_ix2: false,
            });
        }
        mat.push(row);
    }

    // Initialize first row and column
    mat[0][0].ix = -h - g;
    mat[0][0].iy = -h - g - f64::INFINITY;  // no gap in target
    mat[0][0].ix2 = -h2;
    mat[0][0].m = mat[0][0].ix.max(mat[0][0].iy).max(mat[0][0].ix2 - p);

    // initialize first row and column are also follow the formula
    for j in 1..q_len + 1 {
        mat[0][j].ix = -f64::INFINITY;
        mat[0][j].iy = (mat[0][j - 1].m - h - g - f64::INFINITY).max(mat[0][j - 1].iy - g - f64::INFINITY);
        mat[0][j].ix2 = -f64::INFINITY;
        mat[0][j].m = mat[0][j].ix.max(mat[0][j].iy).max(mat[0][j].ix2 - p);
    }

    for i in 1..t_len + 1 {
        mat[i][0].ix = (mat[i - 1][0].m - h - g).max(mat[i - 1][0].ix - g);
        mat[i][0].iy = -f64::INFINITY;
        mat[i][0].ix2 = (mat[i - 1][0].m - p - h2).max(mat[i - 1][0].ix2);
        mat[i][0].m = mat[i][0].ix.max(mat[i][0].iy).max(mat[i][0].ix2 - p);
    }

    // Fill in matrices
    for i in 1..t_len + 1 {
        for j in 1..q_len + 1 {
            let qbase = query.as_bytes()[j - 1];
            // let tbase = target.chars().nth(i - 1).unwrap();
            let col = &profile[i - 1];
            let sij = 2.0 - 4.0 * col.get_score(&qbase);

            // if target is dash, the cost of gap open and gap extension is 0
            if col.get_major_base() == b'-' {
                mat[i][j].ix = mat[i - 1][j].ix.max(mat[i - 1][j].m);
                if mat[i][j].ix == mat[i - 1][j].m {
                    mat[i][j].ix_prev_m = true;
                } else if mat[i][j].ix == mat[i - 1][j].ix {
                    mat[i][j].ix_prev_ix = true;
                }
            } else {
                mat[i][j].ix = (mat[i - 1][j].ix - g).max(mat[i - 1][j].m - h - g);
                if mat[i][j].ix == mat[i - 1][j].m - h - g {
                    mat[i][j].ix_prev_m = true;
                } else if mat[i][j].ix == mat[i - 1][j].ix - g {
                    mat[i][j].ix_prev_ix = true;
                }
            }

            mat[i][j].iy = (mat[i][j - 1].iy - g - f64::INFINITY).max(mat[i][j - 1].m - h - g - f64::INFINITY);
            if (mat[i][j].iy == mat[i][j - 1].m - h - g - f64::INFINITY) {
                mat[i][j].iy_prev_m = true;
            } else if (mat[i][j].iy == mat[i][j - 1].iy - g - f64::INFINITY) {
                mat[i][j].iy_prev_iy = true;
            }

            if col.get_major_base() == b'-' {
                mat[i][j].ix2 = mat[i - 1][j].ix2.max(mat[i - 1][j].m);
                if (mat[i][j].ix2 == mat[i - 1][j].m) {
                    mat[i][j].ix2_prev_m = true;
                } else if (mat[i][j].ix2 == mat[i - 1][j].ix2) {
                    mat[i][j].ix2_prev_ix2 = true;
                }
            } else {
                mat[i][j].ix2 = mat[i - 1][j].ix2.max(mat[i - 1][j].m - p - h2);
                if (mat[i][j].ix2 == mat[i - 1][j].m - p - h2) {
                    mat[i][j].ix2_prev_m = true;
                } else if (mat[i][j].ix2 == mat[i - 1][j].ix2) {
                    mat[i][j].ix2_prev_ix2 = true;
                }
            }

            mat[i][j].m = (mat[i - 1][j - 1].m + sij).max(mat[i][j].ix.max(mat[i][j].iy.max(mat[i][j].ix2 - p)));
            if (mat[i][j].m == mat[i][j].ix) {
                mat[i][j].m_prev_ix = true;
            } else if (mat[i][j].m == mat[i][j].iy) {
                mat[i][j].m_prev_iy = true;
            } else if (mat[i][j].m == mat[i][j].ix2 - p) {
                mat[i][j].m_prev_ix2 = true;
            } else if (mat[i][j].m == mat[i - 1][j - 1].m + sij) {
                mat[i][j].m_prev_m = true;
            }
        }
    }

    // // print matrix
    // // println!("m matrix:");
    // for i in 0..t_len + 1 {
    //     print!("M  :\t");
    //     for j in 0..q_len + 1 {
    //         print!("{}  ", mat[i][j].m);
    //     }
    //     println!();
    //     print!("Ix :\t");
    //     for j in 0..q_len + 1 {
    //         print!("{}  ", mat[i][j].ix);
    //     }
    //     println!();
    //     print!("Iy :\t");
    //     for j in 0..q_len + 1 {
    //         print!("{}  ", mat[i][j].iy);
    //     }
    //     println!();
    //     print!("Ix2:\t");
    //     for j in 0..q_len + 1 {
    //         print!("{}  ", mat[i][j].ix2);
    //     }
    //     println!();
    //     println!();
    //     // for _ in 0..q_len + 1 {
    //     //     print!("--------");
    //     // }
    //     // println!();
    // }


    // Trace back
    let mut aligned_query = String::new();
    let mut ref_target = String::new();
    let mut major_target = String::new();

    let mut i = t_len;
    let mut j = q_len;
    let alignment_score = mat[i][j].m.max(mat[i][j].ix.max(mat[i][j].iy.max(mat[i][j].ix2)));

    let mut trace_back_stat;

    if mat[i][j].m_prev_m {
        trace_back_stat = TraceBack::M;
    } else if mat[i][j].m_prev_ix {
        trace_back_stat = TraceBack::IX;
    } else if mat[i][j].m_prev_iy {
        trace_back_stat = TraceBack::IY;
    } else if mat[i][j].m_prev_ix2 {
        trace_back_stat = TraceBack::IX2;
    } else {
        panic!("Error: traceback");
    }

    while i > 0 && j > 0 {
        let qbase = query.chars().nth(j - 1).unwrap();
        let ref_base = profile[i - 1].get_ref_base() as char;
        let major_base = profile[i - 1].get_major_base() as char;
        if trace_back_stat == TraceBack::M {
            if mat[i][j].m_prev_m {
                aligned_query.push(qbase);
                ref_target.push(ref_base);
                major_target.push(major_base);
                i -= 1;
                j -= 1;
                trace_back_stat = TraceBack::M;
            } else if mat[i][j].m_prev_ix {
                trace_back_stat = TraceBack::IX;
            } else if mat[i][j].m_prev_iy {
                trace_back_stat = TraceBack::IY;
            } else if mat[i][j].m_prev_ix2 {
                trace_back_stat = TraceBack::IX2;
            }
        } else if trace_back_stat == TraceBack::IX {
            if mat[i][j].ix_prev_m {
                aligned_query.push('-');
                ref_target.push(ref_base);
                major_target.push(major_base);
                i -= 1;
                trace_back_stat = TraceBack::M;
            } else if (mat[i][j].ix_prev_ix) {
                aligned_query.push('-');
                ref_target.push(ref_base);
                major_target.push(major_base);
                i -= 1;
                trace_back_stat = TraceBack::IX;
            }
        } else if trace_back_stat == TraceBack::IY {
            println!("Error: dash can not appear on target. gap cost on target is infinity.");
            process::exit(1);
            if mat[i][j].iy_prev_m {
                aligned_query.push(qbase);
                ref_target.push('-');
                major_target.push('-');
                j -= 1;
                trace_back_stat = TraceBack::M;
            } else if mat[i][j].iy_prev_iy {
                aligned_query.push(qbase);
                ref_target.push('-');
                major_target.push('-');
                j -= 1;
                trace_back_stat = TraceBack::IY;
            }
        } else if trace_back_stat == TraceBack::IX2 {
            if mat[i][j].ix2_prev_m {
                aligned_query.push('N');
                ref_target.push(ref_base);
                major_target.push(major_base);
                i -= 1;
                trace_back_stat = TraceBack::M;
            } else if mat[i][j].ix2_prev_ix2 {
                aligned_query.push('N');
                ref_target.push(ref_base);
                major_target.push(major_base);
                i -= 1;
                trace_back_stat = TraceBack::IX2;
            }
        }
    }

    while i > 0 {
        let ref_base = profile[i - 1].get_ref_base() as char;
        let major_base = profile[i - 1].get_major_base() as char;
        aligned_query.push(' ');
        ref_target.push(ref_base);
        major_target.push(major_base);
        i -= 1;
    }
    while j > 0 {
        let qbase = query.chars().nth(j - 1).unwrap();
        aligned_query.push(qbase);
        ref_target.push(' ');
        major_target.push(' ');
        j -= 1;
    }
    aligned_query = aligned_query.chars().rev().collect();
    ref_target = ref_target.chars().rev().collect();
    major_target = major_target.chars().rev().collect();


    (alignment_score, aligned_query, ref_target, major_target)
}