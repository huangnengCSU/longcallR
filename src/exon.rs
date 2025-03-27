// use std::collections::HashMap;

#[derive(Debug, Clone, Default, PartialEq, Eq, Hash)]
pub struct Exon {
    pub chr: String,
    pub start: i64,
    // start position on the reference, 0-based, inclusive
    pub end: i64,
    // end position on the reference, 0-based, exclusive
    pub state: u8, // 0: start exon, 1: internal exon, 2: end exon, 3: whole read is a single exon
}

// pub fn exon_cluster(
//     mut exons: Vec<Exon>,
//     smallest_start: i64,
//     largest_end: i64,
//     min_sup: i32,
// ) -> HashMap<Exon, Vec<Exon>> {
//     let mut cover_vec = vec![0; (largest_end - smallest_start) as usize];
//     let mut cover_exon_idx: Vec<Vec<usize>> = vec![Vec::new(); (largest_end - smallest_start) as usize];
//     for idx in 0..exons.len() {
//         let e = &exons[idx];
//         for j in e.start..e.end {
//             cover_vec[j as usize - smallest_start as usize] += 1;
//         }
//         cover_exon_idx[e.start as usize - smallest_start as usize].push(idx);
//     }
//
//     // for overlapped exons, merge them into one cluster
//     let mut clustersI: Vec<Vec<Exon>> = Vec::new();
//     let mut ts = -1; // index in cover_vec
//     let mut te = -1; // index in cover_vec
//     for i in 0..cover_vec.len() {
//         if cover_vec[i] == 0 {
//             if ts >= 0 && te >= 0 {
//                 let mut t_cluster: Vec<Exon> = Vec::new();
//                 for ti in (ts as usize)..=(te as usize) {
//                     for idx in cover_exon_idx[ti].iter() {
//                         t_cluster.push(exons[*idx].clone());
//                     }
//                 }
//                 clustersI.push(t_cluster);
//                 ts = -1;
//                 te = -1;
//             }
//             continue;
//         } else {
//             if ts == -1 && te == -1 {
//                 ts = i as i64;
//                 te = i as i64;
//             } else {
//                 te = i as i64;
//             }
//         }
//     }
//     if ts >= 0 && te >= 0 {
//         let mut t_cluster: Vec<Exon> = Vec::new();
//         for ti in (ts as usize)..=(te as usize) {
//             for idx in cover_exon_idx[ti].iter() {
//                 t_cluster.push(exons[*idx].clone());
//             }
//         }
//         clustersI.push(t_cluster);
//     }
//
//     // for each level I cluster, divide them into level II clusters (as hierarchical clustering)
//     let mut clusters: HashMap<Exon, Vec<Exon>> = HashMap::new();
//     for c in clustersI.iter() {
//         let mut exon_hashmap: HashMap<Exon, Vec<Exon>> = HashMap::new();
//         for e in c.iter() {
//             let mut texon = e.clone();
//             if e.state == 0 {
//                 // start exon ignores the start position when matching exons, because the start position of start exon is not consistent
//                 texon.start = -1;
//             } else if e.state == 2 {
//                 // end exon ignores the end position when matching exons, because the end position of end exon is not consistent
//                 texon.end = -1;
//             }
//             if exon_hashmap.contains_key(&texon) {
//                 let exon_vec = exon_hashmap.get_mut(&texon).unwrap();
//                 exon_vec.push(e.clone());
//             } else {
//                 exon_hashmap.insert(texon.clone(), vec![e.clone()]);
//             }
//         }
//
//         // ignore exons without enough support
//         let mut filter_exons: Vec<Exon> = Vec::new();
//         for (e, v) in exon_hashmap.iter() {
//             if (v.len() as i32) < min_sup {
//                 filter_exons.push(e.clone());
//             }
//         }
//         for e in filter_exons.iter() {
//             exon_hashmap.remove(e);
//         }
//         for (e, v) in exon_hashmap.iter() {
//             clusters.insert(e.clone(), v.clone());
//         }
//     }
//     clusters
// }