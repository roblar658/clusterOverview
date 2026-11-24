use std::collections::{HashMap, HashSet, VecDeque};

use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;

use crate::algorithms::cat1_centroid::{init_kmeans_pp, run_lloyd_from_centers};
use crate::algorithms::cat2_density::run_dbscan_core;
use crate::algorithms::{AlgorithmOutput, RunConfig};
use crate::data::PreparedData;
use crate::math::{
    estimate_eps_kdist, euclidean, euclidean_sq, jaccard_sim_cat, median_pairwise_dist,
};

pub fn run(id: usize, data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    match id {
        56 => lance_williams_run(data, cfg, Linkage::Single, "Single Linkage (SLINK)"),
        57 => lance_williams_run(data, cfg, Linkage::Complete, "Complete Linkage (CLINK)"),
        58 => lance_williams_run(data, cfg, Linkage::Average, "Average Linkage (UPGMA)"),
        59 => lance_williams_run(data, cfg, Linkage::Ward, "Ward's Method"),
        60 => lance_williams_run(data, cfg, Linkage::Wpgma, "WPGMA"),
        61 => lance_williams_run(data, cfg, Linkage::Centroid, "Centroid Linkage (UPGMC)"),
        62 => lance_williams_run(data, cfg, Linkage::Median, "Median Linkage (WPGMC)"),
        63 => diana(data, cfg),
        64 => birch(data, cfg),
        65 => cure(data, cfg),
        66 => chameleon(data, cfg),
        67 => rock(data, cfg),
        68 => mona(data, cfg),
        69 => genie(data, cfg),
        70 => hca_dbscan(data, cfg),
        71 => bubble(data, cfg),
        72 => lance_williams_general(data, cfg),
        73 => minimax_linkage(data, cfg),
        74 => lance_williams_run(data, cfg, Linkage::FlexibleBeta(-0.25), "Flexible Beta (beta=-0.25)"),
        75 => neighbor_joining(data, cfg),
        76 => perch(data, cfg),
        77 => grinch(data, cfg),
        78 => bico(data, cfg),
        79 => scc_reciprocal_nn(data, cfg),
        80 => finch(data, cfg),
        _ => unreachable!(),
    }
}

#[derive(Clone, Copy)]
pub enum Linkage {
    Single,
    Complete,
    Average,
    Ward,
    Wpgma,
    Centroid,
    Median,
    FlexibleBeta(f64),
    Custom { alpha_i: f64, alpha_j: f64, beta: f64, gamma: f64 },
}

/// Unified Lance-Williams agglomerative hierarchical clustering engine.
pub fn agglomerative_lance_williams(
    dist: &[Vec<f64>],
    k: usize,
    linkage: Linkage,
) -> Vec<i32> {
    let n = dist.len();
    let kk = k.clamp(1, n);
    let mut active: Vec<bool> = vec![true; n];
    let mut members: Vec<Vec<usize>> = (0..n).map(|i| vec![i]).collect();
    let mut dmat: Vec<Vec<f64>> = match linkage {
        Linkage::Ward | Linkage::Centroid | Linkage::Median => dist
            .iter()
            .map(|r| r.iter().map(|&v| v * v).collect())
            .collect(),
        _ => dist.to_vec(),
    };

    let mut remaining = n;
    while remaining > kk {
        let mut best_a = 0;
        let mut best_b = 0;
        let mut best_d = f64::INFINITY;
        for i in 0..n {
            if !active[i] {
                continue;
            }
            for j in (i + 1)..n {
                if !active[j] {
                    continue;
                }
                if dmat[i][j] < best_d {
                    best_d = dmat[i][j];
                    best_a = i;
                    best_b = j;
                }
            }
        }
        if !best_d.is_finite() {
            break;
        }

        let ni = members[best_a].len() as f64;
        let nj = members[best_b].len() as f64;
        let d_ij = dmat[best_a][best_b];

        // Update distances d(i∪j, c) via Lance-Williams recurrence
        for c in 0..n {
            if !active[c] || c == best_a || c == best_b {
                continue;
            }
            let nc = members[c].len() as f64;
            let dic = dmat[best_a][c];
            let djc = dmat[best_b][c];
            let new_d = match linkage {
                Linkage::Single => dic.min(djc),
                Linkage::Complete => dic.max(djc),
                Linkage::Average => (ni * dic + nj * djc) / (ni + nj),
                Linkage::Wpgma => 0.5 * (dic + djc),
                Linkage::Ward => {
                    let tot = ni + nj + nc;
                    ((ni + nc) * dic + (nj + nc) * djc - nc * d_ij) / tot
                }
                Linkage::Centroid => {
                    let s = ni + nj;
                    (ni * dic + nj * djc) / s - (ni * nj * d_ij) / (s * s)
                }
                Linkage::Median => 0.5 * dic + 0.5 * djc - 0.25 * d_ij,
                Linkage::FlexibleBeta(beta) => {
                    let alpha = 0.5 * (1.0 - beta);
                    alpha * dic + alpha * djc + beta * d_ij
                }
                Linkage::Custom { alpha_i, alpha_j, beta, gamma } => {
                    alpha_i * dic + alpha_j * djc + beta * d_ij + gamma * (dic - djc).abs()
                }
            }
            .max(0.0);
            dmat[best_a][c] = new_d;
            dmat[c][best_a] = new_d;
        }

        let moved = std::mem::take(&mut members[best_b]);
        members[best_a].extend(moved);
        active[best_b] = false;
        remaining -= 1;
    }

    let mut labels = vec![0i32; n];
    let mut cid = 0i32;
    for i in 0..n {
        if active[i] {
            for &p in &members[i] {
                labels[p] = cid;
            }
            cid += 1;
        }
    }
    labels
}

fn lance_williams_run(
    data: &PreparedData,
    cfg: &RunConfig,
    linkage: Linkage,
    tag: &str,
) -> AlgorithmOutput {
    let k = cfg.effective_k(data.n_samples);
    let labels = agglomerative_lance_williams(&data.dist_matrix, k, linkage);
    AlgorithmOutput::new(
        labels,
        format!("k={}, {}", k, tag),
        format!("Agglomerativt dendrogram kuttet ved k={}", k),
    )
}

// 63. DIANA (Kaufman & Rousseeuw 1990 — Divisive Analysis splitting largest-diameter cluster)
fn diana(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let k = cfg.effective_k(n);
    let mut clusters: Vec<Vec<usize>> = vec![(0..n).collect()];

    while clusters.len() < k {
        // Choose cluster with largest diameter
        let split_idx = clusters
            .iter()
            .enumerate()
            .filter(|(_, m)| m.len() >= 2)
            .max_by(|(_, a), (_, b)| {
                let da = cluster_diameter(&data.dist_matrix, a);
                let db = cluster_diameter(&data.dist_matrix, b);
                da.total_cmp(&db)
            })
            .map(|(idx, _)| idx);

        let Some(idx) = split_idx else { break };
        let mut main_c = clusters.remove(idx);

        // Find seed of splinter group: point with highest average dissimilarity to rest of cluster
        let splinter_seed_pos = (0..main_c.len())
            .max_by(|&a, &b| {
                let avg_a: f64 = main_c.iter().map(|&j| data.dist_matrix[main_c[a]][j]).sum();
                let avg_b: f64 = main_c.iter().map(|&j| data.dist_matrix[main_c[b]][j]).sum();
                avg_a.total_cmp(&avg_b)
            })
            .unwrap_or(0);
        let mut splinter = vec![main_c.remove(splinter_seed_pos)];

        // Iteratively move points closer on average to splinter than to remaining main_c
        loop {
            if main_c.len() <= 1 {
                break;
            }
            let mut best_pos = None;
            let mut best_diff = 0.0;
            for (pos, &p) in main_c.iter().enumerate() {
                let d_main = main_c
                    .iter()
                    .filter(|&&q| q != p)
                    .map(|&q| data.dist_matrix[p][q])
                    .sum::<f64>()
                    / (main_c.len() - 1) as f64;
                let d_spl = splinter
                    .iter()
                    .map(|&q| data.dist_matrix[p][q])
                    .sum::<f64>()
                    / splinter.len() as f64;
                let diff = d_main - d_spl;
                if diff > best_diff {
                    best_diff = diff;
                    best_pos = Some(pos);
                }
            }
            if let Some(pos) = best_pos {
                splinter.push(main_c.remove(pos));
            } else {
                break;
            }
        }
        clusters.push(main_c);
        clusters.push(splinter);
    }

    let mut labels = vec![0i32; n];
    for (cid, mems) in clusters.iter().enumerate() {
        for &p in mems {
            labels[p] = cid as i32;
        }
    }
    AlgorithmOutput::new(
        labels,
        format!("k={}, top-down DIANA", k),
        format!("Divisiv splitt basert på klyngediameter til {} klynger", clusters.len()),
    )
}

fn cluster_diameter(dist: &[Vec<f64>], members: &[usize]) -> f64 {
    let mut diam = 0.0_f64;
    for i in 0..members.len() {
        for j in (i + 1)..members.len() {
            diam = diam.max(dist[members[i]][members[j]]);
        }
    }
    diam
}

// 64. BIRCH (Zhang et al. 1996 — Clustering Feature entries + global agglomerative phase)
struct CfEntry {
    n: usize,
    ls: Vec<f64>,
    ss: f64,
    members: Vec<usize>,
}

impl CfEntry {
    fn new(p: &[f64], idx: usize) -> Self {
        let ss: f64 = p.iter().map(|v| v * v).sum();
        Self { n: 1, ls: p.to_vec(), ss, members: vec![idx] }
    }
    fn centroid(&self) -> Vec<f64> {
        self.ls.iter().map(|&v| v / self.n as f64).collect()
    }
    fn radius_if_added(&self, p: &[f64]) -> f64 {
        let nn = (self.n + 1) as f64;
        let ss_new = self.ss + p.iter().map(|v| v * v).sum::<f64>();
        let ls_sq: f64 = self
            .ls
            .iter()
            .zip(p.iter())
            .map(|(a, b)| (a + b).powi(2))
            .sum();
        ((ss_new / nn) - ls_sq / (nn * nn)).max(0.0).sqrt()
    }
    fn add(&mut self, p: &[f64], idx: usize) {
        self.n += 1;
        for (a, &b) in self.ls.iter_mut().zip(p.iter()) {
            *a += b;
        }
        self.ss += p.iter().map(|v| v * v).sum::<f64>();
        self.members.push(idx);
    }
}

fn birch(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let k = cfg.effective_k(n);
    let threshold = median_pairwise_dist(&data.dist_matrix) * 0.35;

    let mut cfs: Vec<CfEntry> = Vec::new();
    for i in 0..n {
        let p = &data.x[i];
        let nearest = cfs.iter().enumerate().min_by(|(_, a), (_, b)| {
            euclidean_sq(p, &a.centroid()).total_cmp(&euclidean_sq(p, &b.centroid()))
        });
        if let Some((idx, entry)) = nearest {
            if entry.radius_if_added(p) <= threshold {
                cfs[idx].add(p, i);
                continue;
            }
        }
        cfs.push(CfEntry::new(p, i));
    }

    // Global Ward clustering on CF centroids
    let m = cfs.len();
    let cf_centroids: Vec<Vec<f64>> = cfs.iter().map(|c| c.centroid()).collect();
    let mut cf_dist = vec![vec![0.0; m]; m];
    for a in 0..m {
        for b in (a + 1)..m {
            let dd = euclidean(&cf_centroids[a], &cf_centroids[b]);
            cf_dist[a][b] = dd;
            cf_dist[b][a] = dd;
        }
    }
    let cf_labels = if m > k {
        agglomerative_lance_williams(&cf_dist, k, Linkage::Ward)
    } else {
        (0..m as i32).collect()
    };

    let mut labels = vec![0i32; n];
    for (ci, entry) in cfs.iter().enumerate() {
        for &p in &entry.members {
            labels[p] = cf_labels[ci];
        }
    }
    AlgorithmOutput::new(
        labels,
        format!("k={}, T={:.2}, CF-noder={}", k, threshold, m),
        format!("1-pass CF-tre komprimerte {} punkter til {} subklynger", n, m),
    )
}

// 65. CURE (Guha et al. 1998 — well-scattered representatives shrunk toward cluster centroid)
fn cure(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let c_reps = 4usize;
    let alpha = 0.25_f64;

    struct CureCluster {
        members: Vec<usize>,
        reps: Vec<Vec<f64>>,
    }

    let make_reps = |members: &[usize]| -> Vec<Vec<f64>> {
        let mut mean = vec![0.0; d];
        for &m in members {
            for j in 0..d {
                mean[j] += data.x[m][j];
            }
        }
        for v in &mut mean {
            *v /= members.len() as f64;
        }
        if members.len() <= c_reps {
            return members
                .iter()
                .map(|&m| (0..d).map(|j| data.x[m][j] + alpha * (mean[j] - data.x[m][j])).collect())
                .collect();
        }
        // Pick scattered points farthest from mean, then farthest from already chosen
        let mut chosen = Vec::with_capacity(c_reps);
        let first = *members
            .iter()
            .max_by(|&&a, &&b| euclidean_sq(&data.x[a], &mean).total_cmp(&euclidean_sq(&data.x[b], &mean)))
            .unwrap();
        chosen.push(first);
        while chosen.len() < c_reps {
            let next = *members
                .iter()
                .filter(|m| !chosen.contains(m))
                .max_by(|&&a, &&b| {
                    let da = chosen.iter().map(|&c| data.dist_matrix[a][c]).fold(f64::INFINITY, f64::min);
                    let db = chosen.iter().map(|&c| data.dist_matrix[b][c]).fold(f64::INFINITY, f64::min);
                    da.total_cmp(&db)
                })
                .unwrap();
            chosen.push(next);
        }
        chosen
            .into_iter()
            .map(|m| (0..d).map(|j| data.x[m][j] + alpha * (mean[j] - data.x[m][j])).collect())
            .collect()
    };

    let rep_dist = |a: &[Vec<f64>], b: &[Vec<f64>]| -> f64 {
        let mut md = f64::INFINITY;
        for ra in a {
            for rb in b {
                let dd = euclidean(ra, rb);
                if dd < md {
                    md = dd;
                }
            }
        }
        md
    };

    // Start from ~min(n, 35) micro-clusters via Ward for speed on larger datasets, then CURE merge
    let init_k = (k * 6).clamp(12, 32).min(n);
    let init_lbl = if n > init_k {
        agglomerative_lance_williams(&data.dist_matrix, init_k, Linkage::Average)
    } else {
        (0..n as i32).collect()
    };
    let mut map: HashMap<i32, Vec<usize>> = HashMap::new();
    for (i, &l) in init_lbl.iter().enumerate() {
        map.entry(l).or_default().push(i);
    }
    let mut clusters: Vec<CureCluster> = map
        .into_values()
        .map(|m| {
            let reps = make_reps(&m);
            CureCluster { members: m, reps }
        })
        .collect();

    while clusters.len() > k {
        let mut best = (0, 1, f64::INFINITY);
        for i in 0..clusters.len() {
            for j in (i + 1)..clusters.len() {
                let dd = rep_dist(&clusters[i].reps, &clusters[j].reps);
                if dd < best.2 {
                    best = (i, j, dd);
                }
            }
        }
        let b = clusters.remove(best.1);
        clusters[best.0].members.extend(b.members);
        clusters[best.0].reps = make_reps(&clusters[best.0].members);
    }

    let mut labels = vec![0i32; n];
    for (cid, cl) in clusters.iter().enumerate() {
        for &p in &cl.members {
            labels[p] = cid as i32;
        }
    }
    AlgorithmOutput::new(
        labels,
        format!("k={}, reps={}, krymp={:.2}", k, c_reps, alpha),
        "Spredte representanter krympet 25% mot sentroidene".to_string(),
    )
}

// 66. CHAMELEON (Karypis et al. 1999 — k-NN graph partitioning + RI * RC^alpha dynamic merging)
fn chameleon(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let k = cfg.effective_k(n);
    let alpha = 2.0_f64;

    // Phase 1: Fine-grained sub-clusters on k-NN graph via average linkage
    let sub_k = (k * 4).clamp(8, 24).min(n / 2).max(k);
    let init_lbl = agglomerative_lance_williams(&data.dist_matrix, sub_k, Linkage::Ward);
    let mut map: HashMap<i32, Vec<usize>> = HashMap::new();
    for (i, &l) in init_lbl.iter().enumerate() {
        map.entry(l).or_default().push(i);
    }
    let mut clusters: Vec<Vec<usize>> = map.into_values().collect();

    // Build adjacency weight lookup from shared k-NN graph
    let mut wmat = vec![vec![0.0; n]; n];
    for i in 0..n {
        for &(j, w) in &data.knn_graph[i] {
            wmat[i][j] = w;
        }
    }

    let internal_stats = |mems: &[usize]| -> (f64, f64) {
        let mut sum_w = 0.0;
        let mut cnt = 0usize;
        for a in 0..mems.len() {
            for b in (a + 1)..mems.len() {
                let w = wmat[mems[a]][mems[b]];
                if w > 0.0 {
                    sum_w += w;
                    cnt += 1;
                }
            }
        }
        let avg_w = if cnt > 0 { sum_w / cnt as f64 } else { 1e-4 };
        (sum_w.max(1e-4), avg_w)
    };

    // Phase 2: Merge pair maximizing RI(Ci, Cj) * RC(Ci, Cj)^alpha
    while clusters.len() > k {
        let mut best_pair = (0, 1);
        let mut best_score = f64::NEG_INFINITY;

        for a in 0..clusters.len() {
            let (ec_a, avg_a) = internal_stats(&clusters[a]);
            for b in (a + 1)..clusters.len() {
                let (ec_b, avg_b) = internal_stats(&clusters[b]);
                let mut ec_ab = 0.0;
                let mut cnt_ab = 0usize;
                for &u in &clusters[a] {
                    for &v in &clusters[b] {
                        let w = wmat[u][v];
                        if w > 0.0 {
                            ec_ab += w;
                            cnt_ab += 1;
                        }
                    }
                }
                let score = if cnt_ab > 0 {
                    let ri = ec_ab / (0.5 * (ec_a + ec_b));
                    let avg_ab = ec_ab / cnt_ab as f64;
                    let na = clusters[a].len() as f64;
                    let nb = clusters[b].len() as f64;
                    let denom = (na * avg_a + nb * avg_b) / (na + nb);
                    let rc = avg_ab / denom.max(1e-6);
                    ri * rc.powf(alpha)
                } else {
                    // Fallback negative distance if disconnected in sparse kNN graph
                    let mut min_d = f64::INFINITY;
                    for &u in &clusters[a] {
                        for &v in &clusters[b] {
                            min_d = min_d.min(data.dist_matrix[u][v]);
                        }
                    }
                    -min_d * 1e-3
                };
                if score > best_score {
                    best_score = score;
                    best_pair = (a, b);
                }
            }
        }
        let right = clusters.remove(best_pair.1);
        clusters[best_pair.0].extend(right);
    }

    let mut labels = vec![0i32; n];
    for (cid, mems) in clusters.iter().enumerate() {
        for &p in mems {
            labels[p] = cid as i32;
        }
    }
    AlgorithmOutput::new(
        labels,
        format!("k={}, sub={}, alpha={:.1}", k, sub_k, alpha),
        "2-stegs k-NN deling + dynamisk RI*RC^2 sammenslåing".to_string(),
    )
}

// 67. ROCK (Guha et al. 1999 — RObust Clustering using linKs for categorical attributes)
fn rock(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let k = cfg.effective_k(n);
    let theta = 0.35_f64;
    let f_theta = (1.0 - theta) / (1.0 + theta);

    // Adjacency matrix A_ij = 1 if Jaccard(i,j) >= theta
    let mut nbrs: Vec<HashSet<usize>> = vec![HashSet::new(); n];
    for i in 0..n {
        for j in (i + 1)..n {
            if jaccard_sim_cat(&data.categorical[i], &data.categorical[j]) >= theta {
                nbrs[i].insert(j);
                nbrs[j].insert(i);
            }
        }
    }

    // Link matrix: links[i][j] = |nbrs[i] ∩ nbrs[j]|
    let mut links = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let l = nbrs[i].intersection(&nbrs[j]).count() as f64;
            links[i][j] = l;
            links[j][i] = l;
        }
    }

    // Agglomerative merge maximizing ROCK goodness measure g(Ci, Cj)
    let mut active = vec![true; n];
    let mut members: Vec<Vec<usize>> = (0..n).map(|i| vec![i]).collect();
    let mut remaining = n;
    let exp_p = 1.0 + 2.0 * f_theta;

    while remaining > k {
        let mut best_pair = None;
        let mut best_g = f64::NEG_INFINITY;
        for i in 0..n {
            if !active[i] {
                continue;
            }
            let ni = members[i].len() as f64;
            for j in (i + 1)..n {
                if !active[j] {
                    continue;
                }
                let nj = members[j].len() as f64;
                let denom = ((ni + nj).powf(exp_p) - ni.powf(exp_p) - nj.powf(exp_p)).max(1e-6);
                let g = links[i][j] / denom;
                if g > best_g {
                    best_g = g;
                    best_pair = Some((i, j));
                }
            }
        }
        let Some((a, b)) = best_pair else { break };
        for c in 0..n {
            if active[c] && c != a && c != b {
                links[a][c] += links[b][c];
                links[c][a] = links[a][c];
            }
        }
        let moved = std::mem::take(&mut members[b]);
        members[a].extend(moved);
        active[b] = false;
        remaining -= 1;
    }

    let mut labels = vec![0i32; n];
    let mut cid = 0i32;
    for i in 0..n {
        if active[i] {
            for &p in &members[i] {
                labels[p] = cid;
            }
            cid += 1;
        }
    }
    AlgorithmOutput::new(
        labels,
        format!("k={}, theta={:.2}", k, theta),
        "Kategorisk hierarkisk sammenslåing via felles nabo-lenker".to_string(),
    )
}

// 68. MONA (Kaufman & Rousseeuw 1990 — Monothetic Analysis divisive binary clustering)
fn mona(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.x_binary[0].len();
    let k = cfg.effective_k(n);
    let mut clusters: Vec<Vec<usize>> = vec![(0..n).collect()];
    let mut split_vars = Vec::new();

    while clusters.len() < k {
        // Choose largest cluster that can be split on a binary variable
        let split_idx = clusters
            .iter()
            .enumerate()
            .filter(|(_, m)| m.len() >= 4)
            .max_by_key(|(_, m)| m.len())
            .map(|(i, _)| i);
        let Some(idx) = split_idx else { break };
        let members = clusters.remove(idx);

        // Find binary attribute j maximizing total association sum_{f != j} |a*d - b*c|
        let mut best_var = None;
        let mut best_assoc = -1.0_f64;
        for j in 0..d {
            let n1 = members.iter().filter(|&&m| data.x_binary[m][j] == 1).count();
            if n1 == 0 || n1 == members.len() {
                continue;
            }
            let mut assoc = 0.0;
            for f in 0..d {
                if f == j {
                    continue;
                }
                let mut a = 0.0_f64;
                let mut b = 0.0_f64;
                let mut c = 0.0_f64;
                let mut dd = 0.0_f64;
                for &m in &members {
                    match (data.x_binary[m][j], data.x_binary[m][f]) {
                        (1, 1) => a += 1.0,
                        (1, 0) => b += 1.0,
                        (0, 1) => c += 1.0,
                        _ => dd += 1.0,
                    }
                }
                assoc += (a * dd - b * c).abs();
            }
            if assoc > best_assoc {
                best_assoc = assoc;
                best_var = Some(j);
            }
        }

        if let Some(j) = best_var {
            let left: Vec<usize> = members.iter().copied().filter(|&m| data.x_binary[m][j] == 0).collect();
            let right: Vec<usize> = members.iter().copied().filter(|&m| data.x_binary[m][j] == 1).collect();
            split_vars.push(j);
            clusters.push(left);
            clusters.push(right);
        } else {
            clusters.push(members);
            break;
        }
    }

    let mut labels = vec![0i32; n];
    for (cid, mems) in clusters.iter().enumerate() {
        for &p in mems {
            labels[p] = cid as i32;
        }
    }
    AlgorithmOutput::new(
        labels,
        format!("k={}, binære splitter={}", k, split_vars.len()),
        "Monotetisk binær splitt etter maksimal variabelassosiasjon".to_string(),
    )
}

// 69. Genie (Gagolewski et al. 2016 — Single Linkage constrained by Gini index threshold g)
fn genie(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let k = cfg.effective_k(n);
    let gini_thresh = 0.30_f64;

    let mut active = vec![true; n];
    let mut members: Vec<Vec<usize>> = (0..n).map(|i| vec![i]).collect();
    let mut dmat = data.dist_matrix.clone();
    let mut remaining = n;
    let mut gini_triggers = 0usize;

    let compute_gini = |act: &[bool], mems: &[Vec<usize>]| -> f64 {
        let mut sizes: Vec<f64> = act
            .iter()
            .enumerate()
            .filter(|&(_, &a)| a)
            .map(|(i, _)| mems[i].len() as f64)
            .collect();
        let m = sizes.len();
        if m <= 1 {
            return 0.0;
        }
        sizes.sort_by(|a, b| a.total_cmp(b));
        let sum: f64 = sizes.iter().sum();
        let mut num = 0.0;
        for (i, &s) in sizes.iter().enumerate() {
            num += (2.0 * (i + 1) as f64 - m as f64 - 1.0) * s;
        }
        (num / ((m - 1) as f64 * sum)).clamp(0.0, 1.0)
    };

    while remaining > k {
        let g = compute_gini(&active, &members);
        let (a, b) = if g > gini_thresh {
            gini_triggers += 1;
            // Force merging the smallest active cluster with its nearest neighbor
            let min_size = (0..n)
                .filter(|&i| active[i])
                .map(|i| members[i].len())
                .min()
                .unwrap_or(1);
            let mut best = (0, 0, f64::INFINITY);
            for i in 0..n {
                if !active[i] || members[i].len() != min_size {
                    continue;
                }
                for j in 0..n {
                    if !active[j] || i == j {
                        continue;
                    }
                    if dmat[i][j] < best.2 {
                        best = (i.min(j), i.max(j), dmat[i][j]);
                    }
                }
            }
            (best.0, best.1)
        } else {
            let mut best = (0, 0, f64::INFINITY);
            for i in 0..n {
                if !active[i] {
                    continue;
                }
                for j in (i + 1)..n {
                    if !active[j] {
                        continue;
                    }
                    if dmat[i][j] < best.2 {
                        best = (i, j, dmat[i][j]);
                    }
                }
            }
            (best.0, best.1)
        };
        if a == b {
            break;
        }
        for c in 0..n {
            if active[c] && c != a && c != b {
                let nd = dmat[a][c].min(dmat[b][c]);
                dmat[a][c] = nd;
                dmat[c][a] = nd;
            }
        }
        let moved = std::mem::take(&mut members[b]);
        members[a].extend(moved);
        active[b] = false;
        remaining -= 1;
    }

    let mut labels = vec![0i32; n];
    let mut cid = 0i32;
    for i in 0..n {
        if active[i] {
            for &p in &members[i] {
                labels[p] = cid;
            }
            cid += 1;
        }
    }
    AlgorithmOutput::new(
        labels,
        format!("k={}, Gini-terskel={:.2}", k, gini_thresh),
        format!("Single Linkage med {} Gini-balanserte sammenslåinger", gini_triggers),
    )
}

// 70. HCA-DBSCAN (Hybrid Hierarchical Clustering on DBSCAN density cores)
fn hca_dbscan(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let k = cfg.effective_k(n);
    let min_pts = 3usize.min(n / 4).max(2);
    let eps = estimate_eps_kdist(&data.dist_matrix, min_pts, 0.42);
    let core_lbl = run_dbscan_core(&data.dist_matrix, eps, min_pts);

    let mut groups: HashMap<i32, Vec<usize>> = HashMap::new();
    let mut noise = Vec::new();
    for (i, &l) in core_lbl.iter().enumerate() {
        if l >= 0 {
            groups.entry(l).or_default().push(i);
        } else {
            noise.push(i);
        }
    }
    if groups.len() <= 1 {
        let labels = agglomerative_lance_williams(&data.dist_matrix, k, Linkage::Average);
        return AlgorithmOutput::new(
            labels,
            format!("k={}, eps={:.2}", k, eps),
            "Fallback til hierarkisk kobling (få tetthetskjerner)".to_string(),
        );
    }

    let cores: Vec<Vec<usize>> = groups.into_values().collect();
    let m = cores.len();
    let mut cdist = vec![vec![0.0; m]; m];
    for a in 0..m {
        for b in (a + 1)..m {
            let mut sum = 0.0;
            for &u in &cores[a] {
                for &v in &cores[b] {
                    sum += data.dist_matrix[u][v];
                }
            }
            let avg = sum / (cores[a].len() * cores[b].len()) as f64;
            cdist[a][b] = avg;
            cdist[b][a] = avg;
        }
    }
    let merged = if m > k {
        agglomerative_lance_williams(&cdist, k, Linkage::Average)
    } else {
        (0..m as i32).collect()
    };

    let mut labels = vec![-1i32; n];
    for (ci, pts) in cores.iter().enumerate() {
        for &p in pts {
            labels[p] = merged[ci];
        }
    }
    // Attach noise points to nearest merged core if within 2*eps
    for &p in &noise {
        if let Some(nn) = (0..n)
            .filter(|&j| labels[j] >= 0)
            .min_by(|&a, &b| data.dist_matrix[p][a].total_cmp(&data.dist_matrix[p][b]))
        {
            if data.dist_matrix[p][nn] <= eps * 2.2 {
                labels[p] = labels[nn];
            }
        }
    }
    AlgorithmOutput::new(
        labels,
        format!("k={}, {} tetthetskjerner", k, m),
        "Hierarkisk kobling av DBSCAN-tetthetskjerner".to_string(),
    )
}

// 71. BUBBLE (Ganti et al. 1999 — metric-space extension of BIRCH using clustroids & row-sums without coordinates)
fn bubble(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let k = cfg.effective_k(n);
    let threshold = median_pairwise_dist(&data.dist_matrix) * 0.42;

    struct MetricCf {
        members: Vec<usize>,
        clustroid: usize,
        radius: f64,
    }
    let update_clustroid = |mems: &[usize]| -> (usize, f64) {
        let mut best_c = mems[0];
        let mut best_sum = f64::INFINITY;
        for &u in mems {
            let s: f64 = mems.iter().map(|&v| data.dist_matrix[u][v].powi(2)).sum();
            if s < best_sum {
                best_sum = s;
                best_c = u;
            }
        }
        let rad = (best_sum / mems.len().max(1) as f64).sqrt();
        (best_c, rad)
    };

    let mut nodes: Vec<MetricCf> = Vec::new();
    for i in 0..n {
        let nearest = nodes.iter().enumerate().min_by(|(_, a), (_, b)| {
            data.dist_matrix[i][a.clustroid].total_cmp(&data.dist_matrix[i][b.clustroid])
        });
        if let Some((idx, node)) = nearest {
            if data.dist_matrix[i][node.clustroid] <= threshold && node.radius <= threshold {
                nodes[idx].members.push(i);
                let (cl, rad) = update_clustroid(&nodes[idx].members);
                nodes[idx].clustroid = cl;
                nodes[idx].radius = rad;
                continue;
            }
        }
        nodes.push(MetricCf { members: vec![i], clustroid: i, radius: 0.0 });
    }

    let m = nodes.len();
    let mut ndist = vec![vec![0.0; m]; m];
    for a in 0..m {
        for b in (a + 1)..m {
            let dd = data.dist_matrix[nodes[a].clustroid][nodes[b].clustroid];
            ndist[a][b] = dd;
            ndist[b][a] = dd;
        }
    }
    let nlbl = if m > k {
        agglomerative_lance_williams(&ndist, k, Linkage::Average)
    } else {
        (0..m as i32).collect()
    };

    let mut labels = vec![0i32; n];
    for (ci, node) in nodes.iter().enumerate() {
        for &p in &node.members {
            labels[p] = nlbl[ci];
        }
    }
    AlgorithmOutput::new(
        labels,
        format!("k={}, {} metriske CF-klynger", k, m),
        "Koordinatfri BUBBLE med metriske clustroids".to_string(),
    )
}

// 72. Lance-Williams algoritme (general recursive formula with balanced alpha/beta/gamma parameters)
fn lance_williams_general(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let k = cfg.effective_k(data.n_samples);
    let linkage = Linkage::Custom {
        alpha_i: 0.5,
        alpha_j: 0.5,
        beta: -0.1,
        gamma: 0.25, // mildly space-dilating compromise between average and complete
    };
    let labels = agglomerative_lance_williams(&data.dist_matrix, k, linkage);
    AlgorithmOutput::new(
        labels,
        format!("k={}, a=0.5, b=-0.1, g=0.25", k),
        "Generell Lance-Williams rekursiv formel".to_string(),
    )
}

// 73. Minimax Linkage (Bien & Tibshirani 2011 — minimizes cluster minimax radius)
fn minimax_linkage(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let k = cfg.effective_k(n);
    let mut clusters: Vec<Vec<usize>> = (0..n).map(|i| vec![i]).collect();

    let minimax_radius = |a: &[usize], b: &[usize]| -> f64 {
        let mut best_r = f64::INFINITY;
        for &p in a.iter().chain(b.iter()) {
            let mut max_d = 0.0_f64;
            for &q in a.iter().chain(b.iter()) {
                max_d = max_d.max(data.dist_matrix[p][q]);
            }
            if max_d < best_r {
                best_r = max_d;
            }
        }
        best_r
    };

    // To remain super fast for n up to 240, cache pair radii and update only merged cluster's row
    let mut active = vec![true; n];
    let mut rmat = data.dist_matrix.clone();
    let mut remaining = n;

    while remaining > k {
        let mut best = (0, 0, f64::INFINITY);
        for i in 0..n {
            if !active[i] {
                continue;
            }
            for j in (i + 1)..n {
                if !active[j] {
                    continue;
                }
                if rmat[i][j] < best.2 {
                    best = (i, j, rmat[i][j]);
                }
            }
        }
        let (a, b, _) = best;
        if a == b {
            break;
        }
        let moved = std::mem::take(&mut clusters[b]);
        clusters[a].extend(moved);
        active[b] = false;
        remaining -= 1;

        for c in 0..n {
            if active[c] && c != a {
                let r = minimax_radius(&clusters[a], &clusters[c]);
                rmat[a][c] = r;
                rmat[c][a] = r;
            }
        }
    }

    let mut labels = vec![0i32; n];
    let mut cid = 0i32;
    for i in 0..n {
        if active[i] {
            for &p in &clusters[i] {
                labels[p] = cid;
            }
            cid += 1;
        }
    }
    AlgorithmOutput::new(
        labels,
        format!("k={}, prototype-radius", k),
        "Minimax-kobling med garanterte klyngeprototyper".to_string(),
    )
}

// 75. Neighbor-Joining (Saitou & Nei 1987 — Q-matrix bottom-up phylogenetic tree construction)
fn neighbor_joining(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let k = cfg.effective_k(n);
    let mut active = vec![true; n];
    let mut members: Vec<Vec<usize>> = (0..n).map(|i| vec![i]).collect();
    let mut dmat = data.dist_matrix.clone();
    let mut remaining = n;

    while remaining > k && remaining > 2 {
        let active_ids: Vec<usize> = (0..n).filter(|&i| active[i]).collect();
        let r_cnt = active_ids.len() as f64;
        let mut row_sum = vec![0.0; n];
        for &i in &active_ids {
            row_sum[i] = active_ids.iter().map(|&j| dmat[i][j]).sum();
        }

        // Minimize Saitou-Nei Q(i,j) = (r - 2) * d(i,j) - R_i - R_j
        let mut best_a = active_ids[0];
        let mut best_b = active_ids[1];
        let mut best_q = f64::INFINITY;
        for idx_i in 0..active_ids.len() {
            let i = active_ids[idx_i];
            for idx_j in (idx_i + 1)..active_ids.len() {
                let j = active_ids[idx_j];
                let q = (r_cnt - 2.0) * dmat[i][j] - row_sum[i] - row_sum[j];
                if q < best_q {
                    best_q = q;
                    best_a = i;
                    best_b = j;
                }
            }
        }

        let d_ab = dmat[best_a][best_b];
        for &c in &active_ids {
            if c != best_a && c != best_b {
                let new_d = (0.5 * (dmat[best_a][c] + dmat[best_b][c] - d_ab)).max(0.0);
                dmat[best_a][c] = new_d;
                dmat[c][best_a] = new_d;
            }
        }
        let moved = std::mem::take(&mut members[best_b]);
        members[best_a].extend(moved);
        active[best_b] = false;
        remaining -= 1;
    }

    let mut labels = vec![0i32; n];
    let mut cid = 0i32;
    for i in 0..n {
        if active[i] {
            for &p in &members[i] {
                labels[p] = cid;
            }
            cid += 1;
        }
    }
    AlgorithmOutput::new(
        labels,
        format!("k={}, Saitou-Nei Q-matrise", k),
        "Fylogenetisk Neighbor-Joining uten ultrametrisk antagelse".to_string(),
    )
}

// 76. PERCH (Kobren et al. 2017 — Incremental hierarchical bounding-box tree with local rotations)
fn perch(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let max_leaves = (k * 5).clamp(12, 36).min(n);

    struct PerchLeaf {
        members: Vec<usize>,
        bmin: Vec<f64>,
        bmax: Vec<f64>,
    }
    impl PerchLeaf {
        fn new(p: &[f64], idx: usize) -> Self {
            Self { members: vec![idx], bmin: p.to_vec(), bmax: p.to_vec() }
        }
        fn box_dist(&self, p: &[f64]) -> f64 {
            let mut s = 0.0;
            for j in 0..p.len() {
                if p[j] < self.bmin[j] {
                    s += (self.bmin[j] - p[j]).powi(2);
                } else if p[j] > self.bmax[j] {
                    s += (p[j] - self.bmax[j]).powi(2);
                }
            }
            s.sqrt()
        }
        fn center(&self) -> Vec<f64> {
            self.bmin
                .iter()
                .zip(self.bmax.iter())
                .map(|(a, b)| 0.5 * (a + b))
                .collect()
        }
        fn insert(&mut self, p: &[f64], idx: usize) {
            self.members.push(idx);
            for j in 0..p.len() {
                self.bmin[j] = self.bmin[j].min(p[j]);
                self.bmax[j] = self.bmax[j].max(p[j]);
            }
        }
    }

    let mut leaves: Vec<PerchLeaf> = Vec::new();
    for i in 0..n {
        if leaves.len() < max_leaves {
            leaves.push(PerchLeaf::new(&data.x[i], i));
        } else {
            let best = (0..leaves.len())
                .min_by(|&a, &b| {
                    let da = leaves[a].box_dist(&data.x[i]) + 0.1 * euclidean(&leaves[a].center(), &data.x[i]);
                    let db = leaves[b].box_dist(&data.x[i]) + 0.1 * euclidean(&leaves[b].center(), &data.x[i]);
                    da.total_cmp(&db)
                })
                .unwrap_or(0);
            leaves[best].insert(&data.x[i], i);
        }
    }

    // Build hierarchy over bounding-box centers
    let m = leaves.len();
    let centers: Vec<Vec<f64>> = leaves
        .iter()
        .map(|lf| {
            let mut mu = vec![0.0; d];
            for &id in &lf.members {
                for j in 0..d {
                    mu[j] += data.x[id][j];
                }
            }
            for v in &mut mu {
                *v /= lf.members.len() as f64;
            }
            mu
        })
        .collect();
    let mut ldist = vec![vec![0.0; m]; m];
    for a in 0..m {
        for b in (a + 1)..m {
            let dd = euclidean(&centers[a], &centers[b]);
            ldist[a][b] = dd;
            ldist[b][a] = dd;
        }
    }
    let llbl = if m > k {
        agglomerative_lance_williams(&ldist, k, Linkage::Ward)
    } else {
        (0..m as i32).collect()
    };

    let mut labels = vec![0i32; n];
    for (ci, lf) in leaves.iter().enumerate() {
        for &p in &lf.members {
            labels[p] = llbl[ci];
        }
    }
    AlgorithmOutput::new(
        labels,
        format!("k={}, maks_blad={}", k, max_leaves),
        format!("Inkrementelt bounding-box-tre ({} bladnoder)", m),
    )
}

// 77. GRINCH (Monath et al. 2019 — Incremental hierarchical clustering with rotate & graft)
fn grinch(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 77);

    // Maintain incremental sub-clusters; perform online insert + rotate/graft consistency check
    let target_sub = (k * 4).clamp(10, 28).min(n);
    let mut clusters: Vec<Vec<usize>> = Vec::new();
    let mut centroids: Vec<Vec<f64>> = Vec::new();
    let mut grafts = 0usize;

    let recompute_cen = |mems: &[usize]| -> Vec<f64> {
        let mut mu = vec![0.0; d];
        for &p in mems {
            for j in 0..d {
                mu[j] += data.x[p][j];
            }
        }
        for v in &mut mu {
            *v /= mems.len().max(1) as f64;
        }
        mu
    };

    for i in 0..n {
        if clusters.len() < target_sub {
            clusters.push(vec![i]);
            centroids.push(data.x[i].clone());
        } else {
            let best = (0..clusters.len())
                .min_by(|&a, &b| {
                    euclidean_sq(&data.x[i], &centroids[a])
                        .total_cmp(&euclidean_sq(&data.x[i], &centroids[b]))
                })
                .unwrap_or(0);
            clusters[best].push(i);
            centroids[best] = recompute_cen(&clusters[best]);

            // Graft check: if two sub-clusters became closer than internal spread, graft & split worst
            if i % 8 == 7 && clusters.len() >= 3 {
                let (mut ca, mut cb, mut min_d) = (0, 1, f64::INFINITY);
                for a in 0..clusters.len() {
                    for b in (a + 1)..clusters.len() {
                        let dd = euclidean(&centroids[a], &centroids[b]);
                        if dd < min_d {
                            min_d = dd;
                            ca = a;
                            cb = b;
                        }
                    }
                }
                let worst = (0..clusters.len())
                    .filter(|&c| c != ca && c != cb && clusters[c].len() >= 4)
                    .max_by(|&a, &b| {
                        cluster_diameter(&data.dist_matrix, &clusters[a])
                            .total_cmp(&cluster_diameter(&data.dist_matrix, &clusters[b]))
                    });
                if let Some(wc) = worst {
                    if cluster_diameter(&data.dist_matrix, &clusters[wc]) > min_d * 2.2 {
                        // Graft cb into ca and split wc into 2
                        let moved = clusters[cb].clone();
                        clusters[ca].extend(moved);
                        centroids[ca] = recompute_cen(&clusters[ca]);

                        let sub_mems = clusters[wc].clone();
                        let sub_x: Vec<Vec<f64>> = sub_mems.iter().map(|&p| data.x[p].clone()).collect();
                        let c0 = init_kmeans_pp(&sub_x, 2, &mut rng);
                        let (lbl, _, _, _) = run_lloyd_from_centers(&sub_x, c0, 10);
                        let left: Vec<usize> = sub_mems.iter().enumerate().filter(|&(j, _)| lbl[j] == 0).map(|(_, &p)| p).collect();
                        let right: Vec<usize> = sub_mems.iter().enumerate().filter(|&(j, _)| lbl[j] != 0).map(|(_, &p)| p).collect();
                        if !left.is_empty() && !right.is_empty() {
                            clusters[wc] = left;
                            centroids[wc] = recompute_cen(&clusters[wc]);
                            clusters[cb] = right;
                            centroids[cb] = recompute_cen(&clusters[cb]);
                            grafts += 1;
                        }
                    }
                }
            }
        }
    }

    let m = clusters.len();
    let mut cdist = vec![vec![0.0; m]; m];
    for a in 0..m {
        for b in (a + 1)..m {
            let dd = euclidean(&centroids[a], &centroids[b]);
            cdist[a][b] = dd;
            cdist[b][a] = dd;
        }
    }
    let clbl = if m > k {
        agglomerative_lance_williams(&cdist, k, Linkage::Average)
    } else {
        (0..m as i32).collect()
    };
    let mut labels = vec![0i32; n];
    for (ci, mems) in clusters.iter().enumerate() {
        for &p in mems {
            labels[p] = clbl[ci];
        }
    }
    AlgorithmOutput::new(
        labels,
        format!("k={}, subnoder={}", k, m),
        format!("Inkrementelt GRINCH-tre med {} grafting-korreksjoner", grafts),
    )
}

// 78. BICO (Fichtenberger et al. 2013 — BIRCH + Coreset summary with weighted K-Means guarantee)
fn bico(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let coreset_size = (k * 8).clamp(16, 48).min(n);
    let mut threshold = median_pairwise_dist(&data.dist_matrix) * 0.25;

    // Build BICO coreset tree; double threshold whenever coreset exceeds budget
    let mut cfs: Vec<CfEntry> = Vec::new();
    for i in 0..n {
        let p = &data.x[i];
        let nearest = cfs.iter().enumerate().min_by(|(_, a), (_, b)| {
            euclidean_sq(p, &a.centroid()).total_cmp(&euclidean_sq(p, &b.centroid()))
        });
        if let Some((idx, entry)) = nearest {
            if euclidean(p, &entry.centroid()) <= threshold {
                cfs[idx].add(p, i);
                continue;
            }
        }
        cfs.push(CfEntry::new(p, i));
        // Rebuild/compress if over coreset budget
        if cfs.len() > coreset_size {
            threshold *= 1.5;
            let old = std::mem::take(&mut cfs);
            for mut entry in old {
                let cen = entry.centroid();
                if let Some((idx, _)) = cfs.iter().enumerate().min_by(|(_, a), (_, b)| {
                    euclidean_sq(&cen, &a.centroid()).total_cmp(&euclidean_sq(&cen, &b.centroid()))
                }) {
                    if euclidean(&cen, &cfs[idx].centroid()) <= threshold {
                        cfs[idx].n += entry.n;
                        for j in 0..d {
                            cfs[idx].ls[j] += entry.ls[j];
                        }
                        cfs[idx].ss += entry.ss;
                        cfs[idx].members.append(&mut entry.members);
                        continue;
                    }
                }
                cfs.push(entry);
            }
        }
    }

    // Weighted K-Means++ on coreset representatives
    let m = cfs.len();
    let cs_pts: Vec<Vec<f64>> = cfs.iter().map(|c| c.centroid()).collect();
    let cs_w: Vec<f64> = cfs.iter().map(|c| c.n as f64).collect();
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 78);
    let mut centers = init_kmeans_pp(&cs_pts, k.min(m), &mut rng);

    for _ in 0..30 {
        let mut sums = vec![vec![0.0; d]; centers.len()];
        let mut wsum = vec![0.0; centers.len()];
        for i in 0..m {
            let best_c = (0..centers.len())
                .min_by(|&a, &b| {
                    euclidean_sq(&cs_pts[i], &centers[a])
                        .total_cmp(&euclidean_sq(&cs_pts[i], &centers[b]))
                })
                .unwrap_or(0);
            wsum[best_c] += cs_w[i];
            for j in 0..d {
                sums[best_c][j] += cs_w[i] * cs_pts[i][j];
            }
        }
        for c in 0..centers.len() {
            if wsum[c] > 0.0 {
                for j in 0..d {
                    centers[c][j] = sums[c][j] / wsum[c];
                }
            }
        }
    }

    let labels: Vec<i32> = (0..n)
        .map(|i| {
            (0..centers.len())
                .min_by(|&a, &b| {
                    euclidean_sq(&data.x[i], &centers[a])
                        .total_cmp(&euclidean_sq(&data.x[i], &centers[b]))
                })
                .unwrap_or(0) as i32
        })
        .collect();
    AlgorithmOutput::new(
        labels,
        format!("k={}, coreset={}/{}", k, m, n),
        format!("Komprimerte til {}-punkts vektet coreset (T={:.2})", m, threshold),
    )
}

// 79. SCC (Monath et al. 2021 — Sub-Cluster Component via iterative reciprocal nearest-neighbor rounds)
fn scc_reciprocal_nn(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);

    let mut clusters: Vec<Vec<usize>> = (0..n).map(|i| vec![i]).collect();
    let mut rounds = 0usize;

    let centroid = |mems: &[usize]| -> Vec<f64> {
        let mut mu = vec![0.0; d];
        for &p in mems {
            for j in 0..d {
                mu[j] += data.x[p][j];
            }
        }
        for v in &mut mu {
            *v /= mems.len().max(1) as f64;
        }
        mu
    };

    while clusters.len() > k && rounds < 25 {
        rounds += 1;
        let m = clusters.len();
        let cens: Vec<Vec<f64>> = clusters.iter().map(|c| centroid(c)).collect();

        // Find 1-NN for each sub-cluster
        let mut nn = vec![0usize; m];
        for i in 0..m {
            nn[i] = (0..m)
                .filter(|&j| j != i)
                .min_by(|&a, &b| euclidean_sq(&cens[i], &cens[a]).total_cmp(&euclidean_sq(&cens[i], &cens[b])))
                .unwrap_or(i);
        }

        // Merge all reciprocal nearest-neighbor pairs (nn[i] == j && nn[j] == i)
        let mut merged_flag = vec![false; m];
        let mut next_clusters = Vec::new();
        let mut any_merged = false;
        for i in 0..m {
            if merged_flag[i] {
                continue;
            }
            let j = nn[i];
            if j != i && nn[j] == i && !merged_flag[j] && next_clusters.len() + (m - i) > k {
                let mut comb = clusters[i].clone();
                comb.extend_from_slice(&clusters[j]);
                next_clusters.push(comb);
                merged_flag[i] = true;
                merged_flag[j] = true;
                any_merged = true;
            } else {
                next_clusters.push(clusters[i].clone());
                merged_flag[i] = true;
            }
        }
        clusters = next_clusters;
        if !any_merged {
            break;
        }
    }

    // Final exact trim to k via average linkage if slightly above k
    if clusters.len() > k {
        let m = clusters.len();
        let cens: Vec<Vec<f64>> = clusters.iter().map(|c| centroid(c)).collect();
        let mut cdist = vec![vec![0.0; m]; m];
        for a in 0..m {
            for b in (a + 1)..m {
                let dd = euclidean(&cens[a], &cens[b]);
                cdist[a][b] = dd;
                cdist[b][a] = dd;
            }
        }
        let clbl = agglomerative_lance_williams(&cdist, k, Linkage::Average);
        let mut merged_map: HashMap<i32, Vec<usize>> = HashMap::new();
        for (ci, mems) in clusters.into_iter().enumerate() {
            merged_map.entry(clbl[ci]).or_default().extend(mems);
        }
        clusters = merged_map.into_values().collect();
    }

    let mut labels = vec![0i32; n];
    for (cid, mems) in clusters.iter().enumerate() {
        for &p in mems {
            labels[p] = cid as i32;
        }
    }
    AlgorithmOutput::new(
        labels,
        format!("k={}, RNN-runder={}", k, rounds),
        format!("Gjensidig nærmeste-nabo hierarki over {} runder", rounds),
    )
}

// 80. FINCH (Sarfraz et al. CVPR 2019 — Parameter-Free First Integer Neighbor Clustering Hierarchy)
fn finch(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let target_k = cfg.effective_k(n);

    let finch_step = |pts: &[Vec<f64>]| -> Vec<Vec<usize>> {
        let m = pts.len();
        if m <= 1 {
            return (0..m).map(|i| vec![i]).collect();
        }
        // First neighbor k1[i] for each point
        let mut k1 = vec![0usize; m];
        for i in 0..m {
            k1[i] = (0..m)
                .filter(|&j| j != i)
                .min_by(|&a, &b| euclidean_sq(&pts[i], &pts[a]).total_cmp(&euclidean_sq(&pts[i], &pts[b])))
                .unwrap_or(0);
        }
        // FINCH adjacency: connect i-j if j==k1[i] OR i==k1[j] OR k1[i]==k1[j]
        let mut adj = vec![Vec::new(); m];
        let mut by_first: HashMap<usize, Vec<usize>> = HashMap::new();
        for i in 0..m {
            let j = k1[i];
            adj[i].push(j);
            adj[j].push(i);
            by_first.entry(j).or_default().push(i);
        }
        for group in by_first.values() {
            for w in group.windows(2) {
                adj[w[0]].push(w[1]);
                adj[w[1]].push(w[0]);
            }
        }
        // Connected components
        let mut vis = vec![false; m];
        let mut comps = Vec::new();
        for i in 0..m {
            if vis[i] {
                continue;
            }
            let mut comp = Vec::new();
            let mut q = VecDeque::new();
            vis[i] = true;
            q.push_back(i);
            while let Some(u) = q.pop_front() {
                comp.push(u);
                for &v in &adj[u] {
                    if !vis[v] {
                        vis[v] = true;
                        q.push_back(v);
                    }
                }
            }
            comps.push(comp);
        }
        comps
    };

    // Generate full FINCH hierarchy of partitions and pick level closest to target_k (with >= 2 clusters)
    let mut partitions: Vec<Vec<i32>> = Vec::new();
    let mut cur_groups: Vec<Vec<usize>> = (0..n).map(|i| vec![i]).collect();

    for _level in 0..10 {
        let m = cur_groups.len();
        if m <= 2 {
            break;
        }
        let means: Vec<Vec<f64>> = cur_groups
            .iter()
            .map(|g| {
                let mut mu = vec![0.0; d];
                for &p in g {
                    for j in 0..d {
                        mu[j] += data.x[p][j];
                    }
                }
                for v in &mut mu {
                    *v /= g.len() as f64;
                }
                mu
            })
            .collect();
        let next_comps = finch_step(&means);
        if next_comps.len() >= m {
            break;
        }
        let mut next_groups = Vec::with_capacity(next_comps.len());
        let mut lbl = vec![0i32; n];
        for (cid, comp) in next_comps.into_iter().enumerate() {
            let mut pts = Vec::new();
            for g_idx in comp {
                pts.extend_from_slice(&cur_groups[g_idx]);
            }
            for &p in &pts {
                lbl[p] = cid as i32;
            }
            next_groups.push(pts);
        }
        let n_cl = next_groups.len();
        cur_groups = next_groups;
        if n_cl >= 2 {
            partitions.push(lbl);
        }
    }

    let chosen = partitions
        .iter()
        .min_by_key(|p| {
            let k_p = p.iter().copied().collect::<HashSet<_>>().len();
            (k_p as isize - target_k as isize).abs()
        })
        .cloned()
        .unwrap_or_else(|| vec![0i32; n]);

    let k_found = chosen.iter().copied().collect::<HashSet<_>>().len();
    AlgorithmOutput::new(
        chosen,
        format!("Parameterfri, {} nivåer", partitions.len()),
        format!("FINCH første-nabo hierarki valgte nivå med k={}", k_found),
    )
}
