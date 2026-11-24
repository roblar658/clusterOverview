use std::collections::{BinaryHeap, HashMap, HashSet, VecDeque};

use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
use rayon::prelude::*;

use crate::algorithms::cat1_centroid::{init_kmeans_pp, run_lloyd_from_centers};
use crate::algorithms::{AlgorithmOutput, RunConfig};
use crate::data::PreparedData;
use crate::math::{estimate_eps_kdist, euclidean, euclidean_sq, median_pairwise_dist};

pub fn run(id: usize, data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    match id {
        31 => dbscan(data, cfg),
        32 => hdbscan(data, cfg),
        33 => optics(data, cfg),
        34 => mean_shift(data, cfg),
        35 => denclue(data, cfg),
        36 => snn_clustering(data, cfg),
        37 => pdbscan(data, cfg),
        38 => dbclasd(data, cfg),
        39 => gdbscan(data, cfg),
        40 => st_dbscan(data, cfg),
        41 => quickshift(data, cfg),
        42 => medoidshift(data, cfg),
        43 => denclue_2(data, cfg),
        44 => k_dbscan(data, cfg),
        45 => rnn_dbscan(data, cfg),
        46 => debacl(data, cfg),
        47 => dbscan_plus_plus(data, cfg),
        48 => gridbscan(data, cfg),
        49 => mr_dbscan(data, cfg),
        50 => vdbscan(data, cfg),
        51 => density_peaks(data, cfg),
        52 => snn_dpc(data, cfg),
        53 => optics_of(data, cfg),
        54 => blurring_mean_shift(data, cfg),
        55 => mean_shift_plus_plus(data, cfg),
        _ => unreachable!(),
    }
}

pub fn run_dbscan_core(dist: &[Vec<f64>], eps: f64, min_pts: usize) -> Vec<i32> {
    let n = dist.len();
    let mut labels = vec![-2i32; n]; // -2 = unvisited, -1 = noise
    let mut cluster_id = 0i32;

    for i in 0..n {
        if labels[i] != -2 {
            continue;
        }
        let nbrs: Vec<usize> = (0..n).filter(|&j| dist[i][j] <= eps).collect();
        if nbrs.len() < min_pts {
            labels[i] = -1;
            continue;
        }
        labels[i] = cluster_id;
        let mut queue: VecDeque<usize> = nbrs.into_iter().collect();
        while let Some(q) = queue.pop_front() {
            if labels[q] == -1 {
                labels[q] = cluster_id;
            }
            if labels[q] != -2 {
                continue;
            }
            labels[q] = cluster_id;
            let q_nbrs: Vec<usize> = (0..n).filter(|&j| dist[q][j] <= eps).collect();
            if q_nbrs.len() >= min_pts {
                for nb in q_nbrs {
                    if labels[nb] == -2 || labels[nb] == -1 {
                        queue.push_back(nb);
                    }
                }
            }
        }
        cluster_id += 1;
    }
    labels
}

fn auto_tune_dbscan_eps(dist: &[Vec<f64>], min_pts: usize, target_k: usize) -> (f64, Vec<i32>) {
    let quantiles = [0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85];
    let mut best_eps = estimate_eps_kdist(dist, min_pts, 0.6);
    let mut best_labels = run_dbscan_core(dist, best_eps, min_pts);
    let mut best_score = f64::NEG_INFINITY;

    for &q in &quantiles {
        let eps = estimate_eps_kdist(dist, min_pts, q);
        let lbl = run_dbscan_core(dist, eps, min_pts);
        let n_clusters = lbl.iter().filter(|&&l| l >= 0).collect::<HashSet<_>>().len();
        let noise_frac = lbl.iter().filter(|&&l| l < 0).count() as f64 / dist.len().max(1) as f64;
        if n_clusters >= 2 && noise_frac < 0.75 {
            let sil = crate::metrics::silhouette_score(dist, &lbl).unwrap_or(-0.5);
            let k_pen = -0.05 * ((n_clusters as isize - target_k as isize).abs() as f64);
            let score = sil - 0.25 * noise_frac + k_pen;
            if score > best_score {
                best_score = score;
                best_eps = eps;
                best_labels = lbl;
            }
        }
    }
    (best_eps, best_labels)
}

// 31. DBSCAN
fn dbscan(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let min_pts = (data.n_features + 1).clamp(3, 8).min(data.n_samples / 3).max(2);
    let (eps, labels) = auto_tune_dbscan_eps(&data.dist_matrix, min_pts, cfg.effective_k(data.n_samples));
    let noise = labels.iter().filter(|&&l| l < 0).count();
    AlgorithmOutput::new(
        labels,
        format!("eps={:.3}, minPts={}", eps, min_pts),
        format!("Auto-kalibrert eps fra k-dist; {} støypunkter", noise),
    )
}

// 32. HDBSCAN (Mutual Reachability MST + Condensed Hierarchy + Stability Extraction)
fn hdbscan(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let min_samples = (data.n_features + 1).clamp(3, 6).min(n / 4).max(2);
    let min_cluster_size = (n / (cfg.effective_k(n) * 3)).clamp(4, 15).min(n / 2).max(3);

    // Core distances = distance to min_samples-th nearest neighbor
    let mut core_dist = vec![0.0; n];
    for i in 0..n {
        let mut row = data.dist_matrix[i].clone();
        row.sort_by(|a, b| a.total_cmp(b));
        core_dist[i] = row[min_samples.min(n - 1)];
    }

    // Prim's MST on mutual reachability distance
    let mut in_tree = vec![false; n];
    let mut min_edge = vec![f64::INFINITY; n];
    let mut parent = vec![0usize; n];
    let mut mst_edges: Vec<(usize, usize, f64)> = Vec::with_capacity(n.saturating_sub(1));
    min_edge[0] = 0.0;

    for _ in 0..n {
        let u = (0..n)
            .filter(|&i| !in_tree[i])
            .min_by(|&a, &b| min_edge[a].total_cmp(&min_edge[b]))
            .unwrap_or(0);
        in_tree[u] = true;
        if u != 0 && min_edge[u].is_finite() {
            mst_edges.push((parent[u], u, min_edge[u]));
        }
        for v in 0..n {
            if !in_tree[v] {
                let mrd = core_dist[u].max(core_dist[v]).max(data.dist_matrix[u][v]);
                if mrd < min_edge[v] {
                    min_edge[v] = mrd;
                    parent[v] = u;
                }
            }
        }
    }

    // Sort MST edges ascending -> build single-linkage cuts and pick cut maximizing stability / valid clusters
    mst_edges.sort_by(|a, b| a.2.total_cmp(&b.2));

    let mut best_labels = vec![-1i32; n];
    let mut best_score = f64::NEG_INFINITY;
    let mut best_lambda = 0.0;

    // Evaluate candidate cuts along the MST lambda spectrum
    let step = (mst_edges.len() / 24).max(1);
    for idx in (0..mst_edges.len()).step_by(step) {
        let thresh = mst_edges[idx].2;
        let mut adj = vec![Vec::new(); n];
        for &(u, v, w) in &mst_edges {
            if w <= thresh {
                adj[u].push(v);
                adj[v].push(u);
            }
        }
        let mut lbl = vec![-1i32; n];
        let mut cid = 0i32;
        let mut vis = vec![false; n];
        let mut stability = 0.0;

        for start in 0..n {
            if vis[start] {
                continue;
            }
            let mut comp = Vec::new();
            let mut q = VecDeque::new();
            vis[start] = true;
            q.push_back(start);
            while let Some(u) = q.pop_front() {
                comp.push(u);
                for &v in &adj[u] {
                    if !vis[v] {
                        vis[v] = true;
                        q.push_back(v);
                    }
                }
            }
            if comp.len() >= min_cluster_size {
                let lam_birth = 1.0 / thresh.max(1e-5);
                for &u in &comp {
                    lbl[u] = cid;
                    let lam_pt = 1.0 / core_dist[u].max(1e-5);
                    stability += (lam_pt - lam_birth).max(0.0);
                }
                cid += 1;
            }
        }
        if cid >= 2 {
            let sil = crate::metrics::silhouette_score(&data.dist_matrix, &lbl).unwrap_or(0.0);
            let noise_frac = lbl.iter().filter(|&&l| l < 0).count() as f64 / n as f64;
            let score = sil + 0.02 * (stability / n as f64).ln_1p() - 0.35 * noise_frac;
            if score > best_score {
                best_score = score;
                best_labels = lbl;
                best_lambda = 1.0 / thresh.max(1e-5);
            }
        }
    }

    AlgorithmOutput::new(
        best_labels,
        format!("min_cluster={}, min_samples={}", min_cluster_size, min_samples),
        format!("Gjensidig rekkevidde-MST, optimal lambda={:.2}", best_lambda),
    )
}

#[derive(Copy, Clone, PartialEq)]
struct OrdFloat(f64, usize);
impl Eq for OrdFloat {}
impl PartialOrd for OrdFloat {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for OrdFloat {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // Min-heap reversal
        other.0.total_cmp(&self.0).then_with(|| self.1.cmp(&other.1))
    }
}

fn optics_ordering(dist: &[Vec<f64>], min_pts: usize, max_eps: f64) -> (Vec<usize>, Vec<f64>, Vec<f64>) {
    let n = dist.len();
    let mut processed = vec![false; n];
    let mut reach = vec![f64::INFINITY; n];
    let mut core = vec![f64::INFINITY; n];
    let mut order = Vec::with_capacity(n);

    for i in 0..n {
        let mut d_nbrs: Vec<f64> = (0..n)
            .filter(|&j| dist[i][j] <= max_eps)
            .map(|j| dist[i][j])
            .collect();
        d_nbrs.sort_by(|a, b| a.total_cmp(b));
        if d_nbrs.len() >= min_pts {
            core[i] = d_nbrs[min_pts - 1];
        }
    }

    for start in 0..n {
        if processed[start] {
            continue;
        }
        let mut heap = BinaryHeap::new();
        heap.push(OrdFloat(0.0, start));
        while let Some(OrdFloat(_, p)) = heap.pop() {
            if processed[p] {
                continue;
            }
            processed[p] = true;
            order.push(p);
            if core[p].is_finite() {
                for q in 0..n {
                    if !processed[q] && dist[p][q] <= max_eps {
                        let new_reach = core[p].max(dist[p][q]);
                        if new_reach < reach[q] {
                            reach[q] = new_reach;
                            heap.push(OrdFloat(new_reach, q));
                        }
                    }
                }
            }
        }
    }
    (order, reach, core)
}

// 33. OPTICS
fn optics(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let min_pts = (data.n_features + 1).clamp(3, 7).min(n / 4).max(2);
    let max_eps = median_pairwise_dist(&data.dist_matrix) * 2.0;
    let (order, reach, core) = optics_ordering(&data.dist_matrix, min_pts, max_eps);

    // Extract clusters by valley threshold on reachability plot
    let mut finite_r: Vec<f64> = reach.iter().copied().filter(|r| r.is_finite()).collect();
    finite_r.sort_by(|a, b| a.total_cmp(b));
    let eps_cut = if finite_r.is_empty() {
        max_eps * 0.5
    } else {
        finite_r[(finite_r.len() * 72 / 100).min(finite_r.len() - 1)]
    };

    let mut labels = vec![-1i32; n];
    let mut cid = -1i32;
    for &p in &order {
        if reach[p] > eps_cut {
            if core[p] <= eps_cut {
                cid += 1;
                labels[p] = cid;
            } else {
                labels[p] = -1;
            }
        } else {
            if cid < 0 {
                cid = 0;
            }
            labels[p] = cid;
        }
    }

    // Fallback to clean valley cut at target k if single cluster emerged
    let n_found = labels.iter().filter(|&&l| l >= 0).collect::<HashSet<_>>().len();
    if n_found < 2 {
        let k = cfg.effective_k(n);
        let mut peaks: Vec<(usize, f64)> = order
            .iter()
            .enumerate()
            .skip(1)
            .map(|(pos, &p)| (pos, if reach[p].is_finite() { reach[p] } else { max_eps * 3.0 }))
            .collect();
        peaks.sort_by(|a, b| b.1.total_cmp(&a.1));
        let mut cuts: Vec<usize> = peaks.into_iter().take(k - 1).map(|(pos, _)| pos).collect();
        cuts.sort_unstable();
        let mut cur = 0i32;
        let mut ci = 0;
        for (pos, &p) in order.iter().enumerate() {
            if ci < cuts.len() && pos >= cuts[ci] {
                cur += 1;
                ci += 1;
            }
            labels[p] = cur;
        }
    }

    AlgorithmOutput::new(
        labels,
        format!("minPts={}, eps'={:.3}", min_pts, eps_cut),
        format!("Ordnet rekkeviddegraf over {} punkter", order.len()),
    )
}

fn merge_modes_to_labels(modes: &[Vec<f64>], merge_tol: f64) -> Vec<i32> {
    let mut centers: Vec<Vec<f64>> = Vec::new();
    let mut labels = vec![0i32; modes.len()];
    for (i, m) in modes.iter().enumerate() {
        let found = centers
            .iter()
            .position(|c| euclidean(m, c) <= merge_tol);
        if let Some(cid) = found {
            labels[i] = cid as i32;
        } else {
            labels[i] = centers.len() as i32;
            centers.push(m.clone());
        }
    }
    labels
}

// 34. Mean Shift (Gaussian kernel mode-seeking)
fn mean_shift(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let h = estimate_eps_kdist(&data.dist_matrix, (n / 6).clamp(4, 18), 0.5).max(1e-3);
    let h2 = 2.0 * h * h;
    let max_iter = cfg.max_iters(35);

    let mut shifted = data.x.clone();
    for i in 0..n {
        for _ in 0..max_iter {
            let mut num = vec![0.0; d];
            let mut den = 0.0;
            for p in &data.x {
                let w = (-euclidean_sq(&shifted[i], p) / h2).exp();
                den += w;
                for j in 0..d {
                    num[j] += w * p[j];
                }
            }
            if den <= 1e-12 {
                break;
            }
            let next: Vec<f64> = num.iter().map(|&v| v / den).collect();
            let step = euclidean(&shifted[i], &next);
            shifted[i] = next;
            if step < 1e-3 * h {
                break;
            }
        }
    }
    let labels = merge_modes_to_labels(&shifted, 0.5 * h);
    AlgorithmOutput::new(
        labels,
        format!("båndbredde h={:.3}", h),
        "Gaussisk kjerne-vindu mot lokale tetthetstopper".to_string(),
    )
}

// 35. DENCLUE (Hinneburg & Keim 1998 — analytical KDE gradient hill-climbing)
fn denclue(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let sigma = estimate_eps_kdist(&data.dist_matrix, (n / 7).clamp(4, 16), 0.5).max(1e-3);
    let s2 = 2.0 * sigma * sigma;
    let step_size = 0.25 * sigma;
    let max_iter = cfg.max_iters(40);

    let mut attractors = data.x.clone();
    let mut densities = vec![0.0; n];

    for i in 0..n {
        let mut x = data.x[i].clone();
        for _ in 0..max_iter {
            let mut grad = vec![0.0; d];
            let mut dens = 0.0;
            for p in &data.x {
                let k = (-euclidean_sq(&x, p) / s2).exp();
                dens += k;
                for j in 0..d {
                    grad[j] += k * (p[j] - x[j]);
                }
            }
            densities[i] = dens / n as f64;
            let gnorm: f64 = grad.iter().map(|v| v * v).sum::<f64>().sqrt();
            if gnorm < 1e-4 {
                break;
            }
            for j in 0..d {
                x[j] += step_size * grad[j] / gnorm;
            }
        }
        attractors[i] = x;
    }

    let mut labels = merge_modes_to_labels(&attractors, 0.6 * sigma);
    let mut sorted_d = densities.clone();
    sorted_d.sort_by(|a, b| a.total_cmp(b));
    let xi = sorted_d[n / 20]; // bottom 5% low-density noise threshold
    for i in 0..n {
        if densities[i] < xi {
            labels[i] = -1;
        }
    }
    AlgorithmOutput::new(
        labels,
        format!("sigma={:.3}, xi={:.4}", sigma, xi),
        "Analytisk tetthetsestimat med fast-steg bakkeklatring".to_string(),
    )
}

// 36. SNN (Shared Nearest Neighbor — Ertöz et al. 2003)
fn snn_clustering(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let k_nn = ((n as f64).sqrt().round() as usize).clamp(6, 18).min(n - 1);
    let mut knn_sets: Vec<HashSet<usize>> = Vec::with_capacity(n);
    for i in 0..n {
        let mut nbrs: Vec<(usize, f64)> = (0..n)
            .filter(|&j| j != i)
            .map(|j| (j, data.dist_matrix[i][j]))
            .collect();
        nbrs.sort_by(|a, b| a.1.total_cmp(&b.1));
        knn_sets.push(nbrs.into_iter().take(k_nn).map(|(j, _)| j).collect());
    }

    // SNN dissimilarity = 1 - |N(i) ∩ N(j)| / k_nn
    let mut snn_dist = vec![vec![1.0; n]; n];
    for i in 0..n {
        snn_dist[i][i] = 0.0;
        for j in (i + 1)..n {
            if knn_sets[i].contains(&j) && knn_sets[j].contains(&i) {
                let shared = knn_sets[i].intersection(&knn_sets[j]).count();
                let d_val = 1.0 - (shared as f64 / k_nn as f64);
                snn_dist[i][j] = d_val;
                snn_dist[j][i] = d_val;
            }
        }
    }

    let min_pts = (k_nn / 3).max(3);
    let (eps, labels) = auto_tune_dbscan_eps(&snn_dist, min_pts, cfg.effective_k(n));
    AlgorithmOutput::new(
        labels,
        format!("k_NN={}, SNN-eps={:.2}, minPts={}", k_nn, eps, min_pts),
        "Tetthetsclustering på delte nærmeste naboer (SNN-graf)".to_string(),
    )
}

// 37. PDBSCAN (Parallelized DBSCAN using Rayon parallel region queries & parallel component merge)
fn pdbscan(data: &PreparedData, _cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let min_pts = (data.n_features + 1).clamp(3, 8).min(n / 3).max(2);
    let eps = estimate_eps_kdist(&data.dist_matrix, min_pts, 0.62);

    // Parallel neighborhood computation across CPU cores
    let nbr_lists: Vec<Vec<usize>> = (0..n)
        .into_par_iter()
        .map(|i| {
            (0..n)
                .filter(|&j| data.dist_matrix[i][j] <= eps)
                .collect()
        })
        .collect();
    let is_core: Vec<bool> = nbr_lists.iter().map(|v| v.len() >= min_pts).collect();

    // Union-Find merge over core edges
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(p: &mut [usize], mut x: usize) -> usize {
        while p[x] != x {
            p[x] = p[p[x]];
            x = p[x];
        }
        x
    }
    for i in 0..n {
        if !is_core[i] {
            continue;
        }
        for &j in &nbr_lists[i] {
            if is_core[j] {
                let ri = find(&mut parent, i);
                let rj = find(&mut parent, j);
                if ri != rj {
                    parent[rj] = ri;
                }
            }
        }
    }

    let mut labels = vec![-1i32; n];
    for i in 0..n {
        if is_core[i] {
            labels[i] = find(&mut parent, i) as i32;
        }
    }
    // Assign border points to nearest core neighbor
    for i in 0..n {
        if !is_core[i] {
            if let Some(&core_nb) = nbr_lists[i]
                .iter()
                .filter(|&&j| is_core[j])
                .min_by(|&&a, &&b| data.dist_matrix[i][a].total_cmp(&data.dist_matrix[i][b]))
            {
                labels[i] = labels[core_nb];
            }
        }
    }
    let cores = is_core.iter().filter(|&&b| b).count();
    AlgorithmOutput::new(
        labels,
        format!("eps={:.3}, minPts={}, parallell", eps, min_pts),
        format!("Parallell Rayon-nabolag + Union-Find ({} kjernepunkter)", cores),
    )
}

// 38. DBCLASD (Xu et al. 1998 — nonparametric NN-distance distribution test)
fn dbclasd(data: &PreparedData, _cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let mut nn1 = vec![0.0; n];
    for i in 0..n {
        nn1[i] = (0..n)
            .filter(|&j| j != i)
            .map(|j| data.dist_matrix[i][j])
            .fold(f64::INFINITY, f64::min);
    }

    // Order candidates starting from densest points; grow cluster while NN-distance stays within chi-sq/Poisson expectation
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| nn1[a].total_cmp(&nn1[b]));

    let mut labels = vec![-1i32; n];
    let mut cid = 0i32;

    for &seed in &order {
        if labels[seed] >= 0 {
            continue;
        }
        // Seed initial cluster with seed + 3 nearest unassigned neighbors
        let mut init_nbrs: Vec<usize> = (0..n)
            .filter(|&j| j != seed && labels[j] < 0)
            .collect();
        init_nbrs.sort_by(|&a, &b| data.dist_matrix[seed][a].total_cmp(&data.dist_matrix[seed][b]));
        if init_nbrs.len() < 3 {
            continue;
        }
        let mut cluster = vec![seed, init_nbrs[0], init_nbrs[1], init_nbrs[2]];
        for &u in &cluster {
            labels[u] = cid;
        }

        // Iteratively expand with candidates whose NN distance fits cluster's expected NN distribution
        loop {
            let mut mean_nn = cluster.iter().map(|&u| nn1[u]).sum::<f64>() / cluster.len() as f64;
            let var_nn = cluster
                .iter()
                .map(|&u| (nn1[u] - mean_nn).powi(2))
                .sum::<f64>()
                / cluster.len().max(1) as f64;
            let std_nn = var_nn.sqrt().max(mean_nn * 0.25);
            mean_nn = mean_nn.max(1e-5);
            let upper_tol = mean_nn + 2.2 * std_nn;

            let mut added = false;
            let mut best_cand = None;
            let mut best_d = upper_tol;
            for cand in 0..n {
                if labels[cand] >= 0 {
                    continue;
                }
                let d_to_c = cluster
                    .iter()
                    .map(|&u| data.dist_matrix[cand][u])
                    .fold(f64::INFINITY, f64::min);
                if d_to_c <= best_d {
                    best_d = d_to_c;
                    best_cand = Some(cand);
                }
            }
            if let Some(c) = best_cand {
                labels[c] = cid;
                cluster.push(c);
                added = true;
            }
            if !added {
                break;
            }
        }
        cid += 1;
    }

    AlgorithmOutput::new(
        labels,
        "Parameterfri (NN-fordeling)".to_string(),
        format!("Vokste {} klynger via forventet 1-NN-avstandsfordeling", cid),
    )
}

// 39. GDBSCAN (Sander et al. 1998 — Generalized DBSCAN with weighted points & cosine+euclidean predicate)
fn gdbscan(data: &PreparedData, _cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let eps = estimate_eps_kdist(&data.dist_matrix, 5.min(n - 1), 0.62);
    // Point weight = local density inverse variance proxy (weighted cardinality)
    let weights: Vec<f64> = (0..n)
        .map(|i| {
            let cnt = (0..n).filter(|&j| data.dist_matrix[i][j] <= eps).count();
            1.0 + 0.1 * (cnt as f64).ln_1p()
        })
        .collect();
    let min_weight = 4.5_f64;

    let mut labels = vec![-2i32; n];
    let mut cid = 0i32;

    let n_pred = |i: usize, j: usize| -> bool {
        data.dist_matrix[i][j] <= eps && crate::math::cosine_sim(&data.x[i], &data.x[j]) > -0.2
    };

    for i in 0..n {
        if labels[i] != -2 {
            continue;
        }
        let nbrs: Vec<usize> = (0..n).filter(|&j| n_pred(i, j)).collect();
        let w_sum: f64 = nbrs.iter().map(|&j| weights[j]).sum();
        if w_sum < min_weight {
            labels[i] = -1;
            continue;
        }
        labels[i] = cid;
        let mut q: VecDeque<usize> = nbrs.into_iter().collect();
        while let Some(u) = q.pop_front() {
            if labels[u] == -1 {
                labels[u] = cid;
            }
            if labels[u] != -2 {
                continue;
            }
            labels[u] = cid;
            let unbrs: Vec<usize> = (0..n).filter(|&j| n_pred(u, j)).collect();
            let uw: f64 = unbrs.iter().map(|&j| weights[j]).sum();
            if uw >= min_weight {
                for nb in unbrs {
                    if labels[nb] == -2 || labels[nb] == -1 {
                        q.push_back(nb);
                    }
                }
            }
        }
        cid += 1;
    }

    AlgorithmOutput::new(
        labels,
        format!("eps={:.3}, minVekt={:.1}", eps, min_weight),
        "Generalisert nabolagspredikat (euklidsk + vinkel) med punktvekter".to_string(),
    )
}

// 40. ST-DBSCAN (Birant & Kut 2007 — Spatiotemporal DBSCAN with separate eps1 spatial and eps2 temporal)
fn st_dbscan(data: &PreparedData, _cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let eps1 = estimate_eps_kdist(&data.dist_matrix, 5.min(n - 1), 0.65);
    let eps2 = 0.45_f64; // normalized temporal window in [0, 1]
    let min_pts = 4usize.min(n / 3).max(2);

    let mut labels = vec![-2i32; n];
    let mut cid = 0i32;
    let st_nbrs = |i: usize| -> Vec<usize> {
        (0..n)
            .filter(|&j| {
                data.dist_matrix[i][j] <= eps1
                    && (data.time_values[i] - data.time_values[j]).abs() <= eps2
            })
            .collect()
    };

    for i in 0..n {
        if labels[i] != -2 {
            continue;
        }
        let nbrs = st_nbrs(i);
        if nbrs.len() < min_pts {
            labels[i] = -1;
            continue;
        }
        labels[i] = cid;
        let mut q: VecDeque<usize> = nbrs.into_iter().collect();
        while let Some(u) = q.pop_front() {
            if labels[u] == -1 {
                labels[u] = cid;
            }
            if labels[u] != -2 {
                continue;
            }
            labels[u] = cid;
            let un = st_nbrs(u);
            if un.len() >= min_pts {
                for nb in un {
                    if labels[nb] == -2 || labels[nb] == -1 {
                        q.push_back(nb);
                    }
                }
            }
        }
        cid += 1;
    }

    AlgorithmOutput::new(
        labels,
        format!("eps_rom={:.2}, eps_tid={:.2}, minPts={}", eps1, eps2, min_pts),
        "Rom-tid-tetthet med separate terskler for rom og tidsakse".to_string(),
    )
}

// 41. Quickshift (Vedaldi & Soatto 2008 — connect each point to nearest neighbor with strictly higher KDE)
fn quickshift(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let sigma = estimate_eps_kdist(&data.dist_matrix, (n / 8).clamp(4, 15), 0.5).max(1e-3);
    let s2 = 2.0 * sigma * sigma;
    let max_dist = sigma * 2.2;

    let density: Vec<f64> = (0..n)
        .map(|i| {
            (0..n)
                .map(|j| (-data.dist_matrix[i][j].powi(2) / s2).exp())
                .sum::<f64>()
        })
        .collect();

    let mut parent: Vec<usize> = (0..n).collect();
    let mut link_dist = vec![f64::INFINITY; n];
    for i in 0..n {
        let mut best_j = i;
        let mut best_d = f64::INFINITY;
        for j in 0..n {
            if density[j] > density[i] && data.dist_matrix[i][j] < best_d {
                best_d = data.dist_matrix[i][j];
                best_j = j;
            }
        }
        if best_d <= max_dist {
            parent[i] = best_j;
            link_dist[i] = best_d;
        }
    }

    // If only 1 root, cut the (k-1) longest parent links
    let mut roots: HashSet<usize> = (0..n)
        .map(|mut u| {
            while parent[u] != u {
                u = parent[u];
            }
            u
        })
        .collect();
    if roots.len() < 2 {
        let k = cfg.effective_k(n);
        let mut edges: Vec<(usize, f64)> = (0..n)
            .filter(|&i| parent[i] != i)
            .map(|i| (i, link_dist[i]))
            .collect();
        edges.sort_by(|a, b| b.1.total_cmp(&a.1));
        for &(i, _) in edges.iter().take(k - 1) {
            parent[i] = i;
        }
        roots = (0..n)
            .map(|mut u| {
                while parent[u] != u {
                    u = parent[u];
                }
                u
            })
            .collect();
    }

    let labels: Vec<i32> = (0..n)
        .map(|mut u| {
            while parent[u] != u {
                u = parent[u];
            }
            u as i32
        })
        .collect();
    AlgorithmOutput::new(
        labels,
        format!("sigma={:.3}, max_dist={:.2}", sigma, max_dist),
        format!("Koblet til nærmeste høyere tetthet ({} modustopper)", roots.len()),
    )
}

// 42. Medoidshift (Sheikh et al. 2007 — mode-seeking constrained to existing data points)
fn medoidshift(data: &PreparedData, _cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let sigma = estimate_eps_kdist(&data.dist_matrix, (n / 7).clamp(4, 15), 0.5).max(1e-3);
    let s2 = 2.0 * sigma * sigma;

    // Precompute Gaussian kernel weights
    let mut w = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            w[i][j] = (-data.dist_matrix[i][j].powi(2) / s2).exp();
        }
    }

    // For each point i, next medoid y(i) = argmin_k sum_j w(i,j) * d(k, j)^2
    let mut next_med = vec![0usize; n];
    for i in 0..n {
        let mut best_k = i;
        let mut best_val = f64::INFINITY;
        for k in 0..n {
            let mut val = 0.0;
            for j in 0..n {
                val += w[i][j] * data.dist_matrix[k][j].powi(2);
            }
            if val < best_val {
                best_val = val;
                best_k = k;
            }
        }
        next_med[i] = best_k;
    }

    // Follow shift pointers to fixed points (modes)
    let mut labels = vec![0i32; n];
    for i in 0..n {
        let mut cur = i;
        for _ in 0..n {
            let nxt = next_med[cur];
            if nxt == cur {
                break;
            }
            cur = nxt;
        }
        labels[i] = cur as i32;
    }
    AlgorithmOutput::new(
        labels,
        format!("sigma={:.3}, datapunkt-medoider", sigma),
        "Begrenset modusskift til eksisterende datapunkter".to_string(),
    )
}

// 43. DENCLUE 2.0 (Hinneburg & Gabriel 2007 — adaptive EM-like step-size hill climbing)
fn denclue_2(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let h = estimate_eps_kdist(&data.dist_matrix, (n / 6).clamp(4, 16), 0.55).max(1e-3);
    let h2 = 2.0 * h * h;
    let max_iter = cfg.max_iters(30);

    let mut attractors = data.x.clone();
    let mut total_steps = 0usize;

    for i in 0..n {
        let mut x = data.x[i].clone();
        for _ in 0..max_iter {
            total_steps += 1;
            let mut num = vec![0.0; d];
            let mut den = 0.0;
            for p in &data.x {
                let k = (-euclidean_sq(&x, p) / h2).exp();
                den += k;
                for j in 0..d {
                    num[j] += k * p[j];
                }
            }
            if den <= 1e-12 {
                break;
            }
            let next: Vec<f64> = num.iter().map(|&v| v / den).collect();
            let diff = euclidean(&x, &next);
            x = next;
            if diff < 1e-3 * h {
                break;
            }
        }
        attractors[i] = x;
    }

    let labels = merge_modes_to_labels(&attractors, 0.55 * h);
    AlgorithmOutput::new(
        labels,
        format!("h={:.3}, adaptiv EM-klatring", h),
        format!("Snitt {:.1} adaptive steg per punkt", total_steps as f64 / n as f64),
    )
}

// 44. K-DBSCAN (Coarse K-Means partitioning + local DBSCAN refinement inside partitions)
fn k_dbscan(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let k_coarse = cfg.effective_k(n).clamp(2, 6);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 44);
    let c0 = init_kmeans_pp(&data.x, k_coarse, &mut rng);
    let (coarse_lbl, _, _, _) = run_lloyd_from_centers(&data.x, c0, 25);

    let mut labels = vec![-1i32; n];
    let mut next_cid = 0i32;

    for part in 0..k_coarse {
        let members: Vec<usize> = (0..n).filter(|&i| coarse_lbl[i] == part as i32).collect();
        if members.len() < 4 {
            for &m in &members {
                labels[m] = next_cid;
            }
            next_cid += 1;
            continue;
        }
        let m_len = members.len();
        let mut sub_dist = vec![vec![0.0; m_len]; m_len];
        for a in 0..m_len {
            for b in (a + 1)..m_len {
                let dd = data.dist_matrix[members[a]][members[b]];
                sub_dist[a][b] = dd;
                sub_dist[b][a] = dd;
            }
        }
        let min_pts = 3usize.min(m_len / 2).max(2);
        let eps = estimate_eps_kdist(&sub_dist, min_pts, 0.72);
        let sub_lbl = run_dbscan_core(&sub_dist, eps, min_pts);
        let mut local_map: HashMap<i32, i32> = HashMap::new();
        for (idx, &sl) in sub_lbl.iter().enumerate() {
            if sl >= 0 {
                let global_id = *local_map.entry(sl).or_insert_with(|| {
                    let id = next_cid;
                    next_cid += 1;
                    id
                });
                labels[members[idx]] = global_id;
            }
        }
    }

    AlgorithmOutput::new(
        labels,
        format!("grov_k={}, lokal DBSCAN", k_coarse),
        format!("K-Means grovdeling ({}) + lokal tetthetsdeteksjon", k_coarse),
    )
}

// 45. RNN-DBSCAN (Bryant & Cios 2018 — Reverse Nearest Neighbors for varying densities)
fn rnn_dbscan(data: &PreparedData, _cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let k_nn = ((n as f64).sqrt().round() as usize).clamp(4, 12).min(n - 1);

    let mut knn: Vec<Vec<usize>> = Vec::with_capacity(n);
    let mut rnn: Vec<Vec<usize>> = vec![Vec::new(); n];
    for i in 0..n {
        let mut nbrs: Vec<(usize, f64)> = (0..n)
            .filter(|&j| j != i)
            .map(|j| (j, data.dist_matrix[i][j]))
            .collect();
        nbrs.sort_by(|a, b| a.1.total_cmp(&b.1));
        let top: Vec<usize> = nbrs.into_iter().take(k_nn).map(|(j, _)| j).collect();
        for &j in &top {
            rnn[j].push(i);
        }
        knn.push(top);
    }

    let is_core: Vec<bool> = (0..n).map(|i| rnn[i].len() >= k_nn).collect();
    let mut labels = vec![-1i32; n];
    let mut cid = 0i32;

    for i in 0..n {
        if labels[i] >= 0 || !is_core[i] {
            continue;
        }
        labels[i] = cid;
        let mut q = VecDeque::new();
        q.push_back(i);
        while let Some(u) = q.pop_front() {
            // Neighborhood = kNN(u) ∪ {y in RNN(u) | is_core(y)}
            for &v in &knn[u] {
                if labels[v] < 0 {
                    labels[v] = cid;
                    if is_core[v] {
                        q.push_back(v);
                    }
                }
            }
            for &v in &rnn[u] {
                if is_core[v] && labels[v] < 0 {
                    labels[v] = cid;
                    q.push_back(v);
                }
            }
        }
        cid += 1;
    }

    // Assign remaining border points to nearest core in their kNN
    for i in 0..n {
        if labels[i] < 0 {
            if let Some(&nb) = knn[i].iter().find(|&&j| is_core[j] && labels[j] >= 0) {
                labels[i] = labels[nb];
            }
        }
    }
    let cores = is_core.iter().filter(|&&b| b).count();
    AlgorithmOutput::new(
        labels,
        format!("k_RNN={}", k_nn),
        format!("Reverse-NN kjernepunkter: {}/{}", cores, n),
    )
}

// 46. DeBaCl (Kent et al. 2013 — Density-Based Clustering via Level Set Trees)
fn debacl(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let k_nn = ((n as f64).sqrt().round() as usize).clamp(4, 12).min(n - 1);
    let gamma = (n / (cfg.effective_k(n) * 4)).clamp(3, 12);

    // k-NN density estimate: f(x) ~ 1 / d_k(x)
    let mut dens = vec![0.0; n];
    let mut knn_adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for i in 0..n {
        let mut nbrs: Vec<(usize, f64)> = (0..n)
            .filter(|&j| j != i)
            .map(|j| (j, data.dist_matrix[i][j]))
            .collect();
        nbrs.sort_by(|a, b| a.1.total_cmp(&b.1));
        let dk = nbrs[k_nn - 1].1.max(1e-6);
        dens[i] = 1.0 / dk;
        for &(j, _) in nbrs.iter().take(k_nn) {
            knn_adj[i].push(j);
            knn_adj[j].push(i);
        }
    }

    // Scan upper-level sets {x: f(x) >= lambda} from high density down to find leaves of level-set tree
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| dens[b].total_cmp(&dens[a]));

    let mut best_labels = vec![-1i32; n];
    let mut best_score = f64::NEG_INFINITY;
    for frac in [0.45, 0.55, 0.65, 0.75, 0.85] {
        let active_n = ((n as f64) * frac).round() as usize;
        let active: HashSet<usize> = order.iter().take(active_n).copied().collect();
        let mut vis = HashSet::new();
        let mut lbl = vec![-1i32; n];
        let mut cid = 0i32;
        for &u in &active {
            if vis.contains(&u) {
                continue;
            }
            let mut comp = Vec::new();
            let mut q = VecDeque::new();
            vis.insert(u);
            q.push_back(u);
            while let Some(v) = q.pop_front() {
                comp.push(v);
                for &w in &knn_adj[v] {
                    if active.contains(&w) && !vis.contains(&w) {
                        vis.insert(w);
                        q.push_back(w);
                    }
                }
            }
            if comp.len() >= gamma {
                for x in comp {
                    lbl[x] = cid;
                }
                cid += 1;
            }
        }
        if cid >= 2 {
            let sil = crate::metrics::silhouette_score(&data.dist_matrix, &lbl).unwrap_or(-0.5);
            if sil > best_score {
                best_score = sil;
                best_labels = lbl;
            }
        }
    }

    // Assign remaining lower-density points to nearest assigned point (first-kNN assignment)
    for i in 0..n {
        if best_labels[i] < 0 {
            if let Some(nn) = (0..n)
                .filter(|&j| best_labels[j] >= 0)
                .min_by(|&a, &b| data.dist_matrix[i][a].total_cmp(&data.dist_matrix[i][b]))
            {
                best_labels[i] = best_labels[nn];
            }
        }
    }
    AlgorithmOutput::new(
        best_labels,
        format!("k_NN={}, pruning_gamma={}", k_nn, gamma),
        "Nivåsett-tre (upper level set) med k-NN tetthetsestimat".to_string(),
    )
}

// 47. DBSCAN++ (Jang & Jiang 2019 — compute core density only on greedy k-center subset m << n)
fn dbscan_plus_plus(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let m = ((n as f64) * 0.35).round().clamp(12.0, n as f64) as usize;
    let min_pts = (data.n_features + 1).clamp(3, 7).min(n / 4).max(2);
    let eps = estimate_eps_kdist(&data.dist_matrix, min_pts, 0.65);

    // Greedy k-center subset selection
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 47);
    let mut subset = vec![(rng.gen_range(0..n)) as usize];
    let mut min_d: Vec<f64> = (0..n).map(|i| data.dist_matrix[i][subset[0]]).collect();
    while subset.len() < m {
        let far = (0..n)
            .max_by(|&a, &b| min_d[a].total_cmp(&min_d[b]))
            .unwrap_or(0);
        subset.push(far);
        for i in 0..n {
            min_d[i] = min_d[i].min(data.dist_matrix[i][far]);
        }
    }

    // Only subset points are evaluated for core status
    let core_subset: Vec<usize> = subset
        .iter()
        .copied()
        .filter(|&s| (0..n).filter(|&j| data.dist_matrix[s][j] <= eps).count() >= min_pts)
        .collect();

    // Connected components of core_subset under 2*eps / eps_db
    let edge_eps = eps * 1.35;
    let mut core_lbl = vec![-1i32; core_subset.len()];
    let mut cid = 0i32;
    for i in 0..core_subset.len() {
        if core_lbl[i] >= 0 {
            continue;
        }
        core_lbl[i] = cid;
        let mut q = VecDeque::new();
        q.push_back(i);
        while let Some(u) = q.pop_front() {
            for v in 0..core_subset.len() {
                if core_lbl[v] < 0 && data.dist_matrix[core_subset[u]][core_subset[v]] <= edge_eps {
                    core_lbl[v] = cid;
                    q.push_back(v);
                }
            }
        }
        cid += 1;
    }

    // Assign every point to nearest core in core_subset within eps
    let mut labels = vec![-1i32; n];
    for i in 0..n {
        if let Some((ci, _)) = core_subset
            .iter()
            .enumerate()
            .min_by(|&(_, &a), &(_, &b)| data.dist_matrix[i][a].total_cmp(&data.dist_matrix[i][b]))
        {
            if data.dist_matrix[i][core_subset[ci]] <= edge_eps {
                labels[i] = core_lbl[ci];
            }
        }
    }
    use rand::Rng;
    AlgorithmOutput::new(
        labels,
        format!("utvalg m={}/{}, eps={:.3}", m, n, eps),
        format!("k-center kjerne-delutvalg fant {} kjernepunkter", core_subset.len()),
    )
}

// 48. GriDBSCAN (Mahran & Mahar 2008 — regular grid cell spatial index + local DBSCAN + boundary merge)
fn gridbscan(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let min_pts = (data.n_features + 1).clamp(3, 7).min(n / 4).max(2);
    let (eps, _) = auto_tune_dbscan_eps(&data.dist_matrix, min_pts, cfg.effective_k(n));

    // Partition in 2D PCA space into eps-sized cells for spatial indexing
    let cell_size = (eps * 1.2).max(1e-3);
    let mut grid: HashMap<(i32, i32), Vec<usize>> = HashMap::new();
    for i in 0..n {
        let cx = (data.pca_2d[i][0] / cell_size).floor() as i32;
        let cy = (data.pca_2d[i][1] / cell_size).floor() as i32;
        grid.entry((cx, cy)).or_default().push(i);
    }

    // Query only neighboring cells (3x3 halo) to form neighborhoods
    let mut nbrs = vec![Vec::new(); n];
    for (&(cx, cy), pts) in &grid {
        let mut cand = Vec::new();
        for dx in -1..=1 {
            for dy in -1..=1 {
                if let Some(v) = grid.get(&(cx + dx, cy + dy)) {
                    cand.extend_from_slice(v);
                }
            }
        }
        for &i in pts {
            for &j in &cand {
                if data.dist_matrix[i][j] <= eps {
                    nbrs[i].push(j);
                }
            }
        }
    }

    let mut labels = vec![-2i32; n];
    let mut cid = 0i32;
    for i in 0..n {
        if labels[i] != -2 {
            continue;
        }
        if nbrs[i].len() < min_pts {
            labels[i] = -1;
            continue;
        }
        labels[i] = cid;
        let mut q: VecDeque<usize> = nbrs[i].iter().copied().collect();
        while let Some(u) = q.pop_front() {
            if labels[u] == -1 {
                labels[u] = cid;
            }
            if labels[u] != -2 {
                continue;
            }
            labels[u] = cid;
            if nbrs[u].len() >= min_pts {
                for &nb in &nbrs[u] {
                    if labels[nb] == -2 || labels[nb] == -1 {
                        q.push_back(nb);
                    }
                }
            }
        }
        cid += 1;
    }

    AlgorithmOutput::new(
        labels,
        format!("celler={}, eps={:.3}", grid.len(), eps),
        "Rutenettpartisjonert 3x3 nabosøk + lokal DBSCAN".to_string(),
    )
}

// 49. MR-DBSCAN (He et al. 2011 — MapReduce spatial splits with eps-margin + Reduce boundary relabeling)
fn mr_dbscan(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let min_pts = (data.n_features + 1).clamp(3, 7).min(n / 4).max(2);
    let (eps, _) = auto_tune_dbscan_eps(&data.dist_matrix, min_pts, cfg.effective_k(n));

    // Split along PC1 into 4 overlapping Map partitions (with eps margin)
    let mut pc1: Vec<f64> = data.pca_2d.iter().map(|p| p[0]).collect();
    pc1.sort_by(|a, b| a.total_cmp(b));
    let splits = [pc1[n / 4], pc1[n / 2], pc1[(3 * n) / 4]];

    let n_parts = 4;
    let partitions: Vec<Vec<usize>> = (0..n_parts)
        .map(|p| {
            let lo = if p == 0 { f64::NEG_INFINITY } else { splits[p - 1] - eps };
            let hi = if p == n_parts - 1 { f64::INFINITY } else { splits[p] + eps };
            (0..n)
                .filter(|&i| data.pca_2d[i][0] >= lo && data.pca_2d[i][0] <= hi)
                .collect()
        })
        .collect();

    // Map phase: parallel local DBSCAN per partition
    let local_results: Vec<Vec<(usize, i32)>> = partitions
        .par_iter()
        .enumerate()
        .map(|(p_id, pts)| {
            let m = pts.len();
            let mut sub_d = vec![vec![0.0; m]; m];
            for a in 0..m {
                for b in (a + 1)..m {
                    let dd = data.dist_matrix[pts[a]][pts[b]];
                    sub_d[a][b] = dd;
                    sub_d[b][a] = dd;
                }
            }
            let lbl = run_dbscan_core(&sub_d, eps, min_pts);
            pts.iter()
                .enumerate()
                .map(|(i, &gid)| {
                    let cid = if lbl[i] >= 0 { (p_id as i32) * 1000 + lbl[i] } else { -1 };
                    (gid, cid)
                })
                .collect()
        })
        .collect();

    // Reduce phase: merge overlapping boundary cluster IDs via Union-Find
    let mut uf: HashMap<i32, i32> = HashMap::new();
    fn uf_find(uf: &mut HashMap<i32, i32>, x: i32) -> i32 {
        let mut cur = x;
        while let Some(&p) = uf.get(&cur) {
            if p == cur {
                break;
            }
            cur = p;
        }
        cur
    }
    let mut pt_assigned: Vec<Vec<i32>> = vec![Vec::new(); n];
    for part_res in &local_results {
        for &(gid, cid) in part_res {
            if cid >= 0 {
                uf.entry(cid).or_insert(cid);
                pt_assigned[gid].push(cid);
            }
        }
    }
    for cids in &pt_assigned {
        if cids.len() >= 2 {
            let r0 = uf_find(&mut uf, cids[0]);
            for &c in &cids[1..] {
                let rc = uf_find(&mut uf, c);
                if rc != r0 {
                    uf.insert(rc, r0);
                }
            }
        }
    }

    let mut labels = vec![-1i32; n];
    for i in 0..n {
        if let Some(&cid) = pt_assigned[i].first() {
            labels[i] = uf_find(&mut uf, cid);
        }
    }
    AlgorithmOutput::new(
        labels,
        format!("4 Map-partisjoner, eps={:.3}", eps),
        "MapReduce lokal DBSCAN + grense-sammenslåing i Reduce-fase".to_string(),
    )
}

// 50. VDBSCAN (Liu et al. 2007 — Varied-Density DBSCAN with multi-level eps from k-dist plot)
fn vdbscan(data: &PreparedData, _cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let min_pts = (data.n_features + 1).clamp(3, 6).min(n / 4).max(2);
    let eps_levels = [
        estimate_eps_kdist(&data.dist_matrix, min_pts, 0.35),
        estimate_eps_kdist(&data.dist_matrix, min_pts, 0.60),
        estimate_eps_kdist(&data.dist_matrix, min_pts, 0.85),
    ];

    let mut labels = vec![-1i32; n];
    let mut next_cid = 0i32;

    for &eps in &eps_levels {
        let unassigned: Vec<usize> = (0..n).filter(|&i| labels[i] < 0).collect();
        if unassigned.len() < min_pts {
            break;
        }
        let m = unassigned.len();
        let mut sub_d = vec![vec![0.0; m]; m];
        for a in 0..m {
            for b in (a + 1)..m {
                let dd = data.dist_matrix[unassigned[a]][unassigned[b]];
                sub_d[a][b] = dd;
                sub_d[b][a] = dd;
            }
        }
        let sub_lbl = run_dbscan_core(&sub_d, eps, min_pts);
        let mut map = HashMap::new();
        for (i, &sl) in sub_lbl.iter().enumerate() {
            if sl >= 0 {
                let gid = *map.entry(sl).or_insert_with(|| {
                    let id = next_cid;
                    next_cid += 1;
                    id
                });
                labels[unassigned[i]] = gid;
            }
        }
    }

    AlgorithmOutput::new(
        labels,
        format!("eps=[{:.2}, {:.2}, {:.2}]", eps_levels[0], eps_levels[1], eps_levels[2]),
        "Fler-tetthets DBSCAN i 3 k-dist-nivåer".to_string(),
    )
}

// 51. DPC (Density Peaks Clustering — Rodriguez & Laio, Science 2014)
fn density_peaks(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let k = cfg.effective_k(n);
    let dc = estimate_eps_kdist(&data.dist_matrix, (n / 12).clamp(3, 12), 0.5).max(1e-4);

    // Gaussian kernel local density rho_i
    let rho: Vec<f64> = (0..n)
        .map(|i| {
            (0..n)
                .filter(|&j| j != i)
                .map(|j| (-(data.dist_matrix[i][j] / dc).powi(2)).exp())
                .sum::<f64>()
        })
        .collect();

    // Order points by descending density
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| rho[b].total_cmp(&rho[a]));

    // delta_i = min distance to any point with higher density; nearest_higher[i]
    let mut delta = vec![0.0; n];
    let mut nearest_higher = vec![order[0]; n];
    let mut max_delta = 0.0_f64;

    for idx in 1..n {
        let i = order[idx];
        let mut best_d = f64::INFINITY;
        let mut best_j = order[0];
        for prev in 0..idx {
            let j = order[prev];
            if data.dist_matrix[i][j] < best_d {
                best_d = data.dist_matrix[i][j];
                best_j = j;
            }
        }
        delta[i] = best_d;
        nearest_higher[i] = best_j;
        if best_d > max_delta {
            max_delta = best_d;
        }
    }
    delta[order[0]] = max_delta.max(1e-4);

    // Decision graph score gamma_i = rho_i * delta_i -> pick top k centers
    let mut gamma: Vec<(usize, f64)> = (0..n).map(|i| (i, rho[i] * delta[i])).collect();
    gamma.sort_by(|a, b| b.1.total_cmp(&a.1));

    let mut labels = vec![-1i32; n];
    for (cid, &(center, _)) in gamma.iter().take(k).enumerate() {
        labels[center] = cid as i32;
    }
    // Single-pass assignment in descending density order
    for &i in &order {
        if labels[i] < 0 {
            labels[i] = labels[nearest_higher[i]];
        }
    }

    AlgorithmOutput::new(
        labels,
        format!("k={}, dc={:.3}", k, dc),
        "Beslutningsgraf (rho * delta) med enkelt-pass tildeling".to_string(),
    )
}

// 52. SNN-DPC (Liu et al. 2018 — Shared Nearest Neighbor Density Peaks)
fn snn_dpc(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let k = cfg.effective_k(n);
    let k_nn = ((n as f64).sqrt().round() as usize).clamp(5, 14).min(n - 1);

    let mut knn_sets: Vec<HashSet<usize>> = Vec::with_capacity(n);
    for i in 0..n {
        let mut nbrs: Vec<(usize, f64)> = (0..n)
            .filter(|&j| j != i)
            .map(|j| (j, data.dist_matrix[i][j]))
            .collect();
        nbrs.sort_by(|a, b| a.1.total_cmp(&b.1));
        knn_sets.push(nbrs.into_iter().take(k_nn).map(|(j, _)| j).collect());
    }

    // SNN similarity Sim(i,j) = |SNN|^2 / sum_{p in SNN} (d(i,p) + d(j,p))
    let mut sim = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            if knn_sets[i].contains(&j) && knn_sets[j].contains(&i) {
                let shared: Vec<usize> = knn_sets[i].intersection(&knn_sets[j]).copied().collect();
                if !shared.is_empty() {
                    let denom: f64 = shared
                        .iter()
                        .map(|&p| data.dist_matrix[i][p] + data.dist_matrix[j][p])
                        .sum::<f64>()
                        .max(1e-6);
                    let s = (shared.len() * shared.len()) as f64 / denom;
                    sim[i][j] = s;
                    sim[j][i] = s;
                }
            }
        }
    }

    // Local SNN density = sum of top k_nn similarities
    let rho: Vec<f64> = (0..n)
        .map(|i| {
            let mut vals = sim[i].clone();
            vals.sort_by(|a, b| b.total_cmp(a));
            vals.into_iter().take(k_nn).sum()
        })
        .collect();

    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| rho[b].total_cmp(&rho[a]));

    let mut delta = vec![0.0; n];
    let mut nearest_higher = vec![order[0]; n];
    let mut max_d = 0.0_f64;
    for idx in 1..n {
        let i = order[idx];
        let mut bd = f64::INFINITY;
        let mut bj = order[0];
        for prev in 0..idx {
            let j = order[prev];
            if data.dist_matrix[i][j] < bd {
                bd = data.dist_matrix[i][j];
                bj = j;
            }
        }
        delta[i] = bd;
        nearest_higher[i] = bj;
        max_d = max_d.max(bd);
    }
    delta[order[0]] = max_d.max(1e-4);

    let mut gamma: Vec<(usize, f64)> = (0..n).map(|i| (i, rho[i] * delta[i])).collect();
    gamma.sort_by(|a, b| b.1.total_cmp(&a.1));

    let mut labels = vec![-1i32; n];
    for (cid, &(cen, _)) in gamma.iter().take(k).enumerate() {
        labels[cen] = cid as i32;
    }
    for &i in &order {
        if labels[i] < 0 {
            labels[i] = labels[nearest_higher[i]];
        }
    }
    AlgorithmOutput::new(
        labels,
        format!("k={}, k_NN={}", k, k_nn),
        "SNN-basert lokal tetthet for robustere DPC-sentre".to_string(),
    )
}

// 53. OPTICS-OF (Breunig et al. 1999 — OPTICS + Local Outlier Factor anomaly flagging)
fn optics_of(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let min_pts = (data.n_features + 1).clamp(3, 7).min(n / 4).max(2);
    let max_eps = median_pairwise_dist(&data.dist_matrix) * 2.0;
    let (order, reach, core) = optics_ordering(&data.dist_matrix, min_pts, max_eps);

    // Compute Local Reachability Density (LRD) and OPTICS-OF (LOF) score
    let mut lrd = vec![0.0; n];
    let mut knn_idx: Vec<Vec<usize>> = Vec::with_capacity(n);
    for i in 0..n {
        let mut nbrs: Vec<(usize, f64)> = (0..n)
            .filter(|&j| j != i)
            .map(|j| (j, data.dist_matrix[i][j]))
            .collect();
        nbrs.sort_by(|a, b| a.1.total_cmp(&b.1));
        let top: Vec<usize> = nbrs.into_iter().take(min_pts).map(|(j, _)| j).collect();
        let sum_rd: f64 = top
            .iter()
            .map(|&j| core[j].min(max_eps).max(data.dist_matrix[i][j]))
            .sum::<f64>()
            .max(1e-8);
        lrd[i] = top.len() as f64 / sum_rd;
        knn_idx.push(top);
    }

    let of_scores: Vec<f64> = (0..n)
        .map(|i| {
            let avg_nbr_lrd =
                knn_idx[i].iter().map(|&j| lrd[j]).sum::<f64>() / knn_idx[i].len().max(1) as f64;
            avg_nbr_lrd / lrd[i].max(1e-8)
        })
        .collect();

    // Cluster via OPTICS ordering then flag top outlier-factor points (> 1.45) as anomalies (-1)
    let k = cfg.effective_k(n);
    let mut peaks: Vec<(usize, f64)> = order
        .iter()
        .enumerate()
        .skip(1)
        .map(|(pos, &p)| (pos, if reach[p].is_finite() { reach[p] } else { max_eps * 3.0 }))
        .collect();
    peaks.sort_by(|a, b| b.1.total_cmp(&a.1));
    let mut cuts: Vec<usize> = peaks.into_iter().take(k - 1).map(|(pos, _)| pos).collect();
    cuts.sort_unstable();

    let mut labels = vec![0i32; n];
    let mut cur = 0i32;
    let mut ci = 0;
    for (pos, &p) in order.iter().enumerate() {
        if ci < cuts.len() && pos >= cuts[ci] {
            cur += 1;
            ci += 1;
        }
        labels[p] = cur;
    }

    let mut flagged = 0usize;
    for i in 0..n {
        if of_scores[i] > 1.45 {
            labels[i] = -1;
            flagged += 1;
        }
    }
    AlgorithmOutput::new(
        labels,
        format!("minPts={}, OF-terskel=1.45", min_pts),
        format!("Detekterte {} lokale uteliggere (anomalier) via OPTICS-OF", flagged),
    )
}

// 54. Blurring Mean Shift (Carreira-Perpiñán 2006 — updates the entire dataset in-place each step)
fn blurring_mean_shift(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let h = estimate_eps_kdist(&data.dist_matrix, (n / 7).clamp(4, 14), 0.5).max(1e-3);
    let h2 = 2.0 * h * h;
    let max_iter = cfg.max_iters(15).min(15);

    let mut pts = data.x.clone();
    let mut iters = 0;
    for it in 0..max_iter {
        iters = it + 1;
        let mut next_pts = vec![vec![0.0; d]; n];
        let mut max_shift = 0.0_f64;
        for i in 0..n {
            let mut num = vec![0.0; d];
            let mut den = 0.0;
            for j in 0..n {
                let w = (-euclidean_sq(&pts[i], &pts[j]) / h2).exp();
                den += w;
                for m in 0..d {
                    num[m] += w * pts[j][m];
                }
            }
            for m in 0..d {
                next_pts[i][m] = num[m] / den.max(1e-12);
            }
            max_shift = max_shift.max(euclidean(&pts[i], &next_pts[i]));
        }
        pts = next_pts;
        // Stop before total collapse to single mode if clusters already separated sharply
        let unique_est = merge_modes_to_labels(&pts, 0.35 * h)
            .into_iter()
            .collect::<HashSet<_>>()
            .len();
        if max_shift < 1e-3 * h || unique_est <= cfg.effective_k(n) {
            break;
        }
    }

    let labels = merge_modes_to_labels(&pts, 0.35 * h);
    AlgorithmOutput::new(
        labels,
        format!("h={:.3}, in-place blur", h),
        format!("Datasett oppdatert synkront over {} blurring-iterasjoner", iters),
    )
}

// 55. Mean Shift++ (Jang & Jiang 2021 — linear-time grid-accelerated Mean Shift)
fn mean_shift_plus_plus(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let h = estimate_eps_kdist(&data.dist_matrix, (n / 6).clamp(4, 15), 0.52).max(1e-3);
    let max_iter = cfg.max_iters(25);
    // Project/hash on top up to 4 dimensions for grid binning (keeps 3^dim neighbor lookup O(1)!)
    let g_dim = d.min(4);
    let mut pts = data.x.clone();

    for _ in 0..max_iter {
        let mut cell_sum: HashMap<Vec<i32>, (Vec<f64>, usize)> = HashMap::new();
        for p in &pts {
            let key: Vec<i32> = (0..g_dim).map(|j| (p[j] / h).floor() as i32).collect();
            let entry = cell_sum.entry(key).or_insert_with(|| (vec![0.0; d], 0));
            for j in 0..d {
                entry.0[j] += p[j];
            }
            entry.1 += 1;
        }

        let mut max_shift = 0.0_f64;
        let mut next_pts = vec![vec![0.0; d]; n];
        for i in 0..n {
            let base: Vec<i32> = (0..g_dim).map(|j| (pts[i][j] / h).floor() as i32).collect();
            let mut sum = vec![0.0; d];
            let mut cnt = 0usize;
            // Enumerate 3^g_dim adjacent grid cells
            let total_offsets = 3usize.pow(g_dim as u32);
            let mut probe = vec![0i32; g_dim];
            for code in 0..total_offsets {
                let mut tmp = code;
                for j in 0..g_dim {
                    let off = (tmp % 3) as i32 - 1;
                    tmp /= 3;
                    probe[j] = base[j] + off;
                }
                if let Some((csum, ccnt)) = cell_sum.get(&probe) {
                    cnt += *ccnt;
                    for j in 0..d {
                        sum[j] += csum[j];
                    }
                }
            }
            if cnt > 0 {
                for j in 0..d {
                    next_pts[i][j] = sum[j] / cnt as f64;
                }
            } else {
                next_pts[i] = pts[i].clone();
            }
            max_shift = max_shift.max(euclidean(&pts[i], &next_pts[i]));
        }
        pts = next_pts;
        if max_shift < 1e-3 * h {
            break;
        }
    }

    let labels = merge_modes_to_labels(&pts, 0.6 * h);
    AlgorithmOutput::new(
        labels,
        format!("rutenett h={:.3}, d_grid={}", h, g_dim),
        "Lineær rutenettakselerert Mean Shift++".to_string(),
    )
}
