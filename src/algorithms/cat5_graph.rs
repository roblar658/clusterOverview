use std::collections::{HashMap, HashSet, VecDeque};

use rand::seq::SliceRandom;
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;

use crate::algorithms::cat1_centroid::{init_kmeans_pp, run_lloyd_from_centers};
use crate::algorithms::cat3_hierarchical::{agglomerative_lance_williams, Linkage};
use crate::algorithms::{AlgorithmOutput, RunConfig};
use crate::data::PreparedData;
use crate::math::{
    bottom_k_symmetric_eigen, median_pairwise_dist, normalize_vec, rbf_affinity_matrix,
    top_k_symmetric_eigen,
};

pub fn run(id: usize, data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    match id {
        103 => spectral_njw(data, cfg),
        104 => normalized_cuts(data, cfg),
        105 => ratio_cut(data, cfg),
        106 => affinity_propagation(data, cfg),
        107 => louvain(data, cfg),
        108 => leiden(data, cfg),
        109 => markov_clustering(data, cfg),
        110 => label_propagation(data, cfg),
        111 => chinese_whispers(data, cfg),
        112 => girvan_newman(data, cfg),
        113 => infomap(data, cfg),
        114 => walktrap(data, cfg),
        115 => fast_greedy_cnm(data, cfg),
        116 => leading_eigenvector(data, cfg),
        117 => spinglass(data, cfg),
        118 => hcs_mincut(data, cfg),
        119 => fluidc(data, cfg),
        120 => kernighan_lin(data, cfg),
        121 => metis_multilevel(data, cfg),
        122 => chaco_spectral_multilevel(data, cfg),
        123 => spectral_biclustering(data, cfg),
        124 => oslom(data, cfg),
        125 => bigclam(data, cfg),
        126 => slpa(data, cfg),
        127 => clique_percolation(data, cfg),
        128 => demon(data, cfg),
        129 => fast_unfolding_blondel(data, cfg),
        130 => sbm(data, cfg, false),
        131 => sbm(data, cfg, true),
        132 => power_iteration_clustering(data, cfg),
        _ => unreachable!(),
    }
}

fn dense_graph_from_knn(data: &PreparedData) -> (Vec<Vec<f64>>, Vec<f64>, f64) {
    let n = data.n_samples;
    let mut w = vec![vec![0.0; n]; n];
    let mut deg = vec![0.0; n];
    let mut two_m = 0.0;
    for i in 0..n {
        for &(j, wt) in &data.knn_graph[i] {
            w[i][j] = wt;
            deg[i] += wt;
            two_m += wt;
        }
    }
    (w, deg, two_m.max(1e-12))
}

fn graph_modularity(adj: &[Vec<(usize, f64)>], labels: &[i32]) -> f64 {
    let n = adj.len();
    let deg: Vec<f64> = adj.iter().map(|nbrs| nbrs.iter().map(|&(_, w)| w).sum()).collect();
    let two_m: f64 = deg.iter().sum::<f64>().max(1e-12);
    let mut q = 0.0;
    for i in 0..n {
        for &(j, w) in &adj[i] {
            if labels[i] == labels[j] {
                q += w;
            }
        }
    }
    let mut c_deg: HashMap<i32, f64> = HashMap::new();
    for i in 0..n {
        *c_deg.entry(labels[i]).or_insert(0.0) += deg[i];
    }
    let expected: f64 = c_deg.values().map(|&s| (s / two_m).powi(2)).sum();
    (q / two_m) - expected
}

// 103. Spectral Clustering (Ng-Jordan-Weiss 2001 — Symmetric Normalized Laplacian L_sym + row normalize + K-Means)
fn spectral_njw(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let k = cfg.effective_k(n);
    let sigma = median_pairwise_dist(&data.dist_matrix).max(1e-4);
    let aff = rbf_affinity_matrix(&data.dist_matrix, sigma);

    let deg: Vec<f64> = (0..n).map(|i| aff[i].iter().sum::<f64>().max(1e-12)).collect();
    // Normalized affinity M = D^-1/2 A D^-1/2; top-k eigenvectors = bottom-k of L_sym
    let mut norm_a = vec![vec![0.0; n]; n];
    for i in 0..n {
        let di = 1.0 / deg[i].sqrt();
        for j in 0..n {
            norm_a[i][j] = aff[i][j] * di / deg[j].sqrt();
        }
    }

    let (_, evecs) = top_k_symmetric_eigen(&norm_a, k, cfg.max_iters(35));
    let mut embed = vec![vec![0.0; k]; n];
    for i in 0..n {
        for c in 0..k {
            embed[i][c] = evecs[c][i];
        }
        embed[i] = normalize_vec(&embed[i]);
    }
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 103);
    let c0 = init_kmeans_pp(&embed, k, &mut rng);
    let (labels, _, _, _) = run_lloyd_from_centers(&embed, c0, 40);
    AlgorithmOutput::new(
        labels,
        format!("k={}, RBF sigma={:.2}", k, sigma),
        "Normalisert Laplace-matrise (NJW) + radnormalisering + K-Means".to_string(),
    )
}

// 104. Normalized Cuts (Shi & Malik 2000 — recursive 2-way generalized Fiedler eigenvector D^-1/2 (D-W) D^-1/2)
fn normalized_cuts(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let k = cfg.effective_k(n);
    let (w, _, _) = dense_graph_from_knn(data);

    let mut groups: Vec<Vec<usize>> = vec![(0..n).collect()];
    while groups.len() < k {
        let idx = groups
            .iter()
            .enumerate()
            .filter(|(_, g)| g.len() >= 4)
            .max_by_key(|(_, g)| g.len())
            .map(|(i, _)| i);
        let Some(gi) = idx else { break };
        let nodes = groups.remove(gi);
        let m = nodes.len();

        let mut lsym = vec![vec![0.0; m]; m];
        let mut deg = vec![0.0; m];
        for a in 0..m {
            for b in 0..m {
                if a != b {
                    deg[a] += w[nodes[a]][nodes[b]];
                }
            }
            deg[a] = deg[a].max(1e-6);
        }
        for a in 0..m {
            lsym[a][a] = 1.0;
            for b in 0..m {
                if a != b {
                    lsym[a][b] = -w[nodes[a]][nodes[b]] / (deg[a] * deg[b]).sqrt();
                }
            }
        }
        let (_, evecs) = bottom_k_symmetric_eigen(&lsym, 2.min(m), 30);
        let fiedler = evecs.get(1).unwrap_or(&evecs[0]);
        // Scale by D^-1/2 and sweep threshold minimizing Ncut(A,B)
        let mut vals: Vec<(usize, f64)> = (0..m)
            .map(|a| (a, fiedler[a] / deg[a].sqrt()))
            .collect();
        vals.sort_by(|a, b| a.1.total_cmp(&b.1));

        let mut best_cut = m / 2;
        let mut best_ncut = f64::INFINITY;
        let step = (m / 16).max(1);
        for split in (1..(m - 1)).step_by(step) {
            let mut in_left = vec![false; m];
            for &(a, _) in &vals[..split] {
                in_left[a] = true;
            }
            let mut cut_w = 0.0;
            let mut vol_a = 0.0;
            let mut vol_b = 0.0;
            for a in 0..m {
                if in_left[a] {
                    vol_a += deg[a];
                    for b in 0..m {
                        if !in_left[b] {
                            cut_w += w[nodes[a]][nodes[b]];
                        }
                    }
                } else {
                    vol_b += deg[a];
                }
            }
            let ncut = cut_w / vol_a.max(1e-8) + cut_w / vol_b.max(1e-8);
            if ncut < best_ncut {
                best_ncut = ncut;
                best_cut = split;
            }
        }
        let left: Vec<usize> = vals[..best_cut].iter().map(|&(a, _)| nodes[a]).collect();
        let right: Vec<usize> = vals[best_cut..].iter().map(|&(a, _)| nodes[a]).collect();
        groups.push(left);
        groups.push(right);
    }

    let mut labels = vec![0i32; n];
    for (cid, g) in groups.iter().enumerate() {
        for &p in g {
            labels[p] = cid as i32;
        }
    }
    AlgorithmOutput::new(
        labels,
        format!("k={}, rekursiv Shi-Malik Ncut", k),
        "Minimerte normalisert kutt via Fiedler-egenvektor-sveip".to_string(),
    )
}

// 105. Ratio Cut (Unnormalized Laplacian L = D - W -> bottom-k eigenvectors -> K-Means)
fn ratio_cut(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let k = cfg.effective_k(n);
    let (w, deg, _) = dense_graph_from_knn(data);

    let mut lap = vec![vec![0.0; n]; n];
    for i in 0..n {
        lap[i][i] = deg[i];
        for j in 0..n {
            if i != j {
                lap[i][j] = -w[i][j];
            }
        }
    }

    let (_, evecs) = bottom_k_symmetric_eigen(&lap, k, cfg.max_iters(35));
    let mut embed = vec![vec![0.0; k]; n];
    for i in 0..n {
        for c in 0..k {
            embed[i][c] = evecs[c][i];
        }
    }
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 105);
    let c0 = init_kmeans_pp(&embed, k, &mut rng);
    let (labels, _, _, _) = run_lloyd_from_centers(&embed, c0, 40);
    AlgorithmOutput::new(
        labels,
        format!("k={}, unormalisert L=D-W", k),
        "Spektral Ratio Cut-relaksering via bunn-k Laplace-egenvektorer".to_string(),
    )
}

// 106. Affinity Propagation (Frey & Dueck, Science 2007 — Responsibility & Availability message passing)
fn affinity_propagation(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let damping = 0.7_f64;
    let max_iter = cfg.max_iters(50);

    // Similarity s(i,k) = -||xi - xk||^2; preference s(k,k) = median similarity * 1.4
    let mut sim = vec![vec![0.0; n]; n];
    let mut vals = Vec::with_capacity(n * (n - 1) / 2);
    for i in 0..n {
        for j in (i + 1)..n {
            let s = -data.dist_matrix[i][j].powi(2);
            sim[i][j] = s;
            sim[j][i] = s;
            vals.push(s);
        }
    }
    vals.sort_by(|a, b| a.total_cmp(b));
    let pref = if vals.is_empty() { -1.0 } else { vals[vals.len() / 2] * 2.2 };
    for i in 0..n {
        sim[i][i] = pref;
    }

    let mut r = vec![vec![0.0; n]; n];
    let mut a = vec![vec![0.0; n]; n];

    for _ in 0..max_iter {
        // Update responsibilities R(i,k) = s(i,k) - max_{k' != k} (A(i,k') + s(i,k'))
        for i in 0..n {
            let mut max1 = f64::NEG_INFINITY;
            let mut max2 = f64::NEG_INFINITY;
            let mut arg1 = 0;
            for k in 0..n {
                let v = a[i][k] + sim[i][k];
                if v > max1 {
                    max2 = max1;
                    max1 = v;
                    arg1 = k;
                } else if v > max2 {
                    max2 = v;
                }
            }
            for k in 0..n {
                let mx = if k == arg1 { max2 } else { max1 };
                let r_new = sim[i][k] - mx;
                r[i][k] = damping * r[i][k] + (1.0 - damping) * r_new;
            }
        }
        // Update availabilities A(i,k)
        for k in 0..n {
            let mut pos_sum = 0.0;
            for i in 0..n {
                if i != k {
                    pos_sum += r[i][k].max(0.0);
                }
            }
            for i in 0..n {
                let a_new = if i == k {
                    pos_sum
                } else {
                    (r[k][k] + pos_sum - r[i][k].max(0.0)).min(0.0)
                };
                a[i][k] = damping * a[i][k] + (1.0 - damping) * a_new;
            }
        }
    }

    let exemplars: Vec<usize> = (0..n).filter(|&i| a[i][i] + r[i][i] > 0.0).collect();
    let labels: Vec<i32> = if !exemplars.is_empty() {
        (0..n)
            .map(|i| {
                exemplars
                    .iter()
                    .enumerate()
                    .max_by(|(_, &e1), (_, &e2)| sim[i][e1].total_cmp(&sim[i][e2]))
                    .map(|(cid, _)| cid as i32)
                    .unwrap_or(0)
            })
            .collect()
    } else {
        (0..n)
            .map(|i| {
                (0..n)
                    .max_by(|&c1, &c2| (a[i][c1] + r[i][c1]).total_cmp(&(a[i][c2] + r[i][c2])))
                    .unwrap_or(0) as i32
            })
            .collect()
    };

    AlgorithmOutput::new(
        labels,
        format!("demping={:.1}, pref={:.2}", damping, pref),
        format!("Meldingsutveksling valgte {} exemplars", exemplars.len().max(1)),
    )
}

fn louvain_core(adj: &[Vec<(usize, f64)>], resolution: f64, seed: u64) -> (Vec<i32>, f64) {
    let n = adj.len();
    let deg: Vec<f64> = adj.iter().map(|nbrs| nbrs.iter().map(|&(_, w)| w).sum()).collect();
    let two_m: f64 = deg.iter().sum::<f64>().max(1e-12);

    let mut comm: Vec<usize> = (0..n).collect();
    let mut tot_deg = deg.clone();
    let mut order: Vec<usize> = (0..n).collect();
    let mut rng = ChaCha8Rng::seed_from_u64(seed);

    for _pass in 0..20 {
        order.shuffle(&mut rng);
        let mut moved = false;
        for &i in &order {
            let cur_c = comm[i];
            let ki = deg[i];
            tot_deg[cur_c] -= ki;

            let mut comm_weights: HashMap<usize, f64> = HashMap::new();
            for &(j, w) in &adj[i] {
                if j != i {
                    *comm_weights.entry(comm[j]).or_insert(0.0) += w;
                }
            }
            let mut best_c = cur_c;
            let mut best_gain = comm_weights.get(&cur_c).copied().unwrap_or(0.0)
                - resolution * tot_deg[cur_c] * ki / two_m;

            for (&cand_c, &k_in) in &comm_weights {
                let gain = k_in - resolution * tot_deg[cand_c] * ki / two_m;
                if gain > best_gain + 1e-10 {
                    best_gain = gain;
                    best_c = cand_c;
                }
            }
            comm[i] = best_c;
            tot_deg[best_c] += ki;
            if best_c != cur_c {
                moved = true;
            }
        }
        if !moved {
            break;
        }
    }

    let lbl: Vec<i32> = comm.into_iter().map(|c| c as i32).collect();
    let q = graph_modularity(adj, &lbl);
    (lbl, q)
}

// 107. Louvain Method
fn louvain(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let (labels, q) = louvain_core(&data.knn_graph, 1.0, cfg.seed + 107);
    AlgorithmOutput::new(
        labels,
        "oppløsning=1.0, k-NN graf".to_string(),
        format!("Grådig modularitetsoptimering (Q={:.3})", q),
    )
}

// 108. Leiden Algorithm (Traag et al. 2019 — Louvain + connectivity refinement pass guaranteeing well-connected communities)
fn leiden(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let (mut labels, _) = louvain_core(&data.knn_graph, 1.05, cfg.seed + 108);

    // Refinement step: split any internally disconnected community into its connected sub-components
    let mut comm_nodes: HashMap<i32, Vec<usize>> = HashMap::new();
    for (i, &c) in labels.iter().enumerate() {
        comm_nodes.entry(c).or_default().push(i);
    }
    let mut refined = vec![-1i32; n];
    let mut next_cid = 0i32;
    let mut splits = 0usize;

    for nodes in comm_nodes.values() {
        let node_set: HashSet<usize> = nodes.iter().copied().collect();
        let mut vis = HashSet::new();
        let mut sub_count = 0usize;
        for &u in nodes {
            if vis.contains(&u) {
                continue;
            }
            sub_count += 1;
            let mut q = VecDeque::new();
            vis.insert(u);
            q.push_back(u);
            while let Some(cur) = q.pop_front() {
                refined[cur] = next_cid;
                for &(nb, _) in &data.knn_graph[cur] {
                    if node_set.contains(&nb) && !vis.contains(&nb) {
                        vis.insert(nb);
                        q.push_back(nb);
                    }
                }
            }
            next_cid += 1;
        }
        if sub_count > 1 {
            splits += sub_count - 1;
        }
    }
    labels = refined;
    let q = graph_modularity(&data.knn_graph, &labels);
    AlgorithmOutput::new(
        labels,
        "gamma=1.05 + sammenhengs-refinement".to_string(),
        format!("Leiden garantert sammenheng ({} splitt, Q={:.3})", splits, q),
    )
}

// 109. Markov Clustering (MCL — van Dongen 2000 — Expand M=M*M + Inflate M_ij^r)
fn markov_clustering(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let inflation = 2.0_f64;
    let max_iter = cfg.max_iters(25).min(25);

    let (mut m, _, _) = dense_graph_from_knn(data);
    // Add self-loops & column-normalize to stochastic transition matrix
    for i in 0..n {
        m[i][i] = 1.0;
    }
    let col_norm = |mat: &mut [Vec<f64>]| {
        let sz = mat.len();
        for j in 0..sz {
            let s: f64 = (0..sz).map(|i| mat[i][j]).sum::<f64>().max(1e-12);
            for i in 0..sz {
                mat[i][j] /= s;
            }
        }
    };
    col_norm(&mut m);

    for _ in 0..max_iter {
        // Expansion: M2 = M * M (skip near-zeros for speed!)
        let mut m2 = vec![vec![0.0; n]; n];
        for i in 0..n {
            for k in 0..n {
                let mik = m[i][k];
                if mik < 1e-4 {
                    continue;
                }
                for j in 0..n {
                    if m[k][j] > 1e-4 {
                        m2[i][j] += mik * m[k][j];
                    }
                }
            }
        }
        // Inflation: element-wise power r + prune + normalize
        let mut diff = 0.0_f64;
        for i in 0..n {
            for j in 0..n {
                let v = if m2[i][j] < 1e-4 { 0.0 } else { m2[i][j].powf(inflation) };
                m2[i][j] = v;
            }
        }
        col_norm(&mut m2);
        for i in 0..n {
            for j in 0..n {
                diff = diff.max((m2[i][j] - m[i][j]).abs());
            }
        }
        m = m2;
        if diff < 1e-3 {
            break;
        }
    }

    // Extract attractors: each column j belongs to row argmax_i M[i][j]
    let labels: Vec<i32> = (0..n)
        .map(|j| {
            (0..n)
                .max_by(|&a, &b| m[a][j].total_cmp(&m[b][j]))
                .unwrap_or(j) as i32
        })
        .collect();
    AlgorithmOutput::new(
        labels,
        format!("inflasjon r={:.1}, ekspansjon=2", inflation),
        "Markov-vandring ekspansjon/inflasjon konvergerte til attraktorer".to_string(),
    )
}

// 110. Label Propagation (Raghavan et al. 2007 — synchronous/asynchronous neighbor majority voting)
fn label_propagation(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let max_iter = cfg.max_iters(40);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 110);
    let mut labels: Vec<i32> = (0..n as i32).collect();
    let mut order: Vec<usize> = (0..n).collect();
    let mut iters = 0;

    for it in 0..max_iter {
        iters = it + 1;
        order.shuffle(&mut rng);
        let mut changed = false;
        for &i in &order {
            if data.knn_graph[i].is_empty() {
                continue;
            }
            let mut votes: HashMap<i32, f64> = HashMap::new();
            for &(j, w) in &data.knn_graph[i] {
                *votes.entry(labels[j]).or_insert(0.0) += w;
            }
            if let Some((&best_l, _)) = votes.iter().max_by(|a, b| a.1.total_cmp(b.1)) {
                if labels[i] != best_l {
                    labels[i] = best_l;
                    changed = true;
                }
            }
        }
        if !changed {
            break;
        }
    }
    AlgorithmOutput::new(
        labels,
        "k-NN graf, asynkron LPA".to_string(),
        format!("Konvergerte på {} iterasjoner (lineær tid)", iters),
    )
}

// 111. Chinese Whispers (Biemann 2006 — randomized linear-weighted graph clustering with top-class selection)
fn chinese_whispers(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let max_iter = cfg.max_iters(25);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 111);
    let mut labels: Vec<i32> = (0..n as i32).collect();
    let mut order: Vec<usize> = (0..n).collect();

    // Chinese Whispers uses inverse-degree/log-weighted voting in randomized passes
    for _ in 0..max_iter {
        order.shuffle(&mut rng);
        let mut changed = false;
        for &i in &order {
            let mut scores: HashMap<i32, f64> = HashMap::new();
            for &(j, w) in &data.knn_graph[i] {
                let deg_j = (data.knn_graph[j].len() as f64).ln_1p().max(1.0);
                *scores.entry(labels[j]).or_insert(0.0) += w / deg_j;
            }
            if let Some((&best_l, _)) = scores.iter().max_by(|a, b| a.1.total_cmp(b.1)) {
                if labels[i] != best_l {
                    labels[i] = best_l;
                    changed = true;
                }
            }
        }
        if !changed {
            break;
        }
    }
    AlgorithmOutput::new(
        labels,
        "log-gradvektet Chinese Whispers".to_string(),
        "Randomisert nabovalg over k-NN-nettverk".to_string(),
    )
}

// 112. Girvan-Newman (2002 — Brandes edge betweenness centrality divisive removal)
fn girvan_newman(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let k = cfg.effective_k(n);
    // Start with sparse symmetric 5-NN unweighted graph so edge removal splits cleanly & quickly
    let mut adj: Vec<HashSet<usize>> = vec![HashSet::new(); n];
    let deg_k = 5usize.min(n.saturating_sub(1).max(1));
    for i in 0..n {
        let mut nbrs: Vec<(usize, f64)> = (0..n)
            .filter(|&j| j != i)
            .map(|j| (j, data.dist_matrix[i][j]))
            .collect();
        nbrs.sort_by(|a, b| a.1.total_cmp(&b.1));
        for &(j, _) in nbrs.iter().take(deg_k) {
            adj[i].insert(j);
            adj[j].insert(i);
        }
    }

    let connected_components = |g: &[HashSet<usize>]| -> Vec<i32> {
        let mut lbl = vec![-1i32; g.len()];
        let mut cid = 0i32;
        for s in 0..g.len() {
            if lbl[s] >= 0 {
                continue;
            }
            let mut q = VecDeque::new();
            lbl[s] = cid;
            q.push_back(s);
            while let Some(u) = q.pop_front() {
                for &v in &g[u] {
                    if lbl[v] < 0 {
                        lbl[v] = cid;
                        q.push_back(v);
                    }
                }
            }
            cid += 1;
        }
        lbl
    };

    let mut removed = 0usize;
    let mut labels = connected_components(&adj);
    let max_steps = (n * 2).min(160);

    for _ in 0..max_steps {
        let n_comp = labels.iter().copied().collect::<HashSet<_>>().len();
        if n_comp >= k {
            break;
        }
        // Brandes edge betweenness (subsampled sources if n > 80 for snappy speed)
        let mut edge_betw: HashMap<(usize, usize), f64> = HashMap::new();
        let stride = (n / 45).max(1);
        for s in (0..n).step_by(stride) {
            let mut stack = Vec::new();
            let mut pred: Vec<Vec<usize>> = vec![Vec::new(); n];
            let mut sigma = vec![0.0_f64; n];
            let mut dist = vec![-1i32; n];
            sigma[s] = 1.0;
            dist[s] = 0;
            let mut q = VecDeque::new();
            q.push_back(s);
            while let Some(v) = q.pop_front() {
                stack.push(v);
                for &w in &adj[v] {
                    if dist[w] < 0 {
                        dist[w] = dist[v] + 1;
                        q.push_back(w);
                    }
                    if dist[w] == dist[v] + 1 {
                        sigma[w] += sigma[v];
                        pred[w].push(v);
                    }
                }
            }
            let mut delta = vec![0.0; n];
            while let Some(w) = stack.pop() {
                for &v in &pred[w] {
                    let c = (sigma[v] / sigma[w].max(1e-12)) * (1.0 + delta[w]);
                    let e = (v.min(w), v.max(w));
                    *edge_betw.entry(e).or_insert(0.0) += c;
                    delta[v] += c;
                }
            }
        }
        // Remove top highest-betweenness edge(s)
        let Some((&(u, v), _)) = edge_betw.iter().max_by(|a, b| a.1.total_cmp(b.1)) else {
            break;
        };
        adj[u].remove(&v);
        adj[v].remove(&u);
        removed += 1;
        labels = connected_components(&adj);
    }

    AlgorithmOutput::new(
        labels,
        format!("mål_k={}, 5-NN graf", k),
        format!("Fjernet {} kanter med høyest edge betweenness", removed),
    )
}

// 113. Infomap (Rosvall & Bergstrom 2008 — Greedy minimization of the Map Equation MDL L(M))
fn infomap(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let deg: Vec<f64> = data.knn_graph.iter().map(|v| v.iter().map(|&(_, w)| w).sum()).collect();
    let two_m: f64 = deg.iter().sum::<f64>().max(1e-12);
    let p_node: Vec<f64> = deg.iter().map(|&d| d / two_m).collect();

    let h_term = |p: f64| -> f64 {
        if p > 1e-12 {
            p * p.log2()
        } else {
            0.0
        }
    };

    let map_equation_l = |comm: &[usize]| -> f64 {
        let mut p_m: HashMap<usize, f64> = HashMap::new();
        let mut q_m: HashMap<usize, f64> = HashMap::new();
        for i in 0..n {
            *p_m.entry(comm[i]).or_insert(0.0) += p_node[i];
            for &(j, w) in &data.knn_graph[i] {
                if comm[j] != comm[i] {
                    *q_m.entry(comm[i]).or_insert(0.0) += w / two_m;
                }
            }
        }
        let q_tot: f64 = q_m.values().sum();
        let term_index = h_term(q_tot);
        let term_exit: f64 = q_m.values().map(|&q| h_term(q)).sum();
        let term_node: f64 = p_node.iter().map(|&p| h_term(p)).sum();
        let term_mod: f64 = p_m
            .iter()
            .map(|(m, &pm)| h_term(pm + q_m.get(m).copied().unwrap_or(0.0)))
            .sum();
        term_index - 2.0 * term_exit - term_node + term_mod
    };

    let mut comm: Vec<usize> = (0..n).collect();
    let mut cur_l = map_equation_l(&comm);
    let mut order: Vec<usize> = (0..n).collect();
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 113);

    for _ in 0..12 {
        order.shuffle(&mut rng);
        let mut moved = false;
        for &i in &order {
            let old_c = comm[i];
            let cand_comms: HashSet<usize> = data.knn_graph[i].iter().map(|&(j, _)| comm[j]).collect();
            let mut best_c = old_c;
            let mut best_l = cur_l;
            for c in cand_comms {
                if c == old_c {
                    continue;
                }
                comm[i] = c;
                let l_new = map_equation_l(&comm);
                if l_new + 1e-6 < best_l {
                    best_l = l_new;
                    best_c = c;
                }
            }
            comm[i] = best_c;
            if best_c != old_c {
                cur_l = best_l;
                moved = true;
            }
        }
        if !moved {
            break;
        }
    }

    AlgorithmOutput::new(
        comm.into_iter().map(|c| c as i32).collect(),
        "Map Equation MDL".to_string(),
        format!("Minimerte tilfeldig vandring-kodelengde til {:.3} bits", cur_l),
    )
}

// 114. Walktrap (Pons & Latapy 2005 — t-step random walk probability distance + Ward linkage)
fn walktrap(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let k = cfg.effective_k(n);
    let steps = 4usize;
    let (w, deg, _) = dense_graph_from_knn(data);

    // Transition matrix P_ij = W_ij / deg[i]
    let mut p = vec![vec![0.0; n]; n];
    for i in 0..n {
        let di = deg[i].max(1e-12);
        for j in 0..n {
            p[i][j] = w[i][j] / di;
        }
    }
    // Compute P^t by repeated squaring (t=4 -> 2 matrix squarings)
    let mat_mul = |a: &[Vec<f64>], b: &[Vec<f64>]| -> Vec<Vec<f64>> {
        let mut c = vec![vec![0.0; n]; n];
        for i in 0..n {
            for kk in 0..n {
                let aik = a[i][kk];
                if aik < 1e-6 {
                    continue;
                }
                for j in 0..n {
                    c[i][j] += aik * b[kk][j];
                }
            }
        }
        c
    };
    let p2 = mat_mul(&p, &p);
    let pt = mat_mul(&p2, &p2);

    // Walktrap distance r_ij = sqrt(sum_k (P^t_ik - P^t_jk)^2 / deg[k])
    let mut rdist = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let s: f64 = (0..n)
                .map(|m| (pt[i][m] - pt[j][m]).powi(2) / deg[m].max(1e-8))
                .sum();
            let d = s.sqrt();
            rdist[i][j] = d;
            rdist[j][i] = d;
        }
    }
    let labels = agglomerative_lance_williams(&rdist, k, Linkage::Ward);
    AlgorithmOutput::new(
        labels,
        format!("k={}, steg t={}", k, steps),
        "Ward-hierarki på 4-stegs tilfeldig-vandring-avstander".to_string(),
    )
}

// 115. Fast Greedy Modularity (Clauset-Newman-Moore 2004 — agglomerative delta-Q merging)
fn fast_greedy_cnm(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let target_k = cfg.effective_k(n);
    let (w, deg, two_m) = dense_graph_from_knn(data);

    let mut active = vec![true; n];
    let mut members: Vec<Vec<usize>> = (0..n).map(|i| vec![i]).collect();
    let mut e_mat = vec![vec![0.0; n]; n];
    let mut a_vec = vec![0.0; n];
    for i in 0..n {
        a_vec[i] = deg[i] / two_m;
        for j in 0..n {
            e_mat[i][j] = w[i][j] / two_m;
        }
    }

    let mut cur_q: f64 = (0..n).map(|i| e_mat[i][i] - a_vec[i] * a_vec[i]).sum();
    let mut best_q = cur_q;
    let mut best_members = members.clone();
    let mut best_active = active.clone();
    let mut remaining = n;

    while remaining > 2 {
        let mut best_pair = None;
        let mut best_dq = f64::NEG_INFINITY;
        for i in 0..n {
            if !active[i] {
                continue;
            }
            for j in (i + 1)..n {
                if !active[j] || e_mat[i][j] <= 0.0 {
                    continue;
                }
                let dq = 2.0 * (e_mat[i][j] - a_vec[i] * a_vec[j]);
                if dq > best_dq {
                    best_dq = dq;
                    best_pair = Some((i, j));
                }
            }
        }
        let Some((a, b)) = best_pair else { break };
        if best_dq <= 0.0 && remaining <= target_k {
            break;
        }
        cur_q += best_dq;
        e_mat[a][a] += e_mat[b][b] + 2.0 * e_mat[a][b];
        a_vec[a] += a_vec[b];
        for c in 0..n {
            if active[c] && c != a && c != b {
                e_mat[a][c] += e_mat[b][c];
                e_mat[c][a] = e_mat[a][c];
            }
        }
        let moved = std::mem::take(&mut members[b]);
        members[a].extend(moved);
        active[b] = false;
        remaining -= 1;

        if cur_q > best_q {
            best_q = cur_q;
            best_members = members.clone();
            best_active = active.clone();
        }
    }

    let mut labels = vec![0i32; n];
    let mut cid = 0i32;
    for i in 0..n {
        if best_active[i] {
            for &p in &best_members[i] {
                labels[p] = cid;
            }
            cid += 1;
        }
    }
    AlgorithmOutput::new(
        labels,
        "Clauset-Newman-Moore delta-Q".to_string(),
        format!("Agglomerativ modularitetstopp ved Q={:.3} (k={})", best_q, cid),
    )
}

// 116. Leading Eigenvector (Newman 2006 — recursive split by top eigenvector of Modularity Matrix B)
fn leading_eigenvector(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let k = cfg.effective_k(n);
    let (w, deg, two_m) = dense_graph_from_knn(data);

    let mut groups: Vec<Vec<usize>> = vec![(0..n).collect()];
    while groups.len() < k {
        let idx = groups
            .iter()
            .enumerate()
            .filter(|(_, g)| g.len() >= 4)
            .max_by_key(|(_, g)| g.len())
            .map(|(i, _)| i);
        let Some(gi) = idx else { break };
        let nodes = groups.remove(gi);
        let m = nodes.len();

        // Subgraph modularity matrix B^(g)_ij = B_ij - delta_ij * sum_{l in g} B_il
        let mut bmat = vec![vec![0.0; m]; m];
        for a in 0..m {
            let mut row_sum = 0.0;
            for b in 0..m {
                let val = w[nodes[a]][nodes[b]] - (deg[nodes[a]] * deg[nodes[b]]) / two_m;
                bmat[a][b] = val;
                row_sum += val;
            }
            bmat[a][a] -= row_sum;
        }

        let (evals, evecs) = top_k_symmetric_eigen(&bmat, 1, 30);
        if evals.first().copied().unwrap_or(0.0) <= 1e-6 {
            groups.push(nodes);
            break;
        }
        let u = &evecs[0];
        let left: Vec<usize> = (0..m).filter(|&a| u[a] >= 0.0).map(|a| nodes[a]).collect();
        let right: Vec<usize> = (0..m).filter(|&a| u[a] < 0.0).map(|a| nodes[a]).collect();
        if left.is_empty() || right.is_empty() {
            groups.push(nodes);
            break;
        }
        groups.push(left);
        groups.push(right);
    }

    let mut labels = vec![0i32; n];
    for (cid, g) in groups.iter().enumerate() {
        for &p in g {
            labels[p] = cid as i32;
        }
    }
    let q = graph_modularity(&data.knn_graph, &labels);
    AlgorithmOutput::new(
        labels,
        format!("mål_k={}, Newman B-matrise", k),
        format!("Spektral modularitetssplitt (k={}, Q={:.3})", groups.len(), q),
    )
}

// 117. Spinglass Community Detection (Reichardt & Bornholdt 2006 — Potts Hamiltonian simulated annealing)
fn spinglass(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let q_spins = cfg.effective_k(n);
    let gamma = 1.0_f64;
    let (_, deg, two_m) = dense_graph_from_knn(data);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 117);

    let mut spins: Vec<usize> = (0..n).map(|_| rng.gen_range(0..q_spins)).collect();
    let mut spin_deg = vec![0.0; q_spins];
    for i in 0..n {
        spin_deg[spins[i]] += deg[i];
    }

    let sweeps = cfg.max_iters(45);
    for step in 0..sweeps {
        let temp = 1.0 * (0.02_f64).powf(step as f64 / sweeps.max(1) as f64);
        for i in 0..n {
            let old_s = spins[i];
            let cand_s = rng.gen_range(0..q_spins);
            if cand_s == old_s {
                continue;
            }
            let ki = deg[i];
            let mut w_old = 0.0;
            let mut w_new = 0.0;
            for &(j, w) in &data.knn_graph[i] {
                if spins[j] == old_s {
                    w_old += w;
                } else if spins[j] == cand_s {
                    w_new += w;
                }
            }
            // Change in Potts Hamiltonian H = -sum (A_ij - gamma * p_ij) delta(s_i, s_j)
            let delta_h = -(w_new - w_old)
                + gamma * (ki / two_m) * (spin_deg[cand_s] - (spin_deg[old_s] - ki));
            if delta_h <= 0.0 || rng.gen::<f64>() < (-delta_h / temp.max(1e-4)).exp() {
                spin_deg[old_s] -= ki;
                spin_deg[cand_s] += ki;
                spins[i] = cand_s;
            }
        }
    }

    let labels: Vec<i32> = spins.into_iter().map(|s| s as i32).collect();
    let q = graph_modularity(&data.knn_graph, &labels);
    AlgorithmOutput::new(
        labels,
        format!("spinn q={}, gamma={:.1}", q_spins, gamma),
        format!("Potts-modell simulert avkjøling ({} sveip, Q={:.3})", sweeps, q),
    )
}

// 118. HCS (Hartuv & Shamir 2000 — Highly Connected Subgraphs via recursive Stoer-Wagner minimum cut)
fn hcs_mincut(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let target_k = cfg.effective_k(n);
    let (w, _, _) = dense_graph_from_knn(data);

    // Stoer-Wagner minimum cut on a subgraph
    let stoer_wagner = |nodes: &[usize]| -> (f64, Vec<usize>, Vec<usize>) {
        let m = nodes.len();
        let mut supernodes: Vec<Vec<usize>> = (0..m).map(|i| vec![nodes[i]]).collect();
        let mut active = vec![true; m];
        let mut gw = vec![vec![0.0; m]; m];
        for a in 0..m {
            for b in 0..m {
                if a != b {
                    gw[a][b] = w[nodes[a]][nodes[b]];
                }
            }
        }

        let mut best_cut = f64::INFINITY;
        let mut best_part = vec![nodes[0]];
        let mut rem = m;

        while rem > 1 {
            let active_list: Vec<usize> = (0..m).filter(|&i| active[i]).collect();
            let mut in_a = vec![false; m];
            let mut key = vec![0.0_f64; m];
            let mut prev = active_list[0];
            let mut last = active_list[0];

            for step in 0..active_list.len() {
                let u = if step == 0 {
                    active_list[0]
                } else {
                    *active_list
                        .iter()
                        .filter(|&&x| !in_a[x])
                        .max_by(|&&a, &&b| key[a].total_cmp(&key[b]))
                        .unwrap()
                };
                in_a[u] = true;
                prev = last;
                last = u;
                for &v in &active_list {
                    if !in_a[v] {
                        key[v] += gw[u][v];
                    }
                }
            }
            if key[last] < best_cut {
                best_cut = key[last];
                best_part = supernodes[last].clone();
            }
            // Contract last into prev
            let moved = std::mem::take(&mut supernodes[last]);
            supernodes[prev].extend(moved);
            for &v in &active_list {
                if v != prev && v != last {
                    gw[prev][v] += gw[last][v];
                    gw[v][prev] = gw[prev][v];
                }
            }
            active[last] = false;
            rem -= 1;
        }
        let part_set: HashSet<usize> = best_part.iter().copied().collect();
        let other: Vec<usize> = nodes.iter().copied().filter(|x| !part_set.contains(x)).collect();
        (best_cut, best_part, other)
    };

    let mut groups: Vec<Vec<usize>> = vec![(0..n).collect()];
    let mut cuts = 0usize;
    while groups.len() < target_k {
        let Some(gi) = groups
            .iter()
            .enumerate()
            .filter(|(_, g)| g.len() >= 6)
            .max_by_key(|(_, g)| g.len())
            .map(|(i, _)| i)
        else {
            break;
        };
        let nodes = groups.remove(gi);
        let (_, left, right) = stoer_wagner(&nodes);
        if left.is_empty() || right.is_empty() {
            groups.push(nodes);
            break;
        }
        cuts += 1;
        groups.push(left);
        groups.push(right);
    }

    let mut labels = vec![0i32; n];
    for (cid, g) in groups.iter().enumerate() {
        for &p in g {
            labels[p] = cid as i32;
        }
    }
    AlgorithmOutput::new(
        labels,
        format!("k={}, Stoer-Wagner min-cut", target_k),
        format!("Rekursiv minimum-kutt deling ({} kutt)", cuts),
    )
}

// 119. FLUIDC (Parés et al. 2018 — Fluid Communities with inverse-size fluid density propagation)
fn fluidc(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let k = cfg.effective_k(n);
    let max_iter = cfg.max_iters(40);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 119);

    // Seed k fluids on well-spaced points
    let mut seeds = vec![rng.gen_range(0..n)];
    while seeds.len() < k {
        let far = (0..n)
            .filter(|i| !seeds.contains(i))
            .max_by(|&a, &b| {
                let da = seeds.iter().map(|&s| data.dist_matrix[a][s]).fold(f64::INFINITY, f64::min);
                let db = seeds.iter().map(|&s| data.dist_matrix[b][s]).fold(f64::INFINITY, f64::min);
                da.total_cmp(&db)
            })
            .unwrap_or(0);
        seeds.push(far);
    }

    let mut comm = vec![-1i32; n];
    let mut counts = vec![0usize; k];
    for (c, &s) in seeds.iter().enumerate() {
        comm[s] = c as i32;
        counts[c] = 1;
    }

    let mut order: Vec<usize> = (0..n).collect();
    for _ in 0..max_iter {
        order.shuffle(&mut rng);
        let mut changed = false;
        for &u in &order {
            let mut fluid_score = vec![0.0; k];
            if comm[u] >= 0 {
                let c = comm[u] as usize;
                fluid_score[c] += 1.0 / counts[c].max(1) as f64;
            }
            for &(v, w) in &data.knn_graph[u] {
                if comm[v] >= 0 {
                    let c = comm[v] as usize;
                    fluid_score[c] += w / counts[c].max(1) as f64;
                }
            }
            let (best_c, best_s) = fluid_score
                .iter()
                .enumerate()
                .max_by(|a, b| a.1.total_cmp(b.1))
                .map(|(c, &s)| (c, s))
                .unwrap_or((0, 0.0));
            if best_s > 0.0 && comm[u] != best_c as i32 {
                if comm[u] >= 0 {
                    let old = comm[u] as usize;
                    if counts[old] <= 1 {
                        continue; // preserve every fluid
                    }
                    counts[old] -= 1;
                }
                comm[u] = best_c as i32;
                counts[best_c] += 1;
                changed = true;
            }
        }
        if !changed && comm.iter().all(|&c| c >= 0) {
            break;
        }
    }

    AlgorithmOutput::new(
        comm,
        format!("k={} væsker, tetthet=1/|C|", k),
        "Ekspanderende/sammentrekkende væskesimulering på graf".to_string(),
    )
}

fn kl_bisect(w: &[Vec<f64>], nodes: &[usize]) -> (Vec<usize>, Vec<usize>, f64) {
    let m = nodes.len();
    let half = m / 2;
    // Initialize balanced partition along PC1-like spectral or index split
    let mut part = vec![false; m];
    for i in half..m {
        part[i] = true;
    }

    for _pass in 0..8 {
        let mut locked = vec![false; m];
        let mut d_val = vec![0.0; m];
        for i in 0..m {
            let mut ext = 0.0;
            let mut intl = 0.0;
            for j in 0..m {
                if i == j {
                    continue;
                }
                let wt = w[nodes[i]][nodes[j]];
                if part[i] == part[j] {
                    intl += wt;
                } else {
                    ext += wt;
                }
            }
            d_val[i] = ext - intl;
        }

        let mut swaps = Vec::new();
        let mut gains = Vec::new();
        for _step in 0..half.min(m - half) {
            let mut best_a = None;
            let mut best_b = None;
            let mut best_g = f64::NEG_INFINITY;
            for a in 0..m {
                if locked[a] || part[a] {
                    continue;
                }
                for b in 0..m {
                    if locked[b] || !part[b] {
                        continue;
                    }
                    let g = d_val[a] + d_val[b] - 2.0 * w[nodes[a]][nodes[b]];
                    if g > best_g {
                        best_g = g;
                        best_a = Some(a);
                        best_b = Some(b);
                    }
                }
            }
            let (Some(a), Some(b)) = (best_a, best_b) else { break };
            locked[a] = true;
            locked[b] = true;
            swaps.push((a, b));
            gains.push(best_g);

            for x in 0..m {
                if locked[x] {
                    continue;
                }
                let sign = if !part[x] { 1.0 } else { -1.0 };
                d_val[x] += 2.0 * sign * (w[nodes[x]][nodes[a]] - w[nodes[x]][nodes[b]]);
            }
        }

        let mut best_pref = 0usize;
        let mut best_sum = 0.0;
        let mut acc = 0.0;
        for (i, &g) in gains.iter().enumerate() {
            acc += g;
            if acc > best_sum + 1e-9 {
                best_sum = acc;
                best_pref = i + 1;
            }
        }
        if best_pref == 0 {
            break;
        }
        for &(a, b) in &swaps[..best_pref] {
            part[a] = true;
            part[b] = false;
        }
    }

    let mut cut = 0.0;
    for i in 0..m {
        if !part[i] {
            for j in 0..m {
                if part[j] {
                    cut += w[nodes[i]][nodes[j]];
                }
            }
        }
    }
    let left: Vec<usize> = (0..m).filter(|&i| !part[i]).map(|i| nodes[i]).collect();
    let right: Vec<usize> = (0..m).filter(|&i| part[i]).map(|i| nodes[i]).collect();
    (left, right, cut)
}

// 120. Kernighan-Lin algoritme (1970 — recursive balanced bisection minimizing cut cost via pair swaps)
fn kernighan_lin(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let k = cfg.effective_k(n);
    let (w, _, _) = dense_graph_from_knn(data);

    let mut groups: Vec<Vec<usize>> = vec![(0..n).collect()];
    let mut total_cut = 0.0;
    while groups.len() < k {
        let Some(gi) = groups
            .iter()
            .enumerate()
            .filter(|(_, g)| g.len() >= 4)
            .max_by_key(|(_, g)| g.len())
            .map(|(i, _)| i)
        else {
            break;
        };
        let nodes = groups.remove(gi);
        let (l, r, c) = kl_bisect(&w, &nodes);
        total_cut += c;
        groups.push(l);
        groups.push(r);
    }

    let mut labels = vec![0i32; n];
    for (cid, g) in groups.iter().enumerate() {
        for &p in g {
            labels[p] = cid as i32;
        }
    }
    AlgorithmOutput::new(
        labels,
        format!("k={}, rekursiv KL-bytte", k),
        format!("Minimerte kuttkostnad via nodepar-bytting (kutt={:.2})", total_cut),
    )
}

// 121. METIS (Karypis & Kumar 1998 — Multilevel Heavy-Edge Matching Coarsening -> Partition -> FM Uncoarsening)
fn metis_multilevel(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let k = cfg.effective_k(n);
    let (mut w, _, _) = dense_graph_from_knn(data);

    // Phase 1: Multilevel Heavy-Edge Matching (HEM) coarsening stack
    let mut levels: Vec<Vec<usize>> = Vec::new(); // fine -> coarse mapping per level
    let mut cur_n = n;
    while cur_n > (k * 4).max(16) {
        let mut matched = vec![false; cur_n];
        let mut map = vec![0usize; cur_n];
        let mut next_id = 0usize;
        for u in 0..cur_n {
            if matched[u] {
                continue;
            }
            matched[u] = true;
            map[u] = next_id;
            if let Some(v) = (0..cur_n)
                .filter(|&x| !matched[x] && w[u][x] > 0.0)
                .max_by(|&a, &b| w[u][a].total_cmp(&w[u][b]))
            {
                matched[v] = true;
                map[v] = next_id;
            }
            next_id += 1;
        }
        if next_id >= cur_n {
            break;
        }
        let mut w_next = vec![vec![0.0; next_id]; next_id];
        for i in 0..cur_n {
            for j in 0..cur_n {
                if map[i] != map[j] {
                    w_next[map[i]][map[j]] += w[i][j];
                }
            }
        }
        levels.push(map);
        w = w_next;
        cur_n = next_id;
    }

    // Phase 2: Initial k-way partition on coarsest graph via Kernighan-Lin bisections
    let mut groups: Vec<Vec<usize>> = vec![(0..cur_n).collect()];
    while groups.len() < k {
        let Some(gi) = groups
            .iter()
            .enumerate()
            .filter(|(_, g)| g.len() >= 2)
            .max_by_key(|(_, g)| g.len())
            .map(|(i, _)| i)
        else {
            break;
        };
        let nodes = groups.remove(gi);
        let (l, r, _) = kl_bisect(&w, &nodes);
        groups.push(l);
        groups.push(r);
    }
    let mut part = vec![0usize; cur_n];
    for (cid, g) in groups.iter().enumerate() {
        for &u in g {
            part[u] = cid;
        }
    }

    // Phase 3: Uncoarsen back to original graph + boundary Fiduccia-Mattheyses refinement
    for map in levels.into_iter().rev() {
        let fine_n = map.len();
        part = (0..fine_n).map(|u| part[map[u]]).collect();
    }

    // Final boundary FM refinement on original graph
    let (w_orig, _, _) = dense_graph_from_knn(data);
    let max_size = ((n as f64 / k as f64) * 1.25).ceil() as usize;
    let mut sizes = vec![0usize; k];
    for &c in &part {
        sizes[c] += 1;
    }
    for _ in 0..6 {
        let mut improved = false;
        for u in 0..n {
            let cur = part[u];
            if sizes[cur] <= 2 {
                continue;
            }
            let cur_in: f64 = data.knn_graph[u]
                .iter()
                .filter(|&&(v, _)| part[v] == cur)
                .map(|&(_, wt)| wt)
                .sum();
            let mut best_c = cur;
            let mut best_gain = 0.0;
            for c in 0..k {
                if c == cur || sizes[c] >= max_size {
                    continue;
                }
                let cand_in: f64 = data.knn_graph[u]
                    .iter()
                    .filter(|&&(v, _)| part[v] == c)
                    .map(|&(_, wt)| wt)
                    .sum();
                if cand_in - cur_in > best_gain + 1e-8 {
                    best_gain = cand_in - cur_in;
                    best_c = c;
                }
            }
            if best_c != cur {
                sizes[cur] -= 1;
                sizes[best_c] += 1;
                part[u] = best_c;
                improved = true;
            }
        }
        if !improved {
            break;
        }
    }
    let _ = w_orig;

    AlgorithmOutput::new(
        part.into_iter().map(|c| c as i32).collect(),
        format!("k={}, HEM flernivå + FM-refinement", k),
        "Flernivå coarsening -> partisjonering -> FM uncoarsening".to_string(),
    )
}

// 122. Chaco (Hendrickson & Leland 1995 — Spectral Quadrisection + Multilevel refinement)
fn chaco_spectral_multilevel(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let k = cfg.effective_k(n);
    let (w, deg, _) = dense_graph_from_knn(data);

    // Spectral coordinates using 2nd and 3rd smallest Laplacian eigenvectors (Chaco quadrisection)
    let mut lap = vec![vec![0.0; n]; n];
    for i in 0..n {
        lap[i][i] = deg[i];
        for j in 0..n {
            if i != j {
                lap[i][j] = -w[i][j];
            }
        }
    }
    let (_, evecs) = bottom_k_symmetric_eigen(&lap, 3.min(n), 30);
    let dim_e = evecs.len().saturating_sub(1).max(1);
    let mut coords = vec![vec![0.0; dim_e]; n];
    for i in 0..n {
        for c in 0..dim_e {
            coords[i][c] = evecs[(c + 1).min(evecs.len() - 1)][i];
        }
    }
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 122);
    let c0 = init_kmeans_pp(&coords, k, &mut rng);
    let (mut labels, _, _, _) = run_lloyd_from_centers(&coords, c0, 25);

    // Local Kernighan-Lin style boundary polish
    for _ in 0..5 {
        let mut changed = false;
        for u in 0..n {
            let mut score = vec![0.0; k];
            for &(v, wt) in &data.knn_graph[u] {
                score[labels[v] as usize] += wt;
            }
            let best_c = (0..k).max_by(|&a, &b| score[a].total_cmp(&score[b])).unwrap_or(0) as i32;
            if best_c != labels[u] && score[best_c as usize] > score[labels[u] as usize] * 1.1 {
                labels[u] = best_c;
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }
    AlgorithmOutput::new(
        labels,
        format!("k={}, spektral + lokal polering", k),
        "Chaco spektral partisjonering med grense-refinement".to_string(),
    )
}

// 123. Spectral Biclustering (Kluger et al. 2003 — Bistochastic normalization + SVD checkerboard)
fn spectral_biclustering(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);

    // Normalize non-negative matrix A_n = R^-1/2 A C^-1/2
    let r_sum: Vec<f64> = (0..n).map(|i| data.x_nonneg[i].iter().sum::<f64>().max(1e-8)).collect();
    let c_sum: Vec<f64> = (0..d)
        .map(|j| (0..n).map(|i| data.x_nonneg[i][j]).sum::<f64>().max(1e-8))
        .collect();
    let mut an = vec![vec![0.0; d]; n];
    for i in 0..n {
        for j in 0..d {
            an[i][j] = data.x_nonneg[i][j] / (r_sum[i] * c_sum[j]).sqrt();
        }
    }

    // Left singular vectors = eigenvectors of A_n * A_n^T
    let mut aat = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in i..n {
            let val: f64 = (0..d).map(|m| an[i][m] * an[j][m]).sum();
            aat[i][j] = val;
            aat[j][i] = val;
        }
    }
    let m_sv = (k.min(d)).max(2);
    let (_, evecs) = top_k_symmetric_eigen(&aat, m_sv, 30);
    let mut embed = vec![vec![0.0; m_sv]; n];
    for i in 0..n {
        for c in 0..m_sv {
            embed[i][c] = evecs[c][i] / r_sum[i].sqrt();
        }
    }

    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 123);
    let c0 = init_kmeans_pp(&embed, k, &mut rng);
    let (labels, _, _, _) = run_lloyd_from_centers(&embed, c0, 35);
    AlgorithmOutput::new(
        labels,
        format!("k={}, {} singulærvektorer", k, m_sv),
        "Kluger spektral sjakkbrett-biklynging på R^-1/2 A C^-1/2".to_string(),
    )
}

// 124. OSLOM (Lancichinetti et al. 2011 — Statistically significant local fitness expansion)
fn oslom(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let deg: Vec<f64> = data.knn_graph.iter().map(|v| v.iter().map(|&(_, w)| w).sum()).collect();
    let two_m: f64 = deg.iter().sum::<f64>().max(1e-12);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 124);

    let mut visited = vec![false; n];
    let mut memberships: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    let mut cid = 0usize;
    let mut order: Vec<usize> = (0..n).collect();
    order.shuffle(&mut rng);

    for &seed in &order {
        if visited[seed] {
            continue;
        }
        let mut cluster: HashSet<usize> = HashSet::new();
        cluster.insert(seed);
        let mut c_vol = deg[seed];

        // Expand cluster by adding neighbors with statistically significant internal degree vs null model
        for _ in 0..(n / 2).max(5) {
            let mut best_cand = None;
            let mut best_z = 0.8_f64; // significance threshold
            let mut boundary: HashSet<usize> = HashSet::new();
            for &u in &cluster {
                for &(v, _) in &data.knn_graph[u] {
                    if !cluster.contains(&v) {
                        boundary.insert(v);
                    }
                }
            }
            for v in boundary {
                let k_in: f64 = data.knn_graph[v]
                    .iter()
                    .filter(|&&(nb, _)| cluster.contains(&nb))
                    .map(|&(_, w)| w)
                    .sum();
                let exp_in = deg[v] * c_vol / two_m;
                let z = (k_in - exp_in) / (exp_in.sqrt().max(0.25));
                if z > best_z {
                    best_z = z;
                    best_cand = Some((v, z));
                }
            }
            if let Some((v, z)) = best_cand {
                cluster.insert(v);
                c_vol += deg[v];
                memberships[v].push((cid, z));
            } else {
                break;
            }
        }
        for &u in &cluster {
            visited[u] = true;
            if memberships[u].iter().all(|&(c, _)| c != cid) {
                memberships[u].push((cid, 1.0));
            }
        }
        cid += 1;
    }

    let overlap_pts = memberships.iter().filter(|m| m.len() > 1).count();
    let labels: Vec<i32> = memberships
        .iter()
        .map(|m| {
            m.iter()
                .max_by(|a, b| a.1.total_cmp(&b.1))
                .map(|&(c, _)| c as i32)
                .unwrap_or(-1)
        })
        .collect();
    AlgorithmOutput::new(
        labels,
        format!("z-terskel=0.8, {} signifikante klynger", cid),
        format!("Fant {} overlappende noder (tilordnet sterkeste klynge)", overlap_pts),
    )
}

// 125. BigCLAM (Yang & Leskovec WSDM 2013 — Nonnegative Matrix Factorization community affiliation F >= 0)
fn bigclam(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let k = cfg.effective_k(n);
    let max_iter = cfg.max_iters(40);
    let lr = 0.03_f64;
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 125);

    // Initialize F_uc >= 0
    let mut f = vec![vec![0.0; k]; n];
    for u in 0..n {
        for c in 0..k {
            f[u][c] = rng.gen_range(0.05..0.5);
        }
    }

    // Maximize log-likelihood sum_{(u,v) in E} log(1 - exp(-F_u . F_v)) - sum_{(u,v) not in E} F_u . F_v
    for _ in 0..max_iter {
        let mut f_sum = vec![0.0; k];
        for u in 0..n {
            for c in 0..k {
                f_sum[c] += f[u][c];
            }
        }
        for u in 0..n {
            let mut grad = vec![0.0; k];
            let mut nbr_sum = vec![0.0; k];
            for &(v, _) in &data.knn_graph[u] {
                let dp: f64 = (0..k).map(|c| f[u][c] * f[v][c]).sum::<f64>().clamp(1e-4, 20.0);
                let coeff = (-dp).exp() / (1.0 - (-dp).exp()).max(1e-6);
                for c in 0..k {
                    grad[c] += f[v][c] * coeff;
                    nbr_sum[c] += f[v][c];
                }
            }
            for c in 0..k {
                let non_nbr = (f_sum[c] - f[u][c] - nbr_sum[c]).max(0.0);
                grad[c] -= non_nbr / (n as f64).max(1.0);
                f[u][c] = (f[u][c] + lr * grad[c]).clamp(0.0, 5.0);
            }
        }
    }

    let delta = (1.0 / n as f64).sqrt();
    let overlap = f.iter().filter(|row| row.iter().filter(|&&v| v > delta).count() > 1).count();
    let labels: Vec<i32> = (0..n)
        .map(|u| {
            (0..k)
                .max_by(|&a, &b| f[u][a].total_cmp(&f[u][b]))
                .unwrap_or(0) as i32
        })
        .collect();
    AlgorithmOutput::new(
        labels,
        format!("k={}, F>=0 gradientoppstigning", k),
        format!("BigCLAM affiliasjonsmatrise ({} overlappende noder)", overlap),
    )
}

// 126. SLPA (Xie et al. 2011 — Speaker-Listener Label Propagation with per-node memory buffer)
fn slpa(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let t_steps = cfg.max_iters(25);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 126);

    let mut memory: Vec<Vec<i32>> = (0..n).map(|i| vec![i as i32]).collect();
    let mut order: Vec<usize> = (0..n).collect();

    for _ in 0..t_steps {
        order.shuffle(&mut rng);
        for &listener in &order {
            if data.knn_graph[listener].is_empty() {
                continue;
            }
            let mut spoken: HashMap<i32, f64> = HashMap::new();
            for &(speaker, w) in &data.knn_graph[listener] {
                let mem = &memory[speaker];
                let lbl = mem[rng.gen_range(0..mem.len())];
                *spoken.entry(lbl).or_insert(0.0) += w;
            }
            if let Some((&most_pop, _)) = spoken.iter().max_by(|a, b| a.1.total_cmp(b.1)) {
                memory[listener].push(most_pop);
            }
        }
    }

    let r_thresh = 0.25;
    let mut overlap_cnt = 0usize;
    let mut labels = vec![0i32; n];
    for i in 0..n {
        let mut freq: HashMap<i32, usize> = HashMap::new();
        for &l in &memory[i] {
            *freq.entry(l).or_insert(0) += 1;
        }
        let tot = memory[i].len() as f64;
        let above = freq.values().filter(|&&c| (c as f64 / tot) >= r_thresh).count();
        if above > 1 {
            overlap_cnt += 1;
        }
        labels[i] = freq
            .into_iter()
            .max_by_key(|&(_, c)| c)
            .map(|(l, _)| l)
            .unwrap_or(0);
    }

    AlgorithmOutput::new(
        labels,
        format!("T={}, terskel r={:.2}", t_steps, r_thresh),
        format!("SLPA minne-propagering ({} overlappende noder)", overlap_cnt),
    )
}

// 127. CPM (Palla et al. Nature 2005 — Clique Percolation Method rolling 3-cliques/triangles)
fn clique_percolation(data: &PreparedData, _cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let adj_set: Vec<HashSet<usize>> = data
        .knn_graph
        .iter()
        .map(|nbrs| nbrs.iter().map(|&(j, _)| j).collect())
        .collect();

    // Find all 3-cliques (triangles)
    let mut cliques: Vec<[usize; 3]> = Vec::new();
    for u in 0..n {
        for &v in &adj_set[u] {
            if v <= u {
                continue;
            }
            for &w in adj_set[u].intersection(&adj_set[v]) {
                if w > v {
                    cliques.push([u, v, w]);
                }
            }
        }
    }

    if cliques.is_empty() {
        let (lbl, _) = louvain_core(&data.knn_graph, 1.0, 127);
        return AlgorithmOutput::new(lbl, "3-klikk (ingen triangler)".to_string(), "Fallback til Louvain".to_string());
    }

    // Two 3-cliques are adjacent iff they share 2 nodes (an edge)
    let mut edge_to_cliques: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
    for (ci, &[a, b, c]) in cliques.iter().enumerate() {
        edge_to_cliques.entry((a, b)).or_default().push(ci);
        edge_to_cliques.entry((a, c)).or_default().push(ci);
        edge_to_cliques.entry((b, c)).or_default().push(ci);
    }

    let m = cliques.len();
    let mut clique_comm = vec![-1i32; m];
    let mut cid = 0i32;
    for i in 0..m {
        if clique_comm[i] >= 0 {
            continue;
        }
        let mut q = VecDeque::new();
        clique_comm[i] = cid;
        q.push_back(i);
        while let Some(cur) = q.pop_front() {
            let [a, b, c] = cliques[cur];
            for edge in [(a, b), (a, c), (b, c)] {
                if let Some(nbr_cliques) = edge_to_cliques.get(&edge) {
                    for &nc in nbr_cliques {
                        if clique_comm[nc] < 0 {
                            clique_comm[nc] = cid;
                            q.push_back(nc);
                        }
                    }
                }
            }
        }
        cid += 1;
    }

    // Assign each node to its most frequent percolated clique community (unassigned -> -1)
    let mut node_votes: Vec<HashMap<i32, usize>> = vec![HashMap::new(); n];
    for (ci, &[a, b, c]) in cliques.iter().enumerate() {
        let cc = clique_comm[ci];
        *node_votes[a].entry(cc).or_insert(0) += 1;
        *node_votes[b].entry(cc).or_insert(0) += 1;
        *node_votes[c].entry(cc).or_insert(0) += 1;
    }
    let labels: Vec<i32> = node_votes
        .into_iter()
        .map(|v| v.into_iter().max_by_key(|&(_, cnt)| cnt).map(|(c, _)| c).unwrap_or(-1))
        .collect();
    AlgorithmOutput::new(
        labels,
        format!("k-klikk=3, {} triangler", m),
        format!("Perkolerte {} klikksamfunn gjennom k-NN-grafen", cid),
    )
}

// 128. DEMON (Coscia et al. KDD 2012 — Local ego-network label propagation + democratic merge)
fn demon(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let eps_merge = 0.35_f64;
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 128);
    let mut local_comms: Vec<HashSet<usize>> = Vec::new();

    // Sample ego-networks and run local LPA on each ego-minus-ego graph
    let step = (n / 35).max(1);
    for ego in (0..n).step_by(step) {
        let ego_nbrs: Vec<usize> = data.knn_graph[ego].iter().map(|&(v, _)| v).collect();
        if ego_nbrs.len() < 3 {
            continue;
        }
        let ego_set: HashSet<usize> = ego_nbrs.iter().copied().collect();
        let mut lbl: HashMap<usize, usize> = ego_nbrs.iter().map(|&v| (v, v)).collect();
        for _ in 0..6 {
            let mut ord = ego_nbrs.clone();
            ord.shuffle(&mut rng);
            for u in ord {
                let mut votes: HashMap<usize, f64> = HashMap::new();
                for &(v, w) in &data.knn_graph[u] {
                    if ego_set.contains(&v) {
                        *votes.entry(lbl[&v]).or_insert(0.0) += w;
                    }
                }
                if let Some((&best, _)) = votes.iter().max_by(|a, b| a.1.total_cmp(b.1)) {
                    lbl.insert(u, best);
                }
            }
        }
        let mut groups: HashMap<usize, HashSet<usize>> = HashMap::new();
        for (&u, &l) in &lbl {
            groups.entry(l).or_default().insert(u);
        }
        for mut g in groups.into_values() {
            g.insert(ego);
            if g.len() >= 3 {
                local_comms.push(g);
            }
        }
    }

    // Merge overlapping ego communities with Jaccard/subset overlap >= eps_merge
    let mut merged: Vec<HashSet<usize>> = Vec::new();
    for comm in local_comms {
        let mut did_merge = false;
        for existing in &mut merged {
            let inter = existing.intersection(&comm).count() as f64;
            let denom = existing.len().min(comm.len()).max(1) as f64;
            if inter / denom >= eps_merge {
                existing.extend(comm.iter().copied());
                did_merge = true;
                break;
            }
        }
        if !did_merge {
            merged.push(comm);
        }
    }

    let mut labels = vec![-1i32; n];
    for i in 0..n {
        if let Some((cid, _)) = merged
            .iter()
            .enumerate()
            .filter(|(_, g)| g.contains(&i))
            .max_by_key(|(_, g)| g.len())
        {
            labels[i] = cid as i32;
        }
    }
    AlgorithmOutput::new(
        labels,
        format!("eps_merge={:.2}, {} ego-samfunn", eps_merge, merged.len()),
        "Demokratisk lokal ego-nettverk LPA + sammenslåing".to_string(),
    )
}

// 129. Fast Unfolding (Blondel et al. 2008 — full two-phase iterative aggregation Louvain)
fn fast_unfolding_blondel(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let mut cur_adj = data.knn_graph.clone();
    let mut node_to_orig: Vec<Vec<usize>> = (0..n).map(|i| vec![i]).collect();
    let mut levels = 0usize;

    for lvl in 0..5 {
        let (lbl, _) = louvain_core(&cur_adj, 1.0, cfg.seed + 129 + lvl as u64);
        let mut groups: HashMap<i32, Vec<usize>> = HashMap::new();
        for (u, &c) in lbl.iter().enumerate() {
            groups.entry(c).or_default().push(u);
        }
        let m = groups.len();
        if m >= cur_adj.len() || m <= 1 {
            break;
        }
        levels += 1;
        let group_list: Vec<Vec<usize>> = groups.into_values().collect();
        let mut super_id = vec![0usize; cur_adj.len()];
        let mut next_orig = vec![Vec::new(); m];
        for (sid, members) in group_list.iter().enumerate() {
            for &u in members {
                super_id[u] = sid;
                next_orig[sid].extend_from_slice(&node_to_orig[u]);
            }
        }
        // Build aggregated super-node weighted graph
        let mut edge_w: HashMap<(usize, usize), f64> = HashMap::new();
        for u in 0..cur_adj.len() {
            for &(v, w) in &cur_adj[u] {
                let su = super_id[u];
                let sv = super_id[v];
                if su != sv {
                    *edge_w.entry((su, sv)).or_insert(0.0) += w;
                }
            }
        }
        let mut next_adj = vec![Vec::new(); m];
        for ((su, sv), w) in edge_w {
            next_adj[su].push((sv, w));
        }
        cur_adj = next_adj;
        node_to_orig = next_orig;
    }

    let mut labels = vec![0i32; n];
    for (cid, members) in node_to_orig.iter().enumerate() {
        for &p in members {
            labels[p] = cid as i32;
        }
    }
    let q = graph_modularity(&data.knn_graph, &labels);
    AlgorithmOutput::new(
        labels,
        format!("{} hierarkiske foldingsnivåer", levels.max(1)),
        format!("Blondel 2-fase modularitetsfolding (Q={:.3})", q),
    )
}

// 130 & 131. SBM & Degree-Corrected SBM (Karrer & Newman 2011 — Variational EM blockmodel)
fn sbm(data: &PreparedData, cfg: &RunConfig, degree_corrected: bool) -> AlgorithmOutput {
    let n = data.n_samples;
    let k = cfg.effective_k(n);
    let max_iter = cfg.max_iters(35);
    let (w, deg, _) = dense_graph_from_knn(data);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + if degree_corrected { 131 } else { 130 });

    // Initialize via K-Means++ soft membership
    let c0 = init_kmeans_pp(&data.x, k, &mut rng);
    let (lbl0, _, _, _) = run_lloyd_from_centers(&data.x, c0, 10);
    let mut q_mat = vec![vec![0.05 / k as f64; k]; n];
    for i in 0..n {
        q_mat[i][lbl0[i] as usize] += 0.95;
    }

    let mut omega = vec![vec![0.1; k]; k];
    let mut pi = vec![1.0 / k as f64; k];

    for _ in 0..max_iter {
        // M-step: estimate block affinity matrix omega[r][s]
        let mut nk = vec![1e-4; k];
        let mut kappa = vec![1e-4; k];
        for i in 0..n {
            for r in 0..k {
                nk[r] += q_mat[i][r];
                kappa[r] += q_mat[i][r] * deg[i];
            }
        }
        for r in 0..k {
            pi[r] = nk[r] / n as f64;
            for s in r..k {
                let mut e_rs = 1e-4;
                for i in 0..n {
                    for &(j, wt) in &data.knn_graph[i] {
                        e_rs += q_mat[i][r] * q_mat[j][s] * wt;
                    }
                }
                let val = if degree_corrected {
                    e_rs / (kappa[r] * kappa[s]).max(1e-6)
                } else {
                    (e_rs / (nk[r] * nk[s]).max(1e-6)).clamp(1e-4, 0.999)
                };
                omega[r][s] = val.max(1e-5);
                omega[s][r] = omega[r][s];
            }
        }

        // E-step: variational mean-field update for node block posteriors
        let mut next_q = vec![vec![0.0; k]; n];
        for i in 0..n {
            let mut lp = vec![0.0; k];
            let mut mx = f64::NEG_INFINITY;
            for r in 0..k {
                let mut l = pi[r].max(1e-8).ln();
                if degree_corrected {
                    // Poisson DC-SBM log-likelihood contribution
                    let mut expected = 0.0;
                    for s in 0..k {
                        expected += deg[i] * kappa[s] * omega[r][s];
                    }
                    let mut obs = 0.0;
                    for &(j, wt) in &data.knn_graph[i] {
                        for s in 0..k {
                            obs += wt * q_mat[j][s] * (omega[r][s] * deg[i].max(1e-4)).ln();
                        }
                    }
                    l += obs - expected;
                } else {
                    for j in 0..n {
                        if i == j {
                            continue;
                        }
                        let a_ij = if w[i][j] > 0.0 { 1.0 } else { 0.0 };
                        for s in 0..k {
                            let p = omega[r][s].clamp(1e-4, 0.999);
                            l += q_mat[j][s] * (a_ij * p.ln() + (1.0 - a_ij) * (1.0 - p).ln());
                        }
                    }
                }
                lp[r] = l;
                mx = mx.max(l);
            }
            let sp: f64 = lp.iter().map(|&v| (v - mx).exp()).sum::<f64>().max(1e-15);
            for r in 0..k {
                next_q[i][r] = 0.5 * q_mat[i][r] + 0.5 * ((lp[r] - mx).exp() / sp);
            }
        }
        q_mat = next_q;
    }

    let labels: Vec<i32> = q_mat
        .iter()
        .map(|row| {
            row.iter()
                .enumerate()
                .max_by(|a, b| a.1.total_cmp(b.1))
                .map(|(c, _)| c as i32)
                .unwrap_or(0)
        })
        .collect();
    let tag = if degree_corrected { "DC-SBM (gradkorrigert)" } else { "Standard SBM" };
    AlgorithmOutput::new(
        labels,
        format!("k={}, {}", k, tag),
        format!("Variasjons-EM på stokastisk blokkmodell ({})", tag),
    )
}

// 132. Power Iteration Clustering (Lin & Cohen ICML 2010 — truncated power iteration v^(t) = W v^(t-1))
fn power_iteration_clustering(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let k = cfg.effective_k(n);
    let sigma = median_pairwise_dist(&data.dist_matrix).max(1e-4);
    let aff = rbf_affinity_matrix(&data.dist_matrix, sigma);

    // Row-normalized affinity matrix W = D^-1 A
    let mut w = vec![vec![0.0; n]; n];
    for i in 0..n {
        let s: f64 = aff[i].iter().sum::<f64>().max(1e-12);
        for j in 0..n {
            w[i][j] = aff[i][j] / s;
        }
    }

    // Initialize v^0 proportional to degree (or random) and iterate until velocity delta converges
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 132);
    let mut v: Vec<f64> = (0..n).map(|_| rng.gen_range(0.1..1.0)).collect();
    let v_sum: f64 = v.iter().sum();
    for x in &mut v {
        *x /= v_sum;
    }

    let mut prev_delta = vec![0.0; n];
    let eps = 1e-5 / n as f64;
    let mut iters = 0;
    for it in 0..cfg.max_iters(60) {
        iters = it + 1;
        let mut next = vec![0.0; n];
        for i in 0..n {
            let mut s = 0.0;
            for j in 0..n {
                s += w[i][j] * v[j];
            }
            next[i] = s;
        }
        let l1: f64 = next.iter().map(|x| x.abs()).sum::<f64>().max(1e-12);
        for x in &mut next {
            *x /= l1;
        }
        let delta: Vec<f64> = (0..n).map(|i| (next[i] - v[i]).abs()).collect();
        let accel: f64 = (0..n).map(|i| (delta[i] - prev_delta[i]).abs()).fold(0.0, f64::max);
        v = next;
        if it > 3 && accel < eps {
            break;
        }
        prev_delta = delta;
    }

    // 1D embedding clustering via K-Means++
    let embed: Vec<Vec<f64>> = v.iter().map(|&val| vec![val]).collect();
    let c0 = init_kmeans_pp(&embed, k, &mut rng);
    let (labels, _, _, _) = run_lloyd_from_centers(&embed, c0, 40);
    AlgorithmOutput::new(
        labels,
        format!("k={}, avkortet potensiterasjon", k),
        format!("1D PIC-embedding stoppet etter {} iterasjoner", iters),
    )
}
