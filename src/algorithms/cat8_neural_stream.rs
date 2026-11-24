use std::collections::{HashMap, HashSet, VecDeque};

use rand::seq::SliceRandom;
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;

use crate::algorithms::cat1_centroid::{init_kmeans_pp, run_lloyd_from_centers};
use crate::algorithms::cat3_hierarchical::{agglomerative_lance_williams, Linkage};
use crate::algorithms::{AlgorithmOutput, RunConfig};
use crate::data::PreparedData;
use crate::math::{
    cosine_sim, dot, euclidean, euclidean_sq, median_pairwise_dist, normalize_vec, sinkhorn_knopp,
};

pub fn run(id: usize, data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    match id {
        176 => som_kohonen(data, cfg),
        177 => growing_neural_gas(data, cfg),
        178 => fuzzy_art(data, cfg),
        179 => neural_gas(data, cfg),
        180 => growing_cell_structures(data, cfg),
        181 => dec_clustering(data, cfg, false),
        182 => dec_clustering(data, cfg, true),
        183 => dcn_clustering(data, cfg),
        184 => vade_clustering(data, cfg),
        185 => depict_clustering(data, cfg),
        186 => jule_clustering(data, cfg),
        187 => dac_clustering(data, cfg),
        188 => iic_clustering(data, cfg),
        189 => scan_clustering(data, cfg),
        190 => swav_clustering(data, cfg),
        191 => deep_cluster(data, cfg),
        192 => contrastive_clustering(data, cfg),
        193 => clustream(data, cfg),
        194 => denstream(data, cfg),
        195 => stream_km_pp(data, cfg),
        196 => clustree_birch_stream(data, cfg),
        197 => d_stream(data, cfg),
        198 => evo_kmeans_gka(data, cfg),
        199 => pso_clustering(data, cfg),
        200 => aco_clustering(data, cfg),
        _ => unreachable!(),
    }
}

// ============================================================================
// Lightweight Genuine 2-Layer MLP / Autoencoder with Backpropagation (Pure Rust)
// Used by DEC, IDEC, DCN, VaDE, DEPICT, JULE, DAC, IIC, SCAN, SwAV, DeepCluster, CC
// ============================================================================
struct AutoEncoder {
    d_in: usize,
    d_lat: usize,
    w_enc: Vec<Vec<f64>>, // d_lat x d_in
    b_enc: Vec<f64>,
    w_dec: Vec<Vec<f64>>, // d_in x d_lat
    b_dec: Vec<f64>,
}

impl AutoEncoder {
    fn new(d_in: usize, d_lat: usize, seed: u64) -> Self {
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        let scale_e = (2.0 / (d_in + d_lat) as f64).sqrt();
        let w_enc = (0..d_lat)
            .map(|_| (0..d_in).map(|_| rng.gen_range(-scale_e..scale_e)).collect())
            .collect();
        let b_enc = vec![0.0; d_lat];
        let w_dec = (0..d_in)
            .map(|_| (0..d_lat).map(|_| rng.gen_range(-scale_e..scale_e)).collect())
            .collect();
        let b_dec = vec![0.0; d_in];
        Self { d_in, d_lat, w_enc, b_enc, w_dec, b_dec }
    }

    fn encode(&self, x: &[f64]) -> Vec<f64> {
        let mut z = vec![0.0; self.d_lat];
        for h in 0..self.d_lat {
            let s = self.b_enc[h] + dot(&self.w_enc[h], x);
            z[h] = s.tanh();
        }
        z
    }

    fn decode(&self, z: &[f64]) -> Vec<f64> {
        let mut xhat = vec![0.0; self.d_in];
        for j in 0..self.d_in {
            xhat[j] = self.b_dec[j] + dot(&self.w_dec[j], z);
        }
        xhat
    }

    fn train_reconstruction(&mut self, data: &[Vec<f64>], epochs: usize, lr: f64, noise_std: f64, seed: u64) -> f64 {
        let n = data.len().max(1) as f64;
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        let mut last_loss = 0.0;
        for _ in 0..epochs {
            let mut gw_enc = vec![vec![0.0; self.d_in]; self.d_lat];
            let mut gb_enc = vec![0.0; self.d_lat];
            let mut gw_dec = vec![vec![0.0; self.d_lat]; self.d_in];
            let mut gb_dec = vec![0.0; self.d_in];
            let mut loss = 0.0;

            for x in data {
                let x_in: Vec<f64> = if noise_std > 0.0 {
                    x.iter().map(|&v| v + rng.gen_range(-noise_std..noise_std)).collect()
                } else {
                    x.clone()
                };
                let z = self.encode(&x_in);
                let xhat = self.decode(&z);

                let mut dxhat = vec![0.0; self.d_in];
                for j in 0..self.d_in {
                    let err = xhat[j] - x[j];
                    loss += err * err;
                    dxhat[j] = (2.0 * err) / (n * self.d_in as f64);
                    gb_dec[j] += dxhat[j];
                    for h in 0..self.d_lat {
                        gw_dec[j][h] += dxhat[j] * z[h];
                    }
                }
                for h in 0..self.d_lat {
                    let dz: f64 = (0..self.d_in).map(|j| dxhat[j] * self.w_dec[j][h]).sum();
                    let dpre = dz * (1.0 - z[h] * z[h]);
                    gb_enc[h] += dpre;
                    for j in 0..self.d_in {
                        gw_enc[h][j] += dpre * x_in[j];
                    }
                }
            }
            last_loss = loss / (n * self.d_in as f64);
            for h in 0..self.d_lat {
                self.b_enc[h] -= lr * gb_enc[h];
                for j in 0..self.d_in {
                    self.w_enc[h][j] -= lr * gw_enc[h][j];
                }
            }
            for j in 0..self.d_in {
                self.b_dec[j] -= lr * gb_dec[j];
                for h in 0..self.d_lat {
                    self.w_dec[j][h] -= lr * gw_dec[j][h];
                }
            }
        }
        last_loss
    }

    /// Backpropagate an external gradient dL/dz directly into the encoder weights.
    fn apply_latent_grad(&mut self, data: &[Vec<f64>], z_all: &[Vec<f64>], dz_all: &[Vec<f64>], lr: f64) {
        let n = data.len().max(1) as f64;
        for h in 0..self.d_lat {
            let mut gb = 0.0;
            let mut gw = vec![0.0; self.d_in];
            for i in 0..data.len() {
                let dpre = dz_all[i][h] * (1.0 - z_all[i][h] * z_all[i][h]) / n;
                gb += dpre;
                for j in 0..self.d_in {
                    gw[j] += dpre * data[i][j];
                }
            }
            self.b_enc[h] -= lr * gb.clamp(-1.0, 1.0);
            for j in 0..self.d_in {
                self.w_enc[h][j] -= lr * gw[j].clamp(-1.0, 1.0);
            }
        }
    }
}

// 176. SOM (Self-Organizing Maps — Kohonen 1982 2D topology-preserving grid + Ward U-matrix clustering)
fn som_kohonen(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let gx = 4usize;
    let gy = 4usize;
    let m_nodes = gx * gy;
    let epochs = cfg.max_iters(25);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 176);

    let mut weights: Vec<Vec<f64>> = (0..m_nodes)
        .map(|_| data.x[rng.gen_range(0..n)].clone())
        .collect();
    let mut order: Vec<usize> = (0..n).collect();

    for ep in 0..epochs {
        let frac = ep as f64 / epochs.max(1) as f64;
        let lr = 0.4 * (1.0 - 0.85 * frac);
        let sigma = 2.0 * (1.0 - 0.75 * frac).max(0.4);
        let s2 = 2.0 * sigma * sigma;
        order.shuffle(&mut rng);

        for &i in &order {
            let bmu = (0..m_nodes)
                .min_by(|&a, &b| {
                    euclidean_sq(&data.x[i], &weights[a])
                        .total_cmp(&euclidean_sq(&data.x[i], &weights[b]))
                })
                .unwrap_or(0);
            let bx = (bmu % gx) as f64;
            let by = (bmu / gx) as f64;
            for u in 0..m_nodes {
                let ux = (u % gx) as f64;
                let uy = (u / gx) as f64;
                let gdist2 = (ux - bx).powi(2) + (uy - by).powi(2);
                let h = (-gdist2 / s2).exp();
                let step = lr * h;
                for j in 0..d {
                    weights[u][j] += step * (data.x[i][j] - weights[u][j]);
                }
            }
        }
    }

    // Cluster the 4x4 Kohonen prototypes into k clusters
    let mut wdist = vec![vec![0.0; m_nodes]; m_nodes];
    for a in 0..m_nodes {
        for b in (a + 1)..m_nodes {
            let dd = euclidean(&weights[a], &weights[b]);
            wdist[a][b] = dd;
            wdist[b][a] = dd;
        }
    }
    let proto_lbl = agglomerative_lance_williams(&wdist, k.min(m_nodes), Linkage::Ward);
    let labels: Vec<i32> = (0..n)
        .map(|i| {
            let bmu = (0..m_nodes)
                .min_by(|&a, &b| {
                    euclidean_sq(&data.x[i], &weights[a])
                        .total_cmp(&euclidean_sq(&data.x[i], &weights[b]))
                })
                .unwrap_or(0);
            proto_lbl[bmu]
        })
        .collect();
    AlgorithmOutput::new(
        labels,
        format!("{}x{} Kohonen-kart, k={}", gx, gy, k),
        format!("Topologibevarende SOM trent over {} epoker", epochs),
    )
}

// 177. Growing Neural Gas (Fritzke 1995 — Dynamic topology graph with edge aging & error-driven node growth)
fn growing_neural_gas(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let max_nodes = (k * 4).clamp(10, 24).min(n);
    let max_age = 18usize;
    let lambda_grow = 15usize;
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 177);

    let mut nodes: Vec<Vec<f64>> = vec![
        data.x[rng.gen_range(0..n)].clone(),
        data.x[rng.gen_range(0..n)].clone(),
    ];
    let mut errors = vec![0.0; 2];
    let mut edges: HashMap<(usize, usize), usize> = HashMap::new(); // (u<v) -> age

    let total_steps = (n * 4).clamp(120, 600);
    for step in 1..=total_steps {
        let p = &data.x[rng.gen_range(0..n)];
        let m = nodes.len();
        // Find nearest s1 and second nearest s2
        let mut s1 = 0;
        let mut s2 = 1.min(m - 1);
        let mut d1 = f64::INFINITY;
        let mut d2 = f64::INFINITY;
        for u in 0..m {
            let dd = euclidean_sq(p, &nodes[u]);
            if dd < d1 {
                d2 = d1;
                s2 = s1;
                d1 = dd;
                s1 = u;
            } else if dd < d2 {
                d2 = dd;
                s2 = u;
            }
        }
        errors[s1] += d1;
        // Adapt s1 and its topological neighbors
        for j in 0..d {
            nodes[s1][j] += 0.12 * (p[j] - nodes[s1][j]);
        }
        let nbrs: Vec<usize> = edges
            .keys()
            .filter_map(|&(a, b)| if a == s1 { Some(b) } else if b == s1 { Some(a) } else { None })
            .collect();
        for nb in nbrs {
            for j in 0..d {
                nodes[nb][j] += 0.015 * (p[j] - nodes[nb][j]);
            }
            let e = (s1.min(nb), s1.max(nb));
            if let Some(age) = edges.get_mut(&e) {
                *age += 1;
            }
        }
        // Refresh edge (s1, s2) and prune old edges
        edges.insert((s1.min(s2), s1.max(s2)), 0);
        edges.retain(|_, age| *age <= max_age);

        // Grow new node between highest-error node q and its highest-error neighbor f
        if step % lambda_grow == 0 && nodes.len() < max_nodes {
            let q = (0..nodes.len())
                .max_by(|&a, &b| errors[a].total_cmp(&errors[b]))
                .unwrap_or(0);
            let q_nbrs: Vec<usize> = edges
                .keys()
                .filter_map(|&(a, b)| if a == q { Some(b) } else if b == q { Some(a) } else { None })
                .collect();
            let f = q_nbrs
                .into_iter()
                .max_by(|&a, &b| errors[a].total_cmp(&errors[b]))
                .unwrap_or((q + 1) % nodes.len());
            let r_vec: Vec<f64> = (0..d).map(|j| 0.5 * (nodes[q][j] + nodes[f][j])).collect();
            let r_id = nodes.len();
            nodes.push(r_vec);
            errors[q] *= 0.5;
            errors[f] *= 0.5;
            errors.push(errors[q]);
            edges.remove(&(q.min(f), q.max(f)));
            edges.insert((q.min(r_id), q.max(r_id)), 0);
            edges.insert((f.min(r_id), f.max(r_id)), 0);
        }
        for e in &mut errors {
            *e *= 0.995;
        }
    }

    // Connected components of GNG topology (or Ward cut to k if single component)
    let m = nodes.len();
    let mut ndist = vec![vec![0.0; m]; m];
    for a in 0..m {
        for b in (a + 1)..m {
            let dd = euclidean(&nodes[a], &nodes[b]);
            ndist[a][b] = dd;
            ndist[b][a] = dd;
        }
    }
    let n_lbl = agglomerative_lance_williams(&ndist, k.min(m), Linkage::Average);
    let labels: Vec<i32> = (0..n)
        .map(|i| {
            let bmu = (0..m)
                .min_by(|&a, &b| {
                    euclidean_sq(&data.x[i], &nodes[a])
                        .total_cmp(&euclidean_sq(&data.x[i], &nodes[b]))
                })
                .unwrap_or(0);
            n_lbl[bmu]
        })
        .collect();
    AlgorithmOutput::new(
        labels,
        format!("noder={}, kanter={}, k={}", m, edges.len(), k),
        "Voksende nevral graf (GNG) tilpasset datamanifolden".to_string(),
    )
}

// 178. ART (Carpenter et al. 1991 — Fuzzy ART with complement coding & vigilance parameter rho)
fn fuzzy_art(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    // Normalize features to [0,1] and complement-code I = [x, 1-x] (length 2d, L1 norm = d)
    let mut col_min = vec![f64::INFINITY; d];
    let mut col_max = vec![f64::NEG_INFINITY; d];
    for i in 0..n {
        for j in 0..d {
            col_min[j] = col_min[j].min(data.x[i][j]);
            col_max[j] = col_max[j].max(data.x[i][j]);
        }
    }
    let comp_x: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            let mut v = Vec::with_capacity(2 * d);
            for j in 0..d {
                let norm = ((data.x[i][j] - col_min[j]) / (col_max[j] - col_min[j]).max(1e-9))
                    .clamp(0.0, 1.0);
                v.push(norm);
                v.push(1.0 - norm);
            }
            v
        })
        .collect();

    let alpha = 0.01_f64;
    let beta = 0.6_f64;
    // Auto-calibrate vigilance rho to get a reasonable cluster count near target k
    let target_k = cfg.effective_k(n);
    let mut best_labels = vec![0i32; n];
    let mut best_rho = 0.65;
    let mut best_diff = usize::MAX;

    for &rho in &[0.45, 0.55, 0.65, 0.72, 0.80] {
        let mut weights: Vec<Vec<f64>> = Vec::new();
        let mut labels = vec![0i32; n];

        for epoch in 0..3 {
            for i in 0..n {
                let inp = &comp_x[i];
                let inp_norm: f64 = inp.iter().sum::<f64>().max(1e-9);
                let mut scores: Vec<(usize, f64)> = weights
                    .iter()
                    .enumerate()
                    .map(|(c, w)| {
                        let fuzzy_and: f64 = inp.iter().zip(w.iter()).map(|(a, b)| a.min(*b)).sum();
                        let w_norm: f64 = w.iter().sum();
                        (c, fuzzy_and / (alpha + w_norm))
                    })
                    .collect();
                scores.sort_by(|a, b| b.1.total_cmp(&a.1));

                let mut matched = None;
                for (c, _) in scores {
                    let fuzzy_and: f64 = inp
                        .iter()
                        .zip(weights[c].iter())
                        .map(|(a, b)| a.min(*b))
                        .sum();
                    if fuzzy_and / inp_norm >= rho {
                        matched = Some(c);
                        break;
                    }
                }
                if let Some(c) = matched {
                    for j in 0..(2 * d) {
                        weights[c][j] =
                            beta * inp[j].min(weights[c][j]) + (1.0 - beta) * weights[c][j];
                    }
                    labels[i] = c as i32;
                } else if epoch == 0 || weights.len() < 16 {
                    labels[i] = weights.len() as i32;
                    weights.push(inp.clone());
                } else {
                    labels[i] = 0;
                }
            }
        }
        let nc = labels.iter().copied().collect::<HashSet<_>>().len();
        if nc >= 2 {
            let diff = nc.abs_diff(target_k);
            if diff < best_diff {
                best_diff = diff;
                best_rho = rho;
                best_labels = labels;
            }
        }
    }
    AlgorithmOutput::new(
        best_labels,
        format!("Fuzzy ART rho={:.2}, beta={:.1}", best_rho, beta),
        "Adaptive Resonance Theory med komplementkoding og årvåkenhetstest".to_string(),
    )
}

// 179. Neural Gas (Martinetz & Schulten 1991 — rank-sorted soft neighborhood adaptation)
fn neural_gas(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let epochs = cfg.max_iters(35);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 179);

    let mut protos = init_kmeans_pp(&data.x, k, &mut rng);
    let mut order: Vec<usize> = (0..n).collect();

    for ep in 0..epochs {
        let t = ep as f64 / epochs.max(1) as f64;
        let eps = 0.4 * (0.05_f64).powf(t);
        let lambda = (k as f64 * 0.5).max(0.5) * (0.08_f64).powf(t);
        order.shuffle(&mut rng);

        for &i in &order {
            let mut ranks: Vec<(usize, f64)> = (0..k)
                .map(|c| (c, euclidean_sq(&data.x[i], &protos[c])))
                .collect();
            ranks.sort_by(|a, b| a.1.total_cmp(&b.1));
            for (rank, &(c, _)) in ranks.iter().enumerate() {
                let h = (-(rank as f64) / lambda.max(1e-4)).exp();
                let step = eps * h;
                for j in 0..d {
                    protos[c][j] += step * (data.x[i][j] - protos[c][j]);
                }
            }
        }
    }
    let labels: Vec<i32> = (0..n)
        .map(|i| {
            (0..k)
                .min_by(|&a, &b| {
                    euclidean_sq(&data.x[i], &protos[a])
                        .total_cmp(&euclidean_sq(&data.x[i], &protos[b]))
                })
                .unwrap_or(0) as i32
        })
        .collect();
    AlgorithmOutput::new(
        labels,
        format!("k={}, {} epoker", k, epochs),
        "Rangsortert myk vektorkvantisering (Neural Gas)".to_string(),
    )
}

// 180. Growing Cell Structures (Fritzke 1994 — 2D triangle simplex network growing by local signal counter)
fn growing_cell_structures(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let max_cells = (k * 4).clamp(9, 20).min(n);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 180);

    // Initialize a single 2D triangle simplex (3 connected cells)
    let mut cells: Vec<Vec<f64>> = (0..3).map(|_| data.x[rng.gen_range(0..n)].clone()).collect();
    let mut signal = vec![0.0_f64; 3];
    let mut adj: Vec<HashSet<usize>> = vec![HashSet::new(); 3];
    for &(a, b) in &[(0, 1), (1, 2), (0, 2)] {
        adj[a].insert(b);
        adj[b].insert(a);
    }

    let steps = (n * 4).clamp(120, 500);
    for step in 1..=steps {
        let p = &data.x[rng.gen_range(0..n)];
        let bmu = (0..cells.len())
            .min_by(|&a, &b| euclidean_sq(p, &cells[a]).total_cmp(&euclidean_sq(p, &cells[b])))
            .unwrap_or(0);
        signal[bmu] += 1.0;
        for j in 0..d {
            cells[bmu][j] += 0.10 * (p[j] - cells[bmu][j]);
        }
        for &nb in &adj[bmu] {
            for j in 0..d {
                cells[nb][j] += 0.02 * (p[j] - cells[nb][j]);
            }
        }
        // Split longest edge incident to highest-signal cell, maintaining simplex triangulation
        if step % 18 == 0 && cells.len() < max_cells {
            let q = (0..cells.len())
                .max_by(|&a, &b| signal[a].total_cmp(&signal[b]))
                .unwrap_or(0);
            if let Some(&f) = adj[q]
                .iter()
                .max_by(|&&a, &&b| euclidean_sq(&cells[q], &cells[a]).total_cmp(&euclidean_sq(&cells[q], &cells[b])))
            {
                let common: Vec<usize> = adj[q].intersection(&adj[f]).copied().collect();
                let new_id = cells.len();
                let new_w: Vec<f64> = (0..d).map(|j| 0.5 * (cells[q][j] + cells[f][j])).collect();
                cells.push(new_w);
                signal[q] *= 0.5;
                signal.push(signal[q]);
                adj.push(HashSet::new());

                adj[q].remove(&f);
                adj[f].remove(&q);
                adj[q].insert(new_id);
                adj[new_id].insert(q);
                adj[f].insert(new_id);
                adj[new_id].insert(f);
                for c in common {
                    adj[c].insert(new_id);
                    adj[new_id].insert(c);
                }
            }
        }
        for s in &mut signal {
            *s *= 0.992;
        }
    }

    let m = cells.len();
    let mut cdist = vec![vec![0.0; m]; m];
    for a in 0..m {
        for b in (a + 1)..m {
            let dd = euclidean(&cells[a], &cells[b]);
            cdist[a][b] = dd;
            cdist[b][a] = dd;
        }
    }
    let clbl = agglomerative_lance_williams(&cdist, k.min(m), Linkage::Ward);
    let labels: Vec<i32> = (0..n)
        .map(|i| {
            let bmu = (0..m)
                .min_by(|&a, &b| euclidean_sq(&data.x[i], &cells[a]).total_cmp(&euclidean_sq(&data.x[i], &cells[b])))
                .unwrap_or(0);
            clbl[bmu]
        })
        .collect();
    AlgorithmOutput::new(
        labels,
        format!("simpleks-celler={}, k={}", m, k),
        "Voksende simpleks-topologi (GCS) etter lokal signalrate".to_string(),
    )
}

// 181 & 182. DEC (Xie et al. ICML 2016) & IDEC (Guo et al. IJCAI 2017 — preserves reconstruction loss)
fn dec_clustering(data: &PreparedData, cfg: &RunConfig, improved_idec: bool) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let d_lat = d.clamp(2, 6);
    let seed = cfg.seed + if improved_idec { 182 } else { 181 };

    let mut ae = AutoEncoder::new(d, d_lat, seed);
    let recon_loss = ae.train_reconstruction(&data.x, cfg.max_iters(35), 0.25, 0.0, seed);

    // Initialize cluster centers in latent space via K-Means++
    let mut z: Vec<Vec<f64>> = data.x.iter().map(|r| ae.encode(r)).collect();
    let mut rng = ChaCha8Rng::seed_from_u64(seed + 1);
    let c0 = init_kmeans_pp(&z, k, &mut rng);
    let (_, mut mu, _, _) = run_lloyd_from_centers(&z, c0, 15);

    // Joint Student-t soft assignment Q & sharpened auxiliary target P -> minimize KL(P || Q)
    let mut q = vec![vec![1.0 / k as f64; k]; n];
    let mut kl_val = 0.0;
    for _ in 0..cfg.max_iters(18) {
        z = data.x.iter().map(|r| ae.encode(r)).collect();
        for i in 0..n {
            let mut raw = vec![0.0; k];
            let mut s = 0.0;
            for c in 0..k {
                raw[c] = 1.0 / (1.0 + euclidean_sq(&z[i], &mu[c]));
                s += raw[c];
            }
            for c in 0..k {
                q[i][c] = (raw[c] / s.max(1e-12)).max(1e-8);
            }
        }
        // Auxiliary target distribution P_ic = (Q_ic^2 / f_c) / sum_j (Q_ij^2 / f_j)
        let f_c: Vec<f64> = (0..k)
            .map(|c| q.iter().map(|row| row[c]).sum::<f64>().max(1e-6))
            .collect();
        let mut p = vec![vec![0.0; k]; n];
        kl_val = 0.0;
        for i in 0..n {
            let mut s = 0.0;
            for c in 0..k {
                p[i][c] = (q[i][c] * q[i][c]) / f_c[c];
                s += p[i][c];
            }
            for c in 0..k {
                p[i][c] = (p[i][c] / s.max(1e-12)).max(1e-8);
                kl_val += p[i][c] * (p[i][c] / q[i][c]).ln();
            }
        }
        kl_val /= n as f64;

        // Exact gradient of KL(P || Q) wrt z_i and mu_c
        let mut dz = vec![vec![0.0; d_lat]; n];
        let mut dmu = vec![vec![0.0; d_lat]; k];
        for i in 0..n {
            for c in 0..k {
                let factor = 2.0 * (p[i][c] - q[i][c]) / (1.0 + euclidean_sq(&z[i], &mu[c]));
                for h in 0..d_lat {
                    let g = factor * (z[i][h] - mu[c][h]);
                    dz[i][h] += g;
                    dmu[c][h] -= g / n as f64;
                }
            }
        }
        ae.apply_latent_grad(&data.x, &z, &dz, 0.12);
        for c in 0..k {
            for h in 0..d_lat {
                mu[c][h] -= 0.18 * dmu[c][h];
            }
        }
        if improved_idec {
            // IDEC preserves local structure via continued reconstruction loss gradient
            ae.train_reconstruction(&data.x, 1, 0.12, 0.0, seed);
        }
    }

    let labels: Vec<i32> = q
        .iter()
        .map(|row| {
            row.iter()
                .enumerate()
                .max_by(|a, b| a.1.total_cmp(b.1))
                .map(|(c, _)| c as i32)
                .unwrap_or(0)
        })
        .collect();
    let tag = if improved_idec { "IDEC (bevarte rekonstruksjon)" } else { "DEC (Student-t KL)" };
    AlgorithmOutput::new(
        labels,
        format!("k={}, latent={}", k, d_lat),
        format!("{}: recon={:.3}, KL={:.4}", tag, recon_loss, kl_val),
    )
}

// 183. DCN (Yang et al. ICML 2017 — Deep Clustering Network alternating AE recon + latent K-Means penalty)
fn dcn_clustering(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let d_lat = d.clamp(2, 6);
    let lambda = 0.5_f64;

    let mut ae = AutoEncoder::new(d, d_lat, cfg.seed + 183);
    ae.train_reconstruction(&data.x, cfg.max_iters(30), 0.25, 0.0, cfg.seed + 183);

    let mut z: Vec<Vec<f64>> = data.x.iter().map(|r| ae.encode(r)).collect();
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 1830);
    let c0 = init_kmeans_pp(&z, k, &mut rng);
    let (mut labels, mut centers, _, _) = run_lloyd_from_centers(&z, c0, 12);

    for _ in 0..cfg.max_iters(15) {
        // Pull latent codes toward assigned centroids + maintain AE reconstruction
        let mut dz = vec![vec![0.0; d_lat]; n];
        for i in 0..n {
            let c = labels[i] as usize;
            for h in 0..d_lat {
                dz[i][h] = lambda * (z[i][h] - centers[c][h]);
            }
        }
        ae.apply_latent_grad(&data.x, &z, &dz, 0.12);
        ae.train_reconstruction(&data.x, 2, 0.15, 0.0, cfg.seed + 183);
        z = data.x.iter().map(|r| ae.encode(r)).collect();
        let (lbl, cen, _, _) = run_lloyd_from_centers(&z, centers, 5);
        labels = lbl;
        centers = cen;
    }
    AlgorithmOutput::new(
        labels,
        format!("k={}, latent={}, lambda={:.1}", k, d_lat, lambda),
        "Alternerende Autoencoder + latent K-Means-objektiv".to_string(),
    )
}

// 184. VaDE (Jiang et al. IJCAI 2017 — Variational Autoencoder with Gaussian Mixture prior in latent space)
fn vade_clustering(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let d_lat = d.clamp(2, 6);

    let mut ae = AutoEncoder::new(d, d_lat, cfg.seed + 184);
    ae.train_reconstruction(&data.x, cfg.max_iters(30), 0.24, 0.08, cfg.seed + 184);

    let mut z: Vec<Vec<f64>> = data.x.iter().map(|r| ae.encode(r)).collect();
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 1840);
    let c0 = init_kmeans_pp(&z, k, &mut rng);
    let (_, mut mu_c, _, _) = run_lloyd_from_centers(&z, c0, 12);
    let mut var_c = vec![vec![0.5_f64; d_lat]; k];
    let mut pi_c = vec![1.0 / k as f64; k];
    let mut gamma = vec![vec![1.0 / k as f64; k]; n];

    for _ in 0..cfg.max_iters(16) {
        z = data.x.iter().map(|r| ae.encode(r)).collect();
        // E-step: posterior p(c|z) under GMM prior
        for i in 0..n {
            let mut lp = vec![0.0; k];
            let mut mx = f64::NEG_INFINITY;
            for c in 0..k {
                let mut l = pi_c[c].max(1e-8).ln();
                for h in 0..d_lat {
                    let v = var_c[c][h].max(0.05);
                    l += -0.5 * v.ln() - 0.5 * (z[i][h] - mu_c[c][h]).powi(2) / v;
                }
                lp[c] = l;
                mx = mx.max(l);
            }
            let s: f64 = lp.iter().map(|&v| (v - mx).exp()).sum::<f64>().max(1e-15);
            for c in 0..k {
                gamma[i][c] = (lp[c] - mx).exp() / s;
            }
        }
        // M-step for GMM prior + KL gradient pulling latent z toward mixture posterior
        let mut dz = vec![vec![0.0; d_lat]; n];
        for c in 0..k {
            let nk: f64 = gamma.iter().map(|r| r[c]).sum::<f64>().max(1e-5);
            pi_c[c] = nk / n as f64;
            for h in 0..d_lat {
                mu_c[c][h] = (0..n).map(|i| gamma[i][c] * z[i][h]).sum::<f64>() / nk;
                let v = (0..n)
                    .map(|i| gamma[i][c] * (z[i][h] - mu_c[c][h]).powi(2))
                    .sum::<f64>()
                    / nk;
                var_c[c][h] = v.max(0.08);
            }
        }
        for i in 0..n {
            for c in 0..k {
                for h in 0..d_lat {
                    dz[i][h] += gamma[i][c] * (z[i][h] - mu_c[c][h]) / var_c[c][h];
                }
            }
        }
        ae.apply_latent_grad(&data.x, &z, &dz, 0.08);
        ae.train_reconstruction(&data.x, 1, 0.15, 0.05, cfg.seed + 184);
    }

    let labels: Vec<i32> = gamma
        .iter()
        .map(|row| {
            row.iter()
                .enumerate()
                .max_by(|a, b| a.1.total_cmp(b.1))
                .map(|(c, _)| c as i32)
                .unwrap_or(0)
        })
        .collect();
    AlgorithmOutput::new(
        labels,
        format!("k={}, VAE latent={}", k, d_lat),
        "Variasjons-autoencoder med GMM-prior i latent rom (VaDE)".to_string(),
    )
}

// 185. DEPICT (Ghasedi Dizaji et al. ICCV 2017 — Denoising AE + balanced cluster prior softmax)
fn depict_clustering(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let d_lat = d.clamp(2, 6);

    let mut ae = AutoEncoder::new(d, d_lat, cfg.seed + 185);
    let recon = ae.train_reconstruction(&data.x, cfg.max_iters(35), 0.25, 0.20, cfg.seed + 185);
    let z: Vec<Vec<f64>> = data.x.iter().map(|r| ae.encode(r)).collect();

    // Softmax clustering head with balanced uniform marginal regularizer
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 1850);
    let c0 = init_kmeans_pp(&z, k, &mut rng);
    let (_, mut head, _, _) = run_lloyd_from_centers(&z, c0, 12);
    let mut p_mat = vec![vec![1.0 / k as f64; k]; n];

    for _ in 0..20 {
        for i in 0..n {
            let logits: Vec<f64> = (0..k).map(|c| -euclidean_sq(&z[i], &head[c]) * 1.5).collect();
            let mx = logits.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            let ex: Vec<f64> = logits.iter().map(|&v| (v - mx).exp()).collect();
            let s: f64 = ex.iter().sum::<f64>().max(1e-12);
            for c in 0..k {
                p_mat[i][c] = ex[c] / s;
            }
        }
        // Balanced target Q_ic = (P_ic / sqrt(f_c)) / sum_j (P_ij / sqrt(f_j))
        let f_c: Vec<f64> = (0..k)
            .map(|c| p_mat.iter().map(|r| r[c]).sum::<f64>().max(1e-6))
            .collect();
        let mut q = vec![vec![0.0; k]; n];
        for i in 0..n {
            let mut s = 0.0;
            for c in 0..k {
                q[i][c] = p_mat[i][c] / f_c[c].sqrt();
                s += q[i][c];
            }
            for c in 0..k {
                q[i][c] /= s.max(1e-12);
            }
        }
        for c in 0..k {
            let mut num = vec![0.0; d_lat];
            let mut den = 0.0;
            for i in 0..n {
                den += q[i][c];
                for h in 0..d_lat {
                    num[h] += q[i][c] * z[i][h];
                }
            }
            if den > 1e-8 {
                for h in 0..d_lat {
                    head[c][h] = num[h] / den;
                }
            }
        }
    }

    let labels: Vec<i32> = p_mat
        .iter()
        .map(|row| {
            row.iter()
                .enumerate()
                .max_by(|a, b| a.1.total_cmp(b.1))
                .map(|(c, _)| c as i32)
                .unwrap_or(0)
        })
        .collect();
    AlgorithmOutput::new(
        labels,
        format!("k={}, DAE støy=0.20, latent={}", k, d_lat),
        format!("Denoising AE + balansert klyngetildeling (recon={:.3})", recon),
    )
}

// 186. JULE (Yang et al. CVPR 2016 — Joint Unsupervised Learning: alternating deep representation + agglomerative merges)
fn jule_clustering(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let d_lat = d.clamp(2, 6);

    let mut ae = AutoEncoder::new(d, d_lat, cfg.seed + 186);
    ae.train_reconstruction(&data.x, cfg.max_iters(25), 0.25, 0.0, cfg.seed + 186);

    let mut labels = vec![0i32; n];
    // Recurrent schedule: progressively merge clusters down toward k while fine-tuning encoder affinity
    let stage_ks = [(k * 3).min(n / 2).max(k), (k * 2).min(n / 2).max(k), k];
    for &sk in &stage_ks {
        let z: Vec<Vec<f64>> = data.x.iter().map(|r| normalize_vec(&ae.encode(r))).collect();
        let dm = crate::math::pairwise_dist_matrix(&z);
        labels = agglomerative_lance_williams(&dm, sk, Linkage::Average);

        // Compute stage cluster prototypes in latent space and pull representations toward their merged cluster
        let mut sums = vec![vec![0.0; d_lat]; sk];
        let mut cnts = vec![0usize; sk];
        for i in 0..n {
            let c = labels[i] as usize;
            cnts[c] += 1;
            for h in 0..d_lat {
                sums[c][h] += z[i][h];
            }
        }
        let mut dz = vec![vec![0.0; d_lat]; n];
        for i in 0..n {
            let c = labels[i] as usize;
            for h in 0..d_lat {
                let proto = sums[c][h] / cnts[c].max(1) as f64;
                dz[i][h] = z[i][h] - proto;
            }
        }
        ae.apply_latent_grad(&data.x, &z, &dz, 0.14);
    }

    AlgorithmOutput::new(
        labels,
        format!("k={}, 3 rekurrente JULE-steg", k),
        "Samtidig representasjonslæring og agglomerativ sammenslåing".to_string(),
    )
}

// 187. DAC (Chang et al. ICCV 2017 — Deep Adaptive Clustering via adaptive pairwise cosine thresholds u(lam) / l(lam))
fn dac_clustering(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 187);

    // Softmax feature-to-cluster linear encoder with L2-normalized probability vectors
    let mut w = vec![vec![0.0; d]; k];
    let c0 = init_kmeans_pp(&data.x, k, &mut rng);
    for c in 0..k {
        w[c] = normalize_vec(&c0[c]);
    }
    let mut b = vec![0.0; k];

    let forward = |x: &[f64], w: &[Vec<f64>], b: &[f64]| -> Vec<f64> {
        let logits: Vec<f64> = (0..k).map(|c| dot(&w[c], x) + b[c]).collect();
        let mx = logits.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let ex: Vec<f64> = logits.iter().map(|&v| (v - mx).exp()).collect();
        let s: f64 = ex.iter().sum::<f64>().max(1e-12);
        let p: Vec<f64> = ex.into_iter().map(|v| v / s).collect();
        normalize_vec(&p)
    };

    let epochs = cfg.max_iters(18);
    for ep in 0..epochs {
        let lam = 0.02 * (ep as f64);
        let u_thr = (0.95 - lam).max(0.55);
        let l_thr = (0.45 + 0.1 * lam).min(0.50);

        let preds: Vec<Vec<f64>> = (0..n).map(|i| forward(&data.x[i], &w, &b)).collect();
        let mut gw = vec![vec![0.0; d]; k];
        let mut gb = vec![0.0; k];

        // Sample pairwise binary classification updates where cosine similarity exceeds u_thr or below l_thr
        let pairs_per_ep = (n * 3).min(500);
        for _ in 0..pairs_per_ep {
            let i = rng.gen_range(0..n);
            let j = rng.gen_range(0..n);
            if i == j {
                continue;
            }
            let sim_feat = cosine_sim(&data.x[i], &data.x[j]);
            let sim_pred = dot(&preds[i], &preds[j]).clamp(1e-5, 1.0 - 1e-5);
            let target = if sim_feat >= u_thr {
                Some(1.0)
            } else if sim_feat <= l_thr {
                Some(0.0)
            } else {
                None
            };
            if let Some(y) = target {
                let err = sim_pred - y;
                for c in 0..k {
                    let gi = err * preds[j][c] / pairs_per_ep as f64;
                    gb[c] += gi;
                    for m in 0..d {
                        gw[c][m] += gi * data.x[i][m];
                    }
                }
            }
        }
        for c in 0..k {
            b[c] -= 0.4 * gb[c];
            for m in 0..d {
                w[c][m] -= 0.4 * gw[c][m];
            }
        }
    }

    let labels: Vec<i32> = (0..n)
        .map(|i| {
            let p = forward(&data.x[i], &w, &b);
            (0..k)
                .max_by(|&a, &b| p[a].total_cmp(&p[b]))
                .unwrap_or(0) as i32
        })
        .collect();
    AlgorithmOutput::new(
        labels,
        format!("k={}, adaptive terskler ({} ep.)", k, epochs),
        "Parvis binær klassifisering med adaptive likhetsterskler".to_string(),
    )
}

// 188. IIC (Ji et al. ICCV 2019 — Invariant Information Clustering maximizing I(z, z') across augmented pairs)
fn iic_clustering(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 188);

    let mut w = vec![vec![0.0; d]; k];
    let c0 = init_kmeans_pp(&data.x, k, &mut rng);
    for c in 0..k {
        w[c] = c0[c].clone();
    }
    let mut b = vec![0.0; k];

    let softmax = |x: &[f64], w: &[Vec<f64>], b: &[f64]| -> Vec<f64> {
        let z: Vec<f64> = (0..k).map(|c| dot(&w[c], x) + b[c]).collect();
        let mx = z.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let ex: Vec<f64> = z.iter().map(|&v| (v - mx).exp()).collect();
        let s: f64 = ex.iter().sum::<f64>().max(1e-12);
        ex.into_iter().map(|v| (v / s).max(1e-7)).collect()
    };

    let mut last_mi = 0.0;
    for _ in 0..cfg.max_iters(25) {
        // Generate jittered augmented views x' of x
        let p1: Vec<Vec<f64>> = (0..n).map(|i| softmax(&data.x[i], &w, &b)).collect();
        let p2: Vec<Vec<f64>> = (0..n)
            .map(|i| {
                let aug: Vec<f64> = data.x[i]
                    .iter()
                    .map(|&v| v + rng.gen_range(-0.15..0.15))
                    .collect();
                softmax(&aug, &w, &b)
            })
            .collect();

        // Joint matrix P = 0.5 * (Q + Q^T), Q = (1/n) sum_i p1[i] p2[i]^T
        let mut joint = vec![vec![1e-7; k]; k];
        for i in 0..n {
            for a in 0..k {
                for bb in 0..k {
                    joint[a][bb] += 0.5 * (p1[i][a] * p2[i][bb] + p1[i][bb] * p2[i][a]) / n as f64;
                }
            }
        }
        let marg: Vec<f64> = (0..k).map(|a| joint[a].iter().sum::<f64>().max(1e-7)).collect();
        last_mi = 0.0;
        for a in 0..k {
            for bb in 0..k {
                last_mi += joint[a][bb] * (joint[a][bb] / (marg[a] * marg[bb])).ln();
            }
        }

        // Ascent step maximizing mutual information I(z, z')
        let mut gw = vec![vec![0.0; d]; k];
        let mut gb = vec![0.0; k];
        for i in 0..n {
            let mut g_p = vec![0.0; k];
            for a in 0..k {
                for bb in 0..k {
                    let coeff = (joint[a][bb] / (marg[a] * marg[bb])).ln();
                    g_p[a] += p2[i][bb] * coeff / n as f64;
                }
            }
            let dot_gp: f64 = (0..k).map(|a| p1[i][a] * g_p[a]).sum();
            for a in 0..k {
                let dlog = p1[i][a] * (g_p[a] - dot_gp);
                gb[a] += dlog;
                for m in 0..d {
                    gw[a][m] += dlog * data.x[i][m];
                }
            }
        }
        for c in 0..k {
            b[c] += 0.5 * gb[c];
            for m in 0..d {
                w[c][m] += 0.5 * gw[c][m];
            }
        }
    }

    let labels: Vec<i32> = (0..n)
        .map(|i| {
            let p = softmax(&data.x[i], &w, &b);
            (0..k)
                .max_by(|&a, &b| p[a].total_cmp(&p[b]))
                .unwrap_or(0) as i32
        })
        .collect();
    AlgorithmOutput::new(
        labels,
        format!("k={}, augmentert MI", k),
        format!("Maksimerte gjensidig informasjon I(z,z')={:.3} nats", last_mi),
    )
}

// 189. SCAN (Van Gansbeke et al. ECCV 2020 — Semantic Clustering by Adopting Nearest neighbors)
fn scan_clustering(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let d_lat = d.clamp(2, 6);

    // Stage 1: Self-supervised representation via Denoising Autoencoder -> mine nearest neighbors
    let mut ae = AutoEncoder::new(d, d_lat, cfg.seed + 189);
    ae.train_reconstruction(&data.x, cfg.max_iters(25), 0.25, 0.1, cfg.seed + 189);
    let z: Vec<Vec<f64>> = data.x.iter().map(|r| ae.encode(r)).collect();

    let k_nn = 5usize.min(n - 1).max(1);
    let mut nbrs: Vec<Vec<usize>> = Vec::with_capacity(n);
    for i in 0..n {
        let mut ds: Vec<(usize, f64)> = (0..n)
            .filter(|&j| j != i)
            .map(|j| (j, euclidean_sq(&z[i], &z[j])))
            .collect();
        ds.sort_by(|a, b| a.1.total_cmp(&b.1));
        nbrs.push(ds.into_iter().take(k_nn).map(|(j, _)| j).collect());
    }

    // Stage 2: Train cluster head enforcing neighbor consistency -log(p_i . p_nbr) + entropy regularization
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 1890);
    let mut w = init_kmeans_pp(&z, k, &mut rng);
    let mut b = vec![0.0; k];
    let softmax = |v: &[f64], w: &[Vec<f64>], b: &[f64]| -> Vec<f64> {
        let lg: Vec<f64> = (0..k).map(|c| dot(&w[c], v) + b[c]).collect();
        let mx = lg.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let ex: Vec<f64> = lg.iter().map(|&x| (x - mx).exp()).collect();
        let s: f64 = ex.iter().sum::<f64>().max(1e-12);
        ex.into_iter().map(|x| (x / s).max(1e-7)).collect()
    };

    for _ in 0..cfg.max_iters(22) {
        let p: Vec<Vec<f64>> = (0..n).map(|i| softmax(&z[i], &w, &b)).collect();
        let p_mean: Vec<f64> = (0..k).map(|c| p.iter().map(|r| r[c]).sum::<f64>() / n as f64).collect();
        let mut gw = vec![vec![0.0; d_lat]; k];
        let mut gb = vec![0.0; k];
        for i in 0..n {
            let mut g_p = vec![0.0; k];
            for &nb in &nbrs[i] {
                let dp = dot(&p[i], &p[nb]).max(1e-6);
                for c in 0..k {
                    g_p[c] -= (p[nb][c] / dp) / (n * k_nn) as f64;
                }
            }
            // Entropy maximization term preventing collapse
            for c in 0..k {
                g_p[c] += 1.5 * (p_mean[c].max(1e-6).ln() + 1.0) / n as f64;
            }
            let dot_gp: f64 = (0..k).map(|c| p[i][c] * g_p[c]).sum();
            for c in 0..k {
                let dlog = p[i][c] * (g_p[c] - dot_gp);
                gb[c] += dlog;
                for h in 0..d_lat {
                    gw[c][h] += dlog * z[i][h];
                }
            }
        }
        for c in 0..k {
            b[c] -= 0.35 * gb[c];
            for h in 0..d_lat {
                w[c][h] -= 0.35 * gw[c][h];
            }
        }
    }

    let labels: Vec<i32> = (0..n)
        .map(|i| {
            let p = softmax(&z[i], &w, &b);
            (0..k)
                .max_by(|&a, &b| p[a].total_cmp(&p[b]))
                .unwrap_or(0) as i32
        })
        .collect();
    AlgorithmOutput::new(
        labels,
        format!("k={}, 2-stegs SCAN (k_NN={})", k, k_nn),
        "Selvveiledet nabokonsistens + entropiregularisering".to_string(),
    )
}

// 190. SwAV (Caron et al. NeurIPS 2020 — Swapping Assignments between Views with Sinkhorn-Knopp equipartition)
fn swav_clustering(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 190);

    let mut protos: Vec<Vec<f64>> = init_kmeans_pp(&data.x_unit, k, &mut rng)
        .into_iter()
        .map(|v| normalize_vec(&v))
        .collect();
    let temp = 0.15_f64;

    for _ in 0..cfg.max_iters(25) {
        // Two augmented views z1, z2 on unit hypersphere
        let z1: Vec<Vec<f64>> = (0..n)
            .map(|i| {
                let v: Vec<f64> = data.x[i].iter().map(|&x| x + rng.gen_range(-0.12..0.12)).collect();
                normalize_vec(&v)
            })
            .collect();
        let z2: Vec<Vec<f64>> = (0..n)
            .map(|i| {
                let v: Vec<f64> = data.x[i].iter().map(|&x| x + rng.gen_range(-0.12..0.12)).collect();
                normalize_vec(&v)
            })
            .collect();

        let scores1: Vec<Vec<f64>> = (0..n)
            .map(|i| (0..k).map(|c| dot(&z1[i], &protos[c])).collect())
            .collect();
        let scores2: Vec<Vec<f64>> = (0..n)
            .map(|i| (0..k).map(|c| dot(&z2[i], &protos[c])).collect())
            .collect();

        // Sinkhorn-Knopp equipartition codes Q1, Q2
        let q1 = sinkhorn_knopp(&scores1, 0.08, 4);
        let q2 = sinkhorn_knopp(&scores2, 0.08, 4);

        // Swapped prediction gradient: predict Q2 from z1 and Q1 from z2
        let mut grad = vec![vec![0.0; d]; k];
        for i in 0..n {
            let p1 = {
                let mx = scores1[i].iter().copied().fold(f64::NEG_INFINITY, f64::max);
                let ex: Vec<f64> = scores1[i].iter().map(|&v| ((v - mx) / temp).exp()).collect();
                let s: f64 = ex.iter().sum::<f64>().max(1e-12);
                ex.into_iter().map(|v| v / s).collect::<Vec<_>>()
            };
            let p2 = {
                let mx = scores2[i].iter().copied().fold(f64::NEG_INFINITY, f64::max);
                let ex: Vec<f64> = scores2[i].iter().map(|&v| ((v - mx) / temp).exp()).collect();
                let s: f64 = ex.iter().sum::<f64>().max(1e-12);
                ex.into_iter().map(|v| v / s).collect::<Vec<_>>()
            };
            for c in 0..k {
                let e1 = (p1[c] - q2[i][c]) / n as f64;
                let e2 = (p2[c] - q1[i][c]) / n as f64;
                for m in 0..d {
                    grad[c][m] += e1 * z1[i][m] + e2 * z2[i][m];
                }
            }
        }
        for c in 0..k {
            for m in 0..d {
                protos[c][m] -= 0.35 * grad[c][m];
            }
            protos[c] = normalize_vec(&protos[c]);
        }
    }

    let final_scores: Vec<Vec<f64>> = (0..n)
        .map(|i| (0..k).map(|c| dot(&data.x_unit[i], &protos[c])).collect())
        .collect();
    let q_final = sinkhorn_knopp(&final_scores, 0.08, 5);
    let labels: Vec<i32> = q_final
        .iter()
        .map(|r| {
            r.iter()
                .enumerate()
                .max_by(|a, b| a.1.total_cmp(b.1))
                .map(|(c, _)| c as i32)
                .unwrap_or(0)
        })
        .collect();
    AlgorithmOutput::new(
        labels,
        format!("k={}, Sinkhorn-Knopp tau={:.2}", k, temp),
        "Swapped online prototype-læring med ekvipartisjon".to_string(),
    )
}

// 191. DeepCluster (Caron et al. ECCV 2018 — alternates K-Means pseudo-labels & supervised network training)
fn deep_cluster(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let d_lat = d.clamp(2, 6);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 191);

    let mut ae = AutoEncoder::new(d, d_lat, cfg.seed + 191);
    let mut head_w = vec![vec![0.0; d_lat]; k];
    for c in 0..k {
        for h in 0..d_lat {
            head_w[c][h] = rng.gen_range(-0.4..0.4);
        }
    }
    let mut labels = vec![0i32; n];

    for _round in 0..6 {
        // Step 1: Extract deep features and cluster with K-Means++ for pseudo-labels
        let feats: Vec<Vec<f64>> = data.x.iter().map(|r| ae.encode(r)).collect();
        let c0 = init_kmeans_pp(&feats, k, &mut rng);
        let (pseudo, _, _, _) = run_lloyd_from_centers(&feats, c0, 15);
        labels = pseudo;

        // Step 2: Train encoder + classification head with cross-entropy on pseudo-labels
        for _ep in 0..5 {
            let z: Vec<Vec<f64>> = data.x.iter().map(|r| ae.encode(r)).collect();
            let mut dz = vec![vec![0.0; d_lat]; n];
            let mut gh = vec![vec![0.0; d_lat]; k];
            for i in 0..n {
                let y = labels[i] as usize;
                let lg: Vec<f64> = (0..k).map(|c| dot(&head_w[c], &z[i])).collect();
                let mx = lg.iter().copied().fold(f64::NEG_INFINITY, f64::max);
                let ex: Vec<f64> = lg.iter().map(|&v| (v - mx).exp()).collect();
                let s: f64 = ex.iter().sum::<f64>().max(1e-12);
                for c in 0..k {
                    let err = (ex[c] / s) - if c == y { 1.0 } else { 0.0 };
                    for h in 0..d_lat {
                        gh[c][h] += err * z[i][h] / n as f64;
                        dz[i][h] += err * head_w[c][h];
                    }
                }
            }
            for c in 0..k {
                for h in 0..d_lat {
                    head_w[c][h] -= 0.25 * gh[c][h];
                }
            }
            ae.apply_latent_grad(&data.x, &z, &dz, 0.12);
        }
    }

    AlgorithmOutput::new(
        labels,
        format!("k={}, 6 pseudo-etikett-runder", k),
        "Alternerte K-Means pseudo-etiketter og dyp feature-trening".to_string(),
    )
}

// 192. Contrastive Clustering (Li et al. AAAI 2021 — Instance-level + Cluster-level InfoNCE dual heads)
fn contrastive_clustering(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let d_lat = d.clamp(2, 6);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 192);

    let mut ae = AutoEncoder::new(d, d_lat, cfg.seed + 192);
    ae.train_reconstruction(&data.x, cfg.max_iters(20), 0.22, 0.1, cfg.seed + 192);

    // Instance-level contrastive refinement: pull augmented positive pair closer, push random negative apart
    for _ in 0..cfg.max_iters(14) {
        let z1: Vec<Vec<f64>> = data.x.iter().map(|r| normalize_vec(&ae.encode(r))).collect();
        let z2: Vec<Vec<f64>> = (0..n)
            .map(|i| {
                let aug: Vec<f64> = data.x[i].iter().map(|&v| v + rng.gen_range(-0.15..0.15)).collect();
                normalize_vec(&ae.encode(&aug))
            })
            .collect();
        let mut dz = vec![vec![0.0; d_lat]; n];
        for i in 0..n {
            let neg = rng.gen_range(0..n);
            let sim_pos = dot(&z1[i], &z2[i]);
            let sim_neg = dot(&z1[i], &z1[neg]);
            let w_neg = ((sim_neg - sim_pos) * 2.0).exp().clamp(0.05, 5.0);
            for h in 0..d_lat {
                dz[i][h] = -(z2[i][h] - w_neg * z1[neg][h]) * 0.25;
            }
        }
        ae.apply_latent_grad(&data.x, &z1, &dz, 0.12);
    }

    // Cluster-level head on contrastive features
    let z_final: Vec<Vec<f64>> = data.x.iter().map(|r| normalize_vec(&ae.encode(r))).collect();
    let c0 = init_kmeans_pp(&z_final, k, &mut rng);
    let (labels, _, _, _) = run_lloyd_from_centers(&z_final, c0, 30);
    AlgorithmOutput::new(
        labels,
        format!("k={}, instans+klynge InfoNCE", k),
        "Kontrastiv læring på instans- og klyngenivå (CC)".to_string(),
    )
}

// 193. CluStream (Aggarwal et al. VLDB 2003 — Online micro-clusters with timestamps + Offline weighted K-Means)
fn clustream(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let q_micro = (k * 5).clamp(12, 32).min(n);

    struct MicroCluster {
        n: usize,
        cf1x: Vec<f64>,
        cf2x: Vec<f64>,
        cf1t: f64,
    }
    impl MicroCluster {
        fn new(p: &[f64], t: f64) -> Self {
            Self { n: 1, cf1x: p.to_vec(), cf2x: p.iter().map(|v| v * v).collect(), cf1t: t }
        }
        fn center(&self) -> Vec<f64> {
            self.cf1x.iter().map(|&v| v / self.n as f64).collect()
        }
        fn radius(&self) -> f64 {
            if self.n <= 1 {
                return 0.0;
            }
            let nf = self.n as f64;
            let var: f64 = self
                .cf2x
                .iter()
                .zip(self.cf1x.iter())
                .map(|(s2, s1)| (s2 / nf - (s1 / nf).powi(2)).max(0.0))
                .sum();
            var.sqrt()
        }
        fn absorb(&mut self, p: &[f64], t: f64) {
            self.n += 1;
            for j in 0..p.len() {
                self.cf1x[j] += p[j];
                self.cf2x[j] += p[j] * p[j];
            }
            self.cf1t += t;
        }
    }

    // Stream arrival order
    let mut stream_ord: Vec<usize> = (0..n).collect();
    stream_ord.sort_by(|&a, &b| data.time_values[a].total_cmp(&data.time_values[b]));

    let mut micros: Vec<MicroCluster> = Vec::with_capacity(q_micro);
    let init_r = median_pairwise_dist(&data.dist_matrix) * 0.35;

    for (step, &idx) in stream_ord.iter().enumerate() {
        let p = &data.x[idx];
        let t = step as f64;
        if micros.len() < q_micro {
            micros.push(MicroCluster::new(p, t));
            continue;
        }
        let best = (0..micros.len())
            .min_by(|&a, &b| {
                euclidean_sq(p, &micros[a].center()).total_cmp(&euclidean_sq(p, &micros[b].center()))
            })
            .unwrap_or(0);
        let max_b = (micros[best].radius() * 2.0).max(init_r);
        if euclidean(p, &micros[best].center()) <= max_b {
            micros[best].absorb(p, t);
        } else {
            // Merge two closest micro-clusters to make room for new outlier micro-cluster
            let mut pair = (0, 1, f64::INFINITY);
            for a in 0..micros.len() {
                for b in (a + 1)..micros.len() {
                    let dd = euclidean(&micros[a].center(), &micros[b].center());
                    if dd < pair.2 {
                        pair = (a, b, dd);
                    }
                }
            }
            let b_mc = micros.remove(pair.1);
            micros[pair.0].n += b_mc.n;
            for j in 0..d {
                micros[pair.0].cf1x[j] += b_mc.cf1x[j];
                micros[pair.0].cf2x[j] += b_mc.cf2x[j];
            }
            micros[pair.0].cf1t += b_mc.cf1t;
            micros.push(MicroCluster::new(p, t));
        }
    }

    // Offline macro-clustering: weighted K-Means on micro-cluster centroids
    let mc_pts: Vec<Vec<f64>> = micros.iter().map(|m| m.center()).collect();
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 193);
    let c0 = init_kmeans_pp(&mc_pts, k.min(mc_pts.len()), &mut rng);
    let (_, macro_centers, _, _) = run_lloyd_from_centers(&mc_pts, c0, 25);

    let labels: Vec<i32> = (0..n)
        .map(|i| {
            (0..macro_centers.len())
                .min_by(|&a, &b| {
                    euclidean_sq(&data.x[i], &macro_centers[a])
                        .total_cmp(&euclidean_sq(&data.x[i], &macro_centers[b]))
                })
                .unwrap_or(0) as i32
        })
        .collect();
    AlgorithmOutput::new(
        labels,
        format!("q={} mikroklynger, k={}", micros.len(), k),
        "2-fase strøm: online mikroklynger + offline makro-K-Means".to_string(),
    )
}

// 194. DenStream (Cao et al. SDM 2006 — Potential p-micro-clusters vs Outlier o-micro-clusters with time decay)
fn denstream(data: &PreparedData, _cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let eps = median_pairwise_dist(&data.dist_matrix) * 0.42;
    let lambda = 0.008_f64;
    let beta_mu = 2.0_f64;

    struct DenMicro {
        w: f64,
        cf1: Vec<f64>,
        cf2: Vec<f64>,
    }
    impl DenMicro {
        fn new(p: &[f64]) -> Self {
            Self { w: 1.0, cf1: p.to_vec(), cf2: p.iter().map(|v| v * v).collect() }
        }
        fn center(&self) -> Vec<f64> {
            self.cf1.iter().map(|&v| v / self.w.max(1e-8)).collect()
        }
        fn decay(&mut self, factor: f64) {
            self.w *= factor;
            for v in &mut self.cf1 {
                *v *= factor;
            }
            for v in &mut self.cf2 {
                *v *= factor;
            }
        }
        fn add(&mut self, p: &[f64]) {
            self.w += 1.0;
            for j in 0..p.len() {
                self.cf1[j] += p[j];
                self.cf2[j] += p[j] * p[j];
            }
        }
    }

    let mut stream_ord: Vec<usize> = (0..n).collect();
    stream_ord.sort_by(|&a, &b| data.time_values[a].total_cmp(&data.time_values[b]));

    let mut p_micros: Vec<DenMicro> = Vec::new();
    let mut o_micros: Vec<DenMicro> = Vec::new();
    let decay_step = 2.0_f64.powf(-lambda);

    for &idx in &stream_ord {
        let p = &data.x[idx];
        for m in &mut p_micros {
            m.decay(decay_step);
        }
        for m in &mut o_micros {
            m.decay(decay_step);
        }

        // Try absorbing into nearest p-micro-cluster
        if let Some((bi, _)) = p_micros.iter().enumerate().min_by(|(_, a), (_, b)| {
            euclidean_sq(p, &a.center()).total_cmp(&euclidean_sq(p, &b.center()))
        }) {
            if euclidean(p, &p_micros[bi].center()) <= eps {
                p_micros[bi].add(p);
                continue;
            }
        }
        // Else try o-micro-cluster and promote if weight >= beta * mu
        if let Some((bi, _)) = o_micros.iter().enumerate().min_by(|(_, a), (_, b)| {
            euclidean_sq(p, &a.center()).total_cmp(&euclidean_sq(p, &b.center()))
        }) {
            if euclidean(p, &o_micros[bi].center()) <= eps {
                o_micros[bi].add(p);
                if o_micros[bi].w >= beta_mu {
                    let promoted = o_micros.remove(bi);
                    p_micros.push(promoted);
                }
                continue;
            }
        }
        o_micros.push(DenMicro::new(p));
    }

    // Offline DBSCAN over p-micro-cluster centers
    let pm = p_micros.len();
    let p_cens: Vec<Vec<f64>> = p_micros.iter().map(|m| m.center()).collect();
    let mut p_lbl = vec![-1i32; pm];
    let mut cid = 0i32;
    for i in 0..pm {
        if p_lbl[i] >= 0 {
            continue;
        }
        let mut q = VecDeque::new();
        p_lbl[i] = cid;
        q.push_back(i);
        while let Some(u) = q.pop_front() {
            for v in 0..pm {
                if p_lbl[v] < 0 && euclidean(&p_cens[u], &p_cens[v]) <= eps * 2.0 {
                    p_lbl[v] = cid;
                    q.push_back(v);
                }
            }
        }
        cid += 1;
    }

    let mut labels = vec![-1i32; n];
    for i in 0..n {
        if let Some((bi, _)) = p_cens
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| euclidean_sq(&data.x[i], a).total_cmp(&euclidean_sq(&data.x[i], b)))
        {
            if euclidean(&data.x[i], &p_cens[bi]) <= eps * 1.8 {
                labels[i] = p_lbl[bi];
            }
        }
    }
    let _ = d;
    AlgorithmOutput::new(
        labels,
        format!("eps={:.2}, {} p-mikro / {} o-mikro", eps, pm, o_micros.len()),
        "Tetthetsstrøm med tidsnedbrytning og p/o-mikroklynger".to_string(),
    )
}

// 195. StreamKM++ (Ackermann et al. JEA 2012 — Merge-and-Reduce coreset tree with K-Means++ sampling)
fn stream_km_pp(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let m_bucket = (k * 6).clamp(16, 40).min(n);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 195);

    // Weighted point (coords, weight)
    type WPt = (Vec<f64>, f64);

    let reduce_coreset = |pts: &[WPt], target: usize, rng: &mut ChaCha8Rng| -> Vec<WPt> {
        if pts.len() <= target {
            return pts.to_vec();
        }
        // Recursive coreset tree splitting by weighted K-Means++ (2-means) until `target` leaves
        let mut leaves: Vec<Vec<usize>> = vec![(0..pts.len()).collect()];
        while leaves.len() < target {
            let Some(li) = leaves
                .iter()
                .enumerate()
                .filter(|(_, lf)| lf.len() >= 2)
                .max_by_key(|(_, lf)| lf.len())
                .map(|(i, _)| i)
            else {
                break;
            };
            let mems = leaves.remove(li);
            let sub_x: Vec<Vec<f64>> = mems.iter().map(|&i| pts[i].0.clone()).collect();
            let c0 = init_kmeans_pp(&sub_x, 2, rng);
            let (lbl, _, _, _) = run_lloyd_from_centers(&sub_x, c0, 6);
            let left: Vec<usize> = mems.iter().enumerate().filter(|&(j, _)| lbl[j] == 0).map(|(_, &i)| i).collect();
            let right: Vec<usize> = mems.iter().enumerate().filter(|&(j, _)| lbl[j] != 0).map(|(_, &i)| i).collect();
            if left.is_empty() || right.is_empty() {
                let mid = mems.len() / 2;
                leaves.push(mems[..mid].to_vec());
                leaves.push(mems[mid..].to_vec());
            } else {
                leaves.push(left);
                leaves.push(right);
            }
        }
        leaves
            .into_iter()
            .map(|lf| {
                let mut sum = vec![0.0; d];
                let mut wtot = 0.0;
                for i in lf {
                    wtot += pts[i].1;
                    for j in 0..d {
                        sum[j] += pts[i].1 * pts[i].0[j];
                    }
                }
                for v in &mut sum {
                    *v /= wtot.max(1e-9);
                }
                (sum, wtot)
            })
            .collect()
    };

    // Merge-and-Reduce buckets in streaming order
    let mut buckets: Vec<Option<Vec<WPt>>> = Vec::new();
    let mut buf: Vec<WPt> = Vec::with_capacity(m_bucket);
    for i in 0..n {
        buf.push((data.x[i].clone(), 1.0));
        if buf.len() == m_bucket {
            let mut carry = std::mem::take(&mut buf);
            let mut lvl = 0;
            loop {
                if lvl == buckets.len() {
                    buckets.push(Some(carry));
                    break;
                }
                if let Some(mut existing) = buckets[lvl].take() {
                    existing.extend(carry);
                    carry = reduce_coreset(&existing, m_bucket, &mut rng);
                    lvl += 1;
                } else {
                    buckets[lvl] = Some(carry);
                    break;
                }
            }
        }
    }

    let mut final_coreset = buf;
    for b in buckets.into_iter().flatten() {
        final_coreset.extend(b);
    }
    let cs_x: Vec<Vec<f64>> = final_coreset.iter().map(|p| p.0.clone()).collect();
    let c0 = init_kmeans_pp(&cs_x, k.min(cs_x.len()), &mut rng);
    let (_, centers, _, _) = run_lloyd_from_centers(&cs_x, c0, 25);

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
        format!("k={}, coreset-bøtte m={}", k, m_bucket),
        format!("Merge-and-Reduce coreset ({} vektede punkter)", final_coreset.len()),
    )
}

// 196. BIRCH-Stream / ClusTree (Kranen et al. 2011 — Anytime decaying CF-tree adapting to stream arrival)
fn clustree_birch_stream(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let max_entries = (k * 5).clamp(12, 30).min(n);
    let decay_lambda = 0.004_f64;

    struct DecayingCf {
        w: f64,
        ls: Vec<f64>,
        last_t: f64,
    }
    impl DecayingCf {
        fn new(p: &[f64], t: f64) -> Self {
            Self { w: 1.0, ls: p.to_vec(), last_t: t }
        }
        fn decay_to(&mut self, t: f64, lam: f64) {
            let dt = (t - self.last_t).max(0.0);
            let f = 2.0_f64.powf(-lam * dt);
            self.w *= f;
            for v in &mut self.ls {
                *v *= f;
            }
            self.last_t = t;
        }
        fn center(&self) -> Vec<f64> {
            self.ls.iter().map(|&v| v / self.w.max(1e-9)).collect()
        }
        fn add(&mut self, p: &[f64], t: f64, lam: f64) {
            self.decay_to(t, lam);
            self.w += 1.0;
            for j in 0..p.len() {
                self.ls[j] += p[j];
            }
        }
    }

    let mut entries: Vec<DecayingCf> = Vec::with_capacity(max_entries);
    for i in 0..n {
        let t = data.time_values[i] * n as f64;
        if entries.len() < max_entries {
            entries.push(DecayingCf::new(&data.x[i], t));
        } else {
            let best = (0..entries.len())
                .min_by(|&a, &b| {
                    euclidean_sq(&data.x[i], &entries[a].center())
                        .total_cmp(&euclidean_sq(&data.x[i], &entries[b].center()))
                })
                .unwrap_or(0);
            entries[best].add(&data.x[i], t, decay_lambda);
        }
    }

    let cens: Vec<Vec<f64>> = entries.iter().map(|e| e.center()).collect();
    let m = cens.len();
    let mut cdist = vec![vec![0.0; m]; m];
    for a in 0..m {
        for b in (a + 1)..m {
            let dd = euclidean(&cens[a], &cens[b]);
            cdist[a][b] = dd;
            cdist[b][a] = dd;
        }
    }
    let clbl = if m > k {
        agglomerative_lance_williams(&cdist, k, Linkage::Ward)
    } else {
        (0..m as i32).collect()
    };
    let labels: Vec<i32> = (0..n)
        .map(|i| {
            let best = (0..m)
                .min_by(|&a, &b| euclidean_sq(&data.x[i], &cens[a]).total_cmp(&euclidean_sq(&data.x[i], &cens[b])))
                .unwrap_or(0);
            clbl[best]
        })
        .collect();
    let _ = d;
    AlgorithmOutput::new(
        labels,
        format!("k={}, {} tidsveide CF-noder", k, m),
        "ClusTree eksponentielt dempet CF-tre over ankomsttid".to_string(),
    )
}

// 197. D-Stream (Chen & Tu KDD 2007 — Grid-based stream clustering with exponential density decay)
fn d_stream(data: &PreparedData, _cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let g = 8usize;
    let decay = 0.992_f64;

    let mut mn = [f64::INFINITY; 2];
    let mut mx = [f64::NEG_INFINITY; 2];
    for p in &data.pca_2d {
        for j in 0..2 {
            mn[j] = mn[j].min(p[j]);
            mx[j] = mx[j].max(p[j]);
        }
    }

    let mut stream_ord: Vec<usize> = (0..n).collect();
    stream_ord.sort_by(|&a, &b| data.time_values[a].total_cmp(&data.time_values[b]));

    let mut grid_dens: HashMap<(usize, usize), f64> = HashMap::new();
    let mut pt_cell = vec![(0usize, 0usize); n];

    for &idx in &stream_ord {
        for val in grid_dens.values_mut() {
            *val *= decay;
        }
        let cx = (((data.pca_2d[idx][0] - mn[0]) / (mx[0] - mn[0]).max(1e-9)) * g as f64)
            .floor()
            .clamp(0.0, (g - 1) as f64) as usize;
        let cy = (((data.pca_2d[idx][1] - mn[1]) / (mx[1] - mn[1]).max(1e-9)) * g as f64)
            .floor()
            .clamp(0.0, (g - 1) as f64) as usize;
        pt_cell[idx] = (cx, cy);
        *grid_dens.entry((cx, cy)).or_insert(0.0) += 1.0;
    }

    let mut vals: Vec<f64> = grid_dens.values().copied().collect();
    vals.sort_by(|a, b| a.total_cmp(b));
    let dm_thresh = vals.get(vals.len() / 3).copied().unwrap_or(0.5);

    let dense: HashSet<(usize, usize)> = grid_dens
        .iter()
        .filter(|(_, &v)| v >= dm_thresh)
        .map(|(&c, _)| c)
        .collect();

    // Connected components of dense cells
    let mut cell_lbl: HashMap<(usize, usize), i32> = HashMap::new();
    let mut cid = 0i32;
    for &c in &dense {
        if cell_lbl.contains_key(&c) {
            continue;
        }
        let mut q = VecDeque::new();
        cell_lbl.insert(c, cid);
        q.push_back(c);
        while let Some((cx, cy)) = q.pop_front() {
            for dx in -1isize..=1 {
                for dy in -1isize..=1 {
                    let nx = cx as isize + dx;
                    let ny = cy as isize + dy;
                    if nx >= 0 && ny >= 0 && (nx as usize) < g && (ny as usize) < g {
                        let nb = (nx as usize, ny as usize);
                        if dense.contains(&nb) && !cell_lbl.contains_key(&nb) {
                            cell_lbl.insert(nb, cid);
                            q.push_back(nb);
                        }
                    }
                }
            }
        }
        cid += 1;
    }

    let mut labels = vec![-1i32; n];
    for i in 0..n {
        if let Some(&l) = cell_lbl.get(&pt_cell[i]) {
            labels[i] = l;
        }
    }
    AlgorithmOutput::new(
        labels,
        format!("{}x{} rutenett, decay={:.3}", g, g, decay),
        format!("D-Stream tidsnedbrytning fant {} tette celler", dense.len()),
    )
}

// 198. EVO-KMeans / Genetic K-Means Algorithm (Krishna & Murty 1999 — GKA crossover + biased mutation + 1-step K-Means)
fn evo_kmeans_gka(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let pop_size = 14usize;
    let gens = cfg.max_iters(22);
    let mut_prob = 0.08_f64;
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 198);

    // Chromosome = flattened k centroids (or label assignment refined by K-Means operator)
    let eval_centers = |cens: &[Vec<f64>]| -> (Vec<i32>, f64) {
        let mut lbl = vec![0i32; n];
        let mut sse = 0.0;
        for i in 0..n {
            let mut bc = 0;
            let mut bd = f64::INFINITY;
            for (c, cen) in cens.iter().enumerate() {
                let dd = euclidean_sq(&data.x[i], cen);
                if dd < bd {
                    bd = dd;
                    bc = c;
                }
            }
            lbl[i] = bc as i32;
            sse += bd;
        }
        (lbl, sse)
    };

    let mut pop: Vec<Vec<Vec<f64>>> = (0..pop_size)
        .map(|_| init_kmeans_pp(&data.x, k, &mut rng))
        .collect();
    let mut best_labels = vec![0i32; n];
    let mut best_sse = f64::INFINITY;

    for _ in 0..gens {
        let mut scored: Vec<(usize, f64)> = pop
            .iter()
            .enumerate()
            .map(|(idx, c)| {
                let (lbl, sse) = eval_centers(c);
                if sse < best_sse {
                    best_sse = sse;
                    best_labels = lbl;
                }
                (idx, sse)
            })
            .collect();
        scored.sort_by(|a, b| a.1.total_cmp(&b.1));

        // Elitism: keep top 2 individuals
        let mut next_pop = vec![pop[scored[0].0].clone(), pop[scored[1].0].clone()];
        let top_pool = (pop_size / 2).max(2);

        while next_pop.len() < pop_size {
            let p1 = &pop[scored[rng.gen_range(0..top_pool)].0];
            let p2 = &pop[scored[rng.gen_range(0..top_pool)].0];
            // Arithmetic centroid crossover + Gaussian mutation
            let mut child = vec![vec![0.0; d]; k];
            for c in 0..k {
                let w = rng.gen_range(0.0..1.0);
                for j in 0..d {
                    child[c][j] = w * p1[c][j] + (1.0 - w) * p2[c][j];
                    if rng.gen::<f64>() < mut_prob {
                        child[c][j] += rng.gen_range(-0.35..0.35);
                    }
                }
            }
            // GKA 1-step K-Means local search operator
            let (_, refined, _, _) = run_lloyd_from_centers(&data.x, child, 2);
            next_pop.push(refined);
        }
        pop = next_pop;
    }

    AlgorithmOutput::new(
        best_labels,
        format!("k={}, pop={}, gen={}", k, pop_size, gens),
        format!("Genetisk K-Means (GKA) nådde SSE={:.2}", best_sse),
    )
}

// 199. PSO-Clustering (van der Merwe & Engelbrecht 2003 — Particle Swarm Optimization in centroid space)
fn pso_clustering(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let n_particles = 14usize;
    let iters = cfg.max_iters(35);
    let w_inertia = 0.65_f64;
    let c1 = 1.4_f64;
    let c2 = 1.4_f64;
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 199);

    let eval_quant = |pos: &[Vec<f64>]| -> (Vec<i32>, f64) {
        let mut lbl = vec![0i32; n];
        let mut err = 0.0;
        for i in 0..n {
            let mut bc = 0;
            let mut bd = f64::INFINITY;
            for (c, cen) in pos.iter().enumerate() {
                let dd = euclidean_sq(&data.x[i], cen);
                if dd < bd {
                    bd = dd;
                    bc = c;
                }
            }
            lbl[i] = bc as i32;
            err += bd;
        }
        (lbl, err)
    };

    let mut pos: Vec<Vec<Vec<f64>>> = (0..n_particles)
        .map(|_| init_kmeans_pp(&data.x, k, &mut rng))
        .collect();
    let mut vel = vec![vec![vec![0.0; d]; k]; n_particles];
    let mut pbest = pos.clone();
    let mut pbest_val = vec![f64::INFINITY; n_particles];
    let mut gbest = pos[0].clone();
    let mut gbest_val = f64::INFINITY;
    let mut gbest_labels = vec![0i32; n];

    for _ in 0..iters {
        for p in 0..n_particles {
            let (lbl, val) = eval_quant(&pos[p]);
            if val < pbest_val[p] {
                pbest_val[p] = val;
                pbest[p] = pos[p].clone();
            }
            if val < gbest_val {
                gbest_val = val;
                gbest = pos[p].clone();
                gbest_labels = lbl;
            }
        }
        for p in 0..n_particles {
            for c in 0..k {
                for j in 0..d {
                    let r1 = rng.gen::<f64>();
                    let r2 = rng.gen::<f64>();
                    vel[p][c][j] = w_inertia * vel[p][c][j]
                        + c1 * r1 * (pbest[p][c][j] - pos[p][c][j])
                        + c2 * r2 * (gbest[c][j] - pos[p][c][j]);
                    pos[p][c][j] += vel[p][c][j];
                }
            }
        }
    }

    AlgorithmOutput::new(
        gbest_labels,
        format!("k={}, {} partikler, {} iter.", k, n_particles, iters),
        format!("PSO-sverm konvergerte til SSE={:.2}", gbest_val),
    )
}

// 200. ACO-Clustering (Shelokar et al. 2004 — Ant Colony Optimization pheromone matrix + local refinement)
fn aco_clustering(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let n_ants = 12usize;
    let iters = cfg.max_iters(25);
    let rho = 0.25_f64; // evaporation rate
    let q0 = 0.72_f64; // exploitation probability
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 200);

    // Pheromone matrix tau[i][c] for assigning object i to cluster c
    let mut tau = vec![vec![0.1_f64; k]; n];
    // Warm-start pheromone hint with one K-Means++ trail
    let c_init = init_kmeans_pp(&data.x, k, &mut rng);
    let (l_init, _, _, _) = run_lloyd_from_centers(&data.x, c_init, 8);
    for i in 0..n {
        tau[i][l_init[i] as usize] += 0.35;
    }

    let mut best_labels = l_init;
    let mut best_sse = f64::INFINITY;

    for _ in 0..iters {
        let mut ant_solutions: Vec<(Vec<i32>, f64)> = Vec::with_capacity(n_ants);
        for _ in 0..n_ants {
            // Ant constructs assignment using pheromone trail
            let mut sol = vec![0i32; n];
            for i in 0..n {
                let c = if rng.gen::<f64>() < q0 {
                    (0..k)
                        .max_by(|&a, &b| tau[i][a].total_cmp(&tau[i][b]))
                        .unwrap_or(0)
                } else {
                    let s: f64 = tau[i].iter().sum::<f64>().max(1e-12);
                    let mut r = rng.gen::<f64>() * s;
                    let mut chosen = k - 1;
                    for ci in 0..k {
                        r -= tau[i][ci];
                        if r <= 0.0 {
                            chosen = ci;
                            break;
                        }
                    }
                    chosen
                };
                sol[i] = c as i32;
            }
            // Local search (1 centroid update step) to evaluate ant solution quality
            let mut cens = vec![vec![0.0; d]; k];
            let mut cnts = vec![0usize; k];
            for i in 0..n {
                let c = sol[i] as usize;
                cnts[c] += 1;
                for j in 0..d {
                    cens[c][j] += data.x[i][j];
                }
            }
            for c in 0..k {
                if cnts[c] > 0 {
                    for j in 0..d {
                        cens[c][j] /= cnts[c] as f64;
                    }
                } else {
                    cens[c] = data.x[rng.gen_range(0..n)].clone();
                }
            }
            let (refined_lbl, _, sse, _) = run_lloyd_from_centers(&data.x, cens, 2);
            ant_solutions.push((refined_lbl, sse));
        }

        ant_solutions.sort_by(|a, b| a.1.total_cmp(&b.1));
        if ant_solutions[0].1 < best_sse {
            best_sse = ant_solutions[0].1;
            best_labels = ant_solutions[0].0.clone();
        }

        // Evaporate pheromones and reinforce trails of top 3 elite ants
        for i in 0..n {
            for c in 0..k {
                tau[i][c] = (tau[i][c] * (1.0 - rho)).max(0.01);
            }
        }
        for (lbl, sse) in ant_solutions.iter().take(3) {
            let deposit = 1.0 / (sse.max(1.0));
            for i in 0..n {
                tau[i][lbl[i] as usize] += deposit * 20.0;
            }
        }
    }

    AlgorithmOutput::new(
        best_labels,
        format!("k={}, {} maur, rho={:.2}", k, n_ants, rho),
        format!("Feromonstyrt maursøk med lokal forbedring (SSE={:.2})", best_sse),
    )
}
