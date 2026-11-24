use std::collections::HashMap;

use rand::seq::SliceRandom;
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;

use crate::algorithms::{AlgorithmOutput, RunConfig};
use crate::data::PreparedData;
use crate::math::{
    anderson_darling_normal, dot, euclidean, euclidean_sq, hamming_dist, jacobi_eigen, manhattan,
    median_of_slice, median_pairwise_dist, normalize_vec, weighted_minkowski,
};

pub fn run(id: usize, data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    match id {
        1 => kmeans_lloyd(data, cfg),
        2 => kmeans_plus_plus(data, cfg),
        3 => mini_batch_kmeans(data, cfg),
        4 => kmedoids_pam(data, cfg),
        5 => clara(data, cfg),
        6 => clarans(data, cfg),
        7 => kmodes(data, cfg),
        8 => kprototypes(data, cfg),
        9 => bisecting_kmeans(data, cfg),
        10 => xmeans_bic(data, cfg),
        11 => gmeans_anderson_darling(data, cfg),
        12 => spherical_kmeans(data, cfg),
        13 => elkan_kmeans(data, cfg),
        14 => hamerly_kmeans(data, cfg),
        15 => kmedians(data, cfg),
        16 => kernel_kmeans(data, cfg),
        17 => global_kmeans(data, cfg),
        18 => isodata(data, cfg),
        19 => hartigan_wong(data, cfg),
        20 => cop_kmeans(data, cfg),
        21 => macqueen_kmeans(data, cfg),
        22 => forgy_kmeans(data, cfg),
        23 => kmeans_parallel(data, cfg),
        24 => faster_pam(data, cfg),
        25 => clara_med(data, cfg),
        26 => trimmed_kmeans(data, cfg),
        27 => power_kmeans(data, cfg),
        28 => weighted_kmeans(data, cfg),
        29 => minkowski_weighted_kmeans(data, cfg),
        30 => lbg_vector_quantization(data, cfg),
        _ => unreachable!(),
    }
}

pub fn init_kmeans_pp(x: &[Vec<f64>], k: usize, rng: &mut ChaCha8Rng) -> Vec<Vec<f64>> {
    let n = x.len();
    let mut centers = Vec::with_capacity(k);
    centers.push(x[rng.gen_range(0..n)].clone());
    for _ in 1..k {
        let dists: Vec<f64> = x
            .iter()
            .map(|p| {
                centers
                    .iter()
                    .map(|c| euclidean_sq(p, c))
                    .fold(f64::INFINITY, f64::min)
            })
            .collect();
        let sum: f64 = dists.iter().sum();
        if sum <= 1e-12 {
            centers.push(x[rng.gen_range(0..n)].clone());
        } else {
            let mut r = rng.gen::<f64>() * sum;
            let mut chosen = n - 1;
            for (i, &w) in dists.iter().enumerate() {
                r -= w;
                if r <= 0.0 {
                    chosen = i;
                    break;
                }
            }
            centers.push(x[chosen].clone());
        }
    }
    centers
}

pub fn run_lloyd_from_centers(
    x: &[Vec<f64>],
    mut centers: Vec<Vec<f64>>,
    max_iter: usize,
) -> (Vec<i32>, Vec<Vec<f64>>, f64, usize) {
    let n = x.len();
    let d = x[0].len();
    let k = centers.len();
    let mut labels = vec![0i32; n];
    let mut iters = 0;

    for it in 0..max_iter {
        iters = it + 1;
        let mut changed = false;
        for i in 0..n {
            let mut best_c = 0;
            let mut best_d = f64::INFINITY;
            for (c, cen) in centers.iter().enumerate() {
                let dist = euclidean_sq(&x[i], cen);
                if dist < best_d {
                    best_d = dist;
                    best_c = c;
                }
            }
            if labels[i] != best_c as i32 {
                labels[i] = best_c as i32;
                changed = true;
            }
        }
        if !changed && it > 0 {
            break;
        }
        let mut sums = vec![vec![0.0; d]; k];
        let mut cnts = vec![0usize; k];
        for i in 0..n {
            let c = labels[i] as usize;
            cnts[c] += 1;
            for j in 0..d {
                sums[c][j] += x[i][j];
            }
        }
        for c in 0..k {
            if cnts[c] > 0 {
                for j in 0..d {
                    centers[c][j] = sums[c][j] / cnts[c] as f64;
                }
            } else {
                // Re-seed empty cluster to farthest point from its assigned centroid
                let mut far_idx = 0;
                let mut far_d = -1.0;
                for i in 0..n {
                    let dd = euclidean_sq(&x[i], &centers[labels[i] as usize]);
                    if dd > far_d {
                        far_d = dd;
                        far_idx = i;
                    }
                }
                centers[c] = x[far_idx].clone();
            }
        }
    }

    let sse: f64 = (0..n)
        .map(|i| euclidean_sq(&x[i], &centers[labels[i] as usize]))
        .sum();
    (labels, centers, sse, iters)
}

// 1. K-Means (Lloyd's algoritme) - Random partition initialization + Lloyd updates
fn kmeans_lloyd(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let max_iter = cfg.max_iters(60);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed);

    // Random Partition initialization (classic Lloyd)
    let mut sums = vec![vec![0.0; d]; k];
    let mut cnts = vec![0usize; k];
    for i in 0..n {
        let c = rng.gen_range(0..k);
        cnts[c] += 1;
        for j in 0..d {
            sums[c][j] += data.x[i][j];
        }
    }
    let centers: Vec<Vec<f64>> = (0..k)
        .map(|c| {
            if cnts[c] > 0 {
                sums[c].iter().map(|&v| v / cnts[c] as f64).collect()
            } else {
                data.x[rng.gen_range(0..n)].clone()
            }
        })
        .collect();

    let (labels, _, sse, iters) = run_lloyd_from_centers(&data.x, centers, max_iter);
    AlgorithmOutput::new(
        labels,
        format!("k={}, init=Tilfeldig partisjon", k),
        format!("Konvergerte på {} iterasjoner (SSE={:.2})", iters, sse),
    )
}

// 2. K-Means++
fn kmeans_plus_plus(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let k = cfg.effective_k(data.n_samples);
    let max_iter = cfg.max_iters(60);
    let restarts = 5;
    let mut best_labels = vec![0i32; data.n_samples];
    let mut best_sse = f64::INFINITY;
    let mut best_iters = 0;

    for r in 0..restarts {
        let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed.wrapping_add(r as u64 * 97 + 1));
        let centers = init_kmeans_pp(&data.x, k, &mut rng);
        let (labels, _, sse, iters) = run_lloyd_from_centers(&data.x, centers, max_iter);
        if sse < best_sse {
            best_sse = sse;
            best_labels = labels;
            best_iters = iters;
        }
    }
    AlgorithmOutput::new(
        best_labels,
        format!("k={}, n_init={}", k, restarts),
        format!("Beste av {} kjøringer: {} iter., SSE={:.2}", restarts, best_iters, best_sse),
    )
}

// 3. Mini-Batch K-Means
fn mini_batch_kmeans(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let batch_size = (n / 3).clamp(8, 64).min(n);
    let max_iter = cfg.max_iters(80);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 3);

    let mut centers = init_kmeans_pp(&data.x, k, &mut rng);
    let mut counts = vec![0usize; k];
    let mut indices: Vec<usize> = (0..n).collect();

    for _ in 0..max_iter {
        indices.shuffle(&mut rng);
        let batch = &indices[..batch_size];
        let assigned: Vec<usize> = batch
            .iter()
            .map(|&idx| {
                (0..k)
                    .min_by(|&a, &b| {
                        euclidean_sq(&data.x[idx], &centers[a])
                            .total_cmp(&euclidean_sq(&data.x[idx], &centers[b]))
                    })
                    .unwrap_or(0)
            })
            .collect();

        for (&idx, &c) in batch.iter().zip(assigned.iter()) {
            counts[c] += 1;
            let eta = 1.0 / counts[c] as f64;
            for j in 0..d {
                centers[c][j] = (1.0 - eta) * centers[c][j] + eta * data.x[idx][j];
            }
        }
    }

    let labels: Vec<i32> = (0..n)
        .map(|i| {
            (0..k)
                .min_by(|&a, &b| {
                    euclidean_sq(&data.x[i], &centers[a])
                        .total_cmp(&euclidean_sq(&data.x[i], &centers[b]))
                })
                .unwrap_or(0) as i32
        })
        .collect();

    AlgorithmOutput::new(
        labels,
        format!("k={}, batch={}", k, batch_size),
        format!("{} mini-batch oppdateringer fullført", max_iter),
    )
}

fn pam_on_indices(dist: &[Vec<f64>], subset: &[usize], k: usize, max_iter: usize) -> (Vec<usize>, f64) {
    let m = subset.len();
    let kk = k.min(m).max(1);
    // BUILD phase: greedy medoid selection
    let mut medoids = Vec::with_capacity(kk);
    let first = (0..m)
        .min_by(|&a, &b| {
            let sa: f64 = (0..m).map(|j| dist[subset[a]][subset[j]]).sum();
            let sb: f64 = (0..m).map(|j| dist[subset[b]][subset[j]]).sum();
            sa.total_cmp(&sb)
        })
        .unwrap_or(0);
    medoids.push(first);

    while medoids.len() < kk {
        let mut best_cand = 0;
        let mut best_gain = f64::NEG_INFINITY;
        for c in 0..m {
            if medoids.contains(&c) {
                continue;
            }
            let mut gain = 0.0;
            for j in 0..m {
                let d_cur = medoids
                    .iter()
                    .map(|&med| dist[subset[j]][subset[med]])
                    .fold(f64::INFINITY, f64::min);
                let d_new = dist[subset[j]][subset[c]];
                if d_cur > d_new {
                    gain += d_cur - d_new;
                }
            }
            if gain > best_gain {
                best_gain = gain;
                best_cand = c;
            }
        }
        medoids.push(best_cand);
    }

    // SWAP phase
    let eval_cost = |meds: &[usize]| -> f64 {
        (0..m)
            .map(|j| {
                meds.iter()
                    .map(|&med| dist[subset[j]][subset[med]])
                    .fold(f64::INFINITY, f64::min)
            })
            .sum()
    };

    let mut cur_cost = eval_cost(&medoids);
    for _ in 0..max_iter {
        let mut best_swap: Option<(usize, usize, f64)> = None;
        for mi in 0..kk {
            for cand in 0..m {
                if medoids.contains(&cand) {
                    continue;
                }
                let old = medoids[mi];
                medoids[mi] = cand;
                let cost = eval_cost(&medoids);
                medoids[mi] = old;
                if cost + 1e-9 < best_swap.map(|s| s.2).unwrap_or(cur_cost) {
                    best_swap = Some((mi, cand, cost));
                }
            }
        }
        if let Some((mi, cand, new_cost)) = best_swap {
            medoids[mi] = cand;
            cur_cost = new_cost;
        } else {
            break;
        }
    }

    let global_meds: Vec<usize> = medoids.into_iter().map(|i| subset[i]).collect();
    (global_meds, cur_cost)
}

fn assign_to_medoids(dist: &[Vec<f64>], medoids: &[usize]) -> (Vec<i32>, f64) {
    let n = dist.len();
    let mut labels = vec![0i32; n];
    let mut cost = 0.0;
    for i in 0..n {
        let mut best_c = 0;
        let mut best_d = f64::INFINITY;
        for (c, &m) in medoids.iter().enumerate() {
            if dist[i][m] < best_d {
                best_d = dist[i][m];
                best_c = c;
            }
        }
        labels[i] = best_c as i32;
        cost += best_d;
    }
    (labels, cost)
}

// 4. K-Medoids (PAM)
fn kmedoids_pam(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let k = cfg.effective_k(n);
    let max_iter = cfg.max_iters(15);
    let all: Vec<usize> = (0..n).collect();
    let (medoids, _) = pam_on_indices(&data.dist_matrix, &all, k, max_iter);
    let (labels, total_cost) = assign_to_medoids(&data.dist_matrix, &medoids);
    AlgorithmOutput::new(
        labels,
        format!("k={}, BUILD+SWAP", k),
        format!("Total medoidkostnad={:.2}, medoider={:?}", total_cost, &medoids[..medoids.len().min(4)]),
    )
}

// 5. CLARA
fn clara(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let k = cfg.effective_k(n);
    let sample_size = (40 + 2 * k).min(n);
    let n_samples_runs = 5;
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 5);
    let mut all_idx: Vec<usize> = (0..n).collect();

    let mut best_labels = vec![0i32; n];
    let mut best_cost = f64::INFINITY;

    for _ in 0..n_samples_runs {
        all_idx.shuffle(&mut rng);
        let sub = &all_idx[..sample_size];
        let (meds, _) = pam_on_indices(&data.dist_matrix, sub, k, 10);
        let (labels, cost) = assign_to_medoids(&data.dist_matrix, &meds);
        if cost < best_cost {
            best_cost = cost;
            best_labels = labels;
        }
    }
    AlgorithmOutput::new(
        best_labels,
        format!("k={}, utvalg={}, runder={}", k, sample_size, n_samples_runs),
        format!("Beste totalkostnad over {} delutvalg: {:.2}", n_samples_runs, best_cost),
    )
}

// 6. CLARANS
fn clarans(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let k = cfg.effective_k(n);
    let num_local = 3;
    let max_neighbor = ((k * (n.saturating_sub(k))) as f64 * 0.05).round().clamp(15.0, 120.0) as usize;
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 6);
    let mut all_idx: Vec<usize> = (0..n).collect();

    let mut global_best_meds = Vec::new();
    let mut global_best_cost = f64::INFINITY;

    for _ in 0..num_local {
        all_idx.shuffle(&mut rng);
        let mut current: Vec<usize> = all_idx[..k].to_vec();
        let (_, mut cur_cost) = assign_to_medoids(&data.dist_matrix, &current);
        let mut tries = 0;
        while tries < max_neighbor {
            let slot = rng.gen_range(0..k);
            let cand = rng.gen_range(0..n);
            if current.contains(&cand) {
                tries += 1;
                continue;
            }
            let old = current[slot];
            current[slot] = cand;
            let (_, new_cost) = assign_to_medoids(&data.dist_matrix, &current);
            if new_cost + 1e-9 < cur_cost {
                cur_cost = new_cost;
                tries = 0;
            } else {
                current[slot] = old;
                tries += 1;
            }
        }
        if cur_cost < global_best_cost {
            global_best_cost = cur_cost;
            global_best_meds = current;
        }
    }
    let (labels, _) = assign_to_medoids(&data.dist_matrix, &global_best_meds);
    AlgorithmOutput::new(
        labels,
        format!("k={}, numlocal={}, maxnbr={}", k, num_local, max_neighbor),
        format!("Randomisert grafsøk fant kostnad={:.2}", global_best_cost),
    )
}

// 7. K-Modes
fn kmodes(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let c_dim = data.categorical[0].len();
    let k = cfg.effective_k(n);
    let max_iter = cfg.max_iters(50);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 7);

    let mut idxs: Vec<usize> = (0..n).collect();
    idxs.shuffle(&mut rng);
    let mut modes: Vec<Vec<usize>> = idxs[..k].iter().map(|&i| data.categorical[i].clone()).collect();
    let mut labels = vec![0i32; n];
    let mut iters = 0;

    for it in 0..max_iter {
        iters = it + 1;
        let mut changed = false;
        for i in 0..n {
            let mut best_c = 0;
            let mut best_d = usize::MAX;
            for (c, mode) in modes.iter().enumerate() {
                let hd = hamming_dist(&data.categorical[i], mode);
                if hd < best_d {
                    best_d = hd;
                    best_c = c;
                }
            }
            if labels[i] != best_c as i32 {
                labels[i] = best_c as i32;
                changed = true;
            }
        }
        if !changed && it > 0 {
            break;
        }
        for c in 0..k {
            let members: Vec<usize> = (0..n).filter(|&i| labels[i] == c as i32).collect();
            if members.is_empty() {
                modes[c] = data.categorical[rng.gen_range(0..n)].clone();
                continue;
            }
            for j in 0..c_dim {
                let mut freq: HashMap<usize, usize> = HashMap::new();
                for &m in &members {
                    *freq.entry(data.categorical[m][j]).or_insert(0) += 1;
                }
                if let Some((&mode_val, _)) = freq.iter().max_by_key(|(_, &cnt)| cnt) {
                    modes[c][j] = mode_val;
                }
            }
        }
    }

    let source_note = if data.has_real_categorical {
        "brukte kategoriske kolonner"
    } else {
        "brukte kvantil-diskretiserte kategorier"
    };
    AlgorithmOutput::new(
        labels,
        format!("k={}, Hamming-avstand", k),
        format!("Konvergerte på {} iter. ({})", iters, source_note),
    )
}

// 8. K-Prototypes
fn kprototypes(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let c_dim = data.categorical[0].len();
    let k = cfg.effective_k(n);
    let gamma = 0.5;
    let max_iter = cfg.max_iters(50);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 8);

    let mut idxs: Vec<usize> = (0..n).collect();
    idxs.shuffle(&mut rng);
    let mut num_centers: Vec<Vec<f64>> = idxs[..k].iter().map(|&i| data.x[i].clone()).collect();
    let mut cat_modes: Vec<Vec<usize>> = idxs[..k].iter().map(|&i| data.categorical[i].clone()).collect();
    let mut labels = vec![0i32; n];

    for it in 0..max_iter {
        let mut changed = false;
        for i in 0..n {
            let mut best_c = 0;
            let mut best_cost = f64::INFINITY;
            for c in 0..k {
                let cost = euclidean_sq(&data.x[i], &num_centers[c])
                    + gamma * hamming_dist(&data.categorical[i], &cat_modes[c]) as f64;
                if cost < best_cost {
                    best_cost = cost;
                    best_c = c;
                }
            }
            if labels[i] != best_c as i32 {
                labels[i] = best_c as i32;
                changed = true;
            }
        }
        if !changed && it > 0 {
            break;
        }
        for c in 0..k {
            let members: Vec<usize> = (0..n).filter(|&i| labels[i] == c as i32).collect();
            if members.is_empty() {
                continue;
            }
            for j in 0..d {
                num_centers[c][j] = members.iter().map(|&m| data.x[m][j]).sum::<f64>() / members.len() as f64;
            }
            for j in 0..c_dim {
                let mut freq: HashMap<usize, usize> = HashMap::new();
                for &m in &members {
                    *freq.entry(data.categorical[m][j]).or_insert(0) += 1;
                }
                if let Some((&mv, _)) = freq.iter().max_by_key(|(_, &cnt)| cnt) {
                    cat_modes[c][j] = mv;
                }
            }
        }
    }
    AlgorithmOutput::new(
        labels,
        format!("k={}, gamma={:.2}", k, gamma),
        "Kombinerte euklidsk avstand og Hamming-modus".to_string(),
    )
}

// 9. Bisecting K-Means
fn bisecting_kmeans(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let k = cfg.effective_k(n);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 9);

    let mut clusters: Vec<Vec<usize>> = vec![(0..n).collect()];
    while clusters.len() < k {
        // Choose cluster with largest SSE
        let split_idx = clusters
            .iter()
            .enumerate()
            .filter(|(_, members)| members.len() >= 2)
            .max_by(|(_, a), (_, b)| cluster_sse(&data.x, a).total_cmp(&cluster_sse(&data.x, b)))
            .map(|(idx, _)| idx);

        let Some(idx) = split_idx else { break };
        let members = clusters.remove(idx);
        let sub_x: Vec<Vec<f64>> = members.iter().map(|&i| data.x[i].clone()).collect();

        let mut best_split = (Vec::new(), Vec::new());
        let mut best_sse = f64::INFINITY;
        for _ in 0..4 {
            let c0 = init_kmeans_pp(&sub_x, 2, &mut rng);
            let (lbl, _, sse, _) = run_lloyd_from_centers(&sub_x, c0, 30);
            if sse < best_sse {
                best_sse = sse;
                let left: Vec<usize> = members
                    .iter()
                    .enumerate()
                    .filter(|&(j, _)| lbl[j] == 0)
                    .map(|(_, &id)| id)
                    .collect();
                let right: Vec<usize> = members
                    .iter()
                    .enumerate()
                    .filter(|&(j, _)| lbl[j] != 0)
                    .map(|(_, &id)| id)
                    .collect();
                best_split = (left, right);
            }
        }
        if best_split.0.is_empty() || best_split.1.is_empty() {
            let mid = members.len() / 2;
            clusters.push(members[..mid].to_vec());
            clusters.push(members[mid..].to_vec());
        } else {
            clusters.push(best_split.0);
            clusters.push(best_split.1);
        }
    }

    let mut labels = vec![0i32; n];
    for (c, members) in clusters.iter().enumerate() {
        for &i in members {
            labels[i] = c as i32;
        }
    }
    AlgorithmOutput::new(
        labels,
        format!("k={}, 4 bisections/steg", k),
        format!("Top-down rekursiv splitting til {} klynger", clusters.len()),
    )
}

fn cluster_sse(x: &[Vec<f64>], members: &[usize]) -> f64 {
    if members.is_empty() {
        return 0.0;
    }
    let d = x[0].len();
    let mut mean = vec![0.0; d];
    for &i in members {
        for j in 0..d {
            mean[j] += x[i][j];
        }
    }
    for v in &mut mean {
        *v /= members.len() as f64;
    }
    members.iter().map(|&i| euclidean_sq(&x[i], &mean)).sum()
}

fn bic_for_partition(x: &[Vec<f64>], clusters: &[Vec<usize>]) -> f64 {
    let n: usize = clusters.iter().map(|c| c.len()).sum();
    let k = clusters.len();
    let d = x[0].len();
    if n <= k || k == 0 {
        return f64::NEG_INFINITY;
    }
    let total_sse: f64 = clusters.iter().map(|c| cluster_sse(x, c)).sum();
    let var = (total_sse / ((n - k) as f64 * d as f64)).max(1e-8);
    let mut log_lik = 0.0;
    for c in clusters {
        let rn = c.len() as f64;
        if rn > 0.0 {
            log_lik += rn * (rn / n as f64).ln()
                - 0.5 * rn * d as f64 * (2.0 * std::f64::consts::PI * var).ln()
                - 0.5 * cluster_sse(x, c) / var;
        }
    }
    let p_params = (k * (d + 1)) as f64;
    log_lik - 0.5 * p_params * (n as f64).ln()
}

// 10. X-Means (Pelleg & Moore BIC automatic k)
fn xmeans_bic(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let k_max = (cfg.effective_k(n) + 5).min(12).min(n / 3).max(2);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 10);

    let mut clusters: Vec<Vec<usize>> = vec![(0..n).collect()];
    loop {
        if clusters.len() >= k_max {
            break;
        }
        let mut next_clusters = Vec::new();
        let mut any_split = false;

        for c in clusters {
            if c.len() < 6 || next_clusters.len() + 1 >= k_max {
                next_clusters.push(c);
                continue;
            }
            let parent_bic = bic_for_partition(&data.x, std::slice::from_ref(&c));
            let sub_x: Vec<Vec<f64>> = c.iter().map(|&i| data.x[i].clone()).collect();
            let c0 = init_kmeans_pp(&sub_x, 2, &mut rng);
            let (lbl, _, _, _) = run_lloyd_from_centers(&sub_x, c0, 25);
            let left: Vec<usize> = c
                .iter()
                .enumerate()
                .filter(|&(j, _)| lbl[j] == 0)
                .map(|(_, &id)| id)
                .collect();
            let right: Vec<usize> = c
                .iter()
                .enumerate()
                .filter(|&(j, _)| lbl[j] != 0)
                .map(|(_, &id)| id)
                .collect();

            if left.len() >= 2 && right.len() >= 2 {
                let child_bic = bic_for_partition(&data.x, &[left.clone(), right.clone()]);
                if child_bic > parent_bic {
                    next_clusters.push(left);
                    next_clusters.push(right);
                    any_split = true;
                } else {
                    next_clusters.push(c);
                }
            } else {
                next_clusters.push(c);
            }
        }
        clusters = next_clusters;
        if !any_split {
            break;
        }
    }

    let final_bic = bic_for_partition(&data.x, &clusters);
    let mut labels = vec![0i32; n];
    for (cid, members) in clusters.iter().enumerate() {
        for &i in members {
            labels[i] = cid as i32;
        }
    }
    AlgorithmOutput::new(
        labels,
        format!("Auto-k (maks={}), BIC", k_max),
        format!("BIC valgte k={} (BIC={:.1})", clusters.len(), final_bic),
    )
}

// 11. G-Means (Hamerly & Elkan 2003 — Anderson-Darling normality test along split axis)
fn gmeans_anderson_darling(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k_max = (cfg.effective_k(n) + 5).min(12).min(n / 4).max(2);
    let crit_ad = 1.8692; // alpha = 0.0001 critical value (standard G-Means default)
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 11);

    let mut clusters: Vec<Vec<usize>> = vec![(0..n).collect()];
    let mut splits_performed = 0;

    for _round in 0..8 {
        if clusters.len() >= k_max {
            break;
        }
        let mut next = Vec::new();
        let mut changed = false;

        for c in clusters {
            if c.len() < 8 || next.len() + 1 >= k_max {
                next.push(c);
                continue;
            }
            let sub_x: Vec<Vec<f64>> = c.iter().map(|&i| data.x[i].clone()).collect();
            let c0 = init_kmeans_pp(&sub_x, 2, &mut rng);
            let (lbl, centers, _, _) = run_lloyd_from_centers(&sub_x, c0, 25);

            // Project onto vector v = c1 - c2
            let mut v: Vec<f64> = (0..d).map(|j| centers[0][j] - centers[1][j]).collect();
            let vn = dot(&v, &v).sqrt();
            if vn < 1e-9 {
                next.push(c);
                continue;
            }
            for x in &mut v {
                *x /= vn;
            }
            let proj: Vec<f64> = sub_x.iter().map(|row| dot(row, &v)).collect();
            let ad = anderson_darling_normal(&proj);

            let left: Vec<usize> = c
                .iter()
                .enumerate()
                .filter(|&(j, _)| lbl[j] == 0)
                .map(|(_, &id)| id)
                .collect();
            let right: Vec<usize> = c
                .iter()
                .enumerate()
                .filter(|&(j, _)| lbl[j] != 0)
                .map(|(_, &id)| id)
                .collect();

            if ad > crit_ad && left.len() >= 3 && right.len() >= 3 {
                next.push(left);
                next.push(right);
                splits_performed += 1;
                changed = true;
            } else {
                next.push(c);
            }
        }
        clusters = next;
        if !changed {
            break;
        }
    }

    let mut labels = vec![0i32; n];
    for (cid, members) in clusters.iter().enumerate() {
        for &i in members {
            labels[i] = cid as i32;
        }
    }
    AlgorithmOutput::new(
        labels,
        format!("AD-terskel={:.2}, maks_k={}", crit_ad, k_max),
        format!("Utledet k={} etter {} normalitetssplitter", clusters.len(), splits_performed),
    )
}

// 12. Spherical K-Means
fn spherical_kmeans(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let max_iter = cfg.max_iters(60);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 12);

    let mut centers = init_kmeans_pp(&data.x_unit, k, &mut rng)
        .into_iter()
        .map(|c| normalize_vec(&c))
        .collect::<Vec<_>>();
    let mut labels = vec![0i32; n];
    let mut iters = 0;

    for it in 0..max_iter {
        iters = it + 1;
        let mut changed = false;
        for i in 0..n {
            let mut best_c = 0;
            let mut best_sim = f64::NEG_INFINITY;
            for (c, cen) in centers.iter().enumerate() {
                let sim = dot(&data.x_unit[i], cen);
                if sim > best_sim {
                    best_sim = sim;
                    best_c = c;
                }
            }
            if labels[i] != best_c as i32 {
                labels[i] = best_c as i32;
                changed = true;
            }
        }
        if !changed && it > 0 {
            break;
        }
        let mut sums = vec![vec![0.0; d]; k];
        let mut cnts = vec![0usize; k];
        for i in 0..n {
            let c = labels[i] as usize;
            cnts[c] += 1;
            for j in 0..d {
                sums[c][j] += data.x_unit[i][j];
            }
        }
        for c in 0..k {
            if cnts[c] > 0 {
                centers[c] = normalize_vec(&sums[c]);
            }
        }
    }

    let avg_cos: f64 = (0..n)
        .map(|i| dot(&data.x_unit[i], &centers[labels[i] as usize]))
        .sum::<f64>()
        / n as f64;
    AlgorithmOutput::new(
        labels,
        format!("k={}, cosinus-likhet", k),
        format!("Konvergerte på {} iter., snitt cosinus={:.3}", iters, avg_cos),
    )
}

// 13. Elkan's K-Means (exact triangle inequality acceleration with per-point-per-center lower bounds)
fn elkan_kmeans(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let max_iter = cfg.max_iters(60);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 13);

    let mut centers = init_kmeans_pp(&data.x, k, &mut rng);
    let mut lower = vec![vec![0.0; k]; n];
    let mut upper = vec![0.0; n];
    let mut labels = vec![0usize; n];
    let mut dist_evals = 0usize;
    let mut brute_evals = 0usize;

    // Initial assignment
    for i in 0..n {
        let mut best_c = 0;
        let mut best_d = f64::INFINITY;
        for c in 0..k {
            let dist = euclidean(&data.x[i], &centers[c]);
            dist_evals += 1;
            lower[i][c] = dist;
            if dist < best_d {
                best_d = dist;
                best_c = c;
            }
        }
        labels[i] = best_c;
        upper[i] = best_d;
    }

    let mut iters = 0;
    for it in 0..max_iter {
        iters = it + 1;
        brute_evals += n * k;

        let mut cdist = vec![vec![0.0; k]; k];
        let mut s_half = vec![f64::INFINITY; k];
        for a in 0..k {
            for b in (a + 1)..k {
                let dd = euclidean(&centers[a], &centers[b]);
                cdist[a][b] = dd;
                cdist[b][a] = dd;
                let half = 0.5 * dd;
                if half < s_half[a] {
                    s_half[a] = half;
                }
                if half < s_half[b] {
                    s_half[b] = half;
                }
            }
        }

        let mut changed = false;
        for i in 0..n {
            let mut c_cur = labels[i];
            if upper[i] <= s_half[c_cur] {
                continue;
            }
            let mut r_outdated = true;
            for c in 0..k {
                if c == c_cur {
                    continue;
                }
                if upper[i] <= lower[i][c] || upper[i] <= 0.5 * cdist[c_cur][c] {
                    continue;
                }
                if r_outdated {
                    upper[i] = euclidean(&data.x[i], &centers[c_cur]);
                    lower[i][c_cur] = upper[i];
                    dist_evals += 1;
                    r_outdated = false;
                    if upper[i] <= lower[i][c] || upper[i] <= 0.5 * cdist[c_cur][c] {
                        continue;
                    }
                }
                let d_ic = euclidean(&data.x[i], &centers[c]);
                lower[i][c] = d_ic;
                dist_evals += 1;
                if d_ic < upper[i] {
                    c_cur = c;
                    labels[i] = c;
                    upper[i] = d_ic;
                    changed = true;
                }
            }
        }

        // Update centers and drift
        let old_centers = centers.clone();
        let mut sums = vec![vec![0.0; d]; k];
        let mut cnts = vec![0usize; k];
        for i in 0..n {
            let c = labels[i];
            cnts[c] += 1;
            for j in 0..d {
                sums[c][j] += data.x[i][j];
            }
        }
        let mut drift = vec![0.0; k];
        for c in 0..k {
            if cnts[c] > 0 {
                for j in 0..d {
                    centers[c][j] = sums[c][j] / cnts[c] as f64;
                }
                drift[c] = euclidean(&centers[c], &old_centers[c]);
            }
        }
        for i in 0..n {
            for c in 0..k {
                lower[i][c] = (lower[i][c] - drift[c]).max(0.0);
            }
            upper[i] += drift[labels[i]];
        }
        if !changed && it > 0 {
            break;
        }
    }

    let saved_pct = if brute_evals > 0 {
        (1.0 - (dist_evals as f64 / (brute_evals + n * k) as f64)).clamp(0.0, 1.0) * 100.0
    } else {
        0.0
    };
    let out_labels: Vec<i32> = labels.into_iter().map(|c| c as i32).collect();
    AlgorithmOutput::new(
        out_labels,
        format!("k={}, Elkan-grenser", k),
        format!("{} iter., sparte {:.1}% avstandsberegninger", iters, saved_pct),
    )
}

// 14. Hamerly's K-Means (single upper & single lower bound per point)
fn hamerly_kmeans(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let max_iter = cfg.max_iters(60);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 14);

    let mut centers = init_kmeans_pp(&data.x, k, &mut rng);
    let mut upper = vec![f64::INFINITY; n];
    let mut lower = vec![0.0; n];
    let mut labels = vec![0usize; n];
    let mut skipped = 0usize;
    let mut total_checks = 0usize;

    for i in 0..n {
        let mut d1 = f64::INFINITY;
        let mut d2 = f64::INFINITY;
        let mut best_c = 0;
        for c in 0..k {
            let dd = euclidean(&data.x[i], &centers[c]);
            if dd < d1 {
                d2 = d1;
                d1 = dd;
                best_c = c;
            } else if dd < d2 {
                d2 = dd;
            }
        }
        labels[i] = best_c;
        upper[i] = d1;
        lower[i] = d2;
    }

    for it in 0..max_iter {
        let mut s_half = vec![f64::INFINITY; k];
        for a in 0..k {
            for b in (a + 1)..k {
                let half = 0.5 * euclidean(&centers[a], &centers[b]);
                if half < s_half[a] {
                    s_half[a] = half;
                }
                if half < s_half[b] {
                    s_half[b] = half;
                }
            }
        }

        let mut changed = false;
        for i in 0..n {
            total_checks += 1;
            let m = s_half[labels[i]].max(lower[i]);
            if upper[i] <= m {
                skipped += 1;
                continue;
            }
            upper[i] = euclidean(&data.x[i], &centers[labels[i]]);
            if upper[i] <= m {
                skipped += 1;
                continue;
            }
            let old_c = labels[i];
            let mut d1 = f64::INFINITY;
            let mut d2 = f64::INFINITY;
            let mut best_c = 0;
            for c in 0..k {
                let dd = euclidean(&data.x[i], &centers[c]);
                if dd < d1 {
                    d2 = d1;
                    d1 = dd;
                    best_c = c;
                } else if dd < d2 {
                    d2 = dd;
                }
            }
            labels[i] = best_c;
            upper[i] = d1;
            lower[i] = d2;
            if best_c != old_c {
                changed = true;
            }
        }

        let old_centers = centers.clone();
        let mut sums = vec![vec![0.0; d]; k];
        let mut cnts = vec![0usize; k];
        for i in 0..n {
            cnts[labels[i]] += 1;
            for j in 0..d {
                sums[labels[i]][j] += data.x[i][j];
            }
        }
        let mut drift = vec![0.0; k];
        for c in 0..k {
            if cnts[c] > 0 {
                for j in 0..d {
                    centers[c][j] = sums[c][j] / cnts[c] as f64;
                }
                drift[c] = euclidean(&centers[c], &old_centers[c]);
            }
        }
        // Largest and second largest center drift
        let mut max_drift = 0.0;
        let mut max_c = 0;
        let mut second_drift = 0.0;
        for c in 0..k {
            if drift[c] > max_drift {
                second_drift = max_drift;
                max_drift = drift[c];
                max_c = c;
            } else if drift[c] > second_drift {
                second_drift = drift[c];
            }
        }
        for i in 0..n {
            upper[i] += drift[labels[i]];
            let sub = if labels[i] == max_c { second_drift } else { max_drift };
            lower[i] = (lower[i] - sub).max(0.0);
        }
        if !changed && it > 0 {
            break;
        }
    }

    let skip_pct = if total_checks > 0 {
        (skipped as f64 / total_checks as f64) * 100.0
    } else {
        0.0
    };
    AlgorithmOutput::new(
        labels.into_iter().map(|c| c as i32).collect(),
        format!("k={}, 2 grenser/punkt", k),
        format!("Hoppet over {:.1}% fulle punkt-søk via Hamerly-grenser", skip_pct),
    )
}

// 15. K-Medians (L1 Manhattan distance + component-wise median)
fn kmedians(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let max_iter = cfg.max_iters(60);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 15);

    let mut centers = init_kmeans_pp(&data.x, k, &mut rng);
    let mut labels = vec![0i32; n];
    let mut iters = 0;

    for it in 0..max_iter {
        iters = it + 1;
        let mut changed = false;
        for i in 0..n {
            let mut best_c = 0;
            let mut best_d = f64::INFINITY;
            for (c, cen) in centers.iter().enumerate() {
                let dd = manhattan(&data.x[i], cen);
                if dd < best_d {
                    best_d = dd;
                    best_c = c;
                }
            }
            if labels[i] != best_c as i32 {
                labels[i] = best_c as i32;
                changed = true;
            }
        }
        if !changed && it > 0 {
            break;
        }
        for c in 0..k {
            let members: Vec<usize> = (0..n).filter(|&i| labels[i] == c as i32).collect();
            if members.is_empty() {
                continue;
            }
            for j in 0..d {
                let mut vals: Vec<f64> = members.iter().map(|&m| data.x[m][j]).collect();
                centers[c][j] = median_of_slice(&mut vals);
            }
        }
    }

    let l1_loss: f64 = (0..n)
        .map(|i| manhattan(&data.x[i], &centers[labels[i] as usize]))
        .sum();
    AlgorithmOutput::new(
        labels,
        format!("k={}, L1 (Manhattan)", k),
        format!("Konvergerte på {} iter., L1-tap={:.2}", iters, l1_loss),
    )
}

// 16. Kernel K-Means (RBF Gaussian Kernel in RKHS)
fn kernel_kmeans(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let k = cfg.effective_k(n);
    let max_iter = cfg.max_iters(40);
    let sigma = median_pairwise_dist(&data.dist_matrix).max(1e-4);
    let gamma = 1.0 / (2.0 * sigma * sigma);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 16);

    // Precompute kernel matrix K_ij = exp(-gamma * ||xi - xj||^2)
    let mut kmat = vec![vec![0.0; n]; n];
    for i in 0..n {
        kmat[i][i] = 1.0;
        for j in (i + 1)..n {
            let val = (-gamma * data.dist_matrix[i][j].powi(2)).exp();
            kmat[i][j] = val;
            kmat[j][i] = val;
        }
    }

    // Initialize via K-Means++ labels
    let c0 = init_kmeans_pp(&data.x, k, &mut rng);
    let (mut labels, _, _, _) = run_lloyd_from_centers(&data.x, c0, 10);

    let mut iters = 0;
    for it in 0..max_iter {
        iters = it + 1;
        let mut members: Vec<Vec<usize>> = vec![Vec::new(); k];
        for i in 0..n {
            members[labels[i] as usize].push(i);
        }
        // Third term per cluster: (1 / |C|^2) * sum_{p,q in C} K(p,q)
        let mut term3 = vec![0.0; k];
        for c in 0..k {
            let sz = members[c].len();
            if sz == 0 {
                continue;
            }
            let mut s = 0.0;
            for &p in &members[c] {
                for &q in &members[c] {
                    s += kmat[p][q];
                }
            }
            term3[c] = s / (sz * sz) as f64;
        }

        let mut changed = false;
        for i in 0..n {
            let mut best_c = 0;
            let mut best_d = f64::INFINITY;
            for c in 0..k {
                let sz = members[c].len();
                if sz == 0 {
                    continue;
                }
                let term2: f64 = members[c].iter().map(|&p| kmat[i][p]).sum::<f64>() * (2.0 / sz as f64);
                let dist_rkhs = kmat[i][i] - term2 + term3[c];
                if dist_rkhs < best_d {
                    best_d = dist_rkhs;
                    best_c = c;
                }
            }
            if labels[i] != best_c as i32 {
                labels[i] = best_c as i32;
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }

    AlgorithmOutput::new(
        labels,
        format!("k={}, RBF sigma={:.2}", k, sigma),
        format!("RBF-kjerne i RKHS konvergerte på {} iterasjoner", iters),
    )
}

// 17. Global K-Means (Likas et al. 2003 — deterministic incremental center search)
fn global_kmeans(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);

    // Start with 1-cluster centroid = global mean
    let mut mean = vec![0.0; d];
    for row in &data.x {
        for j in 0..d {
            mean[j] += row[j];
        }
    }
    for v in &mut mean {
        *v /= n as f64;
    }
    let mut centers = vec![mean];
    let mut best_labels = vec![0i32; n];
    let mut best_sse = f64::INFINITY;

    // Fast Global K-Means upper-bound candidate scoring (deterministic!)
    for _step in 1..k {
        let cur_dists: Vec<f64> = (0..n)
            .map(|i| {
                centers
                    .iter()
                    .map(|c| euclidean_sq(&data.x[i], c))
                    .fold(f64::INFINITY, f64::min)
            })
            .collect();

        // Pick top candidates maximizing guaranteed reduction b_j = sum_i max(d_i - ||xi - xj||^2, 0)
        let mut scores: Vec<(usize, f64)> = (0..n)
            .map(|cand| {
                let b: f64 = (0..n)
                    .map(|i| (cur_dists[i] - data.dist_matrix[i][cand].powi(2)).max(0.0))
                    .sum();
                (cand, b)
            })
            .collect();
        scores.sort_by(|a, b| b.1.total_cmp(&a.1));

        let mut step_best_centers = centers.clone();
        let mut step_best_sse = f64::INFINITY;
        let mut step_best_labels = vec![0i32; n];

        for &(cand, _) in scores.iter().take(6.min(n)) {
            let mut trial = centers.clone();
            trial.push(data.x[cand].clone());
            let (lbl, cen, sse, _) = run_lloyd_from_centers(&data.x, trial, 25);
            if sse < step_best_sse {
                step_best_sse = sse;
                step_best_centers = cen;
                step_best_labels = lbl;
            }
        }
        centers = step_best_centers;
        best_labels = step_best_labels;
        best_sse = step_best_sse;
    }

    AlgorithmOutput::new(
        best_labels,
        format!("k={}, deterministisk", k),
        format!("Inkrementelt deterministisk søk, SSE={:.2}", best_sse),
    )
}

// 18. ISODATA (Iterative Self-Organizing Data Analysis Technique — dynamic split & merge)
fn isodata(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k_target = cfg.effective_k(n);
    let min_pts = (n / (k_target * 4)).max(3);
    let max_iter = cfg.max_iters(35);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 18);

    // Compute global std to calibrate split/merge thresholds
    let mut global_var = 0.0;
    for j in 0..d {
        let mu = data.x.iter().map(|r| r[j]).sum::<f64>() / n as f64;
        global_var += data.x.iter().map(|r| (r[j] - mu).powi(2)).sum::<f64>() / n as f64;
    }
    let std_split = (global_var / d as f64).sqrt() * 0.85;
    let merge_dist = median_pairwise_dist(&data.dist_matrix) * 0.35;

    let mut centers = init_kmeans_pp(&data.x, k_target, &mut rng);
    let mut labels = vec![0i32; n];
    let mut splits = 0usize;
    let mut merges = 0usize;

    for it in 0..max_iter {
        let k = centers.len();
        for i in 0..n {
            labels[i] = (0..k)
                .min_by(|&a, &b| {
                    euclidean_sq(&data.x[i], &centers[a])
                        .total_cmp(&euclidean_sq(&data.x[i], &centers[b]))
                })
                .unwrap_or(0) as i32;
        }
        // Recompute centers & drop tiny clusters
        let mut next_centers = Vec::new();
        let mut cluster_members = Vec::new();
        for c in 0..k {
            let mems: Vec<usize> = (0..n).filter(|&i| labels[i] == c as i32).collect();
            if mems.len() >= min_pts || k <= 2 {
                if !mems.is_empty() {
                    let mut mean = vec![0.0; d];
                    for &m in &mems {
                        for j in 0..d {
                            mean[j] += data.x[m][j];
                        }
                    }
                    for v in &mut mean {
                        *v /= mems.len() as f64;
                    }
                    next_centers.push(mean);
                    cluster_members.push(mems);
                }
            }
        }
        if next_centers.is_empty() {
            next_centers = init_kmeans_pp(&data.x, 2, &mut rng);
            centers = next_centers;
            continue;
        }
        centers = next_centers;

        // Split on odd iterations if k < 2*k_target
        if it % 2 == 1 && centers.len() < (k_target * 2).min(12) {
            let mut split_added = Vec::new();
            let mut remove_idx = None;
            for (c, mems) in cluster_members.iter().enumerate() {
                if mems.len() < 2 * min_pts {
                    continue;
                }
                let mut max_dim = 0;
                let mut max_std = 0.0;
                for j in 0..d {
                    let var = mems
                        .iter()
                        .map(|&m| (data.x[m][j] - centers[c][j]).powi(2))
                        .sum::<f64>()
                        / mems.len() as f64;
                    let s = var.sqrt();
                    if s > max_std {
                        max_std = s;
                        max_dim = j;
                    }
                }
                if max_std > std_split {
                    let mut c1 = centers[c].clone();
                    let mut c2 = centers[c].clone();
                    c1[max_dim] += 0.5 * max_std;
                    c2[max_dim] -= 0.5 * max_std;
                    split_added.push(c1);
                    split_added.push(c2);
                    remove_idx = Some(c);
                    splits += 1;
                    break;
                }
            }
            if let Some(ri) = remove_idx {
                centers.remove(ri);
                centers.extend(split_added);
            }
        } else if it % 2 == 0 && centers.len() > 2 {
            // Merge closest pair below threshold
            let mut best_pair = None;
            let mut best_d = merge_dist;
            for a in 0..centers.len() {
                for b in (a + 1)..centers.len() {
                    let dd = euclidean(&centers[a], &centers[b]);
                    if dd < best_d {
                        best_d = dd;
                        best_pair = Some((a, b));
                    }
                }
            }
            if let Some((a, b)) = best_pair {
                let merged: Vec<f64> = (0..d)
                    .map(|j| 0.5 * (centers[a][j] + centers[b][j]))
                    .collect();
                centers.remove(b);
                centers[a] = merged;
                merges += 1;
            }
        }
    }

    for i in 0..n {
        labels[i] = (0..centers.len())
            .min_by(|&a, &b| {
                euclidean_sq(&data.x[i], &centers[a])
                    .total_cmp(&euclidean_sq(&data.x[i], &centers[b]))
            })
            .unwrap_or(0) as i32;
    }

    AlgorithmOutput::new(
        labels,
        format!("mål_k={}, minPts={}", k_target, min_pts),
        format!("Endte på k={} ({} splitter, {} sammenslåinger)", centers.len(), splits, merges),
    )
}

// 19. Hartigan-Wong K-Means (exact point-by-point variance-change transfer)
fn hartigan_wong(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let max_iter = cfg.max_iters(50);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 19);

    let init_centers = init_kmeans_pp(&data.x, k, &mut rng);
    let (lbl0, mut centers, _, _) = run_lloyd_from_centers(&data.x, init_centers, 3);
    let mut labels: Vec<usize> = lbl0.into_iter().map(|l| l as usize).collect();
    let mut counts = vec![0usize; k];
    for &c in &labels {
        counts[c] += 1;
    }

    let mut transfers = 0usize;
    for _ in 0..max_iter {
        let mut changed = false;
        for i in 0..n {
            let cur = labels[i];
            if counts[cur] <= 1 {
                continue;
            }
            let na = counts[cur] as f64;
            let r_leave = (na / (na - 1.0)) * euclidean_sq(&data.x[i], &centers[cur]);

            let mut best_c = cur;
            let mut best_r = r_leave;
            for c in 0..k {
                if c == cur {
                    continue;
                }
                let nb = counts[c] as f64;
                let r_enter = (nb / (nb + 1.0)) * euclidean_sq(&data.x[i], &centers[c]);
                if r_enter + 1e-10 < best_r {
                    best_r = r_enter;
                    best_c = c;
                }
            }
            if best_c != cur {
                // Immediate incremental centroid update
                for j in 0..d {
                    centers[cur][j] = (centers[cur][j] * na - data.x[i][j]) / (na - 1.0);
                    let nb = counts[best_c] as f64;
                    centers[best_c][j] = (centers[best_c][j] * nb + data.x[i][j]) / (nb + 1.0);
                }
                counts[cur] -= 1;
                counts[best_c] += 1;
                labels[i] = best_c;
                transfers += 1;
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }

    AlgorithmOutput::new(
        labels.into_iter().map(|c| c as i32).collect(),
        format!("k={}, Hartigan-Wong", k),
        format!("{} direkte punkt-overføringer med umiddelbar oppdatering", transfers),
    )
}

// 20. COP-KMeans (Constrained K-Means with auto-derived high-confidence must-link / cannot-link constraints)
fn cop_kmeans(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let max_iter = cfg.max_iters(40);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 20);

    // Derive pairwise constraints: if ground truth exists, sample 8% semi-supervised pairs;
    // otherwise derive natural mutual-nearest-neighbor must-links & farthest-pair cannot-links.
    let mut must_link: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut cannot_link: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut n_ml = 0usize;
    let mut n_cl = 0usize;

    if let Some(gt) = &data.ground_truth {
        let n_pairs = (n / 5).clamp(4, 25);
        for _ in 0..n_pairs {
            let a = rng.gen_range(0..n);
            let b = rng.gen_range(0..n);
            if a == b {
                continue;
            }
            if gt[a] == gt[b] {
                must_link[a].push(b);
                must_link[b].push(a);
                n_ml += 1;
            } else {
                cannot_link[a].push(b);
                cannot_link[b].push(a);
                n_cl += 1;
            }
        }
    } else {
        // Mutual 1-NN must-link (top 8% tightest pairs) and farthest pairs cannot-link
        let mut pairs = Vec::new();
        for i in 0..n {
            if let Some(nn) = (0..n)
                .filter(|&j| j != i)
                .min_by(|&a, &b| data.dist_matrix[i][a].total_cmp(&data.dist_matrix[i][b]))
            {
                if i < nn {
                    pairs.push((i, nn, data.dist_matrix[i][nn]));
                }
            }
        }
        pairs.sort_by(|a, b| a.2.total_cmp(&b.2));
        for &(a, b, _) in pairs.iter().take((n / 12).clamp(2, 12)) {
            must_link[a].push(b);
            must_link[b].push(a);
            n_ml += 1;
        }
        let mut far = Vec::new();
        for i in 0..n {
            for j in (i + 1)..n {
                far.push((i, j, data.dist_matrix[i][j]));
            }
        }
        far.sort_by(|a, b| b.2.total_cmp(&a.2));
        for &(a, b, _) in far.iter().take((n / 15).clamp(2, 10)) {
            cannot_link[a].push(b);
            cannot_link[b].push(a);
            n_cl += 1;
        }
    }

    let mut centers = init_kmeans_pp(&data.x, k, &mut rng);
    let mut labels = vec![-1i32; n];

    for _ in 0..max_iter {
        let mut new_labels = vec![-1i32; n];
        for i in 0..n {
            let mut order: Vec<usize> = (0..k).collect();
            order.sort_by(|&a, &b| {
                euclidean_sq(&data.x[i], &centers[a])
                    .total_cmp(&euclidean_sq(&data.x[i], &centers[b]))
            });
            let mut assigned = order[0] as i32;
            for &c in &order {
                let ci = c as i32;
                let violates_ml = must_link[i]
                    .iter()
                    .any(|&m| new_labels[m] >= 0 && new_labels[m] != ci);
                let violates_cl = cannot_link[i].iter().any(|&m| new_labels[m] == ci);
                if !violates_ml && !violates_cl {
                    assigned = ci;
                    break;
                }
            }
            new_labels[i] = assigned;
        }
        if new_labels == labels {
            break;
        }
        labels = new_labels;
        let mut sums = vec![vec![0.0; d]; k];
        let mut cnts = vec![0usize; k];
        for i in 0..n {
            let c = labels[i] as usize;
            cnts[c] += 1;
            for j in 0..d {
                sums[c][j] += data.x[i][j];
            }
        }
        for c in 0..k {
            if cnts[c] > 0 {
                for j in 0..d {
                    centers[c][j] = sums[c][j] / cnts[c] as f64;
                }
            }
        }
    }

    AlgorithmOutput::new(
        labels,
        format!("k={}, ML={}, CL={}", k, n_ml, n_cl),
        format!("Semi-veiledet COP-KMeans med {} ML / {} CL begrensninger", n_ml, n_cl),
    )
}

// 21. MacQueen's K-Means (sequential / online immediate centroid update)
fn macqueen_kmeans(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let passes = cfg.max_iters(12);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 21);

    let mut order: Vec<usize> = (0..n).collect();
    order.shuffle(&mut rng);
    let mut centers: Vec<Vec<f64>> = order[..k].iter().map(|&i| data.x[i].clone()).collect();
    let mut counts = vec![1usize; k];
    let mut labels = vec![0i32; n];

    for _ in 0..passes {
        order.shuffle(&mut rng);
        for &i in &order {
            let best_c = (0..k)
                .min_by(|&a, &b| {
                    euclidean_sq(&data.x[i], &centers[a])
                        .total_cmp(&euclidean_sq(&data.x[i], &centers[b]))
                })
                .unwrap_or(0);
            counts[best_c] += 1;
            let eta = 1.0 / counts[best_c] as f64;
            for j in 0..d {
                centers[best_c][j] += eta * (data.x[i][j] - centers[best_c][j]);
            }
            labels[i] = best_c as i32;
        }
    }
    // Final crisp assignment
    for i in 0..n {
        labels[i] = (0..k)
            .min_by(|&a, &b| {
                euclidean_sq(&data.x[i], &centers[a])
                    .total_cmp(&euclidean_sq(&data.x[i], &centers[b]))
            })
            .unwrap_or(0) as i32;
    }
    AlgorithmOutput::new(
        labels,
        format!("k={}, online-pass={}", k, passes),
        format!("Sekvensiell oppdatering etter hvert datapunkt ({} gjennomløp)", passes),
    )
}

// 22. Forgy's K-Means (random data-point seeds + batch Lloyd)
fn forgy_kmeans(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let k = cfg.effective_k(n);
    let max_iter = cfg.max_iters(60);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 22);

    let mut idxs: Vec<usize> = (0..n).collect();
    idxs.shuffle(&mut rng);
    let centers: Vec<Vec<f64>> = idxs[..k].iter().map(|&i| data.x[i].clone()).collect();
    let (labels, _, sse, iters) = run_lloyd_from_centers(&data.x, centers, max_iter);
    AlgorithmOutput::new(
        labels,
        format!("k={}, Forgy-frø", k),
        format!("Batch-variant konvergerte på {} iter. (SSE={:.2})", iters, sse),
    )
}

// 23. K-Means|| (Bahmani et al. 2012 — scalable oversampling rounds + weighted re-clustering)
fn kmeans_parallel(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let k = cfg.effective_k(n);
    let oversample_l = (2 * k) as f64;
    let rounds = 5;
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 23);

    let mut c_idx = vec![rng.gen_range(0..n)];
    for _ in 0..rounds {
        let dists: Vec<f64> = (0..n)
            .map(|i| {
                c_idx
                    .iter()
                    .map(|&ci| data.dist_matrix[i][ci].powi(2))
                    .fold(f64::INFINITY, f64::min)
            })
            .collect();
        let phi: f64 = dists.iter().sum::<f64>().max(1e-12);
        for i in 0..n {
            let prob = (oversample_l * dists[i] / phi).min(1.0);
            if rng.gen::<f64>() < prob && !c_idx.contains(&i) {
                c_idx.push(i);
            }
        }
    }

    // Weight each candidate by number of points closest to it, then cluster candidates to k centers
    let m = c_idx.len();
    let mut weights = vec![0.0_f64; m];
    for i in 0..n {
        let best = (0..m)
            .min_by(|&a, &b| data.dist_matrix[i][c_idx[a]].total_cmp(&data.dist_matrix[i][c_idx[b]]))
            .unwrap_or(0);
        weights[best] += 1.0;
    }

    let mut cand_order: Vec<usize> = (0..m).collect();
    cand_order.sort_by(|&a, &b| weights[b].total_cmp(&weights[a]));
    let mut centers: Vec<Vec<f64>> = cand_order
        .into_iter()
        .take(k)
        .map(|ci| data.x[c_idx[ci]].clone())
        .collect();
    while centers.len() < k {
        centers.push(data.x[rng.gen_range(0..n)].clone());
    }

    let (labels, _, sse, iters) = run_lloyd_from_centers(&data.x, centers, cfg.max_iters(50));
    AlgorithmOutput::new(
        labels,
        format!("k={}, l={:.0}, runder={}", k, oversample_l, rounds),
        format!("Samplet {} kandidater -> {} iter., SSE={:.2}", m, iters, sse),
    )
}

// 24. FastPAM / FasterPAM (Schubert & Rousseeuw 2021 — eager O((n-k)^2) delta-TD swaps)
fn faster_pam(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let k = cfg.effective_k(n);
    let max_iter = cfg.max_iters(25);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 24);

    // Initialize medoids via K-Means++ style farthest seeding
    let mut medoids = vec![rng.gen_range(0..n)];
    while medoids.len() < k {
        let next = (0..n)
            .filter(|i| !medoids.contains(i))
            .max_by(|&a, &b| {
                let da = medoids.iter().map(|&m| data.dist_matrix[a][m]).fold(f64::INFINITY, f64::min);
                let db = medoids.iter().map(|&m| data.dist_matrix[b][m]).fold(f64::INFINITY, f64::min);
                da.total_cmp(&db)
            })
            .unwrap_or(0);
        medoids.push(next);
    }

    // Maintain nearest & second-nearest medoid per point
    let mut nearest = vec![0usize; n];
    let mut d_near = vec![0.0; n];
    let mut d_second = vec![0.0; n];
    let recompute_cache = |meds: &[usize], nearest: &mut [usize], d1: &mut [f64], d2: &mut [f64]| {
        for i in 0..n {
            let mut b1 = f64::INFINITY;
            let mut b2 = f64::INFINITY;
            let mut m1 = 0;
            for (mi, &m) in meds.iter().enumerate() {
                let dd = data.dist_matrix[i][m];
                if dd < b1 {
                    b2 = b1;
                    b1 = dd;
                    m1 = mi;
                } else if dd < b2 {
                    b2 = dd;
                }
            }
            nearest[i] = m1;
            d1[i] = b1;
            d2[i] = b2;
        }
    };
    recompute_cache(&medoids, &mut nearest, &mut d_near, &mut d_second);

    let mut swaps = 0usize;
    for _ in 0..max_iter {
        let mut improved = false;
        for cand in 0..n {
            if medoids.contains(&cand) {
                continue;
            }
            // Compute delta_TD(m_i, cand) for all i in O(n) simultaneously (FastPAM core insight!)
            let mut delta_remove = vec![0.0; k];
            let mut delta_add = 0.0;
            for o in 0..n {
                let d_oc = data.dist_matrix[o][cand];
                let mi = nearest[o];
                if d_oc < d_near[o] {
                    delta_add += d_oc - d_near[o];
                    delta_remove[mi] += d_near[o] - d_oc; // offset double counting if mi removed
                }
                delta_remove[mi] += d_oc.min(d_second[o]) - d_near[o];
            }
            // Find best medoid to replace eagerly (FasterPAM)
            let mut best_mi = 0;
            let mut best_delta = 0.0;
            for mi in 0..k {
                let total_delta = delta_add + delta_remove[mi];
                if total_delta < best_delta - 1e-9 {
                    best_delta = total_delta;
                    best_mi = mi;
                }
            }
            if best_delta < -1e-9 {
                medoids[best_mi] = cand;
                recompute_cache(&medoids, &mut nearest, &mut d_near, &mut d_second);
                swaps += 1;
                improved = true;
            }
        }
        if !improved {
            break;
        }
    }

    let total_td: f64 = d_near.iter().sum();
    let labels: Vec<i32> = nearest.into_iter().map(|c| c as i32).collect();
    AlgorithmOutput::new(
        labels,
        format!("k={}, FasterPAM O((n-k)²)", k),
        format!("Utførte {} raske delta-TD-bytter (TD={:.2})", swaps, total_td),
    )
}

// 25. CLARA-MED (Median-cost objective variant of CLARA for outlier robustness)
fn clara_med(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let k = cfg.effective_k(n);
    let sample_size = (40 + 2 * k).min(n);
    let runs = 6;
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 25);
    let mut all_idx: Vec<usize> = (0..n).collect();

    let mut best_labels = vec![0i32; n];
    let mut best_med_cost = f64::INFINITY;

    for _ in 0..runs {
        all_idx.shuffle(&mut rng);
        let sub = &all_idx[..sample_size];
        let (meds, _) = pam_on_indices(&data.dist_matrix, sub, k, 8);
        let (labels, _) = assign_to_medoids(&data.dist_matrix, &meds);
        let mut pt_dists: Vec<f64> = (0..n)
            .map(|i| data.dist_matrix[i][meds[labels[i] as usize]])
            .collect();
        let med_cost = median_of_slice(&mut pt_dists);
        if med_cost < best_med_cost {
            best_med_cost = med_cost;
            best_labels = labels;
        }
    }
    AlgorithmOutput::new(
        best_labels,
        format!("k={}, utvalg={}, median-kriterium", k, sample_size),
        format!("Laveste medianavstand til medoid: {:.3}", best_med_cost),
    )
}

// 26. Trimmed K-Means (Cuesta-Albertos et al. 1997 — excludes alpha fraction of outliers as noise -1)
fn trimmed_kmeans(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let alpha = 0.08;
    let trim_count = ((n as f64) * alpha).floor() as usize;
    let max_iter = cfg.max_iters(50);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 26);

    let mut centers = init_kmeans_pp(&data.x, k, &mut rng);
    let mut labels = vec![0i32; n];

    for _ in 0..max_iter {
        let mut pt_best: Vec<(usize, usize, f64)> = (0..n)
            .map(|i| {
                let mut bc = 0;
                let mut bd = f64::INFINITY;
                for (c, cen) in centers.iter().enumerate() {
                    let dd = euclidean_sq(&data.x[i], cen);
                    if dd < bd {
                        bd = dd;
                        bc = c;
                    }
                }
                (i, bc, bd)
            })
            .collect();
        pt_best.sort_by(|a, b| a.2.total_cmp(&b.2));

        let keep_n = n.saturating_sub(trim_count).max(k);
        let mut new_labels = vec![-1i32; n];
        let mut sums = vec![vec![0.0; d]; k];
        let mut cnts = vec![0usize; k];

        for &(i, c, _) in pt_best.iter().take(keep_n) {
            new_labels[i] = c as i32;
            cnts[c] += 1;
            for j in 0..d {
                sums[c][j] += data.x[i][j];
            }
        }
        if new_labels == labels {
            break;
        }
        labels = new_labels;
        for c in 0..k {
            if cnts[c] > 0 {
                for j in 0..d {
                    centers[c][j] = sums[c][j] / cnts[c] as f64;
                }
            }
        }
    }

    AlgorithmOutput::new(
        labels,
        format!("k={}, alpha={:.0}% trim", k, alpha * 100.0),
        format!("Trimmet {} mest avvikende punkter som støy", trim_count),
    )
}

// 27. Power K-Means (Xu & Lange 2019 — annealing power mean s -> -infinity via MM updates)
fn power_kmeans(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let max_iter = cfg.max_iters(50);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 27);

    let mut centers = init_kmeans_pp(&data.x, k, &mut rng);
    let mut s = -1.0_f64;

    for it in 0..max_iter {
        let mut num = vec![vec![0.0; d]; k];
        let mut den = vec![0.0; k];
        for i in 0..n {
            let d2s: Vec<f64> = (0..k)
                .map(|c| euclidean_sq(&data.x[i], &centers[c]).max(1e-8))
                .collect();
            // Compute soft power weights w_ic = d_ic^(2*(s-1)) / (sum_j d_ij^(2*s))^(1 - 1/s)
            // Numerically stabilized using min distance
            let min_d = d2s.iter().copied().fold(f64::INFINITY, f64::min);
            let scaled: Vec<f64> = d2s.iter().map(|&v| (v / min_d).powf(s)).collect();
            let sum_scaled: f64 = scaled.iter().sum::<f64>().max(1e-12);
            let mean_term = (sum_scaled / k as f64).powf(1.0 / s - 1.0);
            for c in 0..k {
                let w = ((d2s[c] / min_d).powf(s - 1.0) * mean_term / k as f64).min(1e6);
                den[c] += w;
                for j in 0..d {
                    num[c][j] += w * data.x[i][j];
                }
            }
        }
        for c in 0..k {
            if den[c] > 1e-12 {
                for j in 0..d {
                    centers[c][j] = num[c][j] / den[c];
                }
            }
        }
        if it % 4 == 3 {
            s = (s * 1.25).max(-20.0);
        }
    }

    let labels: Vec<i32> = (0..n)
        .map(|i| {
            (0..k)
                .min_by(|&a, &b| {
                    euclidean_sq(&data.x[i], &centers[a])
                        .total_cmp(&euclidean_sq(&data.x[i], &centers[b]))
                })
                .unwrap_or(0) as i32
        })
        .collect();
    AlgorithmOutput::new(
        labels,
        format!("k={}, s=-1 -> {:.1}", k, s),
        "Annealing med potensgjennomsnitt (MM-optimering)".to_string(),
    )
}

// 28. W-K-Means (Huang et al. 2005 — automatic variable weight learning with beta exponent)
fn weighted_kmeans(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let beta = 2.0_f64;
    let max_iter = cfg.max_iters(50);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 28);

    let mut centers = init_kmeans_pp(&data.x, k, &mut rng);
    let mut weights = vec![1.0 / d as f64; d];
    let mut labels = vec![0i32; n];

    for it in 0..max_iter {
        let mut changed = false;
        for i in 0..n {
            let mut best_c = 0;
            let mut best_d = f64::INFINITY;
            for c in 0..k {
                let dist: f64 = (0..d)
                    .map(|j| weights[j].powf(beta) * (data.x[i][j] - centers[c][j]).powi(2))
                    .sum();
                if dist < best_d {
                    best_d = dist;
                    best_c = c;
                }
            }
            if labels[i] != best_c as i32 {
                labels[i] = best_c as i32;
                changed = true;
            }
        }
        // Update centers
        let mut sums = vec![vec![0.0; d]; k];
        let mut cnts = vec![0usize; k];
        for i in 0..n {
            let c = labels[i] as usize;
            cnts[c] += 1;
            for j in 0..d {
                sums[c][j] += data.x[i][j];
            }
        }
        for c in 0..k {
            if cnts[c] > 0 {
                for j in 0..d {
                    centers[c][j] = sums[c][j] / cnts[c] as f64;
                }
            }
        }
        // Update weights inversely proportional to within-cluster dispersion D_j
        let mut disp = vec![1e-6; d];
        for i in 0..n {
            let c = labels[i] as usize;
            for j in 0..d {
                disp[j] += (data.x[i][j] - centers[c][j]).powi(2);
            }
        }
        let inv_exp = 1.0 / (beta - 1.0);
        let raw: Vec<f64> = disp.iter().map(|&dj| (1.0 / dj).powf(inv_exp)).collect();
        let sum_raw: f64 = raw.iter().sum::<f64>().max(1e-12);
        for j in 0..d {
            weights[j] = raw[j] / sum_raw;
        }
        if !changed && it > 2 {
            break;
        }
    }

    let top_feat = (0..d)
        .max_by(|&a, &b| weights[a].total_cmp(&weights[b]))
        .unwrap_or(0);
    AlgorithmOutput::new(
        labels,
        format!("k={}, beta={:.1}", k, beta),
        format!(
            "Høyeste lærte vekt: {} ({:.2})",
            data.feature_names.get(top_feat).map(|s| s.as_str()).unwrap_or("x0"),
            weights[top_feat]
        ),
    )
}

// 29. Minkowski Weighted K-Means (de Amorim & Mirkin 2012 — cluster-specific Minkowski weights)
fn minkowski_weighted_kmeans(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let p = 1.5_f64;
    let max_iter = cfg.max_iters(45);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 29);

    let mut centers = init_kmeans_pp(&data.x, k, &mut rng);
    let mut weights = vec![vec![1.0 / d as f64; d]; k];
    let mut labels = vec![0i32; n];

    for it in 0..max_iter {
        let mut changed = false;
        for i in 0..n {
            let mut best_c = 0;
            let mut best_d = f64::INFINITY;
            for c in 0..k {
                let dist = weighted_minkowski(&data.x[i], &centers[c], &weights[c], p);
                if dist < best_d {
                    best_d = dist;
                    best_c = c;
                }
            }
            if labels[i] != best_c as i32 {
                labels[i] = best_c as i32;
                changed = true;
            }
        }
        // Update centers (steepest descent / weighted median-mean blend for L_p) & per-cluster weights
        for c in 0..k {
            let mems: Vec<usize> = (0..n).filter(|&i| labels[i] == c as i32).collect();
            if mems.is_empty() {
                continue;
            }
            let mut disp = vec![1e-4; d];
            for j in 0..d {
                let mu = mems.iter().map(|&m| data.x[m][j]).sum::<f64>() / mems.len() as f64;
                let mut vals: Vec<f64> = mems.iter().map(|&m| data.x[m][j]).collect();
                let med = median_of_slice(&mut vals);
                centers[c][j] = 0.5 * (mu + med);
                disp[j] += mems
                    .iter()
                    .map(|&m| (data.x[m][j] - centers[c][j]).abs().powf(p))
                    .sum::<f64>();
            }
            let inv_exp = 1.0 / (p - 1.0);
            let raw: Vec<f64> = disp.iter().map(|&v| (1.0 / v).powf(inv_exp)).collect();
            let s: f64 = raw.iter().sum::<f64>().max(1e-12);
            for j in 0..d {
                weights[c][j] = raw[j] / s;
            }
        }
        if !changed && it > 2 {
            break;
        }
    }

    AlgorithmOutput::new(
        labels,
        format!("k={}, Minkowski p={:.1}", k, p),
        "Lærte klyngespesifikke Minkowski-attributtvekter".to_string(),
    )
}

// 30. LBG-algoritmen (Linde-Buzo-Gray 1980 — iterative codebook splitting by +/- epsilon perturbation)
fn lbg_vector_quantization(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let eps = 0.02;

    // Start with 1 codebook vector = dataset centroid
    let mut mean = vec![0.0; d];
    for row in &data.x {
        for j in 0..d {
            mean[j] += row[j];
        }
    }
    for v in &mut mean {
        *v /= n as f64;
    }

    // Principal direction for clean deterministic perturbation
    let mut cov = vec![vec![0.0; d]; d];
    for row in &data.x {
        for a in 0..d {
            for b in a..d {
                cov[a][b] += (row[a] - mean[a]) * (row[b] - mean[b]) / n as f64;
            }
        }
    }
    for a in 0..d {
        for b in a..d {
            cov[b][a] = cov[a][b];
        }
    }
    let (_, eigvecs) = jacobi_eigen(&cov, 20);
    let dir = eigvecs.first().cloned().unwrap_or_else(|| vec![1.0; d]);

    let mut codebook = vec![mean];
    let mut split_rounds = 0;

    while codebook.len() < k {
        split_rounds += 1;
        let mut expanded = Vec::with_capacity((codebook.len() * 2).min(k));
        for cw in &codebook {
            if expanded.len() + (codebook.len() - expanded.len() / 2) < k || expanded.len() + 2 <= k {
                let c_plus: Vec<f64> = (0..d).map(|j| cw[j] + eps * dir[j]).collect();
                let c_minus: Vec<f64> = (0..d).map(|j| cw[j] - eps * dir[j]).collect();
                expanded.push(c_plus);
                if expanded.len() < k {
                    expanded.push(c_minus);
                }
            } else {
                expanded.push(cw.clone());
            }
        }
        expanded.truncate(k);
        let (_, refined, _, _) = run_lloyd_from_centers(&data.x, expanded, 25);
        codebook = refined;
    }

    let (labels, _, distortion, _) = run_lloyd_from_centers(&data.x, codebook, 30);
    AlgorithmOutput::new(
        labels,
        format!("kodebok={}, eps={}", k, eps),
        format!("{} kodebok-splitt, kvantiseringsfeil={:.2}", split_rounds, distortion),
    )
}
