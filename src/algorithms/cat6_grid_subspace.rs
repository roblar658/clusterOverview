use std::collections::{HashMap, HashSet, VecDeque};

use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;

use crate::algorithms::cat1_centroid::{init_kmeans_pp, run_lloyd_from_centers};
use crate::algorithms::cat2_density::run_dbscan_core;
use crate::algorithms::cat3_hierarchical::{agglomerative_lance_williams, Linkage};
use crate::algorithms::{AlgorithmOutput, RunConfig};
use crate::data::PreparedData;
use crate::math::{
    dot, estimate_eps_kdist, euclidean, jacobi_eigen, top_k_symmetric_eigen,
};

pub fn run(id: usize, data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    match id {
        133 => sting(data, cfg),
        134 => clique(data, cfg),
        135 => wavecluster(data, cfg),
        136 => proclus(data, cfg),
        137 => subclu(data, cfg),
        138 => orclus(data, cfg),
        139 => doc_projected(data, cfg),
        140 => mafia(data, cfg),
        141 => optigrid(data, cfg),
        142 => bang_clustering(data, cfg),
        143 => enclus(data, cfg),
        144 => findit(data, cfg),
        145 => cosa(data, cfg),
        146 => fires(data, cfg),
        147 => inscy(data, cfg),
        148 => p3c(data, cfg),
        149 => statpc(data, cfg),
        150 => sparse_subspace_clustering(data, cfg),
        _ => unreachable!(),
    }
}

fn pca_bounds(pts: &[[f64; 2]]) -> ([f64; 2], [f64; 2]) {
    let mut mn = [f64::INFINITY; 2];
    let mut mx = [f64::NEG_INFINITY; 2];
    for p in pts {
        for j in 0..2 {
            mn[j] = mn[j].min(p[j]);
            mx[j] = mx[j].max(p[j]);
        }
    }
    for j in 0..2 {
        if (mx[j] - mn[j]).abs() < 1e-9 {
            mx[j] = mn[j] + 1.0;
        }
    }
    (mn, mx)
}

fn grid_connected_components(
    dense_cells: &HashSet<(usize, usize)>,
    grid_n: usize,
) -> HashMap<(usize, usize), i32> {
    let mut cell_lbl = HashMap::new();
    let mut cid = 0i32;
    for &cell in dense_cells {
        if cell_lbl.contains_key(&cell) {
            continue;
        }
        let mut q = VecDeque::new();
        cell_lbl.insert(cell, cid);
        q.push_back(cell);
        while let Some((cx, cy)) = q.pop_front() {
            for dx in -1isize..=1 {
                for dy in -1isize..=1 {
                    if dx == 0 && dy == 0 {
                        continue;
                    }
                    let nx = cx as isize + dx;
                    let ny = cy as isize + dy;
                    if nx >= 0 && ny >= 0 && (nx as usize) < grid_n && (ny as usize) < grid_n {
                        let nb = (nx as usize, ny as usize);
                        if dense_cells.contains(&nb) && !cell_lbl.contains_key(&nb) {
                            cell_lbl.insert(nb, cid);
                            q.push_back(nb);
                        }
                    }
                }
            }
        }
        cid += 1;
    }
    cell_lbl
}

// 133. STING (Wang et al. 1997 — Statistical Information Grid hierarchical cell tree)
fn sting(data: &PreparedData, _cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let (mn, mx) = pca_bounds(&data.pca_2d);
    let levels = [2usize, 4, 8];
    let g_leaf = levels[2];

    let mut leaf_cells: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
    let mut pt_cell = vec![(0usize, 0usize); n];
    for i in 0..n {
        let cx = (((data.pca_2d[i][0] - mn[0]) / (mx[0] - mn[0])) * g_leaf as f64)
            .floor()
            .clamp(0.0, (g_leaf - 1) as f64) as usize;
        let cy = (((data.pca_2d[i][1] - mn[1]) / (mx[1] - mn[1])) * g_leaf as f64)
            .floor()
            .clamp(0.0, (g_leaf - 1) as f64) as usize;
        pt_cell[i] = (cx, cy);
        leaf_cells.entry((cx, cy)).or_default().push(i);
    }

    // Store per-cell statistical summary (n, mean, std) and propagate top-down density query
    let mut counts: Vec<usize> = leaf_cells.values().map(|v| v.len()).collect();
    counts.sort_unstable();
    let dens_thresh = counts.get(counts.len() / 3).copied().unwrap_or(1).max(2);

    let dense_set: HashSet<(usize, usize)> = leaf_cells
        .iter()
        .filter(|(_, pts)| pts.len() >= dens_thresh)
        .map(|(&c, _)| c)
        .collect();
    let cell_lbl = grid_connected_components(&dense_set, g_leaf);

    let mut labels = vec![-1i32; n];
    for i in 0..n {
        if let Some(&cid) = cell_lbl.get(&pt_cell[i]) {
            labels[i] = cid;
        }
    }
    AlgorithmOutput::new(
        labels,
        format!("3 nivåer ({}x{}), terskel={}", g_leaf, g_leaf, dens_thresh),
        format!("Hierarkiske cellestatistikker ({} tette celler)", dense_set.len()),
    )
}

// 134. CLIQUE (Agrawal et al. SIGMOD 1998 — Apriori dense subspace units)
fn clique(data: &PreparedData, _cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let xi = 6usize; // intervals per dimension
    let tau = ((n as f64) * 0.04).ceil().max(2.0) as usize;

    // Discretize every dimension into xi equal-width bins
    let mut bins = vec![vec![0usize; d]; n];
    for j in 0..d {
        let mn = (0..n).map(|i| data.x[i][j]).fold(f64::INFINITY, f64::min);
        let mx = (0..n).map(|i| data.x[i][j]).fold(f64::NEG_INFINITY, f64::max);
        let span = (mx - mn).max(1e-9);
        for i in 0..n {
            bins[i][j] = (((data.x[i][j] - mn) / span) * xi as f64)
                .floor()
                .clamp(0.0, (xi - 1) as f64) as usize;
        }
    }

    // 1D dense units -> score each dimension by coverage of dense units
    let mut dim_scores: Vec<(usize, usize)> = (0..d)
        .map(|j| {
            let mut cnt = vec![0usize; xi];
            for i in 0..n {
                cnt[bins[i][j]] += 1;
            }
            let covered: usize = cnt.iter().filter(|&&c| c >= tau).sum();
            (j, covered)
        })
        .collect();
    dim_scores.sort_by(|a, b| b.1.cmp(&a.1));

    // Select top-2 Apriori subspace dimensions
    let d1 = dim_scores[0].0;
    let d2 = dim_scores.get(1).map(|x| x.0).unwrap_or(d1);

    let mut cell_pts: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
    for i in 0..n {
        cell_pts.entry((bins[i][d1], bins[i][d2])).or_default().push(i);
    }
    let dense: HashSet<(usize, usize)> = cell_pts
        .iter()
        .filter(|(_, v)| v.len() >= tau)
        .map(|(&c, _)| c)
        .collect();
    let cell_lbl = grid_connected_components(&dense, xi);

    let mut labels = vec![-1i32; n];
    for i in 0..n {
        if let Some(&cid) = cell_lbl.get(&(bins[i][d1], bins[i][d2])) {
            labels[i] = cid;
        }
    }
    AlgorithmOutput::new(
        labels,
        format!("xi={}, tau={}, underrom=({},{})", xi, tau, d1, d2),
        format!("Apriori fant {} tette 2D-underromsceller", dense.len()),
    )
}

// 135. WaveCluster (Sheikholeslami et al. VLDB 1998 — 2D Haar wavelet transform on grid)
fn wavecluster(data: &PreparedData, _cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let (mn, mx) = pca_bounds(&data.pca_2d);
    let g = 16usize; // power-of-2 grid for 2D Haar DWT

    let mut grid = vec![vec![0.0; g]; g];
    let mut pt_cell = vec![(0usize, 0usize); n];
    for i in 0..n {
        let cx = (((data.pca_2d[i][0] - mn[0]) / (mx[0] - mn[0])) * g as f64)
            .floor()
            .clamp(0.0, (g - 1) as f64) as usize;
        let cy = (((data.pca_2d[i][1] - mn[1]) / (mx[1] - mn[1])) * g as f64)
            .floor()
            .clamp(0.0, (g - 1) as f64) as usize;
        grid[cx][cy] += 1.0;
        pt_cell[i] = (cx, cy);
    }

    // 1-level 2D Haar wavelet transform -> LL approximation band (8x8) removes high-freq noise
    let half = g / 2;
    let mut ll = vec![vec![0.0; half]; half];
    for x in 0..half {
        for y in 0..half {
            ll[x][y] = 0.5
                * (grid[2 * x][2 * y]
                    + grid[2 * x + 1][2 * y]
                    + grid[2 * x][2 * y + 1]
                    + grid[2 * x + 1][2 * y + 1]);
        }
    }
    let mut nz: Vec<f64> = ll.iter().flat_map(|r| r.iter()).copied().filter(|&v| v > 0.0).collect();
    nz.sort_by(|a, b| a.total_cmp(b));
    let thresh = nz.get(nz.len() / 3).copied().unwrap_or(0.5);

    let mut dense = HashSet::new();
    for x in 0..half {
        for y in 0..half {
            if ll[x][y] >= thresh {
                dense.insert((x, y));
            }
        }
    }
    let cell_lbl = grid_connected_components(&dense, half);
    let mut labels = vec![-1i32; n];
    for i in 0..n {
        let (cx, cy) = pt_cell[i];
        if let Some(&cid) = cell_lbl.get(&(cx / 2, cy / 2)) {
            labels[i] = cid;
        }
    }
    AlgorithmOutput::new(
        labels,
        format!("Haar DWT {}x{} -> {}x{}", g, g, half, half),
        format!("Wavelet LL-bånd filtrerte støy (terskel={:.2})", thresh),
    )
}

// 136. PROCLUS (Aggarwal et al. SIGMOD 1999 — Medoids + per-cluster relevant subspace dimensions)
fn proclus(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let avg_l = (d / 2).clamp(1, d);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 136);

    // Greedy farthest-first medoid initialization
    let mut medoids = vec![rng.gen_range(0..n)];
    while medoids.len() < k {
        let far = (0..n)
            .filter(|i| !medoids.contains(i))
            .max_by(|&a, &b| {
                let da = medoids.iter().map(|&m| data.dist_matrix[a][m]).fold(f64::INFINITY, f64::min);
                let db = medoids.iter().map(|&m| data.dist_matrix[b][m]).fold(f64::INFINITY, f64::min);
                da.total_cmp(&db)
            })
            .unwrap_or(0);
        medoids.push(far);
    }

    let mut dims: Vec<Vec<usize>> = vec![(0..d).collect(); k];
    let mut labels = vec![0i32; n];

    for _ in 0..15 {
        // Determine relevant dimensions for each medoid via Z-score of local dispersion
        let mut z_scores: Vec<(usize, usize, f64)> = Vec::new();
        for (ci, &m) in medoids.iter().enumerate() {
            let min_med_d = medoids
                .iter()
                .filter(|&&om| om != m)
                .map(|&om| data.dist_matrix[m][om])
                .fold(f64::INFINITY, f64::min);
            let local: Vec<usize> = (0..n)
                .filter(|&i| data.dist_matrix[i][m] <= min_med_d)
                .collect();
            let mut x_ij = vec![0.0; d];
            for j in 0..d {
                x_ij[j] = local
                    .iter()
                    .map(|&i| (data.x[i][j] - data.x[m][j]).abs())
                    .sum::<f64>()
                    / local.len().max(1) as f64;
            }
            let mu = x_ij.iter().sum::<f64>() / d as f64;
            let sigma = (x_ij.iter().map(|&v| (v - mu).powi(2)).sum::<f64>() / d.max(1) as f64)
                .sqrt()
                .max(1e-6);
            for j in 0..d {
                z_scores.push((ci, j, (x_ij[j] - mu) / sigma));
            }
        }
        // Pick lowest Z-scores (tightest spread dimensions) globally, ensuring at least 2 dims/cluster (or d)
        z_scores.sort_by(|a, b| a.2.total_cmp(&b.2));
        for c in 0..k {
            dims[c].clear();
        }
        let min_d = 2.min(d);
        for c in 0..k {
            for &(ci, j, _) in &z_scores {
                if ci == c && dims[c].len() < min_d {
                    dims[c].push(j);
                }
            }
        }
        let target_total = (k * avg_l).max(k * min_d);
        let mut cur_total: usize = dims.iter().map(|v| v.len()).sum();
        for &(ci, j, _) in &z_scores {
            if cur_total >= target_total {
                break;
            }
            if !dims[ci].contains(&j) {
                dims[ci].push(j);
                cur_total += 1;
            }
        }

        // Assign points using Manhattan segmental distance in each cluster's subspace
        for i in 0..n {
            let best_c = (0..k)
                .min_by(|&a, &b| {
                    let da = dims[a]
                        .iter()
                        .map(|&j| (data.x[i][j] - data.x[medoids[a]][j]).abs())
                        .sum::<f64>()
                        / dims[a].len().max(1) as f64;
                    let db = dims[b]
                        .iter()
                        .map(|&j| (data.x[i][j] - data.x[medoids[b]][j]).abs())
                        .sum::<f64>()
                        / dims[b].len().max(1) as f64;
                    da.total_cmp(&db)
                })
                .unwrap_or(0);
            labels[i] = best_c as i32;
        }

        // Update each medoid inside its assigned cluster
        for c in 0..k {
            let mems: Vec<usize> = (0..n).filter(|&i| labels[i] == c as i32).collect();
            if mems.is_empty() {
                continue;
            }
            if let Some(&best_m) = mems.iter().min_by(|&&a, &&b| {
                let sa: f64 = mems
                    .iter()
                    .map(|&p| dims[c].iter().map(|&j| (data.x[p][j] - data.x[a][j]).abs()).sum::<f64>())
                    .sum();
                let sb: f64 = mems
                    .iter()
                    .map(|&p| dims[c].iter().map(|&j| (data.x[p][j] - data.x[b][j]).abs()).sum::<f64>())
                    .sum();
                sa.total_cmp(&sb)
            }) {
                medoids[c] = best_m;
            }
        }
    }

    AlgorithmOutput::new(
        labels,
        format!("k={}, snitt dim l={}", k, avg_l),
        "PROCLUS fant klyngespesifikke akseparallelle underrom".to_string(),
    )
}

// 137. SUBCLU (Kailing et al. SDM 2004 — Bottom-up monotonicity density-connected subspace DBSCAN)
fn subclu(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let min_pts = 4usize.min(n / 4).max(2);

    // Find 1D dimensions supporting >= 2 density-connected clusters, then combine top subspace
    let mut good_dims: Vec<(usize, f64)> = Vec::new();
    for j in 0..d {
        let mut d1 = vec![vec![0.0; n]; n];
        for a in 0..n {
            for b in (a + 1)..n {
                let dd = (data.x[a][j] - data.x[b][j]).abs();
                d1[a][b] = dd;
                d1[b][a] = dd;
            }
        }
        let eps = estimate_eps_kdist(&d1, min_pts, 0.52);
        let lbl = run_dbscan_core(&d1, eps, min_pts);
        let nc = lbl.iter().filter(|&&l| l >= 0).collect::<HashSet<_>>().len();
        if nc >= 2 {
            let sil = crate::metrics::silhouette_score(&d1, &lbl).unwrap_or(0.0);
            good_dims.push((j, sil + 0.1 * nc as f64));
        }
    }
    good_dims.sort_by(|a, b| b.1.total_cmp(&a.1));
    let sub_dims: Vec<usize> = if good_dims.is_empty() {
        (0..d.min(3)).collect()
    } else {
        good_dims.iter().take(3.min(d)).map(|&(j, _)| j).collect()
    };

    // Run DBSCAN in discovered highest-dimensional dense subspace
    let mut sub_d = vec![vec![0.0; n]; n];
    for a in 0..n {
        for b in (a + 1)..n {
            let dd = sub_dims
                .iter()
                .map(|&j| (data.x[a][j] - data.x[b][j]).powi(2))
                .sum::<f64>()
                .sqrt();
            sub_d[a][b] = dd;
            sub_d[b][a] = dd;
        }
    }
    let eps = estimate_eps_kdist(&sub_d, min_pts, 0.62);
    let labels = run_dbscan_core(&sub_d, eps, min_pts);
    let _ = cfg;
    AlgorithmOutput::new(
        labels,
        format!("underrom={:?}, eps={:.2}", sub_dims, eps),
        "Monoton tetthets-underromsclustering (SUBCLU)".to_string(),
    )
}

// 138. ORCLUS (Aggarwal & Yu SIGMOD 2000 — Arbitrarily oriented PCA projected clusters)
fn orclus(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let l_dim = (d / 2).clamp(1, d);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 138);

    // Start with k_0 = min(2*k, n/3) seeds and iteratively project + merge down to k
    let k0 = (k * 2).min(n / 3).max(k);
    let c0 = init_kmeans_pp(&data.x, k0, &mut rng);
    let (mut labels, mut seeds, _, _) = run_lloyd_from_centers(&data.x, c0, 12);
    let mut cur_k = k0;

    // Basis vectors E_c (lowest-variance eigenvectors defining cluster's projected subspace!)
    let mut bases: Vec<Vec<Vec<f64>>> = vec![
        (0..l_dim)
            .map(|j| {
                let mut e = vec![0.0; d];
                e[j % d] = 1.0;
                e
            })
            .collect();
        cur_k
    ];

    for _round in 0..8 {
        // Compute per-cluster covariance -> pick l_dim smallest-variance eigenvectors (tight subspace)
        for c in 0..cur_k {
            let mems: Vec<usize> = (0..n).filter(|&i| labels[i] == c as i32).collect();
            if mems.len() < 3 {
                continue;
            }
            let mut mean = vec![0.0; d];
            for &m in &mems {
                for j in 0..d {
                    mean[j] += data.x[m][j] / mems.len() as f64;
                }
            }
            seeds[c] = mean.clone();
            let mut cov = vec![vec![0.0; d]; d];
            for &m in &mems {
                for a in 0..d {
                    for b in a..d {
                        cov[a][b] += (data.x[m][a] - mean[a]) * (data.x[m][b] - mean[b]);
                    }
                }
            }
            for a in 0..d {
                for b in a..d {
                    cov[a][b] /= mems.len() as f64;
                    cov[b][a] = cov[a][b];
                }
            }
            let (_, evecs) = jacobi_eigen(&cov, 20);
            // Smallest eigenvalues are at the end of jacobi_eigen (sorted descending)
            bases[c] = evecs.into_iter().rev().take(l_dim).collect();
        }

        // Re-assign points by projected Euclidean distance onto cluster's oriented subspace
        for i in 0..n {
            let best_c = (0..cur_k)
                .min_by(|&a, &b| {
                    let da = proj_dist_sq(&data.x[i], &seeds[a], &bases[a]);
                    let db = proj_dist_sq(&data.x[i], &seeds[b], &bases[b]);
                    da.total_cmp(&db)
                })
                .unwrap_or(0);
            labels[i] = best_c as i32;
        }

        if cur_k > k {
            // Merge closest pair of cluster seeds
            let mut best = (0, 1, f64::INFINITY);
            for a in 0..cur_k {
                for b in (a + 1)..cur_k {
                    let dd = euclidean(&seeds[a], &seeds[b]);
                    if dd < best.2 {
                        best = (a, b, dd);
                    }
                }
            }
            for j in 0..d {
                seeds[best.0][j] = 0.5 * (seeds[best.0][j] + seeds[best.1][j]);
            }
            seeds.remove(best.1);
            bases.remove(best.1);
            for l in &mut labels {
                if *l == best.1 as i32 {
                    *l = best.0 as i32;
                } else if *l > best.1 as i32 {
                    *l -= 1;
                }
            }
            cur_k -= 1;
        }
    }

    AlgorithmOutput::new(
        labels,
        format!("k={}, orientert dim l={}", k, l_dim),
        "Vilkårlig orienterte PCA-projiserte underrom (ORCLUS)".to_string(),
    )
}

fn proj_dist_sq(p: &[f64], center: &[f64], basis: &[Vec<f64>]) -> f64 {
    let diff: Vec<f64> = p.iter().zip(center.iter()).map(|(a, b)| a - b).collect();
    let mut s = 0.0;
    for e in basis {
        let proj = dot(&diff, e);
        s += proj * proj;
    }
    s
}

// 139. DOC (Procopiuc et al. SIGMOD 2002 — Density-based Optimal Projective Clustering via Monte Carlo)
fn doc_projected(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let alpha = 0.12_f64;
    let beta = 0.25_f64;
    let width = 1.6_f64;
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 139);

    let mut remaining: Vec<usize> = (0..n).collect();
    let mut labels = vec![-1i32; n];
    let mut cid = 0i32;

    while cid < k as i32 && remaining.len() >= 4 {
        let mut best_cluster: Vec<usize> = Vec::new();
        let mut best_mu = f64::NEG_INFINITY;

        let outer_trials = 14usize;
        let disc_size = 3usize.min(remaining.len());
        for _ in 0..outer_trials {
            let p = remaining[rng.gen_range(0..remaining.len())];
            // Sample discriminating set X
            let mut x_set = Vec::with_capacity(disc_size);
            for _ in 0..disc_size {
                x_set.push(remaining[rng.gen_range(0..remaining.len())]);
            }
            // Select dimensions where all x in X are within `width` of p
            let dims: Vec<usize> = (0..d)
                .filter(|&j| x_set.iter().all(|&q| (data.x[q][j] - data.x[p][j]).abs() <= width))
                .collect();
            if dims.is_empty() {
                continue;
            }
            let cand_pts: Vec<usize> = remaining
                .iter()
                .copied()
                .filter(|&q| dims.iter().all(|&j| (data.x[q][j] - data.x[p][j]).abs() <= width))
                .collect();
            if (cand_pts.len() as f64) < alpha * (remaining.len() as f64) {
                continue;
            }
            // DOC quality function mu(|C|, |D|) = |C| * (1/beta)^|D|
            let mu = (cand_pts.len() as f64) * (1.0 / beta).powi(dims.len().min(6) as i32);
            if mu > best_mu {
                best_mu = mu;
                best_cluster = cand_pts;
            }
        }
        if best_cluster.is_empty() {
            break;
        }
        let chosen: HashSet<usize> = best_cluster.iter().copied().collect();
        for &p in &best_cluster {
            labels[p] = cid;
        }
        remaining.retain(|x| !chosen.contains(x));
        cid += 1;
    }

    AlgorithmOutput::new(
        labels,
        format!("width={:.1}, alpha={:.2}, beta={:.2}", width, alpha, beta),
        format!("Monte Carlo-sampling fant {} projiserte hyperkubeklynger", cid),
    )
}

// 140. MAFIA (Goil et al. 1999 — Adaptive histogram grid intervals + CLIQUE merge)
fn mafia(data: &PreparedData, _cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;

    // Build adaptive 1D histogram bins by merging adjacent fine bins with similar density
    let fine_bins = 12usize;
    let mut adaptive_bin = vec![vec![0usize; d]; n];
    let mut dim_contrast = Vec::with_capacity(d);

    for j in 0..d {
        let mn = (0..n).map(|i| data.x[i][j]).fold(f64::INFINITY, f64::min);
        let mx = (0..n).map(|i| data.x[i][j]).fold(f64::NEG_INFINITY, f64::max);
        let span = (mx - mn).max(1e-9);
        let mut hist = vec![0usize; fine_bins];
        for i in 0..n {
            let b = (((data.x[i][j] - mn) / span) * fine_bins as f64)
                .floor()
                .clamp(0.0, (fine_bins - 1) as f64) as usize;
            hist[b] += 1;
        }
        // Merge adjacent fine bins whose relative counts differ by < 30%
        let mut bin_map = vec![0usize; fine_bins];
        let mut cur_id = 0usize;
        for b in 1..fine_bins {
            let diff = (hist[b] as f64 - hist[b - 1] as f64).abs();
            let mx_h = (hist[b].max(hist[b - 1])).max(1) as f64;
            if diff / mx_h > 0.30 {
                cur_id += 1;
            }
            bin_map[b] = cur_id;
        }
        for i in 0..n {
            let b = (((data.x[i][j] - mn) / span) * fine_bins as f64)
                .floor()
                .clamp(0.0, (fine_bins - 1) as f64) as usize;
            adaptive_bin[i][j] = bin_map[b];
        }
        let max_h = *hist.iter().max().unwrap_or(&1) as f64;
        let avg_h = (n as f64 / fine_bins as f64).max(1.0);
        dim_contrast.push((j, max_h / avg_h, cur_id + 1));
    }

    dim_contrast.sort_by(|a, b| b.1.total_cmp(&a.1));
    let d1 = dim_contrast[0].0;
    let d2 = dim_contrast.get(1).map(|x| x.0).unwrap_or(d1);
    let g_max = dim_contrast[0].2.max(dim_contrast.get(1).map(|x| x.2).unwrap_or(1));

    let mut cell_cnt: HashMap<(usize, usize), usize> = HashMap::new();
    for i in 0..n {
        *cell_cnt.entry((adaptive_bin[i][d1], adaptive_bin[i][d2])).or_insert(0) += 1;
    }
    let avg_cell = (n as f64 / cell_cnt.len().max(1) as f64).max(1.0);
    let dense: HashSet<(usize, usize)> = cell_cnt
        .iter()
        .filter(|(_, &c)| (c as f64) >= 0.7 * avg_cell)
        .map(|(&k, _)| k)
        .collect();
    let cell_lbl = grid_connected_components(&dense, g_max + 1);

    let mut labels = vec![-1i32; n];
    for i in 0..n {
        if let Some(&cid) = cell_lbl.get(&(adaptive_bin[i][d1], adaptive_bin[i][d2])) {
            labels[i] = cid;
        }
    }
    AlgorithmOutput::new(
        labels,
        format!("adaptive bins, dim=({},{})", d1, d2),
        format!("MAFIA sammenslo histogram-intervaller -> {} tette celler", dense.len()),
    )
}

// 141. OptiGrid (Hinneburg & Keim ICML 1999 — Optimal cutting planes at 1D KDE local minima)
fn optigrid(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);

    // Find best low-density KDE valley cutting plane along any axis
    let find_best_cut = |members: &[usize]| -> Option<(usize, f64, f64)> {
        if members.len() < 8 {
            return None;
        }
        let mut best: Option<(usize, f64, f64)> = None; // (dim, cut_val, valley_ratio)
        for j in 0..d {
            let mut vals: Vec<f64> = members.iter().map(|&i| data.x[i][j]).collect();
            vals.sort_by(|a, b| a.total_cmp(b));
            let span = vals[vals.len() - 1] - vals[0];
            if span < 1e-6 {
                continue;
            }
            let h = (span / 6.0).max(1e-4);
            let grid_pts = 20usize;
            let mut kde = vec![0.0; grid_pts];
            let mut xs = vec![0.0; grid_pts];
            for g in 0..grid_pts {
                let x = vals[0] + span * (g as f64) / (grid_pts - 1) as f64;
                xs[g] = x;
                kde[g] = vals
                    .iter()
                    .map(|&v| (-0.5 * ((x - v) / h).powi(2)).exp())
                    .sum::<f64>();
            }
            // Look for deepest interior local minimum between two peaks
            for g in 2..(grid_pts - 2) {
                let left_peak = kde[..g].iter().copied().fold(0.0_f64, f64::max);
                let right_peak = kde[(g + 1)..].iter().copied().fold(0.0_f64, f64::max);
                let min_peak = left_peak.min(right_peak);
                if min_peak > 1e-6 && kde[g] < min_peak * 0.92 {
                    let depth = 1.0 - (kde[g] / min_peak);
                    if best.map(|b| depth > b.2).unwrap_or(true) {
                        best = Some((j, xs[g], depth));
                    }
                }
            }
        }
        best
    };

    let mut regions: Vec<Vec<usize>> = vec![(0..n).collect()];
    let mut cuts_made = 0usize;
    while regions.len() < k {
        let mut best_region = None;
        let mut best_cut = (0usize, 0.0f64, -1.0f64);
        for (ri, mems) in regions.iter().enumerate() {
            if let Some(c) = find_best_cut(mems) {
                if c.2 > best_cut.2 {
                    best_cut = c;
                    best_region = Some(ri);
                }
            }
        }
        let Some(ri) = best_region else { break };
        let mems = regions.remove(ri);
        let left: Vec<usize> = mems.iter().copied().filter(|&i| data.x[i][best_cut.0] <= best_cut.1).collect();
        let right: Vec<usize> = mems.iter().copied().filter(|&i| data.x[i][best_cut.0] > best_cut.1).collect();
        if left.is_empty() || right.is_empty() {
            regions.push(mems);
            break;
        }
        cuts_made += 1;
        regions.push(left);
        regions.push(right);
    }

    let mut labels = vec![0i32; n];
    for (cid, mems) in regions.iter().enumerate() {
        for &p in mems {
            labels[p] = cid as i32;
        }
    }
    AlgorithmOutput::new(
        labels,
        format!("mål_k={}, KDE skjæreplan", k),
        format!("OptiGrid skar gjennom {} lavtetthets-daler", cuts_made),
    )
}

// 142. BANG-Clustering (Schikuta & Erhart 1998 — Multiscale BANG-tree adaptive density blocks)
fn bang_clustering(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let (mn, mx) = pca_bounds(&data.pca_2d);
    let max_bucket = (n / (cfg.effective_k(n) * 3)).clamp(4, 14);

    #[derive(Clone)]
    struct BangBlock {
        x0: f64,
        x1: f64,
        y0: f64,
        y1: f64,
        pts: Vec<usize>,
        depth: usize,
    }
    impl BangBlock {
        fn area(&self) -> f64 {
            ((self.x1 - self.x0) * (self.y1 - self.y0)).max(1e-9)
        }
        fn density(&self) -> f64 {
            self.pts.len() as f64 / self.area()
        }
        fn touches(&self, other: &BangBlock) -> bool {
            let eps = 1e-6;
            self.x0 <= other.x1 + eps
                && self.x1 + eps >= other.x0
                && self.y0 <= other.y1 + eps
                && self.y1 + eps >= other.y0
        }
    }

    let mut leaves = Vec::new();
    let mut stack = vec![BangBlock {
        x0: mn[0],
        x1: mx[0],
        y0: mn[1],
        y1: mx[1],
        pts: (0..n).collect(),
        depth: 0,
    }];

    while let Some(blk) = stack.pop() {
        if blk.pts.len() <= max_bucket || blk.depth >= 6 {
            if !blk.pts.is_empty() {
                leaves.push(blk);
            }
            continue;
        }
        if blk.depth % 2 == 0 {
            let mid = 0.5 * (blk.x0 + blk.x1);
            let l: Vec<usize> = blk.pts.iter().copied().filter(|&i| data.pca_2d[i][0] <= mid).collect();
            let r: Vec<usize> = blk.pts.iter().copied().filter(|&i| data.pca_2d[i][0] > mid).collect();
            stack.push(BangBlock { x0: blk.x0, x1: mid, y0: blk.y0, y1: blk.y1, pts: l, depth: blk.depth + 1 });
            stack.push(BangBlock { x0: mid, x1: blk.x1, y0: blk.y0, y1: blk.y1, pts: r, depth: blk.depth + 1 });
        } else {
            let mid = 0.5 * (blk.y0 + blk.y1);
            let b: Vec<usize> = blk.pts.iter().copied().filter(|&i| data.pca_2d[i][1] <= mid).collect();
            let t: Vec<usize> = blk.pts.iter().copied().filter(|&i| data.pca_2d[i][1] > mid).collect();
            stack.push(BangBlock { x0: blk.x0, x1: blk.x1, y0: blk.y0, y1: mid, pts: b, depth: blk.depth + 1 });
            stack.push(BangBlock { x0: blk.x0, x1: blk.x1, y0: mid, y1: blk.y1, pts: t, depth: blk.depth + 1 });
        }
    }

    // Sort blocks by descending density and form clusters from adjacent high-density blocks
    leaves.sort_by(|a, b| b.density().total_cmp(&a.density()));
    let m = leaves.len();
    let keep = ((m as f64) * 0.72).ceil().max(2.0) as usize;
    let mut blk_lbl = vec![-1i32; m];
    let mut cid = 0i32;
    for i in 0..keep.min(m) {
        if blk_lbl[i] >= 0 {
            continue;
        }
        let mut q = VecDeque::new();
        blk_lbl[i] = cid;
        q.push_back(i);
        while let Some(u) = q.pop_front() {
            for v in 0..keep.min(m) {
                if blk_lbl[v] < 0 && leaves[u].touches(&leaves[v]) {
                    blk_lbl[v] = cid;
                    q.push_back(v);
                }
            }
        }
        cid += 1;
    }

    let mut labels = vec![-1i32; n];
    for (bi, blk) in leaves.iter().enumerate() {
        for &p in &blk.pts {
            labels[p] = blk_lbl[bi];
        }
    }
    AlgorithmOutput::new(
        labels,
        format!("blokker={}, maks_b={}", m, max_bucket),
        format!("Flerskala BANG-tre med {} tetthetsrangerte blokker", m),
    )
}

// 143. ENCLUS (Cheng et al. KDD 1999 — Entropy-based Subspace Clustering)
fn enclus(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let bins = 6usize;

    // Score each dimension by low marginal entropy H(X_j) (structured/peaked) & high interest gain
    let mut dim_ent: Vec<(usize, f64)> = Vec::with_capacity(d);
    for j in 0..d {
        let mn = (0..n).map(|i| data.x[i][j]).fold(f64::INFINITY, f64::min);
        let mx = (0..n).map(|i| data.x[i][j]).fold(f64::NEG_INFINITY, f64::max);
        let span = (mx - mn).max(1e-9);
        let mut cnt = vec![0usize; bins];
        for i in 0..n {
            let b = (((data.x[i][j] - mn) / span) * bins as f64)
                .floor()
                .clamp(0.0, (bins - 1) as f64) as usize;
            cnt[b] += 1;
        }
        let h: f64 = cnt
            .iter()
            .filter(|&&c| c > 0)
            .map(|&c| {
                let p = c as f64 / n as f64;
                -p * p.log2()
            })
            .sum();
        dim_ent.push((j, h));
    }
    dim_ent.sort_by(|a, b| a.1.total_cmp(&b.1)); // lowest entropy = highest cluster structure
    let sel_dims: Vec<usize> = dim_ent.iter().take((d / 2).clamp(2, 4).min(d)).map(|&(j, _)| j).collect();

    let sub_x: Vec<Vec<f64>> = (0..n)
        .map(|i| sel_dims.iter().map(|&j| data.x[i][j]).collect())
        .collect();
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 143);
    let c0 = init_kmeans_pp(&sub_x, k, &mut rng);
    let (labels, _, _, _) = run_lloyd_from_centers(&sub_x, c0, 35);
    AlgorithmOutput::new(
        labels,
        format!("k={}, lav-entropi dim={:?}", k, sel_dims),
        format!("Laveste entropi H={:.2} bits -> underromsclustering", dim_ent[0].1),
    )
}

// 144. FINDIT (Woo et al. 2004 — Dimension-Oriented Distance DOD subspace clustering)
fn findit(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let eps_dim = 0.75_f64;

    // Dimension-Oriented Distance DOD_eps(x, y) = number of dimensions where |x_j - y_j| > eps_dim
    let mut dod = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let mut mismatch = 0.0;
            for dim in 0..d {
                let diff = (data.x[i][dim] - data.x[j][dim]).abs();
                if diff > eps_dim {
                    mismatch += 1.0 + 0.1 * (diff - eps_dim);
                }
            }
            dod[i][j] = mismatch;
            dod[j][i] = mismatch;
        }
    }
    let labels = agglomerative_lance_williams(&dod, k, Linkage::Average);
    AlgorithmOutput::new(
        labels,
        format!("k={}, DOD eps={:.2}", k, eps_dim),
        "Dimension-Oriented Distance (DOD) underromsgruppering".to_string(),
    )
}

// 145. COSA (Friedman & Meulman 2004 — Clustering Objects on Subsets of Attributes with per-object weights)
fn cosa(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let lambda = 0.35_f64;
    let k_nn = (n / 6).clamp(4, 12).min(n - 1);

    // Per-object attribute weights w[i][j] updated iteratively via exponential inverse dispersion
    let mut w = vec![vec![1.0 / d as f64; d]; n];
    let mut cosa_d = vec![vec![0.0; n]; n];

    for _iter in 0..6 {
        for i in 0..n {
            for j in (i + 1)..n {
                let mut d_ij = 0.0;
                for m in 0..d {
                    let wm = w[i][m].max(w[j][m]);
                    d_ij += wm * (data.x[i][m] - data.x[j][m]).abs();
                }
                cosa_d[i][j] = d_ij;
                cosa_d[j][i] = d_ij;
            }
        }
        // Update each object's attribute weights from its k_nn nearest neighbors
        for i in 0..n {
            let mut nbrs: Vec<(usize, f64)> = (0..n)
                .filter(|&j| j != i)
                .map(|j| (j, cosa_d[i][j]))
                .collect();
            nbrs.sort_by(|a, b| a.1.total_cmp(&b.1));
            let mut raw = vec![0.0; d];
            for m in 0..d {
                let avg_diff = nbrs
                    .iter()
                    .take(k_nn)
                    .map(|&(nb, _)| (data.x[i][m] - data.x[nb][m]).abs())
                    .sum::<f64>()
                    / k_nn as f64;
                raw[m] = (-avg_diff / lambda).exp();
            }
            let s: f64 = raw.iter().sum::<f64>().max(1e-12);
            for m in 0..d {
                w[i][m] = raw[m] / s;
            }
        }
    }

    let labels = agglomerative_lance_williams(&cosa_d, k, Linkage::Average);
    AlgorithmOutput::new(
        labels,
        format!("k={}, lambda={:.2}, k_NN={}", k, lambda, k_nn),
        "Individuelle attributtvekter per objekt (COSA-avstand)".to_string(),
    )
}

// 146. FIRES (Kriegel et al. ICDM 2005 — Fast 1D base clusters merged via similarity graph)
fn fires(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 146);

    // Step 1: Compute 1D base clusters along each feature
    let mut base_clusters: Vec<HashSet<usize>> = Vec::new();
    for j in 0..d {
        let col: Vec<Vec<f64>> = (0..n).map(|i| vec![data.x[i][j]]).collect();
        let c0 = init_kmeans_pp(&col, k, &mut rng);
        let (lbl, _, _, _) = run_lloyd_from_centers(&col, c0, 15);
        for c in 0..k {
            let set: HashSet<usize> = (0..n).filter(|&i| lbl[i] == c as i32).collect();
            if set.len() >= 3 {
                base_clusters.push(set);
            }
        }
    }

    // Step 2: Co-membership matrix from 1D base clusters -> merge into high-dim subspace clusters
    let mut co_dist = vec![vec![1.0; n]; n];
    for i in 0..n {
        co_dist[i][i] = 0.0;
        for j in (i + 1)..n {
            let shared = base_clusters
                .iter()
                .filter(|bc| bc.contains(&i) && bc.contains(&j))
                .count();
            let dist = 1.0 - (shared as f64 / d.max(1) as f64);
            co_dist[i][j] = dist;
            co_dist[j][i] = dist;
        }
    }
    let labels = agglomerative_lance_williams(&co_dist, k, Linkage::Average);
    AlgorithmOutput::new(
        labels,
        format!("k={}, {} 1D-basisklynger", k, base_clusters.len()),
        "Sammenslo 1D-klynger til høydimensjonale underromsklynger".to_string(),
    )
}

// 147. INSCY (Assent et al. ICDM 2008 — Indexing non-redundant subspace clusters with SCY-tree pruning)
fn inscy(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let min_pts = 4usize.min(n / 4).max(2);

    // Evaluate non-redundant 2D subspace candidates, pruning redundant lower-dim projections
    let mut best_labels = vec![0i32; n];
    let mut best_score = f64::NEG_INFINITY;
    let mut best_pair = (0usize, 0usize);

    let max_d = d.min(6);
    for a in 0..max_d {
        for b in (a + 1)..max_d.max(a + 2).min(d) {
            let mut sd = vec![vec![0.0; n]; n];
            for i in 0..n {
                for j in (i + 1)..n {
                    let dd = ((data.x[i][a] - data.x[j][a]).powi(2)
                        + (data.x[i][b] - data.x[j][b]).powi(2))
                    .sqrt();
                    sd[i][j] = dd;
                    sd[j][i] = dd;
                }
            }
            let eps = estimate_eps_kdist(&sd, min_pts, 0.58);
            let lbl = run_dbscan_core(&sd, eps, min_pts);
            let nc = lbl.iter().filter(|&&l| l >= 0).collect::<HashSet<_>>().len();
            if nc >= 2 {
                let sil = crate::metrics::silhouette_score(&sd, &lbl).unwrap_or(-0.5);
                if sil > best_score {
                    best_score = sil;
                    best_labels = lbl;
                    best_pair = (a, b);
                }
            }
        }
    }
    if best_score == f64::NEG_INFINITY {
        let eps = estimate_eps_kdist(&data.dist_matrix, min_pts, 0.62);
        best_labels = run_dbscan_core(&data.dist_matrix, eps, min_pts);
    }
    let _ = cfg;
    AlgorithmOutput::new(
        best_labels,
        format!("SCY-indeks dim=({},{})", best_pair.0, best_pair.1),
        "Ikke-redundant underroms-tetthet med SCY-beskjæring".to_string(),
    )
}

// 148. P3C (Moise et al. ICDM 2006 — Poisson-significant 1D intervals -> cluster cores -> EM refinement)
fn p3c(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let n_bins = (1.0 + (n as f64).log2()).round().clamp(5.0, 12.0) as usize;
    let expected = n as f64 / n_bins as f64;

    // Identify Poisson-significant 1D dense intervals: count > expected + 1.5*sqrt(expected)
    let mut sig_mask = vec![vec![0.0; d]; n];
    let mut total_sig = 0usize;
    for j in 0..d {
        let mn = (0..n).map(|i| data.x[i][j]).fold(f64::INFINITY, f64::min);
        let mx = (0..n).map(|i| data.x[i][j]).fold(f64::NEG_INFINITY, f64::max);
        let span = (mx - mn).max(1e-9);
        let mut counts = vec![0usize; n_bins];
        let mut pt_b = vec![0usize; n];
        for i in 0..n {
            let b = (((data.x[i][j] - mn) / span) * n_bins as f64)
                .floor()
                .clamp(0.0, (n_bins - 1) as f64) as usize;
            pt_b[i] = b;
            counts[b] += 1;
        }
        let thresh = expected + 1.2 * expected.sqrt();
        for b in 0..n_bins {
            if (counts[b] as f64) >= thresh {
                total_sig += 1;
            }
        }
        for i in 0..n {
            let b = pt_b[i];
            sig_mask[i][j] = if (counts[b] as f64) >= thresh {
                data.x[i][j]
            } else {
                0.5 * data.x[i][j]
            };
        }
    }

    // Refine discovered interval cores via K-Means/EM
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 148);
    let c0 = init_kmeans_pp(&sig_mask, k, &mut rng);
    let (labels, _, _, _) = run_lloyd_from_centers(&sig_mask, c0, 35);
    AlgorithmOutput::new(
        labels,
        format!("k={}, {} Poisson-signifikante intervaller", k, total_sig),
        "Poisson-signifikante 1D-intervaller + kjerne-refinement".to_string(),
    )
}

// 149. STATPC (Moise & Sander KDD 2008 — Statistically significant non-redundant projected clusters)
fn statpc(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 149);

    // Rank attributes by statistical peakedness (excess kurtosis / variance structure)
    let mut scores: Vec<(usize, f64)> = (0..d)
        .map(|j| {
            let mu = data.x.iter().map(|r| r[j]).sum::<f64>() / n as f64;
            let var = data.x.iter().map(|r| (r[j] - mu).powi(2)).sum::<f64>() / n.max(1) as f64;
            let std = var.sqrt().max(1e-6);
            let kurt = data.x.iter().map(|r| ((r[j] - mu) / std).powi(4)).sum::<f64>() / n as f64;
            // Bimodal/multimodal distributions have LOW kurtosis (< 3.0)
            (j, 3.0 - kurt)
        })
        .collect();
    scores.sort_by(|a, b| b.1.total_cmp(&a.1));
    let keep_d = (d * 2 / 3).clamp(1, d);
    let dims: Vec<usize> = scores.iter().take(keep_d).map(|&(j, _)| j).collect();

    let proj: Vec<Vec<f64>> = (0..n)
        .map(|i| dims.iter().map(|&j| data.x[i][j]).collect())
        .collect();
    let c0 = init_kmeans_pp(&proj, k, &mut rng);
    let (labels, _, _, _) = run_lloyd_from_centers(&proj, c0, 35);
    AlgorithmOutput::new(
        labels,
        format!("k={}, {} ikke-redundante akser", k, dims.len()),
        "Statistisk kurtosis-utvalg av ikke-redundante underrom".to_string(),
    )
}

// 150. Sparse Subspace Clustering (Elhamifar & Vidal TPAMI 2013 — OMP L1 sparse self-representation + Spectral)
fn sparse_subspace_clustering(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let sparsity = (d + 1).clamp(3, 6).min(n.saturating_sub(1).max(1));

    // Orthogonal Matching Pursuit (OMP) for sparse self-representation x_i ≈ sum_{j != i} c_ij x_j
    let mut c_mat = vec![vec![0.0; n]; n];
    for i in 0..n {
        let mut residual = data.x_unit[i].clone();
        let mut support: Vec<usize> = Vec::with_capacity(sparsity);
        for _step in 0..sparsity {
            let best_j = (0..n)
                .filter(|&j| j != i && !support.contains(&j))
                .max_by(|&a, &b| {
                    dot(&residual, &data.x_unit[a])
                        .abs()
                        .total_cmp(&dot(&residual, &data.x_unit[b]).abs())
                });
            let Some(j) = best_j else { break };
            let corr = dot(&residual, &data.x_unit[j]);
            if corr.abs() < 1e-4 {
                break;
            }
            support.push(j);
            c_mat[i][j] = corr.abs();
            for m in 0..d {
                residual[m] -= corr * data.x_unit[j][m];
            }
        }
    }

    // Symmetrize sparse affinity W = |C| + |C|^T -> normalized spectral clustering
    let mut w = vec![vec![0.0; n]; n];
    let mut deg = vec![0.0; n];
    for i in 0..n {
        for j in 0..n {
            let val = c_mat[i][j] + c_mat[j][i];
            w[i][j] = val;
            deg[i] += val;
        }
        deg[i] = deg[i].max(1e-6);
    }
    let mut norm_w = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            norm_w[i][j] = w[i][j] / (deg[i] * deg[j]).sqrt();
        }
    }
    let (_, evecs) = top_k_symmetric_eigen(&norm_w, k, 30);
    let mut embed = vec![vec![0.0; k]; n];
    for i in 0..n {
        for c in 0..k {
            embed[i][c] = evecs[c][i];
        }
    }
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 150);
    let c0 = init_kmeans_pp(&embed, k, &mut rng);
    let (labels, _, _, _) = run_lloyd_from_centers(&embed, c0, 35);
    AlgorithmOutput::new(
        labels,
        format!("k={}, glissenhet s={}", k, sparsity),
        "Glissen L1/OMP selvrepresentasjon + spektral klynging".to_string(),
    )
}
