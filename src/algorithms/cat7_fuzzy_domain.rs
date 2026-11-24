use std::collections::{HashMap, HashSet};

use rand::seq::SliceRandom;
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;

use crate::algorithms::cat1_centroid::{init_kmeans_pp, run_lloyd_from_centers};
use crate::algorithms::cat2_density::run_dbscan_core;
use crate::algorithms::cat3_hierarchical::{agglomerative_lance_williams, Linkage};
use crate::algorithms::{AlgorithmOutput, RunConfig};
use crate::data::PreparedData;
use crate::math::{
    dtw_distance_and_path, estimate_eps_kdist, euclidean, euclidean_sq, inv_and_logdet,
    mahalanobis_sq, median_pairwise_dist, normalize_vec, rbf_affinity_matrix, sbd_distance,
    soft_dtw, top_k_symmetric_eigen,
};

pub fn run(id: usize, data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    match id {
        151 => fuzzy_c_means(data, cfg),
        152 => possibilistic_c_means(data, cfg),
        153 => gustafson_kessel(data, cfg),
        154 => gath_geva(data, cfg),
        155 => kernel_fcm(data, cfg),
        156 => spatial_fcm(data, cfg),
        157 => cheng_church(data, cfg),
        158 => plaid_biclustering(data, cfg),
        159 => bimax(data, cfg),
        160 => isa_biclustering(data, cfg),
        161 => fabia(data, cfg),
        162 => k_shape(data, cfg),
        163 => dba_kmeans(data, cfg),
        164 => soft_dtw_clustering(data, cfg),
        165 => rocket_clustering(data, cfg),
        166 => slic_superpixel(data, cfg),
        167 => watershed_clustering(data, cfg),
        168 => dbscan_3d_pcl(data, cfg),
        169 => traclus(data, cfg),
        170 => consensus_clustering(data, cfg),
        171 => cspa_ensemble(data, cfg),
        172 => hgpa_ensemble(data, cfg),
        173 => mcla_ensemble(data, cfg),
        174 => multi_view_kmeans(data, cfg),
        175 => coreg_multi_view_spectral(data, cfg),
        _ => unreachable!(),
    }
}

fn run_fcm_core(
    x: &[Vec<f64>],
    k: usize,
    m: f64,
    max_iter: usize,
    seed: u64,
) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let n = x.len();
    let d = x[0].len();
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    let mut centers = init_kmeans_pp(x, k, &mut rng);
    let mut u = vec![vec![1.0 / k as f64; k]; n];
    let exp_p = 2.0 / (m - 1.0);

    for _ in 0..max_iter {
        for i in 0..n {
            let dists: Vec<f64> = (0..k)
                .map(|c| euclidean(&x[i], &centers[c]).max(1e-8))
                .collect();
            for c in 0..k {
                let denom: f64 = (0..k).map(|j| (dists[c] / dists[j]).powf(exp_p)).sum();
                u[i][c] = 1.0 / denom.max(1e-12);
            }
        }
        for c in 0..k {
            let mut num = vec![0.0; d];
            let mut den = 0.0;
            for i in 0..n {
                let w = u[i][c].powf(m);
                den += w;
                for j in 0..d {
                    num[j] += w * x[i][j];
                }
            }
            if den > 1e-12 {
                for j in 0..d {
                    centers[c][j] = num[j] / den;
                }
            }
        }
    }
    (u, centers)
}

fn defuzzify(u: &[Vec<f64>]) -> (Vec<i32>, f64) {
    let mut max_sum = 0.0;
    let labels: Vec<i32> = u
        .iter()
        .map(|row| {
            let (best_c, &best_u) = row
                .iter()
                .enumerate()
                .max_by(|a, b| a.1.total_cmp(b.1))
                .unwrap_or((0, &0.0));
            max_sum += best_u;
            best_c as i32
        })
        .collect();
    (labels, max_sum / u.len().max(1) as f64)
}

// 151. Fuzzy C-Means (Bezdek 1981)
fn fuzzy_c_means(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let k = cfg.effective_k(data.n_samples);
    let m = 2.0;
    let (u, _) = run_fcm_core(&data.x, k, m, cfg.max_iters(40), cfg.seed + 151);
    let (labels, avg_u) = defuzzify(&u);
    AlgorithmOutput::new(
        labels,
        format!("k={}, fuzzifier m={:.1}", k, m),
        format!("Snitt maksimalt medlemskap={:.3}", avg_u),
    )
}

// 152. Possibilistic C-Means (Krishnapuram & Keller 1993 — relaxes sum-to-1 constraint for noise robustness)
fn possibilistic_c_means(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let m = 2.0_f64;
    let (mut u, mut centers) = run_fcm_core(&data.x, k, m, 15, cfg.seed + 152);

    // Estimate per-cluster bandwidth eta_c
    let mut eta = vec![1.0; k];
    for c in 0..k {
        let num: f64 = (0..n)
            .map(|i| u[i][c].powf(m) * euclidean_sq(&data.x[i], &centers[c]))
            .sum();
        let den: f64 = (0..n).map(|i| u[i][c].powf(m)).sum::<f64>().max(1e-6);
        eta[c] = (num / den).max(0.05);
    }

    for _ in 0..cfg.max_iters(30) {
        for i in 0..n {
            for c in 0..k {
                let d2 = euclidean_sq(&data.x[i], &centers[c]);
                u[i][c] = 1.0 / (1.0 + (d2 / eta[c]).powf(1.0 / (m - 1.0)));
            }
        }
        for c in 0..k {
            let mut num = vec![0.0; d];
            let mut den = 0.0;
            for i in 0..n {
                let w = u[i][c].powf(m);
                den += w;
                for j in 0..d {
                    num[j] += w * data.x[i][j];
                }
            }
            if den > 1e-12 {
                for j in 0..d {
                    centers[c][j] = num[j] / den;
                }
            }
        }
    }

    // Points with max typicality < 0.15 are flagged as outliers (-1)
    let mut labels = vec![-1i32; n];
    let mut noise = 0usize;
    for i in 0..n {
        let (best_c, &best_t) = u[i]
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.total_cmp(b.1))
            .unwrap_or((0, &0.0));
        if best_t < 0.15 {
            noise += 1;
        } else {
            labels[i] = best_c as i32;
        }
    }
    AlgorithmOutput::new(
        labels,
        format!("k={}, m={:.1}, typikalitet", k, m),
        format!("Possibilistisk medlemskap isolerte {} uteliggere", noise),
    )
}

// 153. Gustafson-Kessel (1979 — Adaptive fuzzy covariance matrix with volume constraint det(A_c) = 1)
fn gustafson_kessel(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let m = 2.0_f64;
    let reg = 0.05;
    let (mut u, mut centers) = run_fcm_core(&data.x, k, m, 10, cfg.seed + 153);

    for _ in 0..cfg.max_iters(25) {
        let mut norm_invs = Vec::with_capacity(k);
        for c in 0..k {
            let mut fcov = vec![vec![0.0; d]; d];
            let den: f64 = (0..n).map(|i| u[i][c].powf(m)).sum::<f64>().max(1e-6);
            for i in 0..n {
                let w = u[i][c].powf(m);
                for a in 0..d {
                    for b in 0..d {
                        fcov[a][b] += w * (data.x[i][a] - centers[c][a]) * (data.x[i][b] - centers[c][b]);
                    }
                }
            }
            for a in 0..d {
                for b in 0..d {
                    fcov[a][b] = fcov[a][b] / den + if a == b { reg } else { 0.0 };
                }
            }
            let (inv, logdet) = inv_and_logdet(&fcov, reg);
            // Scale by det(F)^(1/d) for Gustafson-Kessel unit-volume norm-inducing matrix
            let det_root = (logdet / d as f64).exp().clamp(0.05, 50.0);
            let a_mat: Vec<Vec<f64>> = inv
                .into_iter()
                .map(|row| row.into_iter().map(|v| v * det_root).collect())
                .collect();
            norm_invs.push(a_mat);
        }

        // Update fuzzy memberships with adaptive Mahalanobis-GK distance
        for i in 0..n {
            let ds: Vec<f64> = (0..k)
                .map(|c| mahalanobis_sq(&data.x[i], &centers[c], &norm_invs[c]).max(1e-8))
                .collect();
            for c in 0..k {
                let denom: f64 = (0..k).map(|j| (ds[c] / ds[j]).powf(1.0 / (m - 1.0))).sum();
                u[i][c] = 1.0 / denom.max(1e-12);
            }
        }
        // Update centers
        for c in 0..k {
            let mut num = vec![0.0; d];
            let mut den = 0.0;
            for i in 0..n {
                let w = u[i][c].powf(m);
                den += w;
                for j in 0..d {
                    num[j] += w * data.x[i][j];
                }
            }
            if den > 1e-12 {
                for j in 0..d {
                    centers[c][j] = num[j] / den;
                }
            }
        }
    }

    let (labels, avg_u) = defuzzify(&u);
    AlgorithmOutput::new(
        labels,
        format!("k={}, adaptiv kovarians", k),
        format!("Ellipsoide GK-klynger (snitt maks-u={:.3})", avg_u),
    )
}

// 154. Gath-Geva (1989 — Fuzzy Maximum Likelihood Estimation exponential distance)
fn gath_geva(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let m = 2.0_f64;
    let reg = 0.08;
    let (mut u, mut centers) = run_fcm_core(&data.x, k, m, 10, cfg.seed + 154);

    for _ in 0..cfg.max_iters(20) {
        let mut inv_lds = Vec::with_capacity(k);
        let mut priors = vec![1.0 / k as f64; k];
        for c in 0..k {
            let den: f64 = (0..n).map(|i| u[i][c].powf(m)).sum::<f64>().max(1e-6);
            priors[c] = (den / n as f64).max(1e-3);
            let mut fcov = vec![vec![0.0; d]; d];
            for i in 0..n {
                let w = u[i][c].powf(m);
                for a in 0..d {
                    for b in 0..d {
                        fcov[a][b] += w * (data.x[i][a] - centers[c][a]) * (data.x[i][b] - centers[c][b]);
                    }
                }
            }
            for a in 0..d {
                for b in 0..d {
                    fcov[a][b] = fcov[a][b] / den + if a == b { reg } else { 0.0 };
                }
            }
            inv_lds.push(inv_and_logdet(&fcov, reg));
        }

        // Exponential FMLE distance: D_FMLE(x, c) = sqrt(det F_c) / pi_c * exp(0.5 * mahalanobis^2)
        for i in 0..n {
            let mut log_d = vec![0.0; k];
            let mut mn = f64::INFINITY;
            for c in 0..k {
                let m2 = mahalanobis_sq(&data.x[i], &centers[c], &inv_lds[c].0).min(40.0);
                let ld = 0.5 * inv_lds[c].1 - priors[c].ln() + 0.5 * m2;
                log_d[c] = ld;
                mn = mn.min(ld);
            }
            let ds: Vec<f64> = log_d.iter().map(|&v| (v - mn).exp().max(1e-8)).collect();
            for c in 0..k {
                let denom: f64 = (0..k).map(|j| ds[c] / ds[j]).sum();
                u[i][c] = 1.0 / denom.max(1e-12);
            }
        }
        for c in 0..k {
            let mut num = vec![0.0; d];
            let mut den = 0.0;
            for i in 0..n {
                let w = u[i][c].powf(m);
                den += w;
                for j in 0..d {
                    num[j] += w * data.x[i][j];
                }
            }
            if den > 1e-12 {
                for j in 0..d {
                    centers[c][j] = num[j] / den;
                }
            }
        }
    }

    let (labels, avg_u) = defuzzify(&u);
    AlgorithmOutput::new(
        labels,
        format!("k={}, Gaussisk FMLE-avstand", k),
        format!("Gath-Geva eksponentiell fuzzy-MLE (snitt u={:.3})", avg_u),
    )
}

// 155. Kernel Fuzzy C-Means (Zhang & Chen 2004 — RBF kernel distance 2*(1 - K(x_i, v_c)))
fn kernel_fcm(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let m = 2.0_f64;
    let sigma = median_pairwise_dist(&data.dist_matrix).max(1e-3);
    let s2 = 2.0 * sigma * sigma;
    let (mut u, mut centers) = run_fcm_core(&data.x, k, m, 8, cfg.seed + 155);

    for _ in 0..cfg.max_iters(30) {
        let mut kval = vec![vec![0.0; k]; n];
        for i in 0..n {
            for c in 0..k {
                kval[i][c] = (-euclidean_sq(&data.x[i], &centers[c]) / s2).exp();
            }
            let dists: Vec<f64> = (0..k).map(|c| (1.0 - kval[i][c]).max(1e-7)).collect();
            for c in 0..k {
                let denom: f64 = (0..k).map(|j| (dists[c] / dists[j]).powf(1.0 / (m - 1.0))).sum();
                u[i][c] = 1.0 / denom.max(1e-12);
            }
        }
        // Kernel-weighted prototype update v_c = sum_i u_ic^m K(x_i, v_c) x_i / sum_i u_ic^m K(x_i, v_c)
        for c in 0..k {
            let mut num = vec![0.0; d];
            let mut den = 0.0;
            for i in 0..n {
                let w = u[i][c].powf(m) * kval[i][c];
                den += w;
                for j in 0..d {
                    num[j] += w * data.x[i][j];
                }
            }
            if den > 1e-12 {
                for j in 0..d {
                    centers[c][j] = num[j] / den;
                }
            }
        }
    }
    let (labels, avg_u) = defuzzify(&u);
    AlgorithmOutput::new(
        labels,
        format!("k={}, RBF sigma={:.2}", k, sigma),
        format!("Kjernebasert KFCM (snitt maks-u={:.3})", avg_u),
    )
}

// 156. Spatial Fuzzy C-Means (Chuang et al. 2006 — incorporates spatial k-NN neighborhood membership h_ic)
fn spatial_fcm(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let m = 2.0_f64;
    let (mut u, mut centers) = run_fcm_core(&data.x, k, m, 10, cfg.seed + 156);

    for _ in 0..cfg.max_iters(25) {
        // Compute standard FCM memberships
        for i in 0..n {
            let dists: Vec<f64> = (0..k)
                .map(|c| euclidean(&data.x[i], &centers[c]).max(1e-8))
                .collect();
            for c in 0..k {
                let denom: f64 = (0..k).map(|j| (dists[c] / dists[j]).powf(2.0 / (m - 1.0))).sum();
                u[i][c] = 1.0 / denom.max(1e-12);
            }
        }
        // Spatial function h_ic = weighted average of neighbors' membership in k-NN graph
        let mut u_sp = vec![vec![0.0; k]; n];
        for i in 0..n {
            let mut h = vec![0.0; k];
            let mut w_tot = 0.0;
            for &(nb, wt) in &data.knn_graph[i] {
                w_tot += wt;
                for c in 0..k {
                    h[c] += wt * u[nb][c];
                }
            }
            if w_tot > 0.0 {
                for c in 0..k {
                    h[c] /= w_tot;
                }
                let sum_uh: f64 = (0..k).map(|c| u[i][c] * h[c]).sum::<f64>().max(1e-12);
                for c in 0..k {
                    u_sp[i][c] = (u[i][c] * h[c]) / sum_uh;
                }
            } else {
                u_sp[i] = u[i].clone();
            }
        }
        u = u_sp;
        for c in 0..k {
            let mut num = vec![0.0; d];
            let mut den = 0.0;
            for i in 0..n {
                let w = u[i][c].powf(m);
                den += w;
                for j in 0..d {
                    num[j] += w * data.x[i][j];
                }
            }
            if den > 1e-12 {
                for j in 0..d {
                    centers[c][j] = num[j] / den;
                }
            }
        }
    }
    let (labels, avg_u) = defuzzify(&u);
    AlgorithmOutput::new(
        labels,
        format!("k={}, romlig p=1, q=1", k),
        format!("sFCM med naboglattet medlemskap (snitt u={:.3})", avg_u),
    )
}

// 157. Cheng-Church Biclustering (ISMB 2000 — iterative node deletion minimizing Mean Squared Residue MSR)
fn cheng_church(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 157);
    let mut mat = data.x.clone();

    let compute_msr = |m: &[Vec<f64>], rows: &[usize], cols: &[usize]| -> (f64, Vec<f64>, Vec<f64>) {
        let nr = rows.len().max(1) as f64;
        let nc = cols.len().max(1) as f64;
        let mut row_mean = vec![0.0; rows.len()];
        let mut col_mean = vec![0.0; cols.len()];
        let mut grand = 0.0;
        for (ri, &r) in rows.iter().enumerate() {
            for (ci, &c) in cols.iter().enumerate() {
                let v = m[r][c];
                row_mean[ri] += v;
                col_mean[ci] += v;
                grand += v;
            }
        }
        grand /= nr * nc;
        for v in &mut row_mean {
            *v /= nc;
        }
        for v in &mut col_mean {
            *v /= nr;
        }
        let mut row_res = vec![0.0; rows.len()];
        let mut col_res = vec![0.0; cols.len()];
        let mut msr = 0.0;
        for (ri, &r) in rows.iter().enumerate() {
            for (ci, &c) in cols.iter().enumerate() {
                let res = (m[r][c] - row_mean[ri] - col_mean[ci] + grand).powi(2);
                row_res[ri] += res / nc;
                col_res[ci] += res / nr;
                msr += res;
            }
        }
        (msr / (nr * nc), row_res, col_res)
    };

    let mut labels = vec![-1i32; n];
    let mut available_rows: Vec<usize> = (0..n).collect();
    let delta = 0.28_f64;

    for b_id in 0..k {
        if available_rows.len() < 3 {
            break;
        }
        let mut rows = available_rows.clone();
        let mut cols: Vec<usize> = (0..d).collect();
        let target_min_r = (available_rows.len() / (k - b_id).max(1)).max(3);

        for _ in 0..60 {
            let (msr, r_res, c_res) = compute_msr(&mat, &rows, &cols);
            if msr <= delta || rows.len() <= target_min_r {
                break;
            }
            let (max_ri, &mr) = r_res
                .iter()
                .enumerate()
                .max_by(|a, b| a.1.total_cmp(b.1))
                .unwrap();
            let (max_ci, &mc) = c_res
                .iter()
                .enumerate()
                .max_by(|a, b| a.1.total_cmp(b.1))
                .unwrap();
            if mc > mr && cols.len() > 2.min(d) {
                cols.remove(max_ci);
            } else {
                rows.remove(max_ri);
            }
        }
        let chosen: HashSet<usize> = rows.iter().copied().collect();
        for &r in &rows {
            labels[r] = b_id as i32;
            // Mask discovered bicluster entries with random noise (classic Cheng-Church masking)
            for &c in &cols {
                mat[r][c] = rng.gen_range(-1.0..1.0);
            }
        }
        available_rows.retain(|r| !chosen.contains(r));
    }
    // Assign any remaining rows to last bicluster
    for &r in &available_rows {
        labels[r] = (k - 1) as i32;
    }

    AlgorithmOutput::new(
        labels,
        format!("k={}, MSR delta={:.2}", k, delta),
        "Iterativ Cheng-Church biklynge-ekstraksjon med MSR-reduksjon".to_string(),
    )
}

// 158. PLAID (Lazzeroni & Owen 2002 — Additive layer biclustering X ≈ mu_0 + sum_l theta_l r_l c_l^T)
fn plaid_biclustering(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);

    // Subtract global mean layer
    let grand = data.x.iter().flat_map(|r| r.iter()).sum::<f64>() / (n * d).max(1) as f64;
    let mut resid: Vec<Vec<f64>> = data
        .x
        .iter()
        .map(|row| row.iter().map(|&v| v - grand).collect())
        .collect();

    let mut row_scores = vec![vec![0.0; k]; n];
    for layer in 0..k {
        // Rank-1 alternating least squares layer fit on residual matrix
        let mut r_vec = vec![1.0 / (n as f64).sqrt(); n];
        let mut c_vec = vec![1.0 / (d as f64).sqrt(); d];
        for i in 0..n {
            r_vec[i] = ((i + layer + 1) as f64).sin();
        }
        r_vec = normalize_vec(&r_vec);
        for _ in 0..12 {
            // Update c_vec = R^T r
            for j in 0..d {
                c_vec[j] = (0..n).map(|i| resid[i][j] * r_vec[i]).sum();
            }
            c_vec = normalize_vec(&c_vec);
            // Update r_vec = R c
            for i in 0..n {
                r_vec[i] = (0..d).map(|j| resid[i][j] * c_vec[j]).sum();
            }
            let rn: f64 = r_vec.iter().map(|v| v * v).sum::<f64>().sqrt().max(1e-8);
            for i in 0..n {
                r_vec[i] /= rn;
            }
        }
        // Estimate layer amplitude and peel off additive layer
        let mut mu_l = 0.0;
        for i in 0..n {
            for j in 0..d {
                mu_l += resid[i][j] * r_vec[i] * c_vec[j];
            }
        }
        for i in 0..n {
            row_scores[i][layer] = (mu_l * r_vec[i]).abs();
            for j in 0..d {
                resid[i][j] -= mu_l * r_vec[i] * c_vec[j];
            }
        }
    }

    let labels: Vec<i32> = (0..n)
        .map(|i| {
            (0..k)
                .max_by(|&a, &b| row_scores[i][a].total_cmp(&row_scores[i][b]))
                .unwrap_or(0) as i32
        })
        .collect();
    AlgorithmOutput::new(
        labels,
        format!("{} additive PLAID-lag", k),
        "Sekvensiell lagvis dekomponering av residualmatrise".to_string(),
    )
}

// 159. Bimax (Prelic et al. Bioinformatics 2006 — Divide-and-conquer maximal all-1s binary biclusters)
fn bimax(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.x_binary[0].len();
    let k = cfg.effective_k(n);

    // Divide-and-conquer binary submatrix search for dense all-1 blocks
    let mut biclusters: Vec<(Vec<usize>, Vec<usize>)> = Vec::new();
    let mut stack: Vec<(Vec<usize>, Vec<usize>)> = vec![((0..n).collect(), (0..d).collect())];

    while let Some((rows, cols)) = stack.pop() {
        if rows.len() < 3 || cols.is_empty() || biclusters.len() >= 24 {
            continue;
        }
        let all_ones = rows.iter().all(|&r| cols.iter().all(|&c| data.x_binary[r][c] == 1));
        if all_ones {
            biclusters.push((rows, cols));
            continue;
        }
        // Pick row with both 0s and 1s to partition columns into C_u (1s) and C_v (0s)
        let pivot = rows.iter().copied().find(|&r| {
            let ones = cols.iter().filter(|&&c| data.x_binary[r][c] == 1).count();
            ones > 0 && ones < cols.len()
        });
        let Some(pr) = pivot else {
            // Every row is either all-0 or all-1
            let ones_rows: Vec<usize> = rows
                .into_iter()
                .filter(|&r| data.x_binary[r][cols[0]] == 1)
                .collect();
            if ones_rows.len() >= 3 {
                biclusters.push((ones_rows, cols));
            }
            continue;
        };
        let c_u: Vec<usize> = cols.iter().copied().filter(|&c| data.x_binary[pr][c] == 1).collect();
        let c_v: Vec<usize> = cols.iter().copied().filter(|&c| data.x_binary[pr][c] == 0).collect();

        let r_u: Vec<usize> = rows
            .iter()
            .copied()
            .filter(|&r| c_v.iter().all(|&c| data.x_binary[r][c] == 0))
            .collect();
        let r_w: Vec<usize> = rows
            .iter()
            .copied()
            .filter(|&r| {
                c_u.iter().any(|&c| data.x_binary[r][c] == 1)
                    && c_v.iter().any(|&c| data.x_binary[r][c] == 1)
            })
            .collect();
        let r_v: Vec<usize> = rows
            .iter()
            .copied()
            .filter(|&r| c_u.iter().all(|&c| data.x_binary[r][c] == 0))
            .collect();

        let mut ru_rw = r_u;
        ru_rw.extend_from_slice(&r_w);
        if !ru_rw.is_empty() && !c_u.is_empty() {
            stack.push((ru_rw, c_u));
        }
        let mut rw_rv = r_w;
        rw_rv.extend_from_slice(&r_v);
        if !rw_rv.is_empty() && !c_v.is_empty() {
            stack.push((rw_rv, cols));
        }
    }

    biclusters.sort_by_key(|(r, c)| std::cmp::Reverse(r.len() * c.len()));
    let mut labels = vec![-1i32; n];
    for (cid, (rows, _)) in biclusters.iter().take(k).enumerate() {
        for &r in rows {
            if labels[r] < 0 {
                labels[r] = cid as i32;
            }
        }
    }
    AlgorithmOutput::new(
        labels,
        format!("binær splitt, {} maksimale blokker", biclusters.len()),
        "Bimax fant maksimale binære all-1 delmatriser".to_string(),
    )
}

// 160. ISA (Bergmann et al. 2003 — Iterative Signature Algorithm alternating row/col thresholding)
fn isa_biclustering(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 160);

    // Column-standardized & row-standardized matrices
    let mut ec = data.x.clone();
    let mut er = data.x.clone();
    for i in 0..n {
        let mu = er[i].iter().sum::<f64>() / d as f64;
        let std = (er[i].iter().map(|v| (v - mu).powi(2)).sum::<f64>() / d.max(1) as f64)
            .sqrt()
            .max(1e-6);
        for j in 0..d {
            er[i][j] = (er[i][j] - mu) / std;
        }
    }
    for j in 0..d {
        let mu = ec.iter().map(|r| r[j]).sum::<f64>() / n as f64;
        let std = (ec.iter().map(|r| (r[j] - mu).powi(2)).sum::<f64>() / n.max(1) as f64)
            .sqrt()
            .max(1e-6);
        for i in 0..n {
            ec[i][j] = (ec[i][j] - mu) / std;
        }
    }

    let mut sig_scores = vec![vec![0.0; k]; n];
    let t_r = 0.6;
    let t_c = 0.5;

    for s_id in 0..k {
        let mut r_w: Vec<f64> = (0..n).map(|_| rng.gen_range(0.0..1.0)).collect();
        for _ in 0..18 {
            // c = threshold(E_c^T r)
            let mut c_w = vec![0.0; d];
            for j in 0..d {
                c_w[j] = (0..n).map(|i| ec[i][j] * r_w[i]).sum();
            }
            let c_mu = c_w.iter().sum::<f64>() / d as f64;
            let c_std = (c_w.iter().map(|v| (v - c_mu).powi(2)).sum::<f64>() / d.max(1) as f64)
                .sqrt()
                .max(1e-6);
            for j in 0..d {
                if (c_w[j] - c_mu).abs() < t_c * c_std {
                    c_w[j] = 0.0;
                }
            }
            c_w = normalize_vec(&c_w);
            // r = threshold(E_r c)
            for i in 0..n {
                r_w[i] = (0..d).map(|j| er[i][j] * c_w[j]).sum();
            }
            let r_mu = r_w.iter().sum::<f64>() / n as f64;
            let r_std = (r_w.iter().map(|v| (v - r_mu).powi(2)).sum::<f64>() / n.max(1) as f64)
                .sqrt()
                .max(1e-6);
            for i in 0..n {
                if (r_w[i] - r_mu).abs() < t_r * r_std {
                    r_w[i] = 0.0;
                }
            }
            let rn: f64 = r_w.iter().map(|v| v * v).sum::<f64>().sqrt();
            if rn > 1e-8 {
                for v in &mut r_w {
                    *v /= rn;
                }
            }
        }
        for i in 0..n {
            sig_scores[i][s_id] = r_w[i].abs();
        }
    }

    let labels: Vec<i32> = (0..n)
        .map(|i| {
            let (best_c, &val) = sig_scores[i]
                .iter()
                .enumerate()
                .max_by(|a, b| a.1.total_cmp(b.1))
                .unwrap_or((0, &0.0));
            if val <= 1e-9 { -1 } else { best_c as i32 }
        })
        .collect();
    AlgorithmOutput::new(
        labels,
        format!("k={}, terskler t_r={:.1}, t_c={:.1}", k, t_r, t_c),
        "Iterativ rad/kolonne-terskling fant ko-regulerte moduler".to_string(),
    )
}

// 161. FABIA (Hochreiter et al. Bioinformatics 2010 — Sparse Laplacian Factor Analysis X = L Z + eps)
fn fabia(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let alpha_l1 = 0.12_f64;
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 161);

    // Factors Z (n x k) and loadings L (k x d) with L1 soft-thresholding (Laplacian prior)
    let mut z = vec![vec![0.0; k]; n];
    let mut l = vec![vec![0.0; d]; k];
    for i in 0..n {
        for c in 0..k {
            z[i][c] = rng.gen_range(-0.5..0.5);
        }
    }
    for c in 0..k {
        for j in 0..d {
            l[c][j] = rng.gen_range(-0.5..0.5);
        }
    }
    let soft_thresh = |v: f64, lam: f64| -> f64 {
        if v > lam {
            v - lam
        } else if v < -lam {
            v + lam
        } else {
            0.0
        }
    };

    for _ in 0..cfg.max_iters(35) {
        // Update Z given L via gradient + L1 proximal step
        for i in 0..n {
            for c in 0..k {
                let mut pred_err = 0.0;
                let mut l_norm = 1e-6;
                for j in 0..d {
                    let pred: f64 = (0..k).map(|m| z[i][m] * l[m][j]).sum();
                    pred_err += (data.x[i][j] - pred + z[i][c] * l[c][j]) * l[c][j];
                    l_norm += l[c][j] * l[c][j];
                }
                z[i][c] = soft_thresh(pred_err / l_norm, alpha_l1);
            }
        }
        // Update L given Z with L1 sparsity
        for c in 0..k {
            let z_norm: f64 = (0..n).map(|i| z[i][c] * z[i][c]).sum::<f64>().max(1e-6);
            for j in 0..d {
                let mut num = 0.0;
                for i in 0..n {
                    let pred: f64 = (0..k).map(|m| z[i][m] * l[m][j]).sum();
                    num += (data.x[i][j] - pred + z[i][c] * l[c][j]) * z[i][c];
                }
                l[c][j] = soft_thresh(num / z_norm, alpha_l1 * 0.5);
            }
        }
    }

    let labels: Vec<i32> = (0..n)
        .map(|i| {
            (0..k)
                .max_by(|&a, &b| z[i][a].abs().total_cmp(&z[i][b].abs()))
                .unwrap_or(0) as i32
        })
        .collect();
    AlgorithmOutput::new(
        labels,
        format!("k={}, Laplace L1={:.2}", k, alpha_l1),
        "Glissen faktoranalyse-biklynging X ≈ Z*L".to_string(),
    )
}

// 162. k-Shape (Paparrizos & Gravano SIGMOD 2015 — Shape-Based Distance via normalized cross-correlation)
fn k_shape(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 162);

    // Z-normalize each observation treated as a 1D ordered profile/series
    let series: Vec<Vec<f64>> = data
        .x
        .iter()
        .map(|row| {
            let mu = row.iter().sum::<f64>() / d.max(1) as f64;
            let std = (row.iter().map(|v| (v - mu).powi(2)).sum::<f64>() / d.max(1) as f64)
                .sqrt()
                .max(1e-6);
            row.iter().map(|&v| (v - mu) / std).collect()
        })
        .collect();

    let mut labels: Vec<usize> = (0..n).map(|i| i % k).collect();
    labels.shuffle(&mut rng);
    let mut centroids: Vec<Vec<f64>> = (0..k).map(|c| series[c % n].clone()).collect();

    for _ in 0..cfg.max_iters(20) {
        // Shape extraction: shift-align members to centroid & compute top eigenvector of S = X^T X
        for c in 0..k {
            let mems: Vec<usize> = (0..n).filter(|&i| labels[i] == c).collect();
            if mems.is_empty() {
                continue;
            }
            let mut aligned = vec![vec![0.0; d]; mems.len()];
            for (idx, &m) in mems.iter().enumerate() {
                let (_, shift) = sbd_distance(&series[m], &centroids[c]);
                for t in 0..d {
                    let src = t as isize - shift;
                    if src >= 0 && (src as usize) < d {
                        aligned[idx][t] = series[m][src as usize];
                    }
                }
            }
            let mut s_mat = vec![vec![0.0; d]; d];
            for row in &aligned {
                for a in 0..d {
                    for b in a..d {
                        s_mat[a][b] += row[a] * row[b];
                    }
                }
            }
            for a in 0..d {
                for b in a..d {
                    s_mat[b][a] = s_mat[a][b];
                }
            }
            let (_, evecs) = top_k_symmetric_eigen(&s_mat, 1, 15);
            if let Some(v) = evecs.first() {
                centroids[c] = v.clone();
            }
        }

        // Assignment via Shape-Based Distance
        let mut changed = false;
        for i in 0..n {
            let best_c = (0..k)
                .min_by(|&a, &b| {
                    sbd_distance(&series[i], &centroids[a])
                        .0
                        .total_cmp(&sbd_distance(&series[i], &centroids[b]).0)
                })
                .unwrap_or(0);
            if labels[i] != best_c {
                labels[i] = best_c;
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }

    AlgorithmOutput::new(
        labels.into_iter().map(|c| c as i32).collect(),
        format!("k={}, SBD krysskorrelasjon", k),
        "Fase-invariant tidsserie-clustering med normalisert krysskorrelasjon".to_string(),
    )
}

// 163. DBA (Petitjean et al. 2011 — DTW Barycenter Averaging K-Means)
fn dba_kmeans(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 163);
    let mut centers = init_kmeans_pp(&data.x, k, &mut rng);
    let mut labels = vec![0i32; n];

    for _ in 0..cfg.max_iters(15) {
        let mut changed = false;
        for i in 0..n {
            let best_c = (0..k)
                .min_by(|&a, &b| {
                    dtw_distance_and_path(&data.x[i], &centers[a])
                        .0
                        .total_cmp(&dtw_distance_and_path(&data.x[i], &centers[b]).0)
                })
                .unwrap_or(0);
            if labels[i] != best_c as i32 {
                labels[i] = best_c as i32;
                changed = true;
            }
        }
        if !changed {
            break;
        }
        // DTW Barycenter Averaging update for each cluster centroid
        for c in 0..k {
            let mems: Vec<usize> = (0..n).filter(|&i| labels[i] == c as i32).collect();
            if mems.is_empty() {
                continue;
            }
            let mut sums = vec![0.0; d];
            let mut cnts = vec![0usize; d];
            for &m in &mems {
                let (_, path) = dtw_distance_and_path(&centers[c], &data.x[m]);
                for (ci, mi) in path {
                    sums[ci] += data.x[m][mi];
                    cnts[ci] += 1;
                }
            }
            for t in 0..d {
                if cnts[t] > 0 {
                    centers[c][t] = sums[t] / cnts[t] as f64;
                }
            }
        }
    }
    AlgorithmOutput::new(
        labels,
        format!("k={}, Dynamic Time Warping", k),
        "DTW Barycenter Averaging (DBA) over sekvensdimensjoner".to_string(),
    )
}

// 164. Soft-DTW Clustering (Cuturi & Blondel ICML 2017)
fn soft_dtw_clustering(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let gamma = 0.5_f64;
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 164);
    let mut centers = init_kmeans_pp(&data.x, k, &mut rng);
    let mut labels = vec![0i32; n];

    for _ in 0..cfg.max_iters(14) {
        let mut changed = false;
        for i in 0..n {
            let best_c = (0..k)
                .min_by(|&a, &b| {
                    soft_dtw(&data.x[i], &centers[a], gamma)
                        .total_cmp(&soft_dtw(&data.x[i], &centers[b], gamma))
                })
                .unwrap_or(0);
            if labels[i] != best_c as i32 {
                labels[i] = best_c as i32;
                changed = true;
            }
        }
        if !changed {
            break;
        }
        // Finite-difference gradient step on Soft-DTW barycenter objective
        let eps = 1e-3;
        for c in 0..k {
            let mems: Vec<usize> = (0..n).filter(|&i| labels[i] == c as i32).collect();
            if mems.is_empty() {
                continue;
            }
            // Start from Euclidean mean then refine 2 gradient steps on Soft-DTW
            let mut bary = vec![0.0; d];
            for &m in &mems {
                for j in 0..d {
                    bary[j] += data.x[m][j] / mems.len() as f64;
                }
            }
            for _step in 0..2 {
                let mut grad = vec![0.0; d];
                let base_loss: f64 = mems.iter().map(|&m| soft_dtw(&bary, &data.x[m], gamma)).sum();
                for j in 0..d {
                    bary[j] += eps;
                    let plus_loss: f64 = mems.iter().map(|&m| soft_dtw(&bary, &data.x[m], gamma)).sum();
                    bary[j] -= eps;
                    grad[j] = (plus_loss - base_loss) / (eps * mems.len() as f64);
                }
                for j in 0..d {
                    bary[j] -= 0.1 * grad[j];
                }
            }
            centers[c] = bary;
        }
    }
    AlgorithmOutput::new(
        labels,
        format!("k={}, Soft-DTW gamma={:.1}", k, gamma),
        "Deriverbar glattet Soft-DTW tidsserie-barysenter".to_string(),
    )
}

// 165. ROCKET-Clustering (Dempster et al. 2020 — Random Convolutional Kernel Transform -> PPV & Max -> K-Means)
fn rocket_clustering(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let n_kernels = 32usize;
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 165);

    // Generate random 1D kernels (length 2..=min(5, d), weights zero-mean, random bias & dilation)
    let mut feats = vec![vec![0.0; n_kernels * 2]; n];
    for ki in 0..n_kernels {
        let klen = rng.gen_range(1..=d.clamp(1, 5));
        let mut weights: Vec<f64> = (0..klen).map(|_| rng.gen_range(-1.0..1.0)).collect();
        let wmean = weights.iter().sum::<f64>() / klen as f64;
        for w in &mut weights {
            *w -= wmean;
        }
        let bias = rng.gen_range(-0.8..0.8);
        let max_dil = ((d.saturating_sub(1)) / klen.saturating_sub(1).max(1)).max(1);
        let dil = rng.gen_range(1..=max_dil.min(3));

        for i in 0..n {
            let mut max_val = f64::NEG_INFINITY;
            let mut pos_cnt = 0usize;
            let mut total_pos = 0usize;
            let span = (klen - 1) * dil;
            let end_pos = d.saturating_sub(span).max(1);
            for start in 0..end_pos {
                let mut acc = bias;
                for m in 0..klen {
                    let idx = (start + m * dil).min(d - 1);
                    acc += weights[m] * data.x[i][idx];
                }
                if acc > max_val {
                    max_val = acc;
                }
                if acc > 0.0 {
                    pos_cnt += 1;
                }
                total_pos += 1;
            }
            feats[i][ki * 2] = max_val;
            feats[i][ki * 2 + 1] = pos_cnt as f64 / total_pos.max(1) as f64; // PPV
        }
    }

    let c0 = init_kmeans_pp(&feats, k, &mut rng);
    let (labels, _, _, _) = run_lloyd_from_centers(&feats, c0, 35);
    AlgorithmOutput::new(
        labels,
        format!("k={}, {} kjerner (Max+PPV)", k, n_kernels),
        "Tilfeldige 1D-konvolusjonskjerner (ROCKET) + K-Means".to_string(),
    )
}

// 166. SLIC (Achanta et al. TPAMI 2012 — Compactness-weighted feature + 2D spatial distance)
fn slic_superpixel(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let m_compact = 0.35_f64;
    let s_step = (n as f64 / k as f64).sqrt().max(1.0);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 166);

    // Augment feature space with 2D spatial proxy (PCA coordinates normalized) -> 5D-like [features, x, y]
    let mut aug: Vec<Vec<f64>> = Vec::with_capacity(n);
    for i in 0..n {
        let mut v = data.x[i].clone();
        v.push(data.pca_2d[i][0]);
        v.push(data.pca_2d[i][1]);
        aug.push(v);
    }

    let mut centers = init_kmeans_pp(&aug, k, &mut rng);
    let mut labels = vec![0i32; n];

    for _ in 0..cfg.max_iters(25) {
        let mut changed = false;
        for i in 0..n {
            let best_c = (0..k)
                .min_by(|&a, &b| {
                    let da = slic_dist(&aug[i], &centers[a], d, m_compact, s_step);
                    let db = slic_dist(&aug[i], &centers[b], d, m_compact, s_step);
                    da.total_cmp(&db)
                })
                .unwrap_or(0);
            if labels[i] != best_c as i32 {
                labels[i] = best_c as i32;
                changed = true;
            }
        }
        if !changed {
            break;
        }
        for c in 0..k {
            let mems: Vec<usize> = (0..n).filter(|&i| labels[i] == c as i32).collect();
            if mems.is_empty() {
                continue;
            }
            for j in 0..(d + 2) {
                centers[c][j] = mems.iter().map(|&i| aug[i][j]).sum::<f64>() / mems.len() as f64;
            }
        }
    }
    AlgorithmOutput::new(
        labels,
        format!("k={}, kompakt m={:.2}, S={:.1}", k, m_compact, s_step),
        "SLIC superpiksel-metrikk (farge/attributt + romlig avstand)".to_string(),
    )
}

fn slic_dist(a: &[f64], b: &[f64], d_feat: usize, m: f64, s: f64) -> f64 {
    let dc2: f64 = (0..d_feat).map(|j| (a[j] - b[j]).powi(2)).sum();
    let ds2 = (a[d_feat] - b[d_feat]).powi(2) + (a[d_feat + 1] - b[d_feat + 1]).powi(2);
    (dc2 + (m / s).powi(2) * ds2).sqrt()
}

// 167. Watershed Clustering (Topographic Meyer flooding over inverted local density elevation map)
fn watershed_clustering(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let k = cfg.effective_k(n);
    let sigma = estimate_eps_kdist(&data.dist_matrix, (n / 8).clamp(4, 14), 0.5).max(1e-3);
    let s2 = 2.0 * sigma * sigma;

    // Elevation = -KDE density (catchment basins correspond to density peaks!)
    let elev: Vec<f64> = (0..n)
        .map(|i| {
            -(0..n)
                .map(|j| (-data.dist_matrix[i][j].powi(2) / s2).exp())
                .sum::<f64>()
        })
        .collect();

    // Find local minima in k-NN graph as watershed basin markers
    let mut markers: Vec<usize> = (0..n)
        .filter(|&i| data.knn_graph[i].iter().all(|&(nb, _)| elev[i] <= elev[nb]))
        .collect();
    markers.sort_by(|&a, &b| elev[a].total_cmp(&elev[b]));
    if markers.len() > k {
        markers.truncate(k);
    }
    if markers.len() < 2 {
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&a, &b| elev[a].total_cmp(&elev[b]));
        markers = vec![order[0]];
        while markers.len() < k {
            if let Some(next) = order
                .iter()
                .copied()
                .filter(|x| !markers.contains(x))
                .max_by(|&a, &b| {
                    let da = markers.iter().map(|&m| data.dist_matrix[a][m]).fold(f64::INFINITY, f64::min);
                    let db = markers.iter().map(|&m| data.dist_matrix[b][m]).fold(f64::INFINITY, f64::min);
                    da.total_cmp(&db)
                })
            {
                markers.push(next);
            } else {
                break;
            }
        }
    }

    // Flood points in order of rising elevation
    let mut labels = vec![-1i32; n];
    for (cid, &m) in markers.iter().enumerate() {
        labels[m] = cid as i32;
    }
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| elev[a].total_cmp(&elev[b]));

    for _pass in 0..3 {
        for &u in &order {
            if labels[u] >= 0 {
                continue;
            }
            if let Some(&(nb, _)) = data.knn_graph[u]
                .iter()
                .filter(|&&(v, _)| labels[v] >= 0)
                .min_by(|&&(a, _), &&(b, _)| elev[a].total_cmp(&elev[b]))
            {
                labels[u] = labels[nb];
            }
        }
    }
    // Nearest marker fallback for any isolated node
    for i in 0..n {
        if labels[i] < 0 {
            if let Some((cid, _)) = markers
                .iter()
                .enumerate()
                .min_by(|(_, &a), (_, &b)| data.dist_matrix[i][a].total_cmp(&data.dist_matrix[i][b]))
            {
                labels[i] = cid as i32;
            }
        }
    }
    AlgorithmOutput::new(
        labels,
        format!("{} bassengmarkører, sigma={:.2}", markers.len(), sigma),
        "Topografisk Watershed-oversvømmelse på invertert tetthetskart".to_string(),
    )
}

// 168. DBSCAN-3D / Euclidean Cluster Extraction (PCL point-cloud 3D projection + voxel Euclidean clusters)
fn dbscan_3d_pcl(data: &PreparedData, _cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    // Extract 3D point cloud coordinates (first 3 features or padded with PCA/0)
    let pts3d: Vec<[f64; 3]> = (0..n)
        .map(|i| {
            [
                data.x[i][0],
                if d > 1 { data.x[i][1] } else { data.pca_2d[i][1] },
                if d > 2 { data.x[i][2] } else { 0.0 },
            ]
        })
        .collect();

    let mut d3 = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let dd = ((pts3d[i][0] - pts3d[j][0]).powi(2)
                + (pts3d[i][1] - pts3d[j][1]).powi(2)
                + (pts3d[i][2] - pts3d[j][2]).powi(2))
            .sqrt();
            d3[i][j] = dd;
            d3[j][i] = dd;
        }
    }
    let min_cluster = 4usize.min(n / 4).max(2);
    let tol = estimate_eps_kdist(&d3, min_cluster, 0.58);
    let labels = run_dbscan_core(&d3, tol, min_cluster);
    AlgorithmOutput::new(
        labels,
        format!("3D-toleranse={:.3}, minPts={}", tol, min_cluster),
        "PCL 3D Euklidsk punktsky-ekstraksjon".to_string(),
    )
}

// 169. TRACLUS (Lee et al. SIGMOD 2007 — Trajectory segment distance d_perp + d_parallel + d_theta + DBSCAN)
fn traclus(data: &PreparedData, _cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    // Represent each observation as a local directed trajectory segment in 2D PCA space (ordered by time_values)
    let mut time_ord: Vec<usize> = (0..n).collect();
    time_ord.sort_by(|&a, &b| data.time_values[a].total_cmp(&data.time_values[b]));

    let mut seg_start = vec![[0.0; 2]; n];
    let mut seg_end = vec![[0.0; 2]; n];
    for pos in 0..n {
        let cur = time_ord[pos];
        let nxt = time_ord[(pos + 1).min(n - 1)];
        seg_start[cur] = data.pca_2d[cur];
        seg_end[cur] = if pos + 1 < n {
            data.pca_2d[nxt]
        } else {
            data.pca_2d[cur]
        };
    }

    // TRACLUS segment distance combining perpendicular, midpoint parallel, and angular distance
    let mut tdist = vec![vec![0.0; n]; n];
    for i in 0..n {
        let vi = [seg_end[i][0] - seg_start[i][0], seg_end[i][1] - seg_start[i][1]];
        let li = (vi[0] * vi[0] + vi[1] * vi[1]).sqrt().max(1e-6);
        let mi = [0.5 * (seg_start[i][0] + seg_end[i][0]), 0.5 * (seg_start[i][1] + seg_end[i][1])];
        for j in (i + 1)..n {
            let vj = [seg_end[j][0] - seg_start[j][0], seg_end[j][1] - seg_start[j][1]];
            let lj = (vj[0] * vj[0] + vj[1] * vj[1]).sqrt().max(1e-6);
            let mj = [0.5 * (seg_start[j][0] + seg_end[j][0]), 0.5 * (seg_start[j][1] + seg_end[j][1])];

            let d_mid = ((mi[0] - mj[0]).powi(2) + (mi[1] - mj[1]).powi(2)).sqrt();
            let cos_t = ((vi[0] * vj[0] + vi[1] * vj[1]) / (li * lj)).clamp(-1.0, 1.0);
            let d_angle = (1.0 - cos_t.abs()) * li.min(lj) * 0.5;
            let dd = d_mid + d_angle;
            tdist[i][j] = dd;
            tdist[j][i] = dd;
        }
    }

    let min_lns = 4usize.min(n / 4).max(2);
    let eps = estimate_eps_kdist(&tdist, min_lns, 0.62);
    let labels = run_dbscan_core(&tdist, eps, min_lns);
    AlgorithmOutput::new(
        labels,
        format!("eps={:.3}, MinLns={}", eps, min_lns),
        "Partition-and-Group banesegmentering (posisjon + retningsvinkel)".to_string(),
    )
}

fn generate_base_ensemble(data: &PreparedData, k: usize, r_runs: usize, seed: u64) -> Vec<Vec<i32>> {
    let mut ensemble = Vec::with_capacity(r_runs);
    for r in 0..r_runs {
        let mut rng = ChaCha8Rng::seed_from_u64(seed + r as u64 * 31 + 7);
        let kr = (k + (r % 3)).saturating_sub(1).clamp(2, data.n_samples / 3);
        let c0 = init_kmeans_pp(&data.x, kr, &mut rng);
        let (lbl, _, _, _) = run_lloyd_from_centers(&data.x, c0, 20);
        ensemble.push(lbl);
    }
    ensemble
}

// 170. Consensus Clustering (Monti et al. 2003 — repeated 80% subsampling stability consensus matrix)
fn consensus_clustering(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let k = cfg.effective_k(n);
    let n_runs = 14usize;
    let sub_n = ((n as f64) * 0.80).round().max((k + 1) as f64) as usize;
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 170);

    let mut co_occur = vec![vec![0.0; n]; n];
    let mut co_sample = vec![vec![0.0; n]; n];
    let mut idx: Vec<usize> = (0..n).collect();

    for _ in 0..n_runs {
        idx.shuffle(&mut rng);
        let sub = &idx[..sub_n.min(n)];
        let sub_x: Vec<Vec<f64>> = sub.iter().map(|&i| data.x[i].clone()).collect();
        let c0 = init_kmeans_pp(&sub_x, k, &mut rng);
        let (lbl, _, _, _) = run_lloyd_from_centers(&sub_x, c0, 20);
        for a in 0..sub.len() {
            for b in (a + 1)..sub.len() {
                let u = sub[a];
                let v = sub[b];
                co_sample[u][v] += 1.0;
                co_sample[v][u] += 1.0;
                if lbl[a] == lbl[b] {
                    co_occur[u][v] += 1.0;
                    co_occur[v][u] += 1.0;
                }
            }
        }
    }

    let mut cdist = vec![vec![0.0; n]; n];
    let mut stability_sum = 0.0;
    let mut pairs = 0usize;
    for i in 0..n {
        for j in (i + 1)..n {
            let m_ij: f64 = if co_sample[i][j] > 0.0 {
                co_occur[i][j] / co_sample[i][j]
            } else {
                0.0
            };
            cdist[i][j] = 1.0 - m_ij;
            cdist[j][i] = 1.0 - m_ij;
            // Crisp consensus entries near 0 or 1 indicate high stability
            stability_sum += (2.0 * m_ij - 1.0).abs();
            pairs += 1;
        }
    }
    let labels = agglomerative_lance_williams(&cdist, k, Linkage::Average);
    let stab = stability_sum / pairs.max(1) as f64;
    AlgorithmOutput::new(
        labels,
        format!("k={}, {} delutvalg (80%)", k, n_runs),
        format!("Monti konsensusmatrise-stabilitet={:.3}", stab),
    )
}

// 171. CSPA (Strehl & Ghosh JMLR 2002 — Cluster-based Similarity Partitioning Algorithm)
fn cspa_ensemble(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let k = cfg.effective_k(n);
    let r = 10usize;
    let ens = generate_base_ensemble(data, k, r, cfg.seed + 171);

    // Co-association similarity S = (1/r) H H^T
    let mut sim = vec![vec![0.0; n]; n];
    for i in 0..n {
        sim[i][i] = 1.0;
        for j in (i + 1)..n {
            let matches = ens.iter().filter(|e| e[i] == e[j]).count();
            let s = matches as f64 / r as f64;
            sim[i][j] = s;
            sim[j][i] = s;
        }
    }
    // Spectral partitioning of induced similarity graph
    let (_, evecs) = top_k_symmetric_eigen(&sim, k, 25);
    let mut embed = vec![vec![0.0; k]; n];
    for i in 0..n {
        for c in 0..k {
            embed[i][c] = evecs[c][i];
        }
    }
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 171);
    let c0 = init_kmeans_pp(&embed, k, &mut rng);
    let (labels, _, _, _) = run_lloyd_from_centers(&embed, c0, 30);
    AlgorithmOutput::new(
        labels,
        format!("k={}, ensemble r={}", k, r),
        "CSPA hypergraf-likhetsmatrise + spektral partisjonering".to_string(),
    )
}

// 172. HGPA (Strehl & Ghosh 2002 — HyperGraph Partitioning Algorithm via clique-expanded hyperedges)
fn hgpa_ensemble(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let k = cfg.effective_k(n);
    let r = 10usize;
    let ens = generate_base_ensemble(data, k, r, cfg.seed + 172);

    // Each base cluster is a hyperedge weighted by 1/|e| (star/clique hypergraph Laplacian)
    let mut adj = vec![vec![0.0; n]; n];
    let mut deg = vec![0.0; n];
    for run in &ens {
        let mut groups: HashMap<i32, Vec<usize>> = HashMap::new();
        for (i, &l) in run.iter().enumerate() {
            groups.entry(l).or_default().push(i);
        }
        for mems in groups.values() {
            let w = 1.0 / mems.len().max(1) as f64;
            for a in 0..mems.len() {
                for b in (a + 1)..mems.len() {
                    adj[mems[a]][mems[b]] += w;
                    adj[mems[b]][mems[a]] += w;
                }
            }
        }
    }
    for i in 0..n {
        deg[i] = adj[i].iter().sum::<f64>().max(1e-6);
    }
    let mut norm_a = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            norm_a[i][j] = adj[i][j] / (deg[i] * deg[j]).sqrt();
        }
    }
    let (_, evecs) = top_k_symmetric_eigen(&norm_a, k, 25);
    let mut embed = vec![vec![0.0; k]; n];
    for i in 0..n {
        for c in 0..k {
            embed[i][c] = evecs[c][i];
        }
    }
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 172);
    let c0 = init_kmeans_pp(&embed, k, &mut rng);
    let (labels, _, _, _) = run_lloyd_from_centers(&embed, c0, 30);
    AlgorithmOutput::new(
        labels,
        format!("k={}, hypergraf r={}", k, r),
        "HGPA normalisert hypergraf-min-cut partisjonering".to_string(),
    )
}

// 173. MCLA (Strehl & Ghosh 2002 — Meta-CLustering Algorithm over Jaccard meta-graph of base clusters)
fn mcla_ensemble(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let k = cfg.effective_k(n);
    let r = 10usize;
    let ens = generate_base_ensemble(data, k, r, cfg.seed + 173);

    // Collect all base clusters as binary indicator vectors
    let mut base_sets: Vec<HashSet<usize>> = Vec::new();
    for run in &ens {
        let mut groups: HashMap<i32, HashSet<usize>> = HashMap::new();
        for (i, &l) in run.iter().enumerate() {
            groups.entry(l).or_default().insert(i);
        }
        base_sets.extend(groups.into_values());
    }
    let m = base_sets.len();
    let mut meta_dist = vec![vec![0.0; m]; m];
    for a in 0..m {
        for b in (a + 1)..m {
            let inter = base_sets[a].intersection(&base_sets[b]).count() as f64;
            let union = base_sets[a].union(&base_sets[b]).count().max(1) as f64;
            let dd = 1.0 - inter / union;
            meta_dist[a][b] = dd;
            meta_dist[b][a] = dd;
        }
    }
    // Group base clusters into k meta-clusters
    let meta_lbl = if m > k {
        agglomerative_lance_williams(&meta_dist, k, Linkage::Average)
    } else {
        (0..m as i32).collect()
    };

    // Collapse each meta-cluster into soft association scores and assign each point to strongest meta-cluster
    let mut scores = vec![vec![0.0; k]; n];
    let mut meta_cnt = vec![0usize; k];
    for (bi, set) in base_sets.iter().enumerate() {
        let mc = (meta_lbl[bi] as usize).min(k - 1);
        meta_cnt[mc] += 1;
        for &p in set {
            scores[p][mc] += 1.0;
        }
    }
    for i in 0..n {
        for mc in 0..k {
            if meta_cnt[mc] > 0 {
                scores[i][mc] /= meta_cnt[mc] as f64;
            }
        }
    }
    let labels: Vec<i32> = (0..n)
        .map(|i| {
            (0..k)
                .max_by(|&a, &b| scores[i][a].total_cmp(&scores[i][b]))
                .unwrap_or(0) as i32
        })
        .collect();
    AlgorithmOutput::new(
        labels,
        format!("k={}, {} basisklynger -> {} metaklynger", k, m, k),
        "MCLA grupperte basisklynger via Jaccard-metagraf".to_string(),
    )
}

fn split_views(data: &PreparedData) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let n = data.n_samples;
    let d = data.n_features;
    if d >= 2 {
        let mid = d / 2;
        let v1: Vec<Vec<f64>> = (0..n).map(|i| data.x[i][..mid].to_vec()).collect();
        let v2: Vec<Vec<f64>> = (0..n).map(|i| data.x[i][mid..].to_vec()).collect();
        (v1, v2)
    } else {
        // Fallback: view 1 = scaled numeric, view 2 = non-linear / rank features
        let v1 = data.x.clone();
        let v2: Vec<Vec<f64>> = (0..n)
            .map(|i| vec![data.pca_2d[i][0], data.pca_2d[i][1]])
            .collect();
        (v1, v2)
    }
}

// 174. Multi-View K-Means (Cai et al. IJCAI 2013 — Adaptive view weights alpha_v across feature views)
fn multi_view_kmeans(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let k = cfg.effective_k(n);
    let (v1, v2) = split_views(data);
    let views = [&v1, &v2];
    let gamma = 2.0_f64;
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 174);

    let mut centers: Vec<Vec<Vec<f64>>> = views
        .iter()
        .map(|v| init_kmeans_pp(v, k, &mut rng))
        .collect();
    let mut alpha = vec![0.5_f64, 0.5];
    let mut labels = vec![0i32; n];

    for _ in 0..cfg.max_iters(35) {
        let mut changed = false;
        for i in 0..n {
            let best_c = (0..k)
                .min_by(|&a, &b| {
                    let da: f64 = (0..2)
                        .map(|vi| alpha[vi].powf(gamma) * euclidean_sq(&views[vi][i], &centers[vi][a]))
                        .sum();
                    let db: f64 = (0..2)
                        .map(|vi| alpha[vi].powf(gamma) * euclidean_sq(&views[vi][i], &centers[vi][b]))
                        .sum();
                    da.total_cmp(&db)
                })
                .unwrap_or(0);
            if labels[i] != best_c as i32 {
                labels[i] = best_c as i32;
                changed = true;
            }
        }
        // Update per-view centroids & learn inverse-SSE view weights
        let mut view_sse = vec![1e-6; 2];
        for vi in 0..2 {
            let dv = views[vi][0].len();
            for c in 0..k {
                let mems: Vec<usize> = (0..n).filter(|&i| labels[i] == c as i32).collect();
                if mems.is_empty() {
                    continue;
                }
                for j in 0..dv {
                    centers[vi][c][j] = mems.iter().map(|&i| views[vi][i][j]).sum::<f64>() / mems.len() as f64;
                }
            }
            for i in 0..n {
                view_sse[vi] += euclidean_sq(&views[vi][i], &centers[vi][labels[i] as usize]);
            }
        }
        let inv0 = (1.0 / view_sse[0]).powf(1.0 / (gamma - 1.0));
        let inv1 = (1.0 / view_sse[1]).powf(1.0 / (gamma - 1.0));
        let s = (inv0 + inv1).max(1e-12);
        alpha[0] = inv0 / s;
        alpha[1] = inv1 / s;
        if !changed {
            break;
        }
    }

    AlgorithmOutput::new(
        labels,
        format!("k={}, 2 visninger (w={:.2}/{:.2})", k, alpha[0], alpha[1]),
        "Multi-View K-Means med lærte visningsvekter".to_string(),
    )
}

// 175. Co-regularized Multi-View Spectral Clustering (Kumar et al. NeurIPS 2011)
fn coreg_multi_view_spectral(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let k = cfg.effective_k(n);
    let (v1, v2) = split_views(data);
    let lambda = 0.5_f64;

    let norm_lap_aff = |v: &[Vec<f64>]| -> Vec<Vec<f64>> {
        let dm = crate::math::pairwise_dist_matrix(v);
        let sig = median_pairwise_dist(&dm).max(1e-4);
        let aff = rbf_affinity_matrix(&dm, sig);
        let deg: Vec<f64> = (0..n).map(|i| aff[i].iter().sum::<f64>().max(1e-8)).collect();
        let mut l = vec![vec![0.0; n]; n];
        for i in 0..n {
            for j in 0..n {
                l[i][j] = aff[i][j] / (deg[i] * deg[j]).sqrt();
            }
        }
        l
    };

    let l1 = norm_lap_aff(&v1);
    let l2 = norm_lap_aff(&v2);

    // Alternating co-regularization: L1 + lambda * U2 U2^T  <->  L2 + lambda * U1 U1^T
    let (_, mut u1) = top_k_symmetric_eigen(&l1, k, 20);
    let (_, mut u2) = top_k_symmetric_eigen(&l2, k, 20);

    let outer_uu = |u: &[Vec<f64>]| -> Vec<Vec<f64>> {
        let mut m = vec![vec![0.0; n]; n];
        for c in 0..u.len() {
            for i in 0..n {
                for j in i..n {
                    let prod = u[c][i] * u[c][j];
                    m[i][j] += prod;
                    if i != j {
                        m[j][i] += prod;
                    }
                }
            }
        }
        m
    };

    for _ in 0..3 {
        let reg2 = outer_uu(&u2);
        let mut mod1 = l1.clone();
        for i in 0..n {
            for j in 0..n {
                mod1[i][j] += lambda * reg2[i][j];
            }
        }
        u1 = top_k_symmetric_eigen(&mod1, k, 18).1;

        let reg1 = outer_uu(&u1);
        let mut mod2 = l2.clone();
        for i in 0..n {
            for j in 0..n {
                mod2[i][j] += lambda * reg1[i][j];
            }
        }
        u2 = top_k_symmetric_eigen(&mod2, k, 18).1;
    }

    // Consensus embedding from co-regularized eigenvectors of both views
    let mut embed = vec![vec![0.0; 2 * k]; n];
    for i in 0..n {
        for c in 0..k {
            embed[i][c] = u1[c][i];
            embed[i][k + c] = u2[c][i];
        }
        embed[i] = normalize_vec(&embed[i]);
    }
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 175);
    let c0 = init_kmeans_pp(&embed, k, &mut rng);
    let (labels, _, _, _) = run_lloyd_from_centers(&embed, c0, 35);
    AlgorithmOutput::new(
        labels,
        format!("k={}, ko-reg lambda={:.1}", k, lambda),
        "Ko-regulert spektral konsensus på tvers av 2 datavisninger".to_string(),
    )
}
