use std::f64::consts::PI;

#[inline]
pub fn euclidean_sq(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| {
            let d = x - y;
            d * d
        })
        .sum()
}

#[inline]
pub fn euclidean(a: &[f64], b: &[f64]) -> f64 {
    euclidean_sq(a, b).sqrt()
}

#[inline]
pub fn manhattan(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| (x - y).abs()).sum()
}

#[inline]
pub fn weighted_minkowski(a: &[f64], b: &[f64], weights: &[f64], p: f64) -> f64 {
    let sum: f64 = a
        .iter()
        .zip(b.iter())
        .zip(weights.iter())
        .map(|((x, y), &w)| (w * (x - y).abs()).powf(p))
        .sum();
    sum.powf(1.0 / p)
}

#[inline]
pub fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

#[inline]
pub fn norm(a: &[f64]) -> f64 {
    dot(a, a).sqrt()
}

#[inline]
pub fn normalize_vec(a: &[f64]) -> Vec<f64> {
    let n = norm(a);
    if n > 1e-12 {
        a.iter().map(|&v| v / n).collect()
    } else {
        vec![0.0; a.len()]
    }
}

#[inline]
pub fn cosine_sim(a: &[f64], b: &[f64]) -> f64 {
    let na = norm(a);
    let nb = norm(b);
    if na > 1e-12 && nb > 1e-12 {
        (dot(a, b) / (na * nb)).clamp(-1.0, 1.0)
    } else {
        0.0
    }
}

#[allow(dead_code)]
#[inline]
pub fn cosine_dist(a: &[f64], b: &[f64]) -> f64 {
    (1.0 - cosine_sim(a, b)).max(0.0)
}

#[inline]
pub fn hamming_dist(a: &[usize], b: &[usize]) -> usize {
    a.iter().zip(b.iter()).filter(|(x, y)| x != y).count()
}

#[inline]
pub fn jaccard_sim_cat(a: &[usize], b: &[usize]) -> f64 {
    if a.is_empty() {
        return 0.0;
    }
    let matches = a.iter().zip(b.iter()).filter(|(x, y)| x == y).count();
    let union = 2 * a.len() - matches;
    if union == 0 {
        0.0
    } else {
        matches as f64 / union as f64
    }
}

pub fn pairwise_dist_matrix(x: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = x.len();
    let mut dist = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let d = euclidean(&x[i], &x[j]);
            dist[i][j] = d;
            dist[j][i] = d;
        }
    }
    dist
}

pub fn median_pairwise_dist(dist: &[Vec<f64>]) -> f64 {
    let n = dist.len();
    if n < 2 {
        return 1.0;
    }
    let mut vals = Vec::with_capacity(n * (n - 1) / 2);
    for i in 0..n {
        for j in (i + 1)..n {
            if dist[i][j] > 0.0 {
                vals.push(dist[i][j]);
            }
        }
    }
    if vals.is_empty() {
        return 1.0;
    }
    vals.sort_by(|a, b| a.total_cmp(b));
    vals[vals.len() / 2].max(1e-6)
}

pub fn estimate_eps_kdist(dist: &[Vec<f64>], k: usize, quantile: f64) -> f64 {
    let n = dist.len();
    if n < 2 {
        return 1.0;
    }
    let kk = k.clamp(1, n - 1);
    let mut kdists = Vec::with_capacity(n);
    for row in dist {
        let mut sorted = row.clone();
        sorted.sort_by(|a, b| a.total_cmp(b));
        kdists.push(sorted[kk]);
    }
    kdists.sort_by(|a, b| a.total_cmp(b));
    let idx = ((kdists.len() as f64 - 1.0) * quantile.clamp(0.05, 0.98)).round() as usize;
    kdists[idx.min(kdists.len() - 1)].max(1e-5)
}

/// Symmetric Gaussian-weighted k-NN graph adjacency list: adj[i] = Vec<(neighbor, weight)>
pub fn build_knn_graph(dist: &[Vec<f64>], k: usize) -> Vec<Vec<(usize, f64)>> {
    let n = dist.len();
    if n == 0 {
        return vec![];
    }
    let kk = k.clamp(1, n.saturating_sub(1).max(1));
    let sigma = estimate_eps_kdist(dist, kk, 0.5).max(1e-5);
    let denom = 2.0 * sigma * sigma;

    let mut weights = vec![vec![0.0; n]; n];
    for i in 0..n {
        let mut nbrs: Vec<(usize, f64)> = (0..n)
            .filter(|&j| j != i)
            .map(|j| (j, dist[i][j]))
            .collect();
        nbrs.sort_by(|a, b| a.1.total_cmp(&b.1));
        for &(j, d) in nbrs.iter().take(kk) {
            let w = (-(d * d) / denom).exp().max(1e-4);
            if w > weights[i][j] {
                weights[i][j] = w;
                weights[j][i] = w;
            }
        }
    }

    let mut adj = vec![Vec::new(); n];
    for i in 0..n {
        for j in 0..n {
            if weights[i][j] > 0.0 {
                adj[i].push((j, weights[i][j]));
            }
        }
    }
    adj
}

pub fn rbf_affinity_matrix(dist: &[Vec<f64>], sigma: f64) -> Vec<Vec<f64>> {
    let n = dist.len();
    let s2 = 2.0 * sigma.max(1e-5).powi(2);
    let mut aff = vec![vec![0.0; n]; n];
    for i in 0..n {
        aff[i][i] = 1.0;
        for j in (i + 1)..n {
            let w = (-(dist[i][j] * dist[i][j]) / s2).exp();
            aff[i][j] = w;
            aff[j][i] = w;
        }
    }
    aff
}

/// Exact symmetric Jacobi eigendecomposition for small/medium symmetric matrices (d <= 64).
/// Returns (eigenvalues sorted descending, eigenvectors as columns: vecs[eig_idx][dim]).
pub fn jacobi_eigen(mat: &[Vec<f64>], max_sweeps: usize) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = mat.len();
    if n == 0 {
        return (vec![], vec![]);
    }
    let mut a = mat.to_vec();
    let mut v = vec![vec![0.0; n]; n];
    for i in 0..n {
        v[i][i] = 1.0;
    }

    let max_iter = max_sweeps * n * n.max(1);
    for _ in 0..max_iter.min(50_000) {
        let mut p = 0;
        let mut q = 0;
        let mut max_off = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                let val = a[i][j].abs();
                if val > max_off {
                    max_off = val;
                    p = i;
                    q = j;
                }
            }
        }
        if max_off < 1e-10 {
            break;
        }
        let app = a[p][p];
        let aqq = a[q][q];
        let apq = a[p][q];
        let tau = (aqq - app) / (2.0 * apq);
        let t = if tau >= 0.0 {
            1.0 / (tau + (1.0 + tau * tau).sqrt())
        } else {
            -1.0 / (-tau + (1.0 + tau * tau).sqrt())
        };
        let c = 1.0 / (1.0 + t * t).sqrt();
        let s = t * c;

        a[p][p] = app - t * apq;
        a[q][q] = aqq + t * apq;
        a[p][q] = 0.0;
        a[q][p] = 0.0;

        for r in 0..n {
            if r != p && r != q {
                let arp = a[r][p];
                let arq = a[r][q];
                a[r][p] = c * arp - s * arq;
                a[p][r] = a[r][p];
                a[r][q] = s * arp + c * arq;
                a[q][r] = a[r][q];
            }
        }
        for r in 0..n {
            let vrp = v[r][p];
            let vrq = v[r][q];
            v[r][p] = c * vrp - s * vrq;
            v[r][q] = s * vrp + c * vrq;
        }
    }

    let mut pairs: Vec<(f64, Vec<f64>)> = (0..n)
        .map(|i| {
            let col: Vec<f64> = (0..n).map(|r| v[r][i]).collect();
            (a[i][i], col)
        })
        .collect();
    pairs.sort_by(|x, y| y.0.total_cmp(&x.0));
    let vals = pairs.iter().map(|p| p.0).collect();
    let vecs = pairs.into_iter().map(|p| p.1).collect();
    (vals, vecs)
}

/// Fast orthogonal subspace iteration for top-k eigenvectors of an n x n symmetric matrix.
/// Returns (eigenvalues descending, eigenvectors[eig_idx][row]).
pub fn top_k_symmetric_eigen(
    mat: &[Vec<f64>],
    k: usize,
    max_iter: usize,
) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = mat.len();
    if n == 0 || k == 0 {
        return (vec![], vec![]);
    }
    let kk = k.min(n);
    if n <= 24 {
        let (vals, vecs) = jacobi_eigen(mat, 30);
        return (vals.into_iter().take(kk).collect(), vecs.into_iter().take(kk).collect());
    }

    // Shift matrix so eigenvalues are positive-biased for reliable top-k convergence
    let gershgorin_shift: f64 = mat
        .iter()
        .enumerate()
        .map(|(i, row)| {
            let off: f64 = row
                .iter()
                .enumerate()
                .filter(|&(j, _)| j != i)
                .map(|(_, v)| v.abs())
                .sum();
            (off - row[i]).max(0.0)
        })
        .fold(0.0_f64, f64::max)
        + 1.0;

    // Deterministic quasi-random initial basis
    let mut basis = vec![vec![0.0; n]; kk];
    for b in 0..kk {
        for i in 0..n {
            let angle = ((i + 1) as f64) * ((b + 1) as f64) * 0.61803398875;
            basis[b][i] = (angle * PI).sin() + 0.13 * ((i * 7 + b * 13) % 19) as f64;
        }
    }
    orthonormalize(&mut basis);

    let mut next = vec![vec![0.0; n]; kk];
    for _ in 0..max_iter {
        for b in 0..kk {
            for i in 0..n {
                let mut sum = gershgorin_shift * basis[b][i];
                let row = &mat[i];
                for j in 0..n {
                    sum += row[j] * basis[b][j];
                }
                next[b][i] = sum;
            }
        }
        orthonormalize(&mut next);
        std::mem::swap(&mut basis, &mut next);
    }

    // Rayleigh-Ritz projection onto the kk-dimensional subspace
    let mut proj = vec![vec![0.0; kk]; kk];
    let mut av = vec![vec![0.0; n]; kk];
    for b in 0..kk {
        for i in 0..n {
            let mut sum = 0.0;
            let row = &mat[i];
            for j in 0..n {
                sum += row[j] * basis[b][j];
            }
            av[b][i] = sum;
        }
    }
    for i in 0..kk {
        for j in i..kk {
            let val = dot(&basis[i], &av[j]);
            proj[i][j] = val;
            proj[j][i] = val;
        }
    }
    let (eigvals, eigvecs_small) = jacobi_eigen(&proj, 30);
    let mut final_vecs = vec![vec![0.0; n]; kk];
    for c in 0..kk {
        for b in 0..kk {
            let coeff = eigvecs_small[c][b];
            for i in 0..n {
                final_vecs[c][i] += coeff * basis[b][i];
            }
        }
        let nn = norm(&final_vecs[c]);
        if nn > 1e-12 {
            for v in &mut final_vecs[c] {
                *v /= nn;
            }
        }
    }
    (eigvals, final_vecs)
}

/// Compute the k smallest eigenvalues & eigenvectors of symmetric matrix mat.
pub fn bottom_k_symmetric_eigen(
    mat: &[Vec<f64>],
    k: usize,
    max_iter: usize,
) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = mat.len();
    if n == 0 || k == 0 {
        return (vec![], vec![]);
    }
    let mut neg = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            neg[i][j] = -mat[i][j];
        }
    }
    let (vals, vecs) = top_k_symmetric_eigen(&neg, k, max_iter);
    let pos_vals = vals.into_iter().map(|v| -v).collect();
    (pos_vals, vecs)
}

pub fn orthonormalize(vecs: &mut [Vec<f64>]) {
    let m = vecs.len();
    if m == 0 {
        return;
    }
    let n = vecs[0].len();
    for i in 0..m {
        // Modified Gram-Schmidt twice for numerical stability
        for _ in 0..2 {
            for j in 0..i {
                let r = dot(&vecs[i], &vecs[j]);
                for k in 0..n {
                    vecs[i][k] -= r * vecs[j][k];
                }
            }
        }
        let nn = norm(&vecs[i]);
        if nn > 1e-12 {
            for k in 0..n {
                vecs[i][k] /= nn;
            }
        } else {
            // Fallback canonical direction
            vecs[i].fill(0.0);
            vecs[i][i % n] = 1.0;
        }
    }
}

/// 2D PCA projection returning coordinates and explained variance ratio [pc1, pc2].
pub fn pca_2d(data: &[Vec<f64>]) -> (Vec<[f64; 2]>, [f64; 2]) {
    let n = data.len();
    if n == 0 {
        return (vec![], [0.0, 0.0]);
    }
    let d = data[0].len();
    if d == 0 {
        return (vec![[0.0, 0.0]; n], [0.0, 0.0]);
    }
    if d == 1 {
        let mean: f64 = data.iter().map(|r| r[0]).sum::<f64>() / n as f64;
        let coords = data
            .iter()
            .enumerate()
            .map(|(i, r)| [r[0] - mean, (i as f64) / (n.max(1) as f64) - 0.5])
            .collect();
        return (coords, [1.0, 0.0]);
    }

    let mut mean = vec![0.0; d];
    for row in data {
        for j in 0..d {
            mean[j] += row[j];
        }
    }
    for v in &mut mean {
        *v /= n as f64;
    }

    let centered: Vec<Vec<f64>> = data
        .iter()
        .map(|r| r.iter().enumerate().map(|(j, &v)| v - mean[j]).collect())
        .collect();

    let mut cov = vec![vec![0.0; d]; d];
    let denom = (n.saturating_sub(1)).max(1) as f64;
    for row in &centered {
        for i in 0..d {
            for j in i..d {
                cov[i][j] += row[i] * row[j];
            }
        }
    }
    for i in 0..d {
        for j in i..d {
            cov[i][j] /= denom;
            cov[j][i] = cov[i][j];
        }
    }

    let (vals, vecs) = if d <= 40 {
        jacobi_eigen(&cov, 40)
    } else {
        top_k_symmetric_eigen(&cov, 2, 40)
    };

    let total_var: f64 = (0..d).map(|i| cov[i][i].max(0.0)).sum::<f64>().max(1e-12);
    let v1 = vecs.first().cloned().unwrap_or_else(|| {
        let mut e = vec![0.0; d];
        e[0] = 1.0;
        e
    });
    let v2 = vecs.get(1).cloned().unwrap_or_else(|| {
        let mut e = vec![0.0; d];
        if d > 1 {
            e[1] = 1.0;
        }
        e
    });

    let coords = centered
        .iter()
        .map(|row| [dot(row, &v1), dot(row, &v2)])
        .collect();

    let ev1 = (vals.first().copied().unwrap_or(0.0).max(0.0) / total_var).clamp(0.0, 1.0);
    let ev2 = (vals.get(1).copied().unwrap_or(0.0).max(0.0) / total_var).clamp(0.0, 1.0);
    (coords, [ev1, ev2])
}

/// Regularized Cholesky inverse and log-determinant for symmetric PSD covariance matrix.
pub fn inv_and_logdet(cov: &[Vec<f64>], reg: f64) -> (Vec<Vec<f64>>, f64) {
    let d = cov.len();
    let mut a = cov.to_vec();
    for i in 0..d {
        a[i][i] += reg;
    }

    // Cholesky L * L^T = A
    let mut l = vec![vec![0.0; d]; d];
    for i in 0..d {
        for j in 0..=i {
            let mut sum = a[i][j];
            for k in 0..j {
                sum -= l[i][k] * l[j][k];
            }
            if i == j {
                if sum <= 1e-10 {
                    // Diagonal fallback
                    let mut inv = vec![vec![0.0; d]; d];
                    let mut logdet = 0.0;
                    for m in 0..d {
                        let diag = (cov[m][m] + reg).max(1e-4);
                        inv[m][m] = 1.0 / diag;
                        logdet += diag.ln();
                    }
                    return (inv, logdet);
                }
                l[i][j] = sum.sqrt();
            } else {
                l[i][j] = sum / l[j][j];
            }
        }
    }

    let logdet: f64 = 2.0 * (0..d).map(|i| l[i][i].ln()).sum::<f64>();

    // Invert lower-triangular L
    let mut linv = vec![vec![0.0; d]; d];
    for i in 0..d {
        linv[i][i] = 1.0 / l[i][i];
        for j in (i + 1)..d {
            let mut sum = 0.0;
            for k in i..j {
                sum -= l[j][k] * linv[k][i];
            }
            linv[j][i] = sum / l[j][j];
        }
    }

    // A^-1 = L^-T * L^-1
    let mut inv = vec![vec![0.0; d]; d];
    for i in 0..d {
        for j in i..d {
            let mut sum = 0.0;
            for k in j..d {
                sum += linv[k][i] * linv[k][j];
            }
            inv[i][j] = sum;
            inv[j][i] = sum;
        }
    }
    (inv, logdet)
}

pub fn mahalanobis_sq(x: &[f64], mean: &[f64], inv_cov: &[Vec<f64>]) -> f64 {
    let d = x.len();
    let diff: Vec<f64> = (0..d).map(|i| x[i] - mean[i]).collect();
    let mut sum = 0.0;
    for i in 0..d {
        let mut row_s = 0.0;
        for j in 0..d {
            row_s += inv_cov[i][j] * diff[j];
        }
        sum += diff[i] * row_s;
    }
    sum.max(0.0)
}

/// Dynamic Time Warping distance and optimal warping path.
pub fn dtw_distance_and_path(a: &[f64], b: &[f64]) -> (f64, Vec<(usize, usize)>) {
    let m = a.len();
    let n = b.len();
    if m == 0 || n == 0 {
        return (0.0, vec![]);
    }
    let mut dp = vec![vec![f64::INFINITY; n + 1]; m + 1];
    dp[0][0] = 0.0;
    for i in 1..=m {
        for j in 1..=n {
            let cost = (a[i - 1] - b[j - 1]).powi(2);
            let prev = dp[i - 1][j].min(dp[i][j - 1]).min(dp[i - 1][j - 1]);
            dp[i][j] = cost + prev;
        }
    }
    let mut path = Vec::new();
    let mut i = m;
    let mut j = n;
    while i > 0 && j > 0 {
        path.push((i - 1, j - 1));
        let diag = dp[i - 1][j - 1];
        let up = dp[i - 1][j];
        let left = dp[i][j - 1];
        if diag <= up && diag <= left {
            i -= 1;
            j -= 1;
        } else if up <= left {
            i -= 1;
        } else {
            j -= 1;
        }
    }
    path.reverse();
    (dp[m][n].sqrt(), path)
}

/// Soft-DTW distance (Cuturi & Blondel 2017) with smoothing gamma > 0.
pub fn soft_dtw(a: &[f64], b: &[f64], gamma: f64) -> f64 {
    let m = a.len();
    let n = b.len();
    if m == 0 || n == 0 {
        return 0.0;
    }
    let g = gamma.max(1e-3);
    let mut r = vec![vec![f64::INFINITY; n + 1]; m + 1];
    r[0][0] = 0.0;
    for i in 1..=m {
        for j in 1..=n {
            let cost = (a[i - 1] - b[j - 1]).powi(2);
            let r0 = -r[i - 1][j - 1] / g;
            let r1 = -r[i - 1][j] / g;
            let r2 = -r[i][j - 1] / g;
            let rmax = r0.max(r1).max(r2);
            let softmin = if rmax.is_finite() {
                -g * (rmax + ((r0 - rmax).exp() + (r1 - rmax).exp() + (r2 - rmax).exp()).ln())
            } else {
                0.0
            };
            r[i][j] = cost + softmin;
        }
    }
    r[m][n]
}

/// Shape-Based Distance (Paparrizos & Gravano 2015) for k-Shape.
pub fn sbd_distance(a: &[f64], b: &[f64]) -> (f64, isize) {
    let m = a.len();
    if m == 0 {
        return (0.0, 0);
    }
    let na = norm(a);
    let nb = norm(b);
    if na < 1e-12 || nb < 1e-12 {
        return (1.0, 0);
    }
    let mut best_ncc = f64::NEG_INFINITY;
    let mut best_shift = 0isize;
    let max_shift = (m as isize) - 1;
    for shift in -max_shift..=max_shift {
        let mut cc = 0.0;
        for i in 0..m {
            let j = i as isize - shift;
            if j >= 0 && (j as usize) < m {
                cc += a[i] * b[j as usize];
            }
        }
        let ncc = cc / (na * nb);
        if ncc > best_ncc {
            best_ncc = ncc;
            best_shift = shift;
        }
    }
    ((1.0 - best_ncc).clamp(0.0, 2.0), best_shift)
}

/// Sinkhorn-Knopp optimal transport normalization (used by SwAV).
pub fn sinkhorn_knopp(scores: &[Vec<f64>], eps: f64, iters: usize) -> Vec<Vec<f64>> {
    let n = scores.len();
    if n == 0 {
        return vec![];
    }
    let k = scores[0].len();
    if k == 0 {
        return vec![vec![]; n];
    }
    let max_s = scores
        .iter()
        .flat_map(|r| r.iter())
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let mut q: Vec<Vec<f64>> = scores
        .iter()
        .map(|row| {
            row.iter()
                .map(|&v| ((v - max_s) / eps.max(1e-3)).exp().max(1e-15))
                .collect()
        })
        .collect();

    for _ in 0..iters {
        // Normalize columns to sum to n/k
        for c in 0..k {
            let col_sum: f64 = (0..n).map(|i| q[i][c]).sum::<f64>().max(1e-15);
            let target = (n as f64) / (k as f64);
            for i in 0..n {
                q[i][c] *= target / col_sum;
            }
        }
        // Normalize rows to sum to 1
        for i in 0..n {
            let row_sum: f64 = q[i].iter().sum::<f64>().max(1e-15);
            for c in 0..k {
                q[i][c] /= row_sum;
            }
        }
    }
    q
}

/// Error function approximation (Abramowitz & Stegun 7.1.26) for standard normal CDF.
pub fn normal_cdf(x: f64) -> f64 {
    let t = 1.0 / (1.0 + 0.2316419 * x.abs() / std::f64::consts::SQRT_2);
    let poly = t
        * (0.254829592
            + t * (-0.284496736
                + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    let erf = 1.0 - poly * (-(x * x) / 2.0).exp();
    let signed_erf = if x >= 0.0 { erf } else { -erf };
    0.5 * (1.0 + signed_erf)
}

/// Anderson-Darling test statistic A*^2 for normality (used by G-Means).
pub fn anderson_darling_normal(vals: &[f64]) -> f64 {
    let n = vals.len();
    if n < 5 {
        return 0.0;
    }
    let mean = vals.iter().sum::<f64>() / n as f64;
    let var = vals.iter().map(|&v| (v - mean).powi(2)).sum::<f64>() / (n - 1) as f64;
    let std = var.sqrt();
    if std < 1e-12 {
        return 0.0;
    }
    let mut z: Vec<f64> = vals.iter().map(|&v| (v - mean) / std).collect();
    z.sort_by(|a, b| a.total_cmp(b));

    let nf = n as f64;
    let mut sum = 0.0;
    for i in 0..n {
        let cdf_i = normal_cdf(z[i]).clamp(1e-12, 1.0 - 1e-12);
        let cdf_rev = normal_cdf(z[n - 1 - i]).clamp(1e-12, 1.0 - 1e-12);
        sum += (2 * (i + 1) - 1) as f64 * (cdf_i.ln() + (1.0 - cdf_rev).ln());
    }
    let a2 = -nf - sum / nf;
    a2 * (1.0 + 4.0 / nf - 25.0 / (nf * nf))
}

/// Relabels cluster IDs cleanly to 0..k-1 ordered by descending cluster size, preserving -1 noise.
pub fn relabel_contiguous(labels: &[i32]) -> Vec<i32> {
    use std::collections::HashMap;
    let mut counts: HashMap<i32, usize> = HashMap::new();
    for &l in labels {
        if l >= 0 {
            *counts.entry(l).or_insert(0) += 1;
        }
    }
    let mut sorted: Vec<(i32, usize)> = counts.into_iter().collect();
    sorted.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));

    let mut map = HashMap::new();
    for (new_id, (old_id, _)) in sorted.into_iter().enumerate() {
        map.insert(old_id, new_id as i32);
    }

    labels
        .iter()
        .map(|&l| if l < 0 { -1 } else { *map.get(&l).unwrap_or(&-1) })
        .collect()
}

pub fn compute_centroids(x: &[Vec<f64>], labels: &[i32], k: usize) -> Vec<Vec<f64>> {
    if x.is_empty() || k == 0 {
        return vec![];
    }
    let d = x[0].len();
    let mut sums = vec![vec![0.0; d]; k];
    let mut counts = vec![0usize; k];
    for (i, &l) in labels.iter().enumerate() {
        if l >= 0 && (l as usize) < k {
            let c = l as usize;
            counts[c] += 1;
            for j in 0..d {
                sums[c][j] += x[i][j];
            }
        }
    }
    for c in 0..k {
        if counts[c] > 0 {
            for j in 0..d {
                sums[c][j] /= counts[c] as f64;
            }
        }
    }
    sums
}

pub fn median_of_slice(vals: &mut [f64]) -> f64 {
    if vals.is_empty() {
        return 0.0;
    }
    vals.sort_by(|a, b| a.total_cmp(b));
    let mid = vals.len() / 2;
    if vals.len() % 2 == 1 {
        vals[mid]
    } else {
        0.5 * (vals[mid - 1] + vals[mid])
    }
}
