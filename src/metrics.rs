use std::collections::HashMap;

use crate::data::PreparedData;
use crate::math::{compute_centroids, euclidean, euclidean_sq};

#[derive(Debug, Clone)]
pub struct ClusterMetrics {
    pub n_clusters: usize,
    pub n_noise: usize,
    pub noise_ratio: f64,
    pub silhouette: Option<f64>,
    pub davies_bouldin: Option<f64>,
    pub calinski_harabasz: Option<f64>,
    pub dunn_index: Option<f64>,
    pub inertia_sse: f64,
    pub balance_entropy: f64,
    pub adjusted_rand_index: Option<f64>,
    pub cluster_sizes: Vec<usize>,
}

pub fn evaluate_clustering(data: &PreparedData, labels: &[i32]) -> ClusterMetrics {
    let n = labels.len();
    let mut n_noise = 0usize;
    let mut counts_map: HashMap<i32, usize> = HashMap::new();
    for &l in labels {
        if l < 0 {
            n_noise += 1;
        } else {
            *counts_map.entry(l).or_insert(0) += 1;
        }
    }

    let n_clusters = counts_map.len();
    let mut cluster_sizes = vec![0usize; n_clusters];
    for (&cid, &cnt) in &counts_map {
        if (cid as usize) < n_clusters {
            cluster_sizes[cid as usize] = cnt;
        }
    }
    if cluster_sizes.iter().sum::<usize>() != n - n_noise {
        let mut vals: Vec<usize> = counts_map.values().copied().collect();
        vals.sort_by(|a, b| b.cmp(a));
        cluster_sizes = vals;
    }

    let noise_ratio = if n > 0 { n_noise as f64 / n as f64 } else { 0.0 };

    // Normalized cluster size balance entropy in [0, 1]
    let valid_n = n.saturating_sub(n_noise);
    let balance_entropy = if n_clusters >= 2 && valid_n > 0 {
        let h: f64 = cluster_sizes
            .iter()
            .filter(|&&c| c > 0)
            .map(|&c| {
                let p = c as f64 / valid_n as f64;
                -p * p.ln()
            })
            .sum();
        (h / (n_clusters as f64).ln()).clamp(0.0, 1.0)
    } else {
        0.0
    };

    let centroids = compute_centroids(&data.x, labels, n_clusters);
    let mut inertia_sse = 0.0;
    for (i, &l) in labels.iter().enumerate() {
        if l >= 0 && (l as usize) < centroids.len() {
            inertia_sse += euclidean_sq(&data.x[i], &centroids[l as usize]);
        }
    }

    let silhouette = silhouette_score(&data.dist_matrix, labels);
    let davies_bouldin = davies_bouldin_index(&data.x, labels, n_clusters, &centroids);
    let calinski_harabasz = calinski_harabasz_index(&data.x, labels, n_clusters, &centroids, inertia_sse);
    let dunn_index = dunn_index(&data.dist_matrix, labels, n_clusters);
    let adjusted_rand_index = data
        .ground_truth
        .as_ref()
        .map(|gt| adjusted_rand_index(gt, labels));

    ClusterMetrics {
        n_clusters,
        n_noise,
        noise_ratio,
        silhouette,
        davies_bouldin,
        calinski_harabasz,
        dunn_index,
        inertia_sse,
        balance_entropy,
        adjusted_rand_index,
        cluster_sizes,
    }
}

pub fn silhouette_score(dist: &[Vec<f64>], labels: &[i32]) -> Option<f64> {
    let mut members: HashMap<i32, Vec<usize>> = HashMap::new();
    for (i, &l) in labels.iter().enumerate() {
        if l >= 0 {
            members.entry(l).or_default().push(i);
        }
    }
    if members.len() < 2 {
        return None;
    }

    let mut total_s = 0.0;
    let mut valid_pts = 0usize;
    for (i, &li) in labels.iter().enumerate() {
        if li < 0 {
            continue;
        }
        let own = &members[&li];
        if own.len() <= 1 {
            valid_pts += 1;
            continue; // s(i) = 0 for singleton
        }
        let a = own
            .iter()
            .filter(|&&j| j != i)
            .map(|&j| dist[i][j])
            .sum::<f64>()
            / (own.len() - 1) as f64;

        let mut b = f64::INFINITY;
        for (&other_l, other_pts) in &members {
            if other_l == li || other_pts.is_empty() {
                continue;
            }
            let avg = other_pts.iter().map(|&j| dist[i][j]).sum::<f64>() / other_pts.len() as f64;
            if avg < b {
                b = avg;
            }
        }
        let denom = a.max(b);
        let s = if denom > 1e-12 { (b - a) / denom } else { 0.0 };
        total_s += s;
        valid_pts += 1;
    }

    if valid_pts > 0 {
        Some((total_s / valid_pts as f64).clamp(-1.0, 1.0))
    } else {
        None
    }
}

fn davies_bouldin_index(
    x: &[Vec<f64>],
    labels: &[i32],
    k: usize,
    centroids: &[Vec<f64>],
) -> Option<f64> {
    if k < 2 || centroids.len() < 2 {
        return None;
    }
    let mut scatter = vec![0.0; k];
    let mut counts = vec![0usize; k];
    for (i, &l) in labels.iter().enumerate() {
        if l >= 0 && (l as usize) < k {
            let c = l as usize;
            scatter[c] += euclidean(&x[i], &centroids[c]);
            counts[c] += 1;
        }
    }
    let active: Vec<usize> = (0..k).filter(|&c| counts[c] > 0).collect();
    if active.len() < 2 {
        return None;
    }
    for &c in &active {
        scatter[c] /= counts[c] as f64;
    }

    let mut sum = 0.0;
    for &i in &active {
        let mut max_rij = 0.0_f64;
        for &j in &active {
            if i == j {
                continue;
            }
            let m_ij = euclidean(&centroids[i], &centroids[j]);
            if m_ij > 1e-12 {
                let rij = (scatter[i] + scatter[j]) / m_ij;
                if rij > max_rij {
                    max_rij = rij;
                }
            }
        }
        sum += max_rij;
    }
    Some(sum / active.len() as f64)
}

fn calinski_harabasz_index(
    x: &[Vec<f64>],
    labels: &[i32],
    k: usize,
    centroids: &[Vec<f64>],
    wcss: f64,
) -> Option<f64> {
    if k < 2 || x.is_empty() {
        return None;
    }
    let d = x[0].len();
    let mut global_mean = vec![0.0; d];
    let mut n_valid = 0usize;
    let mut counts = vec![0usize; k];
    for (i, &l) in labels.iter().enumerate() {
        if l >= 0 && (l as usize) < k {
            counts[l as usize] += 1;
            n_valid += 1;
            for j in 0..d {
                global_mean[j] += x[i][j];
            }
        }
    }
    let active_k = counts.iter().filter(|&&c| c > 0).count();
    if active_k < 2 || n_valid <= active_k || wcss <= 1e-12 {
        return None;
    }
    for v in &mut global_mean {
        *v /= n_valid as f64;
    }

    let mut bcss = 0.0;
    for c in 0..k {
        if counts[c] > 0 {
            bcss += (counts[c] as f64) * euclidean_sq(&centroids[c], &global_mean);
        }
    }
    let ch = (bcss / (active_k - 1) as f64) / (wcss / (n_valid - active_k) as f64);
    if ch.is_finite() {
        Some(ch)
    } else {
        None
    }
}

fn dunn_index(dist: &[Vec<f64>], labels: &[i32], k: usize) -> Option<f64> {
    if k < 2 {
        return None;
    }
    let mut groups: Vec<Vec<usize>> = vec![Vec::new(); k];
    for (i, &l) in labels.iter().enumerate() {
        if l >= 0 && (l as usize) < k {
            groups[l as usize].push(i);
        }
    }
    let groups: Vec<Vec<usize>> = groups.into_iter().filter(|g| !g.is_empty()).collect();
    if groups.len() < 2 {
        return None;
    }

    let mut max_diam = 0.0_f64;
    for g in &groups {
        for i in 0..g.len() {
            for j in (i + 1)..g.len() {
                let d = dist[g[i]][g[j]];
                if d > max_diam {
                    max_diam = d;
                }
            }
        }
    }
    if max_diam <= 1e-12 {
        return None;
    }

    let mut min_inter = f64::INFINITY;
    for a in 0..groups.len() {
        for b in (a + 1)..groups.len() {
            for &i in &groups[a] {
                for &j in &groups[b] {
                    let d = dist[i][j];
                    if d < min_inter {
                        min_inter = d;
                    }
                }
            }
        }
    }
    if min_inter.is_finite() {
        Some(min_inter / max_diam)
    } else {
        None
    }
}

/// Hubert & Arabie Adjusted Rand Index (ARI) between two labelings.
pub fn adjusted_rand_index(a: &[i32], b: &[i32]) -> f64 {
    let n = a.len().min(b.len());
    if n < 2 {
        return 0.0;
    }
    let mut map_a: HashMap<i32, usize> = HashMap::new();
    let mut map_b: HashMap<i32, usize> = HashMap::new();
    let mut pairs: HashMap<(usize, usize), usize> = HashMap::new();
    for i in 0..n {
        let next_a = map_a.len();
        let ia = *map_a.entry(a[i]).or_insert(next_a);
        let next_b = map_b.len();
        let ib = *map_b.entry(b[i]).or_insert(next_b);
        *pairs.entry((ia, ib)).or_insert(0) += 1;
    }

    let mut row_sums = vec![0usize; map_a.len()];
    let mut col_sums = vec![0usize; map_b.len()];
    for (&(ia, ib), &cnt) in &pairs {
        row_sums[ia] += cnt;
        col_sums[ib] += cnt;
    }

    let comb2 = |x: usize| -> f64 {
        if x < 2 {
            0.0
        } else {
            (x as f64) * ((x - 1) as f64) / 2.0
        }
    };

    let sum_nij: f64 = pairs.values().map(|&c| comb2(c)).sum();
    let sum_a: f64 = row_sums.iter().map(|&c| comb2(c)).sum();
    let sum_b: f64 = col_sums.iter().map(|&c| comb2(c)).sum();
    let total = comb2(n);
    if total <= 0.0 {
        return 0.0;
    }
    let expected = (sum_a * sum_b) / total;
    let max_idx = 0.5 * (sum_a + sum_b);
    let denom = max_idx - expected;
    if denom.abs() < 1e-12 {
        if (sum_nij - expected).abs() < 1e-12 {
            1.0
        } else {
            0.0
        }
    } else {
        ((sum_nij - expected) / denom).clamp(-1.0, 1.0)
    }
}
