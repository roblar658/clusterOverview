use std::collections::HashMap;
use std::f64::consts::PI;

use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;

use crate::algorithms::cat1_centroid::{init_kmeans_pp, run_lloyd_from_centers};
use crate::algorithms::{AlgorithmOutput, RunConfig};
use crate::data::PreparedData;
use crate::math::{dot, inv_and_logdet, mahalanobis_sq, normalize_vec};

pub fn run(id: usize, data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    match id {
        81 => gmm_em(data, cfg),
        82 => bayesian_dpgmm(data, cfg),
        83 => lda_gibbs(data, cfg),
        84 => plsa_em(data, cfg),
        85 => cobweb(data, cfg),
        86 => classit(data, cfg),
        87 => autoclass(data, cfg),
        88 => student_t_mixture(data, cfg),
        89 => vmf_mixture(data, cfg),
        90 => latent_class_analysis(data, cfg),
        91 => hdp_topic(data, cfg),
        92 => nmf_clustering(data, cfg),
        93 => mclust(data, cfg),
        94 => mixmod(data, cfg),
        95 => poisson_mixture(data, cfg),
        96 => bernoulli_mixture(data, cfg),
        97 => dmm_short_text(data, cfg),
        98 => skew_normal_mixture(data, cfg),
        99 => igmm_crp(data, cfg),
        100 => correlated_topic_model(data, cfg),
        101 => biterm_topic_model(data, cfg),
        102 => snob_mml(data, cfg),
        _ => unreachable!(),
    }
}

#[derive(Clone, Copy)]
enum CovModel {
    Eii, // Spherical equal volume
    Vii, // Spherical varying volume
    Eei, // Diagonal equal
    Vvi, // Diagonal varying
    Eee, // Full shared covariance
    Vvv, // Full unconstrained covariance
}

impl CovModel {
    fn name(&self) -> &'static str {
        match self {
            CovModel::Eii => "EII",
            CovModel::Vii => "VII",
            CovModel::Eei => "EEI",
            CovModel::Vvi => "VVI",
            CovModel::Eee => "EEE",
            CovModel::Vvv => "VVV",
        }
    }
    fn n_params(&self, k: usize, d: usize) -> usize {
        let base = (k - 1) + k * d;
        let cov_p = match self {
            CovModel::Eii => 1,
            CovModel::Vii => k,
            CovModel::Eei => d,
            CovModel::Vvi => k * d,
            CovModel::Eee => d * (d + 1) / 2,
            CovModel::Vvv => k * d * (d + 1) / 2,
        };
        base + cov_p
    }
}

fn fit_gmm_model(
    x: &[Vec<f64>],
    k: usize,
    model: CovModel,
    max_iter: usize,
    seed: u64,
) -> (Vec<i32>, f64, f64) {
    let n = x.len();
    let d = x[0].len();
    let reg = 1e-3;
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    let c0 = init_kmeans_pp(x, k, &mut rng);
    let (init_lbl, mut means, _, _) = run_lloyd_from_centers(x, c0, 12);

    let mut weights = vec![1.0 / k as f64; k];
    let mut covs = vec![vec![vec![0.0; d]; d]; k];
    for c in 0..k {
        let mems: Vec<usize> = (0..n).filter(|&i| init_lbl[i] == c as i32).collect();
        weights[c] = (mems.len() as f64 / n as f64).max(1e-3);
        for &i in &mems {
            for a in 0..d {
                for b in 0..d {
                    covs[c][a][b] += (x[i][a] - means[c][a]) * (x[i][b] - means[c][b]);
                }
            }
        }
        let denom = mems.len().max(1) as f64;
        for a in 0..d {
            for b in 0..d {
                covs[c][a][b] = covs[c][a][b] / denom + if a == b { reg } else { 0.0 };
            }
        }
    }

    let mut resp = vec![vec![1.0 / k as f64; k]; n];
    let mut log_lik = f64::NEG_INFINITY;
    let const_term = -0.5 * (d as f64) * (2.0 * PI).ln();

    for _ in 0..max_iter {
        // Enforce covariance structure
        apply_cov_constraint(&mut covs, &weights, model, reg);

        let inv_and_ld: Vec<(Vec<Vec<f64>>, f64)> =
            covs.iter().map(|cv| inv_and_logdet(cv, reg)).collect();

        // E-step
        let mut new_ll = 0.0;
        for i in 0..n {
            let mut log_p = vec![0.0; k];
            let mut max_lp = f64::NEG_INFINITY;
            for c in 0..k {
                let m2 = mahalanobis_sq(&x[i], &means[c], &inv_and_ld[c].0);
                let lp = weights[c].max(1e-12).ln() + const_term - 0.5 * inv_and_ld[c].1 - 0.5 * m2;
                log_p[c] = lp;
                if lp > max_lp {
                    max_lp = lp;
                }
            }
            let sum_exp: f64 = log_p.iter().map(|&lp| (lp - max_lp).exp()).sum::<f64>().max(1e-15);
            let lse = max_lp + sum_exp.ln();
            new_ll += lse;
            for c in 0..k {
                resp[i][c] = (log_p[c] - lse).exp();
            }
        }

        // M-step
        for c in 0..k {
            let nk: f64 = resp.iter().map(|r| r[c]).sum::<f64>().max(1e-6);
            weights[c] = (nk / n as f64).max(1e-4);
            for j in 0..d {
                means[c][j] = (0..n).map(|i| resp[i][c] * x[i][j]).sum::<f64>() / nk;
            }
            for a in 0..d {
                for b in 0..d {
                    let mut s = 0.0;
                    for i in 0..n {
                        s += resp[i][c] * (x[i][a] - means[c][a]) * (x[i][b] - means[c][b]);
                    }
                    covs[c][a][b] = s / nk + if a == b { reg } else { 0.0 };
                }
            }
        }
        let wsum: f64 = weights.iter().sum();
        for w in &mut weights {
            *w /= wsum;
        }

        if (new_ll - log_lik).abs() < 1e-4 {
            log_lik = new_ll;
            break;
        }
        log_lik = new_ll;
    }

    let labels: Vec<i32> = resp
        .iter()
        .map(|r| {
            r.iter()
                .enumerate()
                .max_by(|a, b| a.1.total_cmp(b.1))
                .map(|(c, _)| c as i32)
                .unwrap_or(0)
        })
        .collect();
    let bic = 2.0 * log_lik - (model.n_params(k, d) as f64) * (n as f64).ln();
    (labels, log_lik, bic)
}

fn apply_cov_constraint(covs: &mut [Vec<Vec<f64>>], weights: &[f64], model: CovModel, reg: f64) {
    let k = covs.len();
    let d = covs[0].len();
    match model {
        CovModel::Vvv => {}
        CovModel::Eee => {
            let mut pooled = vec![vec![0.0; d]; d];
            for c in 0..k {
                for a in 0..d {
                    for b in 0..d {
                        pooled[a][b] += weights[c] * covs[c][a][b];
                    }
                }
            }
            for c in 0..k {
                covs[c] = pooled.clone();
            }
        }
        CovModel::Vvi => {
            for c in 0..k {
                for a in 0..d {
                    for b in 0..d {
                        if a != b {
                            covs[c][a][b] = 0.0;
                        }
                    }
                }
            }
        }
        CovModel::Eei => {
            let mut diag = vec![0.0; d];
            for c in 0..k {
                for a in 0..d {
                    diag[a] += weights[c] * covs[c][a][a];
                }
            }
            for c in 0..k {
                for a in 0..d {
                    for b in 0..d {
                        covs[c][a][b] = if a == b { diag[a].max(reg) } else { 0.0 };
                    }
                }
            }
        }
        CovModel::Vii => {
            for c in 0..k {
                let s = ((0..d).map(|a| covs[c][a][a]).sum::<f64>() / d as f64).max(reg);
                for a in 0..d {
                    for b in 0..d {
                        covs[c][a][b] = if a == b { s } else { 0.0 };
                    }
                }
            }
        }
        CovModel::Eii => {
            let mut s = 0.0;
            for c in 0..k {
                s += weights[c] * (0..d).map(|a| covs[c][a][a]).sum::<f64>() / d as f64;
            }
            s = s.max(reg);
            for c in 0..k {
                for a in 0..d {
                    for b in 0..d {
                        covs[c][a][b] = if a == b { s } else { 0.0 };
                    }
                }
            }
        }
    }
}

// 81. Gaussian Mixture Models (GMM)
fn gmm_em(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let k = cfg.effective_k(data.n_samples);
    let (labels, ll, bic) = fit_gmm_model(&data.x, k, CovModel::Vvv, cfg.max_iters(40), cfg.seed + 81);
    AlgorithmOutput::new(
        labels,
        format!("k={}, full kovarians", k),
        format!("EM konvergerte: logL={:.1}, BIC={:.1}", ll, bic),
    )
}

// 82. Bayesian GMM (DPGMM — Variational Dirichlet Process Stick-Breaking)
fn bayesian_dpgmm(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let t_cap = (cfg.effective_k(n) + 4).min(10).min(n / 3).max(3);
    let alpha_dp = 0.6_f64;
    let reg = 0.02;
    let max_iter = cfg.max_iters(35);

    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 82);
    let mut means = init_kmeans_pp(&data.x, t_cap, &mut rng);
    let mut vars = vec![vec![1.0_f64; d]; t_cap];
    let mut stick_weights = vec![1.0 / t_cap as f64; t_cap];
    let mut resp = vec![vec![1.0 / t_cap as f64; t_cap]; n];

    // Global prior mean
    let mut prior_mean = vec![0.0; d];
    for row in &data.x {
        for j in 0..d {
            prior_mean[j] += row[j] / n as f64;
        }
    }

    for _ in 0..max_iter {
        // Variational E-step
        for i in 0..n {
            let mut log_r = vec![0.0; t_cap];
            let mut mx = f64::NEG_INFINITY;
            for c in 0..t_cap {
                let mut quad = 0.0;
                let mut ldet = 0.0;
                for j in 0..d {
                    let v = vars[c][j].max(reg);
                    quad += (data.x[i][j] - means[c][j]).powi(2) / v;
                    ldet += v.ln();
                }
                let lp = stick_weights[c].max(1e-12).ln() - 0.5 * (ldet + quad);
                log_r[c] = lp;
                if lp > mx {
                    mx = lp;
                }
            }
            let denom: f64 = log_r.iter().map(|&v| (v - mx).exp()).sum::<f64>().max(1e-15);
            for c in 0..t_cap {
                resp[i][c] = ((log_r[c] - mx).exp() / denom).max(1e-15);
            }
        }

        // Variational M-step with Stick-Breaking GEM(alpha_dp) prior
        let mut nk = vec![0.0; t_cap];
        for c in 0..t_cap {
            nk[c] = resp.iter().map(|r| r[c]).sum();
        }
        let mut rem_stick = 1.0;
        for c in 0..t_cap {
            let tail: f64 = nk[(c + 1)..].iter().sum();
            let gamma1 = 1.0 + nk[c];
            let gamma2 = alpha_dp + tail;
            let v_c = gamma1 / (gamma1 + gamma2);
            stick_weights[c] = (v_c * rem_stick).max(1e-10);
            rem_stick *= 1.0 - v_c;

            // MAP shrinkage toward prior_mean
            let kappa0 = 0.5;
            for j in 0..d {
                let s: f64 = (0..n).map(|i| resp[i][c] * data.x[i][j]).sum();
                means[c][j] = (kappa0 * prior_mean[j] + s) / (kappa0 + nk[c]);
            }
            for j in 0..d {
                let s: f64 = (0..n)
                    .map(|i| resp[i][c] * (data.x[i][j] - means[c][j]).powi(2))
                    .sum();
                vars[c][j] = ((s + 0.5) / (nk[c] + 1.0)).max(reg);
            }
        }
    }

    let labels: Vec<i32> = resp
        .iter()
        .map(|r| {
            r.iter()
                .enumerate()
                .max_by(|a, b| a.1.total_cmp(b.1))
                .map(|(c, _)| c as i32)
                .unwrap_or(0)
        })
        .collect();
    let active_k = labels.iter().copied().collect::<std::collections::HashSet<_>>().len();
    AlgorithmOutput::new(
        labels,
        format!("T_maks={}, alpha_DP={:.1}", t_cap, alpha_dp),
        format!("Stick-breaking Dirichlet-prosess aktiverte {} av {} komponenter", active_k, t_cap),
    )
}

// 83. Latent Dirichlet Allocation (LDA — Collapsed Gibbs Sampling on count matrix)
fn lda_gibbs(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let v = data.n_features;
    let k = cfg.effective_k(n);
    let alpha = 0.5_f64;
    let beta = 0.1_f64;
    let n_iters = cfg.max_iters(35);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 83);

    // Expand each row's feature counts into token lists (doc, word)
    let mut doc_tokens: Vec<Vec<(usize, usize)>> = Vec::with_capacity(n); // (word_id, topic)
    let mut n_dk = vec![vec![0usize; k]; n];
    let mut n_kw = vec![vec![0usize; v]; k];
    let mut n_k = vec![0usize; k];

    for doc in 0..n {
        let mut toks = Vec::new();
        for w in 0..v {
            let cnt = data.counts[doc][w].min(12);
            for _ in 0..cnt {
                let t = rng.gen_range(0..k);
                toks.push((w, t));
                n_dk[doc][t] += 1;
                n_kw[t][w] += 1;
                n_k[t] += 1;
            }
        }
        if toks.is_empty() {
            let t = rng.gen_range(0..k);
            toks.push((0, t));
            n_dk[doc][t] += 1;
            n_kw[t][0] += 1;
            n_k[t] += 1;
        }
        doc_tokens.push(toks);
    }

    let v_beta = v as f64 * beta;
    for _ in 0..n_iters {
        for doc in 0..n {
            for idx in 0..doc_tokens[doc].len() {
                let (w, old_t) = doc_tokens[doc][idx];
                n_dk[doc][old_t] -= 1;
                n_kw[old_t][w] -= 1;
                n_k[old_t] -= 1;

                let mut probs = vec![0.0; k];
                let mut sum_p = 0.0;
                for t in 0..k {
                    let p = (n_dk[doc][t] as f64 + alpha)
                        * ((n_kw[t][w] as f64 + beta) / (n_k[t] as f64 + v_beta));
                    probs[t] = p;
                    sum_p += p;
                }
                let mut r = rng.gen::<f64>() * sum_p;
                let mut new_t = k - 1;
                for t in 0..k {
                    r -= probs[t];
                    if r <= 0.0 {
                        new_t = t;
                        break;
                    }
                }
                doc_tokens[doc][idx].1 = new_t;
                n_dk[doc][new_t] += 1;
                n_kw[new_t][w] += 1;
                n_k[new_t] += 1;
            }
        }
    }

    let labels: Vec<i32> = (0..n)
        .map(|doc| {
            (0..k)
                .max_by_key(|&t| n_dk[doc][t])
                .unwrap_or(0) as i32
        })
        .collect();
    AlgorithmOutput::new(
        labels,
        format!("emner={}, alpha={:.1}, beta={:.1}", k, alpha, beta),
        format!("Collapsed Gibbs-sampling ({} iterasjoner)", n_iters),
    )
}

// 84. pLSA (Hofmann 1999 — Probabilistic Latent Semantic Analysis via EM)
fn plsa_em(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let v = data.n_features;
    let k = cfg.effective_k(n);
    let max_iter = cfg.max_iters(40);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 84);

    // P(z|d) and P(w|z)
    let mut p_zd = vec![vec![0.0; k]; n];
    for d in 0..n {
        let mut s = 0.0;
        for z in 0..k {
            p_zd[d][z] = rng.gen_range(0.1..1.0);
            s += p_zd[d][z];
        }
        for z in 0..k {
            p_zd[d][z] /= s;
        }
    }
    let mut p_wz = vec![vec![0.0; v]; k];
    for z in 0..k {
        let mut s = 0.0;
        for w in 0..v {
            p_wz[z][w] = rng.gen_range(0.1..1.0);
            s += p_wz[z][w];
        }
        for w in 0..v {
            p_wz[z][w] /= s;
        }
    }

    for _ in 0..max_iter {
        let mut new_zd = vec![vec![1e-8; k]; n];
        let mut new_wz = vec![vec![1e-8; v]; k];
        for d in 0..n {
            for w in 0..v {
                let x_dw = data.x_nonneg[d][w];
                let denom: f64 = (0..k).map(|z| p_zd[d][z] * p_wz[z][w]).sum::<f64>().max(1e-12);
                for z in 0..k {
                    let post = (p_zd[d][z] * p_wz[z][w]) / denom;
                    let contrib = x_dw * post;
                    new_zd[d][z] += contrib;
                    new_wz[z][w] += contrib;
                }
            }
        }
        for d in 0..n {
            let s: f64 = new_zd[d].iter().sum::<f64>().max(1e-12);
            for z in 0..k {
                p_zd[d][z] = new_zd[d][z] / s;
            }
        }
        for z in 0..k {
            let s: f64 = new_wz[z].iter().sum::<f64>().max(1e-12);
            for w in 0..v {
                p_wz[z][w] = new_wz[z][w] / s;
            }
        }
    }

    let labels: Vec<i32> = (0..n)
        .map(|d| {
            (0..k)
                .max_by(|&a, &b| p_zd[d][a].total_cmp(&p_zd[d][b]))
                .unwrap_or(0) as i32
        })
        .collect();
    AlgorithmOutput::new(
        labels,
        format!("k={}, EM-iter={}", k, max_iter),
        "Sannsynlighetsbasert matrisefaktorisering P(w|z)P(z|d)".to_string(),
    )
}

// 85. COBWEB (Fisher 1987 — Incremental conceptual clustering maximizing Category Utility)
fn cobweb(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let c_dim = data.categorical[0].len();
    let max_k = (cfg.effective_k(n) + 2).min(10);

    // Compute category utility CU = (1/K) sum_k P(C_k) [ sum_i sum_v P(A_i=v|C_k)^2 - P(A_i=v)^2 ]
    let mut base_sq = 0.0;
    for j in 0..c_dim {
        let mut freq: HashMap<usize, usize> = HashMap::new();
        for i in 0..n {
            *freq.entry(data.categorical[i][j]).or_insert(0) += 1;
        }
        for &cnt in freq.values() {
            base_sq += (cnt as f64 / n as f64).powi(2);
        }
    }

    let eval_cu = |clusters: &[Vec<usize>]| -> f64 {
        let total: usize = clusters.iter().map(|c| c.len()).sum();
        if total == 0 || clusters.is_empty() {
            return 0.0;
        }
        let mut sum = 0.0;
        for cl in clusters {
            if cl.is_empty() {
                continue;
            }
            let p_c = cl.len() as f64 / total as f64;
            let mut sq = 0.0;
            for j in 0..c_dim {
                let mut freq: HashMap<usize, usize> = HashMap::new();
                for &i in cl {
                    *freq.entry(data.categorical[i][j]).or_insert(0) += 1;
                }
                for &cnt in freq.values() {
                    sq += (cnt as f64 / cl.len() as f64).powi(2);
                }
            }
            sum += p_c * (sq - base_sq);
        }
        sum / clusters.len() as f64
    };

    let mut clusters: Vec<Vec<usize>> = Vec::new();
    for i in 0..n {
        if clusters.is_empty() {
            clusters.push(vec![i]);
            continue;
        }
        // Try inserting into each existing concept OR creating a new concept
        let mut best_cu = f64::NEG_INFINITY;
        let mut best_choice: Option<usize> = None;

        for c in 0..clusters.len() {
            clusters[c].push(i);
            let cu = eval_cu(&clusters);
            clusters[c].pop();
            if cu > best_cu {
                best_cu = cu;
                best_choice = Some(c);
            }
        }
        if clusters.len() < max_k {
            clusters.push(vec![i]);
            let cu_new = eval_cu(&clusters);
            clusters.pop();
            if cu_new > best_cu {
                best_choice = None;
            }
        }
        match best_choice {
            Some(c) => clusters[c].push(i),
            None => clusters.push(vec![i]),
        }
    }

    let final_cu = eval_cu(&clusters);
    let mut labels = vec![0i32; n];
    for (cid, cl) in clusters.iter().enumerate() {
        for &p in cl {
            labels[p] = cid as i32;
        }
    }
    AlgorithmOutput::new(
        labels,
        format!("maks_k={}, Category Utility", max_k),
        format!("Inkrementell COBWEB: {} konsepter (CU={:.3})", clusters.len(), final_cu),
    )
}

// 86. CLASSIT (Gennari et al. 1989 — Continuous Gaussian extension of COBWEB with acuity)
fn classit(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let max_k = (cfg.effective_k(n) + 2).min(10);
    let acuity = 0.25_f64; // minimum std floor

    // Continuous CU = (1/K) * sum_k P(C_k) * sum_j (1/(2*sqrt(pi))) * (1/sigma_kj - 1/sigma_j)
    let mut global_inv_std = 0.0;
    for j in 0..d {
        let mu = data.x.iter().map(|r| r[j]).sum::<f64>() / n as f64;
        let var = data.x.iter().map(|r| (r[j] - mu).powi(2)).sum::<f64>() / n.max(1) as f64;
        global_inv_std += 1.0 / var.sqrt().max(acuity);
    }

    let eval_cu = |clusters: &[Vec<usize>]| -> f64 {
        let total: usize = clusters.iter().map(|c| c.len()).sum();
        if total == 0 || clusters.is_empty() {
            return 0.0;
        }
        let mut s = 0.0;
        for cl in clusters {
            if cl.is_empty() {
                continue;
            }
            let pc = cl.len() as f64 / total as f64;
            let mut inv_s = 0.0;
            for j in 0..d {
                let mu = cl.iter().map(|&i| data.x[i][j]).sum::<f64>() / cl.len() as f64;
                let var = cl.iter().map(|&i| (data.x[i][j] - mu).powi(2)).sum::<f64>() / cl.len() as f64;
                inv_s += 1.0 / var.sqrt().max(acuity);
            }
            s += pc * (inv_s - global_inv_std);
        }
        s / (2.0 * PI.sqrt() * clusters.len() as f64)
    };

    let mut clusters: Vec<Vec<usize>> = Vec::new();
    for i in 0..n {
        if clusters.is_empty() {
            clusters.push(vec![i]);
            continue;
        }
        let mut best_cu = f64::NEG_INFINITY;
        let mut best_c = None;
        for c in 0..clusters.len() {
            clusters[c].push(i);
            let cu = eval_cu(&clusters);
            clusters[c].pop();
            if cu > best_cu {
                best_cu = cu;
                best_c = Some(c);
            }
        }
        if clusters.len() < max_k {
            clusters.push(vec![i]);
            let cu_new = eval_cu(&clusters);
            clusters.pop();
            if cu_new > best_cu {
                best_c = None;
            }
        }
        match best_c {
            Some(c) => clusters[c].push(i),
            None => clusters.push(vec![i]),
        }
    }

    let final_cu = eval_cu(&clusters);
    let mut labels = vec![0i32; n];
    for (cid, cl) in clusters.iter().enumerate() {
        for &p in cl {
            labels[p] = cid as i32;
        }
    }
    AlgorithmOutput::new(
        labels,
        format!("acuity={:.2}, maks_k={}", acuity, max_k),
        format!("Kontinuerlig CLASSIT: {} klynger (CU={:.3})", clusters.len(), final_cu),
    )
}

// 87. AutoClass (Cheeseman & Stutz 1996 — Bayesian mixed Gaussian + Categorical Naive-Bayes mixture)
fn autoclass(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let c_dim = data.categorical[0].len();
    let max_k = (cfg.effective_k(n) + 2).min(8).max(2);
    let max_cat = data.categorical.iter().flat_map(|r| r.iter()).copied().max().unwrap_or(0) + 1;

    let mut best_labels = vec![0i32; n];
    let mut best_score = f64::NEG_INFINITY;
    let mut best_k = 2;

    for k in 2..=max_k {
        let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 87 + k as u64);
        let c0 = init_kmeans_pp(&data.x, k, &mut rng);
        let (lbl0, mut means, _, _) = run_lloyd_from_centers(&data.x, c0, 10);
        let mut vars = vec![vec![1.0; d]; k];
        let mut pi = vec![1.0 / k as f64; k];
        let mut cat_p = vec![vec![vec![1.0 / max_cat as f64; max_cat]; c_dim]; k];
        let mut resp = vec![vec![1.0 / k as f64; k]; n];

        for i in 0..n {
            for c in 0..k {
                resp[i][c] = if lbl0[i] == c as i32 { 0.85 } else { 0.15 / (k - 1).max(1) as f64 };
            }
        }

        let mut ll = f64::NEG_INFINITY;
        for _ in 0..25 {
            // M-step (MAP with weak conjugate priors)
            for c in 0..k {
                let nk: f64 = resp.iter().map(|r| r[c]).sum::<f64>().max(1e-5);
                pi[c] = (nk + 1.0) / (n as f64 + k as f64);
                for j in 0..d {
                    means[c][j] = (0..n).map(|i| resp[i][c] * data.x[i][j]).sum::<f64>() / nk;
                    let v = (0..n)
                        .map(|i| resp[i][c] * (data.x[i][j] - means[c][j]).powi(2))
                        .sum::<f64>()
                        / nk;
                    vars[c][j] = v.max(0.05);
                }
                for j in 0..c_dim {
                    let mut cnts = vec![0.5; max_cat]; // Dirichlet(0.5) prior
                    for i in 0..n {
                        cnts[data.categorical[i][j]] += resp[i][c];
                    }
                    let tot: f64 = cnts.iter().sum();
                    for v in 0..max_cat {
                        cat_p[c][j][v] = cnts[v] / tot;
                    }
                }
            }
            // E-step
            ll = 0.0;
            for i in 0..n {
                let mut lp = vec![0.0; k];
                let mut mx = f64::NEG_INFINITY;
                for c in 0..k {
                    let mut l = pi[c].ln();
                    for j in 0..d {
                        l += -0.5 * (2.0 * PI * vars[c][j]).ln()
                            - 0.5 * (data.x[i][j] - means[c][j]).powi(2) / vars[c][j];
                    }
                    for j in 0..c_dim {
                        l += 0.35 * cat_p[c][j][data.categorical[i][j]].max(1e-8).ln();
                    }
                    lp[c] = l;
                    mx = mx.max(l);
                }
                let sum_e: f64 = lp.iter().map(|&v| (v - mx).exp()).sum::<f64>().max(1e-15);
                let lse = mx + sum_e.ln();
                ll += lse;
                for c in 0..k {
                    resp[i][c] = (lp[c] - lse).exp();
                }
            }
        }
        let n_params = k * (2 * d + c_dim * (max_cat - 1) + 1);
        let bic_approx = ll - 0.5 * (n_params as f64) * (n as f64).ln();
        if bic_approx > best_score {
            best_score = bic_approx;
            best_k = k;
            best_labels = resp
                .iter()
                .map(|r| {
                    r.iter()
                        .enumerate()
                        .max_by(|a, b| a.1.total_cmp(b.1))
                        .map(|(c, _)| c as i32)
                        .unwrap_or(0)
                })
                .collect();
        }
    }

    AlgorithmOutput::new(
        best_labels,
        format!("Bayesiansk miks, søk k=2..{}", max_k),
        format!("Valgte k={} (marginal log-evidens={:.1})", best_k, best_score),
    )
}

// 88. Student's t-Mixture Models (Peel & McLachlan 2000 — heavy-tailed robust EM with gamma weights u_ic)
fn student_t_mixture(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let nu = 4.0_f64; // degrees of freedom (heavy tails!)
    let reg = 1e-2;
    let max_iter = cfg.max_iters(35);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 88);

    let c0 = init_kmeans_pp(&data.x, k, &mut rng);
    let (_, mut means, _, _) = run_lloyd_from_centers(&data.x, c0, 10);
    let mut covs = vec![vec![vec![0.0; d]; d]; k];
    for c in 0..k {
        for j in 0..d {
            covs[c][j][j] = 1.0;
        }
    }
    let mut pi = vec![1.0 / k as f64; k];
    let mut resp = vec![vec![1.0 / k as f64; k]; n];
    let mut u_weights = vec![vec![1.0; k]; n];

    for _ in 0..max_iter {
        let inv_ld: Vec<(Vec<Vec<f64>>, f64)> =
            covs.iter().map(|cv| inv_and_logdet(cv, reg)).collect();

        // E-step: Student-t posterior + robustness scale weight u_ic = (nu + d) / (nu + delta^2)
        for i in 0..n {
            let mut lp = vec![0.0; k];
            let mut mx = f64::NEG_INFINITY;
            for c in 0..k {
                let m2 = mahalanobis_sq(&data.x[i], &means[c], &inv_ld[c].0);
                u_weights[i][c] = (nu + d as f64) / (nu + m2);
                let l = pi[c].max(1e-12).ln()
                    - 0.5 * inv_ld[c].1
                    - 0.5 * (nu + d as f64) * (1.0 + m2 / nu).ln();
                lp[c] = l;
                mx = mx.max(l);
            }
            let s: f64 = lp.iter().map(|&v| (v - mx).exp()).sum::<f64>().max(1e-15);
            for c in 0..k {
                resp[i][c] = ((lp[c] - mx).exp() / s).max(1e-15);
            }
        }

        // M-step weighted by resp[i][c] * u_weights[i][c] (automatically downweights outliers!)
        for c in 0..k {
            let nk: f64 = resp.iter().map(|r| r[c]).sum::<f64>().max(1e-6);
            pi[c] = nk / n as f64;
            let wu_sum: f64 = (0..n).map(|i| resp[i][c] * u_weights[i][c]).sum::<f64>().max(1e-6);
            for j in 0..d {
                means[c][j] =
                    (0..n).map(|i| resp[i][c] * u_weights[i][c] * data.x[i][j]).sum::<f64>() / wu_sum;
            }
            for a in 0..d {
                for b in 0..d {
                    let mut s = 0.0;
                    for i in 0..n {
                        s += resp[i][c]
                            * u_weights[i][c]
                            * (data.x[i][a] - means[c][a])
                            * (data.x[i][b] - means[c][b]);
                    }
                    covs[c][a][b] = s / nk + if a == b { reg } else { 0.0 };
                }
            }
        }
    }

    let labels: Vec<i32> = resp
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
        format!("k={}, nu={:.1} (t-fordeling)", k, nu),
        "Tunghalet EM med automatisk nedvekting av uteliggere".to_string(),
    )
}

// 89. von Mises-Fisher Mixture Models (Banerjee et al. 2005 — spherical directional EM)
fn vmf_mixture(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let max_iter = cfg.max_iters(40);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 89);

    let mut mus: Vec<Vec<f64>> = init_kmeans_pp(&data.x_unit, k, &mut rng)
        .into_iter()
        .map(|v| normalize_vec(&v))
        .collect();
    let mut kappas = vec![4.0_f64; k];
    let mut pi = vec![1.0 / k as f64; k];
    let mut resp = vec![vec![1.0 / k as f64; k]; n];

    for _ in 0..max_iter {
        // E-step on unit hypersphere
        for i in 0..n {
            let mut lp = vec![0.0; k];
            let mut mx = f64::NEG_INFINITY;
            for c in 0..k {
                // Asymptotic log C_d(kappa) approximation for stability across dimensions
                let kap = kappas[c].clamp(0.1, 40.0);
                let log_cd = (d as f64 * 0.5 - 1.0) * kap.ln() - kap;
                let l = pi[c].max(1e-12).ln() + log_cd + kap * dot(&data.x_unit[i], &mus[c]);
                lp[c] = l;
                mx = mx.max(l);
            }
            let s: f64 = lp.iter().map(|&v| (v - mx).exp()).sum::<f64>().max(1e-15);
            for c in 0..k {
                resp[i][c] = ((lp[c] - mx).exp() / s).max(1e-15);
            }
        }
        // M-step with Banerjee et al. kappa approximation: kappa = r_bar*(d - r_bar^2)/(1 - r_bar^2)
        for c in 0..k {
            let nk: f64 = resp.iter().map(|r| r[c]).sum::<f64>().max(1e-6);
            pi[c] = nk / n as f64;
            let mut r_vec = vec![0.0; d];
            for i in 0..n {
                for j in 0..d {
                    r_vec[j] += resp[i][c] * data.x_unit[i][j];
                }
            }
            let r_norm = dot(&r_vec, &r_vec).sqrt();
            if r_norm > 1e-10 {
                mus[c] = r_vec.iter().map(|&v| v / r_norm).collect();
                let r_bar = (r_norm / nk).clamp(0.05, 0.96);
                kappas[c] = (r_bar * (d as f64 - r_bar * r_bar) / (1.0 - r_bar * r_bar)).clamp(0.5, 35.0);
            }
        }
    }

    let avg_kappa = kappas.iter().sum::<f64>() / k as f64;
    let labels: Vec<i32> = resp
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
        format!("k={}, snitt kappa={:.2}", k, avg_kappa),
        "Retningsfordeling på enhetshypersfære (vMF EM)".to_string(),
    )
}

// 90. Latent Class Analysis (LCA — categorical latent-class EM)
fn latent_class_analysis(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let c_dim = data.categorical[0].len();
    let k = cfg.effective_k(n);
    let max_cat = data.categorical.iter().flat_map(|r| r.iter()).copied().max().unwrap_or(0) + 1;
    let max_iter = cfg.max_iters(45);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 90);

    let mut pi = vec![1.0 / k as f64; k];
    let mut theta = vec![vec![vec![0.0_f64; max_cat]; c_dim]; k];
    for c in 0..k {
        for j in 0..c_dim {
            let mut s = 0.0;
            for v in 0..max_cat {
                theta[c][j][v] = rng.gen_range(0.1..1.0);
                s += theta[c][j][v];
            }
            for v in 0..max_cat {
                theta[c][j][v] /= s;
            }
        }
    }

    let mut resp = vec![vec![1.0 / k as f64; k]; n];
    for _ in 0..max_iter {
        for i in 0..n {
            let mut lp = vec![0.0; k];
            let mut mx = f64::NEG_INFINITY;
            for c in 0..k {
                let mut l = pi[c].max(1e-12).ln();
                for j in 0..c_dim {
                    l += theta[c][j][data.categorical[i][j]].max(1e-8).ln();
                }
                lp[c] = l;
                mx = mx.max(l);
            }
            let s: f64 = lp.iter().map(|&v| (v - mx).exp()).sum::<f64>().max(1e-15);
            for c in 0..k {
                resp[i][c] = (lp[c] - mx).exp() / s;
            }
        }
        for c in 0..k {
            let nk: f64 = resp.iter().map(|r| r[c]).sum::<f64>().max(1e-6);
            pi[c] = nk / n as f64;
            for j in 0..c_dim {
                let mut cnt = vec![0.2; max_cat];
                for i in 0..n {
                    cnt[data.categorical[i][j]] += resp[i][c];
                }
                let tot: f64 = cnt.iter().sum();
                for v in 0..max_cat {
                    theta[c][j][v] = cnt[v] / tot;
                }
            }
        }
    }

    let labels: Vec<i32> = resp
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
        format!("k={}, {} kategoriske var.", k, c_dim),
        "Latent Class Analysis (EM på betinget uavhengige kategorier)".to_string(),
    )
}

// 91. Hierarchical Dirichlet Process (HDP — two-level stick-breaking topic model)
fn hdp_topic(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let v = data.n_features;
    let t_cap = (cfg.effective_k(n) + 4).min(10).max(3);
    let gamma_top = 1.0_f64;
    let alpha_local = 0.8_f64;
    let eta = 0.15_f64;
    let max_iter = cfg.max_iters(30);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 91);

    // Global topic weights beta (top-level stick breaking) & word-topic phi
    let mut beta_top = vec![1.0 / t_cap as f64; t_cap];
    let mut phi = vec![vec![0.0; v]; t_cap];
    for t in 0..t_cap {
        let mut s = 0.0;
        for w in 0..v {
            phi[t][w] = rng.gen_range(0.1..1.0);
            s += phi[t][w];
        }
        for w in 0..v {
            phi[t][w] /= s;
        }
    }
    let mut theta = vec![vec![1.0 / t_cap as f64; t_cap]; n];

    for _ in 0..max_iter {
        let mut topic_word = vec![vec![eta; v]; t_cap];
        let mut topic_total = vec![0.0; t_cap];

        for d in 0..n {
            let mut doc_topic = vec![0.0; t_cap];
            for w in 0..v {
                let xdw = data.x_nonneg[d][w];
                let mut probs = vec![0.0; t_cap];
                let mut sp = 0.0;
                for t in 0..t_cap {
                    let p = theta[d][t] * phi[t][w];
                    probs[t] = p;
                    sp += p;
                }
                let sp = sp.max(1e-12);
                for t in 0..t_cap {
                    let c = xdw * (probs[t] / sp);
                    doc_topic[t] += c;
                    topic_word[t][w] += c;
                    topic_total[t] += c;
                }
            }
            // Local Dirichlet posterior centered on alpha_local * beta_top
            let denom: f64 = doc_topic.iter().sum::<f64>() + alpha_local;
            for t in 0..t_cap {
                theta[d][t] = (doc_topic[t] + alpha_local * beta_top[t]) / denom.max(1e-12);
            }
        }
        // Top-level stick-breaking update for beta_top
        let mut rem = 1.0;
        for t in 0..t_cap {
            let tail: f64 = topic_total[(t + 1)..].iter().sum();
            let v_t = (1.0 + topic_total[t]) / (1.0 + gamma_top + topic_total[t] + tail);
            beta_top[t] = (v_t * rem).max(1e-6);
            rem *= 1.0 - v_t;
            let sw: f64 = topic_word[t].iter().sum::<f64>().max(1e-12);
            for w in 0..v {
                phi[t][w] = topic_word[t][w] / sw;
            }
        }
    }

    let labels: Vec<i32> = (0..n)
        .map(|d| {
            (0..t_cap)
                .max_by(|&a, &b| theta[d][a].total_cmp(&theta[d][b]))
                .unwrap_or(0) as i32
        })
        .collect();
    let active = labels.iter().copied().collect::<std::collections::HashSet<_>>().len();
    AlgorithmOutput::new(
        labels,
        format!("T_maks={}, gamma={:.1}, alpha={:.1}", t_cap, gamma_top, alpha_local),
        format!("Hierarkisk Dirichlet-prosess utledet {} aktive emner", active),
    )
}

// 92. NMF (Lee & Seung 1999 — Multiplicative Updates V ≈ W * H)
fn nmf_clustering(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let max_iter = cfg.max_iters(60);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 92);

    let mut w = vec![vec![0.0; k]; n];
    let mut h = vec![vec![0.0; d]; k];
    for i in 0..n {
        for c in 0..k {
            w[i][c] = rng.gen_range(0.05..1.0);
        }
    }
    for c in 0..k {
        for j in 0..d {
            h[c][j] = rng.gen_range(0.05..1.0);
        }
    }

    for _ in 0..max_iter {
        // Update H: H <- H * (W^T V) / (W^T W H)
        let mut wt_v = vec![vec![0.0; d]; k];
        let mut wt_w = vec![vec![0.0; k]; k];
        for i in 0..n {
            for a in 0..k {
                for j in 0..d {
                    wt_v[a][j] += w[i][a] * data.x_nonneg[i][j];
                }
                for b in 0..k {
                    wt_w[a][b] += w[i][a] * w[i][b];
                }
            }
        }
        for a in 0..k {
            for j in 0..d {
                let denom: f64 = (0..k).map(|b| wt_w[a][b] * h[b][j]).sum::<f64>().max(1e-10);
                h[a][j] *= wt_v[a][j] / denom;
            }
        }

        // Update W: W <- W * (V H^T) / (W H H^T)
        let mut h_ht = vec![vec![0.0; k]; k];
        for a in 0..k {
            for b in 0..k {
                for j in 0..d {
                    h_ht[a][b] += h[a][j] * h[b][j];
                }
            }
        }
        for i in 0..n {
            for a in 0..k {
                let num: f64 = (0..d).map(|j| data.x_nonneg[i][j] * h[a][j]).sum();
                let den: f64 = (0..k).map(|b| w[i][b] * h_ht[b][a]).sum::<f64>().max(1e-10);
                w[i][a] *= num / den;
            }
        }
    }

    let labels: Vec<i32> = (0..n)
        .map(|i| {
            (0..k)
                .max_by(|&a, &b| w[i][a].total_cmp(&w[i][b]))
                .unwrap_or(0) as i32
        })
        .collect();
    AlgorithmOutput::new(
        labels,
        format!("k={}, multiplikativ NMF", k),
        format!("Faktoriserte V ≈ W*H over {} iterasjoner", max_iter),
    )
}

// 93. MCLUST (Fraley & Raftery 2002 — BIC selection across EII, VII, EEI, VVI, EEE, VVV models)
fn mclust(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let k_lo = 2usize;
    let k_hi = (cfg.effective_k(n) + 2).min(7).min(n / 4).max(2);
    let models = [
        CovModel::Eii,
        CovModel::Vii,
        CovModel::Eei,
        CovModel::Vvi,
        CovModel::Eee,
        CovModel::Vvv,
    ];

    let mut best_bic = f64::NEG_INFINITY;
    let mut best_labels = vec![0i32; n];
    let mut best_name = "VVV";
    let mut best_k = k_lo;

    for k in k_lo..=k_hi {
        for &m in &models {
            let (lbl, _, bic) = fit_gmm_model(&data.x, k, m, 22, cfg.seed + 93 + k as u64);
            if bic > best_bic {
                best_bic = bic;
                best_labels = lbl;
                best_name = m.name();
                best_k = k;
            }
        }
    }
    AlgorithmOutput::new(
        best_labels,
        format!("valgt={} (k={}), søk 6 modeller", best_name, best_k),
        format!("BIC valgte kovariansstruktur {} med k={} (BIC={:.1})", best_name, best_k, best_bic),
    )
}

// 94. Mixmod (Biernacki et al. 2006 — Parsimonious Gaussian + Multinomial model selection via ICL)
fn mixmod(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let k = cfg.effective_k(n);
    let models = [CovModel::Eii, CovModel::Vvi, CovModel::Eee, CovModel::Vvv];
    let mut best_icl = f64::NEG_INFINITY;
    let mut best_labels = vec![0i32; n];
    let mut best_m = "VVI";

    for &m in &models {
        let (lbl, _, bic) = fit_gmm_model(&data.x, k, m, 25, cfg.seed + 94);
        // ICL penalizes BIC by cluster overlap entropy proxy
        let mut counts = vec![0usize; k];
        for &l in &lbl {
            if l >= 0 && (l as usize) < k {
                counts[l as usize] += 1;
            }
        }
        let entropy: f64 = counts
            .iter()
            .filter(|&&c| c > 0)
            .map(|&c| {
                let p = c as f64 / n as f64;
                -p * p.ln()
            })
            .sum();
        let icl = bic - entropy;
        if icl > best_icl {
            best_icl = icl;
            best_labels = lbl;
            best_m = m.name();
        }
    }
    AlgorithmOutput::new(
        best_labels,
        format!("k={}, Mixmod ({})", k, best_m),
        format!("ICL-valgt parsimonisk modell {} (ICL={:.1})", best_m, best_icl),
    )
}

// 95. Poisson Mixture Models (EM for discrete count rates lambda_cj)
fn poisson_mixture(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let max_iter = cfg.max_iters(40);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 95);

    let mut lambdas = vec![vec![1.0; d]; k];
    for c in 0..k {
        let seed_i = rng.gen_range(0..n);
        for j in 0..d {
            lambdas[c][j] = (data.counts[seed_i][j] as f64 + 0.5).max(0.2);
        }
    }
    let mut pi = vec![1.0 / k as f64; k];
    let mut resp = vec![vec![1.0 / k as f64; k]; n];

    for _ in 0..max_iter {
        for i in 0..n {
            let mut lp = vec![0.0; k];
            let mut mx = f64::NEG_INFINITY;
            for c in 0..k {
                let mut l = pi[c].max(1e-12).ln();
                for j in 0..d {
                    let x = data.counts[i][j] as f64;
                    let lam = lambdas[c][j].max(1e-4);
                    // log Poisson up to constant -ln(x!)
                    l += x * lam.ln() - lam;
                }
                lp[c] = l;
                mx = mx.max(l);
            }
            let s: f64 = lp.iter().map(|&v| (v - mx).exp()).sum::<f64>().max(1e-15);
            for c in 0..k {
                resp[i][c] = (lp[c] - mx).exp() / s;
            }
        }
        for c in 0..k {
            let nk: f64 = resp.iter().map(|r| r[c]).sum::<f64>().max(1e-6);
            pi[c] = nk / n as f64;
            for j in 0..d {
                let s: f64 = (0..n).map(|i| resp[i][c] * data.counts[i][j] as f64).sum();
                lambdas[c][j] = ((s + 0.1) / (nk + 0.1)).max(0.05);
            }
        }
    }

    let labels: Vec<i32> = resp
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
        format!("k={}, Poisson-rater", k),
        "EM på diskrete telledata (Poisson-blanding)".to_string(),
    )
}

// 96. Bernoulli Mixture Models (EM on binary feature representation)
fn bernoulli_mixture(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.x_binary[0].len();
    let k = cfg.effective_k(n);
    let max_iter = cfg.max_iters(40);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 96);

    let mut mu = vec![vec![0.5_f64; d]; k];
    for c in 0..k {
        for j in 0..d {
            mu[c][j] = rng.gen_range(0.2..0.8);
        }
    }
    let mut pi = vec![1.0 / k as f64; k];
    let mut resp = vec![vec![1.0 / k as f64; k]; n];

    for _ in 0..max_iter {
        for i in 0..n {
            let mut lp = vec![0.0; k];
            let mut mx = f64::NEG_INFINITY;
            for c in 0..k {
                let mut l = pi[c].max(1e-12).ln();
                for j in 0..d {
                    let p = mu[c][j].clamp(1e-4, 1.0 - 1e-4);
                    l += if data.x_binary[i][j] == 1 { p.ln() } else { (1.0 - p).ln() };
                }
                lp[c] = l;
                mx = mx.max(l);
            }
            let s: f64 = lp.iter().map(|&v| (v - mx).exp()).sum::<f64>().max(1e-15);
            for c in 0..k {
                resp[i][c] = (lp[c] - mx).exp() / s;
            }
        }
        for c in 0..k {
            let nk: f64 = resp.iter().map(|r| r[c]).sum::<f64>().max(1e-6);
            pi[c] = nk / n as f64;
            for j in 0..d {
                let s: f64 = (0..n)
                    .map(|i| resp[i][c] * data.x_binary[i][j] as f64)
                    .sum();
                mu[c][j] = ((s + 0.5) / (nk + 1.0)).clamp(1e-4, 1.0 - 1e-4);
            }
        }
    }

    let labels: Vec<i32> = resp
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
        format!("k={}, binær Bernoulli EM", k),
        "Blandingsmodell for binære attributtforekomster".to_string(),
    )
}

// 97. Dirichlet Multinomial Mixture (Yin & Wang KDD 2014 — GSDMM Gibbs Sampling)
fn dmm_short_text(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let v = data.n_features;
    let k = cfg.effective_k(n);
    let alpha = 0.1_f64;
    let beta = 0.1_f64;
    let max_iter = cfg.max_iters(35);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 97);

    let doc_len: Vec<usize> = data.counts.iter().map(|r| r.iter().sum()).collect();
    let mut z: Vec<usize> = (0..n).map(|_| rng.gen_range(0..k)).collect();
    let mut m_z = vec![0usize; k];
    let mut n_z = vec![0usize; k];
    let mut n_zw = vec![vec![0usize; v]; k];

    for d in 0..n {
        let c = z[d];
        m_z[c] += 1;
        n_z[c] += doc_len[d];
        for w in 0..v {
            n_zw[c][w] += data.counts[d][w];
        }
    }

    let v_beta = v as f64 * beta;
    for _ in 0..max_iter {
        for d in 0..n {
            let old = z[d];
            m_z[old] -= 1;
            n_z[old] -= doc_len[d];
            for w in 0..v {
                n_zw[old][w] -= data.counts[d][w];
            }

            let mut log_p = vec![0.0; k];
            let mut mx = f64::NEG_INFINITY;
            for c in 0..k {
                let mut lp = (m_z[c] as f64 + alpha).ln();
                // Product over word occurrences in document d
                let mut denom_acc = 0.0;
                for j in 0..doc_len[d] {
                    denom_acc += (n_z[c] as f64 + v_beta + j as f64).ln();
                }
                let mut num_acc = 0.0;
                for w in 0..v {
                    for j in 0..data.counts[d][w] {
                        num_acc += (n_zw[c][w] as f64 + beta + j as f64).ln();
                    }
                }
                lp += num_acc - denom_acc;
                log_p[c] = lp;
                mx = mx.max(lp);
            }
            let probs: Vec<f64> = log_p.iter().map(|&v| (v - mx).exp()).collect();
            let sp: f64 = probs.iter().sum::<f64>().max(1e-12);
            let mut r = rng.gen::<f64>() * sp;
            let mut new_c = k - 1;
            for c in 0..k {
                r -= probs[c];
                if r <= 0.0 {
                    new_c = c;
                    break;
                }
            }
            z[d] = new_c;
            m_z[new_c] += 1;
            n_z[new_c] += doc_len[d];
            for w in 0..v {
                n_zw[new_c][w] += data.counts[d][w];
            }
        }
    }

    AlgorithmOutput::new(
        z.into_iter().map(|c| c as i32).collect(),
        format!("k={}, GSDMM a={:.1}, b={:.1}", k, alpha, beta),
        "Dirichlet Multinomial Mixture (1-klynge-per-observasjon Gibbs)".to_string(),
    )
}

// 98. Skew-Normal Mixture Models (Lin 2009 — EM with skewness shape parameter delta_c)
fn skew_normal_mixture(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let k = cfg.effective_k(n);
    let max_iter = cfg.max_iters(35);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 98);

    let c0 = init_kmeans_pp(&data.x, k, &mut rng);
    let (_, mut xi, _, _) = run_lloyd_from_centers(&data.x, c0, 10);
    let mut omega = vec![vec![1.0_f64; d]; k]; // diagonal scale
    let mut delta = vec![vec![0.25; d]; k]; // skewness vector
    let mut pi = vec![1.0 / k as f64; k];
    let mut resp = vec![vec![1.0 / k as f64; k]; n];

    for _ in 0..max_iter {
        // E-step with Azzalini skew-normal density 2 * phi(z) * Phi(alpha * z)
        let mut tau = vec![vec![0.0; k]; n];
        for i in 0..n {
            let mut lp = vec![0.0; k];
            let mut mx = f64::NEG_INFINITY;
            for c in 0..k {
                let mut l = pi[c].max(1e-12).ln();
                let mut proj = 0.0;
                for j in 0..d {
                    let s2 = omega[c][j].max(0.05);
                    let diff = data.x[i][j] - xi[c][j];
                    l += -0.5 * (2.0 * PI * s2).ln() - 0.5 * diff * diff / s2;
                    proj += delta[c][j] * diff / s2.sqrt();
                }
                let cdf = crate::math::normal_cdf(proj).max(1e-8);
                let pdf = (-0.5 * proj * proj).exp() / (2.0 * PI).sqrt();
                tau[i][c] = (proj + pdf / cdf).clamp(-3.0, 3.0);
                l += (2.0 * cdf).ln();
                lp[c] = l;
                mx = mx.max(l);
            }
            let s: f64 = lp.iter().map(|&v| (v - mx).exp()).sum::<f64>().max(1e-15);
            for c in 0..k {
                resp[i][c] = (lp[c] - mx).exp() / s;
            }
        }
        // M-step updating location xi, skewness delta, and scale omega
        for c in 0..k {
            let nk: f64 = resp.iter().map(|r| r[c]).sum::<f64>().max(1e-6);
            pi[c] = nk / n as f64;
            for j in 0..d {
                let num_xi: f64 = (0..n)
                    .map(|i| resp[i][c] * (data.x[i][j] - delta[c][j] * tau[i][c]))
                    .sum();
                xi[c][j] = num_xi / nk;
                let num_del: f64 = (0..n)
                    .map(|i| resp[i][c] * tau[i][c] * (data.x[i][j] - xi[c][j]))
                    .sum();
                let den_del: f64 = (0..n)
                    .map(|i| resp[i][c] * (tau[i][c] * tau[i][c] + 1.0))
                    .sum::<f64>()
                    .max(1e-6);
                delta[c][j] = (num_del / den_del).clamp(-1.5, 1.5);
                let v: f64 = (0..n)
                    .map(|i| {
                        resp[i][c] * (data.x[i][j] - xi[c][j] - delta[c][j] * tau[i][c]).powi(2)
                    })
                    .sum::<f64>()
                    / nk;
                omega[c][j] = v.max(0.05);
            }
        }
    }

    let labels: Vec<i32> = resp
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
        format!("k={}, skjev-normal EM", k),
        "Asymmetrisk Azzalini skjev-normalfordelingsblanding".to_string(),
    )
}

// 99. Infinite Gaussian Mixture Model (Rasmussen 2000 — Chinese Restaurant Process Collapsed Gibbs)
fn igmm_crp(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let alpha = 0.8_f64;
    let max_iter = cfg.max_iters(25);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 99);

    let mut global_mean = vec![0.0; d];
    for row in &data.x {
        for j in 0..d {
            global_mean[j] += row[j] / n as f64;
        }
    }
    let sigma2_prior = 1.5_f64;
    let sigma2_like = 0.45_f64;

    // Start with 2 tables
    let mut assign: Vec<usize> = (0..n).map(|i| i % 2).collect();
    let mut tables: Vec<Vec<usize>> = vec![Vec::new(); 2];
    for (i, &t) in assign.iter().enumerate() {
        tables[t].push(i);
    }

    for _ in 0..max_iter {
        for i in 0..n {
            let cur_t = assign[i];
            if let Some(pos) = tables[cur_t].iter().position(|&x| x == i) {
                tables[cur_t].swap_remove(pos);
            }

            // Compute CRP + conjugate Gaussian predictive log-likelihood for existing tables + 1 new table
            let mut log_p = Vec::with_capacity(tables.len() + 1);
            let mut active_t = Vec::new();
            for (t, mems) in tables.iter().enumerate() {
                if mems.is_empty() {
                    continue;
                }
                let nt = mems.len() as f64;
                let post_var = 1.0 / (1.0 / sigma2_prior + nt / sigma2_like);
                let pred_var = sigma2_like + post_var;
                let mut quad = 0.0;
                for j in 0..d {
                    let s: f64 = mems.iter().map(|&m| data.x[m][j]).sum();
                    let post_mu = post_var * (global_mean[j] / sigma2_prior + s / sigma2_like);
                    quad += (data.x[i][j] - post_mu).powi(2) / pred_var;
                }
                let lp = nt.ln() - 0.5 * (d as f64) * (2.0 * PI * pred_var).ln() - 0.5 * quad;
                log_p.push(lp);
                active_t.push(Some(t));
            }
            // New table probability (alpha * prior predictive)
            if active_t.len() < 10 {
                let pred_var0 = sigma2_like + sigma2_prior;
                let quad0: f64 = (0..d)
                    .map(|j| (data.x[i][j] - global_mean[j]).powi(2) / pred_var0)
                    .sum();
                let lp_new = alpha.ln() - 0.5 * (d as f64) * (2.0 * PI * pred_var0).ln() - 0.5 * quad0;
                log_p.push(lp_new);
                active_t.push(None);
            }

            let mx = log_p.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            let probs: Vec<f64> = log_p.iter().map(|&v| (v - mx).exp()).collect();
            let sp: f64 = probs.iter().sum::<f64>().max(1e-12);
            let mut r = rng.gen::<f64>() * sp;
            let mut pick = active_t.len() - 1;
            for (idx, &p) in probs.iter().enumerate() {
                r -= p;
                if r <= 0.0 {
                    pick = idx;
                    break;
                }
            }

            match active_t[pick] {
                Some(t) => {
                    tables[t].push(i);
                    assign[i] = t;
                }
                None => {
                    let new_id = tables.len();
                    tables.push(vec![i]);
                    assign[i] = new_id;
                }
            }
        }
    }

    let n_tables = tables.iter().filter(|t| !t.is_empty()).count();
    AlgorithmOutput::new(
        assign.into_iter().map(|t| t as i32).collect(),
        format!("CRP alpha={:.1}, Gibbs={}", alpha, max_iter),
        format!("Chinese Restaurant Process konvergerte til {} bord/klynger", n_tables),
    )
}

// 100. Correlated Topic Model (Blei & Lafferty 2007 — Logistic-Normal correlated topic proportions)
fn correlated_topic_model(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let v = data.n_features;
    let k = cfg.effective_k(n);
    let max_iter = cfg.max_iters(30);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 100);

    // Latent Gaussian log-odds eta[d][k] -> softmax -> topic proportions theta[d][k]
    let mut eta = vec![vec![0.0; k]; n];
    for d in 0..n {
        for c in 0..k {
            eta[d][c] = rng.gen_range(-0.5..0.5);
        }
    }
    let mut mu_prior = vec![0.0; k];
    let mut cov_prior = vec![vec![0.0; k]; k];
    for c in 0..k {
        cov_prior[c][c] = 1.0;
    }

    let mut beta_kw = vec![vec![1.0 / v as f64; v]; k];
    for c in 0..k {
        let mut s = 0.0;
        for w in 0..v {
            beta_kw[c][w] = rng.gen_range(0.1..1.0);
            s += beta_kw[c][w];
        }
        for w in 0..v {
            beta_kw[c][w] /= s;
        }
    }

    let softmax = |z: &[f64]| -> Vec<f64> {
        let mx = z.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let e: Vec<f64> = z.iter().map(|&v| (v - mx).exp()).collect();
        let s: f64 = e.iter().sum::<f64>().max(1e-12);
        e.into_iter().map(|v| v / s).collect()
    };

    for _ in 0..max_iter {
        let (inv_cov, _) = inv_and_logdet(&cov_prior, 0.1);
        let mut new_beta = vec![vec![1e-4; v]; k];

        for d in 0..n {
            let theta = softmax(&eta[d]);
            let mut phi_sum = vec![0.0; k];
            let doc_tot: f64 = data.x_nonneg[d].iter().sum();
            for w in 0..v {
                let xdw = data.x_nonneg[d][w];
                let denom: f64 = (0..k).map(|c| theta[c] * beta_kw[c][w]).sum::<f64>().max(1e-12);
                for c in 0..k {
                    let post = xdw * theta[c] * beta_kw[c][w] / denom;
                    phi_sum[c] += post;
                    new_beta[c][w] += post;
                }
            }
            // Gradient ascent step on eta[d] under logistic-normal prior N(mu_prior, cov_prior)
            for c in 0..k {
                let prior_grad: f64 = (0..k)
                    .map(|b| inv_cov[c][b] * (eta[d][b] - mu_prior[b]))
                    .sum();
                let grad = phi_sum[c] - doc_tot * theta[c] - 0.15 * prior_grad;
                eta[d][c] = (eta[d][c] + 0.12 * grad).clamp(-5.0, 5.0);
            }
        }
        // Update mu_prior, correlated covariance cov_prior, and topic-word beta
        for c in 0..k {
            mu_prior[c] = (0..n).map(|d| eta[d][c]).sum::<f64>() / n as f64;
            let sb: f64 = new_beta[c].iter().sum::<f64>().max(1e-12);
            for w in 0..v {
                beta_kw[c][w] = new_beta[c][w] / sb;
            }
        }
        for a in 0..k {
            for b in 0..k {
                let mut s = 0.0;
                for d in 0..n {
                    s += (eta[d][a] - mu_prior[a]) * (eta[d][b] - mu_prior[b]);
                }
                cov_prior[a][b] = s / n as f64 + if a == b { 0.1 } else { 0.0 };
            }
        }
    }

    let labels: Vec<i32> = (0..n)
        .map(|d| {
            (0..k)
                .max_by(|&a, &b| eta[d][a].total_cmp(&eta[d][b]))
                .unwrap_or(0) as i32
        })
        .collect();
    AlgorithmOutput::new(
        labels,
        format!("k={}, logistisk-normal prior", k),
        "Modellerte emnekorrelasjoner via kovariansmatrise i logit-rom".to_string(),
    )
}

// 101. Biterm Topic Model (Yan et al. WWW 2013 — corpus-wide feature-pair biterm Gibbs)
fn biterm_topic_model(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let v = data.n_features;
    let k = cfg.effective_k(n);
    let alpha = 0.5_f64;
    let beta = 0.1_f64;
    let max_iter = cfg.max_iters(25);
    let mut rng = ChaCha8Rng::seed_from_u64(cfg.seed + 101);

    // Extract corpus biterms (w1, w2) from active features in each observation
    let mut biterms: Vec<(usize, usize, usize)> = Vec::new(); // (w1, w2, topic)
    let mut doc_biterms: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n];
    for d in 0..n {
        let active: Vec<usize> = (0..v).filter(|&w| data.counts[d][w] > 0).collect();
        for i in 0..active.len() {
            for j in (i + 1)..active.len().min(i + 5) {
                let t = rng.gen_range(0..k);
                biterms.push((active[i], active[j], t));
                doc_biterms[d].push((active[i], active[j]));
            }
        }
        if active.len() == 1 {
            let t = rng.gen_range(0..k);
            biterms.push((active[0], active[0], t));
            doc_biterms[d].push((active[0], active[0]));
        }
    }
    if biterms.is_empty() {
        return AlgorithmOutput::new(vec![0; n], format!("k={}", k), "Ingen aktive biterms".to_string());
    }

    let mut n_z = vec![0usize; k];
    let mut n_zw = vec![vec![0usize; v]; k];
    for &(w1, w2, t) in &biterms {
        n_z[t] += 1;
        n_zw[t][w1] += 1;
        n_zw[t][w2] += 1;
    }

    let vb = v as f64 * beta;
    for _ in 0..max_iter {
        for b in 0..biterms.len() {
            let (w1, w2, old_t) = biterms[b];
            n_z[old_t] -= 1;
            n_zw[old_t][w1] -= 1;
            n_zw[old_t][w2] -= 1;

            let mut probs = vec![0.0; k];
            let mut sp = 0.0;
            for t in 0..k {
                let denom = (2 * n_z[t]) as f64 + vb;
                let p = (n_z[t] as f64 + alpha)
                    * ((n_zw[t][w1] as f64 + beta) * (n_zw[t][w2] as f64 + beta))
                    / (denom * (denom + 1.0));
                probs[t] = p;
                sp += p;
            }
            let mut r = rng.gen::<f64>() * sp;
            let mut new_t = k - 1;
            for t in 0..k {
                r -= probs[t];
                if r <= 0.0 {
                    new_t = t;
                    break;
                }
            }
            biterms[b].2 = new_t;
            n_z[new_t] += 1;
            n_zw[new_t][w1] += 1;
            n_zw[new_t][w2] += 1;
        }
    }

    // Infer P(z|d) from learned P(z) and P(w|z)
    let tot_b = biterms.len() as f64 + k as f64 * alpha;
    let p_z: Vec<f64> = (0..k).map(|t| (n_z[t] as f64 + alpha) / tot_b).collect();
    let p_wz: Vec<Vec<f64>> = (0..k)
        .map(|t| {
            let den = (2 * n_z[t]) as f64 + vb;
            (0..v).map(|w| (n_zw[t][w] as f64 + beta) / den).collect()
        })
        .collect();

    let labels: Vec<i32> = (0..n)
        .map(|d| {
            let mut scores = vec![0.0; k];
            for &(w1, w2) in &doc_biterms[d] {
                let den: f64 = (0..k).map(|t| p_z[t] * p_wz[t][w1] * p_wz[t][w2]).sum::<f64>().max(1e-12);
                for t in 0..k {
                    scores[t] += p_z[t] * p_wz[t][w1] * p_wz[t][w2] / den;
                }
            }
            (0..k)
                .max_by(|&a, &b| scores[a].total_cmp(&scores[b]))
                .unwrap_or(0) as i32
        })
        .collect();
    AlgorithmOutput::new(
        labels,
        format!("k={}, {} biterms", k, biterms.len()),
        "Korpus-globale ordpar (biterms) via Gibbs-sampling".to_string(),
    )
}

// 102. SNOB (Wallace & Boulton 1968 — Minimum Message Length (MML) mixture model selection)
fn snob_mml(data: &PreparedData, cfg: &RunConfig) -> AlgorithmOutput {
    let n = data.n_samples;
    let d = data.n_features;
    let max_k = (cfg.effective_k(n) + 3).min(8).min(n / 4).max(2);

    let mut best_mml = f64::INFINITY;
    let mut best_labels = vec![0i32; n];
    let mut best_k = 2;

    // Wallace-Freeman MML87 message length approximation for diagonal Gaussian mixture
    for k in 2..=max_k {
        let (lbl, log_lik, _) = fit_gmm_model(&data.x, k, CovModel::Vvi, 25, cfg.seed + 102 + k as u64);
        let p = (k - 1 + 2 * k * d) as f64;
        // Lattice constant c_p ≈ 1 / (2*pi*e)
        let fisher_penalty = 0.5 * p * ((n as f64 / (2.0 * PI * std::f64::consts::E)).max(1.0).ln());
        // Assignment & parameter code length (nats)
        let mml_len = -log_lik + fisher_penalty + 0.5 * (k as f64) * (n as f64).ln();
        if mml_len < best_mml {
            best_mml = mml_len;
            best_labels = lbl;
            best_k = k;
        }
    }
    AlgorithmOutput::new(
        best_labels,
        format!("MML87 søk k=2..{}", max_k),
        format!("Minimum Message Length valgte k={} (MML={:.1} nats)", best_k, best_mml),
    )
}
