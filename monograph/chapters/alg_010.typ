#import "../template/ieee.typ": *

#alg_entry(
  id: 10,
  name: "X-Means",
  category: "Cat. I: Centroid & Partitioning",
  year: "2000",
  paradigm: "BIC-Guided Divisive Model-Selection K-Means",
  complexity: $O(n k_max d T)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

A central limitation of Lloyd's K-Means (Alg. 1) and Bisecting K-Means (Alg. 9) is that the cluster count $k$ must be supplied a priori. *X-Means*, introduced by Pelleg and Moore @pelleg2000xmeans, automates model-order selection across an admissible range $[k_min, k_max]$ by embedding K-Means inside a regularized statistical model-selection loop driven by Schwarz's *Bayesian Information Criterion (BIC)* @schwarz1978bic.

Given a dataset $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$, X-Means associates each candidate $k$-partition $cal(C) = {C_1, dots, C_k}$ (with centroids $bold(mu) = (mu_1, dots, mu_k)$ and sizes $n_j = |C_j|$) with an *identical-spherical Gaussian generative surrogate* $cal(M)_k$: each cluster $C_j$ is modeled by $cal(N)(mu_j, hat(sigma)^2 I_d)$ with empirical mixing proportion $hat(pi)_j = n_j / n$ and pooled isotropic variance estimate. Under the exact maximum-likelihood variance estimator
$ hat(sigma)_"MLE"^2 = 1 / (n d) sum_(j=1)^k sum_(i in C_j) norm(x_i - mu_j)_2^2 = (J^*(cal(C))) / (n d), $
the hard-assignment (classification) log-likelihood of $cal(X)$ under $cal(M)_k$ evaluates in closed form to
$ hat(ell)(cal(X) | cal(M)_k) &= sum_(j=1)^k sum_(i in C_j) [ log hat(pi)_j - d / 2 log(2 pi hat(sigma)_"MLE"^2) \
  &quad - (norm(x_i - mu_j)_2^2) / (2 hat(sigma)_"MLE"^2) ] \
  &= sum_(j=1)^k n_j log(n_j / n) - (n d) / 2 log(2 pi hat(sigma)_"MLE"^2) - (n d) / 2. $
(Pelleg and Moore @pelleg2000xmeans alternatively employ the degrees-of-freedom-adjusted variance $hat(sigma)^2 = J^*(cal(C)) / (d (n - k))$, shifting the final constant to $-d(n - k) / 2$ while preserving the same asymptotic BIC scaling.) Counting $k - 1$ independent mixing weights, $k d$ centroid coordinates, and $1$ shared isotropic variance gives $p_k = k(d + 1)$ free parameters. In maximization sign convention @pelleg2000xmeans, the *Bayesian Information Criterion* is
$ "BIC"(cal(M)_k) = hat(ell)(cal(X) | cal(M)_k) - (p_k) / 2 log n. $

#definition(title: "Local BIC Split Criterion")[
  Let $C subset.eq {1, dots, n}$ be a parent cluster of size $n_C = |C| >= 3$, and let $(A, B)$ be a candidate 2-means bisection of $C$. Denote by $"BIC"_1 (C)$ the 1-component isotropic Gaussian BIC on $C$ ($p_1 = d + 1$) and by $"BIC"_2 (A, B)$ the 2-component shared-variance BIC on $C$ ($p_2 = 2(d + 1)$). The local split is *BIC-justified* if and only if
  $ Delta "BIC"(C -> A, B) := "BIC"_2 (A, B) - "BIC"_1 (C) > 0. $
]

X-Means alternates two phases—*Improve-Params* (global Lloyd refinement) and *Improve-Structure* (local BIC-gated binary splitting)—until no cluster justifies splitting or $k = k_max$ (Algorithm 1).

#algbox(
  title: "Pelleg–Moore X-Means Algorithm",
  inputs: [Dataset $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$, bounds $1 <= k_min <= k_max$, Lloyd iteration cap $T$.],
  outputs: [Selected cluster count $k^*$, partition $cal(C)^*$, centroids $bold(mu)^*$.],
)[
  1. Initialize $k <- k_min$ and seed $bold(mu)$ (e.g., via K-Means++ @arthur2007kmeanspp); set $"best_bic" <- -infinity$.
  2. *Repeat:*
  3. #h(1em) *Improve-Params:* Run Lloyd's K-Means (Alg. 1) globally from $bold(mu)$ for at most $T$ iterations to obtain $(cal(C), bold(mu))$.
  4. #h(1em) Evaluate global $"BIC"(cal(M)_k)$; if $"BIC"(cal(M)_k) > "best_bic"$, update $(k^*, cal(C)^*, bold(mu)^*)$ and $"best_bic"$.
  5. #h(1em) *If* $k = k_max$, *break*.
  6. #h(1em) *Improve-Structure:* For each cluster $C_j$ with $|C_j| >= 3$:
     - Spawn child seeds $mu_j plus.minus delta v_j$ along a direction $v_j$ (principal axis or random unit vector, $delta prop sqrt("SSE"(C_j) / (n_j d))$).
     - Run local 2-means restricted to $C_j$ yielding $(A_j, B_j)$ and compute $Delta "BIC"(C_j -> A_j, B_j)$.
  7. #h(1em) Accept splits with $Delta "BIC" > 0$ in descending order of $Delta "BIC"$ up to budget $k_max - k$, replacing $mu_j$ with $(mu_(A_j), mu_(B_j))$ and updating $k$.
  8. *Until* no cluster splits.
  9. *Return* $(k^*, cal(C)^*, bold(mu)^*)$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five formal theoretical results for X-Means: (1) Schwarz's Laplace asymptotic derivation of BIC and its mixture-regularity caveat; (2) an exact closed-form variance-reduction threshold for local splitting proving asymptotic protection against over-splitting pure isotropic Gaussians; (3) an original analytical counterexample proving *statistical inconsistency* via anisotropic over-splitting even at modest $2 : 1$ aspect ratios; (4) exact termination and doubling-versus-incremental complexity bounds; and (5) exact similarity-group invariance alongside a $1 / n$ model-order outlier breakdown proof.

#subsec("A", "Laplace Evidence Asymptotics of BIC")

#theorem(title: "Laplace Asymptotic Justification of BIC", origin: "Classical — Cited")[
  Let $cal(M)$ be a regular $p$-parameter statistical model with twice continuously differentiable log-likelihood $ell_n (theta)$ on $n$ i.i.d. observations and smooth positive prior $pi(theta)$ at the interior MLE $hat(theta)_n$. If the observed Fisher information $H_n = -nabla^2 ell_n (hat(theta)_n)$ satisfies $H_n / n -> cal(I)(theta^*) succ 0$, then the log marginal evidence obeys @schwarz1978bic @bishop2006pattern
  $ log p(cal(X) | cal(M)) = ell_n (hat(theta)_n) - p / 2 log n + O(1) quad (n -> infinity). $
]

#proof(label: "Proof Synthesis (Cited — Schwarz 1978)")[
  Expanding $g_n (theta) := ell_n (theta) + log pi(theta)$ to second order around the MLE $hat(theta)_n$ (where $nabla ell_n (hat(theta)_n) = bold(0)$) gives
  $ g_n (theta) &approx ell_n (hat(theta)_n) + log pi(hat(theta)_n) \
    &quad - 1 / 2 (theta - hat(theta)_n)^T H_n (theta - hat(theta)_n). $
  Integrating $exp(g_n (theta))$ over $bb(R)^p$ via the multivariate Gaussian integral yields Laplace's approximation:
  $ p(cal(X) | cal(M)) approx e^(ell_n (hat(theta)_n)) pi(hat(theta)_n) (2 pi)^(p / 2) (det H_n)^(-1 / 2). $
  Taking logarithms and substituting $det H_n = det(n cal(I)(theta^*) + o(n)) = n^p det cal(I)(theta^*) (1 + o(1))$ gives
  $ log p(cal(X) | cal(M)) &= ell_n (hat(theta)_n) - p / 2 log n + log pi(hat(theta)_n) \
    &quad + p / 2 log(2 pi) - 1 / 2 log det cal(I)(theta^*) + o(1), $
  where all terms after $-p / 2 log n$ are $O(1)$ as $n -> infinity$. Dropping $O(1)$ constants recovers $"BIC"(cal(M))$. As noted by McLachlan and Peel @mclachlan2000finite, finite mixtures violate strict Fisher regularity on boundary strata where $pi_j = 0$ or $mu_j = mu_ell$, rendering X-Means' hard-assignment BIC a tractable asymptotic surrogate rather than an exact Bayes factor.
]

#subsec("B", "Exact Local Split Threshold and Isotropic Stability")

#theorem(title: "Exact Split-BIC Threshold and Isotropic Non-Splitting", origin: "Original Derivation")[
  For any parent cluster $C$ of size $n_C >= 3$ partitioned into $(A, B)$ with $alpha = |A| / n_C in (0, 1)$ and retained variance ratio $r := ("SSE"(A) + "SSE"(B)) / "SSE"(C) in (0, 1]$, let $H_b (alpha) = -alpha log alpha - (1 - alpha) log(1 - alpha)$ denote binary Shannon entropy (in nats). Then under exact MLE variances:
  1. $Delta "BIC"(C -> A, B) = -(n_C d) / 2 log r - n_C H_b (alpha) - (d + 1) / 2 log n_C$.
  2. Consequently, $Delta "BIC"(C -> A, B) > 0$ if and only if
     $ r < exp( -(2 H_b (alpha)) / d ) dot n_C^(-(d + 1) / (n_C d)). $
  3. *Isotropic Gaussian Stability:* For i.i.d. draws from a single isotropic Gaussian $cal(N)(mu, sigma^2 I_d)$ in *any* dimension $d >= 1$, symmetric bisection along any unit axis yields asymptotic retained ratio $r_"iso" (d) = 1 - 2 / (pi d) > 4^(-1 / d)$, whence $n_C^(-1) Delta "BIC" -> -d / 2 log(1 - 2 / (pi d)) - log 2 <= -0.186 < 0$ almost surely.
]

#proof[
  (1) Under exact MLE variances $hat(sigma)_1^2 = "SSE"(C) / (n_C d)$ and $hat(sigma)_2^2 = ("SSE"(A) + "SSE"(B)) / (n_C d) = r hat(sigma)_1^2$, the 1-cluster and 2-cluster hard-assignment BIC scores on $C$ are
  $ "BIC"_1 (C) &= -(n_C d) / 2 log(2 pi hat(sigma)_1^2) - (n_C d) / 2 - (d + 1) / 2 log n_C, \
    "BIC"_2 (A, B) &= -n_C H_b (alpha) - (n_C d) / 2 log(2 pi r hat(sigma)_1^2) \
    &quad - (n_C d) / 2 - (d + 1) log n_C, $
  where $|A| log alpha + |B| log(1 - alpha) = -n_C H_b (alpha)$. Subtracting $"BIC"_1 (C)$ from $"BIC"_2 (A, B)$ cancels $-(n_C d) / 2 log(2 pi hat(sigma)_1^2) - (n_C d) / 2$ identically, yielding Claim 1.

  (2) Rearranging $-(n_C d) / 2 log r > n_C H_b (alpha) + (d + 1) / 2 log n_C$ gives $-log r > (2 H_b (alpha)) / d + ((d + 1) / (n_C d)) log n_C$; exponentiating both sides proves Claim 2.

  (3) For $X tilde.op cal(N)(0, I_d)$ bisected symmetrically at $chevron.l x, u chevron.r = 0$ ($norm(u)_2 = 1$), the half-normal mean along $u$ is $bb(E)[Z | Z > 0] = sqrt(2 / pi)$, reducing variance along $u$ from $1$ to $1 - 2 / pi$ while leaving all $d - 1$ orthogonal coordinates unchanged. Thus $r_"iso" (d) = ((d - 1) + (1 - 2 / pi)) / d = 1 - 2 / (pi d)$ and $alpha -> 1 / 2$ a.s. Setting $u = 1 / d in (0, 1]$ and $c = 2 / pi in (0, 1)$, the Mercator expansion gives
  $ -d / 2 log(1 - 2 / (pi d)) &= 1 / 2 sum_(m=1)^infinity (c^m) / m u^(m - 1) \
    &<= 1 / 2 sum_(m=1)^infinity (c^m) / m = -1 / 2 log(1 - 2 / pi). $
  Because $pi > 8 / 3$, we have $1 - 2 / pi > 1 / 4$, hence $-1 / 2 log(1 - 2 / pi) < 1 / 2 log 4 = log 2$. By SLLN and $((d + 1) / (2 n_C)) log n_C -> 0$, $n_C^(-1) Delta "BIC" -> -d / 2 log r_"iso" (d) - log 2 <= -1 / 2 log(1 - 2 / pi) - log 2 < -0.186 < 0$ a.s.
]

#subsec("C", "Anisotropic Over-Splitting Instability (Inconsistency of k)")

While Theorem 2(3) protects pure spherical Gaussians, real clusters are rarely isotropic. The following original result proves that even modest covariance elongation causes X-Means to shatter a *single unimodal Gaussian* with probability tending to $1$.

#theorem(title: "Anisotropic Over-Splitting Inconsistency of X-Means", origin: "Original Counterexample")[
  Let $X_1, dots, X_(n_C)$ be i.i.d. draws from a *single unimodal* elliptical Gaussian $cal(N)(bold(0), "diag"(kappa, 1, dots, 1))$ in dimension $d >= 2$ with variance aspect ratio $kappa > 1$. Define the critical threshold
  $ kappa^*(d) := ((d - 1)(1 - 4^(-1 / d))) / (2 / pi - (1 - 4^(-1 / d))). $
  Then $kappa^*(d)$ is finite and positive for every $d >= 2$ (specifically $kappa^*(2) = pi / (4 - pi) approx 3.66$, decreasing monotonically to $lim_(d -> infinity) kappa^*(d) = pi log 2 approx 2.18$). Whenever $kappa > kappa^*(d)$, bisecting the single cluster along its principal axis $x_1 = 0$ yields $Delta "BIC" -> +infinity$ almost surely as $n_C -> infinity$; hence X-Means is *statistically inconsistent* for $k$ on elliptical clusters.
]

#proof[
  First verify positivity of the denominator: because $4^(-1 / d)$ increases strictly in $d$, $1 - 4^(-1 / d)$ is maximized over $d >= 2$ at $d = 2$, where $1 - 4^(-1 / 2) = 1 / 2 < 2 / pi$ (as $pi < 4$). Hence $2 / pi - (1 - 4^(-1 / d)) > 0$ for all $d >= 2$.

  Now bisect the population along the stationary 2-means CVT hyperplane $x_1 = 0$, which separates symmetric population centroids $mu_plus.minus = (plus.minus sqrt(2 kappa / pi), 0, dots, 0)$. The population variance before splitting is $"Tr"(Sigma) = kappa + (d - 1)$, and the variance removed along axis $1$ is $(2 / pi) kappa$. Therefore, by SLLN, the empirical split proportions obey $hat(alpha) -> 1 / 2$ a.s. and the empirical retained variance ratio satisfies
  $ hat(r)_(n_C) -> r_"aniso" (kappa) := 1 - 2 / pi dot kappa / (kappa + d - 1) quad "a.s." $
  By Theorem 2(1), the per-sample split gain converges almost surely to $-d / 2 log r_"aniso" (kappa) - log 2$, which is strictly positive if and only if $r_"aniso" (kappa) < 4^(-1 / d)$. Substituting $r_"aniso" (kappa)$ gives the equivalent condition
  $ 2 / pi dot kappa / (kappa + d - 1) > 1 - 4^(-1 / d) <==> kappa > kappa^*(d). $
  Thus for any $kappa > kappa^*(d)$, there exists $epsilon > 0$ such that $Delta "BIC" >= epsilon n_C - ((d + 1) / 2) log n_C -> +infinity$ almost surely. Because local 2-means initialized along the principal axis converges to this CVT (or achieves an even lower retained ratio), X-Means accepts the spurious split with probability tending to $1$, overestimating $k = 1$ as $hat(k) >= 2$. Finally, as $d -> infinity$, $1 - 4^(-1 / d) = (log 4) / d + O(d^(-2))$, giving $kappa^*(d) -> (pi log 4) / 2 = pi log 2 approx 2.18$ (a standard-deviation elongation of merely $sqrt(kappa^*) approx 1.48 times$).
]

#subsec("D", "Termination, Complexity Regimes, Invariance, and Breakdown")

#theorem(title: "Finite Termination and Doubling vs. Linear Cost", origin: "Original Proof")[
  Algorithm 1 halts in at most $k_max - k_min + 1$ outer iterations (and at most $ceil(log_2 (k_max / k_min))$ rounds under full doubling). Moreover:
  1. Each Improve-Structure sweep costs $O(n d T)$ flops *independent of $k$*.
  2. Under geometric doubling ($k_r = k_min 2^r$), total runtime is $O(n k_max d T)$—within a constant factor of a single flat $k_max$-means call.
  3. Under single-cluster incremental splitting ($k_r = k_min + r$), total runtime is $Theta(n (k_max^2 - k_min^2) d T)$.
]

#proof[
  Every non-terminal round increases $k$ by at least $1$ up to $k_max$, bounding outer iterations by $k_max - k_min + 1$. During Improve-Structure at round $r$, local 2-means executes on pairwise disjoint clusters $C_1, dots, C_(k_r)$ with $sum_(j=1)^(k_r) |C_j| = n$, costing $sum_j O(2 |C_j| d T) = O(n d T)$ regardless of $k_r$. Improve-Params runs global $k_r$-means in $O(n k_r d T)$ flops. Under doubling, the geometric sum $sum_(r=0)^(L) n d T (k_min 2^r + 1) <= 2 n k_max d T + O(n d T log k_max) = O(n k_max d T)$. Under single-cluster splitting, $sum_(k = k_min)^(k_max) O(n k d T) = Theta(n (k_max^2 - k_min^2) d T)$.
]

#proposition(title: "Similarity Invariance and 1/n Model-Order Breakdown", origin: "Original Proof")[
  1. Under any uniform similarity map $T(x) = s Q x + b$ ($s > 0$, $Q in O(d)$, $b in bb(R)^d$), global $"BIC"(cal(M)_k)$ shifts by $-n d log s$ (independent of $k$) and every local score $Delta "BIC"(C -> A, B)$ is *strictly invariant*.
  2. The model-order breakdown point of X-Means is $epsilon^* = 1 / n$: corrupting a single point $x_1 in C$ ($|C| >= 3$, $"SSE"(C backslash {x_1}) > 0$) by sending $norm(x_1)_2 -> infinity$ forces $Delta "BIC"(C -> {x_1}, C backslash {x_1}) -> +infinity$, guaranteeing a spurious split.
]

#proof[
  (1) Under $T$, every squared distance scales by $s^2$, so $hat(sigma)_"MLE"^2 |-> s^2 hat(sigma)_"MLE"^2$ and $-(n d / 2) log(2 pi s^2 hat(sigma)_"MLE"^2)$ shifts by $-n d log s$, independent of $k$. Furthermore, $r = ("SSE"(A) + "SSE"(B)) / "SSE"(C)$ is invariant to $s > 0$, leaving $Delta "BIC"$ unchanged by Theorem 2(1).

  (2) Let $B = C backslash {x_1}$ with $S_B = "SSE"(B) > 0$. By Lemma 1 of Alg. 9, $"SSE"(C) = S_B + ((n_C - 1) / n_C) norm(x_1 - overline(x)_B)_2^2 -> infinity$ as $norm(x_1)_2 -> infinity$. Peeling $A = {x_1}$ gives $alpha = 1 / n_C$ (fixed) and $r = S_B / "SSE"(C) -> 0^+$. By Theorem 2(1), $Delta "BIC" = -(n_C d / 2) log r - n_C H_b (1 / n_C) - ((d + 1) / 2) log n_C -> +infinity$. Thus a single outlier forces an irreversible false-positive split.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Because Improve-Structure partitions $cal(X)$ across disjoint clusters (Theorem 4), the dominant cost of X-Means lies in the Improve-Params global Lloyd sweeps—which Pelleg and Moore @pelleg2000xmeans accelerate further via multiresolution kd-trees (MRKD-trees) in low-to-moderate dimensions ($d <= 12$). Memory footprint remains strictly linear at $O((n + k_max) d)$. Unlike naive grid search over $k in [k_min, k_max]$ using Silhouette @rousseeuw1987silhouettes ($O(n^2)$), Caliński–Harabasz @calinski1974dendrite, or Davies–Bouldin @davies1979cluster, X-Means reuses warm-started split centroids across rounds and evaluates BIC in $O(k)$ time directly from cluster counts and SSEs.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for X-Means],
    [Time / Space Complexity], [$O(n k_max d T)$ doubling, $O(n k_max^2 d T)$ incremental; $O((n+k_max)d)$ space.],
    [Inductive Bias], [Spherical equal-variance Gaussian mixtures; Voronoi tessellations.],
    [Model-Order Consistency], [Protected on isotropic Gaussians (Theorem 2); inconsistent on ellipses $kappa > kappa^*(d)$ (Theorem 3).],
    [Invariances], [Strictly invariant to similarity group $s Q x + b$; sensitive to anisotropic scaling.],
    [Outlier Robustness], [Model-order breakdown $epsilon^* = 1 / n$ (single outlier triggers spurious split).],
    [Recommended Upgrades], [G-Means @hamerly2003gmeans for elliptical clusters; Mclust/GMM (Alg. 81 @fraley2002mclust) for full covariance BIC.],
  ),
  caption: [Summary of theoretical guarantees and operational profile of X-Means.],
)

#sec("IV", "Engineering and Scientific Applications")

X-Means is widely deployed as an unsupervised exploratory discovery engine whenever cluster cardinality is unknown beforehand and approximate spherical compactness holds after feature whitening.

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Automated Customer & Cohort Discovery:* In CRM analytics where marketing teams lack a pre-specified segment count $k$, X-Means automatically discovers natural coarse-to-medium cohorts within $[k_min, k_max]$ without manual elbow inspection, running in linear memory and millisecond latency.
2. *Autonomous Sensor & Astronomical Catalog Partitioning:* Originally tested on Sloan Digital Sky Survey (SDSS) astrophysical catalogs @pelleg2000xmeans, MRKD-accelerated X-Means scales to millions of multi-band photometric objects where manual $k$-tuning across sky tiles is infeasible.
3. *Adaptive Quantization & Codebook Sizing:* In audio, image patch, and IoT telemetry compression, X-Means dynamically adapts dictionary size $k^*$ to signal complexity rather than over-allocating codewords to homogeneous signals or under-allocating to heterogeneous ones.
4. *Cybersecurity Baseline & Behavioral Profiling:* Enterprise user-behavior analytics (UEBA) profiles distinct operational roles across thousands of hosts whose true mode count evolves over time; similarity invariance (Proposition 1) ensures unit changes across telemetry pipelines never alter split decisions.
5. *Rapid Preprocessing for Downstream Mixture Models:* Because full-covariance EM (Alg. 81 @dempster1977em @fraley2002mclust) is expensive across large $k$ ranges, X-Means provides a fast $O(n k_max d T)$ structural initialization supplying candidate $k^*$ and warm centroids for subsequent full-covariance refinement.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Theorem 3 and Proposition 1 pinpoint three acute failure regimes: (a) *elliptical or correlated clusters*, where aspect ratios exceeding just $sqrt(kappa^*) approx 1.5$–$1.9$ trigger systematic over-splitting into chains of spherical fragments (resolved by G-Means @hamerly2003gmeans or full-covariance GMM BIC @fraley2002mclust); (b) *non-convex shapes* (crescents, rings), which X-Means approximates by tiling many small balls; and (c) *heavy-tailed outliers or noisy subspaces*, where a single extreme observation triggers a false-positive split (Proposition 1) and high noise dimensions weaken variance-drop signals.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Why X-Means Fits / Prerequisite], [Caution / Alternative],
    [CRM Discovery], [Automatic $k$ without manual elbow], [Whiten/PCA correlated features first],
    [Astronomy / SDSS], [MRKD-tree acceleration on compact blobs], [Guard $k_max$ against outlier splits],
    [Telemetry / UEBA], [Scale-invariant BIC split thresholds], [Winsorize extreme spikes (Prop. 1)],
    [Elliptical Biology], [Fast coarse initialization only], [Use G-Means @hamerly2003gmeans or Mclust @fraley2002mclust],
    [Non-Convex Manifolds], [Not structurally matched], [Use DBSCAN (Alg. 31) or Spectral (Alg. 103)],
  ),
  caption: [Application suitability and domain engineering guide for X-Means.],
)

#sec("V", "Structural Relations and Empirical Behavior")

X-Means unifies divisive K-Means (Alg. 9) and penalized information-theoretic model selection (Part V) through an explicit log-distortion and entropy functional.

#proposition(title: "Penalized Log-SSE Equivalence and Structural Position", origin: "Original Derivation")[
  Under the exact MLE isotropic Gaussian likelihood, maximizing global $"BIC"(cal(M)_k)$ over partitions $cal(C)$ and $k in [k_min, k_max]$ is strictly equivalent to minimizing the penalized log-distortion objective
  $ cal(L)_"BIC" (cal(C)) = (n d) / 2 log J^*(cal(C)) + n H(hat(pi)) + (k (d + 1)) / 2 log n, $
  where $H(hat(pi)) = -sum_(j=1)^k (n_j / n) log(n_j / n)$ is the cluster-size Shannon entropy.
]

#proof[
  Substituting $hat(sigma)_"MLE"^2 = J^*(cal(C)) / (n d)$ and $sum_j n_j log(n_j / n) = -n H(hat(pi))$ into $"BIC"(cal(M)_k)$ yields
  $ "BIC"(cal(M)_k) &= -n H(hat(pi)) - (n d) / 2 log( (2 pi) / (n d) J^*(cal(C)) ) \
    &quad - (n d) / 2 - (k (d + 1)) / 2 log n \
    &= C_(n, d) - [ (n d) / 2 log J^*(cal(C)) \
    &quad + n H(hat(pi)) + (k (d + 1)) / 2 log n ], $
  where $C_(n, d) = -(n d / 2)(log(2 pi / (n d)) + 1)$ depends solely on $(n, d)$ and is independent of $(cal(C), k)$. Hence maximizing $"BIC"(cal(M)_k)$ is identically equivalent to minimizing $cal(L)_"BIC" (cal(C))$.
]

Proposition 2 clarifies X-Means' exact structural lineage: (i) compared to Bisecting K-Means (Alg. 9 @steinbach2000comparison), X-Means replaces fixed $k$ with local BIC gating and inserts global Improve-Params Lloyd passes after every splitting round, eliminating the non-Voronoi trapping of Alg. 9's Theorem 2; (ii) compared to ISODATA @ball1965isodata, X-Means replaces ad-hoc variance split thresholds with a theoretically calibrated likelihood-ratio penalty (though omitting ISODATA's merge step); and (iii) compared to G-Means @hamerly2003gmeans (Anderson–Darling normality test @anderson1952asymptotic) and Mclust/ICL @fraley2002mclust @biernacki2006mixmod, X-Means pays for its $O(k d)$ parameter parsimony with Theorem 3's spherical-assumption bias.

On the ClusterOverview 4-dataset empirical benchmark, X-Means achieves rapid execution (mean runtime $0.341$ ms) and delivers its *highest* comparative score on Customers ($"Sil" = 0.393$, outperforming both K-Prototypes at $0.342$ and Bisecting K-Means at $0.351$) by adapting $k^*$ directly to the customer distribution's coarse spherical modes. Conversely, it records lower scores on Iris ($"Sil" = 0.298$, $"ARI" = 0.306$), Multi-Shape ($"Sil" = 0.580$, $"ARI" = "N/A"$), and High-D Subspace+Noise ($"Sil" = 0.253$, $"ARI" = "N/A"$). This empirical split directly validates Theorems 2–3 and Proposition 1: BIC excels when coarse isotropic geometry matches the model surrogate, yet over-splits elongated Iris species, curved Multi-Shape manifolds, and noise-perturbed high-dimensional clusters.
