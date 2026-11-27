#import "../template/ieee.typ": *

#alg_entry(
  id: 165,
  name: "ROCKET-Clustering",
  category: "Cat. VII: Fuzzy, Biclustering & Domain",
  year: "2020",
  paradigm: "Random Convolutional Feature Partitioning",
  complexity: $O(n M m + T K n M)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Prototype-based elastic clustering methods such as DBA (Alg. 163) and Soft-DTW (Alg. 164) achieve strong alignment fidelity but incur $O(m^2)$ per-pair dynamic programming costs and struggle to capture discriminative local motifs occurring at multiple temporal scales simultaneously. Dempster, Petitjean, and Webb @dempster2020rocket demonstrated that convolving time series with a large bank of *Random Convolutional Kernels* (*ROCKET*) followed by global max-pooling and Proportion of Positive Values (PPV) pooling extracts highly expressive, multi-scale representations in strictly linear $O(n M m)$ time without any gradient-based filter training. *ROCKET-Clustering* couples this randomized convolutional feature map $Phi_M : bb(R)^m -> bb(R)^(2 M)$ with standard Euclidean prototype clustering (such as $k$-Means++ @arthur2007kmeanspp @lloyd1982least, Mini-Batch $k$-Means @sculley2010webscale, or spectral clustering @ng2002spectral) in feature space.

Let $X = {x_1, dots, x_n} subset bb(R)^m$ be univariate time series of length $m$. Each random convolutional kernel $theta_j = (ell_j, w_j, b_j, d_j, p_j)$ ($j = 1, dots, M$) is drawn independently from the ROCKET distribution $cal(P)_"ROCKET"$ @dempster2020rocket:
1. *Kernel length*: $ell_j tilde.op "Unif"({7, 9, 11})$;
2. *Zero-mean Gaussian weights*: draw $tilde(w)_j tilde.op cal(N)(0, I_(ell_j))$ and center $w_j = Q_(ell_j) tilde(w)_j = tilde(w)_j - bar(tilde(w))_j bold(1)_(ell_j)$, ensuring $bold(1)_(ell_j)^top w_j = 0$;
3. *Bias*: $b_j tilde.op "Unif"([-1, 1])$;
4. *Exponential dilation*: $d_j = floor(2^(u_j))$ where $u_j tilde.op "Unif"([0, log_2 ((m - 1) / (ell_j - 1))])$ (clamped to $d_j = 1$ if $m <= ell_j$);
5. *Padding flag*: $p_j tilde.op "Bern"(1 / 2)$, appending $P_j = p_j floor((ell_j - 1) d_j / 2)$ zeros to both ends of $x$ to form $tilde(x)^((j)) in bb(R)^(m + 2 P_j)$.

#definition(title: "ROCKET Embedding and ROCKET-Clustering Objective")[
  For kernel $theta_j$ with receptive span $R_j = (ell_j - 1) d_j + 1$ and output length $L_j = m + 2 P_j - R_j + 1 >= 1$, define the dilated convolution feature map $z_j (x) in bb(R)^(L_j)$ by
  $ (z_j (x))_t = b_j + sum_(r = 0)^(ell_j - 1) w_(j, r) tilde(x)_(t + r d_j)^((j)), quad t in {0, dots, L_j - 1}. $
  The two global shift-pooling statistics for kernel $j$ @dempster2020rocket are
  $ phi_j^"max" (x) & = max_(0 <= t < L_j) (z_j (x))_t, \
    phi_j^"PPV" (x) & = 1 / (L_j) sum_(t = 0)^(L_j - 1) bb(I)[(z_j (x))_t > 0]. $
  Stacking across $M$ i.i.d. kernels defines the normalized *ROCKET embedding* $Phi_M : bb(R)^m -> bb(R)^(2 M)$:
  $ Phi_M (x) & = 1 / sqrt(M) (phi_1^"max" (x), phi_1^"PPV" (x), \
    & quad dots, phi_M^"max" (x), phi_M^"PPV" (x))^top. $
  *ROCKET-Clustering* optionally standardizes each feature coordinate across $X$ (to balance unbounded max-pool responses against bounded $[0, 1]$ PPV features) and partitions ${Phi_M (x_i)}_(i = 1)^n$ by minimizing the $k$-Means objective
  $ J_"ROCKET" (cal(C), {mu_k}) = sum_(k = 1)^K sum_(i in C_k) norm(Phi_M (x_i) - mu_k)_2^2. $
]

#algbox(
  title: "ROCKET-Clustering",
  inputs: [Series $X = {x_1, dots, x_n} subset bb(R)^m$, kernel count $M$, clusters $K$, iterations $T$],
  outputs: [Partition $cal(C) = {C_1, dots, C_K}$ and feature-space centroids ${mu_1, dots, mu_K} subset bb(R)^(2 M)$],
)[
  1. *Sample Kernel Bank*: Draw $M$ i.i.d. kernels ${theta_j}_(j = 1)^M tilde.op cal(P)_"ROCKET"$ (lengths, zero-mean Gaussian weights, biases, dilations, and padding flags).
  2. *Convolutional Feature Extraction*: For every $i in {1, dots, n}$ and $j in {1, dots, M}$:
     (a) Compute dilated convolution $z_j (x_i) in bb(R)^(L_j)$ in $O(ell_j L_j)$ flops.
     (b) Extract $phi_j^"max" (x_i) = max_t (z_j (x_i))_t$ and $phi_j^"PPV" (x_i) = L_j^(-1) sum_t bb(I)[(z_j (x_i))_t > 0]$.
  3. *Feature Scaling*: Form $u_i = Phi_M (x_i) in bb(R)^(2 M)$ (and optionally standardize feature-wise to zero mean and unit variance across $X$).
  4. *Prototype Partitioning*: Initialize ${mu_k^((0))}_(k = 1)^K$ via $k$-Means++ seeding on ${u_i}_(i = 1)^n$ and run Lloyd's $k$-Means for at most $T$ iterations until convergence.
  5. *Return*: Output final cluster partition $cal(C)$ and centroids ${mu_k}$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five formal analytical properties: exact DC-offset invariance and bounded shift stability of PPV pooling, uniform Hoeffding concentration of the empirical ROCKET kernel toward a positive-definite expected convolutional kernel, exact arc-cosine geometry of single-scale PPV co-activation, equivalence to randomized Kernel $k$-Means with an $O(log K)$ approximation bound, and an exact $R_j$-dependent variance collapse theorem explaining high-dimensional non-sequential degeneracy.

#subsec("A", "DC-Offset Invariance and Shift-Pooling Stability")

#theorem(title: "Exact Additive Offset Invariance and Shift Stability of PPV", origin: "Original Proof")[
  Let $theta_j = (ell_j, w_j, b_j, d_j, p_j)$ be any ROCKET kernel with receptive span $R_j = (ell_j - 1) d_j + 1 <= m$ and unpadded flag $p_j = 0$ (so $L_j = m - R_j + 1$). Then:
  (i) *Exact DC invariance*: For every $x in bb(R)^m$ and every constant shift $c in bb(R)$, $z_j (x + c bold(1)_m) = z_j (x)$ identically, whence $phi_j^"max" (x + c bold(1)_m) = phi_j^"max" (x)$ and $phi_j^"PPV" (x + c bold(1)_m) = phi_j^"PPV" (x)$.
  (ii) *Bounded shift sensitivity*: For any zero-padded lag $s in {-(m - 1), dots, m - 1}$ with shifted series $(T_s (x))_t = x_(t - s)$ on $0 <= t - s < m$ and $0$ otherwise,
  $ |phi_j^"PPV" (x) - phi_j^"PPV" (T_s (x))| <= min(1, (|s|) / (L_j)). $
]

#proof[
  *(i)* Because $p_j = 0$, every convolution window lies entirely within ${0, dots, m - 1}$. Expanding $z_j (x + c bold(1)_m)$ at any $t in {0, dots, L_j - 1}$ gives
  $ (z_j (x + c bold(1)_m))_t & = b_j + sum_(r = 0)^(ell_j - 1) w_(j, r) (x_(t + r d_j) + c) \
    & = (z_j (x))_t + c (bold(1)_(ell_j)^top w_j). $
  Because $w_j = Q_(ell_j) tilde(w)_j$ satisfies $bold(1)_(ell_j)^top w_j = 0$, the additive term vanishes identically at every position $t$. Taking maximum and positive-fraction averages preserves equality. (For padded kernels $p_j = 1$, exact invariance holds on all interior positions $P_j <= t < L_j - P_j$, while boundary windows overlapping zero-padding reflect signal level relative to zero.)

  *(ii)* Assume $s >= 0$ without loss of generality (the case $s < 0$ is symmetric by swapping $x$ and $T_s (x)$). If $s >= L_j$, the bound $1$ holds trivially since $phi_j^"PPV" in [0, 1]$. For $0 <= s < L_j$, consider any output index $t in {s, dots, L_j - 1}$. Every receptive-field sample $t + r d_j$ ($0 <= r <= ell_j - 1$) satisfies $s <= t + r d_j <= (L_j - 1) + (R_j - 1) = m - 1$, so $(T_s (x))_(t + r d_j) = x_(t - s + r d_j)$ identically. Consequently $(z_j (T_s (x)))_t = (z_j (x))_(t - s)$ for all $t in {s, dots, L_j - 1}$. Canceling the $L_j - s$ matching indicators in the PPV sums yields
  $ phi_j^"PPV" (T_s (x)) - phi_j^"PPV" (x) & = 1 / (L_j) (sum_(t = 0)^(s - 1) bb(I)[(z_j (T_s (x)))_t > 0] \
    & quad - sum_(u = L_j - s)^(L_j - 1) bb(I)[(z_j (x))_u > 0]). $
  Because each sum lies in $[0, s]$, their difference lies in $[-s, s]$, proving $|phi_j^"PPV" (x) - phi_j^"PPV" (T_s (x))| <= s / L_j$.
]

#subsec("B", "Uniform Kernel Concentration of Random Convolutional Features")

#theorem(title: "Uniform Convergence to the Expected Convolutional Kernel", origin: "Original Proof")[
  Define the empirical PPV kernel $k_M^"PPV" (x, y) = M^(-1) sum_(j = 1)^M phi_j^"PPV" (x) phi_j^"PPV" (y)$ and the limiting expected kernel $k_infinity^"PPV" (x, y) = bb(E)_(theta tilde.op cal(P)_"ROCKET") [phi_theta^"PPV" (x) phi_theta^"PPV" (y)]$. Then:
  (i) $k_infinity^"PPV"$ is a symmetric positive semidefinite kernel bounded in $[0, 1]$;
  (ii) For any fixed dataset $X = {x_1, dots, x_n} subset bb(R)^m$ and $epsilon > 0$,
  $ "Pr"(max_(1 <= a, b <= n) & |k_M^"PPV" (x_a, x_b) - k_infinity^"PPV" (x_a, x_b)| >= epsilon) \
    & <= n (n + 1) exp(-2 M epsilon^2); $
  (iii) Whenever $M >= (2 epsilon^2)^(-1) log(n (n + 1) / delta)$, with probability at least $1 - delta$ all pairwise squared PPV distances satisfy $|norm(Phi_M^"PPV" (x_a) - Phi_M^"PPV" (x_b))_2^2 - d_infinity^2 (x_a, x_b)| <= 4 epsilon$, where $d_infinity^2 (x_a, x_b) = k_infinity^"PPV" (x_a, x_a) + k_infinity^"PPV" (x_b, x_b) - 2 k_infinity^"PPV" (x_a, x_b)$.
]

#proof[
  *(i)* Symmetry is immediate, and for any $c_1, dots, c_n in bb(R)$,
  $ sum_(a, b = 1)^n c_a c_b k_infinity^"PPV" (x_a, x_b) & = bb(E)_theta [(sum_(a = 1)^n c_a phi_theta^"PPV" (x_a))^2] \
    & >= 0, $
  so $k_infinity^"PPV"$ is positive semidefinite (Mercer-type @scholkopf1998nonlinear), with values in $[0, 1]$ since $phi_theta^"PPV" in [0, 1]$.

  *(ii)* For any fixed pair $(a, b)$, the summands $Y_j = phi_j^"PPV" (x_a) phi_j^"PPV" (x_b) in [0, 1]$ are i.i.d. across $j = 1, dots, M$ with mean $k_infinity^"PPV" (x_a, x_b)$. By Hoeffding's inequality for $[0, 1]$-bounded variables,
  $ "Pr"(|k_M^"PPV" (x_a, x_b) - k_infinity^"PPV" (x_a, x_b)| >= epsilon) <= 2 exp(-2 M epsilon^2). $
  There are $n (n + 1) / 2$ distinct unordered pairs ${a, b}$ (including diagonals $a = b$). Applying the union bound gives failure probability at most $2 dot.op (n (n + 1) / 2) exp(-2 M epsilon^2) = n (n + 1) exp(-2 M epsilon^2)$.

  *(iii)* Setting $n (n + 1) exp(-2 M epsilon^2) <= delta$ yields $M >= (2 epsilon^2)^(-1) log(n (n + 1) / delta)$. On the complementary high-probability event, expanding $norm(Phi_M^"PPV" (x_a) - Phi_M^"PPV" (x_b))_2^2 = k_M^"PPV" (x_a, x_a) + k_M^"PPV" (x_b, x_b) - 2 k_M^"PPV" (x_a, x_b)$ and applying the triangle inequality bounds the squared-distance distortion uniformly by $epsilon + epsilon + 2 epsilon = 4 epsilon$.
]

#subsec("C", "Arc-Cosine Local Correlation Geometry of PPV Features")

#lemma(title: "Exact Bias-Averaged Ramp and Zero-Bias Arc-Cosine Kernel", origin: "Original Derivation")[
  Fix kernel length $ell$, dilation $d$, and two dilated local patches $u, v in bb(R)^ell$ with centered projections $tilde(u) = Q_ell u$, $tilde(v) = Q_ell v$ ($tilde(u), tilde(v) != 0$). For centered Gaussian weights $w = Q_ell tilde(w)$ ($tilde(w) tilde.op cal(N)(0, I_ell)$) and $b tilde.op "Unif"([-1, 1])$:
  (i) $bb(E)_b [bb(I)[w^top u + b > 0] | w] = "clip"_([0, 1]) ((1 + w^top tilde(u)) / 2)$, and $bb(E)_(w, b) [bb(I)[w^top u + b > 0]] = 1 / 2$;
  (ii) At threshold $b = 0$, the expected patch co-activation probability equals the exact arc-cosine similarity of local Pearson correlation $rho(tilde(u), tilde(v)) = chevron.l tilde(u), tilde(v) chevron.r / (norm(tilde(u))_2 norm(tilde(v))_2)$:
  $ "Pr"_w (w^top u > 0 "and" w^top v > 0) = 1 / 2 - (arccos(rho(tilde(u), tilde(v)))) / (2 pi). $
]

#proof[
  *(i)* Because $Q_ell^top = Q_ell$, we have $w^top u = tilde(w)^top Q_ell u = tilde(w)^top tilde(u) = w^top tilde(u)$. For $b tilde.op "Unif"([-1, 1])$, $"Pr"_b (b > -w^top tilde(u)) = (1 - "clip"_([-1, 1]) (-w^top tilde(u))) / 2 = "clip"_([0, 1]) ((1 + w^top tilde(u)) / 2)$. Since the distribution of $w$ is symmetric under $w |-> -w$ and $"clip"_([0, 1]) ((1 + t) / 2) + "clip"_([0, 1]) ((1 - t) / 2) = 1$ for all $t in bb(R)$, taking expectation over $w$ yields $1 / 2$.

  *(ii)* Because $tilde(w) tilde.op cal(N)(0, I_ell)$ is rotationally invariant, projecting $tilde(w)$ onto $"span"{tilde(u), tilde(v)}$ yields a standard bivariate Gaussian $(g_1, g_2)$ whose signed inner products with unit vectors $hat(u) = tilde(u) / norm(tilde(u))_2$ and $hat(v) = tilde(v) / norm(tilde(v))_2$ are simultaneously positive iff the angular coordinate of $(g_1, g_2)$ falls in the open intersection of the two half-planes $chevron.l g, hat(u) chevron.r > 0$ and $chevron.l g, hat(v) chevron.r > 0$. With inter-vector angle $alpha = arccos(chevron.l hat(u), hat(v) chevron.r) = arccos(rho(tilde(u), tilde(v))) in [0, pi]$, this angular wedge has Lebesgue measure $pi - alpha$, giving probability $(pi - alpha) / (2 pi) = 1 / 2 - arccos(rho(tilde(u), tilde(v))) / (2 pi)$.
]

Lemma 1 uncovers the mathematical mechanism behind ROCKET-Clustering: PPV aggregates local dilated-patch Pearson correlations via an arc-cosine kernel across exponentially spaced scales, while random biases $b in [-1, 1]$ inject amplitude sensitivity through the ramp $(1 + w^top tilde(u)) / 2$.

#subsec("D", "Equivalence to Kernel k-Means and Approximation Bound")

#theorem(title: "Randomized Kernel k-Means Equivalence and O(log K) Bound", origin: "Classical — Cited")[
  Conditional on the drawn kernel bank ${theta_j}_(j = 1)^M$, let $K_M (x_a, x_b) = chevron.l Phi_M (x_a), Phi_M (x_b) chevron.r$ be the empirical ROCKET kernel. Then:
  (i) Minimizing $J_"ROCKET" (cal(C), {mu_k})$ is algebraically identical to Kernel $k$-Means @dhillon2004kernel with kernel matrix $K_M$;
  (ii) Initializing with $k$-Means++ @arthur2007kmeanspp and running Lloyd refinement @lloyd1982least decreases $J_"ROCKET"$ monotonically at every step and satisfies
  $ bb(E)[J_"ROCKET" (cal(C)_"final")] <= 8 (ln K + 2) min_(cal(C), {mu_k}) J_"ROCKET" (cal(C), {mu_k}). $
]

#proof[
  *(i)* For any fixed partition $cal(C)$, the optimal feature-space centroid is $mu_k = |C_k|^(-1) sum_(i in C_k) Phi_M (x_i)$. Substituting into $J_"ROCKET"$ yields
  $ sum_(k = 1)^K sum_(i in C_k) norm(Phi_M (x_i) - mu_k)_2^2 & = sum_(i = 1)^n K_M (x_i, x_i) \
    & quad - sum_(k = 1)^K 1 / (|C_k|) sum_(i, j in C_k) K_M (x_i, x_j), $
  which is precisely the Kernel $k$-Means objective @dhillon2004kernel. *(ii)* Because ${Phi_M (x_i)}_(i = 1)^n$ lies in Euclidean space $bb(R)^(2 M)$, the $k$-Means++ seeding theorem @arthur2007kmeanspp guarantees $bb(E)[J_"ROCKET" (cal(C)_"init")] <= 8 (ln K + 2) J_"ROCKET"^ast$. Since each subsequent Lloyd assignment and centroid update satisfies $J_"ROCKET" (cal(C)_"final") <= J_"ROCKET" (cal(C)_"init")$ deterministically, taking expectations completes the proof.
]

#subsec("E", "High-Dimensional Concentration on Non-Sequential Data")

#proposition(title: "Variance Collapse of PPV under Isotropic Non-Sequential Inputs", origin: "Original Derivation")[
  Suppose $x = (x_0, dots, x_(m - 1))^top tilde.op cal(N)(0, I_m)$ has i.i.d. coordinates without temporal autocorrelation @beyer1999nearest. Fix any unpadded kernel $theta_j$ with receptive span $R_j = (ell_j - 1) d_j + 1 <= m$ and output length $L_j = m - R_j + 1$. Let $Phi_cal(N)$ be the standard normal CDF. Then:
  (i) $bb(E)[phi_j^"PPV" (x) | theta_j] = Phi_cal(N) (b_j / norm(w_j)_2)$, independent of $x$;
  (ii) $"Var"(phi_j^"PPV" (x) | theta_j) <= (2 R_j - 1) / (4 L_j)$, so for independent draws $x_a, x_b tilde.op cal(N)(0, I_m)$,
  $ bb(E)[(phi_j^"PPV" (x_a) - phi_j^"PPV" (x_b))^2 | theta_j] <= (2 R_j - 1) / (2 L_j) -> 0 $
  as $L_j / R_j -> infinity$.
]

#proof[
  *(i)* For each $t in {0, dots, L_j - 1}$, $(z_j (x))_t = b_j + sum_(r = 0)^(ell_j - 1) w_(j, r) x_(t + r d_j) tilde.op cal(N)(b_j, norm(w_j)_2^2)$ because coordinates of $x$ are i.i.d. standard normal. Hence $Y_t = bb(I)[(z_j (x))_t > 0]$ is Bernoulli with success probability $"Pr"((z_j (x))_t > 0) = Phi_cal(N) (b_j / norm(w_j)_2)$ for all $t$, giving $bb(E)[phi_j^"PPV" (x) | theta_j] = L_j^(-1) sum_t bb(E)[Y_t] = Phi_cal(N) (b_j / norm(w_j)_2)$.

  *(ii)* Two indicators $Y_t$ and $Y_(t')$ depend on disjoint coordinate sets whenever $|t - t'| >= R_j$, so $"Cov"(Y_t, Y_(t')) = 0$ for $|t - t'| >= R_j$. For the at most $2 R_j - 1$ pairs with $|t - t'| < R_j$ per index $t$, Cauchy-Schwarz and $"Var"(Y_t) <= 1 / 4$ give $|"Cov"(Y_t, Y_(t'))| <= 1 / 4$. Therefore
  $ "Var"(phi_j^"PPV" (x) | theta_j) & = 1 / (L_j^2) sum_(t = 0)^(L_j - 1) sum_(|t' - t| < R_j) "Cov"(Y_t, Y_(t')) \
    & <= (L_j (2 R_j - 1)) / (4 L_j^2) = (2 R_j - 1) / (4 L_j). $
  For independent $x_a, x_b$, $bb(E)[(phi_j^"PPV" (x_a) - phi_j^"PPV" (x_b))^2 | theta_j] = 2 "Var"(phi_j^"PPV" (x) | theta_j) <= (2 R_j - 1) / (2 L_j)$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Because every ROCKET kernel has bounded length $ell_j <= 11$, applying one kernel to one series of length $m$ requires at most $11 m$ multiply-add operations, so extracting all $2 M$ features for $n$ series costs $O(n M m)$ time — strictly linear in *both* dataset size $n$ and series length $m$, with zero dependence on $m^2$ or $m log m$! Subsequent $k$-Means++ clustering in $bb(R)^(2 M)$ costs $O(T K n M)$ time, yielding total complexity $O(n M m + T K n M)$. Even more importantly, feature extraction is embarrassingly parallel across both $n$ series and $M$ kernels, vectorizes cleanly into SIMD fused multiply-adds, and avoids dynamic-programming loop dependencies entirely.

Memory footprint is $O(n M + K M + n m)$ when materializing the dense $n times 2 M$ feature matrix (or $O(n r + K r + n m)$ when streaming features through an online $r$-dimensional randomized PCA / Johnson-Lindenstrauss sketch prior to clustering). Geometrically, ROCKET-Clustering replaces pairwise elastic alignment with a fixed random-feature linear embedding; clusters correspond to Euclidean Voronoi polyhedra in $bb(R)^(2 M)$, which pull back to rich, multi-scale nonlinear decision regions in raw time-series space $bb(R)^m$ without constructing an explicit $O(n^2)$ Gram matrix.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Representation], [Random dilated 1D convolutions + global max & PPV pooling],
    [Clustering Engine], [$k$-Means++ / Lloyd (or Mini-Batch / Spectral) in $bb(R)^(2 M)$],
    [Time Complexity], [$O(n M m + T K n M)$ — linear in $n$ and $m$],
    [Space Complexity], [$O(n M + K M + n m)$ (reducible via streaming PCA)],
    [Invariances], [Exact DC-offset invariance ($p_j = 0$); bounded shift stability],
    [Approximation], [$8 (ln K + 2)$ expected $k$-Means++ factor in RKHS (Theorem 3)],
    [Interpretability], [Feature-space centroids lack direct raw time-domain waveforms],
  ),
  caption: [Operational and geometric profile of ROCKET-Clustering (Alg. 165).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *High-Frequency Industrial IoT and Predictive Maintenance*: Accelerometer, current, and acoustic streams sampled at $10$--$50$ kHz produce series lengths $m >= 10^4$ where quadratic DTW methods (Algs. 163--164) are infeasible. ROCKET's exponentially spaced dilations simultaneously detect narrow transient impulses (small $d_j$) and low-frequency envelope modulations (large $d_j$) in $O(n M m)$ time.
2. *Large-Scale Wearable Sensor Activity Discovery*: Continuous smartwatch and IMU recordings contain recurring micro-motifs (steps, arm swings, tremors) appearing at arbitrary phases within long windows. Global PPV pooling encodes motif frequency and duty cycle robustly across phase shifts (Theorem 1), outperforming rigid whole-series metrics when activities repeat asynchronously.
3. *Wearable ECG / PPG Arrhythmia and Sleep-Stage Triage*: Multi-scale dilated kernels capture both fine QRS spike slopes and broader P-R / Q-T intervals simultaneously while ignoring baseline DC drift ($bold(1)^top w_j = 0$), enabling millisecond-scale unsupervised triage over millions of segments.
4. *Acoustic Event, Bioacoustics, and Ultrasonic Inspection*: Short transient calls or structural flaw echoes embedded in long noisy background recordings are poorly matched by whole-series Euclidean or DTW alignment, whereas ROCKET's global max-pooling activates strongly whenever a matching local pattern appears anywhere in the window.
5. *Streaming Network Telemetry and Cloud Microservice Anomaly Clustering*: CPU, latency, and packet-rate time series exhibit bursty multi-scale shapes (spikes, step changes, periodic sawtooths). Because the random kernel bank is data-independent, new streaming series embed in a single forward pass without refitting an alignment dictionary.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

ROCKET-Clustering has four principal limitations: (i) *non-sequential high-dimensional tabular data*, where coordinates lack local neighborhood autocorrelation and PPV features collapse toward deterministic constants $Phi_cal(N) (b_j / norm(w_j)_2)$ (Proposition 1), destroying cluster contrast; (ii) *very short sequences ($m < 7$)*, which fall below the standard ROCKET kernel lengths ${7, 9, 11}$ and force heavy zero-padding or kernel truncation; (iii) *absence of native time-domain prototypes*, since cluster centroids $mu_k in bb(R)^(2 M)$ are not directly invertible to interpretable waveform shapes (requiring post-hoc medoid extraction or DBA refinement within each discovered cluster); and (iv) *unsupervised feature dilution*, because data-independent random kernels weight irrelevant noise frequencies equally with discriminative motifs unless paired with variance filtering or PCA.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [High-Rate IoT / Vibration], [Very High: $O(n M m)$ multi-scale motifs], [Extract medoids for visual prototypes],
    [Wearable Activity Streams], [Very High: phase-robust PPV duty cycles], [Standardize max vs. PPV scales],
    [Transient Acoustic / Echo], [High: global max-pool local detection], [Use MiniRocket for deterministic speed],
    [Short / Tabular Vectors], [Poor: $m < 7$ padding + PPV collapse], [Use $k$-Means++, DBA, or Soft-DTW],
    [Prototype-Driven Science], [Moderate: centroids live in $bb(R)^(2 M)$], [Post-process with DBA (Alg. 163)],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for ROCKET-Clustering.],
)

#sec("V", "Structural Relations and Empirical Behavior")

ROCKET-Clustering bridges classical random Fourier / random feature approximations of kernel machines @scholkopf1998nonlinear @dhillon2004kernel with domain-specific time-series clustering (Algs. 162--164).

#proposition(title: "Spectral Reduction to Multi-Scale Random-Feature Spectral Clustering", origin: "Original Derivation")[
  Replacing the Lloyd $k$-Means step in ROCKET-Clustering with normalized Laplacian eigendecomposition of the Nyström-free affinity matrix $W = Phi_M (X) Phi_M (X)^top$ (projected to nonnegative entries, e.g. using PPV blocks $Phi_M^"PPV" (X) in [0, 1]^(n times M)$) reduces ROCKET-Clustering to $O(n M K)$ linear-time spectral clustering @ng2002spectral whose top $K$ left singular vectors are computed directly via thin SVD of $D^(-1 / 2) Phi_M^"PPV" (X)$ without ever forming the $n times n$ similarity matrix $W$.
]

#proof[
  Because every entry of $Phi_M^"PPV" (X) in [0, 1]^(n times M)$ is nonnegative, $W = Phi_M^"PPV" (X) Phi_M^"PPV" (X)^top$ is entrywise nonnegative and symmetric positive semidefinite. Its degree vector $d = W bold(1)_n = Phi_M^"PPV" (X) (Phi_M^"PPV" (X)^top bold(1)_n)$ is computed in $O(n M)$ time. Letting $B = "diag"(d)^(-1 / 2) Phi_M^"PPV" (X) in bb(R)^(n times M)$, the normalized affinity satisfies $D^(-1 / 2) W D^(-1 / 2) = B B^top$, whose leading $K$ eigenvectors coincide identically with the leading $K$ left singular vectors of $B$ computable in $O(n M K)$ time via randomized SVD.
]

The ClusterOverview four-dataset benchmark provides a textbook empirical validation of our theoretical analysis @rousseeuw1987silhouettes @hubert1985comparing. ROCKET-Clustering attains Iris $"Sil" = 0.315$ ($"ARI" = 0.515$), Customers $"Sil" = 0.143$, Multi-Shape $"Sil" = 0.609$ ($"ARI" = "N/A"$), and High-D $"Sil" = 0.002$ ($"ARI" = "N/A"$), with an extraordinary mean runtime of just $1.091$ ms — $13.5 times$ faster than k-Shape ($14.685$ ms), $24 times$ faster than DBA ($26.159$ ms), and $135 times$ faster than Soft-DTW ($147.236$ ms). On Multi-Shape, its nonlinear random-feature embedding achieves the *highest Silhouette* ($0.609$) among all four time-series algorithms (Algs. 162--165), surpassing Soft-DTW ($0.603$), DBA ($0.459$), and k-Shape ($0.200$). Conversely, on the High-D benchmark, its Silhouette collapses to $0.002$ — precisely matching Proposition 1's exact variance-collapse bound $"Var"(phi_j^"PPV" (x)) <= (2 R_j - 1) / (4 L_j)$ when convolutional kernels scan unstructured high-dimensional coordinates lacking local temporal correlation.
