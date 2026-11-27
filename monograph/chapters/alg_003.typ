#import "../template/ieee.typ": *

#alg_entry(
  id: 3,
  name: "Mini-Batch K-Means",
  category: "Cat. I: Centroid & Partitioning",
  year: "2010",
  paradigm: "Per-Center Stochastic Gradient Descent",
  complexity: $O(b k d T)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Let $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$ be a large dataset with $n >> 1$, and let $hat(P)_n = n^(-1) sum_(i=1)^n delta_(x_i)$ denote its empirical measure. Express the Euclidean $k$-means objective in normalized expected-risk form:
$ F(bold(mu)) &:= bb(E)_(x ~ hat(P)_n) [f(x; bold(mu))] \
              &= 1 / (2 n) sum_(i=1)^n min_(1 <= j <= k) norm(x_i - mu_j)_2^2, $
where the per-sample loss is $f(x; bold(mu)) = 1 / 2 norm(x - mu_(c(x; bold(mu))))_2^2$ and $c(x; bold(mu)) = arg min_(1 <= j <= k) norm(x - mu_j)_2^2$ (breaking ties by smallest index). Computing a single exact Lloyd step (Alg. 1) @lloyd1982least requires $Theta(n k d)$ operations over the entire dataset, which becomes prohibitive when $n$ exceeds main-memory capacity or streaming latency budgets.

*Mini-Batch K-Means* @sculley2010webscale replaces full-dataset passes with stochastic mini-batch updates of size $b << n$, combining the low per-step cost of online MacQueen K-Means (Alg. 21) @macqueen1967some with the variance reduction and SIMD vectorization of batch updates, governed by an adaptive per-center learning rate $eta_j = 1 / v_j$.

#definition(title: "Per-Center Stochastic K-Means Operator")[
  Given mini-batch size $b in {1, dots, n}$, cumulative per-center hit counts $v^((t)) = (v_1^((t)), dots, v_k^((t))) in bb(Z)_(>= 0)^k$, and centroids $bold(mu)^((t)) in (bb(R)^d)^k$, round $t >= 1$ draws a uniform random mini-batch $M_t subset.eq {1, dots, n}$ of cardinality $b$, caches assignments $c_i^((t)) = c(x_i; bold(mu)^((t-1)))$ for all $i in M_t$, and partitions $M_t$ into batch slices $M_(j,t) := {i in M_t : c_i^((t)) = j}$ of size $m_(j,t) = |M_(j,t)|$. For every $j$ with $m_(j,t) >= 1$, Sculley's update iterates over $i in M_(j,t)$:
  $ v_j <- v_j + 1, quad eta_j <- 1 / v_j, quad mu_j <- (1 - eta_j) mu_j + eta_j x_i. $
]

Unlike fixed global step sizes $eta_t = 1 / t$ that under-step rare clusters and over-step dense clusters, maintaining independent counters $v_j$ adapts the step schedule to each cluster's empirical visitation rate.

#algbox(
  title: "Sculley's Mini-Batch K-Means Algorithm",
  inputs: [Dataset $cal(X) subset.eq bb(R)^d$, clusters $k >= 1$, batch size $b$, iterations $T$, initial $bold(mu)^((0))$ (via K-Means++ on pilot sample).],
  outputs: [Centroids $bold(mu)^((T))$ and induced Voronoi partition $cal(C)$.],
)[
  1. Initialize per-center counts $v_j <- 0$ for all $j in {1, dots, k}$.
  2. *For* $t = 1, dots, T$ *do*:
  3. #h(1em) Draw a uniform random mini-batch $M_t$ of $b$ indices from ${1, dots, n}$.
  4. #h(1em) *Cache assignments:* For each $i in M_t$, compute
     $ c_i^((t)) <- arg min_(j in {1, dots, k}) norm(x_i - mu_j^((t-1)))_2^2. $
  5. #h(1em) *Per-center stochastic update:* For each $i in M_t$ with $j = c_i^((t))$:
  6. #h(2em) $v_j <- v_j + 1$; #h(0.8em) $eta <- 1 / v_j$;
  7. #h(2em) $mu_j <- (1 - eta) mu_j + eta x_i$.
  8. #h(1em) *If* relative centroid shift $max_j norm(mu_j^((t)) - mu_j^((t-1)))_2 < epsilon$ for consecutive checks, *break*.
  9. *Return* $bold(mu)^((t))$ and final Voronoi assignments $c(x_i; bold(mu)^((t)))$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five formal theoretical properties of Mini-Batch K-Means: (1) exact order-invariance and algebraic equivalence between sequential and vectorized batch updates; (2) exact equivalence to a block-diagonally preconditioned stochastic Newton step; (3) an explicit counterexample proving strict non-monotonicity of full-dataset inertia; (4) an exact $O(1 / (b t))$ local mean-squared error rate in which cluster-imbalance factors cancel globally; and (5) almost-sure convergence under Robbins–Monro conditions.

#subsec("A", "Order-Invariance and Stochastic Preconditioned Newton Identity")

#lemma(title: "Exact Equivalence of Sequential and Aggregated Batch Updates", origin: "Original Proof")[
  Fix round $t >= 1$ and cluster $j$ with prior count $v_j^((t-1)) >= 0$, prior centroid $mu_j^((t-1))$, and cached batch slice $M_(j,t)$ of size $m_(j,t) >= 1$ and slice mean $overline(x)_(M_(j,t)) = m_(j,t)^(-1) sum_(i in M_(j,t)) x_i$. Applying lines 6–7 of Algorithm 1 sequentially across $M_(j,t)$ in *any* permutation yields the exact order-invariant closed form
  $ mu_j^((t)) = (1 - (m_(j,t)) / (v_j^((t)))) mu_j^((t-1)) + (m_(j,t)) / (v_j^((t))) overline(x)_(M_(j,t)), $
  where $v_j^((t)) = v_j^((t-1)) + m_(j,t)$. Moreover, starting from $v_j^((0)) = 0$, $mu_j^((t))$ equals the unweighted arithmetic mean of the multiset $cal(H)_j^((t)) = union.sq_(s=1)^t M_(j,s)$ of all historical samples ever assigned to center $j$.
]

#proof[
  Enumerate $M_(j,t) = {i_1, dots, i_m}$ with $m = m_(j,t)$, and let $u_r, w_r$ denote the count and centroid after processing the $r$-th element ($u_0 = v_j^((t-1))$, $w_0 = mu_j^((t-1))$). We prove by induction on $r in {0, dots, m}$ that
  $ u_r w_r = u_0 w_0 + sum_(ell=1)^r x_(i_ell). $
  For $r = 0$, both sides equal $u_0 w_0$. Assuming the identity at $r - 1$, step $r$ sets $u_r = u_(r-1) + 1$ and $eta_r = 1 / u_r$, whence
  $ u_r w_r &= u_r ((1 - 1 / u_r) w_(r-1) + 1 / u_r x_(i_r)) \
            &= (u_r - 1) w_(r-1) + x_(i_r) = u_(r-1) w_(r-1) + x_(i_r) \
            &= u_0 w_0 + sum_(ell=1)^r x_(i_ell). $
  Setting $r = m$ and dividing by $u_m = v_j^((t)) = v_j^((t-1)) + m_(j,t)$ gives
  $ mu_j^((t)) = (v_j^((t-1)) mu_j^((t-1)) + m_(j,t) overline(x)_(M_(j,t))) / (v_j^((t))), $
  which is strictly symmetric in the elements of $M_(j,t)$ and therefore permutation-invariant. Inducting across rounds $s = 1, dots, t$ with $v_j^((0)) = 0$ yields $v_j^((t)) mu_j^((t)) = sum_(s=1)^t sum_(i in M_(j,s)) x_i$, proving that $mu_j^((t))$ is the exact running arithmetic average over $cal(H)_j^((t))$.
]

#theorem(title: "Preconditioned Stochastic Newton Interpretation", origin: "Original Derivation")[
  Let $cal(U)$ be an open neighborhood of a strict centroidal Voronoi equilibrium $bold(mu)^*$ on which all assignments $c(x_i; bold(mu)) equiv c_i^*$ are locally constant with cluster frequencies $p_j = |C_j^*| / n > 0$. Define the mini-batch stochastic gradient block $g_j^((t)) := b^(-1) sum_(i in M_(j,t)) (mu_j^((t-1)) - x_i)$ and empirical curvature block $H_j^((t)) := (v_j^((t)) / (b t)) I_d$. Then whenever $v_j^((t)) >= 1$, Algorithm 1 satisfies
  $ mu_j^((t)) = mu_j^((t-1)) - 1 / t (H_j^((t)))^(-1) g_j^((t)), $
  where $bb(E)[g_j^((t)) | cal(F)_(t-1)] = nabla_(mu_j) F(bold(mu)^((t-1)))$ and $H_j^((t)) -> p_j I_d = nabla^2_(mu_j, mu_j) F(bold(mu)^*)$ almost surely as $t -> infinity$.
]

#proof[
  Substituting the definitions of $H_j^((t))$ and $g_j^((t))$ directly gives
  $ &1 / t (H_j^((t)))^(-1) g_j^((t)) \
    &quad = 1 / t dot (b t) / (v_j^((t))) dot 1 / b sum_(i in M_(j,t)) (mu_j^((t-1)) - x_i) \
    &quad = (m_(j,t)) / (v_j^((t))) (mu_j^((t-1)) - overline(x)_(M_(j,t))). $
  Subtracting this expression from $mu_j^((t-1))$ yields $(1 - m_(j,t)/v_j^((t))) mu_j^((t-1)) + (m_(j,t)/v_j^((t))) overline(x)_(M_(j,t))$, recovering Lemma 1 identically. Because each index in $M_t$ is drawn uniformly from ${1, dots, n}$,
  $ bb(E)[g_j^((t)) | cal(F)_(t-1)] &= 1 / n sum_(i in C_j^*) (mu_j^((t-1)) - x_i) \
                                    &= nabla_(mu_j) F(bold(mu)^((t-1))). $
  Finally, on $cal(U)$ the batch counts $m_(j,s)$ are i.i.d. with mean $b p_j$. By the Strong Law of Large Numbers,
  $ (v_j^((t))) / (b t) = 1 / t sum_(s=1)^t (m_(j,s)) / b -> p_j quad "a.s." $
  Since $F(bold(mu)) = (2 n)^(-1) sum_(j=1)^k sum_(i in C_j^*) norm(x_i - mu_j)_2^2$ on $cal(U)$, its Hessian is $nabla^2_(mu_j, mu_ell) F = p_j delta_(j ell) I_d$. Thus Sculley's per-center step size $1 / v_j^((t))$ automatically applies asymptotically exact inverse-Hessian preconditioning across unbalanced clusters.
]

#subsec("B", "Non-Monotonicity of Full-Dataset Inertia")

#proposition(title: "Strict Non-Monotonicity of Mini-Batch K-Means", origin: "Original Counterexample")[
  Unlike Lloyd's algorithm (Alg. 1, Theorem 1), Mini-Batch K-Means does not guarantee monotonic decrease of $F(bold(mu))$: even when initialized at the unique global optimum $bold(mu)^*$, for any non-trivial batch size $1 <= b < n$ and any initial count $v_0 >= 0$, there exists an instance $cal(X) subset.eq bb(R)$ on which $bb(P)[F(bold(mu)^((1))) > F(bold(mu)^*)] > 0$.
]

#proof[
  Consider $n = 4$ points in $bb(R)$ with $k = 2$: $x_1 = -1$, $x_2 = 1$, $x_3 = 9$, $x_4 = 11$. The unique global optimum is $bold(mu)^* = (0, 10)$ with clusters $C_1^* = {x_1, x_2}$, $C_2^* = {x_3, x_4}$, and objective $F(bold(mu)^*) = 1 / 2$. Initialize $bold(mu)^((0)) = bold(mu)^*$ with equal prior counts $v_1^((0)) = v_2^((0)) = v_0 >= 0$. For $b = 2$, the mini-batch $M_1 = {x_2, x_4}$ occurs with strictly positive probability ($1 / 6$ without replacement, $1 / 8$ with replacement), yielding $m_(1,1) = m_(2,1) = 1$. By Lemma 1,
  $ mu_1^((1)) = 1 / (v_0 + 1), quad mu_2^((1)) = 10 + 1 / (v_0 + 1). $
  The updated Voronoi midpoint is $5 + (v_0 + 1)^(-1) in (5, 6]$, so all four points retain their optimal assignments $C_1^*, C_2^*$. Applying the parallel-axis identity (Alg. 1, Lemma 1) yields
  $ F(bold(mu)^((1))) - F(bold(mu)^*) &= 1 / 4 (mu_1^((1)))^2 + 1 / 4 (mu_2^((1)) - 10)^2 \
                                      &= 1 / (2 (v_0 + 1)^2) > 0. $
  Appending duplicate copies of each cluster extends the construction to arbitrary $n$ and $1 <= b < n$.
]

#subsec("C", "Finite-Sample MSE Rate and Almost-Sure Convergence")

#theorem(title: "Exact O(1 / (b t)) Local MSE and Imbalance-Free Excess Risk", origin: "Original Proof")[
  Suppose Algorithm 1 operates inside a fixed Voronoi basin $cal(U)$ around $bold(mu)^*$ with i.i.d. mini-batches of size $b$ drawn with replacement, $v_j^((0)) = 0$, cluster probabilities $p_j = |C_j^*| / n > 0$, and intra-cluster trace variances $sigma_j^2 = |C_j^*|^(-1) sum_(i in C_j^*) norm(x_i - mu_j^*)_2^2$. Then for every $t >= 1$ and cluster $j$:
  1. $bb(E)[norm(mu_j^((t)) - mu_j^*)_2^2 | v_j^((t)) >= 1] <= (2 sigma_j^2) / ((b t + 1) p_j)$.
  2. On the full-coverage event $E_t = {forall j : v_j^((t)) >= 1}$ (whose probability satisfies $bb(P)(E_t) >= 1 - sum_(j=1)^k (1 - p_j)^(b t)$), the excess normalized inertia satisfies the *cluster-imbalance-free* bounds
     $ bb(E)[(F(bold(mu)^((t))) - F(bold(mu)^*)) bold(1)_(E_t)] &<= 1 / (b t + 1) sum_(j=1)^k sigma_j^2, \
       bb(E)[F(bold(mu)^((t))) - F(bold(mu)^*) | E_t] &<= 1 / (bb(P)(E_t) (b t + 1)) sum_(j=1)^k sigma_j^2. $
]

#proof[
  Across $t$ rounds of $b$ i.i.d. draws from $hat(P)_n$, the total number of samples $N_j := v_j^((t))$ falling in $C_j^*$ follows $"Bin"(M, p_j)$ where $M = b t$. Conditional on $N_j = r >= 1$, the $r$ points $Z_1, dots, Z_r$ assigned to cluster $j$ are i.i.d. uniform draws from $C_j^*$ with mean $mu_j^*$ and $bb(E)[norm(Z_ell - mu_j^*)_2^2] = sigma_j^2$. By Lemma 1, $mu_j^((t)) = r^(-1) sum_(ell=1)^r Z_ell$, so independence of cross-terms yields the exact conditional MSE
  $ bb(E)[norm(mu_j^((t)) - mu_j^*)_2^2 | N_j = r] = (sigma_j^2) / r. $
  For every integer $r >= 1$, $1 / r <= 2 / (r + 1)$. Using the exact binomial identity for $N ~ "Bin"(M, p)$ with $q = 1 - p$,
  $ bb(E)[1 / (N + 1)] &= sum_(r=0)^M 1 / (r + 1) binom(M, r) p^r q^(M - r) \
                       &= 1 / ((M + 1) p) sum_(r=0)^M binom(M + 1, r + 1) p^(r+1) q^(M - r) \
                       &= (1 - q^(M+1)) / ((M + 1) p). $
  Subtracting the $r = 0$ atom $q^M$ and conditioning on $N >= 1$ gives
  $ bb(E)[1 / (N + 1) | N >= 1] = (1 - q^(M+1) - (M + 1) p q^M) / ((M + 1) p (1 - q^M)). $
  Observe that $q^M - q^(M+1) - (M + 1) p q^M = p q^M - (M + 1) p q^M = -M p q^M <= 0$, so the numerator is bounded above by $1 - q^M$, canceling the denominator and yielding $bb(E)[(N + 1)^(-1) | N >= 1] <= ((M + 1) p)^(-1)$. Multiplying by $2 sigma_j^2$ proves Claim 1.

  For Claim 2, because $bold(mu)^((t)) in cal(U)$ shares the partition $cal(C)^*$, the parallel-axis identity gives
  $ F(bold(mu)^((t))) - F(bold(mu)^*) = 1 / 2 sum_(j=1)^k p_j norm(mu_j^((t)) - mu_j^*)_2^2. $
  Since $bold(1)_(E_t) <= bold(1)_(N_j >= 1)$ pointwise for every $j$ and each squared norm is non-negative,
  $ &bb(E)[(F(bold(mu)^((t))) - F(bold(mu)^*)) bold(1)_(E_t)] \
    &quad <= 1 / 2 sum_(j=1)^k p_j bb(E)[norm(mu_j^((t)) - mu_j^*)_2^2 bold(1)_(N_j >= 1)] \
    &quad <= 1 / 2 sum_(j=1)^k p_j dot 2 sigma_j^2 bb(E)[1 / (N_j + 1)] \
    &quad <= 1 / 2 sum_(j=1)^k p_j dot (2 sigma_j^2) / ((b t + 1) p_j) \
    &quad = 1 / (b t + 1) sum_(j=1)^k sigma_j^2, $
  where the cluster proportions $p_j$ cancel identically. Dividing by $bb(P)(E_t) >= 1 - sum_(j=1)^k (1 - p_j)^(b t)$ establishes the conditional bound.
]

#theorem(title: "Almost-Sure Convergence via Robbins–Monro Decay", origin: "Classical — Cited")[
  For i.i.d. sampling from a bounded-support distribution $P$ with isolated centroidal Voronoi equilibria, the per-center step sizes $eta_(j,t) = 1 / v_j^((t))$ satisfy $sum_(t=1)^infinity eta_(j,t) = infinity$ and $sum_(t=1)^infinity eta_(j,t)^2 < infinity$ almost surely for every persistently active cluster ($p_j > 0$), guaranteeing almost-sure convergence of $bold(mu)^((t))$ to a stationary CVT @bottou1995convergence @sculley2010webscale.
]

#proof(label: "Proof Synthesis (Cited)")[
  For any cluster $j$ visited with asymptotic frequency $lim inf_(t -> infinity) v_j^((t)) / (b t) >= p_"min" > 0$, there exist almost surely finite random constants $0 < c_1 <= c_2$ such that $c_1 t <= v_j^((t)) <= c_2 t$ for all sufficiently large $t$. Consequently,
  $ sum_(t=1)^infinity 1 / (v_j^((t))) &>= 1 / c_2 sum_(t=t_0)^infinity 1 / t = infinity, \
    sum_(t=1)^infinity 1 / ((v_j^((t)))^2) &<= 1 / (c_1^2) sum_(t=t_0)^infinity 1 / (t^2) < infinity. $
  Verifying bounded stochastic gradient variance on the compact support of $P$, Bottou and Bengio @bottou1995convergence establish via supermartingale Lyapunov convergence on the piecewise-smooth risk $F(bold(mu))$ that $nabla F(bold(mu)^((t))) -> bold(0)$ almost surely, matching Sculley's empirical stabilization guarantee @sculley2010webscale.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Each iteration of Algorithm 1 evaluates $b k$ distances in $bb(R)^d$, costing $Theta(b k d)$ arithmetic operations—a factor $n / b$ reduction relative to Lloyd's $Theta(n k d)$. Crucially, working memory is $O((b + k) d)$, independent of total dataset size $n$: choosing $b in [256, 4096]$ keeps the active batch and centroid matrices resident inside CPU L2/L3 cache or GPU SRAM, enabling out-of-core streaming over terabyte corpora. Because Lemma 1 proves exact equivalence to vectorized slice averaging, high-performance implementations replace the inner loop (lines 5–7) with a single BLAS-3 distance matrix product $X_(M_t) M^T$ followed by scatter-add centroid updates.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for Mini-Batch K-Means],
    [Time / Space Complexity], [$O(b k d T)$ time; $O((b+k)d)$ streaming space (independent of $n$!).],
    [Convergence Behavior], [Non-monotonic $F(bold(mu))$ (Prop. 1); $O(1/(b t))$ local MSE (Theorem 2); a.s. CVT limit (Theorem 3).],
    [Geometric Inductive Bias], [Convex isotropic Euclidean Voronoi cells (identical to Lloyd's limit).],
    [Batch-Size Trade-Off], [Small $b$: fast steps, noisy centroids; large $b$: lower variance, approaches Lloyd.],
    [Outlier Breakdown Point], [$epsilon^* = 1/n$ (inherited from squared Euclidean mean; damped over time by $1/v_j$).],
    [Recommended Upgrades], [K-Means++ (Alg. 2) pilot seeding; Elkan/Hamerly (Algs. 13–14) for exact batch; BIRCH @zhang1996birch for trees.],
  ),
  caption: [Summary of theoretical guarantees and operational profile of Mini-Batch K-Means.],
)

#sec("IV", "Engineering and Scientific Applications")

Mini-Batch K-Means was engineered specifically for web-scale regimes ($n in [10^6, 10^9]$) where full-batch Lloyd passes saturate disk I/O and memory bandwidth, trading a slight ($1\%$–$3\%$) steady-state inertia gap for $10^2$–$10^3 times$ wall-clock speedups @sculley2010webscale.

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Web-Scale Document & Image Bag-of-Features Codebooks:* Originally validated by Sculley @sculley2010webscale on sparse TF-IDF web corpora and millions of SIFT descriptors, where $O(b k d)$ mini-batch iterations converge before Lloyd completes two full epochs while fitting sparse batches in cache.
2. *Out-of-Core Embedding Index Training:* Training $k = 10^4$–$10^5$ coarse quantizer centroids on billion-scale neural embedding datasets streams mini-batches directly from SSD/NVMe without loading the $n times d$ matrix into RAM.
3. *High-Velocity Clickstream & Ad-Tech User Segmentation:* Continuous behavioral event streams can be partitioned in bounded $O((b+k)d)$ memory, while per-center counts $v_j$ naturally stabilize mature cohorts without starving niche audience segments (Theorem 1).
4. *Large-Scale Single-Cell Transcriptomics Preprocessing:* On atlas-scale datasets ($n > 10^6$ cells in PCA space), Mini-Batch K-Means rapidly constructs coarse metacells prior to graph-based Leiden/Louvain community refinement.
5. *Edge & IoT Sensor Quantization:* Resource-constrained gateways with strict SRAM ceilings deploy Algorithm 1 to maintain online codebooks of vibration spectra and telemetry vectors without retaining historical buffers.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Mini-Batch K-Means exhibits four characteristic failure modes: (a) *minority cluster starvation* when $b p_j << 1$, leaving rare clusters unvisited across early rounds and prone to reassignment (mitigate by enlarging $b >= 10 k$ or reinitializing idle centers); (b) *step-size lock-in under non-stationary concept drift*, since $1 / v_j -> 0$ freezes centroids permanently as $t -> infinity$ (mitigate via exponential count forgetting); (c) *high-dimensional variance inflation*, where large $sum_j sigma_j^2$ amplifies the $(b t + 1)^(-1)$ MSE floor of Theorem 2; and (d) *non-convex shapes and extreme outliers* inherited from $L_2$ centroids.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Why Mini-Batch Fits / Prerequisite], [Caution / Alternative],
    [Web / Text Corpora], [Cache-resident $O((b+k)d)$ sparse steps], [Use Spherical (Alg. 12) for cosine],
    [Billion-Vector ANN], [Out-of-core streaming without full RAM], [Seed with K-Means++ (Alg. 2)],
    [Streaming Telemetry], [Constant memory; adaptive $1/v_j$ rates], [Use CluStream @aggarwal2003clustream for drift],
    [High-D Noisy Data], [Fast baseline after PCA denoising], [Noise inflates MSE (Thm. 2); denoise first],
    [Small Clean Data], [Negligible gain when $n < 10^4$], [Use exact K-Means++ (Alg. 2) / Elkan (Alg. 13)],
  ),
  caption: [Application suitability and engineering matrix for Mini-Batch K-Means.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Mini-Batch K-Means sits on an exact algebraic bridge linking sequential online clustering to full-batch Lloyd refinement:

#proposition(title: "Exact Parametric Continuum from MacQueen to Lloyd K-Means", origin: "Original Unified Proof")[
  Parameterize Algorithm 1 by batch size $b in {1, dots, n}$ and a per-round count-retention factor $gamma in [0, 1]$ modifying line 6 to $v_j^((t)) = gamma v_j^((t-1)) + m_(j,t)$:
  1. Setting $(b, gamma) = (1, 1)$ with $v_j^((0)) = 0$ recovers *MacQueen's sequential online K-Means* (Alg. 21) @macqueen1967some identically.
  2. Setting $(b, gamma) = (n, 0)$ (sampling without replacement, so $M_t = cal(X)$) recovers *Lloyd's batch K-Means* (Alg. 1) @lloyd1982least identically.
  3. Setting $gamma in (0, 1)$ produces an exponentially weighted moving-average (EWMA) tracker with constant effective memory window $(1 - gamma)^(-1)$, resolving concept-drift lock-in while preserving order-invariance inside every mini-batch.
]

#proof[
  By Lemma 1 with retention factor $gamma$, every non-empty slice $M_(j,t)$ updates
  $ mu_j^((t)) = (gamma v_j^((t-1)) mu_j^((t-1)) + m_(j,t) overline(x)_(M_(j,t))) / (gamma v_j^((t-1)) + m_(j,t)). $
  (1) At $b = 1, gamma = 1$, exactly one cluster $j = c(x_i; bold(mu)^((t-1)))$ has $m_(j,t) = 1$, giving $mu_j^((t)) = (v_j^((t-1)) mu_j^((t-1)) + x_i) / (v_j^((t-1)) + 1)$, which is MacQueen's exact incremental mean rule. (2) At $b = n$ without replacement and $gamma = 0$, $M_(j,t) = C_j^((t))$ and the prior term vanishes identically, leaving $mu_j^((t)) = overline(x)_(C_j^((t)))$, i.e., Lloyd's exact centroid update. (3) For $gamma in (0, 1)$, unrolling the recursion weights round-$s$ samples by $gamma^(t - s)$, bounding the effective denominator by $b p_j / (1 - gamma)$ and preventing the learning rate $m_(j,t) / v_j^((t))$ from vanishing under distribution drift.
]

This continuum also links Mini-Batch K-Means to streaming summary structures such as BIRCH @zhang1996birch, BICO (Alg. 78) @fichtenberger2013bico, and StreamKM++ (Alg. 195) @ackermann2012streamkmpp, which compress data geometrically rather than subsampling uniformly.

On the ClusterOverview 4-dataset empirical benchmark, Mini-Batch K-Means records *Iris Silhouette = 0.447* and *Iris ARI = 0.594*, *Customers Silhouette = 0.365*, *Multi-Shape Silhouette = 0.611* (ARI N/A), and *High-D Subspace+Noise Silhouette = 0.179* (ARI N/A), with *mean runtime = 0.387 ms*. Notably, its stochastic perturbations yield strong Iris ARI (0.594) and higher Customers Silhouette (0.365 vs. 0.334 for K-Means++), yet on High-D Subspace+Noise its Silhouette drops sharply to 0.179 (vs. 0.303 for K-Means++), directly illustrating the noise-variance amplification $sum_j sigma_j^2 / (b t + 1)$ proven in Theorem 2. Likewise, on small benchmark datasets ($n <= 240$), batch-sampling overhead yields 0.387 ms vs. 0.229 ms for K-Means++, confirming that Mini-Batch's $n / b$ computational advantage takes hold strictly in the large-$n$ regime ($n >> b$).
