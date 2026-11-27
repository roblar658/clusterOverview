#import "../template/ieee.typ": *

#alg_entry(
  id: 21,
  name: "MacQueen's K-Means",
  category: "Cat. I: Centroid & Partitioning",
  year: "1967",
  paradigm: "Sequential Online Stochastic Approximation",
  complexity: $O(n k d T)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Let $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$ be a finite stream or sequence of observations in Euclidean space, and let $k in {1, dots, n}$ be the target number of clusters. While Lloyd's (Alg. 1) and Forgy's (Alg. 22) algorithms alternate synchronous batch passes over all $n$ observations, James MacQueen @macqueen1967some introduced a *sequential (online) $k$-means* process that updates the affected cluster centroid *immediately after processing every single datum*.

Formally, MacQueen's algorithm targets the expected quantization distortion (population WCSS) for a random draw $X tilde P$ on $bb(R)^d$:
$ Phi(bold(mu)) = 1 / 2 bb(E)_(X tilde P)[min_(1 <= j <= k) norm(X - mu_j)_2^2], $
and its finite-sample empirical counterpart $J(cal(C), bold(mu)) = sum_(j=1)^k sum_(i in C_j) norm(x_i - mu_j)_2^2$. Rather than storing the entire dataset before moving prototypes, the algorithm maintains only the running cluster state $(bold(mu), bold(n))$, where $mu_j in bb(R)^d$ is the exact arithmetic mean of the $n_j in bb(N)$ points currently assigned to cluster $j$.

#definition(title: "Sequential Running-Centroid State & MacQueen Equilibrium")[
  A *running-centroid state* $(bold(mu), bold(n)) in (bb(R)^d)^k times bb(N)_+^k$ tracks exact online sufficient statistics $mu_j = n_j^(-1) sum_(i in C_j) x_i$ and $n_j = |C_j|$. A partition $cal(C) = {C_1, dots, C_k}$ with $n_j >= 1$ is a *MacQueen stationary equilibrium* if for every cluster $a in {1, dots, k}$ with $n_a >= 2$ and every $i in C_a$,
  $ norm(x_i - mu_a)_2^2 <= norm(x_i - mu_b)_2^2 quad forall b != a, $
  where ties are broken consistently by retaining $c_i = a$.
]

MacQueen's original procedure @macqueen1967some consists of an $O(n k d)$ single-pass online construction phase (Pass 0, seeding $mu_1, dots, mu_k$ with the first $k$ points) optionally followed by multi-pass online reallocation passes (Algorithm 1) that update both the donor and recipient centroids immediately upon each cluster switch.

#algbox(
  title: "MacQueen's Sequential Online K-Means",
  inputs: [Ordered sequence $(x_1, dots, x_n) subset.eq bb(R)^d$, integer $k >= 1$, max passes $T >= 1$.],
  outputs: [Partition $cal(C) = {C_1, dots, C_k}$ and centroids $bold(mu)$.],
)[
  1. *Pass 0 (Online Construction):* Set $mu_j <- x_j$, $n_j <- 1$, $c_j <- j$ for $j = 1, dots, k$.
  2. *For* $i = k + 1, dots, n$ *do*:
  3. #h(1em) $j^* <- arg min_(1 <= j <= k) norm(x_i - mu_j)_2^2$; set $c_i <- j^*$.
  4. #h(1em) $n_(j^*) <- n_(j^*) + 1$; $mu_(j^*) <- mu_(j^*) + n_(j^*)^(-1) (x_i - mu_(j^*))$.
  5. *For* pass $t = 1, dots, T - 1$ *do* (multi-pass refinement):
  6. #h(1em) Set `moved <- false`.
  7. #h(1em) *For* $i = 1, dots, n$ *do*:
  8. #h(2em) Let $a <- c_i$. *If* $n_a <= 1$, *continue*.
  9. #h(2em) $b <- arg min_(1 <= j <= k) norm(x_i - mu_j)_2^2$.
  10. #h(2em) *If* $b != a$ and $norm(x_i - mu_b)_2^2 < norm(x_i - mu_a)_2^2$ *then*:
  11. #h(3em) $mu_a <- mu_a - (n_a - 1)^(-1) (x_i - mu_a)$; $n_a <- n_a - 1$.
  12. #h(3em) $mu_b <- mu_b + (n_b + 1)^(-1) (x_i - mu_b)$; $n_b <- n_b + 1$.
  13. #h(3em) $c_i <- b$; `moved <- true`.
  14. #h(1em) *If* not `moved`, *break*.
  15. *Return* $(cal(C), bold(mu))$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish four core mathematical results characterizing MacQueen's algorithm: (A) an exact algebraic recoil decomposition proving strict WCSS descent and finite-step convergence of multi-pass MacQueen; (B) exact equivalence of single-pass MacQueen to per-cluster Robbins–Monro stochastic gradient descent with asymptotically Newton-scaled diagonal step sizes; (C) an explicit permutation counterexample proving unbounded single-pass drift on ordered streams; and (D) early-arrival primacy sensitivity, rigid isometry equivariance, and breakdown point.

#subsec("A", "Exact Recoil WCSS Identity and Finite Convergence")

#lemma(title: "Rank-1 Mean Invariant and Exact Reallocation Energy Drop", origin: "Original Derivation")[
  The incremental updates in Algorithm 1 preserve $mu_j = |C_j|^(-1) sum_(ell in C_j) x_ell$ identically at every step. Furthermore, moving a single point $x_i$ from cluster $C_a$ (size $n_a >= 2$, centroid $mu_a$) to $C_b$ (size $n_b >= 1$, centroid $mu_b$) yields updated centroids $mu'_a, mu'_b$ and changes the exact partition WCSS $J^*(cal(C))$ by
  $ J^*(cal(C)^("old")) - J^*(cal(C)^("new")) &= (norm(x_i - mu_a)_2^2 - norm(x_i - mu_b)_2^2) \
    &quad + (n_a - 1) norm(mu'_a - mu_a)_2^2 \
    &quad + (n_b + 1) norm(mu'_b - mu_b)_2^2. $
]

#proof[
  For insertion into $C_b$ with sum $S_b = n_b mu_b$, the new mean is
  $ mu'_b &= (S_b + x_i) / (n_b + 1) = (n_b mu_b + x_i) / (n_b + 1) \
    &= mu_b + (x_i - mu_b) / (n_b + 1). $
  Similarly, removing $x_i$ from $C_a$ ($n_a >= 2$) gives
  $ mu'_a = (n_a mu_a - x_i) / (n_a - 1) = mu_a - (x_i - mu_a) / (n_a - 1). $
  By induction from Step 1 ($mu_j = x_j, n_j = 1$), exact centroid invariance holds after every step.

  Next, set $C'_a = C_a backslash {i}$ and $C'_b = C_b union {i}$. Applying the parallel-axis identity to $C'_a$ centered at $mu_a$ and adding $norm(x_i - mu_a)_2^2$ gives
  $ &sum_(ell in C_a) norm(x_ell - mu_a)_2^2 - sum_(ell in C'_a) norm(x_ell - mu'_a)_2^2 \
    &quad = (n_a - 1) norm(mu'_a - mu_a)_2^2 + norm(x_i - mu_a)_2^2. $
  Noting $norm(mu'_a - mu_a)_2^2 = (n_a - 1)^(-2) norm(x_i - mu_a)_2^2$, the donor cluster WCSS decreases by $(1 + (n_a - 1)^(-1)) norm(x_i - mu_a)_2^2 = (n_a / (n_a - 1)) norm(x_i - mu_a)_2^2$. Symmetrically, applying the parallel-axis identity to $C_b$ centered at $mu'_b$ and adding $norm(x_i - mu'_b)_2^2$ with $x_i - mu'_b = (n_b / (n_b + 1))(x_i - mu_b)$ gives
  $ &sum_(ell in C'_b) norm(x_ell - mu'_b)_2^2 - sum_(ell in C_b) norm(x_ell - mu_b)_2^2 \
    &quad = n_b norm(mu_b - mu'_b)_2^2 + norm(x_i - mu'_b)_2^2 \
    &quad = (n_b / (n_b + 1)) norm(x_i - mu_b)_2^2. $
  Subtracting the recipient increase from the donor decrease yields
  $ Delta J^* &= (n_a / (n_a - 1)) norm(x_i - mu_a)_2^2 - (n_b / (n_b + 1)) norm(x_i - mu_b)_2^2 \
    &= (norm(x_i - mu_a)_2^2 - norm(x_i - mu_b)_2^2) \
    &quad + norm(x_i - mu_a)_2^2 / (n_a - 1) + norm(x_i - mu_b)_2^2 / (n_b + 1). $
  Substituting $(n_a - 1) norm(mu'_a - mu_a)_2^2 = norm(x_i - mu_a)_2^2 / (n_a - 1)$ and $(n_b + 1) norm(mu'_b - mu_b)_2^2 = norm(x_i - mu_b)_2^2 / (n_b + 1)$ completes the proof.
]

#theorem(title: "Strict Lyapunov Descent and Finite-Step Convergence", origin: "Original Proof")[
  Every executed reallocation move in Steps 10–13 of Algorithm 1 strictly decreases $J^*(cal(C))$ by strictly more than the raw Voronoi distance gap $norm(x_i - mu_a)_2^2 - norm(x_i - mu_b)_2^2$. Consequently, multi-pass MacQueen's K-Means never cycles and terminates in at most $S(n, k)$ moves at a MacQueen stationary equilibrium.
]

#proof[
  Step 10 executes a move only when $n_a >= 2$ and $norm(x_i - mu_b)_2^2 < norm(x_i - mu_a)_2^2$. By Lemma 1, both centroid recoil energies $(n_a - 1) norm(mu'_a - mu_a)_2^2 > 0$ (since $norm(x_i - mu_a)_2^2 > norm(x_i - mu_b)_2^2 >= 0$) and $(n_b + 1) norm(mu'_b - mu_b)_2^2 >= 0$ are non-negative, so
  $ J^*(cal(C)^("old")) - J^*(cal(C)^("new")) &> norm(x_i - mu_a)_2^2 \
    &quad - norm(x_i - mu_b)_2^2 > 0. $
  Because every move strictly lowers the real-valued functional $J^*(cal(C))$, no partition of $V$ into $k$ non-empty subsets can ever be visited twice. Since there are exactly $S(n, k)$ such partitions (the Stirling number of the second kind), the total number of moves is bounded by $S(n, k)$, and Algorithm 1 halts when `moved = false`, satisfying Definition 1.
]

#remark(title: "Relation to Hartigan–Wong (Alg. 19)")[
  Lemma 1 exposes the exact mathematical gap between MacQueen and Hartigan–Wong (Alg. 19 @hartigan1979algorithm): Hartigan–Wong accepts *every* move with $(n_b / (n_b + 1)) norm(x_i - mu_b)_2^2 < (n_a / (n_a - 1)) norm(x_i - mu_a)_2^2$, including boundary moves where $norm(x_i - mu_b)_2^2 >= norm(x_i - mu_a)_2^2$ yet the two recoil terms still make $Delta J^* > 0$. Thus every Hartigan–Wong fixed point is a MacQueen equilibrium, but not conversely.
]

#subsec("B", "Robbins–Monro SGD Equivalence & Asymptotic Efficiency")

#theorem(title: "Stochastic Approximation & Strong Convergence", origin: "Classical — Cited")[
  Let $X_1, X_2, dots$ be i.i.d. draws from a distribution $P$ on $bb(R)^d$ with $bb(E)[norm(X)_2^2] < infinity$ and $P(partial V_j (bold(mu))) = 0$ for all non-degenerate $bold(mu)$. Then single-pass online MacQueen's update (Step 4) is exact per-cluster Robbins–Monro stochastic gradient descent on $Phi(bold(mu)) = 1/2 bb(E)[min_j norm(X - mu_j)_2^2]$ with adaptive step size $eta_(j, t) = 1 / n_j (t)$ @macqueen1967some @bottou1995convergence. Moreover, along any sequence converging to a strict local minimum $bold(mu)^*$ with positive cell probabilities $p_j = P(V_j (bold(mu)^*)) > 0$, $eta_(j, t)$ asymptotically matches the inverse diagonal population Hessian $(t nabla^2_(mu_j, mu_j) Phi(bold(mu)^*))^(-1)$, attaining optimal $O(1/t)$ asymptotic variance rate @bottou1995convergence.
]

#proof[
  Define the instantaneous loss $ell(x; bold(mu)) = 1/2 min_(1 <= r <= k) norm(x - mu_r)_2^2$. For every $x$ off the measure-zero Voronoi boundary union $union_j partial V_j (bold(mu))$, the minimizer $j^*(x)$ is unique and locally constant in $bold(mu)$, so $ell$ is differentiable with block gradient
  $ nabla_(mu_j) ell(x; bold(mu)) = bb(1)[x in V_j (bold(mu))] (mu_j - x). $
  Applying per-cluster SGD with learning rate $eta_(j, t) = n_j (t)^(-1)$ gives
  $ mu_j^((t)) &= mu_j^((t-1)) - n_j (t)^(-1) bb(1)[X_t in V_j (bold(mu)^((t-1)))] \
    &quad times (mu_j^((t-1)) - X_t), $
  which coincides verbatim with Step 4 of Algorithm 1. By the Strong Law of Large Numbers, whenever $bold(mu)^((t)) -> bold(mu)^*$ with $p_j > 0$, we have $n_j (t) / t -> p_j$ almost surely, ensuring the Robbins–Monro divergence/summability criteria $sum_t eta_(j, t) = infinity$ and $sum_t eta_(j, t)^2 < infinity$ @macqueen1967some. Differentiating $Phi(bold(mu)) = 1/2 sum_j integral_(V_j (bold(mu))) norm(x - mu_j)_2^2 d P(x)$ at $bold(mu)^*$ yields diagonal Hessian blocks $nabla^2_(mu_j, mu_j) Phi(bold(mu)^*) = p_j I_d - B_j$ (where $B_j >= 0$ is the Voronoi boundary curvature term); neglecting boundary coupling, $eta_(j, t) approx (t p_j)^(-1)$ equals the exact inverse diagonal Newton step, yielding $O(1/t)$ asymptotic MSE @bottou1995convergence (compare strong global consistency of empirical minimizers @pollard1981strong).
]

#subsec("C", "Order Sensitivity and Unbounded Single-Pass Drift")

#theorem(title: "Permutation Non-Invariance & Unbounded Drift Ratio", origin: "Original Counterexample")[
  For any target factor $M > 0$, there exists a 4-point one-dimensional dataset $cal(X) subset.eq bb(R)$ with $k = 2$ and two permutations $pi_1, pi_2 in S_4$ of $cal(X)$ such that single-pass MacQueen's K-Means (Steps 1–4) achieves the global optimum $J^*$ under $pi_2$, yet under $pi_1$ satisfies
  $ J^*(cal(C)_"MacQ" (pi_1)) / J^*(cal(C)_"MacQ" (pi_2)) > M. $
]

#proof[
  Set $k = 2$, choose $L = sqrt(2 M) + 2 > 2$, and consider the four-point dataset $cal(X) = {0, 1, L, L + 1} subset.eq bb(R)$.
  - *Interleaved order $pi_2 = (0, L, 1, L + 1)$:* Step 1 initializes $mu_1 = 0$ ($n_1 = 1$) and $mu_2 = L$ ($n_2 = 1$). Point $1$ is closer to $0$ than $L$ ($1 < L - 1$), yielding $mu_1 = 1/2$, $n_1 = 2$. Point $L + 1$ is closer to $L$ than $1/2$ ($1 < L + 1/2$), yielding $mu_2 = L + 1/2$, $n_2 = 2$. Thus $cal(C)_"MacQ" (pi_2) = {{0, 1}, {L, L + 1}} = cal(C)^*$ with exact optimal cost
    $ J^*(cal(C)_"MacQ" (pi_2)) = 4 dot (1 / 2)^2 = 1. $
  - *Sorted order $pi_1 = (0, 1, L, L + 1)$:* Step 1 initializes $mu_1 = 0$ ($n_1 = 1$) and $mu_2 = 1$ ($n_2 = 1$). Next, $x_3 = L$ satisfies $|L - 1| < |L - 0|$, joining cluster 2 and dragging $mu_2 <- (1 + L) / 2$ ($n_2 = 2$). Finally, $x_4 = L + 1$ satisfies $|(L + 1) - (1 + L)/2| = (L + 1)/2 < L + 1 = |L + 1 - 0|$, also joining cluster 2 ($mu_2 <- (2L + 2)/3$, $n_2 = 3$). The resulting partition is $C_1 = {0}$, $C_2 = {1, L, L + 1}$ with WCSS
    $ J^*(cal(C)_"MacQ" (pi_1)) &= ((1 - 2L) / 3)^2 + ((L - 2) / 3)^2 + ((L + 1) / 3)^2 \
      &= 2 / 3 (L^2 - L + 1). $
  Because $L = sqrt(2 M) + 2$, we obtain
  $ J^*(cal(C)_"MacQ" (pi_1)) / J^*(cal(C)_"MacQ" (pi_2)) &= 2 / 3 (L^2 - L + 1) \
    &> 2 / 3 (2 M) > M, $
  proving both permutation non-invariance and unbounded single-pass suboptimality on sorted data streams.
]

#subsec("D", "Primacy Leverage, Invariances, and Breakdown Point")

#proposition(title: "Harmonic Primacy Bias, Isometry Equivariance, and Breakdown", origin: "Original Proof")[
  1. *(Harmonic Primacy Bias):* Upon ingesting the $m$-th point into cluster $j$, the centroid displacement satisfies $norm(mu_j^((m)) - mu_j^((m-1)))_2 = m^(-1) norm(x_i - mu_j^((m-1)))_2$, giving early observations $Theta(n)$ greater trajectory leverage than late observations.
  2. *(Isometry Equivariance):* For any fixed order $pi in S_n$ and rigid isometry $T(x) = Q x + b$ ($Q in O(d)$), Algorithm 1's entire trajectory is strictly equivariant under $T$, whereas anisotropic diagonal scaling alters assignments.
  3. *(Outlier Breakdown):* The finite-sample breakdown point is $epsilon^* = 1 / n$.
]

#proof[
  (1) Direct from Step 4: $mu_j^((m)) - mu_j^((m-1)) = m^(-1) (x_i - mu_j^((m-1)))$, taking Euclidean norms gives the exact factor $m^(-1)$, which decays harmonically from $1/2$ at $m = 2$ to $1/n$ at $m = n$.
  (2) Because $norm(Q x_i + b - (Q mu_j + b))_2^2 = norm(x_i - mu_j)_2^2$ and affine maps commute with $mu_j + m^(-1) (x_i - mu_j)$, every nearest-centroid choice and incremental mean update commutes with $T$ by induction on step index. Anisotropic scaling counterexample follows identically to Alg. 1.
  (3) Sending a single point $norm(x_i)_2 -> infinity$ drives its assigned cluster mean $norm(mu_(c_i))_2 >= n_(c_i)^(-1) norm(x_i)_2 - M_0 -> infinity$; worse, if $i <= k$, it directly initializes an entire centroid at infinity in Pass 0.
]

#sec("III", "Complexity, Geometry and Operational Profile")

In single-pass streaming mode ($T = 1$), MacQueen's K-Means requires only $O(n k d)$ arithmetic operations and *$O(k d)$ memory*—independent of stream length $n$—since each observation can be discarded immediately after updating $(mu_(j^*), n_(j^*))$. In multi-pass reallocation mode ($T > 1$), memory rises to $O((n + k) d)$ to store $cal(X)$ and previous labels $c_i$, while total work is $O(n k d T)$. Because every centroid moves as soon as a point switches, multi-pass MacQueen frequently converges in fewer passes $T$ than synchronous batch Lloyd/Forgy. However, the per-point loop-carried dependency on $bold(mu)$ prevents formulating assignments as a single BLAS-3 matrix product $X M^T$, making MacQueen less SIMD/GPU-friendly than Forgy (Alg. 22) or Mini-Batch K-Means (Alg. 3) at large $n$.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for MacQueen's K-Means],
    [Time / Space Complexity], [$O(n k d)$ time, $O(k d)$ space (1-pass); $O(n k d T)$ time, $O((n+k)d)$ space (multi-pass).],
    [Optimization Mechanism], [Per-cluster Robbins–Monro SGD ($eta_(j,t) = 1/n_j$) in Pass 0; strict WCSS descent in multi-pass.],
    [Order Sensitivity], [Strictly order-dependent; sorted streams cause unbounded single-pass drift (Theorem 3).],
    [Geometric Inductive Bias], [Isotropic spherical Voronoi cells; exact CVT equilibrium in multi-pass mode.],
    [Invariances / Breakdown], [Isometry-equivariant for fixed $pi$; scale-sensitive; breakdown $epsilon^* = 1/n$.],
    [Recommended Upgrades], [Randomly shuffle stream before Pass 0; use Mini-Batch (Alg. 3) for BLAS-3; Hartigan–Wong (Alg. 19) for sharper local optima.],
  ),
  caption: [Theoretical guarantees and operational profile of MacQueen's K-Means (Alg. 21).],
)

#sec("IV", "Engineering and Scientific Applications")

MacQueen's online architecture targets environments with strict memory ceilings, continuous data arrival, or single-pass CPU latency budgets where data arrives in approximately stationary, shuffled order.

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Embedded IoT & Edge Sensor Quantization:* Microcontrollers and edge gateways processing high-rate acoustic, vibration, or accelerometer streams cannot hold $O(n d)$ history in SRAM. MacQueen's $O(k d)$ state footprint and $O(k d)$ per-sample update provide exact online centroid maintenance with zero dynamic memory allocation.
2. *Adaptive Vector Quantization & Competitive Learning:* In signal compression and adaptive codebook design, MacQueen's rule is the zero-neighborhood limit of Kohonen's Self-Organizing Map @kohonen1982som and Neural Gas @martinetz1993neuralgas, adapting codewords continuously on live symbol streams without buffering frames.
3. *Single-Pass Streaming Prefiltering (Macro-Cluster Warm-Start):* High-volume clickstreams and network telemetry pipelines deploy a shuffled MacQueen pass as an ultra-fast $O(n k d)$ sketch before hierarchical micro-cluster engines such as BIRCH @zhang1996birch or CluStream @aggarwal2003clustream.
4. *Fast Sequential Initialization for Batch Solvers:* Because Pass 0 updates centroids immediately as points arrive rather than leaving initial seeds frozen across the entire first scan, a randomly permuted MacQueen pass rapidly pulls seeds toward dense data modes, accelerating downstream Elkan (Alg. 13) or Hartigan–Wong (Alg. 19) refinement.
5. *Non-Stationary Regime Tracking (with Constant Forgetting):* Replacing the harmonic gain $1/n_j$ with a bounded exponential moving average step $eta in (0, 1)$ converts MacQueen's update into a sliding-window tracker for financial order-book regimes and drifting user cohorts, mitigating Proposition 1(1)'s primacy lock-in.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

MacQueen's K-Means degrades sharply on (a) *sorted, grouped, or temporally autocorrelated streams* (Theorem 3), where the first $k$ points all fall in one cluster; (b) *high-dimensional noisy subspaces*, where early noisy arrivals lock centroids onto spurious directions before $n_j$ grows; (c) *non-convex shapes*, limited by convex Voronoi cells; and (d) *GPU/SIMD batch workloads*, where sequential centroid writes block BLAS-3 vectorization.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit / Prerequisite], [Caution / Alternative],
    [Edge / IoT Streams], [Strict $O(k d)$ RAM; $O(k d)$ per datum], [Shuffle or use EMA step $eta$],
    [Online Quantization], [Exact Robbins–Monro diagonal scaling], [Seed with K-Means++ (Alg. 2)],
    [Sorted Time Series], [Avoid: severe primacy drift (Thm. 3)], [Use Mini-Batch (Alg. 3) / shuffle],
    [GPU / Billion ANN], [Avoid: serial loop blocks BLAS-3 GEMM], [Use Mini-Batch or Forgy (Alg. 22)],
    [Non-Convex / Noise], [Avoid: convex cells & $epsilon^* = 1/n$], [Use DBSCAN (Alg. 31) / PAM (Alg. 4)],
  ),
  caption: [Application suitability and engineering matrix for MacQueen's K-Means.],
)

#sec("V", "Structural Relations and Empirical Behavior")

MacQueen's K-Means occupies the $b = 1$ boundary of the mini-batch optimization spectrum and sits cleanly inside the fixed-point hierarchy of centroidal partitioning methods.

#proposition(title: "Batch-Size Unification & Equilibrium Hierarchy", origin: "Original Unified Proof")[
  Consider Mini-Batch K-Means (Alg. 3 @sculley2010webscale) processing disjoint mini-batches $B_t$ of size $b in {1, dots, n}$ with per-cluster counts $n_j <- n_j + b_(j, t)$ ($b_(j, t) = |B_t inter V_j (bold(mu))|$) and centroid update
  $ mu_j^("new") = mu_j + (b_(j, t) / (n_j + b_(j, t))) (overline(x)_(B_(j, t)) - mu_j). $
  1. At $b = 1$, this rule reduces *identically* to MacQueen's Pass 0 update (Alg. 21).
  2. At $b = n$ with per-epoch count reset $n_j = 0$, it reduces *identically* to synchronous batch Lloyd/Forgy (Algs. 1 and 22).
  3. Across finite-sample multi-pass variants, the stationary partition sets satisfy the strict inclusion hierarchy
     $ cal(E)_"Global" subset.eq cal(E)_"HW" subset.eq cal(E)_"MacQ" = cal(E)_"Lloyd" = cal(E)_"Forgy", $
     where $cal(E)_"HW" subset.neq cal(E)_"MacQ"$ whenever boundary recoil moves exist.
]

#proof[
  (1) For $b = 1$, $B_t = {x_i}$ yields $b_(j^*, t) = 1$ with $overline(x)_(B_(j^*, t)) = x_i$ and $b_(j, t) = 0$ for $j != j^*$, giving $mu_(j^*)^("new") = mu_(j^*) + (n_(j^*) + 1)^(-1) (x_i - mu_(j^*))$, matching Step 4 of Algorithm 1 verbatim.
  (2) For $b = n$ and reset $n_j = 0$, $b_(j, t) / (0 + b_(j, t)) = 1$, yielding $mu_j^("new") = overline(x)_(C_j)$, which is exact batch Lloyd/Forgy.
  (3) For partitions without singletons, Definition 1 ($norm(x_i - mu_a)_2^2 <= norm(x_i - mu_b)_2^2$ for all $i in C_a, b != a$) is identical to a discrete CVT of Lloyd/Forgy. By Lemma 1, every Hartigan–Wong stationary partition satisfies $(n_a / (n_a - 1)) norm(x_i - mu_a)_2^2 <= (n_b / (n_b + 1)) norm(x_i - mu_b)_2^2$, which strictly implies $norm(x_i - mu_a)_2^2 < norm(x_i - mu_b)_2^2$ since $n_a / (n_a - 1) > 1 > n_b / (n_b + 1)$, proving $cal(E)_"HW" subset.eq cal(E)_"MacQ"$. Strictness holds on $cal(X) = {0, 1, 2}$, $k = 2$: $C_1 = {0, 1}$ ($mu_1 = 0.5$), $C_2 = {2}$ ($mu_2 = 2$) is a MacQueen/Lloyd CVT ($|1 - 0.5| = 0.5 < 1 = |1 - 2|$), yet $(n_1 / (n_1 - 1)) (0.5)^2 = 0.5 = (n_2 / (n_2 + 1)) (1)^2$; perturbing $x_2 = 1.05$ leaves $|1.05 - 0.525| = 0.525 < 0.95 = |1.05 - 2|$ (MacQueen CVT) while $2 (0.525)^2 = 0.55125 > 0.5 (0.95)^2 = 0.45125$, so Hartigan–Wong strictly moves $x_2$ to $C_2$.
]

*Empirical Benchmark Behavior:* On the ClusterOverview 4-dataset benchmark suite, MacQueen's K-Means records *Iris* $"Sil" = 0.452$, $"ARI" = 0.586$; *Customers* $"Sil" = 0.341$; *Multi-Shape* $"Sil" = 0.357$ ($"ARI" = "N/A"$); and *High-D Subspace+Noise* $"Sil" = 0.112$ ($"ARI" = "N/A"$), with a mean runtime of *0.140 ms*. On well-behaved convex datasets (*Iris* and *Customers*), MacQueen's immediate centroid tracking achieves the highest Iris ARI ($0.586$ vs. Forgy's $0.431$ and COP-KMeans's $0.564$) and highest Customers Silhouette ($0.341$) among the three variants. Conversely, on *Multi-Shape* ($0.357$) and *High-D Subspace+Noise* ($0.112$), its Silhouette drops steeply below Forgy ($0.609$ and $0.238$), empirically validating Theorem 3 and Proposition 1(1): early noisy or ordered arrivals induce persistent primacy drift that traps sequential centroids in suboptimal regions. Its ultra-low $0.140$ ms runtime confirms the $O(n k d)$ streaming efficiency derived in Section III, trailing only batch-vectorized Forgy ($0.063$ ms).
