#import "../template/ieee.typ": *

#alg_entry(
  id: 47,
  name: "DBSCAN++ (Fast Subsampled Density Clustering)",
  category: "Cat. II: Density & Mode-Seeking",
  year: "2019",
  paradigm: "Subsampled Core-Point Density Clustering",
  complexity: $O(m n d)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Standard DBSCAN (Alg. 31) @ester1996dbscan requires evaluating an $epsilon$-range density query for *every* observation $x_i in cal(X)$, incurring worst-case $Theta(n^2 d)$ arithmetic complexity and degrading sharply when spatial indices ($k$-d trees, $R^*$-trees) lose pruning efficiency in moderate-to-high dimensions @schubert2017dbscanrevisited @beyer1999nearest. Introduced by Jang and Jiang @jang2019dbscanpp, *DBSCAN++* overcomes this quadratic bottleneck through a sharp geometric insight: to estimate connected components of an upper density level set $L_f (lambda) = {x in bb(R)^d : f(x) >= lambda}$, it is statistically and geometrically redundant to evaluate empirical density at all $n$ points—a carefully chosen subset $S subset.eq cal(X)$ of size $m = ceil(p n) << n$ suffices to cover the core manifold, identify core points, and reconstruct cluster connectivity.

Let $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$ be equipped with the Euclidean norm $norm(dot)_2$, and fix a subsampling budget $m in {1, dots, n}$ (sampling fraction $p = m / n in (0, 1]$), minimum point threshold $"minPts" in bb(N)$, density estimation radius $epsilon_d > 0$, and core connectivity radius $epsilon_c > 0$ (canonically $epsilon_d = epsilon_c = epsilon$, or decoupled with $epsilon_c >= epsilon_d$ to compensate for subsampling gaps). DBSCAN++ selects a representative subset $S subset.eq cal(X)$ of size $|S| = m$ via either *uniform random subsampling* or *greedy $k$-center initialization* (farthest-first traversal), which iteratively appends
$ s_(t+1) = arg max_(x in cal(X)) min_(1 <= j <= t) norm(x - s_j)_2, quad t = 1, dots, m - 1. $

#definition(title: "Subsampled Core Graph & Nearest-Core Partition")[
  Given a sampled subset $S subset.eq cal(X)$ of size $m$, the *sampled core set* at parameters $(epsilon_d, "minPts")$ is
  $ C_S := { s in S : |B(s, epsilon_d) inter cal(X)| >= "minPts" }, $
  where $B(s, r) := {x in bb(R)^d : norm(x - s)_2 <= r}$ (equivalently, $r_k (s) <= epsilon_d$ for $k = "minPts" - 1$ excluding $s$). The *core-point cluster graph* is the undirected geometric graph $G_C = (C_S, E_C)$ with edges $E_C := {{s, s'} subset.eq C_S : norm(s - s')_2 <= epsilon_c}$. Let $cal(K)_1, dots, cal(K)_K$ be the connected components of $G_C$. Each observation $x in cal(X)$ is assigned to the cluster of its *nearest core point* $s^*(x) := arg min_(s in C_S) norm(x - s)_2$ if $norm(x - s^*(x))_2 <= epsilon_c$, and declared noise ($c(x) = -1$) otherwise.
]

Crucially, density counts $|B(s, epsilon_d) inter cal(X)|$ are computed against the *full* dataset $cal(X)$—preserving full statistical power of sample size $n$—while only evaluating queries centered at the $m$ points in $S$ (Algorithm 1).

#algbox(
  title: "DBSCAN++ Subsampled Core-Point Clustering",
  inputs: [Dataset $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$, subsample size $m = ceil(p n)$, seeding mode (`uniform` or `k-center`), radii $epsilon_d, epsilon_c > 0$, threshold $"minPts" >= 1$.],
  outputs: [Cluster labels $c: cal(X) -> {1, dots, K} union {-1}$.],
)[
  1. *Subset Selection ($|S| = m$):*
  2. #h(1em) *If* `uniform`: draw $S subset.eq cal(X)$ of size $m$ uniformly without replacement.
  3. #h(1em) *If* `k-center`: initialize $s_1 in cal(X)$, set $d_i <- norm(x_i - s_1)_2$; for $t = 2, dots, m$, choose $s_t = arg max_(x_i in cal(X)) d_i$ and update $d_i <- min(d_i, norm(x_i - s_t)_2)$.
  4. *Subsampled Core Detection:* Compute distances $D_(s, x) = norm(s - x)_2$ for all $(s, x) in S times cal(X)$ (or issue $m$ spatial range queries on $cal(X)$); form
     $ C_S = { s in S : |{x in cal(X) : D_(s, x) <= epsilon_d}| >= "minPts" }. $
  5. *Core Graph Construction:* Construct $G_C = (C_S, E_C)$ with edge ${s, s'} in E_C$ iff $D_(s, s') <= epsilon_c$.
  6. *Connected Components:* Compute components $cal(K)_1, dots, cal(K)_K$ of $G_C$ via BFS/Union-Find; set $c(s) = j$ for all $s in cal(K)_j$.
  7. *Nearest-Core Assignment:* For each $x in cal(X) backslash C_S$:
  8. #h(1em) Compute $s^*(x) = arg min_(s in C_S) D_(s, x)$ (breaking ties by smallest core index).
  9. #h(1em) *If* $D_(s^*(x), x) <= epsilon_c$, set $c(x) = c(s^*(x))$; *else* set $c(x) = -1$ (noise).
  10. *Return* cluster assignments $c$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five formal theoretical results governing DBSCAN++: (A) the exact 2-approximation covering bound of $k$-center seeding and deterministic preservation of global DBSCAN core connectivity under $delta$-covering; (B) classical minimax-optimal Hausdorff level-set consistency at subquadratic sample size $m = o(n)$; (C) explicit finite-sample counterexamples proving minority-cluster miss probability under uniform seeding, core-chain fracture, and radius-inflation collapse; and (D) border-point permutation determinism, isometry covariance, and exact pairwise distance reduction.

#subsec("A", "K-Center Covering Radius & Core-Graph Preservation")

#lemma(title: "2-Approximation Covering Bound of K-Center Seeding", origin: "Original Proof")[
  Let $S_m = {s_1, dots, s_m} subset.eq cal(X)$ be generated by greedy $k$-center seeding (Step 3 of Algorithm 1), and let
  $ R^*(m) := min_(S subset.eq cal(X), space |S| = m) max_(x in cal(X)) "dist"(x, S) $
  be the optimal $m$-center covering radius of $cal(X)$. Then the achieved covering radius $R(S_m) := max_(x in cal(X)) "dist"(x, S_m)$ satisfies
  $ R(S_m) <= 2 R^*(m). $
  Moreover, if $cal(X) subset.eq [0, L]^d$, then $R(S_m) <= (sqrt(d) L) / floor(m^(1 / d))$.
]

#proof[
  Let $s_(m+1) := arg max_(x in cal(X)) "dist"(x, S_m)$ so that $R(S_m) = min_(1 <= j <= m) norm(s_(m+1) - s_j)_2$. Because each greedy step $t in {2, dots, m}$ chooses $s_t$ maximizing the minimum distance to ${s_1, dots, s_(t-1)}$, the greedy distances are non-increasing: for all $1 <= a < b <= m + 1$,
  $ norm(s_b - s_a)_2 >= min_(j < b) norm(s_b - s_j)_2 >= min_(j <= m) norm(s_(m+1) - s_j)_2 = R(S_m). $
  Now let $S^* = {c_1^*, dots, c_m^*} subset.eq cal(X)$ attain $R^*(m)$. Since $|S_m union {s_(m+1)}| = m + 1 > m$, by the pigeonhole principle two distinct points $s_a, s_b$ fall inside the same optimal ball $B(c_r^*, R^*(m))$. By the triangle inequality,
  $ R(S_m) <= norm(s_a - s_b)_2 <= norm(s_a - c_r^*)_2 + norm(c_r^* - s_b)_2 <= 2 R^*(m). $
  For $cal(X) subset.eq [0, L]^d$, set $q = floor(m^(1 / d)) >= 1$ and partition $[0, L]^d$ into $q^d <= m$ axis-aligned subcubes of side $L / q$ and Euclidean diameter $sqrt(d) L / q$. Selecting one dataset point from each non-empty subcube yields a subset $tilde(S) subset.eq cal(X)$ of size $<= m$ covering $cal(X)$ within $sqrt(d) L / q$, hence $R^*(m) <= sqrt(d) L / q$. Wait—even tighter, every point in a subcube lies within $sqrt(d) L / (2 q)$ of its geometric center, so two points in the same subcube lie within $sqrt(d) L / q$; substituting $R(S_m) <= 2 R^*(m)$? Even simpler: by the pigeonhole principle on the $q^d <= m$ subcubes, two of the $m + 1$ points in $S_m union {s_(m+1)}$ must lie in the *same* subcube of diameter $sqrt(d) L / q$, directly giving $R(S_m) <= norm(s_a - s_b)_2 <= (sqrt(d) L) / floor(m^(1 / d))$!
]

#theorem(title: "Deterministic Core Connectivity Preservation under Covering", origin: "Original Proof")[
  Let $C_"full" := {x in cal(X) : |B(x, epsilon_0) inter cal(X)| >= "minPts"}$ be the global DBSCAN core set at parameters $(epsilon_0, "minPts")$, and suppose the sample $S subset.eq cal(X)$ is a $delta$-cover of $C_"full"$ (i.e., $max_(x in C_"full") "dist"(x, S) <= delta$, guaranteed with $delta = R(S_m)$ under $k$-center seeding). Run DBSCAN++ on $S$ with inflated radii $epsilon_d >= epsilon_0 + delta$ and $epsilon_c >= epsilon_0 + 2 delta$. Then:
  1. *Zero Core Loss:* Every global core point $x in C_"full"$ is assigned to a non-noise cluster of DBSCAN++.
  2. *Full Path Preservation:* Any two global core points $x, y in C_"full"$ connected by an $epsilon_0$-chain in $C_"full"$ are assigned to the *exact same* DBSCAN++ cluster.
]

#proof[
  (1) Fix any $x in C_"full"$. Because $S$ is a $delta$-cover of $C_"full"$, there exists $s_x in S$ with $norm(x - s_x)_2 <= delta$. For every $z in B(x, epsilon_0)$, the triangle inequality gives
  $ norm(z - s_x)_2 <= norm(z - x)_2 + norm(x - s_x)_2 <= epsilon_0 + delta <= epsilon_d. $
  Hence $B(x, epsilon_0) subset.eq B(s_x, epsilon_d)$, which implies
  $ |B(s_x, epsilon_d) inter cal(X)| >= |B(x, epsilon_0) inter cal(X)| >= "minPts". $
  Therefore $s_x in C_S$. Consequently, $x$'s nearest core point $s^*(x) in C_S$ obeys $norm(x - s^*(x))_2 <= norm(x - s_x)_2 <= delta < epsilon_c$, so $x$ is guaranteed to receive a non-noise cluster label.

  (2) Suppose $x = u_0, u_1, dots, u_L = y$ is a path in $C_"full"$ with $norm(u_t - u_(t+1))_2 <= epsilon_0$ for $t = 0, dots, L - 1$. By part (1), each $u_t$ admits a witness core point $s_t in C_S$ with $norm(u_t - s_t)_2 <= delta$. Applying the triangle inequality along each step gives
  $ norm(s_t - s_(t+1))_2 &<= norm(s_t - u_t)_2 + norm(u_t - u_(t+1))_2 + norm(u_(t+1) - s_(t+1))_2 \
                          &<= delta + epsilon_0 + delta = epsilon_0 + 2 delta <= epsilon_c. $
  Thus ${s_t, s_(t+1)} in E_C$ for all $t$, placing $s_0$ and $s_L$ in the same connected component of $G_C$. Finally, for the actual nearest core point $s^*(x) in C_S$ assigned to $x$,
  $ norm(s^*(x) - s_0)_2 <= norm(s^*(x) - x)_2 + norm(x - s_0)_2 <= 2 delta <= epsilon_c, $
  so ${s^*(x), s_0} in E_C$ (and analogously ${s^*(y), s_L} in E_C$). Hence $c(x) = c(s^*(x)) = c(s_0) = c(s_L) = c(s^*(y)) = c(y)$.
]

#subsec("B", "Minimax-Optimal Level-Set Consistency & Subquadratic Runtime")

#theorem(title: "Statistical Consistency & Subquadratic Sample Scaling", origin: "Classical — Cited")[
  Let $X_1, dots, X_n$ be i.i.d. draws from a twice-continuously differentiable density $f$ on compact support $cal(M) subset.eq bb(R)^d$ whose upper level set $L_f (lambda) = union.big_(j=1)^K C_j$ consists of $K$ separated connected components with non-vanishing boundary gradient $norm(nabla f(x))_2 >= g_0 > 0$ on $partial L_f (lambda)$ @jang2019dbscanpp. Configure DBSCAN++ with $k$-NN density cutoff ($k = "minPts" asymp n^(4 / (d+4)) log n$, corresponding to $epsilon_d asymp (k / (n v_d lambda))^(1 / d)$), connectivity radius $epsilon_c asymp n^(-1 / (d+4)) (log n)^(1 / d)$, and $k$-center (or uniform) subsample size
  $ m = Theta(n^(d / (d+4)) log n). $
  Then with probability tending to $1$ as $n -> infinity$:
  1. DBSCAN++ recovers the exact number of components $K$ and achieves the minimax-optimal Hausdorff level-set error $d_H (hat(C)_j, C_j) = tilde(O)(n^(-1 / (d+4)))$ for every $j = 1, dots, K$.
  2. Exact brute-force computation runs in strictly subquadratic time
     $ O(m n d) = tilde(O)(d dot n^(2 - 4 / (d+4))) = o(n^2) $
     for every fixed dimension $d >= 1$ @jang2019dbscanpp.
]

#proof(label: "Proof Synthesis (Cited — Jang & Jiang 2019)")[
  By second-order Taylor expansion of $f$ and VC-dimension uniform ball concentration @chaudhuri2014consistent, the empirical ball density $hat(f)_(epsilon_d) (x)$ achieves uniform error $O(epsilon_d^2 + sqrt((d log n) / (n epsilon_d^d)))$ over $cal(M)$. Balancing squared bias $epsilon_d^2$ against stochastic variance $(n epsilon_d^d)^(-1/2)$ fixes $epsilon_d asymp n^(-1 / (d+4))$ (up to log factors), yielding density error $O(n^(-2 / (d+4)))$ and, via the non-degenerate gradient lower bound $norm(nabla f)_2 >= g_0 > 0$, Hausdorff boundary perturbation $O(n^(-1 / (d+4)))$ for the full core set. By Lemma 1, $k$-center seeding with $m = Theta(n^(d / (d+4)))$ attains covering radius $delta = R(S_m) = O(m^(-1/d)) = O(n^(-1 / (d+4)))$, matching the statistical boundary uncertainty without inflating its order! Applying Theorem 1 with $epsilon_c asymp delta$ guarantees each true component $C_j$ forms a single connected component in $G_C$, while fixed inter-component separation $Delta > 2 epsilon_c$ for large $n$ prevents cross-valley merges. Substituting $m = tilde(Theta)(n^(d / (d+4)))$ into $O(m n d)$ gives $tilde(O)(d n^(2 - 4/(d+4)))$ @jang2019dbscanpp.
]

#subsec("C", "Subsampling Failure Modes: Minority Erasure & Collapse")

While asymptotically optimal as $n -> infinity$, finite-sample subsampling introduces sharp failure regimes absent in full DBSCAN.

#theorem(title: "Minority-Cluster Erasure & Radius-Inflation Collapse", origin: "Original Counterexample")[
  Consider DBSCAN++ with subsample size $m = ceil(p n) < n$:
  1. *Uniform-Seeding Minority Erasure:* Suppose $cal(X)$ contains an isolated true cluster $C_"min"$ of size $n_0$ ($"minPts" <= n_0 < n - m$) separated from $cal(X) backslash C_"min"$ by distance $> epsilon_c$. Under uniform random subsampling without replacement, the exact probability that DBSCAN++ completely erases $C_"min"$ (labeling all $n_0$ points as noise) obeys
     $ (1 - m / (n - n_0 + 1))^(n_0) <= bb(P)(S inter C_"min" = emptyset) = binom(n - n_0, m) / binom(n, m) <= e^(-p n_0). $
     Consequently, missing $C_"min"$ occurs with constant probability $Omega(1)$ whenever $n_0 = O(1 / p)$.
  2. *Chain Fracture vs. Single-Cluster Collapse:* There exist 1D datasets where setting $epsilon_c = epsilon_d = epsilon$ fractures a connected DBSCAN cluster under subsampling, whereas inflating $epsilon_c = epsilon + 2 delta$ (Theorem 1) bridges an inter-cluster valley and collapses multiple clusters into $K = 1$.
]

#proof[
  (1) By counting subsets of size $m$ disjoint from $C_"min"$,
  $ bb(P)(S inter C_"min" = emptyset) = binom(n - n_0, m) / binom(n, m) = product_(j=0)^(n_0 - 1) (n - m - j) / (n - j) = product_(j=0)^(n_0 - 1) (1 - m / (n - j)). $
  For every $j in {0, dots, n_0 - 1}$, $m / n <= m / (n - j) <= m / (n - n_0 + 1)$, so
  $ (1 - m / (n - n_0 + 1))^(n_0) <= bb(P)(S inter C_"min" = emptyset) <= (1 - m / n)^(n_0) <= e^(-p n_0), $
  using $1 - u <= e^(-u)$ and $m / n >= p$. When $S inter C_"min" = emptyset$, no point of $C_"min"$ enters $C_S$; since $"dist"(C_"min", cal(X) backslash C_"min") > epsilon_c$, every $x in C_"min"$ satisfies $"dist"(x, C_S) > epsilon_c$ and is labeled noise. If $p = 0.05$ and $n_0 = 10$, then $bb(P) approx 0.60$.

  (2) *Fracture:* Place $2 N$ points along $bb(R)$ at $x_i = i (epsilon / 2)$ ($i = 1, dots, 2N$, $"minPts" = 3$). Every point is a global core point forming one DBSCAN cluster. If any two consecutive interior points $x_r, x_(r+1)$ are excluded from $S$, the gap between adjacent sampled points jumps to $x_(r+2) - x_(r-1) = 3 epsilon / 2 > epsilon$, fracturing the cluster into two disconnected pieces when $epsilon_c = epsilon$. *Collapse:* Conversely, take two dense intervals $[0, 1]$ and $[1 + w, 2 + w]$ separated by a valley of width $w = (3 / 2) epsilon$. Global DBSCAN ($epsilon_0 = epsilon$) cleanly separates the two clusters ($K = 2$). However, inflating $epsilon_c = epsilon + 2 delta$ to guarantee zero fracture for a subsample with covering radius $delta = epsilon / 2$ yields $epsilon_c = 2 epsilon > w$, bridging the valley and collapsing the entire dataset into a single cluster ($K = 1$).
]

#subsec("D", "Border Determinism, Invariance & Exact Distance Savings")

#proposition(title: "Order-Independence, Covariance & Distance Reduction", origin: "Original Proof")[
  For any fixed subsample $S subset.eq cal(X)$ of size $m = p n$ and radii $(epsilon_d, epsilon_c)$:
  1. *Border-Point Permutation Determinism:* Assuming unique distances to $C_S$, Algorithm 1's output labels are strictly *invariant* to arbitrary permutations of $cal(X) backslash S$ (eliminating standard DBSCAN's order-dependent border-point ambiguity @schubert2017dbscanrevisited).
  2. *Similarity Covariance:* Under any Euclidean similarity $T(x) = c Q x + b$ ($c > 0$, $Q in O(d)$), rescaling $(epsilon_d, epsilon_c) |-> (c epsilon_d, c epsilon_c)$ preserves the partition identically for fixed sample indices $S$.
  3. *Exact Distance-Matrix Reduction:* Computing the single $m times n$ cross-distance matrix $D_(S, cal(X))$ suffices for seeding ($k$-center), density counting, core-graph construction, and nearest-core assignment, evaluating exactly
     $ N_"dist" = m(n - m) + binom(m, 2) = (2 p - p^2) binom(n, 2) - O(m) $
     unique pairs—an exact speedup factor of $(2 p - p^2)^(-1) > 1 / (2 p)$ over symmetric all-pairs DBSCAN (and $1 / p$ over naive $n$-query DBSCAN) with $O(m n)$ memory instead of $O(n^2)$.
]

#proof[
  (1) The core set $C_S$ depends only on cardinality counts $|B(s, epsilon_d) inter cal(X)|$, which are invariant to point ordering; connected components of $G_C$ are graph-theoretic invariants; and each $x in cal(X) backslash C_S$ maps to $arg min_(s in C_S) norm(x - s)_2$, which depends solely on $x$ and $C_S$ rather than traversal order. (2) Because $norm(T(x) - T(y))_2 = c norm(x - y)_2$, all ball inclusions $norm(T(x) - T(s))_2 <= c epsilon$ and nearest-core minimizers commute identically with $T$. (3) Counting unique unordered pairs incident to $S$ gives $m(n - m)$ cross-pairs between $S$ and $cal(X) backslash S$ plus $m(m - 1) / 2$ internal pairs within $S$, totaling $m n - m(m + 1) / 2$. Dividing by $binom(n, 2) = n(n - 1) / 2$ with $m = p n$ yields $(2 p n^2 - p^2 n^2 - p n) / (n^2 - n) = 2 p - p^2 + O(1 / n)$. The single matrix $D_(S, cal(X))$ directly furnishes Steps 3, 4, 5 (via its $C_S times C_S$ submatrix), and 8 without extra distance evaluations.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Under exact Euclidean evaluation, DBSCAN++ requires $O(m n d)$ time and $O(m n)$ temporary space (or $O((n + m) d)$ streaming memory by processing $S$ in blocks), reducing runtime by an order of magnitude at $p in [0.05, 0.15]$. When paired with a $k$-d or ball tree in low dimensions, issuing only $m$ range queries reduces average-case cost to $O(m log n + n log m)$, while avoiding high-dimensional index breakdown because the $m times n$ distance matrix computes cleanly via dense BLAS-3 GEMM ($S X^T$). Note that $k$-center seeding adds only $O(m n d)$ work (zero asymptotic overhead) and substantially lowers variance over uniform seeding by preventing cluster gaps—though it can be drawn toward isolated extreme outliers unless paired with outlier-trimmed initialization @jang2019dbscanpp.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for DBSCAN++],
    [Time / Space Complexity], [$O(m n d)$ exact ($tilde(O)(d n^(2 - 4/(d+4)))$ optimal); $O(m n)$ or $O(n d)$ streaming space.],
    [Geometric Inductive Bias], [Arbitrary non-convex density level sets covered by $m$ representative core balls.],
    [Determinism / Seeding], [Border-deterministic given $S$; stochastic under `uniform`, seed-dependent under `k-center`.],
    [Invariances], [Covariant under Euclidean similarities $c Q x + b$; sensitive to anisotropic feature scaling.],
    [Outlier Robustness], [Noise points excluded via $norm(x - s^*)_2 > epsilon_c$; `k-center` wastes seeds on extreme outliers.],
    [Recommended Companions], [DBSCAN (Alg. 31) for small $n$; HDBSCAN (Alg. 32) for multi-density; Mean Shift++ (Alg. 55) for mode-seeking.],
  ),
  caption: [Summary of theoretical guarantees and operational properties of DBSCAN++.],
)

#sec("IV", "Engineering and Scientific Applications")

DBSCAN++ targets large-scale metric datasets where standard DBSCAN's $n$ neighborhood queries exceed latency or memory budgets, yet non-convex cluster shapes and noise rejection remain essential.

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Large-Scale LiDAR & 3D Point-Cloud Segmentation:* Autonomous-driving scans yield $10^5$–$10^6$ 3D returns per frame containing non-convex vehicles, pedestrians, and road barriers amid sparse dust noise. Because smooth surfaces have low intrinsic dimension $d_"int" approx 2$, $k$-center seeding with $p = 0.05$ covers surface manifolds with tiny radius $delta = O(m^(-1/2))$ (Lemma 1), achieving real-time $10 times$–$20 times$ speedup over DBSCAN.
2. *High-Volume Geospatial Mobility & Ride-Hailing Hotspots:* Millions of GPS pickup/dropoff coordinates form irregular road-network corridors and dense urban hubs. DBSCAN++ preserves arbitrary street geometries while querying only $m = O(n^(1/3))$ candidate cores, enabling interactive city-scale clustering.
3. *Large Embedding Deduplication & Anomaly Screening:* Dense neural embeddings ($d = 16$–$128$) break $k$-d trees @beyer1999nearest, forcing brute-force distance scans. Formulating DBSCAN++'s $m times n$ core query as a single GPU matrix multiplication $S X^T$ accelerates density clustering by $> 15 times$ over all-pairs DBSCAN.
4. *Astronomical Catalog & High-Energy Particle Clustering:* Sky surveys and calorimeter hits exhibit dense non-spherical filaments embedded in uniform background radiation. DBSCAN++'s full-$n$ density counting around subsampled cores preserves statistical signal-to-noise ratio while scaling to millions of events.
5. *Interactive Hyperparameter Sweeps for Density Discovery:* Because $D_(S, cal(X))$ is computed once in $O(m n d)$ time (Proposition 1(3)), sweeping dozens of $(epsilon_d, epsilon_c, "minPts")$ pairs requires only fast $O(m n)$ thresholding and $O(m^2)$ connected-component passes, accelerating model selection dramatically.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

DBSCAN++ exhibits four concrete failure modes: (a) *rare minority-cluster erasure* under uniform seeding when $n_0 < 1 / p$ (Theorem 3(1)); (b) *outlier-starved $k$-center budgets*, where $k$-center seeding spends many of its $m$ centers on isolated extreme outliers before covering dense interior clusters; (c) *over-merging collapse ($K = 1$)* on small or tightly spaced datasets when $epsilon_c$ exceeds inter-cluster valleys (Theorem 3(2)); and (d) *multi-density blind spots* inherited from fixed global $(epsilon_d, "minPts")$ (resolved by HDBSCAN Alg. 32 or RNN-DBSCAN Alg. 45 @bryant2018rnndbscan).

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Why DBSCAN++ Fits / Prerequisite], [Caution / Alternative],
    [3D LiDAR / Robotics], [Low intrinsic $d$; $10 times$+ speedup via $k$-center], [Pre-filter extreme far outliers],
    [GPS / Mobility ($n > 10^5$)], [Arbitrary road shapes; subquadratic queries], [Avoid `uniform` if rare zones matter],
    [Medium-D Embeddings], [BLAS-3 $S X^T$ beats broken $k$-d trees], [Use SNN (Alg. 36) if $d > 50$],
    [Small/Dense Tabular], [Subsampling gains minimal; collapse risk], [Use exact DBSCAN (Alg. 31) or GriDBSCAN (Alg. 48)],
    [Varying-Density Data], [Single global $epsilon_d$ limits resolution], [Use HDBSCAN (Alg. 32) or DeBaCl (Alg. 46)],
  ),
  caption: [Application suitability and engineering trade-off matrix for DBSCAN++.],
)

#sec("V", "Structural Relations and Empirical Behavior")

DBSCAN++ bridges classical density connectivity (DBSCAN, HDBSCAN) with metric covering coresets (K-Means++ Alg. 2 @arthur2007kmeanspp and $k$-center seeding) and grid-accelerated mode-seeking (Mean Shift++ Alg. 55 @jang2021meanshiftpp).

#proposition(title: "Exact Recovery of Core-Deterministic DBSCAN at p = 1", origin: "Original Unified Proof")[
  When $m = n$ ($p = 1$, so $S = cal(X)$) and $epsilon_d = epsilon_c = epsilon$, DBSCAN++ produces:
  1. The *exact* global core set $C_S = C_"full"$ and the *exact* core-point cluster equivalence classes of classical DBSCAN$(epsilon, "minPts")$ (Alg. 31) @ester1996dbscan; and
  2. A deterministic, permutation-invariant border assignment identical to DBSCAN\* @campello2013hdbscan augmented with nearest-core projection—eliminating classical DBSCAN's traversal-order non-determinism @schubert2017dbscanrevisited.
]

#proof[
  When $S = cal(X)$ and $epsilon_d = epsilon_c = epsilon$, $C_S = {x in cal(X) : |B(x, epsilon) inter cal(X)| >= "minPts"} = C_"full"$, and $G_C$ has an edge between $u, v in C_"full"$ iff $norm(u - v)_2 <= epsilon$. By Ester et al.'s density-connectivity lemma @ester1996dbscan, two core points are density-connected in DBSCAN if and only if they lie in the same connected component of $G_C$, establishing (1). Every non-core point $x in cal(X) backslash C_"full"$ has $norm(x - s^*(x))_2 <= epsilon$ iff $x$ lies in the $epsilon$-ball of at least one core point (i.e., $x$ is a DBSCAN border point), and assigning $x$ to $c(s^*(x))$ resolves multi-core ties purely geometrically rather than by scan order, proving (2).
]

While GriDBSCAN (Alg. 48) @mahran2008gridbscan and MR-DBSCAN (Alg. 49) @he2011mrdbscan accelerate DBSCAN via *spatial space partitioning* (suffering exponential $2^d$–$3^d$ boundary replication in high dimensions), DBSCAN++ accelerates via *data-manifold subsampling*, whose covering complexity depends on intrinsic rather than ambient grid faces. Its authors later extended the same subsampling-and-grid philosophy to mode-seeking in Mean Shift++ (Alg. 55) @jang2021meanshiftpp.

*Empirical Behavior on ClusterOverview Benchmarks:* On the 4-dataset ClusterOverview benchmark suite @rousseeuw1987silhouettes @hubert1985comparing, DBSCAN++ achieves the fastest runtime among the three algorithms in this group with a mean execution time of *1.346 ms*, empirically confirming Proposition 1(3)'s distance-reduction advantage even at modest sample sizes. On *Multi-Shape*, DBSCAN++ scores *Sil = 0.670* (ARI = N/A), cleanly capturing non-convex geometry while isolating background noise; on *High-D Subspace+Noise* ($d = 12$), it achieves a strong *Sil = 0.348* (ARI = N/A), substantially outperforming DeBaCl (`0.217`) because subsampled $epsilon$-balls resist full-assignment noise contamination. Conversely, on *Iris* (*Sil = 0.452*, *ARI = 0.360*) subsampling slightly fragments overlapping clusters, and on the compact *Customers* benchmark (*Sil = N/A*) DBSCAN++ collapses the dataset into a single cluster ($K = 1$), furnishing a textbook empirical illustration of Theorem 3(2)'s radius/subsampling over-merging failure mode on small, low-margin datasets.
