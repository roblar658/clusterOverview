#import "../template/ieee.typ": *

#alg_entry(
  id: 18,
  name: "ISODATA",
  category: "Cat. I: Centroid & Partitioning",
  year: "1965",
  paradigm: "Dynamic Split-and-Merge Centroid Partitioning",
  complexity: $O((n k + k^2) d T)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Classical partitioning methods such as Forgy's and Lloyd's $k$-means (Alg. 1) @forgy1965cluster @lloyd1982least require the cluster cardinality $k$ to be fixed a priori. *ISODATA* (*Iterative Self-Organizing Data Analysis Technique*), introduced by Ball and Hall @ball1965isodata in 1965, generalizes prototype clustering to a *variable-cardinality* dynamical system over $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$ by interleaving Voronoi–centroid updates with three heuristic topology operators: *pruning* sparse clusters, *splitting* high-dispersion clusters, and *merging* closely spaced centroids.

ISODATA is governed by a six-parameter control vector $Theta = (k_0, theta_N, theta_S, theta_C, L, T)$ together with a split-offset factor $gamma in (0, 1]$ (classically $gamma = 1$):
- $k_0 >= 1$: target nominal cluster count defining the forced-split regime ($k_t <= k_0 / 2$) and forced-merge regime ($k_t >= 2 k_0$);
- $theta_N in bb(N)$: minimum admissible cluster size (pruning threshold);
- $theta_S > 0$: maximum admissible coordinate-wise standard deviation (splitting threshold);
- $theta_C > 0$: minimum admissible inter-centroid Euclidean distance (merging threshold);
- $L in bb(N)$: maximum number of centroid pairs merged per iteration;
- $T in bb(N)$: hard iteration budget.

At iteration $t$ with $k_t$ non-empty clusters $cal(C)^((t)) = {C_1, dots, C_(k_t)}$ and centroids $mu_j = |C_j|^(-1) sum_(i in C_j) x_i$, ISODATA computes three diagnostic statistics. First, the per-cluster and global mean Euclidean radii:
$ overline(D)_j &:= 1 / (|C_j|) sum_(i in C_j) norm(x_i - mu_j)_2, \
  overline(D) &:= 1 / (n_"act") sum_(j=1)^(k_t) |C_j| overline(D)_j, $
where $n_"act" = sum_(j=1)^(k_t) |C_j|$ counts retained active points. Second, the coordinate-wise standard deviation vector $sigma_j in bb(R)_+^d$ with components (for $r in {1, dots, d}$)
$ sigma_(j, r) := sqrt(1 / (|C_j|) sum_(i in C_j) (x_(i, r) - mu_(j, r))^2), $
principal coordinate $r_j^* := arg max_(1 <= r <= d) sigma_(j, r)$, and maximum spread $sigma_(j, max) := sigma_(j, r_j^*)$. Third, the pairwise centroid distance matrix $D_(j ell) := norm(mu_j - mu_ell)_2$ for $1 <= j < ell <= k_t$.

#definition(title: "ISODATA Structural Equilibrium")[
  A configuration $(cal(C), bold(mu), k)$ on active subset $cal(X)_"act" subset.eq cal(X)$ is an *ISODATA structural equilibrium* under $Theta$ if: (i) $(cal(C), bold(mu))$ is a discrete centroidal Voronoi tessellation (CVT) on $cal(X)_"act"$; (ii) every cluster satisfies $|C_j| >= theta_N$; (iii) no cluster triggers a split (either $sigma_(j, max) <= theta_S$, or both $k > k_0 / 2$ and $(overline(D)_j <= overline(D) or |C_j| <= 2(theta_N + 1))$); and (iv) no pair satisfies $D_(j ell) < theta_C$ whenever merging is active.
]

To prevent simultaneous conflicting topology edits in one pass, Ball and Hall @ball1965isodata schedule splitting on odd iterations (or whenever $k_t <= k_0 / 2$) and merging on even iterations (or whenever $k_t >= 2 k_0$), as formalized in Algorithm 1.

#algbox(
  title: "Ball–Hall ISODATA Algorithm",
  inputs: [Dataset $cal(X) subset.eq bb(R)^d$, parameters $Theta = (k_0, theta_N, theta_S, theta_C, L, T)$, offset $gamma in (0, 1]$, initial centroids $bold(mu)^((1))$ of size $k_1$.],
  outputs: [Dynamic partition $(cal(C), bold(mu), k)$ with variable $k$.],
)[
  1. *For* $t = 1, dots, T$ *do*:
  2. #h(1em) *Voronoi assignment:* Assign each $x_i in cal(X)$ to $c_i = arg min_(1 <= j <= k_t) norm(x_i - mu_j)_2^2$ to form $C_1, dots, C_(k_t)$.
  3. #h(1em) *Cardinality pruning:* Discard every $C_j$ with $|C_j| < theta_N$ and decrement $k_t$ accordingly.
  4. #h(1em) *Centroid & radius update:* Recompute $mu_j = |C_j|^(-1) sum_(i in C_j) x_i$, $overline(D)_j$, and global average $overline(D)$.
  5. #h(1em) *Branch selection:* If $t = T$, set $theta_C = 0$ and *break*. Go to *Split* if $k_t <= k_0 / 2$, or if ($t$ is odd and $k_t < 2 k_0$); else go to *Merge*.
  6. #h(1em) *Split step:* For each $C_j$, compute $sigma_(j, r)$, $r_j^* = arg max_r sigma_(j, r)$, and $sigma_(j, max) = sigma_(j, r_j^*)$. If $sigma_(j, max) > theta_S$ and (($overline(D)_j > overline(D)$ and $|C_j| > 2(theta_N + 1)$) or $k_t <= k_0 / 2$), replace $mu_j$ by two centroids
     $ mu_j^plus.minus = mu_j plus.minus gamma sigma_(j, max) e_(r_j^*), $
     increment $k_t <- k_t + 1$, and *continue* to next $t$.
  7. #h(1em) *Merge step:* Compute $D_(j ell) = norm(mu_j - mu_ell)_2$ for $j < ell$. Sort all pairs with $D_(j ell) < theta_C$ in ascending order. Greedily merge up to $L$ vertex-disjoint pairs $(j, ell)$ into
     $ mu_(j ell) = (|C_j| mu_j + |C_ell| mu_ell) / (|C_j| + |C_ell|), $
     decrementing $k_t$ per merge.
  8. *Return* final partition $(cal(C), bold(mu), k_T)$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

Unlike fixed-$k$ Lloyd's algorithm, ISODATA's split and merge operators act in opposite directions on the within-cluster sum of squares $J$. We rigorously analyze: (1) exact WCSS variations under weighted merging and axis-aligned splitting (uncovering an optimal split-offset formula); (2) an explicit $8$-point counterexample proving perpetual period-$2$ oscillation; (3) a hysteresis separation condition preventing immediate re-merging and ensuring finite convergence once edits cease; (4) exact proof of rotational non-invariance by a factor up to $sqrt(d)$ and partial outlier pruning; and (5) classical statistical contrast with consistent model-order criteria.

#subsec("A", "Exact Energy Variation Under Merge and Axis-Aligned Split")

#lemma(title: "Exact WCSS Identities for Merge and Axis Split", origin: "Original Derivation")[
  Let $C_j, C_ell subset.eq cal(X)$ be disjoint non-empty clusters with sizes $n_j, n_ell$ and centroids $mu_j, mu_ell$.
  1. *Merge Penalty:* Merging $C_j$ and $C_ell$ into $C_(j ell) = C_j union C_ell$ at barycenter $mu_(j ell)$ increases WCSS by the exact amount
     $ Delta J_"merge" (j, ell) &= (n_j n_ell) / (n_j + n_ell) norm(mu_j - mu_ell)_2^2 \
       &<= (n_j n_ell) / (n_j + n_ell) theta_C^2. $
  2. *Pre-Update Split Change:* Splitting $C_j$ into $mu_j^plus.minus = mu_j plus.minus gamma sigma_(j, max) e_(r^*)$ ($gamma > 0$) and partitioning $C_j$ by proximity to $mu_j^plus.minus$ changes its WCSS by
     $ Delta J_"split"^((0)) (j) = n_j (2 gamma sigma_(j, max) nu_1 - gamma^2 sigma_(j, max)^2), $
     where $nu_1 := n_j^(-1) sum_(i in C_j) |x_(i, r^*) - mu_(j, r^*)|$ is the mean absolute deviation along axis $r^*$. This drop is uniquely maximized at $gamma^* = nu_1 / sigma_(j, max) in (0, 1]$ and can be *negative* if $gamma > 2 nu_1 / sigma_(j, max)$.
  3. *Post-Update Split Guarantee:* Regardless of $gamma > 0$, after updating the two half-cluster means $mu_+$ and $mu_-$ on $C_j^+ = {i in C_j : x_(i, r^*) >= mu_(j, r^*)}$ and $C_j^- = C_j backslash C_j^+$, the WCSS of $C_j$ *strictly decreases* whenever $sigma_(j, max) > 0$ by
     $ Delta J_"split"^"upd" (j) &= (n_+ n_-) / n_j norm(mu_+ - mu_-)_2^2 \
       &>= (n_j^3 nu_1^2) / (4 n_+ n_-) >= n_j nu_1^2 > 0. $
]

#proof[
  (1) Applying the parallel-axis identity to $C_j$ and $C_ell$ centered at $mu_(j ell)$ gives
  $ Delta J_"merge" (j, ell) &= n_j norm(mu_j - mu_(j ell))_2^2 \
    &quad + n_ell norm(mu_ell - mu_(j ell))_2^2. $
  Since $mu_j - mu_(j ell) = n_ell / (n_j + n_ell) (mu_j - mu_ell)$ and $mu_ell - mu_(j ell) = - n_j / (n_j + n_ell) (mu_j - mu_ell)$, substituting yields
  $ Delta J_"merge" (j, ell) &= (n_j n_ell^2 + n_ell n_j^2) / (n_j + n_ell)^2 norm(mu_j - mu_ell)_2^2 \
    &= (n_j n_ell) / (n_j + n_ell) norm(mu_j - mu_ell)_2^2, $
  and the merge condition $norm(mu_j - mu_ell)_2 < theta_C$ yields the upper bound.

  (2) Decompose each $x_i - mu_j$ into its scalar coordinate $z_i := x_(i, r^*) - mu_(j, r^*)$ along $e_(r^*)$ and orthogonal residual $u_i perp e_(r^*)$, so $norm(x_i - mu_j)_2^2 = norm(u_i)_2^2 + z_i^2$. With $delta := gamma sigma_(j, max) > 0$, point $x_i$ is closer to $mu_j^+$ iff $z_i >= 0$ and to $mu_j^-$ iff $z_i < 0$, yielding squared distance $norm(u_i)_2^2 + (|z_i| - delta)^2$. Subtracting from $norm(u_i)_2^2 + z_i^2$ and summing over $i in C_j$ gives
  $ Delta J_"split"^((0)) (j) &= sum_(i in C_j) (2 delta |z_i| - delta^2) \
    &= n_j (2 delta nu_1 - delta^2). $
  This concave quadratic in $delta$ achieves its unique maximum at $delta^* = nu_1$, i.e., $gamma^* = nu_1 / sigma_(j, max) <= 1$ (by Cauchy–Schwarz, $nu_1 <= sigma_(j, max)$), and turns negative whenever $delta > 2 nu_1$.

  (3) Because $sum_(i in C_j) z_i = 0$ and $sigma_(j, max) > 0$, both half-clusters $C_j^+$ and $C_j^-$ are non-empty ($n_+, n_- >= 1$, $n_+ + n_- = n_j$) and satisfy
  $ sum_(i in C_j^+) z_i &= sum_(i in C_j^-) (-z_i) \
    &= 1 / 2 sum_(i in C_j) |z_i| = (n_j nu_1) / 2. $
  Consequently, the $r^*$-coordinate separation between updated means is
  $ mu_(+, r^*) - mu_(-, r^*) &= (n_j nu_1) / 2 (1 / n_+ + 1 / n_-) \
    &= (n_j^2 nu_1) / (2 n_+ n_-). $
  Applying Claim 1 in reverse to the partition $C_j = C_j^+ union C_j^-$ and dropping the non-negative orthogonal components of $norm(mu_+ - mu_-)_2^2$ gives
  $ Delta J_"split"^"upd" (j) &= (n_+ n_-) / n_j norm(mu_+ - mu_-)_2^2 \
    &>= (n_+ n_-) / n_j (mu_(+, r^*) - mu_(-, r^*))^2 \
    &= (n_j^3 nu_1^2) / (4 n_+ n_-) >= n_j nu_1^2 > 0, $
  where the final step uses the AM–GM bound $4 n_+ n_- <= (n_+ + n_-)^2 = n_j^2$.
]

#subsec("B", "Non-Monotonicity and Limit-Cycle Oscillations")

Because merging increases WCSS while splitting decreases it, ISODATA possesses *no* global Lyapunov function and can cycle perpetually unless capped by $T$.

#theorem(title: "Existence of Infinite Period-2 Split–Merge Limit Cycles", origin: "Original Counterexample")[
  There exists a strictly distinct $8$-point dataset $cal(X) subset.eq bb(R)$ and an admissible parameter tuple $Theta = (k_0, theta_N, theta_S, theta_C, L, T)$ with $gamma = 1$ on which Algorithm 1 enters an infinite period-$2$ limit cycle oscillating between $k = 1$ and $k = 2$ for every $T >= 2$.
]

#proof[
  Define $cal(X) = {-1.1, -1.05, -0.95, -0.9, 0.9, 0.95, 1.05, 1.1} subset.eq bb(R)$, and set $k_0 = 2$, $theta_N = 2$, $theta_S = 0.5$, $theta_C = 2.5$, $L = 1$, $gamma = 1$, initialized at $k_1 = 1$, $mu_1^((1)) = 0$.
  - *Odd iteration $t = 2 m - 1$ (from $k = 1, mu_1 = 0$):* All $8$ points assign to $C_1$, satisfying $|C_1| = 8 >= theta_N = 2$ and $|C_1| > 2(theta_N + 1) = 6$. The centroid remains $mu_1 = 0$ with variance
    $ sigma_1^2 = 1 / 4 (1.1^2 + 1.05^2 + 0.95^2 + 0.9^2) = 1.00625, $
    so $sigma_1 approx 1.0031 > theta_S = 0.5$. Because $k = 1 <= k_0 / 2$ (and $t$ is odd), Algorithm 1 enters *Split* and replaces $mu_1 = 0$ with $mu^plus.minus = plus.minus sigma_1 approx plus.minus 1.0031$, giving $k = 2$.
  - *Even iteration $t = 2 m$ (from $k = 2, mu = (-sigma_1, +sigma_1)$):* Voronoi assignment places the four negative points in $C_1$ and four positive points in $C_2$. Both satisfy $|C_j| = 4 >= theta_N = 2$. Mean updates yield exact centroids $mu_1 = -1$ and $mu_2 = +1$. Because $t$ is even and $k_0 / 2 < k = 2 < 2 k_0$, Algorithm 1 enters *Merge*. The pairwise distance is $D_(1 2) = |1 - (-1)| = 2.0 < theta_C = 2.5$, so $(C_1, C_2)$ merge back to $mu_(1 2) = (4(-1) + 4(1)) / 8 = 0$ with $k = 1$.

  Iteration $t = 2 m + 1$ recovers the identical state $(k = 1, mu_1 = 0)$, establishing an infinite period-$2$ limit cycle.
]

#subsec("C", "Hysteresis Stability and Finite-Step Convergence")

#theorem(title: "Hysteresis Deadband and Eventual-Quiescence Convergence", origin: "Original Proof")[
  1. *Minimum Post-Split Separation:* Whenever a cluster $C_j$ with $sigma_(j, max) > theta_S$ splits along axis $r^*$ and updates its half-cluster means to $mu_+, mu_-$, their Euclidean separation obeys
     $ norm(mu_+ - mu_-)_2 >= (n_j^2 nu_1) / (2 n_+ n_-) >= 2 nu_1. $
     Consequently, for any cluster family with MAD-to-std ratio $nu_1 / sigma_(j, max) >= c > 0$ (e.g., $c = sqrt(2 / pi) approx 0.798$ for Gaussian, $c = sqrt(3) / 2 approx 0.866$ for uniform), immediate re-merging of freshly split siblings is impossible whenever the *hysteresis separation condition* $theta_C < 2 c theta_S$ holds.
  2. *Finite-Step Convergence:* If topological edits cease after some finite iteration $t_0$ (for example, by imposing a finite budget $M_"op"$ on split/merge/prune operations), then Algorithm 1 converges in at most $S(n_"act", k_(t_0))$ subsequent iterations to a discrete CVT on $cal(X)_"act"$.
]

#proof[
  (1) From the proof of Lemma 1(3) and AM–GM ($n_+ n_- <= n_j^2 / 4$),
  $ norm(mu_+ - mu_-)_2 &>= |mu_(+, r^*) - mu_(-, r^*)| \
    &= (n_j^2 nu_1) / (2 n_+ n_-) >= 2 nu_1. $
  If $nu_1 >= c sigma_(j, max)$, then $sigma_(j, max) > theta_S$ implies $norm(mu_+ - mu_-)_2 > 2 c theta_S > theta_C$, so the sibling pair violates $D_(j ell) < theta_C$ at the subsequent merge step.

  (2) Once no pruning, splitting, or merging fires for $t >= t_0$, the active point set $cal(X)_"act"$ and cluster count $k_(t_0)$ remain constant, and Steps 2 and 4 execute exact Lloyd iterations on $cal(X)_"act"$. By strict WCSS Lyapunov descent (Alg. 1, Theorem 1), no partition can recur and the algorithm terminates in at most $S(n_"act", k_(t_0))$ steps.
]

#subsec("D", "Rotational Non-Invariance and Partial Outlier Pruning")

#proposition(title: "Coordinate-Frame Non-Invariance and Cardinality Pruning", origin: "Original Proof")[
  Let $Sigma_j in bb(R)^(d times d)$ be the empirical covariance matrix of cluster $C_j$ with largest eigenvalue $lambda_"max" (Sigma_j)$.
  1. *Rotational Sensitivity up to Factor $sqrt(d)$:* Under orthogonal rotations $Q in O(d)$, the splitting statistic $sigma_(j, max) (Q) := max_(1 <= r <= d) sqrt((Q Sigma_j Q^T)_(r r))$ satisfies the tight bounds
     $ sqrt(lambda_"max" (Sigma_j) / d) <= sigma_(j, max) (Q) <= sqrt(lambda_"max" (Sigma_j)). $
     Hence ISODATA is *not* invariant under $O(d)$: rotating an elongated cluster away from coordinate axes can suppress splitting whenever $theta_S in (sqrt(lambda_"max" / d), sqrt(lambda_"max"))$.
  2. *Opportunistic Outlier Pruning:* Any isolated outlier group of size $s < theta_N$ that splits into its own Voronoi cell is permanently pruned by Step 3; however, if $s = 1$ outlier remains absorbed in a cluster of size $>= theta_N$, the cluster centroid still diverges with Donoho–Huber breakdown point $epsilon^* = 1 / n$.
]

#proof[
  (1) Let $v_1$ ($norm(v_1)_2 = 1$) be a top eigenvector of $Sigma_j$ so $Sigma_j succ.eq lambda_"max" v_1 v_1^T$. For any $Q in O(d)$, $w := Q v_1$ satisfies $norm(w)_2 = 1$, hence
  $ (Q Sigma_j Q^T)_(r r) = e_r^T Q Sigma_j Q^T e_r <= lambda_"max" (Sigma_j), $
  giving the upper bound (attained when $w = e_1$). Conversely,
  $ max_(1 <= r <= d) (Q Sigma_j Q^T)_(r r) &>= 1 / d "Tr"(Q Sigma_j Q^T) \
    &>= lambda_"max" (Sigma_j) / d, $
  which is attained with equality when $Sigma_j = lambda_"max" v_1 v_1^T$ (rank-1 elongated cluster) and $Q$ rotates $v_1$ to $d^(-1/2) bold(1)$. Choosing $theta_S$ strictly inside $(sqrt(lambda_"max" / d), sqrt(lambda_"max"))$ causes the split predicate $sigma_(j, max) > theta_S$ to evaluate true before rotation and false after rotation.

  (2) If an extreme outlier $x_1$ ($norm(x_1)_2 -> infinity$) occupies a cluster $C_j$ with $|C_j| < theta_N$, Step 3 removes $C_j$ completely, leaving remaining centroids bounded. Conversely, if $x_1$ stays inside a cluster $C_j$ with $|C_j| >= theta_N$ (for instance when $theta_S$ is not exceeded or $t$ is on a merge cycle), $norm(mu_j)_2 >= |C_j|^(-1) norm(x_1)_2 - M_0 -> infinity$.
]

#subsec("E", "Classical Statistical Perspective on Automatic k Selection")

#theorem(title: "Fixed Thresholds vs. Consistent Information Criteria", origin: "Classical — Cited")[
  Fixed scale-dependent thresholds $(theta_S, theta_C)$ in ISODATA @ball1965isodata do not adapt to sample size $n$ or dimension $d$: under high-dimensional isotropic distributions, pairwise Euclidean distances concentrate around $Theta(sqrt(d))$ @beyer1999nearest while axis-aligned spreads remain $O(1)$, systematically decoupling $theta_S$ from $theta_C$. By contrast, penalized likelihood splitting via BIC @schwarz1978bic in X-Means @pelleg2000xmeans and Gaussianity testing via Anderson–Darling @anderson1952asymptotic in G-Means @hamerly2003gmeans scale automatically with $(n, d)$, achieving asymptotically consistent model-order selection under Gaussian mixtures @fraley2002mclust.
]

#proof(label: "Proof Synthesis (Cited)")[
  For independent coordinates with variance $sigma^2$, inter-centroid separation between distinct components grows as $Omega(sqrt(d) sigma)$ while individual axis variances stay $sigma^2$ (and shrink by $d^(-1/2)$ along rotated principal axes by Proposition 1). Thus no dimension-free fixed ratio $theta_C / theta_S$ remains calibrated across $d$. Pelleg and Moore @pelleg2000xmeans replace $(theta_S, theta_C)$ with the BIC difference $Delta "BIC" = ell(k+1) - ell(k) - (d + 1) / 2 ln n$, whose penalty grows logarithmically with $n$ and linearly with $d$, while Hamerly and Elkan @hamerly2003gmeans project each cluster onto its principal eigenvector (restoring $O(d)$ invariance) and test 1D normality at fixed significance level $alpha$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

At iteration $t$ with $k_t <= k_"max"$ clusters, Voronoi assignment costs $O(n k_t d)$, centroid and coordinate-wise variance evaluation costs $O(n d)$, and pairwise centroid distance computation plus sorting for merges costs $O(k_t^2 (d + log k_t))$. Across $T$ iterations, total runtime is $O((n k_"max" + k_"max"^2) d T)$ (simplifying to $O(n k_"max" d T)$ when $k_"max" <= n$) with $O((n + k_"max") d + k_"max"^2)$ memory.

In practice, ISODATA's operational bottleneck is *not* asymptotic per-iteration cost, but rather (1) slower convergence caused by repeated non-monotonic split–merge re-equilibration, and (2) tuning six coupled, scale-dependent thresholds $(k_0, theta_N, theta_S, theta_C, L, T)$. Best engineering practice standardizes features to unit variance, enforces the hysteresis rule $theta_C <= 1.5 theta_S$ (Theorem 2), caps total topology edits, and evaluates final partitions using internal indices such as Davies–Bouldin @davies1979cluster, Caliński–Harabasz @calinski1974dendrite, or Silhouette @rousseeuw1987silhouettes.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for ISODATA],
    [Time / Space Complexity], [$O((n k + k^2) d T)$ time; $O((n + k) d + k^2)$ space (with dynamic $k <= k_"max"$).],
    [Geometric Inductive Bias], [Variable-$k$ convex Voronoi cells with axis-aligned split geometry.],
    [Convergence Behavior], [Non-monotonic in $J$; can enter period-$2$ split–merge limit cycles (Thm. 1); requires $T$ or hysteresis (Thm. 2).],
    [Invariances], [Translation invariant; *not* invariant to rotations $Q in O(d)$ (up to $sqrt(d)$ factor, Prop. 1) or feature scaling.],
    [Outlier Handling], [Prunes isolated small clusters ($|C_j| < theta_N$), but retains $epsilon^* = 1/n$ for absorbed outliers.],
    [Recommended Upgrades], [X-Means / G-Means for parameter-free $k$; PCA-axis splitting to restore $O(d)$ invariance; Trimmed K-Means (Alg. 26) for robust pruning.],
  ),
  caption: [Summary of theoretical guarantees and operational profile of ISODATA.],
)

#sec("IV", "Engineering and Scientific Applications")

Despite predating modern model-selection theory, ISODATA remains deeply entrenched in geospatial, remote-sensing, and industrial inspection systems where domain engineers possess physical intuition for acceptable cluster spreads ($theta_S$) and minimum class separations ($theta_C$).

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Multispectral Remote Sensing & Land-Cover Mapping:* ISODATA is the historic workhorse of unsupervised satellite classification (Landsat, Sentinel-2, MODIS) in GIS suites (ENVI, ERDAS IMAGINE, ArcGIS). Because radiometric reflectance bands have physically meaningful calibration units and natural landscapes contain an unknown number of spectral sub-classes (water, forest, urban, cropland), setting $(theta_S, theta_C, theta_N)$ directly encodes sensor noise floors and minimum parcel areas.
2. *Precision Agriculture & Aerial Parcel Zoning:* Multispectral UAV imagery maps soil moisture and vegetation indices (NDVI) into dynamic management zones; cardinality pruning ($theta_N$) automatically suppresses single-pixel shadow or glint artifacts while splitting separates heterogeneous soil blocks.
3. *Medical & Histopathology Multispectral Segmentation:* In fluorescence microscopy and multi-channel MRI, anatomical structures vary in count across slices (e.g., presence or absence of lesions). Dynamic splitting accommodates unexpected pathological classes without refitting a fixed-$k$ model.
4. *Adaptive Vector Quantization & Radar Signal Sorting:* Electronic-support radar deinterleaving clusters pulse descriptors (carrier frequency, pulse width, PRI) from an unknown number of emitters; $theta_S$ bounds intra-emitter jitter while $theta_C$ merges duplicated track hypotheses.
5. *Exploratory Customer & Telemetry Cohort Discovery:* When an analyst seeks approximately $k_0$ operational segments (e.g., $k_0 in [4, 10]$) without committing rigidly to one count, ISODATA automatically eliminates micro-segments smaller than $theta_N$ and consolidates redundant overlapping cohorts.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

ISODATA degrades severely under four conditions: (a) *high-dimensional correlated subspaces*, where axis-aligned $sigma_(j, max)$ underestimates rotated variance by up to $sqrt(d)$ (Proposition 1) while inter-centroid distances inflate as $sqrt(d)$; (b) *misconfigured thresholds* violating $theta_C < 2 c theta_S$, inducing oscillatory split–merge thrashing (Theorem 1); (c) *non-convex shapes*, inherited from Euclidean Voronoi cells; and (d) *unstandardized heterogeneous units*, where a single large-scale feature dominates both $sigma_(j, max)$ and $D_(j ell)$.

#figure(
  table(
    columns: (1.05fr, 1.15fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Why ISODATA Fits], [Caution / Alternative],
    [Remote Sensing / GIS], [Calibrated bands match $(theta_S, theta_C, theta_N)$], [Normalize; cap $T$ to avoid cycles],
    [UAV Agriculture], [Pruning rejects glint/shadow pixels], [Use GMM (Alg. 81) for covariance],
    [Radar / Signal Sort.], [Dynamic emitter count with jitter bounds], [Use DBSCAN (Alg. 31) for drift],
    [Customer Cohorts], [Enforces min segment size $theta_N$], [Tune $theta_C < 1.5 theta_S$ after Z-score],
    [High-D Omics / Text], [Avoid: axis split loses $sqrt(d)$ invariance], [Use G-Means or Spherical (Alg. 12)],
  ),
  caption: [Application suitability and domain engineering guide for ISODATA.],
)

#sec("V", "Structural Relations and Empirical Behavior")

ISODATA bridges partitional $k$-means (Alg. 1), hierarchical Ward linkage (Alg. 59) @ward1963hierarchical, and divisive codebook splitting (LBG, Alg. 30) @linde1980algorithm—yet differs subtly in how it prioritizes merges.

#proposition(title: "Unweighted Merge Suboptimality Relative to Ward Linkage", origin: "Original Proof")[
  ISODATA prioritizes merges by unweighted centroid distance $D_(j ell) = norm(mu_j - mu_ell)_2$, whereas Ward's minimum-variance criterion @ward1963hierarchical orders merges by exact WCSS increase $Delta J_"merge" (j, ell) = (n_j n_ell / (n_j + n_ell)) D_(j ell)^2$ (Lemma 1). Consequently, for any $M > 1$, there exist three clusters where ISODATA greedily merges a large pair incurring WCSS increase $> M$ times larger than an available alternative merge.
]

#proof[
  Consider three clusters in $bb(R)$: two large clusters $C_1, C_2$ of size $m$ separated by $D_(1 2) = theta_C - epsilon$ ($0 < epsilon < theta_C / 2$), and a small cluster $C_3$ of size $1$ separated from $C_2$ by $D_(2 3) = theta_C - epsilon / 2$ (with $L = 1$). Because $D_(1 2) < D_(2 3) < theta_C$, ISODATA merges $(C_1, C_2)$ first, incurring WCSS increase
  $ Delta J_"merge" (1, 2) = m / 2 (theta_C - epsilon)^2. $
  Merging $(C_2, C_3)$ instead would incur only
  $ Delta J_"merge" (2, 3) = m / (m + 1) (theta_C - epsilon / 2)^2 < theta_C^2. $
  The ratio $Delta J_"merge" (1, 2) / Delta J_"merge" (2, 3) > (m / 2) (1 - epsilon / theta_C)^2 -> infinity$ as $m -> infinity$, exceeding any prescribed $M > 1$.
]

This exact contrast explains both ISODATA's historical legacy—inspiring statistically principled successors like X-Means @pelleg2000xmeans, G-Means @hamerly2003gmeans, and Trimmed K-Means (Alg. 26) @cuesta1997trimmed—and its distinct empirical signature.

On the ClusterOverview 4-dataset benchmark suite @rousseeuw1987silhouettes @jain2010data, ISODATA records: *Iris* $"Sil" = 0.447$, $"ARI" = 0.583$; *Customers* $"Sil" = 0.427$; *Multi-Shape* $"Sil" = 0.713$ ($"ARI" = "N/A"$); and *High-D Subspace+Noise* $"Sil" = 0.181$ ($"ARI" = "N/A"$), with a mean runtime of $1.432$ ms. Two empirical phenomena directly validate our theorems: (1) on low-dimensional *Customers* data, ISODATA achieves the highest Silhouette among the three methods ($0.427$ vs. $0.365$ for Global K-Means and $0.332$ for Hartigan–Wong), because dynamic split/merge and cardinality pruning adapt cluster cardinality to the data's natural density pockets; (2) conversely, on *High-D Subspace+Noise*, ISODATA drops sharply to $"Sil" = 0.181$ (well below Global K-Means's $0.303$ and Hartigan–Wong's $0.259$), confirming Proposition 1 and Theorem 3 that axis-aligned $sigma_(j, max)$ splitting and fixed $theta_C$ thresholds deteriorate under high-dimensional rotated subspaces and noise concentration, while split–merge re-equilibration inflates runtime to $1.432$ ms ($28.6 times$ slower than Hartigan–Wong).
