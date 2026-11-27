#import "../template/ieee.typ": *

#alg_entry(
  id: 44,
  name: "K-DBSCAN (Hybrid K-Means–DBSCAN)",
  category: "Cat. II: Density & Mode-Seeking",
  year: "2010–2015",
  paradigm: "Two-Stage Voronoi-Partitioned Density Clustering",
  complexity: $O(n k d T + sum_(r=1)^k n_r^2 d)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Monolithic DBSCAN (Alg. 31, @ester1996dbscan) suffers from two fundamental bottlenecks: (i) $O(n^2 d)$ pairwise distance cost whenever spatial indices degrade in moderate-to-high dimensions @schubert2017dbscanrevisited, and (ii) inability to separate coexisting clusters with markedly different local densities under a single global radius $epsilon$ @liu2007vdbscan. *K-DBSCAN* addresses both bottlenecks via a two-stage *coarse-to-fine* decomposition combining prototype partitioning with localized density connectivity.

Let $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$ and choose a coarse partition count $k in {2, dots, n - 1}$. In *Stage 1 (Coarse Voronoi Partitioning)*, K-DBSCAN seeds prototypes $bold(mu)^((0))$ via K-Means++ @arthur2007kmeanspp and runs Lloyd's algorithm @lloyd1982least (Algs. 1–2) to minimize within-cluster sum of squares (WCSS),
$ J(cal(P), bold(mu)) = sum_(r=1)^k sum_(x in P_r) norm(x - mu_r)_2^2, $
yielding $k$ disjoint coarse cells $cal(P) = {P_1, dots, P_k}$ of sizes $n_r = |P_r|$ ($sum_(r=1)^k n_r = n$), where each cell lies inside the convex Voronoi polyhedron
$ V_r (bold(mu)) = {x in bb(R)^d : norm(x - mu_r)_2 <= norm(x - mu_s)_2 space forall s != r}. $
In *Stage 2 (Localized Density Detection)*, each non-trivial cell $P_r$ is assigned a *cell-local* radius $epsilon_r > 0$ and threshold $"minPts"_r >= 2$ (typically estimated from a quantile $q in (0, 1)$ of the intra-cell $m$-nearest-neighbor distances $d_(m, r) (x)$ within $P_r$), and DBSCAN is executed independently on $(P_r, epsilon_r, "minPts"_r)$.

#definition(title: "Voronoi-Localized Density Connectivity and Collar")[
  Fix coarse cells $cal(P) = {P_1, dots, P_k}$ and local parameters ${(epsilon_r, "minPts"_r)}_(r=1)^k$. For $x in P_r$:
  1. Its *local $epsilon_r$-neighborhood* is $N_(epsilon_r)^((r)) (x) := {y in P_r : norm(x - y)_2 <= epsilon_r}$, and $x$ is a *local core point* ($x in cal(C)_r$) iff $|N_(epsilon_r)^((r)) (x)| >= "minPts"_r$.
  2. Points $u, v in P_r$ are *locally density-connected* if a chain in $cal(C)_r$ links them via steps of length $<= epsilon_r$. Maximal local connected subsets form the cell subclusters $cal(K)_r$, while $cal(N)_r = P_r backslash union.big_(C in cal(K)_r) C$ is local noise.
  3. For homogeneous $epsilon_r equiv epsilon$, the *Voronoi $epsilon$-collar* of $P_r$ is $tilde(P)_r := P_r union {y in cal(X) backslash P_r : "dist"(y, P_r) <= epsilon}$, bounded by separating bisector hyperplanes $H_(r s) = {z : chevron.l z - (mu_r + mu_s)/2, mu_s - mu_r chevron.r = 0}$.
]

The core disjoint K-DBSCAN formulation (Algorithm 1) enables heterogeneous regional scales $epsilon_r$ with zero inter-partition communication; an optional boundary-reconciliation step (analyzed in Proposition 1) merges clusters split across bisectors $H_(r s)$ when global homogeneous density connectivity is required.

#algbox(
  title: "K-DBSCAN Two-Stage Hybrid Clustering",
  inputs: [Dataset $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$, coarse $k >= 2$, Lloyd cap $T$, local quantile $q in (0, 1)$, base $"minPts"$.],
  outputs: [Global cluster partition $cal(C)$ and noise set $cal(N)$.],
)[
  1. *Stage 1 (Coarse Partition):* Initialize $bold(mu)^((0))$ with K-Means++; execute up to $T$ Lloyd iterations to obtain disjoint Voronoi cells $P_1, dots, P_k$.
  2. Initialize global cluster counter $"cid" <- 0$, labels $ell_i <- -1$.
  3. *Stage 2 (Fine-Grained Local DBSCAN):* *For* $r = 1, dots, k$ (in parallel) *do*:
  4. #h(1em) If $n_r = |P_r| < "minPts" + 1$, assign $P_r$ as a small cluster or noise and *continue*.
  5. #h(1em) Set $m_r <- min("minPts", floor(n_r / 2))$ and compute local $m_r$-NN distances inside $P_r$.
  6. #h(1em) Estimate local radius $epsilon_r <- "Quantile"_q ({d_(m_r, r) (x) : x in P_r})$.
  7. #h(1em) Run DBSCAN on $(P_r, epsilon_r, m_r)$ to obtain local clusters $cal(K)_r$ and noise $cal(N)_r$.
  8. #h(1em) Relabel each $C in cal(K)_r$ with a unique global ID $"cid" <- "cid" + 1$.
  9. *(Optional)* If homogeneous $epsilon$ is used, reconcile shared core points across $epsilon$-collars $tilde(P)_r inter tilde(P)_s$ via Union-Find.
  10. *Return* $(cal(C), cal(N))$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five formal results governing K-DBSCAN: (1) exact pairwise-work reduction as a function of cell-size imbalance and the closed-form optimal coarse count $k^*$; (2) constructive proof of exact multi-density recovery where monolithic DBSCAN provably fails; (3) an explicit counterexample proving Voronoi boundary fracture and spurious noise without reconciliation; (4) exact equivalence to global DBSCAN under $epsilon$-collar reconciliation; and (5) classical seeding, invariance, and high-dimensional indexing bounds.

#subsec("A", "Quadratic Speedup and Optimal Coarse Partitioning")

#theorem(title: "Exact Pairwise Work Reduction and Optimal Coarse Count k*", origin: "Original Proof")[
  Let $P_1, dots, P_k$ be any non-empty partition of $n$ points with cell sizes $n_r >= 1$ ($sum_(r=1)^k n_r = n$), mean cell size $overline(n) = n / k$, and squared coefficient of variation $"CV"^2 = k^(-1) sum_(r=1)^k (n_r / overline(n) - 1)^2$. Then:
  1. The exact Stage-2 pairwise distance count $W_"loc" = sum_(r=1)^k binom(n_r, 2)$ satisfies
     $ W_"loc" = n^2 / (2 k) (1 + "CV"^2) - n / 2, $
     bounded tightly by
     $ n^2 / (2 k) - n / 2 <= W_"loc" <= binom(n - k + 1, 2). $
  2. Modeling total distance evaluations across Stages 1 and 2 as $F(k) = c_1 n k T + c_2 n^2 (1 + "CV"^2) / (2 k)$ (with constants $c_1, c_2 > 0$ and fixed $"CV"^2$), $F$ is strictly convex on $(0, infinity)$ with unique minimizer
     $ k^* = sqrt((c_2 (1 + "CV"^2) n) / (2 c_1 T)) = Theta(sqrt(n / T)), $
     reducing unindexed complexity from $O(n^2 d)$ to $O(n^(3/2) d sqrt(T))$.
]

#proof[
  (1) Expand the variance of the cell sizes $n_r$:
  $ sum_(r=1)^k n_r^2 &= overline(n)^2 sum_(r=1)^k ((n_r / overline(n) - 1) + 1)^2 \
                      &= overline(n)^2 (k "CV"^2 + 2 sum_(r=1)^k (n_r / overline(n) - 1) + k) \
                      &= k overline(n)^2 (1 + "CV"^2) = n^2 / k (1 + "CV"^2), $
  since $sum_(r=1)^k (n_r / overline(n) - 1) = 0$. Therefore,
  $ W_"loc" = 1 / 2 sum_(r=1)^k n_r (n_r - 1) = n^2 / (2 k) (1 + "CV"^2) - n / 2. $
  Because $"CV"^2 >= 0$ with equality iff $n_1 = dots = n_k = n / k$, the lower bound follows immediately (and coincides with Cauchy–Schwarz). For the upper bound, $g(u) = u(u - 1) / 2$ is convex on $[1, infinity)$, so the separable convex sum $sum_(r=1)^k g(n_r)$ over the simplex ${n_r >= 1, sum_r n_r = n}$ attains its maximum at an extreme point $(n - k + 1, 1, dots, 1)$, where $W_"loc" = binom(n - k + 1, 2)$.

  (2) Differentiating $F(k)$ with respect to continuous $k > 0$:
  $ F'(k) = c_1 n T - (c_2 (1 + "CV"^2) n^2) / (2 k^2), quad F''(k) = (c_2 (1 + "CV"^2) n^2) / k^3 > 0. $
  Strict convexity implies the unique stationary point $F'(k^*) = 0$ is the global minimizer $k^* = sqrt(c_2 (1 + "CV"^2) n / (2 c_1 T))$. Evaluating $F(k^*)$ yields
  $ F(k^*) = sqrt(2 c_1 c_2 (1 + "CV"^2) T) n^(3/2), $
  and multiplying by $O(d)$ per Euclidean distance establishes $O(n^(3/2) d sqrt(T))$ total arithmetic work without relying on spatial trees.
]

#subsec("B", "Exact Multi-Density Separation Beyond Monolithic DBSCAN")

#theorem(title: "Exact Recovery Under Arbitrarily Large Density Contrast", origin: "Original Proof")[
  Fix any integer $m >= 2$ and any target density-scale contrast ratio $Lambda > 1$. There exists a dataset $cal(X) = C_(1 a) union C_(1 b) union C_2 subset.eq bb(R)^d$ of three separated clusters such that:
  1. *Monolithic DBSCAN impossibility:* For *every* global radius $epsilon > 0$ (at $"minPts" = m$), standard DBSCAN @ester1996dbscan fails to recover ${C_(1 a), C_(1 b), C_2}$—either merging $C_(1 a)$ and $C_(1 b)$ or dissolving $C_2$ into noise.
  2. *Exact K-DBSCAN recovery:* K-DBSCAN with $k = 2$ and local scale estimation recovers the exact three-cluster partition ${C_(1 a), C_(1 b), C_2}$ with zero noise.
]

#proof[
  Fix base spacing $delta_1 > 0$, gap $g = 2 delta_1$, and sparse spacing $delta_2 = (Lambda + 2) delta_1 > g$. Construct:
  - $C_(1 a)$ and $C_(1 b)$: two dense clusters of $2 m$ points each, internally $m$-connected at pairwise scale $delta_1$, separated by exact minimum inter-cluster distance $"dist"(C_(1 a), C_(1 b)) = g$, both contained inside the ball $B_R (bold(0))$ of radius $R = 3 m delta_1$.
  - $C_2$: a sparse regular chain/lattice of $2 m$ points whose nearest-neighbor spacing is $delta_2$ (so every point's $m$-th neighbor distance satisfies $d_m (x) >= delta_2$), contained inside $B_(R_2) (L e_1)$ with $R_2 = m delta_2$ and macro-separation $L = 6 max(R, R_2)$.

  (1) Consider monolithic DBSCAN with fixed $"minPts" = m$ and any global $epsilon > 0$. If $epsilon < delta_2$, then every $x in C_2$ has $|B_epsilon (x) inter cal(X)| < m$ (as even its $m$-th nearest neighbor lies at distance $>= delta_2 > epsilon$), so all of $C_2$ is classified as noise. Conversely, if $epsilon >= delta_2$, then $epsilon > g = "dist"(C_(1 a), C_(1 b))$, so facing core points in $C_(1 a)$ and $C_(1 b)$ become directly density-reachable and merge into a single cluster. Thus no global $epsilon$ recovers all three clusters.

  (2) Consider K-DBSCAN with $k = 2$. Because $L = 6 max(R, R_2)$, any partition placing $P_1 = C_(1 a) union C_(1 b)$ and $P_2 = C_2$ has centroids $mu_1 in B_R (bold(0))$ and $mu_2 in B_(R_2) (L e_1)$. For every $x in P_1$,
  $ norm(x - mu_1)_2 <= 2 R < 4 max(R, R_2) <= L - R - R_2 <= norm(x - mu_2)_2, $
  and symmetrically for every $x in P_2$, $norm(x - mu_2)_2 <= 2 R_2 < norm(x - mu_1)_2$. Thus $(P_1, P_2)$ is a strict global minimizer of WCSS (and the unique CVT whenever K-Means++ draws one seed per macro-ball, which occurs with probability $1 - O((R_2 / L)^2)$). Inside $P_1$, every point's intra-cell $m$-NN distance is $<= delta_1$, so local scale selection yields $epsilon_1 in [delta_1, g)$ and local DBSCAN cleanly separates $C_(1 a)$ and $C_(1 b)$. Inside $P_2$, local $m$-NN distances equal $Theta(delta_2)$, giving $epsilon_2 >= delta_2$ and recovering $C_2$ intact.
]

#subsec("C", "Voronoi Boundary Fracture and Exact Collar Reconciliation")

Without cross-cell communication, Stage 1's convex Voronoi bisectors can slice continuous non-convex or elongated density manifolds.

#theorem(title: "Voronoi Over-Segmentation and Boundary Noise Without Reconciliation", origin: "Original Counterexample")[
  Let $cal(X) = {x_j = j delta : j = 1, dots, n} subset.eq bb(R)$ be a single uniform density cluster of $n = 2 q m$ points ($q >= 3$, $m >= 2$) with homogeneous DBSCAN parameters $epsilon = m delta$ and $"minPts" = m$. Then:
  1. Monolithic DBSCAN recovers the single cluster $cal(X)$ with zero noise.
  2. Disjoint K-DBSCAN with $k = 2$ (without collar reconciliation) deterministically bisects $cal(X)$ into *two separate clusters* $P_1 = {x_1, dots, x_(q m)}$ and $P_2 = {x_(q m + 1), dots, x_(2 q m)}$.
  3. Moreover, if a 1-D cluster has local spacing such that a boundary point $x in P_1$ has only $a < m$ neighbors in $P_1$ within $epsilon$ and $b >= m - a$ neighbors in $P_2$ within $epsilon$, disjoint K-DBSCAN demotes $x$ from a true global core point to border or spurious noise.
]

#proof[
  (1) On $cal(X)$, every point $x_j$ has at least $m$ points (including itself) within radius $epsilon = m delta$ and consecutive points satisfy $|x_(j+1) - x_j| = delta <= epsilon$; hence all $n$ points are core and form one connected component.

  (2) For $k = 2$ on symmetric equispaced points, the unique optimal 2-means CVT splits $cal(X)$ at the midpoint $(n + 1) delta / 2$, producing equal halves $P_1$ and $P_2$ with centroids $mu_1 = (q m + 1) delta / 2$ and $mu_2 = (3 q m + 1) delta / 2$. Because Stage 2 executes DBSCAN separately on $P_1$ and $P_2$ and assigns distinct global IDs, no edge ever crosses the bisector $H_(1 2)$, splitting the single true cluster into two fragments (for general $k$, up to $k$ fragments).

  (3) For $x_j$ adjacent to $H_(1 2)$, $N_epsilon^((1)) (x_j) = B_epsilon (x_j) inter P_1$ strictly excludes $B_epsilon (x_j) inter P_2$. Setting $epsilon = (m - 1) delta$ with $"minPts" = m + 1$ makes $x_(q m)$ a global core point ($|B_epsilon (x_(q m)) inter cal(X)| = 2 m - 1 >= m + 1$ for $m >= 2$), whereas $|N_epsilon^((1)) (x_(q m))| = m < m + 1$, strictly demoting $x_(q m)$ out of the core set; at an isolated neck of width $< m$, such points become false noise.
]

#proposition(title: "Exact Global DBSCAN Equivalence via One-Epsilon Collar Merge", origin: "Original Proof")[
  Suppose all cells use a common homogeneous threshold $(epsilon, "minPts")$. Augment each Voronoi cell $P_r$ with its $epsilon$-collar $tilde(P)_r = P_r union {y in cal(X) backslash P_r : "dist"(y, P_r) <= epsilon}$ (filtered in $O(d)$ via the hyperplane distance $|chevron.l y - (mu_r + mu_s)/2, mu_s - mu_r chevron.r| / norm(mu_s - mu_r)_2 <= epsilon$). Evaluate each $x in P_r$'s core status using $N_epsilon^(tilde(P)_r) (x)$, run local core expansion on $P_r$, and merge the cluster IDs of any core pair $(u, v)$ with $u in P_r$, $v in tilde(P)_r inter P_s$, and $norm(u - v)_2 <= epsilon$. Then the resulting core-point partition is *identically equal* to monolithic DBSCAN's core-point partition.
]

#proof[
  First, for every $x in P_r$ and any $y in cal(X)$ with $norm(y - x)_2 <= epsilon$, we have $"dist"(y, P_r) <= norm(y - x)_2 <= epsilon$, so $y in tilde(P)_r$. Hence $N_epsilon^(tilde(P)_r) (x) = B_epsilon (x) inter cal(X)$ identically, guaranteeing that every $x in cal(X)$ receives its exact global core classification $x in cal(C)$.

  Second, let $G_"core" = (cal(C), E_epsilon)$ be monolithic DBSCAN's undirected core adjacency graph, where ${u, v} in E_epsilon$ iff $norm(u - v)_2 <= epsilon$. For any edge ${u, v} in E_epsilon$ with $u in P_r$ and $v in P_s$: if $r = s$, the edge is traversed during local DBSCAN inside $P_r$; if $r != s$, then $v in tilde(P)_r inter P_s$ and $norm(u - v)_2 <= epsilon$, so the collar reconciliation pass triggers a Union-Find merge between the local components containing $u$ and $v$. Conversely, no merge ever occurs without an edge in $E_epsilon$. Therefore the merged components match the connected components of $G_"core"$ identically (and border points attach to any core neighbor within $epsilon$, matching global DBSCAN up to standard border tie-breaking @schubert2017dbscanrevisited).
]

#subsec("D", "Seeding Quality, Invariances, and High-D Indexing Bounds")

#theorem(title: "Coarse Seeding Guarantee, Breakdown, and Indexing Context", origin: "Classical — Cited")[
  1. (*K-Means++ Coarse Quality:*) Stage-1 K-Means++ seeding guarantees $bb(E)[J(cal(P)^((0)))] <= 8 (ln k + 2) J^*$ @arthur2007kmeanspp, preventing degenerate empty cells and curbing $"CV"^2$ in Theorem 1.
  2. (*Invariances & Outlier Decoupling:*) Both stages are strictly invariant to rigid Euclidean isometries $Q x + b$ ($Q in O(d)$) and sensitive to anisotropic scaling. While Stage 1 inherits K-Means' $1/n$ centroid breakdown point (a single extreme outlier can hijack one coarse cell), Stage-2 DBSCAN isolates that outlier as noise inside its cell without corrupting core clusters in the remaining $k - 1$ cells.
  3. (*High-Dimensional DBSCAN Bottleneck:*) Schubert et al. @schubert2017dbscanrevisited show that spatial-tree DBSCAN acceleration collapses to $O(n^2 d)$ in high dimensions due to nearest-neighbor distance concentration @beyer1999nearest, whereas coarse Voronoi partitioning retains its $O(k)$ work factor without the $g^d$ exponential cell explosion of axis-aligned grids @mahran2008gridbscan.
]

#proof(label: "Proof Synthesis (Cited)")[
  Arthur and Vassilvitskii @arthur2007kmeanspp bound expected $D^2$-seeding WCSS by induction on covered optimal clusters, ensuring prototypes spread proportionally to macro-cluster mass. Rigid motions preserve both squared distances (hence Voronoi bisectors) and local radii $epsilon_r$. In high dimension $d$, axis-aligned grid partitioning (GriDBSCAN @mahran2008gridbscan) requires $O(3^d)$ neighbor-cell probes per query, whereas Voronoi partitioning maintains exactly $k$ data-adaptive cells regardless of $d$, with hyperplane collar screening costing only $O(k d)$ per point and zero dimensional explosion.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Stage 1 executes in $O(n k d T)$ time (further accelerable via Elkan's triangle-inequality pruning @elkan2003using, Alg. 13) and $O((n + k) d)$ memory. Stage 2 executes in $O(d sum_(r=1)^k n_r^2)$ exact time (or $O(sum_(r=1)^k n_r log n_r)$ using per-cell $k$-d trees in low $d$) and requires only $O(max_r n_r^2)$ temporary scratch memory if local distance matrices are materialized (or $O(n)$ streaming). Crucially, each cell $P_r$ of size $approx n / k$ fits cleanly in L2/L3 CPU cache and runs independently across threads or distributed workers, yielding wall-clock speedups that frequently exceed the theoretical factor $k$ due to reduced cache misses.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for K-DBSCAN],
    [Time / Space Complexity], [$O(n k d T + sum_r n_r^2 d)$; optimal $k^* = Theta(sqrt(n / T))$ gives $O(n^(3/2) d sqrt(T))$.],
    [Geometric Inductive Bias], [Arbitrary-shape local density clusters inside convex Voronoi macro-zones.],
    [Multi-Density Capability], [High via cell-local $epsilon_r$ estimation (Theorem 2); exact global match with $epsilon$-collar (Prop. 1).],
    [Determinism / Seeding], [Stochastic via Stage-1 K-Means++ seeding; deterministic given coarse partition $cal(P)$.],
    [Outlier Robustness], [High in Stage 2 (noise label $-1$), though extreme outliers can occupy a Stage-1 cell.],
    [Recommended Companions], [HDBSCAN (Alg. 32) or RNN-DBSCAN (Alg. 45) for fracture-free multi-density; MR-DBSCAN (Alg. 49) for distributed exactness.],
  ),
  caption: [Operational and theoretical profile of K-DBSCAN.],
)

#sec("IV", "Engineering and Scientific Applications")

K-DBSCAN is engineered for large-scale or multi-zone datasets where monolithic DBSCAN is either too slow or unable to accommodate regional density contrasts.

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Urban Geospatial & LiDAR Point-Cloud Segmentation:* Multi-million-point airborne LiDAR scans and GPS trajectory logs exhibit dense urban cores alongside sparse suburban fringes. Coarse Voronoi zoning isolates cache-resident spatial tiles in $O(n k d T)$, while localized DBSCAN extracts non-convex road/building footprints and strips sensor noise.
2. *Multi-Region Telecom & IoT Traffic Analytics:* Cell-tower load and sensor telemetry vary by orders of magnitude across metropolitan and rural zones. By fitting a separate local radius $epsilon_r$ inside each coarse macro-region (Theorem 2), K-DBSCAN detects localized outage bursts and irregular traffic clusters without wiping out low-density regions as noise.
3. *Distributed & Sharded Big-Data Clustering:* Coarse K-Means centroids act as compact $O(k d)$ routing tables across compute nodes—mirroring MR-DBSCAN @he2011mrdbscan and PDBSCAN @xu1999pdbscan—allowing embarrassing parallelism in Stage 2 with optional lightweight collar-point exchange (Proposition 1).
4. *Moderate-Dimensional Customer & Behavioral Segmentation:* In $10$–$50$ standardized business features where spatial trees degenerate and grid methods suffer $g^d$ explosion, K-DBSCAN delivers rapid sub-millisecond partitioning while rejecting anomalous accounts inside each coarse cohort.
5. *Astronomical Survey Tiling & Spatial Epidemiology:* Regional background densities shift smoothly across galactic latitude or population density; localized DBSCAN inside coarse spatial partitions identifies irregular filaments or disease clusters relative to regional baselines.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

K-DBSCAN exhibits three primary failure regimes: (a) *manifold over-segmentation* (Theorem 3), when disjoint Stage-1 Voronoi boundaries bisect continuous arcs, rings, or large homogeneous blobs—mitigated by collar reconciliation (Proposition 1) or by switching to OPTICS (Alg. 33, @ankerst1999optics) / HDBSCAN (Alg. 32, @campello2013hdbscan); (b) *cell-size skew or outlier hijacking*, where extreme Stage-1 outliers consume coarse prototypes (mitigated by K-Medoids/Trimmed coarse seeding); and (c) *under-populated local cells*, where small $n_r$ causes noisy quantile estimates for $epsilon_r$.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Why K-DBSCAN Fits], [Caution / Alternative],
    [Geospatial / LiDAR], [$k$-fold speedup; regional $epsilon_r$], [Enable $epsilon$-collar merge (Prop. 1)],
    [Customer Cohorts], [Fast local outlier filtering], [Standardize features first],
    [Distributed Logs], [Natural Voronoi shard routing], [Use MR-DBSCAN (Alg. 49) for exactness],
    [Interlocking Manifolds], [Risk of Voronoi bisection], [Use HDBSCAN (Alg. 32) or RNN-DBSCAN (Alg. 45)],
    [Ultra-High-D Sparse], [Euclidean $d$-concentration], [Use SNN (Alg. 36) or Spherical K-Means (Alg. 12)],
  ),
  caption: [Domain suitability and engineering guide for K-DBSCAN.],
)

#sec("V", "Structural Relations and Empirical Behavior")

K-DBSCAN continuously interpolates between centroid partitioning (Category I) and density-based connectivity (Category II).

#proposition(title: "Exact Boundary Reductions to DBSCAN and Lloyd's K-Means", origin: "Original Derivation")[
  Let $cal(C)_"K-DB"(k, {epsilon_r}, {"minPts"_r})$ denote the partition of disjoint K-DBSCAN (Algorithm 1). Then:
  1. (*Monolithic DBSCAN Limit:*) Setting $k = 1$ yields $P_1 = cal(X)$, so $cal(C)_"K-DB"(1, {epsilon}, {"minPts"})$ reduces *identically* to standard DBSCAN (Alg. 31, @ester1996dbscan).
  2. (*Lloyd's K-Means Limit:*) Setting $"minPts"_r = 1$ and $epsilon_r >= "diam"(P_r)$ for every $r in {1, dots, k}$ makes every non-empty Voronoi cell $P_r$ a single connected clique with zero noise, so $cal(C)_"K-DB"$ reduces *identically* to K-Means++ / Lloyd's K-Means (Algs. 1–2, @lloyd1982least, @arthur2007kmeanspp).
]

#proof[
  (1) When $k = 1$, Stage 1 assigns every point to $P_1 = cal(X)$ in one trivial step, and Stage 2 invokes DBSCAN verbatim on $(cal(X), epsilon, "minPts")$.

  (2) When $"minPts"_r = 1$ and $epsilon_r >= "diam"(P_r) = max_(u, v in P_r) norm(u - v)_2$, every $x in P_r$ satisfies $N_(epsilon_r)^((r)) (x) = P_r$, hence $|N_(epsilon_r)^((r)) (x)| = |P_r| >= 1 = "minPts"_r$. Therefore every point in $P_r$ is core and directly density-reachable from every other point in $P_r$. Stage 2 outputs exactly one cluster per non-empty Voronoi cell $P_r$ with $cal(N) = emptyset$, coinciding identically with the Stage-1 K-Means partition.
]

Within Category II, K-DBSCAN complements GriDBSCAN (Alg. 48, @mahran2008gridbscan, which partitions via axis-aligned grids rather than data-adaptive Voronoi cells), MR-DBSCAN (Alg. 49, @he2011mrdbscan, distributed balanced partitioning with collar merging), VDBSCAN (Alg. 50, @liu2007vdbscan, global $k$-dist stratification across densities), DBSCAN++ (Alg. 47, @jang2019dbscanpp, subsampling candidate cores), and RNN-DBSCAN (Alg. 45, @bryant2018rnndbscan, parameter-free multi-density via reverse nearest neighbors).

Empirically on the ClusterOverview 4-dataset benchmark, K-DBSCAN records a mean runtime of just $0.351$ ms—over $20 times$ faster than DENCLUE 2.0 ($7.234$ ms) and $4 times$ faster than RNN-DBSCAN ($1.456$ ms)—directly reflecting Theorem 1's $k$-fold distance reduction and cache locality. On heterogeneous and noisy benchmarks it outperforms both neighbors, achieving *Customers* $"Sil" = 0.526$ and *High-D Subspace+Noise* $"Sil" = 0.358$ ($"ARI" = "N/A"$), validating Theorem 2's regional scale adaptation and Stage-2 noise rejection. Conversely, its lower scores on *Iris* ($"Sil" = 0.392$, $"ARI" = 0.333$) and *Multi-Shape* ($"Sil" = 0.556$, $"ARI" = "N/A"$) empirically confirm Theorem 3: without collar reconciliation, coarse Voronoi bisectors split continuous manifolds into over-segmented local clusters.
