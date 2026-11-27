#import "../template/ieee.typ": *

#alg_entry(
  id: 31,
  name: "DBSCAN (Density-Based Spatial Clustering of Applications with Noise)",
  category: "Cat. II: Density & Mode-Seeking",
  year: "1996",
  paradigm: "Fixed-Radius Density-Level Connected Components",
  complexity: $O(n log n) "/" O(n^2)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Let $(cal(X), d)$ be a finite metric space with dataset $cal(X) = {x_1, dots, x_n}$ (canonically $cal(X) subset.eq bb(R)^d$ with Euclidean metric $d(p, q) = norm(p - q)_2$). Rather than optimizing a global parametric distortion over a pre-specified cluster count $k$, DBSCAN @ester1996dbscan defines clusters topologically as maximal connected components of an empirical upper density level set, parameterized by a neighborhood radius $epsilon > 0$ and a minimum point count $m = "MinPts" in bb(N)_(>= 2)$.

For every point $p in cal(X)$, define its closed $epsilon$-neighborhood by
$ N_epsilon (p) = { q in cal(X) : d(p, q) <= epsilon }. $
In $bb(R)^d$, letting $V_d = pi^(d/2) / Gamma(d/2 + 1)$ denote the volume of the unit Euclidean ball, the spherical uniform-kernel density estimator at $p$ is $hat(f)_epsilon (p) = |N_epsilon (p)| / (n V_d epsilon^d)$. Thresholding $hat(f)_epsilon (p) >= lambda$ at density level $lambda = m / (n V_d epsilon^d)$ induces a canonical partition of $cal(X)$ into three disjoint strata:
$ cal(X)_"core" &= { p in cal(X) : |N_epsilon (p)| >= m }, \
  cal(X)_"border" &= { q in cal(X) backslash cal(X)_"core" : exists p in cal(X)_"core" "s.t." q in N_epsilon (p) }, \
  cal(X)_"noise" &= cal(X) backslash (cal(X)_"core" union cal(X)_"border"). $

#definition(title: "Density-Reachability, Connectivity, and DBSCAN Clusters")[
  Fix $epsilon > 0$ and $m >= 2$.
  1. *Directly density-reachable ($p ->_(epsilon, m) q$):* $q in N_epsilon (p)$ and $p in cal(X)_"core"$.
  2. *Density-reachable ($p ->_(epsilon, m)^* q$):* There exists a finite chain $p = p_0, p_1, dots, p_r = q$ in $cal(X)$ ($r >= 1$) with $p_(j-1) ->_(epsilon, m) p_j$ for all $j = 1, dots, r$.
  3. *Density-connected ($p <->_(epsilon, m) q$):* There exists $o in cal(X)_"core"$ such that $o ->_(epsilon, m)^* p$ and $o ->_(epsilon, m)^* q$.
  4. *DBSCAN Cluster:* A non-empty subset $C subset.eq cal(X)$ satisfying *maximality* ($forall p in C, space p ->_(epsilon, m)^* q => q in C$) and *connectivity* ($forall p, q in C, space p <->_(epsilon, m) q$). All remaining points constitute noise.
]

Because $->_(epsilon, m)$ is asymmetric between core and border points, DBSCAN discovers clusters via breadth-first graph expansion initiated exclusively from unvisited core points (Algorithm 1), while leaving unexpanded non-core points in $cal(X)_"noise"$ @ester1996dbscan @schubert2017dbscanrevisited.

#algbox(
  title: "DBSCAN Algorithm (Ester et al., 1996)",
  inputs: [Metric dataset $(cal(X), d)$ of size $n$, radius $epsilon > 0$, threshold $m = "MinPts" >= 2$.],
  outputs: [Cluster labels $L(p) in {1, dots, k} union {"noise"}$ for all $p in cal(X)$.],
)[
  1. Initialize $L(p) <- "unvisited"$ for all $p in cal(X)$, and cluster counter $c <- 0$.
  2. *For* each $p in cal(X)$ *do*:
  3. #h(1em) *If* $L(p) != "unvisited"$, *continue*.
  4. #h(1em) Query $N <- N_epsilon (p) = {q in cal(X) : d(p, q) <= epsilon}$.
  5. #h(1em) *If* $|N| < m$, set $L(p) <- "noise"$ and *continue*.
  6. #h(1em) Increment $c <- c + 1$; set $L(p) <- c$; initialize queue $Q <- N backslash {p}$.
  7. #h(1em) *While* $Q != emptyset$ *do*:
  8. #h(2em) Pop $q$ from $Q$.
  9. #h(2em) *If* $L(q) = "noise"$, set $L(q) <- c$ (border point absorption).
  10. #h(2em) *If* $L(q) != "unvisited"$, *continue*.
  11. #h(2em) Set $L(q) <- c$ and query $N' <- N_epsilon (q)$.
  12. #h(2em) *If* $|N'| >= m$, append $N' backslash {q}$ to $Q$.
  13. *Return* label map $L$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five foundational results governing DBSCAN: (1) core-set equivalence and the exact threshold $m = 4$ for permutation sensitivity; (2) an explicit impossibility proof showing fixed-$epsilon$ failure on multi-density clusters; (3) quantitative adversarial robustness bounds against outlier injection and bridge merging; (4) classical asymptotic level-set consistency; and (5) high-dimensional radius concentration.

#subsec("A", "Core Equivalence Relation and Exact Order-Sensitivity Threshold")

#lemma(title: "Core-Graph Equivalence and Order-Invariant Core Components", origin: "Original Proof")[
  Restricted to $cal(X)_"core"$, density-reachability $->_(epsilon, m)^*$ is an equivalence relation and coincides identically with graph connectivity in the undirected core $epsilon$-graph $G_"core" = (cal(X)_"core", E_"core")$, where ${p, q} in E_"core" <=> d(p, q) <= epsilon$. Consequently, the core partition ${C inter cal(X)_"core"}$ and the noise set $cal(X)_"noise"$ produced by Algorithm 1 are strictly unique and invariant to input ordering.
]

#proof[
  Let $p, q, w in cal(X)_"core"$. *Reflexivity:* Since $d(p, p) = 0 <= epsilon$ and $p in cal(X)_"core"$, $p ->_(epsilon, m) p$, hence $p ->_(epsilon, m)^* p$. *Transitivity:* Concatenating chains $p ->_(epsilon, m)^* q$ and $q ->_(epsilon, m)^* w$ yields $p ->_(epsilon, m)^* w$. *Symmetry:* Suppose $p = p_0, p_1, dots, p_r = q$ satisfies $p_(j-1) ->_(epsilon, m) p_j$ for $1 <= j <= r$. For every $j in {0, dots, r - 1}$, $p_j in cal(X)_"core"$ by definition of direct density-reachability, and $p_r = q in cal(X)_"core"$ by hypothesis. Since $d$ is symmetric, $d(p_j, p_(j-1)) = d(p_(j-1), p_j) <= epsilon$ implies $p_j ->_(epsilon, m) p_(j-1)$ for all $j = r, dots, 1$. Reversing the chain yields $q ->_(epsilon, m)^* p$, and every edge ${p_(j-1), p_j}$ lies in $E_"core"$.

  Conversely, every path in $G_"core"$ induces a valid reachability chain. In Algorithm 1, a label $c$ expands through Step 12 exclusively across edges incident to core points, visiting the entire connected component of $G_"core"$ containing its seed before terminating the loop. Because connected components of an undirected graph and the sets $cal(X)_"core"$, $cal(X)_"noise"$ depend solely on the unordered metric space $(cal(X), d)$, they are strictly invariant to data permutation.
]

#proposition(title: "Exact Border-Point Ambiguity Threshold at m = 4", origin: "Original Counterexample")[
  The full output partition of Algorithm 1 is strictly permutation-invariant for all $m <= 3$, but fails permutation invariance for every $m >= 4$: for every $m >= 4$, there exists a one-dimensional dataset of size $n = 2m + 1$ whose DBSCAN border assignment depends strictly on the traversal order.
]

#proof[
  *Part 1 ($m <= 3$):* Suppose a border point $q in cal(X)_"border"$ could be absorbed into two distinct core components $K_1 != K_2$ depending on ordering. Then there exist core points $p_1 in K_1$ and $p_2 in K_2$ with $q in N_epsilon (p_1) inter N_epsilon (p_2)$. Since $K_1 inter K_2 = emptyset$ and $q in.not cal(X)_"core"$, the points $p_1, p_2, q$ are pairwise distinct, and symmetry of $d$ gives ${p_1, p_2, q} subset.eq N_epsilon (q)$. Thus $|N_epsilon (q)| >= 3 >= m$, contradicting $q in.not cal(X)_"core"$. Every border point therefore connects to a unique core component, making Algorithm 1 strictly order-invariant for $m <= 3$.

  *Part 2 ($m >= 4$):* Fix $epsilon = 1$ and any $m >= 4$. Construct $2m + 1$ distinct points in $bb(R)$:
  $ A &= {0} union { -0.1 - 0.05 j / m : j = 1, dots, m - 1 }, \
    B &= {2} union { 2.1 + 0.05 j / m : j = 1, dots, m - 1 }, \
    q &= 1, quad cal(X) = A union B union {q}. $
  Every pair inside $A$ has distance $<= 0.15 < 1$, so $|N_1 (a)| >= m$ for all $a in A$; symmetrically, every $b in B$ has $|N_1 (b)| >= m$. Because $d(0, 2) = 2 > 1$, $A$ and $B$ form two disjoint connected components of $G_"core"$. However, $q = 1$ satisfies $d(q, 0) = d(q, 2) = 1$, while every other point of $A union B$ lies at distance $>= 1.1 > 1$ from $q$. Hence $N_1 (q) = {0, 1, 2}$ has cardinality $3 < m$, proving $q in cal(X)_"border"$ with $0 ->_(1, m) q$ and $2 ->_(1, m) q$. If Algorithm 1 encounters $0$ before $2$, Step 11 permanently assigns $q$ to $A$'s cluster; under the reversed input order, $q$ is permanently assigned to $B$'s cluster.
]

#remark(title: "Order-Invariant DBSCAN\* Variant")[
  Campello et al. @campello2015hierarchical eliminate border ambiguity altogether in *DBSCAN\** by classifying all non-core points $cal(X) backslash cal(X)_"core"$ as noise, rendering the clustering strictly identical to the connected components of $G_"core"$ for every $m >= 2$.
]

#subsec("B", "Single-Scale Multi-Density Impossibility")

#theorem(title: "Impossibility of Fixed-Epsilon Multi-Density Recovery", origin: "Original Counterexample")[
  For every $m >= 2$, there exists a one-dimensional dataset $cal(X) = C_1 union C_2 union C_3$ consisting of three well-separated uniform-grid clusters (each of size $N = 2m$) such that for *every* radius $epsilon > 0$, DBSCAN($epsilon, m$) either classifies the entire sparse cluster $C_3$ as noise or merges the two dense clusters $C_1$ and $C_2$ into a single component.
]

#proof[
  Set dense step $delta = 1$, inter-cluster gap $g = m$, sparse step $Delta = m + 1$, and offset $R = 10 m^2$. Define
  $ C_1 &= { 0, 1, dots, 2m - 1 }, \
    C_2 &= { 2m - 1 + g + j : j = 0, dots, 2m - 1 }, \
    C_3 &= { R + j Delta : j = 0, dots, 2m - 1 }. $
  Let $r = floor(m / 2) >= 1$, and consider any $epsilon > 0$:

  *Case 1 ($epsilon < r Delta$):* Every point $x in C_3$ is at distance $>= R - (4m - 1 + m) > r Delta$ from $C_1 union C_2$. Inside $C_3$, adjacent points are spaced by $Delta$, so a closed ball of radius $epsilon < r Delta$ centered at $x in C_3$ contains at most $r - 1$ points to the left and $r - 1$ points to the right of $x$. Thus
  $ |N_epsilon (x)| <= 1 + 2(r - 1) = 2 floor(m / 2) - 1 <= m - 1 < m. $
  Consequently no point in $C_3$ is core, and since no point of $C_1 union C_2$ lies within $epsilon$ of $C_3$, every element of $C_3$ is labeled noise.

  *Case 2 ($epsilon >= r Delta$):* Because $r >= 1$, we have $epsilon >= Delta = m + 1 > m = g > m - 1$. Consider the facing endpoints $u = 2m - 1 in C_1$ and $v = 2m - 1 + g in C_2$. The $m$ points ${m, m+1, dots, 2m-1} subset.eq C_1$ all lie within distance $m - 1 < epsilon$ of $u$, so $|N_epsilon (u)| >= m$; symmetrically, $|N_epsilon (v)| >= m$. Hence $u, v in cal(X)_"core"$ with $d(u, v) = g < epsilon$, so ${u, v} in E_"core"$. By Lemma 1, $C_1$ and $C_2$ merge into a single cluster. No $epsilon > 0$ can separate $C_1, C_2$ while retaining $C_3$.
]

#subsec("C", "Geometric Invariance and Adversarial Contamination Bounds")

#proposition(title: "Outlier Fabrication and Bridge-Merging Lower Bounds", origin: "Original Proof")[
  Let $cal(X) = C_1 union C_2$ consist of two core components of DBSCAN($epsilon, m$) separated by inter-component distance $D = min_(x in C_1, y in C_2) d(x, y) > epsilon$, and let an adversary inject a set $Z$ of arbitrary points into $(bb(R)^d, d)$. Then:
  1. *Spurious-Cluster Threshold:* If $d(Z, cal(X)) > epsilon$, injecting $Z$ creates a spurious cluster if and only if $|Z| >= m$; for $|Z| < m$, all points in $Z$ are isolated as noise and $cal(X)_"core"$ is unchanged.
  2. *Bridge-Merging Lower Bound:* Merging $C_1$ and $C_2$ into a single core component of $cal(X) union Z$ requires injecting at least
     $ |Z| >= max( ceil(D / epsilon) - 1, space m (ceil(D / (3 epsilon)) - 1) ) $
     adversarial points.
]

#proof[
  (1) When $d(Z, cal(X)) > epsilon$, neighborhoods of $cal(X)$ and $Z$ are disjoint: $N_epsilon (z) subset.eq Z$ for all $z in Z$. Thus $|N_epsilon (z)| <= |Z| < m$ precludes any core point in $Z$, whereas placing $m$ points inside a single ball of diameter $<= epsilon$ yields a valid core cluster of size $m$.

  (2) Let $p_0, z_1, dots, z_L, p_(L+1)$ be a shortest path in $G_"core" (cal(X) union Z)$ with $p_0 in C_1$, $p_(L+1) in C_2$, and $z_1, dots, z_L in Z$. By the triangle inequality, $D <= d(p_0, p_(L+1)) <= (L + 1) epsilon$, giving $|Z| >= L >= ceil(D / epsilon) - 1$.

  For the density-weighted bound, let $k = ceil(D / (3 epsilon)) - 1$. If $k <= 0$ the bound is trivial; assume $k >= 1$, so $3 k epsilon < D$. Define the 1-Lipschitz distance potential $phi(x) = d(x, C_1)$, satisfying $|phi(x) - phi(y)| <= d(x, y)$. Along the path, $phi(p_0) = 0$, $phi(p_(L+1)) >= D > 3 k epsilon$, and each step increases $phi$ by at most $epsilon$. For each $r in {1, dots, k}$, let $u_r$ be the first path vertex satisfying $phi(u_r) > (3r - 2) epsilon$. Its predecessor has $phi <= (3r - 2) epsilon$, whence
  $ (3r - 2) epsilon < phi(u_r) <= (3r - 1) epsilon. $
  Because $phi(u_r) > epsilon$, we have $N_epsilon (u_r) inter C_1 = emptyset$. For every $y in C_2$, the triangle inequality yields
  $ d(u_r, y) >= d(C_1, y) - phi(u_r) >= D - (3r - 1) epsilon > epsilon, $
  so $N_epsilon (u_r) inter C_2 = emptyset$. Thus $u_r in Z$ and $N_epsilon (u_r) subset.eq Z$. Moreover, for any $1 <= r < s <= k$,
  $ d(u_s, u_r) >= phi(u_s) - phi(u_r) > (3(s - r) - 1) epsilon >= 2 epsilon. $
  Hence the closed $epsilon$-balls $N_epsilon (u_1), dots, N_epsilon (u_k)$ are pairwise disjoint subsets of $Z$. Because every $u_r$ is a core vertex in $cal(X) union Z$, each ball satisfies $|N_epsilon (u_r)| >= m$, yielding $|Z| >= sum_(r=1)^k |N_epsilon (u_r)| >= k m$.
]

Rigid Euclidean isometries $T(x) = Q x + b$ ($Q in O(d)$) preserve all pairwise distances and thus leave $cal(X)_"core"$, $cal(X)_"noise"$, and $G_"core"$ strictly invariant; isotropic scaling $x |-> alpha x$ requires co-scaling $epsilon |-> alpha epsilon$, while anisotropic feature scaling alters neighborhoods non-trivially.

#subsec("D", "Asymptotic Level-Set Consistency and High-Dimensional Collapse")

#theorem(title: "Density Level-Set Consistency of DBSCAN", origin: "Classical — Cited")[
  Let $X_1, dots, X_n$ be i.i.d. draws from a continuous density $f$ on $bb(R)^d$ whose upper level set $L_f (lambda) = {x : f(x) >= lambda}$ has finitely many compact connected components separated by a positive margin with regular boundaries. If $epsilon_n -> 0$ and $m_n -> infinity$ satisfy $m_n / (n V_d epsilon_n^d) -> lambda$ and $m_n / log n -> infinity$, then with probability tending to $1$ as $n -> infinity$, the core components of DBSCAN($epsilon_n, m_n$) stand in one-to-one correspondence with the connected components of $L_f (lambda)$ @hartigan1981consistency @chaudhuri2014consistent.
]

#proof(label: "Proof Synthesis (Cited)")[
  Euclidean balls in $bb(R)^d$ form a Vapnik–Chervonenkis class of dimension $d + 2$. Applying uniform Bernstein/VC deviation bounds to the empirical density $hat(f)_(epsilon_n) (x) = |N_(epsilon_n)(x)| / (n V_d epsilon_n^d)$ yields $sup_x |hat(f)_(epsilon_n)(x) - f(x)| -> 0$ in probability whenever $n epsilon_n^d / log n -> infinity$ and $epsilon_n -> 0$ @chaudhuri2014consistent. For any tolerance $eta_n -> 0$, the sandwich inclusion $L_f (lambda + eta_n) inter cal(X) subset.eq cal(X)_"core" subset.eq L_f (lambda - eta_n)$ holds with high probability. Regularity of $partial L_f (lambda)$ ensures that distinct components of $L_f (lambda - eta_n)$ remain separated by distance $> epsilon_n$ (preventing false merges), while an $epsilon_n / 2$-covering net inside each component of $L_f (lambda + eta_n)$ receives at least one sample per ball w.h.p., ensuring internal connectivity in $G_"core"$ @hartigan1981consistency.
]

#proposition(title: "High-Dimensional Concentration of Euclidean Pairwise Distances", origin: "Original Derivation")[
  Let $X_1, dots, X_n ~^(i.i.d.) cal(N)(0, I_d)$ in $bb(R)^d$. For every pair $i != j$ and any $t > 0$,
  $ bb(P)( | norm(X_i - X_j)_2 / sqrt(2 d) - 1 | >= t ) <= 2 / (d t^2). $
  In particular, as $d -> infinity$ with $n = o(exp(c d))$, all pairwise distances concentrate tightly around $sqrt(2d)$, collapsing the contrast ratio $(D_"max" - D_"min") / D_"min" ->^(P) 0$ @beyer1999nearest.
]

#proof[
  Because $X_i - X_j ~ cal(N)(0, 2 I_d)$, the normalized squared distance $Z = norm(X_i - X_j)_2^2 / (2 d)$ is distributed as $chi_d^2 / d$, with mean $bb(E)[Z] = 1$ and variance $"Var"(Z) = 2 / d$. Since $sqrt(Z) >= 0$, the elementary inequality $|sqrt(Z) - 1| = |Z - 1| / (sqrt(Z) + 1) <= |Z - 1|$ holds pointwise. Applying Chebyshev's inequality gives
  $ bb(P)(|sqrt(Z) - 1| >= t) <= bb(P)(|Z - 1| >= t) <= ("Var"(Z)) / t^2 = 2 / (d t^2). $
  Tightening via Laurent–Massart sub-exponential tails for $chi_d^2$ and taking a union bound over all $binom(n, 2)$ pairs implies that every distance lies inside $sqrt(2d) (1 plus.minus O(sqrt((log n) / d)))$ w.h.p. Thus a tiny perturbation of $epsilon$ around $sqrt(2d)$ transitions abruptly from empty neighborhoods (all noise) to a single percolated giant component.
]

#sec("III", "Complexity, Geometry and Operational Profile")

A crucial algorithmic property emphasized by Schubert et al. @schubert2017dbscanrevisited is that Algorithm 1 never requires storing an $O(n^2)$ distance matrix: streaming range queries consume strictly $O(n)$ auxiliary memory for labels and the BFS queue. When $d$ is low ($d <= 8$) and $epsilon$ yields $O(1)$ average neighborhood cardinality, spatial indices (R\*-trees, KD-trees, Ball-trees) answer each query in $O(log n)$ amortized time, giving $O(n log n)$ overall runtime. In high dimensions (Proposition 3) or dense regimes where $|N_epsilon (p)| = Theta(n)$, spatial pruning collapses and Algorithm 1 reverts to $Theta(n^2 d)$ time. Grid-hashing variants (GriDBSCAN Alg. 48 @mahran2008gridbscan), subsampled core selection (DBSCAN++ Alg. 47 @jang2019dbscanpp), and distributed partition-merge schemes (PDBSCAN Alg. 37 @xu1999pdbscan, MR-DBSCAN Alg. 49 @he2011mrdbscan) address massive-$n$ scaling.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for DBSCAN],
    [Time / Space Complexity], [$O(n log n)$ indexed (low $d$), $O(n^2 d)$ unindexed; $O(n)$ workspace without distance matrix @schubert2017dbscanrevisited.],
    [Geometric Inductive Bias], [Arbitrary non-convex, curved, or nested manifolds of *uniform* local density $lambda = m / (n V_d epsilon^d)$.],
    [Determinism / Ordering], [Core sets & noise strictly unique (Lemma 1); strictly order-invariant iff $m <= 3$ (Prop. 1).],
    [Invariances], [Strictly invariant under $Q in O(d)$ & translation; requires co-scaling $epsilon$ under isotropic scaling.],
    [Outlier Robustness], [Explicit noise rejection; $>= m$ points required to fabricate a cluster (Prop. 2).],
    [Recommended Upgrades], [HDBSCAN (Alg. 32) / OPTICS (Alg. 33) for multi-density; SNN (Alg. 36) for high-$d$; ST-DBSCAN (Alg. 40) for space-time.],
  ),
  caption: [Theoretical guarantees and operational profile of DBSCAN.],
)

#sec("IV", "Engineering and Scientific Applications")

DBSCAN's ability to discover an unknown number of arbitrarily shaped clusters while isolating background clutter makes it the default spatial clustering primitive across physical sciences and engineering.

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Autonomous Driving 3D LiDAR Point-Cloud Segmentation:* Vehicle hulls, cyclists, and roadside structures form non-convex 3D surfaces embedded in sparse atmospheric/ground returns @rusu2011pcl. In $bb(R)^3$, KD-tree range queries execute in $O(n log n)$ real time, and Proposition 2 guarantees that isolated sensor glitches never distort object bounding boxes.
2. *Geospatial Mobility, Urban Crime, and Trajectory Corridors:* Road networks, coastal settlements, and GPS stop clusters follow non-convex ribbons rather than circular Voronoi cells. DBSCAN isolates rare off-route excursions as noise and extends cleanly to spatiotemporal windows via ST-DBSCAN (Alg. 40 @birant2007stdbscan) and trajectory segments @lee2007traclus.
3. *Astrophysics and High-Energy Particle Calorimetry:* Galactic star streams and electromagnetic shower hits form filamentary over-densities atop uniform cosmic-ray or detector noise, directly matching the upper density level-set model $L_f (lambda)$ of Theorem 2.
4. *Unsupervised Industrial Telemetry and Fraud Anomaly Isolation:* Because $cal(X)_"noise"$ requires no prior cluster count $k$ and tolerates arbitrary healthy operating manifolds (multi-modal steady states), isolated failure bursts surface immediately as noise points.
5. *Spatial Transcriptomics and Fluorescence Microscopy:* Tissue micro-niches and subcellular puncta exhibit irregular morphologies separated by low-density extracellular regions, where non-parametric level-set connectivity avoids splitting elongated cell colonies.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

DBSCAN degrades severely under: (a) *heterogeneous cluster densities* (Theorem 1), where no global $epsilon$ suffices; (b) *high ambient dimensionality* ($d > 15$–$20$), where Proposition 3 collapses $epsilon$ sensitivity; (c) *thin high-density bridges* connecting distinct clusters ("chaining"), unless $m$ is raised; and (d) *uncalibrated feature units*, where anisotropic axes distort spherical $epsilon$-balls.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Why DBSCAN Fits / Prerequisite], [Caution / Alternative],
    [3D LiDAR / Robotics], [Non-convex surfaces; fast 3D KD-tree], [Adapt $epsilon$ with range; HDBSCAN (Alg. 32)],
    [Geospatial Hotspots], [Irregular street arcs; noise rejection], [Use ST-DBSCAN (Alg. 40) for time],
    [Multi-Density Retail], [Fails global $epsilon$ if densities vary], [Prefer HDBSCAN (Alg. 32) or VDBSCAN (Alg. 50)],
    [High-D Embeddings], [Distance concentration (Prop. 3)], [Use UMAP + HDBSCAN or SNN (Alg. 36)],
    [Streaming IoT], [Static batch range queries], [Use DenStream or DBSCAN++ (Alg. 47)],
  ),
  caption: [Application suitability and engineering matrix for DBSCAN.],
)

#sec("V", "Structural Relations and Empirical Behavior")

DBSCAN anchors the entire density-based paradigm (Algorithms 31–55) and connects directly to hierarchical linkage via mutual reachability.

#proposition(title: "Exact Cut Equivalence Between DBSCAN\* and Mutual-Reachability Single-Linkage", origin: "Original Unified Proof")[
  Let $d_"core" (p)$ denote the distance from $p in cal(X)$ to its $m$-th nearest neighbor (including $p$), and define the symmetric mutual reachability distance $d_"mreach" (p, q) = max{d_"core" (p), space d_"core" (q), space d(p, q)}$ for $p != q$. For every $epsilon > 0$, removing all edges of weight $> epsilon$ from the complete weighted graph $(cal(X), d_"mreach")$—equivalently cutting the HDBSCAN single-linkage dendrogram at height $epsilon$—yields non-noise components (vertices with $d_"core" (p) <= epsilon$) identical to the core components of DBSCAN\*($epsilon, m$) @campello2015hierarchical.
]

#proof[
  Observe that $|N_epsilon (p)| >= m$ holds if and only if the $m$-th nearest neighbor of $p$ lies at distance $<= epsilon$, i.e., $d_"core" (p) <= epsilon$. Hence for any distinct $p, q in cal(X)$,
  $ d_"mreach" (p, q) <= epsilon &<=> d_"core" (p) <= epsilon and d_"core" (q) <= epsilon and d(p, q) <= epsilon \
    &<=> p in cal(X)_"core" and q in cal(X)_"core" and d(p, q) <= epsilon \
    &<=> {p, q} in E_"core". $
  Thus the subgraph of $(cal(X), d_"mreach")$ retaining edges of weight $<= epsilon$, restricted to vertices satisfying $d_"core" (p) <= epsilon$, is edge-for-edge identical to $G_"core"$. By Lemma 1, its connected components coincide identically with the clusters of DBSCAN\*($epsilon, m$).
]

This exact cut duality explains why HDBSCAN (Alg. 32) and OPTICS (Alg. 33) eliminate DBSCAN's single-scale rigidity by integrating over all $epsilon > 0$, why GDBSCAN (Alg. 39 @sander1998gdbscan) generalizes $N_epsilon$ to arbitrary symmetric predicates, and why RNN-DBSCAN (Alg. 45 @bryant2018rnndbscan) and SNN (Alg. 36 @ertoz2003snn) replace fixed-radius balls with scale-adaptive neighbor counts.

*Empirical Behavior on the ClusterOverview 4-Dataset Suite:* Across the standardized benchmark, DBSCAN achieves its strongest performance on *Multi-Shape* ($"Sil" = 0.803$, $"ARI" = "N/A"$), confirming its geometric superiority over convex centroidal methods (Alg. 1) on interlocking arcs and irregular manifolds. On *Iris* ($"Sil" = 0.618$, $"ARI" = 0.534$), DBSCAN cleanly isolates *Iris setosa* with high Silhouette while merging the density-bridged *versicolor* / *virginica* pair, capping 3-class ARI. On *Customers* ($"Sil" = 0.450$), heterogeneous spending densities directly expose Theorem 1's single-$epsilon$ bottleneck—trailing HDBSCAN ($0.777$) substantially. On *High-D Subspace+Noise* ($"Sil" = 0.352$, $"ARI" = "N/A"$), distance concentration (Proposition 3) limits separation, though explicit noise isolation avoids centroid drift. Its mean runtime of $14.321$ ms reflects unindexed neighborhood scans and scale-search overhead relative to single-pass MST-based HDBSCAN ($2.469$ ms) and OPTICS ($2.182$ ms).
