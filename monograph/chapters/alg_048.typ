#import "../template/ieee.typ": *

#alg_entry(
  id: 48,
  name: "GriDBSCAN (Grid-Partitioned Accelerated DBSCAN)",
  category: "Cat. II: Density & Mode-Seeking",
  year: "2008",
  paradigm: "Grid-Partitioned Exact Local-DBSCAN Acceleration",
  complexity: $O(n d + sum_g |tilde(P)_g|^2 d)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

A central obstacle to scaling DBSCAN (Alg. 31) @ester1996dbscan is that global $epsilon$-neighborhood queries must inspect points across the entire domain unless mediated by spatial indices whose pointer-chasing and bounding-box overlaps degrade under large memory footprints or parallel execution @schubert2017dbscanrevisited. Introduced by Mahran and Mahar @mahran2008gridbscan, *GriDBSCAN* accelerates DBSCAN via *divide-and-conquer spatial grid partitioning*: it decomposes the data domain into disjoint axis-aligned hyperrectangular grid cells, enlarges each cell by an exact *$epsilon$-halo margin* so that local density calculations are self-contained, executes independent local DBSCAN instances inside each enlarged box, and deterministically merges local clusters across shared boundary points.

Let $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$ lie within the axis-aligned bounding box $Omega = product_(r=1)^d [a_r, b_r]$, and fix DBSCAN parameters $(epsilon, "minPts")$ with $epsilon > 0$ and $"minPts" in bb(N)$. Partition each axis $r in {1, dots, d}$ into $g_r >= 1$ intervals of width $w_r >= epsilon$ (canonically uniform $w_r = (b_r - a_r) / g_r$, or adaptive cost-balanced cuts @mahran2008gridbscan), yielding $M = product_(r=1)^d g_r$ pairwise disjoint cells $G_1, dots, G_M$ covering $Omega$, where
$ G_g = product_(r=1)^d [u_(g, r), space v_(g, r)), quad v_(g, r) - u_(g, r) = w_r >= epsilon. $
To ensure every interior point sees its complete Euclidean $epsilon$-ball locally, each cell $G_g$ is expanded into its *$epsilon$-haloed box*
$ tilde(G)_g := product_(r=1)^d [u_(g, r) - epsilon, space v_(g, r) + epsilon]. $

#definition(title: "Interior Set, Boundary Halo & Local Core Equivalence")[
  For each grid cell $G_g$ ($g = 1, dots, M$), define its *interior point set* $P_g := cal(X) inter G_g$, its *halo point set* $H_g := cal(X) inter (tilde(G)_g backslash G_g)$, and its *enlarged cell dataset* $tilde(P)_g := P_g union H_g = cal(X) inter tilde(G)_g$. Every $x in cal(X)$ owns a unique *home cell* index $h(x) in {1, dots, M}$ satisfying $x in P_(h(x))$. Running local DBSCAN$(epsilon, "minPts")$ on $tilde(P)_g$ yields local core set $C_g^"loc" := {x in tilde(P)_g : |B(x, epsilon) inter tilde(P)_g| >= "minPts"}$ and local core clusters ${L_(g, 1), dots, L_(g, m_g)}$.
]

Because $G_g$ and $tilde(G)_g$ differ strictly by a margin of width $epsilon$ along every coordinate face, local DBSCAN on each non-empty $tilde(P)_g$ can execute completely independently across CPU cores or distributed nodes, followed by a single Union-Find merge pass across overlapping boundary points (Algorithm 1).

#algbox(
  title: "GriDBSCAN (Grid-Partitioned Exact DBSCAN)",
  inputs: [Dataset $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$, radius $epsilon > 0$, threshold $"minPts" >= 1$, grid widths $w_r >= epsilon$ (optimal $w_r = 2 epsilon$).],
  outputs: [Global cluster labels $c: cal(X) -> {1, dots, K} union {-1}$.],
)[
  1. *Grid Partitioning:* Compute bounding box $Omega = product_(r=1)^d [a_r, b_r]$; map each $x in cal(X)$ to home cell $h(x) = (floor((x_1 - a_1) / w_1), dots, floor((x_d - a_d) / w_d))$ to form $P_g$.
  2. *Adaptive Cell Merging (Optional):* Merge adjacent sparse cells whose combined size stays below budget to reduce overhead @mahran2008gridbscan.
  3. *Halo Replication:* For each $x in cal(X)$ and each cell $G_g$ with $x in tilde(G)_g backslash G_g$, append $x$ to halo $H_g$; set $tilde(P)_g = P_g union H_g$.
  4. *Parallel Local DBSCAN:* *For each* non-empty cell $G_g$ *in parallel do*:
  5. #h(1em) Run DBSCAN$(epsilon, "minPts")$ on $tilde(P)_g$, identifying interior local cores $C_g^"int" := C_g^"loc" inter P_g$ and local core clusters $L_(g, 1), dots, L_(g, m_g)$.
  6. #h(1em) For every halo point $y in H_g$, record the *full set* $cal(R)_g (y)$ of local core clusters $L_(g, r)$ containing an interior core $x in C_g^"int"$ with $norm(x - y)_2 <= epsilon$.
  7. *Cross-Cell Union-Find Merging:* Initialize DSU over all local core clusters ${L_(g, r)}$. For every $g$ and every $y in H_g$: if $y in C_(h(y))^"int"$ (i.e., $y$ is core in its home cell), merge $L_(h(y))(y)$ with *every* $L_(g, r) in cal(R)_g (y)$.
  8. *Global Labeling:* Assign each global core $x in union.big_g C_g^"int"$ the DSU root label of $L_(h(x))(x)$. Assign each non-core point $z$ to a merged cluster touching $z$ within $epsilon$ (or nearest core for determinism), else noise ($-1$).
  9. *Return* global labels $c$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five formal theoretical results characterizing GriDBSCAN: (A) exact equivalence of multi-label halo GriDBSCAN to global DBSCAN on all core clusters; (B) an explicit 16-point counterexample proving that naive *single-label* halo merging fails to merge globally density-connected clusters; (C) an exact closed-form expected-work theorem proving $w^* = 2 epsilon$ uniquely minimizes pairwise distance work while exposing the $2^d$–$3^d$ dimensional halo barrier; (D) classical parallel scalability under skew-mitigated partitioning; and (E) exact invariance of core clusters under grid-origin shifts and Euclidean rotations.

#subsec("A", "Exact Core-Graph Equivalence Theorem")

#theorem(title: "Exact Global Core Equivalence of GriDBSCAN", origin: "Original Proof")[
  Let $C_"full" := {x in cal(X) : |B(x, epsilon) inter cal(X)| >= "minPts"}$ and let $G_"core" = (C_"full", E_"core")$ be the global DBSCAN core graph with ${x, y} in E_"core"$ iff $norm(x - y)_2 <= epsilon$. Under Algorithm 1 with $epsilon$-halo boxes $tilde(G)_g$:
  1. *Exact Interior Core Identification:* Every interior point $x in P_g$ satisfies $B(x, epsilon) inter tilde(P)_g = B(x, epsilon) inter cal(X)$, hence $C_g^"int" = C_"full" inter P_g$ and $union.big_(g=1)^M C_g^"int" = C_"full"$.
  2. *Exact Global Component Recovery:* The merged DSU equivalence classes produced by Steps 5–7 induce the *exact* connected components of $G_"core"$, identically matching Ester et al.'s global DBSCAN @ester1996dbscan on all core points.
]

#proof[
  (1) Fix $g$ and $x in P_g = cal(X) inter G_g$, so $u_(g, r) <= x_r < v_(g, r)$ for all $r = 1, dots, d$. Take any $y in B(x, epsilon) inter cal(X)$. Then for every coordinate $r$,
  $ |y_r - x_r| <= norm(y - x)_2 <= epsilon ==> u_(g, r) - epsilon <= y_r <= v_(g, r) + epsilon. $
  Therefore $y in tilde(G)_g inter cal(X) = tilde(P)_g$. Conversely, $tilde(P)_g subset.eq cal(X)$ gives the reverse inclusion, so $B(x, epsilon) inter tilde(P)_g = B(x, epsilon) inter cal(X)$ identically. Counting cardinalities yields $x in C_g^"loc" <==> x in C_"full"$. Because ${P_g}_(g=1)^M$ partitions $cal(X)$, $union.big_g C_g^"int" = C_"full"$.

  (2) Let $tilde.op$ denote the DSU equivalence relation on $C_"full"$ produced by Steps 5–7. We prove $x tilde.op y <==> x, y$ lie in the same connected component of $G_"core"$.
  - *Soundness ($==>$):* Within any local cell $tilde(P)_g$, every local core point $z in C_g^"loc"$ satisfies $|B(z, epsilon) inter cal(X)| >= |B(z, epsilon) inter tilde(P)_g| >= "minPts"$, so $C_g^"loc" subset.eq C_"full"$. Hence every local core-to-core path inside $tilde(P)_g$ is a valid path in $G_"core"$. Likewise, every cross-cell merge in Step 7 joins $L_(h(y))(y)$ (with $y in C_(h(y))^"int" subset.eq C_"full"$) to $L_(g, r)$ only when there exists $x in C_g^"int" subset.eq C_"full"$ with $norm(x - y)_2 <= epsilon$, i.e., ${x, y} in E_"core"$. By transitivity, $x tilde.op y$ implies $x$ and $y$ are connected in $G_"core"$.
  - *Completeness ($<==$):* Take any edge ${x, y} in E_"core"$ with $x, y in C_"full"$ and $norm(x - y)_2 <= epsilon$. Let $g = h(x)$ and $g' = h(y)$. If $g = g'$, both $x, y in C_g^"int"$ are adjacent inside $tilde(P)_g$, hence share the same local cluster $L_g (x) = L_g (y)$. If $g != g'$, then by part (1), $y in B(x, epsilon) inter cal(X) subset.eq tilde(P)_g$, and since $y in.not P_g$ we have $y in H_g$. Because $x in C_g^"int"$ and $norm(x - y)_2 <= epsilon$, Step 6 places $L_g (x) in cal(R)_g (y)$. Since $y in C_(g')^"int"$, Step 7 explicitly unites $L_(g')(y)$ with $L_g (x)$, giving $x tilde.op y$. By induction along edges of $G_"core"$, every globally connected pair satisfies $x tilde.op y$.
]

#subsec("B", "Necessity of Multi-Label Halo Merging")

A subtle implementation flaw arises if local DBSCAN stores only a *single* cluster label for each halo border point $y in H_g$ (as in standard single-pass DBSCAN) and merges cells solely via matched single labels.

#proposition(title: "Failure Counterexample for Naive Single-Label Halo Merging", origin: "Original Counterexample")[
  There exists a 16-point planar dataset ($d = 2$, $epsilon = 1$, $"minPts" = 4$) split by a single grid hyperplane $x_1 = 0$ on which global DBSCAN produces $K = 1$ connected cluster, whereas naive single-label halo merging splits the dataset into $K = 2$ disconnected clusters.
]

#proof[
  Set $d = 2$, $epsilon = 1$, $"minPts" = 4$, left cell $G_L = {x_1 < 0}$, right cell $G_R = {x_1 >= 0}$, with enlarged boxes $tilde(G)_L = {x_1 <= 1}$ and $tilde(G)_R = {x_1 >= -1}$. Construct four 4-point groups:
  - $A_L$: $x = (-0.4, 0)$ plus 3 companions $(-1.2, 0 plus.minus 0.05), (-1.2, 0)$;
  - $B_L$: $z = (-0.4, 1.15)$ plus 3 companions $(-1.2, 1.15 plus.minus 0.05), (-1.2, 1.15)$;
  - $A_R$: $y = (0.4, 0.575)$ plus 3 companions $(1.2, 0.575 plus.minus 0.05), (1.2, 0.575)$;
  - $B_R$: $w = (0.4, 1.725)$ plus 3 companions $(1.2, 1.725 plus.minus 0.05), (1.2, 1.725)$.
  Every companion lies within $sqrt(0.8^2 + 0.05^2) approx 0.802 < 1$ of its anchor ($x, z, y, w$), so all four anchors are global core points. Meanwhile, $norm(x - z)_2 = norm(y - w)_2 = 1.15 > 1$, inter-group companion distances exceed $1.05 > 1$, and all cross-boundary anchor distances satisfy
  $ norm(x - y)_2 = norm(z - y)_2 = norm(z - w)_2 = sqrt(0.8^2 + 0.575^2) approx 0.985 < 1. $
  Thus global DBSCAN merges $x - y - z - w$ into a *single* cluster ($K = 1$). Now inspect both local boxes:
  - In $tilde(P)_L = A_L union B_L union {y, w}$ (all $A_R, B_R$ companions lie at $x_1 = 1.2 > 1$): $x$ (core, cluster $L_1$) and $z$ (core, cluster $L_2$) are disconnected locally. Halo point $y$ has only 3 local neighbors ${y, x, z} < 4$, making it a local *border* point; if $L_1$ expands first, $y$ receives single label $L_1$ (while $w$ receives $L_2$).
  - In $tilde(P)_R = A_R union B_R union {x, z}$ (all $A_L, B_L$ companions lie at $x_1 = -1.2 < -1$): $y$ (core, cluster $R_1$) and $w$ (core, cluster $R_2$) are disconnected locally. Halo point $z$ has only 3 local neighbors ${z, y, w} < 4$, making it a local *border* point; if $R_2$ expands first, $z$ receives single label $R_2$ (while $x$ receives $R_1$).
  Hence the four shared points yield only pairs $x |-> (L_1, R_1)$, $y |-> (L_1, R_1)$, $z |-> (L_2, R_2)$, $w |-> (L_2, R_2)$, merging into *two* disjoint components ${L_1, R_1}$ and ${L_2, R_2}$—missing the global core edge ${z, y}$! Recording the full set $cal(R)_g (y) = {L_1, L_2}$ (Step 6 of Algorithm 1) merges $R_1$ with both $L_1$ and $L_2$, restoring exactness.
]

#subsec("C", "Optimal Grid Width w* = 2*eps & Dimensional Halo Barrier")

#theorem(title: "Exact Expected Distance Work & Optimal Grid Spacing", origin: "Original Derivation")[
  Let $cal(X) = {x_1, dots, x_n}$ be i.i.d. uniform on the flat $d$-torus $(bb(R) / L bb(Z))^d$ (or an upper bound on $[0, L]^d$), partitioned into $M = (L / w)^d$ regular grid cells of side length $w >= epsilon$ ($L / w in bb(N)$). Then:
  1. *Pointwise Halo Replication:* Every point $x in cal(X)$ belongs to at most $3^d$ enlarged boxes $tilde(G)_g$, and at most $2^d$ enlarged boxes whenever $w > 2 epsilon$ (or $w >= 2 epsilon$ for almost every $x$).
  2. *Exact Expected Pairwise Work:* The expected total number of local pairwise distance evaluations across all $M$ cells is exactly
     $ bb(E)[W(w)] = sum_(g=1)^M bb(E)[binom(|tilde(P)_g|, 2)] = binom(n, 2) 1 / L^d (w + 4 epsilon + (4 epsilon^2) / w)^d. $
  3. *Unique Optimal Cell Width:* Minimizing $bb(E)[W(w)]$ over $w > 0$ yields the *unique* optimal grid spacing $w^* = 2 epsilon$, attaining
     $ bb(E)[W(w^*)] = binom(n, 2) ((8 epsilon) / L)^d. $
     Consequently, uniform GriDBSCAN achieves exponential speedup $(L / (8 epsilon))^d$ when $8 epsilon < L$, but suffers *negative speedup* (more work than unpartitioned DBSCAN) whenever $8 epsilon >= L$ or $2^d$ dominates cell partitioning.
]

#proof[
  (1) Along any single axis $r$, a coordinate $x_r$ lies in the enlarged interval $[k w - epsilon, (k + 1) w + epsilon]$ iff $(x_r - epsilon) / w - 1 <= k <= (x_r + epsilon) / w$. This interval of real numbers has length $1 + 2 epsilon / w$. If $w > 2 epsilon$, the length is strictly less than $2$, so at most $2$ integers $k$ lie inside it; taking the Cartesian product over $d$ axes bounds replication by $2^d$. If $epsilon <= w <= 2 epsilon$, the length lies in $[2, 3]$, admitting at most $3$ cells per axis (for interior coordinates) and at most $3^d$ enlarged boxes total.

  (2) Each enlarged box $tilde(G)_g$ has side length $w + 2 epsilon$ and volume fraction $p_g = ((w + 2 epsilon) / L)^d$ (assuming $w + 2 epsilon <= L$). Under i.i.d. uniform draws, $|tilde(P)_g| tilde.op "Binomial"(n, p_g)$, whose second factorial moment is $bb(E)[|tilde(P)_g|(|tilde(P)_g| - 1)] = n(n - 1) p_g^2$. Summing $1 / 2 bb(E)[|tilde(P)_g|(|tilde(P)_g| - 1)]$ over all $M = (L / w)^d$ boxes gives
  $ bb(E)[W(w)] &= M dot binom(n, 2) ((w + 2 epsilon) / L)^(2 d) = binom(n, 2) 1 / L^d ((w + 2 epsilon)^2 / w)^d \
                &= binom(n, 2) 1 / L^d (w + 4 epsilon + (4 epsilon^2) / w)^d. $

  (3) Since $u |-> u^d$ is strictly increasing on $(0, infinity)$, minimizing $bb(E)[W(w)]$ is equivalent to minimizing $h(w) := w + 4 epsilon + 4 epsilon^2 / w$. By the AM–GM inequality,
  $ w + (4 epsilon^2) / w >= 2 sqrt(w dot (4 epsilon^2) / w) = 4 epsilon, $
  with unique equality at $w = 4 epsilon^2 / w <==> w^* = 2 epsilon$. Substituting $w^* = 2 epsilon$ yields $h(w^*) = 8 epsilon$, so $bb(E)[W(w^*)] = binom(n, 2) (8 epsilon / L)^d$. Whenever $8 epsilon >= L$, the ratio $(8 epsilon / L)^d >= 1$, proving that halo duplication exceeds full unpartitioned pairwise work.
]

#subsec("D", "Adaptive Partitioning, Parallel Span & Geometric Invariance")

#theorem(title: "Skew Bottleneck & Parallel Scalability of Grid DBSCAN", origin: "Classical — Cited")[
  Let $P$ parallel workers execute Step 4 of Algorithm 1 over $M >= P$ grid boxes with maximum box occupancy $n_"max" := max_g |tilde(P)_g|$ and total replicated point count $N_"rep" := sum_g |tilde(P)_g|$ @mahran2008gridbscan @he2011mrdbscan. Then:
  1. The parallel wall-clock time satisfies $T_"par" = O(n_"max"^2 d + N_"rep" alpha(n))$ using brute-force local scans (or $O(n_"max" log n_"max" + N_"rep" alpha(n))$ with local $R^*$-trees).
  2. Under uniform grids on skewed densities where a single dense cluster concentrates $Theta(n)$ points in one cell, $n_"max" = Theta(n)$ and parallel speedup collapses to $O(1)$; adaptive cost-balanced grid splitting @mahran2008gridbscan @he2011mrdbscan restores near-linear speedup whenever $n_"max" = O(n / P)$ and boundary replication $N_"rep" / n = O(1)$.
]

#proof(label: "Proof Synthesis (Cited — Mahran & Mahar 2008; He et al. 2011)")[
  Because local DBSCAN executions across boxes are mutually independent, the parallel span of Step 4 is governed strictly by the heaviest box $max_g T_"loc" (|tilde(P)_g|)$, while the sequential Union-Find merge pass processes at most $O(N_"rep")$ halo-incidence records in $O(N_"rep" alpha(n))$ time @mahran2008gridbscan. If a single cell contains $c n$ points ($c > 0$), that single worker incurs $Theta(c^2 n^2 d)$ time regardless of $P$, bounding speedup by $1 / c^2 = O(1)$. Recursively splitting overloaded cells and merging adjacent sparse cells via marginal projection histograms @mahran2008gridbscan or cost-based $k$-d partitioning @he2011mrdbscan equalizes $n_"max"$ while keeping boundary surface-to-volume ratios small in low $d$, recovering near-linear scalability.
]

#proposition(title: "Grid-Shift & Rotation Invariance of Exact Core Clusters", origin: "Original Proof")[
  Unlike lossy grid clustering algorithms (STING Alg. 133 @wang1997sting; CLIQUE Alg. 134 @agrawal1998clique), GriDBSCAN's core-point partition under Algorithm 1 is *strictly invariant* under both (i) arbitrary grid-origin shifts $t in bb(R)^d$, and (ii) arbitrary rigid Euclidean isometries $T(x) = Q x + b$ ($Q in O(d)$)—even though the cell occupancies $|tilde(P)_g|$ and execution runtime vary with $t$ and $Q$.
]

#proof[
  By Theorem 1(2), Algorithm 1 recovers the exact connected components of the global core graph $G_"core" (cal(X), epsilon, "minPts")$, which depends solely on pairwise Euclidean distances $norm(x_i - x_j)_2$. Rigid isometries satisfy $norm(T(x_i) - T(x_j))_2 = norm(x_i - x_j)_2$, and shifting the auxiliary grid origin by $t$ does not modify $cal(X)$ or $G_"core"$. Conversely, rotating a 1D line of $n$ points along the $x_1$-axis by $45 degree$ in $bb(R)^2$ spreads it across diagonal cells and increases the number of halo-crossing faces, demonstrating that $|tilde(P)_g|$ and runtime are rotation-dependent while the core output remains invariant.
]

#sec("III", "Complexity, Geometry and Operational Profile")

GriDBSCAN replaces global $O(n^2 d)$ scans with $M$ smaller local problems costing $O(sum_g |tilde(P)_g|^2 d)$ (or $O(sum_g |tilde(P)_g| log |tilde(P)_g|)$ with local spatial trees), plus $O(n d)$ grid hashing and $O(N_"rep" alpha(n))$ merging. In 2D–3D spatial workloads where $L / (8 epsilon) >> 1$ and $2^d <= 8$, Theorem 2 delivers $10^2 times$–$10^4 times$ distance reductions with trivial cache locality—each box $tilde(P)_g$ fits cleanly inside L2/L3 cache! However, as $d$ grows beyond $6$–$10$, the $2^d$–$3^d$ halo replication factor causes $N_"rep"$ to explode exponentially, making grid indexing slower than brute-force or subsampled DBSCAN++.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for GriDBSCAN],
    [Time / Space Complexity], [$O(n d + sum_g |tilde(P)_g|^2 d)$; optimal $w^* = 2 epsilon$ gives $(8 epsilon / L)^d$ ratio; $O(2^d n d)$ halo space.],
    [Geometric Inductive Bias], [Exact arbitrary non-convex DBSCAN density level sets with explicit noise isolation.],
    [Determinism / Seeding], [Core clusters strictly deterministic and grid-invariant (Prop. 2); zero random seeding.],
    [Invariances], [Core clusters invariant to rigid isometries & grid shifts; runtime varies with rotation & scale.],
    [Outlier Robustness], [High: sparse outliers isolated in low-occupancy cells and labeled noise ($-1$).],
    [Recommended Companions], [MR-DBSCAN (Alg. 49) for distributed clusters; DBSCAN++ (Alg. 47) for high-$d$; HDBSCAN (Alg. 32) for multi-density.],
  ),
  caption: [Summary of theoretical guarantees and operational properties of GriDBSCAN.],
)

#sec("IV", "Engineering and Scientific Applications")

GriDBSCAN is ideally suited to low-dimensional ($d in {2, 3, 4}$), large-scale physical and geospatial workloads that demand *exact* DBSCAN semantics alongside multi-core or shared-nothing parallelism.

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Geographic Information Systems (GIS) & Earth Observation:* Satellite wildfire detections, seismic catalogs, and urban POI databases live in $d = 2$ or $3$ with $L / (8 epsilon) in [10, 10^3]$. By Theorem 2, setting $w^* = 2 epsilon$ yields $10^2 times$–$10^6 times$ pairwise distance reduction while preserving exact coastal and fault-line cluster boundaries without grid-quantization artifacts (Proposition 2).
2. *Multi-Core Real-Time 3D Point-Cloud & Lidar Processing:* Autonomous perception pipelines partition 3D voxel spaces ($d = 3$, replication $<= 2^3 = 8$) across multi-core CPUs/GPUs; each enlarged tile runs cache-resident local DBSCAN independently, meeting strict sub-50 ms frame budgets with exact obstacle clustering.
3. *Materials Science & Molecular Dynamics Defect Detection:* Identifying vacancy clusters, dislocation loops, and crystal grain boundaries in 3D atomic lattices involves $10^6$–$10^8$ atoms with uniform short-range interaction cutoff $epsilon$. Regular $w = 2 epsilon$ spatial decomposition mirrors classic Verlet cell lists, achieving $O(n)$ wall-clock scaling.
4. *Maritime AIS & Aviation Trajectory Corridor Mining:* Millions of 2D/3D vessel or aircraft positions form dense shipping lanes amid sparse open-ocean outliers. Adaptive cell merging (Step 2) collapses empty ocean tiles while keeping dense harbor cells small, achieving balanced multi-threaded throughput.
5. *Biomedical Spatial Transcriptomics & Histopathology Nuclei Clustering:* Whole-slide pathology images contain $10^6+$ cell centroids in 2D tissue coordinates. GriDBSCAN's exact halo merging guarantees that tumor microenvironments crossing tile boundaries are never artificially severed—solving a notorious tiling artifact of naive sliding-window image analysis.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

GriDBSCAN encounters four primary failure modes: (a) *exponential $2^d$–$3^d$ halo explosion* in high dimensions ($d >= 8$), where Theorem 2 proves halo duplication destroys speedup; (b) *small-$n$ indexing overhead*, where grid hashing, halo copying, and merge tables cost more than direct $O(n^2)$ scans; (c) *severe spatial skew* under non-adaptive uniform grids (Theorem 3); and (d) *naive single-label merge bugs* (Proposition 1) if multi-label halo reachability is omitted.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Why GriDBSCAN Fits / Prerequisite], [Caution / Alternative],
    [GIS / Spatial DBs ($d <= 3$)], [Exact DBSCAN; massive $(L/(8 epsilon))^d$ speedup], [Use adaptive merge for dense cities],
    [3D LiDAR / Atomic Lattices], [Verlet-like $w = 2 epsilon$ cache-resident tiles], [Avoid $w < epsilon$ (breaks halo)],
    [Whole-Slide Pathology], [Seamless cross-tile exact core stitching], [Must use multi-label halo merge],
    [High-D Features ($d >= 10$)], [$2^d$–$3^d$ halo blow-up negates speedup], [Use DBSCAN++ (Alg. 47) or SNN (Alg. 36)],
    [Petabyte Distributed Data], [Single-node memory cap on merge map], [Use MR-DBSCAN (Alg. 49) or PDBSCAN (Alg. 37)],
  ),
  caption: [Application suitability and engineering trade-off matrix for GriDBSCAN.],
)

#sec("V", "Structural Relations and Empirical Behavior")

GriDBSCAN occupies a crucial bridge between *exact* density clustering (DBSCAN Alg. 31 @ester1996dbscan, PDBSCAN Alg. 37 @xu1999pdbscan, MR-DBSCAN Alg. 49 @he2011mrdbscan) and *approximate* grid-based methods (STING Alg. 133 @wang1997sting, CLIQUE Alg. 134 @agrawal1998clique, WaveCluster Alg. 135 @sheikholeslami1998wavecluster, OptiGrid Alg. 141 @hinneburg1999optigrid).

#proposition(title: "Halo-Zero Reduction to Grid Clustering & Distributed Generalization", origin: "Original Unified Proof")[
  Let $eta in [0, epsilon]$ parameterize the halo margin $tilde(G)_g (eta) := product_(r=1)^d [u_(g, r) - eta, v_(g, r) + eta]$ in GriDBSCAN:
  1. At $eta = epsilon$, Algorithm 1 achieves zero core-classification error relative to global DBSCAN (Theorem 1).
  2. Shrinking $eta -> 0$ eliminates cross-cell replication ($N_"rep" = n$) but misclassifies interior boundary strips of volume fraction $1 - (1 - 2 epsilon / w)^d$; replacing local DBSCAN in $G_g$ by cell-count thresholding $|P_g| >= tau$ with adjacent-cell union reduces directly to axis-aligned grid density clustering (CLIQUE / STING paradigm) @agrawal1998clique @wang1997sting, sacrificing rotation and shift invariance.
  3. Replacing uniform grid cells with cost-balanced $k$-d/BSP partitions and MapReduce shuffle-based halo replication generalizes Algorithm 1 directly to MR-DBSCAN (Alg. 49) @he2011mrdbscan.
]

#proof[
  (1) Follows directly from Theorem 1. (2) At $eta = 0$, $H_g = emptyset$ and $tilde(G)_g (0) = G_g$. Only points in the shrunken inner core $product_(r=1)^d [u_(g, r) + epsilon, v_(g, r) - epsilon]$—occupying volume fraction $((w - 2 epsilon) / w)^d = (1 - 2 epsilon / w)^d$ (for $w >= 2 epsilon$)—retain their complete Euclidean $epsilon$-ball inside $G_g$; points in the complementary boundary shell of fraction $1 - (1 - 2 epsilon / w)^d$ suffer truncated $epsilon$-balls, causing false core dropouts and severed clusters unless cell-level adjacency replaces pointwise $epsilon$-balls (as in CLIQUE/STING). (3) Every step of Theorem 1 requires only that each partition region $G_g$ be convex/axis-aligned and enlarged by Minkowski sum $G_g plus.o [-epsilon, epsilon]^d$, which holds verbatim for any cost-balanced $k$-d/BSP partition used in MR-DBSCAN @he2011mrdbscan.
]

GriDBSCAN also composes naturally with generalized neighborhood predicates (GDBSCAN Alg. 39 @sander1998gdbscan) and spatiotemporal halos (ST-DBSCAN Alg. 40 @birant2007stdbscan) by expanding spatial and temporal box faces by $(epsilon_s, epsilon_t)$ respectively.

*Empirical Behavior on ClusterOverview Benchmarks:* On the 4-dataset ClusterOverview benchmark suite @rousseeuw1987silhouettes @hubert1985comparing, GriDBSCAN delivers *exceptionally high clustering quality* across every geometry—precisely because exact halo stitching preserves true DBSCAN density connectivity without subsampling artifacts. On *Multi-Shape*, GriDBSCAN attains an outstanding *Sil = 0.806* (ARI = N/A), cleanly recovering all four non-convex and elongated structures while rejecting peripheral outliers; on *Customers*, it achieves *Sil = 0.663* by isolating sparse peripheral records as noise; on *Iris*, it scores *Sil = 0.618* and *ARI = 0.534*; and on *High-D Subspace+Noise* ($d = 12$), it reaches *Sil = 0.352* (ARI = N/A), edging out both DBSCAN++ (`0.348`) and DeBaCl (`0.217`). Tellingly, however, GriDBSCAN's mean runtime is *9.027 ms*—nearly $5 times$ slower than DeBaCl (`1.836 ms`) and $6.7 times$ slower than DBSCAN++ (`1.346 ms`)—empirically validating Theorem 2 and Section IV-B: on modest dataset sizes ($n <= 240$) and 12-dimensional spaces, grid construction, halo replication, and cross-cell merging overhead dominate local distance savings, underscoring that GriDBSCAN's runtime advantage materializes at large $n$ in low spatial dimensions ($d <= 4$).
