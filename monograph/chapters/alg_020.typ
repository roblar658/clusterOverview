#import "../template/ieee.typ": *

#alg_entry(
  id: 20,
  name: "COP-KMeans (Constrained K-Means)",
  category: "Cat. I: Centroid & Partitioning",
  year: "2001",
  paradigm: "Instance-Constrained Alternating Partitioning",
  complexity: $O(T k (n d + |E_"ML"| + |E_"CL"|))$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Let $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$ be a finite dataset with index set $V = {1, dots, n}$, and let $k in {2, dots, n}$ denote the target number of clusters. In semi-supervised clustering, domain knowledge is encoded through two undirected graphs on $V$: a *must-link (ML)* edge set $E_"ML" subset.eq binom(V, 2)$ requiring co-assignment, and a *cannot-link (CL)* edge set $E_"CL" subset.eq binom(V, 2)$ forbidding co-assignment @wagstaff2001constrained. Define the feasible partition family
$ cal(F)(E_"ML", E_"CL") := { cal(C) &= {C_1, dots, C_k} : \
  c_i &= c_ell space forall (i, ell) in E_"ML", \
  c_i &!= c_ell space forall (i, ell) in E_"CL" }, $
where $c_i in {1, dots, k}$ denotes the cluster index containing $i$. The *constrained $k$-means problem* seeks $(cal(C), bold(mu)) in cal(F)(E_"ML", E_"CL") times (bb(R)^d)^k$ minimizing the within-cluster sum of squares (WCSS):
$ &min_(cal(C) in cal(F)(E_"ML", E_"CL"), space bold(mu) in (bb(R)^d)^k) J(cal(C), bold(mu)) \
  &quad = min_(cal(C), bold(mu)) sum_(j=1)^k sum_(i in C_j) norm(x_i - mu_j)_2^2. $

#definition(title: "ML Transitive Closure & Quotient Constraint Graph")[
  Because co-membership is an equivalence relation, $G_"ML" = (V, E_"ML")$ partitions $V$ into $m <= n$ connected components (super-instances) $cal(K) = {K_1, dots, K_m}$. The *quotient cannot-link graph* is $H = (cal(K), overline(E)_"CL")$, where $(K_a, K_b) in overline(E)_"CL"$ if and only if there exist $i in K_a$ and $ell in K_b$ with $(i, ell) in E_"CL"$. The constraint system is *logically consistent* if $H$ contains no self-loops, and *$k$-feasible* ($cal(F) != emptyset$) if and only if $H$ is $k$-colorable, i.e., its chromatic number satisfies $chi(H) <= k$.
]

Wagstaff et al. @wagstaff2001constrained extend Lloyd's algorithm (Alg. 1) by replacing unconstrained Voronoi projection with a sequential greedy constraint-checking pass (Algorithm 1), while leaving the Euclidean centroid update step unchanged.

#algbox(
  title: "COP-KMeans (Constrained K-Means)",
  inputs: [Dataset $cal(X) subset.eq bb(R)^d$, integer $k$, constraint sets $(E_"ML", E_"CL")$, initial centroids $bold(mu)^((0))$, max iterations $T$.],
  outputs: [Feasible partition $cal(C) in cal(F)$ and centroids $bold(mu)$, or `INFEASIBLE`.],
)[
  1. *For* $t = 0, 1, dots, T - 1$ *do*:
  2. #h(1em) Initialize empty sets $C_1^((t)), dots, C_k^((t)) <- emptyset$ and labels $c_i <- bot$ for all $i in V$.
  3. #h(1em) *For* each $i = 1, dots, n$ (in visitation order $pi$) *do*:
  4. #h(2em) Sort clusters $j_1, dots, j_k$ by increasing $norm(x_i - mu_j^((t)))_2^2$.
  5. #h(2em) Find first $j in {j_1, dots, j_k}$ passing `Violate(i, j) = false`:
     - No assigned $ell$ with $(i, ell) in E_"ML"$ has $c_ell != j$;
     - No assigned $ell$ with $(i, ell) in E_"CL"$ has $c_ell = j$.
  6. #h(2em) *If* no such $j$ exists, *return* `INFEASIBLE`.
  7. #h(2em) Assign $c_i <- j$ and $C_j^((t)) <- C_j^((t)) union {i}$.
  8. #h(1em) *For* each non-empty $C_j^((t))$, set $mu_j^((t+1)) = |C_j^((t))|^(-1) sum_(i in C_j^((t))) x_i$.
  9. #h(1em) *If* $cal(C)^((t)) = cal(C)^((t-1))$, *break*.
  10. *Return* $(cal(C)^((t)), bold(mu)^((t+1)))$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five foundational results governing COP-KMeans: (A) exact algebraic reduction of must-link components to weighted $k$-means; (B) NP-completeness of constraint feasibility and an explicit dead-end failure counterexample on bipartite graphs; (C) monotone convergence under pure must-link constraints versus strict WCSS non-monotonicity under cannot-link constraints; and (D) loss of Voronoi convexity, permutation dependence, and outlier breakdown.

#subsec("A", "Super-Instance Reduction for Must-Link Components")

#lemma(title: "Exact Must-Link Aggregation Identity", origin: "Original Proof")[
  Let $K_1, dots, K_m$ be the connected components of $G_"ML"$, with weights $w_a = |K_a|$ and barycenters $overline(x)_a = w_a^(-1) sum_(i in K_a) x_i$. For every feasible partition $cal(C) in cal(F)(E_"ML", E_"CL")$ inducing component labels $tilde(c)_a in {1, dots, k}$ and any centroids $bold(mu) in (bb(R)^d)^k$,
  $ J(cal(C), bold(mu)) = sum_(a=1)^m w_a norm(overline(x)_a - mu_(tilde(c)_a))_2^2 + Delta_"ML", $
  where $Delta_"ML" = sum_(a=1)^m sum_(i in K_a) norm(x_i - overline(x)_a)_2^2 >= 0$ is independent of $(cal(C), bold(mu))$.
]

#proof[
  Every $cal(C) in cal(F)(E_"ML", E_"CL")$ is constant on each connected component $K_a$ by induction along ML paths; denote its common cluster label by $tilde(c)_a$. Partitioning the outer sum over $V$ into components $K_1, dots, K_m$ gives
  $ J(cal(C), bold(mu)) &= sum_(a=1)^m sum_(i in K_a) norm((x_i - overline(x)_a) \
    &quad + (overline(x)_a - mu_(tilde(c)_a)))_2^2. $
  Expanding each squared Euclidean norm yields
  $ &sum_(i in K_a) norm(x_i - overline(x)_a)_2^2 + w_a norm(overline(x)_a - mu_(tilde(c)_a))_2^2 \
    &quad + 2 chevron.l sum_(i in K_a) (x_i - overline(x)_a), space overline(x)_a - mu_(tilde(c)_a) chevron.r. $
  Since $sum_(i in K_a) (x_i - overline(x)_a) = sum_(i in K_a) x_i - w_a overline(x)_a = bold(0)$, every cross-term vanishes identically. Summing over $a = 1, dots, m$ produces the stated identity, proving that must-link clustering reduces *without approximation* to weighted $k$-means on ${(overline(x)_a, w_a)}_(a=1)^m$ constrained only by $overline(E)_"CL"$.
]

#subsec("B", "Feasibility Complexity and Greedy Dead-End Failure")

#theorem(title: "NP-Completeness of Constrained Feasibility", origin: "Classical — Cited")[
  Deciding whether $cal(F)(E_"ML", E_"CL") != emptyset$ is NP-complete for every fixed $k >= 3$, even when $E_"ML" = emptyset$ and $d = 1$ @wagstaff2001constrained. Conversely, for $k = 2$, feasibility is decidable in linear time $O(n + |E_"ML"| + |E_"CL"|)$.
]

#proof[
  Membership in NP holds because any proposed label vector $c in {1, dots, k}^n$ is verified against $E_"ML"$ and $E_"CL"$ in $O(|E_"ML"| + |E_"CL"|)$ operations. For NP-hardness when $k >= 3$, reduce from Graph $k$-Colorability: given an undirected graph $G = (V, E)$, place arbitrary coordinates $x_1, dots, x_n in bb(R)$ and set $E_"ML" = emptyset$, $E_"CL" = E$. A cluster assignment $c : V -> {1, dots, k}$ satisfies all constraints if and only if $c(i) != c(ell)$ for every $(i, ell) in E$, i.e., $c$ is a proper $k$-coloring of $G$. Because $k$-Colorability is NP-complete for every $k >= 3$, feasibility of COP-KMeans is NP-complete (independent of the NP-hardness of WCSS optimization itself @aloise2009np @mahajan2012planar). For $k = 2$, contracting $G_"ML"$ into $H = (cal(K), overline(E)_"CL")$ via DFS and testing bipartiteness of $H$ via BFS checks self-loops and odd cycles in $O(n + |E_"ML"| + |E_"CL"|)$ time.
]

Even when $H$ is trivially 2-colorable, Algorithm 1's myopic greedy pass can still fail catastrophically.

#theorem(title: "Greedy Sequential Infeasibility on Bipartite Graphs", origin: "Original Counterexample")[
  For every $k >= 2$, there exists a one-dimensional instance $cal(X) subset.eq bb(R)$ of size $n = 2k - 1$ with $E_"ML" = emptyset$ whose constraint graph $H$ is bipartite ($chi(H) = 2 <= k$) and admits a feasible non-empty $k$-partition, yet Algorithm 1 aborts with `INFEASIBLE` on iteration $t = 0$.
]

#proof[
  Construct $cal(X) subset.eq bb(R)$ with $n = 2k - 1$ distinct points:
  $ x_j = j space (1 <= j <= k), quad x_(k+1) = 0, \
    x_(k+1+r) = 2k + r space (1 <= r <= k - 2). $
  Set $E_"ML" = emptyset$ and define the star graph $E_"CL" = {(j, k + 1) : 1 <= j <= k}$ centered at $k + 1$ (with vertices $k + 2, dots, 2k - 1$ isolated). The constraint graph is a forest of a star $K_(1, k)$ and $k - 2$ isolated vertices, hence strictly bipartite ($chi(H) = 2$). Indeed, the partition
  $ C_1 = {1, dots, k}, quad C_2 = {k + 1}, \
    C_(r+2) = {k + 1 + r} space (1 <= r <= k - 2) $
  is non-empty in every cluster and satisfies every CL edge. Now initialize $mu_j^((0)) = j$ for $j = 1, dots, k$ and execute Algorithm 1 in canonical order $(1, 2, dots, n)$ at $t = 0$:
  - For each $j in {1, dots, k}$, the sole CL neighbor $k + 1$ is still unassigned ($c_(k+1) = bot$). Since $norm(x_j - mu_j^((0)))_2 = 0$, point $x_j$ is assigned to cluster $c_j = j$.
  - At step $i = k + 1$, every cluster $j in {1, dots, k}$ already contains vertex $j$ with $(j, k + 1) in E_"CL"$. Hence `Violate(k + 1, j) = true` for *all* $j in {1, dots, k}$, forcing Algorithm 1 to return `INFEASIBLE` immediately.
]

#subsec("C", "Monotonicity Under ML-Only vs. Failure Under CL Constraints")

#theorem(title: "ML Monotonicity and Strict CL Non-Monotonicity", origin: "Original Proof")[
  1. *(ML-Only Convergence):* If $E_"CL" = emptyset$ and Algorithm 1 assigns each ML super-instance $K_a$ atomically (equivalently via transitive closure), then $J(cal(C)^((t)), bold(mu)^((t)))$ is monotonically non-increasing and terminates in at most $S(m, k)$ iterations.
  2. *(CL Non-Monotonicity):* When $E_"CL" != emptyset$, even along feasible iterations, Algorithm 1 can *strictly increase* both fixed-centroid WCSS and true partition WCSS: $J^*(cal(C)^((t+1))) > J^*(cal(C)^((t)))$, converging to a strictly worse partition.
]

#proof[
  (1) By Lemma 1, for $E_"CL" = emptyset$ the optimal assignment of component $K_a$ minimizes $w_a norm(overline(x)_a - mu_j)_2^2$ independently across $a = 1, dots, m$, while the centroid update sets $mu_j^((t+1)) = (sum_(a in tilde(C)_j) w_a)^(-1) sum_(a in tilde(C)_j) w_a overline(x)_a$, which uniquely minimizes the weighted sum of squares. Every non-trivial step strictly lowers $J^*(cal(C))$, ruling out cycles and capping iterations at the Stirling number $S(m, k)$.

  (2) Take $k = 2$, $d = 1$, $n = 4$ points $x_1 = 3, x_2 = 0, x_3 = 10, x_4 = 11$, single constraint $E_"CL" = {(1, 2)}$, $E_"ML" = emptyset$, and initial centroids $mu_1^((0)) = -2$, $mu_2^((0)) = 6$ under order $(1, 2, 3, 4)$:
  - *Iteration $t = 0$:* $|3 - 6| = 3 < 5 = |3 - (-2)| => c_1^((0)) = 2$. Then $x_2 = 0$ is blocked from cluster 2, giving $c_2^((0)) = 1$. Both $x_3 = 10$ and $x_4 = 11$ lie nearer to $6$ than $-2$, giving $c_3^((0)) = c_4^((0)) = 2$. Thus $C_1^((0)) = {2}$, $C_2^((0)) = {1, 3, 4}$. Updating centroids gives $mu_1^((1)) = 0$ and $mu_2^((1)) = (3 + 10 + 11)/3 = 8$, with partition WCSS
    $ J^*(cal(C)^((0))) = 0^2 + (3 - 8)^2 + (10 - 8)^2 + (11 - 8)^2 = 38. $
  - *Iteration $t = 1$:* Using $bold(mu)^((1)) = (0, 8)$, point $x_1 = 3$ satisfies $|3 - 0| = 3 < 5 = |3 - 8|$, so $c_1^((1)) = 1$. This blocks $x_2 = 0$ from cluster 1, forcing $c_2^((1)) = 2$, while $x_3, x_4$ remain in cluster 2. Thus $C_1^((1)) = {1}$, $C_2^((1)) = {2, 3, 4}$, incurring fixed-centroid cost
    $ J(cal(C)^((1)), bold(mu)^((1))) = 3^2 + 8^2 + 2^2 + 3^2 = 86 > 38. $
    Recomputing centroids gives $mu_1^((2)) = 3$, $mu_2^((2)) = (0 + 10 + 11)/3 = 7$, with new partition WCSS
    $ J^*(cal(C)^((1))) = 0^2 + (0 - 7)^2 + (10 - 7)^2 + (11 - 7)^2 = 74 > 38. $
  - *Iteration $t = 2$:* At $bold(mu)^((2)) = (3, 7)$, $x_1 = 3$ stays in cluster 1, blocking $x_2 = 0$ into cluster 2, and $x_3, x_4$ stay in cluster 2. Algorithm 1 halts at $cal(C)^((1))$ with WCSS $74$—nearly double its earlier feasible cost $38$.
]

#subsec("D", "Geometry, Permutation Dependence, and Breakdown")

#proposition(title: "Non-Convex Cells, Order Non-Invariance, and Breakdown", origin: "Original Proof")[
  1. Feasible optimal partitions in $cal(F)(E_"ML", E_"CL")$ can violate Voronoi convexity: $"conv"(C_1) inter "conv"(C_2) != emptyset$.
  2. For fixed visitation order $pi in S_n$, Algorithm 1 is strictly equivariant under rigid isometries $T(x) = Q x + b$ ($Q in O(d)$), yet is *not* invariant under permutations $pi$ when $E_"CL" != emptyset$.
  3. The finite-sample centroid breakdown point is $epsilon^* = 1 / n$.
]

#proof[
  (1) Place four points $x_1 = 0, x_2 = 1, x_3 = 2, x_4 = 3$ in $bb(R)$ with $E_"ML" = {(1, 3), (2, 4)}$ and $E_"CL" = {(1, 2)}$. The unique feasible 2-partition is $C_1 = {1, 3}$, $C_2 = {2, 4}$, whose convex hulls $[0, 2]$ and $[1, 3]$ intersect on $[1, 2]$.

  (2) Rigid isometries preserve squared norms $norm(T(x_i) - T(mu_j))_2^2 = norm(x_i - mu_j)_2^2$ and commute with cluster means, leaving sorted cluster rankings and constraint checks unchanged for fixed $pi$. For order dependence, revisit Theorem 3(2) at $bold(mu)^((1)) = (0, 8)$: order $pi = (1, 2, 3, 4)$ produces $C_1 = {1}, C_2 = {2, 3, 4}$ (cost 74), whereas $pi' = (2, 1, 3, 4)$ assigns $x_2 = 0$ to cluster 1 first, forces $x_1 = 3$ into cluster 2, and preserves $C_1 = {2}, C_2 = {1, 3, 4}$ (cost 38).

  (3) Sending $norm(x_1)_2 -> infinity$ drives the centroid of $x_1$'s cluster $C_j$ (size $m <= n$) to satisfy $norm(mu_j)_2 >= m^(-1) norm(x_1)_2 - max_(ell >= 2) norm(x_ell)_2 -> infinity$; thus corrupting one point ($1/n$) breaks the estimator.
]

#sec("III", "Complexity, Geometry and Operational Profile")

In each iteration of Algorithm 1, computing all $n k$ squared Euclidean distances costs $O(n k d)$ floating-point operations (vectorizable via BLAS-3 GEMM $X M^T$), and sorting $k$ candidate distances per point costs $O(n k log k)$. Storing constraints as adjacency lists, testing `Violate(i, j)` scans $"deg"_"ML" (i) + "deg"_"CL" (i)$ neighbors; across at most $k$ candidate clusters per point and all $i in V$, constraint checking costs $O(k (|E_"ML"| + |E_"CL"|))$ per iteration. Precomputing ML closure via Union-Find in $O(n + |E_"ML"| alpha(n))$ reduces $n$ to $m <= n$ super-instances (Lemma 1), lowering per-iteration work to $O(m k (d + log k) + k |overline(E)_"CL"|)$. Memory usage is $O((n + k) d + |E_"ML"| + |E_"CL"|)$. Crucially, unlike unconstrained Lloyd/Forgy, Step 3's sequential `Violate` check creates a strict loop-carried dependency across $i = 1, dots, n$, preventing parallel assignment over points.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for COP-KMeans],
    [Time / Space Complexity], [$O(T k (n d + |E_"ML"| + |E_"CL"|))$ time; $O((n+k)d + |E_"ML"| + |E_"CL"|)$ space.],
    [Feasibility Complexity], [NP-complete for $k >= 3$; $O(n + |E_"ML"| + |E_"CL"|)$ for $k = 2$ (Theorem 1).],
    [Descent / Termination], [Monotone for ML-only ($<= S(m, k)$); non-monotone & dead-end prone with CL (Theorems 2–3).],
    [Geometric Inductive Bias], [Piecewise Voronoi warped by hard ML equivalence and CL exclusions (non-convex possible).],
    [Invariances / Breakdown], [Isometry-equivariant for fixed $pi$; permutation-dependent; breakdown $epsilon^* = 1/n$.],
    [Recommended Upgrades], [Precompute ML closure (Lemma 1); random order restarts; soft HMRF / Kernel K-Means (Alg. 16) when CL is noisy.],
  ),
  caption: [Theoretical guarantees and operational profile of COP-KMeans (Alg. 20).],
)

#sec("IV", "Engineering and Scientific Applications")

COP-KMeans excels when sparse, high-confidence domain rules must strictly override raw Euclidean proximity—especially where true semantic classes overlap or violate isotropic spherical geometry.

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *GPS Lane-Finding and Road-Network Trajectory Clustering:* In the original benchmark of Wagstaff et al. @wagstaff2001constrained, sequential GPS samples from a single vehicle pass must belong to the same lane (ML), while simultaneous parallel passes across adjacent lanes cannot merge (CL). Because highway lanes are long, narrow, closely spaced ribbons, unconstrained Voronoi bisectors slice lanes transversely; ML/CL constraints realign partitions longitudinally along physical lanes.
2. *Interactive Medical Image and Pathology Segmentation:* In MRI and histopathology, overlapping voxel intensity distributions cause unconstrained K-Means to conflate adjacent anatomical tissues. Clinician scribble annotations supply sparse ML pairs within a lesion and CL pairs across tissue borders, enabling non-convex separation (Proposition 1(1)) with minimal annotation effort.
3. *Enterprise Entity Resolution and Master-Data Deduplication:* Customer or bibliographic records frequently contain noisy numeric attributes alongside deterministic business rules: shared verified tax/account tokens impose exact ML constraints, while mutually exclusive legal identifiers or incompatible timestamps enforce CL rules that must never be violated.
4. *Multi-Speaker Audio Diarization:* Short contiguous acoustic frames within a single uninterrupted speech turn furnish high-precision ML links, whereas overlapping-speech detectors provide CL edges between simultaneous speakers, stabilizing centroid tracking under reverberation and speaker timbre overlap.
5. *Biological Perturbation and Marker-Guided Cell Typing:* Single-cell expression profiles often exhibit continuous developmental gradients where Euclidean boundaries are arbitrary. Biological replicate batches or curated marker-gene exclusivities supply ML and CL graphs that anchor clusters to verified cell phenotypes.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

COP-KMeans fails severely under (a) *noisy or contradictory constraints*—a single erroneous CL edge inside an ML component makes $H$ self-looped and forces immediate `INFEASIBLE` failure; (b) *dense CL graphs* with $chi(H)$ close to $k$, where greedy sequential assignment triggers Theorem 2 dead-ends; (c) *non-shuffled ordered streams*, where sequential order bias magnifies Theorem 3(2) suboptimality; and (d) *heavy outlier contamination* ($epsilon^* = 1/n$).

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit / Prerequisite], [Caution / Alternative],
    [GPS Lane Mapping], [ML/CL overrides narrow ribbon geometry], [Pre-close ML components],
    [Medical Scribbles], [Sparse human ML/CL lifts Voronoi convexity], [Avoid noisy clicks],
    [Entity Resolution], [Hard regulatory CL/ML compliance], [Use graph solver if dense CL],
    [Noisy Crowdsourcing], [Avoid: hard CL aborts on label noise], [Use soft HMRF / Alg. 16],
    [High-Density Arcs], [Works only if ML spans non-convex arcs], [Use DBSCAN (Alg. 31)],
  ),
  caption: [Domain suitability matrix and failure-mode mitigations for COP-KMeans.],
)

#sec("V", "Structural Relations and Empirical Behavior")

COP-KMeans bridges classical prototype partitioning (Alg. 1) and graph-constrained combinatorial optimization, and arises as the exact infinite-penalty limit of soft-constrained and kernelized clustering.

#proposition(title: "Exact Penalty Equivalence to Soft-Constrained K-Means", origin: "Original Proof")[
  Let $J_1 = sum_(i=1)^n norm(x_i - overline(x))_2^2$ be the global total scatter, and define the unconstrained penalized objective
  $ L_lambda (cal(C)) &:= J^*(cal(C)) + lambda sum_((i, ell) in E_"ML") bb(1)[c_i != c_ell] \
    &quad + lambda sum_((i, ell) in E_"CL") bb(1)[c_i = c_ell]. $
  Whenever $cal(F)(E_"ML", E_"CL") != emptyset$ and $lambda > J_1$, the set of global minimizers of $L_lambda$ over all $k$-partitions coincides *identically* with the set of global minimizers of $J^*$ over $cal(F)(E_"ML", E_"CL")$.
]

#proof[
  For any partition $cal(C) = {C_1, dots, C_k}$, applying the parallel-axis identity to each cluster $C_j$ around the global mean $overline(x)$ yields
  $ J_1 &= sum_(j=1)^k sum_(i in C_j) norm(x_i - overline(x))_2^2 \
    &= J^*(cal(C)) + sum_(j=1)^k |C_j| norm(overline(x)_(C_j) - overline(x))_2^2 >= J^*(cal(C)). $
  Hence $0 <= J^*(cal(C)) <= J_1$ for every partition $cal(C)$. Now let $cal(C)_"feas" in cal(F)(E_"ML", E_"CL")$ be any feasible partition (zero violations) and $cal(C)_"infeas" in.not cal(F)$ be any infeasible partition (at least one violation). Because $lambda > J_1$,
  $ L_lambda (cal(C)_"infeas") &>= 0 + lambda > J_1 \
    &>= J^*(cal(C)_"feas") = L_lambda (cal(C)_"feas"). $
  Thus no infeasible partition can minimize $L_lambda$, and on $cal(F)$ we have $L_lambda (cal(C)) = J^*(cal(C))$, proving exact equivalence.
]

For finite $lambda > 0$, relaxing indicator violations into quadratic Gram adjustments links constrained partitioning directly to weighted Kernel K-Means (Alg. 16 @dhillon2004kernel) and Spectral Clustering (Algs. 103–105), avoiding Theorem 2's hard dead-ends when constraints are noisy.

*Empirical Benchmark Behavior:* On the ClusterOverview 4-dataset benchmark suite (measured using Silhouette @rousseeuw1987silhouettes and Adjusted Rand Index @hubert1985comparing), COP-KMeans attains *Iris* $"Sil" = 0.426$, $"ARI" = 0.564$; *Customers* $"Sil" = 0.300$; *Multi-Shape* $"Sil" = 0.712$ ($"ARI" = "N/A"$); and *High-D Subspace+Noise* $"Sil" = 0.303$ ($"ARI" = "N/A"$), with a mean runtime of *1.031 ms*. These measurements corroborate our theory directly: injecting pairwise constraints yields a decisive advantage over unconstrained batch Forgy (Alg. 22, Multi-Shape $0.609$, High-D $0.238$) and online MacQueen (Alg. 21, Multi-Shape $0.357$, High-D $0.112$) by breaking spherical Voronoi restrictions (Proposition 1(1)), while trading a modest reduction in Iris Euclidean compactness ($0.426$ vs. Forgy's $0.473$) for substantially stronger label alignment ($"ARI" = 0.564$ vs. Forgy's $0.431$). Finally, its $1.031$ ms latency ($7.4 times$ slower than MacQueen and $16.4 times$ slower than Forgy) reflects the serial per-point constraint verification overhead established in Section III.
