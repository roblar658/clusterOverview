#import "../template/ieee.typ": *

#alg_entry(
  id: 49,
  name: "MR-DBSCAN: MapReduce-Based Density-Based Spatial Clustering",
  category: "Cat. II: Density & Mode-Seeking",
  year: "2011",
  paradigm: "Distributed Cost-Partitioned Density Reachability",
  complexity: $cal(O)((n / P) log (n / P) + m_"ov" alpha(m_"ov"))$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Standard DBSCAN @ester1996dbscan discovers arbitrarily shaped clusters as maximal density-connected components of a finite dataset $cal(X) = {x_1, dots, x_n} subset bb(R)^d$ under a metric $d(dot, dot)$, a neighborhood radius $epsilon > 0$, and a minimum point threshold $"MinPts" in bb(N)$ with $"MinPts" >= 2$. While low-dimensional spatial indices such as R\*-trees reduce sequential runtime from $cal(O)(n^2)$ to $cal(O)(n log n)$ @schubert2017dbscanrevisited, DBSCAN remains inherently sequential because cluster expansion performs a breadth-first traversal of density-reachable core points across a global spatial index that cannot fit into single-node memory at web scale.

Early parallelizations such as PDBSCAN @xu1999pdbscan relied on distributed R\*-trees with master-worker RPC communication, which suffered from synchronization bottlenecks and poor fault tolerance. Uniform grid-partitioning schemes such as GriDBSCAN @mahran2008gridbscan decoupled local clustering across workers, but equal-area spatial cuts severely underperform on skewed spatial distributions where a single urban core or dense manifold concentrates most computational load on one reducer. MR-DBSCAN @he2011mrdbscan resolves both bottlenecks within the MapReduce paradigm through a four-stage architecture: (i) fine-grid census and Cost-Based Partitioning (CBP) via Partition with Reduced Boundary Points (PRBP), (ii) $epsilon$-margin data replication across $P$ reducers, (iii) embarrassingly parallel local DBSCAN execution on each extended subspace, and (iv) global core-gated equivalence merging and distributed relabeling.

Let $N_epsilon (x) = {y in cal(X) : d(x, y) <= epsilon}$ denote the closed $epsilon$-neighborhood of $x in cal(X)$. A point $x$ is a *core point* (written $"core"(x) = 1$) if $|N_epsilon (x)| >= "MinPts"$, and we write $cal(X)_"core" = {x in cal(X) : "core"(x) = 1}$. The global *core graph* $G_"core" (cal(X)) = (cal(X)_"core", E_"core")$ joins $x, y in cal(X)_"core"$ iff $d(x, y) <= epsilon$. A non-core point $b in cal(X) without cal(X)_"core"$ is a *border point* if $N_epsilon (b) inter cal(X)_"core" != emptyset$, and *noise* otherwise.

#definition(title: "Margin-Extended Partition and Core-Bridge Graph")[
  Let $Omega = product_(j=1)^d [a_j, b_j] supset cal(X)$ be a bounding hyperrectangle partitioned into $P$ pairwise interior-disjoint axis-aligned hyperrectangles $cal(R) = {R_1, dots, R_P}$ with half-open ownership rules so that every $x in cal(X)$ has a unique *home partition* $h(x) in {1, dots, P}$ with $x in R_(h(x))$. For each $r in {1, dots, P}$ with $R_r = product_(j=1)^d [a_(r,j), b_(r,j)]$, define the *$epsilon$-extended region* and *$epsilon$-margin* by
  $
    R_r^(+epsilon) &= product_(j=1)^d [a_(r,j) - epsilon, b_(r,j) + epsilon], \
    M_r &= R_r^(+epsilon) without R_r.
  $
  The local dataset routed to reducer $r$ is $cal(X)_r = cal(X) inter R_r^(+epsilon)$, split into inner points $cal(X)_r^"in" = cal(X) inter R_r$ and margin points $cal(X)_r^"mg" = cal(X) inter M_r$. Let $C_(r, c) subset cal(X)_r$ denote local cluster $c$ produced by running DBSCAN on $cal(X)_r$ with local core indicator $"core"_r (x) = upright(I) [|N_epsilon^(r) (x)| >= "MinPts"]$, where $N_epsilon^(r) (x) = N_epsilon (x) inter cal(X)_r$. The *global core-bridge graph* $cal(H)_"merge"$ has vertices $(r, c)$ and an undirected edge between $(r, c)$ and $(r', c')$ iff there exists a shared point $z in C_(r, c) inter C_(r', c')$ satisfying $"core"_r (z) = 1$ or $"core"_(r') (z) = 1$.
]

To balance skewed datasets before launching local clustering, MR-DBSCAN discretizes $Omega$ into uniform grid cells $cal(G)$ of side length $w >= 2 epsilon$ and assigns each cell $g in cal(G)$ an empirical computation cost $c_g = phi(n_g)$ (e.g., $n_g log n_g$ with R\*-tree indexing or $n_g^2$ for brute-force local search, where $n_g = |cal(X) inter g|$). PRBP recursively bisects a region $R$ along the grid-aligned hyperplane $H$ orthogonal to coordinate $j in {1, dots, d}$ that minimizes the cost-imbalance objective
$
  J(R; j, t) = |Phi(R_L (j, t)) - Phi(R_R (j, t))| + lambda B(R; j, t),
$
where $Phi(R) = sum_(g subset R) c_g$ is total estimated compute work, $B(R; j, t)$ counts points within distance $epsilon$ of the split hyperplane $x_j = t$, and $lambda >= 0$ trades boundary replication against load balance.

#algbox(
  title: "MR-DBSCAN: Distributed Cost-Partitioned DBSCAN",
  inputs: [Dataset $cal(X) subset bb(R)^d$, radius $epsilon > 0$, threshold $"MinPts"$, partition count $P = 2^k$, grid cell width $w >= 2 epsilon$.],
  outputs: [Global cluster labels $ell: cal(X) -> {0, 1, dots, K}$ (with $0$ denoting noise).],
)[
  1. *MapReduce Job 1 (Grid Census & PRBP):* Map each $x_i$ to its grid cell $g(x_i) in cal(G)$; reduce to aggregate cell counts $n_g$ and costs $c_g = phi(n_g)$. Run recursive PRBP bisection on the compact grid summary to construct $P$ hyperrectangles ${R_1, dots, R_P}$ and their $epsilon$-extensions ${R_r^(+epsilon)}$.
  2. *MapReduce Job 2 (Margin Replication & Local DBSCAN):* In the Map phase, emit key-value pair $(r, (x_i, upright(I) [x_i in R_r]))$ for every $r$ such that $x_i in R_r^(+epsilon)$. In Reducer $r$, build a local R\*-tree over $cal(X)_r$ and execute exact sequential DBSCAN$(cal(X)_r, epsilon, "MinPts")$ to obtain local labels $(r, c_r (x_i))$ and flags $"core"_r (x_i)$. Emit boundary records $(x_i, (r, c_r (x_i), "core"_r (x_i)))$ for every $x_i$ falling inside any overlapping margin.
  3. *Global Core-Bridge Merge:* Group emitted boundary records by point ID $x_i$. Whenever $x_i$ appears in local clusters $(r, c)$ and $(r', c')$ with $"core"_r (x_i) = 1$ or $"core"_(r') (x_i) = 1$, union $(r, c)$ and $(r', c')$ in a Disjoint-Set Union (DSU) structure to build the global relabeling map $pi(r, c) in {1, dots, K}$.
  4. *MapReduce Job 3 (Global Relabeling):* Map each home-partition record $(x_i, (h(x_i), c))$ through $pi$ to assign $ell(x_i) = pi(h(x_i), c)$ if clustered, or $ell(x_i) = 0$ if noise.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish four core mathematical results governing MR-DBSCAN: exact interior neighborhood preservation and global core-cluster isomorphism, an explicit counterexample proving why core-gating during merging is indispensable, tight worst-case and expected boundary-replication communication bounds, and a provable parallel speedup bound for PRBP cost bisection.

#subsec("A", "Exact Core-Component Isomorphism")

#lemma(title: "Exact Interior Neighborhood Preservation", origin: "Original Proof")[
  Let $d(x, y) = ||x - y||_p$ for any $1 <= p <= infinity$. For every partition $r in {1, dots, P}$ and every inner point $x in cal(X)_r^"in"$, the local neighborhood in $cal(X)_r$ coincides identically with the global neighborhood: $N_epsilon^(r) (x) = N_epsilon (x)$ and $"core"_r (x) = "core"(x)$. Conversely, for every margin point $y in cal(X)_r^"mg"$, $N_epsilon^(r) (y) subset.eq N_epsilon (y)$ and $"core"_r (y) <= "core"(y)$.
]
#proof[
  Fix $x in cal(X)_r^"in" = cal(X) inter R_r$, so $a_(r,j) <= x_j <= b_(r,j)$ for all $j in {1, dots, d}$. Take any global neighbor $y in N_epsilon (x)$, meaning $||y - x||_p <= epsilon$. By monotonicity of $ell_p$ norms, for every coordinate $j in {1, dots, d}$ we have
  $
    |y_j - x_j| <= ||y - x||_infinity <= ||y - x||_p <= epsilon.
  $
  Consequently, $a_(r,j) - epsilon <= x_j - epsilon <= y_j <= x_j + epsilon <= b_(r,j) + epsilon$ for every $j$, which implies $y in R_r^(+epsilon)$ and thus $y in cal(X)_r$. Hence $N_epsilon (x) subset.eq N_epsilon^(r) (x)$. The reverse inclusion $N_epsilon^(r) (x) subset.eq N_epsilon (x)$ holds trivially since $cal(X)_r subset.eq cal(X)$. Equality of neighborhoods forces $|N_epsilon^(r) (x)| = |N_epsilon (x)|$, so $"core"_r (x) = "core"(x)$. For any margin point $y in cal(X)_r^"mg"$, $cal(X)_r subset.eq cal(X)$ gives $N_epsilon^(r) (y) subset.eq N_epsilon (y)$, hence $|N_epsilon^(r) (y)| <= |N_epsilon (y)|$ and $"core"_r (y) <= "core"(y)$.
]

#theorem(title: "Core-Cluster Isomorphism With Sequential DBSCAN", origin: "Classical — Cited")[
  Under any $ell_p$ metric ($1 <= p <= infinity$), the partition of $cal(X)_"core"$ induced by MR-DBSCAN with core-gated merging is identical to the connected components of the global DBSCAN core graph $G_"core" (cal(X))$ @he2011mrdbscan @schubert2017dbscanrevisited. Moreover, every global noise point receives label $0$, and every global border point receives the label of some global core component within distance $epsilon$.
]
#proof[
  We prove that two global core points $u, v in cal(X)_"core"$ receive the same global MR-DBSCAN label iff they belong to the same connected component of $G_"core" (cal(X))$.

  *($=>$) Soundness:* Suppose local clusters $(r, c)$ and $(r', c')$ are connected by an edge in $cal(H)_"merge"$, witnessed by a shared point $z in C_(r, c) inter C_(r', c')$ with (without loss of generality) $"core"_r (z) = 1$. By Lemma 1, $"core"_r (z) <= "core"(z)$, so $z in cal(X)_"core"$. In partition $r$, every local core point in $C_(r, c)$ is locally core-connected to $z$, and since $"core"_r <= "core"$, every local core-edge in $C_(r, c)$ is a valid edge of $G_"core" (cal(X))$. In partition $r'$, $z$ belongs to $C_(r', c')$, so either $z$ is locally core in $r'$ or $z$ is a local border point directly reachable within $epsilon$ from some local core point $w in C_(r', c')$; since $z, w in cal(X)_"core"$ and $d(z, w) <= epsilon$, $(z, w) in E_"core"$. Hence every local core point of $C_(r, c)$ and $C_(r', c')$ lies in the same connected component of $G_"core" (cal(X))$. By induction along paths in $cal(H)_"merge"$, merged global clusters never unite distinct components of $G_"core" (cal(X))$.

  *($<=$) Completeness:* It suffices to show that every single edge $(u, v) in E_"core"$ maps $u$ and $v$ to the same global MR-DBSCAN cluster. Let $r = h(u)$ and $r' = h(v)$ be their unique home partitions. By Lemma 1, $"core"_r (u) = "core"(u) = 1$ and $"core"_(r') (v) = "core"(v) = 1$. Since $u in R_r$ and $d(u, v) <= epsilon$, Lemma 1 also guarantees $v in N_epsilon^(r) (u) subset.eq cal(X)_r$. Therefore, inside reducer $r$, $v$ is directly density-reachable from the local core point $u$, so $v$ is assigned to the same local cluster $(r, c)$ as $u$ (either as a local core point or as a local border point). Meanwhile, in reducer $r'$, $v$ is a local core point belonging to some local cluster $(r', c')$. Thus $v$ is shared by $(r, c)$ and $(r', c')$ with $"core"_(r') (v) = 1$, creating the edge $((r, c), (r', c')) in cal(H)_"merge"$ and uniting $u$ and $v$. By induction on path length in $G_"core" (cal(X))$, every global core component is intact.

  Finally, any global border point $b$ has $h(b) = r$ and by Lemma 1 possesses at least one global core point $u in N_epsilon^(r) (b)$; either $u$ is locally core in $r$ or $b$ appears in $h(u)$'s margin, ensuring $b$ is attached to a valid core component within $epsilon$ (exact label matching sequential DBSCAN up to border tie-breaking @ester1996dbscan). Any global noise point $x$ has no global core neighbor within $epsilon$; by $"core"_r <= "core"$, it has no local core neighbor in any partition and remains label $0$.
]

#subsec("B", "Failure of Unconstrained Border-Point Merging")

#proposition(title: "Spurious Cluster Fusion Without Core-Gated Merging", origin: "Original Counterexample")[
  Suppose the global merge step links local clusters $(r, c)$ and $(r', c')$ whenever $C_(r, c) inter C_(r', c') != emptyset$ without verifying $"core"_r (z) = 1$ or $"core"_(r') (z) = 1$. Then for every $"MinPts" >= 4$ and $epsilon > 0$, there exists a finite 1D dataset $cal(X)$ and a 2-way partition ${R_1, R_2}$ on which unconstrained MR-DBSCAN merges two disjoint DBSCAN core components separated by distance $Delta in (epsilon, 2 epsilon)$.
]
#proof[
  Fix $"MinPts" = m >= 4$ and $epsilon > 0$. Note that $m >= 4$ is necessary for a point $b$ adjacent to two distinct clusters to remain non-core (since $b$ plus one neighbor from each cluster already gives $|N_epsilon (b)| >= 3$). Construct $cal(X) subset bb(R)$ of size $2 m + 1$ comprising:
  1. Left clique $A$: $m - 1$ points at $-5 epsilon / 4$ and one frontier point $a_+ = -3 epsilon / 4$.
  2. Single bridge point $b = 0$.
  3. Right clique $B$: one frontier point $b_- = 3 epsilon / 4$ and $m - 1$ points at $5 epsilon / 4$.

  Under exact sequential DBSCAN, every $x in A$ satisfies $A subset.eq N_epsilon (x)$ (diameter $epsilon / 2 < epsilon$), so $|N_epsilon (x)| >= m$ and $A subset cal(X)_"core"$; symmetrically $B subset cal(X)_"core"$. The inter-core gap is $Delta = b_- - a_+ = 3 epsilon / 2 in (epsilon, 2 epsilon)$, so no core-core edge connects $A$ and $B$. Meanwhile, $N_epsilon (b) = {a_+, b, b_-}$ has cardinality $3 < m$, confirming that $b in.not cal(X)_"core"$ is strictly a shared border point. Thus sequential DBSCAN yields *two* disjoint core clusters $A$ and $B$ (with $b$ arbitrarily attached to one of them depending on scan order).

  Now partition $bb(R)$ at $x = 0$ into $R_1 = (-infinity, 0]$ and $R_2 = (0, infinity)$. Because $b = 0 in R_1 inter R_2^(+epsilon)$, point $b$ enters both $cal(X)_1$ and $cal(X)_2$. Reducer $1$ assigns $b$ to $A$'s local cluster $(1, 1)$ as a border point of $a_+$; Reducer $2$ assigns $b$ to $B$'s local cluster $(2, 1)$ as a border point of $b_-$. Unconstrained merging observes $b in C_(1, 1) inter C_(2, 1)$ and spuriously fuses $A$ and $B$ into a single cluster, whereas core-gated merging checks $"core"_1 (b) = "core"_2 (b) = 0$ and preserves the two separate clusters.
]

#subsec("C", "Boundary Replication and Communication Complexity")

#theorem(title: "Dimension-Dependent Boundary Replication Bounds", origin: "Original Derivation")[
  Let $cal(R) = {R_1, dots, R_P}$ partition $[0, L]^d$ into axis-aligned hyperrectangles whose side lengths satisfy $ell_(r,j) = b_(r,j) - a_(r,j) >= 2 epsilon$ for all $r, j$. Then:
  1. *(Deterministic Multiplicity Bound)* Every point $x in cal(X)$ is replicated to at most $2^d$ partitions, so the global replication factor $rho_"rep" = (sum_(r=1)^P |cal(X)_r|) / n$ obeys $1 <= rho_"rep" <= 2^d$.
  2. *(Expected Replication Under Uniform Grid)* If $P = m^d$ cubic partitions of side $s = L / m >= 2 epsilon$ cover a dataset drawn i.i.d. from a uniform distribution on $[0, L]^d$, then
  $
    bb(E)[rho_"rep"] <= (1 + (2 epsilon (m - 1)) / L)^d < (1 + (2 epsilon P^(1/d)) / L)^d.
  $
]
#proof[
  *(1)* Fix any $x = (x_1, dots, x_d) in cal(X)$ and coordinate $j in {1, dots, d}$. Suppose $x in R_r^(+epsilon)$, meaning $a_(r,j) - epsilon <= x_j <= b_(r,j) + epsilon$. Along axis $j$, any interval $[a_(r,j), b_(r,j)]$ has length at least $2 epsilon$. Consider the open ball $(x_j - epsilon, x_j + epsilon)$ of length $2 epsilon$: it can overlap at most two disjoint intervals of length $>= 2 epsilon$ (the interval containing $x_j$ and at most one adjacent interval on the nearer side of $x_j$). Taking the Cartesian product across all $d$ coordinates, any point $x$ lies in at most $2^d$ extended regions $R_r^(+epsilon)$. Summing over all $n$ points yields $sum_r |cal(X)_r| <= 2^d n$.

  *(2)* For a regular $m times dots times m$ grid of cell width $s = L / m$, each interior face along axis $j$ is shared by two adjacent cells and adds an $epsilon$-slab of measure $2 epsilon$ to the sum of extended 1D interval lengths. Across $m - 1$ internal split hyperplanes per axis, the sum of 1D extended lengths clipped to $[0, L]$ is exactly $L + 2 epsilon (m - 1)$. Because the $d$ coordinates of a uniform random point $X ~ "Unif"([0, L]^d)$ are independent, linearity of expectation over the $m^d$ cells yields
  $
    bb(E)[rho_"rep"] &= sum_(r=1)^(m^d) "Pr"(X in R_r^(+epsilon)) \
    &= ( (L + 2 epsilon (m - 1)) / L )^d \
    &= (1 + (2 epsilon (m - 1)) / L)^d.
  $
  Substituting $m = P^(1/d)$ yields the strict upper bound.
]

Theorem 3 quantifies the fundamental geometric limitation of spatial MapReduce DBSCAN: keeping communication overhead bounded requires $2 d epsilon P^(1/d) / L = cal(O)(1)$, which breaks down rapidly in high ambient dimension $d$ where the surface-to-volume ratio of hyperrectangles explodes exponentially and nearest-neighbor distances concentrate @beyer1999nearest.

#subsec("D", "PRBP Cost-Bisection Load-Balance Guarantee")

#theorem(title: "Load-Balance and Parallel Speedup Bound for PRBP", origin: "Original Proof")[
  Suppose PRBP recursively partitions a grid with total computational cost $W_"tot" = sum_(g in cal(G)) c_g$ into $P = 2^k$ regions by minimizing $|Phi(R_L) - Phi(R_R)|$ (setting $lambda = 0$) along a chosen axis at each internal node. Let $s_"max"$ be the maximum cost of any single grid hyperplane slab orthogonal to a chosen split axis across all recursion levels. Then every final partition $R_r$ satisfies
  $
    max_(1 <= r <= P) Phi(R_r) <= W_"tot" / P + s_"max" (1 - 1 / P),
  $
  and the resulting parallel speedup $S_P = W_"tot" / max_r Phi(R_r)$ satisfies
  $
    S_P >= P / (1 + (P - 1) (s_"max" / W_"tot")).
  $
]
#proof[
  We prove by induction on recursion depth $t in {0, 1, dots, k}$ that every region $R$ at depth $t$ satisfies $Phi(R) <= W_"tot" / 2^t + s_"max" (1 - 2^(-t))$. The base case $t = 0$ holds with equality since $Phi(Omega) = W_"tot"$.

  For the inductive step, consider any region $R$ at depth $t - 1$ with total cost $W = Phi(R)$. Projecting the grid cells of $R$ along the split axis yields an ordered sequence of nonnegative slab costs $w_1, dots, w_m$ satisfying $sum_(i=1)^m w_i = W$ and $w_i <= s_"max"$ for all $i$. Let $S_j = sum_(i=1)^j w_i$ with $S_0 = 0$, and let $q in {1, dots, m}$ be the smallest index such that $S_q >= W / 2$. Define the two deviations $delta_1 = W / 2 - S_(q-1) > 0$ and $delta_2 = S_q - W / 2 >= 0$. Since
  $
    delta_1 + delta_2 = S_q - S_(q-1) = w_q <= s_"max",
  $
  we have $min(delta_1, delta_2) <= w_q / 2 <= s_"max" / 2$. Cutting before slab $q$ yields child costs $(W / 2 - delta_1, W / 2 + delta_1)$, while cutting after slab $q$ yields $(W / 2 + delta_2, W / 2 - delta_2)$. Because PRBP selects the cut minimizing $|Phi(R_L) - Phi(R_R)|$, the heavier child satisfies
  $
    max(Phi(R_L), Phi(R_R)) <= W / 2 + min(delta_1, delta_2) <= (W + s_"max") / 2.
  $
  Applying the inductive hypothesis to $W$ yields
  $
    max(Phi(R_L), Phi(R_R)) &<= 1 / 2 ( W_"tot" / 2^(t-1) + s_"max" (1 - 2^(-(t-1))) + s_"max" ) \
    &= W_"tot" / 2^t + s_"max" (1 - 2^(-t)).
  $
  Setting $t = k$ with $P = 2^k$ establishes the partition bound. Dividing $W_"tot"$ by the upper bound immediately yields the speedup lower bound.
]

Theorem 4 gives an exact mathematical criterion for MapReduce scalability: PRBP achieves near-linear speedup $S_P approx P$ whenever the heaviest dense slab satisfies $s_"max" << W_"tot" / P$, whereas extreme single-hotspot skew ($s_"max" = Theta(W_"tot")$) caps parallel speedup at $W_"tot" / s_"max"$ regardless of cluster size.

#sec("III", "Complexity, Geometry and Operational Profile")

Let $n_"max" = max_r |cal(X)_r|$ denote the heaviest extended partition and $m_"ov" = sum_r |cal(X)_r^"mg"|$ denote the total boundary-margin records emitted across all reducers.

- *Grid Census & PRBP (Job 1):* Mapping points to grid cells requires $cal(O)(n d)$ total work ($cal(O)(n d / P)$ per mapper). Because $|cal(G)| << n$ on large datasets, centralized PRBP bisection takes $cal(O)(|cal(G)| d log P)$ time and negligible memory.
- *Local Indexed DBSCAN (Job 2):* On reducer $r$ with $n_r = |cal(X)_r|$ points, building an R\*-tree and answering $n_r$ range queries costs $cal(O)(n_r log n_r)$ in low dimension $d$, or $cal(O)(n_r^2 d)$ in worst-case high dimension where spatial trees degenerate to linear scans. Thus parallel wall-clock time is governed by $cal(O)(n_"max" log n_"max")$ (low $d$) or $cal(O)(n_"max"^2 d)$ (high $d$).
- *Global Core-Bridge Merge (DSU):* Processing $m_"ov"$ boundary records with path-compressed, union-by-rank Disjoint-Set Union requires $cal(O)(m_"ov" alpha(m_"ov"))$ time, where $alpha(dot)$ is the inverse Ackermann function.
- *Relabeling (Job 3) & Space:* Final relabeling is a map-only linear pass of $cal(O)(n)$ time. Peak per-reducer memory is $cal(O)(n_"max" d)$, and total Map-to-Reduce network shuffle scales as $cal(O)(n d dot rho_"rep")$.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$cal(O)(n_"max" log n_"max" + m_"ov" alpha(m_"ov"))$ low-$d$ indexed; $cal(O)(n_"max"^2 d)$ worst-case unindexed.],
    [Space / Shuffle], [$cal(O)(n_"max" d)$ RAM per reducer; $cal(O)(n d rho_"rep")$ shuffle volume with $rho_"rep" <= 2^d$.],
    [Inductive Bias], [Exact global single-scale density level sets; arbitrary non-convex manifolds with noise rejection.],
    [Hyperparameters], [Radius $epsilon > 0$, $"MinPts" >= 2$, partition count $P$, grid width $w >= 2 epsilon$, balance weight $lambda$.],
    [Exactness], [Exact core-component isomorphism with sequential DBSCAN (Theorem 2); border order-dependent.],
    [Parallel Bottleneck], [Heaviest grid slab $s_"max"$ (Theorem 4) and MapReduce job-barrier I/O latency on small $n$.],
  ),
  caption: [Operational and complexity profile of MR-DBSCAN.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Web-Scale GPS Mobility & Geosocial Hotspots:* Ride-hailing, fleet telematics, and geotagged social check-ins generate $10^8$–$10^9$ 2D/3D spatial coordinates with severe urban-rural density skew @he2011mrdbscan. Because $d = 2$ keeps boundary replication tiny ($rho_"rep" <= 4$) while PRBP subdivides dense city cores finely, MR-DBSCAN extracts exact irregular road-corridor and venue clusters without memory overflow.
2. *Large-Scale Astronomical Sky Surveys:* Catalogs such as Gaia and SDSS map billions of stars across celestial coordinates and proper-motion subspaces. Filamentary stellar streams and open clusters are strongly non-convex and immersed in sparse field-star noise; MR-DBSCAN preserves exact global core connectivity across sky tiles while discarding background noise.
3. *Distributed 3D LiDAR Point-Cloud Segmentation:* Autonomous-driving map pipelines process terabytes of 3D LiDAR scans across highway corridors @rusu2011pcl. With $d = 3$, $2^3 = 8$ bounds point multiplicity tightly, and exact $epsilon$-margin preservation guarantees that elongated guardrails, curbs, and building facades are never artificially fractured at tile boundaries.
4. *Telecom Call-Detail Records & Base-Station Anomaly Detection:* Spatial aggregation of millions of dropped-call and handover events exhibits sharp localized bursts during outages or mass gatherings. Cost-based partitioning prevents a single stadium-scale outage from stalling the cluster, while noise filtering isolates coherent failure footprints.
5. *Distributed Environmental & IoT Sensor Grids:* Satellite precipitation footprints, wildfire thermal detections, and marine AIS vessel-density logs arrive in massive spatial batches. MR-DBSCAN's exact equivalence to DBSCAN ensures scientifically reproducible perimeter definitions across distributed cloud infrastructure.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

MR-DBSCAN inherits single-scale DBSCAN's structural limitations while adding distributed-systems failure modes: (i) *high-dimensional replication collapse*, where $rho_"rep" approx (1 + 2 epsilon m / L)^d$ multiplies shuffle traffic exponentially for $d > 10$ and R\*-trees degrade to $cal(O)(n_"max"^2)$; (ii) *multi-density blindness*, where a single global $epsilon$ cannot simultaneously resolve dense urban cores and sparse rural clusters without either merging the former or shattering the latter; (iii) *unbreakable single-slab skew*, where a point cloud concentrated inside a single $2 epsilon$-wide slab forces $s_"max" approx W_"tot"$ (Theorem 4); and (iv) *MapReduce synchronization overhead*, where multi-stage disk/shuffle barriers make MR-DBSCAN slower than in-memory clustering on small or medium datasets ($n < 10^5$).

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [2D/3D Geospatial Big Data ($n > 10^6$)], [Ideal: low $d$, small $rho_"rep"$, PRBP handles urban skew.], [Use VDBSCAN / HDBSCAN if densities vary widely.],
    [3D LiDAR / Sky Surveys], [Strong: exact non-convex stitching across tiles.], [Ensure $w >= 2 epsilon$ to cap multiplicity at $8$.],
    [Moderate-Scale In-Memory Data ($n < 10^5$)], [Exact output, but shuffle and merge add overhead.], [Prefer sequential DBSCAN @ester1996dbscan or DBSCAN++ @jang2019dbscanpp.],
    [High-Dimensional Embeddings ($d > 15$)], [Poor: $2^d$ boundary blow-up and distance concentration.], [Use SNN-DPC @liu2018snndpc or spectral / projection methods.],
    [Multi-Density Mixtures], [Weak: single global $epsilon$ forces split-vs-merge trade-off.], [Use VDBSCAN (Alg. 50) or HDBSCAN @campello2013hdbscan.],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for MR-DBSCAN.],
)

#sec("V", "Structural Relations and Empirical Behavior")

#proposition(title: "Exact Reduction Hierarchy Across Distributed Density Algorithms", origin: "Original Derivation")[
  MR-DBSCAN unifies three classical density clustering regimes as strict special cases of its partition operator $cal(R)$: (i) setting $P = 1$ yields $M_1 = emptyset$ and reduces MR-DBSCAN identically to sequential DBSCAN @ester1996dbscan; (ii) replacing PRBP with a static uniform grid of cell width $w >= 2 epsilon$ recovers GriDBSCAN's partitioning topology @mahran2008gridbscan (suffering unbounded load imbalance when $max_r Phi(R_r) >> W_"tot" / P$); and (iii) relaxing Euclidean $ell_p$ neighborhoods to symmetric generalized predicates satisfying a spatial-range envelope generalizes directly to distributed GDBSCAN @sander1998gdbscan.
]
#proof[
  For (i), $P = 1$ sets $R_1 = Omega$, so $cal(X)_1^"in" = cal(X)$, $M_1 = emptyset$, $m_"ov" = 0$, and Job 2 executes a single exact DBSCAN pass with trivial identity relabeling $pi$. For (ii), constraining PRBP split hyperplanes to equal-width coordinates $t_j = a_j + k (b_j - a_j) / m$ produces the identical uniform box partition and $epsilon$-halo replication of GriDBSCAN, forfeiting Theorem 4's cost balance. For (iii), Lemma 1 and Theorem 2 require only the spatial containment property $N_"pred" (x) subset.eq B_epsilon (x)$ for a bounding $ell_p$ radius $epsilon$, under which every inner point's predicate neighborhood remains exact inside $R_r^(+epsilon)$.
]

Compared to adjacent Category II algorithms, MR-DBSCAN focuses exclusively on *horizontal scalability* of fixed-scale density reachability rather than *multi-scale adaptivity*. Whereas VDBSCAN (Alg. 50, @liu2007vdbscan), OPTICS @ankerst1999optics, and HDBSCAN @campello2013hdbscan relax the single-$epsilon$ assumption to capture hierarchical or multi-density clusters, MR-DBSCAN preserves exact single-scale DBSCAN semantics across distributed clusters at petabyte scale.

*Empirical Benchmark Analysis:* On the standardized ClusterOverview four-dataset benchmark (evaluated via Silhouette score @rousseeuw1987silhouettes and Adjusted Rand Index @hubert1985comparing), MR-DBSCAN achieves:
- *Iris:* Silhouette $= 0.618$, $"ARI" = 0.534$;
- *Customers:* Silhouette $= 0.663$;
- *Multi-Shape:* Silhouette $= 0.806$ ($"ARI" =$ N/A);
- *High-D:* Silhouette $= 0.352$ ($"ARI" =$ N/A);
- *Mean Runtime:* $9.373$ ms.

These empirical measurements mirror the theoretical profile in Sections II and III with striking precision. First, MR-DBSCAN attains the highest Silhouette scores among Algs. 49–52 on *Multi-Shape* ($0.806$), *Customers* ($0.663$), and *Iris* ($0.618$), confirming Theorem 2: core-gated boundary stitching faithfully preserves non-convex geometries (moons, rings, elongated arcs) without slicing manifolds along grid cuts, while noise rejection cleanly isolates compact cores. On *Iris*, the strong Silhouette ($0.618$) paired with moderate $"ARI"$ ($0.534$) reflects classical DBSCAN behavior: Setosa separates cleanly while overlapping Versicolor and Virginica merge across a single density bridge. Second, on *High-D*, Silhouette drops sharply to $0.352$, validating Theorem 3 and high-dimensional distance concentration @beyer1999nearest, where fixed-$epsilon$ balls lose contrast and boundary replication inflates partition overlap. Finally, MR-DBSCAN records the longest mean runtime ($9.373$ ms) in the group—more than $2 times$ slower than VDBSCAN ($4.156$ ms) and over $10 times$ slower than DPC ($0.888$ ms)—illustrating the fixed overhead of grid construction, margin duplication, local clustering, and DSU merge-relabeling on small-to-medium benchmarks where distributed parallelism has not yet amortized its coordination cost.
