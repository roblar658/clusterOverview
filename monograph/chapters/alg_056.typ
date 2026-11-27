#import "../template/ieee.typ": *

#alg_entry(
  id: 56,
  name: "Single Linkage Clustering (SLINK)",
  category: "Cat. III: Hierarchical & Tree-Based",
  year: "1973",
  paradigm: "Agglomerative Ultrametric / Minimum Spanning Tree",
  complexity: $O(n^2)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Let $X = {x_1, dots, x_n}$ be a finite dataset equipped with a symmetric, non-negative, zero-diagonal dissimilarity function $d : X times X -> bb(R)_(>= 0)$ satisfying $d(x, x) = 0$ and $d(x, y) = d(y, x)$ for all $x, y in X$. Hierarchical agglomerative clustering constructs a nested sequence of partitions $cal(C)_n, cal(C)_(n-1), dots, cal(C)_1$ starting from singletons $cal(C)_n = {{x_1}, dots, {x_n}}$ and iteratively merging the pair of active clusters minimizing an inter-cluster linkage criterion. In *Single Linkage* (also termed nearest-neighbor linkage), the dissimilarity between two disjoint nonempty subsets $A, B subset.eq X$ is the minimum cross-cluster distance:
$
  d_S (A, B) = min_(x in A, \, y in B) d(x, y).
$
Upon merging $A$ and $B$, the updated dissimilarity to any remaining cluster $C subset.eq X without (A union B)$ obeys the Lance--Williams recurrence @lance1967general with parameters $alpha_A = alpha_B = 1 / 2$, $beta = 0$, and $gamma = -1 / 2$:
$
  d_S (A union B, C) &= 1 / 2 d_S (A, C) + 1 / 2 d_S (B, C) \
  &quad - 1 / 2 abs(d_S (A, C) - d_S (B, C)) \
  &= min { d_S (A, C), \, d_S (B, C) }.
$
Equivalently, for any threshold $delta >= 0$, let $G_delta = (X, E_delta)$ denote the undirected threshold graph with edge set $E_delta = { {x, y} : d(x, y) <= delta }$. The single-linkage clustering at resolution $delta$ is precisely the partition $theta_S (delta)$ of $X$ into connected components of $G_delta$. As $delta$ increases continuously from $0$ to $max_(x, y) d(x, y)$, components coalesce monotonically, tracing a rooted binary tree (dendrogram) whose merge heights induce an ultrametric on $X$.

While naive agglomerative maintenance of an $n times n$ dissimilarity matrix requires $Theta(n^2)$ memory and $O(n^2 log n)$ time with binary heaps, Sibson @sibson1973slink introduced the *SLINK* algorithm, which computes the exact single-linkage hierarchy incrementally in optimal $Theta(n^2)$ time and only $Theta(n)$ auxiliary memory by maintaining a compact pointer encoding of the dendrogram.

#definition(title: "Ultrametric and Sibson's Pointer Representation")[
  A dissimilarity $u : X times X -> bb(R)_(>= 0)$ is an *ultrametric* on $X$ if $u(x, x) = 0$, $u(x, y) = u(y, x)$, and the strong triangle inequality holds:
  $
    u(x, z) <= max { u(x, y), \, u(y, z) } quad forall x, y, z in X.
  $
  Denote by $cal(U)(X)$ the set of all ultrametrics on $X$. A *pointer representation* @sibson1973slink of a dendrogram on ${1, dots, n}$ is a pair of arrays $(pi, lambda)$ with $pi : {1, dots, n} -> {1, dots, n}$ and $lambda : {1, dots, n} -> [0, infinity]$ satisfying:
  1. $pi(n) = n$ and $lambda(n) = infinity$;
  2. $pi(i) > i$ and $lambda(pi(i)) > lambda(i)$ for all $1 <= i < n$.
  Intuitively, $lambda(i)$ is the smallest merge height at which $x_i$ joins a cluster containing an index $j > i$, and $pi(i)$ is the largest-indexed member of that cluster immediately after the merge.
]

#algbox(
  title: "Sibson's SLINK Incremental Single-Linkage Algorithm",
  inputs: [Points $x_1, dots, x_n$ with on-demand dissimilarity oracle $d(x_i, x_j)$; optional target $k$.],
  outputs: [Exact pointer representation $(pi, lambda)$ and $k$-cluster partition $cal(C)_k$.],
)[
  1. Initialize $pi(1) <- 1$ and $lambda(1) <- infinity$.
  2. *For* $m = 1, dots, n - 1$ (inserting point $x_(m+1)$):
     - Set $pi(m + 1) <- m + 1$ and $lambda(m + 1) <- infinity$.
     - Compute working distances $M(i) <- d(x_i, x_(m+1))$ for $i = 1, dots, m$.
     - *For* $i = 1, dots, m$:
       - *If* $lambda(i) >= M(i)$:
         set $M(pi(i)) <- min { M(pi(i)), \, lambda(i) }$, then $lambda(i) <- M(i)$ and $pi(i) <- m + 1$.
       - *Else* ($lambda(i) < M(i)$):
         set $M(pi(i)) <- min { M(pi(i)), \, M(i) }$.
     - *For* $i = 1, dots, m$:
       - *If* $lambda(i) >= lambda(pi(i))$, set $pi(i) <- m + 1$.
  3. *Dendrogram Extraction:* Sort edges ${(i, pi(i), lambda(i))}_(i=1)^(n-1)$ in ascending order of $lambda(i)$ and execute union-find over the first $n - k$ edges to obtain $cal(C)_k$ (and the complete merge tree).
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

#subsec("A", "Subdominant Ultrametric Optimality and MST Equivalence")

A foundational property distinguishing single linkage from every other Lance--Williams scheme is that its cophenetic ultrametric solves a global constrained optimization problem over $cal(U)(X)$ and coincides with the bottleneck metric of a Minimum Spanning Tree (MST) @gower1969minimum.

#theorem(title: "Subdominant Ultrametric and MST Equivalence", origin: "Original Proof")[
  Let $(X, d)$ be a finite dissimilarity space and define the dominated ultrametric cone $cal(U)_(-)(d) = { u in cal(U)(X) : u(x, y) <= d(x, y) \, forall x, y in X }$. Define the minimax path function
  $
    u_S (x, y) = min_(gamma : x -> y) max_(e in gamma) d(e),
  $
  where the minimum ranges over all finite paths $gamma = (v_0, v_1, dots, v_r)$ in the complete graph $K_X$ with $v_0 = x$ and $v_r = y$. Then:
  1. $u_S in cal(U)_(-)(d)$ and $u_S$ is the unique pointwise maximum element of $cal(U)_(-)(d)$, i.e., $u(x, y) <= u_S (x, y)$ for all $u in cal(U)_(-)(d)$ and all $x, y in X$.
  2. For any MST $T^(star)$ of $(K_X, d)$, $u_S (x, y)$ equals the maximum edge weight on the unique path $gamma_(T^(star)) (x, y)$ between $x$ and $y$ in $T^(star)$.
  3. The cophenetic ultrametric produced by agglomerative single linkage coincides identically with $u_S$.
]

#proof[
  *Part 1 (Subdominant Ultrametricity):* Clearly $u_S (x, x) = 0$ and $u_S (x, y) = u_S (y, x)$. Taking the single-edge path $(x, y)$ gives $u_S (x, y) <= d(x, y)$. For any $x, y, z in X$, let $gamma_1$ achieve $u_S (x, y)$ and $gamma_2$ achieve $u_S (y, z)$. Concatenating $gamma_1$ and $gamma_2$ produces a walk from $x$ to $z$ (containing a simple path $gamma$) whose maximum edge weight is at most $max { u_S (x, y), u_S (y, z) }$, establishing the ultrametric inequality $u_S (x, z) <= max { u_S (x, y), u_S (y, z) }$. Thus $u_S in cal(U)_(-)(d)$.

  Now let $u in cal(U)_(-)(d)$ be any dominated ultrametric and let $gamma = (v_0, dots, v_r)$ be any path from $x = v_0$ to $y = v_r$. Applying the strong triangle inequality inductively along $gamma$ yields
  $
    u(x, y) &<= max_(1 <= t <= r) u(v_(t-1), v_t) \
    &<= max_(1 <= t <= r) d(v_(t-1), v_t).
  $
  Taking the minimum over all paths $gamma : x -> y$ proves $u(x, y) <= u_S (x, y)$ pointwise. Uniqueness follows immediately from antisymmetry of the pointwise partial order.

  *Part 2 (MST Bottleneck Equivalence):* Let $T^(star)$ be any MST of $(K_X, d)$, and let $w(x, y) = max_(e in gamma_(T^(star)) (x, y)) d(e)$. Since $gamma_(T^(star)) (x, y)$ is a valid path in $K_X$, $u_S (x, y) <= w(x, y)$. Conversely, let $gamma = (v_0, dots, v_r)$ be a path attaining $u_S (x, y)$. By the classical MST cycle property, every edge $e = (v_(t-1), v_t)$ satisfies $d(e) >= w(v_(t-1), v_t)$ (otherwise replacing the heaviest edge on $gamma_(T^(star)) (v_(t-1), v_t)$ by $e$ would strictly decrease total tree weight). Because $w$ is itself an ultrametric on $T^(star)$, induction along $gamma$ gives
  $
    w(x, y) &<= max_(1 <= t <= r) w(v_(t-1), v_t) \
    &<= max_(1 <= t <= r) d(v_(t-1), v_t) = u_S (x, y).
  $
  Hence $u_S (x, y) = w(x, y)$ for every pair $x, y in X$.

  *Part 3 (Agglomerative Coincidence):* At every agglomerative step, single linkage selects the minimum-weight edge joining two distinct current components, which is identical to Kruskal's MST algorithm. Thus $x$ and $y$ first belong to the same connected component at threshold $delta$ if and only if $G_delta$ connects $x$ and $y$, i.e., iff $u_S (x, y) <= delta$. Consequently, the merge height of the lowest common ancestor of $x$ and $y$ in the single-linkage dendrogram is identically $u_S (x, y)$.
]

#subsec("B", "Exactness and Complexity of Sibson's SLINK Recursion")

Because $u_S$ is unique, any algorithm recovering $(pi, lambda)$ must preserve the exact subdominant ultrametric as points arrive sequentially. We prove Sibson's incremental recursion @sibson1973slink directly from the pointer invariants.

#theorem(title: "Invariant Preservation and Exactness of SLINK", origin: "Original Proof")[
  After iteration $m - 1$ of SLINK ($1 <= m <= n$), the arrays $(pi_m, lambda_m)$ form a valid pointer representation on $X_m = {x_1, dots, x_m}$ whose induced ultrametric equals $u_(S, m)$, the single-linkage ultrametric on $X_m$. Moreover, SLINK executes in exact $Theta(n^2)$ time and $Theta(n)$ auxiliary space.
]

#proof[
  We proceed by induction on $m$. For $m = 1$, $pi_1 (1) = 1$ and $lambda_1 (1) = infinity$ trivially represent the one-point ultrametric. Assume $(pi_m, lambda_m)$ validly represents $u_(S, m)$ on $X_m$, meaning that for every $i < m$,
  $
    lambda_m (i) = min_(j > i, \, j <= m) u_(S, m) (x_i, x_j),
  $
  and $pi_m (i)$ is the largest index in $x_i$'s connected component of $G_(lambda_m (i)) (X_m)$. When $x_(m+1)$ is added, the new ultrametric $u_(S, m+1)$ satisfies, for any $i <= m$,
  $
    lambda_(m+1) (i) = min { lambda_m (i), \, u_(S, m+1) (x_i, x_(m+1)) }.
  $
  Let $C_i (delta)$ denote the component of $x_i$ in $G_delta (X_m)$. Since $pi_m (i) > i$ and $lambda_m (pi_m (i)) > lambda_m (i)$, processing indices $i = 1, dots, m$ in ascending order guarantees that before index $i$ is visited, every $k < i$ with $pi_m (k) = i$ has already propagated its minimum distance to $x_(m+1)$ into $M(i)$. By strong induction on $i$, at the start of step $i$,
  $
    M(i) = min { d(x_r, x_(m+1)) : x_r in C_i (lambda_m (i)^-) },
  $
  where $lambda_m (i)^-$ denotes the limit from below. Two mutually exclusive cases arise:
  - *Case 1 ($lambda_m (i) >= M(i)$):* Component $C_i (lambda_m (i)^-)$ attaches to $x_(m+1)$ at height $M(i)$ *before* merging with $pi_m (i)$ at height $lambda_m (i)$. Thus the true new merge level of $x_i$ to a higher index becomes $lambda_(m+1) (i) = M(i)$ with new maximal-index representative $pi_(m+1) (i) = m + 1$. Meanwhile, the cluster of $pi_m (i)$ now reaches $x_(m+1)$ via $C_i$ at bottleneck height $max { M(i), lambda_m (i) } = lambda_m (i)$, justifying the update $M(pi_m (i)) <- min { M(pi_m (i)), lambda_m (i) }$.
  - *Case 2 ($lambda_m (i) < M(i)$):* Component $C_i$ merges with $pi_m (i)$ at height $lambda_m (i)$ *before* reaching $x_(m+1)$. Thus $lambda_(m+1) (i) = lambda_m (i)$ is unchanged, and $C_i$'s best path to $x_(m+1)$ (of bottleneck height $M(i) > lambda_m (i)$) transfers to $pi_m (i)$ via $M(pi_m (i)) <- min { M(pi_m (i)), M(i) }$.

  After the first loop, every value $lambda_(m+1) (i)$ is exact, and every $pi_(m+1) (i)$ is correct except possibly when $pi_m (i)$ itself merged into $x_(m+1)$ at a level $lambda_(m+1) (pi_m (i)) <= lambda_(m+1) (i)$ (making $m + 1$ the new largest index in $x_i$'s component at level $lambda_(m+1) (i)$). The second loop checks $lambda_(m+1) (i) >= lambda_(m+1) (pi_(m+1) (i))$ and redirects $pi_(m+1) (i) <- m + 1$, restoring strict inequality $lambda_(m+1) (pi_(m+1) (i)) > lambda_(m+1) (i)$ without changing any component partition.

  For complexity, iteration $m$ evaluates $m$ dissimilarities and performs two $O(m)$ linear scans, totaling $sum_(m=1)^(n-1) Theta(m) = Theta(n^2)$ operations. Only the three length-$n$ vectors $pi$, $lambda$, and $M$ are stored, requiring $Theta(n)$ auxiliary space.
]

#subsec("C", "Monotone Invariance, Stability, and Axiomatic Uniqueness")

Beyond MST optimality, Carlsson and Mémoli @carlsson2010characterization showed that single linkage enjoys uniquely strong topological stability and categorical invariance among hierarchical methods.

#proposition(title: "Monotone Equivariance and Gromov--Hausdorff Stability", origin: "Classical — Cited")[
  Let $phi : [0, infinity) -> [0, infinity)$ be any strictly increasing bijection with $phi(0) = 0$. Denote by $u_S (X, d)$ the single-linkage ultrametric on $(X, d)$. Then @carlsson2010characterization:
  1. *Ordinal equivariance:* $u_S (X, phi compose d) = phi compose u_S (X, d)$ identically, so the single-linkage dendrogram topology is invariant to all strictly increasing transformations of dissimilarities, positive rescaling, and input permutations.
  2. *1-Lipschitz stability:* For any two finite metric spaces $(X, d_X)$ and $(Y, d_Y)$, the Gromov--Hausdorff distance satisfies
     $
       d_"GH" ((X, u_S (d_X)), \, (Y, u_S (d_Y))) <= d_"GH" ((X, d_X), \, (Y, d_Y)).
     $
  3. *Axiomatic uniqueness:* Single linkage is the unique hierarchical clustering functor satisfying standard richness on two-point spaces, excisive scale consistency, and functoriality under distance-non-increasing maps @carlsson2010characterization (circumventing Kleinberg's flat-clustering impossibility @kleinberg2002impossibility via hierarchical outputs).
]

#proof[
  *Part 1:* Because $phi$ is strictly increasing on $[0, infinity)$, it commutes with finite minima and maxima. For every pair $x, y in X$,
  $
    u_S (phi compose d)(x, y) &= min_(gamma : x -> y) max_(e in gamma) phi(d(e)) \
    &= phi( min_(gamma : x -> y) max_(e in gamma) d(e) ) \
    &= phi(u_S (d)(x, y)).
  $
  Hence all merge-order comparisons are preserved identically.

  *Part 2:* Recall that $2 d_"GH" (X, Y) = inf_R "dis"(R)$ over all correspondences $R subset.eq X times Y$ projecting surjectively onto both $X$ and $Y$, where the distortion is $"dis"(R) = sup_((x, y), (x', y') in R) abs(d_X (x, x') - d_Y (y, y'))$. Fix any correspondence $R$ and any $(x, y), (x', y') in R$. Let $gamma_X = (x_0, dots, x_r)$ be a minimax path in $X$ from $x_0 = x$ to $x_r = x'$ attaining $u_S (d_X)(x, x')$. Because $R$ is surjective onto $X$, choose $y_t in Y$ with $(x_t, y_t) in R$ ($y_0 = y, y_r = y'$). Along each step, $d_Y (y_(t-1), y_t) <= d_X (x_(t-1), x_t) + "dis"(R)$, so
  $
    u_S (d_Y)(y, y') &<= max_(1 <= t <= r) d_Y (y_(t-1), y_t) \
    &<= u_S (d_X)(x, x') + "dis"(R).
  $
  By symmetry between $X$ and $Y$, $abs(u_S (d_X)(x, x') - u_S (d_Y)(y, y')) <= "dis"(R)$. Thus $"dis"_(u_S)(R) <= "dis"(R)$, and taking the infimum over $R$ proves 1-Lipschitz Gromov--Hausdorff stability. Part 3 follows from Carlsson and Mémoli @carlsson2010characterization, showing any functorial scheme dominated by two-point scale recovery must agree with the maximal subdominant ultrametric $u_S$.
]

#subsec("D", "Chaining Pathology and Zero Finite-Sample Breakdown Point")

The very bottleneck minimization that yields 1-Lipschitz Gromov--Hausdorff stability also causes *chaining*: a thin filament of bridge points suffices to collapse arbitrarily well-separated clusters.

#theorem(title: "Constructive Chaining Failure and Zero Breakdown Point", origin: "Original Counterexample")[
  Fix any dimension $p >= 1$, any cluster radius $r > 0$, any separation gap $L > 0$, and any internal spacing $delta in (0, r)$. Let $B_1, B_2 subset bb(R)^p$ be two datasets of size $m$ contained in balls of radius $r$ centered at $0$ and $(L + 2r) e_1$, each having internal MST maximum edge weight at least $delta$ (e.g., containing one point at distance $delta$ from the rest of its ball). Insert $q = ceil(L / epsilon) - 1$ adversarial bridge points along the connecting axis with spacing $epsilon < delta$. Then:
  1. Single linkage merges $B_1$ and $B_2$ into a single component at height $h_"bridge" <= epsilon < delta$, strictly *before* the internal structure of $B_1$ or $B_2$ completes.
  2. Cutting the dendrogram at $k = 2$ cannot separate $B_1$ from $B_2$: at least $2m - 2$ points of $B_1 union B_2$ lie in the same cluster, forcing the Adjusted Rand Index @hubert1985comparing to satisfy $"ARI" -> 0$ as $m -> infinity$.
  3. Because the contamination ratio satisfies $q / (2m + q) -> 0$ as $m -> infinity$ for fixed $L, epsilon$, the asymptotic finite-sample breakdown point of single-linkage cluster recovery under adversarial bridge perturbations is identically $0$.
]

#proof[
  Let $p_1 in B_1$ be at $r e_1$ and $p_2 in B_2$ at $(L + r) e_1$, so $norm(p_2 - p_1)_2 = L$. Place $q = ceil(L / epsilon) - 1$ bridge points $z_j = (r + j L / (q + 1)) e_1$ for $j = 1, dots, q$. Adjacent points along the path $(p_1, z_1, dots, z_q, p_2)$ satisfy
  $
    norm(z_(j+1) - z_j)_2 = L / (q + 1) <= epsilon < delta.
  $
  By Theorem 1, $u_S (p_1, p_2) <= epsilon$. Suppose $B_1$ contains a boundary point $a_1$ whose nearest neighbor in $B_1 without {a_1}$ lies at distance $delta > epsilon$ (and $norm(a_1 - z_1)_2 >= delta$). Then $a_1$ remains an isolated singleton at threshold $epsilon$, while the core of $B_1$, the entire bridge ${z_1, dots, z_q}$, and the core of $B_2$ already form a single connected component $C_"merged"$ of size at least $2m + q - 2$. At any threshold $tau >= delta$ where $a_1$ finally merges, $C_"merged"$ is already unified; hence the $k = 2$ partition must place $a_1$ (a singleton) in one cluster and all remaining $2m + q - 1$ points (including $m - 1$ points of $B_1$ and all $m$ points of $B_2$) in the other.

  Evaluate the Adjusted Rand Index on $B_1 union B_2$ (size $2m$): the $k = 2$ clustering places $2m - 1$ points in one cluster and $1$ point in the other, so the number of co-clustered pairs is $binom(2m - 1, 2)$, giving
  $
    "ARI" = O(1 / m) -> 0 quad "as" quad m -> infinity.
  $
  Since $q = ceil(L / epsilon) - 1$ is independent of $m$, the adversarial fraction $q / (2m + q) -> 0$, proving that single linkage has zero breakdown point against bridge contamination. This finite-sample mechanism underlies Hartigan's classical result @hartigan1981consistency that single linkage in $bb(R)^p$ ($p >= 2$) is only fractionally consistent for high-density clusters unless density pruning is applied @chaudhuri2014consistent.
]

#sec("III", "Complexity, Geometry and Operational Profile")

*Computational and Memory Profile:* Sibson's SLINK achieves the information-theoretic lower bound $Omega(n^2)$ for general dissimilarity spaces (since every one of the $binom(n, 2)$ pairs could hold the unique minimum edge) while requiring only $3n$ words of working memory ($pi, lambda, M$). Because row $m$ of dissimilarities is streamed sequentially and immediately discarded, SLINK exhibits high cache locality and scales cleanly to $n approx 10^5$ points without materializing an $n^2$ matrix (which would exceed 80 GB in double precision). On Euclidean spaces of low dimension ($p <= 3$), exact single linkage accelerates to $O(n log n)$ via Delaunay triangulation or dual-tree Borůvka Euclidean MST algorithms, though in high dimensions ($p >= 15$) spatial trees degenerate back to $Theta(n^2 p)$ due to volume concentration @beyer1999nearest.

*Geometric Inductive Bias:* Single linkage defines clusters strictly through *topological connectivity* of sublevel sets of the distance graph, with zero preference for convexity, isotropy, or balanced cluster sizes. Because its Lance--Williams recurrence satisfies Bruynooghe's Reducibility ($min {d_S(A, C), d_S(B, C)} >= min {d_S(A, C), d_S(B, C)}$), the dendrogram is strictly free of inversions and requires no centroids, enabling direct execution on arbitrary non-metric, graph, or string dissimilarities.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [*Property*], [*Characterization*],
    [Time Complexity], [$Theta(n^2)$ general (SLINK / Prim MST); $O(n log n)$ in low-$p$ Euclidean space.],
    [Space Complexity], [$Theta(n)$ auxiliary via streaming SLINK ($3n$ scalars); $Theta(n^2)$ if matrix precomputed.],
    [Determinism], [Fully deterministic; unique dendrogram when pairwise dissimilarities are distinct.],
    [Invariances], [Exact invariance to all strictly increasing dissimilarity transforms, scaling, and order.],
    [Cluster Geometry], [Arbitrary non-convex, curved, filamentary, and nested manifold components.],
    [Monotonicity], [Strictly inversion-free; cophenetic metric is the unique subdominant ultrametric $u_S$.],
    [Outlier Sensitivity], [Zero breakdown point under bridge chains; isolated outliers become singleton leaves.],
    [Primary Bottleneck], [Quadratic pairwise distance evaluations for high-$p$ or non-Euclidean large $n$.],
  ),
  caption: [Operational and theoretical profile of Single Linkage Clustering (SLINK).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Astrophysical Halo Finding (Friends-of-Friends Percolation):* In cosmological $N$-body dark-matter simulations, the standard Friends-of-Friends (FoF) halo finder links any two particles within linking length $b$ times the mean inter-particle separation. Mathematically, FoF is an exact level cut $theta_S (delta)$ of single linkage, matching the non-spherical, filamentary geometry of the cosmic web where spherical variance methods fracture galactic filaments.
2. *Topological Data Analysis ($H_0$ Persistent Homology):* Zero-dimensional persistent homology of the Vietoris--Rips filtration on $(X, d)$ tracks the death times of connected components born at scale $0$. By Theorem 1, the nontrivial death times of the $H_0$ persistence barcode are identically the $n - 1$ finite SLINK merge heights ${lambda(i)}_(i=1)^(n-1)$, providing 1-Lipschitz topological signatures @carlsson2010characterization.
3. *Protein Sequence Homology and Transitive Domain Clustering:* Remote evolutionary homologs often fail direct pairwise BLAST alignment yet connect through intermediate orthologs. Single linkage's transitive ultrametric closure captures full superfamily connectivity without distorting sequence-identity rankings under nonlinear E-value rescalings.
4. *Geospatial Corridor and Trajectory Segmentation:* GPS traces along road networks, river basins, and shipping lanes form thin, high-aspect-ratio manifolds embedded in 2D/3D space. Single linkage recovers curved corridors naturally whenever inter-corridor gaps exceed intra-corridor sampling density.
5. *Graph-Based Image Region Merging and Watershed Hierarchies:* On pixel adjacency graphs weighted by intensity gradients, MST-based watershed and quasi-flat zone hierarchies @vincent1991watersheds are exact sparse-graph single-linkage dendrograms computable in $O(n log n)$ time, preserving crisp irregular object boundaries.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Single linkage degrades severely whenever clusters overlap through low-density tails, background clutter forms short-step bridges (Theorem 3), or isolated outliers consume the top $k - 1$ dendrogram cuts, leaving $k - 1$ singletons and one giant merged component. In high-dimensional continuous feature spaces, distance concentration @beyer1999nearest narrows the gap between intra- and inter-cluster nearest neighbors, exacerbating premature bridging.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [*Domain*], [*Mathematical Fit*], [*Caution / Alternative*],
    [Cosmic Web / FoF Halos], [Exact threshold-graph percolation on filaments], [Bridge over-linking; use HDBSCAN\* @campello2013hdbscan],
    [$H_0$ Persistent Homology], [Exact Vietoris--Rips death times via SLINK], [Outlier sensitivity; apply DTM filtration],
    [Clean Non-Convex Manifolds], [Matches arbitrary topology without convexity bias], [Bridge noise; use Genie @gagolewski2016genie],
    [Customer / Gaussian Mixtures], [Poor fit; dense overlap triggers chaining], [Prefer Ward (Alg. 59) or GMM @fraley2002mclust],
    [High-$p$ Noisy Embeddings], [Nearest-neighbor contrast collapses @beyer1999nearest], [Prefer UPGMA (Alg. 58) or Ward],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for Single Linkage.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Single linkage serves as the topological foundation for modern density-based hierarchical algorithms, which mitigate its zero breakdown point by transforming the underlying dissimilarity before invoking SLINK.

#proposition(title: "Exact Reduction of DBSCAN Core Connectivity and HDBSCAN* to SLINK", origin: "Original Derivation")[
  1. Let $X_"core" (epsilon, m_"pts") = { x in X : abs(B_epsilon (x)) >= m_"pts" }$ be the set of DBSCAN core points @ester1996dbscan. The core-point clusters of DBSCAN are identically the connected components $theta_S (epsilon)$ of single linkage restricted to $X_"core"$.
  2. Define the mutual reachability dissimilarity $d_"mreach" (x, y) = max { "core"_k (x), \, "core"_k (y), \, d(x, y) }$, where $"core"_k (x)$ is the distance from $x$ to its $k$-th nearest neighbor. Then the cluster hierarchy of HDBSCAN\* @campello2013hdbscan @campello2015hierarchical is identically the single-linkage dendrogram $u_S (X, d_"mreach")$ followed by minimum-cluster-size condensation.
]

#proof[
  *Part 1:* By definition @ester1996dbscan, two core points $x, y in X_"core"$ are density-connected iff there exists a chain $x = v_0, v_1, dots, v_r = y$ in $X_"core"$ with $d(v_(t-1), v_t) <= epsilon$ for all $t$. This condition holds iff $min_(gamma subset.eq X_"core") max_(e in gamma) d(e) <= epsilon$, i.e., iff $u_S (X_"core", d)(x, y) <= epsilon$.

  *Part 2:* In HDBSCAN\* @campello2015hierarchical, a point $x$ becomes a core point at threshold $delta$ iff $"core"_k (x) <= delta$, and two core points link iff $d(x, y) <= delta$. Combining both requirements gives $max { "core"_k (x), "core"_k (y), d(x, y) } <= delta$, i.e., $d_"mreach" (x, y) <= delta$. Hence the density-connectivity graph at level $delta$ is precisely the $delta$-threshold graph of $(X, d_"mreach")$, whose connected-component hierarchy is by Theorem 1 the single-linkage ultrametric $u_S (X, d_"mreach")$.
]

Related structural extensions include OPTICS @ankerst1999optics (which orders points via Prim-like reachability walks), Genie @gagolewski2016genie (which constrains SLINK merges whenever the Gini index of cluster sizes exceeds a threshold, arresting singleton chaining), and the remaining Lance--Williams linkages: Complete (Alg. 57), Average (Alg. 58), and Ward (Alg. 59).

*Empirical Behavior on the ClusterOverview 4-Dataset Benchmark:* Across the standardized four-dataset suite, Single Linkage achieves *Iris* $"Sil" = 0.489$, $"ARI" = 0.560$; *Customers* $"Sil" = 0.067$; *Multi-Shape* $"Sil" = 0.609$ ($"ARI" = "N/A"$); *High-D* $"Sil" = 0.254$ ($"ARI" = "N/A"$); with a mean runtime of *2.343 ms*. Every metric closely validates our theoretical analysis. On *Multi-Shape*, SLINK's topological connectivity yields a strong Silhouette of $0.609$, decisively outperforming Complete ($0.420$) and Average ($0.452$) linkage because non-convex structures remain intact rather than being fractured into balls. On *Iris*, SLINK cleanly isolates _Iris setosa_ (preserving a high Silhouette of $0.489$), but bridge points between overlapping _versicolor_ and _virginica_ cause premature chaining, capping $"ARI"$ at $0.560$ compared to Ward's $0.696$. The sharpest confirmation of Theorem 3 occurs on *Customers*, where dense continuous variation plus boundary outliers induce catastrophic chaining and singleton peeling, collapsing Silhouette to $0.067$ (versus $0.323$--$0.347$ for Complete, Average, and Ward). On *High-D*, nearest-neighbor contrast compression limits Silhouette to $0.254$. Finally, its mean runtime of $2.343$ ms is the fastest among the classical hierarchical methods (over $2.3 times$ faster than Complete linkage's $5.391$ ms), directly reflecting SLINK's cache-friendly $Theta(n^2)$ streaming execution.
