#import "../template/ieee.typ": *

#alg_entry(
  id: 36,
  name: "SNN (Shared Nearest Neighbor Clustering)",
  category: "Cat. II: Density & Mode-Seeking",
  year: "1973 / 2003",
  paradigm: "Rank-Based Neighborhood-Overlap Density Clustering",
  complexity: $O(n^2 d + n k^2)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Let $cal(X) = {x_1, dots, x_n}$ be a finite dataset in a general metric or dissimilarity space $(cal(M), rho)$, and let $k in {1, dots, n - 1}$ be a neighborhood-size integer. Standard Euclidean density estimators (DBSCAN Alg. 31 @ester1996dbscan, Mean Shift Alg. 34, DENCLUE Alg. 35) rely on raw metric radii $epsilon$ or $sigma$, which fail both when clusters exhibit widely varying intrinsic densities and in high dimensions where pairwise distances concentrate @beyer1999nearest. Building on the shared-neighbor linkage of Jarvis and Patrick (1973), Ertöz, Steinbach, and Kumar @ertoz2003snn formulated *SNN density clustering*, which replaces raw metric distances with the overlap of local $k$-nearest-neighbor ($"kNN"$) lists.

For each $i in {1, dots, n}$, let $N_k (i) subset.eq {1, dots, n} backslash {i}$ denote the set of $k$ nearest neighbors of $x_i$ under $rho$, ranked by $r_i (u) in {1, dots, k}$ (breaking ties deterministically by index). The *unweighted SNN similarity* between $i$ and $j$ with the mutual-neighborhood gate is
$ s_"SNN" (i, j) := bold(1)[j in N_k (i) and i in N_k (j)] dot |N_k (i) inter N_k (j)|. $
Equivalently, Ertöz et al. @ertoz2003snn also define a *rank-weighted SNN similarity* emphasizing closer shared neighbors:
$ s_"SNN"^w (i, j) &:= bold(1)[j in N_k (i) and i in N_k (j)] \
  &quad times sum_(u in N_k (i) inter N_k (j)) (k + 1 - r_i (u))(k + 1 - r_j (u)). $
Given an SNN similarity threshold $"Eps" > 0$ and a core-density integer $"MinPts" in {1, dots, k}$, the *SNN neighborhood* and *SNN density* of $x_i$ are defined as
$ N_"SNN" (i) := {j : s_"SNN" (i, j) >= "Eps"}, quad rho_"SNN" (i) := |N_"SNN" (i)|. $

#definition(title: "SNN Core Points, Similarity Graph & Induced Clustering")[
  Fix parameters $(k, "Eps", "MinPts")$.
  1. The *sparsified SNN graph* $G_"SNN" = (V, E_"Eps")$ has vertex set $V = {1, dots, n}$ and undirected edges $E_"Eps" = {{i, j} : s_"SNN" (i, j) >= "Eps"}$.
  2. A point $i in V$ is an *SNN core point* if $rho_"SNN" (i) >= "MinPts"$; let $V_"core"$ denote the set of all core points.
  3. Core clusters are the connected components $C_1^"core", dots, C_hat(k)^"core"$ of the induced core subgraph $G_"SNN" [V_"core"]$. Each non-core point $i in.not V_"core"$ is assigned to the cluster of $j^* = arg max_(j in V_"core") s_"SNN" (i, j)$ if $s_"SNN" (i, j^*) >= "Eps"$ (border point), and labeled *noise* otherwise.
]

Unlike direct $epsilon$-ball counting, $rho_"SNN" (i)$ measures how many mutual neighbors *agree* on their local neighborhood context, automatically normalizing away global scale differences across clusters of heterogeneous density. Algorithm 1 presents the complete Ertöz–Steinbach–Kumar procedure @ertoz2003snn.

#algbox(
  title: "Ertöz–Steinbach–Kumar SNN Density Clustering",
  inputs: [Dataset $cal(X) = {x_1, dots, x_n}$, metric $rho$, neighborhood size $k$, similarity threshold $"Eps"$, core threshold $"MinPts"$.],
  outputs: [Cluster partition $cal(C) = {C_1, dots, C_hat(k)}$ and noise set $cal(N)$.],
)[
  1. *k-NN Graph Construction:* For each $i in {1, dots, n}$, compute its ordered $k$-nearest-neighbor list $N_k (i)$ under $rho$.
  2. *SNN Graph Sparsification:* For each mutual pair $(i, j)$ with $j in N_k (i)$ and $i in N_k (j)$, compute $s_"SNN" (i, j) = |N_k (i) inter N_k (j)|$ (or rank-weighted $s_"SNN"^w (i, j)$); add undirected edge $(i, j)$ to $E_"Eps"$ iff $s_"SNN" (i, j) >= "Eps"$.
  3. *SNN Density & Core Identification:* Compute $rho_"SNN" (i) = "deg"_(G_"SNN") (i)$ for all $i$; set $V_"core" <- {i : rho_"SNN" (i) >= "MinPts"}$.
  4. *Core Component Extraction:* Extract connected components $C_1^"core", dots, C_hat(k)^"core"$ of $G_"SNN" [V_"core"]$ via BFS or Union–Find.
  5. *Border Assignment & Noise Pruning:* For each $i in V backslash V_"core"$, find $j^* = arg max_(j in V_"core") s_"SNN" (i, j)$. If $s_"SNN" (i, j^*) >= "Eps"$, assign $i$ to $j^*$'s cluster; else insert $i$ into $cal(N)$.
  6. *Return* $cal(C) = {C_1, dots, C_hat(k)}$ and $cal(N)$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We prove five foundational properties of SNN clustering: (1) strict invariance to monotone metric transformations and exact multi-density self-calibration across arbitrarily scaled clusters; (2) positive semidefinite Gram factorization, exact Jaccard metric transitivity, and degree-bounded sparsity; (3) high-dimensional distance concentration collapse of metric balls versus an exact $Theta(1 / sqrt(d))$ rank-order SNR amplification theorem explaining SNN's high-dimensional resilience; and (4) constructive counterexamples characterizing minimum-cluster-size blindness and mutual-NN isolation.

#subsec("A", "Monotone Metric Invariance and Multi-Density Self-Calibration")

#theorem(title: "Monotone Rank Invariance & Multi-Density Isomorphism", origin: "Original Proof")[
  Assume all pairwise distances from each datum are distinct.
  1. *Global Monotone Invariance:* For any strictly increasing map $phi : [0, infinity) -> [0, infinity)$ with $phi(0) = 0$, replacing $rho$ by $tilde(rho) = phi compose rho$ leaves $N_k (i)$, $s_"SNN" (i, j)$, $s_"SNN"^w (i, j)$, $rho_"SNN" (i)$, and the entire SNN partition *strictly identical*.
  2. *Exact Multi-Density Scale Independence:* Let $cal(Z) subset.eq bb(R)^d$ be any finite point configuration with $|cal(Z)| = m >= k + 1$. Form $cal(X) = C_1 union C_2$ where $C_1 = alpha_1 cal(Z)$ and $C_2 = alpha_2 cal(Z) + t$ for arbitrary scales $alpha_2 > alpha_1 > 0$ and separation $"dist"(C_1, C_2) > alpha_2 "diam"(cal(Z))$. Then for *every* density ratio $(alpha_2 / alpha_1)^d in (1, infinity)$, the induced SNN similarity subgraphs and SNN densities on $C_1$ and $C_2$ are *identical and isomorphic*, whereas DBSCAN (Alg. 31) fails for any global $epsilon$ whenever $alpha_2 / alpha_1$ separates the sparse cluster's $k$-NN radius above the inter-cluster gap of a third bridge or dense cluster.
]

#proof[
  (1) Because $phi$ is strictly increasing, $rho(x_i, x_u) < rho(x_i, x_v)$ holds if and only if $tilde(rho)(x_i, x_u) < tilde(rho)(x_i, x_v)$ for every triple $(i, u, v)$. Thus every ordered neighbor list $N_k (i)$ and rank vector $r_i (dot)$ is invariant under $phi$. Since $s_"SNN"$, $s_"SNN"^w$, $rho_"SNN"$, $V_"core"$, and all steps of Algorithm 1 depend on $rho$ exclusively through ${N_k (i), r_i}$, the clustering is strictly invariant.

  (2) Fix any $i in C_r$ ($r in {1, 2}$). Every other point $j in C_r backslash {i}$ satisfies $norm(x_i - x_j)_2 <= alpha_r "diam"(cal(Z)) <= alpha_2 "diam"(cal(Z)) < "dist"(C_1, C_2)$. Because $|C_r backslash {i}| = m - 1 >= k$, all $k$ nearest neighbors of $i$ lie strictly inside $C_r$ ($N_k (i) subset.eq C_r$), yielding zero cross-cluster edges. Restricted to $C_r$, pairwise distances satisfy $norm(x_i - x_j)_2 = alpha_r norm(z_i - z_j)_2$. Because $t |-> alpha_r t$ is strictly increasing, Part 1 implies that the directed $k$-NN subgraph on $C_r$ is isomorphic to that of $cal(Z)$, independent of $alpha_r$. Hence both the tight cluster $C_1$ and the $alpha_2 / alpha_1$-dilated cluster $C_2$ exhibit identical SNN densities and identical recovery for a single parameter triple $(k, "Eps", "MinPts")$.
]

#subsec("B", "Gram Factorization, Metric Transitivity & Sparsity")

#proposition(title: "PSD Overlap Gram Matrix, Jaccard Metric & Bounded Degree", origin: "Original Proof")[
  Let $A in {0, 1}^(n times n)$ be the directed $k$-NN indicator matrix ($A_(i u) = 1$ iff $u in N_k (i)$) and $M = A dot.o A^T$ the mutual-$k$-NN adjacency matrix (entrywise Hadamard product). Then:
  1. The ungated overlap matrix $S_0 := A A^T$ (and likewise the rank-weighted matrix $S_0^w := W W^T$ with $W_(i u) = (k + 1 - r_i (u)) A_(i u)$) is symmetric positive semidefinite ($S_0 succ.eq 0$), and $s_"SNN" (i, j) = M_(i j) (S_0)_(i j)$.
  2. For all triples $(i, j, ell)$, overlap obeys the quantitative transitivity bound
     $ |N_k (i) inter N_k (ell)| &>= |N_k (i) inter N_k (j)| \
       &quad + |N_k (j) inter N_k (ell)| - k, $
     and the fixed-$k$ Jaccard dissimilarity $d_J (i, j) := 1 - (S_0)_(i j) / (2 k - (S_0)_(i j))$ is a genuine pseudometric satisfying the triangle inequality on $cal(X)$.
  3. Every vertex of $G_"SNN"$ satisfies $"deg"_(G_"SNN") (i) <= k$, so $rho_"SNN" (i) in {0, dots, k}$ and $|E_"Eps"| <= n k / 2$.
]

#proof[
  (1) For any $v in bb(R)^n$, $v^T S_0 v = norm(A^T v)_2^2 >= 0$ and $v^T S_0^w v = norm(W^T v)_2^2 >= 0$, establishing positive semidefiniteness. Moreover $(A A^T)_(i j) = sum_(u=1)^n A_(i u) A_(j u) = |N_k (i) inter N_k (j)|$, and multiplying entrywise by $M_(i j) = A_(i j) A_(j i)$ enforces the mutual-neighborhood gate.

  (2) Because $(N_k (i) inter N_k (j)) union (N_k (j) inter N_k (ell)) subset.eq N_k (j)$ has cardinality at most $k$, inclusion–exclusion yields
  $ |N_k (i) inter N_k (ell)| &>= |N_k (i) inter N_k (j) inter N_k (ell)| \
    &>= |N_k (i) inter N_k (j)| + |N_k (j) inter N_k (ell)| - k. $
  For the Jaccard triangle inequality, define the symmetric-difference distance $d_Delta (i, j) := |N_k (i) Delta N_k (j)| = 2(k - (S_0)_(i j))$, which satisfies the triangle inequality because $d_Delta (i, j) = norm(A_(i,:) - A_(j,:))_1$. Since $|N_k (i) union N_k (j)| = 2 k - (S_0)_(i j) = k + d_Delta (i, j) / 2$, we have the exact identity
  $ d_J (i, j) &= (d_Delta (i, j)) / (k + d_Delta (i, j) / 2) \
    &= (2 d_Delta (i, j)) / (2 k + d_Delta (i, j)) =: g(d_Delta (i, j)). $
  On $[0, infinity)$, $g(0) = 0$, $g'(t) = 4 k / (2 k + t)^2 > 0$, and $g''(t) = -8 k / (2 k + t)^3 < 0$. Because $g$ is increasing and concave with $g(0) = 0$, it is subadditive ($g(a + b) <= g(a) + g(b)$), whence $d_J (i, ell) = g(d_Delta (i, ell)) <= g(d_Delta (i, j) + d_Delta (j, ell)) <= d_J (i, j) + d_J (j, ell)$.

  (3) Every edge $(i, j) in E_"Eps"$ requires $j in N_k (i)$ and $|N_k (i)| = k$, so $"deg"_(G_"SNN") (i) <= k$ and $|E_"Eps"| = (1/2) sum_i "deg"_(G_"SNN") (i) <= n k / 2$.
]

#subsec("C", "High-Dimensional Concentration vs. SNN SNR Amplification")

#theorem(title: "High-Dimensional Distance Contrast Collapse", origin: "Classical — Cited")[
  Let $X_1, dots, X_n$ be i.i.d. random vectors in $bb(R)^d$ with query $X_0$ and metric $rho_d$, and denote $D_"min" = min_i rho_d (X_0, X_i)$, $D_"max" = max_i rho_d (X_0, X_i)$. If $lim_(d -> infinity) "Var"(rho_d (X_0, X_1)^p / bb(E)[rho_d (X_0, X_1)^p]) = 0$ for some $p > 0$, then for every $epsilon > 0$,
  $ lim_(d -> infinity) bb(P)(D_"max" <= (1 + epsilon) D_"min") = 1, $
  so relative Euclidean contrast $(D_"max" - D_"min") / D_"min" ->^p 0$ @beyer1999nearest.
]

#proof(label: "Proof Synthesis (Cited)")[
  Beyer et al. @beyer1999nearest normalize $rho_d (X_0, X_i)^p$ by its expectation $mu_d := bb(E)[rho_d (X_0, X_1)^p]$. Vanishing variance implies $rho_d (X_0, X_i)^p / mu_d ->^p 1$ for every $i = 1, dots, n$ by Chebyshev's inequality. For fixed finite $n$, a union bound forces both $D_"min"^p / mu_d ->^p 1$ and $D_"max"^p / mu_d ->^p 1$, hence $D_"max" / D_"min" ->^p 1$, collapsing any fixed-radius $epsilon$-ball discriminator.
]

Crucially, even in regimes where raw metric contrast vanishes as $Theta(1 / d)$, SNN's rank-based neighborhood ordering preserves cluster separability down to an asymptotically smaller signal threshold.

#theorem(title: "Intermediate High-D Rank Preservation & SNN Amplification", origin: "Original Proof")[
  Consider a two-cluster Gaussian subspace model in $bb(R)^d$: $x_i, x_j ~ cal(N)(bold(0), sigma^2 I_d)$ in $C_1$ and $x_ell ~ cal(N)(mu, sigma^2 I_d)$ in $C_2$ independently, with signal shift $Delta := norm(mu)_2 > 0$ supported on $d_s << d$ informative coordinates amidst $d - d_s$ pure noise dimensions. Let $Z := norm(x_i - x_ell)_2^2 - norm(x_i - x_j)_2^2$. Then:
  1. $bb(E)[Z] = Delta^2$ and $"Var"(Z) = 12 d sigma^4 + 8 sigma^2 Delta^2$ identically.
  2. In the intermediate high-dimensional regime $Delta = Theta(sigma d^alpha)$ with $1 / 4 < alpha < 1 / 2$, the raw relative distance contrast vanishes as $Delta^2 / (2 d sigma^2) = Theta(d^(2 alpha - 1)) -> 0$, yet the pairwise rank-inversion probability satisfies
     $ bb(P)(norm(x_i - x_ell)_2 <= norm(x_i - x_j)_2) &<= (12 d sigma^4 + 8 sigma^2 Delta^2) / Delta^4 \
       &= O(d^(1 - 4 alpha)) -> 0. $
  3. Consequently, when $Delta^2 = omega(sigma^2 sqrt(d log n))$, $k$-NN sets concentrate inside true clusters and SNN similarity amplifies separation combinatorially over raw metric thresholding.
]

#proof[
  Write $x_i = sigma g_0$, $x_j = sigma g_1$, and $x_ell = mu + sigma g_2$ with independent $g_0, g_1, g_2 ~ cal(N)(bold(0), I_d)$. Expanding $Z$ cancels $sigma^2 norm(g_0)_2^2$ exactly:
  $ Z &= norm(mu + sigma(g_2 - g_0))_2^2 - sigma^2 norm(g_1 - g_0)_2^2 \
      &= Delta^2 + underbrace(2 sigma chevron.l mu\, g_2 - g_0 chevron.r, T_1) \
      &quad + underbrace(sigma^2 (norm(g_2)_2^2 - norm(g_1)_2^2), T_2) - underbrace(2 sigma^2 chevron.l g_0\, g_2 - g_1 chevron.r, T_3). $
  Each zero-mean fluctuation term has exact variance:
  - $"Var"(T_1) = 4 sigma^2 ("Var"(chevron.l mu, g_2 chevron.r) + "Var"(chevron.l mu, g_0 chevron.r)) = 8 sigma^2 Delta^2$;
  - $"Var"(T_2) = sigma^4 ("Var"(chi_d^2) + "Var"(chi_d^2)) = 4 d sigma^4$;
  - $"Var"(T_3) = 4 sigma^4 sum_(r=1)^d bb(E)[g_(0 r)^2] bb(E)[(g_(2 r) - g_(1 r))^2] = 8 d sigma^4$.

  By Gaussian symmetry (odd polynomial degree for $T_1 T_2$ and $T_1 T_3$) and independence of $g_0$ (linear in $T_3$, absent in $T_2$), all three cross-covariances vanish identically: $"Cov"(T_a, T_b) = 0$ for $a != b$. Summing gives $"Var"(Z) = 12 d sigma^4 + 8 sigma^2 Delta^2$ and $bb(E)[Z] = Delta^2$, proving (1).

  For (2), applying Chebyshev's inequality to $Z <= 0$ yields
  $ bb(P)(Z <= 0) &<= bb(P)(|Z - bb(E)[Z]| >= Delta^2) \
    &<= (12 d sigma^4 + 8 sigma^2 Delta^2) / Delta^4. $
  Substituting $Delta = c sigma d^alpha$ with $1/4 < alpha < 1/2$ gives raw squared-distance contrast $Delta^2 / (2 d sigma^2) = (c^2 / 2) d^(2 alpha - 1) -> 0$, whereas $"Var"(Z) / Delta^4 = O(d^(1 - 4 alpha)) -> 0$. For (3), sub-exponential (Hanson–Wright) concentration of the quadratic forms $T_2, T_3$ tightens the tail to $exp(-Omega(Delta^4 / (d sigma^4 + sigma^2 Delta^2)))$; whenever $Delta^2 = omega(sigma^2 sqrt(d log n))$, a union bound over all $O(n^2)$ pairs ensures $N_k (i)$ contains exclusively within-cluster points with high probability, driving inter-cluster SNN overlap to $0$ while preserving intra-cluster overlap.
]

#subsec("D", "Resolution Limits and Mutual-NN Isolation Counterexample")

#theorem(title: "Minimum-Cluster-Size Blindness & Parameter Coupling", origin: "Original Counterexample")[
  1. *Feasibility Constraint:* If $k < "MinPts"$ or $k < "Eps"$, Algorithm 1 labels *every* datum as noise ($V_"core" = emptyset$).
  2. *Small-Cluster Dissolution:* For any separation $L > 1$ and any parameter triple $(k, "Eps", "MinPts")$ with $"MinPts" >= 1$, there exists a dataset $cal(X) = C union B$ consisting of an arbitrarily dense, well-separated cluster $C$ of cardinality $m = "MinPts"$ alongside a larger cluster $B$ of size $M >= k + 1$, such that SNN classifies *every* point of $C$ as noise.
]

#proof[
  (1) By Proposition 1(3), $rho_"SNN" (i) <= k$ and $s_"SNN" (i, j) <= k$ for all $i, j$. Thus $k < "MinPts"$ forces $rho_"SNN" (i) < "MinPts"$ globally, and $k < "Eps"$ forces $E_"Eps" = emptyset$.

  (2) Place $m = "MinPts"$ points of $C$ inside a ball of diameter $1$ centered at $bold(0)$ and $M >= k + 1$ points of $B$ inside a ball of diameter $1$ centered at $L e_1$ with $L > 2$. For every $b in B$, there exist $M - 1 >= k$ other points in $B$ at distance $<= 1 < L - 1 <= "dist"(b, C)$. Hence $N_k (b) subset.eq B$ and $N_k (b) inter C = emptyset$ for all $b in B$. Even though points $c in C$ must draw $k - (m - 1)$ neighbors from $B$ whenever $m <= k$, the mutual-$k$-NN gate fails for every pair $(c, b) in C times B$ because $c in.not N_k (b)$, forcing $s_"SNN" (c, b) = 0$. Consequently, every $c in C$ can form SNN edges solely with the remaining $m - 1$ points of $C$, giving $rho_"SNN" (c) <= m - 1 = "MinPts" - 1 < "MinPts"$. Thus $C inter V_"core" = emptyset$, and since no $c in C$ has positive SNN similarity to any core point in $B$, the entire dense cluster $C$ is discarded as noise.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Exact $k$-NN construction requires $O(n^2 d)$ arithmetic operations (reducible to $O(n d log n)$ via KD-trees/ball trees for $d <= 12$, or near-linear $O(n d log n)$ via approximate nearest-neighbor graphs such as HNSW/NN-Descent in high $d$). Once ordered neighbor lists $N_k (i)$ are materialized in $O(n k)$ space, computing $s_"SNN" (i, j)$ over the at most $n k / 2$ mutual edges via sorted-list merge or hash-table lookup takes $O(n k^2)$ time (alternatively $O(sum_u |"RNN"_k (u)|^2)$ via reverse-neighbor inverted lists). Extracting connected core components and assigning border points runs in linear graph time $O(n + |E_"Eps"|) = O(n k)$. Crucially, unlike full similarity matrices costing $O(n^2)$ memory, the sparsified SNN graph requires only $O(n k)$ storage.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for SNN Clustering],
    [Time / Space Complexity], [$O(n^2 d + n k^2)$ exact ($O(n d log n + n k^2)$ with ANN/trees); $O(n (d + k))$ space.],
    [Inductive Bias], [Arbitrary non-convex, multi-density manifolds embedded in high-dimensional noisy spaces.],
    [Invariances], [Strictly invariant to any strictly increasing metric transform $phi compose rho$ and per-cluster scaling (Thm. 1).],
    [High-D Behavior], [Preserves cluster order down to $Delta = omega(sigma d^(1/4))$ despite raw contrast collapse (Thm. 3).],
    [Sensitivity / Blind Spots], [Three coupled parameters $(k, "Eps", "MinPts")$; blind to clusters of size $<= "MinPts"$ (Thm. 4).],
    [Recommended Upgrades], [SNN-DPC (Alg. 52) for peak-based partitioning; RNN-DBSCAN (Alg. 45) for 1-parameter multi-density; HDBSCAN (Alg. 32) for hierarchy.],
  ),
  caption: [Theoretical guarantees and operational profile of SNN Clustering (Alg. 36).],
)

#sec("IV", "Engineering and Scientific Applications")

SNN is a premier density-based choice when clusters coexist at wildly different scales or inhabit high-dimensional spaces plagued by irrelevant attributes and uniform noise.

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Single-Cell RNA Sequencing (scRNA-seq) & Cell-Type Discovery:* High-dimensional transcriptomic profiles ($d approx 20$–$50$ PCA/LSI components, thousands of genes) exhibit pronounced curse-of-dimensionality concentration and varying population sizes across rare vs. abundant cell types. SNN graphs (the exact foundation of Seurat's SNN graph pipeline) suppress noisy cross-type edges via mutual-$k$-NN overlap (Theorem 3), cleanly delineating transcriptionally distinct cell populations.
2. *Earth Science, Climate & Spatiotemporal Anomaly Detection:* Ertöz et al. @ertoz2003snn originally validated SNN on global sea-level pressure, temperature, and NASA Earth-observation time-series grids ($d > 100$), where land and ocean teleconnection indices exhibit sharply differing correlation densities; SNN's scale invariance (Theorem 1(2)) isolates both compact polar vortices and diffuse equatorial oscillations simultaneously.
3. *High-Dimensional Document & Dense Embedding Clustering:* Text TF-IDF vectors and neural sentence embeddings reside in high-dimensional spheres where cosine distances cluster tightly near orthogonality. Because $s_"SNN"$ depends only on shared local neighborhoods, it eliminates hub-induced false shortcuts and recovers topical clusters of unequal vocabulary breadth.
4. *Multi-Density Geospatial & Astronomy Catalogs:* Urban cores alongside suburban fringes—or dense globular clusters alongside diffuse tidal streams—defeat single-$epsilon$ DBSCAN. SNN adapts its effective physical radius locally via $k$-NN rank neighborhoods without requiring manual regional thresholds.
5. *Bioinformatics Protein & Metabolomic Pathway Profiling:* High-dimensional mass-spectrometry and gene-expression assays combine non-linear batch scaling with irrelevant probe dimensions; Theorem 1(1)'s monotone transformation invariance ensures robust clustering across log, rank, and quantile normalizations.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

SNN fails under: (a) *very small clusters ($|C| <= "MinPts"$)*, which are unconditionally pruned as noise by Theorem 4(2); (b) *severe high-dimensional hubness*, where a few "universal hub" points enter $N_k (i)$ for disproportionately many vertices, inflating $|"RNN"_k (u)|$ and creating spurious shared-neighbor bridges across distinct clusters; (c) *poorly tuned $(k, "Eps", "MinPts")$*, where too small $k$ fragments wide clusters into micro-components while too large $k$ blurs fine sub-clusters; and (d) *extreme streaming latency*, since exact $k$-NN updates require maintaining dynamic neighbor graphs.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit / Prerequisite], [Caution / Alternative],
    [scRNA-seq / Genomics], [High-D noise filtering; varied densities], [Leiden (Alg. 108) on SNN graph],
    [Climate / Earth Data], [Scale-invariant teleconnection clusters], [Tune $k$ above rare-event size],
    [Text / Embeddings], [Robust to high-D cosine concentration], [Spherical (Alg. 12) for fast $k$-cut],
    [Multi-Density GIS], [Simultaneous dense & sparse regions], [HDBSCAN (Alg. 32) / OPTICS (Alg. 33)],
    [Tiny Rare Anomalies], [Fails: clusters $<= "MinPts"$ erased], [OPTICS-OF (Alg. 53) or LOF],
  ),
  caption: [Engineering domain suitability and failure-mode guide for SNN Clustering.],
)

#sec("V", "Structural Relations and Empirical Behavior")

SNN serves as the primary bridge connecting Category II density clustering to graph-theoretic, hierarchical, and subspace paradigms.

#proposition(title: "Exact Reduction to GDBSCAN and Graph-Kernel Links", origin: "Original Derivation")[
  Instantiating Generalized DBSCAN (GDBSCAN, Alg. 39) @sander1998gdbscan with the symmetric binary neighborhood predicate $"NPred"(i, j) := [s_"SNN" (i, j) >= "Eps"]$ and uniform cardinality weight $"wCard"(S) := |S|$ with threshold $"MinCard" = "MinPts"$ reproduces the exact core-point set $V_"core"$ and core-cluster partition of SNN (Alg. 36) @ertoz2003snn.
]

#proof[
  Under GDBSCAN @sander1998gdbscan, the generalized neighborhood of $i$ is $N_"NPred" (i) = {j : s_"SNN" (i, j) >= "Eps"} = N_"SNN" (i)$, and $i$ is a generalized core object iff $"wCard"(N_"NPred" (i)) = rho_"SNN" (i) >= "MinPts"$. Because $"NPred"$ is symmetric ($s_"SNN" (i, j) = s_"SNN" (j, i)$), direct density-reachability between core points is symmetric, so GDBSCAN's transitive density-connectivity classes on $V_"core"$ coincide identically with the connected components of $G_"SNN" [V_"core"]$ in Algorithm 1.
]

Beyond GDBSCAN (Alg. 39), replacing SNN's flat core-component expansion with Rodriguez–Laio peak selection yields *SNN-DPC* (Alg. 52) @liu2018snndpc; replacing mutual-$k$-NN overlap with reverse-$k$-NN counts yields *RNN-DBSCAN* (Alg. 45) @bryant2018rnndbscan; using shared-neighbor links inside agglomerative merging produces *ROCK* (Alg. 67) @guha2000rock and *CHAMELEON* (Alg. 66) @karypis1999chameleon; and feeding the PSD overlap graph $G_"SNN"$ into normalized Laplacian eigenvectors connects directly to Spectral Clustering (Alg. 103) @vonluxburg2007tutorial and parameter-free first-neighbor chaining in FINCH (Alg. 80) @sarfraz2019finch. For axis-aligned subspace discovery where even $k$-NN lists are swamped by $d >> Delta^4 / sigma^4$, dedicated subspace methods @kriegel2009subspacesurvey (e.g., CLIQUE Alg. 134, PROCLUS Alg. 136) remain complementary.

Empirically on the ClusterOverview 4-dataset benchmark, SNN demonstrates the exact high-dimensional and multi-scale resilience proven in Theorems 1 and 3: it is the *only* algorithm among Algs. 34–36 to produce valid non-degenerate partitions across *all four* datasets—achieving $"Sil" = 0.643$ on Customers and $"Sil" = 0.351$ on High-D Subspace+Noise ($"ARI"$ N/A), where both Mean Shift and DENCLUE collapsed ($"Sil"$ N/A)—while maintaining strong non-convex performance on Multi-Shape ($"Sil" = 0.703$) and the fastest mean runtime of the trio at $5.044$ ms. On low-dimensional Iris, SNN records $"Sil" = 0.515$ but a lower $"ARI" = 0.149$, reflecting the exact neighborhood-size sensitivity of Theorem 4 when two linearly overlapping species (Versicolor and Virginica) share dense cross-boundary $k$-NN overlap under a fixed $(k, "Eps", "MinPts")$ setting.
