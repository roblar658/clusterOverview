#import "../template/ieee.typ": *

#alg_entry(
  id: 45,
  name: "RNN-DBSCAN (Reverse Nearest Neighbor DBSCAN)",
  category: "Cat. II: Density & Mode-Seeking",
  year: "2018",
  paradigm: "Scale-Free Reverse-kNN Density Graph Clustering",
  complexity: $O(n^2 d)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Classical DBSCAN (Alg. 31, @ester1996dbscan) defines local density by counting points inside a fixed metric ball $B_epsilon (x)$, forcing the practitioner to select a scale-dependent radius $epsilon$ that simultaneously fails whenever clusters exhibit distinct densities. Bryant and Cios @bryant2018rnndbscan resolve both limitations in *RNN-DBSCAN* by replacing fixed-$epsilon$ neighborhoods with *reverse $k$-nearest neighbor ($k$-RNN)* density estimates, reducing the entire algorithm to a single scale-free integer hyperparameter $k in {1, dots, n - 1}$.

Let $(cal(X), "dist")$ be a finite metric dataset of cardinality $n > k$ (typically $cal(X) subset.eq bb(R)^d$ with Euclidean distance). For each $x in cal(X)$, let $N_k (x) subset.eq cal(X) backslash {x}$ denote its set of $k$ nearest neighbors (breaking distance ties deterministically by point index), and let $d_k (x) = max_(y in N_k (x)) "dist"(x, y)$ be its $k$-distance. The *reverse $k$-nearest neighbor set* of $x$ is
$ R_k (x) := {y in cal(X) backslash {x} : x in N_k (y)}, $
namely the set of all points that include $x$ among their $k$ closest neighbors. In the directed $k$-NN graph $G_k = (cal(X), E_k)$ with arcs $(u, v) in E_k <=> v in N_k (u)$, every vertex has uniform out-degree $|N_k (u)| = k$, while its in-degree $|R_k (u)|$ acts as a *dimensionless, locally scale-adaptive density proxy*: points deep inside dense or locally homogeneous clusters receive at least $k$ incoming arcs, whereas points in low-density valleys, sparse fringes, or outlier regions receive strictly fewer than $k$ incoming arcs.

#definition(title: "RNN-DBSCAN Core Graph and Cluster Density Radius")[
  Given $k in {1, dots, n - 1}$:
  1. A point $x in cal(X)$ is a *core point* ($x in cal(C)_k$) if $|R_k (x)| >= k$, and *non-core* otherwise.
  2. The *expansion neighborhood* of any core point $u in cal(C)_k$ is
     $ cal(N)_k (u) := N_k (u) union {v in R_k (u) : |R_k (v)| >= k}. $
  3. The *undirected core graph* is $H_k = (cal(C)_k, E_"core")$, where ${u, v} in E_"core"$ iff $u, v in cal(C)_k$ and $(v in N_k (u) or u in N_k (v))$.
  4. For any core cluster $C subset.eq cal(C)_k$ (a connected component of $H_k$), its *cluster-local density radius* is
     $ rho(C) := max_(u in C) max_(v in N_k (u) inter C) "dist"(u, v). $
     A non-core point $z in.not cal(C)_k$ is assigned as a *border point* to the cluster $C$ of the closest core neighbor $y in N_k (z) inter cal(C)_k$ satisfying $"dist"(z, y) <= rho(C)$ (or nearest core in $N_k (z)$ in simplified variants); otherwise $z$ is labeled *noise*.
]

Algorithm 1 details the two-phase procedure of @bryant2018rnndbscan: Phase 1 extracts the exact connected components of $H_k$ via breadth-first search over $cal(N)_k (u)$, and Phase 2 attaches qualifying border points using each cluster's adaptive reachability radius $rho(C)$.

#algbox(
  title: "RNN-DBSCAN Clustering Algorithm",
  inputs: [Dataset $cal(X) = {x_1, dots, x_n}$ with metric $"dist"$, integer $k in {1, dots, n - 1}$.],
  outputs: [Cluster partition $cal(C)$ and noise set $cal(N)$.],
)[
  1. *Graph Construction:* For each $x in cal(X)$, compute $N_k (x)$ and invert edges to populate $R_k (y) = {x : y in N_k (x)}$.
  2. Classify core set $cal(C)_k <- {x in cal(X) : |R_k (x)| >= k}$; initialize labels $ell(x) <- -1$, $"cid" <- 0$.
  3. *Phase 1 (Core Expansion):* *For* each $x in cal(C)_k$ with $ell(x) = -1$ *do*:
  4. #h(1em) Set $ell(x) <- "cid"$, push $x$ onto queue $Q$.
  5. #h(1em) *While* $Q$ is non-empty *do*: pop $u <- Q$, and *for* each $v in cal(N)_k (u)$ *do*:
  6. #h(2em) *If* $v in cal(C)_k$ and $ell(v) = -1$, set $ell(v) <- "cid"$ and push $v$ onto $Q$.
  7. #h(1em) Increment $"cid" <- "cid" + 1$.
  8. *Phase 2 (Border Assignment):* For each cluster $C_r$, compute $rho(C_r) = max_(u in C_r) max_(v in N_k (u) inter C_r) "dist"(u, v)$.
  9. *For* each non-core $z in.not cal(C)_k$ *do*: find closest $y in N_k (z) inter cal(C)_k$ with $"dist"(z, y) <= rho(C_(ell(y)))$; if found set $ell(z) <- ell(y)$, else mark $z in cal(N)$.
  10. *Return* $(cal(C), cal(N))$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five formal theoretical results for RNN-DBSCAN: (1) an exact in-degree conservation law and geometric lower bound on the core-set cardinality $|cal(C)_k|$; (2) permutation invariance of the core and Phase-2 partitions; (3) exact global similarity, monotone metric, and block-wise multi-density scale invariance; (4) an exact analytical calculation proving boundary demotion on uniform clusters and small-cluster dissolution; and (5) asymptotic density-curvature behavior and high-dimensional hubness analysis.

#subsec("A", "In-Degree Conservation and Core-Cardinality Bounds")

#lemma(title: "Exact Sum-Rule Conservation and Geometric Core Lower Bound", origin: "Original Proof")[
  For any $n$-point metric dataset ($n > k$) with maximum in-degree $Delta_"max" = max_(x in cal(X)) |R_k (x)|$:
  1. $sum_(x in cal(X)) |R_k (x)| = n k$, so the core threshold $k$ equals the exact empirical mean of $|R_k (x)|$, and $Delta_"max" >= k$.
  2. The number of core points satisfies the sharp bound
     $ |cal(C)_k| >= ceil(n / (Delta_"max" - k + 1)) >= 1. $
  3. In Euclidean space $bb(R)^d$ (in general position), $Delta_"max" <= tau_d k$, where $tau_d$ is the minimum number of spherical cones of angular diameter $< pi / 3$ covering $bb(S)^(d-1)$ (e.g., $tau_1 = 2$, $tau_2 = 6$, $tau_d = 2^(Theta(d))$), yielding $|cal(C)_k| >= ceil(n / ((tau_d - 1) k + 1))$.
]

#proof[
  (1) Counting the directed edges of $G_k = (cal(X), E_k)$ by source versus target gives
  $ sum_(x in cal(X)) |R_k (x)| = |E_k| = sum_(x in cal(X)) |N_k (x)| = n k. $
  Dividing by $n$ shows the mean in-degree is identically $k$, which immediately forces $Delta_"max" = max_x |R_k (x)| >= k$.

  (2) Let $c = |cal(C)_k|$. Each core point has $|R_k (x)| <= Delta_"max"$, while each of the $n - c$ non-core points has $|R_k (x)| <= k - 1$. Summing over $cal(X)$ gives
  $ n k = sum_(x in cal(C)_k) |R_k (x)| + sum_(x in.not cal(C)_k) |R_k (x)| <= c Delta_"max" + (n - c) (k - 1). $
  Rearranging yields $n <= c (Delta_"max" - k + 1)$. Since $Delta_"max" - k + 1 >= 1$, dividing and taking ceilings proves $c >= ceil(n / (Delta_"max" - k + 1)) >= 1$.

  (3) Partition $bb(S)^(d-1)$ into $tau_d$ open spherical cones centered at $x$ each having angular diameter $< pi / 3$. Suppose by contradiction that some cone $K$ contains $k + 1$ points $y_1, dots, y_(k+1) in R_k (x)$, ordered in general position as $0 < norm(y_1 - x)_2 < dots < norm(y_(k+1) - x)_2$. For every $j in {1, dots, k}$, the angle $theta_j = angle(y_j - x, y_(k+1) - x) < pi / 3$ satisfies $cos theta_j > 1 / 2$, so
  $ norm(y_(k+1) - y_j)_2^2 &= norm(y_(k+1) - x)_2^2 + norm(y_j - x)_2^2 \
    &quad - 2 norm(y_(k+1) - x)_2 norm(y_j - x)_2 cos theta_j \
    &< norm(y_(k+1) - x)_2^2 + norm(y_j - x)_2 (norm(y_j - x)_2 - norm(y_(k+1) - x)_2) \
    &< norm(y_(k+1) - x)_2^2. $
  Thus all $k$ points $y_1, dots, y_k$ are strictly closer to $y_(k+1)$ than $x$ is, contradicting $x in N_k (y_(k+1))$. Hence each cone holds at most $k$ points of $R_k (x)$, giving $|R_k (x)| <= tau_d k$ and completing Claim 3.
]

#theorem(title: "Order-Independence of RNN-DBSCAN Core and Border Partitions", origin: "Original Proof")[
  Assuming strict pairwise distances (no ties), both the Phase-1 core clusters and the full Phase-2 clustering of Algorithm 1 are *strictly invariant* to any permutation of the input order of $cal(X)$.
]

#proof[
  For any two core points $u, v in cal(C)_k$, Algorithm 1 traverses an edge from $u$ to $v$ if and only if $v in cal(N)_k (u) inter cal(C)_k = (N_k (u) inter cal(C)_k) union (R_k (u) inter cal(C)_k)$. Because $v in R_k (u) <=> u in N_k (v)$, this adjacency relation is identically
  $ v in N_k (u) or u in N_k (v), $
  which is strictly symmetric in $(u, v)$. Consequently, Phase-1 BFS extracts precisely the connected components of the undirected graph $H_k = (cal(C)_k, E_"core")$, which depend solely on the unordered set $cal(X)$ and not on loop order. Since each cluster radius $rho(C_r)$ depends only on the component $C_r$, and Phase 2 assigns each non-core $z in.not cal(C)_k$ deterministically to the unique closest qualifying core point $y in N_k (z) inter cal(C)_k$, the entire partition $(cal(C), cal(N))$ is permutation-invariant (eliminating the eager-assignment order dependence noted for standard DBSCAN border points @schubert2017dbscanrevisited).
]

#subsec("B", "Scale Invariance and Arbitrary Multi-Density Resolution")

#theorem(title: "Monotone-Metric and Block-Wise Multi-Density Invariance", origin: "Original Proof")[
  1. (*Monotone Distance & Similarity Invariance:*) Let $"dist"' = psi compose "dist"$ for any strictly increasing function $psi: [0, infinity) -> [0, infinity)$ (including global Euclidean similarities $Phi(x) = alpha Q x + b$, $alpha > 0$, $Q in O(d)$). Then $N_k (x)$, $R_k (x)$, $cal(C)_k$, $H_k$, and the complete RNN-DBSCAN clustering are *identically invariant*.
  2. (*Arbitrary Multi-Density Contrast Invariance:*) Let $cal(X) = union.big_(r=1)^m C_r subset.eq bb(R)^d$ consist of $m$ clusters ($|C_r| > k$) where each cluster $C_r$ is independently rescaled around its centroid by an arbitrary factor $alpha_r > 0$ (producing arbitrary density ratios $(alpha_s / alpha_r)^d$). Provided inter-cluster separations satisfy
     $ "dist"(C_r, C_s) > max_(x in C_r) d_k^((C_r)) (x) quad forall r != s, $
     where $d_k^((C_r)) (x)$ is the intra-cluster $k$-NN distance in $C_r$, the entire RNN-DBSCAN clustering is *strictly invariant* to the choice of $(alpha_1, dots, alpha_m)$.
]

#proof[
  (1) Because $psi$ is strictly increasing, $"dist"(x, y) < "dist"(x, z) <=> "dist"'(x, y) < "dist"'(x, z)$ for all $x, y, z in cal(X)$. Hence every $k$-NN set $N_k (x)$ is unchanged, which fixes $R_k (x)$, $cal(C)_k$, and $H_k$. In Phase 2, for any core cluster $C$, the inequality $"dist"(z, y) <= rho(C) = max_(u in C, v in N_k (u) inter C) "dist"(u, v)$ transforms under $psi$ to $"dist"'(z, y) <= psi(rho(C)) = rho'(C)$, leaving every border/noise decision identically unchanged!

  (2) Under the separation condition $"dist"(C_r, C_s) > max_(x in C_r) d_k^((C_r)) (x)$ for all $s != r$, every point $x in C_r$ has at least $k$ neighbors inside $C_r$ strictly closer than every point in $cal(X) backslash C_r$. Therefore $N_k (x) subset.eq C_r$ for every $x in C_r$, and consequently $R_k (x) subset.eq C_r$ as well. The directed graph $G_k$ decomposes into $m$ completely disconnected subgraphs $G_k [C_r]$. Inside each component $C_r$, applying a cluster-specific dilation $alpha_r > 0$ multiplies all intra-$C_r$ distances by $alpha_r$, which by Claim 1 leaves the internal $N_k$, $R_k$, core graph $H_k [C_r]$, and Phase-2 radius test $"dist"(z, y) <= rho(C_r)$ strictly invariant.
]

#subsec("C", "Exact Boundary Erosion and Small-Cluster Dissolution")

#theorem(title: "Uniform-Cluster Boundary Demotion and Small-Cluster Collapse", origin: "Original Counterexample")[
  1. (*Exact Boundary Demotion:*) Let $cal(X) = {x_j = j delta : j = 1, dots, n} subset.eq bb(R)$ be a uniform 1-D cluster with even $k = 2 m$ and $n >= 4 m + 2$. Then the reverse-neighbor counts satisfy *exactly*:
     $ |R_(2 m) (x_r)| = m + r - 1 < 2 m quad forall r in {1, dots, m}, $
     and symmetrically $|R_(2 m) (x_(n + 1 - r))| = m + r - 1 < 2 m$, while the sub-boundary layer $r in {m + 2, dots, 2 m + 1}$ exhibits elevated $|R_(2 m) (x_r)| = m + r - 1 > 2 m$. Thus the $m$ outermost points on each boundary are systematically demoted to non-core despite having uniform local spacing.
  2. (*Small-Cluster Dissolution:*) Any isolated cluster $C subset.eq cal(X)$ of size $|C| <= k$ separated from the remaining data has $max_(x in C) |R_k (x)| <= |C| - 1 < k$ whenever external points do not reach $C$, causing $C$ to contain *zero core points* and dissolve entirely into noise (or merge spuriously into a neighboring core).
]

#proof[
  (1) For an interior point $x_i$ ($m + 1 <= i <= n - m$), $N_(2 m) (x_i) = {x_(i - m), dots, x_(i - 1), x_(i + 1), dots, x_(i + m)}$. For a left-boundary point $x_i$ ($1 <= i <= m$), there are only $i - 1 < m$ points to its left, so its $2 m$ closest neighbors are ${x_1, dots, x_(2 m + 1)} backslash {x_i}$. Fix any $r in {1, dots, m}$ and count all $i != r$ with $x_r in N_(2 m) (x_i)$:
  - Every $i in {1, dots, m} backslash {r}$ includes $x_r$ (yielding $m - 1$ arcs).
  - Among interior points $i >= m + 1$, $x_r in N_(2 m) (x_i)$ iff $i - r <= m$, i.e., $i in {m + 1, dots, m + r}$ (yielding $r$ arcs).
  Summing gives $|R_(2 m) (x_r)| = (m - 1) + r = m + r - 1 < 2 m$ for all $1 <= r <= m$. For $r in {m + 2, dots, 2 m + 1}$, every left-boundary point $i in {1, dots, m}$ reaches $x_r$, plus interior indices $i in {m + 1, dots, r + m} backslash {r}$ (count $r - 1$), giving $|R_(2 m) (x_r)| = m + r - 1 > 2 m$—balancing Lemma 1's conservation sum.

  (2) If $|C| <= k$ and no point outside $C$ places any $x in C$ in its $k$-NN set, then every $x in C$ can receive incoming arcs only from $C backslash {x}$, so $|R_k (x)| <= |C| - 1 < k$. Hence $C inter cal(C)_k = emptyset$, and $C$ cannot seed a cluster; meanwhile, each $x in C$ must draw $k - (|C| - 1) >= 1$ forward neighbors from outside $C$, risking false attachment if external distances fall below a foreign cluster's $rho(C')$.
]

#subsec("D", "Asymptotic Density Curvature and High-D Hubness")

#theorem(title: "Asymptotic Curvature Bias and High-Dimensional Hubness", origin: "Classical — Cited")[
  1. (*Mode Amplification & Valley Suppression:*) Let $X_1, dots, X_n$ be i.i.d. from a smooth density $f in C^2 (bb(R)^d)$ with $f(x) > 0$, and let $k / n -> 0$ with $k -> infinity$. At any critical point $nabla f(x) = bold(0)$ with characteristic $k$-NN radius $r_0 = (k / (n V_d f(x)))^(1/d)$, the expected reverse-neighbor count satisfies @bryant2018rnndbscan
     $ bb(E)[|R_k (x)|] = k (1 - (r_0^2 Delta f(x)) / (2 d f(x)) + o(r_0^2)). $
     Thus $|R_k (x)|$ exceeds $k$ at local density maxima ($Delta f(x) < 0$) and falls below $k$ in density valleys ($Delta f(x) > 0$).
  2. (*High-Dimensional Hubness Degradation:*) As $d -> infinity$, pairwise distance concentration @beyer1999nearest induces severe right-skewness in $|R_k (x)|$ (hubness), driving $Delta_"max"$ upward toward Lemma 1's $tau_d k = 2^(Theta(d)) k$ ceiling and shrinking the core fraction $|cal(C)_k| / n$.
]

#proof(label: "Proof Synthesis (Cited)")[
  Around a critical point $x$ with Hessian $H = nabla^2 f(x)$, the local $k$-NN ball radius $r_k (y)$ enclosing mass $p = k / n$ satisfies $V_d r_k (y)^d [f(y) + r_k (y)^2 Delta f(y) / (2 (d + 2))] approx p$. Writing $y = x + r theta$ ($theta in bb(S)^(d-1)$) and expanding $f(x + u) = f(x) + 1 / 2 u^T H u + o(norm(u)_2^2)$ gives the boundary radius $R(theta) = r_0 (1 - r_0^2 theta^T H theta / (2 d f(x)) - r_0^2 Delta f(x) / (2 d (d + 2) f(x)))$. Integrating $f(x + r theta) r^(d-1)$ up to $R(theta)$ and averaging $theta^T H theta$ to $Delta f(x) / d$ cancels the $(d + 2)^(-1)$ terms identically, leaving $bb(E)[|R_k (x)|] = n integral f(y) bb(1)_(norm(y - x)_2 <= r_k (y)) d y = k (1 - r_0^2 Delta f(x) / (2 d f(x)) + o(r_0^2))$. In high $d$, concentration of Euclidean distances causes central points near the global mean to dominate $k$-NN lists across the dataset @beyer1999nearest, inflating $Delta_"max"$ and leaving peripheral cluster points with $|R_k| < k$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Computing exact $k$-NN sets takes $O(n^2 d)$ time via brute-force pairwise distances and linear-time selection (or $O(n k log n)$ in low intrinsic dimension via $k$-d, ball, or cover trees). Storing both $N_k$ and $R_k$ requires strictly $O(n k)$ memory by Lemma 1's conservation identity $sum_x |R_k (x)| = n k$. Once $G_k$ is built, Phase-1 BFS visits each edge of $N_k$ and $R_k$ at most twice in $O(n k)$ time, and Phase-2 radius computation plus border assignment also runs in $O(n k)$ time. Thus the entire clustering logic is strictly linear in $n k$, shifting almost all computational effort into $k$-NN graph construction—which also makes RNN-DBSCAN directly compatible with approximate $k$-NN indices (HNSW, NN-Descent) for million-scale datasets.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for RNN-DBSCAN],
    [Time / Space Complexity], [$O(n^2 d)$ exact ($O(n k log n)$ low-$d$ indexed) + $O(n k)$ clustering; $O(n k)$ graph space.],
    [Hyperparameter Burden], [Single scale-free integer $k$ (eliminates metric $epsilon$ tuning completely).],
    [Multi-Density Capability], [Exact invariance under arbitrary block-wise density contrast (Theorem 3).],
    [Determinism / Invariances], [Strictly order-independent (Theorem 2); invariant to any monotone metric warp $psi compose "dist"$.],
    [Boundary / Noise Behavior], [Natural outlier rejection ($|R_k| < k$); erodes $m$-wide cluster edges and clusters of size $<= k$ (Theorem 4).],
    [Recommended Companions], [HDBSCAN (Alg. 32) for hierarchical multi-density; SNN (Alg. 36) for high-$d$ hubness mitigation.],
  ),
  caption: [Operational and theoretical profile of RNN-DBSCAN.],
)

#sec("IV", "Engineering and Scientific Applications")

RNN-DBSCAN excels whenever cluster densities span multiple orders of magnitude, feature scales are non-Euclidean or monotonically calibrated, and operators cannot reliably tune a global distance threshold $epsilon$.

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Heterogeneous-Density Geospatial & Mobility Analytics:* Urban cores, suburban corridors, and rural settlements exhibit $10^2$–$10^3 times$ density contrast. By Theorem 3(2), a single integer $k$ cleanly recovers both dense downtown hotspots and diffuse regional clusters while flagging isolated GPS outliers as noise.
2. *Single-Cell Genomics (scRNA-seq) & Cell-Type Discovery:* Cellular subpopulations in PCA/UMAP embeddings span abundant major lineages and compact rare states with non-convex developmental trajectories. Because RNN-DBSCAN depends purely on neighbor rank orderings (Theorem 3(1)), it is immune to monotone metric distortions introduced by nonlinear manifold embeddings.
3. *Astrophysical Catalogs & Multi-Scale Structure Finding:* Dense globular clusters coexist alongside diffuse stellar associations and tidal streams; reverse-$k$NN in-degree naturally elevates interior cores (Theorem 5) and suppresses inter-cluster field stars without multi-pass $epsilon$-sweeping.
4. *Cybersecurity & Financial Anomaly Isolation:* Normal behavioral clusters vary widely in compactness across enterprise user tiers; genuine anomalies sit in sparse fringes where surrounding points' $k$-NN balls fail to reach them ($|R_k| << k$), yielding clean noise separation.
5. *Non-Metric & Learned-Embedding Retrieval Clustering:* In deep metric learning or bioinformatics sequence databases where similarity scores undergo arbitrary monotone rescaling, rank-invariance guarantees identical clusterings without recalibrating metric radii.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Practitioners must guard against four failure regimes: (a) *micro-clusters of cardinality $<= k$*, which dissolve completely (Theorem 4(2)); (b) *high-dimensional hubness* ($d >= 25$), where right-skewed $|R_k|$ shrinks the core set (Lemma 1, Theorem 5—mitigated by SNN Alg. 36 @ertoz2003snn or SNN-DPC Alg. 52 @liu2018snndpc); (c) *chaining across wide bridges* containing $>= k$ points; and (d) *boundary shrinkage* on thin sheets where outermost points drop to $|R_k| < k$ (Theorem 4(1), largely mitigated by Phase-2 border attachment).

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Why RNN-DBSCAN Fits], [Caution / Alternative],
    [Multi-Density GIS], [Single $k$; arbitrary contrast], [Ensure min cluster size $> k$],
    [scRNA-seq / UMAP], [Monotone-warp rank invariance], [Use Leiden for huge graphs],
    [Astrophysics], [Separates dense cores & streams], [Compare with HDBSCAN (Alg. 32)],
    [Anomaly Detection], [Outliers have $|R_k| << k$], [Check Phase-2 $rho(C)$ cutoff],
    [High-D Subspaces], [Hubness shrinks core fraction], [Use SNN (Alg. 36) or OPTICS (Alg. 33)],
  ),
  caption: [Domain suitability and engineering guide for RNN-DBSCAN.],
)

#sec("V", "Structural Relations and Empirical Behavior")

RNN-DBSCAN bridges generalized density reachability and $k$-NN graph-theoretic clustering.

#proposition(title: "GDBSCAN Specialization and Strict k-NN Graph Hierarchy", origin: "Original Unified Proof")[
  1. (*GDBSCAN Embedding:*) Phase 1 of RNN-DBSCAN is a strict instance of Generalized DBSCAN (Alg. 39, @sander1998gdbscan) with neighborhood predicate $N_"Pred" (u, v) equiv (v in N_k (u)) or (u in N_k (v) and |R_k (v)| >= k)$ and density predicate $"Card"(u) = |R_k (u)| >= k$.
  2. (*Strict Core-Edge Sandwich:*) Let $E_"mut" = {{u, v} : u in N_k (v) and v in N_k (u)}$ and $E_"sym" = {{u, v} : u in N_k (v) or v in N_k (u)}$ denote the mutual and symmetric $k$-NN graphs. Restricted to $cal(C)_k$, RNN-DBSCAN's core graph obeys
     $ E_"mut" |_(cal(C)_k) subset.eq E_"core" = E_"sym" |_(cal(C)_k), $
     with strict inclusion $E_"mut" |_(cal(C)_k) subset.neq E_"core"$ on heterogeneous-density gradients.
]

#proof[
  (1) Substituting $N_"Pred"$ and $"Card"$ directly into GDBSCAN's core-expansion loop @sander1998gdbscan yields core condition $|R_k (u)| >= k$ and expansion set $cal(N)_k (u) = {v : N_"Pred" (u, v)}$, matching Algorithm 1 Steps 2–7 identically.

  (2) By Definition 1(3), ${u, v} in E_"core"$ iff $u, v in cal(C)_k$ and $(v in N_k (u) or u in N_k (v))$, which is precisely $E_"sym" |_(cal(C)_k)$ and manifestly contains $E_"mut" |_(cal(C)_k)$. To prove strictness, place a dense cluster $C_1$ adjacent to a slightly sparser core region $C_2$ so that a boundary core point $u in C_2$ reaches into $C_1$ ($v in N_k (u)$) while $v$'s tighter $k$-NN ball remains inside $C_1$ ($u in.not N_k (v)$); then ${u, v} in E_"core" backslash E_"mut"$.
]

This sandwich explains RNN-DBSCAN's positioning relative to adjacent density methods: it preserves connectivity across gentle density gradients where mutual $k$-NN fractures, yet prunes single-link outlier chains that corrupt full symmetric $k$-NN graphs. Compared to OPTICS (Alg. 33, @ankerst1999optics), HDBSCAN (Alg. 32, @campello2013hdbscan, @campello2015hierarchical), VDBSCAN (Alg. 50, @liu2007vdbscan), and DeBaCl (Alg. 46, @kent2013debacl), RNN-DBSCAN yields a flat multi-density clustering directly in one graph pass without constructing a dendrogram or reachability hierarchy.

Empirically on the ClusterOverview 4-dataset benchmark, RNN-DBSCAN attains *Iris* $"Sil" = 0.586$, $"ARI" = 0.567$ (matching DENCLUE 2.0 and strongly outperforming K-DBSCAN's $0.392 / 0.333$) and *Multi-Shape* $"Sil" = 0.741$ ($"ARI" = "N/A"$), confirming accurate arbitrary-shape manifold tracking without Voronoi fracture. Unlike DENCLUE 2.0, it produces valid multi-cluster partitions across all four datasets—including *Customers* ($"Sil" = 0.364$) and *High-D Subspace+Noise* ($"Sil" = 0.301$, $"ARI" = "N/A"$), where Lemma 1 guarantees non-empty cores even though high-dimensional hubness (Theorem 5) moderates silhouette cohesion. Its mean runtime of $1.456$ ms sits cleanly between K-DBSCAN ($0.351$ ms) and DENCLUE 2.0 ($7.234$ ms), dominated by exact $O(n^2 d)$ $k$-NN construction prior to the linear $O(n k)$ BFS pass.
