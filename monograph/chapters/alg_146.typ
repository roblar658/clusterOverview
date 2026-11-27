#import "../template/ieee.typ": *

#alg_entry(
  id: 146,
  name: "FIRES (Filter-Refinement Subspace Clustering)",
  category: "Cat. VI: Grid, Subspace & Projection",
  year: "2005",
  paradigm: "Filter-Refinement / 1D Base-Cluster Merging",
  complexity: $O(d n log n + n d^2 + m n^2 d_"avg")$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

A central computational bottleneck of bottom-up subspace clustering algorithms such as CLIQUE @agrawal1998clique and SUBCLU @kailing2004subclu is their reliance on Apriori-style level-wise enumeration across the subspace lattice $2^({1, dots, d})$: recovering even a single $s$-dimensional cluster forces the algorithm to materialize and test all $2^s - 1$ non-empty projections, incurring exponential $Omega(2^s)$ runtime. To break this exponential barrier, Kriegel, Kröger, Renz, and Wurst @kriegel2005fires introduced *FIRES* (FIlter REfinement Subspace clustering), a generic three-phase framework that discovers 1D *base clusters* independently on every coordinate, merges mutually overlapping base clusters directly into high-dimensional *subspace cluster approximations* via a similarity graph (filter step), and applies an exact multidimensional clustering algorithm exclusively inside the surviving candidate subspaces (refinement step).

Let $cal(D) = {bold(x)_1, dots, bold(x)_n} subset bb(R)^d$ be a $d$-dimensional dataset, and let $pi_S : bb(R)^d -> bb(R)^(abs(S))$ denote orthogonal projection onto any non-empty attribute subset $S subset.eq {1, dots, d}$. For a 1D clustering operator $cal(A)_1$ (canonically 1D DBSCAN @ester1996dbscan or 1D $k$-means @lloyd1982least) applied to each projection $pi_({k}) (cal(D))$, let $cal(B)_k = {B_(k, 1), dots, B_(k, b_k)}$ be the resulting pairwise-disjoint family of 1D base clusters of cardinality at least $mu_"min"$, and let $cal(B) = union_(k=1)^d cal(B)_k$ with $M = abs(cal(B))$ and dimension map $"dim"(B) = k$ for $B in cal(B)_k$.

#definition(title: "Base-Cluster Similarity Graph and Filter-Refinement Approximation")[
  For two base clusters $B_1, B_2 in cal(B)$ with $"dim"(B_1) != "dim"(B_2)$, define their *overlap similarity* as
  $ s(B_1, B_2) = abs(B_1 inter B_2) $
  (or its normalized variant $tilde(s)(B_1, B_2) = abs(B_1 inter B_2) / min(abs(B_1), abs(B_2))$). For integer pruning parameter $k_"ms" >= 1$ and minimum overlap threshold $tau_"sim" > 0$:
  1. The *$k_"ms"$-most-similar set* $"MS"_(k_"ms") (B)$ contains the $k_"ms"$ base clusters $B' in cal(B)$ ($"dim"(B') != "dim"(B)$) maximizing $s(B, B')$.
  2. The *best-merge set* is $"BM"(B) = {B' in "MS"_(k_"ms") (B) : s(B, B') >= tau_"sim"}$.
  3. The *base-cluster similarity graph* $G_"BC" = (cal(B), E, w)$ places an undirected edge ${B_1, B_2} in E$ with weight $w(B_1, B_2) = s(B_1, B_2)$ iff $B_2 in "BM"(B_1)$ and $B_1 in "BM"(B_2)$ (mutual best-merge).
  A *subspace cluster approximation* is a graph cluster $cal(M)_r subset.eq cal(B)$ of $G_"BC"$ (with distinct coordinate dimensions) inducing candidate subspace $S_r = {"dim"(B) : B in cal(M)_r}$ and candidate support $hat(C)_r = inter_(B in cal(M)_r) B$.
]

Because 1D projections of multiple distinct high-dimensional clusters can collide onto the same coordinate interval, FIRES optionally incorporates a *base-cluster splitting* step prior to graph construction: a base cluster $B in cal(B)_k$ that exhibits bimodal/multi-modal density or high intersection with multiple mutually disjoint clusters on other dimensions is split before building $G_"BC"$. Finally, Phase 3 discards degenerate candidates ($abs(S_r) < d_"min"$ or $abs(hat(C)_r) < mu_"min"$) and runs a multidimensional clusterer $cal(A)_"ref"$ (typically DBSCAN @ester1996dbscan) restricted to each surviving subspace $S_r$.

#algbox(
  title: "FIRES: Three-Phase Filter-Refinement Subspace Clustering",
  inputs: [Dataset $cal(D) subset bb(R)^d$; 1D base clusterer $cal(A)_1$; top-$k$ parameter $k_"ms"$; overlap threshold $tau_"sim"$; min cluster size $mu_"min"$; refinement clusterer $cal(A)_"ref"$.],
  outputs: [Set of refined subspace clusters $cal(R) = {(C_1, S_1), dots, (C_q, S_q)}$.],
)[
  1. *Phase 1 — 1D Base-Cluster Generation:* For each $k = 1, dots, d$, run $cal(A)_1$ on $pi_({k}) (cal(D))$ to obtain pairwise-disjoint 1D base clusters $cal(B)_k$ with $abs(B) >= mu_"min"$. Form $cal(B) = union_(k=1)^d cal(B)_k$ (optionally splitting large overlapping 1D clusters).
  2. *Phase 2 — Filter Step (Subspace Approximation):*
     - Build inverted object-to-base-cluster lists $cal(I)(i) = {B in cal(B) : bold(x)_i in B}$ for $i = 1, dots, n$, and accumulate all nonzero intersections $s(B_1, B_2) = abs(B_1 inter B_2)$.
     - For each $B in cal(B)$, extract its best-merge set $"BM"(B) subset.eq "MS"_(k_"ms") (B)$ satisfying $s(B, B') >= tau_"sim"$, and construct the mutual similarity graph $G_"BC" = (cal(B), E, w)$.
     - Partition $G_"BC"$ into dense/connected base-cluster groups $cal(M)_1, dots, cal(M)_m$ (using greedy high-weight clustering on $G_"BC"$), generating candidate subspaces $S_r = {"dim"(B) : B in cal(M)_r}$ and approximate point sets $hat(C)_r$.
     - Prune candidates failing minimum subspace dimensionality or size thresholds.
  3. *Phase 3 — Refinement Step:* For each surviving candidate $(hat(C)_r, S_r)$, execute $cal(A)_"ref"$ on $pi_(S_r) (cal(D))$ (or on the local neighborhood of $hat(C)_r$ in $S_r$) to extract exact arbitrary-shaped clusters $(C, S_r)$ and add them to $cal(R)$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We prove the projection-containment foundation of 1D base filtering, establish a sharp $O(n d^2)$ inverted-list complexity separation over exponential Apriori lattice search, derive finite-sample hypergeometric Chernoff bounds for spurious-merge suppression, and construct an explicit 1D projection-masking counterexample.

#subsec("A", "Projection Containment and Base-Cluster Recall")

#theorem(title: "Exact 1D Projection Containment for Monotone Density Clusters", origin: "Classical — Cited")[
  Let $cal(A)_1$ and $cal(A)_"ref"$ be $(epsilon, "MinPts")$-DBSCAN @ester1996dbscan under any $L_p$ metric ($p in [1, infinity]$). If $C subset.eq cal(D)$ is a maximal density-connected cluster in subspace $S subset.eq {1, dots, d}$ with $abs(C) >= mu_"min"$, then for every coordinate $k in S$ there exists a unique 1D base cluster $B_k in cal(B)_k$ such that $C subset.eq B_k$. Consequently,
  $ C subset.eq inter_(k in S) B_k, quad "and" quad s(B_(k_1), B_(k_2)) >= abs(C) quad forall k_1, k_2 in S. $
]

#proof[
  For any two objects $bold(x), bold(y) in cal(D)$, any subspace $S$, and any coordinate $k in S$, the monotonicity of $L_p$ norms gives
  $ abs(x_k - y_k) <= ( sum_(j in S) abs(x_j - y_j)^p )^(1/p) = norm(pi_S (bold(x)) - pi_S (bold(y)))_p $
  (with $abs(x_k - y_k) <= max_(j in S) abs(x_j - y_j)$ for $p = infinity$). Therefore every $S$-$epsilon$-neighbor of $bold(x)$ is a 1D-$epsilon$-neighbor on coordinate $k$: $N_epsilon^S (bold(x)) subset.eq N_epsilon^({k}) (bold(x))$. It follows immediately that: (1) every $S$-core object ($abs(N_epsilon^S (bold(x))) >= "MinPts"$) is a 1D-core object on attribute $k$; (2) direct $S$-density-reachability $bold(y) in N_epsilon^S (bold(x))$ implies direct 1D-density-reachability; and (3) by induction along finite chains, $S$-density-connectivity of $C$ implies 1D-density-connectivity on coordinate $k$ @kailing2004subclu @kriegel2005fires. Hence $C$ is wholly contained in a single maximal 1D density-connected component $B_k in cal(B)_k$, which is unique because distinct core-connected components are disjoint (modulo border points assigned to a single base cluster). Taking the intersection over all $k in S$ yields $C subset.eq inter_(k in S) B_k$, so $abs(B_(k_1) inter B_(k_2)) >= abs(C)$ for all $k_1, k_2 in S$.
]

#subsec("B", "Exact Polynomial Separation from Apriori Lattice Enumeration")

A frequently overlooked algorithmic property of FIRES is that computing *all* exact pairwise base-cluster intersections scales strictly quadratically in $d$ and linearly in $n$, independent of the number of base clusters $M$.

#theorem(title: "Exact O(n d^2) Inverted-List Intersection Bound vs. Exponential Apriori Lower Bound", origin: "Original Proof")[
  Suppose each 1D base-cluster family $cal(B)_k$ ($k = 1, dots, d$) consists of pairwise-disjoint subsets of $cal(D)$. Then:
  (i) The sum of all pairwise base-cluster intersection sizes across all distinct-dimension pairs satisfies
  $ sum_(1 <= k_1 < k_2 <= d) sum_(B_1 in cal(B)_(k_1)) sum_(B_2 in cal(B)_(k_2)) abs(B_1 inter B_2) <= n binom(d, 2) = (n d (d - 1)) / 2. $
  Consequently, the exact sparse overlap matrix $s(B_1, B_2)$ has at most $n d (d - 1) / 2$ nonzero entries and is constructed in $O(n d^2)$ worst-case time, independent of $M = abs(cal(B))$.
  (ii) Combined with 1D sorting-sweep DBSCAN in Phase 1 ($O(d n log n)$), Phases 1–2 execute in $O(d n log n + n d^2 + M k_"ms" log M)$ time, whereas any level-wise Apriori subspace algorithm (e.g., CLIQUE @agrawal1998clique or SUBCLU @kailing2004subclu) requires $Omega(2^s)$ time whenever a single true cluster exists in an $s$-dimensional subspace.
]

#proof[
  (i) For each object $bold(x)_i in cal(D)$, let $cal(I)(i) = {B in cal(B) : bold(x)_i in B}$ denote its inverted base-cluster list, and let $d_i = abs(cal(I)(i))$. Because the base clusters in each dimension $k$ are pairwise disjoint, $abs(cal(I)(i) inter cal(B)_k) <= 1$ for every $k in {1, dots, d}$, so $d_i <= d$. Double-counting triples $(i, B_1, B_2)$ with $bold(x)_i in B_1 inter B_2$ and $"dim"(B_1) < "dim"(B_2)$ gives
  $ sum_(k_1 < k_2) sum_(B_1 in cal(B)_(k_1)) sum_(B_2 in cal(B)_(k_2)) abs(B_1 inter B_2) = sum_(i=1)^n binom(d_i, 2) <= sum_(i=1)^n binom(d, 2) = n d (d - 1) / 2. $
  Iterating over $i = 1, dots, n$ and incrementing hash-table accumulators for every unordered pair ${B_1, B_2} subset.eq cal(I)(i)$ executes exactly $sum_(i=1)^n binom(d_i, 2) <= n d (d - 1) / 2$ $O(1)$-expected updates and emits at most that many nonzero edges.
  (ii) In 1D, sorting $pi_({k}) (cal(D))$ takes $O(n log n)$ per attribute, and a two-pointer sliding window of width $2 epsilon$ identifies all 1D core points and connected components in $O(n)$ linear time, totaling $O(d n log n)$ across $d$ attributes. Selecting the top-$k_"ms"$ neighbors per base cluster over at most $min(M^2, n d^2)$ nonzero edges takes $O(n d^2 + M k_"ms" log M)$ using size-$k_"ms"$ min-heaps. By contrast, if $C$ is a dense cluster in subspace $S$ with $abs(S) = s$, monotonicity forces every non-empty subset $S' subset.eq S$ to be dense, so any breadth-first Apriori enumeration must generate and verify all $2^s - 1 = Omega(2^s)$ projections before reaching $S$.
]

#subsec("C", "Hypergeometric Chernoff Calibration of False-Merge Edges")

#theorem(title: "Family-Wise False-Merge Probability Bound under Independent Noise", origin: "Original Derivation")[
  Let $B_1 in cal(B)_(k_1)$ and $B_2 in cal(B)_(k_2)$ ($k_1 != k_2$) be independent noise base clusters of sizes $n_1 = abs(B_1)$ and $n_2 = abs(B_2)$ drawn uniformly without replacement from $cal(D)$, with expected overlap $mu = n_1 n_2 / n$ and $m = min(n_1, n_2)$. Then for every $t > 0$ and $delta > 0$:
  $ "Pr"(abs(B_1 inter B_2) >= mu + t) <= exp( - 2 t^2 / m ), $
  $ "Pr"(abs(B_1 inter B_2) >= (1 + delta) mu) <= exp( - (delta^2 mu) / (2 + delta) ). $
  Consequently, across all at most $binom(M, 2) <= M^2 / 2$ base-cluster pairs, enforcing the adaptive overlap threshold
  $ tau(B_1, B_2) = (n_1 n_2) / n + sqrt( (min(n_1, n_2) / 2) log ( M^2 / (2 alpha) ) ) $
  bounds the probability of *any* spurious noise edge in $G_"BC"$ below $alpha in (0, 1)$.
]

#proof[
  Conditioning on the larger set of size $K = max(n_1, n_2)$, the intersection size $X = abs(B_1 inter B_2)$ follows the hypergeometric distribution $"Hypergeom"(n, K, m)$ with mean $mu = m K / n = n_1 n_2 / n$. By Hoeffding's convex-order theorem for sampling without replacement, $X$ is stochastically dominated in the convex order by a binomial random variable $Y ~ "Bin"(m, K / n)$, meaning $bb(E)[phi(X)] <= bb(E)[phi(Y)]$ for every convex function $phi$; in particular, for all $lambda > 0$, $bb(E)[e^(lambda X)] <= bb(E)[e^(lambda Y)]$. Applying the Chernoff-Hoeffding bound to $Y = sum_(r=1)^m Z_r$ ($Z_r in {0, 1}$ i.i.d.) yields immediately:
  $ "Pr"(X >= mu + t) <= inf_(lambda > 0) e^(-lambda (mu + t)) bb(E)[e^(lambda Y)] <= exp( - 2 t^2 / m ). $
  Similarly, the multiplicative binomial Chernoff bound gives
  $ "Pr"(X >= (1 + delta) mu) <= ( e^delta / (1 + delta)^(1 + delta) )^mu = exp( -mu [ (1 + delta) log(1 + delta) - delta ] ). $
  Define $g(delta) = (1 + delta) log(1 + delta) - delta - delta^2 / (2 + delta)$ on $delta >= 0$. We have $g(0) = 0$, $g'(0) = 0$, and
  $ g''(delta) = 1 / (1 + delta) - 8 / (2 + delta)^3 = (delta^3 + 6 delta^2 + 4 delta) / ((1 + delta) (2 + delta)^3) >= 0, $
  so $g(delta) >= 0$ for all $delta >= 0$, establishing the second tail bound. Setting $exp(-2 t^2 / m) = 2 alpha / M^2$ gives $t = sqrt((m / 2) log(M^2 / (2 alpha)))$, and a union bound over all at most $M^2 / 2$ unordered pairs guarantees $"Pr"(exists "spurious edge") <= alpha$.
]

#subsec("D", "Constructive Failure under 1D Projection Masking")

Theorem II.A guarantees that every true cluster is *contained* in 1D base clusters, yet FIRES can still drop true subspace dimensions in Phase 2 when multiple clusters collapse onto a shared 1D projection.

#proposition(title: "Deterministic Pruning Loss under Co-Linear 1D Projection Masking", origin: "Original Counterexample")[
  For any integer $r >= 2$, any symmetric Jaccard threshold $tau_J > 1 / r$, and any top-$k$ pruning parameter $k_"ms" < r$, there exists a dataset $cal(D) subset bb(R)^2$ containing $r$ pairwise-disjoint $(epsilon, "MinPts")$-DBSCAN clusters $C_1, dots, C_r$ of equal size $m >= "MinPts"$ in full space $S = {1, 2}$ such that:
  (i) Unsplit Phase 1 merges all $r$ clusters into a single base cluster $B_1 = union_(q=1)^r C_q$ on coordinate $1$, reducing every Jaccard similarity $J(B_1, B_(2, q)) = 1 / r < tau_J$ and severing dimension $1$ from every cluster; and
  (ii) Under raw intersection ranking $s(B_1, B_(2, q)) = m$ with top-$k$ pruning $k_"ms" < r$, at least $r - k_"ms"$ true 2D clusters lose their edge to $B_1$ and cannot be recovered in Phase 2.
]

#proof[
  Choose $m >= "MinPts"$ and place $C_q = {(j / (m - 1), 10 q epsilon) : j = 0, dots, m - 1}$ for $q = 1, dots, r$, with $1 / (m - 1) <= epsilon$. Because consecutive points inside $C_q$ are separated by $1 / (m - 1) <= epsilon$ while distinct clusters $C_q, C_(q')$ are separated by at least $10 epsilon > epsilon$, the exact 2D DBSCAN clusters are $C_1, dots, C_r$. Projecting onto coordinate $2$ yields $r$ disjoint 1D base clusters $B_(2, q) = C_q$ of size $m$. Projecting onto coordinate $1$ collapses all $r$ sets onto the identical interval $[0, 1]$, forming one merged base cluster $B_1 = union_(q=1)^r C_q$ of size $r m$. Hence for every $q in {1, dots, r}$:
  $ J(B_1, B_(2, q)) = abs(B_1 inter B_(2, q)) / abs(B_1 union B_(2, q)) = m / (r m) = 1 / r < tau_J. $
  Thus every symmetric Jaccard threshold $tau_J > 1 / r$ eliminates all edges incident to $B_1$, leaving $B_1$ isolated in $G_"BC"$. Even if asymmetric containment $abs(B_1 inter B_(2, q)) / min(abs(B_1), abs(B_(2, q))) = 1$ or raw cardinality $s(B_1, B_(2, q)) = m$ is used, the $k_"ms"$-most-similar set $"MS"_(k_"ms") (B_1)$ has cardinality at most $k_"ms" < r$, so at least $r - k_"ms"$ base clusters $B_(2, q)$ fail the mutual membership test $B_(2, q) in "BM"(B_1)$ and are severed from coordinate $1$.
]

This counterexample provides a rigorous justification for FIRES's optional base-cluster splitting preprocessing and clarifies why $k_"ms"$ must exceed the maximum expected projection multiplicity $r$.

#sec("III", "Complexity, Geometry and Operational Profile")

By Theorem II.B, FIRES decouples subspace dimensionality from exponential enumeration: Phase 1 runs in $O(d n log n)$ using sorted 1D scans, and Phase 2 builds the entire exact base-cluster overlap graph in $O(n d^2)$ time using inverted object lists, followed by $O(M k_"ms" log M)$ pruning and graph clustering. Phase 3 executes the multidimensional refinement algorithm $cal(A)_"ref"$ (e.g., DBSCAN) on $m$ surviving candidate subspaces of average dimension $d_"avg"$, requiring $O(m n^2 d_"avg")$ in worst-case full-dataset refinement (or $O(sum_(r=1)^m abs(hat(C)_r)^2 abs(S_r))$ when localized to candidate supports). Space complexity is $O(n d + M + abs(E))$, avoiding any exponential lattice table.

Geometrically, FIRES is *modular*: its final cluster shapes are inherited directly from the refinement clusterer $cal(A)_"ref"$. When $cal(A)_"ref"$ is DBSCAN, FIRES recovers arbitrary non-convex, elongated, or crescent-shaped clusters within each axis-aligned subspace $S_r$, unlike centroid-based projected methods such as PROCLUS @aggarwal1999proclus or Gaussian-refinement methods like P3C @moise2006p3c.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$O(d n log n + n d^2)$ filter + $O(m n^2 d_"avg")$ refinement],
    [Space Complexity], [$O(n d + M + abs(E))$, strictly polynomial in $d$],
    [Core Mechanism], [1D base clusters $->$ mutual overlap graph $->$ subspace refinement],
    [Invariances], [Coordinate permutation invariant; scale invariant if $cal(A)_1, cal(A)_"ref"$ standardize],
    [Key Hyperparameters], [1D/refinement $(epsilon, "MinPts")$, $k_"ms"$, overlap threshold $tau_"sim"$],
    [Cluster Geometry], [Arbitrary non-convex shapes inside axis-aligned subspaces],
  ),
  caption: [Operational and geometric profile of FIRES (Alg. 146).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *High-Dimensional Functional Genomics and Proteomics:* Expression datasets with hundreds to thousands of features overwhelm exponential Apriori methods (SUBCLU, CLIQUE). FIRES's $O(n d^2)$ filter phase directly merges co-regulated 1D expression clusters into high-dimensional pathways before running exact density refinement @kriegel2005fires.
2. *Multi-Channel Industrial IoT and Sensor-Array Diagnostics:* Modern turbine, battery, and fab telemetry streams produce 50–500 continuous sensor channels where distinct operating regimes or failure modes exhibit non-convex manifold behavior in 5–20 dimensional subspaces. FIRES pairs fast 1D channel segmentation with arbitrary-shape DBSCAN refinement.
3. *Quantitative Multi-Factor Equity and Portfolio Analytics:* Asset universes characterized by dozens of valuation, momentum, liquidity, and volatility factors form regime clusters in changing factor subsets. Mutual best-merge graph clustering identifies coherent multi-factor regimes without requiring a pre-specified cluster count $K$.
4. *High-Dimensional Customer Behavioral Analytics:* E-commerce platforms track dozens of continuous behavioral metrics (recency, category spend, session cadence). FIRES handles heterogeneous background noise dimensions cleanly via Theorem II.C's overlap thresholding while isolating behavioral cohorts in their active feature subsets.
5. *Astronomical Multi-Band Photometric Catalogs:* Sky surveys record fluxes, colors, and morphological indices across dozens of bands where stellar populations and galaxy classes form curved loci in specific color subspaces; density refinement recovers these non-convex structures after 1D base filtering.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Four structural failure modes govern FIRES in practice: (1) *Heavy 1D projection overlap* (Proposition II.D), where $r > k_"ms"$ clusters stack along a coordinate and mask one another unless splitting is invoked; (2) *Global density-scale mismatch*, where a single radius $epsilon$ used for both 1D base discovery and $s$-dimensional DBSCAN refinement under-clusters dense high-$s$ subspaces due to $L_p$ distance inflation ($sqrt(s)$ scaling); (3) *Oblique/correlated subspaces*, where clusters tilted across axes have diffuse 1D marginals lacking clean base clusters; and (4) *Over-pruned sparse clusters*, where small genuine clusters fall below fixed overlap threshold $tau_"sim"$.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [High-$d$ Omics ($d > 50$)], [High: polynomial $O(n d^2)$ bypasses Apriori blowup], [Scale $epsilon$ with $s$; INSCY @assent2008inscy],
    [Non-Convex Subspaces], [High: DBSCAN refinement captures curved clusters], [Watch 1D projection stacking (Prop. II.D)],
    [Ultra-Fast Tabular], [Moderate: Phase-3 refinement costs $O(m n^2 d_"avg")$], [Use P3C @moise2006p3c (sub-ms)],
    [Tilted / Correlated Axes], [Low: 1D marginals smear oblique hyperplanes], [Use ORCLUS @aggarwal2000orclus or SSC @elhamifar2013ssc],
    [Strict Redundancy Control], [Moderate: graph merger can emit overlapping groups], [Use INSCY @assent2008inscy for $R$-pruning],
  ),
  caption: [Domain suitability and failure-mode matrix for FIRES.],
)

#sec("V", "Structural Relations and Empirical Behavior")

FIRES acts as a polynomial-time filter-refinement relaxation of exact bottom-up density subspace clustering (SUBCLU @kailing2004subclu), trading exhaustive lattice completeness for $O(n d^2)$ scalability.

#proposition(title: "Exact Recovery Condition Relative to SUBCLU", origin: "Original Derivation")[
  Let $cal(S)_"SUBCLU"$ be the set of maximal $(epsilon, "MinPts")$-density-connected subspace clusters found by SUBCLU @kailing2004subclu, and suppose every $(C, S) in cal(S)_"SUBCLU"$ satisfies: (a) its 1D projections $B_k supset.eq C$ ($k in S$) are unmasked ($B_k = C$); (b) $abs(C) >= tau_"sim"$; (c) $abs(S) - 1 <= k_"ms"$; and (d) distinct clusters $(C, S), (C', S')$ share strictly fewer than $tau_"sim"$ objects. Then the connected components of $G_"BC"$ coincide bijectively with the maximal subspaces of $cal(S)_"SUBCLU"$, and FIRES Phase 3 recovers $cal(S)_"SUBCLU"$ exactly without Apriori enumeration.
]

#proof[
  Under condition (a), each cluster $(C, S) in cal(S)_"SUBCLU"$ generates an exact family ${B_k : k in S}$ with $B_k = C$, so for every distinct pair $k_1, k_2 in S$, $s(B_(k_1), B_(k_2)) = abs(C) >= tau_"sim"$ by (b). By condition (d), any base cluster $B'$ belonging to a different cluster $(C', S')$ has $s(B_k, B') < tau_"sim"$, so $B' in.not "BM"(B_k)$. Since $B_k$ has exactly $abs(S) - 1 <= k_"ms"$ neighbors within $S$ (condition (c)), $"BM"(B_k) = {B_j : j in S without {k}}$ for all $k in S$. Therefore ${B_k : k in S}$ forms an isolated clique of size $abs(S)$ in $G_"BC"$, yielding candidate subspace $S$ exactly; running $(epsilon, "MinPts")$-DBSCAN on $pi_S (cal(D))$ in Phase 3 reproduces $C$ identically.
]

Compared with CLIQUE @agrawal1998clique (grid-based Apriori) and INSCY (Alg. 147) @assent2008inscy (depth-first SCY-tree with expected-density scaling), FIRES avoids fixed grid discretization and exponential lattice recursion, while differing from COSA (Alg. 145) @friedman2004cosa by identifying crisp discrete subspaces $S_r$ rather than continuous soft weights.

*Empirical Behavior on the ClusterOverview 4-Dataset Benchmark:* On the ClusterOverview evaluation suite (Silhouette @rousseeuw1987silhouettes and ARI @hubert1985comparing), FIRES records *Iris* $"Sil" = 0.419$, $"ARI" = 0.745$; *Customers* $"Sil" = 0.311$; *Multi-Shape* $"Sil" = 0.609$ ($"ARI" = "N/A"$); and *High-D* $"Sil" = 0.259$ ($"ARI" = "N/A"$), at a mean runtime of *15.611 ms*. Notably, FIRES attains the *highest Iris ARI* ($0.745$) among Chapters 145–148—outperforming COSA ($0.727$), P3C ($0.587$), and INSCY ($0.573$)—because merging clean 1D petal-feature base clusters isolates the exact discriminative subspace of *Iris* without distorting cluster boundaries. On *Multi-Shape*, its density-based refinement ties COSA at the top Silhouette ($0.609$), significantly ahead of P3C's elliptical EM refinement ($0.557$). Conversely, on *High-D* ($"Sil" = 0.259$), fixed-scale 1D-to-$s$D density thresholds lose contrast relative to INSCY's dimensionality-normalized expected density ($0.336$), while three-phase graph assembly plus refinement yields the highest mean runtime (15.611 ms) in this group.
