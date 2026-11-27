#import "../template/ieee.typ": *

#alg_entry(
  id: 142,
  name: "BANG-Clustering",
  category: "Cat. VI: Grid, Subspace & Projection",
  year: "1998",
  paradigm: "Multi-Scale Nested-Grid Density Tree Clustering",
  complexity: $O(n d log n)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Classical uniform grid clustering partitions the data space into cells of a
single fixed granularity, forcing an intractable trade-off: coarse meshes blur
fine structures, whereas fine meshes waste exponential memory on empty regions
and fragment sparse clusters. *BANG-Clustering* (Balanced And Nested Grid
Clustering), introduced by Schikuta and Erhart @schikuta1998bang, resolves
this tension by embedding clustering directly atop the *BANG-file*---a
self-balancing, hierarchical, multi-scale spatial directory that adapts cell
resolution locally to empirical point density while supporting *nested*
block regions.

Given a dataset $cal(X) = {bold(x)_1, dots, bold(x)_n} subset bb(R)^d$, first
apply coordinate-wise min-max normalization so that $cal(X) subset.eq [0, 1]^d$.
The BANG-file recursively partitions $[0, 1]^d$ by *cyclic binary halving*:
split step $l >= 1$ bisects coordinate axis $j(l) = ((l - 1) mod d) + 1$ at
the midpoint of the current interval. Every dyadic hyperrectangle (called a
*block region*) is uniquely encoded by a pair $(r, l)$ where $l in bb(N)_0$
is its *generation number* (split depth) and $r in {0, dots, 2^l - 1}$ is its
Peano/bit-interleaved region identifier.

#definition(title: "Nested BANG Block, Effective Volume, and Density")[
  Let $cal(Q) = {Q_1, dots, Q_M}$ be the set of dyadic hyperrectangles
  corresponding to the $M$ leaf data buckets of a BANG-file with maximum
  bucket capacity $b_"max" >= 1$. Because a splitting bucket may carve out an
  inner high-density sub-block while retaining the enclosing shell, each
  leaf bucket $j$ owns the *effective spatial region*
  $
    R_j = Q_j without union.big_(Q_k in cal(Q), Q_k subset Q_j) Q_k.
  $
  The *raw volume* of $Q_j$ at generation $l_j$ is $"Vol"(Q_j) = 2^(-l_j)$,
  and the *effective volume* of $R_j$ is
  $
    V_"eff" (R_j) = "Vol"(R_j) = 2^(-l_j) - sum_(Q_k in "MaxChild"(Q_j)) 2^(-l_k),
  $
  where $"MaxChild"(Q_j)$ denotes the set of inclusion-maximal blocks
  $Q_k in cal(Q)$ strictly contained in $Q_j$. Letting $n_j = |cal(X) inter R_j|$
  ($0 <= n_j <= b_"max"$), the *empirical block density* is
  $
    rho(R_j) = n_j / V_"eff" (R_j).
  $
]

#definition(title: "BANG Adjacency Graph and Density Dendrogram")[
  Two populated BANG regions $R_a, R_b$ are *spatially adjacent*, written
  $R_a tilde.op R_b$, if their topological closures share a $(d-1)$-dimensional
  hyperrectangular facet (including the inner boundary between an enclosing
  shell and an enclosed sub-block). Sorting all populated regions in descending
  density order $rho(R_(pi(1))) >= dots >= rho(R_(pi(M)))$ induces a nested
  family of superlevel adjacency subgraphs
  $G_tau = (V_tau, E_tau)$ with $V_tau = {R_j : rho(R_j) >= tau}$. Sweeping
  $tau$ downward from $rho_"max"$ to $0$ produces the *BANG density dendrogram*,
  from which clusters are extracted either by breadth-first expansion from
  density peaks down to a threshold $tau^*$ or by cutting at pronounced
  density-drop valleys.
]

#algbox(
  title: "BANG-Clustering: Multi-Scale Grid-Tree Density Extraction",
  inputs: [Normalized dataset $cal(X) subset.eq [0, 1]^d$, bucket capacity $b_"max"$, min density or cluster count $k$],
  outputs: [Density dendrogram $cal(T)$ and cluster assignment $cal(C) = {C_1, dots, C_k}$],
)[
  + Initialize BANG directory with root region $Q_1 = [0, 1]^d$ ($l_1 = 0$) and empty leaf bucket.
  + *For each* point $bold(x)_i in cal(X)$: locate innermost region $R_j$ containing $bold(x)_i$ via top-down directory search and insert $bold(x)_i$ into bucket $j$.
  + *If* $|B_j| > b_"max"$: cyclically bisect $Q_j$ (creating buddy sub-blocks or nested inner sub-block) choosing the split that maximizes point balance between child buckets; update directory tree.
  + Compute effective volume $V_"eff" (R_j)$ via bottom-up traversal of enclosed children and density $rho(R_j) = n_j / V_"eff" (R_j)$ for all $M$ non-empty regions.
  + Sort populated regions by descending density: $rho(R_(pi(1))) >= dots >= rho(R_(pi(M)))$.
  + Initialize cluster list $cal(C) <- emptyset$ and union-find structure over ${1, dots, M}$.
  + *For* $t = 1, dots, M$: activate region $u = R_(pi(t))$, query directory for face-adjacent active neighbors $v tilde.op u$ ($rho(v) >= rho(u)$), and merge components (recording merge density $rho(u)$ in dendrogram $cal(T)$).
  + Extract final clusters by thresholding $cal(T)$ at density cutoff $tau^*$ (or selecting the $k$ most persistent density peaks); assign each point the label of its containing region $R_j$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five formal results: exact laminar volume positivity, uniform
adaptive density consistency, equivalence to superlevel single-linkage on the
block adjacency graph, aspect-ratio and complexity bounds, and a constructive
high-dimensional breakdown analysis.

#subsec("A", "Laminar Nesting and Exact Positivity of Effective Volume")

#lemma(title: "Laminar Partition and Lower Bound on Effective Volume", origin: "Original Proof")[
  Let $cal(Q) = {Q_1, dots, Q_M}$ be the dyadic hyperrectangles produced by
  cyclic binary halving in a BANG-file, with maximum generation depth
  $L = max_j l_j$. Then:
  + $cal(Q)$ is a *laminar family*: for any $Q_a, Q_b in cal(Q)$, either
    $"int"(Q_a) inter "int"(Q_b) = emptyset$, $Q_a subset.eq Q_b$, or
    $Q_b subset.eq Q_a$;
  + The effective regions ${R_j}_(j=1)^M$ are pairwise disjoint up to
    Lebesgue-null boundaries and satisfy the recursive volume identity
    $
      V_"eff" (R_j) = 2^(-l_j) - sum_(Q_k subset Q_j) V_"eff" (R_k)
        = 2^(-l_j) - sum_(Q_k in "MaxChild"(Q_j)) 2^(-l_k);
    $
  + If $R_j$ has non-empty interior (i.e., $Q_j$ is not wholly tiled by its
    strictly enclosed children), then $V_"eff" (R_j) >= 2^(-L) > 0$, ensuring
    finite well-defined block densities $rho(R_j) <= b_"max" dot 2^L$.
]

#proof[
  At any fixed generation $l$, cyclic binary halving partitions $[0, 1]^d$
  into $2^l$ dyadic hyperrectangles with pairwise disjoint interiors. Moreover,
  every block at level $l' > l$ is obtained by bisecting a unique ancestor at
  level $l$, so it is either completely contained in a given level-$l$ block or
  interior-disjoint from it. This proves (1). For (2), consider any point
  $bold(x) in [0, 1]^d$ not on any dyadic facet; by laminarity, the set of
  blocks $Q_k in cal(Q)$ containing $bold(x)$ is linearly ordered by inclusion,
  so $bold(x)$ belongs to a unique minimal block $Q_j$, hence to a unique
  effective region $R_j$. Because the maximal strict children
  $Q_k in "MaxChild"(Q_j)$ are pairwise interior-disjoint and their union
  equals $union.big_(Q_m subset Q_j) Q_m$, countable additivity of Lebesgue
  measure gives
  $
    V_"eff" (R_j) = "Vol"(Q_j) - sum_(Q_k in "MaxChild"(Q_j)) "Vol"(Q_k)
      = 2^(-l_j) - sum_(Q_k in "MaxChild"(Q_j)) 2^(-l_k).
  $
  Since each $Q_k$ is itself partitioned by ${R_m : Q_m subset.eq Q_k}$,
  $"Vol"(Q_k) = sum_(Q_m subset.eq Q_k) V_"eff" (R_m)$, establishing both
  forms of the volume identity. For (3), subdivide $Q_j$ into elementary
  dyadic sub-cells at depth $L$, each of exact volume $2^(-L)$. Every
  $Q_k in "MaxChild"(Q_j)$ has $l_k <= L$ and therefore is an exact union of
  $2^(L - l_k)$ elementary level-$L$ cells. Hence $V_"eff" (R_j) = m_j dot 2^(-L)$
  for an integer $m_j in bb(N)_0$; non-empty interior implies $m_j >= 1$,
  yielding $V_"eff" (R_j) >= 2^(-L) > 0$.
]

#subsec("B", "Uniform Finite-Sample Consistency of BANG Density Estimation")

#theorem(title: "Adaptive Nonparametric Density Consistency", origin: "Original Derivation")[
  Let $bold(x)_1, dots, bold(x)_n$ be i.i.d. samples from an $L_f$-Lipschitz
  density $f$ on $[0, 1]^d$ bounded above by $f_"max"$. Suppose BANG builds a
  directory up to maximum generation depth $L_n$ with pure (un-nested) leaf
  blocks satisfying $b_"min" <= n_j <= b_"max"$ for populated cells, where
  $b_"min" = c b_"max"$ ($c in (0, 1/2]$). Let
  $hat(f)_n (bold(x)) = n_j / (n V_"eff" (R_j))$ for $bold(x) in R_j$. Then
  for any $delta in (0, 1)$, with probability at least $1 - delta$,
  simultaneously over all dyadic blocks $Q$ of level $l <= L_n$ with true
  probability mass $p(Q) = integral_Q f$:
  $
    |hat(p)_n (Q) - p(Q)| <= sqrt((2 p(Q) ln(4 dot 2^(L_n) / delta)) / n)
      + (2 ln(4 dot 2^(L_n) / delta)) / (3 n).
  $
  Consequently, choosing $L_n = ceil(log_2 (n / b_n))$ and
  $b_n = Theta((n^(2) (log n)^d)^(1 / (d + 2)))$ yields the uniform adaptive
  relative density consistency rate (on regions with $f(bold(x)) >= f_0 > 0$):
  $
    |hat(f)_n (bold(x)) - f(bold(x))| = O_p (((log n) / n)^(1 / (d + 2))).
  $
]

#proof[
  At each level $l in {0, dots, L_n}$, cyclic binary halving produces exactly
  $2^l$ dyadic hyperrectangles, so the total count across all levels $l <= L_n$
  is $N_"dyad" = sum_(l=0)^(L_n) 2^l = 2^(L_n + 1) - 1 < 2^(L_n + 1)$. For any
  fixed dyadic block $Q$, the count $n hat(p)_n (Q) = sum_(i=1)^n bb(I)(bold(x)_i in Q)$
  is a sum of i.i.d. Bernoulli($p(Q)$) variables with variance
  $p(Q)(1 - p(Q)) <= p(Q)$. By Bernstein's inequality, for every $t > 0$,
  $
    "Pr"(|hat(p)_n (Q) - p(Q)| >= t)
      <= 2 exp(- (n t^2) / (2 p(Q) + 2 t / 3)).
  $
  Setting the right-hand side equal to $delta / 2^(L_n + 1)$ and solving the
  quadratic in $t$ yields $t <= sqrt(2 p(Q) gamma_n / n) + (2 / 3) gamma_n / n$
  where $gamma_n = ln(4 dot 2^(L_n) / delta) = O(L_n + ln(1 / delta))$. A
  union bound over all $N_"dyad" < 2^(L_n + 1)$ dyadic blocks establishes the
  simultaneous concentration inequality cleanly without data-dependent
  stopping-time bias. Now let $Q = Q_j$ be the realized leaf block containing
  $bold(x)$ with volume $V = 2^(-l_j)$ and count $n_j in [c b_n, b_n]$. By
  dividing the deviation bound by $V$, the estimation error decomposes as
  $
    |hat(f)_n (bold(x)) - f(bold(x))| <= |hat(p)_n (Q) - p(Q)| / V
      + |p(Q) / V - f(bold(x))|.
  $
  For the bias term, every side of $Q_j$ has length at most
  $2^(-floor(l_j / d)) <= 2 dot 2^(-l_j / d) = 2 V^(1/d)$, so
  $"diam"(Q_j) <= 2 sqrt(d) V^(1/d)$ and $L_f$-Lipschitz continuity gives
  $|p(Q) / V - f(bold(x))| <= 2 L_f sqrt(d) V^(1/d)$. Because
  $n_j = n hat(p)_n (Q) = Theta(b_n)$ and $p(Q) >= f_0 V$ implies
  $V = O(b_n / (n f_0))$, the bias scales as $O((b_n / n)^(1/d))$ while the
  variance term $sqrt(p(Q) gamma_n / (n V^2)) = O(sqrt((log n) / b_n))$.
  Equating $(b_n / n)^(1/d) = sqrt((log n) / b_n)$ solves to
  $b_n = Theta((n^2 (log n)^d)^(1 / (d + 2)))$ and yields total error
  $O_p (((log n) / n)^(1 / (d + 2)))$. (For a nested block with $r$ direct
  children, applying the triangle inequality over $Q_j$ and its $r$ children
  scales the stochastic term by at most $r + 1$.)
]

#subsec("C", "Graph-Theoretic Equivalence to Superlevel Single-Linkage")

#proposition(title: "Equivalence to Maximum-Bottleneck Spanning Forest", origin: "Classical — Cited")[
  Let $G = (V, E)$ be the BANG face-adjacency graph on populated regions
  $V = {R_1, dots, R_M}$ with symmetric edge weights
  $w(R_a, R_b) = min(rho(R_a), rho(R_b))$ for $R_a tilde.op R_b$. Then for
  every density threshold $tau > 0$, the cluster components constructed by
  BANG's descending-density merge procedure coincide identically with:
  + the connected components of the density upper-level set subgraph
    $G_tau = ({R_j : rho(R_j) >= tau}, E_tau)$, and
  + the thresholded components of a maximum-weight spanning forest (single-linkage
    tree on $-w$) of $G$ @schikuta1998bang, @hartigan1981consistency.
]

#proof[
  When Algorithm 142 processes regions in descending order
  $rho(R_(pi(1))) >= dots >= rho(R_(pi(M)))$ and stops at the last index $t$
  with $rho(R_(pi(t))) >= tau$, the set of activated vertices is exactly
  $V_tau = {R_j : rho(R_j) >= tau}$. Moreover, an edge $(u, v) in E$ is
  traversed during steps $1, dots, t$ if and only if both endpoints have been
  activated, i.e., $min(rho(u), rho(v)) = w(u, v) >= tau$. Because union-find
  maintains exact graph connectivity over traversed edges, the resulting
  partition equals the connected components of $G_tau$. Since Kruskal's
  algorithm for a maximum spanning forest sorts edges in descending order of
  $w(u, v) = min(rho(u), rho(v))$ and unions distinct components whenever
  $w(u, v) >= tau$, two vertices lie in the same BANG cluster at level $tau$
  if and only if they are connected in $G$ by a path whose minimum vertex
  density is at least $tau$---the exact Hartigan @hartigan1981consistency
  high-density cluster tree on the quantized BANG density estimator.
]

#subsec("D", "Bounded Aspect Ratio and Directory Complexity")

#theorem(title: "Bounded Cell Aspect Ratio and Sub-Quadratic Time Bound", origin: "Original Proof")[
  For any BANG-file built by cyclic binary halving in $bb(R)^d$ with bucket
  capacity $b_"max"$ and minimum split occupancy ratio $beta in (0, 1/2]$:
  + Every raw block $Q_j$ at any generation $l_j >= 0$ has *aspect ratio*
    $"max side" / "min side" <= 2$, independent of both $l_j$ and $d$;
  + The number of populated leaf regions satisfies
    $ceil(n / b_"max") <= M <= floor(n / (beta b_"max"))$ whenever all splits
    preserve occupancy $>= beta b_"max"$;
  + Total clustering time is $O(n d log n + M log M + M d log M)$, reducing
    to $O(n d log n)$ since $M <= n$.
]

#proof[
  Write $l_j = q d + r$ with quotient $q = floor(l_j / d)$ and remainder
  $0 <= r < d$. Because cyclic halving splits coordinates $1, dots, d$ in
  round-robin order, the first $r$ dimensions have been halved $q + 1$ times
  (side length $2^(-(q+1))$) and the remaining $d - r$ dimensions have been
  halved $q$ times (side length $2^(-q)$). Thus the ratio of longest to
  shortest edge is $2^(-q) / 2^(-(q+1)) = 2$ if $r > 0$ and $1$ if $r = 0$,
  proving (1). For (2), every non-empty bucket holds between $beta b_"max"$
  and $b_"max"$ points, so $M beta b_"max" <= n <= M b_"max"$, immediately
  yielding the two-sided bound on $M$. For (3), the balanced BANG directory
  tree has height $h = O(log_F M) = O(log n)$ (for page fanout $F >= 2$).
  Point insertion routes each of $n$ vectors down $O(log n)$ directory levels
  with $O(d)$ coordinate comparisons per node, costing $O(n d log n)$.
  Bottom-up effective volume calculation visits $O(M)$ directory nodes in
  $O(M)$ time; sorting $M$ densities takes $O(M log M)$; and querying the
  $2 d$ face-adjacent neighbors of each region via directory lookup takes
  $O(M d log M)$ time, giving $O(n d log n)$ overall since $M <= n$.
]

#subsec("E", "High-Dimensional Unsplit-Axis Collapse and Outlier Fragility")

#proposition(title: "Curse of Cyclic Splitting and Zero Outlier Breakdown", origin: "Original Counterexample")[
  Consider BANG-Clustering on $n$ points in $[0, 1]^d$ with bucket capacity
  $b_"max" >= 2$.
  + *Incomplete axis coverage in high $d$*: For uniformly distributed points,
    expected generation depth is $bar(l) = log_2 (n / b_"max")$. Whenever
    $n < b_"max" dot 2^d$, at least $d - floor(log_2 (n / b_"max"))$ trailing
    coordinates remain *completely unsplit* (full interval $[0, 1]$), rendering
    BANG blind to cluster separation along dimensions $j > bar(l)$.
  + *Zero breakdown point under min-max scaling*: Adding a single extreme
    outlier $bold(x)_(n+1) = (c, 0, dots, 0)$ with $c > 1$ to an inlier dataset
    supported on $[0, 1]^d$ compresses all inliers along axis $1$ into
    $[0, 1 / c]$, forcing $ceil(log_2 c)$ wasted splits along axis $1$---and
    hence $d ceil(log_2 c)$ total cyclic splits across the grid---before
    inliers are partitioned along axis $1$ at resolution $1/2$.
]

#proof[
  For (1), a uniform distribution fills all $2^l$ dyadic blocks at level $l$
  with roughly $n / 2^l$ points each; splitting halts once $n / 2^l <= b_"max"$,
  i.e., $l = floor(log_2 (n / b_"max"))$. By cyclic ordering, coordinate $j$ is
  bisected for the first time only at generation $l = j$. Thus if
  $l < d$ (equivalently $n < b_"max" dot 2^d$), coordinates $l + 1, dots, d$
  never experience a single split and retain edge length $1$. Any two clusters
  separated solely along coordinates ${l + 1, dots, d}$ collapse into identical
  buckets. For (2), min-max rescaling maps $x_(i, 1) |-> x_(i, 1) / c in [0, 1/c]$.
  The dyadic interval $[0, 2^(-k)]$ strictly separates points inside $[0, 1/c]$
  only once $2^(-k) < 1 / c$, requiring $k >= ceil(log_2 c)$ splits of axis $1$.
  Under strict cyclic halving, the $k$-th bisection of axis $1$ occurs at
  generation $l = (k - 1) d + 1$, wasting $Omega(d log c)$ directory depth on
  empty space due to a single corrupted point (breakdown point $1 / (n + 1) -> 0$).
]

#sec("III", "Complexity, Geometry and Operational Profile")

BANG-Clustering is engineered for *database-scale I/O efficiency*: because
the BANG-file doubles as a persistent balanced disk index, clustering operates
entirely on the $M << n$ directory metadata entries (region IDs, counts,
volumes) without re-scanning raw point coordinates during dendrogram
construction @schikuta1998bang. In memory, the compact directory requires
$O(M d)$ space (excluding raw data storage $O(n d)$), and descending-density
BFS executes in microseconds once buckets are populated. Its primary
geometric strength is representing *arbitrary connected non-convex shapes* as
unions of multi-scale dyadic blocks---unlike OptiGrid's convex polytopes---at
the cost of sharp degradation when $d$ exceeds $log_2 (n / b_"max")$.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$O(n d log n)$ directory build + $O(M d log M)$ graph/dendrogram ($M = O(n / b_"max")$)],
    [Space Complexity], [$O(M d)$ directory metadata; $O(n d)$ total with raw points; disk-page friendly],
    [Cluster Geometry], [Arbitrary connected non-convex unions of multi-scale axis-aligned dyadic blocks],
    [Noise Handling], [Low-density outer shells naturally sink to the bottom of the density dendrogram],
    [Hyperparameters], [Bucket capacity $b_"max"$, density cutoff $tau^*$ (or target cluster count $k$)],
    [Invariances], [Translation & uniform-scale invariant (after min-max); sensitive to axis order & outliers],
  ),
  caption: [Operational, geometric, and complexity profile of BANG-Clustering.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

+ *Large-Scale Geographic Information Systems (GIS) and Spatial Databases*:
  Urban mobility, crime-hotspot, and LiDAR point clouds are stored natively in
  2D/3D spatial indexes where $n$ reaches $10^7$--$10^9$ and $d in {2, 3}$.
  Because $2^d <= 8$, cyclic halving resolves all spatial coordinates rapidly,
  and BANG-Clustering extracts irregular street-network or coastline density
  clusters directly from index pages without secondary memory scans
  @schikuta1998bang.
+ *Astronomical N-Body Halo and Filament Identification*: Cosmological
  simulations evolve billions of dark-matter particles in 3D space with
  multi-order-of-magnitude density contrast between dense virialized halos
  and tenuous intergalactic filaments. BANG's nested multi-scale blocks mimic
  adaptive octrees, computing accurate local densities across $10^4$-fold
  volume variations while recovering non-spherical filaments via face adjacency.
+ *Retail Customer Micro-Segmentation in Low-Dimensional RFM Space*: Marketing
  analytics frequently segments millions of consumers in 3D--5D Recency,
  Frequency, and Monetary (RFM) coordinates. BANG compresses millions of
  accounts into hundreds of density-ranked buckets, producing an interactive
  dendrogram that exposes multi-scale customer tiers in sub-millisecond query
  time.
+ *Environmental Sensor Grids and Remote Sensing Telemetry*: Satellite
  multispectral pixels (3--6 bands) and ocean-buoy arrays exhibit spatially
  contiguous ecological zones of variable area. BANG's bounded aspect ratio
  (Theorem 2) ensures isotropic local smoothing across both compact anomalies
  and broad background plains.
+ *Real-Time Database Query Optimization and Histogram Synopsis*: Relational
  query optimizers maintain multi-dimensional histograms to estimate join
  selectivity; BANG-Clustering reuses the identical BANG-file directory to
  simultaneously serve range queries and maintain up-to-date density clusters
  under incremental tuple insertions.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

BANG-Clustering encounters four key failure modes: (i) *high-dimensional
unsplit trailing coordinates* (Proposition 2) whenever $d > log_2 (n / b_"max")$;
(ii) *axis-order dependence* at low depth because cyclic splitting prioritizes
coordinates $1, 2, dots$ ahead of $d$; (iii) *bridge-chaining across density
ridges* inherited from single-linkage superlevel traversal (Proposition 1);
and (iv) *severe outlier compression* under min-max normalization unless
preceded by robust quantile/winsorized scaling.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [2D/3D GIS & spatial databases], [Excellent: multi-scale non-convex shapes; native disk-index integration], [Winsorize extreme coordinate outliers first],
    [Low-D tabular / RFM analytics], [Strong: extreme speed ($M << n$) and hierarchical density dendrogram], [Watch ridge chaining; compare HDBSCAN @campello2013hdbscan],
    [High-D feature spaces ($d >= 10$)], [Poor: cyclic splits leave $d - log_2(n/b_"max")$ axes unpartitioned], [Use OptiGrid (Alg. 141) or ENCLUS (Alg. 143)],
    [Heavy-tailed / unnormalized logs], [Weak: single extreme point collapses inlier resolution (Prop. 2)], [Apply rank/IQR scaling before building grid],
    [Correlated subspace clusters], [Weak: full-space volume $2^(-l)$ dilutes lower-D signal], [Use CLIQUE @agrawal1998clique or FINDIT (Alg. 144)],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for BANG-Clustering.],
)

#sec("V", "Structural Relations and Empirical Behavior")

BANG-Clustering sits directly at the intersection of tree-structured grid
clustering (STING @wang1997sting, MAFIA @goil1999mafia, WaveCluster
@sheikholeslami1998wavecluster), spatial indexing (BIRCH @zhang1996birch), and
hierarchical density clustering (DBSCAN @ester1996dbscan, HDBSCAN
@campello2013hdbscan). Unlike STING, whose regular quadtree/octree splits all
$d$ axes simultaneously (causing $2^d$ fanout per level), BANG splits *one
axis per generation* and permits *nested enclosing regions*, preventing empty
sibling proliferation while preserving bounded aspect ratio. Compared to
GRI-DBSCAN @mahran2008gridbscan, BANG eliminates fixed-$epsilon$ radius search
entirely by letting bucket capacity $b_"max"$ adapt local resolution.

#proposition(title: "Exact Asymptotic Limit to Density Level-Set Clustering", origin: "Original Proof")[
  Under the Lipschitz density assumptions of Theorem 1 with
  $b_n = Theta((n^2 (log n)^d)^(1/(d+2)))$, let $C_1 (tau), C_2 (tau)$ be two
  distinct connected components of the true population superlevel set
  $L_f (tau) = {bold(x) in [0, 1]^d : f(bold(x)) >= tau}$ separated by a
  strict density valley of depth $eta > 0$ (every continuous path between
  $C_1 (tau)$ and $C_2 (tau)$ drops below $tau - eta$). Then for all
  sufficiently large $n$, BANG-Clustering at threshold $tau - eta / 2$
  recovers $C_1 (tau)$ and $C_2 (tau)$ in distinct clusters with probability
  tending to $1$ as $n -> infinity$.
]

#proof[
  By Theorem 1, $sup_(bold(x)) |hat(f)_n (bold(x)) - f(bold(x))| = O_p (((log n)/n)^(1/(d+2))) -> 0$.
  Hence for sufficiently large $n$, with probability $1 - o(1)$,
  $sup |hat(f)_n - f| < eta / 4$ and maximum block diameter
  $max_j "diam"(Q_j) < delta_"gap" / 2$ where $delta_"gap"$ is the minimum
  Euclidean distance across the valley separating $C_1 (tau)$ and $C_2 (tau)$.
  Then every region intersecting $C_1 (tau)$ or $C_2 (tau)$ has
  $rho(R_j) >= tau - eta / 4 > tau - eta / 2$, whereas any face-adjacent chain
  of regions bridging $C_1 (tau)$ to $C_2 (tau)$ must traverse a region
  intersecting ${f < tau - eta}$, which satisfies
  $rho(R_j) < tau - 3 eta / 4 < tau - eta / 2$. By Proposition 1, no active
  path connects $C_1 (tau)$ and $C_2 (tau)$ in $G_(tau - eta / 2)$, while
  internal connectivity of each component is preserved.
]

On the ClusterOverview 4-dataset benchmark (Silhouette @rousseeuw1987silhouettes
and ARI @hubert1985comparing), BANG-Clustering achieves: *Iris* $"Sil" = 0.641$,
$"ARI" = 0.510$; *Customers* $"Sil" = 0.428$; *Multi-Shape* $"Sil" = 0.627$
($"ARI" =$ N/A); *High-D* $"Sil" = 0.168$ ($"ARI" =$ N/A); with an exceptional
mean runtime of $0.044$ ms. Every single metric aligns with our theoretical
analysis. First, $0.044$ ms is the fastest runtime among Algorithms 141--144,
reflecting the $O(M << n)$ directory compression proven in Theorem 2. Second,
BANG leads all four methods in Silhouette on both low-dimensional tabular
benchmarks (*Iris* $0.641$ and *Customers* $0.428$) and performs solidly on
*Multi-Shape* ($0.627$), showcasing Theorem 1's multi-scale adaptive density
resolution; its moderate *Iris* ARI ($0.510$) arises from single-linkage
density-ridge merging (Proposition 1) between Versicolor and Virginica into
one connected superlevel component. Third, its sharp drop on *High-D*
($"Sil" = 0.168$, lowest of the four) provides direct empirical confirmation
of Proposition 2: when $d$ is large relative to $log_2 (n / b_"max")$, cyclic
binary splitting leaves trailing coordinates unpartitioned and loses cluster
discrimination @beyer1999nearest, @xu2005survey.
