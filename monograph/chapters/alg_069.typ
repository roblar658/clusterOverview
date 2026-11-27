#import "../template/ieee.typ": *

#alg_entry(
  id: 69,
  name: "Genie: Inequality-Constrained Agglomerative Clustering",
  category: "Cat. III: Hierarchical & Tree-Based",
  year: "2016",
  paradigm: "Gini-Constrained Single-Linkage Hierarchical Clustering",
  complexity: $O(n^2)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Single Linkage (*SL*) clustering @sibson1973slink @gower1969minimum
recovers arbitrarily shaped non-convex clusters and enjoys an exact
equivalence to Kruskal's Minimum Spanning Tree (*MST*) algorithm, yet
suffers from two notorious pathologies: _chaining_ across thin noise
bridges and _outlier isolation_, wherein cutting an SL dendrogram at
$K$ clusters frequently yields one giant cluster plus $K - 1$ singleton
outliers @carlsson2010characterization @mullner2011modern. Classical
remedies such as Ward's method @ward1963hierarchical or Complete Linkage
@defays1977clink prevent singleton dominance by biasing toward compact
spherical clusters, sacrificing the ability to track elongated manifolds.

*Genie*, introduced by Gagolewski, Bartoszuk, and Cena
@gagolewski2016genie, resolves this tension by augmenting Single Linkage
with an explicit *economic inequality constraint* on the cluster-size
distribution. Let $(X, d)$ be a finite metric or symmetric dissimilarity
space with $X = {x_1, dots, x_n}$ and $d (x_i, x_j) = d (x_j, x_i) >=
0$. At agglomerative step $t in {0, 1, dots, n - 2}$, let
$cal(C)^((t)) = {C_1, dots, C_k}$ denote the active partition of size
$k = n - t >= 2$, with positive integer cardinalities $c_r = abs(C_r)$
satisfying $sum_(r = 1)^k c_r = n$. Let $c_((1)) <= c_((2)) <= dots.c <=
c_((k))$ denote the sorted cluster sizes.

#definition(title: "Normalized Cluster-Size Gini Index and Genie Merge Rule")[
  For $k >= 2$ nonempty clusters with size vector $c = (c_1, dots, c_k)$,
  the *normalized Gini index* @gagolewski2016genie is
  $
    G (c) = (sum_(i = 1)^k sum_(j = 1)^k |c_i - c_j|) / (2 (k - 1) n)
      = (sum_(r = 1)^k (2 r - k - 1) c_((r))) / ((k - 1) n).
  $
  Given a user threshold $g in (0, 1]$ (typically $g = 0.3$) and the
  single-linkage distance $d_"SL" (C_i, C_j) = min_(x in C_i, y in C_j)
  d (x, y)$, let $c_"min" = min_(1 <= r <= k) c_r$. At step $t$, Genie
  merges the pair $(i^*, j^*)$ minimizing $d_"SL" (C_i, C_j)$ over:
  - all distinct pairs $1 <= i < j <= k$ if $G (c) <= g$ (*unconstrained
    SL step*);
  - pairs with $min(c_i, c_j) = c_"min"$ if $G (c) > g$ (*inequality-
    constrained step*).
]

Intuitively, whenever the partition remains reasonably balanced
($G (c) <= g$), Genie tracks arbitrary non-convex manifolds identically
to Single Linkage. The instant outlier isolation or unbalanced growth
pushes size inequality above $g$, Genie forbids merging two already
large clusters and forces a minimum-size cluster to merge with its
nearest neighbor—absorbing isolated noise points before they can
monopolize the top cuts of the dendrogram.

#algbox(
  title: "Genie (MST-Based Agglomerative Algorithm)",
  inputs: [Dataset $(X, d)$ of size $n$; Gini threshold $g in (0, 1]$;
    optional target $K$],
  outputs: [Merge sequence $(i_t, j_t, h_t)_(t = 1)^(n - 1)$ and
    partition $cal(C)^((n - K))$],
)[
  + *MST Construction:* Compute a Minimum Spanning Tree $T^* = (X, E^*)$
    of the complete dissimilarity graph $(X, d)$ (using Prim's algorithm
    in $O(n^2)$ or spatial $k$-d / Delaunay acceleration in low $d$).
  + *Initialization:* Set $cal(C)^((0)) <- {{x_1}, dots, {x_n}}$ with
    sizes $c_r = 1$; initialize the order-statistic size BST (with Gini
    numerator $S = 0$), per-cluster meldable edge heaps $H_r$ of
    incident edges in $T^*$, and a global size-distance BST
    $cal(B)$ keyed by $(c_r, w_"min" (C_r), r)$.
  + *For* $t = 0, 1, dots, n - 2$:
    - Let $k = n - t$ and evaluate $G = S / ((k - 1) n)$ in $O(1)$.
    - *If* $G <= g$: extract the global lightest remaining edge
      $e^* = (C_(i^*), C_(j^*))$ of $T_t$ of weight $h_(t + 1)$.
    - *Else* ($G > g$): query the minimum-key entry in $cal(B)$ to
      obtain $i^*$ with $c_(i^*) = c_"min"$ minimizing $w_"min" (C_i)$,
      and let $e^* = (C_(i^*), C_(j^*))$ be its top alive edge of weight
      $h_(t + 1)$.
    - Contract $e^*$ in $T_t$: merge $C_(i^*)$ and $C_(j^*)$ into
      $C_"new"$ of size $c_(i^*) + c_(j^*)$, mark $e^*$ dead, meld
      $H_(i^*)$ and $H_(j^*)$, lazily pop dead top edges, and update
      $S$ and $cal(B)$ in $O(log n)$.
  + *Return* the merge history (or cut at $k = K$ clusters).
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish four foundational results: the classical MST sufficiency
theorem of Gagolewski et al. @gagolewski2016genie, exact Gini bounds
coupled with a counterexample showing that gradient chaining persists
even at $g = 0$, a universal construction proving merge-height
inversions for every $g in (0, 1)$ alongside monotone distance
invariance, and an original $O(n log n)$ post-MST exact algorithm.

#subsec("A", "Minimum Spanning Tree Edge Sufficiency")

#theorem(
  title: "Exact Reduction of Genie to MST Edge Contractions",
  origin: "Classical — Cited",
)[
  Suppose all pairwise dissimilarities $d (x_u, x_v)$ ($u < v$) are
  distinct (or ties are broken by a fixed global edge order), and let
  $T^* = (X, E^*)$ be the unique Minimum Spanning Tree of $(X, d)$.
  Then for every threshold $g in (0, 1]$ and every step
  $t in {0, dots, n - 2}$, the inter-cluster pair $(C_(i^*), C_(j^*))$
  merged by Genie @gagolewski2016genie is witnessed by an edge
  $e^* = (u^*, v^*) in E^*$ of the MST.
]

#proof[
  We proceed by induction on $t$. At any step $t$, let
  $cal(C)^((t)) = {C_1, dots, C_k}$ be the current partition of $X$, and
  let $e^* = (u^*, v^*)$ realize the Genie merge weight
  $h_(t + 1) = d_"SL" (C_(i^*), C_(j^*))$ with $u^* in C_(i^*)$ and
  $v^* in C_(j^*)$ ($i^* != j^*$).

  In an unconstrained step ($G (c) <= g$), $e^*$ minimizes $d (u, v)$
  over all edges crossing *any* cluster boundary, so in particular $e^*$
  is the unique minimum-weight edge crossing the cut
  $(C_(i^*), X without C_(i^*))$. In a constrained step ($G (c) > g$),
  without loss of generality let $c_(i^*) = c_"min"$. Every edge
  $e' = (u', v')$ crossing the cut $(C_(i^*), X without C_(i^*))$
  connects $C_(i^*)$ to some cluster $C_r$ ($r != i^*$), and hence
  satisfies $min(c_(i^*), c_r) = c_"min"$. Therefore every cut edge of
  $(C_(i^*), X without C_(i^*))$ belongs to the feasible candidate set of
  Case 2, implying $d (e^*) <= d (e')$ for all edges crossing that cut.

  In both cases, $e^*$ is the unique minimum-weight edge crossing the
  nontrivial cut $(C_(i^*), X without C_(i^*))$. By the fundamental Cut
  Property of Minimum Spanning Trees, if $e^* in.not E^*$, then adding
  $e^*$ to $T^*$ creates a unique cycle containing another cut edge
  $e' != e^*$ crossing $(C_(i^*), X without C_(i^*))$ with
  $d (e') > d (e^*)$, so $(T^* without {e'}) union {e^*}$ is a strictly
  lighter spanning tree—a contradiction. Thus $e^* in E^*$ at every step.
]

#subsec("B", "Sharp Gini Extremes and Persistence of Gradient Chaining")

#lemma(
  title: "Sharp Discrete Gini Bounds and Unbounded Chaining at g = 0",
  origin: "Original Counterexample",
)[
  For any partition of $n >= 2$ objects into $k in {2, dots, n}$
  nonempty clusters:
  + $0 <= G (c) <= 1 - k / n$, where $G (c) = 0$ iff $c_1 = dots.c =
    c_k = n / k$, and $G (c) = 1 - k / n$ iff $c$ is a permutation of
    $(n - k + 1, 1, dots, 1)$. Consequently, Genie reduces identically
    to Single Linkage whenever $g >= 1 - 2 / n$.
  + *Chaining Persistence:* Even in the strictest limit $g = 0$, for
    every $n >= 3$ there exists a 1D Euclidean dataset on which Genie
    produces the exact maximally unbalanced Single-Linkage caterpillar
    $(t + 1, 1, dots, 1)$ at every step $t$, attaining the theoretical
    ceiling $G (c^((t))) = t / n -> 1$ as $t -> n - 2$.
]

#proof[
  (1) Non-negativity $G (c) >= 0$ with equality iff all $c_r$ coincide
  follows immediately from $sum_(i, j) |c_i - c_j| >= 0$. For the upper
  bound, write $c_((r)) = 1 + y_r$ with $0 <= y_1 <= dots.c <= y_k$ and
  $sum_(r = 1)^k y_r = n - k$. Since $sum_(r = 1)^k (2 r - k - 1) = 0$
  and $2 r - k - 1 <= k - 1$ for all $r <= k$,
  $
    G (c) = (sum_(r = 1)^k (2 r - k - 1) y_r) / ((k - 1) n)
      <= ((k - 1) sum_(r = 1)^k y_r) / ((k - 1) n) = 1 - k / n,
  $
  with equality iff $y_r = 0$ for all $r < k$ and $y_k = n - k$, i.e.,
  $c_((k)) = n - k + 1$ and all other $c_((r)) = 1$. Since $k >= 2$,
  $G (c) <= 1 - 2 / n$ at all steps, so $g >= 1 - 2 / n$ never triggers
  Case 2 and reproduces Single Linkage identically.

  (2) Choose strictly increasing gaps $0 < w_1 < w_2 < dots.c < w_(n - 1)$
  and define 1D Euclidean points $x_1 = 0$, $x_(r + 1) = x_r + w_r$ for
  $r = 1, dots, n - 1$. The unique MST is the path $x_1 - x_2 - dots.c -
  x_n$. We prove by induction on $t in {0, dots, n - 2}$ that Genie
  merges edge $(x_(t + 1), x_(t + 2))$ at step $t$ regardless of $g >= 0$.
  At $t = 0$, $G = 0 <= g$ and the lightest edge is $(x_1, x_2)$.
  Suppose after $t$ merges the clusters are
  $C_1 = {x_1, dots, x_(t + 1)}$ (size $t + 1$) and singletons
  ${x_(t + 2)}, dots, {x_n}$ (each size $1$). Then $c_"min" = 1$, and
  the globally lightest remaining edge is $(x_(t + 1), x_(t + 2))$ of
  weight $w_(t + 1)$. Because $min(abs(C_1), abs({x_(t + 2)})) =
  min(t + 1, 1) = 1 = c_"min"$, this edge satisfies Case 2 *even when*
  $G (c) > g = 0$. Hence Genie executes the exact Single-Linkage
  caterpillar sequence, hitting $G (c^((t))) = 1 - (n - t) / n = t / n$.
]

Lemma 1 clarifies Genie's exact mechanism: it does *not* force all
clusters to have equal size along monotone density gradients, but rather
blocks premature merges between two *non-minimal* clusters while
singletons remain unabsorbed.

#subsec("C", "Monotone Invariance and Universal Dendrogram Inversions")

#theorem(
  title: "Ordinal Invariance and Arbitrary-Ratio Merge-Height Inversions",
  origin: "Original Proof",
)[
  Let $g in (0, 1)$. Then:
  + *Ordinal Distance Invariance:* For any strictly increasing map
    $phi : bb(R)_(>= 0) -> bb(R)_(>= 0)$ with $phi (0) = 0$, Genie
    generates the identical partition sequence $cal(C)^((0)), dots,
    cal(C)^((n - 1))$ on $(X, phi compose d)$ as on $(X, d)$.
  + *Universal Non-Ultrametricity:* For *every* $g in (0, 1)$ and every
    ratio $eta > 0$, there exists a finite 1D Euclidean dataset $X$ on
    which Genie exhibits a strict merge-height inversion
    $h_(n - 1) / h_(n - 2) < eta$, violating ultrametricity.
]

#proof[
  (1) Because $G (c)$ and $c_"min"$ depend solely on cluster
  cardinalities and never on metric values, and every merge selection
  depends exclusively on comparisons $d (u, v) < d (u', v')$, applying
  any strictly increasing $phi$ preserves the total order of all edge
  weights, leaving every Case 1 and Case 2 decision identical.

  (2) Fix any $g in (0, 1)$ and $eta > 0$. Choose an integer
  $M > (1 + 3 g) / (1 - g)$ (so $M >= 2$), fix $D = 1$ and
  $L > max(1 / eta, 2)$, and pick tiny increments
  $0 < delta_1 < delta_2 < dots.c < delta_M < 1 / (4 M)$. Place $n =
  M + 3$ points on $bb(R)$ in three groups:
  - Group A ($M$ points): $a_1 = 0$, $a_(r + 1) = a_r + delta_r$ for
    $r = 1, dots, M - 1$;
  - Group B ($2$ points): $b_1 = a_M + D$, $b_2 = b_1 + delta_M$;
  - Outlier $z$: $z = b_2 + L$.

  The MST is the sorted path with ordered edge weights
  $delta_1 < dots.c < delta_(M - 1) < delta_M < D < L$. During steps
  $t = 0, dots, M - 2$, by the same induction as Lemma 1(2), each edge
  $(a_(t + 1), a_(t + 2))$ touches singleton ${a_(t + 2)}$, so Genie
  coalesces Group A into a cluster of size $M$ regardless of $G$. At
  step $t = M - 1$, clusters have sizes $(M, 1, 1, 1)$, so $c_"min" = 1$
  and $(b_1, b_2)$ (weight $delta_M$, connecting two singletons) is
  selected, reaching the 3-cluster state $(M, 2, 1)$ with $k = 3$ and
  $n = M + 3$. Evaluating $G (1, 2, M)$ directly gives
  $
    G (1, 2, M) = ((-2)(1) + (0)(2) + (2)(M)) / (2 (M + 3))
      = (M - 1) / (M + 3).
  $
  Our choice $M > (1 + 3 g) / (1 - g)$ implies $(M - 1) / (M + 3) > g$,
  triggering Case 2 with $c_"min" = 1$ (uniquely attained by ${z}$).
  Hence Genie *forbids* merging Group A with Group B (distance $D = 1$)
  and *forces* merging ${z}$ with Group B at height $h_(n - 2) = L$.
  Finally, at $k = 2$, the remaining two clusters merge along
  $(a_M, b_1)$ at height $h_(n - 1) = D = 1$, producing the inversion
  $h_(n - 1) / h_(n - 2) = 1 / L < eta$.
]

#subsec("D", "Exact O(n log n) Post-MST Dynamic Implementation")

#theorem(
  title: "Exact O(n log n) Post-MST Execution via Dynamic Gini Trees",
  origin: "Original Proof",
)[
  Given the Minimum Spanning Tree $T^*$ of $(X, d)$, the entire Genie
  merge hierarchy for any $g in (0, 1]$ can be computed exactly in
  $O(n log n)$ worst-case time and $O(n)$ auxiliary space, tightening
  the $O(n sqrt(n))$ worst-case scan complexity of the baseline
  disjoint-set list implementation @gagolewski2016genie.
]

#proof[
  Maintain three synchronized $O(n)$-space data structures over the
  contracted tree $T_t = (cal(C)^((t)), E_t)$:
  1. *Augmented Order-Statistic BST $cal(S)$* storing the multiset of
     active sizes ${c_1, dots, c_k}$, where every tree node caches its
     subtree count and subtree sum. Let $S = sum_(a < b) |c_a - c_b|$
     be the unordered pair sum, so $G (c) = S / ((k - 1) n)$. When
     deleting a value $v$ from $cal(S)$, the exact decrement is
     $
       Delta S (v) = sum_(x in cal(S)) |x - v|
         = (L_v v - Sigma_(< v)) + (Sigma_(> v) - R_v v),
     $
     where $L_v, R_v$ are the counts and $Sigma_(< v), Sigma_(> v)$ the
     sums of elements strictly below and above $v$. All four quantities
     are queried in $O(log k)$ during a single root-to-leaf walk.
     Deleting $c_i, c_j$ and inserting $c_i + c_j$ updates $G (c)$
     exactly in $O(log n)$ time.
  2. *Per-Cluster Meldable Min-Heaps* $H_r$ (e.g., leftist heaps)
     storing all edges of $T_t$ incident to $C_r$, keyed by weight.
     Initially each MST edge is inserted into its two endpoints' heaps
     ($2 (n - 1)$ entries total). Crucially, because $T_t$ is acyclic,
     contracting an edge $e^* = (C_i, C_j)$ creates neither parallel
     edges nor self-loops outside $e^*$ itself, and leaves the incident
     edge set—and thus the minimum incident weight $w_"min" (C_r)$—of
     *every third cluster* $C_r$ ($r in.not {i, j}$) strictly unchanged!
     Marking $e^*$ dead, melding $H_i$ and $H_j$ into $H_"new"$ in
     $O(log n)$, and lazily popping dead edges from the root of
     $H_"new"$ maintains exact $w_"min" (C_"new")$. Because each of the
     $2 (n - 1)$ heap entries is popped at most once over the entire
     run, total lazy-deletion work across all $n - 1$ steps is
     $O(n log n)$.
  3. *Global Edge Heap and Size-Distance BST $cal(B)$* keyed by
     $(c_r, w_"min" (C_r), r)$. Case 1 pops the lightest alive MST edge
     in $O(log n)$ amortized time; Case 2 reads the minimum element of
     $cal(B)$ in $O(log n)$ time. Each merge removes two entries from
     $cal(B)$ and inserts one updated entry for $C_"new"$ in $O(log n)$.

  Combining all three components yields $O(n log n)$ worst-case
  post-MST time and $O(n)$ total auxiliary space.
]

#sec("III", "Complexity, Geometry and Operational Profile")

End-to-end runtime decomposes cleanly into two stages: (i) *MST
construction*, requiring $O(n^2)$ time and $O(n)$ memory via Prim's
algorithm in arbitrary metric spaces (or $O(n log n)$ in $bb(R)^2$ via
Delaunay triangulation and near-linear expected time in low-dimensional
Euclidean spaces via dual-tree Boruvka @gagolewski2016genie
@mullner2011modern), and (ii) *Gini-constrained contraction*, requiring
only $O(n log n)$ time and $O(n)$ memory by Theorem 3. Because the MST
is independent of $g$, sweeping a grid of $L$ inequality thresholds
$g_1, dots, g_L$ costs only $O(L n log n)$ after a single MST pass—a
decisive operational advantage for model selection.

Geometrically, Genie inherits Single Linkage's topological adaptability
to curved, nested, and filamentary manifolds while truncating the
heavy-tailed singleton fringe induced by low-density outliers. Because
Theorem 2 proves Genie is invariant to monotone distance warps, it does
not assume isotropic covariance or Euclidean coordinates, operating
natively on arbitrary symmetric dissimilarity matrices. Its primary
geometric caveat is *merge-height non-monotonicity* (Theorem 2(2)):
standard dendrogram visualizations require height clamping or rank-order
plotting when $g < 1$.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Input Domain], [Arbitrary symmetric dissimilarity space $(X, d)$],
    [Time (Metric / Post-MST)], [$O(n^2)$ metric MST / $O(n log n)$ contraction],
    [Space Complexity], [$O(n)$ (with on-the-fly Prim MST, no $n^2$ matrix)],
    [Cluster Geometry], [Non-convex manifolds with outlier absorption],
    [Key Hyperparameter], [Gini threshold $g in (0, 1]$ (default $g = 0.3$)],
    [Invariances], [Permutations, arbitrary monotone transforms $phi (d)$],
    [Ultrametricity], [Violated for $g in (0, 1)$ (inversions possible)],
    [Determinism], [Exact up to tied inter-point edge weights],
  ),
  caption: [Operational and theoretical profile of Genie.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

Genie excels in applications combining non-convex cluster geometry with
sparse outlier contamination, where classical Single Linkage breaks down
and Ward's method imposes false sphericity:

+ *Astronomical Filament and Halo Identification:* Galaxy redshift
  surveys exhibit elongated cosmic-web filaments embedded in sparse
  void galaxies. Single Linkage isolates void galaxies as singleton
  clusters at the top of the hierarchy; Genie's Gini constraint absorbs
  isolated field galaxies while preserving MST connectivity along
  non-convex filaments @gagolewski2016genie.
+ *Geospatial Trajectory and Crime-Hotspot Segmentation:* Urban road
  networks and GPS traces lie along 1D non-convex corridors accompanied
  by scattered GPS drift pings. Because Genie is invariant to monotone
  metric transforms and tracks MST paths, it recovers linear transit
  corridors without splitting long avenues into spherical patches.
+ *Single-Cell Transcriptomics and Manifold Phenotyping:* In low-
  dimensional UMAP/diffusion embeddings of single-cell RNA-seq data,
  developmental lineages form curved branching trajectories alongside
  doublet outliers. Setting $g approx 0.2 - 0.4$ prevents doublets from
  consuming cluster budget while respecting continuous branching paths.
+ *Industrial Sensor Anomaly-Resilient Regime Discovery:* Multi-sensor
  telemetry cycles traverse non-convex operating manifolds interrupted
  by transient spike readings. Genie absorbs transient spikes during
  constrained steps while keeping true multi-mode operating regimes
  cleanly separated at target cut $K$.
+ *Fast Multi-Resolution Exploratory Benchmarks:* Because re-running
  Genie for multiple $g$ values reuses the exact same MST in
  $O(n log n)$ time, analysts can inspect stability curves across
  $g in [0.1, 0.9]$ interactively—connecting Single Linkage ($g = 1$)
  to strongly balanced cuts ($g -> 0^+$) in milliseconds.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Four failure regimes require caution. First, *genuinely imbalanced true
clusters*: if a dataset truly contains one large cluster of size $0.95 n$
and a genuine tiny minority class of size $0.05 n$, a small $g$ forces
premature absorption of the valid minority class into the majority
cluster (the flip side of outlier suppression). Second, *monotone
density bridges* (Lemma 1(2)): a dense chain of progressively spaced
points bridging two clusters bypasses Case 2 and still triggers
Single-Linkage chaining. Third, *high-dimensional distance
concentration* @beyer1999nearest: as ambient dimension grows, MST edge
weights concentrate, rendering nearest-neighbor topology noisy. Fourth,
*dendrogram inversions* (Theorem 2(2)) invalidate height-based flat
cutting; cuts must be specified by cluster count $K$ or post-processed
via monotonic height repair.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Non-convex + sparse noise], [Ideal: MST tracks manifolds, Gini absorbs outliers], [Tune $g in [0.2, 0.5]$],
    [Imbalanced rare classes], [Poor for small $g$: forces rare-class absorption], [Raise $g$ or use HDBSCAN @campello2013hdbscan],
    [Dense inter-cluster bridges], [Vulnerable: gradient chaining persists (Lemma 1)], [Use HCA-DBSCAN or Ward],
    [Compact Gaussian blobs], [Moderate: no variance minimization objective], [Ward @ward1963hierarchical or BUBBLE],
    [High-D isotropic spaces], [Degrades as MST hubness/concentration grows], [Reduce dim first; cosine metric],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for Genie.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Genie sits on a continuum between topological linkage and density-based
hierarchies. Unlike Lance--Williams recurrences @lance1967general,
Genie's merge criterion is *state-dependent* on global partition
inequality, trading ultrametricity for outlier robustness
@kleinberg2002impossibility @dasgupta2016cost. It also contrasts sharply
with HDBSCAN @campello2013hdbscan: whereas HDBSCAN modifies edge weights
via mutual reachability distances and discards noise points outside the
final partition, Genie preserves exact raw metric edges and retains a
complete partition of all $n$ points via cardinality-gated MST
contraction.

#proposition(
  title: "Exact Cut-Equivalence Between Genie and Single Linkage Under Balance",
  origin: "Original Proof",
)[
  Suppose an unconstrained Single-Linkage dendrogram on $(X, d)$
  produces partitions $cal(C)_"SL"^((t))$ satisfying
  $G (cal(C)_"SL"^((t))) <= g$ for all steps $t = 0, dots, T$. Then
  Genie with threshold $g$ produces the exact same prefix of merges and
  identical partition $cal(C)_"Genie"^((T)) = cal(C)_"SL"^((T))$ with
  strictly monotone merge heights $h_1 <= h_2 <= dots.c <= h_T$.
]

#proof[
  At $t = 0$, $cal(C)_"Genie"^((0)) = cal(C)_"SL"^((0))$ and $G = 0 <= g$.
  Inductively, if $cal(C)_"Genie"^((t)) = cal(C)_"SL"^((t))$ for some
  $t < T$, then $G (cal(C)_"Genie"^((t))) = G (cal(C)_"SL"^((t))) <= g$,
  so Genie executes an unconstrained Single-Linkage step selecting the
  exact same lightest MST edge as Single Linkage, advancing both
  partitions identically. Since unconstrained MST contraction selects
  edges in non-decreasing weight order, $h_1 <= dots.c <= h_T$.
]

Empirical results on the ClusterOverview 4-dataset benchmark align
closely with these theoretical predictions @rousseeuw1987silhouettes
@hubert1985comparing. On Multi-Shape (curved non-convex manifolds with
scattered geometry), Genie achieves a high Silhouette of *0.718*
($"ARI" = "N/A"$), decisively outperforming MONA (0.304) and BUBBLE
(0.452) by tracking MST connectivity without singleton fragmentation.
On Iris, Genie reaches $"Sil" = 0.447$ and the highest Adjusted Rand
Index (*0.606*) among Algorithms 68--71, confirming that moderate Gini
regularization recovers Iris's three balanced species cleanly without
outlier distortion. Conversely, on Customers ($"Sil" = 0.266$), natural
skew across customer expenditure tiers clashes with uniform size-balance
pressure, forcing uneven sub-segments together, while on High-D
($"Sil" = 0.296$, $"ARI" = "N/A"$), high-dimensional MST edge
concentration attenuates separation. Its mean runtime of *2.781 ms*
reflects full pairwise MST evaluation plus Gini state bookkeeping—
slower than core-compressed HCA-DBSCAN (1.008 ms) and CF-tree BUBBLE
(1.215 ms), yet completely practical for small-to-medium datasets and
rapid multi-$g$ sweeps.
