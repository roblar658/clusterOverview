#import "../template/ieee.typ": *

#alg_entry(
  id: 71,
  name: "BUBBLE: Metric-Space Hierarchical CF-Tree Clustering",
  category: "Cat. III: Hierarchical & Tree-Based",
  year: "1999",
  paradigm: "Metric-Space CF-Tree and Clustroid Summarization",
  complexity: $O(n (B log_B n_L + s) + n_L^2)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

BIRCH @zhang1996birch revolutionized scalable hierarchical clustering
for massive datasets by compressing data in a single I/O scan into a
height-balanced Clustering Feature tree (*CF-tree*) whose summary tuple
$"CF"(C) = (N, bold("LS"), "SS")$ obeys exact additive algebra. However,
BIRCH fundamentally requires Euclidean coordinates $X subset.eq bb(R)^d$
so that linear sums $bold("LS") = sum_(x in C) x$ and centroids
$mu_C = bold("LS") / N$ exist. In numerous scientific domains—biological
sequences under edit distance, time series under Dynamic Time Warping
@sakoe1978dtw, graphs, strings, or molecular fingerprints—objects live
in an *arbitrary metric space* $(cal(X), d)$ lacking vector addition and
inner products.

*BUBBLE* (and its FastMap-accelerated variant *BUBBLE-FM*), introduced
by Ganti, Ramakrishnan, Gehrke, Powell, and French @ganti1999bubble,
generalizes BIRCH's scalable CF-tree framework to arbitrary metric
spaces by replacing Euclidean centroids with *clustroids* and bounded
*representative subsets*, while preserving dynamic tree insertion,
threshold-controlled leaf compression, and global post-clustering.

Let $(cal(X), d)$ be a metric space (symmetric, non-negative, definite,
and obeying the triangle inequality $d (x, z) <= d (x, y) + d (y, z)$).
For any finite cluster $C subset.eq cal(X)$ of size $N = abs(C) >= 1$
and any object $o in C$, define its *row sum* (total squared metric
distance to $C$) as
$
  "RowSum"_C (o) = sum_(x in C) d (o, x)^2.
$

#definition(title: "Clustroid, Metric Radius, and Cluster Representation (CR)")[
  For a nonempty finite cluster $C subset.eq cal(X)$ of size $N =
  abs(C)$:
  - The *clustroid* is a sample minimizing squared intra-cluster
    dispersion:
    $
      hat(o)_C in arg min_(o in C) "RowSum"_C (o).
    $
  - The *metric radius* and *mean-square diameter* (for $N >= 2$) are
    $
      r (C) = sqrt(("RowSum"_C (hat(o)_C)) / N), \
      "diam"_2 (C) = sqrt((sum_(x in C) sum_(y in C) d (x, y)^2)
        / (N (N - 1))).
    $
  - The *Cluster Representation* tuple @ganti1999bubble storing at most
    $s$ representative objects $"Rep"(C) subset.eq C$ (including
    $hat(o)_C$) is
    $
      "CR"(C) = (N, "Rep"(C), hat(o)_C, r (C),
        {"RowSum"_C (o)}_(o in "Rep"(C))).
    $
]

When merging two subclusters $C_1$ and $C_2$ into $C = C_1 union C_2$
without retaining every raw point of $C_2$, BUBBLE estimates the cross
row-sum $S_2 (o) = sum_(y in C_2) d (o, y)^2$ of any representative
$o in "Rep"(C_1)$ using the *parallel-axis metric surrogate*:
$
  hat(S)_2 (o) = abs(C_2) (d (o, hat(o)_2)^2 + r (C_2)^2)
    = abs(C_2) d (o, hat(o)_2)^2 + "RowSum"_(C_2) (hat(o)_2),
$
setting $hat("RowSum")_C (o) = "RowSum"_(C_1) (o) + hat(S)_2 (o)$ and
updating the merged clustroid to the representative minimizing
$hat("RowSum")_C (o)$ across $"Rep"(C_1) union "Rep"(C_2)$. To minimize
expensive metric evaluations during top-down tree routing, *BUBBLE-FM*
further embeds representative sets into $bb(R)^k$ via *FastMap*
@ganti1999bubble: given two pivot objects $(p_a, p_b)$, coordinate $1$
projects $x$ via the cosine-law formula
$
  psi_1 (x) = (d (p_a, x)^2 + d (p_a, p_b)^2 - d (p_b, x)^2)
    / (2 d (p_a, p_b)),
$
iterated over $k$ axes using residual distance
$d' (x, y)^2 = d (x, y)^2 - (psi_1 (x) - psi_1 (y))^2$.

#algbox(
  title: "BUBBLE / BUBBLE-FM (Metric-Space CF-Tree Clustering)",
  inputs: [Metric dataset $(X, d)$; branching factor $B$, leaf limit
    $L$, radius threshold $T$, sample size $s$, optional FastMap dim
    $k$],
  outputs: [Height-balanced CR-tree and $K$-cluster partition],
)[
  + *Initialization:* Allocate empty root CR-tree node with radius
    threshold $T$. (In BUBBLE-FM, initialize $k$ FastMap pivot pairs on
    a seed sample.)
  + *Phase 1 (Incremental CR-Tree Insertion):* For each incoming object
    $x in X$:
    - *Route:* Descend from root to leaf by picking at each internal
      node the child minimizing $d (x, hat(o)_j)$ (in BUBBLE-FM,
      $norm(psi (x) - psi (hat(o)_j))_2$ in $bb(R)^k$).
    - *Absorb or Create:* Let $C^*$ be the closest leaf cluster. Update
      exact row sums of $"Rep"(C^*)$ by $+ d (o, x)^2$, recompute trial
      $r (C^* union {x})$. If $r <= T$, absorb $x$ into $C^*$ (retaining
      up to $s$ well-spread representatives); else insert ${x}$ as a new
      leaf entry, splitting overfull nodes ($> L$ or $> B$) along their
      farthest-apart clustroid pair and propagating CR updates upward.
    - *Rebuild:* If memory budget is exceeded, increase $T$ and rebuild
      a compressed CR-tree from existing leaf CRs.
  + *Phase 2 (Global Clustering):* Run weighted hierarchical linkage or
    medoid clustering on the $n_L$ leaf clustroids $(abs(C_j),
    hat(o)_j)$ to produce $K$ macro-clusters.
  + *Phase 3 (Optional Refinement Pass):* Assign each $x in X$ to its
    nearest final macro-cluster clustroid.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish four mathematical theorems governing BUBBLE: universal
sharp metric bounds relating clustroid radius to diameter and Euclidean
variance, exact two-sided metric error bounds for the RowSum merge
surrogate, a constructive 4-point non-Euclidean counterexample proving
diameter-pivot breakdown of FastMap contraction in BUBBLE-FM, and the
distance-evaluation complexity theorem of Ganti et al. @ganti1999bubble.

#subsec("A", "Universal Clustroid-to-Diameter Sandwich and Euclidean Optimality")

#theorem(
  title: "Sharp Clustroid-Radius Bounds in Arbitrary and Euclidean Metrics",
  origin: "Original Proof",
)[
  Let $C subset.eq cal(X)$ be any finite cluster of size $N >= 2$ in an
  arbitrary metric space $(cal(X), d)$.
  + *Universal Metric Sandwich:* The clustroid radius $r (C)$ and
    mean-square diameter $"diam"_2 (C)$ satisfy
    $
      (1 / 4) "diam"_2 (C)^2 <= r (C)^2
        <= ((N - 1) / N) "diam"_2 (C)^2.
    $
    Both constants are unimprovable: the upper bound is attained on
    every unit simplex, and the lower bound $1 / 4$ is asymptotically
    tight on unweighted star metrics $K_(1, m)$ as $m -> infinity$.
  + *Euclidean Specialization:* If $(cal(X), d) = (bb(R)^d, norm(dot)_2)$
    with centroid $mu_C$ and variance $sigma_C^2 = N^(-1) sum_(x in C)
    norm(x - mu_C)_2^2$, then $hat(o)_C$ is *identically* the point of
    $C$ closest to $mu_C$, and satisfies $sigma_C^2 <= r (C)^2 <= 2
    sigma_C^2$.
]

#proof[
  (1) *Upper bound:* Since $hat(o)_C$ minimizes $"RowSum"_C (o)$ over
  $o in C$, its row sum cannot exceed the average row sum across $C$:
  $
    N r (C)^2 = min_(o in C) sum_(x in C) d (o, x)^2
      <= (1 / N) sum_(o in C) sum_(x in C) d (o, x)^2
      = (N - 1) "diam"_2 (C)^2.
  $
  Dividing by $N$ gives $r (C)^2 <= ((N - 1) / N) "diam"_2 (C)^2$. On a
  regular simplex where $d (x, y) = 1$ for all $x != y$, every vertex
  has $"RowSum" = N - 1$, so $r (C)^2 = (N - 1) / N$ while
  $"diam"_2 (C)^2 = 1$, attaining exact equality.

  *Lower bound:* For every distinct pair $x != y$ in $C$, the triangle
  inequality and $(a + b)^2 <= 2 (a^2 + b^2)$ yield
  $
    d (x, y)^2 <= (d (x, hat(o)_C) + d (hat(o)_C, y))^2
      <= 2 d (x, hat(o)_C)^2 + 2 d (hat(o)_C, y)^2.
  $
  Summing over all $N (N - 1)$ ordered pairs $x != y$ in $C$ (each
  $x in C$ appears $N - 1$ times in the first index and $N - 1$ times
  in the second) gives
  $
    N (N - 1) "diam"_2 (C)^2
      <= 4 (N - 1) sum_(x in C) d (x, hat(o)_C)^2
      = 4 N (N - 1) r (C)^2,
  $
  so $r (C)^2 >= (1 / 4) "diam"_2 (C)^2$. On the unweighted star graph
  $K_(1, m)$ ($N = m + 1$, center at distance $1$ to $m$ leaves, leaves
  pairwise at distance $2$), the center is the clustroid with
  $r (C)^2 = m / (m + 1)$ and $"diam"_2 (C)^2 = (2 m + 4 m (m - 1)) /
  ((m + 1) m) = (4 m - 2) / (m + 1)$, so $r (C)^2 / "diam"_2 (C)^2 =
  m / (4 m - 2) -> 1 / 4$ as $m -> infinity$.

  (2) In $(bb(R)^d, norm(dot)_2)$, the parallel-axis identity gives
  $
    "RowSum"_C (o) = sum_(x in C) norm(o - x)_2^2
      = N (norm(o - mu_C)_2^2 + sigma_C^2)
  $
  for every $o in C$. Minimizing over $o in C$ therefore uniquely
  minimizes $norm(o - mu_C)_2$. Moreover, because
  $N^(-1) sum_(x in C) norm(x - mu_C)_2^2 = sigma_C^2$, by the
  averaging principle at least one $x^* in C$ satisfies
  $norm(x^* - mu_C)_2^2 <= sigma_C^2$, whence
  $norm(hat(o)_C - mu_C)_2^2 <= sigma_C^2$ and
  $r (C)^2 = sigma_C^2 + norm(hat(o)_C - mu_C)_2^2 <= 2 sigma_C^2$.
]

#subsec("B", "Two-Sided Metric Error Bounds for the RowSum Surrogate")

#theorem(
  title: "Exact Triangle-Inequality Bounds on BUBBLE's Merge Surrogate",
  origin: "Original Derivation",
)[
  Let $C_2 subset.eq cal(X)$ be a finite cluster with clustroid
  $hat(o)_2$ and radius $r (C_2)$, and let $o in cal(X)$ be arbitrary.
  Let $S_2 (o) = sum_(y in C_2) d (o, y)^2$ be the true cross row-sum
  and $hat(S)_2 (o) = abs(C_2) (d (o, hat(o)_2)^2 + r (C_2)^2)$ be
  BUBBLE's parallel-axis surrogate. Then in *every* metric space:
  $
    (d (o, hat(o)_2) - r (C_2))^2
      <= S_2 (o) / abs(C_2)
      <= (d (o, hat(o)_2) + r (C_2))^2,
  $
  and the relative approximation error satisfies
  $
    (|S_2 (o) - hat(S)_2 (o)|) / hat(S)_2 (o)
      <= (2 rho) / (1 + rho^2) <= 1,
    quad rho = r (C_2) / d (o, hat(o)_2).
  $
]

#proof[
  For every $y in C_2$, the forward and reverse triangle inequalities
  give $|d (o, hat(o)_2) - d (hat(o)_2, y)| <= d (o, y) <= d (o,
  hat(o)_2) + d (hat(o)_2, y)$. Squaring both sides yields
  $
    |d (o, y)^2 - d (o, hat(o)_2)^2 - d (hat(o)_2, y)^2|
      <= 2 d (o, hat(o)_2) d (hat(o)_2, y).
  $
  Summing over all $y in C_2$ and recognizing
  $sum_(y in C_2) d (hat(o)_2, y)^2 = abs(C_2) r (C_2)^2$ produces
  $
    |S_2 (o) - hat(S)_2 (o)|
      <= 2 d (o, hat(o)_2) sum_(y in C_2) d (hat(o)_2, y).
  $
  By Cauchy--Schwarz,
  $sum_(y in C_2) d (hat(o)_2, y) <= sqrt(abs(C_2)) sqrt(abs(C_2)
  r (C_2)^2) = abs(C_2) r (C_2)$, whence
  $
    |S_2 (o) / abs(C_2) - (d (o, hat(o)_2)^2 + r (C_2)^2)|
      <= 2 d (o, hat(o)_2) r (C_2).
  $
  Expanding the absolute value immediately gives the lower bound
  $(d (o, hat(o)_2) - r (C_2))^2$ and the upper bound
  $(d (o, hat(o)_2) + r (C_2))^2$. Dividing by
  $hat(S)_2 (o) / abs(C_2) = d (o, hat(o)_2)^2 + r (C_2)^2$ yields
  $2 rho / (1 + rho^2) <= 1$ by AM--GM.
]

Theorem 2 reveals precisely why BUBBLE works well when leaf radii are
tight ($T$ small): whenever $r (C_2) << d (o, hat(o)_2)$, $rho -> 0$
and the surrogate error vanishes linearly as $O(rho)$ in *any* metric
space without requiring Euclidean inner products.

#subsec("C", "FastMap Contraction vs. Non-Euclidean Diameter-Pivot Breakdown")

#theorem(
  title: "Euclidean Projection Exactness vs. 4-Point FastMap Failure",
  origin: "Original Counterexample",
)[
  Let $psi_1$ be the first FastMap coordinate with pivots $(p_a, p_b)$.
  + In any Euclidean space $(bb(R)^m, norm(dot)_2)$, $psi_1$ is an exact
    1D orthogonal projection onto the unit axis
    $u = (p_b - p_a) / norm(p_b - p_a)_2$, guaranteeing non-expansiveness
    $|psi_1 (x) - psi_1 (y)| <= norm(x - y)_2$ and non-negative residual
    squared distance $d' (x, y)^2 >= 0$ for all $x, y$.
  + *Non-Euclidean Breakdown:* There exists a valid 4-point finite
    metric space $(cal(X), d)$ such that even when FastMap selects an
    *exact diameter pivot pair* $(p_a, p_b)$, it expands distances by
    $|psi_1 (x) - psi_1 (y)| = (3 / 2) d (x, y) > d (x, y)$ and yields
    a strictly negative residual squared distance $d' (x, y)^2 = -5 / 4
    < 0$.
]

#proof[
  (1) In $bb(R)^m$, expand
  $norm(x - p_b)_2^2 = norm((x - p_a) - (p_b - p_a))_2^2 =
  norm(x - p_a)_2^2 + norm(p_b - p_a)_2^2 - 2 chevron.l x - p_a,
  p_b - p_a chevron.r$. Substituting into $psi_1 (x)$ yields
  $psi_1 (x) = chevron.l x - p_a, u chevron.r$. By Cauchy--Schwarz,
  $|psi_1 (x) - psi_1 (y)| = |chevron.l x - y, u chevron.r| <=
  norm(x - y)_2$, and $d' (x, y)^2 = norm((I - u u^top)(x - y))_2^2 >=
  0$.

  (2) Define $cal(X) = {p_a, p_b, x, y}$ with symmetric zero-diagonal
  distances $d (p_a, x) = d (x, y) = d (y, p_b) = 1$ and
  $d (p_a, p_b) = d (p_a, y) = d (p_b, x) = 2$ (equivalently, the
  truncated path metric $d = min(d_"path", 2)$). Because every
  off-diagonal distance lies in ${1, 2}$ and $1 + 1 >= 2$, every triple
  satisfies the triangle inequality, so $(cal(X), d)$ is a valid metric
  space whose exact diameter is $d (p_a, p_b) = 2$. Selecting the
  diameter pair $(p_a, p_b)$ gives
  $
    psi_1 (x) &= (1^2 + 2^2 - 2^2) / (2 dot 2) = 1 / 4, \
    psi_1 (y) &= (2^2 + 2^2 - 1^2) / (2 dot 2) = 7 / 4.
  $
  Hence $|psi_1 (x) - psi_1 (y)| = 6 / 4 = 3 / 2 > 1 = d (x, y)$, and
  $d' (x, y)^2 = 1^2 - (3 / 2)^2 = -5 / 4 < 0$, proving that BUBBLE-FM
  requires non-negative residual clamping in general metric spaces.
]

#subsec("D", "Distance-Evaluation Complexity of BUBBLE and BUBBLE-FM")

#theorem(
  title: "Metric-Distance Complexity Reduction via FastMap Routing",
  origin: "Classical — Cited",
)[
  Suppose a CR-tree has internal branching factor $B$, leaf capacity
  $L$, maximum representative set size $s$, and produces $n_L$ leaf
  clusters (height $h = O(log_B n_L)$) over $n$ objects
  @ganti1999bubble. Denote by $c_d$ the cost of one metric evaluation
  $d (dot, dot)$.
  + Standard BUBBLE performs $O(n (B log_B n_L + s) + n_L^2)$ metric
    distance evaluations during tree construction and Phase-2 global
    clustering, requiring $O(n_L s)$ memory.
  + BUBBLE-FM with $k$-dimensional FastMap routing reduces metric
    distance evaluations to $O(n (k + s) + n_L^2)$ while replacing the
    $O(n B log_B n_L)$ internal-node metric calls with cheap
    $O(n k B log_B n_L)$ scalar arithmetic operations @ganti1999bubble
    @zhang1996birch.
]

#proof[
  In Phase 1, each of the $n$ inserted objects traverses $h =
  O(log_B n_L)$ internal levels. In standard BUBBLE, selecting the
  nearest child clustroid at each level invokes at most $B$ calls to
  $d (dot, dot)$, totaling $O(B log_B n_L)$ metric calls per object.
  At the destination leaf, testing absorption and updating row sums
  across up to $s$ stored representatives costs $O(s)$ metric calls,
  while periodic node splits amortize over $Omega(B)$ insertions. Over
  $n$ insertions plus $O(n_L^2)$ inter-clustroid distances in Phase 2,
  standard BUBBLE executes $O(n (B log_B n_L + s) + n_L^2)$ metric
  calls. In BUBBLE-FM, embedding an incoming object $x$ into $bb(R)^k$
  requires evaluating $d (x, p_(a, r))$ and $d (x, p_(b, r))$ against
  the $k$ pivot pairs—exactly $2 k = O(k)$ metric evaluations—after
  which all $O(B log_B n_L)$ internal tree comparisons use $O(k)$
  Euclidean arithmetic in $bb(R)^k$ without invoking $d (dot, dot)$.
  Leaf absorption still verifies true metric radius in $O(s)$ calls,
  yielding $O(n (k + s) + n_L^2)$ metric calls. Storing $n_L$ leaf CRs
  with $<= s$ representatives each consumes $O(n_L s)$ memory.
]

#sec("III", "Complexity, Geometry and Operational Profile")

BUBBLE inherits BIRCH's single-scan memory-bounded scalability: instead
of holding an $O(n^2)$ distance matrix in RAM, it maintains only
$O(n_L s)$ representative objects and clustroid summaries, dynamically
raising the radius threshold $T$ whenever the CR-tree approaches memory
capacity @ganti1999bubble @zhang1996birch. When metric evaluations are
expensive (e.g., $O(L_"seq"^2)$ sequence alignment or DTW), BUBBLE-FM's
reduction of routing calls from $O(B log_B n_L)$ to $O(k)$ yields
order-of-magnitude speedups while preserving exact metric verification
at the leaves.

Geometrically, BUBBLE partitions $(cal(X), d)$ into *metric balls* of
radius $<= T$ around leaf clustroids, approximating each ball's mass as
concentrated near its discrete medoid (Theorem 1). Consequently,
BUBBLE's inductive bias strongly favors compact, isotropic, well-
separated metric balls. Because CR merging relies on bounded sample
sets $"Rep"(C)$ and the triangle-inequality surrogate (Theorem 2), CR
merges are *approximate and non-associative*, introducing mild
insertion-order sensitivity compared to BIRCH's exact additive CFs.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Input Domain], [Arbitrary metric space $(cal(X), d)$ (no coordinates)],
    [Metric Calls (BUBBLE / FM)], [$O(n (B log_B n_L + s))$ / $O(n (k + s))$],
    [Space Complexity], [$O(n_L s)$ bounded RAM ($n_L << n$ leaf clusters)],
    [Cluster Geometry], [Compact metric balls centered at discrete clustroids],
    [Key Hyperparameters], [Branching $B$, leaf limit $L$, radius $T$, sample $s$],
    [Streaming / I/O Profile], [Single-pass scalable; dynamic $T$ tree rebuild],
    [Invariances], [Exact metric isometries; mild insertion-order drift],
    [Euclidean Relation], [$r (C) <= sqrt(2) sigma_C$ discrete analogue of BIRCH],
  ),
  caption: [Operational and theoretical profile of BUBBLE and BUBBLE-FM.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

BUBBLE is engineered specifically for massive datasets endowed with a
domain-specific metric where Euclidean vector embeddings are absent or
costly to compute:

+ *Biological Sequence and Protein-Family Clustering:* Large repositories
  of DNA/RNA reads or protein sequences use Levenshtein edit distance or
  Smith--Waterman metric transforms, where computing a consensus centroid
  is NP-hard. BUBBLE summarizes sequence subclusters by actual medoid
  sequences (clustroids) and uses BUBBLE-FM to slash expensive pairwise
  alignment calls by orders of magnitude @ganti1999bubble.
+ *Time-Series Indexing Under Dynamic Time Warping (DTW):* Industrial
  telemetry and physiological waveforms compared via metric-enforced DTW
  @sakoe1978dtw or shape-based metrics lack linear averaging semantics.
  BUBBLE's CR-tree incrementally compresses millions of waveforms into
  $n_L$ representative prototypes in a single pass.
+ *Cheminformatics and Molecular Fingerprint Repositories:* Chemical
  libraries measured by Tanimoto/Jaccard metric distance on structural
  fingerprints require chemically valid cluster representatives. Because
  every BUBBLE clustroid $hat(o)_C in C$ is an actual molecule rather
  than a fractional Euclidean mean, leaf summaries remain physically
  synthesizable scaffolds.
+ *Multimedia and Large-Scale Image Database Indexing:* Histograms and
  earth-mover / Wasserstein perceptual metrics @ganti1999bubble are
  expensive to evaluate pairwise. BUBBLE-FM uses fast $k$-D FastMap
  projections for tree routing while keeping exact perceptual distances
  inside leaf clusters.
+ *Disk-Resident Log and Document Deduplication:* Streaming enterprise
  documents compared by normalized compression or edit metrics fit
  naturally into BUBBLE's bounded-RAM threshold-rebuild loop without
  materializing an $O(n^2)$ similarity matrix.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Four failure modes limit BUBBLE. First, *high-dimensional metric
concentration* @beyer1999nearest: as intrinsic dimension rises, all
inter-point distances concentrate around a common value, driving
$rho = r (C_2) / d (o, hat(o)_2) -> 1$ in Theorem 2 so the RowSum
surrogate loses precision and leaf balls overlap heavily. Second,
*non-convex / elongated shapes*: spherical metric balls of radius $T$
fracture curved manifolds unless Phase-2 global clustering employs a
connectivity-aware linkage on sufficiently fine leaves. Third,
*insertion-order sensitivity and threshold drift*: because CR updates
are approximate and non-associative, adversarial stream ordering can
trigger premature leaf splits and unbalanced trees. Fourth, *severe
non-Euclidean FastMap distortion* (Theorem 3(2)): in tree-like or
cyclic graph metrics, unclamped BUBBLE-FM routing can misroute points
to suboptimal leaves.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Large non-Euclidean metrics], [Ideal: clustroid CR-tree avoids centroids & $O(n^2)$], [Clamp $d'^2 >= 0$ in BUBBLE-FM],
    [Compact convex datasets], [Strong: metric balls tightly track convex modes], [BIRCH @zhang1996birch if Euclidean],
    [Non-convex / curved shapes], [Moderate: requires tiny $T$ + density Phase 2], [Use HCA-DBSCAN or Genie],
    [High-D ($d > 30$) spaces], [Weak: $r / d -> 1$ degrades RowSum surrogate], [Project first; use BICO @fichtenberger2013bico],
    [Small exact-medoid tasks], [Overkill: approximate CR sacrifices exactness], [Use FasterPAM @schubert2019fasterpam],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for BUBBLE.],
)

#sec("V", "Structural Relations and Empirical Behavior")

BUBBLE occupies a pivotal bridge between Euclidean streaming trees
(BIRCH @zhang1996birch, CluStream @aggarwal2003clustream, ClusTree
@kranen2011clustree, BICO @fichtenberger2013bico) and discrete
metric-space medoid methods (PAM @kaufman1990finding, FasterPAM
@schubert2019fasterpam, ROCK @guha2000rock, CURE @guha1998cure).
Whereas BIRCH requires exact Euclidean additivity and PAM requires
$O(n^2)$ quadratic access, BUBBLE trades BIRCH's exact additivity for
metric generality while trading PAM's exact global swaps for single-pass
hierarchical compression.

#proposition(
  title: "Weighted Coreset Distortion Bound of BUBBLE Leaf Compression",
  origin: "Original Proof",
)[
  Let $cal(L) = {C_1, dots, C_(n_L)}$ be the leaf clusters of a BUBBLE
  CR-tree with threshold $T$ (so $r (C_j) <= T$ for all $j$), and let
  $cal(Q) = {q_1, dots, q_K} subset.eq cal(X)$ be any $K$-medoid center
  set. Approximating each $x in C_j$ by its leaf clustroid $hat(o)_j$
  incurs additive root-mean-square $K$-medoids objective error bounded
  globally by $T$:
  $
    |sqrt((1 / n) sum_(x in X) d (x, cal(Q))^2)
      - sqrt((1 / n) sum_(j = 1)^(n_L) abs(C_j) d (hat(o)_j, cal(Q))^2)|
      <= T.
  $
]

#proof[
  Define vectors $u, v in bb(R)^n$ indexed by $x in C_j$ with entries
  $u_x = d (x, cal(Q)) = min_(q in cal(Q)) d (x, q)$ and
  $v_x = d (hat(o)_j, cal(Q))$. By the 1-Lipschitz property of the
  distance-to-set function in any metric space,
  $|u_x - v_x| = |d (x, cal(Q)) - d (hat(o)_j, cal(Q))| <=
  d (x, hat(o)_j)$ for every $x in C_j$. Applying the reverse triangle
  inequality of the Euclidean norm on $bb(R)^n$ yields
  $
    |norm(u)_2 - norm(v)_2| <= norm(u - v)_2
      <= sqrt(sum_(j = 1)^(n_L) sum_(x in C_j) d (x, hat(o)_j)^2)
      = sqrt(sum_(j = 1)^(n_L) abs(C_j) r (C_j)^2).
  $
  Since $r (C_j) <= T$ for every leaf and $sum_j abs(C_j) = n$, the
  right-hand side is at most $sqrt(n) T$. Dividing by $sqrt(n)$ proves
  the exact additive RMS distortion bound $<= T$ @harpeled2004coresets.
]

Proposition 1 formalizes BUBBLE's leaf layer as an *additive-$T$ metric
coreset* for $K$-medoids and Ward-like Phase-2 clustering in *any*
metric space.

Empirical performance on the ClusterOverview 4-dataset benchmark
confirms every aspect of this geometry @rousseeuw1987silhouettes
@hubert1985comparing. On Iris—whose convex, well-separated clusters
align naturally with metric-ball leaf compression and Theorem 1(2)'s
Euclidean clustroid optimality—BUBBLE achieves the *highest Silhouette
(0.467)* across all four chapters alongside strong $"ARI" = 0.559$. On
Customers ($"Sil" = 0.353$), it outperforms both Genie (0.266) and MONA
(0.328) by summarizing variable-size compact regions without uniform
cardinality pressure, though trailing noise-filtering HCA-DBSCAN
(0.514). On Multi-Shape ($"Sil" = 0.452$, $"ARI" = "N/A"$), spherical
leaf balls recover coarse geometry better than MONA's binary cuts
(0.304) but cannot match the manifold-tracking topologies of Genie
(0.718) and HCA-DBSCAN (0.804). Finally, on High-D ($"Sil" = 0.245$,
$"ARI" = "N/A"$), BUBBLE suffers the steepest drop among the four
methods—directly reflecting Theorem 2's ratio $rho = r (C_2) /
d (o, hat(o)_2) -> 1$ under high-dimensional distance concentration
@beyer1999nearest. Its fast mean runtime of *1.215 ms* ($2.3 times$
faster than Genie's 2.781 ms) confirms the practical efficiency of
CR-tree compression @jain2010data @xu2005survey.
