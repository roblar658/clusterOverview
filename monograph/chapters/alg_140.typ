#import "../template/ieee.typ": *

#alg_entry(
  id: 140,
  name: "MAFIA: Merging of Adaptive Finite Intervals",
  category: "Cat. VI: Grid, Subspace & Projection",
  year: "1999",
  paradigm: "Adaptive-Grid Bottom-Up Subspace Clustering",
  complexity: $cal(O)(n (d + |cal(U)_"cand"|))$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

CLIQUE @agrawal1998clique pioneered bottom-up grid-based
subspace clustering by partitioning each coordinate axis into
$xi$ equal-width intervals and applying Apriori pruning to dense
grid cells. While conceptually elegant, uniform quantization
suffers two severe drawbacks: (i) rigid grid boundaries split
dense clusters across neighboring cells, inflating false
negatives exponentially with subspace dimension $k$; and (ii)
allocating equal resolution across empty background regions
inflates the combinatorial cell space $xi^k$ unnecessarily.
MAFIA (Merging of Adaptive Finite Intervals), proposed by Goil,
Nagesh, and Choudhary @goil1999mafia, replaces CLIQUE's fixed
grid with *data-driven variable-width histograms* along every
dimension and couples candidate evaluation with a shared-nothing
data-parallel architecture designed for massive multi-terabyte
datasets.

Let $cal(X) = \{bold(x)_1, ..., bold(x)_n\} subset bb(R)^d$ lie
within the bounding box
$product_(j=1)^d [a_j, b_j]$ of coordinate spans
$L_j = b_j - a_j > 0$. For each axis $j in A = \{1, ..., d\}$,
MAFIA first builds a fine uniform micro-histogram of $B$ bins of
width $h_j = L_j / B$, recording empirical bin counts
$c_(j, t)$ and normalized densities
$rho_(j, t) = c_(j, t) / (n h_j)$ for $t = 1, ..., B$. Adjacent
micro-bins with statistically similar density—within relative
merge tolerance $beta in (0, 1)$—are merged left-to-right into
$m_j << B$ variable-width _Adaptive Finite Intervals_ (AFIs)
$cal(I)_j = \{I_(j, 1), ..., I_(j, m_j)\}$ partitioning
$[a_j, b_j]$, with interval widths $w_(j, r) = |I_(j, r)|$
satisfying $sum_(r=1)^(m_j) w_(j, r) = L_j$.

#definition(title: "Adaptive Grid Unit and Volume-Normalized Density")[
  For any non-empty subspace $S = \{j_1, ..., j_k\} subset.eq A$,
  an _adaptive grid unit_ is a Cartesian product
  $u = I_(j_1, r_1) times ... times I_(j_k, r_k)$ with
  $r_t in \{1, ..., m_(j_t)\}$. Its _normalized volume fraction_
  and _empirical support count_ are
  $
    v_S (u) = product_(t=1)^k (w_(j_t, r_t) / L_(j_t)),
    quad c(u) = |\{bold(x) in cal(X) : pi_S (bold(x)) in u\}|.
  $
  Given a _cluster dominance factor_ $alpha >= 1$, the
  volume-normalized density threshold of $u$ is
  $tau(u) = alpha n v_S (u)$, and $u$ is declared a _dense unit_
  whenever its _density enrichment ratio_ satisfies
  $
    R(u) = (c(u) / n) / v_S (u) >= alpha.
  $
  Two dense units $u, u' in cal(D)_S$ are _face-adjacent_ if they
  share identical interval indices in $k - 1$ coordinates and
  consecutive indices ($|r_j - r'_j| = 1$) in the remaining
  coordinate. Maximal connected components of face-adjacent dense
  units define the subspace clusters of $S$ @goil1999mafia.
]

Unlike CLIQUE's constant count threshold $c(u) >= tau n$,
normalizing by $v_S (u)$ prevents wide, sparse background
intervals merged during the AFI stage from masquerading as dense
clusters solely because they span a large fraction of the range.
To boost recall in high dimensions, MAFIA also relaxes CLIQUE's
strict $(k - 1)$-prefix Apriori join: two $(k - 1)$-dimensional
dense units are joined whenever they share *any* $k - 2$
identical intervals, allowing candidate generation to explore
subspaces more aggressively while parallelizing histogram scans
across $P$ processors.

#algbox(
  title: "MAFIA Adaptive-Grid Parallel Subspace Clustering",
  inputs: [Dataset $cal(X) subset bb(R)^d$ partitioned across $P$
    nodes, fine bin count $B$, merge tolerance $beta in (0, 1)$,
    dominance factor $alpha >= 1$.],
  outputs: [Maximal subspace clusters as unions of connected
    adaptive dense units.],
)[
  1. *Parallel 1D Micro-Histograms*: Each processor $p in \{1, ..., P\}$
     tallies local counts over $B$ bins per axis in
     $cal(O)((n / P) d)$ time; combine via global all-reduce.
  2. *Adaptive Finite Interval (AFI) Merging*: On each axis $j$,
     merge contiguous micro-bins into variable-width intervals
     $I_(j, r)$ whenever density variations stay within relative
     tolerance $beta$; retain 1D dense units
     $cal(D)_1 = \{I_(j, r) : R(I_(j, r)) >= alpha\}$. Set $k <- 1$.
  3. *Level-Wise Loop*: While $cal(D)_k != emptyset$:
  4. *Candidate Dense Unit (CDU) Generation*: Form
     $cal(U)_(k+1)$ by joining pairs $u_1, u_2 in cal(D)_k$
     sharing $k - 1$ identical coordinate intervals (and
     optionally verifying lower-dimensional projections).
  5. *Parallel Population & All-Reduce*: Each processor scans
     its $n / P$ local points to populate counts $c_p (u)$ for
     all $u in cal(U)_(k+1)$; aggregate global counts
     $c(u) = sum_(p=1)^P c_p (u)$ via tree all-reduce.
  6. *Dense Filtering*: Set
     $cal(D)_(k+1) <- \{u in cal(U)_(k+1) : R(u) >= alpha\}$.
     Increment $k <- k + 1$.
  7. *Connected-Component Extraction*: For each maximal subspace
     $S$ with $cal(D)_S != emptyset$, extract connected
     components of the face-adjacency graph $G_S = (cal(D)_S, E_S)$
     and remove redundant lower-dimensional projections.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five rigorous theoretical results for MAFIA: an
original $L_1$ histogram approximation bound for AFI merging, an
original counterexample proving that volume-normalized density
enrichment violates strict Apriori downward closure together with
the exact relaxed screening bound, an original analysis of
boundary-split survival relative to fixed-grid CLIQUE, the
classical parallel speedup and isoefficiency theorem, and a
graph-theoretic characterization of connected-component cluster
extraction and degenerate collapse.

#subsec("A", "Adaptive Interval Merging Accuracy and Compression")

#lemma(
  title: "Sharp L1 Distortion Bound of Adaptive Interval Merging",
  origin: "Original Proof",
)[
  On any coordinate axis $j$, let $hat(f)_B$ be the normalized
  piecewise-constant density of the $B$-bin micro-histogram and
  $hat(f)_"AFI"$ be the merged density obtained by collapsing
  contiguous micro-bins into variable-width intervals $I_(j, r)$
  whose constituent micro-bin densities satisfy
  $(1 - beta) M_r <= rho_t <= M_r$ for $M_r = max_(t in I_(j, r))
  rho_t$ and tolerance $beta in (0, 1)$. Then the global $L_1$
  histogram distortion satisfies
  $
    ||hat(f)_"AFI" - hat(f)_B||_1
    = integral_(a_j)^(b_j) |hat(f)_"AFI" (x) - hat(f)_B (x)| dif x
    <= beta / (2 (1 - beta)).
  $
]

#proof[
  Fix one merged interval $I = I_(j, r)$ spanning $s$ contiguous
  micro-bins of uniform width $h_j$ with densities
  $rho_1, ..., rho_s in [a, b]$, where $a = (1 - beta) M_r$ and
  $b = M_r$. The merged density on $I$ is the arithmetic average
  $overline(rho) = (1 / s) sum_(t=1)^s rho_t in [a, b]$. Define
  $p = (overline(rho) - a) / (b - a) in [0, 1]$. For every
  $t in \{1, ..., s\}$, write $rho_t = (1 - theta_t) a + theta_t b$
  with $theta_t = (rho_t - a) / (b - a) in [0, 1]$, noting that
  $(1 / s) sum_(t=1)^s theta_t = p$. By convexity of the function
  $z |-> |z - overline(rho)|$ on $[a, b]$,
  $
    |rho_t - overline(rho)|
    <= (1 - theta_t) |a - overline(rho)|
       + theta_t |b - overline(rho)|.
  $
  Averaging over $t = 1, ..., s$ and substituting
  $|a - overline(rho)| = p (b - a)$ and
  $|b - overline(rho)| = (1 - p)(b - a)$ yields
  $
    (1 / s) sum_(t=1)^s |rho_t - overline(rho)|
    &<= (1 - p) p (b - a) + p (1 - p)(b - a) \
    &= 2 p (1 - p)(b - a)
    <= (b - a) / 2
    = (beta / 2) M_r,
  $
  since $2 p (1 - p) <= 1 / 2$ for all $p in [0, 1]$. Integrating
  over $I$ (total length $s h_j$) and noting that
  $overline(rho) >= (1 - beta) M_r => M_r <= overline(rho) / (1 - beta)$
  and $s h_j overline(rho) = c(I) / n$, we obtain
  $
    integral_I |hat(f)_"AFI" (x) - hat(f)_B (x)| dif x
    &= h_j sum_(t=1)^s |rho_t - overline(rho)| \
    &<= s h_j dot.op (beta / 2) M_r
    <= (beta / (2 (1 - beta))) dot.op (c(I) / n).
  $
  Summing over all disjoint intervals $I_(j, 1), ..., I_(j, m_j)$
  uses $sum_r c(I_(j, r)) = n$ and yields the global bound
  $||hat(f)_"AFI" - hat(f)_B||_1 <= beta / (2 (1 - beta))$.
]

Lemma 1 establishes a clean quantitative trade-off: for small
$beta$ (e.g., $beta = 0.2$), the AFI histogram alters total
probability mass by at most $12.5\%$ while collapsing flat
background plateaus into single wide intervals, slashing the
combinatorial $k$-D candidate count by a factor of
$product_(j in S) (B / m_j)$ relative to the fine grid.

#subsec("B", "Non-Monotonicity of Volume-Normalized Density Pruning")

A subtle theoretical departure from CLIQUE @agrawal1998clique
arises from MAFIA's volume-normalized threshold $R(u) >= alpha$.
We prove that unlike raw count support, volume-normalized density
enrichment *violates* downward closure, and we derive the exact
relaxed screening threshold required for zero false negatives.

#theorem(
  title: "Failure of Apriori Monotonicity and Exact Relaxed Screen",
  origin: "Original Counterexample",
)[
  1. (_Constructive Counterexample_.) Under MAFIA's
     volume-normalized density criterion $R(u) >= alpha$ with
     $alpha > 1$, downward closure fails: there exists a 2D
     adaptive grid unit $u = I_(1, 1) times I_(2, 1)$ satisfying
     $R(u) >= alpha$ even though its 1D projection $I_(1, 1)$
     fails the threshold ($R(I_(1, 1)) < alpha$).
  2. (_Exact Necessary Screening Bound_.) For any $k$-dimensional
     adaptive unit $u = u_T times I_(j, r)$ with $|T| = k - 1$,
     $R(u) >= alpha$ necessarily implies
     $
       R(u_T)
       >= alpha dot.op (w_(j, r) / L_j)
       >= alpha dot.op min_(1 <= r' <= m_j) (w_(j, r') / L_j).
     $
     Conversely, enforcing an auxiliary fixed support floor
     $c(u) >= tau n$ restores strict downward closure for the
     support condition.
]

#proof[
  (1) Let $d = 2$, $n = 100$, $L_1 = L_2 = 1$, and $alpha = 2$.
  Partition axis $1$ into two equal intervals
  $I_(1, 1) = [0, 0.5]$ ($w_(1, 1) = 0.5$) and
  $I_(1, 2) = (0.5, 1]$ ($w_(1, 2) = 0.5$), and axis $2$ into
  $I_(2, 1) = [0, 0.1]$ ($w_(2, 1) = 0.1$) and
  $I_(2, 2) = (0.1, 1]$ ($w_(2, 2) = 0.9$). Place $40$ points in
  $u = I_(1, 1) times I_(2, 1)$ and $60$ points in
  $I_(1, 2) times I_(2, 2)$ (zero points elsewhere). Check the
  1D and 2D enrichment ratios:
  - In axis $1$: $c(I_(1, 1)) = 40$, $v_(\{1\})(I_(1, 1)) = 0.5$,
    giving $R(I_(1, 1)) = 0.40 / 0.5 = 0.8 < 2 = alpha$
    (pruned as non-dense!).
  - In axis $2$: $c(I_(2, 1)) = 40$, $v_(\{2\})(I_(2, 1)) = 0.1$,
    giving $R(I_(2, 1)) = 0.40 / 0.1 = 4.0 >= alpha$ (dense).
  - In $\{1, 2\}$: $c(u) = 40$,
    $v_(\{1, 2\})(u) = 0.5 dot.op 0.1 = 0.05$, giving
    $R(u) = 0.40 / 0.05 = 8.0 >= 2 = alpha$ (highly dense!).
  Thus strict Apriori pruning at threshold $alpha$ discards
  $I_(1, 1)$ at level $1$ and permanently misses the dense 2D
  unit $u$.

  (2) Because $u subset.eq pi_T^(-1)(u_T)$, point monotonicity
  gives $c(u_T) >= c(u)$, while volume factorizes as
  $v_S (u) = v_T (u_T) dot.op (w_(j, r) / L_j)$. Dividing yields
  $
    R(u_T)
    = (c(u_T) / n) / v_T (u_T)
    >= (c(u) / n) / (v_S (u) dot.op L_j / w_(j, r))
    = R(u) dot.op (w_(j, r) / L_j).
  $
  Substituting $R(u) >= alpha$ yields the tight lower bound.
  Finally, if a fixed count floor $c(u) >= tau n$ is imposed,
  $c(u_T) >= c(u) >= tau n$ preserves strict downward closure.
]

Theorem 1 clarifies why MAFIA's original authors paired volume
normalization with a *relaxed* candidate-joining step and
cautioned against setting $alpha$ aggressively high in 1D: under
pure $R(u) >= alpha$, bottom-up pruning is a fast heuristic
rather than an exact lattice guarantee unless the lower-level
screen is relaxed by $w_(j, "min") / L_j$.

#subsec("C", "Boundary-Split Survival Advantage Over Fixed-Grid CLIQUE")

#proposition(
  title: "Exponential-to-Linear Boundary-Loss Reduction",
  origin: "Original Derivation",
)[
  Consider a $k$-dimensional uniform hypercube cluster
  $C = product_(t=1)^k [theta_t, theta_t + Delta]$ of side $Delta$
  inside a $k$-subspace $S$.
  1. For fixed-grid CLIQUE with cell width $h = Delta$, an
     adversarial shift $theta_t = Delta / 2$ splits $C$ evenly
     across $2^k$ cells, leaving each cell with fraction $2^(-k)$
     of the cluster's points; CLIQUE misses $C$ whenever
     $tau > 2^(-k) |C| / n$.
  2. For MAFIA using fine micro-bins of width $h_"fine" = L / B
     < Delta / 2$, the AFI boundaries align with the density
     discontinuities at $theta_t$ and $theta_t + Delta$ to within
     one micro-bin $h_"fine"$ on each side, ensuring a single
     interior adaptive unit retains at least fraction
     $
       (1 - 2 h_"fine" / Delta)^k
       >= 1 - (2 k h_"fine") / Delta
     $
     of the cluster's volume (linear rather than exponential
     decay in $k$).
]

#proof[
  (1) When $theta_t = Delta / 2$ on a grid with cell boundaries
  at integer multiples of $h = Delta$, each coordinate interval
  $[Delta / 2, 3 Delta / 2]$ intersects $[0, Delta]$ and
  $[Delta, 2 Delta]$ in sub-intervals of length $Delta / 2$.
  Taking the Cartesian product over $k$ dimensions partitions
  $C$ into $2^k$ congruent sub-cubes, each of volume
  $(Delta / 2)^k = 2^(-k) Delta^k$, hence containing fraction
  $2^(-k)$ of a uniformly distributed cluster.

  (2) Along each axis $t$, the step discontinuity of the marginal
  density at $theta_t$ and $theta_t + Delta$ forces AFI boundary
  placement inside the two boundary micro-bins straddling
  $theta_t$ and $theta_t + Delta$, while all interior micro-bins
  have identical expectation and merge into one interval $I_t$ of
  length at least $Delta - 2 h_"fine"$. The $k$-dimensional
  product unit $u = product_(t=1)^k I_t$ therefore covers volume
  at least $(Delta - 2 h_"fine")^k = Delta^k (1 - 2 h_"fine" / Delta)^k$.
  Applying Bernoulli's inequality $(1 - x)^k >= 1 - k x$ (valid
  for $x = 2 h_"fine" / Delta in [0, 1]$ and $k >= 1$) completes
  the proof.
]

#subsec("D", "Parallel Scalability and Isoefficiency Analysis")

#theorem(
  title: "Shared-Nothing Parallel Speedup and Isoefficiency",
  origin: "Classical — Cited",
)[
  Let $n$ records be distributed uniformly ($n / P$ per node)
  across $P$ processors communicating via a binary-tree
  all-reduce, and let $M_"tot" = sum_(k=1)^(k_"max") |cal(U)_k|$
  denote the total number of candidate dense units evaluated
  across all levels @goil1999mafia. Assuming $M_"tot"$ fits in
  each processor's local memory, MAFIA's parallel execution time
  and efficiency satisfy
  $
    T_P = cal(O)((n / P)(d + M_"tot") + M_"tot" log P),
  $
  $
    E_P = T_1 / (P T_P)
    = (1 + Theta((P log P) / n))^(-1).
  $
  Hence MAFIA achieves asymptotically linear speedup ($E_P -> 1$)
  with isoefficiency $n = omega(P log P)$, independent of the
  candidate count $M_"tot"$ when $M_"tot" >= d$.
]

#proof[
  In serial execution ($P = 1$), constructing the 1D histograms
  costs $cal(O)(n d)$ and scanning all $n$ points against the
  $M_"tot"$ candidate units across levels costs
  $cal(O)(n M_"tot")$ (plus $cal(O)(d B + k_"max" |cal(D)|)$
  overhead dominated by $n M_"tot"$ for large $n$), giving
  $T_1 = Theta(n (d + M_"tot"))$.

  Under horizontal data partitioning across $P$ nodes, each node
  independently performs local histogram updates on its $n / P$
  points in $Theta((n / P)(d + M_"tot"))$ time with zero
  inter-node communication during data scans. At the end of each
  level $k$, aggregating the count vector of length $|cal(U)_k|$
  via a binary-tree all-reduce incurs latency and bandwidth cost
  $cal(O)(|cal(U)_k| log P)$ (plus $cal(O)(d B log P)$ at
  level $1$). Summing across levels yields $T_P$. Forming the
  efficiency ratio $E_P = T_1 / (P T_P)$ gives
  $
    E_P
    = (n (d + M_"tot")) /
      (n (d + M_"tot") + Theta(P M_"tot" log P))
    = (1 + Theta((P log P) / n))^(-1),
  $
  since $M_"tot" / (d + M_"tot") in [1 / 2, 1]$ whenever
  $M_"tot" >= d$. Maintaining constant efficiency $E_P >= E_0$
  requires $n = Omega(P log P)$, proving near-linear scalability
  for large databases.
]

#subsec("E", "Face-Adjacency Cluster Graph and Degenerate Collapse")

#proposition(
  title: "Linear-Time Graph Extraction and Coarse-Tolerance Collapse",
  origin: "Original Proof",
)[
  For any subspace $S$ of dimension $k = |S|$ with dense-unit
  set $cal(D)_S$, the undirected face-adjacency graph
  $G_S = (cal(D)_S, E_S)$ satisfies:
  1. Every vertex $u in cal(D)_S$ has degree at most $2 k$, so
     connected-component extraction via BFS/DFS runs in
     $cal(O)(k |cal(D)_S|)$ time and outputs a unique partition
     independent of unit traversal order;
  2. If the relative merge tolerance satisfies
     $beta >= 1 - rho_"min" / rho_"max"$ along every axis
     (equivalently, marginal densities vary by less than a factor
     $(1 - beta)^(-1)$), every axis collapses to $m_j = 1$ single
     interval spanning $[a_j, b_j]$ with $R(u) = 1$, forcing
     MAFIA to return a degenerate partition ($0$ clusters if
     $alpha > 1$, or $1$ trivial all-point cluster if $alpha = 1$).
]

#proof[
  (1) Represent each unit $u in cal(D)_S$ by its integer index
  tuple $(r_1, ..., r_k)$. A face-adjacent neighbor $u'$ differs
  from $u$ by $plus.minus 1$ in exactly one coordinate
  $t in \{1, ..., k\}$ and matches in all other $k - 1$
  coordinates. There are at most $2 k$ such tuples, hence
  $"deg"(u) <= 2 k$ and $|E_S| <= k |cal(D)_S|$. Storing
  $cal(D)_S$ in a hash map permits $cal(O)(k)$ expected neighbor
  lookup per vertex, yielding $cal(O)(k |cal(D)_S|)$ total
  BFS/DFS complexity. Uniqueness follows from the canonical
  equivalence-class decomposition of undirected graph components.

  (2) If $rho_"min" >= (1 - beta) rho_"max"$ on axis $j$, then
  every pair of micro-bins satisfies the merge condition, so the
  left-to-right AFI pass merges all $B$ bins into a single
  interval $I_(j, 1) = [a_j, b_j]$ with $w_(j, 1) = L_j$. For
  any subspace $S$, the sole grid unit $u = product_(j in S)
  [a_j, b_j]$ has $c(u) = n$ and $v_S (u) = 1$, giving
  $R(u) = 1$. Thus $cal(D)_S = emptyset$ whenever $alpha > 1$,
  and $cal(D)_S = \{u\}$ (one cluster containing all $n$ points)
  whenever $alpha = 1$, rendering Silhouette undefined (`N/A`).
]

#sec("III", "Complexity, Geometry and Operational Profile")

MAFIA's serial runtime is
$cal(O)(n d + d B + n M_"tot" + k_"max" |cal(D)|)$, where
$M_"tot" = sum_(k=1)^(k_"max") |cal(U)_k|$ is the total count of
evaluated candidate dense units. Because adaptive merging keeps
$m_j << B$ (often $m_j in [2, 8]$ on smooth marginals), MAFIA
dramatically reduces $M_"tot"$ compared to fine fixed grids
while avoiding SUBCLU's $cal(O)(n^2)$ pairwise neighborhood
scans and ORCLUS's $cal(O)(d^3)$ eigendecompositions. When dense
clusters persist into high subspaces ($k_"max" >= 15$), however,
$M_"tot"$ still scales exponentially in $k_"max"$ as
$cal(O)(c^(k_"max"))$ for average branching factor $c$
@goil1999mafia.

Memory complexity is $cal(O)(n d + d B + max_k |cal(U)_k|)$ (or
$cal(O)((n / P) d + d B + max_k |cal(U)_k|)$ per processor in
distributed mode). Because 1D histograms and CDU lookups rely
exclusively on contiguous array indexing and binary search over
$m_j$ interval endpoints, MAFIA enjoys exceptional CPU cache
locality and vectorization throughput, explaining its ultra-low
empirical runtime.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time complexity],
    [$cal(O)(n (d + M_"tot"))$ serial;
    $cal(O)((n / P)(d + M_"tot") + M_"tot" log P)$ parallel;
    worst-case $cal(O)(c^(k_"max") + n d)$.],
    [Space complexity],
    [$cal(O)((n / P) d + d B + max_k |cal(U)_k|)$ per node;
    avoids $cal(O)(n^2)$ distance storage.],
    [Cluster geometry],
    [Unions of face-adjacent variable-width axis-aligned
    hyper-rectangles in maximal subspaces.],
    [Noise & overlap],
    [Explicit background noise (non-dense cells); overlapping
    clusters permitted across distinct subspaces.],
    [Invariance],
    [Invariant to translations, axis permutations, and
    per-axis positive affine scaling; not rotation-invariant.],
    [Key hyperparameters],
    [Fine bin count $B$, merge tolerance $beta in (0, 1)$,
    cluster dominance factor $alpha >= 1$.],
  ),
  caption: [Operational, geometric, and complexity profile of MAFIA.],
)

Notice a subtle invariance advantage over DOC and SUBCLU: because
MAFIA builds micro-bins over each axis's own span $[a_j, b_j]$
and normalizes volumes by $L_j$, MAFIA is *strictly invariant*
to independent positive affine rescaling of individual features
($x_j |-> s_j x_j + t_j$, $s_j > 0$), eliminating the need for
prior z-score standardization!

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Distributed Scientific Data-Warehouse Mining*: Massive
   petascale climate, genomics, and particle-physics archives
   stored across shared-nothing compute clusters require linear
   horizontal scaling in $n$. MAFIA's tree all-reduce histogram
   architecture (Theorem 2) achieves near-linear speedup across
   dozens to hundreds of nodes without shuffling raw records.
2. *Multi-Unit Industrial Telemetry with Heterogeneous Scales*:
   Process plants log temperatures, pressures, voltages, and flow
   rates spanning six orders of magnitude. Because MAFIA's
   per-axis AFI normalization is invariant to arbitrary positive
   coordinate rescaling, it extracts interpretable multi-variable
   operating boxes directly in native engineering units.
3. *High-Frequency Event and Clickstream Cube Aggregation*:
   Online analytical processing (OLAP) data cubes naturally store
   variable-width range buckets along dimensions such as latency,
   session duration, and revenue. MAFIA's adaptive intervals map
   directly onto OLAP range-partitioned histograms, isolating
   dense multi-attribute traffic regimes in sub-millisecond scans.
4. *Geospatial and Environmental Multi-Layer Grid Analysis*:
   Remote-sensing rasters and ecological surveys combine spatial
   coordinates with elevation, moisture, and spectral indices
   whose marginals exhibit sharp multimodal cliffs (e.g., land
   vs. water). Proposition 1 ensures AFI boundaries snap directly
   to physical transition cliffs without splitting dense biomes.
5. *Fast Exploratory Pre-Screening for High-Volume Pipelines*:
   When exploratory pipelines must inspect tens of millions of
   records interactively to identify candidate subspaces before
   expensive downstream modeling (e.g., ORCLUS or GMM), MAFIA's
   $cal(O)(n)$-linear histogram pass rapidly narrows the
   attribute lattice to high-enrichment subspaces.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

MAFIA breaks down in four regimes: (i) _Coarse-tolerance
degenerate collapse_ (Proposition 2(2)), where overly large
$beta$ or aggressive $alpha$ merges unimodal or smoothly varying
marginals into $0$ or $1$ cluster (`N/A` Silhouette); (ii)
_Volume-normalized false negatives_ (Theorem 1(1)), where 1D
pruning discards genuinely dense high-dimensional cells whose
marginal enrichment falls below $alpha$; (iii) _Curved and
diagonal geometries_, because axis-aligned rectangular unions
approximate oblique manifolds poorly with staircase artifacts;
and (iv) _High intrinsic cluster dimensionality_ ($k_"max" > 15$),
where bottom-up CDU growth explodes combinatorially.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Massive $n$, multimodal marginals],
    [Optimal: $cal(O)(n / P)$ parallel scaling; snaps to cliffs.],
    [Avoid high $beta$ to prevent Prop. 2(2) single-bin collapse.],
    [Mixed-unit tabular features],
    [Strong: exact per-axis positive affine scale invariance.],
    [Combine with support floor $tau$ to curb Thm. 1(1) misses.],
    [Smooth / unimodal marginals],
    [Poor: AFI merges entire axis into 1 interval (degenerate).],
    [Use CLIQUE @agrawal1998clique, ENCLUS @cheng1999enclus, or
    OPTIGRID @hinneburg1999optigrid.],
    [Curved / non-convex manifolds],
    [Weak: staircase grid cells fragment thin curved shapes.],
    [Use SUBCLU @kailing2004subclu (Alg. 137).],
    [Oblique correlated hyperplanes],
    [Poor: axis-aligned boxes smear diagonal correlations.],
    [Use ORCLUS @aggarwal2000orclus (Alg. 138).],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for MAFIA.],
)

#sec("V", "Structural Relations and Empirical Behavior")

MAFIA anchors the adaptive-grid branch of Cat. VI alongside
fixed-grid CLIQUE @agrawal1998clique, entropy-based ENCLUS
@cheng1999enclus, wave-transformed WaveCluster
@sheikholeslami1998wavecluster, and optimal hyperplane-cutting
OPTIGRID @hinneburg1999optigrid. The following proposition
establishes the exact parametric limit reducing MAFIA to CLIQUE.

#proposition(
  title: "Exact Reduction of MAFIA to Fixed-Grid CLIQUE",
  origin: "Original Proof",
)[
  Suppose the AFI merging tolerance is set to $beta -> 0^+$ with
  fine bin count $B = xi$ (so no adjacent bins merge, yielding
  $m_j = xi$ uniform intervals of width $w_(j, r) = L_j / xi$ on
  every axis $j$). Then for every $k$-dimensional grid unit $u$,
  $v_S (u) = xi^(-k)$ is constant across all units of level $k$,
  and MAFIA's dense-unit predicate $R(u) >= alpha$ reduces
  identically to CLIQUE's fixed-count predicate
  $c(u) >= tau_k n$ with level-calibrated selectivity
  $tau_k = alpha xi^(-k)$, matching CLIQUE @agrawal1998clique
  exactly at every single fixed level $k$.
]

#proof[
  When $beta = 0$ and empirical bin densities are distinct (or
  merging is disabled), all $B = xi$ micro-bins remain unmerged,
  so $w_(j, r) / L_j = 1 / xi$ for every $j in A$ and
  $r in \{1, ..., xi\}$. Thus every $k$-dimensional unit $u$ has
  identical normalized volume $v_S (u) = product_(t=1)^k (1 / xi)
  = xi^(-k)$. Substituting into Definition 1 gives
  $R(u) >= alpha <=> c(u) >= alpha n xi^(-k) = tau_k n$. At fixed
  level $k$, the dense unit sets $cal(D)_S$ and their
  face-adjacency graphs $G_S$ coincide identically with CLIQUE.
  Note that CLIQUE keeps $tau$ constant across $k$ (guaranteeing
  strict downward closure), whereas MAFIA's $tau_k = alpha xi^(-k)$
  shrinks geometrically with $k$ to counteract dimensional
  dilution—the exact mechanism behind Theorem 1(1)'s loss of
  strict monotonicity.
]

Empirical results on the ClusterOverview 4-dataset benchmark
(@rousseeuw1987silhouettes, @hubert1985comparing) dramatically
illustrate both MAFIA's speed-precision upside and its
degenerate-collapse vulnerability (Proposition 2(2)). On *Iris*,
MAFIA achieves the highest Silhouette score of the entire quartet
($"Sil" = 0.642$, far ahead of DOC $0.350$, SUBCLU $0.315$, and
ORCLUS $0.292$), paired with $"ARI" = 0.388$: adaptive histogram
cliffs cleanly snap to the well-separated bimodal gap between
_Setosa_ and the combined _Versicolor_ / _Virginica_ cloud,
yielding two high-Silhouette macro-clusters (hence high Sil but
lower 3-class ARI than ORCLUS $0.681$ and SUBCLU $0.516$). On
*Customers* ($"Sil" = "N/A"$) and *Multi-Shape*
($"Sil" = "N/A"$, $"ARI" = "N/A"$), however, smooth unimodal
marginal projections trigger Proposition 2(2) coarse-tolerance
collapse—producing fewer than two clusters so Silhouette cannot
be evaluated—highlighting MAFIA's sensitivity to $(beta, alpha)$
on non-multimodal marginals. On *High-D*, where several axes
retain multimodal density contrast, MAFIA recovers valid
structure ($"Sil" = 0.177$, outperforming DOC $0.068$ while
trailing ORCLUS $0.254$ and SUBCLU $0.212$). Most notably, MAFIA
achieves a blistering mean runtime of $0.058 "ms"$—nearly 2x
faster than DOC ($0.112 "ms"$), 48x faster than ORCLUS
($2.792 "ms"$), and 214x faster than SUBCLU ($12.412 "ms"$)—making
it the clear throughput leader whenever multimodal marginals
support reliable adaptive interval discovery.
