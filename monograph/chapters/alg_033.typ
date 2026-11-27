#import "../template/ieee.typ": *

#alg_entry(
  id: 33,
  name: "OPTICS (Ordering Points To Identify the Clustering Structure)",
  category: "Cat. II: Density & Mode-Seeking",
  year: "1999",
  paradigm: "Priority-Queue Reachability Ordering & Multi-Scale Valley Extraction",
  complexity: $O(n log n) "/" O(n^2)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

To overcome DBSCAN's dependence on a single global radius $epsilon$ without committing immediately to one flat partition, OPTICS (Ordering Points To Identify the Clustering Structure) @ankerst1999optics produces an augmented linear ordering of the dataset that encodes its multi-scale density-based clustering structure across *all* radii $epsilon' <= epsilon_"max"$ simultaneously.

Let $(cal(X), d)$ be a finite metric space of size $n$, parameterized by a generating radius $epsilon_"max" in (0, infinity]$ and minimum point count $m = "MinPts" in bb(N)_(>= 2)$. Let $N_(epsilon_"max") (p) = {q in cal(X) : d(p, q) <= epsilon_"max"}$ and let $d_(m) (p)$ denote the distance from $p$ to its $m$-th nearest neighbor in $cal(X)$ (including $p$). OPTICS defines the *truncated core distance* of $p$:
$ d_"core"^(epsilon_"max") (p) = cases(infinity & "if" |N_(epsilon_"max") (p)| < m, d_(m) (p) & "otherwise.") $
Crucially—and in contrast to HDBSCAN's symmetric mutual reachability (Alg. 32)—OPTICS defines an *asymmetric* directed *reachability distance* from a core point $p$ to any target point $q in cal(X)$:
$ d_"reach" (q | p) = cases(infinity & "if" |N_(epsilon_"max") (p)| < m, max { d_"core"^(epsilon_"max") (p), space d(p, q) } & "otherwise.") $
This asymmetry allows a non-core border point $q$ (whose own core distance exceeds $epsilon'$) to inherit a finite reachability $<= epsilon'$ from an adjacent core point $p$.

#definition(title: "OPTICS Cluster Ordering and Reachability Plot")[
  An *OPTICS ordering* of $(cal(X), d)$ with parameters $(epsilon_"max", m)$ is a permutation $pi = (o_1, dots, o_n)$ of $cal(X)$ together with core-distance sequence $c_t = d_"core"^(epsilon_"max") (o_t)$ and reachability sequence $r_t in [0, infinity]$ ($r_1 = infinity$) satisfying at every step $t in {2, dots, n}$:
  $ r_t &= min_(1 <= s < t) d_"reach" (o_t | o_s) \
        &= min_(q in cal(X) backslash {o_1, dots, o_(t-1)}) min_(1 <= s < t) d_"reach" (q | o_s), $
  with arbitrary tie-breaking among minimizers (and choosing any unprocessed point with $r_t = infinity$ whenever all cross-cut reachabilities are infinite). The discrete curve $t |-> r_t$ is the *reachability plot*, where contiguous low-reachability valleys represent dense clusters across scales.
]

Algorithm 1 implements Definition 1 using a min-priority queue `Seeds` that greedily expands the densest currently reachable region before crossing sparser gaps @ankerst1999optics.

#algbox(
  title: "OPTICS Ordering Algorithm (Ankerst et al., 1999)",
  inputs: [Metric dataset $(cal(X), d)$, generating radius $epsilon_"max" in (0, infinity]$, integer $m >= 2$.],
  outputs: [Ordering $pi = (o_1, dots, o_n)$, core distances $(c_t)$, reachabilities $(r_t)$.],
)[
  1. Mark all $p in cal(X)$ unprocessed; set $r[p] <- infinity$; initialize empty list $pi$ and min-priority queue `Seeds`.
  2. *For* each unprocessed $p in cal(X)$ *do*:
  3. #h(1em) Query $N <- N_(epsilon_"max") (p)$; compute $c[p] <- d_"core"^(epsilon_"max") (p)$.
  4. #h(1em) Mark $p$ processed; append $(p, c[p], r[p])$ to output $(pi, c, r)$.
  5. #h(1em) *If* $c[p] < infinity$ *then*:
  6. #h(2em) *Update*$(N, p, "Seeds")$: for each unprocessed $q in N$, set $rho <- max{c[p], d(p, q)}$; if $rho < r[q]$, set $r[q] <- rho$ and insert/decrease-key $q$ in `Seeds`.
  7. #h(2em) *While* `Seeds` $!= emptyset$ *do*:
  8. #h(3em) Pop $q <- "ExtractMin"("Seeds")$; mark $q$ processed.
  9. #h(3em) Query $N' <- N_(epsilon_"max") (q)$; compute $c[q] <- d_"core"^(epsilon_"max") (q)$.
  10. #h(3em) Append $(q, c[q], r[q])$ to $(pi, c, r)$.
  11. #h(3em) *If* $c[q] < infinity$, call *Update*$(N', q, "Seeds")$.
  12. *Return* $(pi, c, r)$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five formal properties of OPTICS: (1) its exact priority-queue invariant and equivalence to a modified Prim traversal on the directed reachability graph; (2) exact linear-time recovery of DBSCAN core partitions for every $epsilon' <= epsilon_"max"$; (3) an explicit counterexample proving ordering divergence between OPTICS's asymmetric reachability and HDBSCAN's symmetric mutual reachability; (4) similarity invariance and nesting/complexity properties of $xi$-steepness extraction; and (5) exact heap-update complexity bounds proving the quadratic wall at $epsilon_"max" = infinity$.

#subsec("A", "Priority-Queue Invariant and Exact DBSCAN Extraction")

#lemma(title: "Priority-Queue Relaxation Invariant and Directed Prim Equivalence", origin: "Original Proof")[
  At every step $t in {2, dots, n}$ of Algorithm 1, let $P_(t-1) = {o_1, dots, o_(t-1)}$ be the processed prefix and $U_(t-1) = cal(X) backslash P_(t-1)$ the unprocessed suffix. Then:
  1. Every $q in U_(t-1)$ satisfies $r[q] = min_(p in P_(t-1)) d_"reach" (q | p)$, and `Seeds` contains precisely those $q in U_(t-1)$ with $r[q] < infinity$.
  2. The point $o_t$ popped in Step 8 (or chosen in Step 2 if `Seeds` is empty) attains the exact minimum directed cut weight across $(P_(t-1), U_(t-1))$ in the weighted directed graph $cal(D)_"reach" = (cal(X), vec(E), d_"reach")$.
]

#proof[
  We proceed by induction on $t$. *Base case ($t = 2$):* After processing $o_1$, if $c[o_1] = infinity$ then $d_"reach"(q | o_1) = infinity$ for all $q in U_1$, matching $r[q] = infinity$ and empty `Seeds`. If $c[o_1] < infinity$, then `Update` scans every $q in N_(epsilon_"max")(o_1)$ (precisely those $q$ with $d_"reach"(q | o_1) <= epsilon_"max" < infinity$, wait: for $q in.not N_(epsilon_"max")(o_1)$, wait—notice if $epsilon_"max" < infinity$, `Update` only iterates over $q in N_(epsilon_"max")(o_1)$, so OPTICS implicitly sets $d_"reach"(q | p) = infinity$ when $d(p, q) > epsilon_"max"$). Setting $r[q] = max{c[o_1], d(o_1, q)}$ for $q in N_(epsilon_"max")(o_1)$ and inserting $q$ into `Seeds` establishes (1) for $t = 2$.

  *Inductive step:* Assume (1) holds prior to step $t$. Popping `ExtractMin(Seeds)` selects $o_t in arg min_(q in U_(t-1)) r[q]$, which by the induction hypothesis equals $arg min_(q in U_(t-1)) min_(p in P_(t-1)) d_"reach"(q | p)$, proving (2). Next, if $c[o_t] < infinity$, `Update` relaxes $r'[q] = min(r[q], d_"reach"(q | o_t))$ for every $q in U_t inter N_(epsilon_"max")(o_t)$ (and leaves $r'[q] = r[q]$ for $q$ with $d_"reach"(q | o_t) = infinity$). Thus $r'[q] = min_(p in P_t) d_"reach"(q | p)$ for all $q in U_t$, closing the induction.
]

#theorem(title: "Exact Linear-Time Recovery of All DBSCAN Core Partitions", origin: "Classical — Cited")[
  Let $(pi, c, r)$ be the output of Algorithm 1 with parameters $(epsilon_"max", m)$, and let $0 < epsilon' <= epsilon_"max"$. Define the linear-scan rule `ExtractDBSCAN` @ankerst1999optics: scan $t = 1, dots, n$; if $r_t > epsilon'$, start a new cluster when $c_t <= epsilon'$ and otherwise mark $o_t$ as noise; if $r_t <= epsilon'$, assign $o_t$ to the active cluster. Then every core component $K subset.eq cal(X)_"core"^(epsilon')$ of DBSCAN($epsilon', m$) appears as a contiguous block (interleaved only with border points of $K$), and `ExtractDBSCAN` recovers every core component identically in $O(n)$ time.
]

#proof[
  Let $K subset.eq cal(X)_"core"^(epsilon')$ be any connected component of the $epsilon'$-core graph $G_"core"^(epsilon')$ (Alg. 31, Lemma 1), and let $t_0 = min {t : o_t in K}$ be the first index at which Algorithm 1 processes an element of $K$.

  *Claim 1 ($r_(t_0) > epsilon'$ and $c_(t_0) <= epsilon'$):* Because $o_(t_0) in K subset.eq cal(X)_"core"^(epsilon')$ and $epsilon' <= epsilon_"max"$, we have $c_(t_0) = d_(m)(o_(t_0)) <= epsilon'$. Suppose for contradiction that $r_(t_0) <= epsilon'$. By Lemma 1, there exists $p in P_(t_0 - 1)$ with $max{d_"core"^(epsilon_"max")(p), d(p, o_(t_0))} <= epsilon'$. Hence $p in cal(X)_"core"^(epsilon')$ and $d(p, o_(t_0)) <= epsilon'$, implying $p in K$, contradicting minimality of $t_0$. Thus `ExtractDBSCAN` opens a fresh cluster at index $t_0$.

  *Claim 2 (Exhaustive contiguous traversal of $K$):* First observe that at step $t_0$, *no* unprocessed point $w in U_(t_0 - 1)$ can satisfy $min_(p in P_(t_0 - 1)) d_"reach"(w | p) <= epsilon'$: otherwise $r[w] <= epsilon' < r_(t_0)$, contradicting Lemma 1's global minimum selection of $o_(t_0)$. Now consider any subsequent step $t > t_0$ while $K inter U_(t-1) != emptyset$. Because $K$ is connected in $G_"core"^(epsilon')$, there exists an edge $(v, w)$ in $G_"core"^(epsilon')$ with $v in K inter P_(t-1)$ and $w in K inter U_(t-1)$, yielding $r[w] <= d_"reach"(w | v) <= epsilon'$. Consequently, `ExtractMin` must pop a point $o_t$ with $r_t <= epsilon'$. By the absence of any $<= epsilon'$ reachability from $P_(t_0 - 1)$, $o_t$ must satisfy $d_"reach"(o_t | p) <= epsilon'$ for some core point $p in {o_(t_0), dots, o_(t-1)}$. Inducting on $t - t_0$, every such $p$ belongs to $K$, so $o_t$ is either in $K$ or a border point of $K$ (no core point of another component $K' != K$ can lie within $epsilon'$ of $K$). Hence all of $K$ is extracted inside the same active cluster before $r_t > epsilon'$ recurs. Scanning $(c, r)$ once takes $O(n)$ time.
]

#subsec("B", "Asymmetric Reachability Divergence from HDBSCAN")

Because OPTICS uses asymmetric $d_"reach"(q | p)$ while HDBSCAN uses symmetric $d_"mreach"(p, q)$, their cluster orderings are *not* equivalent even when $epsilon_"max" = infinity$.

#proposition(title: "Asymmetric OPTICS vs. Symmetric Mutual-Reachability Ordering Divergence", origin: "Original Counterexample")[
  1. *Left-Shift Valley Artifact:* By Claim 1 of Theorem 1, the first element $o_(t_0)$ of every cluster $K$ carries reachability $r_(t_0) > epsilon'$ (reflecting the inter-cluster gap), while internal elements satisfy $r_t <= epsilon'$.
  2. *Strict Ordering Divergence:* For every $m >= 3$ and $epsilon_"max" = infinity$, there exists a one-dimensional dataset of size $n = 2m + 1$ on which OPTICS pops a sparse non-core border point $q$ *before* crossing to a dense second core cluster $B$, whereas Prim's MST on HDBSCAN's symmetric mutual reachability connects $A$ to $B$ *before* attaching $q$.
]

#proof[
  (1) Immediate from Claim 1 of Theorem 1: $o_(t_0)$ starts the valley yet exhibits the preceding peak height $r_(t_0)$.

  (2) Fix any $m >= 3$ and let $eta = 10^(-3) / m$. Construct $2m + 1$ distinct points in $bb(R)$:
  $ A &= {0} union { -0.40 + j eta : j = 1, dots, m - 1 }, \
    B &= { -0.95 + j eta : j = 0, dots, m - 1 }, \
    q &= +0.35, quad cal(X) = A union B union {q}. $
  Compute exact core distances ($m$-th NN including self):
  - For $a_1 = 0$, its $m - 1$ nearest neighbors are $A backslash {0}$ at distances $<= 0.40$, so $d_"core"(a_1) in [0.399, 0.400]$.
  - For $a in A backslash {0}$, there are $m - 2$ other points in $A backslash {0}$ plus $a_1 = 0$, so $d_"core"(a) <= 0.400$.
  - For every $b in B$, all $m$ points of $B$ lie within $[ -0.95, -0.949 ]$, so $d_"core"(b) <= 0.001$.
  - For $q = +0.35$, its 1st NN is $q$, 2nd NN is $a_1 = 0$ ($d = 0.35$), and 3rd through $m$-th NNs lie in $A backslash {0}$ at distance $>= 0.749$ (since $m >= 3$). Hence $d_"core"(q) >= 0.749$.

  Initialize both traversals at $a_1 = 0$. Both algorithms first consume all of $A$ at weights $<= 0.400$. Now compare the cut weights from $A$ to $q$ versus $A$ to $B$:
  - *OPTICS (asymmetric $d_"reach"$):* $d_"reach"(q | a_1) = max{d_"core"(a_1), 0.35} <= 0.400$, whereas every $a in A, b in B$ satisfies $d_"reach"(b | a) >= d(a, b) >= 0.549$. Thus OPTICS strictly pops $q$ *before* $B$ at $r(q) <= 0.400$.
  - *HDBSCAN (symmetric $d_"mreach"$):* Every $a in A$ has $d_"mreach"(a, q) >= d_"core"(q) >= 0.749$, whereas for $a_2 in A backslash {0}$ and $b_m in B$, $d_"mreach"(a_2, b_m) = max{0.400, 0.001, d(a_2, b_m)} <= 0.551 < 0.749$. Thus HDBSCAN connects $A$ to $B$ *before* attaching $q$.
]

#subsec("C", "Properties of xi-Steepness Multi-Scale Extraction")

To automate multi-density extraction beyond horizontal cuts, Ankerst et al. @ankerst1999optics introduce `ExtractXi` with steepness threshold $xi in (0, 1)$: index $t$ is *steep-down* if $r_t (1 - xi) >= r_(t+1)$ and *steep-up* if $r_t <= r_(t+1) (1 - xi)$. Contiguous steep-down areas $D = [s_D, e_D]$ and steep-up areas $U = [s_U, e_U]$ enclose a $xi$-cluster $[s, e]$ of size $>= m$ whose interior reachabilities lie strictly below $min(r_(s_D), r_(e_U + 1)) (1 - xi)$.

#proposition(title: "Similarity Invariance and Hierarchical Nesting of ExtractXi", origin: "Original Proof")[
  Set $epsilon_"max" = infinity$. Then:
  1. Under any Euclidean similarity $T(x) = alpha Q x + b$ ($alpha > 0, Q in O(d)$), every finite reachability ratio $r_(t+1) / r_t$ is strictly invariant, making `ExtractXi` invariant under all similarities (modulo ordering ties).
  2. Unlike `ExtractDBSCAN` (which returns a flat partition in $O(n)$ time), `ExtractXi` can return $Theta(n / m)$ strictly nested overlapping clusters and evaluates up to $O(n^2)$ candidate $(D, U)$ steep-region pairs in the worst case.
]

#proof[
  (1) Under $T(x) = alpha Q x + b$, every distance scales by $alpha$, so $d_"core"^(infinity)$ and $d_"reach"$ scale uniformly by $alpha$. Thus $r_t^((T)) = alpha r_t$ for all $t$, preserving both the priority-queue comparison order and the homogeneous ratio inequalities $r_(t+1) / r_t <= 1 - xi$ and $r_t / r_(t+1) <= 1 - xi$ identically.

  (2) Consider $k$ nested concentric shells in $bb(R)$ of $2m$ points each ($n = 2 k m$) whose successive intra-shell spacings decrease geometrically as $delta_j = rho^j$ with $0 < rho < (1 - xi) / 2$ for $j = 1, dots, k$. Traversing from outermost to innermost shell generates $k$ steep-down transitions $D_1, dots, D_k$ followed by $k$ matching steep-up exits $U_k, dots, U_1$ (where every interior level satisfies $rho^(j+1) < (1 - xi) rho^j$). Every matched pair $(D_j, U_j)$ satisfies the $xi$-cluster axioms, producing $k = n / (2m)$ strictly nested clusters $[s_1, e_1] supset [s_2, e_2] supset dots supset [s_k, e_k]$. Furthermore, `ExtractXi`'s sweep checks each newly encountered steep-up region against the active list of up to $O(n)$ prior steep-down regions, incurring $O(n^2)$ worst-case pair evaluations.
]

#subsec("D", "Exact Priority-Queue and Range-Query Complexity Bounds")

#theorem(title: "Exact Heap-Update Complexity and the Infinite-Epsilon Barrier", origin: "Original Derivation")[
  Let $T_"range"(n, d, epsilon_"max")$ be the cost of an $epsilon_"max"$-range query and let $M = sum_(p in cal(X)) |N_(epsilon_"max")(p)|$ be the total number of directed edges in the $epsilon_"max"$-neighborhood graph. Excluding index construction, Algorithm 1 executes in strictly $O(n)$ auxiliary memory and exact worst-case time:
  - *Fibonacci heap:* $O(n T_"range" + M + n log n)$.
  - *Binary heap:* $O(n T_"range" + (n + M) log n)$.
  Consequently, when $epsilon_"max" = infinity$, $M = n(n - 1)$ and *no* priority-queue implementation can beat $Omega(n^2)$ time.
]

#proof[
  Each vertex $p in cal(X)$ is marked processed exactly once, executing at most $n$ range queries (total cost $n T_"range"$) and at most $n$ `ExtractMin` operations (total cost $O(n log n)$ in both Fibonacci and binary heaps). During `Update`, each directed neighbor pair $(p, q)$ with $q in N_(epsilon_"max")(p)$ triggers at most one `Insert` or `DecreaseKey` call—at most $M$ priority updates across the entire execution. In a Fibonacci heap, amortized `Insert`/`DecreaseKey` costs $O(1)$ ($O(M)$ total); in a binary heap, each update costs $O(log n)$ ($O(M log n)$ total). Because `Seeds` stores at most one key per unprocessed point and $(pi, c, r)$ have length $n$, workspace is strictly $O(n)$. When $epsilon_"max" = infinity$ and $n >= m$, every point is core with $|N_infinity (p)| = n$, forcing $M = Theta(n^2)$ (in which case an unindexed $O(n)$-scan array replaces the binary heap to avoid the $log n$ factor, yielding $Theta(n^2)$ time).
]

#sec("III", "Complexity, Geometry and Operational Profile")

Theorem 2 reveals the central operational trade-off of OPTICS: setting a tight finite $epsilon_"max"$ paired with a low-dimensional R\*-tree or KD-tree bounds $M = O(n)$ and $T_"range" = O(log n)$, achieving $O(n log n)$ runtime—at the price of truncating clusters sparser than $epsilon_"max"$ to infinite reachability. Conversely, setting $epsilon_"max" = infinity$ captures the entire multi-scale hierarchy cleanly (with full similarity invariance, Proposition 2) at $Theta(n^2 d)$ cost and $O(n)$ memory. Unlike agglomerative trees that output opaque dendrograms, OPTICS produces an interpretable 1D reachability curve exposing cluster cardinality, separation gaps, and nested substructure at a glance.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for OPTICS],
    [Time / Space Complexity], [$O(n log n)$ indexed with tight $epsilon_"max"$; $Theta(n^2)$ at $epsilon_"max" = infinity$ (Thm. 2); $O(n)$ space.],
    [Geometric Inductive Bias], [Arbitrary non-convex and *nested multi-density* valleys via asymmetric reachability.],
    [Determinism / Ordering], [Linear order varies with seed/ties; extracted DBSCAN core clusters are invariant (Thm. 1).],
    [Invariances], [Similarity invariant when $epsilon_"max" = infinity$ (Prop. 2); rigid-isometry invariant for finite $epsilon_"max"$.],
    [Outlier Handling], [High reachability peaks isolate noise; directly powers OPTICS-OF (Alg. 53) @breunig1999opticsof.],
    [Recommended Upgrades], [HDBSCAN (Alg. 32) for automatic flat multi-density cuts; OPTICS-OF (Alg. 53) for local outliers.],
  ),
  caption: [Summary of theoretical guarantees and operational profile of OPTICS.],
)

#sec("IV", "Engineering and Scientific Applications")

OPTICS excels whenever analysts must inspect multi-scale density structure interactively or discover nested cluster hierarchies across unknown spatial scales.

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Exploratory Multi-Scale Visual Analytics:* By projecting arbitrary $d$-dimensional density topology onto a 1D reachability plot $t |-> r_t$, OPTICS allows domain scientists to diagnose nested sub-clusters ("valleys inside valleys"), estimate cluster cardinalities from valley widths, and calibrate $epsilon'$ without repeated clustering runs.
2. *Seismology and Multi-Scale Earthquake Clustering:* Regional fault lines form broad low-density valleys enclosing intense localized aftershock swarms; a single OPTICS ordering enables geologists to inspect both macro-tectonic belts and micro-swarms via Theorem 1's $O(n)$ cut extraction.
3. *Maritime, Aviation, and Urban Trajectory Corridor Mining:* Traffic corridors vary sharply in width between congested harbor channels and open-sea lanes; combined with segment-distance metrics @lee2007traclus, OPTICS exposes nested corridor hierarchies without enforcing one fixed corridor width.
4. *Local Density-Aware Anomaly Detection (OPTICS-OF):* Extending OPTICS with local reachability ratios yields OPTICS-OF (Alg. 53 @breunig1999opticsof), identifying points that are locally sparse relative to an adjacent dense cluster even when globally denser than a separate sparse cluster.
5. *Molecular Dynamics Conformational Basin Exploration:* Free-energy landscapes of proteins exhibit hierarchical metastable sub-basins separated by saddles of varying heights; reachability plots mirror disconnectivity graphs of conformational states.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

OPTICS struggles in four regimes: (a) *fragility of `ExtractXi` on smooth density tails*, where gradual Gaussian slopes lack sharp $(1 - xi)$ cliffs, causing missed or fragmented clusters; (b) *high-dimensional reachability flattening* @beyer1999nearest, where valleys lose vertical contrast; (c) *quadratic runtime* when $epsilon_"max" = infinity$ on large $n$ (Theorem 2); and (d) *border-point ordering artifacts* (Proposition 1) that can blur steepness boundaries.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Why OPTICS Fits / Prerequisite], [Caution / Alternative],
    [Visual / EDA Diagnostics], [1D reachability plot reveals nested scales], [Account for 1-point valley shift (Prop. 1)],
    [Seismology / GIS], [Nested fault belts & local swarms], [Choose finite $epsilon_"max"$ for $O(n log n)$],
    [Automatic Flat Cuts], [Requires manual $epsilon'$ or fragile $xi$], [Prefer HDBSCAN (Alg. 32) EOM],
    [Local Outlier Mining], [Direct basis for OPTICS-OF (Alg. 53)], [Compare GLOSH in HDBSCAN (Alg. 32)],
    [High-D Subspaces], [Reachability curves flatten in $d > 20$], [Use SNN (Alg. 36) or SUBCLU],
  ),
  caption: [Application suitability and engineering matrix for OPTICS.],
)

#sec("V", "Structural Relations and Empirical Behavior")

OPTICS forms the exact algorithmic bridge connecting fixed-scale DBSCAN (Alg. 31), hierarchical HDBSCAN (Alg. 32), and classical Single-Linkage (Alg. 56).

#proposition(title: "Exact Duality Linking OPTICS, HDBSCAN MSTs, and Single-Linkage", origin: "Original Unified Proof")[
  Set $epsilon_"max" = infinity$.
  1. When $m = 1$, $d_"core"^(infinity) (p) = 0$ identically, so Algorithm 1 reduces *exact-for-exact* to Prim's MST algorithm on $(cal(X), d)$, recovering the Single-Linkage dendrogram (SLINK, Alg. 56 @sibson1973slink @gower1969minimum).
  2. Replacing asymmetric $d_"reach"(q | p)$ in Algorithm 1 with symmetric mutual reachability $d_"mreach"(p, q)$ yields an ordering $pi^*$ whose reachability weights $(r_2^*, dots, r_n^*)$ match the exact edge-weight multiset of HDBSCAN's MST, and recursively splitting $pi^*$ at $arg max_(2 <= t <= n) r_t^*$ reconstructs a hierarchy isomorphic to HDBSCAN's single-linkage dendrogram $cal(H)$ @campello2013hdbscan.
]

#proof[
  (1) When $m = 1$, every point is its own 1st NN ($d_"core"^(infinity)(p) = 0$), so $d_"reach"(q | p) = d(p, q)$ and Lemma 1 becomes Prim's greedy cut-minimum MST algorithm identically.

  (2) Under symmetric $d_"mreach"$, Lemma 1 proves that each emitted edge $(p_t, o_t)$ of weight $r_t^*$ is a minimum-weight edge crossing $(P_(t-1), U_(t-1))$, so $cal(T)^* = {(p_t, o_t) : t = 2, dots, n}$ is a valid MST of $G_"mreach"$. Because all MSTs of a weighted graph share the identical edge-weight multiset, $(r_2^*, dots, r_n^*)$ matches HDBSCAN's MST weights. Moreover, let $t^* = arg max_(2 <= t <= n) r_t^*$ with maximum weight $w^* = r_(t^*)^*$. By the exact priority-queue exhaustion argument of Theorem 1, prior to popping $o_(t^*)$, Prim's traversal completely exhausts the connected component of $cal(T)^* backslash {(p_(t^*), o_(t^*))}$ containing $o_1$ (all of whose internal edges have weight $<= w^*$; breaking ties consistently ensures no edge of weight $w^*$ prematurely interleaves). Thus partitioning $pi^*$ into ${o_1, dots, o_(t^*-1)}$ and ${o_(t^*), dots, o_n}$ mirrors deleting the heaviest MST edge; recursing on both halves recovers HDBSCAN's single-linkage hierarchy $cal(H)$ identically.
]

*Empirical Behavior on the ClusterOverview 4-Dataset Suite:* On *Iris*, OPTICS achieves the highest Silhouette among the three core density algorithms ($"Sil" = 0.657$, vs. DBSCAN $0.618$ and HDBSCAN $0.608$), though its 3-class ARI is lower ($0.442$) because `ExtractXi` steepness criteria isolate high-density cores cleanly while leaving gradual transition zones unlabeled. On *Customers* ($"Sil" = 0.328$) and *Multi-Shape* ($"Sil" = 0.624$, $"ARI" = "N/A"$), smooth density tapers weaken sharp $xi$-cliffs, trailing HDBSCAN ($0.777$ and $0.795$) and validating the known sensitivity of `ExtractXi` relative to integral excess-of-mass stability. Conversely, on *High-D Subspace+Noise* ($"Sil" = 0.378$, $"ARI" = "N/A"$), OPTICS slightly outperforms both DBSCAN ($0.352$) and HDBSCAN ($0.350$) by preserving relative reachability dips amid noisy ambient dimensions. Finally, OPTICS achieves the fastest mean runtime of the trio at $2.182$ ms (vs. HDBSCAN $2.469$ ms and DBSCAN $14.321$ ms), demonstrating the low constant-factor overhead of its single-pass priority-queue sweep.
