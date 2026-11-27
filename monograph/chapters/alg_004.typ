#import "../template/ieee.typ": *

#alg_entry(
  id: 4,
  name: "K-Medoids (PAM)",
  category: "Cat. I: Centroid & Partitioning",
  year: "1987 / 1990",
  paradigm: "Discrete Medoid Local-Search Optimization",
  complexity: $O(k (n - k)^2 T)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Let $(cal(X), d)$ be a finite dataset of $n$ objects equipped with an arbitrary non-negative dissimilarity function $d : cal(X) times cal(X) -> bb(R)_(>= 0)$ satisfying $d(x, x) = 0$. Crucially, unlike Euclidean $k$-means (Alg. 1), $cal(X)$ need not reside in a vector space, arithmetic means need not exist, and $d$ may be any domain-specific distance (e.g., $L_1$, Gower, edit, DTW, or shortest-path graph distance).

The *discrete $k$-medoids* (uncapacitated $k$-median) problem @kaufman1990finding selects a subset $M = {m_1, dots, m_k} subset.eq cal(X)$ of $k$ actual dataset exemplars (*medoids*) minimizing the *Total Deviation* (TD):
$ "TD"(M) := sum_(x in cal(X)) min_(m in M) d(x, m). $
For any $M$ with $|M| = k >= 2$ and every $x in cal(X)$, denote its nearest and second-nearest medoid distances by
$ d_1 (x; M) &:= min_(m in M) d(x, m), \
  d_2 (x; M) &:= min_(m in M backslash {m^*(x)}) d(x, m), $
where $m^*(x) in M$ attains $d_1 (x; M)$ (breaking ties deterministically by smallest object index) and induces the cluster $C_m := {x in cal(X) : m^*(x) = m}$.

#definition(title: "1-Swap Medoid Local Equilibrium")[
  A medoid subset $M subset.eq cal(X)$ of cardinality $k$ is a *1-swap local optimum* of $"TD"$ if replacing any current medoid $m in M$ with any non-medoid $h in cal(X) backslash M$ cannot strictly reduce total deviation:
  $ Delta "TD"(m -> h) := "TD"((M backslash {m}) union {h}) - "TD"(M) >= 0 $
  for all $m in M$ and all $h in cal(X) backslash M$.
]

Kaufman and Rousseeuw's *Partitioning Around Medoids (PAM)* algorithm @kaufman1990finding solves this combinatorial optimization in two deterministic phases: greedy *BUILD* initialization followed by steepest-descent *SWAP* refinement (Algorithm 1).

#algbox(
  title: "Partitioning Around Medoids (PAM: BUILD + SWAP)",
  inputs: [Finite set $cal(X)$ ($|cal(X)| = n$), dissimilarity matrix $D_(i j) = d(x_i, x_j)$, cluster count $k in {1, dots, n}$.],
  outputs: [Medoid set $M subset.eq cal(X)$ ($|M| = k$) and induced clusters ${C_m}_(m in M)$.],
)[
  1. *Phase I — BUILD (Greedy Initialization):*
  2. #h(1em) Select $m_1 <- arg min_(x in cal(X)) sum_(y in cal(X)) d(y, x)$; set $M <- {m_1}$.
  3. #h(1em) *For* $r = 2, dots, k$ *do*:
  4. #h(2em) Pick $m_r <- arg max_(h in cal(X) backslash M) sum_(x in cal(X)) max(0, d_1 (x; M) - d(x, h))$.
  5. #h(2em) Set $M <- M union {m_r}$ and update $d_1 (x; M)$ for all $x in cal(X)$.
  6. *Phase II — SWAP (Steepest-Descent 1-Swap Search):*
  7. #h(1em) Cache $(m^*(x), d_1 (x; M), d_2 (x; M))$ for all $x in cal(X)$.
  8. #h(1em) *Repeat*:
  9. #h(2em) Find $(m^*, h^*) <- arg min_(m in M, space h in cal(X) backslash M) Delta "TD"(m -> h)$.
  10. #h(2em) *If* $Delta "TD"(m^* -> h^*) < 0$, set $M <- (M backslash {m^*}) union {h^*}$ and refresh $(m^*, d_1, d_2)$; *else break*.
  11. *Return* $M$ and $C_m = {x in cal(X) : m^*(x) = m}$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five formal guarantees for PAM: (1) an exact additive four-case decomposition of swap gains; (2) unconditional strict descent and finite-step acyclicity for arbitrary dissimilarity matrices; (3) monotone submodularity and a self-contained $(1 - 1/e)$ approximation proof for BUILD; (4) the classical constant-factor $5$-approximation of metric 1-swap equilibria; and (5) exact finite-sample breakdown robustness proving $50\%$ resistance for 1-medoid and fixed clusters alongside the global $k >= 2$ outlier-capture threshold.

#subsec("A", "Exact Swap Decomposition and Finite Termination")

#lemma(title: "Exact Additive Swap-Cost Decomposition", origin: "Original Proof")[
  For any medoid set $M subset.eq cal(X)$ ($|M| = k >= 2$), any $m in M$, and any candidate $h in cal(X) backslash M$, the exact swap gain satisfies $Delta "TD"(m -> h) = sum_(x in cal(X)) delta_x (m, h)$, where
  $ delta_x (m, h) = cases(
      min(d(x, h), d_2 (x; M)) - d_1 (x; M) & "if" m^*(x) = m,
      min(0, d(x, h) - d_1 (x; M)) & "if" m^*(x) != m.
    ) $
]

#proof[
  Let $M' = (M backslash {m}) union {h}$, and partition $cal(X)$ into $C_m = {x : m^*(x) = m}$ and $cal(X) backslash C_m$.
  - If $x in C_m$, removing $m$ leaves $min_(u in M backslash {m}) d(x, u) = d_2 (x; M)$. Inserting $h$ yields new nearest distance $d_1 (x; M') = min(d(x, h), d_2 (x; M))$, so the net change is $d_1 (x; M') - d_1 (x; M) = min(d(x, h), d_2 (x; M)) - d_1 (x; M)$.
  - If $x in.not C_m$, its current nearest medoid $m^*(x)$ remains in $M backslash {m}$, so $min_(u in M backslash {m}) d(x, u) = d_1 (x; M)$. Inserting $h$ yields $d_1 (x; M') = min(d_1 (x; M), d(x, h))$, whose change is $min(0, d(x, h) - d_1 (x; M)) <= 0$.
  Summing $d_1 (x; M') - d_1 (x; M)$ over all $x in cal(X)$ yields $Delta "TD"(m -> h) = sum_(x in cal(X)) delta_x (m, h)$.
]

#remark(title: "Algebraic Basis of FastPAM / FasterPAM")[
  Observe in Lemma 1 that the second branch $min(0, d(x, h) - d_1 (x; M))$ depends on $h$ alone, independent of which $m in M$ is removed! Factoring out $sum_(x in cal(X)) min(0, d(x, h) - d_1 (x; M))$ and correcting only $x in C_m$ collapses each SWAP scan from $O(k (n - k)^2)$ to $O((n - k)^2)$—the core speedup of FastPAM and FasterPAM (Alg. 24) @schubert2019fasterpam @schubert2021fast.
]

#theorem(title: "Unconditional Strict Descent and Finite Acyclicity", origin: "Original Proof")[
  For any finite dissimilarity matrix $D in bb(R)_(>= 0)^(n times n)$—without assuming symmetry or the triangle inequality—PAM's SWAP phase strictly decreases $"TD"(M)$ at every executed move and terminates in at most $T <= binom(n, k) - 1$ iterations at a 1-swap local equilibrium.
]

#proof[
  Line 10 of Algorithm 1 executes a swap $(m^*, h^*)$ if and only if $Delta "TD"(m^* -> h^*) < 0$, which enforces the strict inequality $"TD"(M^((t+1))) < "TD"(M^((t)))$ at every step $t$. Consequently, no medoid subset $M in binom(cal(X), k)$ can ever recur, eliminating cycling unconditionally. Because the feasible space has finite cardinality $binom(n, k)$, at most $binom(n, k) - 1$ strict decreases are possible before the stopping condition $min_(m, h) Delta "TD"(m -> h) >= 0$ holds, satisfying Definition 1.
]

#subsec("B", "Submodularity of BUILD and Metric 5-Approximation of SWAP")

#theorem(title: "Monotone Submodularity and (1 - 1/e)-Optimality of PAM BUILD", origin: "Original Proof")[
  For $k >= 2$, fix the initial 1-medoid $m_1$ from line 2 of Algorithm 1, and define the incremental gain function $L : 2^(cal(X) backslash {m_1}) -> bb(R)_(>= 0)$ by
  $ L(S) := "TD"({m_1}) - "TD"({m_1} union S). $
  Then $L$ is normalized ($L(emptyset) = 0$), monotone non-decreasing, and submodular. Consequently, the set $S_(k-1) = {m_2, dots, m_k}$ produced by BUILD satisfies
  $ L(S_(k-1)) >= (1 - 1 / e) max_(|S| = k - 1) L(S). $
]

#proof[
  Decompose $L(S) = sum_(x in cal(X)) L_x (S)$, where $L_x (S) := d(x, m_1) - d_1 (x; {m_1} union S)$. Clearly $L_x (emptyset) = 0$ and adding elements can only decrease $d_1$, so each $L_x$ is non-negative and monotone. For any $S subset.eq cal(X) backslash {m_1}$ and $h in.not S$, the marginal gain is
  $ L_x (S union {h}) - L_x (S) = max(0, d_1 (x; {m_1} union S) - d(x, h)). $
  Whenever $A subset.eq B$, we have $d_1 (x; {m_1} union A) >= d_1 (x; {m_1} union B)$. Since $u |-> max(0, u - d(x, h))$ is non-decreasing, $L_x (A union {h}) - L_x (A) >= L_x (B union {h}) - L_x (B)$. Summing over $x in cal(X)$ proves submodularity of $L$.

  Notice that line 4 of Algorithm 1 chooses $m_(r+1) = arg max_h [L(S_r union {h}) - L(S_r)]$ where $S_r = {m_2, dots, m_(r+1)}$. Let $S^*$ maximize $L$ among $(k - 1)$-subsets. By monotonicity and submodularity,
  $ L(S^*) - L(S_r) &<= L(S^* union S_r) - L(S_r) \
                    &<= sum_(h in S^* backslash S_r) [L(S_r union {h}) - L(S_r)] \
                    &<= (k - 1) [L(S_(r+1)) - L(S_r)]. $
  Setting $a_r := L(S^*) - L(S_r)$ gives $a_(r+1) <= (1 - (k - 1)^(-1)) a_r$. Inducting from $a_0 = L(S^*)$ yields
  $ L(S_(k-1)) &>= (1 - (1 - 1 / (k - 1))^(k-1)) L(S^*) \
               &>= (1 - 1 / e) L(S^*), $
  using $(1 - 1 / q)^q <= e^(-1)$ for all $q >= 1$.
]

#theorem(title: "Constant-Factor 5-Approximation of Metric 1-Swap Equilibria", origin: "Classical — Cited")[
  Suppose $(cal(X), d)$ is a metric space (symmetric, satisfying the triangle inequality), and let $M^*$ be a globally optimal $k$-medoids set. Then *every* 1-swap local optimum $M$ returned by PAM's SWAP phase satisfies @kanungo2004local
  $ "TD"(M) <= 5 "TD"(M^*). $
]

#proof[
  Map each optimal medoid $o in M^*$ to its closest local medoid $phi(o) := arg min_(m in M) d(o, m)$ (breaking ties with $phi(o) = o$ whenever $o in M$). Partition $M$ by capture multiplicity into $M_0 = {m : |phi^(-1)(m)| = 0}$, $M_1 = {m : |phi^(-1)(m)| = 1}$, and $M_(>= 2) = {m : |phi^(-1)(m)| >= 2}$. Since $|phi^(-1)(M_(>= 2))| = k - |M_1| = |M_0| + |M_(>= 2)|$ and each $m in M_(>= 2)$ captures at least two elements of $M^*$, we obtain $2 |M_(>= 2)| <= |M_0| + |M_(>= 2)|$, hence $|M_(>= 2)| <= |M_0|$ and $|phi^(-1)(M_(>= 2))| <= 2 |M_0|$.

  Construct $k$ test swap pairs $(m, o) in M times M^*$: pair each $m in M_1$ with its unique $o in phi^(-1)(m)$, and assign the remaining $|phi^(-1)(M_(>= 2))|$ optimal medoids to elements of $M_0$ so that every $o in M^*$ enters *exactly once*, every $m in M$ exits *at most twice*, and whenever $(m, o)$ is paired, $phi(o') != m$ for all $o' in M^* backslash {o}$.

  Let $o^*(x)$ denote $x$'s optimal medoid in $M^*$ and $C_o^*$ its optimal cluster. For each test pair $(m, o)$, 1-swap optimality enforces $0 <= "TD"((M backslash {m}) union {o}) - "TD"(M)$ (note that if $o in M$, removing $m$ without adding a new medoid still cannot decrease $"TD"$). Upper-bound this difference by a feasible reassignment under $M' = (M backslash {m}) union {o}$:
  - Send every $x in C_o^*$ to $o$, changing cost by $d(x, o^*(x)) - d(x, m^*(x))$.
  - Send every $x in C_m backslash C_o^*$ to $phi(o^*(x)) in M backslash {m}$ (valid since $o^*(x) != o$). By definition of $phi$ and the triangle inequality,
    $ d(x, phi(o^*(x))) &<= d(x, o^*(x)) + d(o^*(x), phi(o^*(x))) \
                        &<= d(x, o^*(x)) + d(o^*(x), m) \
                        &<= 2 d(x, o^*(x)) + d(x, m), $
    increasing cost by at most $2 d(x, o^*(x))$.
  - Keep all remaining points at their current medoid (zero change).

  Summing across all $k$ pairs, each $x in cal(X)$ appears once in the first sum and at most twice in the second sum (since each $m in M$ is swapped out at most twice):
  $ 0 &<= sum_(x in cal(X)) [d(x, o^*(x)) - d(x, m^*(x))] \
      &quad + 2 sum_(x in cal(X)) 2 d(x, o^*(x)) \
      &= 5 "TD"(M^*) - "TD"(M), $
  which rearranges immediately to $"TD"(M) <= 5 "TD"(M^*)$.
]

#subsec("C", "Exact Finite-Sample Breakdown Point Robustness")

#theorem(title: "50% 1-Medoid Breakdown Point and k >= 2 Outlier Isolation", origin: "Original Proof")[
  Let $(cal(X), d)$ be a metric space of $n$ points initially contained in a ball $B_R (x_0)$.
  1. For $k = 1$, replacing any $q <= floor((n - 1)/2)$ points with *arbitrary* outliers leaves the global 1-medoid strictly inside $B_(D_0) (x_0)$ where $D_0 = ((3 n - 2 q) / (n - 2 q)) R$. Thus the Donoho–Huber breakdown point is $epsilon^* = ceil(n / 2) / n -> 1 / 2$.
  2. For $k >= 2$, holding a cluster $C_j$ fixed inherits breakdown $ceil(|C_j|/2)/|C_j| -> 1/2$; globally, a single outlier at distance $D > (n - 1) dot "diam"(cal(X)_"clean")$ from $cal(X)_"clean"$ captures a dedicated singleton medoid, yet *cannot* pull any remaining clean cluster's medoid out of $cal(X)_"clean"$.
]

#proof[
  (1) Let $cal(X)_"clean"$ ($|cal(X)_"clean"| = n - q$) remain in $B_R (x_0)$ and $cal(X)_"out"$ ($|cal(X)_"out"| = q < n / 2$) be arbitrarily corrupted. Pick any clean point $c in cal(X)_"clean"$ and any corrupted point $z$ with $d(z, x_0) = D > D_0$. For every $y in cal(X)_"clean"$, the triangle inequality gives $d(y, z) >= D - R$ and $d(y, c) <= 2 R$, so $d(y, z) - d(y, c) >= D - 3 R$. For every $w in cal(X)_"out"$, the triangle inequality gives $d(w, c) - d(w, z) <= d(c, z) <= D + R$, hence $d(w, z) - d(w, c) >= -(D + R)$. Summing over all $n$ points:
  $ "TD"({z}) - "TD"({c}) &>= (n - q)(D - 3 R) - q (D + R) \
                         &= (n - 2 q) D - (3 n - 2 q) R > 0 $
  whenever $n - 2 q >= 1$ and $D > D_0$. Hence no point outside $B_(D_0) (x_0)$ can ever be the 1-medoid, establishing $epsilon^* = ceil(n / 2) / n$ (since moving $ceil(n / 2)$ points together to infinity shifts the medoid to infinity).

  (2) Applying (1) to any fixed cluster $C_j$ yields cluster-level breakdown $ceil(|C_j|/2)/|C_j|$. Globally for $k >= 2$, if a single outlier $z$ lies at distance at least $D > (n - 1) dot "diam"(cal(X)_"clean")$ from every clean point, omitting $z$ from $M$ incurs $"TD"(M) >= D$, whereas placing a singleton medoid at $z$ and $k - 1$ medoids on $cal(X)_"clean"$ costs at most $(n - 1) dot "diam"(cal(X)_"clean") < D$. Yet all other $k - 1$ medoids are constrained to actual points in $cal(X)_"clean"$, preventing centroid drift.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Precomputing all pairwise dissimilarities requires $Theta(n^2 d)$ time and $Theta(n^2)$ memory. Phase I (BUILD) performs $O(k n^2)$ distance comparisons, and each iteration of classic Phase II (SWAP) evaluates $k (n - k)$ candidate swaps across $n$ points via Lemma 1 in $O(k (n - k)^2)$ time, giving total runtime $O(k (n - k)^2 T)$. Although $T$ is typically $O(k)$ in practice after submodular BUILD seeding, the quadratic memory and cubic-in-$n$ scaling limit classic PAM to $n <= 5000$ objects, motivating FastPAM/FasterPAM (Alg. 24) @schubert2019fasterpam @schubert2021fast, subsampling CLARA (Alg. 5) @srikantha2011clara, and randomized local-search CLARANS (Alg. 6) @ng2002clarans.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for K-Medoids (PAM)],
    [Time / Space Complexity], [$O(k (n - k)^2 T)$ time; $O(n^2)$ dissimilarity matrix space.],
    [Approximation Guarantee], [BUILD: $(1 - 1/e)$ on gain $L$ (Thm. 2); SWAP: $<= 5 "TD"(M^*)$ in metric spaces (Thm. 3).],
    [Metric Generality], [Operates on *any* dissimilarity $d$ (non-Euclidean, mixed-type, graph, edit, DTW).],
    [Prototype Interpretability], [Every center is a real observed exemplar $m_j in cal(X)$ (medoid).],
    [Outlier Breakdown Point], [$epsilon^* -> 50\%$ for 1-medoid & fixed clusters (Thm. 4); zero clean-centroid drift.],
    [Recommended Upgrades], [FastPAM/FasterPAM (Alg. 24) for $O((n-k)^2)$ exact SWAP; CLARA (Alg. 5) / CLARANS (Alg. 6) for $n > 10^4$.],
  ),
  caption: [Summary of theoretical guarantees and operational profile of K-Medoids (PAM).],
)

#sec("IV", "Engineering and Scientific Applications")

PAM excels whenever *real-object interpretability*, *arbitrary non-Euclidean dissimilarities*, or *heavy-tailed outlier resistance* outweigh raw linear-time scalability.

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Bioinformatics & Biological Sequence Clustering:* Protein and DNA sequences lack a Euclidean vector mean; PAM operates directly on alignment-based edit, BLOSUM, or $k$-mer distances, outputting actual representative sequences as biological consensus exemplars.
2. *Mixed-Type Clinical Phenotyping (EHR Data):* Combining continuous biomarkers, ordinal severity scores, and binary ICD diagnoses via Gower dissimilarity makes arithmetic means meaningless. PAM partitions patients cleanly while surfacing actual representative "prototype patients" for clinical audit.
3. *Discrete Facility Location ($p$-Median Network Siting):* Selecting $k$ warehouses, ambulance depots, or cellular base stations on a road network requires facilities to sit on valid graph vertices under shortest-path travel times—verbatim the $"TD"(M)$ objective with a provable factor-5 metric guarantee (Theorem 3).
4. *Time-Series & Financial Regime Identification:* Under Dynamic Time Warping (DTW) or correlation-based distances, Euclidean averaging smears phase-shifted peaks; PAM selects authentic historical time-series regimes as medoids without distortion.
5. *Robust Fraud, Telemetry, and Exemplar Summarization:* Because within-cluster medoids resist up to $50\%$ corruption (Theorem 4) and use linear $d(x, m)$ rather than squared $d(x, m)^2$ penalties, PAM prevents extreme anomalies from skewing cluster prototypes.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

PAM breaks down on (a) *large-scale datasets* ($n > 10^4$), where the $O(n^2)$ distance matrix and $O(k (n - k)^2)$ SWAP sweeps become infeasible (resolve via FasterPAM Alg. 24, CLARA Alg. 5, or CLARANS Alg. 6); (b) *non-convex manifold structures* (interlocking arcs, filaments), where metric Voronoi partitions slice manifolds (resolve via DBSCAN Alg. 31 or HDBSCAN Alg. 32); (c) *high outlier count with fixed $k$*, where distant outliers consume singleton medoid slots (Theorem 4(2); resolve via Trimmed K-Means Alg. 26 or density filtering); and (d) *high-dimensional noise* causing metric concentration.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Why PAM Fits / Prerequisite], [Caution / Alternative],
    [Bioinformatics / Text], [Arbitrary edit / alignment dissimilarities], [Use FasterPAM (Alg. 24) if $n > 5000$],
    [Clinical / Mixed EHR], [Gower distance + real prototype patients], [Use K-Prototypes (Alg. 8) for huge $n$],
    [Facility Location], [Exact discrete $p$-median on road graphs], [Guaranteed $<= 5 times$ OPT (Thm. 3)],
    [DTW Time Series], [Avoids phase-smearing Euclidean means], [Use CLARANS (Alg. 6) for large sets],
    [Large Numeric Data], [Robust exemplar summary if $n <= 5000$], [Use CLARA (Alg. 5) or Mini-Batch (Alg. 3)],
  ),
  caption: [Application suitability and engineering matrix for K-Medoids (PAM).],
)

#sec("V", "Structural Relations and Empirical Behavior")

Restricting prototypes from an ambient continuous space $Omega$ (such as $L_1$ K-Medians, Alg. 15 @bradley1997clustering) to the finite dataset $cal(X)$ incurs a provably tight cost penalty of at most a factor of $2$:

#proposition(title: "Tight Factor-2 Discrete-vs.-Continuous Medoid Penalty", origin: "Original Derivation")[
  Let $(Omega, d)$ be any ambient metric space containing finite $cal(X) subset.eq Omega$, let $"TD"_Omega^* := min_(C subset.eq Omega, space |C| = k) sum_(x in cal(X)) d_C (x)$ be the continuous $k$-median optimum, and let $"TD"_(cal(X))^* := min_(M subset.eq cal(X), space |M| = k) "TD"(M)$ be the discrete $k$-medoids optimum. Then
  $ "TD"_Omega^* <= "TD"_(cal(X))^* <= 2 "TD"_Omega^*, $
  and the constant $2$ is asymptotically tight.
]

#proof[
  The lower bound is immediate because $cal(X) subset.eq Omega$. For the upper bound, let $C^* = {c_1^*, dots, c_k^*} subset.eq Omega$ attain $"TD"_Omega^*$ with clusters $A_1, dots, A_k$. For each non-empty $A_j$, choose the exemplar $m_j := arg min_(y in A_j) d(y, c_j^*) in cal(X)$. For every $x in A_j$, the triangle inequality and minimality of $m_j$ yield
  $ d(x, m_j) &<= d(x, c_j^*) + d(c_j^*, m_j) \
              &<= d(x, c_j^*) + d(x, c_j^*) = 2 d(x, c_j^*). $
  Summing over all $x in cal(X)$ establishes $"TD"_(cal(X))^* <= "TD"({m_1, dots, m_k}) <= 2 "TD"_Omega^*$. To verify tightness, take $k = 1$ on an $n$-leaf star graph where $cal(X)$ is the set of $n$ leaves at unit distance from an ambient center $c in Omega backslash cal(X)$ and pairwise leaf distance is $2$: then $"TD"_Omega^* = n$, while every discrete medoid $m in cal(X)$ incurs $"TD"({m}) = 2 (n - 1)$, giving ratio $2 (1 - 1 / n) -> 2$ as $n -> infinity$.
]

PAM forms the direct ancestor of FastPAM/FasterPAM (Alg. 24) @schubert2019fasterpam @schubert2021fast, CLARA (Alg. 5) @srikantha2011clara, CLARANS (Alg. 6) @ng2002clarans, CLARA-MED (Alg. 25), and message-passing exemplar clustering via Affinity Propagation @frey2007affinity, while sharing its $L_1$ linear-distance geometry with K-Medians (Alg. 15) and categorical K-Modes (Alg. 7) @huang1998extensions.

On the ClusterOverview 4-dataset empirical benchmark, K-Medoids (PAM) achieves *Iris Silhouette = 0.446* and *Iris ARI = 0.649*, *Customers Silhouette = 0.327*, *Multi-Shape Silhouette = 0.713* (ARI N/A), and *High-D Subspace+Noise Silhouette = 0.302* (ARI N/A), with *mean runtime = 1.588 ms*. Notably, PAM attains the highest *Iris ARI (0.649)* among Algorithms 2–4 (outperforming K-Means++ at 0.583 and Mini-Batch at 0.594) because linear-distance medoids resist boundary leverage points between overlapping *Iris* species, while matching K-Means++ on Multi-Shape (0.713) and High-D Subspace+Noise (0.302). Its 1.588 ms runtime—nearly $7 times$ slower than K-Means++ (0.229 ms)—directly reflects its $O(k (n - k)^2)$ SWAP complexity even on moderate-sized datasets.
