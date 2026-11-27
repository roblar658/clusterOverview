#import "../template/ieee.typ": *

#alg_entry(
  id: 19,
  name: "Hartigan–Wong K-Means",
  category: "Cat. I: Centroid & Partitioning",
  year: "1979",
  paradigm: "Immediate-Update Single-Point Transfer Descent",
  complexity: $O(n k d T)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Let $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$ and $k in {1, dots, n}$. Whereas Lloyd's algorithm (Alg. 1) @lloyd1982least performs batch alternating minimization over $(cal(C), bold(mu))$, Hartigan and Wong @hartigan1979algorithm (Algorithm AS 136) directly minimize the *combinatorial within-cluster sum of squares* over partitions $cal(C) = {C_1, dots, C_k}$ with non-empty clusters ($|C_j| >= 1$):
$ J^* (cal(C)) = sum_(j=1)^k W (C_j), quad W (C_j) := sum_(i in C_j) norm(x_i - mu_j)_2^2, $
where $mu_j := |C_j|^(-1) sum_(i in C_j) x_i$ is always kept exact as the arithmetic centroid of $C_j$. Rather than waiting for a full batch reassignment before moving centroids, Hartigan–Wong evaluates whether relocating a *single* point $x_i$ from its current cluster $C_a$ ($n_a := |C_a| >= 2$) to a candidate cluster $C_b$ ($b != a$, $n_b := |C_b| >= 1$) strictly lowers $J^* (cal(C))$ *after accounting for the simultaneous shift of both centroids $mu_a$ and $mu_b$*.

Define the *exact removal gain* $R_a (x_i)$ and *exact insertion cost* $I_b (x_i)$:
$ R_a (x_i) &:= n_a / (n_a - 1) norm(x_i - mu_a)_2^2 quad (n_a >= 2), \
  I_b (x_i) &:= n_b / (n_b + 1) norm(x_i - mu_b)_2^2. $
As proved in Lemma 1 below, transferring $x_i$ from $C_a$ to $C_b$ reduces $J^* (cal(C))$ by exactly $R_a (x_i) - I_b (x_i)$. Whenever $min_(b != a) I_b (x_i) < R_a (x_i)$, the point is transferred immediately to $b^* = arg min_(b != a) I_b (x_i)$ and both centroids are updated in $O(d)$ time via
$ mu_a' &= mu_a - (x_i - mu_a) / (n_a - 1), \
  mu_b' &= mu_b + (x_i - mu_b) / (n_b + 1). $

#definition(title: "Hartigan–Wong (1-Transfer) Equilibrium")[
  A non-empty partition $cal(C) = {C_1, dots, C_k}$ with exact centroids $mu_j = overline(x)_(C_j)$ is a *Hartigan–Wong equilibrium* (or *$1$-opt partition*) if for every cluster $C_a$ with $|C_a| >= 2$, every point $i in C_a$, and every $b != a$,
  $ (|C_a|) / (|C_a| - 1) norm(x_i - mu_a)_2^2 <= (|C_b|) / (|C_b| + 1) norm(x_i - mu_b)_2^2. $
]

To make single-point transfers computationally competitive with batch Lloyd sweeps, Algorithm AS 136 @hartigan1979algorithm pairs two tightly coupled stages: an *Optimal Transfer Stage (`"OPTRA"`)* guided by a *live set* $cal(L) subset.eq {1, dots, k}$ of recently modified clusters, and a *Quick Transfer Stage (`"QTRAN"`)* that restricts checks to each point's closest cluster $c_1 (i)$ and second-closest cluster $c_2 (i)$ in $O(d)$ time per point (Algorithm 1).

#algbox(
  title: "Hartigan–Wong K-Means (Algorithm AS 136)",
  inputs: [Dataset $cal(X) subset.eq bb(R)^d$, integer $k >= 2$, initial centroids $bold(mu)^((0))$, max outer passes $T$.],
  outputs: [A Hartigan–Wong $1$-opt partition $cal(C)$ and centroids $bold(mu)$.],
)[
  1. *Init:* For each $i in {1, dots, n}$, find its two nearest initial centroids $c_1 (i)$ and $c_2 (i)$. Set $C_j = {i : c_1 (i) = j}$, compute exact means $mu_j$, and set live set $cal(L) = {1, dots, k}$.
  2. *For* pass $t = 1, dots, T$ *do*:
  3. #h(1em) *OPTRA stage:* Initialize $cal(L)_"next" = emptyset$. For $i = 1, dots, n$:
  4. #h(2em) Let $a = c_1 (i)$; if $|C_a| = 1$, *continue*.
  5. #h(2em) Let $cal(B)_i = {1, dots, k} backslash {a}$ if $a in cal(L)$, else $cal(B)_i = cal(L) backslash {a}$. If $cal(B)_i = emptyset$, *continue*.
  6. #h(2em) Compute $R_a (x_i)$ and $b^* = arg min_(b in cal(B)_i) I_b (x_i)$.
  7. #h(2em) *If* $I_(b^*) (x_i) < R_a (x_i)$: transfer $i$ from $C_a$ to $C_(b^*)$; update $mu_a, mu_(b^*)$ in $O(d)$; set $c_2 (i) = a$, $c_1 (i) = b^*$; add $a, b^*$ to $cal(L)_"next"$ and $cal(L)$.
  8. #h(2em) *Else* update $c_2 (i)$ if $b^*$ yields a smaller $I_b (x_i)$ than $c_2 (i)$.
  9. #h(1em) If $cal(L)_"next" = emptyset$, *break* (certified $1$-opt). Set $cal(L) = cal(L)_"next"$.
  10. #h(1em) *QTRAN stage:* Repeat passes over $i = 1, dots, n$ testing only $a = c_1 (i)$ vs. $b = c_2 (i)$ (skipping if $|C_a| = 1$ or both $a, b$ unchanged since last pass); if $I_b (x_i) < R_a (x_i)$, swap $(c_1 (i), c_2 (i))$, update $mu_a, mu_b$, and add $a, b$ to $cal(L)$. Stop QTRAN when a full pass makes zero swaps.
  11. *Return* $(cal(C), bold(mu))$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five rigorous theoretical properties of Hartigan–Wong K-Means: (1) exact single-point transfer gain and $O(d)$ centroid update identities; (2) strict inclusion of Hartigan–Wong equilibria strictly inside the set of Lloyd CVTs, with a $3$-point constructive proof of Lloyd escape; (3) an explicit $4$-point counterexample proving that $1$-opt equilibria still suffer unbounded global approximation ratio; (4) correctness of live-set pruning, non-empty cluster invariance, isometry invariance, and constructive *input-order dependence*; and (5) classical complexity and local-search bounds.

#subsec("A", "Exact Transfer Gain and Rank-1 Update Identities")

#lemma(title: "Exact Hartigan–Wong Single-Point Transfer Formula", origin: "Original Derivation")[
  Let $C_a, C_b$ be disjoint clusters with sizes $n_a >= 2$, $n_b >= 1$ and exact centroids $mu_a, mu_b$. For $x_i in C_a$, set $C_a' = C_a backslash {x_i}$ and $C_b' = C_b union {x_i}$. Then:
  1. The updated centroids satisfy $mu_a' = mu_a - (n_a - 1)^(-1) (x_i - mu_a)$ and $mu_b' = mu_b + (n_b + 1)^(-1) (x_i - mu_b)$.
  2. $W (C_a) - W (C_a') = R_a (x_i) = (n_a / (n_a - 1)) norm(x_i - mu_a)_2^2$.
  3. $W (C_b') - W (C_b) = I_b (x_i) = (n_b / (n_b + 1)) norm(x_i - mu_b)_2^2$.
  4. The exact net WCSS decrease of transferring $x_i$ from $C_a$ to $C_b$ is $J^* (cal(C)) - J^* (cal(C)') = R_a (x_i) - I_b (x_i)$.
]

#proof[
  (1) By definition, $(n_b + 1) mu_b' = n_b mu_b + x_i = (n_b + 1) mu_b + (x_i - mu_b)$, dividing by $n_b + 1$ yields $mu_b'$. Similarly, $(n_a - 1) mu_a' = n_a mu_a - x_i = (n_a - 1) mu_a - (x_i - mu_a)$, dividing by $n_a - 1$ yields $mu_a'$.

  (3) From Claim 1, $mu_b' - mu_b = (x_i - mu_b) / (n_b + 1)$ and $x_i - mu_b' = (n_b / (n_b + 1)) (x_i - mu_b)$. Applying the parallel-axis identity to $C_b$ centered at $mu_b'$ gives
  $ W (C_b') &= sum_(r in C_b) norm(x_r - mu_b')_2^2 + norm(x_i - mu_b')_2^2 \
    &= W (C_b) + n_b norm(mu_b - mu_b')_2^2 + norm(x_i - mu_b')_2^2 \
    &= W (C_b) + (n_b + n_b^2) / (n_b + 1)^2 norm(x_i - mu_b)_2^2 \
    &= W (C_b) + n_b / (n_b + 1) norm(x_i - mu_b)_2^2. $

  (2) Applying Claim 3 in reverse to $C_a = C_a' union {x_i}$ (with $|C_a'| = n_a - 1$) gives $W (C_a) - W (C_a') = ((n_a - 1) / n_a) norm(x_i - mu_a')_2^2$. Since $x_i - mu_a' = (x_i - mu_a) + (n_a - 1)^(-1) (x_i - mu_a) = (n_a / (n_a - 1)) (x_i - mu_a)$, squaring produces
  $ W (C_a) - W (C_a') &= (n_a - 1) / n_a dot n_a^2 / (n_a - 1)^2 norm(x_i - mu_a)_2^2 \
    &= n_a / (n_a - 1) norm(x_i - mu_a)_2^2. $

  (4) Because all clusters $C_ell$ ($ell in.not {a, b}$) are unaffected, $J^* (cal(C)) - J^* (cal(C)') = (W (C_a) - W (C_a')) - (W (C_b') - W (C_b)) = R_a (x_i) - I_b (x_i)$.
]

#subsec("B", "Strict Hierarchy of Equilibria and Finite Convergence")

#theorem(title: "Strict Inclusion Inside Lloyd CVTs & Finite Termination", origin: "Original Proof")[
  Let $cal(E)_"HW" (cal(X), k)$ and $cal(E)_"Lloyd" (cal(X), k)$ denote the sets of non-empty Hartigan–Wong equilibria and discrete Lloyd CVTs on $cal(X)$. Then:
  1. *Strict Inclusion:* $cal(E)_"HW" (cal(X), k) subset.eq cal(E)_"Lloyd" (cal(X), k)$ for every dataset $cal(X) subset.eq bb(R)^d$ and every $k <= n$; moreover, every non-degenerate point $x_i in C_a$ with $|C_a| >= 2$ at a Hartigan–Wong equilibrium lies *strictly inside* the open Voronoi cell $V_a (bold(mu))$.
  2. *Strict Separation ($cal(E)_"HW" subset.neq cal(E)_"Lloyd"$):* For every $k >= 2$ and $alpha in (1, 2)$, there exists a $1$-D dataset ($n = k + 1$) admitting a strict Lloyd CVT that Hartigan–Wong escapes in a *single* transfer, reducing WCSS by factor $4 / alpha^2 -> 4^-$.
  3. *Finite-Step Convergence:* Algorithm 1 preserves $|C_j| >= 1$ for all $j$ and terminates in at most $S(n, k)$ strictly descending transfers at a certified element of $cal(E)_"HW" (cal(X), k)$.
]

#proof[
  (1) Suppose $cal(C) in cal(E)_"HW" (cal(X), k)$, and take $i in C_a$ and $b != a$. If $n_a = 1$, then $x_i = mu_a$, so $norm(x_i - mu_a)_2^2 = 0 <= norm(x_i - mu_b)_2^2$. If $n_a >= 2$, Definition 1 gives
  $ norm(x_i - mu_a)_2^2 &<= n_a / (n_a - 1) norm(x_i - mu_a)_2^2 \
    &<= n_b / (n_b + 1) norm(x_i - mu_b)_2^2 <= norm(x_i - mu_b)_2^2, $
  where the first and third inequalities are *strict* whenever $norm(x_i - mu_a)_2 > 0$ or $norm(x_i - mu_b)_2 > 0$. Hence $cal(C) in cal(E)_"Lloyd" (cal(X), k)$ with strict interior Voronoi separation.

  (2) Set $k = 2$, $alpha in (1, 2)$, and $cal(X) = {0, 2, 2 + alpha} subset.eq bb(R)$. Consider $cal(C)_"Lloyd" = {{0, 2}, {2 + alpha}}$ with centroids $mu_1 = 1$, $mu_2 = 2 + alpha$, and WCSS $J^* = 1^2 + 1^2 + 0 = 2$. Check Lloyd Voronoi stability: $|0 - 1| = 1 < 2 + alpha$, $|2 - 1| = 1 < alpha = |2 - (2 + alpha)|$ (strictly, since $alpha > 1$), and $|(2 + alpha) - (2 + alpha)| = 0$. Thus $cal(C)_"Lloyd"$ is a *strict* Lloyd fixed point. Now test a Hartigan–Wong transfer of $x_2 = 2$ from $C_1$ ($n_1 = 2$) to $C_2$ ($n_2 = 1$):
  $ R_1 (2) &= 2 / 1 (2 - 1)^2 = 2, \
    I_2 (2) &= 1 / 2 (2 - (2 + alpha))^2 = alpha^2 / 2. $
  Because $alpha < 2$, $I_2 (2) = alpha^2 / 2 < 2 = R_1 (2)$. Hartigan–Wong transfers $x_2$ to $C_2$ immediately, reaching ${{0}, {2, 2 + alpha}}$ with optimal WCSS $alpha^2 / 2$ (reduction factor $4 / alpha^2 -> 4^-$ as $alpha -> 1^+$). For $k > 2$, append $k - 2$ isolated singletons at $100, 200, dots$.

  (3) Every accepted transfer in OPTRA or QTRAN requires $n_a >= 2$ (preserving $|C_j| >= 1$ by induction) and satisfies $I_b (x_i) < R_a (x_i)$, strictly decreasing $J^* (cal(C))$ by Lemma 1(4). Because $J^* (cal(C))$ takes strictly decreasing values over the finite set of $S(n, k)$ non-empty $k$-partitions, no partition can repeat and both OPTRA and QTRAN halt finitely. Finally, if a full OPTRA pass makes no move ($cal(L)_"next" = emptyset$), any pair $(a, b)$ skipped by Step 5 satisfied $a, b in.not cal(L)$, meaning neither $C_a$ nor $C_b$ changed since $i$ was last verified against $b$, so $I_b (x_i) >= R_a (x_i)$ still holds for *all* $i, b$. Thus the terminal partition is rigorously certified in $cal(E)_"HW" (cal(X), k)$.
]

#subsec("C", "Unbounded Approximation Ratio of 1-Opt Equilibria")

While Theorem 1 shows Hartigan–Wong eliminates every Voronoi-boundary trap of Lloyd's algorithm, single-point transfers cannot simultaneously merge two nearby clusters and split a distant one.

#theorem(title: "Unbounded Approximation Ratio of Hartigan–Wong", origin: "Original Counterexample")[
  For any target factor $M > 0$ and any $k >= 3$, there exists a $4$-point dataset $cal(X) subset.eq bb(R)$ and a non-empty initial partition $cal(C)^((0))$ that is a *strict Hartigan–Wong equilibrium* satisfying
  $ J^* (cal(C)^((0))) / J^* (cal(C)^"opt") > M. $
]

#proof[
  Set $k = 3$, $D = sqrt(M) + 1$, $L = D + 2$, and construct $cal(X) = {0, 1, L, L + D} subset.eq bb(R)$. Initialize $cal(C)^((0)) = {{0}, {1}, {L, L + D}}$ with centroids $mu_1 = 0$, $mu_2 = 1$, $mu_3 = L + D / 2$ and cost $J^* (cal(C)^((0))) = 2 (D / 2)^2 = D^2 / 2$.
  - Clusters $C_1 = {0}$ and $C_2 = {1}$ are singletons ($n_1 = n_2 = 1$), so neither can transfer a point without emptying a cluster.
  - For $x_3 = L in C_3$ ($n_3 = 2$): $R_3 (L) = 2 (D / 2)^2 = D^2 / 2$, whereas
    $ I_2 (L) = 1 / 2 (L - 1)^2 = (D + 1)^2 / 2 > D^2 / 2 = R_3 (L), $
    and $I_1 (L) = L^2 / 2 > I_2 (L)$.
  - For $x_4 = L + D in C_3$: $R_3 (L + D) = D^2 / 2 < (2 D + 1)^2 / 2 = I_2 (L + D) < I_1 (L + D)$.
  Hence $cal(C)^((0)) in cal(E)_"HW" (cal(X), 3)$ strictly, so Algorithm 1 halts with zero transfers. Meanwhile, $cal(C)^"opt" = {{0, 1}, {L}, {L + D}}$ achieves $J^* (cal(C)^"opt") = 1 / 2$, giving $J^* (cal(C)^((0))) / J^* (cal(C)^"opt") = D^2 > M$. For $k > 3$, append $k - 3$ distant singletons.
]

#subsec("D", "Invariance, Breakdown Point, and Order Dependence")

#proposition(title: "Isometry Invariance, Breakdown 1/n, and Order Sensitivity", origin: "Original Proof")[
  1. *Geometric Invariance:* The set $cal(E)_"HW" (cal(X), k)$ and Algorithm 1's trajectories are strictly invariant under rigid Euclidean isometries $T(x) = Q x + b$ ($Q in O(d)$), sensitive to anisotropic scaling, and have Donoho–Huber breakdown point $epsilon^* = 1 / n$.
  2. *Input-Order Dependence (Constructive Proof):* Unlike batch Lloyd updates, Hartigan–Wong's output depends on the sample scanning order: there exists a $4$-point dataset $cal(X) subset.eq bb(R)$ ($k = 3$) and fixed starting partition $cal(C)^((0))$ from which two point orderings converge to *distinct* strict Hartigan–Wong equilibria whose WCSS values differ by an arbitrarily large factor $D^2$.
]

#proof[
  (1) Rigid isometries preserve squared norms and commute with arithmetic means, leaving $R_a (x_i)$ and $I_b (x_i)$ invariant. Sending a single point $norm(x_1)_2 -> infinity$ drives its cluster centroid to infinity (breakdown point $1 / n$).

  (2) Let $D >= 2$ and define $cal(X) = {0, 1, 10 D, 11 D} subset.eq bb(R)$ with $k = 3$ and initial partition $C_1 = {0, 10 D}$ ($mu_1 = 5 D$), $C_2 = {1}$ ($mu_2 = 1$), $C_3 = {11 D}$ ($mu_3 = 11 D$).
  - *Scan Order A ($0$ visited before $10 D$):* At $x = 0 in C_1$, $R_1 (0) = 2 (5 D)^2 = 50 D^2$, while $I_2 (0) = 1 / 2 (1)^2 = 0.5 < 50 D^2$ and $I_3 (0) = 121 D^2 / 2$. Thus $0$ moves immediately to $C_2$, producing ${{10 D}, {0, 1}, {11 D}}$ with WCSS $0.5$. Because $C_1, C_3$ are now singletons and $R_2 (1) = 0.5 < (10 D - 1)^2 / 2$, Algorithm 1 halts at the global optimum $J^* = 0.5$.
  - *Scan Order B ($10 D$ visited before $0$):* At $x = 10 D in C_1$, $R_1 (10 D) = 50 D^2$, while $I_3 (10 D) = 1 / 2 D^2 < 50 D^2$ and $I_2 (10 D) = (10 D - 1)^2 / 2 > D^2 / 2$. Thus $10 D$ moves immediately to $C_3$, producing ${{0}, {1}, {10 D, 11 D}}$ with WCSS $D^2 / 2$. By Theorem 2, this is a strict Hartigan–Wong equilibrium!

  Hence permuting the scan order alters the terminal WCSS by factor $(D^2 / 2) / (1 / 2) = D^2$.
]

#subsec("E", "Classical Complexity, Consistency, and Multi-Start Bounds")

#theorem(title: "Complexity, K-Means++ Seeding & Statistical Consistency", origin: "Classical — Cited")[
  Global minimization of $J^*$ is NP-hard @aloise2009np @mahajan2012planar, while worst-case local descent can require superpolynomial iterations @vattani2011kmeans yet terminates in smoothed polynomial time under Gaussian perturbations @arthur2011smoothed. Initializing Hartigan–Wong with K-Means++ (Alg. 2) @arthur2007kmeanspp strictly preserves the $O(log k)$ expected approximation bound (since Algorithm 1 is monotonically descending in $J^*$), while multi-centroid swap search @kanungo2004local attains a $(9 + epsilon)$-approximation. As $n -> infinity$ under i.i.d. draws with finite second moment, global empirical minimizers converge almost surely to population minimizers @pollard1981strong.
]

#proof(label: "Proof Synthesis (Cited)")[
  Because Lemma 1 guarantees $J^* (cal(C)_"HW") <= J^* (cal(C)^((0)))$, seeding Algorithm 1 with K-Means++ @arthur2007kmeanspp immediately inherits $bb(E)[J^* (cal(C)_"HW")] <= bb(E)[J^* (cal(C)^((0)))] <= 8 (ln k + 2) J^* (cal(C)^"opt")$ while strictly tightening every local Voronoi equilibrium via Theorem 1(1). Empirical studies in applied statistics @hartigan1979algorithm @hastie2009elements confirm that multi-start Hartigan–Wong consistently outperforms batch Lloyd from identical seeds.
]

#sec("III", "Complexity, Geometry and Operational Profile")

A single full OPTRA pass costs $O(n k d)$ in the worst case when all clusters are live, with $O((n + k) d)$ memory. However, two mechanisms make Hartigan–Wong dramatically faster in practice on CPU workloads: (1) *live-set pruning* restricts distance checks for non-live clusters to $b in cal(L)$, shrinking per-pass work toward $O(n |cal(L)| d)$ as clusters stabilize; and (2) *QTRAN* resolves nearly all boundary oscillations between $c_1 (i)$ and $c_2 (i)$ in $O(d)$ arithmetic per point without scanning all $k$ centroids.

The primary systems trade-off is *sequential data dependency*: because every accepted transfer immediately mutates $mu_a$ and $mu_b$ before point $i + 1$ is evaluated, Hartigan–Wong cannot be expressed as a single dense BLAS-3 matrix product $X M^T$ like batch Lloyd (Alg. 1) or Mini-Batch K-Means @sculley2010webscale. Consequently, Hartigan–Wong is unmatched on single-threaded CPU datasets of small-to-moderate scale ($n <= 10^5$), whereas batch/Elkan/GPU Lloyd dominates billion-scale parallel indexing.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for Hartigan–Wong K-Means],
    [Time / Space Complexity], [$O(n k d T)$ worst-case; $O(d)$ QTRAN inner loop + live-set pruning; $O((n + k) d)$ space.],
    [Local Optimality Strength], [Strictly stronger than Lloyd: $cal(E)_"HW" subset.neq cal(E)_"Lloyd"$ (Thm. 1); never creates empty clusters.],
    [Global Guarantees], [Monotonic descent; unbounded worst-case ratio (Thm. 2); $O(log k)$ when K-Means++ seeded.],
    [Invariances / Order], [Rigid isometry invariant; scale-sensitive; *dependent on input scan order* (Prop. 1).],
    [Outlier Breakdown Point], [$epsilon^* = 1/n$ (standard arithmetic mean vulnerability).],
    [Recommended Pairing], [K-Means++ (Alg. 2) seeding + multi-start shuffling; batch Elkan/Hamerly (Algs. 13–14) for GPU/huge $n$.],
  ),
  caption: [Summary of theoretical guarantees and operational profile of Hartigan–Wong K-Means (AS 136).],
)

#sec("IV", "Engineering and Scientific Applications")

Because Algorithm AS 136 eliminates empty-cluster degeneracies, escapes Voronoi-boundary fixed points, and runs with tiny constant factors on CPU cache hierarchies, it has served as the default `stats::kmeans` algorithm in the R statistical ecosystem for decades @hartigan1979algorithm @hastie2009elements.

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Statistical Computing & Exploratory Tabular Clustering:* In biostatistics, econometrics, and social-science surveys ($n in [50, 10^5]$, $d in [2, 100]$), small cluster sizes ($n_j in [2, 20]$) are common. Because the Hartigan–Wong scaling factor $(n_a / (n_a - 1)) / (n_b / (n_b + 1))$ reaches up to $4$ for tiny clusters, it actively escapes suboptimal small-cluster traps that stall Lloyd's algorithm (Theorem 1(2)).
2. *Local Polishing of K-Means++ & Coreset Partitions:* Running a few passes of Hartigan–Wong as a post-processing refinement step on top of K-Means++ (Alg. 2) or Global K-Means (Alg. 17) guarantees strict improvement (or preservation) of WCSS with negligible latency overhead.
3. *Low-Latency Embedded & Edge CPU Micro-Clustering:* On edge microcontrollers and single-core analytical pipelines without GPU/SIMD matrix accelerators, QTRAN's $O(d)$ transfer test and zero matrix-allocation footprint achieve sub-millisecond convergence.
4. *Chemometrics, Psychometrics & Ecological Ordination:* Standardized spectroscopic, behavioral, and species-abundance tables frequently contain tightly packed overlapping clusters where batch Lloyd can oscillate or collapse a centroid into an empty cluster; Hartigan–Wong's non-empty invariance preserves exact $k$-class structure.
5. *Interactive Data Exploration & Rapid Multi-Start Ensembles:* Because a single Hartigan–Wong run terminates in tens of microseconds on moderate $n$, analysts can execute $50$–$100$ randomized K-Means++ restarts (mitigating Proposition 1(2) order sensitivity and Theorem 2 traps) within milliseconds.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Hartigan–Wong struggles under four regimes: (a) *GPU-accelerated / billion-scale vector search*, where sequential point-by-point centroid updates block SIMD/GPU batching; (b) *unshuffled sorted inputs*, where Proposition 1(2) order dependence can lock into asymmetric local minima (mitigated by random shuffling prior to each restart); (c) *non-convex or heavy-tailed data* (shared with all Euclidean $L_2$ prototype methods); and (d) *very large $k$ with small $d$*, where maintaining exact $(c_1, c_2)$ bounds without spatial indexing loses efficiency to Hamerly/Elkan pruning.

#figure(
  table(
    columns: (1.05fr, 1.15fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Why Hartigan–Wong Fits], [Caution / Alternative],
    [Biostats / Surveys], [Stricter $1$-opt equilibria; no empty clusters], [Shuffle input + K-Means++ multi-start],
    [Edge / CPU Analytics], [Ultra-fast QTRAN $O(d)$ cache-friendly loop], [Standardize heterogeneous units],
    [Post-Seed Polishing], [Strictly tightens any Lloyd/GKM solution], [Costly if $n > 10^6$ without coreset],
    [GPU / Billion ANN], [Avoid: sequential updates block BLAS-3], [Use Mini-Batch or Elkan (Alg. 13)],
    [Outlier / Non-Convex], [Avoid: $epsilon^* = 1/n$ & convex cells], [Use PAM (Alg. 4) or DBSCAN (Alg. 31)],
  ),
  caption: [Application suitability and engineering trade-off matrix for Hartigan–Wong K-Means.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Hartigan–Wong unifies batch Lloyd (Alg. 1), MacQueen's online $k$-means @macqueen1967some, and Ward's hierarchical linkage (Alg. 59) @ward1963hierarchical through exact finite-size weighting factors.

#proposition(title: "Ward Duality, MacQueen Dominance & Large-n Lloyd Limit", origin: "Original Unified Proof")[
  Let $n_"min" := min_(1 <= j <= k) |C_j|$ for a non-empty partition $cal(C)$.
  1. *Exact Ward Duality:* The Hartigan–Wong insertion cost $I_b (x_i)$ equals the exact Ward merge cost @ward1963hierarchical between singleton ${x_i}$ and cluster $C_b$, while $R_a (x_i)$ equals the exact Ward merge cost between ${x_i}$ and $C_a backslash {x_i}$.
  2. *Strict Dominance Over MacQueen's Rule:* MacQueen's sequential update @macqueen1967some moves $x_i$ from $C_a$ to $C_b$ only when $norm(x_i - mu_b)_2^2 < norm(x_i - mu_a)_2^2$. Every MacQueen transfer is strictly accepted by Hartigan–Wong (strictly lowering $J^*$), whereas Hartigan–Wong additionally accepts all improving transfers in the finite-size gap
     $ n_b / (n_b + 1) norm(x_i - mu_b)_2^2 < n_a / (n_a - 1) norm(x_i - mu_a)_2^2 $
     with $norm(x_i - mu_b)_2 >= norm(x_i - mu_a)_2$.
  3. *Asymptotic Lloyd Equivalence:* As $n_"min" -> infinity$, the multiplicative Hartigan–Wong boundary distortion $(n_a / (n_a - 1)) / (n_b / (n_b + 1)) = 1 + O(n_"min"^(-1)) -> 1$, so $cal(E)_"HW"$ approaches $cal(E)_"Lloyd"$ asymptotically on dense balanced clusters.
]

#proof[
  (1) Setting $C_1 = {x_i}$ ($n_1 = 1, mu_1 = x_i$) and $C_2 = C_b$ ($n_2 = n_b, mu_2 = mu_b$) in Ward's merge formula $(n_1 n_2 / (n_1 + n_2)) norm(mu_1 - mu_2)_2^2$ gives $(n_b / (n_b + 1)) norm(x_i - mu_b)_2^2 = I_b (x_i)$. Setting $C_2 = C_a backslash {x_i}$ ($n_2 = n_a - 1$) and noting $x_i - mu_a' = (n_a / (n_a - 1))(x_i - mu_a)$ gives $((n_a - 1) / n_a) norm(x_i - mu_a')_2^2 = R_a (x_i)$.

  (2) Because $n_b / (n_b + 1) < 1 < n_a / (n_a - 1)$ for all $n_a >= 2, n_b >= 1$, any point satisfying MacQueen's criterion $norm(x_i - mu_b)_2^2 < norm(x_i - mu_a)_2^2$ automatically satisfies $I_b (x_i) < norm(x_i - mu_b)_2^2 < norm(x_i - mu_a)_2^2 < R_a (x_i)$, so MacQueen never increases $J^*$ yet halts prematurely on every Lloyd trap of Theorem 1(2).

  (3) Expanding in $n_"min"^(-1)$ gives $n_a / (n_a - 1) = 1 + n_a^(-1) + O(n_a^(-2))$ and $(n_b + 1) / n_b = 1 + n_b^(-1)$, whose product is $1 + n_a^(-1) + n_b^(-1) + O(n_"min"^(-2)) = 1 + O(n_"min"^(-1)) -> 1$.
]

In empirical benchmarks on the ClusterOverview 4-dataset suite @rousseeuw1987silhouettes @hubert1985comparing, Hartigan–Wong K-Means delivers: *Iris* $"Sil" = 0.447$, $"ARI" = 0.583$; *Customers* $"Sil" = 0.332$; *Multi-Shape* $"Sil" = 0.712$ ($"ARI" = "N/A"$); and *High-D Subspace+Noise* $"Sil" = 0.259$ ($"ARI" = "N/A"$), with an extraordinary mean runtime of just *$0.050$ ms*. This runtime is $8.0 times$ faster than Global K-Means ($0.398$ ms, Alg. 17) and $28.6 times$ faster than ISODATA ($1.432$ ms, Alg. 18), vividly validating OPTRA live-set pruning and QTRAN's $O(d)$ inner loop. Meanwhile, while matching Global K-Means identically on Iris and Multi-Shape, its lower Silhouette on Customers ($0.332$ vs. $0.365$ for Alg. 17 and $0.427$ for Alg. 18) and High-D Subspace+Noise ($0.259$ vs. $0.303$ for Alg. 17) empirically mirrors Theorem 2 and Proposition 1(2): single-point local transfer refinement is blisteringly fast and escapes Voronoi-boundary traps, yet still depends on initialization quality and scanning order to select the best global basin.
