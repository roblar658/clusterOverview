#import "../template/ieee.typ": *

#alg_entry(
  id: 24,
  name: "FastPAM / FasterPAM",
  category: "Cat. I: Centroid & Partitioning",
  year: "2019 / 2021",
  paradigm: "Exact O(k)-Accelerated Swap Local Search (k-Medoids)",
  complexity: $O((n - k)^2 T)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Let $(X, d)$ be a finite dissimilarity space of cardinality $|X| = n$ equipped with a symmetric non-negative dissimilarity function $d : X times X -> bb(R)_(>= 0)$ satisfying $d(x, x) = 0$. Unlike centroid-based clustering, $d$ need not arise from an inner product or Euclidean embedding. The *discrete $k$-medoids* (or *discrete $k$-median*) problem @kaufman1990finding seeks a subset $M = {m_1, dots, m_k} subset.eq X$ of $k in {1, dots, n}$ exemplar objects (*medoids*) minimizing the *total deviation* ("TD") functional:
$ "TD"(M) := sum_(o in X) min_(m in M) d(o, m). $
For any $k >= 2$ and fixed medoid set $M$, let each object $o in X$ cache its nearest medoid $n(o) in M$ and second-nearest medoid $s(o) in M backslash {n(o)}$, together with their dissimilarities:
$ d_n (o) := d(o, n(o)) = min_(m in M) d(o, m), \
  d_s (o) := d(o, s(o)) = min_(m in M backslash {n(o)}) d(o, m). $

Classical PAM (Partitioning Around Medoids, Alg. 4 @kaufman1990finding) consists of a greedy `BUILD` phase ($O(k n^2)$) followed by steepest-descent `SWAP` iterations evaluating all $k(n - k)$ medoid–candidate pairs $(m_i, x_c)$ independently in $O(n)$ time each, incurring $O(k (n - k)^2)$ operations per iteration. *FastPAM1* and *FasterPAM*, introduced by Schubert and Rousseeuw @schubert2019fasterpam @schubert2021fast, exploit an exact algebraic loop-inversion identity that evaluates all $k$ possible medoid removals for a candidate $x_c$ *simultaneously* in a single $O(n)$ pass, eliminating the $O(k)$ factor without any approximation loss.

#definition(title: "Swap Delta & Medoid Removal Loss Vector")[
  Given $M = {m_1, dots, m_k}$ ($k >= 2$) and any non-medoid candidate $x_c in X backslash M$, the *exact swap cost change* is
  $ Delta_"TD" (m_i, x_c) := "TD"((M backslash {m_i}) union {x_c}) - "TD"(M). $
  Define the *medoid removal loss vector* $Delta^- in bb(R)_(>= 0)^k$ with components
  $ Delta_i^- := sum_(o : n(o) = m_i) (d_s (o) - d_n (o)), quad i = 1, dots, k, $
  measuring the exact increase in $"TD"$ when medoid $m_i$ is deleted and every point in its cluster falls back to its cached second-nearest medoid $s(o)$.
]

While *FastPAM1* retains steepest-descent selection (executing the globally best improving pair $(m_(i^*), x_(c^*))$ per full pass), *FasterPAM* @schubert2021fast makes two further algorithmic refinements: (i) *eager* (first-improvement) swapping, executing the best medoid swap $m_(i^*)$ for candidate $x_c$ immediately whenever $Delta_"TD" (m_(i^*), x_c) < 0$, and (ii) replacing $O(k n^2)$ `BUILD` initialization with random or subsampled Linear Approximative BUILD (`LAB`) initialization.

#algbox(
  title: "FastPAM1 / FasterPAM Local Search",
  inputs: [Dissimilarity matrix $d in bb(R)_(>= 0)^(n times n)$, integer $k >= 2$, initial medoids $M = {m_1, dots, m_k}$ (via `BUILD`, `LAB`, or random seeding), mode $in {"FastPAM1", "FasterPAM"}$.],
  outputs: [1-swap locally optimal medoid set $M$ and partition $cal(C)$.],
)[
  1. Compute $(n(o), d_n (o), s(o), d_s (o))$ for all $o in X$, and initialize $Delta_i^- <- sum_(o : n(o) = m_i) (d_s (o) - d_n (o))$ for $i = 1, dots, k$.
  2. *Repeat* until a complete pass over $X backslash M$ yields no swap with $Delta_"TD" < 0$:
  3. #h(1em) Set $(Delta^*, i^*, c^*) <- (0, bot, bot)$.
  4. #h(1em) *For each* non-medoid $x_c in X backslash M$ *do*:
  5. #h(2em) Initialize vector $Delta[1..k] <- Delta^-$ and scalar $Delta^+ <- 0$.
  6. #h(2em) *For each* $o in X$ *do*:
  7. #h(3em) Let $d_(o c) <- d(o, x_c)$ and $i <- n(o)$.
  8. #h(3em) *If* $d_(o c) < d_n (o)$ *then*:
  9. #h(4em) $Delta^+ <- Delta^+ + (d_(o c) - d_n (o))$; $Delta[i] <- Delta[i] + (d_n (o) - d_s (o))$.
  10. #h(3em) *Else if* $d_(o c) < d_s (o)$ *then*:
  11. #h(4em) $Delta[i] <- Delta[i] + (d_(o c) - d_s (o))$.
  12. #h(2em) Let $i_c <- arg min_(1 <= i <= k) Delta[i]$ and $Delta_c <- Delta^+ + Delta[i_c]$.
  13. #h(2em) *If* mode $=$ `"FasterPAM"` and $Delta_c < 0$ *then*:
  14. #h(3em) Swap $m_(i_c) <- x_c$; incrementally update $(n, d_n, s, d_s)$ and $Delta^-$.
  15. #h(2em) *Else if* $Delta_c < Delta^*$ *then* $(Delta^*, i^*, c^*) <- (Delta_c, i_c, c)$.
  16. #h(1em) *If* mode $=$ `"FastPAM1"` and $Delta^* < 0$ *then* swap $m_(i^*) <- x_(c^*)$ and update caches.
  17. *Return* $M$ and $C_j = {o in X : n(o) = m_j}$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We prove four core theoretical guarantees for FastPAM1 and FasterPAM: (1) exact algebraic identity showing zero approximation error relative to PAM; (2) finite strict-descent termination and amortized $O(n + k |R_(i^*)|)$ cache-update complexity; (3) the classical metric 1-swap $5$-approximation theorem; and (4) a sharp robustness dichotomy separating $50\%$ intra-cluster breakdown from $1/n$ global medoid hijacking for $k >= 2$.

#subsec("A", "Exact Simultaneous Swap-Delta Identity")

#theorem(title: "Exact Algebraic Equivalence of FastPAM1 and PAM", origin: "Original Proof")[
  For any medoid set $M = {m_1, dots, m_k}$ ($k >= 2$) and any candidate $x_c in X backslash M$, define the shared addition gain
  $ Delta^+ (x_c) := sum_(o in X) min(0, space d(o, x_c) - d_n (o)) <= 0 $
  and, letting $I_i (x_c) := {o in X : n(o) = m_i, space d(o, x_c) < d_s (o)}$, the medoid-specific interaction correction
  $ gamma_i (x_c) := sum_(o in I_i (x_c)) (max(d(o, x_c), d_n (o)) - d_s (o)) <= 0. $
  Then for every $i in {1, dots, k}$ simultaneously,
  $ Delta_"TD" (m_i, x_c) = Delta_i^- + Delta^+ (x_c) + gamma_i (x_c). $
  Consequently, Steps 5–12 of Algorithm 1 compute all $k$ exact swap deltas in $O(n + k)$ time, and FastPAM1 visits the *identical* sequence of medoid sets as classical PAM (Alg. 4) under consistent tie-breaking.
]

#proof[
  Fix $i in {1, dots, k}$ and let $M' = (M backslash {m_i}) union {x_c}$. For each $o in X$, write $d_(o c) := d(o, x_c)$ and let $delta_o (m_i, x_c) := min_(m in M') d(o, m) - d_n (o)$ be the exact contribution of $o$ to $Delta_"TD" (m_i, x_c)$. We verify identity of $delta_o (m_i, x_c)$ with the $o$-summand of $Delta_i^- + Delta^+ (x_c) + gamma_i (x_c)$ across all five exhaustive cases:
  1. *$n(o) != m_i$ and $d_(o c) >= d_n (o)$:* Since $n(o) in M'$ remains nearest, $delta_o = 0$. All three terms contribute $0$.
  2. *$n(o) != m_i$ and $d_(o c) < d_n (o)$:* Now $x_c$ becomes nearest, so $delta_o = d_(o c) - d_n (o)$. Here $Delta_i^-$ and $gamma_i$ contribute $0$, while $Delta^+ (x_c)$ contributes $d_(o c) - d_n (o)$.
  3. *$n(o) = m_i$ and $d_(o c) >= d_s (o)$:* Removing $m_i$ leaves $s(o)$ as the closest element of $M'$ (since $d_(o c) >= d_s (o) >= d_n (o)$), so $delta_o = d_s (o) - d_n (o)$. Here $Delta_i^-$ contributes $d_s (o) - d_n (o)$, while $Delta^+ (x_c)$ and $gamma_i (x_c)$ contribute $0$.
  4. *$n(o) = m_i$ and $d_n (o) <= d_(o c) < d_s (o)$:* Removing $m_i$ makes $x_c$ strictly closer than all remaining medoids in $M backslash {m_i}$, so $delta_o = d_(o c) - d_n (o)$. The formula yields
     $ (d_s (o) - d_n (o)) + 0 + (d_(o c) - d_s (o)) = d_(o c) - d_n (o). $
  5. *$n(o) = m_i$ and $d_(o c) < d_n (o) <= d_s (o)$:* Again $x_c$ is strictly closest in $M'$, so $delta_o = d_(o c) - d_n (o)$. The formula yields
     $ &(d_s (o) - d_n (o)) + (d_(o c) - d_n (o)) \
       &quad + (d_n (o) - d_s (o)) = d_(o c) - d_n (o). $
  Summing over all $o in X$ proves $Delta_"TD" (m_i, x_c) = Delta_i^- + Delta^+ (x_c) + gamma_i (x_c)$ identically. Inspecting Lines 8–11 of Algorithm 1 confirms that Line 9 adds $d_(o c) - d_n (o)$ to $Delta^+$ and $d_n (o) - d_s (o)$ to $Delta[n(o)]$ (Cases 2 and 5), while Line 11 adds $d_(o c) - d_s (o)$ to $Delta[n(o)]$ (Case 4). Because each $o in X$ updates at most one entry $Delta[n(o)]$ and one scalar $Delta^+$, all $k$ exact deltas are evaluated in $O(n + k)$ operations rather than $O(k n)$.
]

#subsec("B", "Strict Descent, Termination & Incremental Cache Complexity")

#theorem(title: "Finite Termination and Exact Cache-Update Complexity", origin: "Original Proof")[
  Let ${M^((t))}$ be the sequence of medoid sets generated by either FastPAM1 or FasterPAM (Algorithm 1). Then:
  1. Every accepted swap strictly decreases $"TD"(M^((t+1))) < "TD"(M^((t)))$, and the algorithm terminates in at most $binom(n, k) - 1$ swaps at a 1-swap local optimum where $Delta_"TD" (m_i, x_c) >= 0$ for all $m_i in M$ and $x_c in X backslash M$.
  2. Each full pass of FastPAM1 costs $O((n - k)(n + k)) = O(n(n - k))$ distance lookups (a strict $Theta(k)$ improvement over PAM's $O(k n (n - k))$).
  3. After an eager swap replacing $m_(i^*)$ with $x_c$ in FasterPAM, let $R_(i^*) := {o in X : n(o) = m_(i^*) "or" s(o) = m_(i^*)}$ denote the set of affected objects. Updating $(n, d_n, s, d_s)$ and $Delta^-$ requires exactly $O(n + k |R_(i^*)|)$ operations; under balanced clusters ($|R_(i^*)| = Theta(n / k)$), the entire eager cache update runs in $O(n)$ time.
]

#proof[
  (1) By Theorem 1, any executed swap satisfies $Delta_"TD" (m_(i_c), x_c) < 0$, so $"TD"(M^((t+1))) < "TD"(M^((t)))$. Hence no $k$-element subset $M subset.eq X$ can be revisited. Because there are exactly $binom(n, k)$ distinct subsets of size $k$, termination occurs in at most $binom(n, k) - 1$ swaps, and the stopping condition guarantees $Delta_"TD" (m_i, x_c) >= 0$ for every pair $(m_i, x_c)$.

  (2) By Theorem 1, evaluating each of the $n - k$ candidates $x_c in X backslash M$ executes $n$ distance lookups plus $k$ comparisons for $arg min_i Delta[i]$, totaling $(n - k)(n + k) = O(n(n - k))$ operations per pass.

  (3) When $m_(i^*)$ is replaced by $x_c$, partition $X$ into $X backslash R_(i^*)$ and $R_(i^*)$. For any $o in.not R_(i^*)$, both $n(o)$ and $s(o)$ remain in $M'$, so the updated nearest and second-nearest distances are simply the two smallest elements among ${d_n (o), d_s (o), d(o, x_c)}$, computed in $O(1)$ time with a single distance lookup $d(o, x_c)$. For any $o in R_(i^*)$, one of the top two medoids was removed, so scanning all $k$ medoids of $M'$ restores $(n(o), d_n (o), s(o), d_s (o))$ in $O(k)$ time. Recomputing $Delta^-$ via a single linear sweep over $o in X$ costs $O(n)$. Summing yields $O(n - |R_(i^*)|) + O(k |R_(i^*)|) + O(n) = O(n + k |R_(i^*)|)$. Under balanced partitions where each medoid is nearest or second-nearest to $O(n / k)$ objects, $k |R_(i^*)| = O(n)$, reducing each eager swap update to $O(n)$—asymptotically equal to scanning a single candidate.
]

#subsec("C", "Metric 1-Swap 5-Approximation Guarantee")

Because FastPAM1 and FasterPAM terminate at exact 1-swap local optima, they inherit the tight constant-factor approximation guarantee of single-swap local search in metric spaces (cf. @kanungo2004local @schubert2021fast).

#theorem(title: "Tight 5-Approximation of Metric 1-Swap Local Optima", origin: "Classical — Cited")[
  Suppose $(X, d)$ is a metric space (symmetric and satisfying the triangle inequality), and let $M subset.eq X$ ($|M| = k$) be any 1-swap local optimum returned by FastPAM1 or FasterPAM. If $M^*$ is a global minimizer of $"TD"$, then
  $ "TD"(M) <= 5 dot "TD"(M^*). $
]

#proof[
  For each $o in X$, let $m(o) in M$ and $m^*(o) in M^*$ denote its nearest medoids in $M$ and $M^*$ respectively, with $a_o := d(o, m(o))$ and $a_o^* := d(o, m^*(o))$. Define the capture map $phi : M^* -> M$ by $phi(m^*) = arg min_(m in M) d(m^*, m)$, and partition $M$ by in-degree $"deg"(m) = |phi^(-1) (m)|$ into $M_0 = {m : "deg"(m) = 0}$, $M_1 = {m : "deg"(m) = 1}$, and $M_(>= 2) = {m : "deg"(m) >= 2}$. Since $sum_(m in M) "deg"(m) = k = |M_0| + |M_1| + |M_(>= 2)|$, we have
  $ |M_0| = sum_(m in M_(>= 2)) ("deg"(m) - 1) >= 1 / 2 sum_(m in M_(>= 2)) "deg"(m). $
  Hence there exists a family $cal(S)$ of $k$ test swaps $(m, m^*) in M times M^*$ such that: (i) every $m^* in M^*$ appears in exactly one swap; (ii) every $m in M_1$ is paired once with its unique $m^* in phi^(-1) (m)$; (iii) no $m in M_(>= 2)$ appears in $cal(S)$; and (iv) every $m in M_0$ appears in at most two swaps.

  Fix any $(m, m^*) in cal(S)$, and let $O^*(m^*) = {o : m^*(o) = m^*}$ and $O(m) = {o : m(o) = m}$. Because $M$ is a 1-swap local optimum, $Delta_"TD" (m, m^*) >= 0$, which we upper-bound by a feasible reassignment: send every $o in O^*(m^*)$ to $m^*$ (cost change $a_o^* - a_o$), and send every $o in O(m) backslash O^*(m^*)$ to $phi(m^*(o))$. Notice $phi(m^*(o)) != m$: if $m in M_0$ then $phi^(-1)(m) = emptyset$; if $m in M_1$ then $m^*(o) != m^*$ implies $phi(m^*(o)) != m$. By the triangle inequality and optimality of $phi(m^*(o))$ over $M$,
  $ d(o, phi(m^*(o))) &<= d(o, m^*(o)) + d(m^*(o), phi(m^*(o))) \
    &<= a_o^* + d(m^*(o), m(o)) <= 2 a_o^* + a_o, $
  so the reassignment cost change for $o in O(m) backslash O^*(m^*)$ is at most $2 a_o^*$. Keeping all other objects at $m(o)$ gives
  $ 0 <= Delta_"TD" (m, m^*) &<= sum_(o in O^*(m^*)) (a_o^* - a_o) \
    &quad + sum_(o in O(m)) 2 a_o^*. $
  Summing over all $(m, m^*) in cal(S)$, each $m^* in M^*$ is counted once and each $m in M$ is counted at most twice:
  $ 0 &<= ("TD"(M^*) - "TD"(M)) + 2 sum_(m in M) sum_(o in O(m)) 2 a_o^* \
    &= 5 "TD"(M^*) - "TD"(M), $
  proving $"TD"(M) <= 5 "TD"(M^*)$.
]

#subsec("D", "Robustness Dichotomy: Intra-Cluster vs. Global Breakdown")

A widespread misconception attributes a $50\%$ global breakdown point to $k$-medoids for all $k$. We establish the exact mathematical dichotomy separating $k = 1$ (and intra-cluster medoids) from global $k >= 2$ medoid hijacking.

#proposition(title: "Breakdown Dichotomy of k-Medoids", origin: "Original Proof")[
  Let $(X, d)$ be a metric space containing a clean subset $X_"clean"$ of diameter $D_0 < infinity$ and $b$ adversarial outliers moved to distance $R -> infinity$ from $X_"clean"$. Then:
  1. For $k = 1$ (or for any fixed cluster $C_j$ of size $n_j$ held stationary), no outlier can become the optimal medoid whenever $b < n / 2$; hence the 1-medoid breakdown point is $epsilon_1^* = ceil(n / 2) / n -> 1/2$.
  2. For any $k >= 2$, corrupting a *single* point ($b = 1$) to distance $R > n D_0$ forces every global minimizer $M^*$ of $"TD"$ to include that outlier ($x_"out" in M^*$), yielding a global medoid-hijacking breakdown point of $epsilon_k^* = 1 / n$—yet every non-hijacked medoid remains strictly inside $X_"clean"$.
]

#proof[
  (1) For $k = 1$ with $b$ outliers at distance $>= R$ from $X_"clean"$, any clean medoid $m in X_"clean"$ incurs $"TD"({m}) <= (n - b) D_0 + b (R + D_0)$ (when outliers lie within distance $R + D_0$ of $X_"clean"$), whereas any outlier $m_"out"$ incurs $"TD"({m_"out"}) >= (n - b) R$. Thus
  $ "TD"({m_"out"}) - "TD"({m}) >= (n - 2 b) R - n D_0. $
  Whenever $b < n / 2$, $n - 2 b >= 1$, so for $R > n D_0$ no outlier can minimize $"TD"$. Conversely, placing $b = ceil(n / 2)$ coincident outliers at distance $R -> infinity$ makes $"TD"({m_"out"}) = (n - b) R < b R <= "TD"({m})$ for $n$ odd (and ties for $n$ even when perturbed by $epsilon$), establishing $epsilon_1^* = ceil(n / 2) / n$.

  (2) Let $k >= 2$ and move a single point $x_"out"$ to distance $R > n D_0$ from $X_"clean"$. Any $k$-medoid set $M subset.eq X_"clean"$ pays at least $R$ for $x_"out"$, so $"TD"(M) >= R > n D_0$. In contrast, selecting $M' = {x_"out"} union M_"clean"$ with any $k - 1$ points $M_"clean" subset.eq X_"clean"$ pays $0$ for $x_"out"$ and at most $D_0$ for each of the $n - 1$ clean points, yielding $"TD"(M') <= (n - 1) D_0 < R <= "TD"(M)$. Hence every global minimizer must place a medoid on $x_"out"$, proving $epsilon_k^* = 1 / n$. Unlike K-Means (Alg. 1), however, every remaining medoid $m in M' backslash {x_"out"}$ is constrained to $X_"clean"$ and cannot drift continuously out of the data support.
]

#sec("III", "Complexity, Geometry and Operational Profile")

By Theorem 2, FastPAM1 reduces every `SWAP` pass from $O(k (n - k)^2)$ to $O(n(n - k))$, yielding speedups of $10 times$–$100 times$ for $k in [20, 500]$. FasterPAM amplifies this gain further: because eager swaps accept multiple improvements *within* a single scan of $X backslash M$, it typically converges in $T_"eager" approx 2$–$5$ passes (vs. $T_"steep" = Omega(k)$ passes for steepest-descent PAM/FastPAM1 when starting from random medoids), while avoiding the $O(k n^2)$ `BUILD` phase entirely @schubert2021fast.

Memory complexity is $O(n^2)$ when precomputing the symmetric dissimilarity matrix $D$ (or $O(n d)$ with on-the-fly distances, trading compute for RAM), plus negligible $O(n + k)$ auxiliary cache storage for $(n, d_n, s, d_s, Delta^-)$. Critically, the inner loop (Lines 6–11) streams sequentially across row $D[c, :]$ of a row-major distance matrix, maximizing L1/L2 CPU cache prefetch utilization.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for FastPAM / FasterPAM],
    [Time / Space Complexity], [$O((n - k)^2 T)$ time ($O(k)$ faster than PAM); $O(n^2)$ matrix RAM (or $O(n d)$ on-the-fly) $+ O(n+k)$ cache.],
    [Metric Generality], [Works on *any* symmetric dissimilarity $d(x, y) >= 0$ (non-Euclidean, graph, edit, Gower).],
    [Approximation Guarantee], [Exact $5$-approximation in any metric space (Theorem 3); FastPAM1 is bit-identical to PAM.],
    [Invariances], [Strictly invariant to isometries and positive scaling $alpha d$; sensitive to nonlinear monotone warps of $d$.],
    [Outlier Behavior], [$50\%$ breakdown for $k = 1$ / intra-cluster; $1/n$ medoid hijacking for global $k >= 2$ (Proposition 1).],
    [Recommended Companions], [CLARA-MED (Alg. 25) or FastCLARA for $n > 5 times 10^4$; Trimmed K-Means (Alg. 26) for heavy outlier rejection.],
  ),
  caption: [Summary of theoretical guarantees and operational profile of FastPAM / FasterPAM.],
)

#sec("IV", "Engineering and Scientific Applications")

FastPAM and FasterPAM are the premier exact local-search solvers whenever prototypes must be *interpretable real data objects* or when the domain geometry is non-Euclidean.

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Cheminformatics & Drug-Discovery Compound Selection:* Chemical libraries encoded as Morgan/ECFP binary fingerprints are compared via Jaccard–Tanimoto dissimilarity (a valid metric). Arithmetic means produce chemically meaningless fractional bit-vectors, whereas FasterPAM returns real synthesizable exemplar compounds with $O(k)$ faster scaling for large $k$ scaffold libraries.
2. *Bioinformatics & Protein Sequence Family Clustering:* Pairwise evolutionary dissimilarities computed from Smith–Waterman alignments or BLOSUM substitution matrices lack a vector space structure. FastPAM operates directly on the precomputed $n times n$ alignment matrix and yields biologically valid representative sequences.
3. *Time-Series Shape Clustering Under DTW:* Dynamic Time Warping @sakoe1978dtw violates Euclidean geometry and makes centroid averaging (DBA @petitjean2011dba) expensive and non-convex. Running FasterPAM on a precomputed DTW or Shape-Based Distance matrix avoids barycenter iteration while scaling smoothly to hundreds of temporal regimes.
4. *Graph Facility Location & Network Edge Placement:* Selecting $k$ logistics hubs, hospital sites, or CDN edge servers on a weighted road or fiber graph $G = (V, E)$ to minimize total shortest-path latency $sum_(v in V) min_(m in M) d_G (v, m)$ is *literally* the metric $k$-medoids objective, inheriting Theorem 3's exact $5$-approximation.
5. *Mixed-Type Clinical & Enterprise Cohorts (Gower Metric):* Electronic health records combine nominal diagnoses, ordinal severity scales, and skewed laboratory markers via Gower dissimilarity. Medoids serve directly as concrete "archetypal patient" case studies for clinical review.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Despite eliminating the $O(k)$ swap factor, FastPAM/FasterPAM still scales quadratically in $n$, becoming memory- and compute-bound beyond $n approx 5 times 10^4$ (resolved by FastCLARA or CLARA-MED, Alg. 25). In addition, it shares the convex-ball bias of metric $k$-median on elongated non-convex manifolds (resolved by DBSCAN, Alg. 31) and can lose a medoid slot to an isolated extreme outlier when $k >= 2$ (Proposition 1(2)).

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Why FastPAM / FasterPAM Fits], [Caution / Alternative],
    [Cheminformatics (Tanimoto)], [Real lead molecules; exact metric 5-approx.], [Use CLARA-MED (Alg. 25) for $n > 10^5$],
    [Time-Series (DTW)], [No expensive DBA barycenter step needed], [Precompute band-limited DTW matrix],
    [Graph Facility Location], [Exact shortest-path $k$-median objective], [Use CLARANS for $n > 10^5$ graphs],
    [Mixed-Type Healthcare], [Native Gower distance; real patient exemplars], [Trim extreme anomalies (Prop. 1)],
    [Curved / Ring Manifolds], [Poor fit (metric balls split arcs)], [Use DBSCAN (Alg. 31) or Spectral (Algs. 103–105)],
  ),
  caption: [Application suitability and engineering trade-off matrix for FastPAM / FasterPAM.],
)

#sec("V", "Structural Relations and Empirical Behavior")

FastPAM/FasterPAM sits at the exact discrete hub connecting PAM (Alg. 4), K-Medians (Alg. 15 @bradley1997clustering), and subsampled medoid solvers (CLARA-MED Alg. 25, CLARANS @ng2002clarans). The next proposition quantifies the exact integrality gap between discrete medoids and continuous geometric medians.

#proposition(title: "Discrete-to-Continuous Integrality Gap of k-Medoids", origin: "Original Derivation")[
  Let $X subset.eq cal(M)$ lie in an ambient metric space $(cal(M), d)$, let $"TD"_"disc"^* = min_(M subset.eq X, |M| = k) "TD"(M)$ be the optimal discrete $k$-medoids cost, and let $"TD"_"cont"^* = min_(C subset.eq cal(M), |C| = k) "TD"(C)$ be the optimal continuous $k$-median cost. Then:
  1. On $(bb(R), |dot|)$, $"TD"_"disc"^* = "TD"_"cont"^*$ (zero integrality gap).
  2. In any metric space $(cal(M), d)$, $"TD"_"disc"^* <= 2 dot "TD"_"cont"^*$, and this factor $2$ is asymptotically tight.
]

#proof[
  (1) For any finite cluster $S subset.eq bb(R)$ of size $m$ sorted as $z_1 <= dots <= z_m$, the piecewise-linear convex function $f(c) = sum_(r=1)^m |z_r - c|$ has subgradient $0 in partial f(c)$ at every median $c in [z_(floor((m+1)/2)), z_(ceil((m+1)/2))]$. Because the endpoints $z_(floor((m+1)/2)), z_(ceil((m+1)/2))$ belong to $S subset.eq X$, every cluster's continuous $L_1$-median is attained by an actual data point, giving $"TD"_"disc"^* = "TD"_"cont"^*$.

  (2) Let $C^* = {c_1^*, dots, c_k^*} subset.eq cal(M)$ attain $"TD"_"cont"^*$ with optimal clusters $X_1, dots, X_k$ and costs $W_j^* = sum_(x in X_j) d(x, c_j^*)$. For each non-empty $X_j$, pick the closest data point $m_j = arg min_(m in X_j) d(m, c_j^*)$. Since $m_j$ minimizes $d(dot, c_j^*)$ over $X_j$, we have $|X_j| d(m_j, c_j^*) <= sum_(x in X_j) d(x, c_j^*) = W_j^*$. By the triangle inequality,
  $ sum_(x in X_j) d(x, m_j) &<= sum_(x in X_j) (d(x, c_j^*) + d(c_j^*, m_j)) \
    &= W_j^* + |X_j| d(m_j, c_j^*) <= 2 W_j^*. $
  Summing over $j = 1, dots, k$ proves $"TD"_"disc"^* <= 2 "TD"_"cont"^*$. For tightness, take $k = 1$ and the $n$ standard basis vectors $X = {e_1, dots, e_n} subset.eq (bb(R)^n, norm(dot)_1)$: the continuous median at $0$ achieves $"TD"_"cont"^* <= sum_i norm(e_i)_1 = n$, while every discrete medoid $e_j$ incurs $"TD"({e_j}) = 2(n - 1)$, giving ratio $2 - 2 / n -> 2$.
]

*Empirical Behavior on ClusterOverview Benchmarks:* On the ClusterOverview 4-dataset suite, FastPAM / FasterPAM delivers a mean runtime of *0.268 ms* ($3.7 times$ faster than multi-subsample CLARA-MED's 0.988 ms while optimizing over the *entire* dataset). On *Iris*, it achieves Silhouette = *0.437* and ARI = *0.663*—surpassing K-Means|| (ARI = 0.431, Alg. 23) by $+0.232$ ARI thanks to exact deterministic swap optimization across overlapping classes. It likewise attains strong cohesion across *Customers* (Sil = *0.327*), *Multi-Shape* (Sil = *0.713*, ARI = N/A), and *High-D Subspace+Noise* (Sil = *0.302*, ARI = N/A), consistently outperforming squared-Euclidean K-Means|| (Sil = 0.301, 0.609, 0.249) because linear dissimilarity penalizes heavy tails and noisy coordinates far less aggressively than squared $L_2$ inertia @rousseeuw1987silhouettes.
