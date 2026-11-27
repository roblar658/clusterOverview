#import "../template/ieee.typ": *

#alg_entry(
  id: 13,
  name: "Elkan's Accelerated K-Means",
  category: "Cat. I: Centroid & Partitioning",
  year: "2003",
  paradigm: "Metric Triangle-Inequality Exact Acceleration",
  complexity: [$O(n k d T)$ (pruned)],
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Standard Lloyd's K-Means (Alg. 1) @lloyd1982least recomputes all $n k$ point-to-centroid distances at every iteration, even though after the first few iterations most points remain firmly interior to their Voronoi cells and centroids shift by vanishingly small increments. *Elkan's Accelerated K-Means* @elkan2003using eliminates redundant distance evaluations *exactly*—producing the identical sequence of partitions and centroids as Lloyd's algorithm—by exploiting the metric triangle inequality to maintain per-point upper and lower distance bounds across iterations.

Let $(cal(M), d)$ be a metric space (canonically Euclidean space $bb(R)^d$ equipped with the true $L_2$ norm $d(x, y) = norm(x - y)_2$, *not* its square, because squared Euclidean distance violates the triangle inequality). Given $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$ and initial centroids $bold(mu)^((0)) = (mu_1^((0)), dots, mu_k^((0)))$, Elkan's algorithm optimizes the exact WCSS objective $J(cal(C), bold(mu)) = sum_(j=1)^k sum_(i in C_j) d(x_i, mu_j)^2$ while maintaining five auxiliary state structures at iteration $t$:
- current cluster assignments $c(i) in {1, dots, k}$ for $i = 1, dots, n$;
- per-point upper bounds $u(i) in bb(R)_(>= 0)$ with staleness flags $r(i) in {"true", "false"}$;
- an $n times k$ matrix of lower bounds $l(i, j) in bb(R)_(>= 0)$;
- pairwise inter-centroid distances $D(j, j') = d(mu_j^((t)), mu_(j')^((t)))$ and half-nearest-center radii
  $ s(j) = 1 / 2 min_(j' != j) D(j, j'), quad j in {1, dots, k}; $
- centroid displacement magnitudes $delta(j) = d(mu_j^((t)), mu_j^((t+1)))$ after each mean update.

#definition(title: "Admissible Metric Bound State")[
  At iteration $t$ with centroid tuple $bold(mu)^((t))$, a state tuple $(c, u, l, r)$ is *admissible* if for every $i in {1, dots, n}$ and $j in {1, dots, k}$:
  $ &0 <= l(i, j) <= d(x_i, mu_j^((t))), \
    &d(x_i, mu_(c(i))^((t))) <= u(i), $
  with exact tightness $u(i) = d(x_i, mu_(c(i))^((t)))$ whenever $r(i) = "false"$.
]

By combining global point-level filtering ($u(i) <= s(c(i))$), pairwise center filtering ($u(i) <= 1/2 D(c(i), j)$), lower-bound filtering ($u(i) <= l(i, j)$), and *lazy* upper-bound refreshment (recomputing $d(x_i, mu_(c(i)))$ only when at least one competitor $j$ survives all cheap tests), Algorithm 1 avoids the vast majority of $d$-dimensional norm evaluations in well-clustered regimes.

#algbox(
  title: "Elkan's Triangle-Inequality K-Means Algorithm",
  inputs: [Dataset $cal(X) subset.eq bb(R)^d$, integer $k >= 1$, initial centroids $bold(mu)^((0))$, max iterations $T$.],
  outputs: [Exact Lloyd partition $cal(C) = {C_1, dots, C_k}$ and centroids $bold(mu)$.],
)[
  1. *Initialize:* Compute $D(j, j') = d(mu_j^((0)), mu_(j')^((0)))$. Set $l(i, j) <- 0$; assign each $x_i$ to nearest center $c(i)$ (using $1/2 D(c(i), j)$ pruning), storing computed distances in $l(i, j)$, setting $u(i) <- d(x_i, mu_(c(i))^((0)))$ and $r(i) <- "false"$.
  2. *For* $t = 0, 1, dots, T - 1$ *do*:
  3. #h(1em) Compute $D(j, j') <- d(mu_j^((t)), mu_(j')^((t)))$ for $j < j'$ and $s(j) <- 1/2 min_(j' != j) D(j, j')$.
  4. #h(1em) *For each* $i in {1, dots, n}$ with $u(i) > s(c(i))$ *do*:
  5. #h(2em) *For each* $j in {1, dots, k} backslash {c(i)}$ *do*:
  6. #h(3em) *If* $u(i) > l(i, j)$ *and* $u(i) > 1/2 D(c(i), j)$ *then*:
  7. #h(4em) *If* $r(i)$ *then* $d_0 <- d(x_i, mu_(c(i))^((t)))$, $u(i) <- d_0$, $l(i, c(i)) <- d_0$, $r(i) <- "false"$.
  8. #h(4em) *If* $u(i) > l(i, j)$ *and* $u(i) > 1/2 D(c(i), j)$ *then*:
  9. #h(5em) $d_1 <- d(x_i, mu_j^((t)))$; $l(i, j) <- d_1$.
  10. #h(5em) *If* $d_1 < u(i)$ *then* $c(i) <- j$, $u(i) <- d_1$.
  11. #h(1em) Update means $mu_j^((t+1)) <- |C_j^((t))|^(-1) sum_(i in C_j^((t))) x_i$ and drifts $delta(j) <- d(mu_j^((t)), mu_j^((t+1)))$.
  12. #h(1em) *If* $max_j delta(j) = 0$ *then break*.
  13. #h(1em) *Propagate bounds:* Set $l(i, j) <- max(0, l(i, j) - delta(j))$, $u(i) <- u(i) + delta(c(i))$, and $r(i) <- "true"$ for all $i, j$.
  14. *Return* $(cal(C)^((t)), bold(mu)^((t+1)))$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish four foundational theoretical results for Elkan's algorithm: (1) formal soundness of the three triangle-inequality pruning filters; (2) inductive invariant preservation proving *exact* trajectory equivalence to Lloyd's K-Means; (3) a quantitative Voronoi-margin theorem proving that interior points collapse to *zero* $d$-dimensional distance evaluations once cumulative centroid drift falls below their margin; and (4) a constructive high-dimensional counterexample proving worst-case pruning collapse and $Omega(n k d)$ overhead under orthogonal noise.

#subsec("A", "Soundness of Metric Pruning Inequalities")

#lemma(title: "Triangle-Inequality Pruning Rules", origin: "Classical — Cited")[
  Let $(cal(M), d)$ be any metric space, let $(c, u, l, r)$ be an admissible bound state for centroids $bold(mu)$, and let $i in {1, dots, n}$ and $j != c(i)$. Then each of the following three conditions implies $d(x_i, mu_j) >= d(x_i, mu_(c(i)))$ @elkan2003using:
  1. *Global Point Filter:* $u(i) <= s(c(i))$ (simultaneously for all $j != c(i)$).
  2. *Pairwise Center Filter:* $u(i) <= 1 / 2 D(c(i), j)$.
  3. *Lower-Bound Filter:* $u(i) <= l(i, j)$.
]

#proof[
  (2) By the triangle inequality and symmetry of $d$,
  $ D(c(i), j) = d(mu_(c(i)), mu_j) <= d(x_i, mu_(c(i))) + d(x_i, mu_j). $
  Rearranging and substituting admissibility $d(x_i, mu_(c(i))) <= u(i)$ together with the hypothesis $D(c(i), j) >= 2 u(i)$ gives
  $ d(x_i, mu_j) &>= D(c(i), j) - d(x_i, mu_(c(i))) \
                 &>= 2 u(i) - u(i) = u(i) >= d(x_i, mu_(c(i))). $

  (1) By definition, $s(c(i)) = 1/2 min_(j' != c(i)) D(c(i), j') <= 1/2 D(c(i), j)$ for every $j != c(i)$. Hence $u(i) <= s(c(i))$ triggers Claim 2 uniformly across all competing centers $j != c(i)$.

  (3) Combining admissibility of both bounds directly yields $d(x_i, mu_j) >= l(i, j) >= u(i) >= d(x_i, mu_(c(i)))$.
]

#remark(title: "Tie-Breaking Exactness")[
  Under strict distance inequalities (generic data in continuous $bb(R)^d$), Lemma 1 guarantees no strictly closer centroid is ever pruned. When exact ties $d(x_i, mu_j) = d(x_i, mu_(c(i)))$ can occur, replacing `>` with `>=` in Steps 4, 6, and 8 whenever $j < c(i)$ preserves deterministic smallest-index tie-breaking identically with Lloyd's algorithm.
]

#subsec("B", "Inductive Bound Soundness and Exact Lloyd Equivalence")

#theorem(title: "Inductive Admissibility and Exact Trajectory Equivalence", origin: "Original Proof")[
  Assume generic non-tied point-to-centroid distances (or tie-safe strict guards). For any initial centroids $bold(mu)^((0))$ and every iteration $t >= 0$:
  1. The state $(c, u, l, r)$ maintained by Algorithm 1 remains *admissible* for $bold(mu)^((t))$ throughout Steps 3–10 and transitions to an admissible state for $bold(mu)^((t+1))$ after Step 13.
  2. Every assignment and centroid generated by Algorithm 1 coincides identically with Lloyd's K-Means (Alg. 1):
     $ cal(C)_"Elkan"^((t)) equiv cal(C)_"Lloyd"^((t)), \
    bold(mu)_"Elkan"^((t)) equiv bold(mu)_"Lloyd"^((t)) quad forall t >= 0. $
  Consequently, Elkan's algorithm inherits exact WCSS Lyapunov descent, finite termination, local Newton–Raphson convergence @bottou1995convergence, and Pollard's strong consistency @pollard1981strong verbatim.
]

#proof[
  We proceed by induction on iteration $t >= 0$.

  *Base case ($t = 0$):* In Step 1, $l(i, j)$ is initialized either to $0 <= d(x_i, mu_j^((0)))$ or to the exact evaluated distance $d(x_i, mu_j^((0)))$, while $u(i) = d(x_i, mu_(c(i))^((0)))$ and $r(i) = "false"$. Any center $j$ skipped during Step 1 satisfies $1/2 D(c(i), j) >= u(i)$, so by Lemma 1(2) $c(i)$ minimizes $d(x_i, mu_j^((0)))$, proving admissibility and Lloyd equivalence at $t = 0$.

  *Inductive step:* Assume $(c, u, l, r)$ is admissible for $bold(mu)^((t)) = bold(mu)_"Lloyd"^((t))$ at the start of iteration $t$. During Steps 4–10, whenever Step 7 executes, $u(i)$ and $l(i, c(i))$ are set to the exact distance $d(x_i, mu_(c(i))^((t)))$ and $r(i) <- "false"$, preserving admissibility. Whenever Step 9 evaluates $d_1 = d(x_i, mu_j^((t)))$, updating $l(i, j) <- d_1$ is exact; if $d_1 < u(i)$, then since Step 7 has already executed (so $r(i) = "false"$ and $u(i) = d(x_i, mu_(c(i))^((t)))$), we have $d_1 < d(x_i, mu_(c(i))^((t)))$, and updating $c(i) <- j$, $u(i) <- d_1$ maintains exact tightness $u(i) = d(x_i, mu_(c(i))^((t)))$ with $r(i) = "false"$.

  Now suppose for contradiction that after Step 10 some point $i$ is not assigned to its true nearest centroid $j^* = arg min_ell d(x_i, mu_ell^((t)))$, so $d(x_i, mu_(j^*)^((t))) < d(x_i, mu_(c_"end" (i))^((t))) <= u_"end" (i)$. Because $u(i)$ only decreases during Steps 4–10, $d(x_i, mu_(j^*)^((t))) < u(i)$ held throughout the entire pass for point $i$. Hence at the moment $j^*$ was tested (or if Step 4 skipped point $i$), Lemma 1 proves none of the three pruning guards $u(i) <= s(c(i))$, $u(i) <= l(i, j^*)$, or $u(i) <= 1/2 D(c(i), j^*)$ could have held. Thus Step 9 necessarily computed $d(x_i, mu_(j^*)^((t)))$ and Step 10 assigned $c(i) <- j^*$, which could only be overwritten by an even closer center—contradicting optimality of $j^*$. Therefore $cal(C)_"Elkan"^((t)) = cal(C)_"Lloyd"^((t))$, and Step 11 computes $bold(mu)_"Elkan"^((t+1)) = bold(mu)_"Lloyd"^((t+1))$.

  Finally, in Step 13, the triangle inequality gives for every $i$:
  $ d(x_i, mu_(c(i))^((t+1))) &<= d(x_i, mu_(c(i))^((t))) + d(mu_(c(i))^((t)), mu_(c(i))^((t+1))) \
                              &<= u(i) + delta(c(i)), $
  while the reverse triangle inequality gives for every $i, j$:
  $ d(x_i, mu_j^((t+1))) &>= d(x_i, mu_j^((t))) - d(mu_j^((t)), mu_j^((t+1))) \
                         &>= l(i, j) - delta(j). $
  Since distances are non-negative, $l'(i, j) = max(0, l(i, j) - delta(j)) <= d(x_i, mu_j^((t+1)))$, and setting $r(i) <- "true"$ restores admissibility for $bold(mu)^((t+1))$, completing the induction.
]

#subsec("C", "Margin-Based Quantitative Pruning Complexity")

#theorem(title: "Cumulative-Drift Margin Criterion for Zero Distance Work", origin: "Original Derivation")[
  For point $x_i$, define its Voronoi separation margin at iteration $t_0$ by
  $ gamma_i^((t_0)) := &min_(j != c^((t_0)) (i)) d(x_i, mu_j^((t_0))) \
                       &- d(x_i, mu_(c^((t_0)) (i))^((t_0))) >= 0, $
  and let $Delta_max^((tau)) := max_(1 <= j <= k) delta^((tau)) (j)$ denote the maximum centroid drift at step $tau$. Suppose at step $t_0$ all bounds for point $i$ are tight ($u^((t_0)) (i) = d(x_i, mu_(c(i))^((t_0)))$ and $l^((t_0)) (i, j) = d(x_i, mu_j^((t_0)))$). Then for every subsequent iteration $t > t_0$ satisfying
  $ 2 sum_(tau = t_0)^(t - 1) Delta_max^((tau)) <= gamma_i^((t_0)), $
  point $i$ triggers *zero* $d$-dimensional distance evaluations in Steps 7 and 9 at iteration $t$.
]

#proof[
  Let $c = c^((t_0)) (i)$ and suppose inductively that point $i$ performs no distance evaluations and does not switch cluster at iterations $t_0 + 1, dots, t - 1$. At each step $tau in {t_0, dots, t - 1}$, Step 13 updates $u^((tau+1)) (i) = u^((tau)) (i) + delta^((tau)) (c)$ and $l^((tau+1)) (i, j) = max(0, l^((tau)) (i, j) - delta^((tau)) (j)) >= l^((tau)) (i, j) - delta^((tau)) (j)$. Telescoping from $t_0$ to $t$ yields, for every $j != c$,
  $ &u^((t)) (i) - l^((t)) (i, j) \
    &<= u^((t_0)) (i) - l^((t_0)) (i, j) \
    &quad + sum_(tau = t_0)^(t - 1) [delta^((tau)) (c) + delta^((tau)) (j)] \
    &<= -gamma_i^((t_0)) + 2 sum_(tau = t_0)^(t - 1) Delta_max^((tau)) <= 0. $
  Consequently, $u^((t)) (i) <= l^((t)) (i, j)$ holds simultaneously for *every* competitor $j != c$. Thus the branch condition $u(i) > l(i, j)$ in Step 6 evaluates to false for all $j != c$, ensuring neither Step 7 (lazy refresh) nor Step 9 (competitor distance) is executed for point $i$ at iteration $t$, and $c(i)$ remains fixed at $c$.
]

Near a strict CVT equilibrium where $Delta_max^((t)) -> 0$ rapidly (Bottou–Bengio superlinear regime @bottou1995convergence), almost all points lie strictly inside their Voronoi cells ($gamma_i > 0$). Theorem 2 proves that once centroid motion drops below the interior margins, the $d$-dimensional arithmetic collapses from $Theta(n k d)$ to $O((k^2 + n_"bdry" k) d)$, where $n_"bdry" << n$ counts only boundary points actively transitioning between cells.

#subsec("D", "High-Dimensional Pruning Collapse and Worst-Case Bounds")

Despite dramatic speedups on structured low-to-moderate dimensional data, Elkan's algorithm offers *no* worst-case asymptotic improvement over Lloyd's algorithm @vattani2011kmeans @aloise2009np and can suffer complete pruning collapse in high dimensions @beyer1999nearest.

#theorem(title: "Deterministic High-Dimensional Zero-Pruning Construction", origin: "Original Counterexample")[
  For every $k >= 2$, $n >= k$, and dimension $d >= k + 1$, there exists a clustered configuration $(cal(X), bold(mu))$ in $bb(R)^d$ with well-separated simplex centroids ($D(j, j') = sqrt(2) r$) on which *every* point $i in {1, dots, n}$ and *every* center pair $(c(i), j)$ violates both the global point filter $u(i) <= s(c(i))$ and the pairwise center filter $u(i) <= 1/2 D(c(i), j)$ simultaneously. Moreover, as orthogonal noise radius $rho -> infinity$, every point's Voronoi margin vanishes as $gamma_i = O(r epsilon / rho)$, forcing $Omega(n k d)$ distance evaluations per iteration whenever centroid drift exceeds $r epsilon / (2 rho)$.
]

#proof[
  Place the $k$ centroids on the first $k$ coordinate axes at $mu_j = r e_j$ ($r > 0$, $j = 1, dots, k$). Then $D(j, j') = norm(r e_j - r e_(j'))_2 = sqrt(2) r$ for all $j != j'$, giving $s(j) = 1/2 D(j, j') = r / sqrt(2)$ uniformly. Next, assign each point $i$ a class $c(i) in {1, dots, k}$ and set
  $ x_i = epsilon e_(c(i)) + rho v_i, $
  where $0 < epsilon < r$, $rho > r / sqrt(2)$, and $v_i in bb(S)^(d-1)$ lies in the orthogonal complement $"span"(e_1, dots, e_k)^perp$ (non-empty since $d >= k + 1$). By orthogonality,
  $ &d(x_i, mu_(c(i)))^2 = (r - epsilon)^2 + rho^2, \
    &d(x_i, mu_j)^2 = r^2 + epsilon^2 + rho^2 $
  for all $j != c(i)$. Hence $c(i)$ is strictly the unique nearest centroid, yet even with an exact tight upper bound $u(i) = d(x_i, mu_(c(i)))$,
  $ u(i) &= sqrt((r - epsilon)^2 + rho^2) > rho \
         &> r / sqrt(2) = s(c(i)) = 1 / 2 D(c(i), j) $
  for *all* $i in {1, dots, n}$ and *all* $j != c(i)$. Thus Steps 4 and 6 achieve *zero* inter-center pruning across the entire dataset. Finally, rationalizing the difference gives the exact Voronoi margin
  $ gamma_i &= sqrt(r^2 + epsilon^2 + rho^2) - sqrt((r - epsilon)^2 + rho^2) \
            &= (2 r epsilon) / (sqrt(r^2 + epsilon^2 + rho^2) + sqrt((r - epsilon)^2 + rho^2)) \
            &<= (r epsilon) / rho. $
  Consequently, in high-dimensional noisy regimes where orthogonal variance $rho^2 = Theta(d sigma^2)$ dwarfs cluster signal $r^2$, $gamma_i -> 0$, any centroid drift exceeding $r epsilon / (2 rho)$ invalidates the lower bounds $l(i, j)$ as well, and Elkan incurs $Theta(n k d)$ distance work *plus* $O(k^2 d)$ pairwise centroid work and $O(n k)$ scalar bound updates.
]

#sec("III", "Complexity, Geometry and Operational Profile")

In worst-case regimes (initial iteration or Theorem 3's high-dimensional noise), each iteration of Algorithm 1 executes $O((n k + k^2) d)$ floating-point operations. In structured regimes where centroids stabilize quickly, empirical distance evaluations often drop by $10 times$–$50 times$, approaching $O((n_"bdry" k + k^2) d)$ distance arithmetic plus $O(n k)$ scalar bound updates in Step 13. Three hardware and scaling factors govern practical throughput:
1. *Auxiliary Memory Wall:* Storing the $n times k$ matrix $l(i, j)$ demands $O((n + k) d + n k + k^2)$ memory; at $n = 10^6, k = 10^4$, $l(i, j)$ alone requires $80$ GB in double precision, exhausting cache and main memory.
2. *Quadratic Center-Matrix Bottleneck:* Computing $D(j, j')$ requires $O(k^2 d)$ time per iteration; when $k = Omega(sqrt(n))$, center-center work can rival point-center work.
3. *Vectorization vs. Branching:* Because Elkan requires true $L_2$ distances (involving `sqrt`) and conditional per-pair branches, dense BLAS-3 GEMM (`X M^T`) implementations of Lloyd's K-Means can outperform Elkan for very small $k$ ($k <= 8$) or on GPU/TPU architectures, whereas Elkan excels on CPU for moderate-to-large $k$ and low-to-moderate $d$.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for Elkan's K-Means],
    [Time Complexity], [$O((n k + k^2) d T)$ worst-case; $O((n_"bdry" k + k^2) d + n k)$ per pruned iteration.],
    [Space Complexity], [$O((n + k) d + n k + k^2)$ due to full $n times k$ lower-bound matrix $l(i, j)$.],
    [Exactness / Guarantees], [Exact Lloyd trajectory (Theorem 1); inherits smoothed poly-time @arthur2011smoothed and $2^(Omega(sqrt(n)))$ worst-case @vattani2011kmeans.],
    [Metric Generality], [Valid in *any* metric space $(cal(M), d)$ satisfying the triangle inequality (not limited to $L_2$).],
    [Outlier Breakdown], [$epsilon^* = 1 / n$ (identical to Lloyd's K-Means).],
    [Recommended Upgrades], [Hamerly (Alg. 14) @hamerly2010making for $O(n)$ memory; K-Means++ @arthur2007kmeanspp seeding; Mini-Batch @sculley2010webscale for huge $n$.],
  ),
  caption: [Summary of theoretical guarantees and operational profile of Elkan's K-Means.],
)

#sec("IV", "Engineering and Scientific Applications")

Because Theorem 1 guarantees exact equivalence to Lloyd's K-Means, Elkan's algorithm serves as a zero-approximation drop-in accelerator across scientific and industrial clustering pipelines whenever $k$ is non-trivial and data exhibit meaningful Voronoi margins.

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *CPU Production K-Means Backends (e.g., scikit-learn):* Scientific software suites deploy Elkan's algorithm as the default exact solver for dense CPU datasets with $k > 1$, delivering deterministic reproducibility alongside multi-fold speedups over naive Lloyd loops once centroids settle.
2. *Large-$k$ Codebook Design and Vector Quantization:* In LBG/Lloyd vector quantizers @linde1980algorithm and IVF coarse quantizers for billion-scale ANN search (with $k in [256, 4096]$), points quickly lock into narrow Voronoi cells with high $gamma_i / Delta_max$, enabling Elkan to skip over $95\%$ of distance calculations in late iterations.
3. *3D LiDAR Point Clouds and Geospatial Zoning:* Spatial point clouds reside in low dimension ($d = 2, 3$), where inter-centroid separations vastly exceed orthogonal noise ($rho approx 0$ in Theorem 3), maximizing the effectiveness of the Step 4 global skip condition $u(i) <= s(c(i))$.
4. *Color Palette Quantization and Superpixel Dictionary Learning:* Image datasets in 3D CIELAB/RGB or low-dimensional patch descriptors combine large $n$ ($10^5$–$10^6$ pixels) with $k in [16, 256]$, matching Elkan's optimal CPU pruning regime.
5. *General Non-Euclidean Metric Partitioning:* When clustering under $L_1$, $L_p$ ($p >= 1$), spherical geodesic, or Wasserstein metrics—where the squared-$L_2$ inner-product GEMM identity $norm(x - mu)_2^2 = norm(x)_2^2 + norm(mu)_2^2 - 2 x^T mu$ does not exist—Elkan provides exact triangle-inequality acceleration directly on the native metric $d$.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Elkan's acceleration becomes counterproductive under four operational regimes: (a) *high-dimensional noisy subspaces* ($d >> 50$ with diffuse noise), where Theorem 3 shows pruning collapses while bound overhead persists; (b) *extreme $n times k$ products*, where the $O(n k)$ lower-bound matrix overflows RAM or L3 cache (resolved by Hamerly @hamerly2010making or Mini-Batch K-Means @sculley2010webscale); (c) *GPU/SIMD-accelerated small-$k$ workloads*, where branchless BLAS-3 GEMM beats scalar pruning; and (d) *non-metric dissimilarities* (such as squared $L_2$, cosine distance $1 - cos$, or KL divergence), where applying Elkan directly violates Lemma 1 and produces *incorrect* partitions!

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Why Elkan Fits / Prerequisite], [Caution / Alternative],
    [CPU Tabular / SciKit], [Exact Lloyd drop-in for moderate $k, d$], [Use BLAS-3 Lloyd if $k <= 4$],
    [3D LiDAR / GIS], [Large Voronoi margins; $>90\%$ pruning], [Use Hamerly (Alg. 14) if $n k$ huge],
    [Color / Patch VQ], [Low $d$, large $k$: ideal pruning ratio], [Pair with K-Means++ @arthur2007kmeanspp],
    [High-D Genomics/NLP], [Noise collapses margins (Theorem 3)], [Use Spherical (Alg. 12) or PCA],
    [Massive Streams], [$O(n k)$ memory wall], [Use Mini-Batch @sculley2010webscale / Scalable++ @bahmani2012scalable],
  ),
  caption: [Application suitability and domain engineering guide for Elkan's K-Means.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Elkan's algorithm initiated an entire lineage of exact metric-bound K-Means accelerators, most notably Hamerly's algorithm (Alg. 14) @hamerly2010making. The following proposition formalizes the exact dominance-versus-memory tradeoff between them.

#proposition(title: "Pruning Dominance vs. Memory Compression: Elkan vs. Hamerly", origin: "Original Derivation")[
  Let $l_"Ham" (i) := min_(j != c(i)) l(i, j)$ denote Hamerly's single per-point lower bound @hamerly2010making summarizing the closest non-assigned center. Under identical initial tight bounds at step $t_0$:
  1. *Strict Bound Dominance:* For every $t >= t_0$, every point $i$, and every $j != c(i)$, $l_"Ham"^((t)) (i) <= l_"Elkan"^((t)) (i, j)$. Consequently, whenever Hamerly prunes point $i$ via $u(i) <= max(s(c(i)), l_"Ham" (i))$, Elkan also prunes point $i$ completely; conversely, when $u(i) > l_"Ham" (i)$, Elkan still skips every center $j$ satisfying $u(i) <= max(l(i, j), 1/2 D(c(i), j))$, whereas Hamerly must evaluate all $k$ distances for point $i$.
  2. *Space-Time Tradeoff:* Hamerly reduces auxiliary memory from $O(n k)$ to $O(n)$ and bound-update time per iteration from $O(n k)$ to $O(n)$, making Hamerly faster in low dimensions ($d <= 10$) where $O(n k)$ memory traffic dominates, while Elkan remains faster in moderate-to-high dimensions ($d >= 20$) where avoiding individual $O(d)$ distance evaluations outweighs $O(n k)$ scalar updates.
]

#proof[
  (1) At $t = t_0$, $l_"Ham"^((t_0)) (i) = min_(j != c(i)) l_"Elkan"^((t_0)) (i, j) <= l_"Elkan"^((t_0)) (i, j)$ holds by definition. At each propagation step, Hamerly subtracts the worst-case drift among all competing centers, $l_"Ham"^((t+1)) (i) = max(0, l_"Ham"^((t)) (i) - max_(q != c(i)) delta(q))$, whereas Elkan subtracts each center's individual drift $delta(j) <= max_(q != c(i)) delta(q)$. By induction, $l_"Ham"^((t)) (i) <= l_"Elkan"^((t)) (i, j)$ for all $t >= t_0$ and all $j != c(i)$. If Hamerly skips point $i$ at step $t$, then either $u(i) <= s(c(i))$ (triggering Elkan's Step 4) or $u(i) <= l_"Ham"^((t)) (i) <= l_"Elkan"^((t)) (i, j)$ for all $j != c(i)$ (failing Step 6 for every $j$). Conversely, if only a single neighbor $j_1$ violates $u(i) <= l(i, j_1)$ while $k - 2$ distant centers satisfy $u(i) <= l(i, j)$, Elkan evaluates at most 2 distances (Steps 7 and 9 for $j_1$) whereas Hamerly triggers a full $k$-center inner loop.

  (2) Elkan stores and touches $n k$ scalar entries in Step 13 ($Theta(n k)$ cache-line accesses per iteration), whereas Hamerly stores only two scalars $(u(i), l_"Ham" (i))$ per point ($Theta(n)$ accesses). Balancing $O(n k)$ scalar overhead against $O(n_"saved" d)$ saved arithmetic determines the crossover dimension $d$.
]

In the ClusterOverview 4-dataset empirical benchmark suite (evaluated via Silhouette @rousseeuw1987silhouettes and Adjusted Rand Index @hubert1985comparing), Elkan's K-Means reflects the exact Euclidean Voronoi geometry guaranteed by Theorem 1: on *Iris*, it achieves *Sil = 0.447* and *ARI = 0.583* (closely tracking Spherical K-Means's 0.444 / 0.589 while maintaining isotropic Euclidean clusters); on *Customers*, it reaches *Sil = 0.334*, outperforming Spherical K-Means (0.319) because raw spend magnitudes remain intact; and on *Multi-Shape*, it records *Sil = 0.712* (ARI = N/A). On *High-D Subspace+Noise*, however, its score drops to *Sil = 0.239* (ARI = N/A)—well below Spherical (0.303) and G-Means (0.309)—confirming both the Euclidean WCSS sensitivity to high-dimensional noise dimensions and the margin-collapse regime formalized in Theorem 3. Its mean runtime of *0.070 ms* is nearly $5 times$ faster than G-Means (0.338 ms) and close to Spherical K-Means (0.048 ms), with the small gap on $n <= 240$ datasets attributable to $O(k^2 d)$ center-matrix and $O(n k)$ bound-array initialization that amortize into large speedups on production-scale datasets.
