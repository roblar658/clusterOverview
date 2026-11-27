#import "../template/ieee.typ": *

#alg_entry(
  id: 53,
  name: "OPTICS-OF",
  category: "Cat. II: Density & Mode-Seeking",
  year: "1999",
  paradigm: "Density-Reachability Ordering & Local Outlier Scoring",
  complexity: $cal(O)(n log n) "–" cal(O)(n^2)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Classical density-based clusterers such as DBSCAN @ester1996dbscan classify any point outside a dense component as undifferentiated binary noise, while hierarchical density ordering in OPTICS @ankerst1999optics represents multi-scale cluster structure via a one-dimensional reachability plot. However, when a dataset simultaneously contains clusters of sharply different local densities, an isolated anomaly hovering near a tight cluster can exhibit a smaller absolute reachability distance than interior members of a sparse cluster. OPTICS-OF (OPTICS with Outlier Factors), introduced by Breunig, Kriegel, Ng, and Sander @breunig1999opticsof, resolves this tension by augmenting the OPTICS reachability traversal with a continuous, scale-relative local outlier score computed from the same neighborhood queries.

Let $cal(D) = {x_1, dots, x_n} subset (cal(X), d)$ be a finite dataset in a metric space $(cal(X), d)$ with radius $epsilon > 0$ and $m = "MinPts" >= 2$. For each $p in cal(D)$, define its punctured $epsilon$-neighborhood $N_epsilon (p) = {q in cal(D) without {p} mid(|) d(p, q) <= epsilon}$ and let $k = m - 1 >= 1$. When $|N_epsilon (p)| >= k$, let $d_k (p)$ denote its $k$-distance (at least $k$ points of $cal(D) without {p}$ lie within $d_k (p)$ and at most $k - 1$ strictly within $d_k (p)$). Following @ankerst1999optics and @breunig1999opticsof, the *core distance* and *reachability distance* are
$
  d_"core" (p) &= cases(
    d_k (p) & "if" |N_epsilon (p)| >= k,
    infinity & "otherwise",
  ), \
  d_"reach" (p, o) &= max { d_"core" (o), d(p, o) }.
$
Whereas OPTICS records only the minimal reachability $r(p) = min_(o "processed") d_"reach" (p, o)$ along a priority-queue traversal, OPTICS-OF averages reachability across $N_epsilon (p)$ to quantify local scale contrast.

#definition(title: "Local Reachability Distance and OPTICS-OF Score")[
  Let $cal(D) subset (cal(X), d)$ satisfy $|N_epsilon (p)| >= k$ for every $p in cal(D)$ (equivalently, $epsilon$ is chosen at least as large as $max_(p in cal(D)) d_k (p)$ so every object has finite core distance). The *local reachability distance* of $p in cal(D)$ is the mean reachability distance of $p$ from its $epsilon$-neighbors @breunig1999opticsof:
  $
    "lrd"_(m, epsilon) (p) = 1 / (|N_epsilon (p)|) sum_(o in N_epsilon (p)) d_"reach" (p, o).
  $
  Equivalently, defining the local reachability density $rho_(m, epsilon) (p) = 1 / "lrd"_(m, epsilon) (p)$ (whenever $"lrd"_(m, epsilon) (p) > 0$), the *OPTICS-OF local outlier factor* of $p$ is
  $
    "OF"_(m, epsilon) (p)
    &= 1 / (|N_epsilon (p)|) sum_(o in N_epsilon (p))
       ("lrd"_(m, epsilon) (p)) / ("lrd"_(m, epsilon) (o)) \
    &= 1 / (|N_epsilon (p)|) sum_(o in N_epsilon (p))
       (rho_(m, epsilon) (o)) / (rho_(m, epsilon) (p)).
  $
  Objects with $"OF"_(m, epsilon) (p) approx 1$ lie in density-homogeneous regions; objects with $"OF"_(m, epsilon) (p) >> 1$ have substantially larger mean reachability distance (lower local density) than their neighbors and are flagged as local outliers at a user threshold $tau > 1$, while clusters are extracted from the OPTICS reachability plot $(pi, r)$ via $xi$-steep valley detection or a level-set cut $epsilon' <= epsilon$ @ankerst1999optics.
]

#algbox(
  title: "OPTICS-OF Joint Ordering and Local Outlier Scoring",
  inputs: [Dataset $cal(D) = {x_1, dots, x_n}$, radius $epsilon > 0$, $m = "MinPts" >= 2$, threshold $tau > 1$],
  outputs: [OPTICS order $pi$, reachabilities $r$, outlier scores ${"OF"(p)}$, outlier set $cal(O)$],
)[
  1. Initialize `processed[p] <- false`, $r(p) <- infinity$, empty order $pi$, and min-heap $Q$.
  2. *OPTICS Traversal & Neighborhood Caching*: For each unprocessed $p in cal(D)$:
     (a) Query $N_epsilon (p)$ and set $d_"core" (p) <- d_(m-1) (p)$ if $|N_epsilon (p)| >= m - 1$, else $infinity$; mark `processed[p] <- true` and append $p$ to $pi$.
     (b) If $d_"core" (p) < infinity$, update $r(q) <- min(r(q), max{d_"core" (p), d(q, p)})$ for each unprocessed $q in N_epsilon (p)$ and insert / decrease-key $q$ in $Q$.
     (c) While $Q$ is non-empty, pop $q = arg min_(u in Q) r(u)$, compute $N_epsilon (q)$ and $d_"core" (q)$, mark processed, append $q$ to $pi$, and update unprocessed neighbors of $q$ in $Q$.
  3. *Local Reachability Aggregation*: For each $p in cal(D)$ with $|N_epsilon (p)| >= m - 1$, compute
     $"lrd"(p) <- |N_epsilon (p)|^(-1) sum_(o in N_epsilon (p)) max{d_"core" (o), d(p, o)}$
     (assigning $"lrd"(p) <- infinity$ if $|N_epsilon (p)| < m - 1$).
  4. *Outlier Factor Evaluation*: For each $p in cal(D)$, compute
     $"OF"(p) <- |N_epsilon (p)|^(-1) sum_(o in N_epsilon (p)) "lrd"(p) / "lrd"(o)$
     (setting $"OF"(p) <- infinity$ for isolated non-core points).
  5. *Joint Output*: Extract clusters from $(pi, r)$ via $xi$-steep valleys or cut $epsilon' <= epsilon$, and return $cal(O) = {p in cal(D) mid(|) "OF"(p) > tau}$.
]

Notice that Step 2 executes the exact priority-queue traversal of OPTICS @ankerst1999optics, while Steps 3--4 reuse the cached neighborhood graph $G_epsilon = (cal(D), E_epsilon)$ and core distances in $cal(O)(|E_epsilon|)$ time without issuing extra range queries.

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five formal guarantees covering interior score stability, finite-sample outlier amplification and self-contamination limits, metric symmetries, separation from global reachability thresholding, and joint algorithmic exactness. Throughout, we assume points in $cal(D)$ are distinct so $d(p, q) > 0$ for $p != q$, ensuring strictly positive $"lrd"(p)$.

#subsec("A", "Interior Stability and Finite-Sample Outlier Bounds")

#lemma(title: "Two-Sided Bounds for Cluster-Interior Points", origin: "Classical — Cited")[
  Let $C subset.eq cal(D)$ be a cluster of size $|C| >= m$ whose distinct points satisfy $0 < d_"min" <= d(p, q) <= d_"max" <= epsilon$ for all $p, q in C$, $p != q$, and suppose $C$ is separated from $cal(D) without C$ by distance strictly greater than $epsilon$, i.e., $d(p, z) > epsilon$ for all $p in C$ and $z in cal(D) without C$ @breunig1999opticsof. Then for every $p in C$:
  $
    d_"min" <= d_"core" (p) <= d_"max",
    quad d_"min" <= "lrd"_(m, epsilon) (p) <= d_"max",
  $
  and consequently
  $
    (d_"min") / (d_"max") <= "OF"_(m, epsilon) (p) <= (d_"max") / (d_"min").
  $
]

#proof[
  Fix any $p in C$. By the separation assumption, every $z in cal(D) without C$ has $d(p, z) > epsilon$, whereas every $q in C without {p}$ satisfies $d(p, q) <= d_"max" <= epsilon$. Hence $N_epsilon (p) = C without {p}$, with cardinality $|N_epsilon (p)| = |C| - 1 >= m - 1 = k$. Because every distance from $p$ to $N_epsilon (p)$ lies in $[d_"min", d_"max"]$, the $k$-th order statistic $d_"core" (p) = d_k (p)$ obeys $d_"min" <= d_"core" (p) <= d_"max"$.

  Since this holds uniformly for every $o in C$, each pair $p, o in C$ ($p != o$) satisfies
  $
    d_"min" <= max { d_"core" (o), d(p, o) } = d_"reach" (p, o) <= d_"max".
  $
  Averaging over the $|N_epsilon (p)|$ neighbors $o in N_epsilon (p) = C without {p}$ preserves both bounds: $d_"min" <= "lrd"_(m, epsilon) (p) <= d_"max"$. Finally, for every $o in N_epsilon (p) subset C$, we obtain the ratio bounds
  $
    (d_"min") / (d_"max") <= ("lrd"_(m, epsilon) (p)) / ("lrd"_(m, epsilon) (o)) <= (d_"max") / (d_"min").
  $
  Averaging over $o in N_epsilon (p)$ yields the claimed two-sided inequality for $"OF"_(m, epsilon) (p)$. In particular, in a homogeneous cluster where $d_"max" / d_"min" -> 1$, $"OF"_(m, epsilon) (p) -> 1$ uniformly across $C$.
]

We next prove an exact, non-asymptotic sandwich theorem for an isolated outlier adjacent to a cluster, revealing both its linear amplification with separation ratio and a sharp *self-contamination ceiling* imposed by cluster size $M$.

#theorem(title: "Tight Finite-Sample Outlier Bounds and Self-Contamination Ceiling", origin: "Original Proof")[
  Let $C subset cal(D)$ be a cluster of cardinality $M = |C| >= m$ with intra-cluster distances in $[d_"min", d_"max"]$ ($0 < d_"min" <= d_"max" <= epsilon$), and let $p in cal(D) without C$ satisfy $d(p, q) in [R_"min", R_"max"]$ for all $q in C$, where $d_"max" < R_"min" <= R_"max" <= epsilon$. Assume no other point of $cal(D) without (C union {p})$ lies within distance $epsilon$ of $C union {p}$. Then:
  $
    (M R_"min") / ((M - 1) d_"max" + R_"max")
    <= "OF"_(m, epsilon) (p)
    <= (M R_"max") / ((M - 1) d_"min" + R_"min").
  $
  In particular, if $R_"min" = R_"max" = R$ and $d_"min" = d_"max" = h$ with $R > h$, then
  $
    "OF"_(m, epsilon) (p) = (M R) / ((M - 1) h + R) < M,
  $
  so for any fixed cluster size $M$, the local outlier factor is strictly bounded above by $M$ regardless of how large $R / h$ grows (provided $R <= epsilon$).
]

#proof[
  By the isolation assumption and $d_"max" < R_"min" <= R_"max" <= epsilon$, the neighborhoods are exact: $N_epsilon (p) = C$ (cardinality $M >= m > k$) and $N_epsilon (q) = (C without {q}) union {p}$ (cardinality $M$) for every $q in C$.

  *Step 1 (Core distances).* For any $q in C$, its $M - 1 >= k$ neighbors in $C without {q}$ all lie at distances in $[d_"min", d_"max"]$, strictly closer than $d(q, p) >= R_"min" > d_"max"$. Thus the $k$-th nearest neighbor of $q$ belongs to $C without {q}$, giving $d_"core" (q) in [d_"min", d_"max"]$. For $p$, all $M > k$ neighbors in $N_epsilon (p) = C$ lie in $[R_"min", R_"max"]$, so $d_"core" (p) in [R_"min", R_"max"]$.

  *Step 2 (Reachability and $"lrd"$ of $p$).* For every $q in C$, $d_"core" (q) <= d_"max" < R_"min" <= d(p, q)$, hence $d_"reach" (p, q) = max { d_"core" (q), d(p, q) } = d(p, q) in [R_"min", R_"max"]$. Averaging over $q in C$ gives
  $
    R_"min" <= "lrd"_(m, epsilon) (p) <= R_"max".
  $

  *Step 3 (Reachability and $"lrd"$ of cluster members $q in C$).* For each $q in C$, partition $N_epsilon (q)$ into $C without {q}$ ($M - 1$ points) and ${p}$ ($1$ point):
  - For every $q' in C without {q}$, both $d_"core" (q')$ and $d(q, q')$ lie in $[d_"min", d_"max"]$, so $d_"reach" (q, q') in [d_"min", d_"max"]$.
  - For $p$, both $d_"core" (p)$ and $d(q, p)$ lie in $[R_"min", R_"max"]$, so $d_"reach" (q, p) = max { d_"core" (p), d(q, p) } in [R_"min", R_"max"]$.

  Summing over the $M$ neighbors of $q$ and dividing by $M$ yields
  $
    ((M - 1) d_"min" + R_"min") / M
    <= "lrd"_(m, epsilon) (q)
    <= ((M - 1) d_"max" + R_"max") / M.
  $

  *Step 4 (Ratio aggregation).* Dividing the bounds for $"lrd"_(m, epsilon) (p)$ by the bounds for $"lrd"_(m, epsilon) (q)$ for each $q in C$ and averaging over all $M$ terms produces
  $
    (M R_"min") / ((M - 1) d_"max" + R_"max")
    &<= 1 / M sum_(q in C) ("lrd"_(m, epsilon) (p)) / ("lrd"_(m, epsilon) (q)) \
    &<= (M R_"max") / ((M - 1) d_"min" + R_"min").
  $
  Setting $d_"min" = d_"max" = h$ and $R_"min" = R_"max" = R$ collapses both sides to $(M R) / ((M - 1) h + R)$, which is strictly less than $M$ because $(M - 1) h > 0$ and approaches $M$ monotonically as $R / h -> infinity$.
]

Theorem 1 clarifies two practical regimes: when $M >> R / h$, $"OF"(p) approx R / h$ tracks the relative separation linearly; when a cluster is small ($M$ comparable to $m$), the outlier's own inclusion in $N_epsilon (q)$ inflates $"lrd"(q)$ and caps $"OF"(p) < M$.

#subsec("B", "Metric Isometry, Scale, and Permutation Invariance")

#proposition(title: "Isometry, Positive Scaling, and Permutation Invariance", origin: "Original Proof")[
  Let $sigma$ be any permutation of ${1, dots, n}$, let $T : cal(X) -> cal(X)'$ be an isometry ($d'(T(x), T(y)) = d(x, y)$), and let $tilde(d) = alpha d$ for any scalar $alpha > 0$ with rescaled radius $tilde(epsilon) = alpha epsilon$. Writing $tilde(cal(D)) = {T(x_(sigma(1))), dots, T(x_(sigma(n)))}$ under $(tilde(d), tilde(epsilon))$:
  1. $tilde(d)_"core" (T(p)) = alpha d_"core" (p)$ and $tilde("lrd")_(m, tilde(epsilon)) (T(p)) = alpha "lrd"_(m, epsilon) (p)$ for all $p in cal(D)$.
  2. The OPTICS-OF score is strictly invariant: $tilde("OF")_(m, tilde(epsilon)) (T(p)) = "OF"_(m, epsilon) (p)$.
  3. Up to tie-breaking on equal distances, the OPTICS cluster ordering $pi$ is preserved and every reachability value scales uniformly by $alpha$.
]

#proof[
  For any $p, q in cal(D)$, $tilde(d)(T(p), T(q)) = alpha d(p, q) <= alpha epsilon$ if and only if $d(p, q) <= epsilon$. Hence $N_(tilde(epsilon)) (T(p)) = T(N_epsilon (p))$ as an unordered set, independent of index permutation $sigma$. Because multiplication by $alpha > 0$ is strictly order-preserving on $[0, infinity)$, the $k$-th order statistic satisfies $tilde(d)_k (T(p)) = alpha d_k (p)$, and therefore $tilde(d)_"core" (T(p)) = alpha d_"core" (p)$.

  Using positive homogeneity of the maximum,
  $
    tilde(d)_"reach" (T(p), T(o))
    = max { alpha d_"core" (o), alpha d(p, o) }
    = alpha d_"reach" (p, o).
  $
  Averaging over $o in N_epsilon (p)$ yields $tilde("lrd")_(m, tilde(epsilon)) (T(p)) = alpha "lrd"_(m, epsilon) (p)$. Substituting into the outlier factor cancels $alpha > 0$ term by term:
  $
    tilde("OF")_(m, tilde(epsilon)) (T(p))
    &= 1 / (|N_epsilon (p)|) sum_(o in N_epsilon (p))
       (alpha "lrd"_(m, epsilon) (p)) / (alpha "lrd"_(m, epsilon) (o)) \
    &= "OF"_(m, epsilon) (p).
  $
  Finally, in Step 2 of the algorithm, every priority-queue comparison $r(u) < r(v)$ scales by $alpha > 0$, preserving strict inequalities in heap extraction and scaling the resulting reachability profile $r$ by $alpha$.
]

#subsec("C", "Density-Contrast Separation and Masking Breakdown")

#proposition(title: "Global Reachability Failure vs. OPTICS-OF and Masking Breakdown", origin: "Original Counterexample")[
  Fix $m = 4$ ($k = 3$) and scale $h > 0$.
  1. *(Multi-density separation)* There exists a 41-point Euclidean dataset $cal(D) = C_1 union C_2 union {p_1} subset bb(R)^(21)$ with generating radius $epsilon = 12 h$ such that:
     (a) No single global reachability threshold $epsilon'$ can separate the anomaly $p_1$ from both clusters $C_1$ and $C_2$ without either merging $p_1$ into $C_1$ or shattering $C_2$ into noise.
     (b) OPTICS-OF cleanly separates $p_1$ from all 40 cluster points with exact scores $"OF"(p_1) = 80 / 23 approx 3.478$, $"OF"(q) = 1543 / 1600 approx 0.964$ for all $q in C_1$, and $"OF"(q) = 1$ for all $q in C_2$.
  2. *(Adversarial masking breakdown)* For any dataset $cal(D)$ and parameters $(m, epsilon)$, inserting an adversarial clump $O = {z_1, dots, z_m}$ of $m$ points with internal diameter $delta <= epsilon$ at distance $> epsilon$ from $cal(D)$ produces $"OF"(z_j) <= d_"max" (O) / d_"min" (O)$ (identically $1$ if $O$ is equidistant). Thus the finite-sample masking breakdown point of OPTICS-OF against cohesive micro-clumps is at most $m / (n + m)$.
]

#proof[
  *Part 1.* Let $e_1, dots, e_(21)$ be the standard orthonormal basis of $bb(R)^(21)$. Define $C_1 = {(h / sqrt(2)) e_j mid(|) j = 1, dots, 20}$ (a 20-vertex regular simplex with exact pairwise Euclidean distance $h$ and centroid $c_1 = (h / (20 sqrt(2))) sum_(j=1)^(20) e_j$). Every $q in C_1$ satisfies $norm(q - c_1)_2^2 = (h^2 / 2)(1 - 1 / 20) = (19 / 40) h^2$. Place the outlier at $p_1 = c_1 + h sqrt(16 - 19 / 40) e_(21)$. By orthogonality of $e_(21)$ to $"span"(C_1)$, $norm(p_1 - q)_2 = 4 h$ for all $q in C_1$. Let $C_2$ be an isometric copy of a 20-vertex regular simplex of edge length $10 h$, translated by $100 h$ along $e_(21)$ so that $d(C_1 union {p_1}, C_2) > 50 h > epsilon = 12 h$.

  Under any global reachability cut $epsilon'$ (equivalently DBSCAN at $(epsilon', m)$): every point in $C_2$ has nearest neighbor at distance $10 h$, so retaining $C_2$ requires $epsilon' >= 10 h$; yet $p_1$ lies at distance $4 h < 10 h$ from every core point of $C_1$, so $epsilon' >= 10 h$ forces $p_1$ into $C_1$. Conversely, any $epsilon' < 4 h$ rejects $p_1$ but shatters $C_2$ into 20 noise points, proving (a).

  For OPTICS-OF with $epsilon = 12 h$ and $m = 4$:
  - Every $q in C_2$ has $N_epsilon (q) = C_2 without {q}$ with all pairwise distances $10 h$, hence $d_"core"(q) = 10 h$, $"lrd"(q) = 10 h$, and $"OF"(q) = 1$.
  - Applying Theorem 1 to $C_1$ ($M = 20$, $d_"min" = d_"max" = h$) and $p_1$ ($R_"min" = R_"max" = 4 h$) gives $"lrd"(p_1) = 4 h$, $"lrd"(q) = (19 h + 4 h) / 20 = (23 / 20) h$ for every $q in C_1$, and
    $
      "OF"(p_1) = (4 h) / ((23 / 20) h) = 80 / 23 approx 3.478.
    $
  - For every $q in C_1$, $N_epsilon (q)$ contains 19 points of $C_1$ (with $"lrd" = (23/20) h$) and $p_1$ (with $"lrd" = 4 h$), yielding
    $
      "OF"(q) &= 1 / 20 (19 dot 1 + ((23 / 20) h) / (4 h)) \
              &= 1 / 20 (19 + 23 / 80) = 1543 / 1600 approx 0.964.
    $
  Thus any threshold $tau in (1, 80 / 23)$ isolates ${p_1}$ with zero false positives or false negatives.

  *Part 2.* For the isolated clump $O$ of size $m$ at distance $> epsilon$ from $cal(D)$, every $z_j in O$ has $N_epsilon (z_j) = O without {z_j}$ of cardinality $m - 1 = k$. Applying Lemma 1 to $O$ yields $"OF"(z_j) <= d_"max"(O) / d_"min"(O)$ (and $"OF"(z_j) = 1$ if $O$ is equidistant), completely masking the $m$ injected anomalies as normal interior points. The contaminated fraction is $m / (n + m)$.
]

#subsec("D", "Complexity and Reachability-Ordering Preservation")

#theorem(title: "Joint Runtime Complexity and OPTICS Equivalence", origin: "Classical — Cited")[
  Let $T_"range" (epsilon)$ be the worst-case time of one $epsilon$-range query on $cal(D)$, and let $E_epsilon = sum_(p in cal(D)) |N_epsilon (p)| <= n (n - 1)$ be the directed edge count of the $epsilon$-neighborhood graph @ankerst1999optics @breunig1999opticsof @schubert2017dbscanrevisited. Then:
  1. Algorithm 1 computes the exact OPTICS cluster ordering $pi$, reachability profile $r$, and all $n$ outlier factors ${"OF"_(m, epsilon) (p)}$ in time
     $
       cal(O)(n T_"range" (epsilon) + E_epsilon + n log n)
     $
     using $cal(O)(n + E_epsilon)$ memory (or $cal(O)(n)$ memory via a second range-query pass).
  2. The ordering $(pi, r)$ is identical to that of OPTICS @ankerst1999optics, so any horizontal cut at $epsilon' <= epsilon$ recovers the DBSCAN density-connected clusters at $(epsilon', m)$ up to border-point assignment.
]

#proof[
  In Step 2, each object $p in cal(D)$ is marked processed once, issuing exactly one range query in total time $cal(O)(n T_"range" (epsilon))$. Finding the $k$-th smallest element in $N_epsilon (p)$ via linear-time selection (`quickselect`) takes $cal(O)(|N_epsilon (p)|)$, summing to $cal(O)(E_epsilon)$. Each directed edge $(p, q) in E_epsilon$ triggers at most one heap insertion or decrease-key on $Q$; with a binary heap this costs $cal(O)(E_epsilon log n)$ worst case (or $cal(O)(E_epsilon + n log n)$ with a Fibonacci heap), while caching adjacency lists takes $cal(O)(n + E_epsilon)$ space.

  In Steps 3 and 4, computing $"lrd"(p)$ and then $"OF"(p)$ traverses each directed edge $(p, o) in E_epsilon$ exactly twice, adding $cal(O)(E_epsilon)$ arithmetic operations without further distance evaluations. Removing caching reduces space to $cal(O)(n)$ at the expense of repeating the $n$ range queries once to accumulate $"OF"(p)$. Finally, because Steps 3--4 never alter the heap keys $r(q)$ or processing order of Step 2, the sequence $(pi, r)$ coincides identically with OPTICS @ankerst1999optics, inheriting its exact level-set correspondence with DBSCAN @ester1996dbscan.
]

#sec("III", "Complexity, Geometry and Operational Profile")

OPTICS-OF inherits the computational duality of spatial density clustering. In low-dimensional Euclidean spaces ($d <= 8$) where R\*-trees or ball trees achieve $T_"range" (epsilon) = cal(O)(log n + |N_epsilon (p)|)$ and $epsilon$ is tuned so average degree $bar("deg")_epsilon = E_epsilon / n$ remains $cal(O)(m)$, the entire joint clustering-and-anomaly pass runs in $cal(O)(n log n)$ time with modest cache footprint. Conversely, in high dimensions or under overly large $epsilon$ (chosen to guarantee $|N_epsilon (p)| >= m - 1$ for sparse peripheral points), $E_epsilon$ can approach $Theta(n^2)$, collapsing spatial index pruning and forcing an $cal(O)(n^2)$ scan @schubert2017dbscanrevisited @beyer1999nearest.

Geometrically, OPTICS-OF operates purely on the metric graph $G_epsilon$ rather than coordinate projections, enabling arbitrary non-convex manifolds, non-Euclidean metrics, and multi-density regimes without parametric cluster-shape assumptions. The critical engineering trade-off lies in memory versus recomputation: caching $E_epsilon$ supports single-query exact evaluation of both $"lrd"$ and $"OF"$ in $cal(O)(n + E_epsilon)$ RAM, whereas streaming implementations bound memory at $cal(O)(n)$ by maintaining a two-phase deferred evaluation buffer @breunig1999opticsof.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time complexity], [$cal(O)(n log n)$ indexed (low $d$, sparse $G_epsilon$); $cal(O)(n^2)$ worst-case or dense $epsilon$],
    [Space complexity], [$cal(O)(n + E_epsilon)$ with graph cache; $cal(O)(n)$ with two-pass range queries],
    [Hyperparameters], [Generating radius $epsilon > 0$, $m = "MinPts" >= 2$, outlier cutoff $tau > 1$, optional $xi$],
    [Cluster geometry], [Arbitrary non-convex shapes and multi-scale local densities via reachability plot],
    [Symmetries], [Exact isometry, permutation, and joint $(d, epsilon)$ positive-scale invariance],
    [Outlier semantics], [Continuous local density-contrast ratio $"OF"(p) in (0, infinity)$ centered near $1$],
    [Main bottleneck], [Choosing $epsilon$ large enough for sparse points without making $E_epsilon = Theta(n^2)$],
  ),
  caption: [Operational and computational profile of OPTICS-OF.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Financial Fraud and Anti-Money-Laundering Analytics*: Transaction embeddings across retail, corporate, and institutional cohorts exhibit order-of-magnitude density disparities. Because OPTICS-OF normalizes each point's reachability distance against its neighbors' $"lrd"$, an anomalous mid-sized transfer orbiting a dense retail cluster achieves a high $"OF"$ score without falsely flagging legitimate sparse high-net-worth clusters (Proposition 2).
2. *Network Intrusion and Cybersecurity Telemetry*: Enterprise packet-flow features mix high-volume routine web protocols (tight clusters) with sparse administrative or backup traffic. Joint density ordering segments protocol regimes while continuous $"OF"$ ranking prioritizes stealthy probing or exfiltration bursts that sit well inside global inter-cluster distances.
3. *Astronomical Survey Outlier Discovery*: In multi-band photometric catalogs (e.g., Gaia / SDSS color-magnitude diagrams), main-sequence stars form dense curved loci while giant branches and halo streams form sparse filaments. OPTICS-OF's non-convex reachability ordering traces curved stellar populations while its scale-invariant ratio flags rare quasars, white-dwarf binaries, and instrument artifacts locally.
4. *Industrial IoT and Rotating-Machinery Condition Monitoring*: Turbines and compressors operate across distinct load setpoints (idle, cruise, peak), each producing a cluster of vibration-spectrum embeddings with different sensor variance. Computing $"OF"$ alongside mode clusters detects incipient bearing wear departing from the active setpoint long before the vibration exceeds global operating bounds.
5. *Geospatial Mobility and Maritime AIS Screening*: Vessel trajectories and urban GPS dwell points have extreme spatial density contrasts between dense ports / city centers and sparse open-sea corridors. Metric-space generality permits geodesic or trajectory distances while ranking localized course deviations relative to surrounding corridor density.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

OPTICS-OF degrades in four identifiable regimes: (i) *high-dimensional distance concentration* @beyer1999nearest, where pairwise distances collapse around their mean ($d_"max" / d_"min" -> 1$), compressing $"OF"(p) -> 1$ by Lemma 1 and rendering $epsilon$ selection brittle; (ii) *cohesive anomaly masking*, where $m$ or more coordinated anomalies form a micro-cluster with $"OF" approx 1$ (Proposition 2); (iii) *unbalanced $epsilon$-radius explosion*, where selecting $epsilon$ large enough to connect sparse regions makes dense clusters nearly complete graphs ($E_epsilon -> n^2$); and (iv) *extreme cluster-size skew*, where small legitimate clusters ($|C| approx m$) cap the outlier factor of nearby anomalies via the Theorem 1 self-contamination ceiling $"OF"(p) < |C|$.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Multi-density fraud & cyber telemetry], [High: local $"lrd"$ ratio isolates anomalies near dense cohorts], [Masking if $m$ bots clump; tune $m$ or use HDBSCAN],
    [Low-D spatial / GIS & astronomy], [High: arbitrary manifolds + $cal(O)(n log n)$ spatial indexing], [Watch $epsilon$ blow-up; consider adaptive $k$-NN LOF],
    [Multi-regime sensor monitoring], [High: separates operating modes and local drift jointly], [Rescale heterogeneous sensors before $d$],
    [High-D embeddings ($d > 20$)], [Poor: distance concentration drives $"OF" -> 1$, $E_epsilon -> n^2$], [Use PCA / UMAP first or subspace / SNN methods],
    [Massive streaming data ($n > 10^6$)], [Low: two-pass $G_epsilon$ and $cal(O)(n^2)$ worst-case scaling], [Use DenStream or approximate $k$-NN LOF],
  ),
  caption: [Domain suitability matrix for OPTICS-OF.],
)

#sec("V", "Structural Relations and Empirical Behavior")

OPTICS-OF occupies the pivotal historical bridge between density-based cluster ordering (OPTICS @ankerst1999optics, Alg. 43) and dedicated local anomaly detection (LOF @breunig1999opticsof). We formalize its exact algebraic connection to LOF and hierarchical density-based outlier scoring.

#proposition(title: "Exact Reduction to LOF Under Adaptive k-Distance Cutoffs", origin: "Original Derivation")[
  Let $N_k (p) = {q in cal(D) without {p} mid(|) d(p, q) <= d_k (p)}$ denote the adaptive $k$-distance neighborhood of $p$ with $k = m - 1$. If the fixed global generating radius $epsilon$ in OPTICS-OF is replaced point-wise by $epsilon(p) = d_k (p)$, then defining the local reachability density $rho_k (p) = ("lrd"_(m, d_k (p)) (p))^(-1)$ yields the exact identity
  $
    "OF"_(m, d_k) (p)
    &= 1 / (|N_k (p)|) sum_(o in N_k (p)) ("lrd"_(m, d_k (p)) (p)) / ("lrd"_(m, d_k (o)) (o)) \
    &= 1 / (|N_k (p)|) sum_(o in N_k (p)) (rho_k (o)) / (rho_k (p))
     = "LOF"_k (p).
  $
]

#proof[
  With $epsilon(p) = d_k (p)$, every $p in cal(D)$ automatically satisfies $|N_(epsilon(p)) (p)| >= k$, so $d_"core" (o) = d_k (o)$ is finite for every $o in cal(D)$ and $d_"reach" (p, o) = max { d_k (o), d(p, o) }$ coincides identically with the $k$-reachability distance of LOF. Substituting $"lrd"(u) = 1 / rho_k (u)$ into every term of the average gives $"lrd"(p) / "lrd"(o) = rho_k (o) / rho_k (p)$, establishing exact algebraic equivalence.
]

Proposition 3 explains why OPTICS-OF inspired both LOF (which decouples outlier scoring from a global $epsilon$ via adaptive $k$-NN neighborhoods) and later hierarchical cluster-outlier frameworks such as HDBSCAN\* / GLOSH @campello2013hdbscan @campello2015hierarchical, which compare each point's density against the peak density of its enclosing cluster component rather than its immediate $epsilon$-neighbors. Relative to GDBSCAN @sander1998gdbscan, VDBSCAN @liu2007vdbscan, SNN @ertoz2003snn, and RNN-DBSCAN @bryant2018rnndbscan, OPTICS-OF remains distinctive in computing a hierarchical reachability ordering and a calibrated local anomaly index in a single unified pass @aggarwal2013data @jain2010data.

Empirically, on the ClusterOverview 4-dataset benchmark evaluated via Silhouette @rousseeuw1987silhouettes and Adjusted Rand Index (ARI) @hubert1985comparing, OPTICS-OF achieves strong low-dimensional separation on Iris (`Sil = 0.616`, `ARI = 0.481`) with a fast mean runtime of `4.412 ms`. The moderate ARI relative to the high Silhouette reflects a characteristic property of density reachability on Iris: the well-separated _Setosa_ cluster is cleanly isolated while the overlapping _Versicolor_ and _Virginica_ clouds merge under a single global valley threshold unless split by fine-grained $xi$-steepness. On the Customers dataset (`Sil = 0.281`), heterogeneous feature densities and peripheral spenders yield a cautious density partition where outlier flagging trims low-density tails at the cost of lower convex compactness. On Multi-Shape (`Sil = 0.446`, `ARI = N/A`), OPTICS-OF successfully tracks non-convex geometry where centroid-based partitioners struggle, although Euclidean Silhouette naturally penalizes elongated curved components compared to compact spherical blobs. Finally, on the High-D benchmark (`Sil = N/A`, `ARI = N/A`), OPTICS-OF fails to produce a non-trivial multi-cluster partition, directly corroborating the theoretical distance-concentration regime of Lemma 1 and @beyer1999nearest, in which a fixed global $epsilon$ either leaves nearly all points as isolated non-core objects or collapses the entire dataset into a single connected component.
