#import "../template/ieee.typ": *

#alg_entry(
  id: 136,
  name: "PROCLUS (Projected Clustering)",
  category: "Cat. VI: Grid, Subspace & Projection",
  year: "1999",
  paradigm: "Medoid-Based Disjoint Projected Subspace Clustering",
  complexity: $cal(O)(I n k d + k d log(k d))$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Bottom-up grid subspace methods such as CLIQUE @agrawal1998clique suffer an unavoidable $Omega(2^k)$ complexity barrier when intrinsic cluster dimensionality $k$ is moderate-to-large and produce overlapping collections of lower-dimensional shadow clusters. PROCLUS (PROjected CLUStering) @aggarwal1999proclus resolves both bottlenecks by seeking a *disjoint partition* of $cal(X) = {bold(x)_1, dots.h, bold(x)_n} subset RR^d$ into $k$ projected clusters plus an explicit outlier set $cal(O)$, associating each cluster $C_i$ with a representative medoid $bold(m)_i in cal(X)$ and a dedicated subset of relevant coordinates $cal(D)_i subset.eq [d] := {1, dots.h, d}$.

#definition(title: "Projected k-Medoids Partition and Segmental Manhattan Dispersion")[
  Given cluster count $k in NN_(>= 2)$ and target average subspace dimensionality $l in [2, d]$, a *projected $k$-clustering* consists of medoids $cal(M) = {bold(m)_1, dots.h, bold(m)_k} subset cal(X)$, coordinate subsets $cal(D)_1, dots.h, cal(D)_k subset.eq [d]$ satisfying the cardinality constraints
  $ |cal(D)_i| >= 2 quad forall i in [k], quad sum_(i=1)^k |cal(D)_i| = k l, $
  and a disjoint partition $cal(X) = C_1 union.sq dots.c union.sq C_k union.sq cal(O)$. The *Manhattan segmental distance* between $bold(x), bold(y) in RR^d$ relative to subspace $cal(D)_i$ is the dimension-normalized $L_1$ metric
  $ d_(cal(D)_i) (bold(x), bold(y)) = 1 / (|cal(D)_i|) sum_(j in cal(D)_i) |x_j - y_j|. $
  Assigning each non-outlier $bold(x)$ to $i^*(bold(x)) = arg min_(1 <= i <= k) d_(cal(D)_i) (bold(x), bold(m)_i)$, PROCLUS evaluates partition quality via the weighted average segmental dispersion objective
  $ J(cal(M), {cal(D)_i}) = 1 / n sum_(i=1)^k |C_i| w_i, $
  where $w_i = (|cal(D)_i| |C_i|)^(-1) sum_(j in cal(D)_i) sum_(bold(x) in C_i) |x_j - m_(i,j)|$ is the average intra-cluster segmental deviation of cluster $i$.
]

To locate $(cal(M), {cal(D)_i})$ without exponential subspace enumeration, PROCLUS executes a three-phase medoid-hill-climbing strategy inspired by CLARANS @ng2002clarans and PAM @kaufman1990finding:

1. *Initialization Phase:* Draw a uniform random sample $S subset cal(X)$ of size $A k$ (typically $A = 50$–$100$), then apply a farthest-point greedy dispersion pass on $S$ to extract a well-separated candidate pool $M subset S$ of size $B k$ (typically $B = 10$). Select $k$ initial medoids $M_"curr" subset M$.
2. *Iterative Phase:* For each medoid $bold(m)_i in M_"curr"$, define its *locality* $cal(L)_i = {bold(x) in cal(X) : norm(bold(x) - bold(m)_i)_1 <= delta_i}$, where $delta_i = min_(r != i) norm(bold(m)_i - bold(m)_r)_1$. Compute per-coordinate average deviations $X_(i,j) = |cal(L)_i|^(-1) sum_(bold(x) in cal(L)_i) |x_j - m_(i,j)|$, row mean $Y_i = d^(-1) sum_(j=1)^d X_(i,j)$, row standard deviation $sigma_i = sqrt(d^(-1) sum_(j=1)^d (X_(i,j) - Y_i)^2)$, and standardized dispersion scores
   $ Z_(i,j) = (X_(i,j) - Y_i) / sigma_i. $
   Select coordinate sets ${cal(D)_i}$ minimizing $sum_(i=1)^k sum_(j in cal(D)_i) Z_(i,j)$ subject to $|cal(D)_i| >= 2$ and $sum_(i=1)^k |cal(D)_i| = k l$. Assign all points via $d_(cal(D)_i)$, evaluate $J$, retain the lowest-$J$ medoid set $M_"best"$, and replace "bad" medoids (the smallest cluster and any cluster with $|C_i| < "min\_dev" dot.op (n / k)$) with random unused candidates from $M without M_"curr"$.
3. *Refinement Phase:* Using $M_"best"$ and its final clusters $C_i$ (instead of full-space localities $cal(L)_i$), recompute $X_(i,j), Z_(i,j)$ and refined subspaces $cal(D)_i$, reassign points via $d_(cal(D)_i)$, and label $bold(x)$ as an outlier ($bold(x) in cal(O)$) if $d_(cal(D)_i) (bold(x), bold(m)_i) > Delta_i$ for all $i in [k]$, where $Delta_i = min_(r != i) d_(cal(D)_i) (bold(m)_i, bold(m)_r)$ is medoid $i$'s *sphere of influence*.

#algbox(
  title: "PROCLUS Three-Phase Projected Clustering",
  inputs: [Dataset $cal(X) subset RR^d$, cluster count $k$, average subspace dimension $l in [2, d]$, pool multipliers $A > B > 1$, minimum deviation ratio $"min\_dev"$.],
  outputs: [Projected clusters ${(C_i, bold(m)_i, cal(D)_i)}_(i=1)^k$ and outlier set $cal(O)$.],
)[
  + *Initialization:* Uniformly sample $S subset cal(X)$ ($|S| = A k$); greedily pick $M subset S$ ($|M| = B k$) maximizing minimum distance to previously selected points; initialize $M_"curr" subset M$ ($|M_"curr"| = k$), $J_"best" <- infinity$.
  + *Iterative Medoid & Subspace Search:* Repeat until $J_"best"$ fails to improve for $T_"max"$ iterations:
    - Form localities $cal(L)_i$ around $bold(m)_i in M_"curr"$ with full-$d$ radius $delta_i = min_(r != i) norm(bold(m)_i - bold(m)_r)_1$.
    - Compute $X_(i,j)$ over $cal(L)_i$ and standardized scores $Z_(i,j) = (X_(i,j) - Y_i) / sigma_i$.
    - *Greedy Subspace ILP:* Pick the 2 smallest $Z_(i,j)$ per cluster $i$, then the $k(l - 2)$ globally smallest remaining $Z_(i,j)$ pairs to define ${cal(D)_i}$.
    - Assign each $bold(x) in cal(X)$ to $arg min_i d_(cal(D)_i) (bold(x), bold(m)_i)$ and compute $J$.
    - If $J < J_"best"$, update $J_"best" <- J$, $M_"best" <- M_"curr"$, and identify bad medoids; else revert to $M_"best"$ and randomly replace bad medoids from $M without M_"best"$.
  + *Refinement & Outlier Pruning:* Recompute $Z_(i,j)$ and ${cal(D)_i}$ from the actual clusters $C_i$ of $M_"best"$; reassign points and move $bold(x)$ to $cal(O)$ if $d_(cal(D)_i) (bold(x), bold(m)_i) > Delta_i$ for all $i$.
  + *Return* ${(C_i, bold(m)_i, cal(D)_i)}_(i=1)^k$ and $cal(O)$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five formal mathematical results for PROCLUS: exact global optimality of its greedy subspace-selection subroutine, probabilistic medoid-pool coverage bounds, exact scale-invariant z-score separation and finite-sample subspace recovery, constructive counterexamples proving ambient-locality pollution and objective non-monotonicity, and the mathematical necessity of $|cal(D)_i|$-normalization.

#subsec("A", "Exact Optimality of Greedy Z-Score Subspace Allocation")

#theorem(title: "Exact ILP Optimality of PROCLUS Dimension Selection", origin: "Original Proof")[
  For any real score matrix $Z in RR^(k times d)$ and integer budget $k l$ with $2 <= l <= d$, consider the 0-1 integer linear program
  $ min_(u in {0, 1}^(k times d)) sum_(i=1)^k sum_(j=1)^d Z_(i,j) u_(i,j) \
    "s.t." quad sum_(j=1)^d u_(i,j) >= 2 quad (forall i in [k]), quad sum_(i=1)^k sum_(j=1)^d u_(i,j) = k l. $
  Then: (i) the constraint matrix is totally unimodular (so the LP relaxation is exact), and (ii) PROCLUS's two-stage greedy algorithm — selecting the two smallest entries $Z_(i,j)$ in each row $i$, followed by the $k(l - 2)$ globally smallest entries among all remaining pairs — computes an *exact global minimizer* in $cal(O)(k d log(k d))$ time.
]

#proof[
  For (i), each variable $u_(i,j)$ appears in exactly one row constraint ($sum_j u_(i,j) >= 2$) with coefficient $+1$ and in the single global cardinality constraint ($sum_(i,j) u_(i,j) = k l$) with coefficient $+1$. Partitioning the constraint rows into the $k$ disjoint row-sum inequalities and the single global equality yields a bipartite vertex-edge incidence structure; by Heller-Tompkins's theorem, every matrix with at most two non-zero entries per column of identical sign partitioned into two disjoint row groups is totally unimodular. Since all right-hand sides ($2$, $k l$, and box bounds $0 <= u_(i,j) <= 1$) are integers, every extreme point of the LP relaxation is integral.

  For (ii), let $u^"gr"$ be the indicator matrix produced by the greedy rule, decomposed as $u^"gr" = u^((1)) + u^((2))$ where $u^((1))$ marks the 2 smallest entries per row and $u^((2))$ marks the $k(l - 2)$ smallest remaining entries globally. Let $u^*$ be any feasible optimal solution minimizing the objective, and assume $u^*$ maximizes overlap $|{(i,j) : u^*_(i,j) = u^"gr"_(i,j) = 1}|$ among all global minimizers.

  First, suppose $u^*_(i, j_1) = 0$ for some entry $(i, j_1)$ selected in Stage 1 ($u^((1))_(i, j_1) = 1$). Because row $i$ satisfies $sum_j u^*_(i,j) >= 2$ and $u^((1))$ selects only the two smallest entries in row $i$, there exists $j_2$ in row $i$ with $u^*_(i, j_2) = 1$, $u^((1))_(i, j_2) = 0$, and $Z_(i, j_1) <= Z_(i, j_2)$. Swapping $u^*_(i, j_2) |-> 0$ and $u^*_(i, j_1) |-> 1$ preserves all row sums and the global sum $k l$, does not increase the objective ($Z_(i, j_1) - Z_(i, j_2) <= 0$), and strictly increases agreement with $u^"gr"$, contradicting maximality of overlap unless $u^* >= u^((1))$ entrywise.

  Second, with $u^* >= u^((1))$ established, both $v^* = u^* - u^((1))$ and $u^((2))$ are binary indicator vectors supported on $cal(E) = {(i,j) : u^((1))_(i,j) = 0}$ with identical cardinality $sum_((i,j) in cal(E)) v^*_(i,j) = k(l - 2)$. Because every $v$ supported on $cal(E)$ of cardinality $k(l - 2)$ automatically preserves $sum_j (u^((1)) + v)_(i,j) >= 2$, minimizing $sum_((i,j) in cal(E)) Z_(i,j) v_(i,j)$ over all subsets of size $k(l - 2)$ is uniquely solved by picking the $k(l - 2)$ smallest values $Z_(i,j)$ in $cal(E)$, which is precisely $u^((2))$. Hence $u^"gr"$ achieves the exact global minimum. Sorting all $k d$ entries takes $cal(O)(k d log(k d))$ time (or $cal(O)(k d)$ via linear-time selection).
]

#subsec("B", "Probabilistic Medoid-Pool Coverage Guarantee")

#theorem(title: "Sample-Size Bound for Complete Cluster Representation", origin: "Classical — Cited")[
  Suppose the dataset $cal(X)$ contains $k$ true clusters $C_1^*, dots.h, C_k^*$ with minimum cluster mass fraction $p_"min" = min_(1 <= i <= k) |C_i^*| / n > 0$. If a uniform random sample $S$ of size $m_S = A k$ is drawn with replacement (or without replacement from $n >> A k$), then $S$ contains at least one representative from every true cluster $C_1^*, dots.h, C_k^*$ with probability at least $1 - delta$ whenever @aggarwal1999proclus
  $ A k >= 1 / p_"min" ln(k / delta). $
  In particular, for balanced clusters ($p_"min" = 1 / k$), $A >= ln(k / delta)$ suffices. Moreover, greedy farthest-point extraction of $B k$ points from $S$ is a $2$-approximation to the $k$-center radius on $S$.
]

#proof[
  Let $E_i$ denote the bad event that sample $S$ misses true cluster $C_i^*$. Under independent uniform draws,
  $ "Pr"(E_i) = (1 - p_i)^(A k) <= (1 - p_"min")^(A k) <= exp(-p_"min" A k), $
  using the elementary inequality $1 - t <= e^(-t)$ for all $t in RR$. (Sampling without replacement yields a strictly smaller failure probability by negative association.) Applying the union bound across all $k$ clusters gives
  $ "Pr"(union.big_(i=1)^k E_i) <= sum_(i=1)^k "Pr"(E_i) <= k exp(-p_"min" A k). $
  Setting $k exp(-p_"min" A k) <= delta$ and taking logarithms yields $p_"min" A k >= ln(k / delta)$, proving the bound.

  For the farthest-point greedy pass choosing $M = {s_1, dots.h, s_(B k)} subset S$, let $r^*$ be the optimal $(B k)$-center radius covering $S$ with $B k$ balls of radius $r^*$. If any point $x in S$ satisfies $min_(s in M) norm(x - s)_1 > 2 r^*$, then every subsequent greedy pick before exhausting $B k$ points was also at distance $> 2 r^*$ from all earlier picks. By the triangle inequality, no two of the $B k + 1$ points $M union {x}$ can lie inside the same optimal ball of radius $r^*$, contradicting the pigeonhole principle. Hence every point of $S$ lies within distance $2 r^*$ of $M$.
]

#subsec("C", "Scale-Invariant Z-Score Separation and Finite-Sample Recovery")

#theorem(title: "Exact Population Z-Scores and High-Probability Dimension Recovery", origin: "Original Derivation")[
  Fix a cluster $C_i$ with true relevant subspace $cal(D)_i^*$ of cardinality $|cal(D)_i^*| = l < d$, and suppose observations in $C_i$ have coordinates in $[0, R]$ with expected absolute deviations $EE|x_j - m_(i,j)| = sigma_"in"$ for $j in cal(D)_i^*$ and $EE|x_j - m_(i,j)| = sigma_"out" > sigma_"in"$ for $j in.not cal(D)_i^*$.
  + *Scale-Invariant Exact Population Scores:* The population standardized scores depend *only* on the dimension ratio $l / d$, independent of the absolute gap $Delta sigma := sigma_"out" - sigma_"in"$:
    $ Z_(i,j)^* = cases(
      -sqrt((d - l) / l) & "if" j in cal(D)_i^*,
      +sqrt(l / (d - l)) & "if" j in.not cal(D)_i^*.
    ) $
  + *Finite-Sample Recovery:* If locality $cal(L)_i$ consists of $N_i$ independent samples of which fraction $eta_i in (0, 1]$ come from $C_i$ and $1 - eta_i$ come from background noise with expected deviation $sigma_"out"$ on all coordinates, then choosing the $l$ smallest empirical scores $Z_(i,j)$ recovers $cal(D)_i^*$ exactly with probability at least $1 - delta$ whenever
    $ N_i >= (2 R^2) / (eta_i^2 (Delta sigma)^2) ln((2 d) / delta). $
]

#proof[
  For (1), compute the population mean across the $d$ coordinates:
  $ Y_i^* = 1 / d (l sigma_"in" + (d - l) sigma_"out") = sigma_"out" - l / d Delta sigma. $
  The coordinate deviations from $Y_i^*$ are $sigma_"in" - Y_i^* = -((d - l) / d) Delta sigma$ for $j in cal(D)_i^*$ and $sigma_"out" - Y_i^* = (l / d) Delta sigma$ for $j in.not cal(D)_i^*$. Hence the population variance across coordinates is
  $ (sigma_i^*)^2 &= 1 / d (l ((d - l) / d)^2 + (d - l) (l / d)^2) (Delta sigma)^2 \
    &= (l (d - l)) / d^2 (Delta sigma)^2. $
  Taking the square root gives $sigma_i^* = (sqrt(l (d - l)) / d) Delta sigma$. Dividing each centered deviation by $sigma_i^*$ cancels $Delta sigma$ and $d$ completely, yielding $Z_(i,j)^* = -sqrt((d - l) / l)$ for $j in cal(D)_i^*$ and $+sqrt(l / (d - l))$ for $j in.not cal(D)_i^*$.

  For (2), under mixture purity $eta_i$, the expected empirical deviation on coordinate $j$ satisfies $EE[X_(i,j)] = eta_i sigma_"in" + (1 - eta_i) sigma_"out" = sigma_"out" - eta_i Delta sigma$ for $j in cal(D)_i^*$, versus $EE[X_(i,j)] = sigma_"out"$ for $j in.not cal(D)_i^*$, creating a separation gap $g = eta_i Delta sigma$. Because $Z_(i,j) = (X_(i,j) - Y_i) / sigma_i$ is a strictly increasing affine function of $X_(i,j)$ within row $i$, the ordering of $Z_(i,j)$ coincides identically with the ordering of $X_(i,j)$. Exact recovery of $cal(D)_i^*$ occurs whenever $max_(j in cal(D)_i^*) X_(i,j) < min_(j in.not cal(D)_i^*) X_(i,j)$, which is implied by $|X_(i,j) - EE[X_(i,j)]| < g / 2$ for every $j in [d]$. Since each $|x_j - m_(i,j)| in [0, R]$, Hoeffding's inequality and a union bound across $d$ coordinates give
  $ "Pr"(exists j : |X_(i,j) - EE[X_(i,j)]| >= g / 2) <= 2 d exp(-(N_i g^2) / (2 R^2)). $
  Bounding the right-hand side by $delta$ and solving for $N_i$ yields $N_i >= (2 R^2 / (eta_i^2 (Delta sigma)^2)) ln(2 d / delta)$.
]

#subsec("D", "Ambient-Locality Pollution and Objective Non-Monotonicity")

#theorem(title: "Locality Pollution and Non-Descent Counterexamples", origin: "Original Counterexample")[
  + *Ambient-Locality Concentration Failure:* Consider two equal-sized projected clusters $C_1^*, C_2^*$ in $[0, 1]^d$ whose relevant $l$-dimensional subspaces are disjoint ($l$ fixed) and whose remaining $d - l$ coordinates are i.i.d. $"Uniform"[0, 1]$. As $d -> infinity$, full-dimensional $L_1$ distances between medoids and points concentrate at $d / 3$ with fluctuations $cal(O)(sqrt(d))$, causing the purity $eta_i$ of PROCLUS's full-$d$ locality $cal(L)_i$ to converge to $1 / 2$ and inflating the sample complexity of Theorem II.C by factor $4$ (or collapsing $eta_i -> 1 / k$ for $k$ clusters).
  + *Non-Monotonicity of $J$:* Unlike Lloyd's $k$-means @lloyd1982least, replacing subspaces ${cal(D)_i}$ by the exact minimizer of the z-score objective $sum_(i,j) Z_(i,j) u_(i,j)$ (Theorem II.A) can strictly *increase* PROCLUS's actual dispersion objective $J(cal(M), {cal(D)_i})$.
]

#proof[
  For (1), let $bold(m)_1 in C_1^*$ have i.i.d. $"Uniform"[0, 1]$ coordinates on its $d - l$ irrelevant axes. For any independent pair $U_1, U_2 tilde.op "Uniform"[0, 1]$, $EE|U_1 - U_2| = 1 / 3$ and $"Var"(|U_1 - U_2|) = 1 / 18$. For $bold(x) in C_1^*$ (zero deviation on its $l$ relevant axes) versus $bold(x)' in C_2^*$ (zero deviation on $l$ disjoint axes), full-$d$ $L_1$ distances obey
  $ norm(bold(x) - bold(m)_1)_1 = S_(d-l), quad norm(bold(x)' - bold(m)_1)_1 = S_(d-2l) + O(l), $
  where $S_m$ is a sum of $m$ noise terms concentrating at $m / 3$ with standard deviation $Theta(sqrt(m))$. Because the $O(l)$ subspace signal is dwarfed by $Theta(sqrt(d))$ noise fluctuations as $d -> infinity$, $"Pr"(norm(bold(x)' - bold(m)_1)_1 < norm(bold(x) - bold(m)_1)_1) -> 1 / 2$ @beyer1999nearest. For $k$ clusters in fixed-$l$ subspaces, the full-$d$ locality $cal(L)_i$ draws nearly uniformly across clusters ($eta_i -> 1 / k$), inflating Theorem II.C's sample bound by factor $k^2$.

  For (2), take a 2-cluster, 4-dimensional instance with $|C_1| = |C_2|$, fixed medoids, and budget $k l = 5$, having deviation rows $X_1 = (0, 0, 2, 100)$ ($Y_1 = 25.5, sigma_1 = 43.02 => Z_(1,3) = -0.546$) and $X_2 = (0, 0, 0.5, 0.6)$ ($Y_2 = 0.275, sigma_2 = 0.277 => Z_(2,3) = +0.812$). After Stage 1 picks ${1, 2}$ in both rows, minimizing $sum Z_(i,j)$ strictly assigns the fifth coordinate to $(1, 3)$ ($Z_(1,3) < Z_(2,3)$), yielding $J = (1 / 2)(2 / 3 + 0) = 1 / 3$. Allocating it to $(2, 3)$ instead yields $J = (1 / 2)(0 + 1 / 6) = 1 / 12 < 1 / 3$. Hence z-score optimization can strictly increase $J$, necessitating best-iterate retention $M_"best"$.
]

#subsec("E", "Necessity of Dimension-Cardinality Normalization")

#proposition(title: "Dimension-Homogeneity of Segmental Manhattan Distance", origin: "Original Proof")[
  Let $bold(x) in RR^d$ have i.i.d. coordinate residuals $|x_j - m_(i,j)|$ with common mean $mu > 0$ across two candidate clusters $C_a, C_b$ having unequal subspace cardinalities $|cal(D)_a| < |cal(D)_b|$. Under unnormalized $L_1$ distance, $EE[norm(pi_(cal(D)_a)(bold(x) - bold(m)_a))_1] < EE[norm(pi_(cal(D)_b)(bold(x) - bold(m)_b))_1]$ strictly, systematically biasing assignments toward lower-dimensional clusters. Under PROCLUS's segmental distance $d_(cal(D)_i)$, expectations are strictly equalized: $EE[d_(cal(D)_a)] = EE[d_(cal(D)_b)] = mu$.
]

#proof[
  By linearity of expectation, unnormalized distance obeys $EE[sum_(j in cal(D)_i) |x_j - m_(i,j)|] = |cal(D)_i| mu$, creating a strict distortion ratio $|cal(D)_b| / |cal(D)_a| > 1$ that causes smaller-subspace clusters to absorb points belonging to higher-dimensional clusters even when per-axis fit is worse. Dividing by $|cal(D)_i|$ yields $EE[d_(cal(D)_i) (bold(x), bold(m)_i)] = |cal(D)_i|^(-1) (|cal(D)_i| mu) = mu$, restoring dimension-homogeneous Voronoi competition across heterogeneous subspace cardinalities.
]

#sec("III", "Complexity, Geometry and Operational Profile")

PROCLUS eliminates CLIQUE's exponential dependence on subspace dimension $k$: initialization requires $cal(O)(A B k^2 d)$ distance evaluations on the $A k$ sample; each iteration computes $k$ localities in $cal(O)(n k d)$ time, solves the greedy subspace ILP in $cal(O)(k d log(k d))$ time (Theorem II.A), and assigns all $n$ points in $cal(O)(n k d)$ time. Across $I$ hill-climbing iterations plus one $cal(O)(n k d)$ refinement pass, total runtime is $cal(O)(I n k d + I k d log(k d) + A B k^2 d)$ — *linear in both $n$ and $d$*! Space complexity is $cal(O)(n d + B k d)$ (or $cal(O)(B k d)$ streaming if data is scanned from disk).

Geometrically, each projected cluster $C_i$ forms a convex $L_1$ hyper-octahedron (cross-polytope) inside its dedicated subspace $cal(D)_i$, extruded as an infinite cylinder across irrelevant axes $[d] without cal(D)_i$, truncated by the sphere-of-influence radius $Delta_i$ for outlier rejection.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [*Property*], [*Characterization*],
    [Time Complexity], [$cal(O)(I n k d + I k d log(k d) + A B k^2 d)$: strictly linear in $n$ and $d$],
    [Space Complexity], [$cal(O)(n d)$ in-memory or $cal(O)(B k d)$ disk-streaming medoid pool],
    [Partition Type], [Disjoint $k$-way partition plus explicit outlier set $cal(O)$],
    [Cluster Geometry], [Convex $L_1$-balls within cluster-specific axis-aligned subspaces $cal(D)_i$],
    [Determinism], [Stochastic (random sample $S$ and bad-medoid swaps); tracks best $J$],
    [Outlier Handling], [Explicit via sphere-of-influence radii $Delta_i = min_(r != i) d_(cal(D)_i)(bold(m)_i, bold(m)_r)$],
    [Primary Hyperparams], [Cluster count $k$, average subspace dimension $l in [2, d]$],
    [Key Limitations], [Requires specifying $k, l$; full-$d$ locality pollution; non-convex blindness],
  ),
  caption: [Computational, geometric, and operational profile of PROCLUS @aggarwal1999proclus.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

*1. Single-Cell and Microarray Transcriptomics Subtyping:* Biological assays measure $d approx 10^3$–$10^4$ genes per cell, where each phenotype coheres along a specific subset of $l approx 20$–$100$ marker genes while others fluctuate noisily. PROCLUS's linear-in-$d$ cost discovers disjoint subtypes alongside their marker subsets, paralleling biclustering @cheng2000biclustering.

*2. E-Commerce Customer Preference Profiling:* Retail tables span hundreds of product categories where shopper personas exhibit tight spending inside small category subsets while varying widely elsewhere. Disjoint projected clustering maps every user to one actionable persona without CLIQUE's overlapping explosion @aggarwal2013data.

*3. Financial Asset Regime and Factor Partitioning:* Multi-asset portfolios track hundreds of fundamental and volatility descriptors reacting to distinct factor subsets. PROCLUS partitions securities into cohesive factor groups while isolating idiosyncratic outliers via the sphere-of-influence radius $Delta_i$.

*4. Industrial Multi-Sensor Fault Isolation:* Semiconductor fabrication and turbine telemetry monitor $d > 100$ sensors where distinct fault modes localize tightly within specific sensor subsets; PROCLUS identifies both the fault cluster and the exact implicated subset $cal(D)_i$ in polynomial time.

*5. High-Dimensional Document and Embedding Partitioning:* High-dimensional text representations concentrate topics within distinct feature subsets; PROCLUS's $L_1$ segmental metric resists single-coordinate outliers better than squared Euclidean $k$-means while avoiding $Omega(2^k)$ grid enumeration.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

PROCLUS degrades in four regimes: (i) *extreme noise dimensionality ($d >> k l$)*, where full-space localities $cal(L)_i$ suffer ambient distance concentration (Theorem II.D(1)); (ii) *arbitrarily oriented (correlated) subspaces*, where axis-aligned coordinates fail to align with oblique hyperplanes; (iii) *non-convex or highly imbalanced clusters*, where $L_1$-medoids split elongated manifolds and the $"min\_dev"$ heuristic prematurely discards small valid clusters; and (iv) *misspecified $(k, l)$*, since the global budget $sum |cal(D)_i| = k l$ forces irrelevant axes into clusters if $l$ is overestimated.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [*Domain*], [*Mathematical Fit*], [*Caution / Alternative*],
    [High-D Tabular & Omics ($d >= 20$)], [Linear in $n, d$; disjoint clusters + explicit feature sets $cal(D)_i$], [Tune $(k, l)$; check locality purity via FINDIT @woo2004findit],
    [Oblique / Correlated Subspaces], [Poor: axis-aligned $cal(D)_i$ misses rotated manifolds], [Use ORCLUS @aggarwal2000orclus or SSC @elhamifar2013ssc],
    [Overlapping Multi-Label Modules], [Poor: hard disjoint partition assigns one cluster per point], [Use CLIQUE @agrawal1998clique or DOC @procopiuc2002doc],
    [Low-D Non-Convex Shapes], [Poor: $L_1$-medoid balls split rings and moons], [Use STING @wang1997sting or WaveCluster @sheikholeslami1998wavecluster],
    [Noisy Unlabeled Outlier Streams], [Explicit sphere-of-influence outlier filtering $Delta_i$], [Small clusters flagged as bad medoids; use FasterPAM @schubert2019fasterpam],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for PROCLUS.],
)

#sec("V", "Structural Relations and Empirical Behavior")

PROCLUS bridges classical prototype partitioning ($k$-medoids / $k$-means) and modern feature-weighted and subspace clustering @kriegel2009subspacesurvey @jain2010data.

#proposition(title: "Algebraic Unification of PAM, ORCLUS, and Weighted Partitioning", origin: "Original Proof")[
  Let $J(cal(M), {cal(D)_i})$ be PROCLUS's segmental dispersion objective over medoids $cal(M)$ and indicator vectors $u_i in {0, 1}^d$ ($cal(D)_i = "supp"(u_i)$).
  + *Full-Space Reduction (PAM / $L_1$ $k$-Medoids):* Setting $l = d$ forces $u_(i,j) = 1$ for all $(i, j)$, reducing $d_(cal(D)_i)$ identically to $d^(-1) norm(bold(x) - bold(m)_i)_1$ and PROCLUS to CLARANS-style $L_1$ $k$-medoids @kaufman1990finding @ng2002clarans.
  + *Arbitrary-Orientation Generalization (ORCLUS @aggarwal2000orclus):* Replacing diagonal coordinate projections $P_i = "diag"(u_i)$ with arbitrary rank-$l_i$ orthogonal projections $P_i = E_i E_i^top$ ($E_i in RR^(d times l_i)$, $E_i^top E_i = I_(l_i)$) and $L_1$ deviation with $L_2$ variance recovers ORCLUS, where selecting the $l_i$ smallest coordinates generalizes to selecting the $l_i$ smallest-eigenvalue eigenvectors of cluster $i$'s covariance matrix.
  + *Continuous-Weight Relaxation (COSA / W-$k$-Means):* Relaxing binary weights $u_i / |cal(D)_i| in {0, 1 / |cal(D)_i|}^d$ to continuous simplex weights $w_i in Delta^(d-1)$ with entropy or $beta$-power regularization yields COSA @friedman2004cosa and automated variable-weighted $k$-means @huang2005automated.
]

#proof[
  For (1), when $l = d$, the global budget $sum_(i=1)^k sum_(j=1)^d u_(i,j) = k d$ with $u_(i,j) <= 1$ admits the unique feasible solution $u_(i,j) equiv 1$. Then $d_(cal(D)_i)(bold(x), bold(m)_i) = d^(-1) norm(bold(x) - bold(m)_i)_1$ and $J$ reduces (up to factor $1 / d$) to the exact $L_1$ $k$-medoids objective (NP-hard in general, akin to $k$-means @aloise2009np).

  For (2), any diagonal projection $P_i = "diag"(u_i)$ minimizes $sum_(x in C_i) norm(P_i (x - mu_i))_2^2 = sum_(j=1)^d u_(i,j) "Var"_i (x_j)$ by picking the $l_i$ axes of minimum variance. Generalizing $P_i$ to any rank-$l_i$ orthogonal projection yields $min_(E_i^top E_i = I_(l_i)) "Tr"(E_i^top Sigma_i E_i)$, whose exact minimizer by the Courant-Fischer / Ky Fan theorem consists of the $l_i$ eigenvectors of cluster covariance $Sigma_i$ with smallest eigenvalues — precisely ORCLUS's subspace update @aggarwal2000orclus.

  For (3), noting that $v_(i,j) := u_(i,j) / |cal(D)_i|$ satisfies $v_(i,j) >= 0$ and $sum_(j=1)^d v_(i,j) = 1$, PROCLUS's segmental distance is a discrete extreme-point restriction of weighted Manhattan distance $sum_(j=1)^d v_(i,j) |x_j - m_(i,j)|$ over the probability simplex $Delta^(d-1)$; smoothing the indicator selection via softmax weights $v_(i,j) prop exp(-X_(i,j) / lambda)$ directly produces the COSA exponential weight update @friedman2004cosa.
]

*Empirical Benchmark Behavior:* On the ClusterOverview 4-dataset benchmark evaluated via Silhouette @rousseeuw1987silhouettes and ARI @hubert1985comparing, PROCLUS achieves *Iris* $"Sil" = 0.344$, $"ARI" = 0.471$; *Customers* $"Sil" = 0.273$; *Multi-Shape* $"Sil" = 0.451$ ($"ARI" = "N/A"$); and *High-D* $"Sil" = 0.294$ ($"ARI" = "N/A"$), with a mean runtime of *$4.477$ ms*. This empirical profile precisely contrasts PROCLUS against the grid methods (Algs. 133–135):
- *High-Dimensional Leadership:* PROCLUS obtains the *highest* Silhouette score on *High-D* ($0.294$) across the entire quartet — surpassing STING ($0.257$), WaveCluster ($0.257$), and CLIQUE ($"N/A"$). By selecting cluster-specific feature subsets in polynomial time (Theorems II.A and II.C), PROCLUS avoids both exponential grid sparsity and fixed-$tau$ selectivity collapse.
- *Convex Inductive Bias on Low-D Datasets:* On *Iris* ($"Sil" = 0.344, "ARI" = 0.471$) and *Multi-Shape* ($"Sil" = 0.451$), PROCLUS trails grid-density methods (STING/WaveCluster $0.593$ on Multi-Shape; CLIQUE $0.611$ on Iris) because $L_1$-medoid spheres split non-convex geometries and fixing a global average dimension $l$ can under- or over-select coordinates on small 4D datasets.
- *Iterative Computational Cost:* Its $4.477$ ms mean runtime is roughly $31 times$ CLIQUE ($0.145$ ms) and $86 times$ STING ($0.052$ ms), directly reflecting the $cal(O)(I n k d)$ iterative locality and medoid-swap loops versus single-pass low-dimensional grid accumulation.
