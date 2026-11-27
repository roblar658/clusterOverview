#import "../template/ieee.typ": *

#alg_entry(
  id: 14,
  name: "Hamerly's K-Means",
  category: "Cat. I: Centroid & Partitioning",
  year: "2010",
  paradigm: "Triangle-Inequality Pruned Lloyd Acceleration",
  complexity: $O(n k d T)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Let $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$ be a finite dataset equipped with the Euclidean metric $d(x, y) = norm(x - y)_2$, and let $k in {1, dots, n}$ be the target cluster count. Hamerly's algorithm @hamerly2010making solves the exact Euclidean $k$-means within-cluster sum of squares (WCSS) problem @lloyd1982least:
$ J(cal(C), bold(mu)) = sum_(j=1)^k sum_(i in C_j) norm(x_i - mu_j)_2^2, $
producing the *exact* same Voronoi assignments and centroid trajectories as Lloyd's algorithm (Alg. 1) while avoiding redundant distance evaluations via metric triangle inequalities. Whereas Elkan's algorithm (Alg. 13) @elkan2003using maintains a full $n times k$ matrix of per-cluster lower bounds requiring $Theta(n k)$ auxiliary storage, Hamerly @hamerly2010making compresses all non-assigned lower bounds for each datum $x_i$ into a *single scalar* $l_i$, reducing per-point bound storage from $O(k)$ to $O(1)$ and total auxiliary state from $O(n k)$ to $O(n + k^2)$.

For each point $i in {1, dots, n}$, the algorithm tracks an assignment index $a_i in {1, dots, k}$, an upper bound $u_i >= 0$ on the distance to its assigned centroid $mu_(a_i)$, and a single uniform lower bound $l_i$ on the distance to its *closest competing* centroid. Across centroids $j in {1, dots, k}$, the algorithm also maintains the half-separation radius
$ s_j (bold(mu)) = 1 / 2 min_(ell != j) d(mu_j, mu_ell), $
together with the per-iteration centroid drift magnitudes $p_j^((t)) = d(mu_j^((t+1)), mu_j^((t)))$ and the indices of the two largest drifts:
$ r = arg max_(1 <= j <= k) p_j^((t)), quad r' = arg max_(j != r) p_j^((t)). $

#definition(title: "Hamerly Bound-Invariant State")[
  A tuple $(bold(a), bold(u), bold(l), bold(mu), bold(s))$ is *Hamerly-admissible* for dataset $cal(X)$ if for every datum $i in {1, dots, n}$ and every cluster $j in {1, dots, k}$:
  1. *(Upper-bound soundness)* $u_i >= d(x_i, mu_(a_i))$;
  2. *(Aggregated lower-bound soundness)* $l_i <= min_(ell != a_i) d(x_i, mu_ell)$;
  3. *(Inter-centroid half-gap consistency)* $s_j = 1 / 2 min_(ell != j) d(mu_j, mu_ell)$.
]

Crucially, Hamerly's algorithm replaces per-cluster bound checks with a fast two-tier scalar filter comparing $u_i$ against $m_i = max(s_(a_i), l_i)$ before ever computing a Euclidean distance in $bb(R)^d$, and updates cluster vector sums $c'_j = sum_(i in C_j) x_i$ and cardinalities $q_j = |C_j|$ incrementally whenever a point switches clusters (Algorithm 1).

#algbox(
  title: "Hamerly's Exact Pruned K-Means Algorithm",
  inputs: [Dataset $cal(X) subset.eq bb(R)^d$, integer $k >= 2$, initial centroids $bold(mu)^((0))$, max iterations $T$.],
  outputs: [Partition $cal(C) = {C_1, dots, C_k}$ and centroids $bold(mu)$.],
)[
  1. *Initialize:* Set $q_j <- 0$, $c'_j <- bold(0)$ for $j = 1, dots, k$. For each $i = 1, dots, n$, compute all $d(x_i, mu_j^((0)))$, set $a_i <- arg min_j d(x_i, mu_j^((0)))$ (smallest index on ties), $u_i <- d(x_i, mu_(a_i)^((0)))$, $l_i <- min_(j != a_i) d(x_i, mu_j^((0)))$, $q_(a_i) <- q_(a_i) + 1$, $c'_(a_i) <- c'_(a_i) + x_i$.
  2. *For* $t = 0, 1, dots, T - 1$ *do*:
  3. #h(1em) Compute $s_j <- 1/2 min_(ell != j) d(mu_j^((t)), mu_ell^((t)))$ for $j = 1, dots, k$.
  4. #h(1em) *For* $i = 1, dots, n$ *do*:
  5. #h(2em) Set $m_i <- max(s_(a_i), l_i)$.
  6. #h(2em) *If* $u_i >= m_i$ *then* *(Gate 1 fails)*:
  7. #h(3em) Tighten $u_i <- d(x_i, mu_(a_i)^((t)))$.
  8. #h(3em) *If* $u_i >= m_i$ *then* *(Gate 2 fails)*:
  9. #h(4em) Compute $d(x_i, mu_j^((t)))$ for all $j = 1, dots, k$; let $a'_i$ be the smallest-index minimizer, $u_i <- d(x_i, mu_(a'_i)^((t)))$, $l_i <- min_(j != a'_i) d(x_i, mu_j^((t)))$.
  10. #h(4em) *If* $a'_i != a_i$ *then*: update $q_(a_i) <- q_(a_i) - 1$, $c'_(a_i) <- c'_(a_i) - x_i$, $q_(a'_i) <- q_(a'_i) + 1$, $c'_(a'_i) <- c'_(a'_i) + x_i$, and set $a_i <- a'_i$.
  11. #h(1em) *Update centroids:* For each $j$ with $q_j > 0$, set $mu_j^((t+1)) <- c'_j / q_j$ and $p_j <- d(mu_j^((t+1)), mu_j^((t)))$.
  12. #h(1em) *If* $max_j p_j = 0$, *break*.
  13. #h(1em) Find $r = arg max_j p_j$ and $r' = arg max_(j != r) p_j$.
  14. #h(1em) *Update bounds:* For each $i = 1, dots, n$, set $u_i <- u_i + p_(a_i)$ and
      $ l_i <- cases(l_i - p_(r') & "if" a_i = r, l_i - p_r & "if" a_i != r). $
  15. *Return* $C_j = {i : a_i = j}$ and $bold(mu)^((t+1))$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish soundness of both scalar gates, inductive equivalence to Lloyd's trajectory, minimax sharpness of the dual-drift update, a drift-spillover counterexample vs. Elkan's algorithm, and high-dimensional pruning collapse under distance concentration.

#subsec("A", "Inductive Bound Invariance and Exact Lloyd Equivalence")

#lemma(title: "Soundness of Hamerly's Two-Gate Scalar Filter", origin: "Original Proof")[
  Let $(bold(a), bold(u), bold(l), bold(mu), bold(s))$ be Hamerly-admissible. For any $i in {1, dots, n}$:
  1. If $u_i <= max(s_(a_i), l_i)$, then $d(x_i, mu_(a_i)) <= d(x_i, mu_j)$ for every $j != a_i$, so $x_i in V_(a_i) (bold(mu))$.
  2. If the inequality is strict, $u_i < max(s_(a_i), l_i)$, then $d(x_i, mu_(a_i)) < d(x_i, mu_j)$ for all $j != a_i$; hence $a_i$ is the *unique* nearest centroid of $x_i$ regardless of tie-breaking.
]

#proof[
  Set $m_i = max(s_(a_i), l_i)$ and consider the two exhaustive branches realizing $m_i$:

  *Case 1 ($m_i = l_i$):* By Definition 1(1)–(2),
  $ d(x_i, mu_(a_i)) <= u_i <= l_i <= min_(j != a_i) d(x_i, mu_j). $
  Thus $d(x_i, mu_(a_i)) <= d(x_i, mu_j)$ for every $j != a_i$, and if $u_i < l_i$ the second inequality is strict.

  *Case 2 ($m_i = s_(a_i)$):* Take any $j != a_i$. By Definition 1(3) and the Euclidean triangle inequality,
  $ 2 s_(a_i) <= d(mu_(a_i), mu_j) <= d(mu_(a_i), x_i) + d(x_i, mu_j). $
  Since $d(mu_(a_i), x_i) = d(x_i, mu_(a_i)) <= u_i <= s_(a_i)$, rearranging gives
  $ d(x_i, mu_j) &>= 2 s_(a_i) - d(x_i, mu_(a_i)) \
                 &>= 2 s_(a_i) - u_i >= u_i >= d(x_i, mu_(a_i)). $
  Whenever $u_i < s_(a_i)$, we have $2 s_(a_i) - u_i > u_i >= d(x_i, mu_(a_i))$, establishing strict inequality $d(x_i, mu_(a_i)) < d(x_i, mu_j)$ for all $j != a_i$. Because tightening $u_i <- d(x_i, mu_(a_i))$ preserves Definition 1(1) with equality, both Gate 1 and Gate 2 are sound.
]

#theorem(title: "Inductive Bound Preservation and Exact Lloyd Equivalence", origin: "Original Proof")[
  For any finite dataset $cal(X) subset.eq bb(R)^d$ and initial centroids $bold(mu)^((0))$ with non-empty clusters, Algorithm 1 satisfies:
  1. At the start of every iteration $t >= 0$, the state $(bold(a)^((t)), bold(u)^((t)), bold(l)^((t)), bold(mu)^((t)), bold(s)^((t)))$ is Hamerly-admissible.
  2. At every iteration $t >= 0$, Algorithm 1 produces the *identical* partition $cal(C)^((t))$ and centroid tuple $bold(mu)^((t+1))$ as Lloyd's algorithm (Alg. 1) with smallest-index tie-breaking, achieves strict monotonic descent of $J(cal(C)^((t)), bold(mu)^((t)))$, and terminates in at most $S(n, k) <= k^n / k!$ iterations at a discrete centroidal Voronoi equilibrium.
]

#proof[
  We proceed by strong induction on $t >= 0$.

  *Base case ($t = 0$):* Step 1 evaluates exact distances $d(x_i, mu_j^((0)))$ for all $i, j$, setting $u_i^((0)) = d(x_i, mu_(a_i)^((0)))$ and $l_i^((0)) = min_(j != a_i) d(x_i, mu_j^((0)))$, and Step 3 computes exact half-gaps $s_j^((0))$. Thus Hamerly-admissibility holds with exact equality at $t = 0$.

  *Inductive step:* Assume admissibility at the start of iteration $t$, with $bold(mu)^((t))$ matching Lloyd's $t$-th iterate. For each $i = 1, dots, n$: if Gate 1 or Gate 2 prunes ($u_i < m_i$), Lemma 1(2) ensures $a_i$ is the unique nearest centroid (matching Lloyd while preserving valid bounds); otherwise Step 9 evaluates all $k$ exact distances, picks the smallest-index minimizer $a'_i$, and resets $u_i, l_i$ with equality. Incremental maintenance of $c'_j$ and $q_j$ yields the exact Lloyd centroid $mu_j^((t+1)) = c'_j / q_j$.

  For Step 14, let $tilde(u)_i, tilde(l)_i$ denote the pre-update valid bounds. By the triangle inequality,
  $ d(x_i, mu_(a_i)^((t+1))) <= d(x_i, mu_(a_i)^((t))) + p_(a_i) <= tilde(u)_i + p_(a_i) = u_i^((t+1)). $
  For any competing cluster $j != a_i$, the reverse triangle inequality gives $d(x_i, mu_j^((t+1))) >= tilde(l)_i - p_j$. Since $max_(j != a_i) p_j$ equals $p_r$ when $a_i != r$ and $p_(r')$ when $a_i = r$, Step 14 yields $l_i^((t+1)) <= min_(j != a_i) d(x_i, mu_j^((t+1)))$ (even if negative, since $d >= 0$). Refreshing $s_j^((t+1))$ completes the induction. Monotonic descent of $J$ and termination in $<= S(n, k) <= k^n / k!$ iterations follow directly from Theorem 1 of Alg. 1.
]

#remark(title: "Tie-Breaking Convention")[
  Hamerly's original pseudocode @hamerly2010making tests $u_i > m_i$, which preserves Voronoi optimality and Lyapunov descent on measure-zero ties; using $u_i >= m_i$ guarantees exact bitwise agreement with smallest-index Lloyd tie-breaking.
]

#subsec("B", "Minimax Sharpness of Dual-Drift Lower-Bound Shifts")

#proposition(title: "Optimality of the Dual-Drift Subtraction Rule", origin: "Original Derivation")[
  Let $Delta(a, bold(p))$ be any scalar shift function updating $l_i^((t+1)) = tilde(l)_i - Delta(a_i, bold(p))$ using solely the assigned cluster $a_i$ and the vector of centroid drift magnitudes $bold(p) = (p_1, dots, p_k) in bb(R)_(>= 0)^k$ (without storing the index of the second-nearest centroid). Then Hamerly's shift
  $ Delta^* (a_i, bold(p)) = cases(p_(r') & "if" a_i = r, p_r & "if" a_i != r) $
  is the *uniformly smallest* shift that preserves Definition 1(2) over all Euclidean configurations.
]

#proof[
  Admissibility of $Delta^*$ was established in Theorem 1. To prove no strictly smaller shift can be universally valid, fix any cluster $a in {1, dots, k}$ and drift vector $bold(p) in bb(R)_(>= 0)^k$, and let $j^* = arg max_(j != a) p_j$, so $p_(j^*) = Delta^* (a, bold(p))$. Choose $L > p_(j^*) + 1$ and construct a collinear configuration along unit vector $e_1 in bb(R)^d$: set $x_i = bold(0)$, $mu_a^((t)) = bold(0)$, $mu_(j^*)^((t)) = L e_1$, and $mu_j^((t)) = 2 L e_1$ for all $j in.not {a, j^*}$. Then $a_i = a$ and the exact pre-update lower bound is $tilde(l)_i = d(x_i, mu_(j^*)^((t))) = L$. Now displace centroid $j^*$ directly toward $x_i$ to $mu_(j^*)^((t+1)) = (L - p_(j^*)) e_1$ while displacing all other centroids orthogonally to $e_1$ by $p_j$. Then
  $ min_(j != a) d(x_i, mu_j^((t+1))) &= d(x_i, mu_(j^*)^((t+1))) \
                                      &= L - p_(j^*) = tilde(l)_i - Delta^* (a, bold(p)). $
  Hence any rule with $Delta(a, bold(p)) < Delta^* (a, bold(p))$ outputs $l_i^((t+1)) > min_(j != a) d(x_i, mu_j^((t+1)))$, violating admissibility.
]

#subsec("C", "Global Drift Spillover Counterexample vs. Elkan's Algorithm")

Compressing $k - 1$ per-centroid lower bounds into one scalar $l_i$ introduces a fundamental structural vulnerability: a single rapidly moving centroid $r$ degrades $l_i$ for *all* points in every other cluster $j != r$, regardless of geometric separation.

#theorem(title: "Single-Cluster Drift Spillover Breakdown", origin: "Original Counterexample")[
  For every integer $k >= 3$ and every $m >= 1$, there exists a one-dimensional dataset $cal(X) subset.eq bb(R)$ (embeddable in $bb(R)^d$ for any $d >= 1$) of size $n = 4 m + k - 1$ and initialization $bold(mu)^((0))$ such that at iteration $t = 1$:
  1. Zero points change cluster membership ($cal(C)^((1)) = cal(C)^((0))$), and $k - 1$ clusters are completely stationary ($p_1 = dots = p_(k-1) = 0$).
  2. Hamerly's algorithm fails *both* Gate 1 and Gate 2 on $m = Theta(n)$ stationary points and executes at least $m (k + 1) = Theta(n k)$ Euclidean distance calculations.
  3. Elkan's algorithm (Alg. 13) prunes all $m$ points using its per-cluster lower bounds and executes *zero* point-to-centroid distance evaluations on those $m$ points.
]

#proof[
  Fix $k >= 3$ and $m >= 1$, and construct $cal(X) subset.eq bb(R)$ ($n = 4 m + k - 1$): cluster $C_1$ has $m$ copies of $-3/2$ and $3 m$ copies of $1/2$; $C_2 = {2}$; $C_j = {10(j - 1)}$ for $3 <= j <= k - 1$; and $C_k = {10 k + 97.5}$. Initialize $mu_1^((0)) = 0$, $mu_2^((0)) = 2$, $mu_j^((0)) = 10(j - 1)$ ($3 <= j <= k - 1$), and $mu_k^((0)) = 10 k + 100$.

  At $t = 0$, every copy of $-3/2$ has $a_i = 1$, $u_i^((0)) = 3/2$, and $l_i^((0)) = d(-3/2, 2) = 7/2$; every copy of $1/2$ has $a_i = 1$ ($u_i = 1/2 < 3/2$); singletons in $C_2, dots, C_(k-1)$ have $u_i = 0$; and $C_k$'s point lies distance $5/2$ from $mu_k^((0))$ vs. $117.5$ from $mu_(k-1)^((0))$. Updating centroids yields
  $ mu_1^((1)) = (m (-3/2) + 3 m (1/2)) / (4 m) = 0 = mu_1^((0)), $
  so $p_1 = dots = p_(k-1) = 0$, while $mu_k^((1)) = 10 k + 97.5$ gives $r = k$, $p_r = 5/2$, and $p_(r') = 0$.

  In Step 14, all $m$ points at $x_i = -3/2$ ($a_i = 1 != r$) receive $u_i^((1)) = 3/2 + 0 = 3/2$ and $l_i^((1)) = 7/2 - 5/2 = 1$. Since $s_1^((1)) = 1/2 d(0, 2) = 1$, iteration $t = 1$ has $m_i = max(1, 1) = 1 < 3/2 = u_i^((1))$. Both Gate 1 and tightened Gate 2 ($d(-3/2, 0) = 3/2 >= 1$) fail, forcing Step 9 to run $k$ distance evaluations on all $m$ points—totaling $m (k + 1) = Theta(n k)$ distances with zero assignment changes. Conversely, Elkan's per-cluster bounds satisfy $l_(i, 2)^((1)) = 7/2$, $l_(i, j)^((1)) >= 43/2$ ($3 <= j <= k - 1$), and $l_(i, k)^((1)) = 10 k + 99$, all strictly exceeding $u_i^((1)) = 3/2$, so Elkan executes zero point-to-centroid distances on all $m$ points.
]

#subsec("D", "High-Dimensional Pruning Collapse and Classical Complexity Bounds")

Even without adversarial drift spillover, Hamerly's two pruning gates undergo a sharp geometric phase transition as ambient dimension $d$ grows, driven by Euclidean distance concentration @beyer1999nearest @hamerly2010making.

#theorem(title: "Asymptotic Failure of Half-Gap Pruning in High Dimensions", origin: "Original Derivation")[
  Suppose data points $x_i$ and centroids $mu_j$ are drawn independently from an isotropic product distribution in $bb(R)^d$ with coordinate mean $0$, variance $sigma^2 > 0$, and finite fourth moment $bb(E)[X_r^4] < infinity$. Then as $d -> infinity$ for fixed $k >= 2$:
  1. The ratio of Hamerly's half-gap bound to the exact assigned distance converges in probability to $1 / 2$:
     $ s_(a_i) / d(x_i, mu_(a_i)) attach(->, t: "p") 1 / 2 < 1. $
  2. The relative margin $(d_2 (x_i) - d_1 (x_i)) / d_1 (x_i)$ between second-nearest and nearest centroid distances satisfies $O_p (d^(-1/2)) -> 0$. Consequently, the probability that Gate 1 prunes via $s_(a_i)$ vanishes as $d -> infinity$, and any centroid drift $p_r = omega(sigma)$ causes $l_i - p_r < u_i$ with probability tending to $1$.
]

#proof[
  For any two independent draws $u, v in bb(R)^d$, write $d(u, v)^2 = sum_(r=1)^d (u_r - v_r)^2$. Since coordinates are i.i.d. with $bb(E)[(u_r - v_r)^2] = 2 sigma^2$ and finite variance $v_0 := "Var"((u_r - v_r)^2) < infinity$, Chebyshev's inequality yields
  $ d(u, v)^2 / (2 d sigma^2) = 1 + O_p (d^(-1/2)). $
  Taking square roots via the delta method gives $d(u, v) = sqrt(2 d) sigma (1 + O_p (d^(-1/2)))$ for each of the $k$ point-to-centroid distances and all $binom(k, 2)$ inter-centroid distances. Taking the minimum over finitely many ($k - 1$) competing pairs preserves the asymptotic scaling:
  $ d_1 (x_i) &= sqrt(2 d) sigma (1 + O_p (d^(-1/2))), \
    d_2 (x_i) &= sqrt(2 d) sigma (1 + O_p (d^(-1/2))), \
    s_(a_i)   &= 1 / 2 min_(ell != a_i) d(mu_(a_i), mu_ell) = sqrt(d / 2) sigma (1 + O_p (d^(-1/2))). $
  Dividing $s_(a_i)$ by $d_1 (x_i) = d(x_i, mu_(a_i))$ cancels $sqrt(2 d) sigma$ and yields $s_(a_i) / d(x_i, mu_(a_i)) = 1/2 + O_p (d^(-1/2)) attach(->, t: "p") 1/2$. Because $u_i >= d(x_i, mu_(a_i))$, the half-gap test $u_i <= s_(a_i)$ holds with probability $-> 0$. Moreover, $l_i - u_i <= d_2 (x_i) - d_1 (x_i) = O_p (sigma)$, so subtracting any drift $p_r$ with $p_r / sigma -> infinity$ drives $l_i - p_r < u_i$ with probability approaching $1$, collapsing Hamerly's work back to $Theta(n k)$ distances per sweep.
]

#corollary(title: "Seeding Guarantees and Worst-Case vs. Smoothed Iterations", origin: "Classical — Cited")[
  Because Theorem 1 guarantees exact trajectory equivalence to Lloyd's algorithm:
  1. Initializing Algorithm 1 with K-Means++ seeding @arthur2007kmeanspp yields expected objective $bb(E)[J(cal(C))] <= 8 (ln k + 2) J^*$ while retaining exact bound-pruned iterations.
  2. Worst-case iteration count remains $2^(Omega(sqrt(n)))$ in the plane @vattani2011kmeans (with NP-hard global optimization @aloise2009np), while Gaussian-perturbed inputs converge in smoothed expected $"poly"(n, k, d, 1/sigma)$ iterations @arthur2011smoothed.
]

#proof(label: "Proof Synthesis (Cited)")[
  By Theorem 1, every iterate $(cal(C)^((t)), bold(mu)^((t)))$ of Hamerly's algorithm coincides with Lloyd's algorithm. Hence the $O(log k)$ expected approximation guarantee of Arthur and Vassilvitskii @arthur2007kmeanspp (which holds immediately after seeding and can only improve under monotonic Lloyd descent), Vattani's $2^(Omega(sqrt(n)))$ planar lower bound @vattani2011kmeans, and the smoothed polynomial bound of Arthur, Manthey, and Röglin @arthur2011smoothed transfer verbatim to Algorithm 1.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Each iteration of Hamerly's algorithm performs two distinct classes of work: (i) $O(k^2 d)$ inter-centroid distance computations to refresh $s_j$ (negligible when $k << n$), and (ii) a point sweep costing $O(n)$ scalar bound checks plus $O(n_"g1" d)$ tightening distances for points failing Gate 1 and $O(n_"g2" k d)$ full inner-loop evaluations for points failing Gate 2. In low-to-moderate dimensions ($d <= 50$) with well-separated clusters, $n_"g2" << n$ after the first few iterations, reducing effective per-iteration complexity from $Theta(n k d)$ toward $O(n + k^2 d)$—often achieving $5 times$ to $50 times$ speedups over unpruned Lloyd @hamerly2010making.

Critically, auxiliary memory is only $O(n + k^2)$ floats (two arrays `u[n]`, `l[n]`, assignment array `a[n]`, and $k$ half-gaps `s[k]`), versus $O(n k)$ for Elkan's algorithm (Alg. 13). When $n = 10^6$ and $k = 10^3$, Elkan requires $8$ GB for double-precision lower bounds alone, causing severe DRAM cache thrashing, whereas Hamerly's $u, l$ arrays fit in $16$ MB of L3 cache. Because every distance evaluation uses exact Euclidean norms, all geometric properties—convex Voronoi polyhedra, $O(d)$ orthogonal invariance, scale sensitivity, and $1/n$ outlier breakdown—are identical to Lloyd's algorithm.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for Hamerly's K-Means],
    [Time Complexity], [Worst-case $O((n k + k^2) d T)$; practical low-$d$ sweep approaches $O(n + k^2 d + n_"g2" k d)$.],
    [Auxiliary Space], [$O(n + k d)$ (just 2 bounds $u_i, l_i$ per point), vs. $O(n k + k d)$ for Elkan (Alg. 13).],
    [Exactness Guarantee], [Bitwise identical trajectory to Lloyd's K-Means (Alg. 1) under $u_i >= m_i$ (Theorem 1).],
    [Optimal Operating Regime], [Low-to-medium $d$ ($d <= 50$), large $n$, moderate-to-large $k$ ($k$ up to $10^3$–$10^4$).],
    [Breakdown / Pruning Risk], [Outlier breakdown $epsilon^* = 1/n$; drift spillover (Thm. 2) and high-$d$ collapse (Thm. 3).],
    [Companion Algorithms], [K-Means++ (Alg. 2) for seeding; Elkan (Alg. 13) for higher $d$ / small $k$; Mini-Batch (Alg. 3) for web-scale.],
  ),
  caption: [Operational and theoretical profile of Hamerly's exact accelerated K-Means.],
)

#sec("IV", "Engineering and Scientific Applications")

Hamerly's algorithm excels wherever *exact* Lloyd convergence is required on large datasets in low-to-moderate dimensions, particularly when memory bandwidth or $n k$ cache footprint makes Elkan's $O(n k)$ matrix prohibitive.

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Image Color Quantization & Superpixel Codebooks:* RGB/CIELAB pixels occupy $d = 3$ with $n = 10^6$–$10^8$ and $k = 16$–$256$. Distance concentration (Theorem 3) is absent, $s_j$ prunes $> 90\%$ of pixels after 3–5 sweeps, and $O(n)$ bounds fit in cache.
2. *Spatial GIS, LiDAR, and Sensor-Grid Partitioning:* In $d = 2$ or $3$, boundary points form an $O(n^((d-1)/d))$ surface fraction of each Voronoi cell, allowing Gate 1 to eliminate nearly all interior points in $O(1)$ scalar tests.
3. *Product-Quantization (PQ) Subspace Codebooks:* PQ decomposes dense embeddings into $M$ subspaces of dimension $d / M in [4, 16]$ with $k = 256$ centroids each—matching Hamerly's exact low-$d$, moderate-$k$ cache sweet spot.
4. *Memory-Constrained Embedded & Edge Analytics:* Where RAM scales linearly with $n$, Elkan's $O(n k)$ lower-bound matrix triggers OOM failures at large $k$, whereas Hamerly sustains exact Lloyd convergence within a $2 n$ scalar budget.
5. *Multi-Start Model Selection & $k$-Sweeps:* Evaluating Elbow/Silhouette curves across $k$ and repeated K-Means++ restarts @arthur2007kmeanspp multiplies work by $10$–$100 times$; Hamerly accelerates every restart with zero WCSS approximation bias.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Hamerly's algorithm degrades in four mathematically predictable regimes: (a) *high intrinsic dimensionality* ($d >> 50$), where Theorem 3 proves $s_(a_i) / u_i -> 1/2$ and relative margins collapse, causing pruning overhead without speedup (favoring BLAS-3 GEMM Lloyd or Elkan for small $k$); (b) *single-cluster drift instability* (Theorem 2), where one noisy outlier cluster repeatedly shifts and invalidates $l_i$ globally (resolved by Elkan Alg. 13 or trimmed/robust updates); (c) *massive $k approx n$*, where the $O(k^2 d)$ half-gap computation dominates; and (d) *non-convex or heavy-tailed cluster geometry*, inherited from the WCSS objective itself (resolved by Kernel K-Means Alg. 16, K-Medians Alg. 15, or DBSCAN Alg. 31).

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Why Hamerly Fits / Key Mechanism], [Caution / Alternative],
    [RGB / CIELAB Quant.], [$d = 3$ yields $>90\%$ pruning; $O(n)$ RAM], [Use Mini-Batch (Alg. 3) if $n > 10^8$],
    [PQ Subspace Codebooks], [Ideal $d in [4,16], k = 256$ cache fit], [Use GEMM Lloyd on GPU],
    [2D/3D GIS & LiDAR], [Interior Voronoi points pruned at Gate 1], [Use DBSCAN (Alg. 31) for arcs],
    [High-D Text / Dense Emb.], [Fails Theorem 3 ($s / u -> 1/2$)], [Use Spherical (Alg. 12) or GEMM],
    [Outlier-Heavy Streams], [Drift spillover (Thm. 2) + $1/n$ breakdown], [Use K-Medians (Alg. 15) / PAM],
  ),
  caption: [Domain suitability matrix and failure-regime guidance for Hamerly's K-Means.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Within Category I, triangle-inequality accelerations form a strict algebraic hierarchy parameterized by the number of stored lower bounds per point $b in {0, 1, dots, k - 1}$.

#proposition(title: "Monotone Memory–Pruning Hierarchy across Lloyd, Hamerly, and Elkan", origin: "Original Unified Proof")[
  For each point $x_i$, let $cal(A)_b$ maintain $u_i >= d(x_i, mu_(a_i))$, the $b - 1$ tightest individual lower bounds $l_(i, j_1), dots, l_(i, j_(b-1))$ for distinct competing clusters, and one residual lower bound $l_(i, "res") <= min_(j in.not {a_i, j_1, dots, j_(b-1)}) d(x_i, mu_j)$. Then:
  1. $b = 0$ recovers unpruned Lloyd (Alg. 1), $b = 1$ recovers Hamerly (Alg. 14), and $b = k - 1$ recovers Elkan (Alg. 13).
  2. At any fixed state $(cal(C)^((t)), bold(mu)^((t)))$ and drift vector $bold(p)$, the post-update effective lower bound $L_i^((b)) (j)$ on every competing cluster $j != a_i$ is monotonically non-decreasing in $b$:
     $ 0 &= L_i^((0)) (j) <= L_i^((1)) (j) <= dots \
         &<= L_i^((k-1)) (j) <= d(x_i, mu_j^((t+1))). $
  3. Auxiliary bound memory grows linearly as $Theta(n (b + 1))$, placing Hamerly ($b = 1$) at the unique *minimum-memory non-trivial vertex* of the exact pruning Pareto frontier.
]

#proof[
  For $b = 0$, no lower bounds are stored ($L_i^((0)) (j) equiv 0$), forcing full $k$-distance evaluation (Lloyd). For $b = 1$, zero individual bounds are stored and $L_i^((1)) (j) = l_(i, "res")^((t+1)) = tilde(l)_i - Delta^* (a_i, bold(p))$ for all $j != a_i$ (Hamerly). For $b = k - 1$, every competing cluster $j != a_i$ has its own bound $L_i^((k-1)) (j) = tilde(l)_(i, j) - p_j$ (Elkan). Passing from $b$ to $b + 1$ promotes one competing cluster $j_b$ out of the residual pool $R_b$ into an individual slot. For $j = j_b$, since $tilde(l)_(i, j_b) >= min_(ell in R_b) tilde(l)_(i, ell) = tilde(l)_(i, "res")$ and $p_(j_b) <= max_(ell in R_b) p_ell$, we have $tilde(l)_(i, j_b) - p_(j_b) >= tilde(l)_(i, "res") - max_(ell in R_b) p_ell$. For every remaining $j in R_(b+1) = R_b backslash {j_b}$, taking the minimum over the strictly smaller set $R_(b+1)$ can only increase $tilde(l)_(i, "res")$ while decreasing $max_(ell in R_(b+1)) p_ell$, so $L_i^((b+1)) (j) >= L_i^((b)) (j)$ uniformly for all $j != a_i$. Memory is $1$ upper bound plus $b$ lower bounds per datum, i.e., $Theta(n (b + 1))$.
]

Beyond Elkan and Lloyd, Hamerly's exact acceleration composes orthogonally with initialization schemes such as K-Means++ (Alg. 2) and parallel K-Means|| @bahmani2012scalable, and stands in contrast to approximate sub-linear accelerators like Mini-Batch K-Means (Alg. 3) @sculley2010webscale and coreset reductions @harpeled2004coresets, which trade exactness for $o(n)$ per-step cost.

*Empirical Behavior on the ClusterOverview Benchmark Suite:* Across the 4-dataset empirical benchmark evaluated by Silhouette @rousseeuw1987silhouettes and Adjusted Rand Index (ARI) @hubert1985comparing, Hamerly's K-Means achieves the fastest mean runtime of all three partitioning variants at *0.049 ms*—roughly $3 times$ faster than K-Medians (0.147 ms) and over $21 times$ faster than Kernel K-Means (1.052 ms)—confirming its cache-resident $O(n)$ pruning efficiency in low-to-moderate dimensions. On cluster quality, it records *Iris Sil = 0.450, ARI = 0.563*, *Customers Sil = 0.341*, and *Multi-Shape Sil = 0.609* (ARI N/A), faithfully matching exact Euclidean Lloyd Voronoi partitions. However, on *High-D Subspace+Noise* its Silhouette falls to *0.158* (ARI N/A), directly corroborating both the isotropic $L_2$ sensitivity to uninformative noise dimensions and the high-dimensional distance-concentration regime formalized in Theorem 3.
