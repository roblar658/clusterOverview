#import "../template/ieee.typ": *

#alg_entry(
  id: 26,
  name: "Trimmed K-Means",
  category: "Cat. I: Centroid & Partitioning",
  year: "1997",
  paradigm: "Impartial Trimming / Concentration-Step Minimization",
  complexity: $O(n (k d + log n) T)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Classical $k$-means clustering (Alg. 1) @lloyd1982least minimizes squared Euclidean distortion across every sample point, rendering its centroids vulnerable to even a single gross outlier (breakdown point $1/n$). *Trimmed $k$-means*, introduced by Cuesta-Albertos, Gordaliza, and Matrán @cuesta1997trimmed, robustifies vector quantization through *impartial trimming*: rather than prespecifying which region of $bb(R)^d$ is aberrant, the procedure jointly selects the most coherent subset of cardinality $h = n - floor(alpha n) = ceil((1 - alpha) n)$ together with $k$ cluster prototypes.

Let $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$ be a finite dataset, let $k in {1, dots, n}$ be the target cluster count, and let $alpha in [0, 1)$ be the trimming proportion satisfying $h := n - floor(alpha n) >= k$ (with $r := floor(alpha n)$ discarded observations). Denote by $cal(A)_h = {A subset.eq {1, dots, n} : |A| = h}$ the family of admissible active index sets. Given $A in cal(A)_h$, a disjoint $k$-partition $cal(C) = {C_1, dots, C_k}$ of $A$, and centroids $bold(mu) = (mu_1, dots, mu_k) in (bb(R)^d)^k$, the *trimmed within-cluster sum of squares* objective @cuesta1997trimmed is
$ J_alpha (A, cal(C), bold(mu)) = 1 / h sum_(j=1)^k sum_(i in C_j) norm(x_i - mu_j)_2^2. $
Minimizing over $(A, cal(C))$ for fixed centroids $bold(mu)$ yields the profile objective
$ F_alpha (bold(mu)) = 1 / h min_(A in cal(A)_h) sum_(i in A) d_i (bold(mu)), $
where $d_i (bold(mu)) := min_(1 <= j <= k) norm(x_i - mu_j)_2^2$. Ordering the residuals as $d_((1)) (bold(mu)) <= d_((2)) (bold(mu)) <= dots <= d_((n)) (bold(mu))$, the optimal active set retains the $h$ smallest residuals, so $F_alpha (bold(mu)) = h^(-1) sum_(s=1)^h d_((s)) (bold(mu))$. Geometrically, letting $rho (bold(mu)) := sqrt(d_((h)) (bold(mu)))$ be the trimming radius, the active region is the union of $k$ closed Euclidean balls of equal radius:
$ B(bold(mu), rho) = union.big_(j=1)^k overline(B) (mu_j, rho), $
intersected with the Voronoi cells $V_j (bold(mu)) = {x in bb(R)^d : norm(x - mu_j)_2^2 <= norm(x - mu_ell)_2^2 space forall ell != j}$.

#definition(title: "Impartially Trimmed Centroidal Equilibrium")[
  A triple $(A, cal(C), bold(mu))$ with $A in cal(A)_h$ and $cal(C)$ partitioning $A$ is a *trimmed centroidal Voronoi equilibrium* if: (i) $max_(i in A) d_i (bold(mu)) <= min_(u in.not A) d_u (bold(mu))$ (optimal trimming set); (ii) every $i in C_j$ satisfies $x_i in V_j (bold(mu))$ (Voronoi optimality); and (iii) $mu_j = |C_j|^(-1) sum_(i in C_j) x_i$ for every non-empty $C_j$ (centroid optimality).
]

Because joint optimization over $(A, cal(C), bold(mu))$ is combinatorial and non-convex, Trimmed $k$-Means is solved via alternating *concentration steps* (C-steps), summarized in Algorithm 1. In practice, multiple random initializations (e.g., $k$-means++ seeding @arthur2007kmeanspp on subsamples) are executed and the equilibrium with lowest $J_alpha$ is retained.

#algbox(
  title: "Concentration-Step Trimmed K-Means",
  inputs: [Dataset $cal(X) subset.eq bb(R)^d$, clusters $k >= 1$, trimming rate $alpha in [0, 1)$ with $h = n - floor(alpha n) >= k$, initial centroids $bold(mu)^((0))$, max iterations $T$.],
  outputs: [Active set $A in cal(A)_h$, partition $cal(C)$ of $A$, trimmed set $A^c$, and centroids $bold(mu)$.],
)[
  1. *For* $t = 0, 1, dots, T - 1$ *do*:
  2. #h(1em) *Distance & Voronoi step:* For each $i in {1, dots, n}$, compute
     $ &c_i^((t)) = arg min_(1 <= j <= k) norm(x_i - mu_j^((t)))_2^2, \
       &d_i^((t)) = norm(x_i - mu_(c_i^((t)))^((t)))_2^2, $
     breaking ties by smallest cluster index $j$.
  3. #h(1em) *Impartial trimming (concentration) step:* Select the $h$ indices with smallest $d_i^((t))$ (breaking ties by smallest point index $i$) to form $A^((t)) in cal(A)_h$; set $C_j^((t)) = {i in A^((t)) : c_i^((t)) = j}$.
  4. #h(1em) *Centroid update:* For each non-empty $C_j^((t))$, compute
     $ mu_j^((t+1)) = 1 / (|C_j^((t))|) sum_(i in C_j^((t))) x_i. $
     (If any $C_j^((t)) = emptyset$, reseed $mu_j^((t+1))$ to $x_(i^*)$ maximizing $d_i^((t))$ over $i in A^((t))$.)
  5. #h(1em) *If* $(A^((t)), cal(C)^((t))) = (A^((t-1)), cal(C)^((t-1)))$ *or* $bold(mu)^((t+1)) = bold(mu)^((t))$, *break*.
  6. *Return* $(A^((t)), cal(C)^((t)), {1, dots, n} backslash A^((t)), bold(mu)^((t+1)))$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five fundamental results for Trimmed $k$-Means: (A) order-statistic optimality and strict C-step Lyapunov convergence in finitely many steps; (B) an exact finite-sample Donoho–Huber breakdown characterization distinguishing scattered from concentrated contamination; (C) existence and strong statistical consistency *without any moment assumptions*; and (D) a constructive counterexample proving catastrophic small-cluster erasure under over-trimming.

#subsec("A", "Concentration-Step Monotonicity and Finite Termination")

#lemma(title: "Order-Statistic Optimality of Impartial Trimming", origin: "Original Proof")[
  For any vector of non-negative residuals $(d_1, dots, d_n) in bb(R)_(>= 0)^n$ and integer $1 <= h <= n$, a subset $A^* in cal(A)_h$ minimizes $A |-> sum_(i in A) d_i$ if and only if
  $ max_(v in A^*) d_v <= min_(u in.not A^*) d_u. $
  Moreover, for any $A in cal(A)_h$, letting $m = |A backslash A^*| = |A^* backslash A|$,
  $ &sum_(i in A) d_i - sum_(i in A^*) d_i \
    &quad >= m ( min_(u in A backslash A^*) d_u - max_(v in A^* backslash A) d_v ) >= 0. $
]

#proof[
  Suppose $max_(v in A^*) d_v <= min_(u in.not A^*) d_u$. For any $A in cal(A)_h$, write $A = (A inter A^*) union (A backslash A^*)$ and $A^* = (A inter A^*) union (A^* backslash A)$. Since $|A| = |A^*| = h$, the symmetric differences have equal size $m = |A backslash A^*| = |A^* backslash A|$. Canceling the shared sum over $A inter A^*$ gives
  $ sum_(i in A) d_i - sum_(i in A^*) d_i &= sum_(u in A backslash A^*) d_u - sum_(v in A^* backslash A) d_v \
    &>= m min_(u in A backslash A^*) d_u - m max_(v in A^* backslash A) d_v >= 0, $
  because $A backslash A^* subset.eq (A^*)^c$ and $A^* backslash A subset.eq A^*$. Conversely, if there existed $v in A^*$ and $u in.not A^*$ with $d_u < d_v$, swapping $v$ for $u$ to form $A' = (A^* backslash {v}) union {u} in cal(A)_h$ yields $sum_(i in A') d_i - sum_(i in A^*) d_i = d_u - d_v < 0$, contradicting optimality of $A^*$.
]

#theorem(title: "Strict C-Step Lyapunov Descent and Finite Convergence", origin: "Original Proof")[
  Let ${(A^((t)), cal(C)^((t)), bold(mu)^((t)))}$ be generated by Algorithm 1 with deterministic lexicographic tie-breaking and non-empty active clusters. Then for every $t >= 0$:
  1. $J_alpha (A^((t)), cal(C)^((t)), bold(mu)^((t+1))) = J_alpha (A^((t)), cal(C)^((t)), bold(mu)^((t))) - h^(-1) sum_(j=1)^k |C_j^((t))| norm(mu_j^((t+1)) - mu_j^((t)))_2^2$.
  2. $J_alpha (A^((t+1)), cal(C)^((t+1)), bold(mu)^((t+1))) <= J_alpha (A^((t)), cal(C)^((t)), bold(mu)^((t+1)))$.
  3. Algorithm 1 halts in at most $T <= binom(n, h) sum_(m=1)^k S(h, m)$ iterations at a trimmed centroidal equilibrium, where $S(h, m)$ is the Stirling number of the second kind.
]

#proof[
  Holding $(A^((t)), cal(C)^((t)))$ fixed, apply the Euclidean parallel-axis identity to each cluster $C_j^((t))$ around its arithmetic centroid $mu_j^((t+1)) = |C_j^((t))|^(-1) sum_(i in C_j^((t))) x_i$:
  $ sum_(i in C_j^((t))) norm(x_i - mu_j^((t)))_2^2 &= sum_(i in C_j^((t))) norm(x_i - mu_j^((t+1)))_2^2 \
    &quad + |C_j^((t))| norm(mu_j^((t+1)) - mu_j^((t)))_2^2. $
  Summing over $j = 1, dots, k$ and dividing by $h$ proves Claim 1 with exact equality.

  For Claim 2, consider iteration $t + 1$ with centroids $bold(mu)^((t+1))$ fixed. First, every point $i in A^((t))$ satisfies $d_i^((t+1)) = min_ell norm(x_i - mu_ell^((t+1)))_2^2 <= norm(x_i - mu_(c_i^((t)))^((t+1)))_2^2$. Second, by Lemma 1 applied to $(d_1^((t+1)), dots, d_n^((t+1)))$, the newly concentrated subset $A^((t+1)) in cal(A)_h$ satisfies $sum_(i in A^((t+1))) d_i^((t+1)) <= sum_(i in A^((t))) d_i^((t+1))$. Chaining both inequalities yields
  $ J_alpha (A^((t+1)), cal(C)^((t+1)), bold(mu)^((t+1))) &= 1 / h sum_(i in A^((t+1))) d_i^((t+1)) \
    &<= 1 / h sum_(i in A^((t))) d_i^((t+1)) \
    &<= J_alpha (A^((t)), cal(C)^((t)), bold(mu)^((t+1))). $

  To establish finite termination without cycling, define the discrete state $Omega^((t)) := (A^((t)), cal(C)^((t)))$ and its reduced energy $J_alpha^* (Omega^((t))) := J_alpha (Omega^((t)), bold(mu)^((t+1)))$. Suppose $J_alpha^* (Omega^((t+1))) = J_alpha^* (Omega^((t)))$. Combining Claims 1 and 2 forces zero drop in every substep: first, $norm(mu_j^((t+2)) - mu_j^((t+1)))_2 = 0$ for all $j$, so $bold(mu)^((t+2)) = bold(mu)^((t+1))$; second, if $bold(mu)^((t+1)) = bold(mu)^((t))$, then the residual vector $d^((t+1)) = d^((t))$ is unchanged, so deterministic tie-breaking by smallest cluster index $j$ and smallest point index $i$ reproduces $Omega^((t+1)) = Omega^((t))$ identically, triggering immediate termination. Otherwise $bold(mu)^((t+1)) != bold(mu)^((t))$ would force a strictly positive drop in Claim 1 at step $t$. Hence every non-terminal iteration strictly decreases $J_alpha^* (Omega)$, so no state $Omega = (A, cal(C))$ can recur. Because there are $binom(n, h)$ choices for $A in cal(A)_h$ and at most $sum_(m=1)^k S(h, m)$ partitions of $A$ into up to $k$ non-empty subsets, Algorithm 1 terminates in finitely many steps at a trimmed centroidal equilibrium.
]

#subsec("B", "Exact Finite-Sample Breakdown Point Characterization")

Unlike univariate location estimators where trimming fraction $alpha$ unconditionally yields breakdown point $(floor(alpha n) + 1) / n$, clustering breakdown depends on *both* $alpha$ and cluster separation, because an adversary can concentrate $m <= floor(alpha n)$ outliers at a distant point to hijack a centroid. For any subset $S subset.eq cal(X)$ and integers $q >= 1$, $s <= |S|$, define the optimal $(q, s)$-trimmed inertia on $S$:
$ &V_S (q, s) := min_(A subset.eq S, |A| = s) min_(|bold(mu)| = q) \
  &quad sum_(i in A) min_(1 <= j <= q) norm(x_i - mu_j)_2^2. $

#theorem(title: "Exact Finite-Sample Breakdown Point of Trimmed K-Means", origin: "Original Proof")[
  Let $cal(X) subset.eq bb(R)^d$ lie in a bounded ball $overline(B)(0, R)$, let $r = floor(alpha n)$, $h = n - r >= k >= 2$, and let $epsilon^* (cal(X))$ denote the finite-sample Donoho–Huber breakdown point of the global minimizer $hat(bold(mu))$ of $F_alpha$ (the smallest fraction $m / n$ of sample replacements able to drive $max_j norm(hat(mu)_j)_2 -> infinity$). Then:
  1. *Universal upper bound:* $epsilon^* (cal(X)) <= (r + 1) / n$.
  2. *Exact threshold for $m <= r$:* The global minimizer resists all $m$-point adversarial corruptions if every $(n - m)$-subset $S subset.eq cal(X)$ satisfies
     $ V_S (k, h) < V_S (k - 1, h - m). $
     Conversely, if $V_S (k - 1, h - m) < V_S (k, h)$ for some $(n - m)$-subset $S subset.eq cal(X)$, then $epsilon^* (cal(X)) <= m / n$.
]

#proof[
  (1) Replace $m = r + 1$ points by identical corrupted copies $x'_1 = dots = x'_(r+1) = M e_1$ with $M -> infinity$, where $e_1$ is the first standard basis vector. Let $cal(O)$ be the corrupted index set ($|cal(O)| = r + 1$). For any active set $A in cal(A)_h$ ($|A| = n - r$), the pigeonhole principle forces
  $ |A inter cal(O)| >= |A| + |cal(O)| - n = (n - r) + (r + 1) - n = 1. $
  Choose $u in A inter cal(O)$ belonging to active cluster $C_j subset.eq A$. Let $p = |C_j inter cal(O)| >= 1$ and $q = |C_j backslash cal(O)| <= h - 1$. Since every clean point satisfies $|chevron.l x_i, e_1 chevron.r| <= R$ and every corrupted point satisfies $chevron.l x'_i, e_1 chevron.r = M$, the optimal centroid $hat(mu)_j = |C_j|^(-1) sum_(i in C_j) tilde(x)_i$ obeys
  $ chevron.l hat(mu)_j, e_1 chevron.r &= (p M + sum_(i in C_j backslash cal(O)) chevron.l x_i, e_1 chevron.r) / (p + q) \
    &>= M / h - R -> infinity. $
  Hence $max_j norm(hat(mu)_j)_2 -> infinity$, proving $epsilon^* (cal(X)) <= (r + 1) / n$.

  (2) Now fix $m <= r$, let $I_m$ be any corrupted index subset of size $m$, and let $S = cal(X) backslash I_m$ ($|S| = n - m >= h$). Trimming all $m$ corrupted points and solving optimal $(k, h)$-trimmed $k$-means on $S$ is always feasible, so the global optimal inertia on the corrupted dataset $tilde(cal(X))$ is bounded above by $V_S (k, h)$.

  Suppose $V_S (k, h) < V_S (k - 1, h - m)$ holds for every such $S$, and assume for contradiction that along some sequence of $m$-point corruptions, a global minimizer $hat(bold(mu))$ has $q >= 1$ diverging centroids ($norm(hat(mu)_j)_2 -> infinity$) while $k - q$ centroids stay bounded. No clean point $x_i in S subset.eq overline(B)(0, R)$ can be active and assigned to a diverging centroid once $norm(hat(mu)_j)_2 > R + sqrt(V_S (k, h))$, as that single point alone would contribute $> V_S (k, h)$ to total inertia. Consequently, all active clean points $A inter S$ must be covered by the $k - q$ bounded centroids. Because $|A inter I_m| <= m$, we have $|A inter S| >= h - m$, so total inertia is at least
  $ V_S (k - q, |A inter S|) >= V_S (k - 1, h - m) > V_S (k, h), $
  contradicting global optimality. Thus all $k$ centroids remain uniformly bounded.

  Conversely, suppose $V_S (k - 1, h - m) < V_S (k, h)$ for some $S = cal(X) backslash I_m$. Place all $m$ corrupted points identically at $z = M e_1$ with $M > 2 R + 2 sqrt(V_S (k, h))$. Placing one centroid at $z$ (zero error on $I_m$) and $k - 1$ centroids optimally on $h - m$ points of $S$ yields total inertia $V_S (k - 1, h - m)$. Yet any candidate keeping all $k$ centroids inside $overline(B)(0, M / 2)$ incurs squared distance $>= (M / 2)^2 > V_S (k, h)$ for any included point of $I_m$, so it must exclude all of $I_m$ and suffer inertia $>= V_S (k, h) > V_S (k - 1, h - m)$. Thus every global minimizer places a centroid outside $overline(B)(0, M / 2)$; sending $M -> infinity$ forces breakdown, so $epsilon^* (cal(X)) <= m / n$.
]

#remark(title: "Maximum Achievable Breakdown")[
  Theorem 2 explains why Trimmed $k$-Means attains $(floor(alpha n) + 1) / n$ on well-separated balanced clusters, whereas unbalanced or weakly separated clusters break down earlier under concentrated outlier clumps.
]

#subsec("C", "Existence and Moment-Free Strong Consistency")

A key statistical advantage of Trimmed $k$-Means over Lloyd's $k$-means is strong consistency even for heavy-tailed laws lacking finite mean or variance (e.g., multivariate Cauchy distributions).

#theorem(title: "Moment-Free Existence and Strong Consistency", origin: "Classical — Cited")[
  Let $P$ be any probability distribution on $bb(R)^d$ (with *no moment assumptions*), and fix $alpha in (0, 1)$ and $k >= 1$. Define the population trimmed objective @cuesta1997trimmed
  $ Phi_alpha (bold(mu)) &= 1 / (1 - alpha) inf_(P(A) >= 1 - alpha) \
    &quad integral_A min_(1 <= j <= k) norm(x - mu_j)_2^2 space d P(x). $
  1. There exists at least one global minimizer $bold(mu)^* in (bb(R)^d)^k$ of $Phi_alpha$, and the optimal trimmed region is a finite union of closed balls $B(bold(mu)^*, rho^*)$.
  2. Suppose $bold(mu)^*$ is unique up to cluster permutation and $P(partial B(bold(mu)^*, rho^*)) = 0$. For i.i.d. draws $X_1, dots, X_n tilde P$, empirical global minimizers $hat(bold(mu))_n$ of $F_alpha$ satisfy $hat(bold(mu))_n -> bold(mu)^*$ almost surely (up to permutation) and $F_alpha (hat(bold(mu))_n) -> Phi_alpha (bold(mu)^*)$ a.s. as $n -> infinity$ @cuesta1997trimmed.
]

#proof(label: "Proof Synthesis (Cited)")[
  Cuesta-Albertos, Gordaliza, and Matrán @cuesta1997trimmed first prove compactness without moment assumptions: because $alpha > 0$, there exists a deterministic radius $R_0 < infinity$ such that $P(overline(B)(0, R_0)) >= 1 - alpha / 2$. Any candidate $bold(mu)$ whose $k$ centroids all lie outside $overline(B)(0, 3 R_0)$ incurs squared distance at least $(2 R_0)^2$ on a set of mass $>= alpha / 2$ inside any $(1 - alpha)$-mass subset, whereas placing centroids inside $overline(B)(0, R_0)$ yields finite cost bounded by $4 R_0^2$. Refining this argument shows that a minimizing sequence admits a subsequence whose non-diverging centroids stay in a compact set and achieve the infimum; continuity of $Phi_alpha$ on compact sets then guarantees existence of $bold(mu)^*$.

  For strong consistency, because trimming restricts integration to the bounded union of balls $B(bold(mu), rho_n (bold(mu)))$ whose radius is uniformly bounded almost surely by the empirical $(1 - alpha / 2)$-quantile of $norm(X)_2$ plus $max_j norm(mu_j)_2$, Pollard's $L^2$ integrable envelope @pollard1981strong is replaced by a *deterministic bounded truncation*. The family of unions of $k$ closed balls in $bb(R)^d$ has finite Vapnik–Chervonenkis dimension, yielding uniform Glivenko–Cantelli convergence of empirical probabilities and truncated second moments over compact centroid domains. Uniqueness of $bold(mu)^*$ and boundary regularity $P(partial B(bold(mu)^*, rho^*)) = 0$ then imply almost-sure convergence of $hat(bold(mu))_n$ to $bold(mu)^*$ up to relabeling.
]

#subsec("D", "Over-Trimming Vulnerability and Minority-Cluster Erasure")

While Theorem 2 favors choosing $alpha$ above the true contamination rate, setting $alpha$ too aggressively can destroy legitimate small clusters regardless of how well separated they are.

#proposition(title: "Small-Cluster Erasure Under Over-Trimming", origin: "Original Counterexample")[
  Let $k = 2$ and let $p in (0, 1 / 2)$. For every separation distance $L > 0$, there exists a dataset $cal(X) subset.eq bb(R)$ of size $n$ containing a zero-variance minority cluster $C_"min"$ of proportion $p$ separated by distance $L$ from a bimodal bulk, such that whenever $alpha >= p$ (i.e., $r = floor(alpha n) >= |C_"min"|$), the unique global minimizer of $J_alpha$ discards $C_"min"$ completely as outliers and splits the bulk instead.
]

#proof[
  Let $n_2 = floor(p n) >= 1$ and $n_1 = n - n_2$, choosing $n_1$ even with $n_1 > r$. Construct $cal(X) subset.eq bb(R)$ by placing $n_1 / 2$ points at $-1$, $n_1 / 2$ points at $+1$, and all $n_2$ minority points at $L > 1$. Assume $alpha >= p$ so $r >= n_2$ and $h = n - r <= n_1$. Consider two centroid configurations:
  - *Bulk-splitting solution:* Set $bold(mu)^((1)) = (-1, +1)$. All $n_1 >= h$ bulk points have zero squared distance to $bold(mu)^((1))$, while all $n_2$ points at $L$ have squared distance $(L - 1)^2 > 0$. Because $h <= n_1$, selecting any $h$ bulk points as active discards all $n_2$ points at $L$ and attains $F_alpha (bold(mu)^((1))) = 0$.
  - *Minority-preserving solution:* Any $bold(mu)$ placing one centroid near $L$ (covering $C_"min"$) leaves only one centroid $mu_1$ to cover at least $h - n_2 >= 1$ bulk points. If $h - n_2 > n_1 / 2$, both $-1$ and $+1$ must appear in $A$, forcing
    $ F_alpha (bold(mu)) >= 1 / h (h - n_2 - n_1 / 2) dot 1^2 > 0, $
    since placing $mu_1$ at either $plus.minus 1$ leaves at least $h - n_2 - n_1 / 2$ active points at squared distance $4$ (and placing $mu_1 = 0$ gives squared distance $1$ on all $h - n_2$ points).
  Thus for *arbitrary* separation $L -> infinity$, $bold(mu)^((1))$ strictly beats every minority-preserving configuration and erases $C_"min"$ entirely.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Each iteration of Algorithm 1 performs three operations: (i) evaluating $n k$ squared Euclidean distances in $bb(R)^d$, costing $Theta(n k d)$ via BLAS-3 matrix products $X M^T$; (ii) identifying the $h$ smallest residuals, costing $O(n log n)$ with full sorting or $O(n)$ expected/worst-case using linear-time selection (Introselect / Quickselect); and (iii) accumulating $k$ active centroids in $O(h d) = O(n d)$ time. Total runtime is therefore $O(n (k d + log n) T)$ with sorting, or $O(n k d T)$ with linear-time order-statistic selection. Memory usage is $O((n + k) d)$—identical to Lloyd's $k$-means plus an $O(n)$ residual buffer.

Geometrically, Trimmed $k$-Means enforces a *shared-radius spherical truncation*: every cluster cell is $overline(B)(mu_j, rho) inter V_j (bold(mu))$ with the *same* global radius $rho = sqrt(d_((h)))$. Because squared Euclidean norms are preserved by orthogonal maps and translations, Trimmed $k$-Means is strictly invariant under rigid isometries $T(x) = Q x + b$ ($Q in O(d)$) and isotropic scaling, but remains sensitive to anisotropic feature scaling just like Alg. 1. Furthermore, minimizing $J_alpha$ is NP-hard even for $alpha = 0$ @aloise2009np @mahajan2012planar.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for Trimmed K-Means],
    [Time / Space Complexity], [$O(n(k d + log n)T)$ with sort, $O(n k d T)$ with Quickselect; $O((n+k)d)$ memory.],
    [Geometric Inductive Bias], [Equal-radius union of $k$ truncated spherical Voronoi balls $overline(B)(mu_j, rho) inter V_j (bold(mu))$.],
    [Determinism / Seeding], [Deterministic given $bold(mu)^((0))$; multi-start concentration steps strongly recommended.],
    [Invariances], [Invariant under $Q in O(d)$, translation, and isotropic scale; sensitive to anisotropic scaling.],
    [Outlier Breakdown Point], [Up to $(floor(alpha n) + 1) / n$ on separated balanced clusters (Theorem 2); moment-free consistency.],
    [Recommended Upgrades], [Student's $t$-Mixtures (Alg. 88) for elliptical/soft tails; PAM (Alg. 4) or K-Medians (Alg. 15) for $L_1$; DBSCAN (Alg. 31) for non-convex shapes.],
  ),
  caption: [Summary of theoretical guarantees and operational profile for Trimmed K-Means (Alg. 26).],
)

#sec("IV", "Engineering and Scientific Applications")

Trimmed $k$-Means excels on approximately spherical clusters contaminated by sensor glitches, heavy-tailed noise, or adversarial anomalies whenever hard outlier flagging and arithmetic prototypes are required.

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Industrial IoT Telemetry & Fault Isolation:* Machinery sensor streams (vibration, pressure, thermal) suffer transient ADC spikes and dropouts. Theorem 3 guarantees moment-free consistency while $O(n)$ selection isolates clean operating modes and quarantines fault spikes into $A^c$.
2. *Financial Regime & Order-Book Clustering:* Asset returns and limit-order snapshots exhibit power-law tails where sample means suffer variance inflation. Impartial trimming sets flash-crash influence to zero outside $B(bold(mu), rho)$, stabilizing regime centroids for risk budgeting.
3. *Astronomical Photometry & Cosmic-Ray Rejection:* Color-magnitude catalogs of star clusters contain background galaxies and cosmic-ray hits. Truncating cells to radius $rho$ recovers stellar centroids while outputting candidate anomalies $A^c$.
4. *Clinical Biomarker & Flow Cytometry Gating:* Single-cell fluorescence panels contain debris and doublets outside biological populations. Impartial trimming removes debris prior to computing mean fluorescence intensity centroids, preserving clinical interpretability.
5. *Geospatial Fleet & Base-Station Zoning:* Urban GPS logs exhibit multipath reflections kilometers outside service corridors. Standard $k$-means drags depot centroids toward outliers, whereas Trimmed $k$-Means locks prototypes onto the core $(1 - alpha)$ demand mass.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Trimmed $k$-Means exhibits four primary failure modes: (a) *over-trimming minority populations*, where $alpha$ exceeds the smallest genuine cluster frequency (Proposition 1); (b) *under-trimming heavy contamination*, where contamination exceeds $floor(alpha n)$ and triggers centroid breakdown (Theorem 2(1)); (c) *heteroscedastic or elongated clusters*, because the shared global radius $rho$ disproportionately trims legitimate boundary points from high-variance clusters while retaining noise near compact clusters (resolved by robust Student's $t$-Mixtures, Alg. 88 @peel2000robust, or TCLUST-style covariance constraints); and (d) *non-convex manifolds*, where convex Voronoi truncation slices arcs and filaments (resolved by DBSCAN Alg. 31 @ester1996dbscan). Selecting $alpha$ in practice relies on inspecting the *trimmed-elbow curve* $alpha |-> F_alpha (hat(bold(mu))_alpha)$ for sharp knees.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Why Trimmed K-Means Fits], [Caution / Alternative],
    [IoT / Telemetry], [Zero influence for extreme spikes; sub-ms speed], [Tune $alpha$ via elbow curve],
    [Quantitative Finance], [Moment-free consistency on fat-tailed returns], [Use Alg. 88 for elliptical tails],
    [Astronomy / Imaging], [Joint core quantization and artifact flagging], [Check small-class erasure],
    [Flow Cytometry], [Clean MFI centroids without debris distortion], [Use Alg. 88 if cluster radii vary],
    [Fraud / Rare Classes], [Flags extreme anomalies cleanly into $A^c$], [Avoid if fraud is target cluster],
  ),
  caption: [Domain suitability matrix and engineering guidance for Trimmed K-Means.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Trimmed $k$-Means bridges hard centroid quantization (Alg. 1) @lloyd1982least, redescending robust M-estimation, and contaminated mixture modeling (Algs. 81 and 88) @peel2000robust.

#proposition(title: "Redescending M-Estimator Equivalence and Zero-Temperature Noise-Mixture Limit", origin: "Original Derivation")[
  Let $(A^*, cal(C)^*, bold(mu)^*)$ be a strict trimmed centroidal equilibrium with unique radius $rho^* = sqrt(d_((h)) (bold(mu)^*)) < sqrt(d_((h+1)) (bold(mu)^*))$. Then:
  1. $bold(mu)^*$ is a stationary point of the *bounded-loss redescending M-estimator*
     $ L_(rho^*) (bold(mu)) = sum_(i=1)^n min ( min_(1 <= j <= k) norm(x_i - mu_j)_2^2, space (rho^*)^2 ), $
     whose influence function vanishes identically outside $B(bold(mu)^*, rho^*)$.
  2. Consider a $(k + 1)$-component mixture of $k$ equal-weight isotropic Gaussians $cal(N)(mu_j, sigma^2 I_d)$ plus a uniform noise component of density $delta_sigma = (2 pi sigma^2)^(-d/2) exp(-(rho^*)^2 / (2 sigma^2))$. As $sigma -> 0^+$, the hard-assignment MAP log-likelihood recovers $L_(rho^*) (bold(mu))$ exactly.
]

#proof[
  (1) Because $d_((h)) (bold(mu)^*) < (rho^*)^2 + epsilon < d_((h+1)) (bold(mu)^*)$ for sufficiently small perturbation, in an open neighborhood $cal(U)$ of $bold(mu)^*$ the active set $A(bold(mu)) equiv A^*$ and Voronoi assignment $cal(C)(bold(mu)) equiv cal(C)^*$ remain constant, with $d_i (bold(mu)) < (rho^*)^2$ for $i in A^*$ and $d_i (bold(mu)) > (rho^*)^2$ for $i in.not A^*$. Thus on $cal(U)$,
  $ L_(rho^*) (bold(mu)) = sum_(j=1)^k sum_(i in C_j^*) norm(x_i - mu_j)_2^2 + (n - h) (rho^*)^2. $
  Differentiating yields $nabla_(mu_j) L_(rho^*) (bold(mu)^*) = 2 sum_(i in C_j^*) (mu_j^* - x_i) = bold(0)$, while points $x_i in.not B(bold(mu)^*, rho^*)$ contribute constant $(rho^*)^2$ and zero gradient (exact redescending influence).

  (2) Taking $-2 sigma^2$ times the maximum component log-density for point $x_i$ gives, after dropping the shared $-d sigma^2 log(2 pi sigma^2)$ term,
  $ &-2 sigma^2 max ( max_j log p_j (x_i), space log delta_sigma ) \
    &quad = min ( min_j norm(x_i - mu_j)_2^2, space (rho^*)^2 ). $
  Summing over $i = 1, dots, n$ recovers $L_(rho^*) (bold(mu))$ identically as $sigma -> 0^+$, proving that Trimmed $k$-Means is the hard zero-temperature limit of a uniform-noise contaminated Gaussian mixture model.
]

At $alpha = 0$, Algorithm 1 collapses identically to Lloyd's $k$-means (Alg. 1). Compared to $L_1$-based robust quantizers such as PAM (Alg. 4) @kaufman1990finding and $k$-Medians (Alg. 15) @bradley1997clustering—whose influence functions are bounded but non-zero at infinity—Trimmed $k$-Means is *strictly redescending*, completely eliminating extreme outlier leverage at the price of tuning $alpha$. Compared to Power $k$-Means (Alg. 27) @xu2019power and W-$k$-Means (Alg. 28) @huang2005automated, Trimmed $k$-Means uniquely targets sample-axis contamination rather than local-minima smoothing or feature-axis noise.

In real empirical benchmarks on the ClusterOverview 4-dataset suite (using Silhouette @rousseeuw1987silhouettes and Adjusted Rand Index @hubert1985comparing), Trimmed $k$-Means attains Iris Sil $= 0.462$ (ARI $= 0.525$), Customers Sil $= 0.374$, Multi-Shape Sil $= 0.649$ (ARI N/A), and High-D Subspace+Noise Sil $= 0.350$ (ARI N/A), with a remarkably fast mean runtime of $0.081$ ms. Notably, it achieves the strongest Silhouette on High-D Subspace+Noise ($0.350$ vs. $0.303$ for Power $k$-Means and $0.243$ for W-$k$-Means) and Customers ($0.374$) because trimming discards ambient noise points that otherwise pull spherical centroids off-center, while its moderate Iris ARI ($0.525$) reflects Proposition 1's trade-off: trimming clean boundary points in overlapping classes slightly perturbs ground-truth class borders @jain2010data.
