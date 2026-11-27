#import "../template/ieee.typ": *

#alg_entry(
  id: 15,
  name: "K-Medians",
  category: "Cat. I: Centroid & Partitioning",
  year: "1997",
  paradigm: "Coordinate-Wise L1 Voronoi Partitioning",
  complexity: $O(n k d T)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Let $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$ be a finite dataset equipped with the Manhattan ($ell_1$) norm $norm(u)_1 = sum_(r=1)^d abs(u_r)$, and let $k in {1, dots, n}$ denote the target number of clusters. The *K-Medians clustering problem* @bradley1997clustering @kaufman1990finding seeks a disjoint partition $cal(C) = {C_1, dots, C_k}$ of ${1, dots, n}$ and prototype vectors $bold(m) = (m_1, dots, m_k) in (bb(R)^d)^k$ minimizing the *within-cluster sum of absolute deviations (WCSAD)*:
$ J_1 (cal(C), bold(m)) &= sum_(j=1)^k sum_(i in C_j) norm(x_i - m_j)_1 \
                        &= sum_(j=1)^k sum_(r=1)^d sum_(i in C_j) abs(x_(i, r) - m_(j, r)). $

Two structural distinctions separate K-Medians from neighboring robust methods. First, unlike K-Medoids / PAM (Alg. 4) @kaufman1990finding, the prototype $m_j in bb(R)^d$ is formed *coordinate-wise* and need not belong to $cal(X)$. Second, unlike the Euclidean geometric median (Fermat–Weber point, minimizing $sum norm(x_i - m_j)_2$ via iterative Weiszfeld updates), the $ell_1$ objective decouples additively across the $d$ coordinates, permitting exact closed-form $O(n d)$ prototype optimization via one-dimensional medians.

#definition(title: "Median Interval and L1-Voronoi Equilibrium")[
  For a non-empty multiset of scalars $S_r = {z_((1)) <= dots <= z_((m))} subset.eq bb(R)$, its *median interval* is the closed segment
  $ "Med"(S_r) := [z_((ceil(m / 2))), space z_((floor(m / 2) + 1))], $
  and its *canonical midpoint median* is $tilde(z) (S_r) = 1/2 (z_((ceil(m / 2))) + z_((floor(m / 2) + 1)))$. A pair $(cal(C), bold(m))$ is an *$L_1$-Voronoi equilibrium* if: (i) for every $i in C_j$, $norm(x_i - m_j)_1 <= norm(x_i - m_ell)_1$ for all $ell != j$, and (ii) for every non-empty cluster $C_j$ and coordinate $r in {1, dots, d}$, $m_(j, r) in "Med"({x_(i, r) : i in C_j})$.
]

When cluster cardinality $m = |C_j|$ is even, $"Med"(S_r)$ can be a non-degenerate interval of constant objective value. To eliminate flat-plateau cycling cleanly (proved in Theorem 1), Algorithm 1 uses an *inertia-preserving* coordinate median update: $m_(j, r)$ retains its previous value whenever it already lies inside $"Med"(S_r)$, and otherwise jumps to the canonical midpoint.

#algbox(
  title: "Alternating K-Medians Algorithm (L1 Partitioning)",
  inputs: [Dataset $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$, integer $k >= 1$, initial prototypes $bold(m)^((0))$, max iterations $T$.],
  outputs: [Partition $cal(C) = {C_1, dots, C_k}$ and $L_1$ prototypes $bold(m)$.],
)[
  1. *For* $t = 0, 1, dots, T - 1$ *do*:
  2. #h(1em) *L1 Assignment step:* For each $i = 1, dots, n$, assign
     $ c_i^((t)) = arg min_(j in {1, dots, k}) sum_(r=1)^d abs(x_(i, r) - m_(j, r)^((t))), $
     breaking ties deterministically by smallest cluster index $j$.
  3. #h(1em) Form $C_j^((t)) = {i : c_i^((t)) = j}$ for $j = 1, dots, k$ (reseed empty clusters to the farthest $L_1$ point).
  4. #h(1em) *Coordinate-wise median update:* For each non-empty $C_j^((t))$ and each axis $r = 1, dots, d$:
     - Compute the order statistics $z_((ceil(m/2)))$ and $z_((floor(m/2)+1))$ of $S_(j, r) = {x_(i, r) : i in C_j^((t))}$ ($m = |C_j^((t))|$) in $O(m)$ time via linear-time selection.
     - *If* $m_(j, r)^((t)) in "Med"(S_(j, r))$ *then* set $m_(j, r)^((t+1)) <- m_(j, r)^((t))$; *else* set $m_(j, r)^((t+1)) <- tilde(z) (S_(j, r))$.
  5. #h(1em) *If* $t >= 1$ and $cal(C)^((t)) = cal(C)^((t-1))$ (and $bold(m)^((t+1)) = bold(m)^((t))$), *break*.
  6. *Return* $(cal(C)^((t)), bold(m)^((t+1)))$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five formal guarantees characterizing K-Medians: (1) exact subdifferential optimality and an integral excess-loss formula for coordinatewise medians; (2) strict Lyapunov non-cycling and finite convergence even under even-cardinality non-uniqueness; (3) exact $50\%$ per-cluster finite-sample breakdown robustness contrasted with global $1/n$ assignment sensitivity; (4) signed-permutation symmetry alongside explicit counterexamples proving rotational non-invariance and non-convexity of $L_1$-Voronoi cells; and (5) exact zero-scale equivalence to factorized Laplace mixtures plus classical NP-hardness, finite-first-moment strong consistency, and high-dimensional $L_1$ contrast.

#subsec("A", "Subdifferential Optimality and Finite-Step Convergence")

#lemma(title: "Subgradient Optimality and Exact Integral Excess of the L1 Median", origin: "Original Proof")[
  For any non-empty finite set $S subset.eq bb(R)^d$, define $f_S (m) = sum_(x in S) norm(x - m)_1$. Then:
  1. $m^* in bb(R)^d$ globally minimizes $f_S$ if and only if $m^*_r in "Med"({x_r : x in S})$ for every $r = 1, dots, d$, equivalently:
     $ abs({x in S : x_r < m^*_r}) &<= |S| / 2, \
       abs({x in S : x_r > m^*_r}) &<= |S| / 2. $
  2. For any $m in bb(R)^d$ and any minimizer $m^*$,
     $ &f_S (m) - f_S (m^*) \
       &= sum_(r=1)^d integral_(m^*_r)^(m_r) (N^<_(S, r) (t) - N^>_(S, r) (t)) d t >= 0, $
     where $N^<_(S, r) (t) = abs({x in S : x_r < t})$ and $N^>_(S, r) (t) = abs({x in S : x_r > t})$. Moreover, $f_S (m) > f_S (m^*)$ strictly whenever $m_r in.not "Med"({x_r : x in S})$ for at least one coordinate $r$.
]

#proof[
  Because $f_S (m) = sum_(r=1)^d g_r (m_r)$ with convex one-dimensional terms $g_r (t) = sum_(x in S) abs(t - x_r)$, global minimization separates coordinate-by-coordinate. At any $t in bb(R)$, let $n_< = N^<_(S, r) (t)$, $n_> = N^>_(S, r) (t)$, and $n_= = abs({x in S : x_r = t})$, so $n_< + n_= + n_> = |S|$. By additivity of the convex subdifferential of $t |-> abs(t - x_r)$,
  $ partial g_r (t) &= sum_(x in S) partial abs(t - x_r) \
                     &= [n_< - n_> - n_=, space n_< - n_> + n_=]. $
  By Fermat's optimality condition for convex functions, $t$ minimizes $g_r$ if and only if $0 in partial g_r (t)$, i.e., $n_< - n_> - n_= <= 0 <= n_< - n_> + n_=$. Substituting $n_= = |S| - n_< - n_>$ transforms these two inequalities identically into
  $ 2 n_< <= |S| quad "and" quad 2 n_> <= |S|, $
  which holds if and only if $t$ lies in the closed interval $[z_((ceil(|S|/2))), z_((floor(|S|/2)+1))] = "Med"({x_r : x in S})$.

  For Claim 2, $g_r$ is Lipschitz continuous (with constant $|S|$) and piecewise linear, hence absolutely continuous with almost-everywhere derivative $g'_r (t) = N^<_(S, r) (t) - N^>_(S, r) (t)$. The Fundamental Theorem of Calculus yields
  $ g_r (m_r) - g_r (m^*_r) = integral_(m^*_r)^(m_r) (N^<_(S, r) (t) - N^>_(S, r) (t)) d t. $
  Let $[a_r, b_r] = "Med"({x_r : x in S})$. For every $t > b_r$, $N^<_(S, r) (t) > |S|/2 >= N^>_(S, r) (t)$, so the integer-valued integrand obeys $N^<_(S, r) (t) - N^>_(S, r) (t) >= 1$; symmetrically, for every $t < a_r$, $N^<_(S, r) (t) - N^>_(S, r) (t) <= -1$. Thus every coordinate integral is non-negative, and if $m_r > b_r$ (or $m_r < a_r$), then $g_r (m_r) - g_r (m^*_r) >= m_r - b_r > 0$ (respectively $a_r - m_r > 0$), establishing strict positivity whenever $m_r in.not "Med"$.
]

#theorem(title: "Strict Lyapunov Non-Cycling and Finite Convergence", origin: "Original Proof")[
  Let ${(cal(C)^((t)), bold(m)^((t)))}$ be generated by Algorithm 1 with smallest-index assignment tie-breaking and inertia-preserving median updates. Then:
  1. $J_1 (cal(C)^((t)), bold(m)^((t+1))) <= J_1 (cal(C)^((t)), bold(m)^((t)))$, with strict inequality whenever $bold(m)^((t+1)) != bold(m)^((t))$.
  2. $J_1 (cal(C)^((t+1)), bold(m)^((t+1))) <= J_1 (cal(C)^((t)), bold(m)^((t+1)))$, with strict inequality whenever any datum moves to a strictly closer $L_1$ prototype.
  3. Algorithm 1 cannot cycle and terminates in at most $S(n, k) <= k^n / k!$ iterations at an $L_1$-Voronoi equilibrium.
]

#proof[
  Applying Lemma 1 independently to each non-empty cluster $C_j^((t))$ proves Claim 1: if $bold(m)^((t+1)) != bold(m)^((t))$, then by Step 4's inertia-preserving rule at least one coordinate satisfies $m_(j, r)^((t)) in.not "Med"(S_(j, r))$, so Lemma 1(2) forces a *strict* decrease $J_1 (cal(C)^((t)), bold(m)^((t+1))) < J_1 (cal(C)^((t)), bold(m)^((t)))$. Claim 2 follows by summing $norm(x_i - m_(c_i^((t+1)))^((t+1)))_1 <= norm(x_i - m_(c_i^((t)))^((t+1)))_1$ over $i = 1, dots, n$.

  To rule out cycling, suppose $J_1 (cal(C)^((t+1)), bold(m)^((t+2))) = J_1 (cal(C)^((t)), bold(m)^((t+1)))$. Then both intermediate inequalities must be equalities. By Claim 1, equality in the update step implies $bold(m)^((t+2)) = bold(m)^((t+1))$. By Claim 2, every point $x_i$ achieves identical $L_1$ distance to $m_(c_i^((t+1)))^((t+1))$ and $m_(c_i^((t)))^((t+1))$. Because Step 2 deterministically chooses the smallest minimizer index under the *same* prototypes $bold(m)^((t+1))$, after one assignment pass no further switch can occur unless prototypes move, so $cal(C)^((t+2)) = cal(C)^((t+1))$ and Step 5 halts. Hence every visited partition prior to termination attains a strictly smaller value of $J_1^*(cal(C)) := min_(bold(m)) J_1 (cal(C), bold(m))$ and cannot repeat, bounding iterations by the number of non-empty $k$-partitions $S(n, k) <= k^n / k!$.
]

#subsec("B", "Finite-Sample Breakdown Robustness Against Outliers")

#theorem(title: "Per-Cluster 50% Breakdown vs. Global Partition Sensitivity", origin: "Original Proof")[
  Let $C_j$ be a cluster of size $m = |C_j| >= 1$ whose uncorrupted observations lie in $[-R, R]^d$ ($R > 0$).
  1. *(Prototype breakdown)* Under arbitrary adversarial corruption of $q$ points in $C_j$, the midpoint coordinatewise median $m_j$ remains inside $[-R, R]^d$ for all $q <= ceil(m / 2) - 1$, whereas $q = ceil(m / 2)$ outliers can drive $norm(m_j)_1 -> infinity$. Hence the exact within-cluster Donoho–Huber breakdown point is
     $ epsilon^*_j = ceil(m / 2) / m = floor((m + 1) / 2) / m -> 1 / 2 quad (m -> infinity). $
  2. *(Global $k$-partition caveat)* For $k >= 2$ without trimming or minimum-size constraints, moving even a *single* datum $x_1$ to $norm(x_1)_1 -> infinity$ causes the global $J_1$ minimizer to assign $x_1$ to a dedicated singleton cluster, so the unconstrained global prototype set still exhibits $1 / n$ breakdown.
]

#proof[
  (1) Fix coordinate $r in {1, dots, d}$ and let $z_((1)) <= dots <= z_((m))$ be the sorted values of $C_j$ after replacing $q$ elements arbitrarily. Since at least $m - q$ values remain in $[-R, R]$, at most $q$ values lie strictly below $-R$ and at most $q$ lie strictly above $R$. The midpoint median averages $z_((a))$ and $z_((b))$ where $a = ceil(m / 2)$ and $b = floor(m / 2) + 1$. Note that $a + b = m + 1$. If $q <= ceil(m / 2) - 1 = a - 1$, then:
  - $a > q$, so at least one value in ${z_((1)), dots, z_((a))}$ is $>= -R$, implying $z_((a)) >= -R$;
  - $m - b = m - (m + 1 - a) = a - 1 >= q$, so at least $q + 1$ values lie at or above $z_((b))$, implying $z_((b)) <= R$.
  Since $z_((a)) <= z_((b))$, both order statistics lie in $[-R, R]$, so $m_(j, r) = 1/2 (z_((a)) + z_((b))) in [-R, R]$ for all $r$. Conversely, setting $q = ceil(m / 2)$ corrupted points to $M e_1$ ($M > R$) forces $z_((b)) = M$ because $m - b + 1 = a = q$, driving $m_(j, 1) >= (M - R) / 2 -> +infinity$ as $M -> infinity$. Thus $epsilon^*_j = ceil(m / 2) / m$.

  (2) Let $cal(X)' = {x_2, dots, x_n} subset.eq [-R, R]^d$ and let $J^*_(k-1)$ be the optimal $(k-1)$-medians cost on $cal(X)'$. Place $x_1 = M e_1$ with $M > J^*_(k-1) + 2 d R$. Any partition placing $x_1$ in a non-singleton cluster $C_j$ alongside at least one $x_ell in cal(X)'$ pays at least $norm(x_1 - m_j)_1 + norm(x_ell - m_j)_1 >= norm(x_1 - x_ell)_1 >= M - d R > J^*_(k-1)$, whereas assigning ${x_1}$ as its own singleton cluster and clustering $cal(X)'$ with $k - 1$ clusters attains total cost $0 + J^*_(k-1) = J^*_(k-1)$. Hence the global minimizer must place a prototype at $m_1 = x_1 = M e_1 -> infinity$.
]

#subsec("C", "Symmetry Group, Rotational Non-Invariance, and Non-Convex Cells")

#proposition(title: "Hyperoctahedral Symmetry, Rotation Sensitivity, and Non-Convexity", origin: "Original Counterexample")[
  Let $B_d subset.eq O(d)$ be the hyperoctahedral group of $2^d d!$ signed permutation matrices.
  1. $J_1$ and its equilibria are strictly invariant under translations, uniform scalings, and every $P in B_d$.
  2. For all $d >= 2$, optimal K-Medians partitions are *not invariant* under generic orthogonal rotations $Q in O(d) backslash B_d$.
  3. For all $d >= 2$, $L_1$-Voronoi cells $V_j (bold(m)) = {x : norm(x - m_j)_1 <= norm(x - m_ell)_1 space forall ell}$ are generally *non-convex*.
]

#proof[
  (1) Any $P in B_d$ permutes and flips coordinate signs, preserving $norm(P u)_1 = norm(u)_1$; moreover, coordinatewise medians commute with sign negation, permutation, translation, and uniform scaling.

  (2) For $d = 2$ (embeddable in $bb(R)^d$ by zero-padding) and $k = 2$, construct the three-point dataset
  $ x_1 = (0, 0), quad x_2 = (6, 0), quad x_3 = (-4, 4). $
  Any 2-partition groups two points together and leaves one singleton; for any pair ${u, v}$, choosing $m = u$ achieves $J_1 ({u, v}) = norm(u - v)_1$. Under $Q = I_2$:
  $ norm(x_2 - x_1)_1 = 6, quad norm(x_3 - x_1)_1 = 8, quad norm(x_3 - x_2)_1 = 14, $
  so the *unique* optimal partition is ${{x_1, x_2}, {x_3}}$ with $J_1^* = 6$. Now apply the $45 degree$ orthogonal rotation $Q = 1/sqrt(2) mat(1, -1; 1, 1)$:
  $ Q x_1 = (0, 0), quad Q x_2 = (3 sqrt(2), 3 sqrt(2)), quad Q x_3 = (-4 sqrt(2), 0). $
  The rotated pairwise $L_1$ costs become:
  $ norm(Q(x_2 - x_1))_1 = 6 sqrt(2), quad norm(Q(x_3 - x_1))_1 = 4 sqrt(2), $
  and $norm(Q(x_3 - x_2))_1 = 10 sqrt(2)$. Thus the unique optimal partition flips strictly to ${{x_1, x_3}, {x_2}}$ ($4 sqrt(2) < 6 sqrt(2)$), proving orthogonal non-invariance.

  (3) Set $m_1 = (0, 0)$ and $m_2 = (6, 3)$ in $bb(R)^2$, and inspect $p = (0, 9)$ and $q = (4, -3)$:
  - $norm(p - m_1)_1 = 9 < 6 + 6 = 12 = norm(p - m_2)_1 => p in V_1 (bold(m))$;
  - $norm(q - m_1)_1 = 7 < 2 + 6 = 8 = norm(q - m_2)_1 => q in V_1 (bold(m))$.
  Yet their midpoint $z = 1/2 (p + q) = (2, 3)$ satisfies
  $ norm(z - m_1)_1 = 2 + 3 = 5 > 4 + 0 = 4 = norm(z - m_2)_1, $
  so $z in V_2 (bold(m)) backslash V_1 (bold(m))$, proving $V_1 (bold(m))$ is non-convex.
]

#subsec("D", "Laplace Mixture Limit, Consistency, and Complexity Guarantees")

#proposition(title: "Zero-Scale Limit of Factorized Laplace Mixtures", origin: "Original Derivation")[
  Consider an equal-weight $k$-component isotropic Laplace mixture $p(x | bold(m), b) = k^(-1) (2 b)^(-d) sum_(j=1)^k exp(- norm(x - m_j)_1 / b)$ with scale $b > 0$. As $b -> 0^+$, the scaled negative log-likelihood $-b sum_(i=1)^n log p(x_i | bold(m), b) - n b log(k (2 b)^d)$ converges uniformly to $min_(cal(C)) J_1 (cal(C), bold(m))$, and the EM algorithm reduces strictly to Algorithm 1.
]

#proof[
  For each $x_i$, let $d_i^* = min_(1 <= j <= k) norm(x_i - m_j)_1$. Bounding the finite exponential sum by its maximum term gives
  $ exp(- d_i^* / b) <= sum_(j=1)^k exp(- norm(x_i - m_j)_1 / b) <= k exp(- d_i^* / b). $
  Taking $-b log(dot)$ and summing over $i = 1, dots, n$ yields
  $ min_(cal(C)) J_1 (cal(C), bold(m)) - n b log k &<= L_b (bold(m)) \
                                                   &<= min_(cal(C)) J_1 (cal(C), bold(m)), $
  where $L_b (bold(m)) = -b sum_(i=1)^n log p(x_i | bold(m), b) - n b log(k (2 b)^d)$. Thus $sup_(bold(m)) abs(L_b (bold(m)) - min_(cal(C)) J_1 (cal(C), bold(m))) <= n b log k -> 0$ as $b -> 0^+$. Whenever $c_i$ is the unique $L_1$-nearest prototype with margin $Delta_i > 0$, the E-step posterior obeys $gamma_(i, c_i) (b) >= (1 + (k-1) exp(-Delta_i / b))^(-1) -> 1$, while the M-step minimizes $sum_i gamma_(i j) norm(x_i - m_j)_1$, recovering the exact coordinatewise median by Lemma 1 in the hard limit.
]

#theorem(title: "NP-Hardness, Finite-First-Moment Consistency, and L1 Contrast", origin: "Classical — Cited")[
  1. *NP-Hardness:* Globally minimizing $J_1 (cal(C), bold(m))$ is NP-hard for general $(k, d)$ @bradley1997clustering @harpeled2004coresets (parallel to Euclidean $k$-means hardness @aloise2009np @mahajan2012planar).
  2. *Strong Consistency under Infinite Variance:* Let $X_1, X_2, dots$ be i.i.d. draws from $P$ on $bb(R)^d$ satisfying only $bb(E)[norm(X)_1] < infinity$. If $Phi_1 (bold(m)) = bb(E)[min_j norm(X - m_j)_1]$ has a unique minimizer $bold(m)^*$ up to permutation, then the empirical minimizer $hat(bold(m))_n -> bold(m)^*$ almost surely as $n -> infinity$ @pollard1981strong.
  3. *High-Dimensional Distance Contrast:* Under independent coordinates, the relative distance contrast $(D_"max" - D_"min") / D_"min"$ decays strictly more slowly for $L_1$ than for $L_2$ or higher $L_p$ norms ($p >= 2$) @beyer1999nearest @deamorim2012minkowski.
]

#proof(label: "Proof Synthesis (Cited)")[
  Bradley, Mangasarian, and Street @bradley1997clustering formulate K-Medians as minimizing a concave piecewise-linear function over a polyhedron, which is globally NP-hard. For strong consistency, Pollard's general M-estimation framework @pollard1981strong requires only an integrable envelope for $g_(bold(m)) (x) = min_j norm(x - m_j)_1$ on compact balls $B_R (0)^k$; by the triangle inequality, $0 <= g_(bold(m)) (x) <= norm(x)_1 + R$, so $bb(E)[norm(X)_1] < infinity$ suffices—allowing heavy-tailed laws (e.g., Student-$t$ with $1 < nu <= 2$) where K-Means fails to be consistent due to infinite variance. Finally, Beyer et al. @beyer1999nearest and de Amorim and Mirkin @deamorim2012minkowski demonstrate that fractional/lower $p$ norms preserve larger normalized spread between nearest and farthest neighbors in high dimensions.
]

#sec("III", "Complexity, Geometry and Operational Profile")

The assignment step of Algorithm 1 evaluates $n k$ Manhattan distances in $O(n k d)$ operations, matching Lloyd's asymptotic assignment cost while avoiding squaring instructions (exploiting hardware absolute-difference SIMD primitives). In the update step, computing exact medians across $d$ axes via naive sorting costs $O(d sum_j |C_j| log |C_j|) <= O(n d log n)$, whereas using randomized Quickselect (expected) or median-of-medians BFPRT (deterministic worst-case) reduces the entire update step to $O(n d)$ total time, yielding overall complexity $O(n k d T)$ with $O((n + k) d)$ memory.

Geometrically, $L_1$ unit balls are cross-polytopes (octahedra / axis-aligned diamonds), and Voronoi boundaries are piecewise-linear bent manifolds aligned with coordinate hyperplanes and $45 degree$ diagonals (non-convex by Proposition 1(3)). This induces a strong *axis-aligned inductive bias*: K-Medians naturally captures sparse, feature-interpretable clusters and resists coordinate-wise spike outliers, at the price of rotational dependence.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for K-Medians],
    [Time / Space Complexity], [$O(n k d T)$ with linear-time selection ($O((n k + n log n) d T)$ with sort); $O((n+k)d)$ space.],
    [Geometric Inductive Bias], [Axis-aligned $L_1$ cross-polytope cells; piecewise-linear non-convex Voronoi regions.],
    [Moment / Tail Requirement], [Strongly consistent under $bb(E)[norm(X)_1] < infinity$ (tolerates infinite variance).],
    [Invariances], [Invariant under translation, uniform scale, and signed permutations $B_d$; *not* rotation-invariant.],
    [Breakdown Robustness], [$epsilon^*_j = ceil(|C_j|/2)/|C_j| -> 50\%$ per cluster; global unconstrained partition still $1/n$ (Thm. 2).],
    [Recommended Companions], [K-Medoids / FasterPAM (Alg. 4) for arbitrary metrics; Trimmed K-Means (Alg. 26) for global outlier rejection.],
  ),
  caption: [Summary of theoretical guarantees and operational profile of K-Medians.],
)

#sec("IV", "Engineering and Scientific Applications")

K-Medians is mathematically suited to domains where features represent distinct, semantically meaningful coordinate axes subject to heavy-tailed spikes, Laplace noise, or high-dimensional sparsity.

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Financial, Retail, and E-Commerce Customer Segmentation:* Transaction features (spend, visit frequency, basket size, return counts) follow skewed, heavy-tailed distributions with whale accounts and logging glitches. By Theorem 2(1), per-cluster medians remain completely anchored by the bulk majority up to $50\%$ contamination, yielding realistic "typical customer" personas rather than whale-distorted means.
2. *Urban Logistics and Manhattan-Grid Facility Location:* In rectilinear road networks, circuit-board routing, and warehouse picking aisles, true travel distance equals the $ell_1$ metric $norm(x - m)_1$. Minimizing $J_1$ directly optimizes total rectilinear delivery distance, making K-Medians the exact continuous relaxation of the $k$-median facility location problem.
3. *High-Dimensional Gene Expression and Proteomics Profiling:* Microarray and single-cell expression matrices exhibit high ambient dimension with technical dropout spikes and heavy-tailed batch artifacts. Both the finite-first-moment consistency (Theorem 3(2)) and superior $L_1$ high-dimensional contrast (Theorem 3(3)) stabilize cluster assignments against noisy gene axes.
4. *Robust Sensor Networks and Impulse-Noise Telemetry:* Industrial IoT and astronomical time-series sensors frequently experience transient bit-flip spikes modeled by Laplace or heavy-tailed noise (Proposition 2). Whereas squared-error K-Means lets a single spike hijack a cluster centroid, K-Medians ignores spike magnitudes entirely once a cluster has at least 3 members.
5. *Sparse Bag-of-Words and Count-Vector Clustering:* On nonnegative frequency tables and count features, coordinatewise medians preserve sparsity (returning exact zeros whenever $> 50\%$ of cluster members lack a feature), producing cleaner, sparser cluster prototypes than dense arithmetic means.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

K-Medians degrades in four settings: (a) *obliquely correlated or rotated clusters*, where Proposition 1(2) shows axis-aligned $L_1$ diamonds misalign with tilted covariances (favoring GMM Alg. 81 or PCA-whitened K-Means); (b) *non-convex manifolds* (arcs, rings), where piecewise-linear $L_1$ boundaries still fail to wrap around curves (favoring Kernel K-Means Alg. 16 or Spectral Clustering); (c) *extreme isolated outliers with unconstrained $k$*, which trigger Theorem 2(2)'s singleton-hijack failure (resolved by Trimmed K-Means Alg. 26 or DBSCAN Alg. 31); and (d) *non-vector/graph dissimilarities*, where coordinatewise medians are undefined (resolved by PAM / FasterPAM Alg. 4 @schubert2019fasterpam).

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Why K-Medians Fits / Mechanism], [Caution / Alternative],
    [Retail / FinTech Cohorts], [Resists heavy-tailed whale outliers (Thm. 2)], [Trim singletons via Alg. 26],
    [Urban / VLSI Routing], [Exact rectilinear $ell_1$ travel objective], [Use PAM (Alg. 4) on graphs],
    [High-D Omics / Noise], [Superior $L_1$ contrast; finite 1st moment], [Rotate/whiten if oblique],
    [IoT / Impulse Telemetry], [Laplace MLE limit (Prop. 2); spike-proof], [Use DTW/k-Shape for lag],
    [Rotated / Curved Shapes], [Fails rotation invariance (Prop. 1(2))], [Use GMM (Alg. 81) / Kernel KM],
  ),
  caption: [Domain suitability and engineering trade-off matrix for K-Medians.],
)

#sec("V", "Structural Relations and Empirical Behavior")

K-Medians bridges Euclidean prototype partitioning (Alg. 1), discrete medoid clustering (Alg. 4), and general Minkowski $L_p$-clustering @deamorim2012minkowski through exact norm-equivalence and factor-2 approximation inequalities.

#proposition(title: "Exact Comparison Bounds Relating K-Medians, K-Means, and K-Medoids", origin: "Original Unified Proof")[
  For any dataset $cal(X) subset.eq bb(R)^d$ and any non-empty $k$-partition $cal(C) = {C_1, dots, C_k}$, let $J_1^* (cal(C))$ denote the optimal K-Medians cost, $J_2^* (cal(C))$ the optimal K-Means WCSS (Alg. 1), and $J_(1, "med")^* (cal(C)) = sum_(j=1)^k min_(y in C_j) sum_(i in C_j) norm(x_i - y)_1$ the optimal $L_1$ K-Medoids cost (Alg. 4). Then:
  1. $J_2^* (cal(C))^(1/2) <= J_1^* (cal(C)) <= sqrt(n d) space J_2^* (cal(C))^(1/2)$.
  2. $J_1^* (cal(C)) <= J_(1, "med")^* (cal(C)) <= 2 J_1^* (cal(C))$.
]

#proof[
  (1) Let $m_j^*$ be optimal $L_1$ medians and $overline(x)_(C_j)$ be arithmetic centroids for $cal(C)$. Because $norm(u)_1 >= norm(u)_2$ and $(sum_i a_i)^2 >= sum_i a_i^2$ for $a_i >= 0$,
  $ (J_1^* (cal(C)))^2 &>= sum_(j=1)^k sum_(i in C_j) norm(x_i - m_j^*)_1^2 \
                       &>= sum_(j=1)^k sum_(i in C_j) norm(x_i - m_j^*)_2^2 >= J_2^* (cal(C)), $
  where the last step uses centroid optimality of $overline(x)_(C_j)$ for $J_2$. Taking square roots proves the lower bound. For the upper bound, suboptimality of $overline(x)_(C_j)$ for $J_1$ plus two applications of Cauchy–Schwarz ($norm(u)_1 <= sqrt(d) norm(u)_2$ and $sum_(i=1)^n b_i <= sqrt(n) (sum b_i^2)^(1/2)$) give
  $ J_1^* (cal(C)) &<= sum_(j=1)^k sum_(i in C_j) norm(x_i - overline(x)_(C_j))_1 \
                   &<= sqrt(d) sum_(j=1)^k sum_(i in C_j) norm(x_i - overline(x)_(C_j))_2 <= sqrt(n d) J_2^* (cal(C))^(1/2). $

  (2) Since $C_j subset.eq bb(R)^d$, $J_1^* (cal(C)) <= J_(1, "med")^* (cal(C))$ is immediate. For the factor-2 upper bound, pick in each cluster $C_j$ the datum $x_(i_j) = arg min_(x in C_j) norm(x - m_j^*)_1$, which satisfies $norm(x_(i_j) - m_j^*)_1 <= |C_j|^(-1) sum_(x in C_j) norm(x - m_j^*)_1$. By the $L_1$ triangle inequality,
  $ sum_(x in C_j) norm(x - x_(i_j))_1 &<= sum_(x in C_j) (norm(x - m_j^*)_1 + norm(m_j^* - x_(i_j))_1) \
                                      &<= 2 sum_(x in C_j) norm(x - m_j^*)_1. $
  Minimizing the left side over all $y in C_j$ and summing over $j = 1, dots, k$ proves $J_(1, "med")^* (cal(C)) <= 2 J_1^* (cal(C))$.
]

In the broader taxonomy, K-Medians corresponds to the $p = 1$ endpoint of Minkowski $L_p$-clustering @deamorim2012minkowski (interpolating smoothly to K-Means at $p = 2$), while Trimmed K-Means (Alg. 26) @cuesta1997trimmed and Student-$t$ mixtures @peel2000robust complement K-Medians by addressing global partition-level outlier breakdown (Theorem 2(2)).

*Empirical Behavior on the ClusterOverview Benchmark Suite:* On the real 4-dataset benchmark @rousseeuw1987silhouettes @hubert1985comparing, K-Medians exhibits a textbook reflection of its theoretical profile. On *High-D Subspace+Noise*, K-Medians achieves *Sil = 0.303* (ARI N/A)—nearly *double* Euclidean Hamerly's K-Means (0.158)—directly validating Theorem 3(3)'s superior $L_1$ high-dimensional contrast and noise resistance. On *Customers*, K-Medians reaches *Sil = 0.357* (improving on Euclidean 0.341), reflecting Theorem 2(1)'s robustness to skewed spending distributions. On *Iris*, K-Medians attains strong geometric compactness (*Sil = 0.478*) yet a lower label alignment (*ARI = 0.440* vs. 0.563 for Euclidean K-Means), consistent with Proposition 1(2)'s rotational sensitivity on obliquely correlated botanical axes, alongside *Multi-Shape Sil = 0.581* (ARI N/A) and fast mean runtime *0.147 ms*.
