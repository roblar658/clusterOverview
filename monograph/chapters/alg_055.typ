#import "../template/ieee.typ": *

#alg_entry(
  id: 55,
  name: "Mean Shift++",
  category: "Cat. II: Density & Mode-Seeking",
  year: "2021",
  paradigm: "Grid-Accelerated Linear-Time Mode-Seeking",
  complexity: $cal(O)(T n d 3^d)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Both standard Mean Shift @comaniciu2002mean and Blurring Mean Shift @cheng1995mean @carreira2006fast require $cal(O)(n^2 d)$ distance evaluations per iteration, creating a severe scalability barrier on large low-dimensional datasets such as megapixel images or dense LiDAR scans. Mean Shift++, proposed by Jang and Jiang @jang2021meanshiftpp, eliminates the quadratic bottleneck by replacing continuous pairwise kernel evaluations with a sparse Cartesian grid partition of side length $h > 0$. By hashing occupied cells and averaging points solely across each cell's $3^d$-cell Moore neighborhood, Mean Shift++ achieves strict linear complexity in $n$ for fixed dimension $d$ while preserving the nonparametric mode-seeking behavior of blurring-style Mean Shift.

Let $X^((0)) = {x_1^((0)), dots, x_n^((0))} subset bb(R)^d$ be a finite dataset and $h > 0$ be the grid bandwidth. Define the coordinate-wise *cell quantization map* $q_h : bb(R)^d -> bb(Z)^d$ by
$
  q_h (x) = (floor(x_1 / h), dots, floor(x_d / h)) in bb(Z)^d.
$
For any integer cell index $c in bb(Z)^d$, denote its half-open hypercube by $Q_h (c) = product_(r=1)^d [c_r h, (c_r + 1) h)$ and its $3^d$-cell Moore neighborhood by
$
  cal(N)_1 (c) = c + {-1, 0, 1}^d = {c' in bb(Z)^d mid(|) norm(c' - c)_infinity <= 1}.
$
The geometric *grid stencil* of a point $x in bb(R)^d$ is the union of its enclosing cell and all adjacent grid cells:
$
  S_h (x) &= union.big_(c' in cal(N)_1 (q_h (x))) Q_h (c') \
          &= {y in bb(R)^d mid(|) norm(q_h (x) - q_h (y))_infinity <= 1}.
$

#definition(title: "Mean Shift++ Hash-Grid Update Operator")[
  At iteration $t >= 0$, let $cal(H)^((t)) = {q_h (x_i^((t))) mid(|) i = 1, dots, n} subset bb(Z)^d$ be the set of occupied grid cells, with $M_t = |cal(H)^((t))| <= n$. For every occupied cell $c in cal(H)^((t))$, define its count and coordinate sum:
  $
    m^((t)) (c) &= |{i in [n] mid(|) q_h (x_i^((t))) = c}|, \
    s^((t)) (c) &= sum_(i : q_h (x_i^((t))) = c) x_i^((t)).
  $
  The *Mean Shift++* cell-level shift target and point update are @jang2021meanshiftpp:
  $
    mu^((t)) (c) &= (sum_(c' in cal(N)_1 (c) inter cal(H)^((t))) s^((t)) (c')) / (sum_(c' in cal(N)_1 (c) inter cal(H)^((t))) m^((t)) (c')), \
    x_i^((t+1)) &= mu^((t)) (q_h (x_i^((t)))).
  $
  Equivalently, defining the symmetric binary cell-adjacency matrix $W_(i j)^((t)) = bb(1){norm(q_h (x_i^((t))) - q_h (x_j^((t))))_infinity <= 1}$ and $bold(D)^((t)) = "diag"(sum_j W_(i j)^((t)))$, the dataset obeys the Markov recursion $bold(X)^((t+1)) = (bold(D)^((t)))^(-1) bold(W)^((t)) bold(X)^((t))$.
]

Because $x_i^((t+1))$ depends on $x_i^((t))$ solely through its discrete cell index $q_h (x_i^((t)))$, every point occupying the same cell $c$ jumps simultaneously to the identical coordinate $mu^((t)) (c)$, enabling rapid cell coalescence across iterations.

#algbox(
  title: "Mean Shift++ (Grid-Accelerated Mode-Seeking)",
  inputs: [Dataset $X^((0)) = {x_1^((0)), dots, x_n^((0))} subset bb(R)^d$, bandwidth $h > 0$, tolerance $epsilon > 0$, max steps $T_"max"$],
  outputs: [Cluster labels $c_1, dots, c_n$ and discrete modes ${m_1, dots, m_K}$],
)[
  1. Set $t <- 0$ and initialize $x_i^((0)) <- x_i$ for $i = 1, dots, n$.
  2. *Repeat* for $t = 0, 1, dots, T_"max" - 1$:
     (a) Initialize an empty hash table $cal(H)^((t))$ mapping $c in bb(Z)^d |-> (m(c), s(c))$.
     (b) *Hash Accumulation*: For each $i = 1, dots, n$, compute $c_i <- floor(x_i^((t)) / h)$ and update $m(c_i) <- m(c_i) + 1$, $s(c_i) <- s(c_i) + x_i^((t))$.
     (c) *Cell-Level Shift*: For each occupied cell $c in cal(H)^((t))$, query all $3^d$ offsets $c' in c + {-1, 0, 1}^d$ present in $cal(H)^((t))$ and compute
         $mu^((t)) (c) <- (sum_(c') s(c')) / (sum_(c') m(c'))$.
     (d) *Point Update*: For each $i = 1, dots, n$, set $x_i^((t+1)) <- mu^((t)) (c_i)$.
     (e) *If* $max_(1 <= i <= n) norm(x_i^((t+1)) - x_i^((t)))_2 <= epsilon$, *break*.
  3. *Mode Partitioning*: Group points sharing identical (or $epsilon$-close) final coordinates $x_i^((t+1))$ into clusters $1, dots, K$ and return labels and mode centroids.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five formal results: an exact geometric sandwich for the grid stencil and cell coalescence, an exact if-and-only-if fixed-point characterization and separation invariant, density-gradient consistency, constructive proof of translation and rotation non-invariance, and sharp complexity bounds quantifying the $3^d$ dimensionality barrier.

#subsec("A", "Geometric Stencil Sandwich and Cell Coalescence")

#lemma(title: "Geometric Inclusion Sandwich and Monotone Cell Coalescence", origin: "Original Proof")[
  For every $x in bb(R)^d$ and $h > 0$, let $z_h (x) = (q_h (x) + 1 / 2 bold(1)_d) h$ be the geometric center of cell $Q_h (q_h (x))$. Then:
  1. *(Exact stencil sandwich)* $S_h (x) = z_h (x) + [-3 h / 2, 3 h / 2)^d$ has exact Lebesgue volume $(3 h)^d$ and satisfies
     $
       B_2 (x, h) &subset.eq B_infinity (x, h) subset.eq S_h (x) \
                  &subset B_infinity (x, 2 h) subset.eq B_2 (x, 2 sqrt(d) h),
     $
     where $B_p (x, r) = {y in bb(R)^d mid(|) norm(y - x)_p < r}$.
  2. *(Cell-coalescence monotonicity)* Let $U_t = |{x_1^((t)), dots, x_n^((t))}|$ be the number of distinct point coordinates at step $t$ and $M_t = |cal(H)^((t))|$ the number of occupied cells. Then for every $t >= 0$,
     $
       1 <= M_(t+1) <= U_(t+1) <= M_t <= U_t <= n.
     $
]

#proof[
  *Part 1.* For any $y in bb(R)^d$, write $x_r = (k_r + u_r) h$ and $y_r = (l_r + v_r) h$ with $k_r, l_r in bb(Z)$ and $u_r, v_r in [0, 1)$.
  - By definition, $y in S_h (x)$ iff $l_r in {k_r - 1, k_r, k_r + 1}$ for every $r$, which holds iff $y_r in [(k_r - 1) h, (k_r + 2) h) = [z_(h, r) (x) - 3 h / 2, z_(h, r) (x) + 3 h / 2)$. Hence $S_h (x) = z_h (x) + [-3 h / 2, 3 h / 2)^d$ with volume $(3 h)^d$.
  - If $y in B_infinity (x, h)$, then $|(k_r - l_r) + (u_r - v_r)| = |x_r - y_r| / h < 1$ for every $r$. Since $|u_r - v_r| < 1$, any $|k_r - l_r| >= 2$ would imply $|(k_r - l_r) + (u_r - v_r)| > 2 - 1 = 1$, a contradiction. Thus $|k_r - l_r| <= 1$, so $B_infinity (x, h) subset.eq S_h (x)$.
  - Conversely, if $y in S_h (x)$, then $|k_r - l_r| <= 1$ and $|u_r - v_r| < 1$, giving $|x_r - y_r| <= (|k_r - l_r| + |u_r - v_r|) h < 2 h$, so $S_h (x) subset B_infinity (x, 2 h)$.
  - Combining with the universal norm inequalities $norm(dot)_infinity <= norm(dot)_2 <= sqrt(d) norm(dot)_infinity$ yields $B_2 (x, h) subset.eq B_infinity (x, h)$ and $B_infinity (x, 2 h) subset.eq B_2 (x, 2 sqrt(d) h)$.

  *Part 2.* Clearly $M_t <= U_t$ because identical coordinates hash to identical cells. By Definition 1, $x_i^((t+1)) = mu^((t)) (q_h (x_i^((t))))$ takes at most one value per occupied cell $c in cal(H)^((t))$, so $U_(t+1) <= M_t$. Applying $M_(t+1) <= U_(t+1)$ closes the chain.
]

#subsec("B", "Exact Fixed-Point Characterization and Separation Invariant")

Unlike Gaussian Blurring Mean Shift (Alg. 54), which globally collapses to a single point as $t -> infinity$, Mean Shift++ naturally preserves multiple separated modes as exact fixed points.

#theorem(title: "Fixed-Point Characterization and Multi-Mode Separation", origin: "Original Proof")[
  Let $G^((t)) = (cal(H)^((t)), E^((t)))$ be the undirected grid-adjacency graph with edges connecting occupied cells $c, c' in cal(H)^((t))$ iff $norm(c - c')_infinity <= 1$.
  1. *(Bounding-box contraction)* Along every coordinate $r in {1, dots, d}$, the bounding interval $[min_i x_(i, r)^((t)), max_i x_(i, r)^((t))]$ is nested across iterations, and $"conv"(X^((t+1))) subset.eq "conv"(X^((t)))$.
  2. *(Exact fixed-point theorem)* A configuration $X^((t))$ is a fixed point of Mean Shift++ ($X^((t+1)) = X^((t))$) *if and only if* every connected component $V$ of $G^((t))$ has collapsed to a single point $z_V in bb(R)^d$. Consequently, any two distinct modes $z_a != z_b$ at a fixed point satisfy the strict separation $norm(z_a - z_b)_infinity > h$.
  3. *(Permanent cluster decoupling)* If at some step $t_0$ the dataset partitions into groups $C_1, dots, C_K$ whose coordinate-aligned bounding boxes satisfy $min_(a != b) "dist"_infinity ("box"(C_a^((t_0))), "box"(C_b^((t_0)))) >= 2 h$, then $C_1, dots, C_K$ never interact for any $t >= t_0$.
]

#proof[
  *Part 1.* Every update $x_i^((t+1)) = sum_(j=1)^n P_(i j)^((t)) x_j^((t))$ is a convex combination of $X^((t))$, immediately implying $"conv"(X^((t+1))) subset.eq "conv"(X^((t)))$ and coordinate-wise bounding-interval nesting.

  *Part 2.* ($arrow.l.double$) Suppose every connected component $V$ of $G^((t))$ consists of points coinciding at a single coordinate $z_V$. Then all points of $V$ occupy a single cell $c_V = q_h (z_V)$, and no occupied cell from any other component lies in $cal(N)_1 (c_V)$ (otherwise the two components would be connected by an edge in $G^((t))$). Thus $mu^((t)) (c_V) = (m(c_V) z_V) / m(c_V) = z_V$, so $x_i^((t+1)) = x_i^((t))$ for all $i$.

  ($arrow.r.double$) Suppose $X^((t+1)) = X^((t))$. Fix any connected component $V$ of $G^((t))$ and any coordinate $r in {1, dots, d}$. Let $u = max_(i : q_h (x_i^((t))) in V) x_(i, r)^((t))$, and let $i^*$ be a point attaining $x_(i^*, r)^((t)) = u$ with cell $c^* = q_h (x_(i^*)^((t))) in V$. By the fixed-point condition,
  $
    u = x_(i^*, r)^((t+1)) = sum_(j : q_h (x_j^((t))) in cal(N)_1 (c^*)) P_(i^* j)^((t)) x_(j, r)^((t)).
  $
  Every $j$ in this sum has $q_h (x_j^((t))) in V$, hence $x_(j, r)^((t)) <= u$, while all weights $P_(i^* j)^((t))$ are strictly positive and sum to $1$. Equality therefore forces $x_(j, r)^((t)) = u$ for every point $j$ in every cell adjacent to $c^*$ in $G^((t))$. Inducting along shortest paths in the connected component $V$ proves $x_(j, r)^((t)) = u$ for every point in $V$. Because $r in {1, dots, d}$ was arbitrary, all points in $V$ coincide at one vector $z_V$.

  For distinct components $V_a != V_b$, their single occupied cells satisfy $norm(q_h (z_a) - q_h (z_b))_infinity >= 2$. Pick a coordinate $r$ with $|k_r - l_r| >= 2$ where $z_(a, r) = (k_r + u_r) h$, $z_(b, r) = (l_r + v_r) h$, and $u_r, v_r in [0, 1)$. Since $|u_r - v_r| < 1$, the reverse triangle inequality gives $|z_(a, r) - z_(b, r)| >= (|k_r - l_r| - |u_r - v_r|) h > (2 - 1) h = h$, hence $norm(z_a - z_b)_infinity > h$.

  *Part 3.* By induction on $t >= t_0$, assume $"box"(C_a^((t))) subset.eq "box"(C_a^((t_0)))$ for all $a$. For $i in C_a$, $j in C_b$ ($a != b$), $norm(x_i^((t)) - x_j^((t)))_infinity >= 2 h$, so by Lemma 1 ($S_h (x) subset B_infinity (x, 2 h)$) we have $x_j^((t)) not in S_h (x_i^((t)))$, i.e., $norm(q_h (x_i^((t))) - q_h (x_j^((t))))_infinity >= 2$. Thus $W_(i j)^((t)) = 0$ across distinct groups, every $C_a$ updates strictly inside its own bounding box, and the invariant holds for $t + 1$.
]

#subsec("C", "Statistical Density-Gradient Consistency")

#theorem(title: "Population Gradient-Ascent Consistency of the Grid Shift", origin: "Classical — Cited")[
  Let $X_1, dots, X_n$ be i.i.d. draws from a density $p in C^3 (bb(R)^d)$ bounded away from zero on a compact domain $Omega subset bb(R)^d$ @jang2021meanshiftpp @ariascastro2016meanshift. For any $x in Omega$ with cell center $z_h (x) = (q_h (x) + 1 / 2 bold(1)_d) h$ whose stencil $S_h (x) subset Omega$, the population grid-stencil expectation satisfies
  $
    bb(E)[Y mid(|) Y in S_h (x)] - z_h (x)
    = (3 h^2) / 4 (nabla p(z_h (x))) / (p(z_h (x))) + cal(O)(h^4).
  $
  Furthermore, if $h = h_n -> 0$ with $n h_n^(d + 2) / log n -> infinity$, then the empirical cell-shift vector $mu^((0)) (q_h (x)) - z_h (x)$ converges uniformly in probability on $Omega$ to $(3 h^2 / 4) nabla log p(z_h (x))$ at rate $cal(O)_p (h^2 + sqrt((log n) / (n h^(d - 2))))$ @jang2021meanshiftpp.
]

#proof[
  By Lemma 1, $S_h (x) = z_0 + [-a, a)^d$ with $z_0 = z_h (x)$ and half-width $a = 3 h / 2$. Changing variables $y = z_0 + u$ and expanding $p(z_0 + u) = p(z_0) + u^top nabla p(z_0) + 1 / 2 u^top bold(H)_p (z_0) u + R_3 (u)$ with $|R_3 (u)| <= C norm(u)_2^3$ on $Omega$:
  - *Denominator*: By symmetry of $[-a, a)^d$, odd monomials integrate to zero, so
    $
      integral_([-a, a]^d) p(z_0 + u) dif u
      = (2 a)^d (p(z_0) + cal(O)(a^2)).
    $
  - *Numerator*: For each coordinate $r in {1, dots, d}$, odd symmetry cancels $integral u_r p(z_0) dif u = 0$, $integral u_r u_s u_t dif u = 0$, and cross-terms $integral u_r u_s dif u = 0$ ($r != s$). The sole leading term is
    $
      integral_([-a, a]^d) u_r^2 partial_r p(z_0) dif u
      &= (2 a)^(d - 1) ((2 a^3) / 3) partial_r p(z_0) \
      &= (2 a)^d (a^2 / 3) partial_r p(z_0).
    $
    The third-order remainder contributes $cal(O)((2 a)^d a^4)$ because $integral norm(u)_2^4 dif u = cal(O)((2 a)^d a^4)$.

  Substituting $a^2 / 3 = (3 h / 2)^2 / 3 = 3 h^2 / 4$ and dividing numerator by denominator yields
  $
    bb(E)[Y - z_0 mid(|) Y in S_h (x)]
    &= ((3 h^2 / 4) nabla p(z_0) + cal(O)(h^4)) / (p(z_0) + cal(O)(h^2)) \
    &= (3 h^2) / 4 (nabla p(z_0)) / (p(z_0)) + cal(O)(h^4).
  $
  In the finite-sample regime, applying Bernstein / Hoeffding concentration over the polynomial number of grid cells intersecting $Omega$ controls the empirical mean deviation around $bb(E)[Y mid(|) Y in S_h (x)]$ at scale $cal(O)_p (h sqrt((log n) / (n h^d)))$ @jang2021meanshiftpp, establishing uniform consistency of the gradient-ascent direction.
]

#subsec("D", "Grid-Phase Aliasing and the Exponential Dimensionality Barrier")

#proposition(title: "Translation and Rotation Non-Invariance of Mean Shift++", origin: "Original Counterexample")[
  Unlike Mean Shift and Blurring Mean Shift (which are strictly invariant under rigid Euclidean motions), Mean Shift++ is invariant under sample permutations, joint isotropic scaling $(x, h) |-> (alpha x, alpha h)$ ($alpha > 0$), and integer grid shifts $x |-> x + h z$ ($z in bb(Z)^d$), but is *neither translation-invariant nor rotation-invariant*:
  1. *(Translation aliasing)* For two points in $bb(R)$ separated by $1.5 h$, ${0.2 h, 1.7 h}$ merges into $1$ cluster in one step, whereas its $+0.5 h$ translate ${0.7 h, 2.2 h}$ remains permanently split into $2$ clusters.
  2. *(Rotation aliasing)* For two points in $bb(R)^2$ separated by Euclidean distance $2.2 h$, the axis-aligned pair ${(0.2 h, 0.2 h), (2.4 h, 0.2 h)}$ forms $2$ fixed clusters, whereas rotating the second point by $45 degree$ around the first merges both points into $1$ cluster in a single iteration.
]

#proof[
  *Invariances.* Permutation invariance holds because cell counts $m(c)$ and sums $s(c)$ are symmetric functions of sample indices. For $alpha > 0$, $floor((alpha x) / (alpha h)) = floor(x / h)$, so cell indices are identical and $mu^((t)) (c)$ scales by $alpha$. For $z in bb(Z)^d$, $floor((x + h z) / h) = floor(x / h) + z$, shifting $cal(H)^((t))$ rigidly by $z$ and translating $mu^((t))$ by $h z$.

  *Part 1.* For $A = {0.2 h, 1.7 h}$, $q_h (0.2 h) = 0$ and $q_h (1.7 h) = 1$, so $|0 - 1| = 1 <= 1$ and both points jump to $0.95 h$ in step $1$. For $A + 0.5 h = {0.7 h, 2.2 h}$, $q_h (0.7 h) = 0$ and $q_h (2.2 h) = 2$, giving $|0 - 2| = 2 > 1$; by Theorem 1(2) the two points form a 2-cluster fixed point.

  *Part 2.* Let $x_1 = (0.2 h, 0.2 h)$ and $x_2 = (2.4 h, 0.2 h)$, so $norm(x_2 - x_1)_2 = 2.2 h$. Here $q_h (x_1) = (0, 0)$ and $q_h (x_2) = (2, 0)$, so $norm(q_h (x_1) - q_h (x_2))_infinity = 2 > 1$, yielding 2 fixed clusters. Rotating $x_2 - x_1$ by $45 degree$ yields $x_2' = x_1 + (1.1 sqrt(2) h, 1.1 sqrt(2) h)$. Since $2.25 < 2.42 = (1.1 sqrt(2))^2 < 2.56$, we have $1.5 < 1.1 sqrt(2) < 1.6$, hence $1.7 < 0.2 + 1.1 sqrt(2) < 1.8 < 2$. Therefore $q_h (x_2') = (1, 1)$ satisfies $norm(q_h (x_1) - q_h (x_2'))_infinity = 1$, collapsing ${x_1, x_2'}$ to their midpoint in one step.
]

#theorem(title: "Linear-in-n Complexity and Exponential-in-d Stencil Barrier", origin: "Original Proof")[
  Using a hash table with $cal(O)(1)$ expected lookup and insertion, iteration $t$ of Mean Shift++ executes in expected time $cal(O)((n + M_t dot 3^d) d)$ and total runtime over $T$ iterations is
  $
    cal(O)(T d (n + bar(M) dot 3^d)) <= cal(O)(T n d dot 3^d),
  $
  where $bar(M) = T^(-1) sum_(t=0)^(T-1) M_t <= M_0 <= n$, with space $cal(O)((n + M_0) d)$. Moreover:
  1. *(Asymptotic crossover)* Whenever $3^d > n$ (i.e., $d > log_3 n$), per-cell Moore-stencil probing is asymptotically slower than brute-force $cal(O)(n^2 d)$ Blurring Mean Shift whenever $M_t = Theta(n)$.
  2. *(High-dimensional singleton collapse)* For $n$ i.i.d. draws from $cal(N)(bold(0), sigma^2 bold(I)_d)$ with fixed ratio $h / sigma = c$, the probability that two distinct draws fall in adjacent cells satisfies
     $
       "Pr"(norm(q_h (X_1) - q_h (X_2))_infinity <= 1) <= ("erf"(c))^d < 1,
     $
     which decays exponentially to $0$ as $d -> infinity$; thus for $d > 2 log n / log(1 / "erf"(c))$, every point becomes an isolated singleton fixed point with high probability.
]

#proof[
  *Complexity.* In Step 2(b), quantizing and inserting $n$ $d$-dimensional vectors into $cal(H)^((t))$ costs $cal(O)(n d)$ expected time. In Step 2(c), iterating over $M_t$ occupied cells and probing $3^d$ neighbor keys of length $d$ takes $cal(O)(M_t d dot 3^d)$ expected time; each matched cell contributes one $d$-dimensional vector addition. Step 2(d) updates $n$ points in $cal(O)(n d)$. Summing over $T$ iterations yields $cal(O)(T d (n + bar(M) 3^d))$, and Lemma 1(2) guarantees $M_t <= M_0 <= n$, bounding memory by $cal(O)((n + M_0) d)$. When $M_t = Theta(n)$ and $3^d > n$, $M_t d dot 3^d = omega(n^2 d)$, proving (1).

  *Singleton collapse.* By Lemma 1(1), $norm(q_h (X_1) - q_h (X_2))_infinity <= 1$ implies $norm(X_1 - X_2)_infinity < 2 h$. Because the coordinates of $Z = X_1 - X_2 tilde cal(N)(bold(0), 2 sigma^2 bold(I)_d)$ are independent with standard deviation $sqrt(2) sigma$,
  $
    "Pr"(norm(Z)_infinity < 2 h)
    &= product_(r=1)^d "Pr"(|Z_r| < 2 h) \
    &= ("erf"((2 h) / (sqrt(2) dot sqrt(2) sigma)))^d
     = ("erf"(c))^d.
  $
  By a union bound over all $binom(n, 2) < n^2 / 2$ unordered pairs, the probability that *any* two points are stencil-adjacent is at most $(n^2 / 2) ("erf"(c))^d$, which tends to $0$ whenever $d > 2 log n / log(1 / "erf"(c))$. When no two points are adjacent, every occupied cell is an isolated singleton component and Mean Shift++ halts immediately at $t = 0$ with $K = n$ trivial clusters (Theorem 1(2)).
]

#sec("III", "Complexity, Geometry and Operational Profile")

In low dimensions ($d in {2, 3, 4, 5}$), Mean Shift++ delivers dramatic speedups—often $100 times$ to $1000 times$ faster than exact Mean Shift on datasets with $n >= 10^5$—because $3^d$ is a modest constant ($9, 27, 81, 243$) and cell coalescence (Lemma 1(2)) rapidly shrinks $M_t << n$ after the very first iteration @jang2021meanshiftpp. Space usage stays strictly linear $cal(O)(n d)$, avoiding the $cal(O)(n^2)$ affinity matrix that blocks dense GBMS at large scale. However, the exponential $3^d$ table-lookup factor and hash-table overhead mean that for small $n$ ($n <= 10^3$), Mean Shift++ can run slightly slower than BLAS-vectorized GBMS, while for $d >= 8$ empty-cell sparsity triggers either prohibitive stencil probes or trivial singleton/single-cluster degeneration (Theorem 3).

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time complexity], [$cal(O)(T d (n + bar(M) dot 3^d))$ expected; linear in $n$ for fixed low $d$],
    [Space complexity], [$cal(O)((n + M_0) d)$ hash table + coordinate array (no $n times n$ matrix)],
    [Hyperparameters], [Single grid bandwidth $h > 0$ and stopping tolerance $epsilon$],
    [Convergence], [Exact multi-mode fixed points separated by $> h$ in $ell_infinity$ (Theorem 1)],
    [Symmetries], [Permutation, isotropic scale, and integer-grid shift invariant; *not* $"SE"(d)$-invariant],
    [Cluster geometry], [Non-convex low-D manifolds; axis-aligned box stencil induces mild anisotropy],
    [Practical sweet spot], [Large $n$ ($10^4$--$10^7$) in low dimensions ($d <= 5$)],
  ),
  caption: [Operational and computational profile of Mean Shift++.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *High-Resolution Image Segmentation and Superpixels*: In 3D $(L^*, a^*, b^*)$ color space or 5D $(L^*, a^*, b^*, x, y)$ joint spatial-color space, $3^3 = 27$ and $3^5 = 243$ are tiny relative to $n = 10^6$ pixels. Mean Shift++ produces segmentation quality matching classical Mean Shift @comaniciu2002mean at interactive frame rates competitive with SLIC superpixels @achanta2012slic @jang2021meanshiftpp.
2. *Real-Time Video Object Tracking*: By tracking color-spatial density modes across consecutive video frames with $cal(O)(n)$ per-frame updates, Mean Shift++ tracks fast non-rigid targets without quadratic kernel bottlenecks @jang2021meanshiftpp.
3. *Large-Scale 2D/3D Geospatial and LiDAR Point Clouds*: Urban GPS pings ($d = 2$) and autonomous-driving LiDAR frames ($d = 3$, $n tilde 10^5$ points per sweep) map naturally onto sparse Cartesian grids ($3^2 = 9$, $3^3 = 27$), enabling real-time obstacle grouping and hotspot discovery with automatic cluster-count selection.
4. *Interactive Exploration of Low-Dimensional Embeddings*: After compressing high-dimensional data into 2D or 3D visualizations via UMAP or t-SNE, millions of projected points can be mode-clustered interactively as users scrub the single bandwidth parameter $h$.
5. *RGB-D Robotics and Volumetric Medical Voxel Grouping*: Intensity-plus-depth or 3D MRI intensity manifolds lie in $d in {3, 4}$, where grid hashing exploits spatial regularity and cell coalescence collapses homogeneous anatomical or scene surfaces in a handful of passes.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Mean Shift++ exhibits four distinct failure regimes: (i) *exponential $3^d$ runtime and singleton collapse* in moderate-to-high dimensions ($d >= 8$, Theorem 3); (ii) *grid-phase aliasing* under arbitrary translations and rotations (Proposition 1), which can split borderline modes or merge near-touching shapes; (iii) *anisotropic feature distortion*, where heterogeneous coordinate variances cause a single scalar $h$ to over-merge low-variance axes while fracturing high-variance axes into singletons; and (iv) *box-stencil artifacts* on diagonally curved manifolds, where the $2 sqrt(d) h$ diagonal reach exceeds the $h$ axis-aligned reach (Lemma 1).

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Image segmentation & video tracking], [Excellent: $d in {3, 5}$, huge $n$, $cal(O)(n)$ speed], [Normalize spatial vs. color scales; cf. SLIC],
    [2D/3D LiDAR & geospatial hotspots], [Excellent: sparse hash grid handles $n > 10^6$], [Jitter/ensemble grid phase if boundary-sensitive],
    [2D/3D UMAP / t-SNE mode discovery], [High: interactive bandwidth sweeps without $K$], [Sensitive to $h$ across multi-scale densities],
    [Heterogeneous tabular data ($d >= 6$)], [Poor: scalar $h$ + $3^d$ sparsity yield collapse/N/A], [Standardize or use Blurring MS / HDBSCAN\*],
    [High-D text / vision embeddings], [Infeasible: $3^d$ blow-up & singleton collapse], [Reduce to $d <= 5$ first or use spherical $k$-means],
  ),
  caption: [Domain suitability matrix for Mean Shift++.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Mean Shift++ sits at the exact intersection of grid-based density clustering (STING @wang1997sting, CLIQUE @agrawal1998clique, WaveCluster @sheikholeslami1998wavecluster) and kernel mode-seeking (Mean Shift @comaniciu2002mean, Blurring Mean Shift @cheng1995mean @carreira2006fast, Quick Shift @vedaldi2008quickshift, Medoid Shift @sheikh2007medoidshift, and DENCLUE 2.0 @hinneburg2007denclue2), echoing the algorithmic philosophy of DBSCAN++ @jang2019dbscanpp in trading exact metric queries for provably consistent fast approximations.

#proposition(title: "Structural Reduction to Quantized Flat-Kernel Blurring Mean Shift", origin: "Original Derivation")[
  Let $bold(P)_"flat"^((t)) (h)$ denote the row-stochastic transition matrix of flat-kernel Blurring Mean Shift with $ell_infinity$ bandwidth $h$, i.e., $W_(i j)^("flat") = bb(1){norm(x_i^((t)) - x_j^((t)))_infinity < h}$. Then the Mean Shift++ transition matrix $bold(P)_"MS++" ^((t)) (h)$ coincides identically with $bold(P)_"flat"^((t)) (3 h / 2)$ applied to the cell-quantized coordinates $z_h (x_i^((t))) = (q_h (x_i^((t))) + 1 / 2 bold(1)_d) h$ for affinity computation while averaging the unquantized coordinates $X^((t))$:
  $
    W_(i j)^("MS++") (h) = bb(1){norm(z_h (x_i^((t))) - z_h (x_j^((t))))_infinity < (3 h) / 2}.
  $
]

#proof[
  For any $x, y in bb(R)^d$, $z_h (x) - z_h (y) = (q_h (x) - q_h (y)) h$. Because $q_h (x) - q_h (y) in bb(Z)^d$ has integer entries, $norm(q_h (x) - q_h (y))_infinity <= 1$ holds if and only if $norm(q_h (x) - q_h (y))_infinity < 3 / 2$, i.e., iff $norm(z_h (x) - z_h (y))_infinity < 3 h / 2$. Thus the binary support of the Mean Shift++ affinity matrix is identically the flat $ell_infinity$-kernel of radius $3 h / 2$ on the quantized lattice centers $z_h (x_i^((t)))$, while the numerator sums the continuous vectors $x_j^((t))$.
]

Proposition 2 illuminates both the speed and the empirical trade-offs observed on the ClusterOverview 4-dataset benchmark evaluated via Silhouette @rousseeuw1987silhouettes and ARI @hubert1985comparing. On Iris (`Sil = 0.586`, `ARI = 0.567`), Mean Shift++ reproduces the *exact* Silhouette and ARI scores of Blurring Mean Shift (Alg. 54), confirming empirically that in low-dimensional, well-structured regimes the quantized grid stencil matches continuous Gaussian mode-seeking. On Multi-Shape (`Sil = 0.567`, `ARI = N/A`), Mean Shift++ outperforms OPTICS-OF (`0.446`) by collapsing curved structures into clean discrete modes, but trails Blurring Mean Shift (`0.787`) because grid-phase aliasing (Proposition 1) and the box-shaped $ell_infinity$ stencil (Lemma 1) blur diagonal curvature more coarsely than an isotropic Gaussian kernel. On Customers (`Sil = N/A`) and High-D (`Sil = N/A`, `ARI = N/A`), Mean Shift++ fails to yield a valid multi-cluster partition: in Customers, heterogeneous tabular scales and moderate dimensionality cause scalar-$h$ grid partitioning to collapse into either a single connected cell cluster or fragmented singletons, while in High-D, Theorem 3's exponential $3^d$ sparsity barrier guarantees singleton isolation @beyer1999nearest. Finally, its benchmark mean runtime of `4.897 ms` slightly exceeds Blurring Mean Shift's `4.568 ms` solely because the benchmark datasets have small $n$ where Python/hash-table overhead outweighs $cal(O)(n^2)$ vectorized BLAS calls; as Theorem 3 guarantees, for $n >> 3^d$ in $d <= 5$, Mean Shift++ achieves orders-of-magnitude speedups over every quadratic mode-seeker.
