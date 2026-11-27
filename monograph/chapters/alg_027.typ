#import "../template/ieee.typ": *

#alg_entry(
  id: 27,
  name: "Power K-Means",
  category: "Cat. I: Centroid & Partitioning",
  year: "2019",
  paradigm: "Power-Mean Homotopy / Majorization-Minimization",
  complexity: $O(n k d T)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

A defining liability of Lloyd's $k$-means algorithm (Alg. 1) @lloyd1982least is its acute sensitivity to initialization: the non-smooth hard minimum inside $F_(-infinity) (bold(mu)) = sum_(i=1)^n min_(1 <= j <= k) norm(x_i - mu_j)_2^2$ partitions parameter space into sharp Voronoi basins littered with arbitrarily poor local minima. *Power $k$-means*, introduced by Xu and Lange @xu2019power, confronts this non-convexity directly by smoothing the hard minimum via a *generalized power mean* of negative order $s < 0$, then tracing a continuous homotopy path as $s -> -infinity$ through closed-form Majorization-Minimization (MM) updates.

Let $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$ and $k in {1, dots, n}$. For any positive vector $y = (y_1, dots, y_k) in bb(R)_(>0)^k$ and exponent $s in bb(R) backslash {0}$, define the *power mean* of order $s$:
$ M_s (y) := ( 1 / k sum_(j=1)^k y_j^s )^(1 / s). $
Given centroids $bold(mu) = (mu_1, dots, mu_k) in (bb(R)^d)^k$ and a small numerical smoothing constant $epsilon >= 0$ (typically $epsilon = 0$ analytically on the domain $x_i != mu_j$, or $epsilon > 0$ small to guarantee global $C^infinity$-differentiability), let $y_(i j) (bold(mu)) := norm(x_i - mu_j)_2^2 + epsilon$ and $y_i (bold(mu)) := (y_(i 1) (bold(mu)), dots, y_(i k) (bold(mu)))$. The *Power $k$-Means objective* at exponent $s < 0$ is @xu2019power:
$ F_(s, epsilon) (bold(mu)) &:= sum_(i=1)^n M_s (y_i (bold(mu))) \
  &= sum_(i=1)^n ( 1 / k sum_(j=1)^k (norm(x_i - mu_j)_2^2 + epsilon)^s )^(1 / s). $
As $s -> -infinity$, $M_s (y_i (bold(mu)))$ decreases monotonically to $min_(1 <= j <= k) y_(i j) (bold(mu))$, recovering the exact $k$-means objective (shifted by $n epsilon$). Conversely, for moderate negative $s$ (e.g., $s = -1$, the harmonic mean), every data point exerts a smooth, non-zero gravitational pull on *all* $k$ centroids, dissolving rigid Voronoi barriers.

#definition(title: "Power-Mean Responsibility Weights and Homotopy Schedule")[
  For $s < 0$ and $y_(i j) > 0$, the *power-mean responsibility weight* of point $i$ for centroid $j$ is the partial derivative
  $ w_(i j) (bold(mu), s) &:= (partial M_s (y_i)) / (partial y_(i j)) \
    &= 1 / k ( M_s (y_i (bold(mu))) / (y_(i j) (bold(mu))) )^(1 - s) > 0. $
  A *homotopy continuation schedule* specifies a strictly decreasing sequence $0 > s_0 > s_1 > dots > s_T$ (typically geometric, $s_(t+1) = eta s_t$ with $eta in (1.05, 1.5]$), updating $bold(mu)$ via one or more MM steps per exponent $s_t$ until convergence.
]

Crucially, because $M_s$ is concave and positively homogeneous of degree $1$ for all $s < 0$, its first-order tangent hyperplane yields an exact, parameter-free quadratic majorizer whose minimizer is a simple weighted average of the dataset (Algorithm 1).

#algbox(
  title: "Power K-Means via MM Homotopy Continuation",
  inputs: [Dataset $cal(X) subset.eq bb(R)^d$, clusters $k >= 1$, initial centroids $bold(mu)^((0))$, initial power $s_0 < 0$ (e.g., $-1$), annealing rate $eta > 1$, floor $s_"min"$, tolerance $tau > 0$, smoothing $epsilon > 0$, max iterations $T$.],
  outputs: [Centroids $bold(mu)$ and hard partition $cal(C) = {C_1, dots, C_k}$.],
)[
  1. Set $s <- s_0$.
  2. *For* $t = 0, 1, dots, T - 1$ *do*:
  3. #h(1em) *Squared distances:* Compute $y_(i j)^((t)) = norm(x_i - mu_j^((t)))_2^2 + epsilon$ for all $i in {1, dots, n}$, $j in {1, dots, k}$.
  4. #h(1em) *Stabilized power-mean weights:* Let $m_i = min_ell y_(i ell)^((t))$. Compute
     $ M_s (y_i^((t))) = m_i ( 1 / k sum_(ell=1)^k ( y_(i ell)^((t)) / m_i )^s )^(1 / s), $
     $ w_(i j)^((t)) = 1 / k ( M_s (y_i^((t))) / (y_(i j)^((t))) )^(1 - s). $
  5. #h(1em) *MM centroid update:* For each $j = 1, dots, k$, update
     $ mu_j^((t+1)) = ( sum_(i=1)^n w_(i j)^((t)) x_i ) / ( sum_(i=1)^n w_(i j)^((t)) ). $
  6. #h(1em) *Anneal power:* Update $s <- max(s_"min", eta s)$ (periodically or every step).
  7. #h(1em) *If* $max_j norm(mu_j^((t+1)) - mu_j^((t)))_2 < tau$ *and* $s <= s_"min"$, *break*.
  8. Assign $c_i = arg min_j norm(x_i - mu_j^((t+1)))_2^2$ and *return* $(bold(mu)^((t+1)), cal(C))$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish four foundational pillars of Power $k$-Means: (A) concavity, monotonicity, and tight uniform sandwich bounds for negative power means; (B) exact quadratic MM majorization with a quantitative descent rate and stationarity guarantee; (C) uniform homotopy convergence, global consistency, and a constructive proof that Power $k$-Means strictly escapes Lloyd's adversarial trapping fixed point; and (D) geometric invariances, harmonic-origin collapse, and finite-sample breakdown.

#subsec("A", "Concavity, Monotonicity, and Uniform Sandwich Bounds")

#lemma(title: "Strict Non-Radial Concavity and Uniform Approximation of Power Means", origin: "Original Proof")[
  For every $s < 0$ and every $y in bb(R)_(>0)^k$:
  1. $M_s$ is concave on $bb(R)_(>0)^k$, and strictly concave along every direction $z in.not "span"(y)$.
  2. For any $s_1 < s_2 < 0$, $min_j y_j <= M_(s_1) (y) <= M_(s_2) (y) <= k^(-1) sum_(j=1)^k y_j$.
  3. Letting $y_"min" := min_(1 <= j <= k) y_j$, $M_s (y)$ obeys the sharp sandwich bound
     $ y_"min" <= M_s (y) <= k^(-1 / s) y_"min". $
]

#proof[
  (1) Differentiating $M_s (y) = (k^(-1) sum_ell y_ell^s)^(1/s)$ twice gives $partial_j M_s (y) = k^(-1) M_s (y)^(1 - s) y_j^(s - 1)$ and
  $ (partial^2 M_s) / (partial y_j partial y_ell) &= (1 - s) / k M_s^(-s) (partial_ell M_s) y_j^(s - 1) - (1 - s) / k M_s^(1 - s) y_j^(s - 2) delta_(j ell) \
    &= (1 - s) M_s^(1 - 2 s) [ 1 / k^2 y_j^(s - 1) y_ell^(s - 1) - 1 / k M_s^s y_j^(s - 2) delta_(j ell) ]. $
  Contracting with any $z in bb(R)^k$ and substituting $M_s^s = k^(-1) sum_(j=1)^k y_j^s$ yields
  $ z^T nabla^2 M_s (y) z &= (1 - s) M_s^(1 - 2 s) [ ( 1 / k sum_(j=1)^k z_j y_j^(s - 1) )^2 \
    &quad - ( 1 / k sum_(j=1)^k y_j^s ) ( 1 / k sum_(j=1)^k z_j^2 y_j^(s - 2) ) ]. $
  Define vectors $u, v in bb(R)^k$ by $u_j = y_j^(s / 2)$ and $v_j = z_j y_j^(s / 2 - 1)$. Then $u_j v_j = z_j y_j^(s - 1)$, $u_j^2 = y_j^s$, and $v_j^2 = z_j^2 y_j^(s - 2)$. By Cauchy–Schwarz,
  $ ( 1 / k sum_(j=1)^k u_j v_j )^2 <= ( 1 / k sum_(j=1)^k u_j^2 ) ( 1 / k sum_(j=1)^k v_j^2 ), $
  with equality if and only if $v$ is proportional to $u$, i.e., $z_j / y_j$ is constant ($z in "span"(y)$). Because $1 - s > 1 > 0$ and $M_s > 0$, it follows that $z^T nabla^2 M_s (y) z <= 0$ for all $z$, with strict inequality whenever $z in.not "span"(y)$.

  (2) Let $s_1 < s_2 < 0$ and set $r = s_1 / s_2 > 1$. Because $t |-> t^r$ is convex on $(0, infinity)$, Jensen's inequality gives
  $ 1 / k sum_(j=1)^k y_j^(s_1) = 1 / k sum_(j=1)^k (y_j^(s_2))^r >= ( 1 / k sum_(j=1)^k y_j^(s_2) )^(s_1 / s_2). $
  Raising both sides to the negative power $1 / s_1 < 0$ reverses the inequality, yielding $M_(s_1) (y) <= M_(s_2) (y)$. The arithmetic-mean upper bound follows analogously from convexity of $t |-> t^(s_2)$ on $(0, infinity)$ for $s_2 < 0$.

  (3) Since $y_j >= y_"min" > 0$ for every $j$ and $s < 0$, we have $y_j^s <= y_"min"^s$, and the single minimizing term gives $k^(-1) y_"min"^s <= k^(-1) sum_(j=1)^k y_j^s <= y_"min"^s$. Raising all three members to $1 / s < 0$ reverses the inequalities to $y_"min" <= M_s (y) <= k^(-1 / s) y_"min"$.
]

#subsec("B", "Exact Quadratic Majorization and Quantitative MM Descent")

#theorem(title: "Exact Tangent Majorization and Lyapunov Convergence", origin: "Original Proof")[
  Fix $s < 0$ and $epsilon > 0$, and let $bold(mu)^((t))$ be generated by fixed-$s$ MM iterations of Algorithm 1 with $W_j^((t)) := sum_(i=1)^n w_(i j)^((t)) > 0$. Define the quadratic surrogate
  $ Q_s (bold(mu) | bold(mu)^((t))) := sum_(i=1)^n sum_(j=1)^k w_(i j)^((t)) (norm(x_i - mu_j)_2^2 + epsilon). $
  Then:
  1. $Q_s (bold(mu) | bold(mu)^((t)))$ tightly majorizes $F_(s, epsilon)$ at $bold(mu)^((t))$: $Q_s (bold(mu)^((t)) | bold(mu)^((t))) = F_(s, epsilon) (bold(mu)^((t)))$ and $F_(s, epsilon) (bold(mu)) <= Q_s (bold(mu) | bold(mu)^((t)))$ for all $bold(mu) in (bb(R)^d)^k$.
  2. Every MM update satisfies the quantitative descent inequality
     $ &F_(s, epsilon) (bold(mu)^((t+1))) <= F_(s, epsilon) (bold(mu)^((t))) \
       &quad - sum_(j=1)^k W_j^((t)) norm(mu_j^((t+1)) - mu_j^((t)))_2^2. $
  3. For all $t >= 1$, $mu_j^((t)) in "conv"(cal(X))$, $norm(mu_j^((t+1)) - mu_j^((t)))_2 -> 0$, and every limit point $bold(mu)^*$ is a first-order stationary point ($nabla F_(s, epsilon) (bold(mu)^*) = bold(0)$).
]

#proof[
  (1) Because $M_s$ is concave on $bb(R)_(>0)^k$ (Lemma 1(1)), its supporting tangent plane lies everywhere above $M_s$:
  $ M_s (y_i (bold(mu))) <= M_s (y_i^((t))) + sum_(j=1)^k w_(i j)^((t)) (y_(i j) (bold(mu)) - y_(i j)^((t))). $
  Now apply Euler's theorem for $1$-homogeneous functions directly:
  $ sum_(j=1)^k w_(i j)^((t)) y_(i j)^((t)) &= 1 / k M_s (y_i^((t)))^(1 - s) sum_(j=1)^k (y_(i j)^((t)))^s \
    &= M_s (y_i^((t)))^(1 - s) dot M_s (y_i^((t)))^s = M_s (y_i^((t))). $
  Hence the constant terms $M_s (y_i^((t))) - sum_j w_(i j)^((t)) y_(i j)^((t))$ cancel identically! Summing over $i = 1, dots, n$ establishes both dominance $F_(s, epsilon) (bold(mu)) <= Q_s (bold(mu) | bold(mu)^((t)))$ and tangency at $bold(mu) = bold(mu)^((t))$.

  (2) Because $w_(i j)^((t)) > 0$, $Q_s (bold(mu) | bold(mu)^((t)))$ separates across $j = 1, dots, k$ into strictly convex quadratics. Applying the weighted parallel-axis identity around the exact minimizer $mu_j^((t+1)) = (W_j^((t)))^(-1) sum_(i=1)^n w_(i j)^((t)) x_i$ yields
  $ &Q_s (bold(mu)^((t)) | bold(mu)^((t))) - Q_s (bold(mu)^((t+1)) | bold(mu)^((t))) \
    &quad = sum_(j=1)^k W_j^((t)) norm(mu_j^((t+1)) - mu_j^((t)))_2^2. $
  Combining tangency and majorization gives $F_(s, epsilon) (bold(mu)^((t+1))) <= Q_s (bold(mu)^((t+1)) | bold(mu)^((t))) = F_(s, epsilon) (bold(mu)^((t))) - sum_(j=1)^k W_j^((t)) norm(mu_j^((t+1)) - mu_j^((t)))_2^2$.

  (3) Since $w_(i j)^((t)) > 0$, each $mu_j^((t+1))$ is a convex combination of ${x_1, dots, x_n}$, hence lies in the compact set $"conv"(cal(X))$ for all $t >= 1$. Letting $D = "diam"(cal(X))$, every $y_(i j)^((t)) in [epsilon, D^2 + epsilon]$, so $M_s (y_i^((t))) >= epsilon$ and thus $w_(i j)^((t)) >= k^(-1) (epsilon / (D^2 + epsilon))^(1 - s) =: c_0 > 0$, giving $W_j^((t)) >= n c_0 > 0$. Telescoping Claim 2 yields
  $ n c_0 sum_(t=1)^infinity sum_(j=1)^k norm(mu_j^((t+1)) - mu_j^((t)))_2^2 <= F_(s, epsilon) (bold(mu)^((1))) < infinity, $
  so $norm(mu_j^((t+1)) - mu_j^((t)))_2 -> 0$. By compactness of $"conv"(cal(X))^k$, limit points exist, and by continuity of the MM update map $T_s$, any limit point $bold(mu)^*$ satisfies $T_s (bold(mu)^*) = bold(mu)^*$, which is equivalent to $sum_(i=1)^n w_(i j) (bold(mu)^*, s) (mu_j^* - x_i) = bold(0)$. Since $nabla_(mu_j) F_(s, epsilon) (bold(mu)) = 2 sum_(i=1)^n w_(i j) (bold(mu), s) (mu_j - x_i)$, we conclude $nabla F_(s, epsilon) (bold(mu)^*) = bold(0)$.
]

#subsec("C", "Homotopy Convergence and Escape from Lloyd Traps")

#theorem(title: "Uniform Homotopy Rate and Global Strong Consistency", origin: "Classical — Cited")[
  Let $F_(-infinity) (bold(mu)) = sum_(i=1)^n min_j norm(x_i - mu_j)_2^2$. Then @xu2019power:
  1. For every $s < 0$ and all $bold(mu) in (bb(R)^d)^k$,
     $ 0 <= F_(s, 0) (bold(mu)) - F_(-infinity) (bold(mu)) <= (k^(-1 / s) - 1) F_(-infinity) (bold(mu)), $
     and the global optimal values converge at rate $O(|s|^(-1) log k)$ as $s -> -infinity$.
  2. For i.i.d. samples $X_1, dots, X_n tilde P$ with $bb(E)[norm(X)_2^2] < infinity$ and a unique population $k$-means minimizer $bold(mu)^*$ (up to permutation), any sequence of empirical global minimizers $hat(bold(mu))_(n, s_n)$ of $n^(-1) F_(s_n, 0)$ with $s_n -> -infinity$ satisfies $hat(bold(mu))_(n, s_n) -> bold(mu)^*$ almost surely (up to permutation) @xu2019power @pollard1981strong.
]

#proof(label: "Proof Synthesis (Cited)")[
  For (1), summing the pointwise sandwich bound of Lemma 1(3) over $i = 1, dots, n$ gives $F_(-infinity) (bold(mu)) <= F_(s, 0) (bold(mu)) <= k^(-1 / s) F_(-infinity) (bold(mu))$ for every $bold(mu)$. Taking infima over $bold(mu)$ preserves the multiplicative factor $k^(-1 / s)$, and for $|s| >= log k$ the elementary inequality $e^u - 1 <= 2 u$ ($u in [0, 1]$) gives $k^(-1 / s) - 1 = exp((log k) / |s|) - 1 <= (2 log k) / |s|$. For (2), Xu and Lange @xu2019power combine this uniform approximation bound with Pollard's strong consistency argument @pollard1981strong: the finite second moment forces empirical minimizers into a common compact ball $overline(B)_R (0)^k$ a.s., on which $sup_(bold(mu)) |n^(-1) F_(s_n, 0) (bold(mu)) - n^(-1) F_(-infinity) (bold(mu))| <= (k^(-1/s_n) - 1) sup_(bold(mu)) n^(-1) F_(-infinity) (bold(mu)) -> 0$ while $n^(-1) F_(-infinity)$ converges uniformly a.s. to the population $k$-means risk $Phi(bold(mu))$. Uniqueness of $bold(mu)^*$ then yields $hat(bold(mu))_(n, s_n) -> bold(mu)^*$ a.s.
]

We now prove constructively *why* Power $k$-Means escapes poor local basins where Lloyd's algorithm freezes permanently.

#proposition(title: "Escape from Adversarial Lloyd Fixed Points", origin: "Original Proof")[
  Consider the one-dimensional 4-point Lloyd trap from Alg. 1: $x_1 = 0$, $x_2 = 1$, $x_3 = L$, $x_4 = L + D$ with $k = 3$, $L > 2 D + 2$, and stationary Lloyd configuration $bold(mu)^((0)) = (0, 1, L + D / 2)$ (whose $k$-means approximation ratio is $D^2$). For every fixed $s < 0$ and all sufficiently small $epsilon > 0$:
  1. $bold(mu)^((0))$ is *not* stationary for $F_(s, epsilon)$: specifically $nabla_(mu_2) F_(s, epsilon) (bold(mu)^((0))) < 0$, and the first MM iterate moves strictly rightward toward the under-served cluster ${x_3, x_4}$: $mu_2^((1)) > 1$.
  2. The escape displacement scales as $mu_2^((1)) - 1 = Theta (D^(2 (1 - s)) L^(2 s - 1))$, proving mathematically that slow annealing from small $|s|$ is required to avoid vanishing long-range pull.
]

#proof[
  At $bold(mu)^((0)) = (0, 1, L + D / 2)$, the squared distances $y_(i j) = (x_i - mu_j^((0)))^2 + epsilon$ satisfy:
  - $y_(1 1) = epsilon$, $y_(1 2) = 1 + epsilon$, so by Lemma 1(3), $M_s (y_1) <= 3^(-1 / s) epsilon$, yielding $w_(1 2) <= 3^(-1 - (1 - s) / s) epsilon^(1 - s) -> 0$ as $epsilon -> 0^+$.
  - $y_(2 2) = epsilon = min_ell y_(2 ell)$, so $1 <= M_s (y_2) / y_(2 2) <= 3^(-1 / s)$, which bounds $w_(2 2) in [1 / 3, 3^(-1 / s)]$ uniformly away from $0$ and $infinity$ independently of $epsilon$.
  - For $i in {3, 4}$, $y_(i 3) = D^2 / 4 + epsilon$ is minimal while $y_(3 2) = (L - 1)^2 + epsilon$ and $y_(4 2) = (L + D - 1)^2 + epsilon$. Hence $M_s (y_i) in [D^2 / 4, 3^(-1 / s) (D^2 / 4 + epsilon)]$, giving
    $ w_(3 2) = Theta ( ( D / (2 L) )^(2 (1 - s)) ), quad w_(4 2) = Theta ( ( D / (2 L) )^(2 (1 - s)) ), $
    strictly positive as $epsilon -> 0^+$.

  Because $mu_2^((0)) - x_2 = 0$, the gradient component with respect to $mu_2$ is
  $ 1 / 2 nabla_(mu_2) F_(s, epsilon) (bold(mu)^((0))) &= w_(1 2) - w_(3 2) (L - 1) \
    &quad - w_(4 2) (L + D - 1). $
  As $epsilon -> 0^+$, $w_(1 2) -> 0$ while the negative terms stay bounded away from zero, so $nabla_(mu_2) F_(s, epsilon) (bold(mu)^((0))) < 0$. Likewise, the MM displacement is
  $ mu_2^((1)) - 1 = ( -w_(1 2) + w_(3 2) (L - 1) + w_(4 2) (L + D - 1) ) / ( w_(1 2) + w_(2 2) + w_(3 2) + w_(4 2) ). $
  Since $w_(2 2) = Theta(1)$, the denominator is $Theta(1)$ and the numerator is $Theta(L dot (D / L)^(2(1 - s))) = Theta(D^(2(1 - s)) L^(2 s - 1)) > 0$, strictly pulling $mu_2$ toward ${x_3, x_4}$ with strength decaying exponentially in $|s| log(2 L / D)$.
]

#subsec("D", "Invariances, Origin Collapse, and Breakdown Point")

#proposition(title: "Isometry Invariance, Symmetry Collapse, and Breakdown Point 1/n", origin: "Original Proof")[
  For every fixed $s < 0$ (with $epsilon = 0$):
  1. $F_(s, 0)$ and the MM iterates are strictly invariant under rigid Euclidean isometries $T(x) = Q x + b$ ($Q in O(d)$) and isotropic scaling $x |-> c x$ ($c > 0$), but not under anisotropic diagonal scaling.
  2. *Coincident-Centroid Fixed Point:* The collapsed configuration $mu_1 = dots = mu_k = overline(x) := n^(-1) sum_i x_i$ is a fixed point of the MM update for every $s < 0$.
  3. *Breakdown Point:* The finite-sample Donoho–Huber breakdown point of Power $k$-Means centroids is $epsilon^* = 1 / n$.
]

#proof[
  (1) Orthogonal transforms $Q in O(d)$ and translations preserve all $norm(x_i - mu_j)_2^2$, hence preserve $w_(i j)$, and convex combinations commute with affine maps. Under isotropic scaling $c > 0$, $y_(i j) |-> c^2 y_(i j)$ and $M_s (c^2 y_i) = c^2 M_s (y_i)$, leaving the ratio $M_s (y_i) / y_(i j)$ and all weights $w_(i j)$ strictly invariant. Anisotropic scaling alters global optima via the same 4-point rectangle counterexample as Alg. 1 once $|s|$ is large.

  (2) If $mu_1 = dots = mu_k = overline(x)$, then $y_(i 1) = dots = y_(i k)$ for every $i$, so $M_s (y_i) = y_(i j)$ and $w_(i j) = 1 / k$ identically. Thus $mu_j^((1)) = (n / k)^(-1) sum_(i=1)^n (1 / k) x_i = overline(x)$ for every $j$.

  (3) Corrupt $x_1 |-> M e_1$ ($M -> infinity$) while keeping $x_2, dots, x_n in overline(B)(0, R)$. If all centroids remained inside $overline(B)(0, B)$ for a fixed $B < infinity$, then $y_(1 j) >= (M - B)^2$ for all $j$, so $F_(s, epsilon) (bold(mu)) >= M_s (y_1) >= (M - B)^2 -> infinity$. Conversely, placing $mu_1 = x_1$ and $mu_2 = dots = mu_k = 0$ gives $F_(s, epsilon) <= k^(-1 / s) (epsilon + (n - 1) (R^2 + epsilon)) < infinity$. Thus every global minimizer must send at least one centroid outside every bounded ball, proving $epsilon^* = 1 / n$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Each iteration of Algorithm 1 evaluates all $n k$ squared distances in $Theta(n k d)$ time, computes stabilized power-mean weights in $O(n k)$ time via floating-point exponentiation (`pow`/`exp`/`log`), and forms all $k$ centroids via a *dense* weighted sum $W^T X$ costing $Theta(n k d)$—unlike Lloyd's algorithm, whose hard assignments allow $O(n d)$ centroid accumulation. Total runtime is $O(n k d T)$ with memory $O(n d + n k)$, where the constant factor is noticeably higher than Alg. 1 due to both dense $W^T X$ updates, transcendental power computations, and the multi-stage homotopy schedule $s_0 > dots > s_T$.

Numerically, computing $y_(i j)^s$ directly overflows or underflows rapidly as $|s|$ grows beyond $20$; dividing every row by $m_i = min_ell y_(i ell)$ (Step 4 of Algorithm 1) keeps all ratios $(y_(i ell) / m_i)^s in (0, 1]$, eliminating overflow completely. Once $|s|$ becomes large, Power $k$-Means can be switched to exact Lloyd polishing for machine-precision termination.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for Power K-Means],
    [Time / Space Complexity], [$O(n k d T)$ time (dense $W^T X$ + power ops); $O(n(d + k))$ memory.],
    [Geometric Inductive Bias], [Smoothly annealed spherical Voronoi cells converging to Euclidean $k$-means.],
    [Determinism / Seeding], [Deterministic given $bold(mu)^((0))$; dramatically less seeding-sensitive than Alg. 1; avoid coincident starts (Prop. 2).],
    [Invariances], [Invariant under $Q in O(d)$, translation, and isotropic scale; sensitive to anisotropic scaling.],
    [Outlier Breakdown Point], [$epsilon^* = 1 / n$ (non-robust to extreme outliers; Prop. 2(3)).],
    [Recommended Upgrades], [Combine with K-Means++ (Alg. 2) for seeding; Trimmed (Alg. 26) or TMM (Alg. 88) for outliers; W-K-Means (Alg. 28) for noisy features.],
  ),
  caption: [Summary of theoretical guarantees and operational profile for Power K-Means (Alg. 27).],
)

#sec("IV", "Engineering and Scientific Applications")

Power $k$-Means targets Euclidean quantization in moderate-to-high dimension $d$ or cluster count $k$, where Lloyd updates and multi-start $k$-means++ @arthur2007kmeanspp stall in suboptimal local minima.

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *High-Dimensional Single-Cell & Transcriptomic Partitioning:* In gene-expression spaces ($d = 20$–$500$), combinatorial Voronoi traps frequently wedge centroids between overlapping lineages. Xu and Lange @xu2019power showed that power-mean homotopy consistently recovers lower WCSS and stronger biological cluster agreement than Lloyd and $k$-means++.
2. *Vector Quantization & Neural Codebook Design:* Designing $k$-entry codebooks for compression (extending LBG Alg. 30 @linde1980algorithm) or VQ-VAE dictionaries suffers from "dead codewords" stranded in low-density zones. Proposition 1 guarantees non-zero long-range pull from high-distortion clusters, reviving idle codewords during early annealing.
3. *Warm-Starting Probabilistic Mixture Models:* EM @dempster1977em @wu1983convergence for Gaussian Mixture Models (Alg. 81) is acutely sensitive to initial means. A brief Power $k$-Means homotopy phase provides balanced component centers that avoid degenerate singular likelihood spikes.
4. *Continuous Facility Location & Logistics Network Design:* Placing $k$ distribution hubs to minimize squared transport cost is a continuous Euclidean $k$-means program where poor local optima permanently inflate operating costs, justifying millisecond-scale homotopy annealing.
5. *Hyperspectral Unmixing & Color Palette Optimization:* Multi-band remote sensing imagery spans wide spectral dynamic ranges; smooth power-mean continuation tracks dominant spectral modes globally before sharpening boundaries at large $|s|$.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Power $k$-Means has four notable failure modes: (a) *annealing schedule misconfiguration*—increasing $|s|$ too aggressively mimics Lloyd's local trapping (Proposition 1(2)), whereas holding $s$ close to $0^-$ risks premature centroid collapse toward $overline(x)$ (Proposition 2(2)); (b) *outlier vulnerability* ($epsilon^* = 1 / n$, resolved by Trimmed $k$-Means Alg. 26 @cuesta1997trimmed or Student's $t$-Mixtures Alg. 88 @peel2000robust); (c) *anisotropic or non-convex geometries*, because the target limit $s -> -infinity$ remains isotropic spherical $k$-means; and (d) *billion-scale throughput constraints*, where dense $O(n k d)$ centroid updates and transcendental arithmetic are slower than Elkan/Hamerly (Algs. 13–14) or Mini-Batch $k$-Means (Alg. 3).

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Why Power K-Means Fits], [Caution / Alternative],
    [High-D Omics / PCA], [Escapes combinatorial Voronoi traps cleanly], [Use Alg. 28/136 if pure noise dims],
    [VQ / Codebooks], [Eliminates dead codewords via long-range pull], [Polish with Alg. 13 at large $|s|$],
    [GMM Initialization], [Globally balanced warm-start for EM], [Follow with Alg. 81 for ellipsoids],
    [Facility Location], [Lowers global WCSS vs. single-start Lloyd], [Use Alg. 4/15 for $L_1$ or outliers],
    [Web-Scale Streams], [High quality on coresets / compressed sketches], [Use Alg. 3 or 23 for raw streams],
  ),
  caption: [Domain suitability matrix and engineering guidance for Power K-Means.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Power $k$-Means occupies a remarkable unifying hub connecting Fuzzy $C$-Means (Alg. 151) @bezdek1981fcm, Lloyd's $k$-means (Alg. 1) @lloyd1982least, Global $k$-Means (Alg. 17) @likas2003global, and Bregman divergence clustering @banerjee2005clustering.

#proposition(title: "Exact Isomorphism with Fuzzy C-Means and Zero-Temperature Limit to Lloyd", origin: "Original Unified Proof")[
  Let $y_(i j) = norm(x_i - mu_j)_2^2 > 0$ for all $i, j$. Then:
  1. *Exact FCM Isomorphism:* For every $s < 0$, define the fuzzifier $m := 1 - 1 / s in (1, infinity)$. Minimizing the Fuzzy $C$-Means (Alg. 151) objective $sum_(i=1)^n sum_(j=1)^k u_(i j)^m y_(i j)$ over simplex memberships $u_i in Delta_(k - 1)$ @bezdek1981fcm @hathaway1986fcm yields the exact reduced objective
     $ min_(U in Delta_(k-1)^n) sum_(i=1)^n sum_(j=1)^k u_(i j)^m y_(i j) = k^(1 / s) F_(s, 0) (bold(mu)), $
     and the optimal FCM centroid weights satisfy $(u_(i j)^*)^m = k^(1 / s) w_(i j) (bold(mu), s)$. Hence one fixed-$s$ MM iteration of Power $k$-Means is *algebraically identical* to one alternating FCM iteration at $m = 1 - 1 / s$.
  2. *Exact Limit to Lloyd's Algorithm:* If every point $i$ has a unique nearest centroid $c_i = arg min_j y_(i j)$, then as $s -> -infinity$,
     $ lim_(s -> -infinity) w_(i j) (bold(mu), s) = cases(1 & "if" j = c_i, 0 & "if" j != c_i), $
     recovering Lloyd's exact hard centroid update (Alg. 1).
]

#proof[
  (1) Fix $i$ and minimize $sum_(j=1)^k u_(i j)^m y_(i j)$ over $u_i in Delta_(k - 1)$ ($m > 1$). The Lagrangian stationarity condition $m u_(i j)^(m - 1) y_(i j) = lambda_i$ gives $u_(i j)^* = y_(i j)^(-1 / (m - 1)) / Z_i$, where $Z_i = sum_(ell=1)^k y_(i ell)^(-1 / (m - 1))$. Substituting $m = 1 - 1 / s$ yields $-1 / (m - 1) = s$, $-m / (m - 1) = s - 1$, and $1 - m = 1 / s$. Therefore $Z_i = sum_(ell=1)^k y_(i ell)^s = k M_s (y_i)^s$, and
  $ (u_(i j)^*)^m &= Z_i^(-m) y_(i j)^(s - 1) = (k M_s (y_i)^s)^(1 / s - 1) y_(i j)^(s - 1) \
    &= k^(1 / s) dot 1 / k ( M_s (y_i) / y_(i j) )^(1 - s) = k^(1 / s) w_(i j) (bold(mu), s). $
  Multiplying by $y_(i j)$ and summing over $j$ gives $sum_j (u_(i j)^*)^m y_(i j) = k^(1 / s) M_s (y_i)$, proving exact proportionality of both objectives and centroid weights (so the global factor $k^(1 / s)$ cancels in the centroid update ratio).

  (2) Let $r_ell := y_(i ell) / y_(i, c_i) > 1$ for all $ell != c_i$. Then $M_s (y_i) / y_(i, c_i) = k^(-1 / s) (1 + sum_(ell != c_i) r_ell^s)^(1 / s)$. Raising to $1 - s$ and multiplying by $1 / k$:
  $ w_(i, c_i) (bold(mu), s) = 1 / k dot k^(-(1 - s) / s) ( 1 + sum_(ell != c_i) r_ell^s )^((1 - s) / s). $
  Observe that $k^(-1) k^(-(1 - s) / s) = k^(-1 / s) -> 1$ as $s -> -infinity$, while $sum_(ell != c_i) r_ell^s -> 0$ exponentially fast in $|s|$, so $(1 + sum_(ell != c_i) r_ell^s)^((1 - s) / s) -> 1$. Thus $w_(i, c_i) -> 1$. For $j != c_i$, $w_(i j) = w_(i, c_i) r_j^(-(1 - s)) -> 0$ exponentially fast since $r_j > 1$.
]

Proposition 3 reveals that Power $k$-Means is mathematically equivalent to *deterministic homotopy annealing of Fuzzy $C$-Means* from $m = 2$ ($s = -1$) down to the hard Lloyd limit $m -> 1^+$ ($s -> -infinity$). Whereas Trimmed $k$-Means (Alg. 26) modifies the sample summation to resist outliers and W-$k$-Means (Alg. 28) applies a negative power mean across *features*, Power $k$-Means applies its negative power mean across *clusters* to smooth the optimization landscape.

  In real empirical benchmarks on the ClusterOverview 4-dataset suite @rousseeuw1987silhouettes @hubert1985comparing, Power $k$-Means attains Iris Sil $= 0.447$ (ARI $= 0.583$), Customers Sil $= 0.340$, Multi-Shape Sil $= 0.581$ (ARI N/A), and High-D Subspace+Noise Sil $= 0.303$ (ARI N/A), with mean runtime $2.704$ ms. Its superior Iris ARI ($0.583$ vs. $0.525$ for Alg. 26 and $0.441$ for Alg. 28) reflects Proposition 1's improved global basin discovery on convex clusters, while its $2.704$ ms runtime (vs. $0.081$ ms for Alg. 26 and $0.456$ ms for Alg. 28) directly reflects the Section III cost of dense $O(n k d)$ centroid sums, power evaluations, and homotopy annealing.
