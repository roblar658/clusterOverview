#import "../template/ieee.typ": *

#alg_entry(
  id: 174,
  name: "Multi-View K-Means",
  category: "Cat. VII: Fuzzy, Biclustering & Domain",
  year: "2013",
  paradigm: "Adaptive-Weight Multi-View Partitional Clustering",
  complexity: $cal(O)(T n k D)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Multi-view datasets describe a common set of $n$ objects through $V >= 2$ heterogeneous feature representations (views) $bold(X)_1, ..., bold(X)_V$, where $bold(X)_v in bb(R)^(n times d_v)$ has $i$-th row $bold(x)_(v, i)^top in bb(R)^(1 times d_v)$ and total ambient dimension $D = sum_(v=1)^V d_v$. Naive feature concatenation treats all views with uniform scale and noise variance, allowing high-dimensional or noisy views to dominate squared Euclidean distances. *Multi-View K-Means* (MVKM) @cai2013multiview resolves this by coupling $V$ view-specific prototype sets to a *single shared* discrete partition matrix while adaptively optimizing a simplex of view weights.

Let $cal(G)_(n, k) = {bold(G) in {0, 1}^(n times k) : bold(G) bold(1)_k = bold(1)_n}$ denote the set of hard $k$-cluster indicator matrices, with cluster counts $n_j = sum_(i=1)^n G_(i j)$ and equivalent assignment map $pi : {1, ..., n} -> {1, ..., k}$ where $G_(i, pi(i)) = 1$. For each view $v in {1, ..., V}$, let $bold(C)_v in bb(R)^(k times d_v)$ store view-specific centroids $bold(c)_(v, 1), ..., bold(c)_(v, k) in bb(R)^(d_v)$, and let $bold(alpha) = (alpha_1, ..., alpha_V)^top$ belong to the probability simplex $Delta_(V-1) = {bold(alpha) in bb(R)_(>= 0)^V : sum_(v=1)^V alpha_v = 1}$.

#definition(title: "Power-Weighted Multi-View K-Means Objective")[
  Given views ${bold(X)_v}_(v=1)^V$, cluster count $k >= 2$, and weighting exponent $gamma > 1$, the *Multi-View K-Means objective* minimizes the power-weighted sum of within-view squared Frobenius distortions
  $
    J_gamma (bold(G), {bold(C)_v}, bold(alpha))
    &= sum_(v=1)^V alpha_v^gamma E_v (bold(G), bold(C)_v), \
    E_v (bold(G), bold(C)_v)
    &= norm(bold(X)_v - bold(G) bold(C)_v)_F^2 \
    &= sum_(i=1)^n sum_(j=1)^k G_(i j) norm(bold(x)_(v, i) - bold(c)_(v, j))_2^2,
  $
  over $bold(G) in cal(G)_(n, k)$, $bold(C)_v in bb(R)^(k times d_v)$, and $bold(alpha) in Delta_(V-1)$. Replacing $E_v$ by the row-wise $ell_(2, 1)$ norm $norm(bold(X)_v - bold(G) bold(C)_v)_(2, 1) = sum_(i=1)^n norm(bold(x)_(v, i) - (bold(G) bold(C)_v)_i)_2$ yields the robust variant RMKMC @cai2013multiview.
]

The exponent $gamma > 1$ prevents trivial degeneracy: at $gamma = 1$, minimizing $sum_v alpha_v E_v$ over $Delta_(V-1)$ assigns mass $1$ exclusively to $arg min_v E_v$ and discards all other views. For $gamma > 1$, the strict convexity of $alpha_v |-> alpha_v^gamma$ balances view reliabilities smoothly, inducing a block-diagonal adaptive metric $d_(bold(alpha))^2 (i, j) = sum_(v=1)^V alpha_v^gamma norm(bold(x)_(v, i) - bold(c)_(v, j))_2^2$ across the product space $bb(R)^D$.

#algbox(
  title: "Alternating Multi-View K-Means (MVKM)",
  inputs: [Views ${bold(X)_v in bb(R)^(n times d_v)}_(v=1)^V$ (Frobenius-normalized), $k >= 2$, exponent $gamma > 1$, tolerance $tau > 0$, max iterations $T$.],
  outputs: [Shared indicator $bold(G) in cal(G)_(n, k)$, centroids ${bold(C)_v}$, weights $bold(alpha) in Delta_(V-1)$.],
)[
  + Initialize $bold(alpha)^((0)) = V^(-1) bold(1)_V$ and initial partition $bold(G)^((0))$ via $k$-means++ @arthur2007kmeanspp on $[alpha_1^(gamma / 2) bold(X)_1, ..., alpha_V^(gamma / 2) bold(X)_V]$; compute initial centroids $bold(C)_v^((0))$ as cluster means of $bold(G)^((0))$.
  + *For* $t = 1, ..., T$:
    - *Assignment step*: For each $i = 1, ..., n$, set
      $pi^((t)) (i) = arg min_(j in {1, ..., k}) sum_(v=1)^V (alpha_v^((t-1)))^gamma norm(bold(x)_(v, i) - bold(c)_(v, j)^((t-1)))_2^2$,
      retaining $pi^((t-1)) (i)$ on exact ties; assemble $bold(G)^((t))$.
    - *Centroid step*: For each view $v = 1, ..., V$ and non-empty cluster $j$, update
      $bold(c)_(v, j)^((t)) = (1 / n_j^((t))) sum_(i : pi^((t)) (i) = j) bold(x)_(v, i)$.
    - *Weight step*: Compute $E_v^((t)) = norm(bold(X)_v - bold(G)^((t)) bold(C)_v^((t)))_F^2$ and update
      $alpha_v^((t)) = (E_v^((t)))^(1 / (1 - gamma)) / sum_(w=1)^V (E_w^((t)))^(1 / (1 - gamma))$.
    - *Stopping check*: If $bold(G)^((t)) = bold(G)^((t-1))$ and $norm(bold(alpha)^((t)) - bold(alpha)^((t-1)))_1 <= tau$, break.
  + Return $(bold(G)^((t)), {bold(C)_v^((t))}, bold(alpha)^((t)))$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish closed-form weight duality and its exact power-mean characterization, Lyapunov monotonicity with finite-step termination, a constructive scale-hijacking counterexample, the majorization-minimization proof of the robust $ell_(2, 1)$ extension, and hardness / seeding guarantees.

#subsec("A", "Closed-Form Weight Optimality and Dual Concave Reduction")

#lemma(title: "Unique Simplex Minimizer and Exact Power-Mean Dual", origin: "Original Proof")[
  Fix $(bold(G), {bold(C)_v})$ such that $E_v = E_v (bold(G), bold(C)_v) > 0$ for all $v = 1, ..., V$, and let $gamma > 1$ with $p = 1 / (gamma - 1) > 0$. Then:
  + The constrained problem $min_(bold(alpha) in Delta_(V-1)) sum_(v=1)^V alpha_v^gamma E_v$ has a unique strictly positive global minimizer
    $
      alpha_v^star = (E_v^(-p)) / (sum_(w=1)^V E_w^(-p)),
      quad v = 1, ..., V.
    $
  + Eliminating $bold(alpha)$ yields the exact reduced objective
    $
      tilde(J)_gamma (bold(G), {bold(C)_v})
      &= min_(bold(alpha) in Delta_(V-1)) J_gamma (bold(G), {bold(C)_v}, bold(alpha)) \
      &= (sum_(v=1)^V E_v^(-p))^(-1 / p)
      = V^(1 - gamma) M_(-p) (E_1, ..., E_V),
    $
    where $M_r (E) = (V^(-1) sum_(v=1)^V E_v^r)^(1 / r)$ is the generalized power mean of order $r$. Moreover, $tilde(J)_gamma$ is concave in $bold(E) = (E_1, ..., E_V) in bb(R)_(> 0)^V$, and the normalized objective $V^(gamma - 1) tilde(J)_gamma$ equals the harmonic mean at $gamma = 2$, converges to $min_v E_v$ as $gamma -> 1^+$, and converges to the geometric mean $(product_(v=1)^V E_v)^(1 / V)$ as $gamma -> infinity$.
]

#proof[
  Consider $f(bold(alpha)) = sum_(v=1)^V alpha_v^gamma E_v$ on the compact convex simplex $Delta_(V-1)$. For $gamma > 1$ and $E_v > 0$, the Hessian is diagonal with entries $partial^2 f / partial alpha_v^2 = gamma (gamma - 1) alpha_v^(gamma - 2) E_v > 0$ on the interior, so $f$ is strictly convex on $Delta_(V-1)$ and attains a unique global minimum. Forming the Lagrangian $cal(L)(bold(alpha), nu) = sum_v alpha_v^gamma E_v - nu (sum_v alpha_v - 1)$, stationarity gives
  $
    gamma alpha_v^(gamma - 1) E_v = nu
    quad => quad
    alpha_v = (nu / (gamma E_v))^p.
  $
  Since $nu > 0$ is required to satisfy $sum_v alpha_v = 1$, every $alpha_v^star > 0$ lies strictly in the interior of $Delta_(V-1)$, so boundary multipliers vanish and normalizing yields $alpha_v^star = E_v^(-p) / Z$ with $Z = sum_(w=1)^V E_w^(-p)$.

  Substituting $alpha_v^star$ back into $J_gamma$ and using $gamma p - 1 = p$:
  $
    tilde(J)_gamma
    = sum_(v=1)^V (E_v^(-gamma p) E_v) / Z^gamma
    = (sum_(v=1)^V E_v^(-p)) / Z^gamma
    = Z^(1 - gamma)
    = Z^(-1 / p).
  $
  Factoring $Z = V (V^(-1) sum_v E_v^(-p))$ gives $Z^(-1 / p) = V^(-1 / p) M_(-p) (E_1, ..., E_V) = V^(1 - gamma) M_(-p) (E_1, ..., E_V)$. Because $tilde(J)_gamma (bold(E)) = min_(bold(alpha) in Delta_(V-1)) sum_v alpha_v^gamma E_v$ is the pointwise infimum of linear functions of $bold(E)$, it is concave on $bb(R)_(> 0)^V$. At $gamma = 2$, $p = 1$ gives the harmonic mean $M_(-1) (bold(E))$. As $gamma -> 1^+$, $p -> infinity$ and $M_(-p) (bold(E)) -> min_v E_v$. As $gamma -> infinity$, $p -> 0^+$; expanding $E_v^(-p) = 1 - p ln E_v + cal(O)(p^2)$ yields
  $
    ln M_(-p) (bold(E))
    = - (1 / p) ln (1 - (p / V) sum_(v=1)^V ln E_v + cal(O)(p^2))
    -> (1 / V) sum_(v=1)^V ln E_v,
  $
  establishing convergence to the geometric mean $(product_(v=1)^V E_v)^(1 / V)$ despite $alpha_v^star -> 1 / V$.
]

#subsec("B", "Monotonic Lyapunov Descent and Finite-Step Convergence")

#theorem(title: "Monotone Lyapunov Descent and Finite Termination", origin: "Original Proof")[
  Suppose $E_v (bold(G), bold(C)_v) > 0$ along all iterates of MVKM with $gamma > 1$ and stay-on-tie assignment breaking. Let $Phi (bold(G)) = tilde(J)_gamma (bold(G), bold(C)^star (bold(G)))$, where $bold(C)_v^star (bold(G))$ contains the exact cluster means of $bold(G)$ in view $v$. Then for every iteration $t >= 1$:
  + The objective satisfies the exact three-stage descent chain
    $
      Phi (bold(G)^((t)))
      &= J_gamma (bold(G)^((t)), bold(C)^((t)), bold(alpha)^((t)))
      <= J_gamma (bold(G)^((t)), bold(C)^((t)), bold(alpha)^((t-1))) \
      &<= J_gamma (bold(G)^((t)), bold(C)^((t-1)), bold(alpha)^((t-1)))
      <= Phi (bold(G)^((t-1))),
    $
    with exact centroid-update gain
    $
      Delta_"cent"^((t))
      = sum_(v=1)^V (alpha_v^((t-1)))^gamma sum_(j=1)^k n_j^((t))
        norm(bold(c)_(v, j)^((t)) - bold(c)_(v, j)^((t-1)))_2^2.
    $
  + Whenever $bold(G)^((t)) != bold(G)^((t-1))$, the descent is strict: $Phi (bold(G)^((t))) < Phi (bold(G)^((t-1)))$. Consequently, MVKM visits no indicator matrix twice and terminates in at most $|cal(G)_(n, k)| = k^n$ iterations at a block-coordinate stationary point.
]

#proof[
  At the end of iteration $t - 1$, both $bold(C)^((t-1)) = bold(C)^star (bold(G)^((t-1)))$ and $bold(alpha)^((t-1)) = bold(alpha)^star (bold(G)^((t-1)), bold(C)^((t-1)))$ are uniquely optimal for $bold(G)^((t-1))$, so $J_gamma (bold(G)^((t-1)), bold(C)^((t-1)), bold(alpha)^((t-1))) = Phi (bold(G)^((t-1)))$.

  *Stage 1 (Assignment).* Because $J_gamma (bold(G), bold(C)^((t-1)), bold(alpha)^((t-1)))$ decomposes additively over $i = 1, ..., n$, choosing $pi^((t)) (i)$ as the minimizer of $d_(bold(alpha)^((t-1)))^2 (i, j)$ guarantees $J_gamma (bold(G)^((t)), bold(C)^((t-1)), bold(alpha)^((t-1))) <= Phi (bold(G)^((t-1)))$. If $bold(G)^((t)) != bold(G)^((t-1))$, the stay-on-tie rule forces at least one $i$ to achieve a strictly smaller weighted distance, so the inequality is strict.

  *Stage 2 (Centroids).* For each view $v$ and cluster $j$ with index set $cal(I)_j = {i : pi^((t)) (i) = j}$ and mean $bold(c)_(v, j)^((t)) = (1 / n_j^((t))) sum_(i in cal(I)_j) bold(x)_(v, i)$, the cross term $sum_(i in cal(I)_j) chevron.l bold(x)_(v, i) - bold(c)_(v, j)^((t)), bold(c)_(v, j)^((t)) - bold(c)_(v, j)^((t-1)) chevron.r$ vanishes identically. Expanding the square gives the parallel-axis identity
  $
    sum_(i in cal(I)_j) norm(bold(x)_(v, i) - bold(c)_(v, j)^((t-1)))_2^2
    &= sum_(i in cal(I)_j) norm(bold(x)_(v, i) - bold(c)_(v, j)^((t)))_2^2 \
    &quad + n_j^((t)) norm(bold(c)_(v, j)^((t)) - bold(c)_(v, j)^((t-1)))_2^2.
  $
  Multiplying by $(alpha_v^((t-1)))^gamma >= 0$ and summing over $j$ and $v$ proves $J_gamma (bold(G)^((t)), bold(C)^((t-1)), bold(alpha)^((t-1))) - J_gamma (bold(G)^((t)), bold(C)^((t)), bold(alpha)^((t-1))) = Delta_"cent"^((t)) >= 0$.

  *Stage 3 (Weights).* By Lemma 1, $bold(alpha)^((t))$ uniquely minimizes $J_gamma (bold(G)^((t)), bold(C)^((t)), bold(alpha))$ over $Delta_(V-1)$, giving $Phi (bold(G)^((t))) <= J_gamma (bold(G)^((t)), bold(C)^((t)), bold(alpha)^((t-1)))$.

  Chaining the three stages yields $Phi (bold(G)^((t))) < Phi (bold(G)^((t-1)))$ whenever $bold(G)^((t)) != bold(G)^((t-1))$. Since $cal(G)_(n, k)$ has finite cardinality $k^n$, the strictly decreasing sequence $Phi (bold(G)^((t)))$ must halt at $bold(G)^((t)) = bold(G)^((t-1))$ within at most $k^n$ steps. At termination, $bold(G)^((t))$, $bold(C)^((t))$, and $bold(alpha)^((t))$ are each globally optimal given the other two blocks.
]

#subsec("C", "Scale Non-Invariance and Constructive View-Collapse Counterexample")

#proposition(title: "Orthogonal Invariance vs. Heterogeneous Scale Hijacking", origin: "Original Counterexample")[
  Under view-wise transformations $tilde(bold(X))_v = s_v bold(X)_v bold(Q)_v + bold(1)_n bold(b)_v^top$ with $bold(Q)_v in cal(O)(d_v)$, $bold(b)_v in bb(R)^(d_v)$, and $s_v > 0$:
  + $tilde(J)_gamma$ and $bold(alpha)^star$ are strictly invariant to arbitrary view-wise rotations, translations, and uniform rescaling ($s_1 = ... = s_V$).
  + Under heterogeneous view scaling ($s_v != s_w$), invariance fails catastrophically: for every $gamma > 1$ and $epsilon in (0, 1)$, there exists a $2$-view, $n = 4$, $k = 2$ instance with clean ground-truth partition $A$ in View 1 and orthogonal conflicting partition $B$ in View 2 (scaled by $s > 0$) whose global MVKM minimizer flips from $A$ to $B$ if and only if $s < 1$, and collapses $alpha_1^star -> 0$ as $s -> 0^+$.
]

#proof[
  For part 1, subtracting $bold(c)_(v, j)$ eliminates $bold(b)_v$, while orthogonal invariance of $norm(dot)_2$ eliminates $bold(Q)_v$. Thus $tilde(E)_v = s_v^2 E_v$. If all $s_v = s$, then $s^(-2 p)$ cancels in $tilde(alpha)_v^star = (s_v^(-2 p) E_v^(-p)) / sum_w (s_w^(-2 p) E_w^(-p))$, leaving $bold(alpha)^star$ unchanged and scaling $tilde(J)_gamma$ by $s^2$.

  For part 2, construct $n = 4$ points in $V = 2$ one-dimensional views ($d_1 = d_2 = 1$):
  $
    bold(x)_1 &= (-1 - epsilon, -1 + epsilon, 1 - epsilon, 1 + epsilon)^top, \
    bold(x)_2 &= s (-1 - epsilon, 1 - epsilon, -1 + epsilon, 1 + epsilon)^top.
  $
  Consider the true View-1 partition $A = {{1, 2}, {3, 4}}$ and the conflicting View-2 partition $B = {{1, 3}, {2, 4}}$. Direct evaluation of optimal within-cluster sums of squares gives:
  $
    E_1 (A) &= 4 epsilon^2, & E_2 (A) &= 4 s^2, \
    E_1 (B) &= 4, & E_2 (B) &= 4 s^2 epsilon^2.
  $
  By Lemma 1 with $p = 1 / (gamma - 1) > 0$,
  $
    tilde(J)_gamma (A) &= 4 (epsilon^(-2 p) + s^(-2 p))^(-1 / p), \
    tilde(J)_gamma (B) &= 4 (1 + s^(-2 p) epsilon^(-2 p))^(-1 / p).
  $
  Because $u |-> u^(-1 / p)$ is strictly decreasing on $(0, infinity)$, we have $tilde(J)_gamma (B) < tilde(J)_gamma (A)$ if and only if
  $
    1 + s^(-2 p) epsilon^(-2 p) - epsilon^(-2 p) - s^(-2 p)
    = (s^(-2 p) - 1) (epsilon^(-2 p) - 1) > 0.
  $
  Since $0 < epsilon < 1$ implies $epsilon^(-2 p) - 1 > 0$, the strict inequality holds if and only if $s^(-2 p) > 1$, i.e., $s < 1$, independently of $gamma > 1$. Moreover, for any fixed $bold(G)$, $alpha_1^star = (1 + (E_1 / (s^2 E_2))^p)^(-1) -> 0$ as $s -> 0^+$, confirming that simply shrinking a conflicting or pure-noise view hijacks the entire multi-view objective unless views are pre-normalized to comparable total variance.
]

#subsec("D", "Robust IRLS Equivalence and Finite-Sample Breakdown")

#theorem(title: "MM Convergence of Robust RMKMC", origin: "Classical — Cited")[
  Consider RMKMC @cai2013multiview minimizing $J_(2, 1) (bold(G), {bold(C)_v}, bold(alpha)) = sum_(v=1)^V alpha_v^gamma norm(bold(E)_v)_(2, 1)$ with residual matrix $bold(E)_v = bold(X)_v - bold(G) bold(C)_v$. Defining diagonal row-reweighting matrices $bold(W)_v^((t))$ with entries $w_(v, i)^((t)) = (2 norm(bold(e)_(v, i)^((t-1)))_2)^(-1)$ for non-zero residual rows $bold(e)_(v, i)^((t-1))$:
  + Updating $(bold(G)^((t)), {bold(C)_v^((t))})$ to decrease the weighted quadratic surrogate $sum_v (alpha_v^((t-1)))^gamma "Tr"((bold(X)_v - bold(G) bold(C)_v)^top bold(W)_v^((t)) (bold(X)_v - bold(G) bold(C)_v))$ followed by updating $bold(alpha)^((t))$ via Lemma 1 with $E_v = norm(bold(E)_v)_(2, 1)$ monotonically decreases $J_(2, 1)$.
  + Whereas squared-error MVKM has a single-view centroid breakdown point of $1 / n$, RMKMC restricts every outlier residual to bounded unit-norm gradient influence on $bold(C)_v$.
]

#proof[
  For any scalars $a >= 0$ and $b > 0$, $(a - b)^2 >= 0$ rearranges to $a - a^2 / (2 b) <= b - b^2 / (2 b)$. Setting $a_i = norm(bold(e)_(v, i)^((t)))_2$ and $b_i = norm(bold(e)_(v, i)^((t-1)))_2 > 0$, summing over $i = 1, ..., n$ yields
  $
    norm(bold(E)_v^((t)))_(2, 1) - norm(bold(E)_v^((t-1)))_(2, 1)
    <= sum_(i=1)^n w_(v, i)^((t)) (a_i^2 - b_i^2).
  $
  Weighting each view by $(alpha_v^((t-1)))^gamma > 0$ and summing over $v = 1, ..., V$ gives
  $
    J_(2, 1) (bold(G)^((t)), bold(C)^((t)), bold(alpha)^((t-1)))
    &- J_(2, 1) (bold(G)^((t-1)), bold(C)^((t-1)), bold(alpha)^((t-1))) \
    &<= sum_(v=1)^V (alpha_v^((t-1)))^gamma
      ("Tr"((bold(E)_v^((t)))^top bold(W)_v^((t)) bold(E)_v^((t))) \
    &quad - "Tr"((bold(E)_v^((t-1)))^top bold(W)_v^((t)) bold(E)_v^((t-1)))) <= 0.
  $
  The subsequent weight update $bold(alpha)^((t))$ minimizes $J_(2, 1)$ by Lemma 1, proving monotone descent of RMKMC @cai2013multiview. For robustness, under $J_gamma$ a single point $bold(x)_(v, 1) -> infinity$ pulls its centroid $bold(c)_(v, j) -> infinity$ (and if all views are corrupted, breaks all centroids at fraction $1 / n$). Under $J_(2, 1)$, the subgradient with respect to $bold(c)_(v, j)$ is $-sum_(i in cal(I)_j) bold(e)_(v, i) / norm(bold(e)_(v, i))_2$, so any single outlier contributes at most unit norm regardless of magnitude.
]

#subsec("E", "NP-Hardness and Multi-View Seeding Approximation")

#theorem(title: "NP-Hardness and Multi-View k-Means++ Approximation", origin: "Classical — Cited")[
  For every $V >= 1$ and $gamma > 1$:
  + Computing the global minimizer of $tilde(J)_gamma (bold(G), {bold(C)_v})$ over $cal(G)_(n, k)$ is NP-hard even for $k = 2$ (general $D$) and for $D = 2 V$ (general $k$).
  + Initializing $bold(G)^((0))$ via $k$-means++ @arthur2007kmeanspp on the weighted concatenation $bold(Z)_(bold(alpha)^((0))) = [(alpha_1^((0)))^(gamma / 2) bold(X)_1, ..., (alpha_V^((0)))^(gamma / 2) bold(X)_V]$ guarantees that the converged MVKM iterate satisfies
    $
      bb(E)[tilde(J)_gamma (bold(G)^((T)), bold(C)^((T)))]
      <= 8 (ln k + 2) min_(bold(G), {bold(C)_v}) J_gamma (bold(G), {bold(C)_v}, bold(alpha)^((0))).
    $
]

#proof[
  Given any single-view Euclidean $k$-means instance $bold(X) in bb(R)^(n times d)$, set $bold(X)_1 = ... = bold(X)_V = bold(X)$. For every partition $bold(G)$ and optimal centroids $bold(C)_v = bold(C)$, all view distortions coincide ($E_1 = ... = E_V = E$), so Lemma 1 yields $alpha_v^star = 1 / V$ and $tilde(J)_gamma = V^(1 - gamma) E$. Thus minimizing $tilde(J)_gamma$ is polynomially equivalent to Euclidean $k$-means, which is NP-hard for $k = 2$ @aloise2009np and for $d = 2$ @mahajan2012planar.

  For the approximation bound, let $bold(C)_(bold(alpha)^((0))) = [(alpha_1^((0)))^(gamma / 2) bold(C)_1, ..., (alpha_V^((0)))^(gamma / 2) bold(C)_V]$; then for fixed $bold(alpha)^((0))$,
  $
    J_gamma (bold(G), {bold(C)_v}, bold(alpha)^((0)))
    = norm(bold(Z)_(bold(alpha)^((0))) - bold(G) bold(C)_(bold(alpha)^((0))))_F^2,
  $
  which is exact Euclidean $k$-means on $bold(Z)_(bold(alpha)^((0)))$. By the $k$-means++ theorem @arthur2007kmeanspp, $bb(E)[J_gamma (bold(G)^((0)), bold(C)^((0)), bold(alpha)^((0)))] <= 8 (ln k + 2) J_gamma^star (bold(alpha)^((0)))$. Combining this with monotone descent (Theorem 2), $tilde(J)_gamma (bold(G)^((T)), bold(C)^((T))) <= J_gamma (bold(G)^((0)), bold(C)^((0)), bold(alpha)^((0)))$, yields the stated expectation bound.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Each iteration of MVKM avoids $n times n$ affinity matrices completely. Evaluating all $n k$ squared distances across $V$ views via the expansion $norm(bold(x)_(v, i) - bold(c)_(v, j))_2^2 = norm(bold(x)_(v, i))_2^2 + norm(bold(c)_(v, j))_2^2 - 2 bold(x)_(v, i)^top bold(c)_(v, j)$ reduces the assignment step to $V$ BLAS-3 matrix products $bold(X)_v bold(C)_v^top$, costing $cal(O)(n k D)$ flops where $D = sum_(v=1)^V d_v$. Updating centroids requires $cal(O)(n D)$ operations, and computing $E_v$ reuses the assigned distances in $cal(O)(n V)$ (or $cal(O)(n D)$ directly), yielding total time $cal(O)(T n k D)$ and linear memory $cal(O)(n D + k D)$. In worst-case adversarial geometries, alternating Lloyd steps can still require super-polynomial $T$ @vattani2011kmeans, though practical convergence typically completes within $10$–$40$ sweeps.

Geometrically, for fixed weights $bold(alpha)$, MVKM induces a convex Voronoi tessellation in $bb(R)^D$ under the block-diagonal Mahalanobis metric $bold(M)_(bold(alpha)) = "diag"(alpha_1^gamma bold(I)_(d_1), ..., alpha_V^gamma bold(I)_(d_V))$. All cluster boundaries are hyperplanes in the product feature space, so MVKM excels on convex, ellipsoidal, or well-separated clusters across views while requiring prior kernel or representation lifting for intertwined manifolds.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$cal(O)(T n k D)$, linear in $n$, $k$, and $D = sum_v d_v$],
    [Space Complexity], [$cal(O)(n D + k D + n + V)$ (no $n times n$ kernel storage)],
    [Decision Geometry], [Block-anisotropic Voronoi polytopes in $bb(R)^D$],
    [Hyperparameters], [Cluster count $k$, view exponent $gamma > 1$],
    [Scale Behavior], [Invariant to $cal(O)(d_v)$; requires view variance normalization],
    [Outlier Robustness], [Breakdown $1 / n$ for $J_gamma$; bounded gradient for $ell_(2, 1)$ RMKMC],
    [Convergence], [Finite-step monotone descent to block stationary point],
    [Parallelism], [Embarrassingly parallel over $n$ and views $v$ via GEMM],
  ),
  caption: [Operational and geometric profile of Multi-View K-Means (Alg. 174).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

+ *Large-Scale Multimodal Web Document Categorization*: Web pages combine lexical TF-IDF / embedding vectors, hyperlink-neighborhood features, and visual page descriptors. Because web corpora contain millions of pages, $cal(O)(n^3)$ spectral fusion is infeasible, whereas MVKM's $cal(O)(T n k D)$ complexity and adaptive $alpha_v^gamma$ weighting downweight noisy visual or OCR views automatically @cai2013multiview.
+ *Multi-Omics Patient Stratification*: Integrating gene expression, DNA methylation, and proteomics matrices over a shared patient cohort requires a single consensus clinical partition $bold(G)$ while keeping assay-specific prototypes $bold(C)_v$ interpretable in their native biological units. After view-wise standardization, $alpha_v^star$ directly ranks platform coherence.
+ *Multimodal Vector Quantization and Retrieval Indexing*: Product-space inverted-file codebooks over paired vision-language or audio-visual embeddings benefit directly from MVKM's block-diagonal Mahalanobis Voronoi cells, allowing sublinear nearest-neighbor routing with closed-form centroid lookup across modalities.
+ *Multi-Sensor Industrial Condition Monitoring*: Rotating machinery monitored simultaneously by vibration spectra, thermal signatures, and stator-current harmonics exhibits convex operating-regime clusters in each sensor space; the $ell_(2, 1)$ RMKMC variant resists impulsive sensor glitches while attenuating drifting channels.
+ *Omnichannel Customer Segmentation*: Combining demographic profiles, RFM transaction summaries, and high-dimensional clickstream histograms causes naive concatenation to be overwhelmed by clickstream dimensionality; Frobenius-normalized MVKM enforces balanced view participation under a unified marketing segmentation.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Four mathematical failure regimes govern deployment. First, *unnormalized heterogeneous view scales* trigger immediate scale-hijacking collapse (Proposition 1), forcing $alpha_v^star -> 1$ on whichever view has the smallest numerical variance regardless of semantic signal. Second, *low-exponent winner-take-all collapse* occurs as $gamma -> 1^+$, reducing MVKM to a single-view $k$-means solution on $arg min_v E_v$ (Lemma 1). Third, *non-convex manifold geometry* (concentric rings, moons, or sparse relational graphs) violates the hyperplane Voronoi bias, demanding spectral co-regularization (Alg. 175) instead. Fourth, *cross-view partition conflict*—where views encode genuinely distinct valid clusterings—cannot be represented by a single shared indicator $bold(G)$.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Multimodal Documents], [Linear $cal(O)(T n k D)$ scaling; noise attenuation], [Normalize TF-IDF & embedding norms],
    [Multi-Omics Cohorts], [Shared cohort $bold(G)$; interpretable $bold(C)_v$], [Use $ell_(2, 1)$ RMKMC for outliers],
    [Sensor Telemetry], [Convex regimes; channel downweighting], [Avoid $gamma approx 1$; verify stationarity],
    [Non-Convex Manifolds], [Poor fit: linear Voronoi boundaries], [Use Alg. 175 (Spectral Co-Reg.)],
    [Conflicting Views], [Single $bold(G)$ forces compromise cut], [Use alternative / subspace clustering],
  ),
  caption: [Domain suitability and failure-mitigation matrix for Alg. 174.],
)

#sec("V", "Structural Relations and Empirical Behavior")

#proposition(title: "Algebraic Reductions Within the Partitional and Multi-View Families", origin: "Original Derivation")[
  Multi-View K-Means embeds three classical algorithms as exact boundary cases:
  + Setting $V = 1$ recovers Lloyd's $k$-means (Alg. 1) @lloyd1982least identically, with strong consistency as $n -> infinity$ @pollard1981strong.
  + Assigning each scalar feature $r in {1, ..., D}$ its own $1$-dimensional view ($V = D$, $d_v = 1$) reduces $J_gamma$ identically to Automated Variable-Weighting $k$-means (W-$k$-means) @huang2005automated and links to Minkowski-weighted $k$-means @deamorim2012minkowski.
  + Fixing $alpha_v = 1 / V$ recovers concatenated Euclidean $k$-means, whereas replacing raw views $bold(X)_v$ by view-specific spectral embeddings $bold(U)_v$ turns MVKM into a linear-space consensus analogue of Co-regularized Multi-View Spectral Clustering (Alg. 175) @kumar2011coregularized and cluster ensembles @strehl2002cluster.
]

#proof[
  For $V = 1$, $Delta_0 = {1}$ fixes $alpha_1 = 1$, leaving $min_(bold(G), bold(C)_1) norm(bold(X)_1 - bold(G) bold(C)_1)_F^2$, which is Lloyd's objective. Setting $V = D$ and $bold(X)_r = bold(x)_(:, r) in bb(R)^(n times 1)$ gives $J_gamma = sum_(r=1)^D alpha_r^gamma sum_(i, j) G_(i j) (x_(i r) - c_(j r))^2$, matching W-$k$-means @huang2005automated verbatim. Finally, fixed $alpha_v = 1 / V$ scales concatenated $k$-means by $V^(-gamma)$, while substituting $bold(X)_v = bold(U)_v$ partitions the concatenated spectral representation under view-adaptive weights.
]

Empirically, on the ClusterOverview four-dataset benchmark @jain2010data, MVKM exhibits strong geometric alignment and sub-millisecond throughput: it attains Silhouette `"Sil" = 0.476` @rousseeuw1987silhouettes and Adjusted Rand Index `"ARI" = 0.450` @hubert1985comparing on Iris, `"Sil" = 0.349` on Customers, `"Sil" = 0.604` (`"ARI" = "N/A"`) on Multi-Shape, and `"Sil" = 0.259` (`"ARI" = "N/A"`) on High-D, with a mean runtime of `0.416 ms`. Two quantitative contrasts with Co-regularized Multi-View Spectral Clustering (Alg. 175) directly validate Sections II–III: first, MVKM runs $43.3 times$ faster (`0.416 ms` vs. `18.015 ms`), reflecting the $cal(O)(T n k D)$ versus $cal(O)(T V n^3)$ complexity separation; second, on Customers (`0.349` vs. `0.032`) and High-D (`0.259` vs. `-0.025`), MVKM avoids the high-dimensional RBF kernel concentration @beyer1999nearest that degrades spectral eigengaps, while adaptive power weighting concentrates weight on compact views and directly supports Euclidean Silhouette separation. On Iris, moderate `"ARI" = 0.450` reflects partial overlap between Versicolor and Virginica when four correlated dimensions are split into separate views.
