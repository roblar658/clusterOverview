#import "../template/ieee.typ": *

#alg_entry(
  id: 150,
  name: "Sparse Subspace Clustering (SSC)",
  category: "Cat. VI: Grid, Subspace & Projection",
  year: "2013",
  paradigm: "Self-Expressive L1-Minimization and Spectral Clustering",
  complexity: $O(T min(n, D) n^2 + n^3)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Many high-dimensional scientific signals---ranging from face images under varying illumination to rigid-body motion trajectories---do not cluster into axis-parallel hyper-rectangles (as in STATPC, Alg. 149) or isotropic Euclidean balls @lloyd1982least, but instead lie near a *union of low-dimensional linear or affine subspaces* embedded in $RR^D$. Sparse Subspace Clustering (SSC) @elhamifar2013ssc exploits the *self-expressiveness property* of linear subspaces: every data point in a union of subspaces can be represented as a linear combination of other points from the *same* subspace, and convex $ell_1$-minimization provably recovers such subspace-preserving representations under broad geometric conditions @soltanolkotabi2012geometric.

Let $bold(X) = [bold(x)_1, dots.c, bold(x)_n] in RR^(D times n)$ be a data matrix whose columns are drawn from a union of $K$ unknown subspaces $union.big_(k=1)^K cal(S)_k subset.eq RR^D$ with intrinsic dimensions $d_k = dim(cal(S)_k) << D$ and cardinalities $n_k >= d_k + 1$ ($sum_(k=1)^K n_k = n$). Because $n > d_k$, each $bold(x)_i in cal(S)_k$ admits infinitely many representations $bold(x)_i = bold(X) bold(c)_i$ with $c_(i i) = 0$. Minimizing the $ell_0$ pseudo-norm $norm(bold(c)_i)_0$ is NP-hard, so noiseless linear SSC solves the tightest convex relaxation, Basis Pursuit:
$ min_(bold(C) in RR^(n times n)) norm(bold(C))_1 quad "s.t." quad bold(X) = bold(X) bold(C), quad "diag"(bold(C)) = bold(0), $
where $norm(bold(C))_1 = sum_(i, j = 1)^n |C_(i j)|$. For *affine* subspaces $cal(A)_k = bold(mu)_k + cal(S)_k$ not passing through the origin, SSC appends the affine normalization constraint $bold(1)^top bold(C) = bold(1)^top$. When observations are corrupted by bounded Gaussian noise $bold(Z)$ and sparse gross outliers $bold(E)$, SSC solves the regularized convex program @elhamifar2013ssc:
$ min_(bold(C), bold(E), bold(Z)) norm(bold(C))_1 + lambda_e norm(bold(E))_1 + frac(lambda_z, 2) norm(bold(Z))_F^2 \ "s.t." quad bold(X) = bold(X) bold(C) + bold(E) + bold(Z), quad "diag"(bold(C)) = bold(0). $
In the pure Gaussian-noise LASSO specialization ($bold(E) = bold(0)$, $lambda = 1 / lambda_z$), the problem decouples across columns $i = 1, dots.c, n$ into
$ min_(bold(c)_i in RR^n : c_(i i) = 0) frac(1, 2) norm(bold(x)_i - bold(X) bold(c)_i)_2^2 + lambda norm(bold(c)_i)_1. $

#definition(title: "Subspace Detection Property (SDP)")[
  Let $kappa : \{1, dots.c, n\} -> \{1, dots.c, K\}$ map each column $bold(x)_i$ to its ground-truth subspace $cal(S)_(kappa(i))$. A self-expression matrix $bold(C) in RR^(n times n)$ with $"diag"(bold(C)) = bold(0)$ satisfies the *Subspace Detection Property* (is *subspace-preserving*) if for every $i in \{1, dots.c, n\}$, its column $bold(c)_i$ is nonzero ($bold(c)_i != bold(0)$) and $C_(j i) = 0$ whenever $kappa(j) != kappa(i)$. Under column permutation by $kappa$, the symmetrized affinity graph $bold(W) = |bold(C)| + |bold(C)|^top$ is strictly block-diagonal: $bold(W) = "diag"(bold(W)_1, dots.c, bold(W)_K)$.
]

Given the optimal coefficient matrix $bold(C)^*$, SSC constructs an undirected nonnegative similarity graph $bold(W) = |bold(C)^*| + |bold(C)^*|^top$ (often after column normalization $tilde(bold(c))_i = bold(c)_i^* / norm(bold(c)_i^*)_infinity$ to equalize scale across subspaces) and extracts $K$ clusters via normalized spectral clustering @ng2002spectral @vonluxburg2007tutorial @shi2000ncut.

#algbox(
  title: "Sparse Subspace Clustering (SSC via ADMM)",
  inputs: [Data matrix $bold(X) in RR^(D times n)$; cluster count $K$; penalty $lambda > 0$ (or $(lambda_e, lambda_z)$); ADMM penalty $rho > 0$; affine flag.],
  outputs: [Subspace partition $ell_1, dots.c, ell_n in \{1, dots.c, K\}$ and sparse affinity graph $bold(W) in RR_(>= 0)^(n times n)$.],
)[
  1. *Column Normalization*: Normalize nonzero columns $bold(x)_i <- bold(x)_i / norm(bold(x)_i)_2$ (for linear SSC; omit or center for affine SSC).
  2. *ADMM Self-Expression Solve*: Split $bold(C) = bold(A) - "diag"(bold(A))$ and iterate for $t = 1, dots.c, T$:
     (a) Solve the ridge system $(bold(X)^top bold(X) + rho bold(I)_n) bold(A)^((t)) = bold(X)^top bold(X) + rho bold(C)^((t-1)) - bold(Lambda)^((t-1))$ (adding Lagrange row $bold(1)^top bold(A) = bold(1)^top$ if affine).
     (b) Soft-threshold $bold(C)^((t)) <- cal(T)_(lambda / rho) (bold(A)^((t)) + bold(Lambda)^((t-1)) / rho)$ and zero its diagonal $"diag"(bold(C)^((t))) <- bold(0)$, where $cal(T)_tau (u) = "sign"(u) max(|u| - tau, 0)$.
     (c) Update dual multiplier $bold(Lambda)^((t)) <- bold(Lambda)^((t-1)) + rho (bold(A)^((t)) - bold(C)^((t)))$.
  3. *Affinity Construction*: Normalize columns $tilde(bold(c))_i <- bold(c)_i^((T)) / norm(bold(c)_i^((T)))_infinity$ (when nonzero) and set $bold(W) <- |tilde(bold(C))| + |tilde(bold(C))|^top$.
  4. *Normalized Spectral Clustering*: Compute $bold(L)_"sym" = bold(I)_n - bold(D)^(-1/2) bold(W) bold(D)^(-1/2)$ with $D_(i i) = sum_j W_(i j)$, extract the $K$ bottom eigenvectors $bold(U) in RR^(n times K)$, row-normalize $bold(U)$, and cluster rows via $k$-means @lloyd1982least.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We present five formal results governing SSC: exact subspace recovery for independent subspaces, the dual-certificate incoherence-versus-inradius theorem for overlapping subspaces, an exact LASSO regularization collapse threshold, a graph-disconnection counterexample showing that SDP alone does not guarantee exact clustering, and orthogonal / affine invariance versus column-scale vulnerability.

#subsec("A", "Exact Recovery Under Linear Independence")

Subspaces $cal(S)_1, dots.c, cal(S)_K subset.eq RR^D$ are *independent* if $dim(cal(S)_1 + dots.c + cal(S)_K) = sum_(k=1)^K dim(cal(S)_k)$, equivalently $cal(S)_k inter sum_(l != k) cal(S)_l = \{bold(0)\}$ for every $k$.

#theorem(title: "Exact Subspace Detection for Independent Subspaces", origin: "Classical — Cited")[
  Suppose $cal(S)_1, dots.c, cal(S)_K$ are independent subspaces and every column $bold(x)_i in cal(S)_k$ lies in $"span"(\{bold(x)_j in cal(S)_k : j != i\})$. Then every optimal solution $bold(C)^*$ of noiseless linear SSC satisfies the Subspace Detection Property @elhamifar2013ssc.
]

#proof[
  Fix any column $bold(x)_i in cal(S)_k$ and let $bold(c)_i^*$ be an optimal solution to $min norm(bold(c)_i)_1$ subject to $bold(X) bold(c)_i = bold(x)_i$, $c_(i i) = 0$. Partition $bold(c)_i^* = bold(c)_i^"in" + bold(c)_i^"out"$, where $bold(c)_i^"in"$ is supported exclusively on $\{j != i : bold(x)_j in cal(S)_k\}$ and $bold(c)_i^"out"$ is supported on $\{j : bold(x)_j in.not cal(S)_k\}$. Rearranging $bold(x)_i = bold(X) bold(c)_i^*$ yields
  $ bold(x)_i - bold(X) bold(c)_i^"in" = bold(X) bold(c)_i^"out". $
  The left-hand side lies in $cal(S)_k$, while the right-hand side lies in $sum_(l != k) cal(S)_l$. By subspace independence, $cal(S)_k inter sum_(l != k) cal(S)_l = \{bold(0)\}$, forcing both sides to vanish identically. Consequently $bold(x)_i = bold(X) bold(c)_i^"in"$, so $bold(c)_i^"in"$ is itself feasible for the $ell_1$ program (and nonzero because $bold(x)_i != bold(0)$). Because $bold(c)_i^"in"$ and $bold(c)_i^"out"$ have disjoint supports, the $ell_1$-norm decomposes additively:
  $ norm(bold(c)_i^*)_1 = norm(bold(c)_i^"in")_1 + norm(bold(c)_i^"out")_1 >= norm(bold(c)_i^"in")_1. $
  If $bold(c)_i^"out" != bold(0)$, the inequality is strict, contradicting optimality of $bold(c)_i^*$. Hence $bold(c)_i^"out" = bold(0)$ for every $i$, proving SDP.
]

#subsec("B", "Geometric Incoherence Bound for Non-Independent Subspaces")

When $sum_(k=1)^K d_k > D$, subspaces cannot be independent and may even share nontrivial intersections. Exact recovery still holds provided the angular separation between subspaces exceeds the inverse *inradius* of the intra-subspace data polytope. Let all columns have unit Euclidean norm ($norm(bold(x)_j)_2 = 1$). Define the first principal-angle cosine between subspaces $cal(S)_k$ and $cal(S)_l$ as
$ cos theta_1 (cal(S)_k, cal(S)_l) = max_(bold(u) in cal(S)_k, bold(w) in cal(S)_l \ norm(bold(u))_2 = norm(bold(w))_2 = 1) |chevron.l bold(u), bold(w) chevron.r|, $
and let $r_i = r(cal(P)_(-i)^((k)))$ denote the *inradius* (radius of the largest Euclidean ball in $cal(S)_k$ centered at $bold(0)$ contained inside the symmetrized convex hull $cal(P)_(-i)^((k)) = "conv"(\{plus.minus bold(x)_j : bold(x)_j in cal(S)_k, j != i\})$).

#theorem(title: "Incoherence vs. Inradius Condition for SDP", origin: "Classical — Cited")[
  Let $bold(X)_(-i)^((k))$ denote the submatrix of unit-norm columns in $cal(S)_k$ excluding $bold(x)_i$, and assume $r_i > 0$. If
  $ max_(l != k) cos theta_1 (cal(S)_k, cal(S)_l) < r_i, $
  then every optimal solution $bold(c)_i^*$ of noiseless SSC satisfies $c_(j i)^* = 0$ for all $bold(x)_j in.not cal(S)_k$ @soltanolkotabi2012geometric @elhamifar2013ssc.
]

#proof[
  Consider the intra-subspace linear program
  $ p_i^* = min_(bold(c)) norm(bold(c))_1 quad "s.t." quad bold(X)_(-i)^((k)) bold(c) = bold(x)_i, $
  which is feasible because $r_i > 0$ implies $"span"(bold(X)_(-i)^((k))) = cal(S)_k$. Its Lagrangian dual restricted to $cal(S)_k$ is
  $ d_i^* = max_(bold(v) in cal(S)_k) chevron.l bold(x)_i, bold(v) chevron.r quad "s.t." quad norm(bold(X)_(-i)^((k) top) bold(v))_infinity <= 1. $
  By strong linear programming duality, $p_i^* = d_i^*$ is attained by an optimal dual vector $bold(v)_i in cal(S)_k$. The dual constraint states that $|chevron.l bold(z), bold(v)_i chevron.r| <= 1$ for every vertex $plus.minus bold(x)_j$ of $cal(P)_(-i)^((k))$, hence for every $bold(z) in cal(P)_(-i)^((k))$ by convexity. Because $cal(P)_(-i)^((k))$ contains the $cal(S)_k$-ball of radius $r_i$, choosing $bold(z) = r_i bold(v)_i / norm(bold(v)_i)_2 in cal(P)_(-i)^((k))$ yields $r_i norm(bold(v)_i)_2 = chevron.l bold(z), bold(v)_i chevron.r <= 1$, i.e., $norm(bold(v)_i)_2 <= 1 / r_i$.

  For every external point $bold(x)_j in cal(S)_l$ ($l != k$, $norm(bold(x)_j)_2 = 1$), the definition of principal angle gives
  $ |chevron.l bold(x)_j, bold(v)_i chevron.r| <= cos theta_1 (cal(S)_k, cal(S)_l) dot.op norm(bold(v)_i)_2 <= frac(cos theta_1 (cal(S)_k, cal(S)_l), r_i) < 1. $
  Set $alpha_"max" = max_(j : bold(x)_j in.not cal(S)_k) |chevron.l bold(x)_j, bold(v)_i chevron.r| < 1$. Now take any feasible $bold(c)$ of the full program ($bold(X) bold(c) = bold(x)_i$, $c_i = 0$) and split $bold(c) = bold(c)^"in" + bold(c)^"out"$. Then
  $ p_i^* = chevron.l bold(x)_i, bold(v)_i chevron.r = sum_(j in cal(S)_k, j != i) c_j chevron.l bold(x)_j, bold(v)_i chevron.r + sum_(j in.not cal(S)_k) c_j chevron.l bold(x)_j, bold(v)_i chevron.r \ <= norm(bold(c)^"in")_1 + alpha_"max" norm(bold(c)^"out")_1 = norm(bold(c))_1 - (1 - alpha_"max") norm(bold(c)^"out")_1. $
  Because the intra-subspace optimum $bold(c)_i^((k) *)$ is feasible for the full program with objective $p_i^*$, any feasible $bold(c)$ with $bold(c)^"out" != bold(0)$ obeys $norm(bold(c))_1 >= p_i^* + (1 - alpha_"max") norm(bold(c)^"out")_1 > p_i^*$, proving that every full-program optimizer has $bold(c)^"out" = bold(0)$.
]

#subsec("C", "Regularization Phase Transition and Energy Identities")

In noisy LASSO-SSC, the regularization parameter $lambda$ governs a sharp phase transition between graph sparsity and complete vertex disconnection.

#proposition(title: "Exact LASSO-SSC Collapse Threshold and Residual Bounds", origin: "Original Derivation")[
  For column $i$, consider the LASSO-SSC objective $J_i (bold(c)) = (1 / 2) norm(bold(x)_i - bold(X)_(-i) bold(c))_2^2 + lambda norm(bold(c))_1$ with $lambda > 0$, optimal solution $bold(c)_i^*$, and residual $bold(r)_i^* = bold(x)_i - bold(X)_(-i) bold(c)_i^*$.
  1. *Exact disconnection threshold*: $bold(c)_i^* = bold(0)$ if and only if
     $ lambda >= lambda_"max" (i) := max_(j != i) |chevron.l bold(x)_j, bold(x)_i chevron.r|. $
  2. *Exact residual identity & energy bound*: Every optimizer obeys
     $ chevron.l bold(x)_i - bold(r)_i^*, bold(r)_i^* chevron.r = lambda norm(bold(c)_i^*)_1, quad norm(bold(c)_i^*)_1 <= frac(norm(bold(x)_i)_2^2 - norm(bold(r)_i^*)_2^2, 2 lambda), $
     and $norm(bold(r)_i^*)_2 < norm(bold(x)_i)_2$ whenever $0 < lambda < lambda_"max" (i)$.
]

#proof[
  Because $J_i$ is convex, $bold(c)_i^*$ is globally optimal iff the first-order KKT subdifferential condition holds:
  $ bold(X)_(-i)^top bold(r)_i^* = lambda bold(g)_i, quad bold(g)_i in partial norm(bold(c)_i^*)_1. $
  Evaluating at $bold(c)_i^* = bold(0)$ gives $bold(r)_i^* = bold(x)_i$ and $partial norm(bold(0))_1 = \{bold(g) : norm(bold(g))_infinity <= 1\}$, so $bold(0)$ is optimal iff $norm(bold(X)_(-i)^top bold(x)_i)_infinity <= lambda$, which is precisely $lambda >= lambda_"max" (i)$, proving Part 1.

  For Part 2, every subgradient $bold(g)_i in partial norm(bold(c)_i^*)_1$ satisfies $chevron.l bold(c)_i^*, bold(g)_i chevron.r = norm(bold(c)_i^*)_1$. Taking the inner product of the KKT condition with $bold(c)_i^*$ yields
  $ lambda norm(bold(c)_i^*)_1 = chevron.l bold(c)_i^*, bold(X)_(-i)^top bold(r)_i^* chevron.r = chevron.l bold(X)_(-i) bold(c)_i^*, bold(r)_i^* chevron.r = chevron.l bold(x)_i - bold(r)_i^*, bold(r)_i^* chevron.r. $
  Next, feasibility of $bold(c) = bold(0)$ implies $J_i (bold(c)_i^*) <= J_i (bold(0))$, i.e., $(1 / 2) norm(bold(r)_i^*)_2^2 + lambda norm(bold(c)_i^*)_1 <= (1 / 2) norm(bold(x)_i)_2^2$, rearranging immediately to $norm(bold(c)_i^*)_1 <= (norm(bold(x)_i)_2^2 - norm(bold(r)_i^*)_2^2) / (2 lambda)$. When $lambda < lambda_"max" (i)$, Part 1 ensures $bold(c)_i^* != bold(0)$, hence $lambda norm(bold(c)_i^*)_1 > 0$, forcing $norm(bold(r)_i^*)_2^2 <= norm(bold(x)_i)_2^2 - 2 lambda norm(bold(c)_i^*)_1 < norm(bold(x)_i)_2^2$.
]

#subsec("D", "Graph Disconnection Counterexample Under Exact SDP")

A subtle theoretical pitfall of SSC is that the Subspace Detection Property prevents *false inter-subspace edges* but does not prevent *missing intra-subspace edges*.

#proposition(title: "Spectral Over-Segmentation Despite Exact SDP", origin: "Original Counterexample")[
  Even when noiseless SSC satisfies the Subspace Detection Property with zero inter-subspace edges, the induced affinity block $bold(W)_k$ of a single $4$-dimensional subspace $cal(S)_k$ containing no collinear point pairs can split into two disconnected components, forcing $"dim ker"(bold(L)_"sym") >= K + 1$ and causing spectral clustering to over-segment $cal(S)_k$.
]

#proof[
  Let $cal(S)_1 = RR^4 = cal(U) plus.o cal(V)$, where $cal(U) = "span"(bold(e)_1, bold(e)_2)$ and $cal(V) = "span"(bold(e)_3, bold(e)_4)$ are orthogonal 2D planes. Fix an integer $m >= 3$ and place $2 m$ distinct unit-norm points in $cal(S)_1$:
  $ bold(u)_r = (cos(r pi / m), sin(r pi / m), 0, 0)^top in cal(U), $
  $ bold(v)_r = (0, 0, cos(r pi / m), sin(r pi / m))^top in cal(V), $
  for $r = 0, dots.c, m - 1$. No two points are collinear because angles lie strictly in $[0, pi)$. Because $m >= 3$, every $bold(u)_r$ lies in $"span"(\{bold(u)_s : s != r\}) = cal(U)$, and similarly for $bold(v)_r$ in $cal(V)$, while all $2 m$ points together span $cal(S)_1 = RR^4$.

  Apply Theorem 1 to the two independent (orthogonal) subspaces $cal(U)$ and $cal(V)$ inside $cal(S)_1$: every $ell_1$-optimal solution for $bold(u)_r$ places strictly zero weight on every $bold(v)_s$, and vice versa. Hence the affinity block $bold(W)_1 = "diag"(bold(W)_cal(U), bold(W)_cal(V))$ has zero cross-edges between $cal(U)$ and $cal(V)$, so the single subspace $cal(S)_1$ contributes at least two independent null eigenvectors to $bold(L)_"sym"$ @vonluxburg2007tutorial @chung1997spectral. In a $K$-subspace mixture containing $cal(S)_1$ plus $K - 1$ other independent subspaces, $"dim ker"(bold(L)_"sym") >= K + 1$, breaking exact recovery of the $K$ clusters by $K$-way spectral clustering.
]

#subsec("E", "Orthogonal Invariance, Scale Hijacking, and Translation")

#proposition(title: "Invariance, Scale Hijacking, and Affine Translation Restoration", origin: "Original Proof")[
  1. *Orthogonal invariance*: Both noiseless and LASSO-SSC produce identical optimal matrices $bold(C)^*$ and affinity graphs $bold(W)$ under any ambient orthogonal transformation $bold(X) |-> bold(Q) bold(X)$ ($bold(Q) in O(D)$).
  2. *Column-scale hijacking in LASSO-SSC*: Fix $lambda > 0$ and let $bold(r)_i^((k)) != bold(0)$ be the unique intra-subspace LASSO residual for $bold(x)_i in cal(S)_k$. Without column normalization, rescaling a single external point $bold(x)_j in cal(S)_l$ ($l != k$) with $chevron.l bold(x)_j, bold(r)_i^((k)) chevron.r != 0$ to $tilde(bold(x))_j = s bold(x)_j$ with $s > lambda / |chevron.l bold(x)_j, bold(r)_i^((k)) chevron.r|$ destroys SDP (forcing $c_(j i)^* != 0$), even when subspaces are linearly independent.
  3. *Affine translation invariance*: Translating $bold(X) |-> tilde(bold(X)) = bold(X) + bold(mu) bold(1)^top$ ($bold(mu) != bold(0)$) breaks linear SSC (driving pairwise cosines to $1$ as $norm(bold(mu))_2 -> infinity$), whereas imposing $bold(1)^top bold(C) = bold(1)^top$ restores exact translation invariance.
]

#proof[
  For Part 1, $bold(Q)^top bold(Q) = bold(I)_D$ implies $bold(Q) bold(x)_i = bold(Q) bold(X) bold(c)_i$ iff $bold(x)_i = bold(X) bold(c)_i$, and $norm(bold(Q)(bold(x)_i - bold(X) bold(c)_i))_2^2 = norm(bold(x)_i - bold(X) bold(c)_i)_2^2$, leaving both feasibility and objective unchanged.

  For Part 2, because the quadratic loss $(1 / 2) norm(bold(x)_i - bold(z))_2^2$ is strictly convex in the fitted vector $bold(z) = bold(X)_(-i)^((k)) bold(c)$, every intra-subspace LASSO-SSC optimizer shares the same unique residual $bold(r)_i^((k))$, which is nonzero since $bold(r)_i^((k)) = bold(0)$ would contradict Proposition 1 ($lambda norm(bold(c)_i^((k) *))_1 = chevron.l bold(x)_i, bold(0) chevron.r = 0 => bold(c)_i^((k) *) = bold(0) => bold(r)_i^((k)) = bold(x)_i != bold(0)$). Let $sigma = "sign"(chevron.l bold(x)_j, bold(r)_i^((k)) chevron.r) in \{-1, +1\}$ and perturb any subspace-preserving optimizer $bold(c)_i^((k) *)$ along coordinate $j$ via $bold(c)(t) = bold(c)_i^((k) *) + t sigma bold(e)_j$ ($t > 0$). Differentiating the full LASSO-SSC objective $J_i (bold(c)(t))$ at $t = 0^+$ yields
  $ frac(dif, dif t) J_i (bold(c)(t)) bar.v_(t = 0^+) = -chevron.l sigma tilde(bold(x))_j, bold(r)_i^((k)) chevron.r + lambda \ = -s |chevron.l bold(x)_j, bold(r)_i^((k)) chevron.r| + lambda < 0 $
  whenever $s > lambda / |chevron.l bold(x)_j, bold(r)_i^((k)) chevron.r|$. Thus $J_i$ strictly decreases for small $t > 0$, proving that no subspace-preserving vector can be globally optimal for the full LASSO-SSC problem and motivating unit column normalization in Step 1.

  For Part 3, under the affine constraint $bold(1)^top bold(c)_i = 1$,
  $ tilde(bold(X)) bold(c)_i - tilde(bold(x))_i = (bold(X) + bold(mu) bold(1)^top) bold(c)_i - (bold(x)_i + bold(mu)) \ = bold(X) bold(c)_i - bold(x)_i + bold(mu)(bold(1)^top bold(c)_i - 1) = bold(X) bold(c)_i - bold(x)_i. $
  Thus both the noiseless constraint and LASSO residual are invariant to $bold(mu)$, whereas without $bold(1)^top bold(c)_i = 1$, $cos angle(tilde(bold(x))_i, tilde(bold(x))_j) -> 1$ for all pairs as $norm(bold(mu))_2 -> infinity$, violating Theorem 2.
]

#sec("III", "Complexity, Geometry and Operational Profile")

*Time Complexity.* Computing the Gram matrix $bold(G) = bold(X)^top bold(X)$ costs $O(D n^2)$. In ADMM Step 2(a), solving $(bold(G) + rho bold(I)_n) bold(A) = bold(R)$ across $T$ iterations via a single Cholesky factorization costs $O(n^3 + T n^3)$ when $D >= n$. Crucially, when $D < n$, the Woodbury matrix identity
$ (bold(X)^top bold(X) + rho bold(I)_n)^(-1) = rho^(-1) bold(I)_n - rho^(-2) bold(X)^top (bold(I)_D + rho^(-1) bold(X) bold(X)^top)^(-1) bold(X) $
factorizes only a $D times D$ matrix in $O(D^3)$ and evaluates each ADMM iteration in $O(D n^2)$, reducing optimization cost to $O(T min(n, D) n^2)$. The final spectral step computes the $K$ bottom eigenvectors of $bold(L)_"sym"$ in $O(n^3)$ dense time (or $O(K "nnz"(bold(W)))$ via Lanczos when $bold(W)$ is sparse), giving overall complexity $O(T min(n, D) n^2 + n^3)$.

*Space and Numerical Profile.* Storing $bold(C), bold(A), bold(Lambda), bold(W) in RR^(n times n)$ requires $O(n^2 + n D)$ memory, becoming the primary bottleneck for $n > 10^4$. ADMM achieves $O(1 / T)$ ergodic objective convergence, while eigengap stability of Step 4 obeys the Davis-Kahan $sin Theta$ theorem @davis1970rotation.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$O(T min(n, D) n^2 + n^3)$ via ADMM + Woodbury + spectral step],
    [Space Complexity], [$O(n^2 + n D)$ for dense self-expression and affinity matrices],
    [Cluster Geometry], [Linear subspaces through origin (or affine flats with $bold(1)^top bold(C) = bold(1)^top$)],
    [Cluster Count $K$], [User-specified (or estimated from zero-eigenvalue multiplicity of $bold(L)_"sym"$)],
    [Noise / Outliers], [Explicitly modeled via $(bold(E), bold(Z))$; robust to sparse corruptions],
    [Rotation Invariance], [Strictly invariant under ambient orthogonal transforms $O(D)$],
    [Scale / Shift Sensitivity], [Requires unit column normalization; linear SSC is shift-sensitive],
    [Determinism], [Convex $ell_1$ solve is deterministic; final $k$-means is seed-dependent],
  ),
  caption: [Operational, geometric, and complexity profile of SSC @elhamifar2013ssc.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *3D Rigid-Body Motion Segmentation in Computer Vision.* Under an affine camera model, $F$-frame feature trajectories $(u_(i f), v_(i f))_(f=1)^F in RR^(2 F)$ of points lying on a single rigid moving object span an affine subspace of dimension at most $4$ (linear dimension at most $4$ in homogeneous coordinates) regardless of frame count $F$. SSC @elhamifar2013ssc recovers exact multi-body segmentations even when motions are partially degenerate or occluded, directly leveraging Theorem 1 and affine translation invariance (Proposition 3).
2. *Face Recognition Under Variable Lambertian Illumination.* By the Basri-Jacobs spherical harmonics theorem, vectorized images of a convex Lambertian face under arbitrary distant lighting conditions reside near a $9$-dimensional linear subspace in pixel space $RR^D$ ($D tilde.op 10^4$). Because distinct subjects' illumination cones have small principal-angle cosines relative to intra-subject inradii, Theorem 2 guarantees subspace-preserving sparse affinities without training labels.
3. *Hyperspectral Remote Sensing Unmixing and Land-Cover Mapping.* Pixels from a single mineral or vegetation endmember class span a low-dimensional linear spectral subspace across hundreds of narrow wavelength bands, even under local illumination scaling. SSC's $ell_1$ self-expressiveness tolerates scale variations along rays and isolates sparse sensor spikes into $bold(E)$, outperforming Euclidean distance clustering under spectral brightness drift.
4. *Hybrid Linear System Identification and Switching Dynamics.* Input-output trajectories of piecewise-linear dynamical systems switch among $K$ unknown ARX hyperplanes in lifted regressor space. Viewing each operating mode as a codimension-$1$ subspace, SSC partitions unlabelled trajectory windows into distinct regimes while handling switching transients via the sparse corruption matrix $bold(E)$.
5. *Handwritten Digit and Appearance Manifold Factorizations.* High-dimensional image patches of handwritten characters under small affine deformations cluster along low-dimensional tangent subspaces. Combining SSC's sparse local-subspace graph with normalized cuts @shi2000ncut separates intersecting digit classes where centroid-based methods @lloyd1982least fail due to overlapping means near the origin.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

SSC fails systematically in four regimes: (i) *compact, non-subspace metric clusters in low ambient dimension* ($D <= 4$ or $sum d_k > D$), where high inter-cluster cosine similarity violates Theorem 2; (ii) *uncentered affine clusters run under linear SSC*, where nonzero means collapse angular separation (Proposition 3); (iii) *mis-calibrated $lambda$*, causing isolated vertices ($lambda >= lambda_"max"$, Proposition 1) or dense false connections ($lambda -> 0$); and (iv) *large sample counts* ($n >> 10^4$), where $O(n^2)$ memory and $O(n^3)$ spectral solves become prohibitive.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Motion Segmentation], [Excellent: trajectories lie on 4D affine subspaces], [Enable affine constraint $bold(1)^top bold(C) = bold(1)^top$],
    [Face Illumination Clustering], [Excellent: 9D Lambertian subspaces in high $D$], [Normalize columns to unit $ell_2$ norm],
    [Hyperspectral Imagery], [Strong: linear spectral mixing & sparse noise $bold(E)$], [Use spatial smoothing for pixel noise],
    [Low-D Tabular / Blobs], [Poor: $sum d_k > D$ & non-subspace cluster geometry], [Use STATPC (Alg. 149) or $k$-means @lloyd1982least],
    [Curved 2D Manifolds], [Poor: global linear combinations cut across curves], [Use HDBSCAN @campello2013hdbscan or kernel spectral],
    [Web-Scale Datasets], [Poor: $O(n^2)$ affinity memory & $O(n^3)$ eigen-solve], [Use landmark Nystrom or PROCLUS @aggarwal1999proclus],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for SSC.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Within Category VI @kriegel2009subspacesurvey, SSC bridges *arbitrarily oriented projected clustering* (such as ORCLUS @aggarwal2000orclus, which alternates local PCA with $k$-medoids assignment) and *graph-theoretic spectral clustering* @ng2002spectral @vonluxburg2007tutorial by replacing Euclidean Gaussian kernels $exp(-norm(bold(x)_i - bold(x)_j)_2^2 / (2 sigma^2))$ with a data-adaptive sparse self-expressive affinity graph.

#proposition(title: "Spectral Reduction to Exact Normalized Cut Under Connected SDP", origin: "Original Derivation")[
  Suppose SSC yields a subspace-preserving matrix $bold(C)^*$ whose symmetrized affinity graph $bold(W) = |bold(C)^*| + |bold(C)^*|^top$ has connected intra-subspace blocks $bold(W)_1, dots.c, bold(W)_K$. Then the normalized cut objective @shi2000ncut achieves its global minimum $"Ncut"(cal(C)_1, dots.c, cal(C)_K) = 0$ uniquely at the ground-truth subspace partition, and the normalized Laplacian $bold(L)_"sym"$ has eigenvalue $0$ with exact multiplicity $K$ whose eigenspace is spanned by the degree-weighted indicator vectors $bold(D)^(1 / 2) bold(1)_(cal(C)_k)$.
]

#proof[
  Because $bold(W) = "diag"(bold(W)_1, dots.c, bold(W)_K)$ has zero inter-block weights, $"cut"(cal(C)_k, overline(cal(C))_k) = 0$ for every $k$, giving $"Ncut" = 0$. The normalized Laplacian likewise decouples as $bold(L)_"sym" = "diag"(bold(L)_"sym"^(1), dots.c, bold(L)_"sym"^(K))$. Because each block $bold(W)_k$ is connected with nonnegative weights and strictly positive degrees (as $bold(c)_i^* != bold(0)$ by SDP), Perron-Frobenius theory @chung1997spectral @vonluxburg2007tutorial guarantees that each $bold(L)_"sym"^(k)$ has simple eigenvalue $0$ with unique (up to scaling) eigenvector $bold(D)_k^(1 / 2) bold(1)_(n_k)$ and strictly positive Fiedler value $lambda_2 (bold(L)_"sym"^(k)) > 0$. Consequently $"dim ker"(bold(L)_"sym") = K$, row-normalizing the $K$ zero eigenvectors maps every point in $cal(S)_k$ identically to the $k$-th standard basis vector $bold(e)_k in RR^K$, and $k$-means recovers the exact partition with zero objective error.
]

*Empirical Behavior on ClusterOverview Benchmarks.* The 4-dataset ClusterOverview benchmark vividly illustrates the boundary between SSC's theoretical subspace regime and generic compact-cluster benchmarks, evaluated via Euclidean Silhouette @rousseeuw1987silhouettes and Adjusted Rand Index @hubert1985comparing:
- *Iris* (`Sil = -0.017`, `ARI = 0.010`): Near-zero ARI and slightly negative Silhouette confirm a sharp inductive-bias mismatch. With ambient dimension $D = 4$ and $K = 3$ compact ellipsoidal species displaced from the origin, pairwise cosines are large, violating Theorem 2's incoherence threshold $cos theta_1 < r_i$ and mixing points across species along shared radial cones (compare STATPC, Alg. 149: `ARI = 0.818`).
- *Customers* (`Sil = -0.078`): Positive-orthant spending vectors cluster in a narrow angular cone around the mean direction; linear self-expression links customers by proportional spending ratios rather than absolute volume, producing overlapping Euclidean clusters with negative Silhouette.
- *Multi-Shape* (`Sil = -0.011`, `ARI = N/A`): In $D = 2$, linear combinations span the entire plane with just two non-collinear points, collapsing subspace separability completely and cutting linearly across curved geometric shapes.
- *High-D* (`Sil = -0.052`, `ARI = N/A`): Even in higher ambient dimension, isotropic Gaussian or axis-aligned clusters embedded in uniform noise lack low-dimensional linear-subspace structure through the origin; combined with finite-sample $lambda$ sensitivity (Proposition 1) and potential intra-block disconnection (Proposition 2), SSC yields poor Euclidean cohesion.
- *Mean Runtime* (`5.253 ms`): Iterative convex optimization plus spectral eigendecomposition makes SSC roughly $27 times$ slower than STATPC (`0.193 ms`), underscoring that SSC should be reserved specifically for genuine union-of-subspaces domains ($D >> d_k$, angularly separated linear or affine flats) where its exact recovery theorems shine.
