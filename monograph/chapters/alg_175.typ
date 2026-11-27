#import "../template/ieee.typ": *

#alg_entry(
  id: 175,
  name: "Co-regularized Multi-View Spectral Clustering",
  category: "Cat. VII: Fuzzy, Biclustering & Domain",
  year: "2011",
  paradigm: "Spectral Multi-View Subspace Co-Regularization",
  complexity: $cal(O)(T V n^3)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

When multi-view data exhibit non-convex manifolds, sparse relational links, or heterogeneous non-Euclidean similarities, Euclidean prototype methods such as Multi-View K-Means (Alg. 174) fail to separate intertwined clusters. *Co-regularized Multi-View Spectral Clustering* @kumar2011coregularized lifts each view $v in {1, ..., V}$ onto a normalized graph Laplacian and couples the resulting spectral subspaces across views on the Grassmann manifold so that view-specific eigenvectors agree on a shared underlying cluster structure.

For each view $v = 1, ..., V$, construct a symmetric non-negative affinity matrix $bold(W)_v in bb(R)_(>= 0)^(n times n)$ (for example via a Gaussian kernel $(bold(W)_v)_(i j) = exp(-norm(bold(x)_(v, i) - bold(x)_(v, j))_2^2 / (2 sigma_v^2))$), degree matrix $bold(D)_v = "diag"(bold(W)_v bold(1)_n)$, and normalized similarity operator
$
  bold(S)_v = bold(D)_v^(-1 / 2) bold(W)_v bold(D)_v^(-1 / 2)
  = bold(I)_n - bold(L)_v,
$
where $bold(L)_v$ is the symmetric normalized graph Laplacian @chung1997spectral @vonluxburg2007tutorial. Single-view normalized spectral clustering @ng2002spectral @shi2000ncut maximizes $"Tr"(bold(U)_v^top bold(S)_v bold(U)_v)$ over the Stiefel manifold $"St"(n, k) = {bold(U) in bb(R)^(n times k) : bold(U)^top bold(U) = bold(I)_k}$. Because any right rotation $bold(U)_v |-> bold(U)_v bold(Q)$ ($bold(Q) in cal(O)(k)$) preserves the cluster subspace, penalizing $norm(bold(U)_v - bold(U)_w)_F^2$ directly is basis-dependent. Instead, Kumar et al. @kumar2011coregularized measure disagreement through the rank-$k$ orthogonal projection matrices $bold(P)_v = bold(U)_v bold(U)_v^top$.

#definition(title: "Grassmannian Co-Regularized Spectral Objectives")[
  Let $bold(P)_v = bold(U)_v bold(U)_v^top$ for $bold(U)_v in "St"(n, k)$ and define the squared chordal Grassmannian disagreement
  $
    D_"Gr" (bold(U)_v, bold(U)_w)
    = (1 / 2) norm(bold(P)_v - bold(P)_w)_F^2
    = k - "Tr"(bold(P)_v bold(P)_w).
  $
  Dropping additive constants, the *pairwise* co-regularized objective with symmetric couplings $lambda_(v w) = lambda_(w v) >= 0$ is
  $
    cal(J)_"pw" ({bold(U)_v})
    = sum_(v=1)^V "Tr"(bold(U)_v^top bold(S)_v bold(U)_v)
      + sum_(1 <= v < w <= V) lambda_(v w) "Tr"(bold(P)_v bold(P)_w),
  $
  and the *centroid-based* objective introducing a consensus basis $bold(U)_star in "St"(n, k)$ ($bold(P)_star = bold(U)_star bold(U)_star^top$) with weights $lambda_v > 0$ is
  $
    cal(J)_"cen" ({bold(U)_v}, bold(U)_star)
    = sum_(v=1)^V "Tr"(bold(U)_v^top bold(S)_v bold(U)_v)
      + sum_(v=1)^V lambda_v "Tr"(bold(P)_v bold(P)_star).
  $
]

Both formulations jointly maximize within-view normalized-cut spectral coherence and cross-view Hilbert-Schmidt similarity $"Tr"(bold(P)_v bold(P)_w)$ between the linear embedding kernels $bold(K)_(bold(U)_v) = bold(U)_v bold(U)_v^top$ @kumar2011coregularized.

#algbox(
  title: "Alternating Co-Regularized Multi-View Spectral Clustering",
  inputs: [Normalized similarity matrices ${bold(S)_v}_(v=1)^V$, clusters $k >= 2$, coupling weights ${lambda_v}$ (centroid) or ${lambda_(v w)}$ (pairwise), tolerance $tau > 0$, max sweeps $T$.],
  outputs: [Consensus spectral embedding $bold(U)_star in "St"(n, k)$ and partition $bold(G) in cal(G)_(n, k)$.],
)[
  + Initialize each $bold(U)_v^((0))$ with the top-$k$ orthonormal eigenvectors of $bold(S)_v$; in the centroid scheme, initialize $bold(U)_star^((0))$ from the top-$k$ eigenvectors of $sum_(v=1)^V lambda_v bold(U)_v^((0)) (bold(U)_v^((0)))^top$.
  + *For* $t = 1, ..., T$:
    - *Centroid scheme*:
      (a) For each $v = 1, ..., V$, set $bold(U)_v^((t))$ to the top-$k$ orthonormal eigenvectors of $tilde(bold(S))_v^((t)) = bold(S)_v + lambda_v bold(P)_star^((t-1))$.
      (b) Form $bold(M)_star^((t)) = sum_(v=1)^V lambda_v bold(P)_v^((t))$ and set $bold(U)_star^((t))$ to its top-$k$ orthonormal eigenvectors.
    - *Pairwise scheme (Gauss-Seidel)*: For $v = 1, ..., V$, set $bold(U)_v^((t))$ to the top-$k$ eigenvectors of $bold(S)_v + sum_(w != v) lambda_(v w) bold(P)_w$ using the latest $bold(P)_w$.
    - Stop if $cal(J)^((t)) - cal(J)^((t-1)) <= tau$.
  + Row-normalize $bold(U)_star^((t))$ (or concatenated / selected view embeddings in the pairwise scheme) to unit $ell_2$ norm @ng2002spectral and run $k$-means @lloyd1982least to obtain $bold(G)$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish principal-angle and Procrustes equivalence of the co-regularization penalty, monotone Ky Fan ascent with an exact eigengap-weighted Grassmannian progress bound, exact strong-coupling asymptotic reduction to summed-Laplacian spectral clustering, Davis-Kahan eigengap amplification under view noise, and a sharp phase-transition counterexample under nuisance-view contamination.

#subsec("A", "Grassmannian Invariance and Principal-Angle Geometry")

#lemma(title: "Orthogonal Invariance and Tight Procrustes Equivalence", origin: "Original Proof")[
  Let $bold(U)_v, bold(U)_w in "St"(n, k)$ span subspaces $cal(U)_v, cal(U)_w subset.eq bb(R)^n$ with canonical principal angles $0 <= theta_1 <= ... <= theta_k <= pi / 2$. Then:
  + $bold(P)_v = bold(U)_v bold(U)_v^top$, $cal(J)_"pw"$, and $cal(J)_"cen"$ are strictly invariant under independent right rotations $bold(U)_v |-> bold(U)_v bold(Q)_v$ with $bold(Q)_v in cal(O)(k)$.
  + The disagreement satisfies the exact principal-angle identity
    $
      D_"Gr" (bold(U)_v, bold(U)_w)
      = k - "Tr"(bold(P)_v bold(P)_w)
      = sum_(r=1)^k sin^2 theta_r in [0, k],
    $
    vanishing if and only if $cal(U)_v = cal(U)_w$.
  + Letting $d_"Proc"^2 (bold(U)_v, bold(U)_w) = (1 / 2) min_(bold(Q) in cal(O)(k)) norm(bold(U)_v - bold(U)_w bold(Q))_F^2$, the two subspace distances obey the tight two-sided equivalence
    $
      d_"Proc"^2 (bold(U)_v, bold(U)_w)
      <= D_"Gr" (bold(U)_v, bold(U)_w)
      <= 2 d_"Proc"^2 (bold(U)_v, bold(U)_w).
    $
]

#proof[
  For any $bold(Q)_v in cal(O)(k)$, $(bold(U)_v bold(Q)_v) (bold(U)_v bold(Q)_v)^top = bold(U)_v bold(Q)_v bold(Q)_v^top bold(U)_v^top = bold(P)_v$. Since $"Tr"(bold(U)_v^top bold(S)_v bold(U)_v) = "Tr"(bold(P)_v bold(S)_v)$, both objectives depend on $bold(U)_v$ solely through $bold(P)_v$, proving part 1.

  For part 2, let $bold(U)_v^top bold(U)_w = bold(A) bold(Sigma) bold(B)^top$ be the singular value decomposition with singular values $sigma_r = cos theta_r in [0, 1]$. Expanding the Frobenius norm and using $bold(P)_v^2 = bold(P)_v$, $"Tr"(bold(P)_v) = k$ gives
  $
    D_"Gr" (bold(U)_v, bold(U)_w)
    &= (1 / 2) "Tr"(bold(P)_v^2 + bold(P)_w^2 - 2 bold(P)_v bold(P)_w) \
    &= k - norm(bold(U)_v^top bold(U)_w)_F^2
    = sum_(r=1)^k (1 - cos^2 theta_r)
    = sum_(r=1)^k sin^2 theta_r.
  $
  Thus $D_"Gr" = 0$ iff $cos theta_r = 1$ for all $r$, i.e., $cal(U)_v = cal(U)_w$.

  For part 3, expanding $norm(bold(U)_v - bold(U)_w bold(Q))_F^2 = 2 k - 2 "Tr"(bold(Q)^top bold(U)_w^top bold(U)_v)$ and maximizing over $bold(Q) in cal(O)(k)$ via von Neumann's trace inequality gives the Procrustes optimum $bold(Q)^star = bold(B) bold(A)^top$ and $d_"Proc"^2 = sum_(r=1)^k (1 - cos theta_r)$. Factoring $sin^2 theta_r = (1 - cos theta_r)(1 + cos theta_r)$ and noting $1 <= 1 + cos theta_r <= 2$ for every $theta_r in [0, pi / 2]$ establishes $d_"Proc"^2 <= D_"Gr" <= 2 d_"Proc"^2$, attained at $theta_r = pi / 2$ (left) and as $theta_r -> 0$ (right).
]

#subsec("B", "Monotone Ky Fan Ascent and Quantitative Progress Bound")

#theorem(title: "Monotone Ky Fan Ascent and Grassmannian Convergence", origin: "Original Proof")[
  Let $mu_1 (bold(A)) >= ... >= mu_n (bold(A))$ denote the eigenvalues of a symmetric matrix $bold(A)$. Along both the pairwise and centroid alternating eigensolvers:
  + Every block update globally maximizes its subproblem via the Ky Fan principle, making $cal(J)^((t))$ monotonically non-decreasing and uniformly bounded above by $k V + k sum_(v < w) lambda_(v w)$ (pairwise) or $k V + k sum_v lambda_v$ (centroid), hence convergent to a finite limit $cal(J)_infinity$.
  + In the centroid scheme, if $bold(M)_star^((t)) = sum_(v=1)^V lambda_v bold(P)_v^((t))$ has $k$-th eigengap $delta_star^((t)) = mu_k (bold(M)_star^((t))) - mu_(k+1) (bold(M)_star^((t))) > 0$, then the consensus update obeys the exact quantitative gain
    $
      cal(J)_"cen" ({bold(U)_v^((t))}, bold(U)_star^((t)))
      &- cal(J)_"cen" ({bold(U)_v^((t))}, bold(U)_star^((t-1))) \
      &>= delta_star^((t)) D_"Gr" (bold(U)_star^((t)), bold(U)_star^((t-1))).
    $
    Consequently, if $inf_t delta_star^((t)) >= delta > 0$, then $sum_(t=1)^infinity D_"Gr" (bold(U)_star^((t)), bold(U)_star^((t-1))) < infinity$ and $D_"Gr" (bold(U)_star^((t)), bold(U)_star^((t-1))) -> 0$.
]

#proof[
  Holding all blocks except one fixed, every subproblem has the canonical form $max_(bold(U) in "St"(n, k)) "Tr"(bold(U)^top bold(A) bold(U))$ for a symmetric operator $bold(A)$ (namely $bold(S)_v + lambda_v bold(P)_star$, $bold(M)_star$, or $bold(S)_v + sum_(w != v) lambda_(v w) bold(P)_w$). By Ky Fan's variational theorem, its global maximum is $sum_(r=1)^k mu_r (bold(A))$, attained by the top-$k$ orthonormal eigenvectors. Hence $cal(J)^((t)) >= cal(J)^((t-1))$. Because every normalized similarity matrix satisfies $-bold(I)_n prec.eq bold(S)_v prec.eq bold(I)_n$ @chung1997spectral, we have $"Tr"(bold(P)_v bold(S)_v) <= k$ and $"Tr"(bold(P)_v bold(P)_w) <= k$, yielding the stated uniform upper bounds and existence of $cal(J)_infinity = lim_(t -> infinity) cal(J)^((t))$.

  For part 2, let $bold(M)_star^((t)) = sum_(r=1)^n mu_r bold(q)_r bold(q)_r^top$ be the spectral decomposition with $bold(P)_+ = bold(P)_star^((t)) = sum_(r=1)^k bold(q)_r bold(q)_r^top$ and $bold(P) = bold(P)_star^((t-1))$. Since $0 prec.eq bold(P) prec.eq bold(I)_n$ and $"Tr"(bold(P)) = k$, the quantities $a_r = bold(q)_r^top bold(P) bold(q)_r in [0, 1]$ satisfy $sum_(r=1)^k a_r = "Tr"(bold(P)_+ bold(P))$ and
  $
    sum_(r=1)^k (1 - a_r)
    = sum_(r=k+1)^n a_r
    = k - "Tr"(bold(P)_+ bold(P))
    = D_"Gr" (bold(U)_star^((t)), bold(U)_star^((t-1))).
  $
  Using $mu_r >= mu_k$ for $r <= k$ and $mu_r <= mu_(k+1)$ for $r >= k + 1$:
  $
    "Tr"((bold(P)_+ - bold(P)) bold(M)_star^((t)))
    &= sum_(r=1)^k mu_r (1 - a_r) - sum_(r=k+1)^n mu_r a_r \
    &>= mu_k sum_(r=1)^k (1 - a_r) - mu_(k+1) sum_(r=k+1)^n a_r \
    &= (mu_k - mu_(k+1)) D_"Gr" (bold(U)_star^((t)), bold(U)_star^((t-1))).
  $
  Telescoping over $t = 1, ..., T$ when $delta_star^((t)) >= delta > 0$ gives $delta sum_(t=1)^T D_"Gr" (bold(U)_star^((t)), bold(U)_star^((t-1))) <= cal(J)_infinity - cal(J)^((0)) < infinity$, proving $D_"Gr" (bold(U)_star^((t)), bold(U)_star^((t-1))) -> 0$.
]

#subsec("C", "Strong-Coupling Asymptotic Reduction to Summed-Laplacian Fusion")

#theorem(title: "Dimension-Free Strong-Coupling Limit", origin: "Original Derivation")[
  Set $lambda_v = lambda w_v$ with fixed positive view weights $w_v > 0$ and coupling strength $lambda > 0$. Define the shifted optimal value $cal(F)(lambda) = max_(bold(U)_v, bold(U)_star) cal(J)_"cen" - lambda k sum_(v=1)^V w_v$ and the fused-graph optimum $Phi_"sum" = max_(bold(U)_star in "St"(n, k)) "Tr"(bold(U)_star^top (sum_(v=1)^V bold(S)_v) bold(U)_star)$. Then every global maximizer $({bold(U)_v (lambda)}, bold(U)_star (lambda))$ satisfies:
  + $sum_(v=1)^V w_v D_"Gr" (bold(U)_v (lambda), bold(U)_star (lambda)) <= (2 k V) / lambda$;
  + $Phi_"sum" <= cal(F)(lambda) <= Phi_"sum" + 2 k sqrt(2 V (sum_(v=1)^V w_v^(-1)) / lambda)$,
    so $lim_(lambda -> infinity) cal(F)(lambda) = Phi_"sum"$ at an $n$-independent $cal(O)(lambda^(-1 / 2))$ rate.
]

#proof[
  Rewrite the shifted objective identically as
  $
    cal(F)(lambda)
    = max_(bold(U)_v, bold(U)_star)
      [sum_(v=1)^V "Tr"(bold(P)_v bold(S)_v)
        - lambda sum_(v=1)^V w_v D_"Gr" (bold(U)_v, bold(U)_star)].
  $
  Evaluating the feasible consensus configuration $bold(U)_1 = ... = bold(U)_V = bold(U)_star$ sets all $D_"Gr" = 0$, proving the lower bound $cal(F)(lambda) >= Phi_"sum"$. Because $-bold(I)_n prec.eq bold(S)_v prec.eq bold(I)_n$, we have $-k V <= Phi_"sum"$ and $sum_v "Tr"(bold(P)_v (lambda) bold(S)_v) <= k V$; substituting into $Phi_"sum" <= cal(F)(lambda)$ yields
  $
    lambda sum_(v=1)^V w_v D_"Gr" (bold(U)_v (lambda), bold(U)_star (lambda))
    <= k V - Phi_"sum" <= 2 k V,
  $
  which proves part 1.

  For part 2, decompose $"Tr"(bold(P)_v (lambda) bold(S)_v) = "Tr"(bold(P)_star (lambda) bold(S)_v) + "Tr"((bold(P)_v (lambda) - bold(P)_star (lambda)) bold(S)_v)$. Since $"rank"(bold(P)_v - bold(P)_star) <= 2 k$ and spectral norm $norm(bold(S)_v)_2 <= 1$, Schatten-$1$ / spectral duality (with nuclear norm $norm(dot)_(S_1)$) and Cauchy-Schwarz on singular values give
  $
    |"Tr"((bold(P)_v - bold(P)_star) bold(S)_v)|
    &<= norm(bold(P)_v - bold(P)_star)_(S_1) norm(bold(S)_v)_2 \
    &<= sqrt(2 k) norm(bold(P)_v - bold(P)_star)_F
    = 2 sqrt(k D_"Gr" (bold(U)_v, bold(U)_star)).
  $
  Dropping the non-positive penalty term and applying Cauchy-Schwarz across $v = 1, ..., V$:
  $
    cal(F)(lambda)
    &<= "Tr"(bold(P)_star (lambda) sum_(v=1)^V bold(S)_v)
      + 2 sqrt(k) sum_(v=1)^V w_v^(-1 / 2) (w_v D_"Gr")^(1 / 2) \
    &<= Phi_"sum"
      + 2 sqrt(k) (sum_(v=1)^V w_v^(-1))^(1 / 2)
        ((2 k V) / lambda)^(1 / 2),
  $
  establishing the $n$-independent $cal(O)(lambda^(-1 / 2))$ convergence to summed normalized-Laplacian spectral clustering.
]

#subsec("D", "Davis-Kahan Eigengap Amplification Under View Noise")

#theorem(title: "Co-Regularized Subspace Perturbation Bound", origin: "Classical — Cited")[
  Let $bold(S)_v$ be a clean normalized similarity matrix with top-$k$ eigenspace $bold(U)_v$ aligned to consensus ($bold(P)_v = bold(P)_star$) and eigengap $delta_v = mu_k (bold(S)_v) - mu_(k+1) (bold(S)_v) > 0$. Suppose view $v$ is observed with symmetric perturbation $hat(bold(S))_v = bold(S)_v + bold(E)_v$ satisfying $norm(bold(E)_v)_2 < delta_v + lambda_v$. Then the perturbed update $hat(bold(U))_v$ (top-$k$ eigenvectors of $hat(bold(S))_v + lambda_v bold(P)_star$) obeys the Davis-Kahan bound @davis1970rotation @vonluxburg2007tutorial
  $
    sqrt(D_"Gr" (hat(bold(U))_v, bold(U)_v))
    <= (min(norm(bold(E)_v)_F, sqrt(k) norm(bold(E)_v)_2))
      / (delta_v + lambda_v - norm(bold(E)_v)_2).
  $
]

#proof[
  Set $bold(A) = bold(S)_v + lambda_v bold(P)_star$ and $hat(bold(A)) = bold(A) + bold(E)_v$. Because $bold(P)_star = bold(U)_v bold(U)_v^top$ commutes with $bold(S)_v$, adding $lambda_v bold(P)_star$ shifts all top-$k$ eigenvalues of $bold(S)_v$ upward by $lambda_v$ while leaving the bottom $n - k$ eigenvalues unchanged, widening the eigengap to $tilde(delta)_v = mu_k (bold(A)) - mu_(k+1) (bold(A)) = delta_v + lambda_v$. Let $hat(bold(U))_(v, perp) in "St"(n, n - k)$ span the orthogonal complement of $hat(bold(U))_v$. Then $norm(hat(bold(U))_(v, perp)^top bold(U)_v)_F^2 = "Tr"(bold(U)_v^top (bold(I)_n - hat(bold(P))_v) bold(U)_v) = D_"Gr" (hat(bold(U))_v, bold(U)_v)$. By Weyl's inequality, $mu_(k+1) (hat(bold(A))) <= mu_(k+1) (bold(A)) + norm(bold(E)_v)_2$, so the spectral separation satisfies $mu_k (bold(A)) - mu_(k+1) (hat(bold(A))) >= delta_v + lambda_v - norm(bold(E)_v)_2 > 0$. Applying the Davis-Kahan $sin Theta$ theorem @davis1970rotation to $bold(A)$ and $hat(bold(A))$ yields
  $
    norm(hat(bold(U))_(v, perp)^top bold(U)_v)_F
    <= (norm(bold(E)_v bold(U)_v)_F)
      / (mu_k (bold(A)) - mu_(k+1) (hat(bold(A)))),
  $
  and $norm(bold(E)_v bold(U)_v)_F <= min(norm(bold(E)_v)_F, sqrt(k) norm(bold(E)_v)_2)$ completes the proof. Thus an accurate consensus $bold(P)_star$ amplifies view $v$'s effective eigengap from $delta_v$ to $delta_v + lambda_v$.
]

#subsec("E", "Phase-Transition Failure Under Nuisance-View Contamination")

#proposition(title: "Exact Hijacking Phase Transition Under Fixed Coupling", origin: "Original Counterexample")[
  Consider $V = 2$ views on $n = 4 m$ vertices ($k = 2$) with balanced true partition $A | B$ and orthogonal nuisance partition $C | D$ ($|A inter C| = |A inter D| = m$). Let $bold(q)_0 = n^(-1 / 2) bold(1)_n$, $bold(q)_1 = n^(-1 / 2) (bold(1)_A - bold(1)_B)$, $bold(q)_2 = n^(-1 / 2) (bold(1)_C - bold(1)_D)$, and define noisy true-view operator $bold(S)_1 = bold(q)_0 bold(q)_0^top + (1 - eta) bold(q)_1 bold(q)_1^top$ ($eta in (0, 1)$) alongside crisp nuisance-view operator $bold(S)_2 = bold(q)_0 bold(q)_0^top + bold(q)_2 bold(q)_2^top$. Under equal-weight centroid co-regularization ($lambda_1 = lambda_2 = lambda > 0$):
  + The true consensus $bold(U)_star = [bold(q)_0, bold(q)_1]$ is *never* globally optimal for any $lambda > 0$.
  + There is a sharp phase transition at $lambda_c = 2 (1 - eta) / eta$: for every $lambda >= lambda_c$, the unique global maximizer collapses completely onto the nuisance subspace $bold(U)_1 = bold(U)_2 = bold(U)_star = [bold(q)_0, bold(q)_2]$, yielding `"ARI" = 0` against $A | B$.
]

#proof[
  Since $bold(q)_0$ is a common top eigenvector of $bold(S)_1$ and $bold(S)_2$, every global maximizer aligns its first basis vector with $bold(q)_0$ and places its second orthonormal vectors $bold(u)_1, bold(u)_2, bold(u)_star$ inside $"span"(bold(q)_1, bold(q)_2)$. Parameterize $bold(u)_star = sqrt(c) bold(q)_1 + sqrt(1 - c) bold(q)_2$ with $c in [0, 1]$, and let $a = 1 - eta in (0, 1)$. Maximizing out $bold(u)_1$ and $bold(u)_2$ for fixed $c$ yields the reduced objective $2 + 2 lambda + g(c)$, where $g(c) = mu_1 (bold(M)_1 (c)) + mu_1 (bold(M)_2 (c))$ for the $2 times 2$ matrices $bold(M)_1 (c) = a bold(e)_1 bold(e)_1^top + lambda bold(u)_star bold(u)_star^top$ and $bold(M)_2 (c) = bold(e)_2 bold(e)_2^top + lambda bold(u)_star bold(u)_star^top$. Computing traces and determinants ($det(bold(M)_1) = a lambda (1 - c)$, $det(bold(M)_2) = lambda c$) gives
  $
    g(c)
    &= (a + 1 + 2 lambda) / 2
      + (1 / 2) sqrt((a - lambda)^2 + 4 a lambda c) \
    &quad + (1 / 2) sqrt((1 + lambda)^2 - 4 lambda c).
  $
  Both square-root terms have strictly negative second derivatives in $c$, so $g$ is strictly concave on $[0, 1]$. For $lambda != 1$, the left derivative at $c = 1$ is
  $
    g'(1^-) = (a lambda) / (a + lambda) - lambda / abs(1 - lambda) < 0
  $
  for all $a in (0, 1)$ (since $a / (a + lambda) < 1$ and $a (lambda - 1) < a + lambda$ when $lambda > 1$), while at $lambda = 1$ we have $g'(1^-) = -infinity$, proving $c = 1$ is never optimal for any $lambda > 0$. At $c = 0$ (for $lambda > a$),
  $
    g'(0)
    = (a lambda) / (lambda - a) - lambda / (1 + lambda)
    = (lambda (2 a - lambda (1 - a))) / ((lambda - a) (1 + lambda)).
  $
  Substituting $a = 1 - eta$ shows $g'(0) <= 0$ iff $lambda >= 2 (1 - eta) / eta = lambda_c$ (which automatically exceeds $1 - eta$). By strict concavity of $g$, $c = 0$ is the unique global maximizer whenever $lambda >= lambda_c$, forcing $bold(U)_1 = bold(U)_2 = bold(U)_star = [bold(q)_0, bold(q)_2]$ and zero overlap with $A | B$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Dense affinity construction requires $cal(O)(n^2 D)$ time and $cal(O)(V n^2)$ memory across $V$ views. In each alternating sweep, exact dense eigendecomposition of $V$ operators $tilde(bold(S))_v in bb(R)^(n times n)$ costs $cal(O)(V n^3)$, giving total dense runtime $cal(O)(T V n^3)$—cubic in sample size $n$. Two algebraic optimizations dramatically improve scalability on sparse graphs: (i) in the centroid update, $bold(M)_star = sum_(v=1)^V lambda_v bold(U)_v bold(U)_v^top = bold(B) bold(B)^top$ with $bold(B) = [sqrt(lambda_1) bold(U)_1, ..., sqrt(lambda_V) bold(U)_V] in bb(R)^(n times k V)$ has rank at most $k V$, so its top-$k$ eigenvectors follow from thin SVD of $bold(B)$ in $cal(O)(n k^2 V^2 + k^3 V^3)$ time; (ii) replacing dense eigensolvers by Lanczos / Arnoldi iterations exploits the low-rank action $(bold(S)_v + lambda_v bold(P)_star) bold(z) = bold(S)_v bold(z) + lambda_v bold(U)_star (bold(U)_star^top bold(z))$, reducing per-view cost to $cal(O)(k ("nnz"(bold(S)_v) + n k))$ per Krylov cycle without ever forming $bold(P)_star$.

Geometrically, co-regularization performs block-coordinate optimization over a product of Grassmannians $"Gr"(n, k)^V$, fusing nonlinear manifold embeddings that relax multi-view normalized graph cuts @shi2000ncut before spherical $k$-means discretization on the unit hypersphere $bb(S)^(k - 1)$ @ng2002spectral.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$cal(O)(T V n^3)$ dense; $cal(O)(T V k ("nnz" + n k))$ via Lanczos],
    [Space Complexity], [$cal(O)(V n^2)$ dense; $cal(O)(sum_v "nnz"(bold(S)_v) + V n k)$ sparse],
    [Decision Geometry], [Non-convex manifold cuts via Grassmannian alignment],
    [Hyperparameters], [Clusters $k$, view kernels $sigma_v$, couplings $lambda_v$ or $lambda_(v w)$],
    [Invariance], [Right-$cal(O)(k)$ invariant; graph-isomorphism invariant],
    [Noise Behavior], [Widens eigengap by $+lambda_v$; vulnerable to fixed-$lambda$ nuisance view],
    [Convergence], [Monotone Ky Fan ascent; $D_"Gr" -> 0$ under positive eigengap],
    [Bottleneck], [Repeated $n times n$ affinity construction & eigendecompositions],
  ),
  caption: [Operational and geometric profile of Co-regularized Multi-View Spectral Clustering (Alg. 175).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

+ *Single-Cell Multi-Omics Manifold Integration*: Joint scRNA-seq and scATAC-seq assays capture cellular differentiation along curved, branching trajectories. Converting each modality into a $k$-NN affinity graph and aligning their Grassmannian subspaces resolves non-convex developmental lineages that Euclidean prototypes cannot separate.
+ *Multilingual and Hyperlinked Document Clustering*: News and scholarly corpora pair lexical cosine similarities with citation or hyperlink graphs @kumar2011coregularized. Normalized Laplacians naturally handle sparse graph neighborhoods and varying degree distributions while co-regularization transfers community signal from dense links to sparse text views.
+ *Multi-Descriptor Computer Vision*: Images represented via color, texture, geometric shape, and deep neural embeddings lie on nonlinear appearance manifolds; projection-kernel alignment fuses complementary invariance properties without requiring feature-space comparability.
+ *Multilayer Social and Biological Network Community Detection*: Protein-protein interaction and co-expression graphs, or multi-platform social graphs, supply native view-wise adjacency matrices $bold(W)_v$ directly; co-regularization recovers consensus modules while allowing bounded layer-specific deviations.
+ *Audio-Visual Speaker Diarization*: Aligning short-time acoustic similarity graphs with facial-track similarity graphs allows visual continuity to bridge acoustic noise and audio continuity to bridge visual occlusion via the eigengap amplification of Theorem 4.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Four core failure regimes arise in practice. First, *cubic sample scaling* ($cal(O)(T V n^3)$) restricts dense implementations to $n <= 10^4$ unless sparse $k$-NN graphs and Lanczos / Nystrom approximations are adopted. Second, *high-dimensional affinity concentration* @beyer1999nearest flattens isotropic Gaussian kernels $bold(W)_v$, shrinking eigengaps $delta_v -> 0$ and triggering Davis-Kahan instability (Theorem 4). Third, *fixed-$lambda$ nuisance-view hijacking* (Proposition 1) lets a crisp but irrelevant view override weaker true signal whenever $lambda >= 2 (1 - eta) / eta$, since Alg. 175 does not adaptively zero out bad views. Fourth, *graph disconnectedness and bandwidth misspecification* ($sigma_v$ too small or too large) introduce spurious unit eigenvalues into $bold(S)_v$, corrupting the top-$k$ subspace before co-regularization begins.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Single-Cell Multi-Omics], [Captures curved lineages & $k$-NN graphs], [Use sparse Lanczos for $n > 10^4$],
    [Multilayer Networks], [Native normalized-cut graph fusion], [Tune $lambda_v$ against noisy layers],
    [Non-Convex Shapes], [Separates rings, arcs, & manifolds], [Sensitive to RBF bandwidth $sigma_v$],
    [High-D Tabular Data], [Poor fit: RBF distance concentration], [Use Alg. 174 (MVKM) or cosine $k$-NN],
    [Large-Scale Web Data], [Poor dense scaling ($cal(O)(T V n^3)$)], [Use Alg. 174 or Nystrom spectral],
  ),
  caption: [Domain suitability and failure-mitigation matrix for Alg. 175.],
)

#sec("V", "Structural Relations and Empirical Behavior")

#proposition(title: "Spectral Reductions to NJW, Graph-Sum Fusion, and Kernel K-Means", origin: "Original Derivation")[
  Co-regularized Multi-View Spectral Clustering connects four foundational paradigms:
  + Setting $lambda_v = 0$ (or $lambda_(v w) = 0$) decouples both $cal(J)_"pw"$ and $cal(J)_"cen"$ into $V$ independent instances of Ng-Jordan-Weiss (NJW) normalized spectral clustering @ng2002spectral @vonluxburg2007tutorial, relaxing Normalized Cut @shi2000ncut (and Ratio Cut @hagen1992ratiocut for unnormalized Laplacians).
  + Letting $lambda -> infinity$ with $lambda_v = lambda w_v$ converges exactly to NJW spectral clustering on the additive normalized similarity matrix $sum_(v=1)^V bold(S)_v$ (Theorem 3).
  + Viewing $bold(P)_v = bold(U)_v bold(U)_v^top$ as a rank-$k$ consensus co-association kernel links the centroid update $bold(M)_star = sum_v lambda_v bold(P)_v$ to spectral consensus clustering @strehl2002cluster @monti2003consensus, weighted kernel $k$-means @dhillon2004kernel, and nonlinearly lifted Multi-View K-Means (Alg. 174) @cai2013multiview.
]

#proof[
  At $lambda_v = 0$, the cross-view trace penalties vanish identically, leaving $V$ decoupled Rayleigh-Ritz trace maximizations $max_(bold(U)_v in "St"(n, k)) "Tr"(bold(U)_v^top bold(S)_v bold(U)_v)$, which is exact NJW spectral clustering @ng2002spectral. Part 2 follows immediately from Theorem 3. For part 3, observe that $bold(M)_star$ is a positive linear combination of positive semidefinite projection kernels $bold(P)_v$, so extracting its top-$k$ eigenvectors solves the Fréchet chordal mean problem $min_(bold(U)_star in "St"(n, k)) sum_v lambda_v D_"Gr" (bold(U)_v, bold(U)_star)$, which parallels spectral ensemble consensus @strehl2002cluster and the trace equivalence with weighted kernel $k$-means @dhillon2004kernel.
]

Empirically, on the ClusterOverview four-dataset benchmark, Co-regularized Multi-View Spectral Clustering achieves Silhouette `"Sil" = 0.185` @rousseeuw1987silhouettes and `"ARI" = 0.408` @hubert1985comparing on Iris, `"Sil" = 0.032` on Customers, `"Sil" = 0.479` (`"ARI" = "N/A"`) on Multi-Shape, and `"Sil" = -0.025` (`"ARI" = "N/A"`) on High-D, with a mean runtime of `18.015 ms`. These metrics directly mirror our theory: the method performs best on Multi-Shape (`"Sil" = 0.479`), where graph affinities resolve non-convex cluster structure, and maintains reasonable external accuracy on Iris (`"ARI" = 0.408`, close to Alg. 174's `0.450`) despite lower Euclidean Silhouette (`0.185` vs. `0.476`) because normalized cuts optimize graph conductance rather than ambient Euclidean variance. Conversely, on Customers (`0.032`) and High-D (`-0.025`), isotropic RBF distance concentration @beyer1999nearest collapses within-view eigengaps $delta_v$, activating the Davis-Kahan perturbation regime (Theorem 4) and fixed-coupling degradation (Proposition 1). Finally, its `18.015 ms` runtime is $43.3 times$ slower than Alg. 174 (`0.416 ms`), quantitatively confirming the $cal(O)(T V n^3)$ eigendecomposition overhead relative to linear-time MVKM.
