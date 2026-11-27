#import "../template/ieee.typ": *

#alg_entry(
  id: 92,
  name: "Non-Negative Matrix Factorization (NMF)",
  category: "Cat. IV: Distribution, Statistical & Topic",
  year: "1999",
  paradigm: "Parts-Based Low-Rank Nonnegative Factorization",
  complexity: $cal(O)(T N D K)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Unconstrained low-rank decompositions such as PCA and SVD allow
arbitrary sign cancellations between positive and negative basis vectors,
yielding holistic components whose coordinates cannot be interpreted as
additive cluster memberships. *Non-Negative Matrix Factorization (NMF)*,
popularized by Lee and Seung @lee1999nmf @lee2001algorithms, enforces
entrywise nonnegativity on both factors so that every observation is
reconstructed strictly as an additive conical combination of $K$ latent
parts or cluster prototypes.

#definition(title: "Non-Negative Matrix Factorization and Simplicial Cone")[
  Given a nonnegative data matrix $bold(X) in bb(R)_+^(N times D)$ and
  target rank $K < min(N, D)$, *NMF* seeks a nonnegative encoding
  (soft-membership) matrix $bold(W) in bb(R)_+^(N times K)$ and a
  nonnegative basis (prototype) matrix $bold(H) in bb(R)_+^(K times D)$
  minimizing a separable Bregman divergence $D_phi (bold(X) || bold(W) bold(H))$
  @banerjee2005clustering. The two canonical objectives @lee2001algorithms
  are the *squared Frobenius loss*
  $ cal(L)_F (bold(W), bold(H)) := frac(1, 2) norm(bold(X) - bold(W) bold(H))_F^2, $
  and the *generalized Kullback-Leibler ($I$-divergence) loss*
  $ D_I (bold(X) || bold(W) bold(H)) := sum_(i=1)^N sum_(j=1)^D (&X_(i j) log frac(X_(i j), (bold(W) bold(H))_(i j)) \
    &- X_(i j) + (bold(W) bold(H))_(i j)). $
  Geometrically, every reconstructed row $hat(bold(x))_i = bold(H)^top bold(w)_i$
  lies inside the $K$-ray *simplicial cone*
  $cal(C)(bold(H)) := {bold(H)^top bold(w) : bold(w) in bb(R)_+^K} subset.eq bb(R)_+^D$.
  Normalizing prototype rows ($norm(bold(h)_k)_1 = 1$ or $norm(bold(h)_k)_2 = 1$)
  absorbs scaling into $bold(W)$, yielding soft cluster weights
  $tilde(W)_(i k) = W_(i k) / sum_(l=1)^K W_(i l)$ and hard labels
  $hat(c)_i = arg max_(1 <= k <= K) W_(i k)$.
]

Although $cal(L)_F$ and $D_I$ are individually convex in $bold(W)$ alone
and in $bold(H)$ alone, both problems are jointly non-convex and NP-hard
in general. Lee and Seung @lee2001algorithms introduced simple
*multiplicative updates (MU)* that preserve nonnegativity automatically
without line search or explicit projection: for Frobenius loss,
$ W_(i k) &<- W_(i k) frac((bold(X) bold(H)^top)_(i k), (bold(W) bold(H) bold(H)^top)_(i k)), \
  H_(k j) &<- H_(k j) frac((bold(W)^top bold(X))_(k j), (bold(W)^top bold(W) bold(H))_(k j)), $
and for generalized KL divergence,
$ W_(i k) &<- W_(i k) frac(sum_(j=1)^D H_(k j) X_(i j) / (bold(W) bold(H))_(i j), sum_(j=1)^D H_(k j)), \
  H_(k j) &<- H_(k j) frac(sum_(i=1)^N W_(i k) X_(i j) / (bold(W) bold(H))_(i j), sum_(i=1)^N W_(i k)). $

#algbox(
  title: "Multiplicative-Update NMF for Soft and Hard Clustering",
  inputs: [Nonnegative matrix $bold(X) in bb(R)_+^(N times D)$; rank $K$; divergence $cal(D) in {cal(L)_F, D_I}$; safeguard $eta > 0$; tolerance $epsilon > 0$; max iterations $T$.],
  outputs: [Factors $(bold(W), bold(H))$; normalized soft memberships $tilde(bold(W))$; hard labels $hat(c)_i = arg max_k W_(i k)$.],
)[
  + *Initialize:* Seed strictly positive $bold(W)^((0)) > 0, bold(H)^((0)) > 0$ via nonnegative SVD (NNDSVDa) or random uniform draws.
  + *For* $t = 0, 1, dots, T - 1$ *do:*
  + #h(1em) *Update $bold(H)$:* Apply the Frobenius or KL multiplicative rule to $bold(H)^((t))$ (adding $eta$ to denominators) to form $bold(H)^((t+1))$.
  + #h(1em) *Update $bold(W)$:* Using $bold(H)^((t+1))$, apply the corresponding multiplicative rule to $bold(W)^((t))$ to form $bold(W)^((t+1))$.
  + #h(1em) *Stop* if $cal(D)(bold(W)^((t)), bold(H)^((t))) - cal(D)(bold(W)^((t+1)), bold(H)^((t+1))) < epsilon$.
  + *Normalize:* Set $s_k = norm(bold(h)_k)_2$, rescale $bold(h)_k <- bold(h)_k / s_k$ and $bold(w)_(dot k) <- s_k bold(w)_(dot k)$; output $tilde(W)_(i k) = W_(i k) / sum_l W_(i l)$ and $hat(c)_i = arg max_k W_(i k)$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five core mathematical results for NMF: auxiliary-function
monotonicity of both multiplicative updates, exact objective equivalence
between KL-NMF and pLSA, an original tight spectral theorem linking
orthogonal NMF to $k$-means, a constructive geometric counterexample
proving non-uniqueness without separability, and exact scaling/zero-locking
properties.

#subsec("A", "Majorization-Minimization Monotonicity of Multiplicative Updates")

#theorem(title: "Auxiliary-Function Monotonicity for Frobenius and KL Updates", origin: "Classical — Cited")[
  For any strictly positive iterates $bold(W)^((t)) > 0, bold(H)^((t)) > 0$
  and $bold(X) in bb(R)_+^(N times D)$ with no zero rows or columns, the
  multiplicative updates satisfy
  $cal(L)_F (bold(W)^((t+1)), bold(H)^((t+1))) <= cal(L)_F (bold(W)^((t)), bold(H)^((t)))$
  and
  $D_I (bold(X) || bold(W)^((t+1)) bold(H)^((t+1))) <= D_I (bold(X) || bold(W)^((t)) bold(H)^((t)))$,
  with equality if and only if the iterate is a fixed point of the update
  map @lee2001algorithms.
]

#proof[
  By symmetry between $(bold(W), bold(H))$ and $(bold(H)^top, bold(W)^top)$,
  it suffices to verify non-increase during the $bold(H)$-update for a
  single column $bold(h) in bb(R)_+^K$ fitting column $bold(x) in bb(R)_+^N$
  of $bold(X)$ at fixed $bold(W) > 0$.

  *Case 1 (Frobenius loss):* Let $F(bold(h)) = (1 / 2) norm(bold(x) - bold(W) bold(h))_2^2$,
  with gradient $nabla F(bold(h)) = bold(W)^top bold(W) bold(h) - bold(W)^top bold(x)$
  and Hessian $bold(W)^top bold(W)$. Define the diagonal matrix
  $bold(K)(bold(h)^((t))) := "diag"(((bold(W)^top bold(W) bold(h)^((t)))_k) / h_k^((t)))$
  and quadratic surrogate
  $ G(bold(h), bold(h)^((t))) &:= F(bold(h)^((t))) + (bold(h) - bold(h)^((t)))^top nabla F(bold(h)^((t))) \
    &+ frac(1, 2) (bold(h) - bold(h)^((t)))^top bold(K)(bold(h)^((t))) (bold(h) - bold(h)^((t))). $
  Clearly $G(bold(h)^((t)), bold(h)^((t))) = F(bold(h)^((t)))$. Moreover,
  $G(bold(h), bold(h)^((t))) - F(bold(h)) = (1 / 2) (bold(h) - bold(h)^((t)))^top bold(M) (bold(h) - bold(h)^((t)))$
  where $bold(M) := bold(K)(bold(h)^((t))) - bold(W)^top bold(W)$. Letting
  $bold(A) := bold(W)^top bold(W) >= 0$ and setting
  $tilde(M)_(a b) := h_a^((t)) M_(a b) h_b^((t))$, for any
  $bold(v) in bb(R)^K$ we obtain
  $ bold(v)^top tilde(bold(M)) bold(v) &= sum_(a, b) v_a^2 h_a^((t)) A_(a b) h_b^((t)) - sum_(a, b) v_a v_b h_a^((t)) A_(a b) h_b^((t)) \
    &= frac(1, 2) sum_(a, b) A_(a b) h_a^((t)) h_b^((t)) (v_a - v_b)^2 >= 0. $
  Because $h_k^((t)) > 0$, $tilde(bold(M)) succ.eq 0$ implies
  $bold(M) succ.eq 0$, so $G(bold(h), bold(h)^((t))) >= F(bold(h))$
  everywhere. Minimizing $G(dot, bold(h)^((t)))$ sets
  $bold(h)^((t+1)) = bold(h)^((t)) - bold(K)(bold(h)^((t)))^(-1) nabla F(bold(h)^((t)))$,
  which algebraically cancels to
  $h_k^((t+1)) = h_k^((t)) (bold(W)^top bold(x))_k / (bold(W)^top bold(W) bold(h)^((t)))_k$,
  guaranteeing
  $F(bold(h)^((t+1))) <= G(bold(h)^((t+1)), bold(h)^((t))) <= G(bold(h)^((t)), bold(h)^((t))) = F(bold(h)^((t)))$.

  *Case 2 (KL divergence):* Define convex weights
  $q_(i k) := W_(i k) h_k^((t)) / sum_(l=1)^K W_(i l) h_l^((t))$
  ($sum_k q_(i k) = 1$). By convexity of $-log$ and Jensen's inequality,
  $ -log(sum_(k=1)^K W_(i k) h_k) <= -sum_(k=1)^K q_(i k) log(frac(W_(i k) h_k, q_(i k))), $
  with equality at $bold(h) = bold(h)^((t))$. Substituting this upper bound
  into $D_I (bold(x) || bold(W) bold(h))$ creates a separable convex
  auxiliary function $G_"KL" (bold(h), bold(h)^((t)))$ whose unique
  stationarity condition
  $partial G_"KL" / partial h_k = sum_(i=1)^N W_(i k) - sum_(i=1)^N X_i q_(i k) / h_k = 0$
  solves directly to
  $h_k^((t+1)) = h_k^((t)) (sum_(i=1)^N W_(i k) X_i / (bold(W) bold(h)^((t)))_i) / sum_(i=1)^N W_(i k)$,
  yielding monotonic descent.
]

#subsec("B", "Exact Equivalence between KL-NMF and pLSA")

#theorem(title: "Equivalence of KL-NMF and Probabilistic Latent Semantic Analysis", origin: "Classical — Cited")[
  Let $bold(X) in bb(R)_+^(N times D)$ be normalized so that
  $sum_(i, j) X_(i j) = 1$. Then minimizing $D_I (bold(X) || bold(W) bold(H))$
  over $bold(W) in bb(R)_+^(N times K), bold(H) in bb(R)_+^(K times D)$ is
  strictly equivalent to maximizing the pLSA log-likelihood
  $ cal(L)_"pLSA" (bold(U), bold(s), bold(V)) = sum_(i=1)^N sum_(j=1)^D X_(i j) log((bold(U) "diag"(bold(s)) bold(V))_(i j)) $
  over column-stochastic $bold(U) in bb(R)_+^(N times K)$
  ($sum_i U_(i k) = 1$), row-stochastic $bold(V) in bb(R)_+^(K times D)$
  ($sum_j V_(k j) = 1$), and $bold(s) in Delta^(K-1)$ @ding2008nmfplsa
  @hofmann1999plsa.
]

#proof[
  For any factor pair $(bold(W), bold(H))$ with positive column sums
  $d_k = sum_i W_(i k) > 0$ and row sums $e_k = sum_j H_(k j) > 0$, define
  $U_(i k) := W_(i k) / d_k$, $V_(k j) := H_(k j) / e_k$, and
  $sigma_k := d_k e_k > 0$. Then
  $bold(W) bold(H) = bold(U) "diag"(bold(sigma)) bold(V)$ and the total
  reconstructed mass is
  $S := sum_(i, j) (bold(W) bold(H))_(i j) = sum_(k=1)^K sigma_k$.
  Consider scaling $(bold(W), bold(H))$ by $c > 0$:
  $ D_I (bold(X) || c bold(W) bold(H)) &= "const" - (sum_(i, j) X_(i j)) log c \
    &- sum_(i, j) X_(i j) log((bold(W) bold(H))_(i j)) + c S. $
  Since $sum_(i, j) X_(i j) = 1$, differentiating with respect to $c$ gives
  $-1 / c + S = 0$, uniquely minimized at $c^* = 1 / S$. Therefore every
  stationary (and globally optimal) factor pair satisfies $S = 1$, whence
  $bold(s) := bold(sigma) in Delta^(K-1)$ and
  $ D_I (bold(X) || bold(W) bold(H)) = sum_(i, j) X_(i j) log X_(i j) - cal(L)_"pLSA" (bold(U), bold(s), bold(V)). $
  Because the entropy term $sum_(i, j) X_(i j) log X_(i j)$ is constant in
  $(bold(W), bold(H))$, KL-NMF and pLSA share the exact same global
  optima and stationary points.
]

#subsec("C", "Orthogonal NMF and Exact Spectral Reduction to k-Means")

A widely cited heuristic equates orthogonal NMF with $k$-means
@lloyd1982least. The following original theorem makes the exact algebraic
relationship precise, isolating the Perron-Frobenius gap between general
orthogonal NMF and binary $k$-means indicators.

#theorem(title: "Orthogonal Frobenius NMF, Hard Partitioning, and k-Means", origin: "Original Derivation")[
  Consider the orthogonal NMF problem
  $ min_(bold(W) >= 0, bold(H) >= 0, bold(W)^top bold(W) = bold(I)_K) norm(bold(X) - bold(W) bold(H))_F^2. $
  Then:
  (i) Every feasible $bold(W) in bb(R)_+^(N times K)$ has at most one
  strictly positive entry per row, inducing a hard disjoint clustering
  $C_1, dots, C_K subset.eq {1, dots, N}$ with $C_k = {i : W_(i k) > 0}$;
  (ii) For any fixed feasible $bold(W)$, the unconstrained minimizer
  $bold(H)^* = bold(W)^top bold(X)$ is automatically nonnegative and
  globally optimal;
  (iii) Letting $cal(F)_"orth" := min_(bold(W) >= 0, bold(W)^top bold(W) = bold(I)_K, bold(H) >= 0) norm(bold(X) - bold(W) bold(H))_F^2$,
  the optimal orthogonal NMF value satisfies the exact spectral identity
  and tight $k$-means bound
  $ cal(F)_"orth" &= norm(bold(X))_F^2 - max_(C_1, dots, C_K "disj.") sum_(k=1)^K lambda_"max" ((bold(X) bold(X)^top)_(C_k, C_k)) \
    &<= min_(C_1, dots, C_K) sum_(k=1)^K sum_(i in C_k) norm(bold(x)_i - bold(mu)_k)_2^2 = "SSE"_(k"-means")^*, $
  with exact equality to $"SSE"_(k"-means")^*$ whenever $bold(W)$ is
  further restricted to cluster-indicator columns
  $W_(i k) = bb(I)(i in C_k) / sqrt(abs(C_k))$.
]

#proof[
  For (i), since $bold(W) >= 0$ and $(bold(W)^top bold(W))_(k l) = sum_(i=1)^N W_(i k) W_(i l) = 0$
  for $k != l$, every summand $W_(i k) W_(i l)$ must vanish. Thus no row
  $i$ can possess two positive entries, so the supports $C_k = {i : W_(i k) > 0}$
  are pairwise disjoint (and nonempty because $norm(bold(w)_(dot k))_2 = 1$).

  For (ii), using $bold(W)^top bold(W) = bold(I)_K$, expand the Frobenius
  norm orthogonally:
  $ norm(bold(X) - bold(W) bold(H))_F^2 &= "Tr"(bold(X)^top bold(X)) - 2 "Tr"(bold(H)^top bold(W)^top bold(X)) \
    &+ "Tr"(bold(H)^top bold(H)) \
    &= norm(bold(X))_F^2 - norm(bold(W)^top bold(X))_F^2 + norm(bold(H) - bold(W)^top bold(X))_F^2. $
  Because $bold(W) >= 0$ and $bold(X) >= 0$, the matrix
  $bold(H)^* = bold(W)^top bold(X)$ satisfies $bold(H)^* >= 0$ and
  annihilates the final squared norm uniquely.

  For (iii), substituting $bold(H)^*$ leaves
  $norm(bold(X))_F^2 - sum_(k=1)^K bold(w)_(dot k)^top (bold(X) bold(X)^top)_(C_k, C_k) bold(w)_(dot k)$.
  Each principal submatrix $(bold(X) bold(X)^top)_(C_k, C_k)$ is symmetric
  and entrywise nonnegative, so by the Perron-Frobenius theorem its
  dominant eigenvalue $lambda_"max" ((bold(X) bold(X)^top)_(C_k, C_k))$ is
  attained by a unit entrywise-nonnegative eigenvector $bold(w)_(dot k) >= 0$.
  Maximizing over all disjoint supports yields the exact equality in (iii).
  Finally, choosing the constant unit vector
  $bold(w)_(dot k) = bold(1)_(C_k) / sqrt(abs(C_k))$ evaluates the
  Rayleigh quotient to $abs(C_k) norm(bold(mu)_k)_2^2$, where
  $bold(mu)_k = abs(C_k)^(-1) sum_(i in C_k) bold(x)_i$, so
  $ norm(bold(X))_F^2 - sum_(k=1)^K abs(C_k) norm(bold(mu)_k)_2^2 = sum_(k=1)^K sum_(i in C_k) norm(bold(x)_i - bold(mu)_k)_2^2. $
  Since $lambda_"max" ((bold(X) bold(X)^top)_(C_k, C_k)) >= abs(C_k) norm(bold(mu)_k)_2^2$,
  orthogonal NMF is a tight rank-1-per-cluster relaxation of $k$-means
  that preserves hard disjoint clusters.
]

#subsec("D", "Simplicial-Cone Non-Uniqueness and Separability")

#proposition(title: "Constructive Non-Uniqueness Counterexample and Separable Rigidity", origin: "Original Counterexample")[
  (i) Exact nonnegative factorization $bold(X) = bold(W) bold(H)$ is
  generically non-unique even modulo scaling and permutation: for the
  exact rank-$3$ matrix $bold(X) = bold(W) bold(I)_3 in bb(R)_+^(4 times 3)$
  with
  $ bold(W) = mat(2, 1, 1; 1, 2, 1; 1, 1, 2; 2, 2, 1), $
  there exists a continuous family $(bold(W)(epsilon), bold(H)(epsilon))$
  for all $epsilon in (0, 1)$ of non-monomial-equivalent nonnegative
  factors satisfying $bold(W)(epsilon) bold(H)(epsilon) = bold(X)$.
  (ii) Conversely, if both $bold(W)$ and $bold(H)^top$ are *separable*
  (each contains a positive diagonal submatrix after row/column
  permutation, i.e., pure anchor rows and anchor features), then every
  exact rank-$K$ factorization $(bold(W) bold(Q), bold(Q)^(-1) bold(H))$
  must have $bold(Q)$ a generalized permutation (monomial) matrix.
]

#proof[
  For (i), set $bold(H)(epsilon) := mat(1, epsilon, 0; 0, 1, epsilon; epsilon, 0, 1) >= 0$
  for $epsilon in (0, 1)$. Its determinant is $1 + epsilon^3 > 0$, with
  exact inverse
  $ bold(Q)(epsilon) = bold(H)(epsilon)^(-1) = frac(1, 1 + epsilon^3) mat(1, -epsilon, epsilon^2; epsilon^2, 1, -epsilon; -epsilon, epsilon^2, 1). $
  Define $bold(W)(epsilon) := bold(W) bold(Q)(epsilon)$. Every row
  $(w_1, w_2, w_3)$ of $bold(W)$ satisfies $1 <= w_r <= 2$, so every entry
  of $(1 + epsilon^3) bold(W)(epsilon)$ has the form
  $w_a - epsilon w_b + epsilon^2 w_c >= 1 - 2 epsilon + epsilon^2 = (1 - epsilon)^2 > 0$.
  Hence $bold(W)(epsilon) > 0$, $bold(H)(epsilon) >= 0$, and
  $bold(W)(epsilon) bold(H)(epsilon) = bold(W) bold(I)_3 = bold(X)$, yet
  $bold(H)(epsilon)$ is not a scaled permutation for any $epsilon in (0, 1)$,
  continuously distorting the soft cluster weights $tilde(bold(W))(epsilon)$.

  For (ii), suppose $bold(W)' = bold(W) bold(Q) >= 0$ and
  $bold(H)' = bold(Q)^(-1) bold(H) >= 0$ with invertible $bold(Q) in bb(R)^(K times K)$.
  Because $bold(W)$ contains a positive diagonal block $bold(D)_1$ as rows,
  $bold(D)_1 bold(Q) >= 0 ==> bold(Q) >= 0$; symmetrically, because
  $bold(H)$ contains a positive diagonal block $bold(D)_2$ as columns,
  $bold(Q)^(-1) bold(D)_2 >= 0 ==> bold(Q)^(-1) >= 0$. If row $k$ of
  $bold(Q)$ had two positive entries $Q_(k l_1), Q_(k l_2) > 0$, then
  $(bold(Q) bold(Q)^(-1))_(k j) = 0$ for all $j != k$ would force
  $(bold(Q)^(-1))_(l_1 j) = (bold(Q)^(-1))_(l_2 j) = 0$ for all $j != k$,
  making rows $l_1$ and $l_2$ of $bold(Q)^(-1)$ linearly dependent
  multiples of $bold(e)_k^top$, contradicting invertibility. Thus $bold(Q)$
  has exactly one positive entry per row and column, proving uniqueness up
  to permutation and positive scaling.
]

#subsec("E", "Scaling Equivariance, Zero-Locking, and Translation Bias")

#proposition(title: "Scaling Homogeneity, MU Zero-Locking, and Translation Bias", origin: "Original Proof")[
  (i) Global positive rescaling $bold(X) |-> c bold(X)$ ($c > 0$) preserves
  the normalized soft memberships $tilde(bold(W))$ identically for both
  Frobenius and KL NMF, whereas non-uniform column scaling
  $bold(X) |-> bold(X) "diag"(bold(s))$ reweights feature $j$ by $s_j^2$
  in $cal(L)_F$ and by $s_j$ in $D_I$.
  (ii) In multiplicative updates, any zero entry $W_(i k)^((t)) = 0$ (or
  $H_(k j)^((t)) = 0$) remains locked at $0$ for all $s >= t$, allowing
  convergence to non-KKT fixed points where
  $nabla_(W_(i k)) cal(L) < 0$.
  (iii) Uniform translation $bold(X) |-> bold(X) + c bold(1) bold(1)^top$
  ($c > 0$) is *not* invariant: as $c -> infinity$, the optimal rank-$K$
  factorization spends one cone ray approximating $c bold(1) bold(1)^top$
  and compresses angular discrimination among rows.
]

#proof[
  Part (i) follows directly from
  $cal(L)_F (c^(1/2) bold(W), c^(1/2) bold(H)) = c^2 cal(L)_F (bold(W), bold(H))$
  and $D_I (c bold(X) || c bold(W) bold(H)) = c D_I (bold(X) || bold(W) bold(H))$,
  while substituting $bold(H)' = bold(H) "diag"(bold(s))$ gives
  $cal(L)_F (bold(W), bold(H)') = (1 / 2) sum_(i, j) s_j^2 (X_(i j) - (bold(W) bold(H))_(i j))^2$
  and $D_I (bold(X) "diag"(bold(s)) || bold(W) bold(H)') = sum_(i, j) s_j d_I (X_(i j) || (bold(W) bold(H))_(i j))$.
  For (ii), the multiplicative factor scales $W_(i k)^((t))$ by a ratio, so
  $W_(i k)^((t)) = 0 ==> W_(i k)^((s)) = 0$ for all $s >= t$ even when the
  KKT dual feasibility condition
  $nabla_(W_(i k)) cal(L)_F = (bold(W) bold(H) bold(H)^top - bold(X) bold(H)^top)_(i k) >= 0$
  is violated (e.g., $(bold(X) bold(H)^top)_(i k) > (bold(W) bold(H) bold(H)^top)_(i k)$).
  For (iii), for any two rows $bold(x)_i, bold(x)_l$, the shifted cosine
  similarity approaches
  $(bold(x)_i + c bold(1))^top (bold(x)_l + c bold(1)) / (norm(bold(x)_i + c bold(1))_2 norm(bold(x)_l + c bold(1))_2) = 1 - cal(O)(c^(-2)) -> 1$
  as $c -> infinity$, collapsing all rows onto the diagonal ray
  $bb(R)_+ bold(1)$ and inflating the cone condition number.
]

#sec("III", "Complexity, Geometry and Operational Profile")

For Frobenius MU, computing $bold(X) bold(H)^top$ costs $cal(O)(N D K)$
dense (or $cal(O)(K "nnz"(bold(X)))$ sparse), while computing
$bold(W) (bold(H) bold(H)^top)$ in optimal associative order costs only
$cal(O)((N + D) K^2)$ rather than forming the $N times D$ product
$bold(W) bold(H)$. Thus each iteration executes in
$cal(O)(K "nnz"(bold(X)) + (N + D) K^2)$ time ($cal(O)(T N D K)$ dense
over $T$ sweeps). In contrast, KL-MU evaluates the ratio
$X_(i j) / (bold(W) bold(H))_(i j)$ strictly on the nonzero support of
$bold(X)$ (requiring $cal(O)(K "nnz"(bold(X)))$), matching Frobenius MU
asymptotically while exhibiting lower BLAS-3 cache reuse on dense inputs.
Memory footprint is $cal(O)("nnz"(bold(X)) + (N + D) K)$. Faster modern
solvers such as Hierarchical Alternating Least Squares (HALS) and
projected-gradient coordinate descent retain the same per-sweep flop
count while eliminating MU zero-locking (Proposition 2) and accelerating
local KKT convergence.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [*Property*], [*Characterization*],
    [Time complexity], [$cal(O)(T N D K)$ dense; $cal(O)(T K "nnz"(bold(X)))$ sparse],
    [Space complexity], [$cal(O)(N D + (N + D) K)$ dense; $cal(O)("nnz"(bold(X)) + (N + D) K)$],
    [Inductive geometry], [$K$-ray nonnegative simplicial cone $cal(C)(bold(H)) subset.eq bb(R)_+^D$],
    [Clustering type], [Additive soft membership $tilde(bold(W))$; hard via $arg max$ or $bold(W)^top bold(W) = bold(I)$],
    [Identifiability], [Unique modulo scaling/permutation under separability (Prop. 1)],
    [Invariance], [Global scale homogeneous; NOT translation or rotation invariant],
    [Optimization trap], [Non-convex NP-hard; MU can lock zeros at non-KKT points],
    [Initialization], [NNDSVDa (deterministic) or multi-start positive random seeds],
  ),
  caption: [Operational and geometric characterization of NMF clustering.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

+ *Document Clustering and Interpretable Topic Discovery:* Nonnegative
  term-frequency or TF-IDF matrices satisfy $bold(X) >= 0$ natively, and
  Theorem 2 proves that KL-NMF recovers the exact pLSA topic-document
  likelihood @ding2008nmfplsa @hofmann1999plsa. Unlike signed LSI/SVD
  axes, every basis row $bold(h)_k >= 0$ lists coherent positive topic
  keywords and every row $bold(w)_i >= 0$ quantifies additive multi-topic
  document membership.
+ *Transcriptomics, Metagenes, and Cancer Subtyping:* Gene-expression
  counts and single-cell UMI matrices are intrinsically nonnegative
  measurements of additive biological pathways @kluger2003spectral.
  Factoring $bold(X) approx bold(W) bold(H)$ simultaneously discovers
  co-regulated *metagene* signatures ($bold(H)$) and soft tumor-subtype
  assignments ($bold(W)$), naturally capturing mixed-lineage biopsies that
  hard partitioning misclassifies.
+ *Audio Source Separation and Spectrogram Factorization:* Magnitude
  short-time Fourier spectrograms $bold(X) in bb(R)_+^(N times D)$ combine
  harmonic instruments additively over time frames. NMF (especially with
  KL or Itakura-Saito divergences) decomposes polyphonic mixtures into
  spectral templates ($bold(H)$) and temporal activation envelopes ($bold(W)$)
  without destructive phase cancellation.
+ *Overlapping Community Detection in Networks:* Factoring symmetric
  nonnegative adjacency matrices $bold(A) approx bold(W) bold(W)^top$
  (Symmetric NMF and its log-likelihood cousin BIGCLAM @yang2013bigclam)
  identifies overlapping social and biological communities where nodes
  participate additively across multiple functional circles.
+ *Hyperspectral Unmixing and Parts-Based Vision:* Reflectance pixels in
  remote sensing are convex combinations of pure mineral endmembers
  ($bold(H)$) weighted by fractional ground abundances ($bold(W)$)
  @lee1999nmf. Here the physical presence of pure pixels satisfies the
  separability condition of Proposition 1, rendering endmember recovery
  mathematically unique.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Four structural failure modes govern NMF in clustering pipelines:
(1) *Signed or shifted Euclidean data:* General real-valued features
must be shifted ($bold(X) - min(bold(X))$) or split into positive/negative
parts before applying NMF; Proposition 2(iii) proves that large shifts
collapse angular ray diversity and degrade cluster separability.
(2) *Identifiability drift without anchor points:* When data sit strictly
inside the interior of $bb(R)_+^D$ without pure-cluster boundary points,
cone ambiguity (Proposition 1) causes different initializations to produce
equally low reconstruction loss with divergent soft memberships.
(3) *MU zero-locking and slow tail convergence:* Initializing with exact
zeros or relying solely on MU can stall at non-KKT boundary points
(Proposition 2(ii)); HALS or NNDSVDa warm-starts mitigate this failure.
(4) *Uncalibrated probability and missing $K$ selection:* Unlike HDP
(Alg. 91) or BIC-LCA (Alg. 90), Frobenius NMF lacks a normalized
predictive density, requiring cophenetic consensus stability or
cross-validated reconstruction error to select $K$.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [*Domain*], [*Mathematical Fit*], [*Caution / Alternative*],
    [Term-count / TF-IDF corpora], [Exact KL-NMF / pLSA equivalence; parts-based], [Use HDP (Alg. 91) / LDA if Bayesian prior needed],
    [Gene expression / spectrograms], [Native $bb(R)_+$; additive metagene / spectral parts], [Use NNDSVDa + HALS to avoid MU zero-locking],
    [Hyperspectral endmembers], [Simplicial cone matches physical linear mixing], [Check anchor separability; add sparsity penalty],
    [Signed Euclidean tabular data], [Shift breaks geometry (Prop. 2(iii))], [Use GMM, vMF (Alg. 89), or $k$-means],
    [Nominal categorical surveys], [One-hot lacks probabilistic local independence], [Use LCA (Alg. 90) for exact discrete mixture],
  ),
  caption: [Domain suitability and failure-mode matrix for NMF.],
)

#sec("V", "Structural Relations and Empirical Behavior")

NMF serves as a unifying algebraic bridge across Cat. IV and matrix/graph
clustering @xu2005survey @jain2010data:
- *Probabilistic topic bridge:* Theorem 2 establishes exact equivalence
  between KL-NMF and pLSA @ding2008nmfplsa @hofmann1999plsa, while adding
  Dirichlet priors over rows of $bold(W)$ and $bold(H)$ bridges directly
  to LDA @blei2003lda and nonparametric HDP (Alg. 91).
- *Partitional and spectral bridge:* Theorem 3 proves that orthogonal
  Frobenius NMF strictly enforces hard disjoint clustering and tightly
  relaxes $k$-means @lloyd1982least, while symmetric NMF
  $norm(bold(S) - bold(W) bold(W)^top)_F^2$ on normalized affinity matrices
  provides a nonnegative relaxation of kernel $k$-means and spectral
  clustering @dhillon2004kernel @cai2013multiview.
- *Biclustering bridge:* Because $bold(W)$ groups rows and $bold(H)$
  simultaneously groups columns via shared rank-$K$ latent factors,
  sparse NMF directly overlaps factor-analytic biclustering methods such
  as FABIA @hochreiter2010fabia and spectral biclustering
  @kluger2003spectral @cheng2000biclustering.

In the ClusterOverview 4-dataset empirical benchmark (evaluated with
Silhouette @rousseeuw1987silhouettes and ARI @hubert1985comparing), NMF
stands out as the *fastest algorithm in the Algs. 89--92 block* (mean
runtime *2.052 ms*, driven by cache-efficient dense matrix products),
while exhibiting performance closely aligned with its simplicial-cone
inductive bias:
- *Multi-Shape:* Silhouette $= 0.710$ ($"ARI" =$ N/A), slightly edging
  movMF ($0.709$) and LCA ($0.709$) as distinct spatial lobes occupy
  well-separated rays after nonnegative alignment.
- *Iris:* Silhouette $= 0.296$, $"ARI" = 0.512$. Because all four raw
  Iris botanical measurements are strictly positive and strongly
  correlated in overall organ size, rows cluster near a narrow interior
  cone; Proposition 1's cone ambiguity and the positive baseline offset
  compress angular separation between _Versicolor_ and _Virginica_,
  trailed behind directional movMF ($"ARI" = 0.577$) and discretized LCA
  ($"ARI" = 0.638$).
- *Customers:* Silhouette $= 0.266$. Wholesale spending channels are
  naturally nonnegative, letting NMF extract interpretable additive spend
  archetypes, though heavy-tailed outliers pull Frobenius prototypes more
  strongly than unit-normalized movMF ($0.351$) or nonparametric HDP
  ($0.380$).
- *High-D:* Silhouette $= 0.287$ ($"ARI" =$ N/A). Low-rank $cal(O)(N D K)$
  factorization stays fast and numerically well-conditioned in high
  dimensions, though translating signed Gaussian features into $bb(R)_+^D$
  incurs the shift compression proved in Proposition 2(iii).

Across scientific and engineering workloads, NMF remains the premier
interpretable, linear-algebraic soft-clustering engine whenever data are
intrinsically nonnegative and generated by additive superposition of
latent parts @lee1999nmf @lee2001algorithms @ding2008nmfplsa.
