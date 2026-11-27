#import "../template/ieee.typ": *

#alg_entry(
  id: 90,
  name: "Latent Class Analysis (LCA)",
  category: "Cat. IV: Distribution, Statistical & Topic",
  year: "1968",
  paradigm: "Finite Categorical Mixture under Local Independence",
  complexity: $cal(O)(T K (N J + M))$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Categorical and survey observations --- such as binary clinical symptoms,
nominal questionnaire responses, or discretized tabular attributes --- lack
an intrinsic Euclidean metric. Applying distance-based clustering to
one-hot encodings conflates arbitrary category cardinality with geometric
separation. *Latent Class Analysis (LCA)*, formalized by Lazarsfeld and
Henry @lazarsfeld1968latent and placed within finite mixture estimation
via EM @dempster1977em @mclachlan2000finite, models multivariate
categorical data as a finite mixture of discrete product distributions
governed by an unobserved nominal class variable $Z in {1, dots, K}$.

#definition(title: "Latent Class Model and Local Independence Tensor")[
  Let $bold(X) = (X_1, dots, X_J)$ be a vector of $J >= 2$ manifest
  categorical variables, where $X_j$ takes values in ${1, dots, M_j}$
  ($M_j >= 2$), $M := sum_(j=1)^J M_j$, and $S := product_(j=1)^J M_j$
  is the contingency-table cell count. A $K$-class *Latent Class Model*
  posits a discrete latent variable $Z in {1, dots, K}$ with mixing
  weights $bold(pi) in Delta^(K-1)$ ($"Pr"(Z = k) = pi_k$) and
  class-conditional item probabilities
  $bold(theta)_(k j) = (theta_(k j 1), dots, theta_(k j M_j)) in Delta^(M_j - 1)$
  such that $X_1, dots, X_J$ are *locally independent* given $Z$:
  $ "Pr"(bold(X) = bold(x) mid(|) Z = k; bold(Theta)) = product_(j=1)^J product_(m=1)^(M_j) theta_(k j m)^(bb(I)(x_j = m)). $
  Marginalizing $Z$ produces a rank-$K$ nonnegative CP decomposition of
  the $J$-way joint probability tensor $cal(P) in [0, 1]^(M_1 times dots.c times M_J)$:
  $ cal(P)_(m_1, dots, m_J) = "Pr"(bold(X) = bold(m) mid(|) bold(Theta)) = sum_(k=1)^K pi_k product_(j=1)^J theta_(k j m_j). $
]

The parameter tuple $bold(Theta) = (bold(pi), {bold(theta)_(k j)})$
contains $nu_K = (K - 1) + K sum_(j=1)^J (M_j - 1) = K(M - J + 1) - 1$
free scalar parameters. Given i.i.d. observations
$cal(X) = {bold(x)_1, dots, bold(x)_N}$, maximum-likelihood estimation
maximizes the incomplete-data log-likelihood
$ cal(L)(bold(Theta)) = sum_(i=1)^N log (sum_(k=1)^K pi_k product_(j=1)^J theta_(k j x_(i j))). $
Notice that if the empirical dataset contains only $U <= min(N, S)$
distinct response patterns with counts $n_u$, the log-likelihood
compresses identically to a weighted sum over those $U$ unique patterns.
To prevent degenerate zero-probability cells in finite samples, Bayesian
MAP-LCA optionally adds a symmetric Dirichlet prior $"Dir"(alpha, dots, alpha)$
($alpha >= 0$, with $alpha = 0$ recovering exact MLE) to each item vector
$bold(theta)_(k j)$ @biernacki2006mixmod.

#algbox(
  title: "Expectation-Maximization for Latent Class Analysis (LCA)",
  inputs: [Categorical records $bold(x)_1, dots, bold(x)_N$ over $J$ items with cardinalities $M_1, dots, M_J$; class count $K$; Dirichlet smoothing $alpha >= 0$; tolerance $epsilon > 0$; max iterations $T$.],
  outputs: [Class priors $bold(pi)$; item-response tables ${bold(theta)_(k j)}$; posterior responsibilities $bold(Gamma) in [0, 1]^(N times K)$; hard labels $hat(c)_i = arg max_k gamma_(i k)$.],
)[
  + *Initialize:* Draw $bold(pi)^((0)) in "int"(Delta^(K-1))$ and strictly positive $bold(theta)_(k j)^((0)) in "int"(Delta^(M_j - 1))$ across multiple random restarts (or $k$-modes seeding).
  + *For* $t = 0, 1, dots, T - 1$ *do:*
  + #h(1em) *E-step:* For each $i in {1, dots, N}$ and $k in {1, dots, K}$, evaluate log-joint probability
    $ ell_(i k) = log pi_k^((t)) + sum_(j=1)^J log theta_(k j x_(i j))^((t)), $
    and compute $gamma_(i k)^((t)) = exp(ell_(i k)) / sum_(r=1)^K exp(ell_(i r))$ via log-sum-exp.
  + #h(1em) *M-step:* For each $k in {1, dots, K}$, compute $N_k = sum_(i=1)^N gamma_(i k)^((t))$, set $pi_k^((t+1)) = N_k / N$, and for each item $j$ and level $m in {1, dots, M_j}$ update
    $ theta_(k j m)^((t+1)) = frac(alpha + sum_(i=1)^N gamma_(i k)^((t)) bb(I)(x_(i j) = m), M_j alpha + N_k). $
  + #h(1em) *Stop* if $cal(L)(bold(Theta)^((t+1))) - cal(L)(bold(Theta)^((t))) < epsilon$.
  + *Return* best restart maximizing $cal(L)(bold(Theta))$ (or selecting $K$ via $"BIC"(K) = -2 cal(L)(hat(bold(Theta))_K) + nu_K log N$).
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five foundational results for LCA: an exact three-term KL
ascent identity giving explicit non-asymptotic stationarity rates, a
constructive spectral proof of generic identifiability for $J >= 3$ items,
an explicit continuum counterexample showing non-identifiability at
$J = 2$, nominal relabeling invariance alongside EM boundary trapping, and
asymptotic consistency of BIC class enumeration.

#subsec("A", "Exact Three-Term KL Ascent Identity and Step-Size Rate")

#theorem(title: "Exact KL Ascent Decomposition and Convergence Rate of LCA-EM", origin: "Original Derivation")[
  Let ${bold(Theta)^((t))}_(t >= 0)$ be generated by unregularized LCA-EM
  ($alpha = 0$) with $N_k^((t)) > 0$ for all $k$. Then for every $t >= 0$:
  (i) the M-step updates are the unique global maximizers of the expected
  complete-data log-likelihood $Q(bold(Theta) mid(|) bold(Theta)^((t)))$;
  (ii) the incomplete-data log-likelihood increment satisfies the exact
  three-term Kullback-Leibler identity
  $ cal(L)(bold(Theta)^((t+1))) - cal(L)(bold(Theta)^((t))) &= N D_"KL" (bold(pi)^((t+1)) || bold(pi)^((t))) \
    &+ sum_(k=1)^K N_k^((t)) sum_(j=1)^J D_"KL" (bold(theta)_(k j)^((t+1)) || bold(theta)_(k j)^((t))) \
    &+ sum_(i=1)^N D_"KL" (bold(gamma)_i^((t)) || bold(gamma)_i^((t+1))) >= 0; $
  (iii) for any horizon $T >= 1$, the minimum weighted $ell_1$ step size
  decays at rate
  $ min_(0 <= t < T) (norm(bold(pi)^((t+1)) - bold(pi)^((t)))_1^2 + sum_(k, j) pi_k^((t+1)) norm(bold(theta)_(k j)^((t+1)) - bold(theta)_(k j)^((t)))_1^2) \
    <= frac(-2 cal(L)(bold(Theta)^((0))), N T). $
]

#proof[
  Write $N_k^((t)) = sum_(i=1)^N gamma_(i k)^((t))$ and
  $N_(k j m)^((t)) = sum_(i=1)^N gamma_(i k)^((t)) bb(I)(x_(i j) = m)$.
  The expected complete-data log-likelihood is
  $ Q(bold(Theta) mid(|) bold(Theta)^((t))) &= sum_(k=1)^K N_k^((t)) log pi_k \
    &+ sum_(k=1)^K sum_(j=1)^J sum_(m=1)^(M_j) N_(k j m)^((t)) log theta_(k j m). $
  Maximizing this strictly concave separable function over
  $Delta^(K-1) times product_(k, j) Delta^(M_j - 1)$ via Lagrange
  multipliers enforces $pi_k = N_k^((t)) / N$ and
  $theta_(k j m) = N_(k j m)^((t)) / N_k^((t))$, proving (i).

  Substituting $N_k^((t)) = N pi_k^((t+1))$ and
  $N_(k j m)^((t)) = N_k^((t)) theta_(k j m)^((t+1))$ into the $Q$-gain
  yields
  $ Q(bold(Theta)^((t+1)) mid(|) bold(Theta)^((t))) - Q(bold(Theta)^((t)) mid(|) bold(Theta)^((t))) &= N sum_(k=1)^K pi_k^((t+1)) log frac(pi_k^((t+1)), pi_k^((t))) \
    &+ sum_(k=1)^K N_k^((t)) sum_(j=1)^J sum_(m=1)^(M_j) theta_(k j m)^((t+1)) log frac(theta_(k j m)^((t+1)), theta_(k j m)^((t))), $
  which equals the first two terms of (ii). By the fundamental EM identity
  @dempster1977em,
  $cal(L)(bold(Theta)) = Q(bold(Theta) mid(|) bold(Theta)^((t))) - H(bold(Theta) mid(|) bold(Theta)^((t)))$
  where $H(bold(Theta) mid(|) bold(Theta)^((t))) = sum_(i, k) gamma_(i k)^((t)) log "Pr"(Z_i = k mid(|) bold(x)_i; bold(Theta))$.
  Hence
  $ H(bold(Theta)^((t)) mid(|) bold(Theta)^((t))) - H(bold(Theta)^((t+1)) mid(|) bold(Theta)^((t))) = sum_(i=1)^N D_"KL" (bold(gamma)_i^((t)) || bold(gamma)_i^((t+1))), $
  proving the exact equality (ii).

  For (iii), since each likelihood factor is a discrete probability in
  $[0, 1]$, we have $cal(L)(bold(Theta)^((T))) <= 0$. Telescoping (ii)
  from $t = 0$ to $T - 1$ and dropping the nonnegative posterior KL term
  gives
  $ sum_(t=0)^(T-1) (&D_"KL" (bold(pi)^((t+1)) || bold(pi)^((t))) \
    &+ sum_(k, j) pi_k^((t+1)) D_"KL" (bold(theta)_(k j)^((t+1)) || bold(theta)_(k j)^((t)))) <= frac(-cal(L)(bold(Theta)^((0))), N). $
  Applying Pinsker's inequality $norm(p - q)_1^2 <= 2 D_"KL" (p || q)$ term
  by term and bounding the sum below by $T$ times its minimum yields (iii).
]

By Gibbs' inequality, equality in Theorem 1(ii) holds if and only if
$bold(Theta)^((t+1)) = bold(Theta)^((t))$, guaranteeing strict monotonic
ascent at every non-stationary iterate and matching Wu's general
convergence theory @wu1983convergence.

#subsec("B", "Constructive Algebraic Identifiability for Three or More Items")

#theorem(title: "Spectral Tripartition Identifiability of Latent Class Models", origin: "Original Proof")[
  Suppose the $J >= 3$ manifest variables admit a partition into three
  nonempty disjoint blocks $A, B, C$ with combined state counts
  $S_A = product_(j in A) M_j$, $S_B = product_(j in B) M_j$, and
  $S_C = product_(j in C) M_j$. Let
  $bold(M)_A in [0, 1]^(S_A times K)$,
  $bold(M)_B in [0, 1]^(S_B times K)$, and
  $bold(M)_C in [0, 1]^(S_C times K)$ denote the block-conditional
  probability matrices (with columns given by Kronecker products of item
  vectors in each block). If $pi_k > 0$ for all $k$,
  $"rank"(bold(M)_A) = "rank"(bold(M)_B) = K$ (necessitating
  $min(S_A, S_B) >= K$), and the columns of $bold(M)_C$ are pairwise
  distinct, then the parameter tuple $(bold(pi), {bold(theta)_(k j)})$ is
  uniquely identifiable from the manifest distribution $cal(P)$ up to a
  global permutation of the $K$ latent classes.
]

#proof[
  By local independence, for every block-$C$ outcome $c in {1, dots, S_C}$,
  the corresponding $S_A times S_B$ slice of the joint probability tensor
  $cal(P)$ factors as
  $ bold(P)_((c)) = bold(M)_A "diag"(bold(pi) dot.o bold(m)_(C, c)) bold(M)_B^top, $
  where $bold(m)_(C, c) = (M_(C, c 1), dots, M_(C, c K))^top$. Summing over
  $c = 1, dots, S_C$ (using $sum_c M_(C, c k) = 1$) yields the marginal
  matrix $bold(P)_0 = bold(M)_A "diag"(bold(pi)) bold(M)_B^top$. Because
  $"rank"(bold(M)_A) = "rank"(bold(M)_B) = K$ and $pi_k > 0$, $bold(P)_0$
  has exact rank $K$; let its compact SVD be
  $bold(P)_0 = bold(U) bold(Sigma) bold(V)^top$ with
  $bold(U) in bb(R)^(S_A times K)$, $bold(V) in bb(R)^(S_B times K)$, and
  invertible $bold(Sigma) in bb(R)^(K times K)$. Then
  $bold(A)_0 := bold(U)^top bold(M)_A$ and
  $bold(B)_0 := bold(V)^top bold(M)_B$ are invertible $K times K$ matrices
  satisfying $bold(A)_0 "diag"(bold(pi)) bold(B)_0^top = bold(Sigma)$.

  Because the columns of $bold(M)_C$ are pairwise distinct, the finite
  union of hyperplanes
  $union.big_(k < l) {bold(w) in bb(R)^(S_C) : bold(w)^top (bold(m)_(C, dot k) - bold(m)_(C, dot l)) = 0}$
  has Lebesgue measure zero; hence there exists $bold(w) in bb(R)^(S_C)$
  such that $v_k := sum_(c=1)^(S_C) w_c M_(C, c k)$ are pairwise distinct
  across $k = 1, dots, K$. Forming the projected slice combination
  $bold(P)_w := sum_(c=1)^(S_C) w_c bold(P)_((c)) = bold(M)_A "diag"(bold(pi) dot.o bold(v)) bold(M)_B^top$,
  we obtain the $K times K$ observable matrix
  $ bold(W) := (bold(U)^top bold(P)_w bold(V)) bold(Sigma)^(-1) = bold(A)_0 "diag"(bold(v)) bold(A)_0^(-1). $
  Since the eigenvalues $v_1, dots, v_K$ of $bold(W)$ are simple, its
  eigenvector matrix $bold(A)_0$ is uniquely determined up to column
  permutation and nonzero column scaling. Because the column space of
  $bold(M)_A$ equals $"span"(bold(U))$, $bold(M)_A = bold(U) bold(A)_0$ is
  likewise unique up to column permutation and scaling. Enforcing the
  probability simplex constraint $bold(1)^top bold(m)_(A, dot k) = 1$
  uniquely pins every column scale. Finally,
  $bold(M)_B "diag"(bold(pi)) = bold(P)_0^top bold(M)_A (bold(M)_A^top bold(M)_A)^(-1)$
  uniquely fixes $bold(pi)$ (via column sums of $bold(M)_B$) and
  $bold(M)_B$, while solving
  $bold(M)_A^+ bold(P)_((c)) (bold(M)_B^+)^top = "diag"(bold(pi) dot.o bold(m)_(C, c))$
  recovers $bold(M)_C$. Individual item distributions $bold(theta)_(k j)$
  follow uniquely by marginalizing each block column.
]

Theorem 2 provides an explicit algebraic realization of Kruskal's tensor
rank identifiability principle for discrete mixtures @lazarsfeld1968latent
@mclachlan2000finite. It also directly certifies why $J >= 3$ binary items
suffice to identify $K = 2$ latent classes ($S_A = S_B = S_C = 2 >= K$),
whereas $J = 2$ items fail generically.

#subsec("C", "Non-Identifiability Counterexample for Two Manifest Items")

#proposition(title: "Continuum Non-Identifiability for Two Binary Variables", origin: "Original Counterexample")[
  For $J = 2$ binary variables ($M_1 = M_2 = 2$) and $K = 2$ latent
  classes, the LCA parameterization is strictly non-identifiable: there
  exists a one-parameter continuum of non-permutation-equivalent interior
  models generating the exact same $2 times 2$ manifest distribution.
]

#proof[
  Consider the symmetric correlated joint distribution
  $ bold(P) = mat(3 / 8, 1 / 8; 1 / 8, 3 / 8). $
  For any scalar $a in (1 / 4, 1 / 2)$, set $b = 1 / (8 a) in (1 / 4, 1 / 2)$
  and define the strictly interior parameter tuple $bold(Theta)(a)$ by
  $pi_1 = pi_2 = 1 / 2$,
  $ theta_(1 1 1) = frac(1, 2) + a, quad theta_(1 2 1) = frac(1, 2) + b, \
    theta_(2 1 1) = frac(1, 2) - a, quad theta_(2 2 1) = frac(1, 2) - b, $
  with $theta_(k j 2) = 1 - theta_(k j 1) in (0, 1)$. Direct expansion
  gives
  $ P_(1 1) (a) &= frac(1, 2) (frac(1, 2) + a)(frac(1, 2) + b) + frac(1, 2) (frac(1, 2) - a)(frac(1, 2) - b) \
    &= frac(1, 4) + a b = frac(1, 4) + frac(1, 8) = frac(3, 8), $
  and identically $P_(2 2) (a) = 1 / 4 + a b = 3 / 8$, while
  $P_(1 2) (a) = P_(2 1) (a) = 1 / 4 - a b = 1 / 8$. For distinct values
  $a_1, a_2 in (1 / 4, 1 / 2)$, the item-1 probabilities
  ${1 / 2 + a_1, 1 / 2 - a_1}$ and ${1 / 2 + a_2, 1 / 2 - a_2}$ are
  disjoint sets, so $bold(Theta)(a_1)$ and $bold(Theta)(a_2)$ cannot be
  related by class permutation yet yield identical likelihoods on all
  datasets.
]

#subsec("D", "Nominal Relabeling Invariance and Boundary Trapping")

#proposition(title: "Category-Permutation Invariance and Zero-Boundary Absorption", origin: "Original Proof")[
  (i) For any product of category permutations
  $(sigma_1, dots, sigma_J) in S_(M_1) times dots.c times S_(M_J)$,
  relabeling $x_(i j) |-> sigma_j (x_(i j))$ and
  $theta_(k j m) |-> theta_(k j, sigma_j^(-1) (m))$ leaves $cal(L)(bold(Theta))$
  and the posterior matrix $bold(Gamma)$ strictly invariant.
  (ii) In unregularized LCA-EM ($alpha = 0$), the boundary faces
  $theta_(k j m) = 0$ and $pi_k = 0$ are forward-invariant absorbing sets:
  if $theta_(k j m)^((t)) = 0$, then $theta_(k j m)^((s)) = 0$ for all
  $s >= t$. Conversely, Dirichlet smoothing $alpha > 0$ eliminates
  boundary trapping by enforcing
  $theta_(k j m)^((t)) >= alpha / (M_j alpha + N) > 0$ for all $t >= 1$.
]

#proof[
  Part (i) follows immediately from
  $theta_(k j, sigma_j^(-1) (sigma_j (x_(i j)))) = theta_(k j x_(i j))$,
  preserving every joint likelihood term $ell_(i k)$ identically.
  For (ii), suppose $theta_(k j m)^((t)) = 0$. For every sample $i$ with
  $x_(i j) = m$, the class-$k$ conditional likelihood vanishes:
  $product_(r=1)^J theta_(k r x_(i r))^((t)) = 0$, forcing
  $gamma_(i k)^((t)) = 0$ whenever $bb(I)(x_(i j) = m) = 1$. Hence the
  M-step numerator satisfies
  $sum_(i=1)^N gamma_(i k)^((t)) bb(I)(x_(i j) = m) = 0$, so
  $theta_(k j m)^((t+1)) = 0$. Induction on $s >= t$ establishes forward
  invariance (and analogously for $pi_k^((t)) = 0$). Under $alpha > 0$,
  because $N_k <= N$ and the numerator is at least $alpha$,
  $theta_(k j m)^((t+1)) >= alpha / (M_j alpha + N) > 0$, keeping every
  iterate strictly inside the interior of the product simplex.
]

#subsec("E", "Asymptotic Consistency of BIC Class Enumeration")

#theorem(title: "BIC Model-Selection Consistency for Identifiable LCA", origin: "Classical — Cited")[
  Suppose data are drawn i.i.d. from an identifiable $K_0$-class LCA
  model $cal(P)_0$ in the interior of the parameter space, and let
  $hat(K)_"BIC" = arg min_(1 <= K <= K_"max") (-2 cal(L)(hat(bold(Theta))_K) + nu_K log N)$
  with $K_"max" >= K_0$ @schwarz1978bic @mclachlan2000finite. Then
  $"Pr"(hat(K)_"BIC" = K_0) -> 1$ as $N -> infinity$.
]

#proof[
  Let $hat(cal(P))_N$ be the empirical cell-frequency tensor over the $S$
  response cells, so $N^(-1) cal(L)(bold(Theta)) = -H(hat(cal(P))_N) - D_"KL" (hat(cal(P))_N || cal(P)(bold(Theta)))$.
  For any under-specified $K < K_0$, because the true rank-$K_0$ model is
  identifiable and not representable in the closed set $cal(M)_K$ of
  $K$-class mixtures, $delta_K := inf_(cal(P) in cal(M)_K) D_"KL" (cal(P)_0 || cal(P)) > 0$.
  By the law of large numbers,
  $2(cal(L)(hat(bold(Theta))_(K_0)) - cal(L)(hat(bold(Theta))_K)) = 2 N delta_K + o_p (N)$,
  which dominates $(nu_(K_0) - nu_K) log N$ so
  $"Pr"("BIC"(K) > "BIC"(K_0)) -> 1$.

  For any over-specified $K > K_0$, the unconstrained multinomial
  log-likelihood $cal(L)_"sat" = -N H(hat(cal(P))_N)$ upper-bounds
  $cal(L)(hat(bold(Theta))_K)$. Because $cal(P)_0$ lies in the interior of
  the $S$-cell simplex, Pearson's chi-square / Wilks' theorem for the
  $K_0$-class model gives
  $ 0 &<= 2(cal(L)(hat(bold(Theta))_K) - cal(L)(hat(bold(Theta))_(K_0))) \
    &<= 2(cal(L)_"sat" - cal(L)(hat(bold(Theta))_(K_0))) ->^d chi^2_(S - 1 - nu_(K_0)). $
  Hence the likelihood ratio gain is $cal(O)_p (1)$ even across the
  non-regular boundary of redundant components, whereas the penalty gap
  $(nu_K - nu_(K_0)) log N = (K - K_0)(M - J + 1) log N -> +infinity$.
  Thus $"Pr"("BIC"(K) > "BIC"(K_0)) -> 1$ for all $K != K_0$.
]

Notice how using the saturated multinomial upper bound $cal(L)_"sat"$
bypasses the singular Fisher-information obstacle of over-fitted mixtures
cleanly and rigorously!

#sec("III", "Complexity, Geometry and Operational Profile")

During each EM iteration, evaluating $ell_(i k)$ visits only the single
observed category $x_(i j)$ per item, requiring $cal(O)(N K J)$ additions.
Accumulating the sufficient counts $N_(k j m)$ likewise touches only
active indicator entries $(i, j, x_(i j))$ in $cal(O)(N K J)$, while
normalizing the parameter tables costs $cal(O)(K M)$ where
$M = sum_(j=1)^J M_j$. Thus total time complexity is
$cal(O)(T K (N J + M))$ (or $cal(O)(T K (U J + M))$ when compressed to
$U <= min(N, S)$ unique response patterns), with space complexity
$cal(O)(N J + N K + K M)$.

Geometrically, LCA does not partition a metric space into convex balls;
instead, it decomposes a $J$-way probability simplex tensor into a convex
combination of $K$ rank-1 product tensors. In log-odds space, posterior
log-ratios between classes $k$ and $l$ are additive linear functions of
one-hot item indicators:
$ log frac(gamma_(i k), gamma_(i l)) = log frac(pi_k, pi_l) + sum_(j=1)^J log frac(theta_(k j x_(i j)), theta_(l j x_(i j))), $
which is structurally identical to a generative Naive Bayes classifier
trained without observed labels @bishop2006pattern.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [*Property*], [*Characterization*],
    [Time complexity], [$cal(O)(T K (N J + M))$; compresses to $cal(O)(T K (U J + M))$],
    [Space complexity], [$cal(O)(N J + N K + K M)$],
    [Free parameters], [$nu_K = K(M - J + 1) - 1$; requires $nu_K <= product_j M_j - 1$],
    [Inductive geometry], [Rank-$K$ nonnegative CP tensor on discrete contingency table],
    [Invariance group], [Arbitrary nominal relabeling $S_(M_1) times dots.c times S_(M_J)$],
    [Missing data], [Native MAR handling by omitting unobserved $j$ in product],
    [Local optima], [Non-concave likelihood; requires multi-start EM + smoothing],
    [Model selection], [Consistent BIC @schwarz1978bic; $G^2$ goodness-of-fit test],
  ),
  caption: [Operational and geometric characterization of Latent Class Analysis.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

+ *Clinical Phenotyping and Diagnostic Nosology:* Electronic health
  records encode comorbidities, binary symptom checklists, and nominal
  medication classes with frequent missingness. LCA uncovers latent
  disease endotypes without a gold-standard diagnostic test, natively
  marginalizing missing-at-random items by skipping unobserved indices in
  the local-independence product @mclachlan2000finite.
+ *Psychometrics, Epidemiology, and Social Survey Research:* Multi-item
  questionnaires measure unobserved behavioral profiles (e.g., consumer
  adoption typologies or risk-factor syndromes) via discrete Likert or
  nominal responses @lazarsfeld1968latent. Category-permutation invariance
  (Proposition 2) prevents arbitrary integer coding of nominal options
  from biasing cluster recovery.
+ *Bioinformatics and Categorical Genotype Clustering:* Single-nucleotide
  polymorphism (SNP) panels take discrete values ${0, 1, 2}$ across loci.
  Under approximate within-population Hardy-Weinberg and linkage
  equilibrium, conditional independence across loci directly mirrors the
  LCA generative factor structure, recovering ancestral subpopulations and
  soft admixture posteriors.
+ *Fraud Detection and Behavioral Event Logs:* Transaction risk flags,
  device-type indicators, merchant categories, and authentication channels
  form heterogeneous nominal tuples. Because LCA outputs calibrated log-odds
  contributions per attribute, investigators obtain transparent, additive
  item-level explanations for each anomalous latent segment.
+ *Customer Segmentation with Mixed Discretized Attributes:* In retail
  CRM databases combining region, channel, loyalty tier, and binned
  purchase frequencies, LCA avoids tuning ad-hoc Gower or Hamming feature
  weights, instead learning per-class item response distributions via
  principled maximum likelihood and BIC selection @biernacki2006mixmod.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Four mathematically distinct failure modes constrain LCA deployment:
(1) *Local-independence violation:* Residual pairwise correlations within
a true subpopulation force LCA to manufacture spurious extra classes to
absorb intra-class dependence (detectable via bivariate residual $G^2$
diagnostics).
(2) *Discretization loss on continuous features:* Binning continuous
measurements discards fine metric ordering and introduces boundary
quantization artifacts, reducing Euclidean compactness.
(3) *High-cardinality sparsity and boundary trapping:* When $M_j$ is
large relative to $N$, empty class-category cells drive $theta_(k j m) -> 0$
(Proposition 2), causing EM to lock onto suboptimal boundary faces unless
Dirichlet smoothing ($alpha > 0$) is applied.
(4) *Under-identification:* Too few items ($J = 2$, Proposition 1) or
excessive $K$ violating $nu_K <= S - 1$ makes class profiles
unrecoverable regardless of sample size.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [*Domain*], [*Mathematical Fit*], [*Caution / Alternative*],
    [Nominal surveys & clinical panels], [Exact discrete likelihood; native MAR handling], [Check bivariate residuals for local dependence],
    [Categorical transaction logs], [Permutation-invariant; interpretable item profiles], [Apply Dirichlet smoothing $alpha > 0$ for rare codes],
    [Continuous Euclidean data], [Requires binning; loses metric geometry], [Use GMM, vMF (Alg. 89), or HDP (Alg. 91)],
    [Two-item contingency tables], [Non-identifiable for $K >= 2$ (Prop. 1)], [Use NMF (Alg. 92) or correspondence analysis],
    [Correlated count corpora], [Multi-word counts exceed fixed-item schema], [Use LDA / HDP (Alg. 91) or pLSA],
  ),
  caption: [Domain suitability and failure-mode matrix for LCA.],
)

#sec("V", "Structural Relations and Empirical Behavior")

LCA serves as the canonical discrete member of the exponential-family
mixture hierarchy @banerjee2005clustering. It directly generalizes
$k$-modes @huang1998extensions (which replaces soft likelihood products
with hard Hamming distances) and sits inside mixed-type Bayesian mixture
frameworks such as AutoClass @cheeseman1996autoclass, Snob
@wallace1968snob, and MIXMOD @biernacki2006mixmod, as well as graph-based
categorical clustering like ROCK @guha2000rock. Formally, LCA also links
directly to $k$-modes via a small-temperature isotropic limit.

#proposition(title: "Isotropic Small-Temperature Reduction to k-Modes", origin: "Original Derivation")[
  Restrict LCA to equal priors $pi_k = 1 / K$ and symmetric per-item
  mode distributions with dominant category $m_(k j) in {1, dots, M_j}$:
  $ theta_(k j m) = cases(1 - epsilon_j & "if" m = m_(k j), epsilon_j / (M_j - 1) & "otherwise",) $
  parameterized by temperature $tau > 0$ via
  $epsilon_j = (M_j - 1) / (M_j - 1 + exp(1 / tau))$. Then as
  $tau -> 0^+$, the scaled log-likelihood satisfies
  $ lim_(tau -> 0^+) tau (cal(L)(bold(Theta)) + N log K) = -sum_(i=1)^N min_(1 <= k <= K) d_H (bold(x)_i, bold(m)_k), $
  where $d_H (bold(x)_i, bold(m)_k) = sum_(j=1)^J bb(I)(x_(i j) != m_(k j))$
  is the Hamming distance of $k$-modes @huang1998extensions.
]

#proof[
  By construction, $1 - epsilon_j = exp(1 / tau) / (M_j - 1 + exp(1 / tau))$
  and $epsilon_j / (M_j - 1) = 1 / (M_j - 1 + exp(1 / tau))$. Hence the
  class-$k$ log-likelihood of observation $bold(x)_i$ equals
  $ sum_(j=1)^J log theta_(k j x_(i j)) &= frac(J - d_H (bold(x)_i, bold(m)_k), tau) \
    &- sum_(j=1)^J log(M_j - 1 + exp(1 / tau)). $
  Observe that $tau log(M_j - 1 + exp(1 / tau)) = 1 + cal(O)(tau exp(-1 / tau))$,
  so $tau sum_(j=1)^J log theta_(k j x_(i j)) = -d_H (bold(x)_i, bold(m)_k) + cal(O)(tau)$.
  Applying the log-sum-exp sandwich across the $K$ classes (with gap at
  most $tau log K -> 0$) yields
  $tau (log p(bold(x)_i mid(|) bold(Theta)) + log K) -> -min_k d_H (bold(x)_i, bold(m)_k)$,
  and summing over $i = 1, dots, N$ completes the proof.
]

In the ClusterOverview 4-dataset empirical benchmark (where continuous
benchmark features are discretized prior to LCA fitting; evaluated with
Silhouette @rousseeuw1987silhouettes and ARI @hubert1985comparing), LCA
runs in a modest *3.064 ms* mean runtime and displays a sharp diagnostic
split between class recovery and Euclidean compactness:
- *Iris:* Silhouette $= 0.386$, $"ARI" = 0.638$. Remarkably, LCA attains
  the *highest Iris ARI* ($0.638$) among Algs. 89--92: discretizing petal
  and sepal features into categorical strata suppresses continuous
  tail-outlier leverage and aligns cleanly with the three botanical
  species via independent per-feature log-odds, even though discretized
  boundaries slightly lower continuous Euclidean Silhouette ($0.386$).
- *Customers:* Silhouette $= 0.182$. Correlated wholesale expenditure
  channels strongly violate local independence after quantization,
  producing coarse axis-aligned bin splits with weak Euclidean cohesion.
- *Multi-Shape:* Silhouette $= 0.709$ ($"ARI" =$ N/A). Disjoint spatial
  regimes separate cleanly into distinct categorical response profiles.
- *High-D:* Silhouette $= 0.302$ ($"ARI" =$ N/A). Additive per-dimension
  log-odds prevent full-covariance singularity, matching movMF ($0.299$)
  in high ambient dimensionality.

When data are genuinely nominal or discretized into meaningful strata,
LCA remains the gold-standard statistically identifiable latent-group
model @lazarsfeld1968latent @mclachlan2000finite.
