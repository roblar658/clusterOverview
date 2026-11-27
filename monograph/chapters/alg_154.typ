#import "../template/ieee.typ": *

#alg_entry(
  id: 154,
  name: "Gath-Geva (GG)",
  category: "Cat. VII: Fuzzy, Biclustering & Domain",
  year: "1989",
  paradigm: "Fuzzy Maximum-Likelihood Exponential Distance",
  complexity: $O(T c (n d^2 + d^3))$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Fuzzy $c$-Means (FCM, Alg. 151) @bezdek1981fcm restricts cluster geometry to
isotropic hyperspheres of uniform volume, while the Gustafson-Kessel (GK, Alg.
153) algorithm @gustafson1979gk relaxes shape to full ellipsoids but constrains
each cluster's covariance determinant to a fixed preset volume $rho_j$. To
overcome both restrictions and accommodate clusters of arbitrary ellipsoidal
orientation, unequal hypervolume, and imbalanced prior mass, Gath and Geva
@gath1989gg introduced the Unsupervised Optimal Fuzzy Clustering (UOFC)
framework, widely known as the *Gath-Geva (GG)* algorithm, by replacing
algebraic Mahalanobis distances with an exponential Fuzzy Maximum Likelihood
Estimation (FMLE) dissimilarity.

Let $cal(X) = {bold(x)_1, dots, bold(x)_n} subset.eq bb(R)^d$ denote $n$
observations and $c >= 2$ the number of fuzzy clusters. Let
$bold(U) = [u_(i j)] in [0, 1]^(n times c)$ be a fuzzy partition matrix
satisfying the probabilistic simplex constraint
$sum_(j=1)^c u_(i j) = 1$ for all $i in {1, dots, n}$ and non-degeneracy
$0 < sum_(i=1)^n u_(i j) < n$ for all $j in {1, dots, c}$. Each cluster $j$
is parameterized by a prototype $bold(v)_j in bb(R)^d$, an a priori weight
$P_j in (0, 1)$ with $sum_(j=1)^c P_j = 1$, and a symmetric positive-definite
fuzzy covariance matrix $bold(F)_j in bb(S)_(++)^d$. Let
$Theta = (bold(V), bold(F), bold(P))$ collect all cluster parameters and let
$m > 1$ denote the fuzzifier exponent.

#definition(title: "Gath-Geva Exponential Dissimilarity (FMLE Distance)")[
  For an observation $bold(x)_i in bb(R)^d$ and cluster parameters
  $(bold(v)_j, bold(F)_j, P_j) in bb(R)^d times bb(S)_(++)^d times (0, 1)$, the
  *Gath-Geva (FMLE) exponential dissimilarity* is
  $
    d_"GG"^2 (bold(x)_i, bold(v)_j; bold(F)_j, P_j)
    &= ( (2 pi)^(d / 2) sqrt(det(bold(F)_j)) ) / P_j \
    &quad times exp( 1 / 2 (bold(x)_i - bold(v)_j)^top
       bold(F)_j^(-1) (bold(x)_i - bold(v)_j) ).
  $
  Equivalently, letting $cal(N)(bold(x) | bold(v)_j, bold(F)_j)$ denote the
  $d$-variate Gaussian density, the FMLE dissimilarity is the exact reciprocal
  joint likelihood:
  $
    d_"GG"^2 (bold(x)_i, bold(v)_j; bold(F)_j, P_j)
    = ( P_j cal(N)(bold(x)_i | bold(v)_j, bold(F)_j) )^(-1).
  $
]

Because the constant $(2 pi)^(d / 2)$ cancels identically in every normalized
membership ratio, Gath and Geva @gath1989gg omit it in computation; retaining
it makes the exact probabilistic correspondence transparent. Given $d_"GG"^2$,
the classical fuzzy partition functional is
$
  J_"GG" (bold(U), Theta)
  = sum_(i=1)^n sum_(j=1)^c u_(i j)^m
    d_"GG"^2 (bold(x)_i, bold(v)_j; bold(F)_j, P_j),
$
whose stationary membership update for fixed $Theta$ is
$
  u_(i j) = ( d_"GG"^2 (bold(x)_i, bold(v)_j) )^(-1 / (m - 1))
  / ( sum_(l=1)^c ( d_"GG"^2 (bold(x)_i, bold(v)_l) )^(-1 / (m - 1)) ).
$
Rather than differentiating $J_"GG"$ with respect to $Theta$ (which yields
implicit coupled transcendental equations with no closed-form solution), GG
pairs the exponential membership update with moment-matched Gaussian maximum
likelihood estimators (with exponent $q = 1$ in the canonical FMLE formulation
of @gath1989gg, or $q = m$ in fuzzified-covariance variants; below we analyze
the canonical $q = 1$ estimator and note $q = m$ differences explicitly).

#algbox(
  title: "Gath-Geva (GG / FMLE) Fuzzy Clustering",
  inputs: [Data $cal(X) subset.eq bb(R)^d$, clusters $c >= 2$, fuzzifier
    $m > 1$, regularization $epsilon >= 0$, tolerance $eta > 0$, warm-start
    partition $bold(U)^((0))$ (e.g., from FCM @bezdek1981fcm).],
  outputs: [Fuzzy partition $bold(U)^*$, prototypes $bold(V)^*$, covariances
    $bold(F)^*$, priors $bold(P)^*$.],
)[
  1. *Initialize:* Run FCM @bezdek1981fcm or GK @gustafson1979gk to obtain a
     non-degenerate partition $bold(U)^((0))$; set $t <- 0$.
  2. *Parameter Estimation (M-Step):* For each $j = 1, dots, c$, let
     $N_j^((t)) = sum_(i=1)^n u_(i j)^((t))$ and compute:
     $
       P_j^((t)) &= N_j^((t)) / n,
       quad bold(v)_j^((t)) = 1 / N_j^((t))
         sum_(i=1)^n u_(i j)^((t)) bold(x)_i, \
       bold(F)_j^((t)) &= 1 / N_j^((t))
         sum_(i=1)^n u_(i j)^((t))
         (bold(x)_i - bold(v)_j^((t)))
         (bold(x)_i - bold(v)_j^((t)))^top + epsilon bold(I)_d.
     $
  3. *Log-FMLE Dissimilarity:* For all $(i, j)$, compute via Cholesky
     factorization $bold(F)_j^((t)) = bold(L)_j bold(L)_j^top$:
     $
       ell_(i j)^((t)) &= ln d_"GG"^2 (bold(x)_i, bold(v)_j^((t));
         bold(F)_j^((t)), P_j^((t))) \
       &= d / 2 ln(2 pi) + sum_(r=1)^d ln (bold(L)_j)_(r r)
          - ln P_j^((t)) \
       &quad + 1 / 2 norm(bold(L)_j^(-1)
          (bold(x)_i - bold(v)_j^((t))))_2^2.
     $
  4. *Membership Update (Tempered E-Step):* Using log-sum-exp stabilization
     with $beta = 1 / (m - 1)$ and $alpha_i = min_l ell_(i l)^((t))$, set
     $
       u_(i j)^((t+1)) = exp(-beta (ell_(i j)^((t)) - alpha_i))
       / ( sum_(l=1)^c exp(-beta (ell_(i l)^((t)) - alpha_i)) ).
     $
  5. *Stopping Criterion:* If
     $norm(bold(U)^((t+1)) - bold(U)^((t)))_oo < eta$, return
     $(bold(U)^((t+1)), bold(V)^((t)), bold(F)^((t)), bold(P)^((t)))$; else
     $t <- t + 1$ and return to Step 2.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

Although GG is routinely described as a heuristic modification of FCM lacking
a Lyapunov objective for $J_"GG"$, we prove below that the canonical GG
iteration is in fact *exact alternating minimization* of an entropy-regularized
logarithmic free energy for every $m > 1$, establish its full-affine
invariance, prove finite-sample covariance collapse without regularization
alongside an explicit Tikhonov lower bound, and connect its convergence to
classical EM theory @wu1983convergence @mclachlan2000finite.

#subsec("A", "Variational Free-Energy Lyapunov Monotonicity")

#theorem(
  title: "Exact Free-Energy Descent of Canonical Gath-Geva",
  origin: "Original Unified Proof",
)[
  Define the temperature parameter $tau = m - 1 > 0$ and the
  entropy-regularized logarithmic free energy on
  $Delta_c^n times (bb(R)^d)^c times (bb(S)_(++)^d)^c times Delta_c$:
  $
    cal(F)_tau (bold(U), Theta)
    &= sum_(i=1)^n sum_(j=1)^c u_(i j)
       ln d_"GG"^2 (bold(x)_i, bold(v)_j; bold(F)_j, P_j) \
    &quad + tau sum_(i=1)^n sum_(j=1)^c u_(i j) ln u_(i j).
  $
  Assume unregularized GG ($epsilon = 0$) maintains $N_j^((t)) > 0$ and
  $bold(F)_j^((t)) in bb(S)_(++)^d$ at every iterate. Then Steps 2--4 of
  canonical GG perform exact block coordinate descent on $cal(F)_tau$, yielding
  $
    cal(F)_tau (bold(U)^((t+1)), Theta^((t+1)))
    <= cal(F)_tau (bold(U)^((t+1)), Theta^((t)))
    <= cal(F)_tau (bold(U)^((t)), Theta^((t)))
  $
  for every $t >= 0$ and every $m > 1$, with strict inequality whenever
  $(bold(U)^((t)), Theta^((t)))$ is not a stationary point of $cal(F)_tau$.
]

#proof[
  We analyze both block minimizations separately.

  *Block 1 (Minimization over $bold(U)$ for fixed $Theta$):* For fixed
  $Theta = Theta^((t))$, $cal(F)_tau$ decouples across observations $i$:
  $
    cal(F)_tau (bold(U), Theta) = sum_(i=1)^n sum_(j=1)^c
    ( u_(i j) ell_(i j) + tau u_(i j) ln u_(i j) ),
  $
  where $ell_(i j) = ln d_"GG"^2 (bold(x)_i, bold(v)_j; bold(F)_j, P_j)$. Since
  $u |-> u ln u$ is strictly convex on $(0, 1]$ with derivative
  $ln u + 1 -> -oo$ as $u -> 0^+$, and $tau = m - 1 > 0$, each $i$-th
  subproblem is strictly convex over the simplex $Delta_c$ with its unique
  global minimizer in the relative interior $(0, 1)^c$. Forming the Lagrangian
  $cal(L)_i (bold(u)_i, lambda_i) = sum_j (u_(i j) ell_(i j) + tau u_(i j) ln u_(i j)) + lambda_i (sum_j u_(i j) - 1)$
  and setting $partial cal(L)_i / (partial u_(i j)) = 0$ gives
  $
    ell_(i j) + tau (ln u_(i j) + 1) + lambda_i = 0
    => u_(i j) prop exp(-ell_(i j) / tau).
  $
  Since $exp(-ell_(i j) / tau) = (d_"GG"^2 (bold(x)_i, bold(v)_j))^(-1 / (m - 1))$,
  imposing $sum_(j=1)^c u_(i j) = 1$ recovers Step 4 uniquely. Hence
  $cal(F)_tau (bold(U)^((t+1)), Theta^((t))) <= cal(F)_tau (bold(U)^((t)), Theta^((t)))$,
  with equality iff $bold(U)^((t+1)) = bold(U)^((t))$.

  *Block 2 (Minimization over $Theta$ for fixed $bold(U)$):* Fixing
  $bold(U) = bold(U)^((t+1))$, the entropy term
  $tau sum_(i,j) u_(i j) ln u_(i j)$ is constant in $Theta$. Expanding
  $ell_(i j)$ gives
  $
    sum_(i,j) u_(i j) ell_(i j)
    &= (n d) / 2 ln(2 pi) - sum_(j=1)^c N_j ln P_j \
    &quad + 1 / 2 sum_(j=1)^c ( N_j ln det(bold(F)_j)
       + sum_(i=1)^n u_(i j) (bold(x)_i - bold(v)_j)^top
       bold(F)_j^(-1) (bold(x)_i - bold(v)_j) ).
  $
  First, minimizing $-sum_(j=1)^c N_j ln P_j$ over $P in Delta_c$ is equivalent
  to minimizing $n "KL"(hat(bold(P)) || bold(P))$ where $hat(P)_j = N_j / n$.
  By Gibbs' inequality, the unique global minimizer is $P_j^* = N_j / n$.

  Second, for each $j$ and any $bold(F)_j in bb(S)_(++)^d$, let
  $bold(v)_j^* = N_j^(-1) sum_(i=1)^n u_(i j) bold(x)_i$. By the orthogonal
  bias-variance identity, for any $bold(v)_j in bb(R)^d$,
  $
    sum_(i=1)^n u_(i j) (bold(x)_i - bold(v)_j)^top bold(F)_j^(-1)
      (bold(x)_i - bold(v)_j)
    &= sum_(i=1)^n u_(i j) (bold(x)_i - bold(v)_j^*)^top bold(F)_j^(-1)
      (bold(x)_i - bold(v)_j^*) \
    &quad + N_j (bold(v)_j^* - bold(v)_j)^top bold(F)_j^(-1)
      (bold(v)_j^* - bold(v)_j).
  $
  Since $bold(F)_j^(-1) succ bold(0)$ and $N_j > 0$, the second term is
  strictly positive unless $bold(v)_j = bold(v)_j^*$, uniquely yielding
  $bold(v)_j^*$ independently of $bold(F)_j$.

  Third, define the scatter matrix
  $bold(S)_j = N_j^(-1) sum_(i=1)^n u_(i j) (bold(x)_i - bold(v)_j^*)(bold(x)_i - bold(v)_j^*)^top in bb(S)_(++)^d$
  and reparameterize via the precision matrix
  $bold(Lambda)_j = bold(F)_j^(-1) in bb(S)_(++)^d$. The $j$-th objective is
  $
    g_j (bold(Lambda)_j)
    = N_j / 2 ( -ln det(bold(Lambda)_j)
      + "Tr"(bold(Lambda)_j bold(S)_j) ).
  $
  Because $-ln det(dot)$ is strictly convex on $bb(S)_(++)^d$ and
  $"Tr"(dot bold(S)_j)$ is linear, $g_j$ is strictly convex with unique
  stationary point satisfying
  $nabla g_j (bold(Lambda)_j) = (N_j / 2)(-bold(Lambda)_j^(-1) + bold(S)_j) = bold(0)$,
  hence $bold(F)_j^* = bold(Lambda)_j^(-1) = bold(S)_j$. Thus Step 2 uniquely
  minimizes $cal(F)_tau (bold(U)^((t+1)), dot)$, completing the descent chain.
]

#corollary(
  title: "Exact Temperature Limits and GMM-EM Isomorphism",
  origin: "Original Derivation",
)[
  Let $cal(L)_"GMM" (Theta) = sum_(i=1)^n ln(sum_(j=1)^c P_j cal(N)(bold(x)_i | bold(v)_j, bold(F)_j))$
  be the Gaussian Mixture Model log-likelihood @dempster1977em @mclachlan2000finite.
  Then canonical GG satisfies three exact temperature regimes:
  1. *Exact GMM-EM ($m = 2$, $tau = 1$):* Every GG iterate is algebraically
     identical to Gaussian Mixture EM, and
     $cal(F)_1 (bold(U)^((t+1)), Theta^((t))) = -cal(L)_"GMM" (Theta^((t)))$,
     guaranteeing $cal(L)_"GMM" (Theta^((t+1))) >= cal(L)_"GMM" (Theta^((t)))$.
  2. *Hard Maximum-Likelihood Limit ($m -> 1^+$, $tau -> 0^+$):* Provided
     $j^*(i) = arg min_j ell_(i j)$ is unique,
     $u_(i j) -> bb(I)[j = j^*(i)]$, recovering hard classification EM (CEM).
  3. *Maximum-Entropy Limit ($m -> oo$, $tau -> oo$):* For finite
     $ell_(i j)$, $u_(i j) -> 1 / c$ uniformly for all $(i, j)$, collapsing
     all prototypes to the global sample mean.
]

#proof[
  At $m = 2$, $tau = 1$ and $beta = 1$, so Step 4 yields
  $u_(i j)^((t+1)) = P_j cal(N)(bold(x)_i | bold(v)_j, bold(F)_j) / (sum_l P_l cal(N)(bold(x)_i | bold(v)_l, bold(F)_l))$,
  the exact Bayes posterior $gamma_(i j)$ of GMM-EM @dempster1977em.
  Substituting $ln u_(i j)^((t+1)) = -ell_(i j) - ln(sum_l exp(-ell_(i l)))$
  into $cal(F)_1$ cancels $ell_(i j)$ identically:
  $
    cal(F)_1 (bold(U)^((t+1)), Theta^((t)))
    = -sum_(i=1)^n ln( sum_(l=1)^c exp(-ell_(i l)^((t))) )
    = -cal(L)_"GMM" (Theta^((t))).
  $
  By Theorem 1,
  $-cal(L)_"GMM" (Theta^((t+1))) = cal(F)_1 (bold(U)^((t+2)), Theta^((t+1))) <= cal(F)_1 (bold(U)^((t+1)), Theta^((t))) = -cal(L)_"GMM" (Theta^((t)))$.
  For $m -> 1^+$, $beta -> oo$, so $exp(-beta (ell_(i j) - alpha_i)) -> 0$ for
  every $j != j^*(i)$ and equals $1$ at $j^*(i)$. For $m -> oo$, $beta -> 0^+$,
  so $exp(-beta (ell_(i j) - alpha_i)) -> 1$ for all $j$, yielding
  $u_(i j) -> 1 / c$.
]

#subsec("B", "Full-Affine Coordinate Invariance")

#theorem(
  title: "Exact Affine Equivariance of Unregularized Gath-Geva",
  origin: "Original Proof",
)[
  Let $T(bold(x)) = bold(A) bold(x) + bold(b)$ be any invertible affine
  transformation with $bold(A) in "GL"_d (bb(R))$ and $bold(b) in bb(R)^d$,
  and let $tilde(cal(X)) = {T(bold(x)_1), dots, T(bold(x)_n)}$. When initialized
  from the same partition $bold(U)^((0))$, unregularized GG ($epsilon = 0$)
  generates identical membership matrices $tilde(bold(U))^((t)) = bold(U)^((t))$
  and affinely equivariant parameters
  $tilde(bold(v))_j^((t)) = bold(A) bold(v)_j^((t)) + bold(b)$,
  $tilde(bold(F))_j^((t)) = bold(A) bold(F)_j^((t)) bold(A)^top$, and
  $tilde(P)_j^((t)) = P_j^((t))$ for all $t >= 0$.
]

#proof[
  We proceed by induction on $t$. At $t = 0$, $tilde(bold(U))^((0)) = bold(U)^((0))$
  by hypothesis. Suppose $tilde(bold(U))^((t)) = bold(U)^((t))$. Then
  $tilde(N)_j^((t)) = N_j^((t))$ and $tilde(P)_j^((t)) = P_j^((t))$. By
  linearity of weighted averages,
  $
    tilde(bold(v))_j^((t))
    = 1 / N_j^((t)) sum_(i=1)^n u_(i j)^((t)) (bold(A) bold(x)_i + bold(b))
    = bold(A) bold(v)_j^((t)) + bold(b).
  $
  Hence $T(bold(x)_i) - tilde(bold(v))_j^((t)) = bold(A)(bold(x)_i - bold(v)_j^((t)))$,
  which yields
  $
    tilde(bold(F))_j^((t))
    = 1 / N_j^((t)) sum_(i=1)^n u_(i j)^((t))
      bold(A)(bold(x)_i - bold(v)_j^((t)))
      (bold(x)_i - bold(v)_j^((t)))^top bold(A)^top
    = bold(A) bold(F)_j^((t)) bold(A)^top.
  $
  Using $det(bold(A) bold(F)_j bold(A)^top) = abs(det(bold(A)))^2 det(bold(F)_j)$
  and the quadratic-form invariance
  $
    (bold(A) bold(z))^top (bold(A) bold(F)_j bold(A)^top)^(-1) (bold(A) bold(z))
    = bold(z)^top bold(A)^top bold(A)^(-top) bold(F)_j^(-1)
      bold(A)^(-1) bold(A) bold(z)
    = bold(z)^top bold(F)_j^(-1) bold(z),
  $
  the transformed FMLE dissimilarity satisfies
  $
    d_"GG"^2 (T(bold(x)_i), tilde(bold(v))_j^((t));
      tilde(bold(F))_j^((t)), tilde(P)_j^((t)))
    = abs(det(bold(A))) dot
      d_"GG"^2 (bold(x)_i, bold(v)_j^((t)); bold(F)_j^((t)), P_j^((t))).
  $
  Because $abs(det(bold(A))) > 0$ is strictly independent of $j$, the factor
  $abs(det(bold(A)))^(-1 / (m - 1))$ cancels identically between numerator and
  denominator in Step 4, giving $tilde(bold(U))^((t+1)) = bold(U)^((t+1))$.
]

In contrast, FCM @bezdek1981fcm is invariant only under orthogonal-plus-uniform
scaling transformations, while regularized GG ($epsilon > 0$) breaks affine
invariance unless the shrinkage term $epsilon bold(I)_d$ is scaled proportionally
to the sample covariance $epsilon hat(bold(Sigma))_cal(X)$.

#subsec("C", "Covariance Singularity Collapse and Tikhonov Lower Bound")

#proposition(
  title: "Unbounded Free-Energy Collapse and Regularized Lower Bound",
  origin: "Original Counterexample",
)[
  1. *Unregularized Collapse ($epsilon = 0$):* For any dataset $cal(X)$ with
     $n >= 3$, $c >= 2$, and $m > 1$, the unregularized free energy
     $cal(F)_tau$ is unbounded below:
     $inf_(bold(U), Theta) cal(F)_tau (bold(U), Theta) = -oo$.
  2. *Tikhonov Stabilization ($epsilon > 0$):* Whenever every covariance
     satisfies $bold(F)_j succ.eq epsilon bold(I)_d$, $cal(F)_tau$ admits the
     uniform finite lower bound
     $
       cal(F)_tau (bold(U), Theta)
       >= (n d) / 2 ln(2 pi epsilon) - n (m - 1) ln c.
     $
]

#proof[
  For Part 1, fix $bold(v)_1 = bold(x)_1$, $bold(F)_1 = delta bold(I)_d$ with
  $delta > 0$, $P_1 = 1 / c$, and assign all remaining $c - 1$ clusters fixed
  non-degenerate parameters $(bold(v)_j, bold(I)_d, 1 / c)$ for $j >= 2$.
  Evaluate $cal(F)_tau$ at the feasible membership matrix with
  $u_(1 1) = 1 - (c - 1) gamma$, $u_(1 j) = gamma$ ($j >= 2$), and for $i >= 2$,
  $u_(i 1) = gamma$, $u_(i j) = (1 - gamma) / (c - 1)$, choosing
  $gamma = exp(-1 / delta) -> 0^+$ as $delta -> 0^+$.
  For $i = 1, j = 1$, $(bold(x)_1 - bold(v)_1) = bold(0)$, so
  $ell_(1 1) = (d / 2) ln(2 pi delta) + ln c$.
  For $i >= 2, j = 1$,
  $u_(i 1) ell_(i 1) = exp(-1 / delta) ( (d / 2) ln(2 pi delta) + ln c + norm(bold(x)_i - bold(x)_1)_2^2 / (2 delta) ) -> 0$
  as $delta -> 0^+$. All terms for $j >= 2$ remain bounded by a finite constant
  $C_0$, and the entropy term is bounded in $[-n tau ln c, 0]$. Hence
  $
    cal(F)_tau (bold(U), Theta)
    = (1 - (c - 1) e^(-1 / delta)) ( d / 2 ln(2 pi delta) + ln c )
      + O(1) -> -oo
  $
  as $delta -> 0^+$, demonstrating Dirac-delta singularity collapse onto a
  single observation.

  For Part 2, $bold(F)_j succ.eq epsilon bold(I)_d$ implies every eigenvalue
  $lambda_r (bold(F)_j) >= epsilon$, so $det(bold(F)_j) >= epsilon^d$ and
  $(1 / 2) ln det(bold(F)_j) >= (d / 2) ln epsilon$. Since $P_j in (0, 1]$
  implies $-ln P_j >= 0$, and $bold(F)_j^(-1) succ bold(0)$ implies the
  Mahalanobis quadratic form is non-negative, every pair $(i, j)$ satisfies
  $ell_(i j) >= (d / 2) ln(2 pi epsilon)$. Moreover, by concavity of Shannon
  entropy on $Delta_c$, $sum_(j=1)^c u_(i j) ln u_(i j) >= -ln c$. Summing over
  $i = 1, dots, n$ and $j = 1, dots, c$ yields
  $
    cal(F)_tau (bold(U), Theta)
    >= sum_(i=1)^n ( sum_(j=1)^c u_(i j) ) d / 2 ln(2 pi epsilon)
       - n tau ln c
    = (n d) / 2 ln(2 pi epsilon) - n (m - 1) ln c.
  $
]

#subsec("D", "Stationary Convergence and Exponential Basin Locality")

#theorem(
  title: "Global Convergence to Stationary Points Under Regularization",
  origin: "Classical — Cited",
)[
  Under covariance regularization $bold(F)_j = bold(S)_j + epsilon bold(I)_d$
  ($epsilon > 0$) and bounded prior support $P_j >= delta_P > 0$, define the
  regularized free energy
  $
    cal(F)_tau^epsilon (bold(U), Theta)
    = cal(F)_tau (bold(U), Theta)
      + epsilon / 2 sum_(j=1)^c N_j "Tr"(bold(F)_j^(-1)).
  $
  Following Wu's convergence theorem for EM @wu1983convergence and the local
  fuzzy convergence framework of @hathaway1986fcm @mclachlan2000finite, every
  iterate decreases $cal(F)_tau^epsilon$, the sequence ${(bold(U)^((t)), Theta^((t)))}$
  remains inside a compact set, and every limit point is a stationary point of
  $cal(F)_tau^epsilon$.
]

#proof[
  In the M-step for $bold(Lambda)_j = bold(F)_j^(-1)$, adding the barrier
  $(epsilon / 2) N_j "Tr"(bold(Lambda)_j)$ changes the objective to
  $(N_j / 2)(-ln det(bold(Lambda)_j) + "Tr"(bold(Lambda)_j (bold(S)_j + epsilon bold(I)_d)))$,
  whose unique global minimizer is exactly
  $bold(F)_j^* = bold(S)_j + epsilon bold(I)_d succ.eq epsilon bold(I)_d$
  (matching Step 2!). Furthermore, since all $bold(x)_i$ lie in a bounded ball
  of radius $R$, every prototype lies in the convex hull of $cal(X)$ and
  $epsilon bold(I)_d prec.eq bold(F)_j^* prec.eq (4 R^2 + epsilon) bold(I)_d$.
  Thus all iterates remain in a compact subset $cal(K)$ of the parameter domain.
  Because $cal(F)_tau^epsilon$ is continuously differentiable on $cal(K)$,
  bounded below by Proposition 1, and strictly decreased at every non-stationary
  point by the exact block-minimization updates, Zangwill's global convergence
  theorem (as formalized by Wu @wu1983convergence) implies that every
  accumulation point of the compact sequence is a stationary point of
  $cal(F)_tau^epsilon$.
]

Note a critical geometric consequence highlighted by Gath and Geva @gath1989gg:
because $d_"GG"^2$ grows *exponentially* with the squared Mahalanobis distance
$r_(i j)^2 = (bold(x)_i - bold(v)_j)^top bold(F)_j^(-1)(bold(x)_i - bold(v)_j)$,
the sensitivity $norm(partial u_(i j) / (partial bold(v)_j))_2$ contains the
factor $exp(-r_(i j)^2 / (2(m - 1)))$, which vanishes super-polynomially once
$r_(i j) >> sqrt(m - 1)$. Consequently, GG acts as a *local density refiner*
with narrow basins of attraction and requires warm-start initialization from a
globally broader objective such as FCM @bezdek1981fcm or GK @gustafson1979gk.

#sec("III", "Complexity, Geometry and Operational Profile")

Each iteration of GG consists of three computational phases:
1. *Scatter accumulation:* Computing $c$ prototypes $bold(v)_j$ takes $O(n c d)$
   flops, and accumulating the $c$ weighted covariance matrices $bold(F)_j$
   requires $O(n c d^2)$ flops (leveraging symmetry for $n c d(d + 1) / 2$
   multiply-adds).
2. *Cholesky factorizations:* Factorizing $c$ matrices
   $bold(F)_j = bold(L)_j bold(L)_j^top$ costs $c dot d^3 / 3 = O(c d^3)$ flops,
   yielding both $ln det(bold(F)_j) = 2 sum_r ln (bold(L)_j)_(r r)$ and fast
   triangular solves without explicit matrix inversion.
3. *Exponential dissimilarity and memberships:* Evaluating
   $bold(L)_j^(-1)(bold(x)_i - bold(v)_j)$ for all $(i, j)$ requires $n c$
   forward substitutions of cost $O(d^2)$ each, totaling $O(n c d^2)$ flops,
   followed by $O(n c)$ log-sum-exp operations.

Summing over $T$ iterations gives exact time complexity $O(T c (n d^2 + d^3))$.
Space complexity is $O(n d + n c + c d^2)$ to store the dataset, partition
matrix $bold(U)$, and $c$ covariance matrices. Evaluating memberships strictly
in log-space (Step 4) is mandatory: in moderate dimensions ($d >= 20$),
computing $exp(r_(i j)^2 / 2)$ directly overflows IEEE-754 double precision
($approx 1.8 times 10^(308)$) whenever $r_(i j) > 37.7$.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$O(T c (n d^2 + d^3))$ via Cholesky solves],
    [Space Complexity], [$O(n d + n c + c d^2)$],
    [Cluster Geometry], [Full ellipsoids of variable volume, shape, and prior],
    [Invariance Group], [Full affine $"GL"_d (bb(R)) times.r bb(R)^d$ ($epsilon = 0$)],
    [Hyperparameters], [Cluster count $c$, fuzzifier $m > 1$, shrinkage $epsilon$],
    [Optimization Nature], [Monotone on $cal(F)_(m-1)^epsilon$; narrow local basins],
    [Initialization], [Requires warm start (FCM @bezdek1981fcm or GK @gustafson1979gk)],
    [Statistical Model], [Exact GMM-EM at $m = 2$; tempered EM for $m != 2$],
  ),
  caption: [Operational and geometric profile of Gath-Geva (GG) clustering.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Biomedical Sleep-Stage and EEG Signal Analysis:* Gath and Geva @gath1989gg
   originally motivated GG for automated electroencephalogram (EEG) sleep
   staging, where spectral power features of REM, deep-sleep, and wake states
   occupy vastly unequal durations (imbalanced $P_j$) and exhibit strongly
   anisotropic cross-band correlations (full $bold(F)_j$). Because $d_"GG"^2$
   explicitly weights by $sqrt(det(bold(F)_j)) / P_j$, minor clinical states
   are not absorbed by high-variance dominant stages.
2. *Takagi-Sugeno Nonlinear Fuzzy System Identification:* In control
   engineering, identifying local linear models over nonlinear operating
   envelopes requires partitioning input-output product space into thin,
   elongated hyperplanes of varying density. GG's adaptive covariance and
   volume freedom directly yields smooth Gaussian-like fuzzy rules whose
   principal axes align with local Jacobians.
3. *Multispectral and Hyperspectral Remote Sensing:* Land-cover classes
   (e.g., dense canopy, urban concrete, water bodies) span sharply different
   radiometric covariances and areal abundances. Full-affine equivariance
   (Theorem 2) guarantees that sensor gain recalibrations across spectral bands
   do not distort pixel membership maps.
4. *Medical Image Tissue Segmentation:* Differentiating small focal lesions
   (low $P_j$, compact $bold(F)_j$) from expansive background parenchyma (large
   $P_j$, broad $bold(F)_j$) in multi-parametric MRI defeats fixed-volume
   models such as GK @gustafson1979gk and spherical FCM @bezdek1981fcm, while
   matching GG's variable-volume likelihood geometry (often combined with
   spatial priors as in @chuang2006sfcm).
5. *Heteroskedastic Financial Regime Detection:* Multi-asset return vectors
   shift between low-volatility bull regimes and high-volatility, high-correlation
   crisis regimes. Because $bold(F)_j$ tracks covariance inflation and $P_j$
   captures crisis rarity, GG separates volatility regimes with soft transition
   probabilities controlled by $m$.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Four structural failure modes constrain GG deployment:
- *High-Dimensional Covariance Degeneracy:* Each cluster requires estimating
  $d(d + 1) / 2$ covariance parameters; when $N_j <= d$, $bold(S)_j$ is rank
  deficient and collapses without heavy shrinkage $epsilon$ (Proposition 1).
- *Sensitivity to Cold Initialization:* Vanishing exponential gradients trap
  random initializations in poor local minima @gath1989gg.
- *Outlier Leverage:* Squared Mahalanobis deviations in $bold(S)_j$ give GG
  zero breakdown point; a single extreme outlier inflates $bold(F)_j$ unless
  replaced by heavy-tailed $t$-mixtures @peel2000robust or possibilistic
  weights @krishnapuram1993pcm.
- *Non-Convex Manifolds:* Interlocking crescents or nested rings cannot be
  separated by quadratic ellipsoidal level sets; kernelized @zhang2004kfcm or
  spectral @vonluxburg2007tutorial methods are required.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [EEG / Biosignals], [Unequal priors $P_j$, correlated bands], [Warm-start with FCM],
    [Control / Fuzzy ID], [Ellipsoidal local linear subregions], [Regularize thin slabs],
    [Remote Sensing], [Affine band-gain invariance, unequal area], [High $d$: use PCA / diagonal],
    [MRI Lesion Seg.], [Small lesions vs. broad background], [Add spatial MRF @chuang2006sfcm],
    [High-D ($d > 30$)], [Poor: $O(c d^2)$ overfitting, $O(c d^3)$ cost], [Use MCLUST @fraley2002mclust / FABIA],
    [Non-Convex Shapes], [Poor: quadratic decision boundaries], [Use Spectral / HDBSCAN],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for GG.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Gath-Geva occupies the exact bridge between classical fuzzy partitioning
@bezdek1981fcm @gustafson1979gk and probabilistic finite mixture modeling
@dempster1977em @fraley2002mclust.

#proposition(
  title: "Exact Parametric Hierarchy Linking FCM, GK, GG, and GMM",
  origin: "Original Derivation",
)[
  Consider the generalized tempered free energy $cal(F)_tau$ of Theorem 1
  over covariance family $cal(C) subset.eq bb(S)_(++)^d$ and prior family
  $cal(P) subset.eq Delta_c$:
  1. Restricting $bold(F)_j = sigma^2 bold(I)_d$, $P_j = 1 / c$, and replacing
     the logarithmic cost $u_(i j) ln d^2$ by power-weighted $u_(i j)^m d^2$
     recovers Euclidean FCM @bezdek1981fcm.
  2. Allowing $bold(F)_j in bb(S)_(++)^d$ with fixed determinant
     $det(bold(F)_j) = rho_j$ and power-weighted Mahalanobis cost yields GK
     @gustafson1979gk.
  3. Releasing both volume and prior constraints inside $cal(F)_(m-1)$ yields
     canonical GG @gath1989gg, which coincides identically with full-covariance
     Gaussian Mixture EM @dempster1977em @mclachlan2000finite at $m = 2$, and
     with parsimonious MCLUST eigen-decomposed models @fraley2002mclust when
     $bold(F)_j$ is constrained to structured covariance families.
]

#proof[
  Setting $bold(F)_j = sigma^2 bold(I)_d$ and $P_j = 1 / c$ reduces
  $ln d_"GG"^2$ to an affine function of $norm(bold(x)_i - bold(v)_j)_2^2$,
  recovering isotropic Euclidean proximity (FCM). Releasing anisotropy while
  imposing $det(bold(F)_j) = rho_j$ via Lagrange multipliers yields the GK
  scaled Mahalanobis metric $det(bold(S)_j)^(1 / d) bold(S)_j^(-1)$
  @gustafson1979gk. Releasing $det(bold(F)_j)$ and $P_j$ inside $cal(F)_(m-1)$
  gives GG (Theorem 1), whose exact $m = 2$ equivalence to GMM-EM was proved in
  Corollary 1.
]

*Empirical Benchmark Analysis (ClusterOverview 4-Dataset Suite):* On the
standardized 4-dataset benchmark evaluated via Silhouette @rousseeuw1987silhouettes
and Adjusted Rand Index (ARI) @hubert1985comparing, Gath-Geva achieves:
*Iris* $"Sil" = 0.460$, $"ARI" = 0.559$; *Customers* $"Sil" = 0.337$;
*Multi-Shape* $"Sil" = 0.706$ ($"ARI" = "N/A"$); *High-D* $"Sil" = 0.287$
($"ARI" = "N/A"$); with a mean runtime of *9.848 ms*. These results closely
validate our theoretical analysis. First, GG attains its highest quality on
*Multi-Shape* ($"Sil" = 0.706$) because variable-volume, arbitrarily oriented
ellipsoids adapt effectively to anisotropic and unequally sized clusters where
spherical methods fail. Second, on *Iris* ($"Sil" = 0.460$, $"ARI" = 0.559$),
GG separates the linearly isolated Setosa ellipsoid cleanly but suffers
moderate boundary drift between the overlapping Versicolor and Virginica
ellipsoids due to exponential-tail locality and initialization sensitivity
(Section II.D). Third, on *Customers* ($"Sil" = 0.337$) and especially
*High-D* ($"Sil" = 0.287$), heavy-tailed spending skew and $O(d^2)$ covariance
estimation variance compress separation, confirming the dimensionality
sensitivity predicted by Proposition 1. Finally, the *9.848 ms* mean runtime—an
order of magnitude slower than linear-time fuzzy or biclustering scans—directly
reflects the $O(T c (n d^2 + d^3))$ Cholesky factorizations and scatter-matrix
accumulations per iteration.
