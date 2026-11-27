#import "../template/ieee.typ": *

#alg_entry(
  id: 89,
  name: "von Mises-Fisher Mixture Models (vMF)",
  category: "Cat. IV: Distribution, Statistical & Topic",
  year: "2005",
  paradigm: "Directional Exponential-Family Mixture via EM",
  complexity: $cal(O)(T N K d)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Directional and high-dimensional angular data reside naturally on the
$(d - 1)$-dimensional unit hypersphere
$bb(S)^(d-1) = {bold(x) in bb(R)^d : norm(bold(x))_2 = 1}$, where $d >= 2$
and Euclidean magnitude is either uninformative or actively confounding
@banerjee2005vmf. Classical Gaussian mixtures implicitly penalize radial
dispersion across $bb(R)^d$, distorting cluster boundaries when
observations represent $ell_2$-normalized TF-IDF document vectors, unit
surface normals, or normalized neural embeddings. The mixture of
von Mises-Fisher distributions (movMF) resolves this mismatch by modeling
observations directly through a $K$-component directional exponential
family on $bb(S)^(d-1)$ parameterized by unit mean directions
$bold(mu)_k in bb(S)^(d-1)$ and scalar concentrations $kappa_k >= 0$
@banerjee2005vmf @mclachlan2000finite.

#definition(title: "von Mises-Fisher Distribution and Mixture Model")[
  A random unit vector $bold(X) in bb(S)^(d-1)$ follows a
  *von Mises-Fisher distribution*, written
  $bold(X) tilde.op "vMF"(bold(mu), kappa)$, with mean direction
  $bold(mu) in bb(S)^(d-1)$ and concentration $kappa >= 0$, if its density
  with respect to the uniform surface measure $sigma$ on $bb(S)^(d-1)$ is
  $ f_d (bold(x) mid(|) bold(mu), kappa) = C_d (kappa) exp(kappa bold(mu)^top bold(x)), $
  where the normalizing constant is given by
  $ C_d (kappa) = frac(kappa^(d/2 - 1), (2 pi)^(d/2) I_(d/2 - 1) (kappa)), $
  and $I_nu (kappa)$ denotes the modified Bessel function of the first kind
  of order $nu = d / 2 - 1$. The *Bessel ratio map* is defined for
  $kappa >= 0$ by
  $ A_d (kappa) := frac(I_(d/2) (kappa), I_(d/2 - 1) (kappa)) in [0, 1). $
  A $K$-component *movMF* model assigns mixture density
  $ p(bold(x) mid(|) bold(Theta)) = sum_(k=1)^K pi_k f_d (bold(x) mid(|) bold(mu)_k, kappa_k), $
  with parameter tuple $bold(Theta) = {(pi_k, bold(mu)_k, kappa_k)}_(k=1)^K$,
  $bold(pi) in Delta^(K-1)$, $bold(mu)_k in bb(S)^(d-1)$, and
  $kappa_k >= 0$.
]

Given an i.i.d. dataset $cal(X) = {bold(x)_1, dots, bold(x)_N} subset.eq bb(S)^(d-1)$,
maximum-likelihood estimation maximizes the incomplete-data log-likelihood
$ cal(L)(bold(Theta)) = sum_(i=1)^N log (sum_(k=1)^K pi_k C_d (kappa_k) exp(kappa_k bold(mu)_k^top bold(x)_i)). $
Introducing latent component indicators $z_i in {1, dots, K}$ yields the
Expectation-Maximization (EM) scheme of Banerjee et al. @banerjee2005vmf
@dempster1977em. In the E-step, posterior responsibilities
$gamma_(i k) = "Pr"(z_i = k mid(|) bold(x)_i; bold(Theta))$ weight each
unit vector to form the component resultant vector
$bold(r)_k = sum_(i=1)^N gamma_(i k) bold(x)_i$. In the M-step, each mean
direction aligns with $bold(r)_k$ and each concentration solves the
implicit Bessel equation $A_d (kappa_k) = bar(r)_k := norm(bold(r)_k)_2 / N_k$,
where $N_k = sum_(i=1)^N gamma_(i k)$. Because inverting $A_d$ directly
via Bessel evaluations is numerically delicate at large $d$, practical
implementations pair Banerjee's closed-form initialization with one or
two safeguarded Newton steps @banerjee2005vmf.

#algbox(
  title: "Expectation-Maximization for movMF (Soft-movMF)",
  inputs: [Unit vectors $bold(x)_1, dots, bold(x)_N in bb(S)^(d-1)$; cluster count $K$; cap $kappa_"max" > 0$; tolerance $epsilon > 0$; max iterations $T$.],
  outputs: [Parameters $bold(Theta) = {(pi_k, bold(mu)_k, kappa_k)}_(k=1)^K$; soft posteriors $bold(Gamma) in [0, 1]^(N times K)$; hard labels $hat(c)_i = arg max_k gamma_(i k)$.],
)[
  + *Initialize:* Seed $bold(mu)_1^((0)), dots, bold(mu)_K^((0)) in bb(S)^(d-1)$ via spherical $k$-means++ or farthest-angle sampling; set $pi_k^((0)) = 1 / K$ and $kappa_k^((0)) = kappa_0 > 0$.
  + *For* $t = 0, 1, dots, T - 1$ *do:*
  + #h(1em) *E-step:* For every $i in {1, dots, N}$ and $k in {1, dots, K}$, compute log-weights
    $ ell_(i k) = log pi_k^((t)) + log C_d (kappa_k^((t))) + kappa_k^((t)) (bold(mu)_k^((t)))^top bold(x)_i, $
    and stabilize via log-sum-exp to form $gamma_(i k)^((t)) = exp(ell_(i k)) / sum_(j=1)^K exp(ell_(i j))$.
  + #h(1em) *M-step:* For each $k in {1, dots, K}$, set $N_k = sum_(i=1)^N gamma_(i k)^((t))$, $pi_k^((t+1)) = N_k / N$, and $bold(r)_k = sum_(i=1)^N gamma_(i k)^((t)) bold(x)_i$.
  + #h(1em) Compute $bar(r)_k = norm(bold(r)_k)_2 / N_k$, $bold(mu)_k^((t+1)) = bold(r)_k / norm(bold(r)_k)_2$, and solve $A_d (kappa_k^((t+1))) = bar(r)_k$ via Newton refinement from $hat(kappa)_"Ban" = bar(r)_k (d - bar(r)_k^2) / (1 - bar(r)_k^2)$, clipping to $[0, kappa_"max"]$.
  + #h(1em) *Stop* if $cal(L)(bold(Theta)^((t+1))) - cal(L)(bold(Theta)^((t))) < epsilon$.
  + *Return* $bold(Theta)^((t+1))$, $bold(Gamma)^((t))$, and $hat(c)_i = arg max_k gamma_(i k)^((t))$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five structural theorems governing the movMF estimator: the
exact Riccati identity and strict uniqueness of the M-step, monotonic EM
ascent, an original two-sided rational envelope proving $1 / (d - 1)$
relative accuracy for Banerjee's estimator, geometric invariance together
with exact Dirac singularity rates, and the zero-temperature reduction to
spherical $k$-means.

#subsec("A", "Bessel Cumulant Structure and M-Step Uniqueness")

#lemma(title: "Riccati Identity and Strict Monotonicity of the Bessel Ratio", origin: "Original Derivation")[
  Let $d >= 2$, $nu = d / 2 - 1$, and $psi_d (kappa) := -log C_d (kappa)$
  be the log-partition function of $"vMF"(bold(mu), kappa)$ for
  $kappa > 0$. Then:
  (i) $psi'_d (kappa) = bb(E)[bold(mu)^top bold(X)] = A_d (kappa)$;
  (ii) $A_d$ obeys the exact Riccati differential equation
  $ A'_d (kappa) = 1 - A_d (kappa)^2 - frac(d - 1, kappa) A_d (kappa) = "Var"(bold(mu)^top bold(X)) > 0; $
  (iii) $A_d : [0, infinity) -> [0, 1)$ is a strictly increasing
  $C^infinity$ diffeomorphism. Consequently, whenever $N_k > 0$ and
  $0 < bar(r)_k < 1$, the M-step subproblem for component $k$ admits a
  unique global maximizer
  $(hat(bold(mu))_k, hat(kappa)_k) = (bold(r)_k / norm(bold(r)_k)_2, A_d^(-1) (bar(r)_k))$.
]

#proof[
  From the normalization condition
  $exp(psi_d (kappa)) = integral_(bb(S)^(d-1)) exp(kappa bold(mu)^top bold(x)) dif sigma(bold(x))$,
  differentiating under the integral sign yields
  $psi'_d (kappa) = bb(E)[bold(mu)^top bold(X)]$ and
  $psi''_d (kappa) = "Var"(bold(mu)^top bold(X))$. On the other hand,
  $psi_d (kappa) = (nu + 1) log(2 pi) + log I_nu (kappa) - nu log kappa$.
  Recall the standard Bessel differential recurrences
  $ I'_nu (kappa) = I_(nu+1) (kappa) + frac(nu, kappa) I_nu (kappa), \
    I'_(nu+1) (kappa) = I_nu (kappa) - frac(nu + 1, kappa) I_(nu+1) (kappa). $
  Applying the first recurrence directly gives
  $ psi'_d (kappa) = frac(I'_nu (kappa), I_nu (kappa)) - frac(nu, kappa) = frac(I_(nu+1) (kappa), I_nu (kappa)) = A_d (kappa), $
  establishing (i). Differentiating $A_d = I_(nu+1) / I_nu$ and invoking
  both recurrences yields
  $ A'_d (kappa) &= frac(I'_(nu+1) I_nu - I_(nu+1) I'_nu, I_nu^2) \
    &= (1 - frac(nu + 1, kappa) A_d (kappa)) - A_d (kappa) (A_d (kappa) + frac(nu, kappa)) \
    &= 1 - A_d (kappa)^2 - frac(2 nu + 1, kappa) A_d (kappa). $
  Since $2 nu + 1 = d - 1$ and $bold(mu)^top bold(X)$ has continuous
  density supported on $[-1, 1]$ for all $kappa >= 0$, its variance is
  strictly positive, proving (ii). Because $A_d (0) = 0$,
  $lim_(kappa -> infinity) A_d (kappa) = 1$, and $A'_d (kappa) > 0$
  everywhere, $A_d$ is a $C^infinity$ diffeomorphism onto $[0, 1)$.

  For the M-step objective of component $k$,
  $ Q_k (bold(mu)_k, kappa_k) = N_k log C_d (kappa_k) + kappa_k bold(mu)_k^top bold(r)_k. $
  For any fixed $kappa_k > 0$, Cauchy-Schwarz gives
  $bold(mu)_k^top bold(r)_k <= norm(bold(r)_k)_2$ with unique equality at
  $hat(bold(mu))_k = bold(r)_k / norm(bold(r)_k)_2$ whenever
  $norm(bold(r)_k)_2 > 0$. Substituting $hat(bold(mu))_k$ produces the
  univariate objective
  $g(kappa) = N_k (kappa bar(r)_k - psi_d (kappa))$, whose derivatives are
  $g'(kappa) = N_k (bar(r)_k - A_d (kappa))$ and
  $g''(kappa) = -N_k A'_d (kappa) < 0$. Hence $g$ is strictly concave on
  $[0, infinity)$ and uniquely maximized at $hat(kappa)_k = A_d^(-1) (bar(r)_k)$.
]

#subsec("B", "Monotonic Variational Ascent and Convergence")

#theorem(title: "Monotonic Ascent and Stationary Convergence of movMF EM", origin: "Classical — Cited")[
  Let $bold(Theta)^((t))$ be generated by Soft-movMF with concentrations
  restricted to the compact interval $[0, kappa_"max"]$ @banerjee2005vmf
  @dempster1977em @wu1983convergence. Then:
  (i) $cal(L)(bold(Theta)^((t+1))) >= cal(L)(bold(Theta)^((t)))$ for all
  $t >= 0$, with equality if and only if $bold(Theta)^((t))$ is a
  stationary point of $cal(L)$ on the constrained parameter space;
  (ii) every limit point of ${bold(Theta)^((t))}$ is a stationary point of
  $cal(L)$, and $cal(L)(bold(Theta)^((t)))$ converges monotonically to a
  finite value $cal(L)^*$.
]

#proof[
  For any distribution $q_i in Delta^(K-1)$ over $z_i$, decompose the
  log-likelihood via the variational free energy (ELBO):
  $ cal(L)(bold(Theta)) = cal(F)(q, bold(Theta)) + sum_(i=1)^N D_"KL" (q_i || p(dot mid(|) bold(x)_i; bold(Theta))), $
  where
  $cal(F)(q, bold(Theta)) = sum_(i, k) q_(i k) log(pi_k f_d (bold(x)_i mid(|) bold(mu)_k, kappa_k) / q_(i k))$.
  The E-step sets $q_i^((t)) = p(dot mid(|) bold(x)_i; bold(Theta)^((t)))$,
  annihilating the KL term so $cal(L)(bold(Theta)^((t))) = cal(F)(q^((t)), bold(Theta)^((t)))$.
  By Lemma 1 and simplex concavity for $bold(pi)$, the M-step uniquely
  maximizes $cal(F)(q^((t)), bold(Theta))$ over the compact domain
  $Delta^(K-1) times (bb(S)^(d-1) times [0, kappa_"max"])^K$, giving
  $ cal(L)(bold(Theta)^((t+1))) &>= cal(F)(q^((t)), bold(Theta)^((t+1))) \
    &>= cal(F)(q^((t)), bold(Theta)^((t))) = cal(L)(bold(Theta)^((t))), $
  where the first inequality uses nonnegativity of KL divergence. If
  equality holds, both $D_"KL" (q_i^((t)) || q_i^((t+1))) = 0$ and the
  unique M-step maximizer satisfies $bold(Theta)^((t+1)) = bold(Theta)^((t))$,
  which implies zero projected gradient of $cal(L)$ at $bold(Theta)^((t))$.
  Because $kappa_k <= kappa_"max"$ bounds $log C_d (kappa_k) + kappa_k$
  uniformly from above, $cal(L)$ is bounded above on the compact domain;
  by Wu's global convergence theorem @wu1983convergence, all limit points
  of the sequence are stationary points of $cal(L)$.
]

#subsec("C", "Sharp Non-Asymptotic Inverse Bessel Bounds")

Exact evaluation of $A_d^(-1) (bar(r))$ requires numerical root-finding.
We now prove a sharp, non-asymptotic two-sided rational envelope that
simultaneously brackets $A_d^(-1) (bar(r))$ for all $d >= 2$ and certifies
the uniform relative accuracy of Banerjee's approximation @banerjee2005vmf.

#theorem(title: "Exact Two-Sided Rational Envelope for the Inverse Bessel Ratio", origin: "Original Proof")[
  For every dimension $d >= 2$ and every mean resultant length
  $bar(r) in (0, 1)$, the exact concentration MLE
  $kappa^* = A_d^(-1) (bar(r))$ satisfies
  $ frac((d - 1) bar(r), 1 - bar(r)^2) < kappa^* < frac(d bar(r), 1 - bar(r)^2). $
  Furthermore, Banerjee's closed-form estimator
  $hat(kappa)_"Ban" := bar(r) (d - bar(r)^2) / (1 - bar(r)^2)$ lies strictly
  inside this interval and obeys the uniform relative error bound
  $ frac(abs(hat(kappa)_"Ban" - kappa^*), kappa^*) < frac(1, d - 1) quad forall bar(r) in (0, 1), thin d >= 2. $
]

#proof[
  Let $kappa = kappa^* > 0$ so that $A_d (kappa) = bar(r) in (0, 1)$. For
  the lower bound, Lemma 1 gives
  $ 0 < A'_d (kappa) = 1 - bar(r)^2 - frac(d - 1, kappa) bar(r). $
  Multiplying by $kappa / (1 - bar(r)^2) > 0$ immediately yields
  $kappa > (d - 1) bar(r) / (1 - bar(r)^2)$.

  For the upper bound, recall the Weierstrass infinite product for
  $I_nu (kappa)$ ($nu = d / 2 - 1 >= 0$):
  $ kappa^(-nu) I_nu (kappa) = frac(2^(-nu), Gamma(nu + 1)) product_(m=1)^infinity (1 + frac(kappa^2, j_(nu, m)^2)), $
  where $0 < j_(nu, 1) < j_(nu, 2) < dots$ are the positive real zeros of
  the Bessel function $J_nu$. Taking the logarithmic derivative and using
  $frac(dif, dif kappa) log(kappa^(-nu) I_nu (kappa)) = A_d (kappa)$ from
  Lemma 1 gives the uniformly convergent series
  $ A_d (kappa) = sum_(m=1)^infinity frac(2 kappa, kappa^2 + j_(nu, m)^2) ==> frac(A_d (kappa), kappa) = sum_(m=1)^infinity frac(2, kappa^2 + j_(nu, m)^2). $
  Term-by-term differentiation on $(0, infinity)$ yields
  $ frac(dif, dif kappa) (frac(A_d (kappa), kappa)) = -sum_(m=1)^infinity frac(4 kappa, (kappa^2 + j_(nu, m)^2)^2) < 0. $
  Expanding the quotient rule gives
  $kappa A'_d (kappa) - A_d (kappa) < 0$, hence
  $A'_d (kappa) < A_d (kappa) / kappa$. Substituting the Riccati identity
  from Lemma 1 into this inequality yields
  $ 1 - bar(r)^2 - frac(d - 1, kappa) bar(r) < frac(bar(r), kappa) ==> 1 - bar(r)^2 < frac(d bar(r), kappa), $
  which rearranges directly to $kappa < d bar(r) / (1 - bar(r)^2)$.

  Finally, for any $bar(r) in (0, 1)$, we have
  $d - 1 < d - bar(r)^2 < d$, so both $kappa^*$ and $hat(kappa)_"Ban"$
  belong to the open interval $(L, U)$ with
  $L = (d - 1) bar(r) / (1 - bar(r)^2)$ and
  $U = d bar(r) / (1 - bar(r)^2)$. Therefore
  $ frac(abs(hat(kappa)_"Ban" - kappa^*), kappa^*) &< frac(U - L, L) \
    &= frac(bar(r) / (1 - bar(r)^2), (d - 1) bar(r) / (1 - bar(r)^2)) = frac(1, d - 1), $
  completing the proof.
]

Theorem 2 reveals why high-dimensional directional clustering becomes
analytically benign: for $d = 101$, a single closed-form evaluation of
$hat(kappa)_"Ban"$ already guarantees $< 1\%$ relative error across all
concentrations, and one Newton update using Lemma 1's exact derivative
$A'_d$ drives the residual to machine precision.

#subsec("D", "Symmetry Group and Dirac Singularity Rate")

#proposition(title: "Orthogonal Invariance and Single-Point Likelihood Singularity", origin: "Original Proof")[
  (i) For any orthogonal matrix $bold(Q) in O(d)$ ($bold(Q)^top bold(Q) = bold(I)_d$)
  and any positive sample scalars $c_1, dots, c_N > 0$, mapping raw vectors
  $bold(u)_i |-> c_i bold(Q) bold(u)_i$ before unit projection
  $bold(x)_i = bold(u)_i / norm(bold(u)_i)_2$ leaves the movMF
  log-likelihood and posterior matrix $bold(Gamma)$ strictly invariant
  under $bold(mu)_k |-> bold(Q) bold(mu)_k$, whereas translation
  $bold(u)_i |-> bold(u)_i + bold(b)$ ($bold(b) != bold(0)$) breaks
  invariance.
  (ii) If $K >= 2$ and $kappa_k$ is unconstrained, the movMF
  log-likelihood is unbounded above: fixing $bold(mu)_1 = bold(x)_1$,
  $pi_1 = 1 / N$, and sending $kappa_1 -> infinity$ diverges at the exact
  asymptotic rate
  $ cal(L)(bold(Theta)) = frac(d - 1, 2) log(frac(kappa_1, 2 pi)) + cal(O)(1). $
]

#proof[
  For (i), $(c_i bold(Q) bold(u)_i) / norm(c_i bold(Q) bold(u)_i)_2 = bold(Q) bold(x)_i$
  because $c_i > 0$ cancels and $norm(bold(Q) bold(u)_i)_2 = norm(bold(u)_i)_2$.
  Inner products satisfy $(bold(Q) bold(mu)_k)^top (bold(Q) bold(x)_i) = bold(mu)_k^top bold(x)_i$,
  leaving every density $f_d$ and responsibility $gamma_(i k)$ unchanged.
  Translating by $bold(b) != bold(0)$ shifts the origin and alters pairwise
  angles $(bold(u)_i + bold(b))^top (bold(u)_j + bold(b)) / (norm(bold(u)_i + bold(b))_2 norm(bold(u)_j + bold(b))_2)$.

  For (ii), keep components $k >= 2$ fixed with $kappa_k = 0$ (uniform
  density $C_d (0) > 0$) and $pi_k = (1 - 1 / N) / (K - 1) > 0$, while
  setting $bold(mu)_1 = bold(x)_1$ and $pi_1 = 1 / N$. Every term $i >= 2$
  is bounded below by $log((1 - 1 / N) C_d (0)) > -infinity$, while term
  $i = 1$ satisfies
  $ log p(bold(x)_1 mid(|) bold(Theta)) >= log(frac(1, N)) + log(C_d (kappa_1) exp(kappa_1)). $
  Using the large-argument Bessel asymptotic
  $I_nu (kappa) = (2 pi kappa)^(-1/2) exp(kappa) (1 + cal(O)(kappa^(-1)))$
  for $nu = d / 2 - 1$, we obtain
  $ C_d (kappa_1) exp(kappa_1) &= frac(kappa_1^(d/2 - 1) sqrt(2 pi kappa_1), (2 pi)^(d/2)) (1 + cal(O)(kappa_1^(-1))) \
    &= (frac(kappa_1, 2 pi))^((d-1)/2) (1 + cal(O)(kappa_1^(-1))). $
  Taking logarithms yields
  $cal(L)(bold(Theta)) = ((d - 1) / 2) log(kappa_1 / (2 pi)) + cal(O)(1) -> +infinity$,
  analogous to the covariance singularity of unconstrained Gaussian
  mixtures @mclachlan2000finite.
]

#subsec("E", "Zero-Temperature Equivalence to Spherical k-Means")

#proposition(title: "Small-Variance Limit and Spherical k-Means Reduction", origin: "Classical — Cited")[
  Consider the restricted movMF family with uniform mixing weights
  $pi_k = 1 / K$ and tied concentration $kappa_k = kappa > 0$
  @banerjee2005vmf @dhillon2001concept. As $kappa -> infinity$:
  (i) for every $bold(x)_i$ with a unique closest prototype
  $c_i^* = arg max_k bold(mu)_k^top bold(x)_i$, the E-step responsibilities
  converge exponentially to the hard assignment
  $gamma_(i k) -> bb(I)(k = c_i^*)$;
  (ii) the scaled log-likelihood converges uniformly to the spherical
  $k$-means objective:
  $ lim_(kappa -> infinity) frac(1, kappa) (cal(L)(bold(Theta)) - N log(frac(C_d (kappa), K))) \
    = sum_(i=1)^N max_(1 <= k <= K) bold(mu)_k^top bold(x)_i; $
  (iii) the M-step direction update reduces to the spherical $k$-means
  centroid $bold(mu)_k = (sum_(i in C_k) bold(x)_i) / norm(sum_(i in C_k) bold(x)_i)_2$.
]

#proof[
  Let $m_i = max_k bold(mu)_k^top bold(x)_i$ and
  $Delta_i = m_i - max_(k != c_i^*) bold(mu)_k^top bold(x)_i > 0$. Because
  $pi_k C_d (kappa)$ cancels identically across components,
  $ gamma_(i k) = frac(exp(kappa (bold(mu)_k^top bold(x)_i - m_i)), sum_(j=1)^K exp(kappa (bold(mu)_j^top bold(x)_i - m_i))). $
  For $k != c_i^*$, $gamma_(i k) <= exp(-kappa Delta_i) -> 0$, so
  $gamma_(i c_i^*) -> 1$. For (ii), observe that
  $ exp(kappa m_i) <= sum_(k=1)^K exp(kappa bold(mu)_k^top bold(x)_i) <= K exp(kappa m_i). $
  Taking logarithms, subtracting $log(C_d (kappa) / K)$ per sample, and
  dividing by $kappa$ yields
  $ sum_(i=1)^N m_i &<= frac(1, kappa) (cal(L)(bold(Theta)) - N log(frac(C_d (kappa), K))) \
    &<= sum_(i=1)^N m_i + frac(N log K, kappa). $
  Letting $kappa -> infinity$ sends $N log K / kappa -> 0$, establishing
  (ii). Substituting $gamma_(i k) -> bb(I)(i in C_k)$ into
  $bold(r)_k = sum_i gamma_(i k) bold(x)_i$ immediately yields (iii).
]

#sec("III", "Complexity, Geometry and Operational Profile")

Each Soft-movMF EM iteration executes two dense matrix operations: the
E-step inner-product matrix $bold(X) bold(M)^top in bb(R)^(N times K)$
and the M-step resultant matrix $bold(Gamma)^top bold(X) in bb(R)^(K times d)$,
both costing exactly $2 N K d$ floating-point operations via BLAS-3 GEMM.
Computing $log C_d (kappa_k)$ and inverting $A_d (bar(r)_k)$ via two
Newton iterations costs $cal(O)(K)$ per step, yielding total time
$cal(O)(T N K d)$ for dense data and $cal(O)(T K "nnz"(bold(X)))$ for
sparse document-term matrices. Memory complexity is $cal(O)(N d + K d + N K)$
(or $cal(O)("nnz"(bold(X)) + K d + N K)$ sparse), and because all
prototype comparisons depend only on $bold(mu)_k^top bold(x)_i$,
cache-blocked matrix multiplication achieves near-peak memory bandwidth.

Geometrically, level sets of $f_d (bold(x) mid(|) bold(mu)_k, kappa_k)$
are *spherical caps* centered at $bold(mu)_k$, and posterior decision
boundaries between components $k$ and $l$ are intersections of
$bb(S)^(d-1)$ with affine hyperplanes
$(kappa_k bold(mu)_k - kappa_l bold(mu)_l)^top bold(x) = log(pi_l C_d (kappa_l) / (pi_k C_d (kappa_k)))$.
Crucially, movMF avoids the $cal(O)(K d^2)$ covariance parameterization
of full GMMs @fraley2002mclust while retaining component-specific angular
spreads $kappa_k$, requiring only $K(d + 1) - 1$ free parameters. In high
dimensions ($d >> 1$), isotropic vMF components naturally mitigate the
Euclidean distance-concentration effect @beyer1999nearest whenever
discriminative signal resides in cosine alignment rather than radial norm.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [*Property*], [*Characterization*],
    [Time complexity], [$cal(O)(T N K d)$ dense; $cal(O)(T K "nnz"(bold(X)))$ sparse],
    [Space complexity], [$cal(O)(N d + K d + N K)$ dense; streaming-friendly],
    [Free parameters], [$K(d + 1) - 1$ (linear in $d$, unlike $cal(O)(K d^2)$ GMM)],
    [Cluster geometry], [Isotropic spherical caps on $bb(S)^(d-1)$; hyperplane boundaries],
    [Invariance group], [Orthogonal group $O(d)$ and positive radial rescaling],
    [Singularity control], [Requires cap $kappa <= kappa_"max"$ or conjugate prior],
    [Hard / soft output], [Exact posterior probabilities $gamma_(i k)$ plus hard MAP],
    [Model selection], [BIC / AIC @schwarz1978bic or held-out directional log-likelihood],
  ),
  caption: [Operational and geometric characterization of movMF clustering.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

+ *High-Dimensional Text and Corpus Clustering:* Sparse TF-IDF document
  representations live in the positive orthant of $bb(R)^d$ with $d ~ 10^4$
  to $10^6$, where document length scales vector norms without altering
  topic identity @banerjee2005vmf @steinbach2000comparison. Projecting to
  $bb(S)^(d-1)$ removes length bias exactly (Proposition 1), while
  per-cluster concentrations $kappa_k$ distinguish tightly focused
  sub-topics from broad thematic clusters without estimating $d times d$
  covariance matrices.
+ *Normalized Deep Representation Clustering:* Contrastive and metric
  learning encoders (e.g., SimCLR, CLIP, ArcFace) explicitly map inputs
  onto the hypersphere $bb(S)^(d-1)$ via cosine temperature objectives.
  Because the induced feature geometry matches the vMF natural parameter
  $kappa bold(mu)^top bold(x)$, movMF acts as the canonical probabilistic
  head for soft clustering, out-of-distribution scoring, and class-prior
  calibration on frozen embeddings.
+ *Bioinformatics and Single-Cell Expression Profiling:* Microarray and
  normalized RNA-seq profiles frequently measure relative co-expression
  angles across thousands of genes where technical library-size factors
  induce multiplicative scaling across samples @banerjee2005vmf. Unit
  normalization cancels sample-specific amplification, and vMF mixtures
  isolate coherent gene-regulation modules with calibrated uncertainty.
+ *Computer Vision and Surface Normal Segmentation:* 3D LiDAR point
  clouds, photometric stereo reconstructions, and optical-flow fields
  directly produce unit directional vectors on $bb(S)^2$. Modeling surface
  normals with a movMF mixture separates planar architectural facets
  (high $kappa_k$) from curved or noisy organic surfaces (low $kappa_k$)
  with exact $O(3)$ rotational equivariance.
+ *Geophysics, Paleomagnetism, and Directional Astrometry:* Remanent
  magnetization vectors in rock strata, wind-direction roses, and
  celestial-sphere event catalogs reside intrinsically on $bb(S)^2$.
  Here the vMF distribution is the maximum-entropy directional law given
  mean resultant vector $bold(r)_k$, supplying physically interpretable
  dispersion parameters $kappa_k$ and exact likelihood-ratio tests for
  multimodal geomagnetic polarity regimes.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

The principal failure regimes of movMF follow directly from its
inductive bias on $bb(S)^(d-1)$:
(1) *Radial magnitude carries class identity:* Concentric shells or
amplitude-separated clusters collapse onto overlapping directions after
$ell_2$-normalization.
(2) *Anisotropic elliptical arcs on the sphere:* Because vMF densities
are rotationally symmetric around $bold(mu)_k$, elongated great-circle
ridges or girdle distributions require Kent (Fisher-Bingham) mixtures or
multiple overlapping vMF components.
(3) *Near-origin and zero vectors:* Points near $bold(0)$ suffer extreme
angular amplification under $bold(x) / norm(bold(x))_2$, and $bold(0)$ is
undefined on $bb(S)^(d-1)$, demanding centering/offset care.
(4) *Component collapse and numerical overflow:* Uncapped $kappa_k -> infinity$
triggers the Dirac singularity (Proposition 1), while naive evaluation of
$I_nu (kappa)$ overflows double precision near $kappa > 700$, requiring
exponentially scaled Bessel routines $e^(-kappa) I_nu (kappa)$ or
asymptotic log-Bessel approximations.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [*Domain*], [*Mathematical Fit*], [*Caution / Alternative*],
    [TF-IDF / neural embeddings], [Exact cosine natural parameter; $cal(O)(K d)$ scale], [Choose $K$ via BIC; HDP (Alg. 91) if $K$ unknown],
    [3D normals & paleomagnetism], [Intrinsic $bb(S)^2$ geometry; $O(3)$ invariant], [Girdle/elliptic arcs: use Kent / Bingham],
    [Concentric / radial regimes], [Poor: normalization erases $norm(bold(u))_2$], [Use GMM or spectral clustering],
    [Near-zero noisy tabular data], [High angular variance near origin], [Center carefully or use Gaussian / $t$-mixture],
    [Nonnegative additive parts], [Captures angles only, not part sums], [Use NMF (Alg. 92) or LDA / HDP],
  ),
  caption: [Domain suitability and failure-mode matrix for movMF.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Within the Bregman-divergence taxonomy of clustering @banerjee2005clustering,
movMF sits at the intersection of exponential-family mixtures and
spherical partitional algorithms. Specifically, because
$norm(bold(x) - bold(mu))_2^2 = 2(1 - bold(mu)^top bold(x))$ on
$bb(S)^(d-1)$, an isotropic Gaussian mixture restricted to unit-norm data
induces the exact same log-density functional form up to the
normalization constant: $C_d (kappa)$ replaces $(2 pi sigma^2)^(-d/2)$ with
$kappa = 1 / sigma^2$, making movMF the proper intrinsic spherical
counterpart of isotropic GMM (as well as the soft probabilistic parent of
spherical $k$-means @dhillon2001concept by Proposition 2).

#proposition(title: "Extrinsic Gaussian vs. Intrinsic vMF Normalization Gap", origin: "Original Derivation")[
  For unit vectors $bold(x) in bb(S)^(d-1)$ and tied variance
  $sigma_k^2 = 1 / kappa_k$, an isotropic Gaussian density
  $g(bold(x) mid(|) bold(mu)_k, sigma_k^2) = (2 pi / kappa_k)^(-d/2) exp(-(kappa_k / 2) norm(bold(x) - bold(mu)_k)_2^2)$
  differs from $f_d (bold(x) mid(|) bold(mu)_k, kappa_k)$ by the exact
  data-independent factor
  $ frac(f_d (bold(x) mid(|) bold(mu)_k, kappa_k), g(bold(x) mid(|) bold(mu)_k, 1 / kappa_k)) = frac(kappa_k^(-1) exp(kappa_k), I_(d/2 - 1) (kappa_k)). $
  Consequently, when concentrations $kappa_k$ are unequal across
  components, isotropic GMM on $bb(S)^(d-1)$ biases posterior odds by
  $cal(O)(sqrt(kappa_k / kappa_l))$ relative to the true intrinsic vMF
  likelihood.
]

#proof[
  Expanding $-(kappa_k / 2) norm(bold(x) - bold(mu)_k)_2^2 = -kappa_k + kappa_k bold(mu)_k^top bold(x)$
  yields $g(bold(x) mid(|) bold(mu)_k, 1 / kappa_k) = (kappa_k / (2 pi))^(d/2) exp(-kappa_k) exp(kappa_k bold(mu)_k^top bold(x))$.
  Dividing $f_d (bold(x) mid(|) bold(mu)_k, kappa_k) = kappa_k^(d/2 - 1) ((2 pi)^(d/2) I_(d/2 - 1) (kappa_k))^(-1) exp(kappa_k bold(mu)_k^top bold(x))$
  by $g$ cancels $(2 pi)^(d/2)$ and $exp(kappa_k bold(mu)_k^top bold(x))$,
  leaving $kappa_k^(-1) exp(kappa_k) / I_(d/2 - 1) (kappa_k)$. Since
  $I_(d/2 - 1) (kappa) ~ (2 pi kappa)^(-1/2) exp(kappa)$ for large $kappa$,
  this ratio scales as $sqrt(2 pi / kappa_k)$, biasing posterior odds
  between components $k$ and $l$ by $sqrt(kappa_l / kappa_k)$ unless the
  intrinsic Bessel normalization $C_d (kappa_k)$ is used.
]

In the ClusterOverview 4-dataset empirical benchmark (evaluated via
Silhouette @rousseeuw1987silhouettes and Adjusted Rand Index
@hubert1985comparing), movMF attains strong directional separation with
very low computational overhead (mean runtime *2.171 ms*):
- *Iris:* Silhouette $= 0.443$, $"ARI" = 0.577$. Among the four
  statistical models in this block (Algs. 89--92), movMF achieves the
  highest Iris Silhouette ($0.443$) because projecting centered/normalized
  sepal-petal profiles onto $bb(S)^3$ sharply separates _Setosa_ and forms
  compact angular cones, while losing partial radial length cues separating
  _Versicolor_ from _Virginica_ (hence moderate $"ARI" = 0.577$).
- *Customers:* Silhouette $= 0.351$. Wholesale spending profiles exhibit
  heavy multiplicative scale variation across business sizes; unit
  projection removes total budget scale and clusters buyers by relative
  category spend mix, outperforming discrete LCA ($0.182$) and NMF ($0.266$).
- *Multi-Shape:* Silhouette $= 0.709$ ($"ARI" =$ N/A). Well-separated
  lobes after angular embedding yield clean spherical-cap partitions.
- *High-D:* Silhouette $= 0.299$ ($"ARI" =$ N/A). Linear parameter scaling
  $cal(O)(K d)$ avoids covariance ill-conditioning and sustains stable
  directional cluster quality in high ambient dimensions.

Overall, movMF provides the mathematically canonical, linear-in-$d$,
$O(d)$-invariant probabilistic clustering model whenever angular
orientation supersedes Euclidean magnitude @banerjee2005vmf
@balakrishnan2017statistical.
