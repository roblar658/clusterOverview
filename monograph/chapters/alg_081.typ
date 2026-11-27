#import "../template/ieee.typ": *

#alg_entry(
  id: 81,
  name: "Gaussian Mixture Models (GMM)",
  category: "Cat. IV: Distribution, Statistical & Topic",
  year: "1977",
  paradigm: "Parametric Finite Mixture / Expectation-Maximization",
  complexity: $cal(O)(T n K d^2 + T K d^3)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

A Gaussian Mixture Model (GMM) models an observed data set $bold(X) = \{bold(x)_1, dots.h, bold(x)_n\} subset.eq bb(R)^d$ as independent realizations drawn from a convex combination of $K$ multivariate normal densities @dempster1977em @mclachlan2000finite @bishop2006pattern. Unlike distance-only partitional algorithms that impose spherical Voronoi partitions, GMM encodes a full generative probability model in which each cluster possesses its own prior weight, mean location, and anisotropic covariance geometry.

Let $Delta^(K-1) = \{bold(pi) in bb(R)_(>= 0)^K : sum_(k=1)^K pi_k = 1\}$ denote the $(K - 1)$-dimensional probability simplex and let $bb(S)_(succ 0)^d$ denote the open cone of real symmetric positive-definite $d times d$ matrices. The parameter tuple is $bold(theta) = \{(pi_k, bold(mu)_k, bold(Sigma)_k)\}_(k=1)^K$ with $bold(pi) in Delta^(K-1)$, $bold(mu)_k in bb(R)^d$, and $bold(Sigma)_k in bb(S)_(succ 0)^d$. The marginal probability density at $bold(x) in bb(R)^d$ is
$
  p(bold(x) | bold(theta)) = sum_(k=1)^K pi_k cal(N)(bold(x) | bold(mu)_k, bold(Sigma)_k) ,
$
where each component density is
$
  cal(N)(bold(x) | bold(mu)_k, bold(Sigma)_k) &= (2 pi)^(-d / 2) det(bold(Sigma)_k)^(-1 / 2) \
  &quad dot exp(-1 / 2 (bold(x) - bold(mu)_k)^top bold(Sigma)_k^(-1) (bold(x) - bold(mu)_k)) .
$
Under the i.i.d. sampling assumption, maximum likelihood estimation (MLE) seeks to maximize the incomplete-data log-likelihood
$
  ell(bold(theta)) = sum_(i=1)^n log (sum_(k=1)^K pi_k cal(N)(bold(x)_i | bold(mu)_k, bold(Sigma)_k)) .
$
Because the logarithm acts outside the sum over components, $ell(bold(theta))$ is non-concave, couples all parameters through posterior responsibilities, and admits no closed-form global maximizer. To make optimization tractable, one introduces one-hot latent assignment vectors $bold(z)_i in \{0, 1\}^K$ with $sum_(k=1)^K z_(i k) = 1$, $z_(i k) tilde.op "Cat"(bold(pi))$, and $bold(x)_i | (z_(i k) = 1) tilde.op cal(N)(bold(mu)_k, bold(Sigma)_k)$. The complete-data log-likelihood factorizes cleanly as
$
  ell_c (bold(theta); bold(X), bold(Z)) = sum_(i=1)^n sum_(k=1)^K z_(i k) [log pi_k + log cal(N)(bold(x)_i | bold(mu)_k, bold(Sigma)_k)] .
$

#definition(title: "Variational Free-Energy Decomposition of GMM")[
  Let $q(bold(Z)) = product_(i=1)^n q_i (bold(z)_i)$ be any factorized distribution over latent assignments with component probabilities $gamma_(i k) = q_i (z_(i k) = 1) in [0, 1]$ and $sum_(k=1)^K gamma_(i k) = 1$. The variational free energy (evidence lower bound) is
  $
    cal(F)(q, bold(theta)) &= bb(E)_q [ell_c (bold(theta); bold(X), bold(Z))] + H(q) \
    &= sum_(i=1)^n sum_(k=1)^K gamma_(i k) log ((pi_k cal(N)(bold(x)_i | bold(mu)_k, bold(Sigma)_k)) / gamma_(i k)) ,
  $
  where $H(q) = -sum_(i=1)^n sum_(k=1)^K gamma_(i k) log gamma_(i k)$ is the discrete Shannon entropy and $0 log 0 := 0$. For every $q$ and $bold(theta)$, the exact log-likelihood identity holds:
  $
    ell(bold(theta)) = cal(F)(q, bold(theta)) + sum_(i=1)^n D_"KL" (q_i || p(bold(z)_i | bold(x)_i, bold(theta))) .
  $
]

Expectation-Maximization (EM) @dempster1977em performs exact block coordinate ascent on $cal(F)(q, bold(theta))$, alternating an E-step over $q$ and an M-step over $bold(theta)$. In practice, a small isotropic covariance ridge $epsilon bold(I)_d$ ($epsilon > 0$) is added during the M-step to prevent numerical degeneracy.

#algbox(
  title: "Expectation-Maximization (EM) for Full-Covariance GMM",
  inputs: [Data $bold(X) = \{bold(x)_i\}_(i=1)^n subset.eq bb(R)^d$, components $K >= 1$, ridge $epsilon >= 0$, tolerance $tau > 0$, max iterations $T$.],
  outputs: [Parameters $hat(bold(theta)) = \{(hat(pi)_k, hat(bold(mu))_k, hat(bold(Sigma))_k)\}_(k=1)^K$, soft responsibilities $bold(Gamma) in [0, 1]^(n times K)$.],
)[
  1. *Initialize* $bold(theta)^((0))$ via $k$-means++ seeding @arthur2007kmeanspp (or multi-start random partitions) with initial weights $pi_k^((0)) > 0$, means $bold(mu)_k^((0))$, and empirical covariances $bold(Sigma)_k^((0)) + epsilon bold(I)_d$.
  2. *For* $t = 0, 1, dots.h, T - 1$:
  3. #h(1em) *E-step:* For each $i in \{1, dots.h, n\}$ and $k in \{1, dots.h, K\}$, compute log-weighted densities $a_(i k) = log pi_k^((t)) + log cal(N)(bold(x)_i | bold(mu)_k^((t)), bold(Sigma)_k^((t)))$ via Cholesky factors of $bold(Sigma)_k^((t))$, and normalize via log-sum-exp:
     $
       gamma_(i k)^((t)) = exp(a_(i k) - log sum_(j=1)^K exp(a_(i j))) .
     $
  4. #h(1em) *M-step:* For each $k in \{1, dots.h, K\}$, set effective count $N_k^((t)) = sum_(i=1)^n gamma_(i k)^((t))$ and update
     $
       pi_k^((t+1)) &= N_k^((t)) / n , quad bold(mu)_k^((t+1)) = 1 / N_k^((t)) sum_(i=1)^n gamma_(i k)^((t)) bold(x)_i , \
       bold(Sigma)_k^((t+1)) &= 1 / N_k^((t)) sum_(i=1)^n gamma_(i k)^((t)) bold(d)_(i k) bold(d)_(i k)^top + epsilon bold(I)_d ,
     $
     where $bold(d)_(i k) = bold(x)_i - bold(mu)_k^((t+1))$.
  5. #h(1em) *Convergence check:* Evaluate $ell(bold(theta)^((t+1)))$; if $|ell(bold(theta)^((t+1))) - ell(bold(theta)^((t)))| < tau$, *break*.
  6. *Return* $hat(bold(theta)) = bold(theta)^((t+1))$ and $bold(Gamma) = (gamma_(i k))$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish the variational ascent geometry, global unboundedness of unregularized MLE, full affine equivariance, classical convergence and local contraction rates, and the exact zero-temperature reduction to $k$-means.

#subsec("A", "Variational Monotonicity and Exact Ascent Gap")

#lemma(title: "Exact Variational Monotonicity and Strict Ascent", origin: "Original Proof")[
  For unregularized EM ($epsilon = 0$) whenever every $bold(Sigma)_k^((t+1)) in bb(S)_(succ 0)^d$, the iterate sequence $\{bold(theta)^((t))\}$ satisfies the exact lower-bound ascent inequality
  $
    ell(bold(theta)^((t+1))) - ell(bold(theta)^((t))) >= sum_(i=1)^n D_"KL" (q_i^((t)) || p(bold(z)_i | bold(x)_i, bold(theta)^((t+1)))) >= 0 ,
  $
  where $q_i^((t)) = p(bold(z)_i | bold(x)_i, bold(theta)^((t)))$. Moreover, $ell(bold(theta)^((t+1))) = ell(bold(theta)^((t)))$ holds if and only if $bold(theta)^((t))$ is a stationary point of $ell(bold(theta))$ on the open parameter interior.
]

#proof[
  By the variational decomposition of Definition 1, for any $bold(theta)$ and any factorized $q$,
  $
    ell(bold(theta)) = cal(F)(q, bold(theta)) + sum_(i=1)^n D_"KL" (q_i || p(bold(z)_i | bold(x)_i, bold(theta))) .
  $
  *Step 1 (E-step optimality).* Holding $bold(theta)^((t))$ fixed, Gibbs' inequality gives $D_"KL" (q_i || p(bold(z)_i | bold(x)_i, bold(theta)^((t)))) >= 0$ with unique zero at $q_i^((t)) (z_(i k) = 1) = gamma_(i k)^((t)) = p(z_(i k) = 1 | bold(x)_i, bold(theta)^((t)))$. Thus $cal(F)(q^((t)), bold(theta)^((t))) = ell(bold(theta)^((t)))$.

  *Step 2 (M-step concavity and unique closed form).* Fixing $q^((t))$, maximizing $cal(F)(q^((t)), bold(theta))$ decouples across $bold(pi)$ and each component $(bold(mu)_k, bold(Lambda)_k)$ with precision matrix $bold(Lambda)_k = bold(Sigma)_k^(-1) in bb(S)_(succ 0)^d$:
  $
    cal(F)(q^((t)), bold(theta)) &= sum_(k=1)^K N_k^((t)) log pi_k + 1 / 2 sum_(k=1)^K N_k^((t)) log det(bold(Lambda)_k) \
    &quad - 1 / 2 sum_(k=1)^K sum_(i=1)^n gamma_(i k)^((t)) (bold(x)_i - bold(mu)_k)^top bold(Lambda)_k (bold(x)_i - bold(mu)_k) + C ,
  $
  where $C$ is independent of $bold(theta)$. Introducing a Lagrange multiplier $lambda$ for $sum_k pi_k = 1$ gives $N_k^((t)) / pi_k = lambda$, so $pi_k^((t+1)) = N_k^((t)) / n$, which is the unique global maximizer of the strictly concave map $bold(pi) |-> sum_k N_k^((t)) log pi_k$. For fixed $bold(Lambda)_k succ bold(0)$, the quadratic form in $bold(mu)_k$ is strictly concave with gradient $bold(Lambda)_k sum_i gamma_(i k)^((t)) (bold(x)_i - bold(mu)_k) = bold(0)$, uniquely yielding $bold(mu)_k^((t+1)) = (1 / N_k^((t))) sum_i gamma_(i k)^((t)) bold(x)_i$. Substituting $bold(mu)_k^((t+1))$ leaves the matrix objective
  $
    g_k (bold(Lambda)_k) = (N_k^((t)) / 2) [log det(bold(Lambda)_k) - "Tr"(bold(S)_k^((t)) bold(Lambda)_k)] ,
  $
  where $bold(S)_k^((t)) = (1 / N_k^((t))) sum_i gamma_(i k)^((t)) bold(d)_(i k) bold(d)_(i k)^top$. Since $log det(bold(Lambda)_k)$ is strictly concave on $bb(S)_(succ 0)^d$ with matrix gradient $bold(Lambda)_k^(-1)$ and $"Tr"(bold(S)_k^((t)) bold(Lambda)_k)$ is linear with gradient $bold(S)_k^((t))$, stationarity $bold(Lambda)_k^(-1) - bold(S)_k^((t)) = bold(0)$ gives the unique global maximizer $bold(Sigma)_k^((t+1)) = bold(S)_k^((t))$. Hence $cal(F)(q^((t)), bold(theta)^((t+1))) >= cal(F)(q^((t)), bold(theta)^((t)))$.

  *Step 3 (Total gain).* Applying Definition 1 at $bold(theta)^((t+1))$ yields
  $
    ell(bold(theta)^((t+1))) - ell(bold(theta)^((t))) &= [cal(F)(q^((t)), bold(theta)^((t+1))) - cal(F)(q^((t)), bold(theta)^((t)))] \
    &quad + sum_(i=1)^n D_"KL" (q_i^((t)) || p(bold(z)_i | bold(x)_i, bold(theta)^((t+1)))) ,
  $
  proving the stated inequality since the bracketed M-step gain is non-negative. Equality forces both terms to vanish; by strict concavity of the M-step in $(bold(pi), bold(mu)_k, bold(Lambda)_k)$, zero M-step gain forces $bold(theta)^((t+1)) = bold(theta)^((t))$ and hence $nabla_(bold(theta)) cal(F)(q^((t)), bold(theta)^((t))) = bold(0)$. Because $D_"KL" (q_i^((t)) || p(bold(z)_i | bold(x)_i, bold(theta)))$ attains its minimum $0$ at $bold(theta) = bold(theta)^((t))$, its gradient with respect to $bold(theta)$ vanishes there, so $nabla_(bold(theta)) ell(bold(theta)^((t))) = nabla_(bold(theta)) cal(F)(q^((t)), bold(theta)^((t))) = bold(0)$.
]

#subsec("B", "Likelihood Singularities and Spectral Regularization")

#theorem(title: "Unboundedness of Unregularized MLE and Boundedness under Spectral Floors", origin: "Original Counterexample")[
  Let $bold(X) = \{bold(x)_1, dots.h, bold(x)_n\} subset.eq bb(R)^d$ with $n >= 2$ and $K >= 2$.
  1. *Singularity:* The unregularized GMM log-likelihood satisfies $sup_(bold(theta)) ell(bold(theta)) = +infinity$, so a global MLE does not exist in the open parameter space.
  2. *Restoration:* Restricting every covariance matrix to $bold(Sigma)_k succ.eq lambda_min bold(I)_d$ for any fixed $lambda_min > 0$ yields the uniform global upper bound
     $
       ell(bold(theta)) <= -(n d / 2) log(2 pi lambda_min) < +infinity .
     $
]

#proof[
  *Part 1.* Fix $bold(mu)_1 = bold(x)_1$, $bold(Sigma)_1 = sigma^2 bold(I)_d$ with $sigma > 0$, $pi_1 = 1 / K$, and for $k >= 2$ fix arbitrary $pi_k = 1 / K$, $bold(mu)_k = bold(0)$, and $bold(Sigma)_k = bold(I)_d$. Denote this parameter curve by $bold(theta)_sigma$. Dropping non-negative terms in the mixture sum gives, for $i = 1$,
  $
    p(bold(x)_1 | bold(theta)_sigma) >= pi_1 cal(N)(bold(x)_1 | bold(x)_1, sigma^2 bold(I)_d) = (1 / K) (2 pi sigma^2)^(-d / 2) ,
  $
  and for every $i in \{2, dots.h, n\}$,
  $
    p(bold(x)_i | bold(theta)_sigma) >= pi_2 cal(N)(bold(x)_i | bold(0), bold(I)_d) = (1 / K) (2 pi)^(-d / 2) exp(-1 / 2 ||bold(x)_i||_2^2) .
  $
  Summing the logarithms over $i = 1, dots.h, n$ yields
  $
    ell(bold(theta)_sigma) >= -d log sigma - n log K - (n d / 2) log(2 pi) - 1 / 2 sum_(i=2)^n ||bold(x)_i||_2^2 .
  $
  As $sigma -> 0^+$, $-d log sigma -> +infinity$ while all remaining terms are finite constants independent of $sigma$. Hence $sup_(bold(theta)) ell(bold(theta)) = +infinity$.

  *Part 2.* Suppose $bold(Sigma)_k succ.eq lambda_min bold(I)_d$ for all $k in \{1, dots.h, K\}$. By Weyl's eigenvalue inequality, every eigenvalue of $bold(Sigma)_k$ satisfies $lambda_j (bold(Sigma)_k) >= lambda_min > 0$, so $det(bold(Sigma)_k) = product_(j=1)^d lambda_j (bold(Sigma)_k) >= lambda_min^d$. Moreover, since $bold(Sigma)_k^(-1) succ bold(0)$, the quadratic form is non-negative, giving $exp(-1 / 2 (bold(x)_i - bold(mu)_k)^top bold(Sigma)_k^(-1) (bold(x)_i - bold(mu)_k)) <= 1$. Therefore, for every $i$ and $k$,
  $
    cal(N)(bold(x)_i | bold(mu)_k, bold(Sigma)_k) <= (2 pi)^(-d / 2) lambda_min^(-d / 2) = (2 pi lambda_min)^(-d / 2) .
  $
  Because $sum_(k=1)^K pi_k = 1$, taking the convex combination preserves this uniform bound: $p(bold(x)_i | bold(theta)) <= (2 pi lambda_min)^(-d / 2)$ for all $i$. Taking logarithms and summing over $i = 1, dots.h, n$ yields $ell(bold(theta)) <= -(n d / 2) log(2 pi lambda_min)$, proving global boundedness.
]

#subsec("C", "Affine Equivariance and Covariance Restrictions")

#proposition(title: "Full Affine Equivariance of Unregularized Full-Covariance EM", origin: "Original Proof")[
  Let $T(bold(x)) = bold(A) bold(x) + bold(b)$ be an invertible affine map with $bold(A) in "GL"(d, bb(R))$ and $bold(b) in bb(R)^d$, and set $bold(x)'_i = T(bold(x)_i)$. Initialize full-covariance EM ($epsilon = 0$) on $bold(X)'$ at $pi_k'^((0)) = pi_k^((0))$, $bold(mu)_k'^((0)) = bold(A) bold(mu)_k^((0)) + bold(b)$, and $bold(Sigma)_k'^((0)) = bold(A) bold(Sigma)_k^((0)) bold(A)^top$. Then for every iteration $t >= 0$:
  1. $gamma_(i k)'^((t)) = gamma_(i k)^((t))$ for all $i, k$;
  2. $pi_k'^((t)) = pi_k^((t))$, $bold(mu)_k'^((t)) = bold(A) bold(mu)_k^((t)) + bold(b)$, and $bold(Sigma)_k'^((t)) = bold(A) bold(Sigma)_k^((t)) bold(A)^top$;
  3. $ell'(bold(theta)'^((t))) = ell(bold(theta)^((t))) - n log |det(bold(A))|$.
  Conversely, diagonal GMM loses rotation equivariance and spherical GMM loses anisotropic scaling equivariance.
]

#proof[
  We proceed by induction on $t$. The base case $t = 0$ holds for the parameters by hypothesis. Assume the parameter relations hold at iteration $t$. For any $i, k$, the transformed Mahalanobis quadratic form is invariant:
  $
    &(bold(x)'_i - bold(mu)_k'^((t)))^top (bold(Sigma)_k'^((t)))^(-1) (bold(x)'_i - bold(mu)_k'^((t))) \
    &= (bold(x)_i - bold(mu)_k^((t)))^top bold(A)^top bold(A)^(-top) (bold(Sigma)_k^((t)))^(-1) \
    &quad dot bold(A)^(-1) bold(A) (bold(x)_i - bold(mu)_k^((t))) \
    &= (bold(x)_i - bold(mu)_k^((t)))^top (bold(Sigma)_k^((t)))^(-1) (bold(x)_i - bold(mu)_k^((t))) .
  $
  Meanwhile, $det(bold(Sigma)_k'^((t))) = det(bold(A))^2 det(bold(Sigma)_k^((t)))$, so every component density scales by the common positive constant $|det(bold(A))|^(-1)$:
  $
    &cal(N)(bold(x)'_i | bold(mu)_k'^((t)), bold(Sigma)_k'^((t))) \
    &= |det(bold(A))|^(-1) cal(N)(bold(x)_i | bold(mu)_k^((t)), bold(Sigma)_k^((t))) .
  $
  This factor cancels identically between numerator and denominator in the E-step ratio, giving $gamma_(i k)'^((t)) = gamma_(i k)^((t))$, and contributes $-log |det(bold(A))|$ to each point's log-likelihood, yielding claim 3. In the M-step, identical responsibilities give $N_k'^((t)) = N_k^((t))$ and $pi_k'^((t+1)) = pi_k^((t+1))$. By linearity of weighted sums,
  $
    bold(mu)_k'^((t+1)) &= 1 / N_k^((t)) sum_(i=1)^n gamma_(i k)^((t)) (bold(A) bold(x)_i + bold(b)) \
    &= bold(A) bold(mu)_k^((t+1)) + bold(b) , \
    bold(Sigma)_k'^((t+1)) &= 1 / N_k^((t)) sum_(i=1)^n gamma_(i k)^((t)) bold(A) bold(d)_(i k) bold(d)_(i k)^top bold(A)^top \
    &= bold(A) bold(Sigma)_k^((t+1)) bold(A)^top ,
  $
  closing the induction. Finally, if $bold(Sigma)_k$ is constrained to diagonal matrices and $bold(A)$ is a $45^degree$ rotation in $bb(R)^2$, then $bold(A) "diag"(2, 1) bold(A)^top$ has nonzero off-diagonal entries $1 / 2$ and leaves the diagonal cone, breaking rotation equivariance; likewise, if $bold(Sigma)_k = sigma_k^2 bold(I)_d$ and $bold(A) = "diag"(2, 1)$, then $bold(A) bold(I)_2 bold(A)^top = "diag"(4, 1)$ leaves the spherical cone, breaking anisotropic scale equivariance.
]

#subsec("D", "Asymptotic Stationarity and Local Linear Contraction")

#theorem(title: "Stationary Convergence and Statistical Contraction of EM", origin: "Classical — Cited")[
  1. *(Global stationarity @wu1983convergence):* Suppose $bold(theta)^((t))$ stays in a compact subset $Theta_0$ of the regularized parameter space ($bold(Sigma)_k succ.eq lambda_min bold(I)_d$) where $cal(F)$ is continuously differentiable in both arguments. Then every limit point of $\{bold(theta)^((t))\}$ is a stationary point of $ell(bold(theta))$, and $ell(bold(theta)^((t)))$ converges monotonically to a stationary value $ell^*$.
  2. *(Local linear contraction @balakrishnan2017statistical):* Consider a balanced two-component isotropic Gaussian mixture $1 / 2 cal(N)(-bold(mu)^*, sigma^2 bold(I)_d) + 1 / 2 cal(N)(bold(mu)^*, sigma^2 bold(I)_d)$ with signal-to-noise ratio $eta = ||bold(mu)^*||_2 / sigma$. Whenever $eta >= eta_0$ (a universal constant) and the initialization lies in the ball $||bold(mu)^((0)) - bold(mu)^*||_2 <= r_0 ||bold(mu)^*||_2$, population EM contracts geometrically with rate $kappa(eta) <= exp(-c eta^2) < 1$:
     $
       ||bold(mu)^((t)) - bold(mu)^*||_2 <= kappa(eta)^t ||bold(mu)^((0)) - bold(mu)^*||_2 ,
     $
     and sample EM achieves statistical error $cal(O)(sigma sqrt(d / n))$ with high probability after $t = cal(O)(log(n / d))$ iterations.
]

#proof[
  For Part 1, by Lemma 1 the objective $ell(bold(theta)^((t)))$ is non-decreasing and by Theorem 1 it is bounded above on $Theta_0$, hence converges to a finite limit $ell^*$. Compactness of $Theta_0$ guarantees convergent subsequences $bold(theta)^((t_m)) -> bar(bold(theta))$. Continuity of the EM map $M(bold(theta))$ implies $ell(M(bar(bold(theta)))) = lim_(m -> infinity) ell(bold(theta)^((t_m + 1))) = ell^* = ell(bar(bold(theta)))$. By the strict-ascent characterization in Lemma 1, $ell(M(bar(bold(theta)))) = ell(bar(bold(theta)))$ forces $bar(bold(theta))$ to be a stationary point of $ell$, reproducing Wu's global convergence criterion @wu1983convergence.

  For Part 2, linearizing the EM operator $M(bold(theta))$ around a non-degenerate interior local maximizer $bold(theta)^*$ yields the classical missing-information identity @dempster1977em @mclachlan2000finite:
  $
    bold(J)_"EM" (bold(theta)^*) := D M(bold(theta)^*) = bold(I)_"com" (bold(theta)^*)^(-1) bold(I)_"mis" (bold(theta)^*) ,
  $
  where $bold(I)_"com" = -bb(E)[nabla^2 ell_c]$ is the complete-data Fisher information and $bold(I)_"mis" = bold(I)_"com" - bold(I)_"obs"$ is the missing information induced by latent labels $bold(Z)$. Because $bold(0) prec.eq bold(I)_"mis" prec bold(I)_"com"$ whenever $bold(theta)^*$ is identifiable, the spectral radius satisfies $rho(bold(J)_"EM" (bold(theta)^*)) in [0, 1)$, governing asymptotic linear convergence. In the symmetric two-component model, the posterior probability $w(bold(x); bold(mu)) = (1 + exp(-2 chevron.l bold(x), bold(mu) chevron.r / sigma^2))^(-1)$ has derivative proportional to $w(1 - w)$, which concentrates exponentially near $0$ except in a thin slab $|chevron.l bold(x), bold(mu) chevron.r| = cal(O)(sigma ||bold(mu)||_2)$ of Gaussian tail mass $exp(-c eta^2)$. Controlling the empirical process deviation between sample and population EM gradients over the contraction ball yields the minimax-optimal $cal(O)(sigma sqrt(d / n))$ error floor @balakrishnan2017statistical.
]

#subsec("E", "Zero-Temperature Isotropic Reduction to k-Means")

#corollary(title: "Exact Asymptotic Limit to Lloyd's k-Means", origin: "Original Derivation")[
  Fix spherical tied covariances $bold(Sigma)_k = sigma^2 bold(I)_d$ and uniform weights $pi_k = 1 / K$ for all $k$. Assume every data point $bold(x)_i$ has a unique nearest mean $k^*(i) = arg min_k ||bold(x)_i - bold(mu)_k||_2^2$. Then as $sigma -> 0^+$:
  1. $gamma_(i k) -> bb(1)[k = k^*(i)]$ exponentially fast in $1 / sigma^2$;
  2. $-2 sigma^2 cal(F)(q, bold(theta)) - n d sigma^2 log(2 pi sigma^2) - 2 n sigma^2 log K -> sum_(i=1)^n ||bold(x)_i - bold(mu)_(k^*(i))||_2^2$.
  Consequently, the E-step and M-step converge identically to the assignment and centroid-update steps of Lloyd's $k$-means @lloyd1982least.
]

#proof[
  Let $Delta_(i k) = ||bold(x)_i - bold(mu)_k||_2^2 - ||bold(x)_i - bold(mu)_(k^*(i))||_2^2 >= 0$, with $Delta_(i k) > 0$ for all $k != k^*(i)$. Under $bold(Sigma)_k = sigma^2 bold(I)_d$ and $pi_k = 1 / K$, the responsibility takes the Gibbs-Boltzmann form
  $
    gamma_(i k) = exp(-Delta_(i k) / (2 sigma^2)) / (1 + sum_(j != k^*(i)) exp(-Delta_(i j) / (2 sigma^2))) .
  $
  As $sigma -> 0^+$, every term with $Delta_(i j) > 0$ vanishes exponentially, so $gamma_(i, k^*(i)) -> 1$ and $gamma_(i k) -> 0$ for $k != k^*(i)$. Substituting this form into $-2 sigma^2 cal(F)(q, bold(theta))$ gives
  $
    -2 sigma^2 cal(F)(q, bold(theta)) &= sum_(i=1)^n sum_(k=1)^K gamma_(i k) ||bold(x)_i - bold(mu)_k||_2^2 \
    &quad + n d sigma^2 log(2 pi sigma^2) + 2 n sigma^2 log K - 2 sigma^2 H(q) .
  $
  Since $0 <= H(q) <= n log K$, the entropy term $2 sigma^2 H(q) -> 0$ as $sigma -> 0^+$. Subtracting the parameter-independent constants leaves the exact within-cluster sum-of-squares objective, and the M-step mean update $bold(mu)_k = (sum_i gamma_(i k) bold(x)_i) / (sum_i gamma_(i k))$ reduces to the unweighted arithmetic centroid of cluster $k$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

*Computational Complexity and Covariance Parameterizations.* Each full-covariance EM iteration requires $K$ Cholesky factorizations $bold(Sigma)_k = bold(L)_k bold(L)_k^top$ costing $cal(O)(K d^3)$ flops, followed by $n K$ triangular solves $bold(L)_k^(-1) (bold(x)_i - bold(mu)_k)$ to compute Mahalanobis quadratic forms in $cal(O)(n K d^2)$ flops. The M-step accumulates $K$ outer-product scatter matrices in $cal(O)(n K d^2)$ flops via level-3 BLAS GEMM kernels. Over $T$ iterations, total time is $cal(O)(T n K d^2 + T K d^3)$ and working memory is $cal(O)(n K + K d^2)$ for responsibilities and covariances (streaming out-of-core variants reduce responsibility storage to $cal(O)(B K)$ per mini-batch $B$). Restricting $bold(Sigma)_k$ across the parsimonious eigenvalue-decomposition families of MCLUST @fraley2002mclust @biernacki2006mixmod trades flexibility for dimensionality scaling: diagonal covariances reduce cost to $cal(O)(T n K d)$ time and $cal(O)(K d)$ covariance memory, while tied covariance shares one Cholesky factorization $cal(O)(d^3)$ across all $K$ clusters.

*Decision-Boundary Geometry.* The Bayes-optimal MAP assignment $hat(z)_i = arg max_k [log pi_k + log cal(N)(bold(x)_i | bold(mu)_k, bold(Sigma)_k)]$ partitions $bb(R)^d$ along quadratic hypersurfaces:
$
  &(bold(x) - bold(mu)_k)^top bold(Sigma)_k^(-1) (bold(x) - bold(mu)_k) \
  &- (bold(x) - bold(mu)_j)^top bold(Sigma)_j^(-1) (bold(x) - bold(mu)_j) = c_(k j) ,
$
where $c_(k j) = 2 log(pi_k / pi_j) - log(det(bold(Sigma)_k) / det(bold(Sigma)_j))$. Thus full GMM produces ellipsoidal, parabolic, or hyperbolic boundaries capable of enclosing one compact cluster completely inside a wider diffuse cluster—a topology impossible under linear Voronoi cells. Under tied covariance ($bold(Sigma)_k = bold(Sigma)$), the quadratic terms cancel and boundaries collapse to Mahalanobis-orthogonal hyperplanes (linear discriminant geometry).

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Objective], [Maximize $ell(bold(theta)) = sum_i log sum_k pi_k cal(N)(bold(x)_i | bold(mu)_k, bold(Sigma)_k)$],
    [Free Parameters], [$K - 1 + K d + K d (d + 1) / 2$ (full); $K - 1 + 2 K d$ (diag)],
    [Time Complexity], [$cal(O)(T n K d^2 + T K d^3)$ full; $cal(O)(T n K d)$ diagonal],
    [Space Complexity], [$cal(O)(n K + K d^2)$ full; $cal(O)(n K + K d)$ diagonal],
    [Cluster Geometry], [Anisotropic ellipsoids; quadratic Bayes decision surfaces],
    [Equivariance], [Full affine equivariant (full cov); permutation invariant],
    [Breakdown Point], [$0\%$ (a single extreme outlier can hijack a component)],
    [Convergence Rate], [Linear; factor $rho(bold(I)_"com"^(-1) bold(I)_"mis")$ slows under overlap],
  ),
  caption: [Operational and theoretical profile of Gaussian Mixture Models (Alg. 81).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Acoustic Modeling and Speaker Verification (GMM-UBM):* Short-time Mel-Frequency Cepstral Coefficients (MFCCs) of speech phonemes form overlapping multimodal densities with strongly distinct vocal-tract formant variances. Diagonal and full GMMs serve as Universal Background Models whose posterior responsibilities provide smooth likelihood-ratio scores and MAP-adapted speaker supervectors directly reflecting phonetic covariance structure.
2. *Computer Vision Foreground Segmentation and Background Subtraction:* Per-pixel RGB or Lab color trajectories under waving foliage, shadows, and illumination drift follow multimodal densities with correlated chromatic channels. Adaptive online GMMs model stationary background illumination modes via multi-component covariances, cleanly flagging moving foreground pixels via Mahalanobis log-likelihood thresholds.
3. *Quantitative Finance and Multi-Asset Regime Modeling:* Empirical asset return vectors alternate between calm low-volatility regimes and crisis high-volatility, high-correlation regimes. A full-covariance GMM directly captures regime-dependent covariance matrices $bold(Sigma)_k$, modeling both excess kurtosis and dynamic cross-asset correlation shifts without requiring arbitrary threshold rules.
4. *Flow Cytometry and Single-Cell Marker Gating:* Multi-parameter fluorescence cytometers measure $d approx 10$--$40$ continuous surface markers across $n approx 10^5$ cells, where biological cell lineages exhibit correlated receptor co-expression (tilted ellipsoidal clouds). Full-covariance GMM replaces manual 2D polygon gating with affine-equivariant probabilistic population assignment and calibrated boundary uncertainty.
5. *Astronomical Photometric and Kinematic Population Separation:* Stellar surveys in color-magnitude and proper-motion spaces contain superimposed galactic disk, thick-disk, and halo populations with distinct velocity dispersion tensors. GMM estimates anisotropic velocity ellipsoids $bold(Sigma)_k$ and membership probabilities $gamma_(i k)$ directly supporting downstream Bayes-weighted astrophysical catalog inference.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

GMM degrades sharply under four structural violations: (i) *heavy-tailed noise and outliers*—exponentially decaying Gaussian tails give GMM a $0\%$ breakdown point, causing single outliers to attract dedicated singleton components or severely inflate $bold(Sigma)_k$ (remedied by Student-$t$ mixtures @peel2000robust); (ii) *high-dimensional small-sample regimes ($d > n / K$)*—sample covariance matrices $bold(S)_k$ become rank-deficient, producing severe overfitting or singularity without Bayesian Wishart shrinkage (Alg. 82) or low-rank factor analyzers; (iii) *non-convex manifold topologies*—interlocking crescents, rings, or thin filaments require many redundant Gaussian patches to approximate curvature, destroying one-to-one cluster interpretability; and (iv) *asymmetric skewness*—skewed clusters are spuriously split into multiple symmetric Gaussians unless skew-normal mixtures @lin2009skewnormal are employed.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Correlated continuous ellipsoids ($n >> K d^2$)], [Exact MLE of anisotropic $bold(Sigma)_k$; affine equivariant], [Use BIC @schwarz1978bic or DPGMM (Alg. 82) for $K$],
    [High-dim sparse counts / text ($d >> n$)], [Poor; $cal(O)(K d^2)$ covariance overparameterizes], [Use LDA (Alg. 83) or pLSA (Alg. 84)],
    [Heavy-tailed / outlier-rich sensors], [Fragile ($0\%$ breakdown); tail distortion], [Student-$t$ mixtures @peel2000robust or trimmed $k$-means],
    [Non-convex manifolds (rings, moons)], [Fragments curved manifolds into multiple Gaussians], [HDBSCAN @campello2013hdbscan or Spectral @ng2002spectral],
    [Skewed biomedical assays], [Spurious extra modes to fit asymmetry], [Skew-normal mixtures @lin2009skewnormal],
  ),
  caption: [Domain suitability matrix and recommended alternatives for GMM.],
)

#sec("V", "Structural Relations and Empirical Behavior")

#proposition(title: "Bregman Exponential-Family Embedding of GMM", origin: "Original Derivation")[
  Let $bold(t)(bold(x)) = (bold(x), bold(x) bold(x)^top)$ be the quadratic sufficient statistic of the multivariate Gaussian family and let $phi$ be its convex conjugate cumulant dual. Then every component density satisfies $log cal(N)(bold(x) | bold(mu)_k, bold(Sigma)_k) = -D_phi (bold(t)(bold(x)) || bold(eta)_k) + c(bold(x))$, where $D_phi$ is the unique Bregman divergence induced by $phi$ and $bold(eta)_k = (bold(mu)_k, bold(Sigma)_k + bold(mu)_k bold(mu)_k^top)$ is the expectation parameter. Consequently, the GMM M-step is the exact Bregman centroid update in the lifted sufficient-statistic space @banerjee2005clustering.
]

#proof[
  Writing the Gaussian density in canonical exponential-family form with natural parameter $bold(theta)_k = (bold(Sigma)_k^(-1) bold(mu)_k, -1 / 2 bold(Sigma)_k^(-1))$ and log-partition function $A(bold(theta)_k) = 1 / 2 bold(mu)_k^top bold(Sigma)_k^(-1) bold(mu)_k + 1 / 2 log det(bold(Sigma)_k) + (d / 2) log(2 pi)$, the Legendre dual $phi(bold(eta)) = sup_(bold(theta)) (chevron.l bold(theta), bold(eta) chevron.r - A(bold(theta)))$ satisfies $A(bold(theta)_k) + phi(bold(t)(bold(x))) - chevron.l bold(theta)_k, bold(t)(bold(x)) chevron.r = D_phi (bold(t)(bold(x)) || bold(eta)_k)$ whenever $bold(eta)_k = nabla A(bold(theta)_k) = (bold(mu)_k, bold(Sigma)_k + bold(mu)_k bold(mu)_k^top)$. Hence maximizing $sum_i gamma_(i k) log cal(N)(bold(x)_i | bold(mu)_k, bold(Sigma)_k)$ is identical to minimizing $sum_i gamma_(i k) D_phi (bold(t)(bold(x)_i) || bold(eta)_k)$, whose unique minimizer by the Bregman mean property @banerjee2005clustering is $bold(eta)_k^((t+1)) = (1 / N_k^((t))) sum_i gamma_(i k)^((t)) bold(t)(bold(x)_i)$, recovering both the mean and second-moment covariance updates in one algebraic identity.
]

Within the 200-algorithm taxonomy, GMM anchors Category IV as the canonical parametric finite mixture: fixing isotropic covariances and taking $sigma -> 0^+$ recovers $k$-means (Alg. 1, Corollary 1); placing conjugate Dirichlet-process and Normal-Inverse-Wishart priors over $(bold(pi), bold(mu)_k, bold(Sigma)_k)$ yields Bayesian GMM / DPGMM (Alg. 82 @rasmussen2000igmm); adding BIC model selection across parsimonious covariance families yields MCLUST @fraley2002mclust; and lifting the generative decoder through a neural network produces Variational Deep Embedding (VaDE) @jiang2017vade.

*Empirical Benchmark Behavior (ClusterOverview 4-Dataset Suite).* On the standardized four-dataset evaluation, GMM achieves: *Iris* Silhouette $= 0.358$, ARI $= 0.950$; *Customers* Silhouette $= 0.332$; *Multi-Shape* Silhouette $= 0.712$ (ARI N/A); *High-D* Silhouette $= 0.181$ (ARI N/A); and a mean runtime of $1.400$ ms. Three theoretical signatures stand out immediately. First, on *Iris*, GMM attains an exceptional ARI of $0.950$ @hubert1985comparing despite a modest Euclidean Silhouette of $0.358$ @rousseeuw1987silhouettes: the *versicolor* and *virginica* species form elongated, strongly correlated ellipsoids whose Mahalanobis separation is far cleaner than their isotropic Euclidean distance, vividly illustrating the advantage of full covariance modeling (Proposition 1) over spherical metrics. Second, on *High-D*, Silhouette drops sharply to $0.181$, reflecting the $cal(O)(K d^2)$ covariance overparameterization diagnosed in Section III—where Bayesian Wishart shrinkage (DPGMM, Alg. 82: $0.208$) and additive low-rank topic models (LDA, Alg. 83: $0.232$; pLSA, Alg. 84: $0.287$) gain a clear dimensional advantage. Third, its $1.400$ ms mean runtime is the fastest among Algorithms 81--84, confirming that vectorized closed-form Cholesky M-steps converge rapidly when $d$ is moderate.
