#import "../template/ieee.typ": *

#alg_entry(
  id: 82,
  name: "Bayesian Gaussian Mixture Models (DPGMM)",
  category: "Cat. IV: Distribution, Statistical & Topic",
  year: "2000",
  paradigm: "Bayesian Nonparametric / Dirichlet Process Variational Inference",
  complexity: $cal(O)(T n K d^2 + T K d^3)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Finite Gaussian Mixture Models (Alg. 81) suffer from two fundamental statistical limitations: the number of components $K$ must be specified a priori (or selected via post-hoc criteria such as BIC @schwarz1978bic @fraley2002mclust), and unregularized maximum likelihood suffers from degenerate covariance singularities where $ell(bold(theta)) -> +infinity$. Bayesian Gaussian Mixture Models with a Dirichlet Process prior—commonly termed Dirichlet Process Gaussian Mixture Models (DPGMM) or Infinite GMMs @rasmussen2000igmm @blei2006variational—resolve both issues simultaneously within a coherent Bayesian nonparametric framework.

Given observations $bold(X) = \{bold(x)_1, dots.h, bold(x)_n\} subset.eq bb(R)^d$, DPGMM places a Dirichlet Process prior $G tilde.op "DP"(alpha_0, G_0)$ over the latent mixing measure, where $alpha_0 > 0$ is the concentration parameter and $G_0$ is a conjugate Normal-Wishart base measure over component means $bold(mu)_k in bb(R)^d$ and precision matrices $bold(Lambda)_k = bold(Sigma)_k^(-1) in bb(S)_(succ 0)^d$:
$
  G_0 (bold(mu), bold(Lambda)) &= cal(N)(bold(mu) | bold(m)_0, (beta_0 bold(Lambda))^(-1)) cal(W)(bold(Lambda) | bold(W)_0, nu_0) ,
$
with prior mean $bold(m)_0 in bb(R)^d$, precision scaling $beta_0 > 0$, positive-definite scale matrix $bold(W)_0 in bb(S)_(succ 0)^d$, and degrees of freedom $nu_0 > d - 1$. Under Sethuraman's stick-breaking representation @blei2006variational @teh2006hdp, a draw $G = sum_(k=1)^infinity pi_k (bold(v)) delta_((bold(mu)_k, bold(Lambda)_k))$ is constructed via independent beta sticks $v_k tilde.op "Beta"(1, alpha_0)$:
$
  pi_k (bold(v)) = v_k product_(j=1)^(k-1) (1 - v_j) , quad k = 1, 2, dots.h ,
$
satisfying $sum_(k=1)^infinity pi_k (bold(v)) = 1$ almost surely. Each observation selects $z_i tilde.op "Cat"(bold(pi)(bold(v)))$ and is generated as $bold(x)_i | (z_i = k) tilde.op cal(N)(bold(mu)_k, bold(Lambda)_k^(-1))$. Because exact posterior integration over the infinite mixture is intractable, practical DPGMM inference employs truncated mean-field Coordinate Ascent Variational Inference (CAVI) @blei2006variational @bishop2006pattern at a conservative upper truncation level $K$, fixing $v_K = 1$ so that $pi_k (bold(v)) = 0$ for all $k > K$.

#definition(title: "Truncated Stick-Breaking Mean-Field ELBO")[
  Let $bold(Theta) = (bold(v), bold(mu), bold(Lambda))$ collect the truncated stick variables $\{v_k\}_(k=1)^(K-1)$ and component parameters $\{(bold(mu)_k, bold(Lambda)_k)\}_(k=1)^K$. Under the mean-field family
  $
    q(bold(Z), bold(Theta)) = (product_(i=1)^n q(bold(z)_i)) (product_(k=1)^(K-1) q(v_k)) (product_(k=1)^K q(bold(mu)_k, bold(Lambda)_k)) ,
  $
  the Evidence Lower Bound (ELBO) is defined as
  $
    cal(L)(q) &= bb(E)_q [log p(bold(X), bold(Z), bold(Theta))] - bb(E)_q [log q(bold(Z), bold(Theta))] \
    &= log p(bold(X)) - D_"KL" (q(bold(Z), bold(Theta)) || p(bold(Z), bold(Theta) | bold(X))) .
  $
  Maximizing $cal(L)(q)$ is strictly equivalent to minimizing the reverse KL divergence from $q$ to the exact nonparametric posterior.
]

#algbox(
  title: "Coordinate Ascent Variational Inference (CAVI) for DPGMM",
  inputs: [Data $bold(X) subset.eq bb(R)^d$, truncation $K$, concentration $alpha_0 > 0$, NIW prior $(bold(m)_0, beta_0, bold(W)_0, nu_0)$, tolerance $tau > 0$.],
  outputs: [Variational factors $\{q(v_k), q(bold(mu)_k, bold(Lambda)_k)\}_(k=1)^K$, soft assignments $bold(R) = (r_(i k))$, active count $K_"eff"$.],
)[
  1. *Initialize* responsibilities $bold(R)^((0))$ (e.g., via $k$-means or random Dirichlet draws) and run one initial M-step to obtain $(gamma_(k,1), gamma_(k,2), beta_k, bold(m)_k, nu_k, bold(W)_k)$.
  2. *For* $t = 0, 1, dots.h, T - 1$:
  3. #h(1em) *Variational E-step:* Using digamma $psi(x) = (dif / dif x) log Gamma(x)$, compute expected log-weights and log-determinants:
     $
       tilde(omega)_k &= psi(gamma_(k,1)) - psi(gamma_(k,1) + gamma_(k,2)) \
       &quad + sum_(j=1)^(k-1) [psi(gamma_(j,2)) - psi(gamma_(j,1) + gamma_(j,2))] , \
       tilde(Lambda)_k &= sum_(j=1)^d psi((nu_k + 1 - j) / 2) + d log 2 + log det(bold(W)_k) .
     $
     For all $i, k$, evaluate expected squared Mahalanobis distance $E_(i k) = d / beta_k + nu_k (bold(x)_i - bold(m)_k)^top bold(W)_k (bold(x)_i - bold(m)_k)$ and normalize responsibilities $r_(i k) prop exp(tilde(omega)_k + 1 / 2 tilde(Lambda)_k - 1 / 2 E_(i k))$.
  4. #h(1em) *Variational M-step:* Compute $N_k = sum_(i=1)^n r_(i k)$, $bar(bold(x))_k = N_k^(-1) sum_i r_(i k) bold(x)_i$, and scatter $bold(S)_k = N_k^(-1) sum_i r_(i k) (bold(x)_i - bar(bold(x))_k)(bold(x)_i - bar(bold(x))_k)^top$. Update stick and Normal-Wishart factors:
     $
       gamma_(k,1) &= 1 + N_k , quad gamma_(k,2) = alpha_0 + sum_(j=k+1)^K N_j , \
       beta_k &= beta_0 + N_k , quad nu_k = nu_0 + N_k , \
       bold(m)_k &= (beta_0 bold(m)_0 + N_k bar(bold(x))_k) / beta_k , \
       bold(W)_k^(-1) &= bold(W)_0^(-1) + N_k bold(S)_k + (beta_0 N_k) / beta_k bold(Delta)_k bold(Delta)_k^top ,
     $
     where $bold(Delta)_k = bar(bold(x))_k - bold(m)_0$.
  5. #h(1em) *Stop* if $|cal(L)(q^((t+1))) - cal(L)(q^((t)))| < tau$.
  6. *Prune* inactive components with $N_k < epsilon_"min"$ to report $K_"eff"$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We prove conjugate optimality and monotonic ELBO ascent, strict singularity immunity via Loewner matrix bounds, automatic Occam pruning and logarithmic CRP cluster growth, exponential stick-breaking truncation decay, and the multivariate Student-$t$ predictive tail property.

#subsec("A", "Exponential-Family Conjugacy and Exact ELBO Ascent")

#theorem(title: "Closed-Form Conjugate Optimality and Strict CAVI Ascent", origin: "Original Proof")[
  Every block coordinate update in Algorithm 82 preserves the exact conjugate family $q^*(v_k) = "Beta"(gamma_(k,1), gamma_(k,2))$, $q^*(bold(mu)_k, bold(Lambda)_k) = cal(N)(bold(m)_k, (beta_k bold(Lambda)_k)^(-1)) cal(W)(bold(W)_k, nu_k)$, and $q^*(bold(z)_i) = "Cat"(bold(r)_i)$. Furthermore, updating any block $q_j$ from $q_j^((t))$ to $q_j^((t+1))$ increases the ELBO by the exact nonnegative gain
  $
    cal(L)(q^((t+1))) - cal(L)(q^((t))) = D_"KL" (q_j^((t)) || q_j^((t+1))) >= 0 ,
  $
  with equality if and only if $q_j^((t+1)) = q_j^((t))$ almost everywhere.
]

#proof[
  Isolate any single factor $q_j (bold(u)_j)$ while holding all remaining factors $q_(-j)$ fixed. Decompose the ELBO as
  $
    cal(L)(q) &= integral q_j (bold(u)_j) log (tilde(p)_j (bold(u)_j) / q_j (bold(u)_j)) dif bold(u)_j \
    &quad + C_(-j) ,
  $
  where $log tilde(p)_j (bold(u)_j) := bb(E)_(q_(-j)) [log p(bold(X), bold(Z), bold(Theta))] - log Z_j$ defines a normalized density $tilde(p)_j$ and $C_(-j)$ is independent of $q_j$. This rewrites the ELBO cleanly as
  $
    cal(L)(q_j, q_(-j)) = log Z_j + C_(-j) - D_"KL" (q_j || tilde(p)_j) .
  $
  By Gibbs' inequality, the unique maximizer is $q_j^((t+1)) = tilde(p)_j$, at which $D_"KL" (q_j^((t+1)) || tilde(p)_j) = 0$. Consequently, the exact step gain equals
  $
    cal(L)(q_j^((t+1)), q_(-j)) - cal(L)(q_j^((t)), q_(-j)) = D_"KL" (q_j^((t)) || q_j^((t+1))) >= 0 ,
  $
  vanishing iff $q_j^((t)) = q_j^((t+1))$ almost everywhere.

  We now verify conjugacy of $tilde(p)_j$ for each block. For stick $v_k$ ($k < K$), the log-joint terms involving $v_k$ are
  $
    &(alpha_0 - 1) log(1 - v_k) \
    &quad + sum_(i=1)^n [z_(i k) log v_k + (sum_(j=k+1)^K z_(i j)) log(1 - v_k)] .
  $
  Taking expectation under $q(bold(Z))$ yields $N_k log v_k + (alpha_0 - 1 + sum_(j=k+1)^K N_j) log(1 - v_k)$, which is the exact log-density of $"Beta"(1 + N_k, alpha_0 + sum_(j=k+1)^K N_j)$. For $(bold(mu)_k, bold(Lambda)_k)$, collecting all prior and likelihood quadratic terms gives
  $
    &1 / 2 (nu_0 + N_k - d) log det(bold(Lambda)_k) - 1 / 2 "Tr"(bold(W)_0^(-1) bold(Lambda)_k) \
    &- 1 / 2 beta_0 (bold(mu)_k - bold(m)_0)^top bold(Lambda)_k (bold(mu)_k - bold(m)_0) \
    &- 1 / 2 sum_(i=1)^n r_(i k) (bold(x)_i - bold(mu)_k)^top bold(Lambda)_k (bold(x)_i - bold(mu)_k) .
  $
  Completing the square in $bold(mu)_k$ around $bold(m)_k = (beta_0 bold(m)_0 + N_k bar(bold(x))_k) / (beta_0 + N_k)$ splits the quadratic sum identically into
  $
    &beta_k (bold(mu)_k - bold(m)_k)^top bold(Lambda)_k (bold(mu)_k - bold(m)_k) \
    &quad + "Tr"((N_k bold(S)_k + (beta_0 N_k) / beta_k bold(Delta)_k bold(Delta)_k^top) bold(Lambda)_k) .
  $
  Matching sufficient statistics immediately recovers $cal(N)(bold(m)_k, (beta_k bold(Lambda)_k)^(-1)) cal(W)(bold(W)_k, nu_k)$ with the exact M-step recursions. Finally, evaluating $bb(E)_q [(bold(x)_i - bold(mu)_k)^top bold(Lambda)_k (bold(x)_i - bold(mu)_k)] = d / beta_k + nu_k (bold(x)_i - bold(m)_k)^top bold(W)_k (bold(x)_i - bold(m)_k)$ yields the E-step responsibilities $r_(i k)$.
]

#subsec("B", "Singularity Immunity and Uniform Posterior Regularization")

#theorem(title: "Global Elimination of Covariance Singularities via Loewner Bounds", origin: "Original Proof")[
  For any finite dataset $bold(X) subset.eq bb(R)^d$, any truncation $K >= 1$, and any Normal-Wishart prior with $bold(W)_0 in bb(S)_(succ 0)^d$ and $nu_0 > d - 1$:
  1. Every variational precision scale satisfies the uniform Loewner inequalities
     $
       bold(W)_k^(-1) &succ.eq bold(W)_0^(-1) , \
       bold(0) prec bb(E)_q [bold(Lambda)_k] &= nu_k bold(W)_k prec.eq (nu_0 + n) bold(W)_0 .
     $
  2. Both the marginal log-evidence $log p(bold(X))$ and the variational ELBO $cal(L)(q)$ are globally bounded above by a finite constant independent of $q$, completely eliminating the unbounded singularities of ML-GMM (Alg. 81).
]

#proof[
  *Part 1.* In the variational M-step, both $N_k bold(S)_k = sum_(i=1)^n r_(i k) (bold(x)_i - bar(bold(x))_k)(bold(x)_i - bar(bold(x))_k)^top$ and $((beta_0 N_k) / beta_k) bold(Delta)_k bold(Delta)_k^top$ are positive semidefinite matrices for every $r_(i k) >= 0$ and $beta_0 > 0$. Therefore,
  $
    bold(W)_k^(-1) &= bold(W)_0^(-1) + N_k bold(S)_k + (beta_0 N_k) / beta_k bold(Delta)_k bold(Delta)_k^top \
    &succ.eq bold(W)_0^(-1) succ bold(0) .
  $
  Because matrix inversion is operator-antitone on $bb(S)_(succ 0)^d$, it follows that $bold(0) prec bold(W)_k prec.eq bold(W)_0$. Since $0 <= N_k <= n$, we have $nu_k = nu_0 + N_k <= nu_0 + n$, so the expected precision matrix satisfies $bb(E)_q [bold(Lambda)_k] = nu_k bold(W)_k prec.eq (nu_0 + n) bold(W)_0$ uniformly across all iterations, initializations, and responsibility configurations—even if a single point occupies component $k$ ($N_k = 1$, $bold(S)_k = bold(0)$).

  *Part 2.* For any $(bold(mu), bold(Lambda))$, integrating the likelihood against the proper Normal-Wishart prior $G_0$ yields a proper marginal density $p(bold(X)) = integral p(bold(X) | bold(Theta)) p(bold(Theta)) dif bold(Theta)$. Specifically, conditioning on any assignment vector $bold(Z) in \{1, dots.h, K\}^n$, each cluster's marginal likelihood is a multivariate non-central Student-$t$ integral whose density is bounded above by the prior-predictive maximum at zero scatter:
  $
    p(bold(X) | bold(Z)) <= product_(k=1)^K &(pi^(-N_k d / 2) (beta_0 / (beta_0 + N_k))^(d / 2) \
    &dot (Gamma_d ((nu_0 + N_k) / 2)) / (Gamma_d (nu_0 / 2)) det(bold(W)_0)^(N_k / 2)) ,
  $
  where $Gamma_d$ is the $d$-variate gamma function and we used $det(bold(W)_k^(-1)) >= det(bold(W)_0^(-1))$. Because $p(bold(X)) = sum_(bold(Z)) p(bold(Z)) p(bold(X) | bold(Z))$ is a convex combination of finitely bounded terms, $log p(bold(X)) <= C(n, d, nu_0, beta_0, bold(W)_0) < +infinity$. Finally, since $D_"KL" >= 0$, Definition 1 gives $cal(L)(q) <= log p(bold(X)) < +infinity$ for all $q$.
]

#subsec("C", "Automatic Occam Pruning and Logarithmic Cluster Growth")

#proposition(title: "Component Deactivation and Exact CRP Cluster-Count Expectation", origin: "Original Derivation")[
  1. *(Automatic pruning):* Under CAVI, if a component $k$ has vanishing responsibility mass $N_k = 0$, its variational parameter posterior collapses identically to the prior $q(bold(mu)_k, bold(Lambda)_k) = G_0$ with zero parameter KL penalty $D_"KL" (q(bold(mu)_k, bold(Lambda)_k) || G_0) = 0$, preventing unused components from distorting $cal(L)(q)$.
  2. *(Logarithmic growth):* Under the exact Chinese Restaurant Process (CRP) induced by marginalizing $G tilde.op "DP"(alpha_0, G_0)$, the number of occupied clusters $K_n$ among $n$ observations satisfies the exact identity
     $
       bb(E)[K_n] &= sum_(i=1)^n alpha_0 / (alpha_0 + i - 1) \
       &= alpha_0 [psi(alpha_0 + n) - psi(alpha_0)] ,
     $
     and exhibits the tight logarithmic asymptotics $bb(E)[K_n] = alpha_0 log(1 + n / alpha_0) + cal(O)(1)$ as $n -> infinity$.
]

#proof[
  *Part 1.* Setting $N_k = 0$ in the M-step immediately yields $beta_k = beta_0$, $bold(m)_k = bold(m)_0$, $nu_k = nu_0$, and $bold(W)_k = bold(W)_0$, hence $q(bold(mu)_k, bold(Lambda)_k) = G_0$ and $D_"KL" (q(bold(mu)_k, bold(Lambda)_k) || G_0) = 0$. By contrast, splitting a homogeneous cluster of size $N$ across two occupied components incurs a double Normal-Wishart KL complexity cost of order $(d (d + 3) / 2) log N$ plus stick-breaking entropy penalty, creating an intrinsic Bayesian Occam's razor that drives redundant responsibilities $r_(i k) -> 0$.

  *Part 2.* Under the Blackwell-MacQueen Pólya urn / CRP scheme obtained by integrating out $G tilde.op "DP"(alpha_0, G_0)$ @rasmussen2000igmm @teh2006hdp, observation $i$ opens a new cluster with conditional probability $p_i = alpha_0 / (alpha_0 + i - 1)$, independently of the existing partition sizes. Writing $K_n = sum_(i=1)^n I_i$ where $I_i = bb(1)["obs. " i " starts new cluster"]$, linearity of expectation gives
  $
    bb(E)[K_n] = sum_(i=1)^n alpha_0 / (alpha_0 + i - 1) = alpha_0 [psi(alpha_0 + n) - psi(alpha_0)] ,
  $
  using the recurrence $psi(x + 1) - psi(x) = 1 / x$. Comparing the monotone decreasing sum with the integral $integral_0^n (alpha_0 / (alpha_0 + t)) dif t = alpha_0 log(1 + n / alpha_0)$ yields the exact two-sided bound
  $
    alpha_0 log(1 + n / alpha_0) &<= bb(E)[K_n] \
    &<= 1 + alpha_0 log(1 + (n - 1) / alpha_0) ,
  $
  confirming sublinear logarithmic growth in $n$ controlled directly by $alpha_0$.
]

#subsec("D", "Stick-Breaking Truncation Error Bound")

#theorem(title: "Exponential Decay of Truncation Approximation Error", origin: "Classical — Cited")[
  Let $p_infinity (bold(X))$ denote the exact marginal likelihood of the infinite DPGMM and $p_K (bold(X))$ denote the marginal likelihood under truncation at level $K$ ($v_K = 1$). Then the expected residual tail weight satisfies $bb(E)[sum_(k=K+1)^infinity pi_k (bold(v))] = (alpha_0 / (1 + alpha_0))^K$, and the total variation distance between the two marginal distributions obeys the exponential bound @blei2006variational
  $
    ||p_K - p_infinity||_"TV" <= 4 n exp(-(K - 1) / alpha_0) .
  $
]

#proof[
  For independent sticks $v_j tilde.op "Beta"(1, alpha_0)$, each factor satisfies $bb(E)[1 - v_j] = alpha_0 / (1 + alpha_0)$. The remaining stick length after $K$ breaks is $R_K := sum_(k=K+1)^infinity pi_k (bold(v)) = product_(j=1)^K (1 - v_j)$, so independence yields
  $
    bb(E)[R_K] = product_(j=1)^K bb(E)[1 - v_j] = (alpha_0 / (1 + alpha_0))^K .
  $
  Now couple $p_K$ and $p_infinity$ on a common probability space sharing $(v_1, dots.h, v_(K-1))$ and component atoms $\{(bold(mu)_k, bold(Lambda)_k)\}$, differing only when at least one of the $n$ observations selects an index $z_i >= K$ in the infinite model. By the union bound over $i = 1, dots.h, n$, the probability of coupling disagreement is bounded by $n bb(E)[R_(K-1)]$, and converting $L_1$ / total variation bounds via the Ishwaran-James inequality cited in @blei2006variational gives
  $
    ||p_K - p_infinity||_"TV" &<= 4 (1 - bb(E)[(1 - R_(K-1))^n]) \
    &<= 4 n bb(E)[R_(K-1)] .
  $
  Using $(alpha_0 / (1 + alpha_0))^(K-1) = (1 + 1 / alpha_0)^(-(K-1))$ and the elementary bound $1 / (1 + u) <= exp(-u / (1 + u))$ (and the tighter continuum tail integral in @blei2006variational), the total variation discrepancy decays exponentially in $K / alpha_0$, so choosing $K = ceil(alpha_0 log(4 n / delta)) + 1$ renders truncation error negligible ($<= delta$).
]

#subsec("E", "Heavy-Tailed Predictive Distribution")

#proposition(title: "Multivariate Student-t Posterior Predictive Density", origin: "Original Derivation")[
  Under the converged variational posterior $q(bold(mu)_k, bold(Lambda)_k) = cal(N)(bold(m)_k, (beta_k bold(Lambda)_k)^(-1)) cal(W)(bold(W)_k, nu_k)$, the posterior predictive density for a new observation $bold(x)_* in bb(R)^d$ given component $k$ is the exact multivariate Student-$t$ distribution
  $
    p(bold(x)_* | z_* = k, bold(X)) approx "St"(bold(x)_* | bold(m)_k, bold(L)_k, nu_k') ,
  $
  with degrees of freedom $nu_k' = nu_k + 1 - d$ and precision scale matrix $bold(L)_k = ((nu_k' beta_k) / (beta_k + 1)) bold(W)_k$.
]

#proof[
  Integrating first over $bold(mu)_k | bold(Lambda)_k tilde.op cal(N)(bold(m)_k, (beta_k bold(Lambda)_k)^(-1))$ convolves two Gaussians, yielding $bold(x)_* | bold(Lambda)_k tilde.op cal(N)(bold(m)_k, ((beta_k + 1) / beta_k) bold(Lambda)_k^(-1))$. Writing $c_k = beta_k / (beta_k + 1)$ and integrating $bold(Lambda)_k$ against the Wishart density $cal(W)(bold(W)_k, nu_k)$ produces
  $
    &integral_(bb(S)_(succ 0)^d) (c_k / (2 pi))^(d / 2) det(bold(Lambda)_k)^((nu_k - d) / 2) \
    &quad dot exp(-1 / 2 "Tr"((bold(W)_k^(-1) + c_k bold(u) bold(u)^top) bold(Lambda)_k)) dif bold(Lambda)_k ,
  $
  where $bold(u) = bold(x)_* - bold(m)_k$. Recognizing the normalizing constant of a Wishart density with updated degrees of freedom $nu_k + 1$ and inverse scale $bold(W)_k^(-1) + c_k bold(u) bold(u)^top$, and applying the matrix determinant lemma $det(bold(W)_k^(-1) + c_k bold(u) bold(u)^top) = det(bold(W)_k^(-1)) (1 + c_k bold(u)^top bold(W)_k bold(u))$, yields
  $
    p(bold(x)_* | z_* = k, bold(X)) prop (1 + (bold(u)^top bold(L)_k bold(u)) / nu_k')^(-(nu_k' + d) / 2) ,
  $
  which is precisely the multivariate Student-$t$ density with $nu_k' = nu_k + 1 - d$ degrees of freedom. For sparse clusters ($N_k$ small), $nu_k'$ is small and tails are polynomially heavy, naturally expressing epistemic uncertainty; as $N_k -> infinity$, $nu_k' -> infinity$ and the predictive density converges to $cal(N)(bold(m)_k, (nu_k bold(W)_k)^(-1))$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

*Computational Cost and Overhead over ML-GMM.* Asymptotically, truncated CAVI for DPGMM matches full-covariance ML-GMM at $cal(O)(T n K d^2 + T K d^3)$ time and $cal(O)(n K + K d^2)$ memory, where $K$ is the truncation cap rather than the final cluster count $K_"eff"$. However, its wall-clock constant factor is typically $3$--$6 times$ larger than ML-GMM for three reasons: (i) $K$ must exceed the anticipated $K_"eff"$ (Theorem 4); (ii) every iteration evaluates $K (d + 2)$ digamma functions $psi(x)$ and explicit inversions or Cholesky solves for $bold(W)_k$; and (iii) the stick-breaking ordering induces sequential coupling across $gamma_(k,2) = alpha_0 + sum_(j > k) N_j$, requiring more iterations to settle component weights. In return, DPGMM eliminates the external model-selection loop over $K in \{1, dots.h, K_"max"\}$ required by BIC-based GMM @fraley2002mclust.

*Geometric and Statistical Profile.* Rather than sharp Gaussian likelihood contours, DPGMM induces posterior-predictive Student-$t$ contours (Proposition 2) whose effective covariance is inflated by $(beta_k + 1) / beta_k$ and shrunk toward $bold(W)_0^(-1)$ whenever $N_k$ is small. This automatic Stein-type covariance shrinkage prevents ill-conditioned needle ellipsoids in high dimensions and keeps decision boundaries smooth even near sparsely populated modes.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Objective], [Maximize stick-breaking ELBO $cal(L)(q) = log p(bold(X)) - D_"KL" (q || p)$],
    [Inference Engine], [Truncated CAVI @blei2006variational (or collapsed Gibbs @rasmussen2000igmm)],
    [Time Complexity], [$cal(O)(T n K d^2 + T K d^3)$ (single run, no $K$-grid search)],
    [Space Complexity], [$cal(O)(n K + K d^2)$ for responsibilities and Wishart scales],
    [Automatic $K$], [Yes; $K_"eff" <= K$ selected via Bayesian Occam pruning],
    [Singularity Safe], [Strictly immune ($bold(W)_k^(-1) succ.eq bold(W)_0^(-1)$, Theorem 2)],
    [Predictive Shape], [Anisotropic multivariate Student-$t$ (converges to Gaussian)],
    [Key Priors], [Concentration $alpha_0$, truncation $K$, NIW $(bold(m)_0, beta_0, bold(W)_0, nu_0)$],
  ),
  caption: [Operational and theoretical profile of Bayesian GMM / DPGMM (Alg. 82).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Single-Cell Mass Cytometry (CyTOF) and Rare Cell Discovery:* Immunological samples contain unknown numbers of rare activated T-cell and stem-cell subpopulations alongside abundant canonical lineages. DPGMM discovers $K_"eff"$ automatically while Wishart shrinkage prevents rare clusters ($N_k approx 20$ cells in $d = 30$ markers) from collapsing to singular covariances.
2. *Zero-Resource Acoustic Unit Discovery:* In undocumented languages, the inventory of phonetic sub-word units is unknown beforehand. DPGMM applied to spectral frame embeddings automatically infers an adaptive phoneme inventory whose size grows sublinearly with corpus size (Proposition 1) without fixed codebook tuning.
3. *Galactic Archeology and Stellar Stream Decomposition:* Identifying accreted dwarf-galaxy remnants within Milky Way phase-space catalogs requires detecting an unknown number of overlapping anisotropic velocity clumps. DPGMM yields calibrated posterior inclusion probabilities and predictive Student-$t$ tails that resist stray halo interlopers.
4. *Cyber-Physical Telemetry and Novel Fault Mode Discovery:* Industrial turbines and network monitors accumulate normal operating regimes while novel fault states emerge unpredictably. Because the Dirichlet Process allocates nonzero prior mass to unseen components, DPGMM accommodates emerging operational regimes without retraining a fixed-$K$ model family.
5. *Multi-Parametric MRI Lesion Segmentation:* Healthy brain tissues (gray matter, white matter, CSF) have known modes, whereas tumors, edema, and necrosis introduce an unknown number of pathological intensity clusters per patient. DPGMM adapts $K_"eff"$ per scan while preserving smooth posterior uncertainty maps along tumor margins.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Four practical failure regimes require careful diagnostics: (i) *prior-scale sensitivity*—an overly diffuse $bold(W)_0^(-1)$ or small $alpha_0$ can aggressively merge moderately separated clusters on small $n$ (over-pruning), whereas overly tight $bold(W)_0^(-1)$ causes spurious over-splitting; (ii) *mean-field local optima and posterior variance underestimation*—factorizing $q(bold(Z)) q(bold(v)) q(bold(mu), bold(Lambda))$ ignores posterior correlations between stick weights and component locations, causing CAVI to converge to suboptimal modes unless paired with split-merge moves or collapsed Gibbs sampling @rasmussen2000igmm; (iii) *truncation saturation*—if all $K$ sticks remain occupied ($K_"eff" = K$), the truncation bound of Theorem 4 is violated and $K$ must be increased; and (iv) *non-convex manifold geometry*—as with GMM, curved manifolds violate the local elliptical base measure.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Unknown-$K$ continuous ellipsoids], [Automatic $K_"eff"$ and singularity-free Wishart shrinkage], [Check $K_"eff" < K$; tune $alpha_0$ and $bold(W)_0$],
    [Small-$n$ closely spaced clusters], [Conservative Occam prior may merge adjacent modes], [Use ML-GMM (Alg. 81) if $K$ is known],
    [Moderate-high $d$ ($d approx 20$--$100$)], [Stronger than ML-GMM via $bold(W)_0^(-1)$ regularization], [Use diagonal prior or LDA/pLSA on counts],
    [Multi-group shared clusters], [Single DP cannot share atoms across groups], [Hierarchical DP (HDP) @teh2006hdp],
    [Curved / density-connected shapes], [Approximates manifolds as disjoint ellipsoids], [HDBSCAN @campello2013hdbscan or Spectral @ng2002spectral],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for DPGMM.],
)

#sec("V", "Structural Relations and Empirical Behavior")

#proposition(title: "Large-Sample BIC Equivalence of the DPGMM Variational Penalty", origin: "Original Derivation")[
  Let $K_"eff"$ be the number of active full-covariance components in a regular non-degenerate mixture with $N_k = Theta(n)$ for every active $k <= K_"eff"$, and let $p_"eff" = K_"eff" - 1 + K_"eff" d + K_"eff" d (d + 1) / 2$ be the active parameter dimension. As $n -> infinity$, the converged DPGMM variational ELBO admits the asymptotic Laplace expansion
  $
    cal(L)(q^*) = ell(hat(bold(theta))_"MLE") - (p_"eff" / 2) log n + cal(O)(1) ,
  $
  establishing that DPGMM's automatic component pruning asymptotically implements the exact Schwarz BIC complexity penalty @schwarz1978bic inside a single continuous optimization.
]

#proof[
  For each active component $k <= K_"eff"$, $N_k = Theta(n)$, so the posterior Normal-Wishart parameters scale as $beta_k = Theta(n)$, $nu_k = Theta(n)$, and $bold(W)_k^(-1) = N_k bold(S)_k + cal(O)(1)$. Evaluating the closed-form KL divergence between $q(bold(mu)_k, bold(Lambda)_k)$ and the fixed prior $G_0$: (i) the conditional Gaussian factor $cal(N)(bold(m)_k, (beta_k bold(Lambda)_k)^(-1))$ over the $d$ mean coordinates has posterior precision scaling as $beta_k = Theta(n)$, contributing $(d / 2) log(beta_k / beta_0) + cal(O)(1) = (d / 2) log n + cal(O)(1)$ to the KL penalty; and (ii) the Wishart factor $cal(W)(bold(W)_k, nu_k)$ concentrates over the $d (d + 1) / 2$ independent symmetric entries of $bold(Lambda)_k$ with posterior curvature scaling as $nu_k = Theta(n)$ per coordinate, so by Stirling's expansion of $log Gamma_d (nu_k / 2)$ its negative differential entropy contributes $(d (d + 1) / 4) log n + cal(O)(1) = (1 / 2) [d (d + 1) / 2] log n + cal(O)(1)$. Together, each active component $k$ contributes $-(p_k / 2) log n + cal(O)(1)$ to the ELBO, where $p_k = d + d (d + 1) / 2$. Similarly, the $K_"eff" - 1$ active Beta stick posteriors have concentration sums $gamma_(k,1) + gamma_(k,2) = Theta(n)$, each contributing $-(1 / 2) log n + cal(O)(1)$. Meanwhile, inactive components ($k > K_"eff"$) contribute $cal(O)(1)$ by Proposition 1. Summing across all active blocks yields $cal(L)(q^*) = ell(hat(bold(theta))_"MLE") - (p_"eff" / 2) log n + cal(O)(1)$.
]

Within Category IV, DPGMM bridges parametric mixtures and nonparametric Bayes: removing the priors recovers ML-GMM (Alg. 81 @dempster1977em); replacing stick-breaking with a BIC outer search yields MCLUST @fraley2002mclust, AutoClass @cheeseman1996autoclass, or Minimum Message Length Snob @wallace1968snob; extending a single DP across grouped corpora yields Hierarchical Dirichlet Processes (HDP) @teh2006hdp; and isotropic automatic-$K$ counterparts include X-means @pelleg2000xmeans and G-means @hamerly2003gmeans.

*Empirical Benchmark Behavior (ClusterOverview 4-Dataset Suite).* Across the four standardized benchmarks, DPGMM records: *Iris* Silhouette $= 0.302$, ARI $= 0.560$; *Customers* Silhouette $= 0.353$; *Multi-Shape* Silhouette $= 0.782$ (ARI N/A); *High-D* Silhouette $= 0.208$ (ARI N/A); and mean runtime $5.162$ ms. Comparing directly against ML-GMM (Alg. 81) reveals a textbook demonstration of the bias-variance trade-off proved in Section II. On *Iris* ($n = 150$), DPGMM's ARI drops to $0.560$ (versus $0.950$ for ML-GMM): the Occam's razor penalty of Proposition 3 merges the closely overlapping *versicolor* and *virginica* ellipsoids into a single broader mode ($K_"eff" = 2$), demonstrating conservative under-splitting when $n$ is small relative to $d (d + 3) / 2$. Conversely, on *every other dataset*, DPGMM strictly outperforms ML-GMM: *Customers* rises from $0.332$ to $0.353$, *Multi-Shape* surges from $0.712$ to $0.782$ (the highest Multi-Shape Silhouette among Algs. 81--84), and *High-D* improves from $0.181$ to $0.208$ thanks to Loewner precision shrinkage (Theorem 2). Its $5.162$ ms mean runtime is $3.7 times$ ML-GMM ($1.400$ ms), accurately matching the digamma and Wishart inversion overhead analyzed in Section III while remaining $4.7 times$ faster than LDA (Alg. 83).
