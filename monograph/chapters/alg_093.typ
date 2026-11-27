#import "../template/ieee.typ": *

#alg_entry(
  id: 93,
  name: "MCLUST",
  category: "Cat. IV: Distribution, Statistical & Topic",
  year: "2002",
  paradigm: "Parsimonious Eigenvalue-Decomposed GMM & BIC Selection",
  complexity: $cal(O)(M T K (n d^2 + d^3))$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Finite Gaussian mixture models (GMMs) provide a principled probabilistic framework for continuous clustering, yet unrestricted covariance estimation requires $K d (d + 1) / 2$ variance parameters, triggering severe overfitting and likelihood singularities whenever $n / K$ is comparable to $d^2$. *MCLUST* @fraley2002mclust resolves this bias-variance trade-off by parameterizing component covariance matrices through an eigenvalue decomposition that decouples geometric volume, shape, and orientation across a hierarchy of 14 parsimonious models, initialized deterministically by model-based hierarchical agglomerative clustering (MBHAC) and selected automatically via the Bayesian Information Criterion (BIC) @schwarz1978bic.

Let $X = {bold(x)_1, dots, bold(x)_n} subset bb(R)^d$ denote an i.i.d. sample drawn from a $K$-component Gaussian mixture density
$ p(bold(x) | Theta) = sum_(k=1)^K pi_k cal(N)(bold(x) | bold(mu)_k, bold(Sigma)_k), $
where $bold(pi) = (pi_1, dots, pi_K) in Delta_(K-1)$ satisfies $pi_k > 0$, $sum_(k=1)^K pi_k = 1$, $bold(mu)_k in bb(R)^d$, and $bold(Sigma)_k in bb(S)_(>0)^d$ is symmetric positive definite. The observed-data log-likelihood is $ell(Theta) = sum_(i=1)^n log p(bold(x)_i | Theta)$, which MCLUST maximizes over a structured family $cal(M)$ of covariance constraints via the Expectation-Maximization (EM) algorithm @dempster1977em.

#definition(title: "Banfield-Raftery-Celeux-Govaert Spectral Covariance Family")[
  Each component covariance matrix $bold(Sigma)_k in bb(S)_(>0)^d$ admits the unique spectral factorization (up to eigenvector sign flips and eigenvalue permutations)
  $ bold(Sigma)_k = lambda_k bold(D)_k bold(A)_k bold(D)_k^top, $
  where $lambda_k = |bold(Sigma)_k|^(1 / d) > 0$ is a scalar controlling cluster *volume*, $bold(D)_k in O(d)$ is an orthogonal matrix of eigenvectors controlling *orientation*, and $bold(A)_k = "diag"(a_(k 1), dots, a_(k d))$ with $a_(k j) > 0$ and $product_(j=1)^d a_(k j) = 1$ controls *shape*. Constraining each factor to be equal across clusters ($"E"$), variable across clusters ($"V"$), or identity ($"I"$) generates 14 parsimonious models: spherical ($"EII"$, $"VII"$), diagonal ($"EEI"$, $"VEI"$, $"EVI"$, $"VVI"$), and ellipsoidal ($"EEE"$, $"EVE"$, $"VEE"$, $"VVE"$, $"EEV"$, $"VEV"$, $"EVV"$, $"VVV"$).
]

For each candidate pair $(cal(M), K)$ in the search grid $cal(G) times {1, dots, K_"max"}$, MCLUST evaluates the penalized log-likelihood under the positive-orientation BIC convention @fraley2002mclust:
$ "BIC"(cal(M), K) = 2 ell(hat(Theta)_(cal(M), K)) - nu_(cal(M), K) log n, $
where $hat(Theta)_(cal(M), K)$ is the constrained maximum-likelihood estimate (MLE) and $nu_(cal(M), K)$ is the exact effective number of free real parameters in model $(cal(M), K)$. Maximizing $"BIC"(cal(M), K)$ is algebraically equivalent to minimizing Schwarz's original penalty $-2 ell + nu log n$ @schwarz1978bic and approximates twice the log marginal likelihood $2 log p(X | cal(M), K)$.

#algbox(
  title: "MCLUST (Model-Based Hierarchical Initialization + Constrained EM + BIC)",
  inputs: [Data matrix $X in bb(R)^(n times d)$, candidate model set $cal(G)$, maximum clusters $K_"max"$, convergence tolerance $epsilon_"tol" > 0$.],
  outputs: [Selected pair $(cal(M)^*, K^*)$, parameter estimate $hat(Theta)^*$, soft responsibilities $tau^*$, MAP partition $hat(bold(c))$.],
)[
  1. *MBHAC Initialization*: Run agglomerative merging maximizing the unconstrained or spherical ($"EII"$) Gaussian classification log-likelihood to build a nested partition tree for $K = 1, dots, K_"max"$.
  2. *Grid Sweep*: For each $K in {1, dots, K_"max"}$ and covariance model $cal(M) in cal(G)$:
     - Initialize responsibilities $tau_(i k)^((0)) = bb(1){c_i^("MBHAC") = k}$ from the $K$-cluster tree cut.
     - *E-step*: For $i = 1, dots, n$ and $k = 1, dots, K$, compute
       $ tau_(i k)^((t)) = (pi_k^((t)) cal(N)(bold(x)_i | bold(mu)_k^((t)), bold(Sigma)_k^((t)))) / (sum_(m=1)^K pi_m^((t)) cal(N)(bold(x)_i | bold(mu)_m^((t)), bold(Sigma)_m^((t)))). $
     - *Constrained M-step*: Update $n_k = sum_(i=1)^n tau_(i k)^((t))$, $pi_k^((t+1)) = n_k / n$, $bold(mu)_k^((t+1)) = (1 / n_k) sum_(i=1)^n tau_(i k)^((t)) bold(x)_i$, scatter matrices
       $ bold(W)_k = sum_(i=1)^n tau_(i k)^((t)) (bold(x)_i - bold(mu)_k^((t+1)))(bold(x)_i - bold(mu)_k^((t+1)))^top, $
       and solve for ${bold(Sigma)_k^((t+1))}_(k=1)^K in cal(M)$ via closed-form spectral updates or alternating MM / orthogonal iterations.
     - Repeat E-step and M-step until $|ell(Theta^((t+1))) - ell(Theta^((t)))| < epsilon_"tol" (1 + |ell(Theta^((t)))|)$.
     - Record $"BIC"(cal(M), K) = 2 ell(hat(Theta)_(cal(M), K)) - nu_(cal(M), K) log n$.
  3. *Model Selection & MAP Assignment*: Choose $(cal(M)^*, K^*) = arg max_(cal(M), K) "BIC"(cal(M), K)$ and assign $hat(c)_i = arg max_k tau_(i k)^*$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

#subsec("A", "Exact Manifold Degrees of Freedom Across All 14 Families")

Because BIC directly penalizes parameter cardinality by $nu_(cal(M), K) log n$, exact dimensional accounting on the orthogonal group $O(d)$ and unit-determinant diagonal manifold is foundational to MCLUST's complexity control.

#proposition(title: "Closed-Form Parameter Complexity of the 14 Covariance Families", origin: "Original Derivation")[
  Let $nu_0 = (K - 1) + K d$ count mixing weights and means, and let $beta = d (d - 1) / 2 = dim O(d)$. The total free parameter count $nu_(cal(M), K) = nu_0 + nu_Sigma (cal(M), K, d)$ decomposes additively into volume ($nu_lambda$), shape ($nu_A$), and orientation ($nu_D$) terms:
  - *Spherical / Diagonal*: $nu_Sigma ("EII") = 1$, $nu_Sigma ("VII") = K$, $nu_Sigma ("EEI") = d$, $nu_Sigma ("VEI") = K + d - 1$, $nu_Sigma ("EVI") = 1 + K (d - 1)$, $nu_Sigma ("VVI") = K d$.
  - *Ellipsoidal (shared $bold(D)$)*: $nu_Sigma ("EEE") = d + beta$, $nu_Sigma ("VEE") = K + d - 1 + beta$, $nu_Sigma ("EVE") = 1 + K (d - 1) + beta$, $nu_Sigma ("VVE") = K d + beta$.
  - *Ellipsoidal (varying $bold(D)_k$)*: $nu_Sigma ("EEV") = d + K beta$, $nu_Sigma ("VEV") = K + d - 1 + K beta$, $nu_Sigma ("EVV") = 1 + K (d - 1) + K beta$, $nu_Sigma ("VVV") = K (d + beta) = K d (d + 1) / 2$.
]

#proof[
  The simplex $Delta_(K-1)$ contributes $K - 1$ parameters and the $K$ mean vectors contribute $K d$, giving $nu_0 = K d + K - 1$. Positive volumes $lambda_k$ add $1$ ($"E"$) or $K$ ($"V"$). Diagonal shapes $bold(A)_k$ obey the smooth unit-determinant constraint $sum_(j=1)^d log a_(k j) = 0$, adding $0$ ($"I"$), $d - 1$ ($"E"$), or $K (d - 1)$ ($"V"$). The orthogonal group $O(d)$ is a compact Lie group whose tangent space at $bold(I)_d$ is the skew-symmetric Lie algebra $frak("so")(d)$ of dimension $beta = d (d - 1) / 2$, adding $0$ ($"I"$), $beta$ ($"E"$), or $K beta$ ($"V"$). Summing $nu_lambda + nu_A + nu_D$ yields all 14 formulas; $"VVV"$ recovers $K d (d + 1) / 2$, whereas diagonal and spherical models reduce covariance degrees of freedom from $Theta(K d^2)$ to $cal(O)(K d)$ and $cal(O)(K)$.
]

#subsec("B", "Variational Optimality of Constrained M-Steps and EM Monotonicity")

We next prove that MCLUST's core spectral M-step updates globally maximize the EM surrogate objective over the constrained covariance manifolds and induce monotonic log-likelihood ascent.

#theorem(title: "Exact Constrained M-Step Optimality and Monotonic Ascent", origin: "Original Proof")[
  At EM iteration $t$, maximizing the expected complete-data log-likelihood $Q(Theta | Theta^((t)))$ over ${bold(Sigma)_k}_(k=1)^K in cal(M)$ is equivalent to minimizing
  $ J(bold(Sigma)_1, dots, bold(Sigma)_K) = sum_(k=1)^K [ n_k log |bold(Sigma)_k| + "Tr"(bold(W)_k bold(Sigma)_k^(-1)) ]. $
  Whenever $bold(W)_k in bb(S)_(>0)^d$, the unique global minimizers for the closed-form families are:
  1. ($"EII"$, $"VII"$): $hat(lambda) = (1 / (n d)) sum_(k=1)^K "Tr"(bold(W)_k)$ and $hat(lambda)_k = "Tr"(bold(W)_k) / (n_k d)$.
  2. ($"EVV"$): $hat(bold(C))_k = bold(W)_k / |bold(W)_k|^(1 / d)$ and $hat(lambda) = (1 / n) sum_(k=1)^K |bold(W)_k|^(1 / d)$.
  3. ($"EEV"$): With spectral decompositions $bold(W)_k = bold(L)_k bold(Omega)_k bold(L)_k^top$ ordered with non-increasing eigenvalues in $bold(Omega)_k$, $hat(bold(D))_k = bold(L)_k$, $hat(bold(A)) = bold(Omega)_"tot" / |bold(Omega)_"tot"|^(1 / d)$, and $hat(lambda) = |bold(Omega)_"tot"|^(1 / d) / n$, where $bold(Omega)_"tot" = sum_(k=1)^K bold(Omega)_k$.

  Moreover, every valid EM update satisfies $ell(Theta^((t+1))) >= ell(Theta^((t)))$.
]

#proof[
  Expanding the Gaussian log-density in $Q(Theta | Theta^((t)))$ and collecting covariance terms yields
  $ Q_Sigma = -1 / 2 sum_(k=1)^K [ n_k log |bold(Sigma)_k| + "Tr"(bold(W)_k bold(Sigma)_k^(-1)) ] + "const", $
  establishing equivalence with minimizing $J$.

  For $"EII"$ ($bold(Sigma)_k = lambda bold(I)_d$), $J(lambda) = n d log lambda + lambda^(-1) sum_k "Tr"(bold(W)_k)$. Differentiating with respect to $lambda > 0$ gives $n d / lambda - lambda^(-2) sum_k "Tr"(bold(W)_k) = 0$, whose unique root $hat(lambda) = (1 / (n d)) sum_k "Tr"(bold(W)_k)$ is strictly global since $J''(hat(lambda)) = n d / hat(lambda)^2 > 0$. Setting $n_k$ per component yields $"VII"$ identically.

  For $"EVV"$ ($bold(Sigma)_k = lambda bold(C)_k$ with $|bold(C)_k| = 1$), fix $lambda > 0$. Each term $"Tr"(bold(W)_k bold(C)_k^(-1))$ satisfies the AM-GM matrix inequality via eigenvalues $gamma_1, dots, gamma_d$ of $bold(W)_k^(1/2) bold(C)_k^(-1) bold(W)_k^(1/2)$:
  $ "Tr"(bold(W)_k bold(C)_k^(-1)) = sum_(j=1)^d gamma_j >= d (product_(j=1)^d gamma_j)^(1 / d) = d |bold(W)_k|^(1 / d), $
  with equality iff $bold(W)_k^(1/2) bold(C)_k^(-1) bold(W)_k^(1/2) = |bold(W)_k|^(1/d) bold(I)_d$, i.e., $hat(bold(C))_k = bold(W)_k / |bold(W)_k|^(1/d)$. Substituting back gives $J(lambda) = n d log lambda + (d / lambda) sum_k |bold(W)_k|^(1/d)$, uniquely minimized at $hat(lambda) = (1 / n) sum_k |bold(W)_k|^(1/d)$.

  For $"EEV"$ ($bold(Sigma)_k = lambda bold(D)_k bold(A) bold(D)_k^top$), fix $lambda$ and $bold(A) = "diag"(a_1, dots, a_d)$ with $a_1 >= dots >= a_d > 0$ and $|bold(A)| = 1$. By von Neumann's trace inequality, for any orthogonal $bold(D)_k in O(d)$,
  $ "Tr"(bold(W)_k bold(D)_k bold(A)^(-1) bold(D)_k^top) >= sum_(j=1)^d omega_(k j) a_j^(-1), $
  where $omega_(k 1) >= dots >= omega_(k d)$ are the ordered eigenvalues of $bold(W)_k$, with global equality achieved at $hat(bold(D))_k = bold(L)_k$. Summing over $k$ reduces the remaining problem to minimizing $n d log lambda + lambda^(-1) sum_(j=1)^d (sum_k omega_(k j)) a_j^(-1)$ subject to $product_j a_j = 1$. Applying AM-GM as above yields $hat(bold(A)) = bold(Omega)_"tot" / |bold(Omega)_"tot"|^(1/d)$ and $hat(lambda) = |bold(Omega)_"tot"|^(1/d) / n$, confirming global optimality.

  Finally, by Jensen's inequality on the posterior $tau^((t))$,
  $ ell(Theta^((t+1))) - ell(Theta^((t))) &= Q(Theta^((t+1)) | Theta^((t))) - Q(Theta^((t)) | Theta^((t))) \
    &quad + sum_(i=1)^n D_"KL" (tau_i^((t)) || tau_i^((t+1))) >= 0, $
  because $Q$ increases at the M-step and $D_"KL" >= 0$ @dempster1977em @wu1983convergence.
]

#subsec("C", "Singularity Immunity of Equal-Volume Models vs. Variable-Volume Blow-Up")

A notorious pathology of unrestricted GMMs ($"VVV"$) is the unboundedness of the likelihood function @mclachlan2000finite. We constructively separate MCLUST's 14 models into unbounded and globally bounded regimes.

#theorem(title: "Unboundedness of V-Volume Models vs. Boundedness of EII/EEE Families", origin: "Original Counterexample")[
  Let $X = {bold(x)_1, dots, bold(x)_n} subset bb(R)^d$ contain at least $K + 1$ distinct points ($K >= 2$).
  1. *Variable-Volume Blow-Up*: For every variable-volume family $cal(M) in {"VII", "VEI", "VVI", "VEE", "VVE", "VEV", "VVV"}$, $sup_(Theta in cal(M)) ell(Theta) = +infinity$.
  2. *Equal-Volume Boundedness*: For $"EII"$ (and any equal-volume family whose shape matrices satisfy bounded condition number $kappa(bold(A)_k) <= C < infinity$), $sup_(Theta in cal(M)) ell(Theta) < +infinity$. Specifically, for $"EII"$,
     $ sup_(Theta in "EII") ell(Theta) <= -(n d) / 2 log(2 pi e (delta_K^2) / (n d)), $
     where $delta_K = min_(|S| <= K) max_(i in [n]) min_(bold(c) in S) norm(bold(x)_i - bold(c))_2 > 0$.
]

#proof[
  For part (1), fix $pi_k = 1 / K$, set $bold(mu)_1 = bold(x)_1$, choose arbitrary fixed $bold(mu)_2, dots, bold(mu)_K in bb(R)^d$, set $bold(A)_k = bold(I)_d$, $bold(D)_k = bold(I)_d$, and fix $lambda_2 = dots = lambda_K = 1$ while letting $lambda_1 -> 0^+$. Notice this trajectory lies inside every variable-volume family because $"VII"$ is nested in all seven of them. Evaluating the log-likelihood gives:
  $ ell(Theta) &= log p(bold(x)_1 | Theta) + sum_(i=2)^n log p(bold(x)_i | Theta) \
    &>= log(pi_1 (2 pi lambda_1)^(-d / 2)) + sum_(i=2)^n log(pi_2 cal(N)(bold(x)_i | bold(mu)_2, bold(I)_d)) \
    &= -(d / 2) log lambda_1 + C_0 -> +infinity quad "as" lambda_1 -> 0^+, $
  since $C_0$ is a finite constant independent of $lambda_1$.

  For part (2), consider $"EII"$ with arbitrary $(bold(pi), bold(mu)_1, dots, bold(mu)_K, lambda)$. Because $X$ has at least $K + 1$ distinct points, the optimal $K$-center radius satisfies $delta_K > 0$: by definition of $delta_K$, for *every* choice of $K$ centers $bold(mu)_1, dots, bold(mu)_K$, at least one point $bold(x)_(i^*)$ satisfies $min_k norm(bold(x)_(i^*) - bold(mu)_k)_2 >= delta_K$. Using $sum_k pi_k = 1$ and $cal(N)(bold(x)_i | bold(mu)_k, lambda bold(I)_d) <= (2 pi lambda)^(-d/2)$ for all $i != i^*$, while for $i^*$ every component satisfies $cal(N)(bold(x)_(i^*) | bold(mu)_k, lambda bold(I)_d) <= (2 pi lambda)^(-d/2) exp(-delta_K^2 / (2 lambda))$, we obtain
  $ ell(Theta) <= -(n d) / 2 log(2 pi lambda) - delta_K^2 / (2 lambda). $
  Maximizing the right-hand side over $lambda > 0$ yields the unique stationary point $lambda^* = delta_K^2 / (n d) > 0$, substituting which produces $ell(Theta) <= -(n d / 2) log(2 pi e delta_K^2 / (n d)) < +infinity$. Under bounded condition number $kappa(bold(A)_k) <= C$, the smallest eigenvalue of $bold(Sigma)_k = lambda bold(D)_k bold(A)_k bold(D)_k^top$ satisfies $lambda_"min"(bold(Sigma)_k) >= lambda C^(-1)$, replacing $delta_K^2$ by $delta_K^2 / C$ and preserving a finite global supremum.
]

#subsec("D", "Laplace Approximation and Asymptotic Consistency of BIC")

#theorem(title: "BIC Laplace Approximation and Model-Order Selection", origin: "Classical — Cited")[
  Let $X$ be generated from a regular parametric mixture model $(cal(M)_0, K_0)$ whose true parameter $Theta_0$ lies in the interior of the parameter space with positive-definite Fisher information $bold(I)(Theta_0) in bb(S)_(>0)^(nu_0)$, and let $pi(Theta)$ be a smooth prior strictly positive at $Theta_0$ @schwarz1978bic @fraley2002mclust. Then as $n -> infinity$, the log marginal likelihood satisfies
  $ 2 log p(X | cal(M), K) = 2 ell(hat(Theta)_(cal(M), K)) - nu_(cal(M), K) log n + cal(O)_p (1). $
  Consequently, for any mis-specified family $(cal(M), K)$ with $inf_(Theta in cal(M)) D_"KL" (p_0 || p_Theta) > 0$, $"BIC"(cal(M)_0, K_0) - "BIC"(cal(M), K) = Omega_p (n) -> +infinity$, while for nested over-parameterized covariance models at $K = K_0$ with $nu_(cal(M), K_0) > nu_(cal(M)_0, K_0)$, $"BIC"(cal(M)_0, K_0) - "BIC"(cal(M), K_0) = (nu_(cal(M), K_0) - nu_(cal(M)_0, K_0)) log n - cal(O)_p (1) -> +infinity$.
]

#proof[
  Expanding $ell(Theta) = log p(X | Theta)$ to second order about the consistent MLE $hat(Theta)$ where $nabla ell(hat(Theta)) = bold(0)$ yields $ell(Theta) = ell(hat(Theta)) - (1 / 2)(Theta - hat(Theta))^top bold(H)_n (Theta - hat(Theta)) + o_p (1)$, where $bold(H)_n = -nabla^2 ell(hat(Theta)) = n bold(I)(Theta_0) + o_p (n)$ is positive definite on a shrinking $cal(O)(n^(-1/2))$ neighborhood. Integrating the Gaussian approximation over $bb(R)^(nu)$ gives
  $ p(X | cal(M), K) approx p(X | hat(Theta)) pi(hat(Theta)) (2 pi)^(nu / 2) |bold(H)_n|^(-1 / 2). $
  Taking $2 log$ and substituting $|bold(H)_n| = n^nu |bold(I)(Theta_0)| (1 + o_p (1))$ yields $2 ell(hat(Theta)) - nu log n + nu log(2 pi) + 2 log pi(hat(Theta)) - log |bold(I)(Theta_0)| + o_p (1) = "BIC"(cal(M), K) + cal(O)_p (1)$ @schwarz1978bic. Under underfitting, the Kullback-Leibler gap produces $2(ell(hat(Theta)_0) - ell(hat(Theta))) = 2 n D_"KL" + cal(O)_p (sqrt(n))$, dominating the $cal(O)(log n)$ penalty. Under nested covariance overfitting at fixed $K_0$, Wilks' theorem yields $2(ell(hat(Theta)) - ell(hat(Theta)_0)) ->^d chi^2_(Delta nu) = cal(O)_p (1)$, which is overwhelmed by $Delta nu log n -> +infinity$. Note that testing $K_0$ versus $K_0 + 1$ involves boundary non-regularity ($pi_(K_0 + 1) = 0$) where the likelihood ratio follows a supremum of a chi-square process rather than standard $chi^2$ @mclachlan2000finite, yet BIC remains empirically and theoretically consistent under mild separation conditions @fraley2002mclust.
]

#subsec("E", "Transformation-Group Equivariance Hierarchy")

#proposition(title: "Affine, Diagonal-Scale, and Similarity Equivariance of MCLUST", origin: "Original Proof")[
  Under an invertible affine transformation $tilde(bold(x))_i = bold(M) bold(x)_i + bold(b)$ ($bold(M) in "GL"(d)$), let $hat(Theta)$ maximize $ell(Theta; X)$ in family $cal(M)$.
  1. *Full Affine Equivariance*: For $cal(M) in {"EEE", "VEE", "EVV", "VVV"}$, the transformed estimate $tilde(bold(mu))_k = bold(M) hat(bold(mu))_k + bold(b)$ and $tilde(bold(Sigma))_k = bold(M) hat(bold(Sigma))_k bold(M)^top$ remains inside $cal(M)$, shifts the log-likelihood by the exact model-independent constant $ell(tilde(Theta); tilde(X)) = ell(hat(Theta); X) - n log |det bold(M)|$, and preserves both MAP clusters and BIC rankings among these four families.
  2. *Axis-Aligned Scale Equivariance*: For diagonal families $cal(M) in {"EEI", "VEI", "EVI", "VVI"}$, invariance holds for all nonsingular diagonal matrices $bold(M) = "diag"(s_1, dots, s_d)$ and coordinate permutations.
  3. *Similarity Equivariance*: All 14 families are strictly equivariant under similarities $bold(M) = c bold(Q)$ ($c > 0$, $bold(Q) in O(d)$), whereas spherical models ($"EII"$, $"VII"$) lose invariance under anisotropic rescaling.
]

#proof[
  By the change-of-variables Jacobian for Gaussian densities, $cal(N)(bold(M) bold(x) + bold(b) | bold(M) bold(mu)_k + bold(b), bold(M) bold(Sigma)_k bold(M)^top) = |det bold(M)|^(-1) cal(N)(bold(x) | bold(mu)_k, bold(Sigma)_k)$, shifting $ell$ by $-n log |det bold(M)|$ for every parameter tuple. It therefore suffices to check closure of the covariance constraint under congruence $bold(Sigma)_k |-> bold(M) bold(Sigma)_k bold(M)^top$. For $"VVV"$ and $"EEE"$, positive definiteness and equality $bold(Sigma)_1 = dots = bold(Sigma)_K$ are trivially preserved. For $"VEE"$ ($bold(Sigma)_k = lambda_k bold(C)$), $bold(M) bold(Sigma)_k bold(M)^top = (lambda_k |det bold(M)|^(2/d)) (|det bold(M)|^(-2/d) bold(M) bold(C) bold(M)^top)$ preserves common shape-and-orientation. For $"EVV"$ ($|bold(Sigma)_k| = lambda^d$), $|bold(M) bold(Sigma)_k bold(M)^top| = (|det bold(M)|^(2/d) lambda)^d$ preserves equal volume across $k$. For diagonal families with $bold(M) = "diag"(s_1, dots, s_d)$, congruence multiplies entry $j$ by $s_j^2$, preserving diagonality, volume equality ($|bold(M) bold(Sigma)_k bold(M)^top| = (product_j s_j^2) lambda_k^d$), and shape equality ($tilde(a)_(k j) = s_j^2 a_(k j) / (product_r s_r^2)^(1/d)$). Under similarity $bold(M) = c bold(Q)$, $bold(M) (lambda_k bold(D)_k bold(A)_k bold(D)_k^top) bold(M)^top = (c^2 lambda_k) (bold(Q) bold(D)_k) bold(A)_k (bold(Q) bold(D)_k)^top$, preserving all 14 constraints because $bold(Q) bold(D)_k in O(d)$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Evaluating all $M <= 14$ covariance models across $K in {1, dots, K_"max"}$ dominates MCLUST's runtime. For a fixed $(cal(M), K)$, each E-step evaluates $n K$ quadratic Mahalanobis forms in $cal(O)(n K d^2)$ time for ellipsoidal models (or $cal(O)(n K d)$ for spherical and diagonal families), while the M-step accumulates $K$ scatter matrices $bold(W)_k$ in $cal(O)(n K d^2)$ and performs $K$ eigendecompositions or Cholesky factorizations in $cal(O)(K d^3)$. Over $T$ EM iterations and $M$ models, the total EM complexity is $cal(O)(M T K (n d^2 + d^3))$. Notice that the MBHAC initialization phase incurs an upfront $cal(O)(n^2 d)$ to $cal(O)(n^2 d^2)$ cost if executed on all $n$ samples, so large-$n$ MCLUST implementations subsample $n_0 << n$ points for MBHAC before running full-sample EM @fraley2002mclust.

Memory usage is $cal(O)(n K + K d^2)$ to hold the $n times K$ responsibility matrix and $K$ covariance / scatter tensors. Geometrically, each ellipsoidal model induces quadratic Mahalanobis decision boundaries (reducing to linear hyperplanes for equal-covariance $"EEE"$ and Euclidean Voronoi cells for $"EII"$). Because covariance matrices remain compact and BLAS-3 matrix products dominate scatter accumulation, cache efficiency is high for moderate $d$, though sweeping all 14 models incurs a constant-factor overhead relative to fixed-structure mixture estimators.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$cal(O)(M T K (n d^2 + d^3))$ EM grid ($cal(O)(M T n K d)$ for spherical / diagonal); MBHAC init $cal(O)(n_0^2 d^2)$.],
    [Space Complexity], [$cal(O)(n d + n K + K d^2)$ for data, responsibilities, and $K$ covariance matrices.],
    [Decision Geometry], [Quadratic Mahalanobis regions ($"VVV"$, $"VEV"$); linear hyperplanes ($"EEE"$); Voronoi cells ($"EII"$).],
    [Invariance], [Similarity-equivariant (all 14); axis-scale (diagonal); full $"GL"(d)$-affine ($"EEE"$, $"VEE"$, $"EVV"$, $"VVV"$).],
    [Model Selection], [Automatic joint $(cal(M), K)$ choice via BIC without requiring manual $K$ specification.],
    [Regularization], [Structural parsimony via spectral tying; optional conjugate Normal-Inverse-Wishart MAP priors.],
  ),
  caption: [Operational and geometric profile of MCLUST across its 14-model covariance hierarchy.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Flow Cytometry and Single-Cell Gating*: Multi-parameter fluorescence profiles of immune subpopulations form overlapping elliptical clouds with varying volumes and tilted spectral spill-over orientations; $"VVV"$, $"VEV"$, and $"EEE"$ recover anisotropic cell populations while BIC detects rare subtypes automatically @fraley2002mclust.
2. *Remote Sensing and Multispectral Segmentation*: Radiance vectors across spectral bands exhibit strong cross-band correlation from illumination geometry; ellipsoidal families ($"VEE"$, $"EEV"$, $"VVV"$) model band covariance directly while equal-volume models stabilize small regions.
3. *Gene Expression & Proteomic Profiling*: In cohorts with $n approx 100$--$500$ and $d approx 20$--$200$, unconstrained $"VVV"$ suffers singular scatter matrices; BIC automatically selects diagonal ($"EVI"$, $"VVI"$) or spherical ($"VII"$) models with $cal(O)(K d)$ parameters.
4. *Astronomical Catalog & Exoplanet Taxonomy*: Stellar kinematics and metallicity-luminosity populations form multi-scale heteroscedastic groups; full affine equivariance of $"VVV"$ / $"EEE"$ / $"EVV"$ guarantees unit-independent clustering across heterogeneous physical scales.
5. *Industrial Telemetry Fault Identification*: Sensor streams operate across nominal regimes plus dispersed fault states; parsimonious GMMs (optionally augmented with a uniform noise component @fraley2002mclust) yield calibrated Mahalanobis anomaly scores.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

MCLUST degrades in four regimes: (i) *non-convex manifolds* (rings, spirals), where Gaussian ellipsoids over-split single connected structures; (ii) *heavy tails or skewness*, where outliers inflate covariances (mitigated by Student-$t$ @peel2000robust or skew-normal @lin2009skewnormal mixtures); (iii) *ultra-high dimensions ($d >> n$)* where even diagonal models accumulate irrelevant coordinate noise; and (iv) *discrete binary or sparse count data*, where Gaussian variance assumptions distort zero-inflated likelihoods.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Flow cytometry & biometrics], [Anisotropic elliptical subpopulations match spectral covariance hierarchy.], [Heavy tails / debris: add noise component or Student-$t$ @peel2000robust.],
    [Continuous tabular segmentation], [Heteroscedastic feature correlations captured with automatic BIC $K$.], [Mixed categorical features: use Mixmod (Alg. 94) @biernacki2006mixmod.],
    [Moderate-$d$ omics ($d approx 50$)], [Diagonal ($"VVI"$, $"EVI"$) avoids $cal(O)(K d^2)$ singular covariance collapse.], [When $d >> n$: apply PCA or sparse subspace clustering.],
    [Curved / manifold data], [Poor fit: ellipsoidal level sets fragment non-convex shapes.], [Use HDBSCAN, Spectral Clustering, or DBSCAN.],
    [Discrete counts / binary text], [Misspecified Gaussian support on $bb(Z)_(>=0)^d$ or ${0, 1}^d$.], [Use Poisson (Alg. 95) or Bernoulli (Alg. 96) mixtures.],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for MCLUST.],
)

#sec("V", "Structural Relations and Empirical Behavior")

MCLUST occupies the central hub of continuous parametric clustering, bridging hard geometric partitional methods, agglomerative dendrograms, and broader finite-mixture frameworks.

#proposition(title: "Exact Reduction of MCLUST-EII to k-Means and Ward's Criterion", origin: "Original Derivation")[
  Consider MCLUST restricted to spherical equal-volume isotropic models ($"EII"$, $bold(Sigma)_k = lambda bold(I)_d$) with equal weights $pi_k = 1 / K$.
  1. As $lambda -> 0^+$, the E-step responsibilities converge pointwise to the hard Voronoi indicator $tau_(i k) -> bb(1){k = arg min_m norm(bold(x)_i - bold(mu)_m)_2^2}$ (off decision boundaries), reducing EM identically to Lloyd's $k$-means @lloyd1982least.
  2. In MBHAC initialization under the $"EII"$ classification log-likelihood $ell_C = -(n d / 2) log(2 pi e "SSE" / (n d))$, merging clusters $C_a$ and $C_b$ increases total within-cluster sum of squares $"SSE"$ by the exact Ward increment
     $ Delta "SSE"(C_a, C_b) = (|C_a| |C_b|) / (|C_a| + |C_b|) norm(bar(bold(x))_a - bar(bold(x))_b)_2^2, $
     making $"EII"$-MBHAC strictly order-equivalent to Ward's minimum-variance hierarchical clustering @ward1963hierarchical.
]

#proof[
  For part (1), $tau_(i k) = (1 + sum_(m != k) exp(-(norm(bold(x)_i - bold(mu)_m)_2^2 - norm(bold(x)_i - bold(mu)_k)_2^2) / (2 lambda)))^(-1)$. If $k$ uniquely minimizes $norm(bold(x)_i - bold(mu)_m)_2^2$, every exponent tends to $-infinity$ as $lambda -> 0^+$, forcing $tau_(i k) -> 1$ and $tau_(i m) -> 0$ for $m != k$, so the M-step mean update coincides with Lloyd's centroid update. For part (2), maximizing the $"EII"$ classification likelihood over $lambda$ gives $hat(lambda) = "SSE" / (n d)$ and profiled likelihood $ell_C = -(n d / 2) log(2 pi e "SSE" / (n d))$, a strictly decreasing function of $"SSE"$. Expanding $sum_(bold(x) in C_a union C_b) norm(bold(x) - bar(bold(x))_(a b))_2^2$ via Huygens' parallel-axis theorem around $bar(bold(x))_(a b) = (|C_a| bar(bold(x))_a + |C_b| bar(bold(x))_b) / (|C_a| + |C_b|)$ yields $Delta "SSE" = (|C_a| |C_b| / (|C_a| + |C_b|)) norm(bar(bold(x))_a - bar(bold(x))_b)_2^2$, proving exact merge-order equivalence with Ward's algorithm.
]

MCLUST generalizes X-means @pelleg2000xmeans (which restricts BIC search to spherical equal-variance models) and is itself embedded within Mixmod (Alg. 94) @biernacki2006mixmod, which adds equal-proportion Gaussian models, qualitative multinomial families, SEM / CEM inference, and ICL selection. Nonparametric Dirichlet Process Gaussian mixtures @rasmussen2000igmm replace discrete BIC model grids with an infinite Bayesian prior.

On the ClusterOverview 4-dataset empirical benchmark (Silhouette @rousseeuw1987silhouettes and Adjusted Rand Index @hubert1985comparing), MCLUST confirms the precision-cost trade-off predicted by our theory:
- *Iris* ($"Sil" = 0.586$, $"ARI" = 0.567$): At small $n = 150$, BIC's $nu log n$ penalty selects a parsimonious low-order fit cleanly separating _Setosa_ from the overlapping _Versicolor_ / _Virginica_ ellipsoid, achieving the highest Iris Silhouette across Algs. 93--96 with conservative label granularity.
- *Customers* ($"Sil" = 0.385$) & *Multi-Shape* ($"Sil" = 0.782$, $"ARI" = "N/A"$): Strongest Silhouette on both datasets among Algs. 93--96, as anisotropic spectral tying captures multi-scale spending correlations and clean convex modes.
- *High-D* ($"Sil" = 0.321$, $"ARI" = "N/A"$): Highest High-D Silhouette in this quartet, directly validating Proposition 1 as BIC prunes $cal(O)(K d^2)$ orientation parameters in favor of diagonal / spherical models.
- *Runtime* ($68.368$ ms mean): Highest runtime across Algs. 93--96 ($5.6 times$ Mixmod's $12.143$ ms and $19.7 times$ Poisson's $3.468$ ms), reflecting exhaustive multi-model EM and BIC grid evaluation.
