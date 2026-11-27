#import "../template/ieee.typ": *

#alg_entry(
  id: 87,
  name: "AutoClass",
  category: "Cat. IV: Distribution, Statistical & Topic",
  year: "1996",
  paradigm: "Bayesian Mixed-Type Finite Mixture with Marginal-Likelihood Selection",
  complexity: $cal(O)(R T n K (d_c + d_m V))$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

AutoClass @cheeseman1996autoclass casts unsupervised classification as full Bayesian model inference over finite mixtures of mixed-type distributions. Rather than fixing the cluster count $K$ a priori or relying on heuristic online tree utility, AutoClass places conjugate priors over all component and mixing parameters, computes maximum a posteriori (MAP) modes via Bayesian Expectation-Maximization (EM) @dempster1977em, and ranks competing model structures $cal(M)_K$ by their marginal likelihood (Bayesian evidence).

Let $cal(D) = {bold(x)_1, dots, bold(x)_n}$ be $n$ observations, each comprising $d_m$ discrete nominal attributes $x_(i, a) in {1, dots, V_a}$ ($V_a <= V$) and $d_c$ continuous real-valued attributes $x_(i, b) in bb(R)$ (or optional covariant Gaussian blocks $bold(x)_(i, c) in bb(R)^(d_b)$), with arbitrary Missing-at-Random (MAR) entries. Given $K$ latent classes with mixing weights $bold(pi) = (pi_1, dots, pi_K) in Delta^(K - 1)$ and latent indicators $z_i in {1, dots, K}$, the conditionally factorized mixture likelihood is
$ p(bold(x)_i | bold(Theta), cal(M)_K) &= sum_(k=1)^K pi_k product_(a=1)^(d_m) theta_(k, a, x_(i, a)) \
  &quad times product_(b=1)^(d_c) cal(N)(x_(i, b) | mu_(k, b), sigma_(k, b)^2), $
where $bold(Theta) = (bold(pi), {bold(theta)_(k, a)}, {mu_(k, b), sigma_(k, b)^2})$. AutoClass specifies weakly informative conjugate priors:
- $bold(pi) ~ "Dir"(alpha_0 / K, dots, alpha_0 / K)$ with $alpha_0 > 0$;
- $bold(theta)_(k, a) ~ "Dir"(beta_a, dots, beta_a)$ for each nominal attribute $a$;
- $(mu_(k, b), sigma_(k, b)^2) ~ "NIG"(mu_(0, b), kappa_(0, b), a_(0, b), b_(0, b))$ (Normal-Inverse-Gamma), or Normal-Inverse-Wishart for covariant continuous blocks.

#definition(title: "AutoClass Bayesian Evidence and Cheeseman-Stutz Score")[
  The exact marginal likelihood (evidence) of a $K$-class model $cal(M)_K$ is
  $ p(cal(D) | cal(M)_K) = integral p(cal(D) | bold(Theta), cal(M)_K) p(bold(Theta) | cal(M)_K) dif bold(Theta). $
  Because the mixture sum renders this integral analytically intractable, AutoClass evaluates the posterior responsibilities $hat(bold(Z)) = (r_(i k))$ at a MAP-EM stationary point $hat(bold(Theta))_"MAP"$, forms the analytic fractional complete-data marginal $p(cal(D), hat(bold(Z)) | cal(M)_K)$, and scores $cal(M)_K$ by the _Cheeseman-Stutz (CS) log-evidence approximation_ @cheeseman1996autoclass:
  $ log hat(p)_"CS" (cal(D) | cal(M)_K) &= log p(cal(D), hat(bold(Z)) | cal(M)_K) + log (K!) \
    &quad + log p(cal(D) | hat(bold(Theta))_"MAP") - log p(cal(D), hat(bold(Z)) | hat(bold(Theta))_"MAP"). $
]

Here the factor $log(K!)$ accounts for the $K!$ identical label-permutation modes of an exchangeable mixture posterior, while the ratio $p(cal(D) | hat(bold(Theta))_"MAP") / p(cal(D), hat(bold(Z)) | hat(bold(Theta))_"MAP")$ adjusts the complete-data evidence for soft assignment entropy.

#algbox(
  title: "AutoClass Bayesian MAP-EM and Model Evidence Search",
  inputs: [Mixed dataset $cal(D)$ (possibly with MAR missing values), prior hyperparameters, candidate cluster counts $cal(K)$, restarts $R$.],
  outputs: [Optimal model $cal(M)_(K^*)$, MAP parameters $hat(bold(Theta))_"MAP"$, soft posteriors $r_(i k)$, and log-evidence $log hat(p)_"CS"$.],
)[
  + For each candidate $K in cal(K)$ and random restart $r = 1, dots, R$:
  + Initialize responsibilities $r_(i k)^((0))$ (or component parameters $bold(Theta)^((0))$) randomly.
  + *Repeat* for $t = 0, 1, dots, T - 1$ until $|cal(J)(bold(Theta)^((t+1))) - cal(J)(bold(Theta)^((t)))| < epsilon$:
  + *E-step*: for each $i, k$, compute unnormalized log-weight summing only over observed attributes $cal(O)_i$:
    $ log tilde(r)_(i k) &= log pi_k^((t)) + sum_(a in cal(O)_(i, m)) log theta_(k, a, x_(i, a))^((t)) \
      &quad + sum_(b in cal(O)_(i, c)) log cal(N)(x_(i, b) | mu_(k, b)^((t)), (sigma_(k, b)^((t)))^2), $
    and normalize $r_(i k)^((t+1)) = tilde(r)_(i k) / sum_(l=1)^K tilde(r)_(i l)$.
  + *MAP M-step*: compute fractional counts $N_k = sum_(i=1)^n r_(i k)$, observed attribute totals, and update $pi_k$, $theta_(k, a, v)$, $mu_(k, b)$, and $(sigma_(k, b))^2$ in closed form via conjugate MAP equations; prune empty classes ($N_k approx 0$).
  + Compute $log hat(p)_"CS" (cal(D) | cal(M)_K)$ at convergence using the exact closed-form conjugate integrals of Lemma 1.
  + Return the model $(cal(M)_(K^*), hat(bold(Theta))_"MAP")$ maximizing $log hat(p)_"CS" (cal(D) | cal(M)_K)$ across all $K$ and restarts.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five formal results: monotonic ascent and strict singularity prevention of conjugate MAP-EM, exact closed-form mixed-type complete-data marginalization, an exact remainder and asymptotic BIC-equivalence theorem for the Cheeseman-Stutz score, consistency under Missing-at-Random observations, and classical local linear convergence coupled to label-switching multiplicity.

#subsec("A", "Monotonic Ascent and Singularity Prevention of Conjugate MAP-EM")

#theorem(title: "Variational Ascent and Variance Lower Bound of Conjugate MAP-EM", origin: "Original Proof")[
  Let $cal(J)(bold(Theta)) = log p(cal(D) | bold(Theta)) + log p(bold(Theta))$ denote the unnormalized log-posterior under conjugate Dirichlet and Normal-Inverse-Gamma priors with $a_(0, b) > 0, b_(0, b) > 0, kappa_(0, b) > 0$.
  1. Every AutoClass MAP-EM iteration satisfies the strict variational ascent identity
     $ cal(J)(bold(Theta)^((t+1))) - cal(J)(bold(Theta)^((t))) &= [ cal(F)(q^((t)), bold(Theta)^((t+1))) - cal(F)(q^((t)), bold(Theta)^((t))) ] \
       &quad + sum_(i=1)^n D_"KL" (q_i^((t)) || q_i^((t+1))) >= 0, $
     where $q_(i k)^((t)) = p(z_i = k | bold(x)_i, bold(Theta)^((t)))$ and $cal(F)(q, bold(Theta))$ is the penalized variational free energy.
  2. Unlike unregularized Gaussian ML-EM, every continuous component variance at every iteration satisfies the strict uniform floor
     $ (hat(sigma)_(k, b)^((t+1)))^2 >= (2 b_(0, b)) / (n + 2 a_(0, b) + 3) > 0. $
     Consequently, $cal(J)(bold(Theta))$ is bounded above on the parameter domain, preventing singleton likelihood explosions and guaranteeing convergence of ${cal(J)(bold(Theta)^((t)))}$ to a finite limit.
]

#proof[
  Define the penalized free energy for any variational distributions $q = {q_i}_(i=1)^n$ on ${1, dots, K}$ by
  $ cal(F)(q, bold(Theta)) &= sum_(i=1)^n sum_(k=1)^K q_(i k) log (p(bold(x)_i, z_i = k | bold(Theta)) / q_(i k)) \
    &quad + log p(bold(Theta)). $
  Using $p(bold(x)_i, z_i = k | bold(Theta)) = p(bold(x)_i | bold(Theta)) p(z_i = k | bold(x)_i, bold(Theta))$, algebraic rearrangement gives
  $ cal(J)(bold(Theta)) &= cal(F)(q, bold(Theta)) \
    &quad + sum_(i=1)^n D_"KL" (q_i || p(z_i = dot.c | bold(x)_i, bold(Theta))). $
  At the E-step, setting $q_i^((t)) = p(z_i = dot.c | bold(x)_i, bold(Theta)^((t)))$ zeroes the KL divergence, so $cal(J)(bold(Theta)^((t))) = cal(F)(q^((t)), bold(Theta)^((t)))$. At the MAP M-step, $bold(Theta)^((t+1)) = arg max_bold(Theta) cal(F)(q^((t)), bold(Theta))$, so $cal(F)(q^((t)), bold(Theta)^((t+1))) >= cal(F)(q^((t)), bold(Theta)^((t)))$. Evaluating the decomposition at $bold(Theta)^((t+1))$ with $q = q^((t))$ and subtracting $cal(J)(bold(Theta)^((t)))$ yields the exact identity in Claim 1; non-negativity follows from M-step optimality and Gibbs' inequality $D_"KL" >= 0$.

  For Claim 2, fix $(k, b)$ and drop the attribute subscript $b$. Let $N_k = sum_(i=1)^n r_(i k) <= n$, $bar(x)_k = N_k^(-1) sum_i r_(i k) x_i$, and $S_k = sum_i r_(i k) (x_i - bar(x)_k)^2 >= 0$. Up to additive constants, the terms of $cal(F)(q^((t)), bold(Theta))$ depending on $(mu_k, v_k)$ with $v_k = sigma_k^2$ are
  $ Phi(mu_k, v_k) &= -(N_k / 2 + a_0 + 3 / 2) log v_k \
    &quad - 1 / (2 v_k) [ 2 b_0 + S_k + N_k (mu_k - bar(x)_k)^2 + kappa_0 (mu_k - mu_0)^2 ]. $
  Completing the square in $mu_k$ yields the unique maximizer $hat(mu)_k = (kappa_0 mu_0 + N_k bar(x)_k) / (kappa_0 + N_k)$ and quadratic remainder
  $ R_k = 2 b_0 + S_k + (kappa_0 N_k) / (kappa_0 + N_k) (bar(x)_k - mu_0)^2 >= 2 b_0. $
  Setting $partial Phi(hat(mu)_k, v_k) / partial v_k = 0$ yields the unique closed-form variance maximizer
  $ hat(v)_k = R_k / (N_k + 2 a_0 + 3) >= (2 b_0) / (n + 2 a_0 + 3) > 0, $
  since $N_k <= n$ and $S_k >= 0$. This uniform lower bound implies $log cal(N)(x | mu_k, hat(v)_k) <= -(1/2) log(2 pi hat(v)_k)$ is uniformly bounded above, while discrete probabilities satisfy $log theta <= 0$ and proper priors are bounded above on their modes. Hence $cal(J)(bold(Theta)^((t)))$ is bounded and monotone increasing, so it converges to a finite limit.
]

#subsec("B", "Exact Closed-Form Mixed-Type Complete-Data Evidence")

#lemma(title: "Analytic Marginalization Under Conjugate Priors", origin: "Original Derivation")[
  For any fractional assignment matrix $bold(Z) = (r_(i k))$ with cluster masses $N_k = sum_(i=1)^n r_(i k)$, categorical counts $N_(k, a, v) = sum_(i : x_(i, a) = v) r_(i k)$, and continuous sufficient statistics $(N_(k, b), bar(x)_(k, b), S_(k, b))$ over observed entries, the conjugate complete-data marginal likelihood factors into exact closed-form Gamma ratios:
  $ p(cal(D), bold(Z) | cal(M)_K) &= p(bold(Z)) product_(k=1)^K [ product_(a=1)^(d_m) p_(k, a) (cal(D) | bold(Z)) \
    &quad times product_(b=1)^(d_c) p_(k, b) (cal(D) | bold(Z)) ], $
  where
  $ p(bold(Z)) = (Gamma(alpha_0)) / (Gamma(n + alpha_0)) product_(k=1)^K (Gamma(N_k + alpha_0 / K)) / (Gamma(alpha_0 / K)), $
  $ p_(k, a) (cal(D) | bold(Z)) = (Gamma(V_a beta_a)) / (Gamma(N_(k, a) + V_a beta_a)) product_(v=1)^(V_a) (Gamma(N_(k, a, v) + beta_a)) / (Gamma(beta_a)), $
  $ p_(k, b) (cal(D) | bold(Z)) &= (2 pi)^(-N_(k, b) / 2) sqrt(kappa_(0, b) / kappa_(k, b)) \
    &quad times (Gamma(a_(k, b)) b_(0, b)^(a_(0, b))) / (Gamma(a_(0, b)) b_(k, b)^(a_(k, b))), $
  with $kappa_(k, b) = kappa_(0, b) + N_(k, b)$, $a_(k, b) = a_(0, b) + N_(k, b) / 2$, and $b_(k, b) = b_(0, b) + S_(k, b) / 2 + (kappa_(0, b) N_(k, b)) / (2 kappa_(k, b)) (bar(x)_(k, b) - mu_(0, b))^2$.
]

#proof[
  Under conditional independence of attributes given $bold(Z)$ and independent factorized priors across $bold(pi)$, $a$, and $b$, the high-dimensional integral over $bold(Theta)$ splits into independent component integrals. For $bold(pi)$ and each $bold(theta)_(k, a)$, integrating a multinomial likelihood against a symmetric Dirichlet prior directly yields the Dirichlet-Multinomial normalizing constant ratios $p(bold(Z))$ and $p_(k, a) (cal(D) | bold(Z))$.

  For continuous attribute $b$ in cluster $k$ (dropping $b$), multiplying the fractional Gaussian likelihood by the $"NIG"(mu_0, kappa_0, a_0, b_0)$ density and completing the square in $mu_k$ yields
  $ p_(k) (cal(D) | bold(Z)) &= (2 pi)^(-N_k / 2) (b_0^(a_0) / Gamma(a_0)) sqrt(kappa_0 / (2 pi)) integral_0^infinity v^(-(a_k + 3 / 2)) e^(-b_k / v) \
    &quad times ( integral_(-infinity)^(+infinity) exp(-(kappa_k (mu - hat(mu)_k)^2) / (2 v)) dif mu ) dif v. $
  The inner Gaussian integral evaluates to $sqrt(2 pi v / kappa_k)$, cancelling $sqrt(2 pi)$ and reducing the exponent of $v$ to $-(a_k + 1)$. Recognizing the remaining integral as the Inverse-Gamma normalizing constant $integral_0^infinity v^(-(a_k + 1)) e^(-b_k / v) dif v = Gamma(a_k) / b_k^(a_k)$ yields the exact formula for $p_(k, b) (cal(D) | bold(Z))$, computable in $cal(O)(K (d_c + d_m V))$ time via log-Gamma evaluations.
]

#subsec("C", "Exact Remainder and BIC Equivalence of Cheeseman-Stutz")

#theorem(title: "Missing-Information Bias and Asymptotic BIC Equivalence of Cheeseman-Stutz", origin: "Original Proof")[
  Let $bold(H)_"comp" = -nabla^2 log p(bold(Theta) | cal(D), hat(bold(Z)))$ and $bold(H)_"obs" = -nabla^2 log p(bold(Theta) | cal(D))$ be the complete-data and observed-data negative log-posterior Hessians at a regular interior MAP mode $hat(bold(Theta))_"MAP"$ of dimension $p_K$, linked by Louis's missing-information identity $bold(H)_"obs" = bold(H)_"comp" - bold(H)_"miss"$ with $0 prec.eq bold(H)_"miss" prec bold(H)_"comp"$. Then:
  1. Up to the $log(K!)$ symmetry factor for exchangeable modes, the exact log-evidence satisfies
     $ log p(cal(D) | cal(M)_K) - log hat(p)_"CS" (cal(D) | cal(M)_K) \
       = 1 / 2 log det(bold(H)_"comp" bold(H)_"obs"^(-1)) + cal(O)(n^(-1)). $
     In particular, because $det(bold(H)_"comp" bold(H)_"obs"^(-1)) = det(bold(I) - bold(H)_"comp"^(-1) bold(H)_"miss")^(-1) >= 1$, $log hat(p)_"CS"$ asymptotically underestimates the true log-evidence by an $cal(O)(1)$ nonnegative missing-information gap that vanishes as cluster separation grows ($bold(H)_"miss" -> bold(0)$).
  2. As $n -> infinity$ with fixed $K$, the Cheeseman-Stutz score agrees with Schwarz's BIC @schwarz1978bic up to $cal(O)(1)$:
     $ log hat(p)_"CS" (cal(D) | cal(M)_K) &= log p(cal(D) | hat(bold(Theta))_"MAP") \
       &quad - (p_K / 2) log n + cal(O)(1). $
]

#proof[
  Applying Bayes' rule to both the observed marginal $p(cal(D) | cal(M)_K)$ and the complete-data marginal $p(cal(D), hat(bold(Z)) | cal(M)_K)$ at any parameter point $bold(Theta)$ and taking their ratio cancels the prior $p(bold(Theta) | cal(M)_K)$ identically:
  $ (p(cal(D) | cal(M)_K)) / (p(cal(D), hat(bold(Z)) | cal(M)_K)) &= (p(cal(D) | bold(Theta))) / (p(cal(D), hat(bold(Z)) | bold(Theta))) \
    &quad times (p(bold(Theta) | cal(D), hat(bold(Z)))) / (p(bold(Theta) | cal(D))). $
  Evaluate this exact identity at $bold(Theta) = hat(bold(Theta))_"MAP"$. Crucially, because $hat(bold(Theta))_"MAP"$ is a fixed point of MAP-EM, it simultaneously maximizes the observed log-posterior $log p(bold(Theta) | cal(D))$ and the fractional complete-data log-posterior $log p(bold(Theta) | cal(D), hat(bold(Z)))$! Applying the Laplace approximation to both posteriors around their shared mode $hat(bold(Theta))_"MAP"$ gives
  $ log p(hat(bold(Theta))_"MAP" | cal(D), hat(bold(Z))) &= 1 / 2 log det(bold(H)_"comp") - (p_K / 2) log(2 pi) + cal(O)(n^(-1)), \
    log p(hat(bold(Theta))_"MAP" | cal(D)) &= 1 / 2 log det(bold(H)_"obs") - (p_K / 2) log(2 pi) + cal(O)(n^(-1)). $
  Subtracting the second line from the first and taking logs in the exact ratio proves Claim 1 (accounting for the $K!$ identical modes of a non-degenerate mixture). Because $bold(H)_"miss" succ.eq 0$, $bold(H)_"obs" prec.eq bold(H)_"comp"$, hence $det(bold(H)_"comp" bold(H)_"obs"^(-1)) >= 1$ and the $cal(O)(1)$ remainder is nonnegative.

  For Claim 2, applying Laplace's method directly to the complete-data integral gives $log p(cal(D), hat(bold(Z)) | cal(M)_K) = log p(cal(D), hat(bold(Z)) | hat(bold(Theta))_"MAP") + log p(hat(bold(Theta))_"MAP") - (1/2) log det(bold(H)_"comp") + (p_K / 2) log(2 pi) + cal(O)(n^(-1))$. Because $bold(H)_"comp" = n bold(I)_"comp" + cal(O)(1)$ for a positive-definite per-sample complete Fisher information matrix $bold(I)_"comp"$, we have $log det(bold(H)_"comp") = p_K log n + cal(O)(1)$. Substituting into Definition 1 cancels $log p(cal(D), hat(bold(Z)) | hat(bold(Theta))_"MAP")$ and leaves $log p(cal(D) | hat(bold(Theta))_"MAP") - (p_K / 2) log n + cal(O)(1)$.
]

Theorem 2 explains why AutoClass avoids overfitting $K$ without ad-hoc pruning thresholds: the conjugate marginal integral automatically embeds an Occam penalty of $-(p_K / 2) log n$ plus finite-sample prior volume corrections, while closely approximating the exact Bayesian evidence.

#subsec("D", "Exact Marginalization Under Missing-at-Random (MAR) Data")

#proposition(title: "Unbiased MAR Marginalization in Mixed-Type MAP-EM", origin: "Original Proof")[
  If instance $bold(x)_i$ has observed indices $cal(O)_i$ and Missing-at-Random (MAR) indices $cal(U)_i$ ($p(bold(m)_i | bold(x)_i) = p(bold(m)_i | bold(x)_(i, cal(O)_i))$), then marginalizing $bold(x)_(i, cal(U)_i)$ from the factorized mixture yields exact E-step responsibilities over $cal(O)_i$ alone and exact M-step updates where missing entries contribute zero weight to attribute sufficient statistics.
]

#proof[
  Under MAR, the missingness mask $bold(m)_i$ factors out of the likelihood of $bold(Theta)$. Integrating out $bold(x)_(i, cal(U)_i)$ inside component $k$ yields
  $ p(bold(x)_(i, cal(O)_i) | z_i = k, bold(Theta)) &= integral product_(j in cal(O)_i) p_k (x_(i, j)) product_(u in cal(U)_i) p_k (x_(i, u)) dif bold(x)_(i, cal(U)_i) \
    &= product_(j in cal(O)_i) p_k (x_(i, j)), $
  as each factorized univariate density integrates to $1$. Hence $p(z_i = k | bold(x)_(i, cal(O)_i), bold(Theta))$ omits all $u in cal(U)_i$, and optimizing the exact observed-data marginal drops missing terms from $N_(k, a)$ and $N_(k, b)$ without mean-imputation variance shrinkage.
]

#subsec("E", "Local Convergence Rate and Mode Multiplicity")

#theorem(title: "Stationary Convergence and Spectral Contraction of MAP-EM", origin: "Classical — Cited")[
  Under Theorem 1's variance floor, every limit point of AutoClass MAP-EM is a stationary point of $cal(J)(bold(Theta))$ @wu1983convergence @dempster1977em. Near a regular maximizer $hat(bold(Theta))_"MAP"$, convergence is locally linear with contraction matrix $bold(M) = bold(H)_"comp"^(-1) bold(H)_"miss"$ and rate $rho(bold(M)) in [0, 1)$ @dempster1977em @mclachlan2000finite. Any solution with $K$ distinct components belongs to an equivalence class of $K!$ permutation-symmetric modes @cheeseman1996autoclass.
]

#proof[
  By Theorem 1, the upper level set ${bold(Theta) : cal(J)(bold(Theta)) >= cal(J)(bold(Theta)^((0)))}$ has strictly positive variance floors and compact simplex weights; continuity of the MAP M-step and strict ascent at non-stationary points satisfy Wu's global convergence conditions @wu1983convergence. Linearizing $bold(Theta)^((t+1)) = cal(M)(bold(Theta)^((t)))$ around $hat(bold(Theta))_"MAP"$ gives Jacobian $D cal(M) = bold(H)_"comp"^(-1) bold(H)_"miss"$ @dempster1977em @bishop2006pattern; because $0 prec.eq bold(H)_"miss" prec bold(H)_"comp"$, all eigenvalues of $bold(H)_"comp"^(-1/2) bold(H)_"miss" bold(H)_"comp"^(-1/2)$ lie in $[0, 1)$, proving $rho(bold(M)) < 1$ and tying slow convergence directly to Theorem 2's overlap gap. Finally, invariance of the mixture and symmetric priors under all $sigma in S_K$ generates $K!$ identically valued modes whenever all $K$ components are distinct.
]

#sec("III", "Complexity, Geometry and Operational Profile")

For conditionally independent attributes, each MAP-EM iteration requires $cal(O)(n K (d_c + d_m))$ log-density evaluations in the E-step and $cal(O)(n K (d_c + d_m) + K d_m V)$ sufficient-statistic accumulations in the M-step. Across $R$ random restarts (and candidate $K$ values) of at most $T$ iterations, total runtime is $cal(O)(R T n K (d_c + d_m V))$, rising to $cal(O)(R T n K d_c^2 + R T K d_c^3)$ when full covariant Gaussian blocks are enabled. Memory scales cleanly as $cal(O)(n K + K (d_c + d_m V))$, storing only the responsibility matrix and compact conjugate statistics.

Because AutoClass searches across multiple restarts and model sizes using Bayesian evidence, it trades single-pass speed (slower than COBWEB and CLASSIT) for global probabilistic calibration, automatic $K$ selection, and smooth soft-membership posteriors.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Input Domain], [Mixed continuous, categorical, and MAR-missing attributes natively],
    [Objective Function], [Bayesian MAP-EM + Cheeseman-Stutz log marginal likelihood $log hat(p)_"CS"$],
    [Time Complexity], [$cal(O)(R T n K (d_c + d_m V))$ diagonal; $cal(O)(R T n K d_c^2)$ full covariant],
    [Space Complexity], [$cal(O)(n K + K (d_c + d_m V))$],
    [Number of Clusters $K$], [Selected automatically via Bayesian evidence (asymptotically BIC-consistent)],
    [Invariances], [Order-invariant in $n$; continuous diagonal terms invariant to affine coordinate rescaling under empirical-Bayes priors],
    [Singularity Control], [Guaranteed positive variance floor via Inverse-Gamma conjugate priors],
    [Parallelism], [Embarrassingly parallel over restarts $R$, candidate $K$, and E-step samples $n$],
  ),
  caption: [Operational and geometric profile of AutoClass.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

*1. Astronomical Spectral Taxonomy (NASA IRAS Catalog):* AutoClass's landmark scientific triumph @cheeseman1996autoclass classified $5{,}425$ Infrared Astronomical Satellite (IRAS) low-resolution spectra into previously unrecognized astrophysical star classes—separating dusty silicate-emission subclasses that astronomers later confirmed optically. Continuous flux intensities across spectral bands with measurement noise matched AutoClass's Gaussian likelihood and automatic Bayesian $K$ discovery directly.

*2. Mixed Clinical-Genomic Patient Stratification:* Precision-medicine cohorts combine continuous biomarkers (lab panels, imaging radiomics) with categorical genotypes, clinical stages, and pervasive missing lab tests. AutoClass's mixed Dirichlet-Normal-Inverse-Gamma model integrates both modalities on a single commensurable log-probability scale without ad-hoc Gower distance weighting, while Proposition 1 ensures unbiased MAR handling.

*3. Bioinformatics and Molecular Sequence Subtyping:* AutoClass has been applied to intron-exon splice sites, protein structural motifs, and multi-omics tables @cheeseman1996autoclass, where discrete residue alphabets coexist with continuous physicochemical descriptors (hydrophobicity, solvent accessibility) and Bayesian evidence guards against splitting scarce biological subtypes into noise.

*4. Industrial Reliability and Mixed Telemetry Forensics:* Aviation and manufacturing failure logs couple continuous sensor readings (temperature, vibration, pressure) with discrete operational modes and fault flags. AutoClass uncovers latent failure regimes with calibrated class probabilities $r_(i k)$, enabling probabilistic root-cause triage across heterogeneous fleet records.

*5. Financial, Insurance, and Demographic Risk Segmentation:* Underwriting datasets combine continuous income, exposure, and claim-amount attributes with discrete occupation, region, and policy-tier indicators. AutoClass yields generative risk archetypes whose conjugate posteriors naturally regularize small niche segments and score novel applicants via exact predictive likelihoods.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

AutoClass encounters four key limitations: (i) *heavy-tailed outliers*, because Gaussian continuous components have quadratic log-loss tails that distort component means and inflate variance estimates under extreme corruptions; (ii) *high-dimensional small-$n$ settings*, where naive independence misses dense correlations while full-covariance blocks suffer $cal(O)(K d_c^2)$ parameter explosion; (iii) *non-convex manifolds*, since unimodal Gaussian components over-split rings, spirals, or crescents into multiple elliptical pieces; and (iv) *local EM traps*, requiring sufficient random restarts $R$ when classes overlap heavily ($rho(bold(M)) -> 1$ in Theorem 3).

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Mixed continuous + categorical tables], [Optimal: principled joint likelihood and exact MAR handling], [Check covariant blocks for correlated features],
    [Scientific discovery with unknown $K$], [High: automatic Bayesian evidence selection without tuning $K$], [Run multiple restarts $R$ to escape local MAP modes],
    [Outlier-contaminated continuous data], [Moderate-Low: Gaussian tails sensitive to extreme glitches], [Use Student's $t$-mixtures (Alg. 88) or trimmed models],
    [Non-convex / manifold geometry], [Low: Gaussian components fragment curved shapes], [Use spectral clustering, HDBSCAN, or density modes],
    [Ultra-high-speed single-pass streams], [Low: batch multi-restart EM unsuitable for microsecond latency], [Use CLASSIT (Alg. 86), BIRCH, or streaming VB],
  ),
  caption: [Domain suitability matrix and failure mitigations for AutoClass.],
)

#sec("V", "Structural Relations and Empirical Behavior")

AutoClass occupies a unifying hub in Category IV statistical clustering.

#proposition(title: "Unification of Classical Latent-Variable Mixture Models", origin: "Original Derivation")[
  AutoClass reduces strictly to: (i) Bayesian Latent Class Analysis @lazarsfeld1968latent when $d_c = 0$; (ii) regularized diagonal or covariant Gaussian Mixture Models (MCLUST @fraley2002mclust) with BIC-consistent model selection when $d_m = 0$ (Theorem 2); and (iii) a finite-dimensional truncation of the Dirichlet Process / Infinite Gaussian Mixture Model @rasmussen2000igmm @blei2006variational as $K$ grows under symmetric prior $"Dir"(alpha_0 / K, dots, alpha_0 / K)$.
]

#proof[
  Setting $d_c = 0$ leaves only the Dirichlet-Multinomial product $product_a theta_(k, a, x_(i, a))$, recovering Bayesian Latent Class Analysis. Setting $d_m = 0$ leaves Gaussian components with conjugate NIG / NIW priors; by Theorem 2, $log hat(p)_"CS"$ matches MCLUST's BIC criterion up to $cal(O)(1)$ while adding finite-sample conjugate regularization (Theorem 1). Finally, Kingman's representation implies that the symmetric Dirichlet prior $"Dir"(alpha_0 / K, dots, alpha_0 / K)$ converges weakly in partition distribution to the Chinese Restaurant Process $"CRP"(alpha_0)$ as $K -> infinity$, linking AutoClass's automatic empty-component pruning directly to nonparametric DP mixtures @rasmussen2000igmm.
]

Compared to COBWEB (Alg. 85) and CLASSIT (Alg. 86), AutoClass eliminates presentation-order dependence and acuity non-convexity while unifying discrete and continuous features in one Bayesian framework; MIXMOD @biernacki2006mixmod and Wallace's MML-based SNOB @wallace1968snob pursue closely related mixed-type penalized-likelihood objectives.

*Empirical Behavior on the ClusterOverview 4-Dataset Benchmark:* Across the 4-dataset benchmark, AutoClass achieves Iris Silhouette $= 0.431$ and $"ARI" = 0.762$ @hubert1985comparing @rousseeuw1987silhouettes, Customers Silhouette $= 0.316$, Multi-Shape Silhouette $= 0.782$ ($"ARI" =$ N/A), and High-D Silhouette $= 0.321$ ($"ARI" =$ N/A), with a mean runtime of $19.192$ ms. The leap over Algs. 85–86 reflects our theoretical guarantees directly:
1. *Global MAP-EM vs. Greedy Trees:* Iris $"ARI"$ surges more than elevenfold from $0.066$ (COBWEB / CLASSIT) to $0.762$ ($"Sil" = 0.431$), while Multi-Shape reaches $"Sil" = 0.782$ and High-D leads the four-algorithm cohort at $"Sil" = 0.321$. Global likelihood optimization, conjugate variance stabilization (Theorem 1), and Bayesian evidence selection (Theorem 2) prevent both greedy fragmentation and variance collapse.
2. *Sensitivity to Heavy Tails and Oblique Overlap:* On heavy-tailed wholesale *Customers* ($"Sil" = 0.316$) and overlapping *Iris* species ($"ARI" = 0.762$), Gaussian tail assumptions remain vulnerable to extreme spenders and borderline points—precisely the regime where Student's $t$-Mixture Models (Alg. 88) achieve higher robustness ($"Sil" = 0.367$ on Customers; $"ARI" = 0.950$ on Iris).
3. *Moderate Runtime Overhead:* At $19.192$ ms, multi-restart MAP-EM costs roughly $8 times$ CLASSIT's runtime while delivering markedly superior clustering fidelity across all four datasets.
