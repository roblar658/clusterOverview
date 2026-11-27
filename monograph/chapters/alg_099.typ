#import "../template/ieee.typ": *

#alg_entry(
  id: 99,
  name: "Infinite Gaussian Mixture Model (IGMM)",
  category: "Cat. IV: Distribution, Statistical & Topic",
  year: "2000",
  paradigm: "Bayesian Nonparametric DP Mixture (CRP)",
  complexity: $cal(O)(T K_"eff" n d^2)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

A central statistical bottleneck of finite Gaussian Mixture Models
@mclachlan2000finite @fraley2002mclust is specifying the cluster count $K$
beforehand or relying on asymptotic criteria such as BIC @schwarz1978bic that
ignore parameter uncertainty across competing cardinalities. The *Infinite
Gaussian Mixture Model (IGMM)*, introduced by Rasmussen @rasmussen2000igmm as
the $K -> infinity$ limit of Bayesian finite Gaussian mixtures and formally
grounded in the Dirichlet Process Mixture Model (DPMM) @blei2006variational
@teh2006hdp, places a nonparametric prior over all countable partitions of
$cal(X) = {bold(x)_1, ..., bold(x)_n} subset.eq bb(R)^d$, letting the active
cluster count $K_"eff"$ adapt organically with the data.

Let $H$ be a conjugate Normal-Inverse-Wishart base measure
$H = "NIW"(bold(mu)_0, kappa_0, nu_0, bold(Psi)_0)$ over Gaussian parameters
$theta = (bold(mu), bold(Sigma)) in bb(R)^d times bb(S)_(++)^d$, with prior
mean $bold(mu)_0 in bb(R)^d$, strength $kappa_0 > 0$, degrees of freedom
$nu_0 > d - 1$, and positive-definite scale matrix $bold(Psi)_0 in bb(S)_(++)^d$:
$
  bold(Sigma) &~ "IW"(nu_0, bold(Psi)_0), \
  bold(mu) | bold(Sigma) &~ cal(N)_d (bold(mu)_0, bold(Sigma) / kappa_0).
$
Given concentration parameter $alpha > 0$, the IGMM generative hierarchy draws
a random discrete probability measure $G$ from a Dirichlet Process
$G ~ "DP"(alpha, H)$ and generates each observation via:
$
  theta_i | G &~ G, \
  bold(x)_i | theta_i &~ cal(N)_d (bold(mu)_i, bold(Sigma)_i), quad i = 1, ..., n.
$
Because draws from $G$ are almost surely discrete ($G = sum_(k=1)^infinity pi_k delta_(theta_k^*)$
with stick-breaking weights $pi_k = v_k product_(j=1)^(k-1) (1 - v_j)$,
$v_k ~ "Beta"(1, alpha)$ @blei2006variational), multiple observations share
identical atoms $theta_k^*$, inducing a random clustering $bold(z)$.

#definition(title: "Chinese Restaurant Process (CRP) Partition Prior")[
  Integrating out the infinite random measure $G ~ "DP"(alpha, H)$ yields the
  *Chinese Restaurant Process* predictive distribution over cluster assignments
  $bold(z) = (z_1, ..., z_n)$: conditioned on $bold(z)_(-i)$ with $K_+$
  occupied clusters of sizes $n_k^(-i) > 0$,
  $
    p(z_i = k | bold(z)_(-i), alpha) = cases(
      n_k^(-i) / (n - 1 + alpha) & "if " k in {1, ..., K_+},
      alpha / (n - 1 + alpha) & "if " k = K_+ + 1 " (new cluster)."
    )
  $
  The induced Exchangeable Partition Probability Function (EPPF) over the
  unlabeled partition $cal(C) = {C_1, ..., C_(K_+)}$ of ${1, ..., n}$ is
  $
    p(cal(C) | alpha) =
    (alpha^(K_+) product_(k=1)^(K_+) (n_k - 1)!)
    / (product_(j=0)^(n-1) (alpha + j)).
  $
]

Conjugacy of $H = "NIW"$ enables integrating out all component means and
covariances analytically, yielding a collapsed Gibbs sampler whose posterior
predictive densities are exact multivariate Student-$t$ distributions.

#algbox(
  title: "Collapsed Gibbs Sampler for Conjugate IGMM",
  inputs: [Data $cal(X) = {bold(x)_i}_(i=1)^n subset.eq bb(R)^d$, prior
    $(alpha, bold(mu)_0, kappa_0, nu_0, bold(Psi)_0)$, sweeps $T$],
  outputs: [Posterior partition $cal(C) = {C_1, ..., C_(K_"eff")}$ and cluster
    parameters],
)[
  + *Initialize:* Assign all points to $K_+ = 1$ cluster (or random $K_0$
    clusters); cache sufficient statistics $(n_k, bar(bold(x))_k, bold(Psi)_k^*)$
    and Cholesky factors $bold(L)_k$ of $bold(Psi)_k^*$.
  + *For* sweep $t = 1, ..., T$ *do:*
    + *For* $i = 1, ..., n$ (in random order) *do:*
      + Remove $bold(x)_i$ from cluster $c = z_i$ via rank-1 Cholesky downdate
        of $bold(Psi)_c^*$; if $n_c$ hits $0$, delete cluster $c$ and decrement $K_+$.
      + For each occupied $k = 1, ..., K_+$, evaluate Student-$t$ predictive
        $p_k = n_k^(-i) dot.op t_(eta_k) (bold(x)_i | bold(mu)_k^*, bold(Lambda)_k^*)$
        (Theorem 2).
      + Evaluate new-cluster weight
        $p_"new" = alpha dot.op t_(eta_0) (bold(x)_i | bold(mu)_0, bold(Lambda)_0)$
        with $eta_0 = nu_0 - d + 1$, $bold(Lambda)_0 = ((kappa_0 + 1) / (kappa_0 eta_0)) bold(Psi)_0$.
      + Normalize and sample $z_i ~ "Cat"(p_1, ..., p_(K_+), p_"new")$.
      + If $z_i = K_+ + 1$, create a new cluster; insert $bold(x)_i$ into $z_i$
        via rank-1 Cholesky update of $bold(Psi)_(z_i)^*$.
    + Optionally resample $alpha$ via Escobar-West auxiliary Beta augmentation.
  + *Return:* MAP / post-burn-in consensus partition $cal(C)$ and $K_"eff" = K_+$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

#subsec("A", "Exact Finite-to-Infinite Limit and Partition Exchangeability")

We prove constructively that the CRP partition prior of Definition 1 arises as
the exact limit of finite symmetric Dirichlet mixtures @rasmussen2000igmm.

#theorem(
  title: "Finite-to-Infinite Limit and Exchangeability of the CRP",
  origin: "Original Proof",
)[
  Consider a finite $K$-component mixture with mixing prior
  $bold(pi) ~ "Dir"_K (alpha / K, ..., alpha / K)$ and assignments
  $z_i | bold(pi) ~ "Cat"(bold(pi))$. Let $p_K (cal(C) | alpha)$ denote the
  exact probability of an unlabeled partition $cal(C) = {C_1, ..., C_(K_+)}$
  with sizes $n_k = |C_k| >= 1$ ($sum_(k=1)^(K_+) n_k = n$). Then for every
  fixed $n$ and $cal(C)$,
  $
    lim_(K -> infinity) p_K (cal(C) | alpha) =
    (alpha^(K_+) product_(k=1)^(K_+) (n_k - 1)!)
    / (product_(j=0)^(n-1) (alpha + j)).
  $
  Consequently, the sequence $(z_1, ..., z_n)$ is infinitely exchangeable and
  admits the exact single-site CRP conditional of Definition 1 regardless of
  observation ordering.
]

#proof[
  For a fixed labeled vector $bold(z) in {1, ..., K}^n$ inducing $K_+ <= K$
  non-empty clusters, integrating out $bold(pi) ~ "Dir"_K (alpha / K, ..., alpha / K)$
  via the multivariate Beta integral yields
  $
    p_K (bold(z) | alpha) =
    Gamma(alpha) / Gamma(n + alpha)
    product_(k: n_k > 0) Gamma(n_k + alpha / K) / Gamma(alpha / K),
  $
  because empty clusters contribute $Gamma(alpha / K) / Gamma(alpha / K) = 1$.
  Since $K! / (K - K_+)!$ distinct ordered assignments of labels in
  ${1, ..., K}$ produce the identical unlabeled partition $cal(C)$ (counting
  each block of $cal(C)$ distinctly by its constituent indices), summing over
  all labelings gives
  $
    p_K (cal(C) | alpha) =
    (K! / (K - K_+)!) (Gamma(alpha) / Gamma(n + alpha))
    product_(k=1)^(K_+) Gamma(n_k + alpha / K) / Gamma(alpha / K).
  $
  Using $Gamma(x + 1) = x Gamma(x)$ repeatedly, write
  $Gamma(alpha / K) = (K / alpha) Gamma(1 + alpha / K)$, so
  $
    Gamma(n_k + alpha / K) / Gamma(alpha / K)
    = (alpha / K) dot.op Gamma(n_k + alpha / K) / Gamma(1 + alpha / K).
  $
  Substituting into $p_K (cal(C) | alpha)$ and pairing $K^(-K_+)$ with the
  falling factorial gives
  $
    p_K (cal(C) | alpha) &=
    (Gamma(alpha) / Gamma(n + alpha)) alpha^(K_+)
    ( product_(j=0)^(K_+ - 1) (1 - j / K) ) \
    &quad times product_(k=1)^(K_+) Gamma(n_k + alpha / K) / Gamma(1 + alpha / K).
  $
  As $K -> infinity$, $1 - j / K -> 1$, $Gamma(1 + alpha / K) -> Gamma(1) = 1$,
  and $Gamma(n_k + alpha / K) -> Gamma(n_k) = (n_k - 1)!$. Finally,
  $Gamma(n + alpha) / Gamma(alpha) = product_(j=0)^(n-1) (alpha + j)$ yields
  the exact EPPF. Since $p(cal(C) | alpha)$ depends only on the multiset of
  cluster sizes ${n_1, ..., n_(K_+)}$, any permutation of ${1, ..., n}$ leaves
  the joint law invariant; treating $bold(x)_i$ as the last customer yields the
  exact CRP conditional $p(z_i = k | bold(z)_(-i), alpha)$ by taking the ratio
  $p(cal(C)) / p(cal(C)_(-i))$.
]

#subsec("B", "Exact Conjugate NIW Marginalization to Student-t Predictive")

#theorem(
  title: "Closed-Form Posterior Predictive of Conjugate IGMM",
  origin: "Original Derivation",
)[
  Let $cal(X)_k = {bold(x)_j : z_j = k}$ be a cluster of $n_k$ points with mean
  $bar(bold(x))_k$ and scatter
  $bold(S)_k = sum_(j in C_k) (bold(x)_j - bar(bold(x))_k)(bold(x)_j - bar(bold(x))_k)^top$.
  Under $H = "NIW"(bold(mu)_0, kappa_0, nu_0, bold(Psi)_0)$, the posterior
  $p(bold(mu)_k, bold(Sigma)_k | cal(X)_k)$ is
  $"NIW"(bold(mu)_k^*, kappa_k^*, nu_k^*, bold(Psi)_k^*)$ with
  $
    kappa_k^* &= kappa_0 + n_k, quad nu_k^* = nu_0 + n_k, \
    bold(mu)_k^* &= (kappa_0 bold(mu)_0 + n_k bar(bold(x))_k) / kappa_k^*, \
    bold(Psi)_k^* &= bold(Psi)_0 + bold(S)_k
    + (kappa_0 n_k) / kappa_k^* (bar(bold(x))_k - bold(mu)_0)
    (bar(bold(x))_k - bold(mu)_0)^top,
  $
  and the collapsed posterior predictive density for a new point $bold(x)$ is
  exact multivariate Student-$t$ with $eta_k = nu_k^* - d + 1$ degrees of freedom:
  $
    p(bold(x) | cal(X)_k) =
    t_(eta_k) (bold(x) | bold(mu)_k^*, bold(Lambda)_k^*), quad
    bold(Lambda)_k^* = (kappa_k^* + 1) / (kappa_k^* eta_k) bold(Psi)_k^*.
  $
]

#proof[
  Expanding the Gaussian log-likelihood sum via the parallel-axis identity
  $sum_(j in C_k) (bold(x)_j - bold(mu))(bold(x)_j - bold(mu))^top = bold(S)_k + n_k (bar(bold(x))_k - bold(mu))(bar(bold(x))_k - bold(mu))^top$
  and completing the quadratic form in $bold(mu)$ together with
  $kappa_0 (bold(mu) - bold(mu)_0)(bold(mu) - bold(mu)_0)^top$ produces
  $kappa_k^* (bold(mu) - bold(mu)_k^*)(bold(mu) - bold(mu)_k^*)^top + (kappa_0 n_k / kappa_k^*)(bar(bold(x))_k - bold(mu)_0)(bar(bold(x))_k - bold(mu)_0)^top$.
  Collecting exponents of $det(bold(Sigma))$ and $"Tr"(bold(Sigma)^(-1) dot.op)$
  confirms the NIW posterior parameters $(bold(mu)_k^*, kappa_k^*, nu_k^*, bold(Psi)_k^*)$.

  The NIW normalizing constant is
  $Z(kappa, nu, bold(Psi)) = (2 pi)^(d/2) kappa^(-d/2) 2^(nu d / 2) Gamma_d (nu / 2) det(bold(Psi))^(-nu / 2)$,
  where $Gamma_d (a) = pi^(d(d-1)/4) product_(j=1)^d Gamma(a + (1 - j) / 2)$.
  Adding one observation $bold(x)$ updates $kappa' = kappa_k^* + 1$,
  $nu' = nu_k^* + 1$, and
  $bold(Psi)' = bold(Psi)_k^* + (kappa_k^* / (kappa_k^* + 1)) (bold(x) - bold(mu)_k^*)(bold(x) - bold(mu)_k^*)^top$.
  By Bayes' rule, $p(bold(x) | cal(X)_k) = (2 pi)^(-d/2) Z(kappa', nu', bold(Psi)') / Z(kappa_k^*, nu_k^*, bold(Psi)_k^*)$.
  Applying the matrix determinant lemma to $det(bold(Psi)')$,
  $
    det(bold(Psi)') = det(bold(Psi)_k^*)
    ( 1 + eta_k^(-1) (bold(x) - bold(mu)_k^*)^top (bold(Lambda)_k^*)^(-1)
    (bold(x) - bold(mu)_k^*) ),
  $
  since $(kappa_k^* / (kappa_k^* + 1)) (bold(Psi)_k^*)^(-1) = eta_k^(-1) (bold(Lambda)_k^*)^(-1)$.
  Meanwhile, the multivariate Gamma ratio telescopes over $j = 1, ..., d$ to
  $
    Gamma_d ((nu_k^* + 1) / 2) / Gamma_d (nu_k^* / 2)
    = Gamma((nu_k^* + 1) / 2) / Gamma((nu_k^* + 1 - d) / 2)
    = Gamma((eta_k + d) / 2) / Gamma(eta_k / 2).
  $
  Gathering the remaining scalar prefactors yields
  $(2 pi)^(-d/2) (kappa_k^* / (kappa_k^* + 1))^(d/2) 2^(d/2) det(bold(Psi)_k^*)^(-1/2) = (eta_k pi)^(-d/2) det(bold(Lambda)_k^*)^(-1/2)$,
  while the determinant ratio enters with exponent $-(nu_k^* + 1) / 2 = -(eta_k + d) / 2$,
  matching the exact $d$-variate Student-$t$ density $t_(eta_k) (bold(x) | bold(mu)_k^*, bold(Lambda)_k^*)$.
]

#subsec("C", "Logarithmic Growth Law and Concentration of Cluster Count")

#theorem(
  title: "Exact Moments and Almost-Sure Logarithmic Cluster Growth",
  origin: "Classical — Cited",
)[
  Under the CRP prior $p(cal(C) | alpha)$ on $n$ observations, the number of
  non-empty clusters $K_+$ decomposes as a sum of independent Bernoulli random
  variables $K_+ = sum_(i=1)^n B_i$ with $B_i ~ "Bernoulli"(alpha / (alpha + i - 1))$.
  Its exact mean and variance satisfy
  $
    bb(E)[K_+] &= sum_(i=1)^n alpha / (alpha + i - 1)
    = alpha log(1 + n / alpha) + cal(O)(1), \
    "Var"(K_+) &= sum_(i=1)^n (alpha (i - 1)) / (alpha + i - 1)^2
    <= bb(E)[K_+],
  $
  and $K_+ / (alpha log n) -> 1$ almost surely as $n -> infinity$
  @rasmussen2000igmm @teh2006hdp.
]

#proof[
  In the sequential CRP construction (Definition 1), customer $i$ starts a new
  table ($B_i = 1$) with probability $alpha / (alpha + i - 1)$, a quantity
  strictly independent of the entire prior seating configuration
  $(z_1, ..., z_(i-1))$ and hence independent of $(B_1, ..., B_(i-1))$.
  Linearity of expectation and independence give the exact sums for
  $bb(E)[K_+]$ and $"Var"(K_+) = sum_i p_i (1 - p_i) <= sum_i p_i = bb(E)[K_+]$.
  Bounding the monotone sum by integrals gives
  $alpha log(1 + n / alpha) <= sum_(i=1)^n alpha / (alpha + i - 1) <= 1 + alpha log(1 + (n - 1) / alpha)$.
  To prove almost-sure convergence, note that $"Var"(B_1) = 0$ and for every
  $i >= 2$, $"Var"(B_i) <= alpha / (alpha + i - 1) <= alpha / (i - 1)$, so
  $sum_(i=2)^infinity "Var"(B_i) / (log i)^2 <= alpha sum_(i=2)^infinity ((i - 1) (log i)^2)^(-1) < infinity$.
  By Kolmogorov's strong law of large numbers for independent non-identically
  distributed variables (via Kronecker's lemma),
  $(K_+ - bb(E)[K_+]) / log n -> 0$ almost surely, whence
  $K_+ / (alpha log n) -> 1$ almost surely.
]

#subsec("D", "Exponential Truncation Bound for Stick-Breaking Variational IGMM")

In large-scale applications where MCMC mixes slowly, Blei and Jordan
@blei2006variational truncate the stick-breaking representation at level $T$
by setting $q(v_T = 1) = 1$ and running mean-field variational inference.

#theorem(
  title: "Exponential Stick-Breaking Truncation Error Bound",
  origin: "Classical — Cited",
)[
  Let $p_infinity (bold(x)_(1:n))$ be the exact IGMM marginal density of $n$
  observations and $p_T (bold(x)_(1:n))$ be the marginal density obtained by
  truncating the stick-breaking prior at $T$ components ($v_T = 1$). Then
  $
    norm(p_T (bold(x)_(1:n)) - p_infinity (bold(x)_(1:n)))_1
    <= 2 n (alpha / (alpha + 1))^T.
  $
  Hence achieving $L_1$ error at most $epsilon$ requires only
  $T >= (log(2 n / epsilon)) / log(1 + 1 / alpha)$ components
  @blei2006variational.
]

#proof[
  Couple the full and truncated processes on the same probability space by
  sharing all atoms ${theta_k^*}_(k >= 1)$, stick lengths ${v_k}_(k=1)^(T-1)$,
  and uniform indicator draws $U_1, ..., U_n ~ "Uniform"(0, 1)$ used to invert
  the cumulative stick weights $F_infinity (k) = sum_(j=1)^k pi_j$. The two
  coupled datasets $bold(X)_(1:n)^((infinity))$ and $bold(X)_(1:n)^((T))$ are
  *identically equal* on the event
  $E = {forall i in {1, ..., n} : U_i <= 1 - R_T}$, where
  $R_T = sum_(k=T+1)^infinity pi_k = product_(k=1)^T (1 - v_k)$ is the
  remaining stick mass beyond level $T$. By the coupling inequality for total
  variation,
  $
    norm(p_T - p_infinity)_1 = 2 d_"TV" (p_T, p_infinity)
    <= 2 bb(P)(E^c) = 2 bb(E)[ 1 - (1 - R_T)^n ].
  $
  Applying Bernoulli's inequality $1 - (1 - R_T)^n <= n R_T$ (for $R_T in [0, 1]$)
  and independence of $v_1, ..., v_T ~ "Beta"(1, alpha)$, each having
  $bb(E)[1 - v_k] = alpha / (alpha + 1)$, yields
  $
    norm(p_T - p_infinity)_1 <= 2 n bb(E)[R_T]
    = 2 n product_(k=1)^T bb(E)[1 - v_k]
    = 2 n (alpha / (alpha + 1))^T.
  $
  Solving $2 n (alpha / (alpha + 1))^T <= epsilon$ for $T$ gives the stated
  logarithmic truncation complexity.
]

#subsec("E", "Polynomial Tail Instability and Outlier Over-Clustering")

#proposition(
  title: "Spontaneous Singleton Creation Under Heavy-Tailed Outliers",
  origin: "Original Proof",
)[
  Suppose an occupied cluster $k$ of size $n_k$ has posterior predictive scale
  $bold(Lambda)_k^*$ while the base prior predictive has scale $bold(Lambda)_0$
  with fixed degrees of freedom $eta_0 = nu_0 - d + 1$. For an observation
  $bold(x)_i$ at large radius $r = norm(bold(x)_i - bold(mu)_k^*) -> infinity$
  along any unit direction $bold(u)$, the posterior odds of opening a brand-new
  singleton cluster versus joining cluster $k$ grow polynomially as
  $
    (p(z_i = K_+ + 1 | bold(z)_(-i), cal(X))) / (p(z_i = k | bold(z)_(-i), cal(X)))
    = Theta( (alpha / n_k) r^(n_k) ) -> infinity.
  $
  Thus even a single heavy-tailed outlier drawn from a non-Gaussian cluster
  forces IGMM to spawn spurious satellite clusters with probability tending to 1.
]

#proof[
  By Theorem 2, as $r = norm(bold(x)_i - bold(mu)_k^*) -> infinity$ along a
  fixed unit ray $bold(u)$, the quadratic forms satisfy
  $(bold(x)_i - bold(mu)_k^*)^top (bold(Lambda)_k^*)^(-1) (bold(x)_i - bold(mu)_k^*) = c_k r^2$
  and
  $(bold(x)_i - bold(mu)_0)^top bold(Lambda)_0^(-1) (bold(x)_i - bold(mu)_0) = c_0 r^2 + cal(O)(r)$
  with $c_k, c_0 > 0$. Consequently, the two Student-$t$ densities decay as
  $t_(eta_k) (bold(x)_i) = Theta(r^(-(eta_k + d)))$ and
  $t_(eta_0) (bold(x)_i) = Theta(r^(-(eta_0 + d)))$. Since
  $eta_k - eta_0 = (nu_0 + n_k - d + 1) - (nu_0 - d + 1) = n_k$, the CRP
  posterior odds ratio satisfies
  $
    p_"new" / p_k = (alpha / n_k) (t_(eta_0) (bold(x)_i) / t_(eta_k) (bold(x)_i))
    = Theta( (alpha / n_k) r^(eta_k - eta_0) )
    = Theta( (alpha / n_k) r^(n_k) ).
  $
  Because $n_k >= 1$, this ratio diverges to $+infinity$ polynomially in $r$,
  driving $p(z_i = K_+ + 1 | bold(z)_(-i), cal(X)) -> 1$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Evaluating all $K_"eff" + 1$ Student-$t$ predictive densities naively requires
$cal(O)(K_"eff" n d^3)$ time per sweep for fresh Cholesky factorizations.
Maintaining Cholesky factors $bold(L)_k$ of $bold(Psi)_k^*$ via rank-1
updates/downdates reduces each point insertion/deletion to $cal(O)(d^2)$ and
each quadratic solve $bold(L)_k^(-1) (bold(x)_i - bold(mu)_k^*)$ to
$cal(O)(d^2)$, yielding $cal(O)(T K_"eff" n d^2)$ total runtime (and
$cal(O)(T K_"eff" n d)$ under diagonal NIW priors). Space complexity is
$cal(O)(n d + K_"eff" d^2)$. Geometrically, IGMM combines full ellipsoidal
cluster shapes with an adaptive heavy-tailed Student-$t$ "background field"
(the new-cluster term) that carves out anomalous or distant regions into new
components automatically.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [*Property*], [*Characterization*],
    [Time Complexity], [$cal(O)(T K_"eff" n d^2)$ with rank-1 Cholesky updates; $cal(O)(T K_"eff" n d^3)$ naive],
    [Space Complexity], [$cal(O)(n d + K_"eff" d^2)$],
    [Cluster Geometry], [Adaptive full-covariance ellipsoids + Student-$t$ novelty envelope],
    [Model Selection], [Automatic Bayesian nonparametric $K_"eff" ~ cal(O)(alpha log n)$ (Thm. 3)],
    [Invariances], [Order-exchangeable (Thm. 1); affine equivariant if $bold(Psi)_0$ transforms],
    [Inference Modes], [Exact collapsed Gibbs @rasmussen2000igmm; truncated VI @blei2006variational],
  ),
  caption: [Operational and complexity profile of the Infinite GMM (IGMM).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

+ *Speaker Diarization and Acoustic Sub-Unit Discovery:* Audio recordings
  contain an unknown number of speakers whose MFCC / embedding feature vectors
  form distinct ellipsoidal clouds; the CRP prior lets IGMM introduce new
  speaker clusters seamlessly whenever a previously unheard voice enters.
+ *Single-Cell Transcriptomics Rare Cell-Type Discovery:* In scRNA-seq (after
  PCA / UMAP reduction), rare stem or malignant subpopulations represent $< 1\%$
  of cells and are routinely absorbed by fixed-$K$ clustering; IGMM's novelty
  term (Proposition 1) naturally isolates low-prevalence distinct clusters.
+ *Astronomical Catalog Classification and Transient Discovery:* Photometric
  color-magnitude catalogs contain both populous main-sequence classes and rare
  unmodeled astrophysical phenomena; IGMM quantifies full posterior uncertainty
  over class count while adapting cluster covariances.
+ *Nonparametric Image Patch and Color-Space Segmentation:* Natural scenes vary
  widely in visual complexity; IGMM adapts the number of Gaussian color/texture
  modes per image automatically without hand-tuning $K$.
+ *Zero-Day Network Anomaly and Behavioral Clustering:* Benign traffic modes
  occupy tight learned ellipsoids while novel intrusion patterns fall outside
  existing $t_(eta_k)$ tails, triggering automatic creation of dedicated
  incident clusters.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

IGMM exhibits four critical failure regimes: (i) *heavy-tailed over-splitting*
(Proposition 1), where mild non-Gaussian tails or outliers trigger logarithmic
proliferation of singleton clusters; (ii) *high-dimensional NIW prior
dominance*, since requiring $nu_0 > d - 1$ injects $d$ pseudo-observations of
$bold(Psi)_0$ into every small cluster, washing out signal when $n_k << d$;
(iii) *MCMC sticky-mode trapping*, where single-point Gibbs moves struggle to
merge or split large clusters simultaneously without split-merge Metropolis
proposals; and (iv) *hyperparameter sensitivity* to $(alpha, bold(Psi)_0)$,
since mis-scaled $bold(Psi)_0$ either merges distinct clusters (too wide) or
shatters every cluster (too narrow).

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [*Domain*], [*Mathematical Fit*], [*Caution / Alternative*],
    [Unknown-$K$ multi-blob data], [Exact Bayesian nonparametric match; full ellipsoids], [Place hyperprior on $alpha$ and $bold(Psi)_0$],
    [Novelty / rare sub-types], [Student-$t$ novelty envelope isolates rare modes], [Filter heavy noise to prevent singletons (Prop. 1)],
    [Very large $n$ ($> 10^5$)], [Gibbs mixes slowly; truncated VI scales well], [Use stick-breaking VI @blei2006variational (Thm. 4)],
    [High-D ($d > 50$, small $n$)], [NIW constraint $nu_0 > d - 1$ over-regularizes], [Use diagonal NIW or PCA pre-projection],
    [Heavy-tailed / skewed data], [Over-splits skewed or $t$-tailed clusters], [Use Skew-Normal @lin2009skewnormal or Student-$t$ @peel2000robust],
  ),
  caption: [Domain suitability matrix for Infinite Gaussian Mixture Model.],
)

#sec("V", "Structural Relations and Empirical Behavior")

IGMM forms the continuous nonparametric bridge across Category IV: it is the
exact $K -> infinity$ limit of finite Bayesian GMMs @mclachlan2000finite
@fraley2002mclust (Theorem 1), the continuous-Gaussian counterpart of a
Dirichlet-Process DMM (Alg. 97 @yin2014dmm), and the single-corpus special
case of the Hierarchical Dirichlet Process (HDP @teh2006hdp), which couples
multiple IGMMs through a shared discrete top-level measure $G_0 ~ "DP"(gamma, H)$.
Compared to heuristic $K$-selection methods such as X-means @pelleg2000xmeans
and G-means @hamerly2003gmeans, which perform greedy binary splits using BIC
or Anderson-Darling normality tests, IGMM maintains a coherent exchangeable
joint distribution over all partitions and propagates full covariance
uncertainty via exact Student-$t$ marginals (Theorem 2).

On the ClusterOverview 4-dataset benchmark, IGMM attains Iris Silhouette
$= 0.415$ and ARI $= 0.563$ @rousseeuw1987silhouettes @hubert1985comparing,
*Customers Silhouette $= 0.381$*, *Multi-Shape Silhouette $= 0.783$* (ARI N/A),
High-D Silhouette $= 0.284$ (ARI N/A), and a mean runtime of $11.610$ ms.
Strikingly, IGMM records the highest Silhouette among Algs. 97--100 on both
*Multi-Shape ($0.783$)* and *Customers ($0.381$)* while nearly matching DMM on
High-D ($0.284$ vs. $0.287$). This empirical advantage stems directly from its
nonparametric CRP adaptation: on Multi-Shape, IGMM allocates the exact number
of full-covariance components needed to tile heterogeneous geometric clusters
without forcing an erroneous fixed $K$, and on Customers it separates main
spending tiers and upper-tail outliers cleanly via its Student-$t$ novelty
envelope. On Iris, its moderate ARI ($0.563$) and Silhouette ($0.415$) reflect
posterior uncertainty (mild splitting or merging) across the overlapping
Versicolor--Virginica boundary under fixed NIW hyperparameters, while its
$11.610$ ms runtime reflects iterative collapsed Gibbs sweeps with $cal(O)(d^2)$
covariance updates—faster than token-level DMM ($16.789$ ms) but slower than
closed-form matrix ECM Skew-Normal ($5.531$ ms).
