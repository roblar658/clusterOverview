#import "../template/ieee.typ": *

#alg_entry(
  id: 95,
  name: "Poisson Mixture Models",
  category: "Cat. IV: Distribution, Statistical & Topic",
  year: "1977",
  paradigm: "Discrete Count-Data Finite Mixture (I-Divergence EM)",
  complexity: $cal(O)(T n K d)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Count-valued datasets—ranging from RNA-sequencing gene counts and document word frequencies to network packet arrivals and photon counts—live on the non-negative integer lattice $bb(Z)_(>= 0)^d$ and exhibit intrinsic heteroscedasticity where variance scales with intensity. Fitting continuous Gaussian mixtures (Algs. 93--94) to count vectors ignores discrete support, wastes probability mass on negative real values, and incurs $cal(O)(K d^2)$ covariance overhead. *Poisson Mixture Models (PMMs)* @dempster1977em @mclachlan2000finite model multivariate count observations directly as a convex combination of conditionally independent Poisson distributions, inducing a natural generalized Kullback-Leibler (I-divergence) geometry @banerjee2005clustering with strictly linear $cal(O)(K d)$ parameterization.

Let $X = {bold(x)_1, dots, bold(x)_n} subset bb(Z)_(>= 0)^d$, where $x_(i j) in {0, 1, 2, dots}$ records the count of feature $j$ in sample $i$. Under a $K$-component PMM with mixing proportions $bold(pi) = (pi_1, dots, pi_K) in Delta_(K-1)$ and component rate vectors $bold(lambda)_k = (lambda_(k 1), dots, lambda_(k d)) in bb(R)_(> 0)^d$, the marginal probability mass function (PMF) is
$ p(bold(x) | Theta) = sum_(k=1)^K pi_k product_(j=1)^d (lambda_(k j)^(x_j) e^(-lambda_(k j))) / (x_j !). $
When sample exposures or library sizes $s_i > 0$ vary across observations, the offset-adjusted extension replaces $lambda_(k j)$ by $s_i lambda_(k j)$ without altering convexity.

#definition(title: "Generalized I-Divergence Representation of Poisson Mixtures")[
  Define the strictly convex negative-entropy Legendre generator $phi(bold(u)) = sum_(j=1)^d (u_j log u_j - u_j)$ on $bb(R)_(>= 0)^d$ (with $0 log 0 = 0$) and its induced Bregman divergence—the *generalized Kullback-Leibler (I-) divergence* @banerjee2005clustering:
  $ D_I (bold(x) || bold(lambda)) = sum_(j=1)^d [ x_j log(x_j / lambda_j) - x_j + lambda_j ]. $
  Each Poisson component log-PMF admits the exact Bregman form
  $ log p_k (bold(x) | bold(lambda)_k) = -D_I (bold(x) || bold(lambda)_k) + c(bold(x)), $
  where $c(bold(x)) = sum_(j=1)^d (x_j log x_j - x_j - log(x_j !))$ depends solely on $bold(x)$ and not on $Theta$.
]

Maximizing the observed log-likelihood $ell(Theta) = sum_(i=1)^n log p(bold(x)_i | Theta)$ (optionally regularized by conjugate $"Gamma"(alpha + 1, beta)$ priors with $alpha, beta >= 0$) proceeds via an exact, closed-form EM iteration @dempster1977em requiring zero matrix inversions. Model order $K$ is selected via BIC, $"BIC"(K) = 2 ell(hat(Theta)_K) - (K d + K - 1) log n$ @schwarz1978bic.

#algbox(
  title: "EM Algorithm for Multivariate Poisson Mixture Models",
  inputs: [Count matrix $X in bb(Z)_(>= 0)^(n times d)$, cluster count $K$, optional Gamma smoothing $(alpha, beta) >= (0, 0)$, tolerance $epsilon_"tol" > 0$.],
  outputs: [Mixing weights $hat(bold(pi))$, rate matrix $hat(bold(Lambda)) in bb(R)_(> 0)^(K times d)$, responsibilities $hat(tau)$, MAP partition $hat(bold(c))$.],
)[
  1. *Initialization*: Seed $bold(pi)^((0)) = (1 / K, dots, 1 / K)$ and positive rates $bold(lambda)_k^((0)) in bb(R)_(> 0)^d$ via $k$-means++ / randomized sample perturbations plus small floor $delta > 0$.
  2. *E-Step (Log-Sum-Exp Stabilized)*: For each $i in {1, dots, n}$ and $k in {1, dots, K}$, compute unnormalized log-posteriors (dropping the $k$-invariant term $-sum_j log(x_(i j)!)$):
     $ eta_(i k)^((t)) = log pi_k^((t)) + sum_(j=1)^d [ x_(i j) log lambda_(k j)^((t)) - lambda_(k j)^((t)) ], $
     and normalize via $tau_(i k)^((t)) = exp(eta_(i k)^((t)) - "lse"(eta_(i 1)^((t)), dots, eta_(i K)^((t))))$.
  3. *M-Step (Closed-Form Weighted Means)*: Compute $n_k^((t)) = sum_(i=1)^n tau_(i k)^((t))$ and update for all $k in [K], j in [d]$:
     $ pi_k^((t+1)) = n_k^((t)) / n, quad lambda_(k j)^((t+1)) = (sum_(i=1)^n tau_(i k)^((t)) x_(i j) + alpha) / (n_k^((t)) + beta). $
  4. *Convergence*: Repeat Steps 2--3 until $|ell(Theta^((t+1))) - ell(Theta^((t)))| < epsilon_"tol" (1 + |ell(Theta^((t)))|)$, then assign $hat(c)_i = arg max_k tau_(i k)^*$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

#subsec("A", "Exact Three-Divergence EM Ascent Identity")

Rather than relying solely on the generic non-negativity inequality of EM @dempster1977em @wu1983convergence, the exponential-family Bregman structure of Poisson mixtures yields an *exact* closed-form Lyapunov ascent identity quantifying every bit of log-likelihood gain.

#theorem(title: "Exact Three-Divergence EM Ascent Identity for PMM", origin: "Original Unified Proof")[
  Let $Theta^((t)) = (bold(pi)^((t)), bold(Lambda)^((t)))$ be any strictly interior iterate of unregularized PMM-EM ($alpha = beta = 0$) with $n_k^((t)) > 0$ and $lambda_(k j)^((t+1)) > 0$. Then:
  1. The weighted sample mean $bold(lambda)_k^((t+1)) = (1 / n_k^((t))) sum_(i=1)^n tau_(i k)^((t)) bold(x)_i$ is the unique global maximizer of $Q(Theta | Theta^((t)))$ over $bb(R)_(> 0)^(K times d)$.
  2. Each EM step increases the observed log-likelihood by the exact non-negative sum of three information divergences:
     $ ell(Theta^((t+1))) - ell(Theta^((t))) &= sum_(k=1)^K n_k^((t)) D_I (bold(lambda)_k^((t+1)) || bold(lambda)_k^((t))) \
       &quad + n D_"KL" (bold(pi)^((t+1)) || bold(pi)^((t))) \
       &quad + sum_(i=1)^n D_"KL" (tau_i^((t)) || tau_i^((t+1))) >= 0. $
  3. Consequently, $ell(Theta^((t+1))) = ell(Theta^((t)))$ if and only if $Theta^((t+1)) = Theta^((t))$, and the cumulative squared step movement converges: $sum_(t=0)^infinity sum_(k=1)^K n_k^((t)) D_I (bold(lambda)_k^((t+1)) || bold(lambda)_k^((t))) < infinity$.
]

#proof[
  Using Definition 1, the expected complete-data log-likelihood decouples completely across components and coordinates:
  $ Q(Theta | Theta^((t))) &= sum_(k=1)^K n_k^((t)) log pi_k - sum_(k=1)^K sum_(i=1)^n tau_(i k)^((t)) D_I (bold(x)_i || bold(lambda)_k) + C(X). $
  For each $(k, j)$, let $S_(k j) = sum_(i=1)^n tau_(i k)^((t)) x_(i j) = n_k^((t)) lambda_(k j)^((t+1))$. Differentiating the $(k, j)$-term $S_(k j) log lambda_(k j) - n_k^((t)) lambda_(k j)$ yields first derivative $S_(k j) / lambda_(k j) - n_k^((t)) = 0$ at $lambda_(k j) = lambda_(k j)^((t+1))$ and second derivative $-S_(k j) / lambda_(k j)^2 < 0$, establishing strict global concavity and proving part (1) (equivalently via Banerjee's Bregman centroid optimality theorem @banerjee2005clustering).

  For part (2), recall the exact EM decomposition $ell(Theta^((t+1))) - ell(Theta^((t))) = Delta Q + sum_(i=1)^n D_"KL" (tau_i^((t)) || tau_i^((t+1)))$, where $Delta Q = Q(Theta^((t+1)) | Theta^((t))) - Q(Theta^((t)) | Theta^((t)))$. Splitting $Delta Q$ into mixing-weight and rate contributions:
  - *Weight gain*: Using $n_k^((t)) = n pi_k^((t+1))$,
    $ sum_(k=1)^K n_k^((t)) log(pi_k^((t+1)) / pi_k^((t))) = n sum_(k=1)^K pi_k^((t+1)) log(pi_k^((t+1)) / pi_k^((t))) = n D_"KL" (bold(pi)^((t+1)) || bold(pi)^((t))). $
  - *Rate gain*: Substituting $S_(k j) = n_k^((t)) lambda_(k j)^((t+1))$,
    $ Delta Q_(k j) &= S_(k j) log(lambda_(k j)^((t+1)) / lambda_(k j)^((t))) - n_k^((t)) (lambda_(k j)^((t+1)) - lambda_(k j)^((t))) \
      &= n_k^((t)) [ lambda_(k j)^((t+1)) log(lambda_(k j)^((t+1)) / lambda_(k j)^((t))) - lambda_(k j)^((t+1)) + lambda_(k j)^((t)) ] \
      &= n_k^((t)) D_I (lambda_(k j)^((t+1)) || lambda_(k j)^((t))). $
  Summing over $k$ and $j$ yields the exact three-divergence identity. Because $D_I$ and $D_"KL"$ are strictly positive whenever their arguments differ, $ell(Theta^((t+1))) = ell(Theta^((t)))$ holds iff $bold(lambda)_k^((t+1)) = bold(lambda)_k^((t))$ and $bold(pi)^((t+1)) = bold(pi)^((t))$. Finally, because $p(bold(x)_i | Theta) <= 1$ for every discrete PMF, $ell(Theta) <= 0$ globally; telescoping the non-negative increments over $t = 0, dots, T$ gives $sum_(t=0)^T sum_k n_k^((t)) D_I (bold(lambda)_k^((t+1)) || bold(lambda)_k^((t))) <= -ell(Theta^((0))) < infinity$.
]

#subsec("B", "Strict Overdispersion and Latent Cross-Feature Covariance")

A single multivariate independent Poisson distribution enforces equidispersion ($"Var"(X_j) = bb(E)[X_j]$) and zero correlation ($"Cov"(X_j, X_m) = 0$). We prove the exact second-moment geometry unlocked by finite mixing.

#theorem(title: "Strict Overdispersion and Induced Covariance of Poisson Mixtures", origin: "Original Proof")[
  Let $bold(X) in bb(Z)_(>= 0)^d$ follow a $K$-component PMM with weights $pi_k > 0$ and mean vector $bar(bold(lambda)) = sum_(k=1)^K pi_k bold(lambda)_k$. Then for all coordinates $j, m in {1, dots, d}$:
  $ "Var"(X_j) = bar(lambda)_j + sum_(k=1)^K pi_k (lambda_(k j) - bar(lambda)_j)^2 >= bar(lambda)_j, $
  $ "Cov"(X_j, X_m) = sum_(k=1)^K pi_k (lambda_(k j) - bar(lambda)_j)(lambda_(k m) - bar(lambda)_m) quad (j != m). $
  Consequently, the marginal Fano factor $"FF"_j = "Var"(X_j) / bb(E)[X_j] >= 1$, with strict overdispersion ($"FF"_j > 1$) whenever component rates vary on coordinate $j$, and the marginal covariance matrix equals $"diag"(bar(bold(lambda))) + bold(B)$ where $"rank"(bold(B)) <= K - 1$.
]

#proof[
  Let $Z in {1, dots, K}$ denote the latent cluster variable with $bb(P)(Z = k) = pi_k$. By the law of iterated expectation, $bb(E)[X_j] = bb(E)[bb(E)[X_j | Z]] = sum_(k=1)^K pi_k lambda_(k j) = bar(lambda)_j$. Applying the law of total variance and using the Poisson identity $"Var"(X_j | Z = k) = lambda_(k j)$,
  $ "Var"(X_j) &= bb(E)["Var"(X_j | Z)] + "Var"(bb(E)[X_j | Z]) \
    &= sum_(k=1)^K pi_k lambda_(k j) + sum_(k=1)^K pi_k (lambda_(k j) - bar(lambda)_j)^2. $
  Dividing by $bar(lambda)_j > 0$ gives $"FF"_j = 1 + (1 / bar(lambda)_j) sum_k pi_k (lambda_(k j) - bar(lambda)_j)^2 >= 1$, strictly exceeding $1$ unless $lambda_(1 j) = dots = lambda_(K j)$. For $j != m$, conditional independence implies $"Cov"(X_j, X_m | Z = k) = 0$; applying the law of total covariance yields
  $ "Cov"(X_j, X_m) &= bb(E)["Cov"(X_j, X_m | Z)] + "Cov"(bb(E)[X_j | Z], bb(E)[X_m | Z]) \
    &= sum_(k=1)^K pi_k (lambda_(k j) - bar(lambda)_j)(lambda_(k m) - bar(lambda)_m). $
  In matrix notation, $"Cov"(bold(X)) = "diag"(bar(bold(lambda))) + bold(B)$ where $bold(B) = sum_(k=1)^K pi_k (bold(lambda)_k - bar(bold(lambda)))(bold(lambda)_k - bar(bold(lambda)))^top$. Since the centered vectors satisfy the linear relation $sum_k pi_k (bold(lambda)_k - bar(bold(lambda))) = bold(0)$, $"rank"(bold(B)) <= K - 1$.
]

Notably, Theorem 2 proves a hard structural impossibility: *no* finite Poisson mixture can represent an underdispersed count variable ($"Var"(X_j) < bb(E)[X_j]$).

#subsec("C", "Global Identifiability in Every Dimension via Generating Functions")

In stark contrast to Bernoulli mixtures (Alg. 96), which fail to be identifiable in low dimensions, finite Poisson mixtures are globally identifiable for *every* $d >= 1$ and *every* finite $K >= 1$.

#theorem(title: "Universal Identifiability of Multivariate Finite Poisson Mixtures", origin: "Classical — Cited")[
  Suppose two finite Poisson mixtures on $bb(Z)_(>= 0)^d$ with distinct component rate vectors ${bold(lambda)_k}_(k=1)^K$ and ${bold(lambda)'_r}_(r=1)^(K')$ and positive weights $pi_k, pi'_r > 0$ induce identical joint PMFs on $bb(Z)_(>= 0)^d$ @mclachlan2000finite. Then $K = K'$ and there exists a permutation $sigma in S_K$ such that $pi_k = pi'_(sigma(k))$ and $bold(lambda)_k = bold(lambda)'_(sigma(k))$ for all $k = 1, dots, K$.
]

#proof[
  Equality of the PMFs implies equality of their multivariate probability generating functions (PGFs) for all $bold(s) in bb(R)^d$:
  $ G(bold(s)) = bb(E)[product_(j=1)^d s_j^(X_j)] = sum_(k=1)^K pi_k exp(chevron.l bold(lambda)_k, bold(s) - bold(1) chevron.r). $
  Setting $bold(u) = bold(s) - bold(1) in bb(R)^d$ gives the entire-function identity
  $ sum_(k=1)^K pi_k e^(chevron.l bold(lambda)_k, bold(u) chevron.r) = sum_(r=1)^(K') pi'_r e^(chevron.l bold(lambda)'_r, bold(u) chevron.r) quad forall bold(u) in bb(R)^d. $
  Let $cal(U) = {bold(lambda)_1, dots, bold(lambda)_K} union {bold(lambda)'_1, dots, bold(lambda)'_(K')}$. The set of vectors $bold(v) in bb(R)^d$ orthogonal to any non-zero difference $bold(a) - bold(b)$ ($bold(a), bold(b) in cal(U), bold(a) != bold(b)$) is a finite union of $(d - 1)$-dimensional hyperplanes of Lebesgue measure zero; hence there exists $bold(v) in bb(R)^d$ such that the scalar projections $gamma(bold(a)) = chevron.l bold(a), bold(v) chevron.r$ are pairwise distinct over $cal(U)$. Evaluating along the ray $bold(u) = t bold(v)$ ($t > 0$) yields $sum_k pi_k e^(t gamma(bold(lambda)_k)) = sum_r pi'_r e^(t gamma(bold(lambda)'_r))$. Let $gamma^* = max_(bold(a) in cal(U)) gamma(bold(a))$. Multiplying both sides by $e^(-t gamma^*)$ and taking $t -> +infinity$ drives every strictly smaller exponent to $0$, forcing $gamma^*$ to appear on *both* sides with equal coefficients: there exist unique $k^*, r^*$ with $bold(lambda)_(k^*) = bold(lambda)'_(r^*)$ and $pi_(k^*) = pi'_(r^*)$. Canceling this matched term and repeating by induction establishes $K = K'$ and exact parameter equality up to permutation.
]

#subsec("D", "Absence of Likelihood Singularities and Zero-Rate Absorbing States")

#proposition(title: "Global Likelihood Boundedness and Boundary Absorption Prevention", origin: "Original Proof")[
  For any count dataset $X in bb(Z)_(>= 0)^(n times d)$, the PMM log-likelihood satisfies the finite universal upper bound
  $ sup_Theta ell(Theta) <= sum_(i=1)^n sum_(j=1)^d [ x_(i j) log x_(i j) - x_(i j) - log(x_(i j) !) ] <= 0. $
  However, under unregularized EM ($alpha = beta = 0$), if $lambda_(k j)^((t)) = 0$ occurs for some $(k, j)$, it forms an *absorbing boundary state*: $tau_(i k)^((t)) = 0$ for every sample with $x_(i j) > 0$, locking $lambda_(k j)^((s)) = 0$ for all $s >= t$. Introducing Gamma MAP smoothing ($alpha > 0, beta >= 0$) enforces $lambda_(k j)^((t)) >= alpha / (n + beta) > 0$, globally eliminating all absorbing faces and guaranteeing finite out-of-sample log-likelihoods.
]

#proof[
  Because $D_I (bold(x)_i || bold(lambda)_k) >= 0$ for all $bold(lambda)_k in bb(R)_(>= 0)^d$, Definition 1 gives $log p_k (bold(x)_i | bold(lambda)_k) <= c(bold(x)_i) <= 0$ (since $p_k <= 1$ for every discrete PMF). Taking the convex combination $sum_k pi_k p_k (bold(x)_i) <= max_k p_k (bold(x)_i) <= exp(c(bold(x)_i))$ and summing over $i$ yields $ell(Theta) <= sum_i c(bold(x)_i) <= 0$, ruling out Gaussian variance blow-up. For boundary absorption, if $lambda_(k j)^((t)) = 0$, the Poisson PMF $0^(x_(i j)) e^0 / x_(i j)!$ vanishes whenever $x_(i j) >= 1$, forcing $tau_(i k)^((t)) = 0$ on every positive-count observation for feature $j$. Thus $sum_(i=1)^n tau_(i k)^((t)) x_(i j) = 0$, yielding $lambda_(k j)^((t+1)) = 0$ inductively for all $s >= t$. Under Gamma MAP smoothing ($alpha > 0, beta >= 0$), $n_k^((t)) <= n$ implies $lambda_(k j)^((t+1)) = (S_(k j) + alpha) / (n_k^((t)) + beta) >= alpha / (n + beta) > 0$, strictly bounding rates away from $0$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

PMM is among the fastest statistical clustering algorithms in existence. Each E-step requires one matrix multiplication $X (log bold(Lambda))^top$ plus a broadcast subtraction of component rate sums $sum_j lambda_(k j)$, costing exactly $cal(O)(n K d)$ floating-point operations (and accelerating to $cal(O)(K ("nnz"(X) + n))$ on sparse count matrices, since $x_(i j) = 0$ contributes zero to $sum_j x_(i j) log lambda_(k j)$!). Each M-step computes $tau^top X$ in $cal(O)(n K d)$ (or $cal(O)(K "nnz"(X))$ sparse) time with *no* Cholesky or eigendecomposition steps. Across $T$ iterations, total runtime is strictly $cal(O)(T n K d)$ dense or $cal(O)(T K ("nnz"(X) + n))$ sparse.

Space complexity is $cal(O)(n K + K d)$ beyond storing $X$. Because both the E-step and M-step reduce to BLAS-3 GEMM calls ($X (log bold(Lambda))^top$ and $tau^top X$), cache locality and SIMD utilization are near-optimal. Geometrically, the log-posterior ratio $log(tau_(i k) / tau_(i m)) = chevron.l log bold(lambda)_k - log bold(lambda)_m, bold(x)_i chevron.r + (log pi_k - sum_j lambda_(k j)) - (log pi_m - sum_j lambda_(m j))$ is *affine* in $bold(x)_i$, so every MAP cluster region is a convex polytope in $bb(R)_(>= 0)^d$—equivalent to an I-divergence Voronoi cell @banerjee2005clustering.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$cal(O)(T n K d)$ dense; $cal(O)(T K ("nnz"(X) + n))$ on sparse count matrices (zero matrix inversions).],
    [Space Complexity], [$cal(O)(n K + K d)$ plus input storage (compact $K times d$ rate matrix).],
    [Decision Geometry], [Linear hyperplanes / convex polytopes in count space (generalized I-divergence Voronoi cells).],
    [Identifiability], [Globally identifiable up to label permutation for *all* $d >= 1$ and *all* finite $K >= 1$.],
    [Dispersion Regime], [Strictly overdispersed marginally ($"FF"_j >= 1$), equidispersed within each cluster.],
    [Regularization], [Gamma conjugate MAP smoothing $(alpha, beta) > 0$ prevents zero-rate absorbing states.],
  ),
  caption: [Operational, computational, and geometric profile of Poisson Mixture Models.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Single-Cell and Bulk RNA-Sequencing Count Clustering*: Gene-by-cell UMI count matrices reside natively on $bb(Z)_(>= 0)^d$ and suffer severe overdispersion across heterogeneous cell types. Offset-adjusted PMMs ($s_i lambda_(k j)$ accounting for cell sequencing depth $s_i$) directly respect discrete sampling noise without arbitrary $log(1 + x)$ variance distortion.
2. *Text Document and Bag-of-Words Corpus Clustering*: Term-frequency matrices contain sparse non-negative word counts across $d approx 10^4$--$10^5$ vocabulary entries. Sparse PMM exploits $cal(O)(T K "nnz"(X))$ updates and I-divergence geometry to group corpora by both lexical profile and document length.
3. *Network Traffic, Cybersecurity, and Event-Rate Telemetry*: Hourly packet counts, API request tallies, and intrusion-sensor hit rates follow inhomogeneous Poisson arrival processes across operational regimes (idle, burst, DDoS attack). PMMs identify latent rate states with exact likelihood-ratio anomaly scoring.
4. *Spatial Ecology and Species Abundance Profiling*: Transect counts of $d$ biological taxa across $n$ habitat plots naturally form non-negative integer vectors with strong zero abundance. Conditional Poisson mixtures uncover ecological niches while quantifying cross-species co-occurrence via Theorem 2.
5. *Astrophysical Photon Counting and Medical PET / SPECT Voxel Segmentation*: Low-count photon imaging is governed physically by quantum Poisson shot noise. PMM segmentation matches the exact physical noise likelihood where Gaussian approximations fail in low-photon regimes.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

PMMs degrade in four well-characterized regimes: (i) *within-cluster overdispersion*, where biological or bursty count variance exceeds the mean even *inside* a single true cluster (causing BIC to split one negative-binomial cluster into multiple Poisson sub-clusters); (ii) *underdispersed or continuous negative-valued data*, where PMM cannot model $"Var" < "Mean"$ and requires clipping / shifting negative inputs; (iii) *extreme zero-inflation*, where structural dropouts dominate sampling zeros (requiring Zero-Inflated Poisson mixtures); and (iv) *uncalibrated exposure differences*, where omitting library-size offsets $s_i$ causes PMM to cluster by total count $sum_j x_(i j)$ rather than relative composition.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [RNA-seq / UMI count profiles], [Exact discrete shot-noise model with $cal(O)(K "nnz"(X))$ sparsity.], [Intra-cluster overdispersion: use Negative Binomial mixtures.],
    [Bag-of-words document clustering], [I-divergence + sparse GEMM scales to large vocabularies.], [Multi-topic documents: use pLSA @hofmann1999plsa or LDA @blei2003lda.],
    [Event arrival / photon telemetry], [Physical Poisson arrival law matches low-count regimes directly.], [Burstiness / autocorrelation: use Hawkes or HMM-Poisson.],
    [Continuous signed tabular data], [Misspecified: requires non-negative rounding / shifting; loses scale.], [Use MCLUST (Alg. 93) or Mixmod (Alg. 94) Gaussian models.],
    [Binary presence / absence data], [Assigns positive mass to counts $x_(i j) >= 2$.], [Use Bernoulli Mixture Models (Alg. 96).],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for Poisson Mixture Models.],
)

#sec("V", "Structural Relations and Empirical Behavior")

PMMs sit at the exact algebraic intersection of probabilistic topic modeling, non-negative matrix factorization (NMF) @lee1999nmf @lee2001algorithms @ding2008nmfplsa, and Bregman divergence clustering @banerjee2005clustering.

#proposition(title: "Exact Equivalence Between Conditioned PMM and Multinomial Topic Mixtures", origin: "Original Derivation")[
  For each component $k$, reparameterize $bold(lambda)_k = Lambda_k bold(theta)_k$, where $Lambda_k = sum_(j=1)^d lambda_(k j) > 0$ is the total component rate and $bold(theta)_k = bold(lambda)_k / Lambda_k in Delta_(d-1)$ is the normalized composition vector. Let $N_i = sum_(j=1)^d x_(i j)$. Then every Poisson component density factorizes identically as
  $ p_k (bold(x)_i | bold(lambda)_k) = "Poisson"(N_i | Lambda_k) dot "Multinomial"(bold(x)_i | N_i, bold(theta)_k). $
  Consequently, whenever total rates are equal ($Lambda_1 = dots = Lambda_K = Lambda$), the posterior responsibilities $tau_(i k)$ of a PMM are *identically equal* to those of a Multinomial Mixture Model (single-topic pLSA / Dirichlet Multinomial Mixture @hofmann1999plsa @yin2014dmm).
]

#proof[
  Substituting $lambda_(k j) = Lambda_k theta_(k j)$ and using $sum_(j=1)^d x_(i j) = N_i$ alongside $sum_(j=1)^d theta_(k j) = 1$,
  $ p_k (bold(x)_i | bold(lambda)_k) &= product_(j=1)^d ((Lambda_k theta_(k j))^(x_(i j)) e^(-Lambda_k theta_(k j))) / (x_(i j) !) \
    &= (Lambda_k^(N_i) e^(-Lambda_k)) / (N_i !) dot (N_i ! product_(j=1)^d theta_(k j)^(x_(i j))) / (product_(j=1)^d x_(i j) !). $
  The first factor is $"Poisson"(N_i | Lambda_k)$ and the second is $"Multinomial"(bold(x)_i | N_i, bold(theta)_k)$. When $Lambda_k = Lambda$ for all $k$, the Poisson factor is independent of $k$ and cancels identically from the numerator and denominator of $tau_(i k) = pi_k p_k (bold(x)_i) / sum_m pi_m p_m (bold(x)_i)$, proving exact posterior equivalence with a Multinomial Mixture Model.
]

Similarly, hard-assignment PMM minimizes $sum_i D_I (bold(x)_i || bold(lambda)_(c_i))$, which is an exact rank-$K$ hard-clustering restriction of KL-divergence NMF $D_I (X || bold(W) bold(H))$ @lee2001algorithms @ding2008nmfplsa.

On the ClusterOverview 4-dataset empirical benchmark (Silhouette @rousseeuw1987silhouettes and ARI @hubert1985comparing), PMM displays extreme computational speed paired with clear sensitivity to distributional dispersion assumptions:
- *Iris*: $"Sil" = 0.435$, $"ARI" = 0.608$. On positive continuous morphological features adapted to rate parameters, PMM outperforms MCLUST ($"ARI" = 0.567$) and Bernoulli ($"ARI" = 0.498$) in external accuracy while trailing Mixmod's flexible elliptical Gaussian fit ($0.950$).
- *Customers*: $"Sil" = 0.189$. Lowest Silhouette among Algs. 93--96: continuous spending features exhibit strong intra-cluster variance-to-mean mismatch (violating within-cluster equidispersion), causing I-divergence hyperplanes to misalign with Euclidean customer cluster boundaries.
- *Multi-Shape*: $"Sil" = 0.712$ ($"ARI" = "N/A"$). Matches Mixmod ($0.712$) cleanly on well-separated non-negative modes, confirming that convex I-divergence polytopes effectively partition isolated positive regions.
- *High-D*: $"Sil" = 0.249$ ($"ARI" = "N/A"$). Linear $cal(O)(K d)$ parameterization prevents singular covariance failure, though continuous standardized high-dimensional features deviate from Poisson variance scaling compared to MCLUST ($0.321$) and Bernoulli ($0.303$).
- *Runtime*: Mean runtime is $3.468$ ms—the *fastest* across all four Category IV distribution algorithms ($19.7 times$ faster than MCLUST's $68.368$ ms and $3.5 times$ faster than Mixmod's $12.143$ ms)—directly reflecting its inversion-free $cal(O)(T n K d)$ closed-form updates.
