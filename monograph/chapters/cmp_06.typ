#import "../template/ieee.typ": *

#cmp_entry(
  roman: "VI",
  title: "Intra-Category Comparative Analysis: Distribution, Statistical, and Topic Models (Algs. 81-102)",
  scope_label: "Head-to-Head Theoretical & Practical Comparison of Algorithms 81-102",
  provenance: "Original Unifying Proofs + Cited Classical Theorems",
)

Part IV of this monograph develops twenty-two probabilistic, statistical, conceptual, and matrix-factorization clustering algorithms (Algs. 81--102) spanning finite Gaussian and parsimonious covariance mixtures, heavy-tailed and asymmetric continuous densities, hyperspherical directional models, Bayesian nonparametric Dirichlet processes, minimum-message-length coding, discrete count and binary latent-class mixtures, incremental conceptual hierarchies, and statistical topic models @mclachlan2000finite @bishop2006pattern @jain2010data @aggarwal2013data. Whereas prototype partitioning (Part I) minimizes deterministic geometric distortions and density clustering (Part II) extracts connected level sets, Category IV posits an explicit *generative statistical mechanism* explaining how observed data $cal(X) = {x_1, dots, x_n}$ arise from latent cluster or topic assignments $Z$. Despite wide variation across continuous Euclidean spaces, unit hyperspheres, contingency tables, sparse word counts, and online streams, all twenty-two methods instantiate a unified *variational penalized-divergence objective*:
$
  cal(F)(q, Theta, K) &= sum_(i=1)^n bb(E)_(q(z_i)) [ log p(x_i, z_i | Theta) ] \
  &quad + H(q) - Omega_"reg" (Theta, K),
$
where $q(z_i)$ is the variational posterior over latent memberships (soft responsibilities in EM/CAVI, hard Dirac indicators in CEM/GSDMM, or tree-path allocations in COBWEB/CLASSIT), $p(x_i, z_i | Theta)$ encodes the emission geometry and latent coupling structure, $H(q) = -sum_i bb(E)_q [log q(z_i)]$ is Shannon assignment entropy, and $Omega_"reg" (Theta, K)$ enforces capacity control via conjugate Bayesian priors, Dirichlet-process stick-breaking, BIC/ICL Laplace penalties, Cheeseman--Stutz marginal approximations, or Fisher-metric MML code lengths.

Four structural mathematical fault lines organize Category IV into coherent subfamilies:
1. *Emission Geometry, Tail Behavior & Covariance Parsimony:* How parametric continuous mixtures balance expressiveness against $O(K d^2)$ sample complexity, likelihood singularities, outlier breakdown, skewness, and manifold curvature (Algs. 81, 88, 89, 93, 94, 98).
2. *Nonparametric Capacity Control & Model-Order Selection:* Whether the component count $K$ is inferred via truncated variational stick-breaking, exact collapsed Gibbs Chinese Restaurant Processes, asymptotic Laplace/entropy penalties (BIC/ICL), approximate complete-data marginals, or optimal Fisher-metric quantization (Algs. 82, 87, 93, 94, 99, 102).
3. *Discrete Local Independence vs. Incremental Conceptual Utility:* How count, binary, and categorical mixtures govern algebraic identifiability and induced cross-feature covariance compared with single-pass Gini/variance conceptual trees (Algs. 85, 86, 90, 95, 96).
4. *Corpus-Level Factorizations, Topic Correlation & Short-Text Sparsity:* How document-level topic mixtures trade parameter compactness, simplex correlation flexibility, cross-document nonparametric sharing, and short-text co-occurrence concentration (Algs. 83, 84, 91, 92, 97, 100, 101).

This chapter conducts a rigorous head-to-head mathematical and algorithmic comparison of *all twenty-two* algorithms across five sections: (I) finite parametric mixtures and covariance parsimony; (II) Bayesian nonparametric mixtures and model selection; (III) discrete, count, binary, and conceptual mixtures; (IV) topic models and nonnegative matrix factorization; and (V) a complete 22-algorithm comparative reference matrix and engineering decision guide.

// =============================================================================
// SECTION I
// =============================================================================
#sec("I", "Finite Parametric Mixtures & Covariance Parsimony")

Six algorithms in Category IV model continuous vector data via finite parametric density mixtures $p(x | Theta) = sum_(k=1)^K pi_k f(x; theta_k)$ optimized by Expectation--Maximization (EM) @dempster1977em @wu1983convergence or its conditional/classification variants: unconstrained Gaussian Mixture Models (GMM, Alg. 81) @mclachlan2000finite, MCLUST (Alg. 93) @fraley2002mclust, Mixmod (Alg. 94) @biernacki2006mixmod, Student's $t$-Mixture Models (TMM, Alg. 88) @peel2000robust, Skew-Normal Mixture Models (SNMM, Alg. 98) @lin2009skewnormal, and von Mises--Fisher Mixture Models (vMF, Alg. 89) @banerjee2005vmf. We contrast these six families along three foundational axes: generic identifiability and unbounded-likelihood singularities, covariance parsimony and EM/CEM/SEM trajectory bias, and heavy-tail influence boundedness versus skewness and spherical curvature.

#subsec("A", "Identifiability & the Exact Volume-Singularity Dichotomy")

Statistical inference in mixture clustering relies on two distinct well-posedness properties: *identifiability* (whether the population density $p(x | Theta)$ uniquely determines the multiset ${(pi_k, theta_k)}_(k=1)^K$ up to label permutation) and *finite-sample MLE boundedness* (whether $sup_Theta ell(Theta) < infinity$ on a finite dataset $cal(X)$). A frequent misconception conflates these properties, yet all six continuous families are generically identifiable while splitting sharply on finite-sample singularity collapse.

Recall the Banfield--Raftery and Celeux--Govaert eigenvalue decomposition employed by MCLUST (Alg. 93) and Mixmod (Alg. 94) @fraley2002mclust @biernacki2006mixmod:
$
  Sigma_k = lambda_k D_k A_k D_k^top, \
  lambda_k = (det Sigma_k)^(1 / d) > 0, quad det A_k = 1, quad D_k in O(d),
$
where $lambda_k$ controls component *volume*, diagonal $A_k$ with unit determinant controls *shape*, and orthogonal $D_k$ controls *orientation*. Constraining each factor to be Equal (`E`), Variable (`V`), or Identity (`I`) generates the 14 parsimonious Gaussian models (`EII` through `VVV`), partitioning cleanly into seven equal-volume (`E**`) and seven variable-volume (`V**`) subfamilies.

#theorem(title: "Continuous-Mixture Identifiability & Exact Singularity Dichotomy", origin: "Original Unified Proof")[
  Let $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$ (or $bb(S)^(d-1)$ for vMF) with $K >= 2$.
  1. *(Unified Permutation Identifiability):* Finite mixtures of non-degenerate multivariate Gaussians (GMM \#81, MCLUST \#93, Mixmod \#94), multivariate Student-$t$ distributions with distinct $(mu_k, Sigma_k, nu_k)$ (TMM \#88), and von Mises--Fisher distributions on $bb(S)^(d-1)$ ($d >= 2$) with distinct natural parameters $eta_k = kappa_k mu_k$ (vMF \#89) are each *strictly identifiable up to component permutation*.
  2. *(Exact Volume-Singularity Dichotomy Across MCLUST/Mixmod & Continuous Families):* On any dataset containing at least $K$ distinct points, the unregularized log-likelihood satisfies $sup_Theta ell(Theta) = +infinity$ for unconstrained GMM (\#81), all seven variable-volume (`V**`) models of MCLUST (\#93) and Mixmod (\#94), unregularized TMM (\#88), unregularized SNMM (\#98), and unregularized vMF (\#89, $d >= 2$). Conversely, for *all seven equal-volume (`E**`) models* (`EII`, `EEI`, `EVI`, `EEE`, `EVE`, `EEV`, `EVV`), whenever no $K$ affine hyperplanes cover $cal(X)$ (holding almost surely for continuous data with $n > K d$, and requiring merely $n > K$ distinct points for `EII`), there exists a geometric separation constant $delta_0 (cal(X), K) > 0$ such that
     $
       sup_Theta ell_("E**") (Theta) <= (n d) / 2 log( (n d) / (2 pi e delta_0^2) ) < infinity.
     $
]

#proof[
  *Part 1 (Identifiability):* Identifiability up to permutation is equivalent to linear independence of distinct component densities in $L^1$. Suppose $sum_(k=1)^M a_k f_k (x) = 0$ identically for distinct component parameters and $a_k != 0$.

  - *Gaussian family (\#81, \#93, \#94):* Taking the moment-generating function $M_k (t) = integral e^(chevron.l t, x chevron.r) f_k (x) dif x = exp(chevron.l t, mu_k chevron.r + (1/2) t^top Sigma_k t)$ along ray $t = r u$ ($r > 0$) gives
    $
      sum_(k=1)^M a_k exp( (r^2 / 2) u^top Sigma_k u + r chevron.l u, mu_k chevron.r ) = 0.
    $
    For any distinct $(mu_k, Sigma_k) != (mu_j, Sigma_j)$, the algebraic set ${u in bb(S)^(d-1) : u^top (Sigma_k - Sigma_j) u = 0 "and" chevron.l u, mu_k - mu_j chevron.r = 0}$ has measure zero on $bb(S)^(d-1)$. Choosing $u$ in the complement of the finite union of these null sets renders the pairs $(u^top Sigma_k u, chevron.l u, mu_k chevron.r)$ strictly lexicographically ordered. Letting $k^*$ achieve the unique lexicographic maximum and dividing by $M_(k^*) (r u)$ sends every other term to $0$ as $r -> +infinity$, forcing $a_(k^*) = 0$, a contradiction.

  - *vMF family (\#89):* Write $f_k (x) = c_d (kappa_k) exp(chevron.l eta_k, x chevron.r)$ on $bb(S)^(d-1)$ ($d >= 2$) with pairwise distinct $eta_k = kappa_k mu_k in bb(R)^d$. The entire holomorphic function $F(z) := sum_(k=1)^M a_k c_d (kappa_k) exp(chevron.l eta_k, z chevron.r)$ vanishes on the real sphere $bb(S)^(d-1)$, which is totally real of maximal real dimension $d - 1$ inside the irreducible complex quadric $S_bb(C) := {z in bb(C)^d : sum_(j=1)^d z_j^2 = 1}$; by the identity theorem for holomorphic functions, $F(z) equiv 0$ on all of $S_bb(C)$. Pick orthonormal real vectors $u, v in bb(R)^d$ such that ${chevron.l eta_k, u chevron.r}$ are pairwise distinct and consider the complex ray $z(r) = (cosh r) u + i (sinh r) v in S_bb(C)$ (satisfying $z(r)^top z(r) = cosh^2 r - sinh^2 r = 1$). Since $|exp(chevron.l eta_k, z(r) chevron.r)| = exp((cosh r) chevron.l eta_k, u chevron.r)$, letting $k^* = arg max_k chevron.l eta_k, u chevron.r$, dividing $F(z(r)) = 0$ by $exp(chevron.l eta_(k^*), z(r) chevron.r)$, and sending $r -> +infinity$ forces $a_(k^*) c_d (kappa_(k^*)) = 0$, hence $a_(k^*) = 0$.

  - *Student-$t$ family (\#88):* Each $t$-density $f_k (x) = C_k Q_k (x)^(-(nu_k + d)/2)$ with $Q_k (x) := nu_k + (x - mu_k)^top Sigma_k^(-1) (x - mu_k)$ extends analytically off the real-symmetric complex quadric hypersurface $cal(Q)_k := {z in bb(C)^d : Q_k (z) = 0}$. Distinct pairs $(mu_k, nu_k^(-1) Sigma_k^(-1))$ yield distinct irreducible quadrics whose pairwise intersections have complex codimension at least $2$. Select a component $k^*$ with maximal exponent $nu_(k^*) + d$ among any group sharing $cal(Q)_(k^*)$ (if a quadric is shared, $(mu_k, nu_k^(-1) Sigma_k^(-1))$ coincide so distinct parameters require strictly distinct $nu_k$), and let $z -> z_0 in cal(Q)_(k^*) backslash union.big_(cal(Q)_j != cal(Q)_(k^*)) cal(Q)_j$. Then $Q_(k^*) (z)^(-(nu_(k^*) + d)/2)$ diverges strictly faster than every other term in $sum_k a_k f_k (z) = 0$, forcing $a_(k^*) = 0$. (For canonical SNMM \#98, generic permutation identifiability up to reflection signs follows from @lin2009skewnormal.)

  *Part 2 (Singularity Dichotomy):*
  - *Unboundedness (`V**`, GMM, TMM, SNMM, vMF):* Fix $pi_k = 1 / K$ for all $k$. For GMM (\#81) and every `V**` model (`VII`, `VEI`, `VVI`, `VEE`, `VVE`, `VEV`, `VVV`), place $mu_1 = x_1$, set $Sigma_1 = lambda_1 I_d$ with $lambda_1 -> 0^+$, and fix $mu_2 = dots = mu_K = x_2$, $Sigma_2 = dots = Sigma_K = I_d$ (feasible in all `V**` families since $A_k = I_d, D_k = I_d$ satisfies both `E` and `V` constraints on shape and orientation). Then $log p(x_1 | Theta) >= log(K^(-1) (2 pi lambda_1)^(-d/2)) = (d / 2) log(lambda_1^(-1)) + C_1$, while every $x_i$ ($i >= 2$) contributes at least $log(K^(-1) cal(N)(x_i; x_2, I_d)) > -infinity$, sending $ell(Theta) -> +infinity$. Identical localization at $x_1$ yields $ell(Theta) >= (d/2) log(lambda_1^(-1)) + C$ for TMM (\#88, fixed $nu_1$) and SNMM (\#98, setting $xi_1 = x_1, lambda_1 = bold(0)$). For vMF (\#89, $d >= 2$), setting $mu_1 = x_1 in bb(S)^(d-1)$, $kappa_1 -> +infinity$, and $kappa_2 = 0$ gives
    $
      log f_"vMF" (x_1; x_1, kappa_1) = log c_d (kappa_1) + kappa_1 \
      = ((d - 1) / 2) log kappa_1 + O(1) -> +infinity.
    $

  - *Global Boundedness of All Seven `E**` Models:* In every `E**` model (`EII`, `EEI`, `EVI`, `EEE`, `EVE`, `EEV`, `EVV`), all components share a common volume parameter $det(Sigma_k) = lambda^d$ ($lambda > 0$). Hence every Gaussian density obeys $cal(N)(x_i; mu_k, Sigma_k) <= (2 pi lambda)^(-d/2)$, yielding $log p(x_i | Theta) <= -(d / 2) log(2 pi lambda)$ for every $i = 1, dots, n$. Moreover, since $product_(j=1)^d sigma_(k, j)^2 = lambda^d$, each covariance satisfies $sigma_(k, "min")^2 <= lambda$. Let $v_k in bb(S)^(d-1)$ be a unit eigenvector of $Sigma_k$ associated with $sigma_(k, "min")^2$; then $Sigma_k^(-1) succ.eq sigma_(k, "min")^(-2) v_k v_k^top succ.eq lambda^(-1) v_k v_k^top$, so for every $i, k$:
    $
      (x_i - mu_k)^top Sigma_k^(-1) (x_i - mu_k) >= lambda^(-1) (chevron.l v_k, x_i chevron.r - c_k)^2,
    $
    where $c_k := chevron.l v_k, mu_k chevron.r in bb(R)$. Let $R := max_i norm(x_i)_2$. Because $chevron.l v_k, x_i chevron.r in [-R, R]$, clamping $tilde(c)_k := "clip"(c_k, -2 R, 2 R)$ satisfies $|chevron.l v_k, x_i chevron.r - c_k| >= |chevron.l v_k, x_i chevron.r - tilde(c)_k|$ for every $i$. Define
    $
      g(v_1, c_1, dots, v_K, c_K) := max_(1 <= i <= n) min_(1 <= k <= K) |chevron.l v_k, x_i chevron.r - c_k|.
    $
    On the compact domain $(bb(S)^(d-1) times [-2 R, 2 R])^K$, $g$ is continuous and strictly positive because no $K$ affine hyperplanes ${x : chevron.l v_k, x chevron.r = c_k}$ cover $cal(X)$. Thus $delta_0 := min g > 0$ is attained, ensuring that for *every* `E**` parameter configuration $Theta$, there exists $i^* in {1, dots, n}$ with $min_k (x_(i^*) - mu_k)^top Sigma_k^(-1) (x_(i^*) - mu_k) >= delta_0^2 / lambda$. Summing upper bounds across $i != i^*$ and the tighter bound at $i^*$ gives
    $
      ell_("E**") (Theta) <= - (n d) / 2 log(2 pi lambda) - (delta_0^2) / (2 lambda).
    $
    Maximizing the right-hand side over $lambda > 0$ attains its unique global maximum at $lambda^* = delta_0^2 / (n d)$, yielding the explicit finite ceiling $((n d) / 2) log((n d) / (2 pi e delta_0^2)) < infinity$. (For `EII`, $Sigma_k = lambda I_d$ replaces hyperplanes by $K$ points ${mu_k}$, so $n > K$ distinct observations already suffice.)
]

Theorem 1 delivers a decisive practical lesson: unconstrained GMM (\#81) and all seven `V**` models of MCLUST/Mixmod have *degenerate global MLEs at $+infinity$*, succeeding only because EM converges to interior local maxima or employs Bayesian/spectral floors (as in DPGMM \#82 or eigenvalue-ratio constraints). By contrast, equal-volume models (`EII`--`EVV`) possess globally well-posed, bounded likelihood surfaces without tuning any variance floor!

#subsec("B", "Covariance Parsimony Compression & EM vs. CEM vs. SEM Trajectories")

Even away from singularities, unconstrained GMM (\#81) estimates $p_"VVV" = (K - 1) + K d + K d (d + 1) / 2$ free parameters, demanding $n_k >= d + 1$ effective observations per component merely to invert $hat(Sigma)_k$. MCLUST (\#93) and Mixmod (\#94) systematically traverse the 14-model parsimony lattice, compressing covariance complexity from $Theta(K d^2)$ down to $Theta(1)$ (isotropic `EII`), $Theta(K + d)$ (diagonal `EEI`/`EVI`), or $Theta(d^2)$ (homoscedastic `EEE`), while Mixmod further expands inference beyond soft EM to Classification EM (CEM) and Stochastic EM (SEM) @biernacki2006mixmod.

#proposition(title: "Parsimony Parameter Compression & Exact CEM Variance-Truncation Bias", origin: "Original Unified Proof")[
  1. *(Asymptotic Parameter Compression Ratio):* Across the Celeux--Govaert covariance lattice for $K >= 2$ clusters in $bb(R)^d$, the free parameter counts satisfy
     $
       p_"EII" &= K d + K, quad p_"VVI" = 2 K d + K - 1, \
       p_"EEE" &= K d + K - 1 + d(d + 1) / 2, \
       p_"VVV" &= K d + K - 1 + K d(d + 1) / 2,
     $
     giving exact asymptotic compression $p_"VVV" / p_"EII" = (d + 3) / 2 - O(d^(-1))$ and $p_"VVV" / p_"EEE" -> K$ as $d -> infinity$.
  2. *(Exact Complete-Data Entropy Sandwich):* For any parameter $Theta$, let $r_(i k) = p(z_i = k | x_i, Theta)$ be the E-step responsibilities, $H(Theta) = -sum_(i=1)^n sum_(k=1)^K r_(i k) log r_(i k)$ the Shannon classification entropy, and $hat(z)_i = arg max_k r_(i k)$ the MAP hard partition. Then the expected complete-data objective $Q(Theta | Theta)$, the CEM complete-data log-likelihood $ell_C (Theta, hat(z))$, and the observed marginal log-likelihood $ell_"EM" (Theta)$ satisfy the exact sandwich
     $
       ell_"EM" (Theta) - H(Theta) = Q(Theta | Theta) <= ell_C (Theta, hat(z)) <= ell_"EM" (Theta),
     $
     with both inequalities strict whenever any observation has non-degenerate posterior uncertainty ($max_k r_(i k) < 1$).
  3. *(Exact Asymptotic Variance-Underestimation & Outward Mean Bias of CEM):* On a balanced symmetric 1D Gaussian mixture $(1/2) cal(N)(-mu, sigma^2) + (1/2) cal(N)(mu, sigma^2)$ with signal-to-noise ratio $s := mu / sigma >= 0$, CEM (Mixmod \#94) converges as $n -> infinity$ to *strictly inconsistent* estimates satisfying
     $
       hat(mu)_"CEM" &= sigma [ s (2 Phi(s) - 1) + 2 phi(s) ] > mu, \
       (hat(sigma)_"CEM"^2) / (sigma^2) &= 1 + s^2 - [ s (2 Phi(s) - 1) + 2 phi(s) ]^2 < 1
     $
     for every finite $s >= 0$ (where $phi, Phi$ are the standard normal PDF and CDF), suffering a $63.7\%$ variance collapse at $s = 0$ and a $36.1\%$ collapse even at moderate separation $s = 1$, whereas EM (\#81/\#93/\#94) is asymptotically consistent ($hat(mu)_"EM" -> mu$, $hat(sigma)_"EM"^2 -> sigma^2$).
]

#proof[
  *Part 1:* Every mixture requires $K - 1$ mixing weights and $K d$ mean coordinates. Covariance parameters add $1$ (`EII`), $K d$ (`VVI`), $d(d+1)/2$ (`EEE`), or $K d(d+1)/2$ (`VVV`), yielding the stated exact counts and ratios immediately.

  *Part 2:* By Bayes' rule, $log(pi_k f_k (x_i)) = log p(x_i | Theta) + log r_(i k)$. Summing over $i$ at $k = hat(z)_i$ yields
  $
    ell_C (Theta, hat(z)) = ell_"EM" (Theta) + sum_(i=1)^n log(max_k r_(i k)) <= ell_"EM" (Theta).
  $
  Conversely, for every $i$ and $k$, $log r_(i k) <= log(max_j r_(i j))$; averaging against $r_(i k)$ gives $sum_k r_(i k) log r_(i k) <= log(max_j r_(i j))$, whence
  $
    Q(Theta | Theta) = ell_"EM" (Theta) - H(Theta) <= ell_C (Theta, hat(z)).
  $

  *Part 3:* By mirror symmetry of $p(x) = (1/2) sigma^(-1) phi((x - mu)/sigma) + (1/2) sigma^(-1) phi((x + mu)/sigma)$, the MAP boundary is $x = 0$, assigning the positive half-line ${X > 0}$ to the right component. Conditioned on $X > 0$, the density is $sigma^(-1)[phi((x - mu)/sigma) + phi((x + mu)/sigma)]$ ($x > 0$), which is *identically* the folded-normal distribution of $|Z|$ for $Z tilde.op cal(N)(mu, sigma^2)$. Consequently, CEM's large-sample M-step updates compute the exact mean and variance of $|Z|$:
  $
    hat(mu)_"CEM" &= bb(E)[abs(Z)] = integral_(-s)^infinity (mu + sigma u) phi(u) dif u \
    &quad + integral_s^infinity (-mu + sigma u) phi(u) dif u \
    &= mu (2 Phi(s) - 1) + 2 sigma phi(s).
  $
  Subtracting $mu$ yields
  $
    hat(mu)_"CEM" - mu = 2 sigma [ phi(s) - s (1 - Phi(s)) ] = 2 sigma integral_s^infinity (u - s) phi(u) dif u > 0
  $
  strictly for every $s >= 0$. Finally, since $bb(E)[abs(Z)^2] = bb(E)[Z^2] = mu^2 + sigma^2$,
  $
    hat(sigma)_"CEM"^2 = "Var"(|Z|) = mu^2 + sigma^2 - hat(mu)_"CEM"^2 = sigma^2 - (hat(mu)_"CEM"^2 - mu^2) < sigma^2.
  $
  At $s = 0$, $hat(sigma)_"CEM"^2 / sigma^2 = 1 - 2 / pi approx 0.3634$; at $s = 1$, $hat(mu)_"CEM" / sigma approx 1.1666$ and $hat(sigma)_"CEM"^2 / sigma^2 approx 0.639$.
]

Proposition 1 clarifies Mixmod's (@biernacki2006mixmod) three-engine strategy: *EM* delivers consistent density estimation but can stall at unstable saddles and slows down on large $n$; *CEM* accelerates finite termination (generalizing Lloyd's $k$-Means, which equals CEM under `EII` with fixed equal weights) at the price of systematic outward mean repulsion and variance shrinkage whenever clusters overlap; and *SEM* injects stochastic Gibbs-like assignment draws $z_i^((t)) tilde.op "Cat"(r_i^((t)))$ to escape saddle points and spurious local traps before final EM polishing.

#subsec("C", "Heavy-Tail Influence Boundedness, Skewness & Directional Curvature")

When data depart from symmetric Euclidean Gaussians—via gross outliers, skewed tails, or unit-norm directional constraints—GMM (\#81) and MCLUST (\#93) fail structurally in three distinct ways addressed respectively by TMM (\#88), SNMM (\#98), and vMF (\#89).

#theorem(title: "Redescending TMM Influence, Skewness Absorption, and vMF Extrinsic Conditioning", origin: "Original Unified Proof")[
  1. *(Redescending Location & Bounded Scatter Influence of TMM \#88 vs. Unbounded GMM \#81):* Consider a single contaminating outlier $x_0$ with component Mahalanobis radius $R := sqrt((x_0 - mu_k)^top Sigma_k^(-1) (x_0 - mu_k)) -> infinity$ assigned non-negligible responsibility $r_(0 k) > 0$. In GMM (\#81), the whitened M-step location and scatter influence terms diverge linearly ($Theta(R)$) and quadratically ($Theta(R^2)$), yielding finite-sample breakdown $1 / n$. In TMM (\#88) with finite degrees of freedom $nu_k > 0$ and Gamma scale weight $u_k (x_0) = (nu_k + d) / (nu_k + R^2)$, the whitened location and scatter influences satisfy uniformly over all $R >= 0$:
     $
       &norm(Sigma_k^(-1/2) u_k (x_0) (x_0 - mu_k))_2 = ((nu_k + d) R) / (nu_k + R^2) <= (nu_k + d) / (2 sqrt(nu_k)), \
       &norm(Sigma_k^(-1/2) u_k (x_0) (x_0 - mu_k)(x_0 - mu_k)^top Sigma_k^(-1/2))_2 \
       &quad = ((nu_k + d) R^2) / (nu_k + R^2) < nu_k + d,
     $
     with the location influence *strictly redescending to $0$* as $R -> infinity$ (peaking at $R^* = sqrt(nu_k)$).
  2. *(Asymmetry Absorption vs. Outlier Vulnerability of SNMM \#98):* In the canonical multivariate Skew-Normal mixture $X | (Z = k, tau) tilde.op cal(N)(xi_k + lambda_k tau, Sigma_k)$ ($tau tilde.op "TN"_+ (0, 1)$) @lin2009skewnormal, any single skewed unimodal component ($lambda_k != bold(0)$) is represented exactly with $K = 1$ (avoiding symmetric GMM's asymptotic over-splitting into multiple Gaussians), yet orthogonal to the skewness vector $lambda_k$ the M-step location and scatter scores remain linearly ($Theta(R)$) and quadratically ($Theta(R^2)$) unbounded.
  3. *(Exact Extrinsic Covariance Ill-Conditioning of GMM on $bb(S)^(d-1)$ vs. Intrinsic vMF \#89):* Let $X tilde.op "vMF"(mu, kappa)$ on $bb(S)^(d-1)$ ($d >= 2$) with concentration $kappa > 0$ and Bessel ratio $A_d (kappa) = (I_(d/2) (kappa)) / (I_(d/2 - 1) (kappa))$. The *exact* extrinsic Euclidean covariance matrix of $X$ is
     $
       Sigma_"ext" = "Cov"(X) = (A_d (kappa)) / kappa (I_d - mu mu^top) + A'_d (kappa) mu mu^top,
     $
     whose exact condition number diverges linearly with concentration $kappa$:
     $
       "cond"(Sigma_"ext") = (A_d (kappa)) / (kappa A'_d (kappa)) = (2 kappa) / (d - 1) + O(1) quad (kappa -> infinity).
     $
     Consequently, fitting extrinsic Euclidean GMM (\#81) to concentrated directional clusters forces near-singular radial collapse ($lambda_parallel = Theta(kappa^(-2))$ vs. $lambda_perp = Theta(kappa^(-1))$) while paying $Theta(K d^2)$ covariance parameters, whereas vMF (\#89) stays intrinsically well-conditioned with only $K d$ parameters and compact bounded score $norm(x_0)_2 = 1$.
]

#proof[
  *Part 1:* In the ECM M-step of TMM (\#88, @peel2000robust), $hat(mu)_k$ and $hat(Sigma)_k$ solve the weighted score equations
  $
    sum_(i=1)^n r_(i k) u_(i k) (x_i - mu_k) = bold(0), \
    Sigma_k = (sum_(i=1)^n r_(i k) u_(i k) (x_i - mu_k)(x_i - mu_k)^top) / (sum_(i=1)^n r_(i k)),
  $
  where $u_(i k) = (nu_k + d) / (nu_k + delta_k (x_i))$ and $delta_k (x_i) = (x_i - mu_k)^top Sigma_k^(-1) (x_i - mu_k)$. Setting $y_0 = Sigma_k^(-1/2) (x_0 - mu_k)$ gives $norm(y_0)_2 = R$ and $delta_k (x_0) = R^2$. In GMM ($nu_k -> infinity$), $u_k (x_0) -> 1$, so $norm(y_0)_2 = R$ and $norm(y_0 y_0^top)_2 = R^2$ diverge without bound. For finite $nu_k > 0$:
  - Location norm: $g(R) = (nu_k + d) R / (nu_k + R^2)$. By AM--GM, $nu_k + R^2 >= 2 sqrt(nu_k) R$ with equality uniquely at $R^* = sqrt(nu_k)$, yielding $g(R) <= (nu_k + d) / (2 sqrt(nu_k))$ and $lim_(R -> infinity) g(R) = 0$.
  - Scatter operator norm: since $y_0 y_0^top$ is rank-1 PSD with spectral norm $R^2$, $norm(u_k (x_0) y_0 y_0^top)_2 = (nu_k + d) R^2 / (nu_k + R^2) < nu_k + d$ strictly for all $R >= 0$.

  *Part 2:* In SNMM (\#98, @lin2009skewnormal), conditioning on $x_0$ gives $tau | (x_0, Z = k) tilde.op "TN"_+ (mu_(tau, 0), sigma_tau^2)$ with
  $
    sigma_tau^2 = (1 + lambda_k^top Sigma_k^(-1) lambda_k)^(-1), \
    mu_(tau, 0) = sigma_tau^2 lambda_k^top Sigma_k^(-1) (x_0 - xi_k).
  $
  Along any perturbation $x_0 = xi_k + R w$ with $w in.not "span"(lambda_k)$ satisfying $lambda_k^top Sigma_k^(-1) w = 0$, we have $mu_(tau, 0) = 0$ identically, so $bb(E)[tau | x_0]$ and $bb(E)[tau^2 | x_0]$ remain fixed $O(1)$ constants independent of $R$. Consequently, the CM-step residual $x_0 - xi_k - lambda_k bb(E)[tau | x_0] = R w + O(1)$ diverges linearly in $R$ and its outer product diverges quadratically in $R^2$, matching GMM's unbounded orthogonal sensitivity.

  *Part 3:* Write the vMF density in canonical exponential-family form $f(x; eta) = exp(chevron.l eta, x chevron.r - A(eta))$ with $eta = kappa mu$ ($kappa = norm(eta)_2$) and log-partition $A(eta) = -log c_d (norm(eta)_2)$. Because $A'(kappa) = A_d (kappa)$ @banerjee2005vmf, differentiating twice with respect to $eta$ yields the exact cumulant identity:
  $
    bb(E)[X] &= nabla_eta A(eta) = A_d (kappa) mu, \
    Sigma_"ext" &= nabla_eta^2 A(eta) = (A_d (kappa)) / kappa (I_d - mu mu^top) + A'_d (kappa) mu mu^top.
  $
  The eigenvalues of $Sigma_"ext"$ are therefore $lambda_perp = (A_d (kappa)) / kappa$ (multiplicity $d - 1$, tangent to $bb(S)^(d-1)$) and $lambda_parallel = A'_d (kappa)$ (multiplicity $1$, radial along $mu$). From the exact Riccati identity for modified Bessel ratios $A'_d (kappa) = 1 - A_d (kappa)^2 - ((d - 1) / kappa) A_d (kappa)$ combined with the large-$kappa$ expansion $A_d (kappa) = 1 - (d - 1) / (2 kappa) + (d - 1)(d - 3) / (8 kappa^2) + O(kappa^(-3))$, we obtain
  $
    lambda_perp &= 1 / kappa - (d - 1) / (2 kappa^2) + O(kappa^(-3)), \
    lambda_parallel &= 1 - (1 - (d - 1) / kappa + ((d - 1)^2 + (d - 1)(d - 3)) / (4 kappa^2)) \
    &quad - (d - 1) / kappa + ((d - 1)^2) / (2 kappa^2) + O(kappa^(-3)) \
    &= (d - 1) / (2 kappa^2) + O(kappa^(-3)).
  $
  Dividing $lambda_perp$ by $lambda_parallel$ gives $"cond"(Sigma_"ext") = (2 kappa) / (d - 1) + O(1)$, diverging linearly with $kappa$.
]

Theorem 2 draws a sharp operational triangle across continuous parametric extensions: *TMM (\#88)* defends against heavy tails and gross outliers via redescending Gamma weights (recovering GMM as $nu_k -> infinity$ and heavy Cauchy mixtures at $nu_k = 1$), *SNMM (\#98)* defends against asymmetric skewness without multi-Gaussian over-splitting (recovering GMM as $lambda_k -> bold(0)$), and *vMF (\#89)* defends against extrinsic spherical curvature and $O(K d^2)$ covariance explosion on normalized embeddings (recovering Spherical $k$-Means, Alg. 12, as $kappa -> infinity$).

#wide_figure(
  caption: [*Table I:* Head-to-Head Structural Comparison of Finite Continuous Parametric Mixture Models (Section I: Algs. 81, 88, 89, 93, 94, 98).],
  table(
    columns: (1.05fr, 1.0fr, 0.92fr, 0.85fr, 0.95fr, 1.23fr),
    align: (left, left, center, center, left, left),
    [Algorithm], [Emission / Structure], [Free Params], [Unreg. MLE], [Outlier / Skew / Sph.], [Primary Failure Mode],
    [GMM (\#81) @dempster1977em], [Full Gaussian $cal(N)(mu_k, Sigma_k)$], [$Theta(K d^2)$], [Unbounded ($+infinity$)], [Unbnd. / Splits / Ill-cond.], [Singularities; $n_k <= d$ collapse],
    [MCLUST (\#93) @fraley2002mclust], [14 models `EII`--`VVV` + BIC], [$Theta(K d)$--$Theta(K d^2)$], [`E**` Bnd. / `V**` Unbnd.], [Unbnd. / Splits / Ill-cond.], [Over-splits skewed/heavy-tailed modes],
    [Mixmod (\#94) @biernacki2006mixmod], [14 Gauss. + Mult.; EM/CEM/SEM], [$Theta(K d)$--$Theta(K d^2)$], [`E**` Bnd. / `V**` Unbnd.], [Unbnd. / Splits / Ill-cond.], [CEM variance shrinkage (Prop. 1.3)],
    [TMM (\#88) @peel2000robust], [Student-$t$ $(mu_k, Sigma_k, nu_k)$], [$Theta(K d^2) + K$], [Unbounded ($+infinity$)], [*Redesc.* / Splits / Ill-cond.], [Slow $nu_k$ CM-step; symmetric tails],
    [SNMM (\#98) @lin2009skewnormal], [Skew-Normal $(xi_k, Sigma_k, lambda_k)$], [$Theta(K d^2) + K d$], [Unbounded ($+infinity$)], [Unbnd. $perp lambda_k$ / *Exact* / No], [Outlier pull $perp lambda_k$; heavy ECM],
    [vMF (\#89) @banerjee2005vmf], [Directional $c_d (kappa_k) e^(kappa_k mu_k^top x)$], [$K d + K - 1$], [Unbounded ($kappa -> infinity$)], [*Compact* / Sym. / *Exact*], [Ignores vector norms; Bessel overflow],
  ),
)

// =============================================================================
// SECTION II
// =============================================================================
#sec("II", "Bayesian Nonparametric Mixtures & Model Selection")

A central liability of Section I's fixed-$K$ finite mixtures is determining the component cardinality $K$ while simultaneously guarding against likelihood singularities and overfitting. Four algorithms in Category IV directly tackle automatic model-order selection and full Bayesian regularization—Bayesian DPGMM (Alg. 82) @blei2006variational, Infinite GMM (IGMM, Alg. 99) @rasmussen2000igmm, AutoClass (Alg. 87) @cheeseman1996autoclass, and SNOB (Alg. 102) @wallace1968snob—forming a natural four-way contest alongside MCLUST's BIC (Alg. 93) @schwarz1978bic and Mixmod's ICL (Alg. 94) @biernacki2006mixmod.

#subsec("A", "Variational Stick-Breaking (DPGMM) vs. Collapsed Gibbs CRP (IGMM)")

Both DPGMM (\#82) and IGMM (\#99) place a Dirichlet Process prior $G tilde.op "DP"(alpha, G_0)$ on the mixing measure over Gaussian component parameters, eliminating manual choice of $K$ and imposing conjugate Normal--Inverse-Wishart (NIW) regularization $G_0 = "NIW"(mu_0, kappa_0, Psi_0, nu_0)$. Yet their inference representations—*truncated stick-breaking mean-field variational Bayes (CAVI)* in DPGMM (\#82) versus *untruncated collapsed Gibbs sampling over the Chinese Restaurant Process (CRP)* in IGMM (\#99)—induce three sharp mathematical separations.

#theorem(title: "DPGMM Truncation Coupling, CRP Growth, and Mean-Field Variance Collapse", origin: "Original Unified Proof")[
  Let $G tilde.op "DP"(alpha, G_0)$ with concentration $alpha > 0$ and conjugate NIW base measure $G_0 = "NIW"(mu_0, kappa_0, Psi_0, nu_0)$ ($Psi_0 succ 0$, $nu_0 > d - 1$).
  1. *(Sharp Coupling Bound for DPGMM Stick-Breaking Truncation):* Let $p_infinity (x_(1:n))$ denote the exact marginal density of $n$ observations under the infinite GEM($alpha$) stick-breaking mixture and $p_T (x_(1:n))$ the $T$-truncated DPGMM (\#82) marginal setting $V_T = 1$ ($pi_T = product_(j=1)^(T-1) (1 - V_j)$). Their total-variation distance obeys the exact finite-sample bound
     $
       "TV"(p_T, p_infinity) <= 1 - bb(E) [ (1 - R_(T-1))^n ] <= n ( alpha / (alpha + 1) )^(T - 1),
     $
     where $R_(T-1) = product_(j=1)^(T-1) (1 - V_j)$, guaranteeing $"TV"(p_T, p_infinity) <= epsilon$ whenever $T >= 1 + (alpha + 1) ln(n / epsilon)$.
  2. *(Exact CRP Cluster Growth of IGMM \#99 & Shared Singularity Immunity):* Under IGMM's untruncated CRP prior, the number of occupied clusters $K_n$ obeys
     $
       bb(E)[K_n] &= sum_(i=1)^n alpha / (alpha + i - 1) = alpha ln(1 + n / alpha) + O(1), \
       "Var"(K_n) &= sum_(i=1)^n (alpha (i - 1)) / (alpha + i - 1)^2 <= bb(E)[K_n].
     $
     Furthermore, both DPGMM (\#82) and IGMM (\#99) globally eliminate GMM's singularity pathology: for every component $k$, the variational/posterior precision matrix obeys the uniform Loewner upper bound
     $
       bb(E)_q [Sigma_k^(-1)] = (nu_0 + N_k) Psi_k^(-1) prec.eq (nu_0 + n) Psi_0^(-1),
     $
     preventing any covariance eigenvalue from collapsing below $(nu_0 + n)^(-1) lambda_"min" (Psi_0) > 0$.
  3. *(Exact Mean-Field Variance Underestimation of DPGMM CAVI vs. Exact IGMM Posterior):* On a correlated bivariate Gaussian target posterior $p(theta_1, theta_2) = cal(N)(bold(0), Sigma)$ with marginal variances $sigma_1^2, sigma_2^2$ and correlation $rho in (-1, 1)$, the unique factorized variational minimizer $q^* (theta_1, theta_2) = q_1^* (theta_1) q_2^* (theta_2) = arg min D_"KL" (q || p)$ satisfies
     $
       "Var"_(q_j^*) (theta_j) = sigma_j^2 (1 - rho^2) quad (j = 1, 2),
     $
     underestimating every marginal posterior variance by the exact factor $1 - rho^2 -> 0$ as $|rho| -> 1$, whereas IGMM's ergodic collapsed Gibbs chain converges to the exact correlated posterior $p$ without factorized bias.
]

#proof[
  *Part 1:* Couple $p_infinity$ and $p_T$ on a shared probability space using identical sticks $(V_1, dots, V_(T-1))$ with $V_j tilde.op "Beta"(1, alpha)$ i.i.d., identical atoms $(theta_1, dots, theta_(T-1)) tilde.op G_0$, and identical allocation uniforms $U_1, dots, U_n tilde.op "Unif"(0, 1)$, where $z_i = k$ iff $U_i in (sum_(j < k) pi_j, sum_(j <= k) pi_j]$. Let $R_(T-1) = 1 - sum_(j=1)^(T-1) pi_j = product_(j=1)^(T-1) (1 - V_j) in (0, 1)$. On the event $E = {max_(1 <= i <= n) U_i <= 1 - R_(T-1)}$, every observation maps identically to $z_i in {1, dots, T - 1}$ in both models, producing identical $x_(1:n)$. By Le Cam's coupling inequality and Bernoulli's inequality $(1 - R_(T-1))^n >= 1 - n R_(T-1)$,
  $
    "TV"(p_T, p_infinity) &<= bb(P)(E^c) = 1 - bb(E)[(1 - R_(T-1))^n] \
    &<= n bb(E)[R_(T-1)] = n product_(j=1)^(T-1) bb(E)[1 - V_j] = n ( alpha / (alpha + 1) )^(T - 1),
  $
  using $bb(E)[1 - V_j] = alpha / (alpha + 1)$. Since $-ln(alpha / (alpha + 1)) = ln(1 + alpha^(-1)) >= (alpha + 1)^(-1)$, choosing $T >= 1 + (alpha + 1) ln(n / epsilon)$ forces $n (alpha / (alpha + 1))^(T-1) <= epsilon$.

  *Part 2:* In the CRP, customer $i$ opens a new table with independent Bernoulli indicator $B_i tilde.op "Bern"(alpha / (alpha + i - 1))$. Summing means and variances of independent $B_i$ yields $bb(E)[K_n] = sum_(i=1)^n alpha / (alpha + i - 1) = alpha ln(1 + n / alpha) + O(1)$ and $"Var"(K_n) = sum_(i=1)^n p_i (1 - p_i) <= sum_i p_i = bb(E)[K_n]$. For singularity immunity, every NIW posterior scale matrix has the PSD sum structure $Psi_k = Psi_0 + S_k + (kappa_0 N_k) / (kappa_0 + N_k) (overline(x)_k - mu_0)(overline(x)_k - mu_0)^top succ.eq Psi_0$ (where $S_k succ.eq bold(0)$ is the within-cluster scatter). By Loewner monotonicity of matrix inversion, $Psi_k^(-1) prec.eq Psi_0^(-1)$, and since $N_k <= n$, $bb(E)_q [Sigma_k^(-1)] = (nu_0 + N_k) Psi_k^(-1) prec.eq (nu_0 + n) Psi_0^(-1)$.

  *Part 3:* Write the precision matrix of $p$:
  $
    Sigma^(-1) = 1 / (1 - rho^2) mat(sigma_1^(-2), -rho (sigma_1 sigma_2)^(-1); -rho (sigma_1 sigma_2)^(-1), sigma_2^(-2)).
  $
  For any factorized $q = q_1 q_2$, the cross-term vanishes at the stationary means $bb(E)_(q_j) [theta_j] = 0$, and Gaussianity of each coordinate maximizer follows from the quadratic exponent of $log p$. Letting $v_j = "Var"_(q_j) (theta_j) > 0$,
  $
    2 D_"KL" (q || p) = v_1 / (sigma_1^2 (1 - rho^2)) + v_2 / (sigma_2^2 (1 - rho^2)) - ln v_1 - ln v_2 + "const".
  $
  Setting $partial / (partial v_j) D_"KL" (q || p) = 0$ uniquely yields $v_j^* = sigma_j^2 (1 - rho^2)$.
]

Theorem 3 exposes the exact engineering duality between DPGMM (\#82) and IGMM (\#99): *DPGMM* delivers deterministic monotonic ELBO ascent in $O(T_"iter" n T d^2)$ time (scaling gracefully to $n >= 10^5$ via stochastic variational inference) with exponentially accurate truncation once $T > 1 + (alpha + 1) ln(n / epsilon)$, at the cost of local-mode trapping and factorized variance shrinkage; *IGMM* preserves exact posterior uncertainty and an unbounded cluster support via collapsed Student-$t$ predictive conditionals, but pays slow MCMC autocorrelation (exponentially slow label-split/merge mixing unless augmented with split--merge proposals) and $O(T_"MCMC" n K_"eff" d^2)$ sequential runtime.

#subsec("B", "Unified Asymptotic Hierarchy: BIC, ICL, Cheeseman–Stutz & SNOB MML")

Rather than placing a nonparametric infinite prior on $K$, MCLUST (\#93), Mixmod (\#94), AutoClass (\#87), and SNOB (\#102) score candidate finite models across $K = 1, 2, dots, K_"max"$ using four distinct penalized-likelihood criteria. Despite arising from different schools—Laplace asymptotics (BIC), complete-data classification entropy (ICL), approximate Bayesian evidence (AutoClass Cheeseman--Stutz), and algorithmic information theory (SNOB Minimum Message Length)—the next unified theorem places all four on a single exact asymptotic scale.

#theorem(title: "Unified O(1) and O(n) Hierarchy of BIC, MML, Cheeseman–Stutz, and ICL", origin: "Original Unified Proof")[
  Let $cal(M)_K$ be a regular identifiable $K$-component mixture with $p_K$ free parameters, MLE $hat(Theta)_K$, smooth positive prior density $h(Theta)$, observed per-sample Fisher information $I_"obs" (Theta) = -n^(-1) nabla^2 ell(Theta)$, complete-data Fisher information $I_"comp" (Theta) = -n^(-1) bb(E)[nabla^2 ell_"comp" (Theta) | cal(X)]$, and posterior classification entropy $H(hat(Theta)_K) = -sum_(i=1)^n sum_(k=1)^K r_(i k) log r_(i k)$. Then, expressed as negative log-evidence code lengths (minimized over $K$):
  1. *(Shared BIC Leading Order):* True Bayesian marginal negative log-evidence $-log p(cal(X) | K)$, MCLUST's $- (1/2) "BIC"(K)$ (\#93), SNOB's two-part MML message length $cal(L)_"MML" (hat(Theta)_K, K)$ (\#102), and AutoClass's negative Cheeseman--Stutz score $-log p_"CS" (cal(X) | K)$ (\#87) all share the identical leading expansion
     $
       -ell(hat(Theta)_K) + (p_K / 2) log n + O(1).
     $
  2. *(Exact Fisher-Metric Refinement of SNOB MML \#102):* SNOB's Wallace--Freeman code length refines BIC at $O(1)$ by retaining the reparameterization-invariant Fisher-volume and lattice-quantization terms @wallace1968snob:
     $
       cal(L)_"MML" (hat(Theta)_K, K) &= -ell(hat(Theta)_K) + (p_K / 2) log n + 1 / 2 log det I_"obs" (hat(Theta)_K) \
       &quad - log h(hat(Theta)_K) + (p_K / 2) (1 + log kappa_(p_K)),
     $
     where $kappa_p$ is the $p$-dimensional optimal quantizer constant ($kappa_1 = 1 / 12$, $kappa_p -> (2 pi e)^(-1)$). Unlike raw BIC, minimizing $cal(L)_"MML"$ jointly over $Theta$ yields the exact finite-sample MML variance estimate $hat(sigma)_"MML"^2 = (N_k - 1)^(-1) S_k$ (for $N_k > 1$) and penalizes singular collapses via $det I_"obs" -> +infinity$.
  3. *(Exact Missing-Information Bias of AutoClass Cheeseman–Stutz \#87):* Defining $log p_"CS" (cal(X) | K) := log p_"comp" (cal(X), hat(Z) | K) + ell(hat(Theta)_K) - ell_"comp" (hat(Theta)_K, hat(Z))$ @cheeseman1996autoclass, its exact asymptotic error relative to true log-evidence is strictly non-positive:
     $
       log p_"CS" (cal(X) | K) - log p(cal(X) | K) &= 1 / 2 log( (det I_"obs" (hat(Theta)_K)) / (det I_"comp" (hat(Theta)_K)) ) \
       &quad + O(n^(-1/2)) <= 0,
     $
     with equality iff clusters have zero overlap ($I_"miss" = I_"comp" - I_"obs" = bold(0)$).
  4. *(O(n) Entropy Separation Gap of Mixmod ICL \#94):* Mixmod's negative half-ICL score satisfies
     $
       - 1 / 2 "ICL"(K) = -ell(hat(Theta)_K) + (p_K / 2) log n + H(hat(Theta)_K).
     $
     Consequently, on any non-Gaussian unimodal cluster (e.g., skewed or heavy-tailed) whose optimal 2-Gaussian approximation gains $n Delta ell > 0$ nats of density fit while inducing $n overline(h) > 0$ nats of posterior overlap entropy, BIC (\#93), MML (\#102), and Cheeseman--Stutz (\#87) *all* over-split to $K >= 2$ asymptotically as $n -> infinity$, whereas ICL (\#94) retains the true single cluster $K = 1$ whenever $overline(h) > Delta ell$.
]

#proof[
  *Part 1 & Part 2:* Under standard regularity conditions around the interior MLE $hat(Theta)_K$, Laplace's method on $p(cal(X) | K) = integral exp(ell(Theta)) h(Theta) dif Theta$ expands the quadratic form $-nabla^2 ell(hat(Theta)_K) = n I_"obs" (hat(Theta)_K)$, yielding
  $
    -log p(cal(X) | K) &= -ell(hat(Theta)_K) + (p_K / 2) log n + 1 / 2 log det I_"obs" (hat(Theta)_K) \
    &quad - log h(hat(Theta)_K) - (p_K / 2) log(2 pi) + O(n^(-1/2)).
  $
  Dropping $O(1)$ terms recovers $-(1/2) "BIC"(K) = -ell(hat(Theta)_K) + (p_K / 2) log n$ @schwarz1978bic. For SNOB (\#102), the two-part Wallace--Freeman message length encodes $Theta$ inside a quantization cell of volume $V(Theta) = (kappa_(p_K)^(p_K) det cal(I)_n (Theta))^(-1/2)$ (with expected Fisher matrix $cal(I)_n (Theta) = n I_"obs" (Theta)$) plus the second-order Taylor expectation $(p_K / 2)$ of the rounding loss in $-ell(Theta)$ @wallace1968snob:
  $
    cal(L)_"MML" (Theta, K) &= -log(h(Theta) V(Theta)) - ell(Theta) + p_K / 2 \
    &= -ell(Theta) + (p_K / 2) log n + 1 / 2 log det I_"obs" (Theta) \
    &quad - log h(Theta) + (p_K / 2) (1 + log kappa_(p_K)).
  $
  Since $1 + log kappa_p -> -log(2 pi)$, $cal(L)_"MML"$ matches the full Laplace negative log-evidence even at $O(1)$ up to $(p_K / 2) log(2 pi e kappa_(p_K)) = o(p_K)$! For a 1D Gaussian component of size $N_k > 1$ and scatter $S_k > 0$, $cal(I)_n (mu, sigma) = "diag"(N_k / sigma^2, 2 N_k / sigma^2)$, so $(1/2) log det cal(I)_n = -2 log sigma + C$; combining with the scale-invariant prior $-log h(sigma) = +log sigma + C'$ ($h(sigma) prop sigma^(-1)$) contributes $-log sigma$ to $-ell = N_k log sigma + S_k / (2 sigma^2)$, shifting the minimizing variance from $S_k / N_k$ to $S_k / (N_k - 1)$ while finite measurement precision $epsilon_"meas"$ blocks degenerate $sigma -> 0^+$ collapse.

  *Part 3:* Applying Laplace's expansion identically to the complete-data marginal $p_"comp" (cal(X), hat(Z) | K)$ around $hat(Theta)_K$ replaces $I_"obs" (hat(Theta)_K)$ by the complete-data Fisher matrix $I_"comp" (hat(Theta)_K)$:
  $
    &log p_"comp" (cal(X), hat(Z) | K) = ell_"comp" (hat(Theta)_K, hat(Z)) + log h(hat(Theta)_K) \
    &quad + (p_K / 2) log((2 pi) / n) - 1 / 2 log det I_"comp" (hat(Theta)_K) + O(n^(-1/2)).
  $
  Adding $ell(hat(Theta)_K) - ell_"comp" (hat(Theta)_K, hat(Z))$ cancels $ell_"comp"$ and yields $log p_"CS"$, whose difference from the Laplace expansion of $log p(cal(X) | K)$ is precisely
  $
    log p_"CS" (cal(X) | K) - log p(cal(X) | K) = 1 / 2 log( (det I_"obs" (hat(Theta)_K)) / (det I_"comp" (hat(Theta)_K)) ) + O(n^(-1/2)).
  $
  By Louis's missing-information identity, $I_"obs" = I_"comp" - I_"miss"$ where $I_"miss" = n^(-1) "Cov"(nabla ell_"comp" | cal(X)) succ.eq bold(0)$, hence $bold(0) prec I_"obs" prec.eq I_"comp"$, forcing $det I_"obs" <= det I_"comp"$ and $log((det I_"obs") / (det I_"comp")) <= 0$ with equality iff $I_"miss" = bold(0)$.

  *Part 4:* By definition @biernacki2006mixmod, $"ICL"(K) = "BIC"(K) - 2 H(hat(Theta)_K)$. On an i.i.d. sample of size $n$ from a mis-specified unimodal density $f_0$, $ell(hat(Theta)_2) - ell(hat(Theta)_1) = n Delta ell + O_p (sqrt(n))$ with $Delta ell = D_"KL" (f_0 || g_1^*) - D_"KL" (f_0 || g_2^*) > 0$, while the complexity penalty $(p_2 - p_1)/2 log n + O(1)$ is strictly $o(n)$ for BIC, MML, and Cheeseman--Stutz—guaranteeing asymptotic selection of $K >= 2$. In contrast, ICL subtracts $H(hat(Theta)_2) = n overline(h) + O_p (sqrt(n))$, giving net score change $n(Delta ell - overline(h)) - ((p_2 - p_1)/2) log n < 0$ for all sufficiently large $n$ whenever $overline(h) > Delta ell$.
]

#remark(title: "Diagonal-Independence Over-Splitting in AutoClass (#87) and SNOB (#102)")[
  Classical AutoClass (\#87, in its independent-attribute mode) and classical SNOB (\#102) assume conditional independence of features within each cluster ($Sigma_k = "diag"(sigma_(k 1)^2, dots, sigma_(k d)^2)$). While this yields $O(n K d)$ scaling and effortless mixed-type (continuous + categorical + MAR missingness) handling in AutoClass @cheeseman1996autoclass, Theorem 4(4) reveals an acute geometric cost: when $cal(X)$ is drawn from a *single* correlated bivariate Gaussian with correlation $rho != 0$, the KL gap between the true correlated ellipse and any single diagonal Gaussian is $Delta ell = -(1/2) log(1 - rho^2) > 0$, driving both SNOB and independent-attribute AutoClass to tile a single tilted cluster into $K -> infinity$ diagonal beads as $n -> infinity$! Full-covariance DPGMM (\#82), IGMM (\#99), and MCLUST `EEE`/`VVV` (\#93) represent the tilted cluster cleanly with $K = 1$.
]

#wide_figure(
  caption: [*Table II:* Head-to-Head Comparison of Bayesian Nonparametric and Model-Selection Mixture Frameworks (Section II: Algs. 82, 87, 93, 94, 99, 102).],
  table(
    columns: (1.02fr, 0.98fr, 0.95fr, 0.95fr, 0.85fr, 1.25fr),
    align: (left, left, left, center, center, left),
    [Algorithm], [Selection Mechanism], [Inference Engine], [Mixed / Missing], [Overlap Penalty], [Primary Strength / Failure Mode],
    [DPGMM (\#82) @blei2006variational], [Truncated GEM($alpha$) DP], [Mean-field CAVI], [Cont. (NIW)], [Prior stick decay], [Fast $O(n T d^2)$; var. shrinkage (Thm. 3.3)],
    [IGMM (\#99) @rasmussen2000igmm], [Exact untruncated CRP], [Collapsed Gibbs], [Cont. (NIW)], [CRP $alpha ln n$], [Exact posterior; slow MCMC split/merge],
    [AutoClass (\#87) @cheeseman1996autoclass], [Cheeseman--Stutz evid.], [Conjugate MAP-EM], [*Mixed + MAR*], [$O(1)$ missing-info], [Best mixed tables; diag. over-splits $rho != 0$],
    [SNOB (\#102) @wallace1968snob], [Two-part Fisher MML], [MML-EM + split/merge], [Mixed (diag.)], [$O(1)$ Fisher vol.], [Singularity-free small-$n$; diag. over-splits],
    [MCLUST (\#93) @fraley2002mclust], [BIC across 14 cov. models], [EM + HC init], [Cont. (`EII`--`VVV`)], [$0$ (pure density)], [Optimal parsimony; over-splits skew/tails],
    [Mixmod (\#94) @biernacki2006mixmod], [BIC, ICL, NEC], [EM / CEM / SEM], [Cont. or Cat.], [*$O(n)$ entropy (ICL)*], [Resists overlap splitting; fixed grid $K$],
  ),
)

// =============================================================================
// SECTION III
// =============================================================================
#sec("III", "Discrete, Count, Binary & Conceptual Mixtures")

Moving beyond continuous Euclidean spaces, five algorithms in Category IV cluster discrete counts, binary vectors, nominal contingency tables, and incremental feature streams: Poisson Mixture Models (PMM, Alg. 95), Bernoulli Mixture Models (BMM, Alg. 96), Latent Class Analysis (LCA, Alg. 90) @lazarsfeld1968latent, COBWEB (Alg. 85) @fisher1987cobweb, and CLASSIT (Alg. 86) @gennari1989classit. All five rely on *conditional feature independence* inside each cluster, yet they split across two fundamental boundaries: (i) an abrupt algebraic identifiability chasm between infinite-support counts (PMM) and finite-support categoricals (BMM/LCA), and (ii) the divergence geometry separating batch Shannon-likelihood EM from incremental Gini/precision Category Utility trees.

#subsec("A", "Identifiability Chasm & Local-Independence Covariance: PMM vs. BMM vs. LCA")

Under local independence, all three parametric discrete mixtures take the product-form density
$
  p(x | Theta) = sum_(k=1)^K pi_k product_(j=1)^d p_(k j) (x_j),
$
where $p_(k j)$ is $"Pois"(lambda_(k j))$ in PMM (\#95), $"Bern"(theta_(k j))$ in BMM (\#96), and $"Cat"(theta_(k j 1), dots, theta_(k j R_j))$ in LCA (\#90). Clearly BMM (\#96) is the exact binary restriction ($R_j = 2$) of LCA (\#90). Despite identical EM monotonicity and product structure, PMM and BMM/LCA behave oppositely regarding algebraic identifiability.

#theorem(title: "Discrete Identifiability Chasm, Overdispersion Law & Multinomial Reduction", origin: "Original Unified Proof")[
  1. *(Unconditional Identifiability of PMM \#95):* Every finite multivariate Poisson mixture on $bb(Z)_(>= 0)^d$ with pairwise distinct rate vectors $lambda_k in bb(R)_(> 0)^d$ and $pi_k > 0$ is *universally identifiable up to permutation* for every dimension $d >= 1$ and every cluster count $K >= 1$.
  2. *(Dimension-Counting Obstruction & Explicit Non-Identifiable Continuum for BMM \#96 / LCA \#90):* For $d$-dimensional binary BMM (\#96), generic identifiability is *impossible* whenever $K (d + 1) > 2^d$. Specifically, at $d = 2, K = 2$, the positively correlated binary table $p(1, 1) = p(0, 0) = 1 / 4 + c$, $p(1, 0) = p(0, 1) = 1 / 4 - c$ ($0 < c < 1 / 4$) is generated *identically* by the entire continuous 1-parameter family
     $
       pi_1 &= t, quad pi_2 = 1 - t, \
       theta_(1 1) &= theta_(1 2) = 1 / 2 + sqrt(c (1 - t) / t), \
       theta_(2 1) &= theta_(2 2) = 1 / 2 - sqrt(c t / (1 - t)),
     $
     for every $t in [4 c / (1 + 4 c), space 1 / (1 + 4 c)]$. Conversely, generic identifiability of LCA/BMM is restored once $d >= 3$ items satisfy Kruskal's tripartition rank bound $r_1 + r_2 + r_3 >= 2 K + 2$.
  3. *(Diagonal-Plus-Low-Rank Covariance & Exact Overdispersion):* Under any local-independence mixture (PMM, BMM, LCA), the population covariance matrix decomposes as $"Cov"(X) = D_"within" + Sigma_"between"$ with $"rank"(Sigma_"between") <= K - 1$. In PMM (\#95), every coordinate exhibits strict overdispersion $"Var"(X_j) = bb(E)[X_j] + "Var"_Z (lambda_(Z j)) > bb(E)[X_j]$ whenever $lambda_(k j)$ is non-constant across $k$; in BMM (\#96), marginal variance stays locked at $overline(theta)_j (1 - overline(theta)_j)$ while all cross-feature correlations reside strictly inside the rank-$<= (K - 1)$ matrix $"Cov"_Z (theta_(Z j), theta_(Z ell))$.
  4. *(Exact Conditioning Bridge from PMM \#95 to Multinomial Topic Mixtures):* Conditioned on the total row sum $L = sum_(j=1)^d X_j$, a multivariate PMM satisfies $X | (L, Z = k) tilde.op "Mult"(L, phi_k)$ with $phi_(k j) = lambda_(k j) / Lambda_k$ ($Lambda_k = sum_j lambda_(k j)$), and $L | (Z = k) tilde.op "Pois"(Lambda_k)$. Whenever component total intensities match ($Lambda_1 = dots = Lambda_K$), PMM posterior responsibilities $p(Z = k | X)$ depend *exclusively* on the normalized multinomial profile $phi_k$, coinciding identically with the multinomial mixture (DMM \#97 / single-topic pLSA \#84).
]

#proof[
  *Part 1:* Suppose $sum_(k=1)^M a_k product_(j=1)^d (e^(-lambda_(k j)) lambda_(k j)^(x_j) / x_j !) = 0$ for all $x in bb(Z)_(>= 0)^d$ with distinct $lambda_k in bb(R)_(> 0)^d$. Multiplying by $product_(j=1)^d s_j^(x_j)$ and summing over $bb(Z)_(>= 0)^d$ yields the exact multivariate probability generating function identity:
  $
    sum_(k=1)^M a_k exp( sum_(j=1)^d lambda_(k j) (s_j - 1) ) = 0 quad forall s in bb(R)^d.
  $
  Set $s_j = 1 + r u_j$ where $u in bb(R)^d$ is chosen outside the finite union of hyperplanes ${u : chevron.l u, lambda_k - lambda_j chevron.r = 0}$ ($k != j$). Then $sum_(k=1)^M a_k exp(r chevron.l u, lambda_k chevron.r) = 0$ for all $r > 0$. Dividing by $exp(r max_k chevron.l u, lambda_k chevron.r)$ and letting $r -> +infinity$ forces the coefficient of the unique maximizer $a_(k^*) = 0$; induction eliminates every $a_k$.

  *Part 2:* A binary distribution on ${0, 1}^d$ lies in the $(2^d - 1)$-simplex while BMM has $K - 1 + K d = K(d + 1) - 1$ free parameters; by Sard's theorem, the polynomial map from parameters to the simplex cannot be locally injective when $K(d + 1) > 2^d$. At $d = 2, K = 2$, set $a(t) := sqrt(c (1 - t) / t)$ and $b(t) := sqrt(c t / (1 - t))$. Feasibility $theta_(k j) in [0, 1]$ holds iff $a(t) <= 1 / 2$ and $b(t) <= 1 / 2$, i.e., $4 c / (1 + 4 c) <= t <= 1 / (1 + 4 c)$, which is a non-degenerate closed interval whenever $0 < c < 1 / 4$. Direct verification gives $t a(t) - (1 - t) b(t) = 0$ and $t a(t)^2 + (1 - t) b(t)^2 = c (1 - t) + c t = c$, whence:
  $
    p(1, 1) &= t (1 / 2 + a)^2 + (1 - t)(1 / 2 - b)^2 = 1 / 4 + c, \
    p(0, 0) &= t (1 / 2 - a)^2 + (1 - t)(1 / 2 + b)^2 = 1 / 4 + c, \
    p(1, 0) &= p(0, 1) = t (1 / 4 - a^2) + (1 - t)(1 / 4 - b^2) = 1 / 4 - c,
  $
  independent of $t$! For $d >= 3$, partitioning variables into three blocks expresses the joint probability tensor as a rank-$K$ 3-way tensor whose unique decomposition follows from Kruskal's theorem when $r_1 + r_2 + r_3 >= 2 K + 2$.

  *Part 3:* By the Law of Total Covariance,
  $
    "Cov"(X) = bb(E)_Z ["Cov"(X | Z)] + "Cov"_Z (bb(E)[X | Z]).
  $
  Conditional independence forces $"Cov"(X | Z = k)$ to be diagonal, so $D_"within" = bb(E)_Z ["Cov"(X | Z)]$ is diagonal. The second matrix $Sigma_"between" = sum_(k=1)^K pi_k (mu_k - overline(mu))(mu_k - overline(mu))^top$ is a sum of $K$ rank-1 terms obeying $sum_k pi_k (mu_k - overline(mu)) = bold(0)$, hence $"rank"(Sigma_"between") <= K - 1$. In PMM, $"Var"(X_j | Z) = lambda_(Z j)$, giving $"Var"(X_j) = bb(E)[lambda_(Z j)] + "Var"_Z (lambda_(Z j)) = bb(E)[X_j] + "Var"_Z (lambda_(Z j))$.

  *Part 4:* Because independent Poisson counts $X_j | (Z = k) tilde.op "Pois"(lambda_(k j))$ sum to $L | (Z = k) tilde.op "Pois"(Lambda_k)$,
  $
    p(X = x | Z = k) = ( e^(-Lambda_k) Lambda_k^L / (L!) ) dot ( L! / (product_j x_j !) product_(j=1)^d (lambda_(k j) / Lambda_k)^(x_j) ).
  $
  When $Lambda_k equiv Lambda$ across all $k$, the Poisson factor $e^(-Lambda) Lambda^L / L!$ cancels identically from Bayes' posterior $p(Z = k | X = x) prop pi_k product_j phi_(k j)^(x_j)$.
]

Notice also an important computational stability advantage of discrete mixtures over Section I's continuous models: because every probability $p_(k j) (x_j) <= 1$, the log-likelihood of PMM (\#95), BMM (\#96), and LCA (\#90) is *globally bounded above by $0$*—completely immune to continuous variance singularities! Instead, their primary numerical pathology is *zero-boundary absorption* ($hat(theta)_(k j) = 0$ or $hat(lambda)_(k j) = 0$ becomes an absorbing fixed point of multiplicative EM updates), effortlessly cured by symmetric Beta/Dirichlet or Gamma Laplace smoothing.

#subsec("B", "Gini vs. Shannon Divergence, Acuity Instability & XOR Blindness: COBWEB & CLASSIT")

COBWEB (\#85, @fisher1987cobweb) and its continuous extension CLASSIT (\#86, @gennari1989classit) depart from batch EM altogether, constructing a probabilistic concept hierarchy *incrementally in a single pass* using four local tree operators (incorporate, create singleton, merge best two, split best) guided by *Category Utility* (CU):
$
  "CU"_"COBWEB" (cal(C)) &= 1 / K sum_(k=1)^K pi_k sum_(j=1)^d sum_(v) [ P(A_j = v | C_k)^2 - P(A_j = v)^2 ], \
  "CU"_"CLASSIT" (cal(C)) &= 1 / (2 sqrt(pi) K) sum_(k=1)^K pi_k sum_(j=1)^d ( 1 / (max(sigma_(k j), sigma_"ac")) - 1 / (max(sigma_(0 j), sigma_"ac")) ).
$
Contrasting COBWEB/CLASSIT directly against LCA (\#90), BMM (\#96), and diagonal GMM (\#81) exposes four sharp theoretical trade-offs.

#theorem(title: "Gini-vs.-Shannon Divergence, Scale Hijacking, Acuity Blinding & Shared Parity Blindness", origin: "Original Unified Proof")[
  1. *(Amortized Gini vs. Unamortized Shannon Divergence):* For any $K$-partition $cal(C) = {C_1, dots, C_K}$ of a categorical dataset with empirical frequencies $P_(k j) (v)$ and root marginals $overline(P)_j (v)$, COBWEB's CU equals amortized squared-$ell_2$ (Gini) divergence, whereas hard-partition LCA's per-sample log-likelihood gain equals unamortized KL (Shannon) mutual information minus cluster entropy:
     $
       "CU"_"COBWEB" (cal(C)) &= 1 / K sum_(k=1)^K pi_k sum_(j=1)^d norm(P_(k j) - overline(P)_j)_2^2, \
       1 / n (ell_"LCA" (cal(C)) - ell_"null") &= sum_(k=1)^K pi_k sum_(j=1)^d D_"KL" (P_(k j) || overline(P)_j) - H(C).
     $
  2. *(Polynomial Scale Hijacking & Acuity Blinding in CLASSIT \#86 vs. Scale-Invariant Diagonal GMM \#81):* Rescaling coordinate $1$ by $x_(i 1) |-> s x_(i 1)$ ($s > 0$) leaves diagonal GMM / MCLUST-`VVI` partitions *identically invariant* (shifting every model's log-likelihood by $-n log s$). In CLASSIT (\#86) with $sigma_(k 1) > sigma_"ac"$, the same rescaling multiplies coordinate $1$'s CU contribution by $s^(-1) -> +infinity$ as $s -> 0^+$, allowing a single tiny-variance noise feature to dominate the entire concept tree! Conversely, whenever a feature's parent standard deviation falls below the acuity floor ($sigma_(0 j) <= sigma_"ac"$), its CU gain vanishes identically ($Delta "CU"_j = 0$), blinding CLASSIT completely to structure below scale $sigma_"ac"$.
  3. *(Minimal Order-Dependence Witness):* Unlike batch EM/BIC, incremental COBWEB (\#85) and CLASSIT (\#86) are *strictly order-dependent*: on the 3-point 1D stream $cal(X) = {0, 1, 2}$ with acuity $sigma_"ac" = 0.55$, presenting $(0, 1, 2)$ yields the strict unique partition ${{0}, {1, 2}}$, whereas the reversed stream $(2, 1, 0)$ yields the distinct partition ${{0, 1}, {2}}$ with zero tie-breaking ambiguity.
  4. *(Universal Higher-Order Parity Blindness of Local-Independence Models):* On the 2-class 3D binary parity dataset $C_1 = {x in {0, 1}^3 : x_1 + x_2 + x_3 equiv 0 space ("mod" 2)}$ ($4$ points) and $C_2 = {0, 1}^3 backslash C_1$ ($4$ points), *every* univariate marginal matches the root ($P(X_j = 1 | C_k) = 1 / 2 = P(X_j = 1)$ for all $j, k$). Consequently, COBWEB assigns $"CU"({C_1, C_2}) = 0$, and $K = 2$ BMM/LCA assigns zero mutual information ($sum_j I(X_j; C) = 0$) and strictly lower BIC than $K = 1$, proving that recovering parity interactions requires either $K >= 4$ local-independence subclusters or non-diagonal covariance (GMM \#81 / MCLUST `VVV` \#93).
]

#proof[
  *Part 1:* Expand $norm(P_(k j) - overline(P)_j)_2^2 = sum_v P_(k j)(v)^2 - 2 sum_v P_(k j)(v) overline(P)_j (v) + sum_v overline(P)_j (v)^2$. Averaging over $k = 1, dots, K$ with weights $pi_k$ uses $sum_k pi_k P_(k j)(v) = overline(P)_j (v)$ to simplify $-2 sum_v overline(P)_j (v)^2 + sum_v overline(P)_j (v)^2 = -sum_v overline(P)_j (v)^2$, giving
  $
    sum_(k=1)^K pi_k norm(P_(k j) - overline(P)_j)_2^2 = sum_(k=1)^K pi_k sum_v [P_(k j)(v)^2 - overline(P)_j (v)^2].
  $
  Dividing by $K$ and summing over $j$ proves the exact COBWEB identity. For hard-partition LCA,
  $
    1 / n ell_"LCA" (cal(C)) &= sum_(k=1)^K pi_k log pi_k + sum_(k=1)^K pi_k sum_(j=1)^d sum_v P_(k j)(v) log P_(k j)(v), \
    1 / n ell_"null" &= sum_(j=1)^d sum_v overline(P)_j (v) log overline(P)_j (v) \
    &= sum_(k=1)^K pi_k sum_(j=1)^d sum_v P_(k j)(v) log overline(P)_j (v),
  $
  whose difference is identically $sum_k pi_k sum_j D_"KL" (P_(k j) || overline(P)_j) - H(C)$.

  *Part 2:* In diagonal GMM, mapping $x_(i 1) |-> s x_(i 1)$ rescales $(mu_(k 1), sigma_(k 1)) |-> (s mu_(k 1), s sigma_(k 1))$, preserving every Mahalanobis term $((x_(i 1) - mu_(k 1)) / sigma_(k 1))^2$ and shifting every log-density by the partition-independent constant $-log s$. In CLASSIT, when $s sigma_(k 1) >= sigma_"ac"$, the term $(2 sqrt(pi))^(-1) ((s sigma_(k 1))^(-1) - (s sigma_(0 1))^(-1)) = s^(-1) Delta_1$ scales as $s^(-1) -> +infinity$ as $s -> 0^+$, overwhelming all unit-scale features $j >= 2$. Conversely, whenever a feature's parent spread obeys $sigma_(0 j) <= sigma_"ac"$ and every child variance obeys $sigma_(k j) <= sigma_"ac"$ (as holds for all variance-reducing splits), both clamped denominators equal $sigma_"ac"$, yielding $sigma_"ac"^(-1) - sigma_"ac"^(-1) = 0$.

  *Part 3:* Run standard CLASSIT with $sigma_"ac" = 0.55$ on $cal(X) = {0, 1, 2}$.
  - *Order $(0, 1, 2)$:* Instances $0$ and $1$ form two singleton leaves ${0}$ and ${1}$ under the root. Upon inserting $x_3 = 2$, the root population standard deviation is $sigma_0 = sqrt(2/3) approx 0.8165$ ($sigma_0^(-1) approx 1.2247$), while $sigma_"ac"^(-1) approx 1.8182$. Evaluating all four CLASSIT operators gives: (a) incorporate $2$ into ${1}$ ($sigma({1, 2}) = 0.5 < 0.55$, both children clamped at $0.55$): $"CU" = (1.8182 - 1.2247) / (2 dot 2 sqrt(pi)) approx 0.2967 / (2 sqrt(pi))$; (b) incorporate $2$ into ${0}$ ($sigma({0, 2}) = 1.0$): $"CU" = ((2/3)(1.0) + (1/3)(1.8182) - 1.2247) / (2 dot 2 sqrt(pi)) approx 0.0240 / (2 sqrt(pi))$; (c) create a third leaf ${2}$: $"CU" = (1.8182 - 1.2247) / (3 dot 2 sqrt(pi)) approx 0.1978 / (2 sqrt(pi))$; (d) merge ${0}, {1}$ and incorporate $2$: $"CU" = 0$. Thus $(0, 1, 2)$ strictly selects ${{0}, {1, 2}}$.
  - *Order $(2, 1, 0)$:* By exact reflection $x |-> 2 - x$, instances $2$ and $1$ form leaves ${2}, {1}$, and inserting $x_3 = 0$ strictly selects ${{0, 1}, {2}}$ with identical margin $0.2967 > 0.1978 > 0.0240 > 0$.

  *Part 4:* Direct counting inside $C_1 = {(0,0,0), (0,1,1), (1,0,1), (1,1,0)}$ and $C_2 = {(1,0,0), (0,1,0), (0,0,1), (1,1,1)}$ yields exactly two $1$s and two $0$s per coordinate in each class. Hence $P_(k j) = overline(P)_j = (1/2, 1/2)$ for every $k in {1, 2}$ and $j in {1, 2, 3}$, forcing $"CU"_"COBWEB" = 0$ and $D_"KL" (P_(k j) || overline(P)_j) = 0$ identically.
]

#wide_figure(
  caption: [*Table III:* Head-to-Head Comparison of Discrete, Count, Binary, and Incremental Conceptual Clustering Algorithms (Section III: Algs. 85, 86, 90, 95, 96).],
  table(
    columns: (1.05fr, 0.95fr, 0.95fr, 0.92fr, 0.88fr, 1.25fr),
    align: (left, left, center, center, center, left),
    [Algorithm], [Data Domain & Model], [Objective / Div.], [Identifiable?], [Order-Inv.?], [Primary Strength / Failure Mode],
    [PMM (\#95)], [Count $bb(Z)_(>= 0)^d$, Pois. mix], [I-div. / Log-lik. EM], [*Always* (Thm. 5.1)], [Batch EM], [Models overdispersion; zero-lock EM],
    [BMM (\#96)], [Binary ${0, 1}^d$, Bern. mix], [KL / Log-lik. EM], [Needs $d >= 3$ (Kruskal)], [Batch EM], [Sparse binary vectors; $d=2$ continuum],
    [LCA (\#90) @lazarsfeld1968latent], [Nominal tables, Cat. mix], [Shannon MI / EM], [Kruskal rank bound], [Batch EM], [Survey/psychometrics; parity blindness],
    [COBWEB (\#85) @fisher1987cobweb], [Nominal stream, tree], [Amort. Gini CU], [Heuristic tree], [*No* (1-pass)], [Incremental anytime tree; order & XOR],
    [CLASSIT (\#86) @gennari1989classit], [Continuous stream, tree], [$sigma^(-1)$ precision CU], [Heuristic tree], [*No* (1-pass)], [Streaming numeric tree; scale & acuity],
  ),
)

// =============================================================================
// SECTION IV
// =============================================================================
#sec("IV", "Topic Models & Nonnegative Matrix Factorization")

Seven algorithms in Category IV perform soft semantic clustering and latent-factor discovery on document--term matrices, text corpora, and nonnegative count tables: pLSA (Alg. 84) @hofmann1999plsa, LDA (Alg. 83) @blei2003lda, HDP (Alg. 91) @teh2006hdp, CTM (Alg. 100) @blei2007ctm, DMM (Alg. 97) @yin2014dmm, BTM (Alg. 101) @yan2013btm, and NMF (Alg. 92) @lee1999nmf @lee2001algorithms. Unlike Section III's instance-level mixtures (which assign one cluster $z_i$ per object), pLSA, LDA, HDP, CTM, and NMF are *mixed-membership (admixture)* models allowing each document $d$ to draw tokens from multiple latent topics simultaneously, whereas DMM restores single-topic document coherence and BTM shifts modeling from document-level tokens to corpus-level word pairs.

#subsec("A", "Exact pLSA–KL-NMF Equivalence & the Generative Overfitting Gap to LDA")

Although NMF (\#92) originated in linear algebra as a parts-based factorization $X approx W H$ ($W in bb(R)_(>= 0)^(V times K)$, $H in bb(R)_(>= 0)^(K times D)$) and pLSA (\#84) originated in statistical NLP as a latent-aspect mixture, the two objectives coincide identically under Kullback--Leibler divergence @ding2008nmfplsa—and both expose the exact statistical defect that motivated LDA (\#83) @blei2003lda.

#theorem(title: "Exact pLSA–KL-NMF Equivalence and Neyman–Scott Overfitting Separation from LDA", origin: "Original Unified Proof")[
  Let $X in bb(R)_(>= 0)^(V times D)$ be a term--document count matrix with total token count $N = sum_(w=1)^V sum_(d=1)^D X_(w d) > 0$.
  1. *(Exact Global Objective & KKT Equivalence of KL-NMF \#92 and pLSA \#84):* For any positive factorization $(W, H) in bb(R)_(> 0)^(V times K) times bb(R)_(> 0)^(K times D)$, decompose $c_k = sum_w W_(w k)$, $r_k = sum_d H_(k d)$, $Phi_(w k) = W_(w k) / c_k$, $Psi_(d k) = H_(k d) / r_k$, $M = sum_k c_k r_k$, and $pi_k = (c_k r_k) / M$. Then the generalized KL-NMF divergence satisfies the exact decomposition
     $
       D_"KL" (X || W H) &= -ell_"pLSA" (Phi, Psi, pi) + (M - N log M - N) \
       &quad + sum_(w, d) X_(w d) log X_(w d).
     $
     Minimizing over the global scale $M > 0$ uniquely fixes $M^* = N$, establishing an exact bijection between global optima (and interior stationary points) of KL-NMF (\#92) and symmetric pLSA (\#84). By contrast, Frobenius NMF ($norm(X - W H)_F^2$) corresponds to homoscedastic Gaussian noise and reduces to hard Euclidean $k$-Means under $H H^top = I_K$.
  2. *(Incidental-Parameter Overfitting of pLSA/NMF vs. Constant-Capacity LDA \#83):* Let a corpus contain $D$ documents of bounded length $L_d <= L_"max"$. The number of free parameters in pLSA (\#84) and NMF (\#92) grows linearly with corpus size:
     $
       p_"pLSA" (D) = K (V - 1) + D (K - 1) + (K - 1) = Theta(D K),
     $
     so the parameter-to-token ratio $(p_"pLSA" (D)) / N >= (K - 1) / L_"max"$ *never vanishes* as $D -> infinity$ (violating classical MLE consistency via the Neyman--Scott incidental-parameter obstruction) and leaving unseen documents outside the probability simplex without ad-hoc fold-in. Integrating out $theta_d tilde.op "Dir"(alpha)$ in LDA (\#83) collapses global parameter count to $p_"LDA" = K V$ independent of $D$ (ratio $(K V) / (D overline(L)) -> 0$ as $D -> infinity$) while defining a coherent exchangeable generative law on unseen documents. Moreover, MAP-pLSA with symmetric Dirichlet priors $"Dir"(alpha), "Dir"(beta)$ recovers the unmarginalized posterior mode of LDA and collapses identically to raw pLSA at $alpha = bold(1), beta = bold(1)$.
]

#proof[
  *Part 1:* By definition of generalized I-divergence,
  $
    D_"KL" (X || W H) = sum_(w, d) [ X_(w d) log(X_(w d) / (W H)_(w d)) - X_(w d) + (W H)_(w d) ].
  $
  Substituting $(W H)_(w d) = sum_(k=1)^K (c_k r_k) Phi_(w k) Psi_(d k) = M sum_(k=1)^K pi_k Phi_(w k) Psi_(d k)$ and noting $sum_(w, d) (W H)_(w d) = M sum_k pi_k (sum_w Phi_(w k))(sum_d Psi_(d k)) = M$ gives
  $
    D_"KL" (X || W H) &= sum_(w, d) X_(w d) log X_(w d) - N + M - N log M \
    &quad - sum_(w, d) X_(w d) log( sum_(k=1)^K pi_k Phi_(w k) Psi_(d k) ).
  $
  The final sum is identically $-ell_"pLSA" (Phi, Psi, pi)$. Differentiating $g(M) = M - N log M$ yields $g'(M) = 1 - N / M = 0 <=> M^* = N$, with $g''(N) = 1 / N > 0$. Because every KKT stationary point of KL-NMF satisfies $sum_(w, d) (W H)_(w d) = N$ (summing the Lee--Seung stationarity condition $sum_w W_(w k) = sum_(w, d) W_(w k) X_(w d) / (W H)_(w d) H_(k d)$ over $k$), the interior stationary sets of KL-NMF and pLSA coincide identically @ding2008nmfplsa.

  *Part 2:* In pLSA, each document introduces an unregularized $(K - 1)$-simplex vector $theta_d$, contributing $D(K - 1)$ incidental parameters. Because total tokens satisfy $N <= D L_"max"$, $(p_"pLSA" (D)) / N >= (K - 1) / L_"max" > 0$ uniformly for all $D >= 1$; on a corpus of length-1 documents ($L_d = 1$), pLSA trivially interpolates every empirical document count with zero perplexity gain on held-out text. In LDA, marginalizing $p(w_d | alpha, Phi) = integral_(Delta^(K-1)) "Dir"(theta_d | alpha) product_(t=1)^(L_d) (sum_k theta_(d k) Phi_(w_(d t), k)) dif theta_d$ eliminates all $D(K - 1)$ per-document parameters from the global model specification, leaving $K(V - 1) + K = K V$ parameters globally and valid out-of-sample likelihoods for any new document $w_"new"$.
]

#subsec("B", "Simplex Correlation Obstruction (LDA vs. CTM) & Cross-Document Sharing (HDP)")

LDA's Dirichlet prior $theta_d tilde.op "Dir"(alpha)$ buys exact conjugacy and $O(K N_"nnz")$ inference, yet imposes two structural rigidities: (i) $K$ must be specified upfront, and (ii) distinct topics cannot positively co-occur in the prior. HDP (\#91, @teh2006hdp) removes the fixed-$K$ rigidity via a two-level Dirichlet process, while CTM (\#100, @blei2007ctm) removes simplex independence via a logistic-normal prior $eta_d tilde.op cal(N)(mu, Sigma)$, $theta_(d k) = e^(eta_(d k)) / sum_j e^(eta_(d j))$.

#theorem(title: "Dirichlet Negative-Correlation Obstruction, CTM Positive Realizability & HDP Atom-Sharing Correlation", origin: "Original Unified Proof")[
  1. *(Strict Negative-Correlation Obstruction of LDA \#83):* Under any Dirichlet prior $theta tilde.op "Dir"(alpha_1, dots, alpha_K)$ ($alpha_k > 0$, $alpha_0 = sum_k alpha_k$), every pair of distinct topics $j != k$ is *strictly negatively correlated*:
     $
       "Cov"(theta_j, theta_k) &= - (alpha_j alpha_k) / (alpha_0^2 (alpha_0 + 1)) < 0, \
       "Corr"(theta_j, theta_k) &= - sqrt( (alpha_j alpha_k) / ((alpha_0 - alpha_j)(alpha_0 - alpha_k)) ) < 0.
     $
  2. *(Constructive Positive-Correlation Realizability & Non-Conjugacy Cost of CTM \#100):* For every $K >= 3$, CTM's logistic-normal prior achieves *strictly positive* pairwise topic correlation $"Corr"(theta_1, theta_2) > 0$ (approaching $+1^-$ as $Sigma_(1 2) / sqrt(Sigma_(1 1) Sigma_(2 2)) -> 1^-$ with $Sigma_(1 1) = Sigma_(2 2) = sigma^2 > 0$). However, loss of Dirichlet--multinomial conjugacy requires per-document variational linearization ($log sum_k e^(eta_(d k)) <= zeta_d^(-1) sum_k e^(eta_(d k)) + log zeta_d - 1$), inflating per-iteration complexity from LDA's $O(K N_"nnz")$ to $O(K N_"nnz" + D K^2 + K^3)$.
  3. *(Exact Cross-Document Measure Correlation & Double-Logarithmic Topic Scaling of HDP \#91):* Under the Hierarchical Dirichlet Process $G_0 tilde.op "DP"(gamma, H)$, $G_d | G_0 tilde.op "DP"(alpha_0, G_0)$ i.i.d., distinct documents $d != d'$ share identical discrete topic atoms almost surely with *exact* cross-document measure correlation for every measurable set $A$ ($0 < H(A) < 1$):
     $
       "Corr"(G_d (A), G_(d') (A)) = (alpha_0 + 1) / (alpha_0 + gamma + 1) in (0, 1).
     $
     Moreover, across $D$ documents of length $L$, the expected number of global topics discovered by the Chinese Restaurant Franchise satisfies the doubly-compressed logarithmic bound
     $
       bb(E)[K_(D, L)] <= gamma ln( 1 + (D alpha_0) / gamma (1 + ln(1 + L / alpha_0)) ) + 1.
     $
]

#proof[
  *Part 1:* Standard Dirichlet moments give $bb(E)[theta_j] = alpha_j / alpha_0$, $"Var"(theta_j) = alpha_j (alpha_0 - alpha_j) / (alpha_0^2 (alpha_0 + 1))$, and $bb(E)[theta_j theta_k] = alpha_j alpha_k / (alpha_0 (alpha_0 + 1))$ ($j != k$). Subtracting $bb(E)[theta_j] bb(E)[theta_k]$ gives $"Cov"(theta_j, theta_k) = -alpha_j alpha_k / (alpha_0^2 (alpha_0 + 1)) < 0$, and dividing by $sqrt("Var"(theta_j) "Var"(theta_k))$ yields the exact negative correlation.

  *Part 2:* Fix $K = 3$ with canonical gauge $eta_3 equiv 0$, $mu = bold(0)$, and bivariate covariance $Sigma(rho) = sigma^2 mat(1, rho; rho, 1)$ ($sigma > 0$, $rho in (0, 1)$). At the degenerate limit $rho = 1$, $eta_1 = eta_2 = Z tilde.op cal(N)(0, sigma^2)$, so $theta_1 = theta_2 = g(Z) := e^Z / (2 e^Z + 1)$. Because $g$ is strictly increasing and bounded ($g in (0, 1/2)$), $"Var"(g(Z)) > 0$, yielding $"Cov"(theta_1, theta_2) = "Var"(g(Z)) > 0$ and $"Corr"(theta_1, theta_2) = +1$. By dominated convergence on bounded $theta_1, theta_2 in (0, 1)$, $"Corr"_(Sigma(rho)) (theta_1, theta_2)$ is continuous in $rho in [0, 1]$ and therefore strictly positive (and arbitrarily close to $+1$) for all $rho$ near $1$ with strictly positive-definite $Sigma(rho) succ 0$. (Embedding $K > 3$ sets $mu_j -> -infinity$ for $j >= 4$.) For computational cost, CTM's M-step estimates and inverts the full $K times K$ covariance $Sigma$ ($O(D K^2 + K^3)$) while each document's variational E-step optimizes $(lambda_d, nu_d, zeta_d)$ via Newton/conjugate-gradient iterations @blei2007ctm.

  *Part 3:* Let $p = H(A) in (0, 1)$ and $q = 1 - p$. By Dirichlet-process marginals, $G_0 (A) tilde.op "Beta"(gamma p, gamma q)$ and $G_d (A) | G_0 tilde.op "Beta"(alpha_0 G_0 (A), alpha_0 (1 - G_0 (A)))$ conditionally independently across $d$. Hence $bb(E)[G_d (A) | G_0] = G_0 (A)$ and $"Var"(G_d (A) | G_0) = G_0 (A)(1 - G_0 (A)) / (alpha_0 + 1)$. For $d != d'$, the Law of Total Covariance gives
  $
    "Cov"(G_d (A), G_(d') (A)) = "Var"(G_0 (A)) = (p q) / (gamma + 1).
  $
  Meanwhile, since $bb(E)[G_0 (A)(1 - G_0 (A))] = p q - "Var"(G_0 (A)) = p q gamma / (gamma + 1)$, the Law of Total Variance gives
  $
    "Var"(G_d (A)) &= "Var"(G_0 (A)) + bb(E)["Var"(G_d (A) | G_0)] \
    &= (p q) / (gamma + 1) ( 1 + gamma / (alpha_0 + 1) ) = (p q (alpha_0 + gamma + 1)) / ((alpha_0 + 1)(gamma + 1)).
  $
  Dividing $"Cov"(G_d (A), G_(d') (A))$ by $"Var"(G_d (A))$ yields $(alpha_0 + 1) / (alpha_0 + gamma + 1)$ identically! For the global topic count $K_(D, L)$, let $M_d$ be the number of local tables in restaurant $d$; by Theorem 3(2), $bb(E)[M_d] = sum_(t=1)^L alpha_0 / (alpha_0 + t - 1) <= 1 + alpha_0 ln(1 + L / alpha_0)$. Conditioned on total tables $M = sum_(d=1)^D M_d$, global dishes obey $bb(E)[K_(D, L) | M] = sum_(m=1)^M gamma / (gamma + m - 1) <= 1 + gamma ln(1 + M / gamma)$. Applying Jensen's inequality to the concave function $m |-> ln(1 + m / gamma)$ gives
  $
    bb(E)[K_(D, L)] <= 1 + gamma ln( 1 + (bb(E)[M]) / gamma ),
  $
  which yields the stated double-logarithmic ceiling upon substituting $bb(E)[M] <= D (1 + alpha_0 ln(1 + L / alpha_0))$.
]

#subsec("C", "Short-Text Sparsity Collapse of LDA vs. DMM and BTM Remedies")

Standard topic models (pLSA \#84, LDA \#83, HDP \#91, CTM \#100) were architected for long documents ($L >= 100$--$1000$ tokens). On short texts—tweets, search queries, chat messages, or clinical notes with $L in [3, 15]$—LDA suffers severe semantic degradation. Two specialized Category IV models resolve this failure through radically different probabilistic architectures: *Dirichlet Multinomial Mixture* (DMM / GSDMM, Alg. 97) @yin2014dmm enforces a single latent topic per document, whereas *Biterm Topic Model* (BTM, Alg. 101) @yan2013btm elevates topic generation from individual document tokens to corpus-wide unordered word pairs (*biterms*). We now prove rigorously *why* LDA collapses on short texts and quantify the exact statistical advantages and trade-offs of DMM and BTM.

#theorem(title: "Short-Text Sparsity Collapse of LDA vs. Exponential DMM Margin & Corpus-Level BTM Variance Reduction", origin: "Original Unified Proof")[
  Consider a corpus of $D$ short documents of length $L >= 2$ generated from $K >= 2$ topics with word distributions $phi_1, dots, phi_K in Delta^(V - 1)$ and symmetric Dirichlet concentration $alpha > 0$.
  1. *(Logarithmic Context Bottleneck & Prior Fragmentation of LDA \#83):*
     - In LDA's collapsed Gibbs conditional for token $w_t$ in document $d$, even when *all* $L - 1$ remaining tokens in $d$ occupy the true topic $k^*$, the log-odds contextual reinforcement favoring $k^*$ over an empty rival topic $k$ grows only *logarithmically* in document length:
       $
         Delta_"LDA" (k^* "vs." k) = log( 1 + (L - 1) / alpha ).
       $
     - Moreover, for a 2-word document ($L = 2$), LDA's prior probability that both words share a common topic is only $(alpha + 1) / (K alpha + 1)$, severing pairwise co-occurrence with probability $((K - 1) alpha) / (K alpha + 1) -> 1$ as $K alpha -> infinity$.
     - Even under oracle token labels ($n_(d, k^*) = L$), LDA's per-document proportion estimate suffers shrinkage bias $(K - 1) alpha / (K alpha + L)$ and estimation variance $"Var"(theta_(d, k^*) | z_d) = Theta((K alpha) / (K alpha + L)^2)$.
  2. *(Linear-in-$L$ Exponential Separation of DMM \#97):* Under DMM's single-topic-per-document assumption $z_d in {1, dots, K}$ @yin2014dmm, all $L$ words pool coherently into one joint posterior log-odds ratio against any rival topic $k != k^*$:
     $
       &bb(E)_(phi_(k^*)) [ log( (p_"DMM" (z_d = k^* | w_d)) / (p_"DMM" (z_d = k | w_d)) ) ] \
       &quad = L dot D_"KL" (phi_(k^*) || phi_k) + log(pi_(k^*) / pi_k),
     $
     growing *linearly* in $L$ (exponentially stronger than LDA's $log(1 + (L-1)/alpha)$ token coupling!) and driving document misclassification probability to $0$ at Chernoff rate $exp(-L C(phi_(k^*), phi_k))$.
  3. *(Exact $O((D L)^(-1))$ Variance Compression & Quadratic Length-Skew Pathology of BTM \#101):* Extracting all corpus biterms $cal(B) = union.big_(d=1)^D {(w_i, w_j) : 1 <= i < j <= L_d}$ ($|cal(B)| = sum_(d=1)^D binom(L_d, 2)$) and modeling a single global topic proportion $theta tilde.op "Dir"(alpha bold(1)_K)$ @yan2013btm reduces mixing-proportion posterior variance relative to LDA's per-document $theta_d$ on uniform length-$L$ corpora by the exact corpus factor:
     $
       ("Var"_"BTM" (theta_k | z_cal(B))) / ("Var"_"LDA" (theta_(d, k) | z_d)) = Theta( 1 / (D L) ) -> 0 quad (D -> infinity),
     $
     while preserving multi-topic soft inference $p(z | d) = sum_(b in cal(B)_d) p(z | b) p(b | d)$. *Conversely*, BTM exhibits *quadratic length-skew dominance* (zero asymptotic breakdown to a single long outlier document): inserting one long document of length $L_"out" = ceil(sqrt(D) L)$ into $D$ short documents of length $L$ captures fraction
     $
       (binom(L_"out", 2)) / (binom(L_"out", 2) + D binom(L, 2)) >= 1 / 2
     $
     of the entire corpus biterm likelihood!
]

#proof[
  *Part 1:* In collapsed LDA, integrating out $theta_d tilde.op "Dir"(alpha bold(1)_K)$ yields the exact token conditional
  $
    p_"LDA" (z_(d, t) = k | z_(d, -t), w_d, Phi) prop (n_(d, -t, k) + alpha) phi_(k, w_t).
  $
  When $n_(d, -t, k^*) = L - 1$ and $n_(d, -t, k) = 0$, the log-ratio of the document-context factors is identically $log((L - 1 + alpha) / alpha) = log(1 + (L - 1) / alpha)$. For $L = 2$, $bb(P)(z_(d, 2) = z_(d, 1)) = bb(E)[sum_(k=1)^K theta_(d, k)^2] = K dot alpha (alpha + 1) / (K alpha (K alpha + 1)) = (alpha + 1) / (K alpha + 1)$, so $bb(P)(z_(d, 2) != z_(d, 1)) = (K - 1) alpha / (K alpha + 1)$. Under oracle counts $(L, 0, dots, 0)$, the conjugate posterior is $"Dir"(L + alpha, alpha, dots, alpha)$ with mean $(L + alpha) / (K alpha + L) = 1 - (K - 1) alpha / (K alpha + L)$ and variance $(L + alpha)(K - 1) alpha / ((K alpha + L)^2 (K alpha + L + 1))$.

  *Part 2:* Under DMM, $p(w_d | z_d = k) = product_(t=1)^L phi_(k, w_t)$, so by Bayes' rule
  $
    log( (p_"DMM" (z_d = k^* | w_d)) / (p_"DMM" (z_d = k | w_d)) ) = log(pi_(k^*) / pi_k) + sum_(t=1)^L log(phi_(k^*, w_t) / phi_(k, w_t)).
  $
  Taking expectation under $w_t tilde.op phi_(k^*)$ i.i.d. yields $L D_"KL" (phi_(k^*) || phi_k) + log(pi_(k^*) / pi_k)$ immediately. By the Chernoff--Stein theorem for hypothesis testing of $L$ i.i.d. draws from $phi_(k^*)$ vs. $phi_k$, the Bayes error decays as $exp(-L C(phi_(k^*), phi_k))$ where $C(phi_(k^*), phi_k) = -min_(lambda in [0, 1]) log sum_w phi_(k^*, w)^lambda phi_(k, w)^(1 - lambda) > 0$ whenever $phi_(k^*) != phi_k$.

  *Part 3:* In BTM with uniform document length $L >= 2$, $|cal(B)| = D L (L - 1) / 2$. For a balanced topic holding fraction $rho_k in (0, 1)$ of biterms, the conjugate posterior $theta | z_cal(B) tilde.op "Dir"(alpha + n_1, dots, alpha + n_K)$ has variance
  $
    "Var"_"BTM" (theta_k | z_cal(B)) &= ((alpha + rho_k |cal(B)|)( (K - 1) alpha + (1 - rho_k)|cal(B)| )) / ((K alpha + |cal(B)|)^2 (K alpha + |cal(B)| + 1)) \
    &= Theta( 1 / (D L^2) ).
  $
  Comparing against LDA's single-document proportion variance $"Var"_"LDA" (theta_(d, k) | z_d) = Theta(1 / L)$ (for fixed $K alpha = O(L)$) gives the exact ratio $Theta(L / (D L^2)) = Theta((D L)^(-1))$. Finally, for length-skew dominance, inserting one outlier document of length $L_"out" = ceil(sqrt(D) L)$ adds $binom(L_"out", 2) = L_"out" (L_"out" - 1) / 2 >= (D L^2 - sqrt(D) L) / 2 >= D L (L - 1) / 2 = D binom(L, 2)$ for all $D >= 1, L >= 2$ (verify: $D L^2 - sqrt(D) L >= D L^2 - D L <=> D >= sqrt(D)$, always true for $D >= 1$!), proving the single outlier accounts for at least $50\%$ of all biterms in $cal(B)$!
]

Theorem 9 clarifies the sharp three-way architectural divide in text clustering: *LDA (\#83)* (and *HDP \#91* / *CTM \#100*) dominate long multi-topic documents ($L >= 50$); *DMM (\#97)* dominates short single-theme texts ($L <= 30$) where fast $O(T K overline(L) D)$ hard document clustering and automatic empty-cluster pruning (GSDMM) are desired; and *BTM (\#101)* dominates short *multi-topic* corpora with bounded document lengths (using a fixed sliding window $|i - j| <= W$ whenever length variance is high to neutralize Theorem 9(3)'s quadratic length-skew dominance!).

#wide_figure(
  caption: [*Table IV:* Head-to-Head Theoretical & Algorithmic Comparison of Topic Models and Matrix Factorization (Section IV: Algs. 83, 84, 91, 92, 97, 100, 101).],
  table(
    columns: (1.02fr, 0.95fr, 0.95fr, 0.98fr, 0.88fr, 1.22fr),
    align: (left, left, center, left, center, left),
    [Algorithm], [Latent Granularity], [Global Params], [Time / Iteration], [Short Text / Corr.], [Primary Strength / Failure Mode],
    [pLSA (\#84) @hofmann1999plsa], [Token-level aspect], [$Theta(D K + K V)$], [$O(K N_"nnz")$], [Poor / Unconstr.], [Simple EM; linear overfit, no fold-in],
    [NMF (\#92) @lee2001algorithms], [Parts-based $W H$], [$Theta(D K + K V)$], [$O(K N_"nnz")$ / $O(D V K)$], [Moderate / Unconstr.], [Sparse interpretable parts; non-unique],
    [LDA (\#83) @blei2003lda], [Token admixture, Dir.], [$K V$ (indep. $D$)], [$O(K N_"nnz")$], [Collapses / Neg. only], [Generative gold standard; short-$L$ fail],
    [HDP (\#91) @teh2006hdp], [Token admixture, 2-DP], [Nonparam. ($K$ auto)], [$O(N K_"eff")$], [Collapses / Neg. only], [Auto-$K$ atom sharing; slow CRF Gibbs],
    [CTM (\#100) @blei2007ctm], [Token admixture, Log-N], [$K V + K^2 / 2$], [$O(K N_"nnz" + D K^2 + K^3)$], [Collapses / *Pos.+Neg.*], [Captures topic graph; $O(K^3)$ non-conj.],
    [DMM (\#97) @yin2014dmm], [*1 topic per doc*], [$K V + K$], [$O(K overline(L) D)$], [*Optimal* / 1-topic], [Fastest short-text clusters; no multi-top.],
    [BTM (\#101) @yan2013btm], [*Corpus biterm pairs*], [$K V + K$], [$O(K |cal(B)|)$], [*Strong* / Global $theta$], [Multi-topic short text; $L^2$ skew vuln.],
  ),
)

// =============================================================================
// SECTION V
// =============================================================================
#sec("V", "Complete 22-Algorithm Head-to-Head Matrix & Practitioner's Decision Guide")

We now synthesize all twenty-two algorithms of Category IV (Algs. 81--102) in a unified head-to-head comparative matrix (Table V), followed by exact parametric reduction laws, a mathematically grounded practitioner's selection guide, and open theoretical frontiers connecting Category IV to the remainder of the monograph.

#wide_figure(
  caption: [*Table V:* Master Head-to-Head Theoretical and Practical Comparison Matrix of All 22 Distribution, Statistical, and Topic Models (Category IV, Algorithms 81--102). Abbreviations: *PM* = Finite Parametric Continuous, *BN* = Bayesian Nonparametric / Model Selection, *DC* = Discrete / Conceptual, *TM* = Topic / Matrix Factorization; *Bnd.* = Bounded Likelihood; *Unb.* = Unbounded Singularity without Reg.; *Auto-$K$* = Automatic Cluster/Topic Cardinality Selection.],
  {
    set text(size: 6.7pt)
    set table(inset: (x: 2.8pt, y: 2.2pt))
    table(
      columns: (1.22fr, 0.36fr, 0.96fr, 1.04fr, 0.72fr, 0.66fr, 0.72fr, 1.28fr),
      align: (left, center, left, left, center, center, center, left),
      [Algorithm], [Grp.], [Data Domain], [Time Complexity], [MLE / Reg.], [Auto-$K$], [Robustness], [Distinctive Strength / Weakness],
      [\#81 GMM @dempster1977em], [PM], [Cont. $bb(R)^d$], [$O(T n K d^2 + T K d^3)$], [Unb. / Floor], [No], [Low ($1/n$)], [Affine-equiv.; singular & $n_k <= d$ fail],
      [\#82 DPGMM @blei2006variational], [BN], [Cont. $bb(R)^d$], [$O(T n K d^2 + T K d^3)$], [*Bnd.* (NIW)], [*Yes* (GEM)], [Med (NIW)], [Fast VI auto-$K$; CAVI var. shrinkage],
      [\#83 LDA @blei2003lda], [TM], [Doc--term counts], [$O(T_"out" (T_"in" K N_"nnz" + K V))$], [*Bnd.* (Dir)], [No], [Long text], [Constant params; short-text collapse],
      [\#84 pLSA @hofmann1999plsa], [TM], [Doc--term counts], [$O(T K N_"nnz")$], [Bnd. / Overfit], [No], [Low], [Exact KL-NMF dual; $Theta(D K)$ overfit],
      [\#85 COBWEB @fisher1987cobweb], [DC], [Nominal stream], [$O(n B^2 d V h)$], [Gini CU], [*Yes* (Tree)], [Order-sens.], [1-pass concept tree; XOR & order bias],
      [\#86 CLASSIT @gennari1989classit], [DC], [Cont. stream], [$O(n B^2 d h)$], [Acuity $sigma_"ac"$], [*Yes* (Tree)], [Scale-sens.], [Streaming numeric tree; $s^(-1)$ hijack],
      [\#87 AutoClass @cheeseman1996autoclass], [BN], [*Mixed + MAR*], [$O(R T n K(d_c + d_m V))$], [*Bnd.* (MAP)], [*Yes* (CS)], [Med (MAP)], [Mixed tabular + MAR; diag. over-split],
      [\#88 TMM @peel2000robust], [PM], [Cont. $bb(R)^d$], [$O(T n K d^2 + T K d^3)$], [Unb. / Floor], [No], [*Redesc.*], [Bounded outlier influence; sym. tails],
      [\#89 vMF @banerjee2005vmf], [PM], [Sphere $bb(S)^(d-1)$], [$O(T N K d)$], [Unb. ($kappa -> infinity$)], [No], [*Compact*], [Intrinsic $bb(S)^(d-1)$ fit; ignores norms],
      [\#90 LCA @lazarsfeld1968latent], [DC], [Nominal tables], [$O(T K (N J + M))$], [*Bnd.* (0-lock)], [Via BIC], [Med], [Psychometric gold std.; Kruskal $d >= 3$],
      [\#91 HDP @teh2006hdp], [TM], [Grouped counts], [$O(T N K_"eff")$], [*Bnd.* (DP)], [*Yes* (CRF)], [Long text], [Auto-$K$ atom sharing; slow MCMC],
      [\#92 NMF @lee1999nmf], [TM], [Nonneg. matrix], [$O(T N D K)$], [Bnd. (KL/Frob)], [No], [Low / Med], [Parts-based; non-unique w/o separab.],
      [\#93 MCLUST @fraley2002mclust], [PM], [Cont. $bb(R)^d$], [$O(M T K(n d^2 + d^3))$], [`E**` *Bnd.*], [*Yes* (BIC)], [Low ($1/n$)], [14-model parsimony; over-splits skew],
      [\#94 Mixmod @biernacki2006mixmod], [PM], [Cont. / Cat.], [$O(M T K(n d^2 + d^3))$], [`E**` *Bnd.*], [*Yes* (ICL)], [Med (ICL)], [EM/CEM/SEM + ICL; CEM var. bias],
      [\#95 PMM], [DC], [Count $bb(Z)_(>= 0)^d$], [$O(T n K d)$], [*Bnd.* (0-lock)], [Via BIC], [Med], [Always identifiable; overdispersion],
      [\#96 BMM], [DC], [Binary ${0, 1}^d$], [$O(T n K d)$], [*Bnd.* (0-lock)], [Via BIC], [Med], [Exact binary LCA; $d = 2$ non-ident.],
      [\#97 DMM @yin2014dmm], [TM], [Short texts], [$O(T K overline(L) n)$], [*Bnd.* (Dir)], [*Soft* (GSDMM)], [*Short text*], [Linear-$L$ margin; 1 topic per doc],
      [\#98 SNMM @lin2009skewnormal], [PM], [Skewed $bb(R)^d$], [$O(T K n d^2 + T K d^3)$], [Unb. / Floor], [No], [Skew *Exact*], [Avoids skew split; unbnd. $perp lambda_k$],
      [\#99 IGMM @rasmussen2000igmm], [BN], [Cont. $bb(R)^d$], [$O(T K_"eff" n d^2)$], [*Bnd.* (NIW)], [*Yes* (CRP)], [Med (NIW)], [Exact untruncated DP; slow Gibbs],
      [\#100 CTM @blei2007ctm], [TM], [Doc--term counts], [$O(T n(K overline(L) + K^2) + T K^3)$], [*Bnd.* (Log-N)], [No], [Long text], [Positive topic corr.; $O(K^3)$ cost],
      [\#101 BTM @yan2013btm], [TM], [Short texts], [$O(T K |cal(B)|)$], [*Bnd.* (Dir)], [No], [Short text], [Corpus biterms; $L^2$ length-skew],
      [\#102 SNOB @wallace1968snob], [BN], [Mixed (diag.)], [$O(T N K d)$], [*Bnd.* (MML)], [*Yes* (MML)], [Small-$n$ safe], [Fisher-metric MML; diag. over-splits],
    )
  },
)

#subsec("A", "Strict Dominance Relations & Exact Cross-Algorithm Reductions")

Theorems 1--9 establish an interlocking lattice of *exact algebraic reductions* and *strict dominance relations* spanning Category IV:
1. *Exact Parametric Limiting Lattice:*
   - *TMM (\#88) $->$ GMM (\#81)* as $nu_k -> infinity$ (and Cauchy mixtures at $nu_k = 1$);
   - *SNMM (\#98) $->$ GMM (\#81)* as $lambda_k -> bold(0)$;
   - *GMM (\#81) $equiv$ MCLUST-`VVV` (\#93) $subset$ Mixmod (\#94)*, while MCLUST-`EII` with equal weights and $lambda -> 0^+$ recovers hard *Lloyd $k$-Means* (Alg. 1, @lloyd1982least) @banerjee2005clustering;
   - *vMF (\#89) $->$ Spherical $k$-Means* (Alg. 12, @dhillon2001concept) as $kappa_k = kappa -> infinity$;
   - *BMM (\#96) $subset$ LCA (\#90)* at $R_j = 2$, reducing to *K-Modes* (Alg. 7, @huang1998extensions) in the zero-temperature symmetric limit;
   - *Conditioned PMM (\#95) $equiv$ Multinomial Mixture* (bridging directly to DMM \#97 and single-topic pLSA \#84);
   - *pLSA (\#84) $equiv$ KL-NMF (\#92)* (Theorem 7.1), while MAP-pLSA recovers the unmarginalized mode of *LDA (\#83)* and *HDP (\#91)* converges to finite *LDA (\#83)* under fixed finite base support.
2. *Strict Practical Dominance Pairs:*
   - *MCLUST (\#93) / Mixmod (\#94) dominate bare unregularized GMM (\#81)* on virtually all continuous clustering problems: they enclose full-covariance GMM as `VVV`, eliminate likelihood singularities automatically across all seven `E**` models (Theorem 1.2), compress parameter complexity by $Theta(d)$ when $n_k$ is modest (Proposition 1.1), and supply principled BIC/ICL model selection.
   - *LDA (\#83) strictly dominates raw pLSA (\#84)* for probabilistic document modeling and out-of-sample generalization by eliminating pLSA's $Theta(D K)$ Neyman--Scott overfitting barrier (Theorem 7.2); retain pLSA/KL-NMF (\#92) only when fast transductive matrix factorization with sparse parts-based interpretability is the explicit target.
   - *DMM (\#97) and windowed BTM (\#101) strictly dominate LDA (\#83) on short texts ($L <= 25$)* by replacing LDA's logarithmic $log(1 + (L-1)/alpha)$ token coupling with linear-$L$ Chernoff document margins (DMM) or $O((D L)^(-1))$-variance corpus biterm aggregation (BTM, Theorem 9).

#subsec("B", "Mathematically Grounded Practitioner's Decision Guide")

Practitioners should navigate Algorithms 81--102 along eight quantitative operational regimes derived directly from Theorems 1--9:

1. *Clean Continuous Tabular Data ($n >= 10^2$, Unknown Covariance Geometry or $K$):* Default to *MCLUST (\#93)* or *Mixmod (\#94)* scanning `EII`--`VVV`. When the primary goal is *density estimation*, select models via *BIC*; when the primary goal is *well-separated cluster discovery* on mildly mis-specified modes, switch to Mixmod's *ICL* (Theorem 4.4) and initialize EM with *SEM* (or hierarchical agglomeration) rather than *CEM* to avoid Proposition 1(3)'s variance-truncation bias. Restrict to `E**` models whenever singular micro-clusters appear without explicit covariance regularization (Theorem 1.2).
2. *Heavy-Tailed Noise, Gross Outliers & Asymmetric Continuous Densities:* Whenever heavy tails or adversarial contamination are present, replace GMM/MCLUST with *TMM (\#88)*, whose Gamma scale weights guarantee redescending location influence and uniformly bounded scatter influence (Theorem 2.1; compare against hard Trimmed $k$-Means, Alg. 26, @cuesta1997trimmed). Conversely, whenever clusters are clean yet noticeably asymmetric (e.g., flow cytometry, financial returns, survival times), choose *SNMM (\#98)* so skewness is absorbed by $lambda_k$ instead of over-splitting a single mode into multiple Gaussians (Theorem 2.2).
3. *High-Dimensional Unit-Norm Embeddings & Directional Data ($x_i in bb(S)^(d-1)$):* Deploy *vMF (\#89)* rather than extrinsic Euclidean GMM (\#81) or MCLUST (\#93): vMF avoids Theorem 2(3)'s linear $(2 kappa) / (d - 1)$ radial covariance ill-conditioning, cuts parameter count from $Theta(K d^2)$ to $K d$, and remains well-posed at high $d$ (using tight rational approximations to $A_d^(-1)$ to prevent Bessel overflow @banerjee2005vmf).
4. *Bayesian Nonparametric Unknown-$K$ Clustering & Small-Sample Rigor:* For continuous datasets with $n >= 5 times 10^3$ where automatic cluster pruning and singularity immunity (NIW Loewner floor) are required, use *DPGMM (\#82)* setting truncation level $T >= 1 + (alpha + 1) ln(n / epsilon)$ (Theorem 3.1), mindful that CAVI underestimates posterior uncertainty on correlated parameters (Theorem 3.3). For smaller datasets ($n <= 5000$) requiring exact posterior uncertainty over $K$ and cluster memberships, deploy *IGMM (\#99)* with split--merge Gibbs moves. For small-$n$ model selection requiring reparameterization-invariant Fisher-volume coding and variance-singularity immunity, use *SNOB (\#102)*—strictly after verifying low within-cluster correlation (or whitening).
5. *Mixed Numeric--Categorical Tables with Missing Data (MAR):* *AutoClass (\#87)* is the premier mixed-type probabilistic engine in Category IV, combining analytic conjugate MAP-EM, principled Missing-At-Random marginalization without imputation bias, and Cheeseman--Stutz model selection (Theorem 4.3). Avoid independent-attribute AutoClass (and diagonal SNOB) on strongly correlated continuous features where Remark 1 guarantees diagonal over-splitting.
6. *Pure Count, Binary, Nominal & Streaming Conceptual Data:* For multivariate overdispersed count data (RNA-seq, web-event frequencies, crime counts), use *PMM (\#95)*—universally identifiable for all $d, K$ (Theorem 5.1)—adding a small Gamma prior to prevent zero-rate absorption. For binary interaction matrices or nominal surveys, use *BMM (\#96)* / *LCA (\#90)* strictly when $d >= 3$ satisfies Kruskal's identifiability bound (Theorem 5.2), recognizing that no local-independence mixture can separate pure parity/XOR clusters at $K = 2$ (Theorem 6.4). Reserve *COBWEB (\#85)* and *CLASSIT (\#86)* exclusively for strict single-pass online conceptual-hierarchy construction where incremental anytime tree updates outweigh batch optimality—always $z$-standardizing numeric attributes and calibrating $sigma_"ac"$ below the minimum meaningful cluster spread before invoking CLASSIT (Theorem 6.2).
7. *Long-Document Topic Modeling & Nonnegative Matrix Factorization ($L >= 50$):* Default to *LDA (\#83)* for scalable generative topic modeling with fixed $K$. Upgrade to *HDP (\#91)* when the corpus topic count $K$ is genuinely unknown and cross-document atom sharing is needed (Theorem 8.3); upgrade to *CTM (\#100)* when mapping positive inter-topic correlations is essential and $K <= 100$ makes $O(K^3)$ covariance updates affordable (Theorem 8.2). Use *NMF (\#92)* (KL divergence for count tables, Frobenius for dense spectrograms/images) when fast parts-based linear factorization without Dirichlet priors is preferred, enforcing separability or sparsity regularizers to resolve rotational non-uniqueness.
8. *Short-Text Clustering ($L <= 25$: Tweets, Queries, Titles, Logs):* *Never* rely on vanilla LDA (\#83), HDP (\#91), or CTM (\#100) without aggregation (Theorem 9.1). Choose *DMM (\#97)* (collapsed GSDMM) when each short document expresses a single dominant theme and fast $O(T K overline(L) D)$ hard document clustering with automatic cluster pruning is desired (Theorem 9.2). Choose *BTM (\#101)* when short documents mix multiple topics and rich corpus-level word co-occurrence topics are needed—always restricting biterms to a bounded window $|i - j| <= W$ on variable-length corpora to eliminate Theorem 9(3)'s quadratic length-skew dominance!

#subsec("C", "Open Theoretical Frontiers & Cross-Paradigm Bridges")

Our comparative synthesis of Algorithms 81--102 also illuminates four deep mathematical frontiers linking Category IV to other parts of this monograph:
1. *Unified Skew-$t$ Parsimonious Nonparametric Mixtures:* Section I showed that TMM (\#88) solves outlier breakdown while SNMM (\#98) solves skewness, and MCLUST (\#93) solves covariance parsimony while DPGMM/IGMM (\#82/\#99) solve unknown $K$. Developing scalable variational and parsimonious eigenvalue-decomposed *infinite skew-$t$* mixtures with provable non-asymptotic EM contraction rates (extending @balakrishnan2017statistical beyond symmetric Gaussians) remains an active theoretical frontier.
2. *Consistent Nonparametric Model-Order Recovery Under Misspecification:* While Theorem 3(2) shows that DP/HDP priors automatically discover clusters logarithmically without pre-setting $K$, Miller--Harrison inconsistency reveals that standard Dirichlet processes asymptotically spawn spurious tiny components under even mild emission misspecification (analogous to BIC over-splitting in Theorem 4.4). Reconciling coarsened posteriors, ICL-style separation penalties, and Mixtures of Finite Mixtures (MFM) with scalable variational inference is a central open challenge.
3. *Non-Euclidean & Neural Amortized Generative Mixtures:* Category IV's parametric exponential-family densities assume linear or spherical emission geometry and shallow latent variables. Extending EM and CAVI to deep nonlinear manifolds bridges directly to Variational Deep Embedding (VaDE, @jiang2017vade @kingma2014vae) in Part VIII, where a neural encoder amortizes the E-step while retaining a GMM prior in latent space; similarly, extending Poisson and multinomial mixtures to relational adjacency matrices connects PMM (\#95) and pLSA/NMF (\#84/\#92) directly to Stochastic Block Models and Degree-Corrected SBMs (@holland1983sbm @karrer2011dcsbm) in Category V.
4. *Optimal Short-Text Topic Inference with Embedded Contextual Priors:* Theorem 9 established the exact boundary separating bag-of-words LDA from DMM and BTM on short texts. Integrating pretrained directional embeddings on $bb(S)^(d-1)$ (vMF \#89) with biterm and single-topic hierarchies (DMM \#97, BTM \#101) to achieve minimax recovery at $L = O(1)$ bridges classical topic modeling and representation learning (Comparison Chapters I, II, XII--XIII).
