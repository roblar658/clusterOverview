#import "../template/ieee.typ": *

#alg_entry(
  id: 94,
  name: "Mixmod",
  category: "Cat. IV: Distribution, Statistical & Topic",
  year: "2006",
  paradigm: "Parsimonious Gaussian & Multinomial Mixture Framework (EM/CEM/SEM)",
  complexity: $cal(O)(M T K (n d^2 + d^3))$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Real-world statistical clustering frequently confronts heterogeneous data modalities—continuous measurements alongside nominal categorical surveys—as well as multi-modal likelihood surfaces where single-start EM stalls in poor local optima and BIC over-selects clusters whenever component densities are slightly misspecified. *Mixmod* (MIXture MODelling) @biernacki2006mixmod addresses all three challenges within one unified probabilistic framework that combines (i) 28 parsimonious Gaussian models for quantitative data and 10 parsimonious multinomial models for qualitative data, (ii) composite estimation strategies chaining EM @dempster1977em, Classification EM (CEM), and Stochastic EM (SEM), and (iii) multiple information criteria including BIC @schwarz1978bic, the Integrated Completed Likelihood (ICL), and the Normalized Entropy Criterion (NEC).

Given observations $X = {bold(x)_1, dots, bold(x)_n}$, Mixmod posits a $K$-component finite mixture density
$ p(bold(x) | Theta) = sum_(k=1)^K pi_k f_k (bold(x) | bold(alpha)_k), $
with mixing weights $bold(pi) = (pi_1, dots, pi_K) in Delta_(K-1)$ constrained either to be *equal* ($pi_k = 1 / K$, denoted $[p]$) or *free* ($pi_k > 0$, denoted $[p_k]$), yielding exact model-complexity granularity across both continuous and discrete domains.

#definition(title: "Mixmod Parsimonious Quantitative and Qualitative Model Families")[
  1. *Quantitative Gaussian Family (28 Models)*: For $bold(x)_i in bb(R)^d$, $f_k (bold(x) | bold(mu)_k, bold(Sigma)_k) = cal(N)(bold(x) | bold(mu)_k, bold(Sigma)_k)$ with spectral covariance decomposition $bold(Sigma)_k = lambda_k bold(D)_k bold(A)_k bold(D)_k^top$ across the 14 volume-shape-orientation constraints @fraley2002mclust, crossed with 2 proportion constraints ($[p]$ vs. $[p_k]$) for $14 times 2 = 28$ models.
  2. *Qualitative Multinomial Family (10 Models)*: For categorical vectors $bold(x)_i = (x_(i 1), dots, x_(i d))$ where attribute $j$ takes values in ${1, dots, m_j}$ (encoded one-hot as $x_(i j h) in {0, 1}$), conditional independence gives $f_k (bold(x)_i) = product_(j=1)^d p_(k j) (x_(i j))$. Under Mixmod's *center-and-dispersion* parameterization, each component-attribute distribution is specified by a modal category $a_(k j) in {1, dots, m_j}$ and dispersion parameter $epsilon_(k j) in (0, 1 - 1 / m_j)$:
     $ p_(k j) (h | a_(k j), epsilon_(k j)) = cases(1 - epsilon_(k j) & "if" h = a_(k j), epsilon_(k j) / (m_j - 1) & "if" h != a_(k j).) $
     Tying dispersion across $k$ and / or $j$ yields five structures ($epsilon$, $epsilon_k$, $epsilon_j$, $epsilon_(k j)$, plus fully unconstrained $alpha_(k j h)$), crossed with $[p]$ / $[p_k]$ for $5 times 2 = 10$ models @biernacki2006mixmod.
]

To pair model flexibility with reliable global optimization, Mixmod allows chaining three iterative engines: *EM* (soft expectation maximization of $ell_"obs" (Theta) = sum_(i=1)^n log p(bold(x)_i | Theta)$), *CEM* (alternating hard MAP classification $z_(i k) = bb(1){k = arg max_m tau_(i m)}$ and complete-data maximization of $ell_C (Theta, bold(z))$), and *SEM* (drawing stochastic latent labels $z_i^((t)) ~ "Categorical"(tau_(i 1)^((t)), dots, tau_(i K)^((t)))$ before updating $Theta$). Candidate models are ranked by BIC, by the entropy-penalized *ICL* criterion
$ "ICL"(cal(M), K) = "BIC"(cal(M), K) - 2 E(hat(tau)), $
where $E(hat(tau)) = -sum_(i=1)^n sum_(k=1)^K hat(tau)_(i k) log hat(tau)_(i k) >= 0$ is the posterior classification entropy, or by *NEC* ($"NEC"(K) = E(hat(tau)_K) / (ell_K - ell_1)$ for $K >= 2$, $"NEC"(1) = 1$).

#algbox(
  title: "Mixmod Composite Estimation and Model Selection",
  inputs: [Data $X$ (quantitative or qualitative), candidate model family $cal(G)$, cluster range ${1, dots, K_"max"}$, strategy sequence (e.g., Small-EM / SEM followed by EM / CEM), criterion $cal(C) in {"BIC", "ICL", "NEC"}$.],
  outputs: [Optimal model $(cal(M)^*, K^*)$, parameter estimate $hat(Theta)^*$, posterior responsibilities $hat(tau)^*$, partition $hat(bold(c))$.],
)[
  1. *Multi-Start Exploration (Small-EM or SEM)*: For each $(cal(M), K)$, draw $R$ random initializations (e.g., $R = 50$ short EM runs of $5$--$10$ steps, or a burn-in SEM Markov chain of $T_"SEM"$ draws) and retain the highest-likelihood seed $Theta^((0))$.
  2. *Core Iterative Refinement*: Starting from $Theta^((0))$, execute the primary engine until convergence:
     - *E-step*: Compute $tau_(i k)^((t)) = pi_k^((t)) f_k (bold(x)_i | bold(alpha)_k^((t))) / sum_m pi_m^((t)) f_m (bold(x)_i | bold(alpha)_m^((t)))$.
     - *Classification / Stochastic Step* (if CEM or SEM): Set hard labels $z_(i k)^((t)) = bb(1){k = arg max_m tau_(i m)^((t))}$ (CEM) or sample $bold(z)_i^((t)) ~ "Categorical"(tau_i^((t)))$ (SEM); skip for pure EM.
     - *Constrained M-step*: Update $Theta^((t+1)) in cal(M)$ maximizing expected (EM) or completed (CEM / SEM) log-likelihood using spectral Gaussian or modal-dispersion multinomial closed forms.
  3. *Criterion Scoring*: Compute $"BIC"$, posterior entropy $E(hat(tau))$, and $"ICL" = "BIC" - 2 E(hat(tau))$ (or $"NEC"$).
  4. *Selection*: Return $(cal(M)^*, K^*)$ optimizing criterion $cal(C)$ and MAP cluster labels $hat(c)_i = arg max_k hat(tau)_(i k)^*$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

#subsec("A", "Entropy Decomposition of Observed Likelihood and ICL Geometry")

A signature theoretical advance of Mixmod is its explicit separation between *density estimation* (favored by BIC) and *well-separated clustering* (favored by ICL) @biernacki2006mixmod. We prove the exact identity governing this gap.

#theorem(title: "Classification-Entropy Identity and ICL Separation Penalty", origin: "Original Proof")[
  For any finite mixture parameter $Theta$ with posterior responsibilities $tau_(i k) = p(Z_i = k | bold(x)_i, Theta)$, let $ell_"comp" (Theta, Z) = sum_(i=1)^n sum_(k=1)^K Z_(i k) log(pi_k f_k (bold(x)_i | bold(alpha)_k))$ denote the complete-data log-likelihood and let $hat(bold(z))$ denote the MAP partition $hat(z)_(i k) = bb(1){k = arg max_m tau_(i m)}$. Then:
  1. The observed log-likelihood satisfies the exact decomposition
     $ ell_"obs" (Theta) = bb(E)_tau [ell_"comp" (Theta, Z)] + E(tau), $
     where $E(tau) = -sum_(i=1)^n sum_(k=1)^K tau_(i k) log tau_(i k) >= 0$.
  2. The hard classification log-likelihood strictly under-bounds $ell_"obs"$ with exact gap
     $ ell_"obs" (Theta) - ell_C (Theta, hat(bold(z))) = sum_(i=1)^n log(1 / (max_k tau_(i k))) >= 0. $
  3. Mixmod's ICL criterion satisfies $"ICL"(cal(M), K) <= "BIC"(cal(M), K)$, with equality if and only if every posterior assignment is deterministic ($tau_(i k) in {0, 1}$ for all $i, k$).
]

#proof[
  By Bayes' rule, $pi_k f_k (bold(x)_i | bold(alpha)_k) = tau_(i k) p(bold(x)_i | Theta)$ whenever $tau_(i k) > 0$. Taking logarithms and averaging against $tau_(i k)$ over $k = 1, dots, K$ gives
  $ sum_(k=1)^K tau_(i k) log(pi_k f_k (bold(x)_i)) &= sum_(k=1)^K tau_(i k) log p(bold(x)_i | Theta) \
    &quad + sum_(k=1)^K tau_(i k) log tau_(i k) \
    &= log p(bold(x)_i | Theta) + sum_(k=1)^K tau_(i k) log tau_(i k). $
  Summing over $i = 1, dots, n$ yields $bb(E)_tau [ell_"comp" (Theta, Z)] = ell_"obs" (Theta) - E(tau)$, proving part (1). For part (2), let $k_i^* = arg max_k tau_(i k)$. Then
  $ ell_C (Theta, hat(bold(z))) &= sum_(i=1)^n log(pi_(k_i^*) f_(k_i^*) (bold(x)_i)) \
    &= sum_(i=1)^n log(tau_(i k_i^*) p(bold(x)_i | Theta)) \
    &= ell_"obs" (Theta) - sum_(i=1)^n log(1 / tau_(i k_i^*)). $
  Because $0 < tau_(i k_i^*) <= 1$, each summand $-log tau_(i k_i^*) >= 0$, establishing the exact gap. For part (3), each term $-tau_(i k) log tau_(i k) >= 0$ with equality iff $tau_(i k) in {0, 1}$; hence $E(hat(tau)) >= 0$ and $"ICL" = "BIC" - 2 E(hat(tau)) <= "BIC"$, directly penalizing overlapping mixture components whose fuzzy entropy grows linearly with the number of ambiguous boundary points.
]

#subsec("B", "Finite Convergence and Exact Asymptotic Bias of Classification EM")

Unlike soft EM, Mixmod's CEM engine terminates in finitely many steps, yet sacrifices statistical consistency whenever mixture components overlap @mclachlan2000finite. We prove both finite termination and an exact closed-form bias formula via a symmetric Gaussian counterexample.

#theorem(title: "Finite Termination and Exact Variance-Truncation Bias of CEM", origin: "Original Counterexample")[
  1. *Finite Monotonic Termination*: Every iteration of CEM satisfies $ell_C (Theta^((t+1)), bold(z)^((t+1))) >= ell_C (Theta^((t)), bold(z)^((t)))$, and CEM terminates at a stationary partition in at most $K^n$ steps.
  2. *Strict Asymptotic Inconsistency*: Let $X ~ (1 / 2) cal(N)(-mu, sigma^2) + (1 / 2) cal(N)(+mu, sigma^2)$ with $mu > 0$, $sigma > 0$, and signal-to-noise ratio $r = mu / sigma$. As $n -> infinity$, the symmetric equal-weight CEM estimators converge almost surely to
     $ hat(mu)_"CEM" = mu (2 Phi(r) - 1) + 2 sigma phi(r) > mu, $
     $ hat(sigma)_"CEM"^2 = sigma^2 - (hat(mu)_"CEM"^2 - mu^2) < sigma^2, $
     where $phi$ and $Phi$ are the standard normal PDF and CDF. In particular, as $r -> 0^+$, $hat(sigma)_"CEM"^2 -> (1 - 2 / pi) sigma^2 approx 0.3634 sigma^2$.
]

#proof[
  For part (1), fixing $Theta^((t))$, the C-step selects $bold(z)^((t+1))$ to maximize $sum_k z_(i k) log(pi_k^((t)) f_k (bold(x)_i | bold(alpha)_k^((t))))$ pointwise over $i$, so $ell_C (Theta^((t)), bold(z)^((t+1))) >= ell_C (Theta^((t)), bold(z)^((t)))$. Fixing $bold(z)^((t+1))$, the M-step maximizes $ell_C (Theta, bold(z)^((t+1)))$ over $Theta in cal(M)$, giving $ell_C (Theta^((t+1)), bold(z)^((t+1))) >= ell_C (Theta^((t)), bold(z)^((t+1)))$. If a partition ever repeated after a strict increase, $ell_C$ would exceed itself—a contradiction; since there are finitely many ($<= K^n$) hard assignments of $n$ points to $K$ labels, CEM terminates in finitely many iterations.

  For part (2), symmetry fixes the optimal Bayes boundary at $x = 0$, assigning $X > 0$ to component $2$ and $X <= 0$ to component $1$, with $bb(P)(X > 0) = 1 / 2$. By the strong law of large numbers, $hat(mu)_"CEM" -> bb(E)[X | X > 0]$. Splitting the expectation across the two components and substituting $u = (x - mu) / sigma$ and $v = (x + mu) / sigma$ yields
  $ bb(E)[X bb(1){X > 0}] &= 1 / 2 integral_(-r)^infinity (mu + sigma u) phi(u) dif u \
    &quad + 1 / 2 integral_r^infinity (-mu + sigma v) phi(v) dif v \
    &= 1 / 2 [ mu Phi(r) + sigma phi(r) - mu (1 - Phi(r)) + sigma phi(r) ]. $
  Dividing by $bb(P)(X > 0) = 1 / 2$ yields $hat(mu)_"CEM" = mu (2 Phi(r) - 1) + 2 sigma phi(r)$. Hence
  $ hat(mu)_"CEM" - mu = 2 sigma [ phi(r) - r (1 - Phi(r)) ] > 0 $
  strictly for every $r > 0$ by Mills' classical Gaussian tail inequality $1 - Phi(r) < phi(r) / r$. Next, because the mixture density is even ($p(-x) = p(x)$), the distribution of $X^2$ conditioned on $X > 0$ equals its unconditional distribution, so $bb(E)[X^2 | X > 0] = bb(E)[X^2] = mu^2 + sigma^2$. Therefore,
  $ hat(sigma)_"CEM"^2 = bb(E)[X^2 | X > 0] - hat(mu)_"CEM"^2 = sigma^2 - (hat(mu)_"CEM"^2 - mu^2) < sigma^2 $
  strictly for all $r > 0$, proving systematic variance underestimation and center outward bias. Taking $r -> 0^+$ gives $hat(mu)_"CEM" -> sigma sqrt(2 / pi)$ and $hat(sigma)_"CEM"^2 -> (1 - 2 / pi) sigma^2$.
]

#subsec("C", "Closed-Form M-Step for Mixmod's Parsimonious Categorical Models")

For qualitative clustering, Mixmod's center-and-dispersion models avoid the $cal(O)(K sum_j m_j)$ parameter explosion of unrestricted multinomial mixtures. We derive their exact global M-step.

#proposition(title: "Global Modal-Dispersion M-Step and Parameter Compression", origin: "Original Derivation")[
  Let $N_(k j h) = sum_(i=1)^n tau_(i k) x_(i j h)$ denote the weighted count of category $h in {1, dots, m_j}$ in cluster $k$ for variable $j$, with $n_k = sum_(h=1)^(m_j) N_(k j h) > 0$. Under the qualitative family $[p_k epsilon_(k j)]$, the unique global maximizer of the M-step surrogate $Q_"cat" (Theta | Theta^((t)))$ satisfies
  $ hat(a)_(k j) = arg max_(h in {1, dots, m_j}) N_(k j h), quad hat(epsilon)_(k j) = 1 - N_(k j hat(a)_(k j)) / n_k. $
  Moreover, $hat(epsilon)_(k j) in [0, 1 - 1 / m_j]$ automatically holds without active box projection. For the fully tied dispersion model $[p epsilon]$, the global M-step reduces to $hat(epsilon) = 1 - (1 / (n d)) sum_(k=1)^K sum_(j=1)^d N_(k j hat(a)_(k j))$, shrinking continuous parameter complexity from $K sum_(j=1)^d (m_j - 1) + K - 1$ down to a single scalar $epsilon$ plus $K d$ discrete modes.
]

#proof[
  Fixing $k$ and $j$, the attribute-$j$ contribution to $Q_"cat"$ for candidate mode $a$ and dispersion $epsilon in (0, 1 - 1 / m_j)$ is
  $ Q_(k j) (a, epsilon) &= N_(k j a) log(1 - epsilon) + (n_k - N_(k j a)) log(epsilon / (m_j - 1)) \
    &= N_(k j a) log(((m_j - 1)(1 - epsilon)) / epsilon) + n_k log(epsilon / (m_j - 1)). $
  Because $epsilon < 1 - 1 / m_j$ is algebraically equivalent to $(m_j - 1)(1 - epsilon) / epsilon > 1$, the logarithmic multiplier of $N_(k j a)$ is strictly positive. Hence for *every* feasible $epsilon$, maximizing over $a in {1, dots, m_j}$ selects $hat(a)_(k j) = arg max_h N_(k j h)$ independently of $epsilon$. Differentiating $Q_(k j) (hat(a)_(k j), epsilon)$ with respect to $epsilon$ yields $-N_(k j hat(a)_(k j)) / (1 - epsilon) + (n_k - N_(k j hat(a)_(k j))) / epsilon = 0$, giving the unique stationary point $hat(epsilon)_(k j) = 1 - N_(k j hat(a)_(k j)) / n_k$. Since $sum_(h=1)^(m_j) N_(k j h) = n_k$, the pigeonhole principle forces $max_h N_(k j h) >= n_k / m_j$, hence $hat(epsilon)_(k j) <= 1 - 1 / m_j$ automatically. Strict concavity follows from $partial^2 Q_(k j) / partial epsilon^2 = -N_(k j hat(a)_(k j)) / (1 - epsilon)^2 - (n_k - N_(k j hat(a)_(k j))) / epsilon^2 < 0$. When $epsilon$ is tied across all $(k, j)$, the monotonicity of the multiplier in $N_(k j a)$ still decouples the mode choice as $hat(a)_(k j) = arg max_h N_(k j h)$, and summing the first-order condition over $(k, j)$ yields $hat(epsilon) = 1 - (1 / (n d)) sum_(k, j) N_(k j hat(a)_(k j))$.
]

#subsec("D", "Stochastic EM Ergodicity and Escape From Saddle Points")

#theorem(title: "Ergodicity of SEM and Elimination of Unstable Saddle Trapping", origin: "Classical — Cited")[
  Let ${Theta^((t))}_(t >= 0)$ be the sequence generated by Mixmod's SEM engine on a compactly regularized parameter space $Theta_"reg"$ with strictly positive component densities @biernacki2006mixmod @mclachlan2000finite. Then ${Theta^((t))}$ defines a time-homogeneous ergodic Markov chain admitting a unique stationary distribution $Psi_n$ whose support concentrates in an $cal(O)(n^(-1/2))$ neighborhood of the consistent MLE. Moreover, unlike deterministic EM—where initializing $bold(mu)_1^((0)) = dots = bold(mu)_K^((0))$ and $bold(Sigma)_1^((0)) = dots = bold(Sigma)_K^((0))$ leaves the iterates trapped permanently on the symmetric saddle manifold—SEM breaks symmetry at step $t = 1$ with probability $1 - K^(1 - n)$.
]

#proof[
  Given $Theta^((t))$, the S-step samples $Z_1^((t)), dots, Z_n^((t))$ independently from $"Categorical"(tau_i^((t)))$, and the subsequent M-step is a deterministic measurable map $Theta^((t+1)) = g(X, Z^((t)))$. Because $tau_(i k)^((t)) > 0$ everywhere on $Theta_"reg"$, every nonempty partition $bold(z)$ receives strictly positive transition probability from every state $Theta^((t))$, establishing irreducibility and aperiodicity on the finite image set $g(X, {1, dots, K}^n)$ and guaranteeing geometric ergodicity to a unique stationary law $Psi_n$ @mclachlan2000finite. Near an isolated local maximizer, linearizing the S-step around $tau(Theta^*)$ adds zero-mean multinomial perturbation noise of covariance $cal(O)(1 / n)$, keeping stationary fluctuations at scale $cal(O)(n^(-1/2))$ while destabilizing every saddle point with a positive Hessian direction @biernacki2006mixmod. At a symmetric initialization $tau_(i k)^((0)) = 1 / K$, the probability that all $n$ i.i.d. draws land in the same cluster is $K (1 / K)^n = K^(1 - n)$; otherwise distinct clusters receive distinct sample subsets and escape the symmetric subspace immediately.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Mixmod's runtime adapts cleanly to data modality and inference strategy. For continuous Gaussian families, each EM / CEM / SEM iteration costs $cal(O)(n K d^2 + K d^3)$ for ellipsoidal covariances and $cal(O)(n K d)$ for diagonal or spherical families; running $R$ small-EM starts of $T_0$ steps followed by $T$ main iterations over $M <= 28$ models totals $cal(O)(M (R T_0 + T) K (n d^2 + d^3))$. For qualitative categorical data with $D = sum_(j=1)^d m_j$ total modality levels, each iteration evaluates $cal(O)(n K d)$ categorical log-probabilities and updates modal counts $N_(k j h)$ in $cal(O)(n K d + K D)$ time—strictly linear in $n, K, d$ and independent of matrix inversions.

Memory usage is $cal(O)(n K + K d^2)$ for Gaussian models and $cal(O)(n K + K D)$ for multinomial models. Because Mixmod's C++ core executes vectorized likelihood accumulation and employs parsimonious small-EM screening before running full convergence only on top seeds, its practical runtime is markedly lower than exhaustive MBHAC + full-grid EM pipelines while remaining highly effective at escaping local maxima @biernacki2006mixmod.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [Gaussian: $cal(O)(M T K (n d^2 + d^3))$; Categorical: $cal(O)(M T K (n d + D))$ with $D = sum_j m_j$.],
    [Space Complexity], [$cal(O)(n d + n K + K d^2)$ (continuous) or $cal(O)(n d + n K + K D)$ (categorical).],
    [Model Families], [28 Gaussian ($14 times 2$ proportions) + 10 multinomial ($5 times 2$ proportions); also mixed-type extensions.],
    [Inference Engines], [Composable pipelines: EM (consistent MLE), CEM (fast hard partitions), SEM (saddle-escaping chain).],
    [Selection Criteria], [BIC (density recovery), ICL (entropy-penalized cluster separation), NEC (normalized entropy).],
    [Initialization], [Small-EM multi-start or SEM burn-in; avoids $cal(O)(n^2)$ agglomerative bottleneck on large $n$.],
  ),
  caption: [Operational profile of the Mixmod statistical clustering framework.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Clinical Phenotyping and Epidemiological Cohort Stratification*: Patient records frequently combine nominal clinical indicators (genotypes, comorbidities, treatment categories) or continuous biomarkers with imbalanced subgroup prevalence. Mixmod's qualitative $[p_k epsilon_(k j)]$ models compress high-cardinality diagnosis tables while ICL ensures clinically distinct phenotypes rather than over-split Gaussian tails @biernacki2006mixmod.
2. *Psychometric Survey and Consumer Market Segmentation*: Likert-scale questionnaires and demographic choice panels contain $d = 20$--$100$ discrete categorical variables where Euclidean distances are meaningless. Center-and-dispersion multinomial models identify prototypical respondent profiles $hat(bold(a))_k$ with interpretable per-item agreement levels $1 - hat(epsilon)_(k j)$.
3. *Population Genetics and Multi-Locus Haplotype Clustering*: Multi-allelic microsatellite or SNP allele profiles across $d$ loci obey conditional independence within panmictic ancestral subpopulations (Hardy-Weinberg / linkage equilibrium), matching Mixmod's qualitative latent-class generative assumption identically.
4. *Semi-Supervised Discriminant Analysis and Industrial Quality Grading*: Because Mixmod unifies unsupervised clustering and supervised mixture discriminant analysis within the same 28 / 10 model families, partially labeled manufacturing inspection datasets can share parsimonious covariance tying across labeled and unlabeled batches.
5. *Ecological Species Assemblage and Environmental Typology*: Continuous physicochemical water or soil descriptors often exhibit mild non-Gaussian skewness that tricks BIC into adding redundant mixture components. Switching to Mixmod's ICL criterion penalizes overlapping entropy $E(hat(tau))$, recovering the true ecological community count.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Mixmod encounters four principal failure regimes: (i) *pure CEM on strongly overlapping clusters*, which incurs the strict variance shrinkage and outward mean shift proven in Theorem 2 (remedied by finishing with EM or SEM-EM); (ii) *conditional dependence within categorical clusters*, where strongly correlated pairs of survey questions violate multinomial local independence and induce spurious extra classes; (iii) *extreme high-dimensionality ($d >> n$)* for continuous ellipsoidal models, where even small-EM seeds struggle unless restricted to diagonal / spherical subsets; and (iv) *over-penalization by ICL on legitimately overlapping true subpopulations*, where BIC retains higher statistical power to detect touching Gaussian modes.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Categorical surveys & genetics], [10 parsimonious multinomial families match nominal attributes naturally.], [Strong intra-cluster item correlations: use log-linear / graphical models.],
    [Continuous overlapping clusters], [28 Gaussian models + Small-EM / SEM locate high-ARI global modes.], [Avoid standalone CEM due to variance truncation; prefer SEM-EM.],
    [Skewed continuous data], [ICL entropy penalty prevents BIC from over-splitting skewed tails.], [Severe asymmetry: use skew-normal @lin2009skewnormal or Student-$t$ @peel2000robust.],
    [Sparse unbounded count data], [Discretizing counts loses Poisson mean-variance structure.], [Use Poisson Mixture Models (Alg. 95) or LDA @blei2003lda.],
    [Non-convex spatial manifolds], [Gaussian ellipsoids cannot trace filaments or nested rings.], [Use HDBSCAN, OPTICS, or Spectral Clustering.],
  ),
  caption: [Domain suitability matrix and mitigation strategies for Mixmod.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Mixmod serves as a unifying algebraic umbrella linking MCLUST (Alg. 93), Lazarsfeld's Latent Class Analysis (LCA) @lazarsfeld1968latent, Bayesian AutoClass @cheeseman1996autoclass / SNOB @wallace1968snob, and hard partitional clustering @lloyd1982least.

#proposition(title: "Exact Embedding of MCLUST and Latent Class Analysis in Mixmod", origin: "Original Derivation")[
  1. *MCLUST Subsumption*: Restricting Mixmod's 28 quantitative models to the 14 free-proportion ($[p_k]$) models with pure EM inference and BIC scoring yields the exact objective space of MCLUST @fraley2002mclust. Including equal-proportion ($[p]$) counterparts strictly lowers parameter complexity by $K - 1$, reducing the BIC penalty by $(K - 1) log n$ whenever cluster sizes are balanced.
  2. *Latent Class Equivalence*: Mixmod's unconstrained qualitative model $[p_k alpha_(k j h)]$ is algebraically identical to Lazarsfeld's polytomous Latent Class Analysis @lazarsfeld1968latent, and when all variables are binary ($m_j = 2$), coincides identically with the $K$-component Bernoulli Mixture Model (Alg. 96).
]

#proof[
  For part (1), setting $bold(Sigma)_k = lambda_k bold(D)_k bold(A)_k bold(D)_k^top$ with unrestricted $bold(pi) in Delta_(K-1)$ yields identical likelihood $ell_"obs" (Theta)$ and identical manifold dimension $nu_(cal(M), K)$ as MCLUST across all 14 covariance structures. Imposing $pi_k = 1 / K$ removes the $(K - 1)$-dimensional simplex parameters, reducing $nu$ by $K - 1$ and increasing BIC by $(K - 1) log n - 2 (ell_"free" - ell_"eq")$, which is positive whenever true weights are uniform (as $2 (ell_"free" - ell_"eq") ->^d chi^2_(K-1) = cal(O)_p (1)$). For part (2), dropping the dispersion tie in Definition 1 lets $p_(k j) (h) = alpha_(k j h)$ range freely over $Delta_(m_j - 1)$, giving $p(bold(x)) = sum_(k=1)^K pi_k product_(j=1)^d product_(h=1)^(m_j) alpha_(k j h)^(x_(j h))$, the defining density of polytomous LCA; for $m_j = 2$, setting $mu_(k j) = alpha_(k j 1)$ and $1 - mu_(k j) = alpha_(k j 0)$ recovers the multivariate Bernoulli mixture identically.
]

On the ClusterOverview 4-dataset benchmark (Silhouette @rousseeuw1987silhouettes and ARI @hubert1985comparing), Mixmod showcases the dramatic value of multi-strategy initialization and equal-proportion parsimony:
- *Iris*: $"Sil" = 0.358$, $"ARI" = 0.950$. Mixmod achieves an extraordinary $"ARI" = 0.950$—far surpassing MCLUST ($0.567$), Poisson ($0.608$), and Bernoulli ($0.498$) mixtures! Because Iris has perfectly balanced classes ($pi_k = 1 / 3$) with overlapping ellipsoidal _Versicolor_ and _Virginica_ clouds, Mixmod's small-EM initialization and equal-proportion Gaussian search recover the true 3-species generative structure almost flawlessly. The lower Euclidean Silhouette ($0.358$ vs. MCLUST's $0.586$) illustrates a fundamental geometric truth: true overlapping elliptical classes have low Euclidean separation even when statistical label recovery is near-exact.
- *Customers*: $"Sil" = 0.328$. Competitive continuous clustering performance slightly behind MCLUST ($0.385$) and ahead of Poisson ($0.189$) and Bernoulli ($0.322$) models.
- *Multi-Shape*: $"Sil" = 0.712$ ($"ARI" = "N/A"$). Strong mode separation matching Poisson mixtures ($0.712$) and second only to MCLUST ($0.782$).
- *High-D*: $"Sil" = 0.295$ ($"ARI" = "N/A"$). Robust high-dimensional behavior enabled by parsimonious diagonal covariance selection, close to Bernoulli ($0.303$) and MCLUST ($0.321$).
- *Runtime*: Mean runtime is $12.143$ ms—$5.6 times$ faster than MCLUST ($68.368$ ms)—demonstrating the computational efficiency of Mixmod's targeted initialization and parsimonious screening while delivering the highest external clustering accuracy on Iris across the entire Category IV quartet.
