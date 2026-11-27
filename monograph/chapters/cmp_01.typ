#import "../template/ieee.typ": *

#cmp_entry(
  roman: "I",
  title: "Grand Unified Variational, Geometric, and Information-Theoretic Framework",
  scope_label: "Unifying All Eight Clustering Families (Algorithms 1-200)",
  provenance: "Original Unifying Proofs + Cited Classical Theorems",
)

#sec("I", "Architectural Unity of the 200-Algorithm Landscape")

Across Parts I–VIII of this monograph, we have formalized, proved, and analyzed *200 distinct clustering algorithms* spanning six decades of applied mathematics, multivariate statistics, spectral graph theory, computational topology, signal processing, and self-supervised deep representation learning @jain2010data, @xu2005survey. At first glance, this vast landscape appears deeply fragmented: Lloyd's Voronoi quantization (Alg. 1), DBSCAN's $epsilon$-reachability components (Alg. 31), Ward's variance-minimizing dendrogram (Alg. 59), Latent Dirichlet Allocation (Alg. 83), Leiden's modularity refinement (Alg. 108), Sparse Subspace Clustering's $ell_1$ self-expression (Alg. 150), Soft-DTW time-series barycenters (Alg. 164), and SwAV's online Sinkhorn prototype swapping (Alg. 190) seem to inhabit disjoint mathematical universes.

The central thesis of Part IX—and the formal objective of this opening comparative synthesis—is that *this apparent algorithmic proliferation collapses into six dual mathematical lenses*. Rather than 200 unrelated heuristics, every algorithm in Categories I–VIII implements a structured variational minimization, a topological functor on level sets, a spectral operator on a weighted diffusion graph, an information-theoretic compression principle, or a constrained optimal-transport coupling. Understanding how these six lenses interlock resolves three foundational questions that haunt clustering theory and engineering alike: (1) _when do two algorithms from different families optimize the exact same functional or converge to one another in an asymptotic limit?_ (proved systematically in Comparison Chapter II); (2) _which geometric, distributional, or regularization assumptions separate success from catastrophic failure across non-convex, heavy-tailed, or high-dimensional data?_ (formalized in Chapter XI and tested empirically in Chapter XIII); and (3) _how can an applied scientist navigate from observable data invariants directly to a mathematically principled algorithm family?_ (codified in Chapter XIV).

#wide_figure(caption: [Mathematical unification map across all eight clustering families.])[#image("../figures/fig_taxonomy_map.svg", width: 96%)]

Figure 1 displays the macro-architecture uniting all eight algorithmic families treated in this compendium. Each family node summarizes its native mathematical object—Bregman/Voronoi distortions in Category I (Algs. 1–30), kernel density level sets and gradient flows in Category II (Algs. 31–55), ultrametric trees and minimum spanning hierarchies in Category III (Algs. 56–80), exponential-family and nonparametric Bayesian mixtures in Category IV (Algs. 81–102), Laplacian cuts and random-walk flows in Category V (Algs. 103–132), quantized cells and sparse low-dimensional subspaces in Category VI (Algs. 133–150), fuzzy memberships, elastic warps, and hypergraph ensembles in Category VII (Algs. 151–175), and deep latent representations, streaming coresets, and stochastic metaheuristics in Category VIII (Algs. 176–200). Crucially, the dashed cross-family edges in Figure 1 are not metaphorical analogies: each edge corresponds to an exact algebraic isomorphism, a zero-temperature thermodynamic limit ($sigma -> 0^+$ or $m -> 1^+$), a metric-space functorial equivalence (MST level sets), or a variational dual identity.

To construct a single analytical scaffold capable of housing all 200 methods without distorting their native semantics, we begin by fixing a general variational-topological setting.

#definition(title: "Generalized Clustering Domain and Assignment Manifold")[
  Let $cal(X) = {x_1, dots, x_n}$ be a finite observation set residing in an ambient space $Omega$ equipped with either a dissimilarity $delta: Omega times Omega -> bb(R)_(>= 0)$, a measurable reference density $nu$, or a weighted affinity graph $G = (V, E, W)$ with $W in bb(R)_(>= 0)^(n times n)$. A *clustering operator* maps $cal(X)$ (optionally through a learned representation $Phi_omega: Omega -> cal(Z)$) into a structured partition representation $U in cal(M)$ together with model parameters $Theta in cal(T)$, where the *assignment manifold* $cal(M)$ takes one of five canonical forms:
  1. *Hard partition vertices:* $cal(M)_"hard" = {U in {0, 1}^(n times k) : U bold(1)_k = bold(1)_n}$ (possibly augmented with a noise column $C_(-1)$).
  2. *Row-stochastic probability simplex:* $cal(M)_"simp" = {U in [0, 1]^(n times k) : U bold(1)_k = bold(1)_n}$.
  3. *Doubly constrained transport polytope:* $Pi(r, c) = {P in bb(R)_(>= 0)^(n times k) : P bold(1)_k = r, P^T bold(1)_n = c}$.
  4. *Possibilistic / overlapping hypercube:* $cal(M)_"poss" = [0, 1]^(n times k)$ (or $bb(R)_(>= 0)^(n times k)$ for matrix factorizations such as NMF and BigCLAM).
  5. *Ultrametric / dendrogram cone:* $cal(U)(cal(X)) = {u: cal(X) times cal(X) -> bb(R)_(>= 0) : u(x, z) <= max(u(x, y), u(y, z))}$, equivalently nested partitions $theta: [0, infinity) -> cal(P)(cal(X))$.
]

Across Sections II–VII, we develop each of the *six master mathematical lenses* in full theoretical rigor, establishing six core unifying theorems—explicitly marking *Original Unified Proofs* versus *Classical — Cited* foundations—before uniting all six lenses in Section VIII via a parametric metric–diffusion–transport simplex (Proposition 1), a canonical six-tuple embedding theorem (Theorem 7), and three exhaustive eight-family comparative matrices (Tables I–III).

#sec("II", "Lens 1: Bregman Divergences and Dual Exponential Families")

The oldest and most pervasive unification in unsupervised learning connects geometric prototype quantization (Category I) with parametric generative mixture modeling (Category IV) and non-negative matrix/topic factorization (Algs. 84, 92). Its engine is convex Legendre duality, formalized by Banerjee, Merugu, Dhillon, and Ghosh @banerjee2005clustering.

Let $Omega subset.eq bb(R)^d$ be a convex domain with non-empty interior $"int"(Omega)$, and let $phi: Omega -> bb(R)$ be a strictly convex, continuously differentiable *Legendre-type* function on $"int"(Omega)$. The *Bregman divergence* $D_phi: Omega times "int"(Omega) -> bb(R)_(>= 0)$ is the first-order Taylor remainder of $phi$:
$ D_phi (x || y) := phi(x) - phi(y) - chevron.l nabla phi(y), x - y chevron.r. $
Strict convexity implies $D_phi (x || y) >= 0$ with equality if and only if $x = y$, yet $D_phi$ is generally neither symmetric nor subject to the triangle inequality. Instead, $D_phi$ obeys an exact *three-point identity* for all $x in Omega$ and $y, z in "int"(Omega)$:
$ D_phi (x || z) = D_phi (x || y) + D_phi (y || z) \
  + chevron.l x - y, nabla phi(y) - nabla phi(z) chevron.r. $
Its convex conjugate $phi^* (theta) := sup_(x in Omega) (chevron.l theta, x chevron.r - phi(x))$ defines the *dual natural coordinate* $theta = nabla phi(y)$ with inverse map $y = nabla phi^* (theta)$ and exact *dual divergence identity*:
$ D_phi (x || y) = D_(phi^*) (nabla phi(y) || nabla phi(x)). $

For a random vector $X$ taking values in $Omega$ with mean $mu = bb(E)[X] in "int"(Omega)$, the *Bregman information* @banerjee2005clustering is the generalized variance
$ I_phi (X) := bb(E)[D_phi (X || bb(E)[X])] = bb(E)[phi(X)] - phi(bb(E)[X]). $
The following landmark theorem establishes why arithmetic centroids dominate hard partitioning, why sufficient statistics compress streaming data without distortion loss, and why every regular exponential-family mixture is secretly a soft Bregman clustering engine.

#theorem(title: "Bregman–Exponential Family Duality & Unique Mean Optimality", origin: "Classical — Cited")[
  Let $phi: Omega -> bb(R)$ be a strictly convex differentiable Legendre-type function with convex conjugate $psi = phi^*$, and let $X$ be a random variable on $Omega$ with finite expectation $mu = bb(E)[X] in "int"(Omega)$ and finite $bb(E)[|phi(X)|]$. Then @banerjee2005clustering:
  1. *(Generalized Pythagorean Decomposition & Centroid Uniqueness):* For every $c in "int"(Omega)$,
     $ bb(E)[D_phi (X || c)] = I_phi (X) + D_phi (mu || c). $
     Consequently, $c^* = bb(E)[X]$ is the *unique* minimizer of $c |-> bb(E)[D_phi (X || c)]$. Conversely, if a smooth non-negative loss $F: Omega times "int"(Omega) -> bb(R)_(>= 0)$ with $F(x, x) = 0$ has $bb(E)[X]$ as the unique minimizer of $bb(E)[F(X, c)]$ for every finitely supported distribution on $Omega$, then $F$ must be a Bregman divergence $D_phi$ for some strictly convex $phi$.
  2. *(Exact Bijection with Regular Exponential Families):* Every regular exponential family $p_(psi, theta) (x) = exp(chevron.l theta, x chevron.r - psi(theta)) nu(x)$ with cumulant function $psi(theta)$ and expectation parameter $mu(theta) = nabla psi(theta)$ admits the exact representation
     $ p_(psi, theta) (x) = exp(-D_phi (x || mu(theta))) dot b_phi (x), $
     where $phi = psi^*$ and $b_phi (x) = exp(phi(x)) nu(x)$ is independent of $theta$. Therefore, maximizing the complete-data log-likelihood of a $k$-component exponential-family mixture equals minimizing a weighted Bregman distortion regularized by prior log-weights.
]

#proof[
  *(1) Pythagorean Decomposition:* Subtracting $D_phi (X || mu)$ from $D_phi (X || c)$ and invoking the three-point identity yields
  $ D_phi (X || c) - D_phi (X || mu) \
    = D_phi (mu || c) + chevron.l X - mu, nabla phi(mu) - nabla phi(c) chevron.r. $
  Taking expectations over $X$ and using linearity of expectation together with $bb(E)[X - mu] = bold(0)$ eliminates the inner product identically:
  $ bb(E)[D_phi (X || c)] = bb(E)[D_phi (X || mu)] + D_phi (mu || c) \
    = I_phi (X) + D_phi (mu || c). $
  Because $D_phi (mu || c) >= 0$ with equality if and only if $c = mu$, the expectation $mu$ uniquely minimizes expected Bregman distortion.

  For the converse (Banerjee et al. @banerjee2005clustering), consider any two-point distribution placing mass $1 - t$ at $x$ and $t$ at $y$ ($t in (0, 1)$), whose mean is $c_t = (1 - t) x + t y$. First-order stationarity of $c |-> (1 - t) F(x, c) + t F(y, c)$ at $c = c_t$ gives
  $ (1 - t) nabla_2 F(x, c_t) + t nabla_2 F(y, c_t) = bold(0), $
  where $nabla_2$ denotes the gradient in the second argument. Setting $y = c_t + ((1 - t) / t) (c_t - x)$ and differentiating with respect to $x$ shows that $nabla_2 F(x, c)$ must be affine in $x$ and vanish at $x = c$, hence $nabla_2 F(x, c) = H(c) (c - x)$ for a matrix-valued function $H(c)$. Integrability ($partial_i partial_j F = partial_j partial_i F$) forces $H(c)$ to be the symmetric Hessian $nabla^2 phi(c)$ of a scalar potential $phi$, and strict minimality forces $nabla^2 phi(c) succ 0$. Integrating $nabla_2 F(x, c) = nabla^2 phi(c)(c - x)$ along the ray from $c = x$ with $F(x, x) = 0$ yields $F(x, c) = phi(x) - phi(c) - chevron.l nabla phi(c), x - c chevron.r = D_phi (x || c)$.

  *(2) Exponential Family Bijection:* Since $phi = psi^*$ is the Legendre conjugate of the strictly convex cumulant $psi$, Fenchel's identity holds at the dual pair $(theta, mu(theta))$ with $mu(theta) = nabla psi(theta)$ and $theta = nabla phi(mu(theta))$:
  $ phi(mu(theta)) + psi(theta) = chevron.l theta, mu(theta) chevron.r. $
  Substituting $psi(theta) = chevron.l theta, mu(theta) chevron.r - phi(mu(theta))$ into the exponent of $p_(psi, theta) (x)$ gives
  $ chevron.l theta, x chevron.r - psi(theta) \
    = chevron.l nabla phi(mu), x - mu chevron.r + phi(mu) \
    = phi(x) - [phi(x) - phi(mu) - chevron.l nabla phi(mu), x - mu chevron.r] \
    = phi(x) - D_phi (x || mu(theta)). $
  Exponentiating both sides produces $p_(psi, theta) (x) = exp(-D_phi (x || mu(theta))) exp(phi(x)) nu(x)$, where $b_phi (x) = exp(phi(x)) nu(x)$ does not depend on $theta$. For a $k$-component mixture with component weights $pi_j$ and assignments $z_i in {1, dots, k}$, the complete-data negative log-likelihood equals
  $ -sum_(i=1)^n log (pi_(z_i) p_(psi, theta_(z_i)) (x_i)) \
    = sum_(i=1)^n D_phi (x_i || mu_(z_i)) - sum_(i=1)^n log pi_(z_i) - sum_(i=1)^n log b_phi (x_i), $
  completing the exact equivalence up to the parameter-free constant $sum_i log b_phi (x_i)$.
]

#subsec("A", "Unification Across Categories I, III, IV, VII, and VIII")

Theorem 1 instantly collapses dozens of seemingly distinct algorithms in Parts I–VIII into special choices of the convex generator $phi$:
- *Quadratic generator $phi(x) = norm(x)_2^2$:* Yields $D_phi (x || y) = norm(x - y)_2^2$, uniting Euclidean K-Means (Alg. 1 @lloyd1982least), K-Means++ (Alg. 2 @arthur2007kmeanspp), Mini-Batch / Hartigan–Wong / Elkan / Hamerly / Global K-Means (Algs. 3, 13, 14, 17, 19), LBG vector quantization (Alg. 30 @linde1980algorithm), Ward's agglomerative criterion (Alg. 59 @ward1963hierarchical), isotropic GMMs (Alg. 81), and coreset-accelerated stream engines BIRCH, BICO, CluStream, and StreamKM++ (Algs. 64, 78, 193, 195 @zhang1996birch, @fichtenberger2013bico, @aggarwal2003clustream, @ackermann2012streamkmpp).
- *Mahalanobis quadratic form $phi_Sigma (x) = x^T Sigma^(-1) x$:* Produces ellipsoidal cluster geometries unifying full-covariance GMMs and MCLUST models EII–VVV (Algs. 81, 93 @fraley2002mclust) with adaptive fuzzy Mahalanobis clusterers Gustafson–Kessel and Gath–Geva (Algs. 153–154 @gustafson1979gk, @gath1989gg).
- *Negative Shannon entropy $phi(x) = sum_(r=1)^d x_r log x_r - x_r$:* Generates the generalized Kullback–Leibler (I-) divergence $D_"KL" (x || y) = sum_r (x_r log(x_r / y_r) - x_r + y_r)$, unifying Poisson Mixtures (Alg. 95), Multinomial / Dirichlet-Multinomial Mixtures (Alg. 97), Bernoulli Mixtures (Alg. 96), pLSA (Alg. 84 @hofmann1999plsa), and KL-divergence NMF (Alg. 92 @lee2001algorithms, @ding2008nmfplsa).
- *Spherical / von Mises–Fisher cumulant:* On the unit hypersphere $bb(S)^(d-1)$, the vMF log-partition $psi_kappa (theta)$ induces cosine distortion $1 - chevron.l x, mu chevron.r$, uniting Spherical K-Means (Alg. 12 @dhillon2001concept), vMF Mixtures (Alg. 89 @banerjee2005vmf), and phase-aligned normalized cross-correlation in k-Shape (Alg. 162 @paparrizos2015kshape).

#remark(title: "Sufficient-Statistic Additivity & The Non-Bregman Frontier")[
  Because the right-hand minimizer of $sum_(i in C_j) D_phi (x_i || mu_j)$ depends *solely* on the additive first-order moment $(n_j, sum_(i in C_j) x_i)$—and the within-cluster Bregman inertia adds $sum_(i in C_j) phi(x_i)$—any Bregman objective admits lossless streaming micro-cluster summaries (the CF-triple $(n_j, "LS"_j, "SS"_j)$ of BIRCH Alg. 64 and CluStream Alg. 193). Conversely, Theorem 1(1) proves why robust non-Bregman objectives—$L_1$ K-Medians (Alg. 15), K-Medoids PAM/FasterPAM (Algs. 4, 24 @kaufman1990finding), Trimmed K-Means (Alg. 26 @cuesta1997trimmed), and heavy-tailed Student-$t$ Mixtures (Alg. 88 @peel2000robust)—cannot be summarized by simple additive arithmetic means: breaking Bregman mean optimality is the exact mathematical price required to lift the finite-sample Donoho–Huber breakdown point above $1/n$.
]

#sec("III", "Lens 2: Variational Free Energy and the Helmholtz Functional")

While Lens 1 unifies *distortion geometry*, it does not yet explain the continuous bridge connecting hard indicator assignments (Category I), probabilistic Expectation–Maximization (Category IV), fuzzy/possibilistic memberships (Category VII), competitive neural vector quantizers (Category VIII), and deep variational autoencoders (Category VIII). All five regimes are governed by a single thermodynamic variational principle: minimizing a *generalized Helmholtz free energy* over soft assignments $U$ and latent/generative parameters $Theta$.

#definition(title: "Generalized Helmholtz Clustering Free Energy")[
  Given per-point cluster discrepancy potentials $ell_Theta (x_i, j) in bb(R)$ (such as Bregman divergence $D_phi (x_i || mu_j)$, negative log-density $-log p(x_i | theta_j)$, or latent reconstruction-plus-prior loss), an inverse temperature $beta > 0$, a prior vector $pi in Delta_k$, a Tsallis-entropic index $gamma > 0$, and a structural parameter regularizer $Omega(Theta)$, the *Generalized Helmholtz Free-Energy Functional* is defined on $cal(M) times cal(T)$ by
  $ cal(F)_(beta, gamma) (U, Theta) := &sum_(i=1)^n sum_(j=1)^k u_(i j)^alpha ell_Theta (x_i, j) \
    &+ 1 / beta sum_(i=1)^n Psi_gamma (u_i || pi) + Omega(Theta), $
  where either: (a) $alpha = 1$ with convex divergence regularizer $Psi_gamma (u_i || pi) = (gamma - 1)^(-1) sum_(j=1)^k u_(i j) ((u_(i j) / pi_j)^(gamma - 1) - 1)$ (whose $gamma -> 1$ limit is the Kullback–Leibler divergence $"KL"(u_i || pi) = sum_j u_(i j) log (u_(i j) / pi_j)$), or (b) *fuzzy power coupling* with $alpha = m > 1$, $Psi_gamma equiv 0$ on $cal(M)_"simp"$ (Bezdek FCM), or $Psi_"PCM" (u_i) = sum_(j=1)^k eta_j (1 - u_(i j))^m$ on $cal(M)_"poss" = [0, 1]^(n times k)$ (Possibilistic C-Means).
]

The next theorem proves that hard partitioning, deterministic annealing, EM / variational Bayes, fuzzy/possibilistic C-Means, and deep amortized clustering all emerge as exact stationary manifolds or asymptotic boundaries of $cal(F)_(beta, gamma)$.

#theorem(title: "Grand Helmholtz Variational Unification & Monotone Block Descent", origin: "Original Unified Proof")[
  Consider alternating block minimization of $cal(F)_(beta, gamma) (U, Theta)$ over $U in cal(M)$ and $Theta in cal(T)$, assuming $Theta |-> sum_(i,j) u_(i j)^alpha ell_Theta (x_i, j) + Omega(Theta)$ admits an exact minimizer at each step. Then:
  1. *(Strict Lyapunov Descent):* Every alternating update $(U^((t)), Theta^((t))) -> (U^((t+1)), Theta^((t))) -> (U^((t+1)), Theta^((t+1)))$ satisfies
     $ cal(F)_(beta, gamma) (U^((t+1)), Theta^((t+1))) <= cal(F)_(beta, gamma) (U^((t+1)), Theta^((t))) \
       <= cal(F)_(beta, gamma) (U^((t)), Theta^((t))), $
     with strict decrease whenever the $U$-block is not already stationary under strict convexity of $u_i |-> cal(F)_(beta, gamma)$.
  2. *(Gibbs–EM–ELBO Equivalence, $alpha = 1, gamma -> 1$):* On $cal(M)_"simp"$, the unique $U$-step minimizer is the Boltzmann–Gibbs posterior
     $ u_(i j)^* = (pi_j exp(-beta ell_Theta (x_i, j))) / (sum_(r=1)^k pi_r exp(-beta ell_Theta (x_i, r))). $
     At $beta = 1$ and $ell_Theta (x_i, j) = -log p_Theta (x_i | c_i = j)$ with $Omega(Theta) = -log p(Theta)$, minimizing $cal(F)_(1, 1) (U, Theta)$ is *identically equivalent* to maximizing the variational Evidence Lower Bound ($"ELBO"$) of the marginal log-posterior:
     $ -cal(F)_(1, 1) (U, Theta) = sum_(i=1)^n log p_Theta (x_i) \
       - sum_(i=1)^n "KL"(u_i || p_Theta (c_i | x_i)) + log p(Theta). $
  3. *(Zero-Temperature Hard Limit, $beta -> infinity$):* For every $i$ with a unique minimizer $j^*(i) = arg min_j ell_Theta (x_i, j)$, the Gibbs minimizer converges exponentially fast to the hard Voronoi indicator:
     $ abs(u_(i, j^*(i))^* - 1) <= (k - 1) (max_r pi_r) / pi_(j^*(i)) exp(-beta Delta_i) -> 0, $
     where $Delta_i = min_(r != j^*(i)) (ell_Theta (x_i, r) - ell_Theta (x_i, j^*(i))) > 0$, recovering hard Bregman / Euclidean $k$-means (Cat. I).
  4. *(Fuzzy & Possibilistic Power Duality, $alpha = m > 1$):* Setting $alpha = m > 1$ with $ell_Theta (x_i, j) = D_phi (x_i || mu_j) > 0$ yields on $cal(M)_"simp"$ the Bezdek Fuzzy C-Means update
     $ u_(i j)^* = (sum_(r=1)^k ((D_phi (x_i || mu_j)) / (D_phi (x_i || mu_r)))^(1 / (m - 1)))^(-1), $
     which converges to hard partitioning as $m -> 1^+$ and to uniform $1/k$ as $m -> infinity$. Dropping the simplex constraint $sum_j u_(i j) = 1$ in favor of $Psi_"PCM" (u_i) = sum_j eta_j (1 - u_(i j))^m$ decouples clusters into the Possibilistic C-Means typicality $u_(i j)^* = (1 + ((D_phi (x_i || mu_j)) / eta_j)^(1 / (m - 1)))^(-1)$.
]

#proof[
  *(1) Lyapunov Descent:* By definition of block coordinate minimization, $U^((t+1)) in arg min_(U in cal(M)) cal(F)_(beta, gamma) (U, Theta^((t)))$ and $Theta^((t+1)) in arg min_(Theta in cal(T)) cal(F)_(beta, gamma) (U^((t+1)), Theta)$. Both inequalities follow immediately. Whenever $u_i |-> cal(F)_(beta, gamma)$ is strictly convex—which holds for $alpha = 1, gamma > 0$ (since $(partial^2 Psi_gamma) / (partial u_(i j)^2) = gamma u_(i j)^(gamma - 2) pi_j^(1 - gamma) > 0$) and for $alpha = m > 1$ with $ell_Theta > 0$ (since $(partial^2 (u_(i j)^m ell)) / (partial u_(i j)^2) = m(m - 1) u_(i j)^(m - 2) ell > 0$)—the $U$-minimizer is unique, forcing strict descent whenever $U^((t+1)) != U^((t))$.

  *(2) Gibbs–EM–ELBO Equivalence:* For $alpha = 1$ and $gamma -> 1$, L'Hôpital's rule gives
  $ lim_(gamma -> 1) (u_(i j) ((u_(i j) / pi_j)^(gamma - 1) - 1)) / (gamma - 1) = u_(i j) log (u_(i j) / pi_j). $
  Introducing Lagrange multipliers $lambda_i$ for $sum_(j=1)^k u_(i j) = 1$, stationarity of the Lagrangian for row $i$ reads:
  $ ell_Theta (x_i, j) + 1 / beta (log (u_(i j) / pi_j) + 1) + lambda_i = 0. $
  Solving for $u_(i j)$ yields $u_(i j) = pi_j exp(-beta ell_Theta (x_i, j) - 1 - beta lambda_i)$, and normalizing over $j = 1, dots, k$ produces the stated Boltzmann–Gibbs softmax formula. Substituting $beta = 1$, $ell_Theta (x_i, j) = -log p_Theta (x_i | c_i = j)$, and Bayes' rule $pi_j p_Theta (x_i | c_i = j) = p_Theta (x_i) p_Theta (c_i = j | x_i)$ into $-cal(F)_(1, 1)$ gives:
  $ -cal(F)_(1, 1) (U, Theta) &= sum_(i,j) u_(i j) log ((pi_j p_Theta (x_i | c_i = j)) / u_(i j)) - Omega(Theta) \
    &= sum_(i,j) u_(i j) log (p_Theta (x_i) dot (p_Theta (c_i = j | x_i)) / u_(i j)) - Omega(Theta) \
    &= sum_(i=1)^n log p_Theta (x_i) - sum_(i=1)^n "KL"(u_i || p_Theta (c_i | x_i)) + log p(Theta). $
  At $u_i = u_i^* = p_Theta (c_i | x_i)$, the KL gap vanishes identically, making the $Theta$-step an exact Majorization–Minimization (EM) lower-bound ascent @dempster1977em, @wu1983convergence, @bishop2006pattern.

  *(3) Zero-Temperature Hard Limit:* Write $j^* = j^*(i)$. Dividing numerator and denominator of $u_(i, j^*)^*$ by $pi_(j^*) exp(-beta ell_Theta (x_i, j^*))$ gives
  $ u_(i, j^*)^* = (1 + sum_(r != j^*) (pi_r / pi_(j^*)) exp(-beta (ell_Theta (x_i, r) - ell_Theta (x_i, j^*))))^(-1). $
  Using $(1 + S)^(-1) >= 1 - S$ for $S >= 0$ and $ell_Theta (x_i, r) - ell_Theta (x_i, j^*) >= Delta_i$ for all $r != j^*$ yields
  $ 1 - u_(i, j^*)^* <= sum_(r != j^*) (pi_r / pi_(j^*)) exp(-beta Delta_i) \
    <= (k - 1) (max_r pi_r) / pi_(j^*) exp(-beta Delta_i), $
  which decays exponentially to $0$ as $beta -> infinity$, reducing the free energy to $min_(cal(C), Theta) sum_i ell_Theta (x_i, c_i) + Omega(Theta)$.

  *(4) Fuzzy & Possibilistic Power Duality:* For $alpha = m > 1$ on $cal(M)_"simp"$, differentiating $sum_j u_(i j)^m d_(i j) + lambda_i (sum_j u_(i j) - 1)$ with $d_(i j) = D_phi (x_i || mu_j) > 0$ yields $m u_(i j)^(m - 1) d_(i j) + lambda_i = 0$, so $u_(i j) prop d_(i j)^(-1 / (m - 1))$. Enforcing $sum_r u_(i r) = 1$ gives Bezdek's formula @bezdek1981fcm, @hathaway1986fcm; as $m -> 1^+$, the exponent $1 / (m - 1) -> infinity$ concentrates weight strictly onto $arg min_j d_(i j)$. Without the simplex constraint, differentiating $u_(i j)^m d_(i j) + eta_j (1 - u_(i j))^m$ with respect to $u_(i j)$ yields $m u_(i j)^(m - 1) d_(i j) = m eta_j (1 - u_(i j))^(m - 1)$, i.e., $(1 - u_(i j)) / u_(i j) = (d_(i j) / eta_j)^(1 / (m - 1))$, giving Krishnapuram–Keller's PCM update @krishnapuram1993pcm.
]

#subsec("A", "Amortized Deep Variational Clustering and Feature-Collapse Dynamics")

Theorem 2 extends seamlessly from finite-dimensional coordinate spaces into deep nonlinear representation clustering (Category VIII, Algs. 181–185). In Variational Deep Embedding (VaDE, Alg. 184 @jiang2017vade), the discrete cluster assignment $c_i$ is coupled with a continuous latent code $z_i in bb(R)^p$ through a generative GMM prior $p_Theta (z_i | c_i = j) = cal(N)(z_i; mu_j, Sigma_j)$ and neural decoder $p_psi (x_i | z_i)$, while an amortized inference network $q_omega (z_i, c_i | x_i) = q_omega (z_i | x_i) u_(i, c_i)$ replaces per-instance variational optimization @kingma2014vae. The resulting negative ELBO is literally an amortized Helmholtz free energy:
$ cal(F)_"VaDE" (omega, psi, Theta, U) &= sum_(i=1)^n bb(E)_(q_omega (z_i | x_i)) [-log p_psi (x_i | z_i)] \
  &quad + sum_(i,j) u_(i j) "KL"(q_omega (z_i | x_i) || p_Theta (z_i | j)) \
  &quad + sum_(i=1)^n "KL"(u_i || pi). $
Analogously, DCN (Alg. 183 @yang2017dcn) corresponds to the deterministic zero-temperature limit ($beta -> infinity$, point-mass $q_omega (z_i | x_i) = delta_(f_omega (x_i))$) with autoencoder reconstruction penalty $Omega(omega, psi) = lambda sum_i norm(x_i - g_psi (f_omega (x_i)))_2^2$, while DEC (Alg. 181 @xie2016dec) drops the reconstruction regularizer $Omega$ after pretraining and minimizes $"KL"(P || Q)$ against a sharpened Student-$t$ auxiliary target—precisely why IDEC (Alg. 182 @guo2017idec) must reinstate $Omega$ to prevent manifold distortion! Even Power K-Means (Alg. 27 @xu2019power), which replaces the hard minimum distance with a power mean $M_s (d_(i 1), dots, d_(i k)) = (k^(-1) sum_j d_(i j)^s)^(1/s)$ annealed from $s -> 0^-$ to $s -> -infinity$, operates as a continuous-homotopy cousin of deterministic annealing ($beta: 0^+ -> infinity$) along the Helmholtz energy landscape. Finally, competitive neural quantizers—SOM (Alg. 176 @kohonen1982som) and Neural Gas (Alg. 179 @martinetz1993neuralgas)—introduce topological neighborhood kernels $h_sigma (j, j^*(i))$ into the $Theta$-update, acting as spatially coupled Gibbs weights across a prototype lattice.

#sec("IV", "Lens 3: Topological Level Sets, Persistent Homology, and Ultrametric Dendrograms")

Convex Bregman balls (Lens 1) and parametric Gaussian hulls (Lens 2) inevitably shatter non-convex, filamentous, crescent, or ring-shaped manifolds into artificial polyhedral fragments. A fundamentally different mathematical principle governs *density-based clustering* (Category II, Algs. 31–55), *hierarchical agglomerative/divisive trees* (Category III, Algs. 56–80), *density-grid and axis-parallel subspace discovery* (Category VI, Algs. 133–142), and *topological image/trajectory segmentation* (Category VII, Algs. 167–169): clusters are defined by *topological connectedness across nested superlevel sets or metric filtrations*.

Let $f: bb(R)^d -> bb(R)_(>= 0)$ be an unknown continuous probability density (estimated empirically by a kernel density estimator $hat(f)_h (x) = (n h^d)^(-1) sum_(i=1)^n K((x - x_i) / h)$ or a $k$-NN density proxy). Following Hartigan @hartigan1981consistency, the *population high-density clusters* at threshold $lambda > 0$ are the connected components $pi_0 (L_lambda (f))$ of the upper level set
$ L_lambda (f) := {x in bb(R)^d : f(x) >= lambda}. $
As $lambda$ descends from $sup f$ to $0$, components are born at local density maxima and merge at saddle points, tracing a canonical *cluster tree* (or *merge tree*) @chaudhuri2014consistent.

On a finite dissimilarity space $(cal(X), d)$, the discrete dual of superlevel connectivity is the *Vietoris–Rips filtration* $"VR"_epsilon (cal(X))$, whose $0$-simplices are $cal(X)$ and whose $1$-simplices are pairs ${x, y}$ with $d(x, y) <= epsilon$. The $0$-dimensional homology group $H_0 ("VR"_epsilon (cal(X)))$ tracks connected components as $epsilon$ increases from $0$ to $infinity$. Remarkably, Carlsson and Mémoli @carlsson2010characterization proved that hierarchical clustering admits a rigid category-theoretic characterization in terms of *ultrametric spaces*. The following unified theorem synthesizes the dendrogram–ultrametric isomorphism, Kruskal MST minimax paths, $H_0$ persistent homology, and density mutual-reachability level-set trees into a single self-contained proof.

#theorem(title: "Isomorphism of Ultrametrics, Single-Linkage Dendrograms & Persistent Level Sets", origin: "Original Unified Proof")[
  Let $(cal(X), d)$ be a finite metric (or symmetric zero-diagonal dissimilarity) space with $|cal(X)| = n$. Let $cal(D)(cal(X))$ denote the set of persistent dendrograms—right-continuous maps $theta: [0, infinity) -> cal(P)(cal(X))$ into partitions of $cal(X)$ such that $theta(0) = {{x} : x in cal(X)}$, $theta(t_1)$ refines $theta(t_2)$ for $t_1 <= t_2$, and $theta(t) = {cal(X)}$ for sufficiently large $t$—and let $cal(U)(cal(X))$ denote the set of ultrametrics on $cal(X)$ satisfying $u(x, z) <= max(u(x, y), u(y, z))$. Then @carlsson2010characterization, @gower1969minimum:
  1. *(Canonical Bijection):* The maps $Upsilon: cal(D)(cal(X)) -> cal(U)(cal(X))$ defined by $Upsilon(theta)(x, y) := inf {t >= 0 : x tilde.op_(theta(t)) y}$ and $Theta: cal(U)(cal(X)) -> cal(D)(cal(X))$ defined by $x tilde.op_(Theta(u)(t)) y <=> u(x, y) <= t$ are mutual inverses.
  2. *(Maximal Sub-Dominant Ultrametric = Single Linkage = MST = $H_0$ Persistence):* There exists a *unique* maximal ultrametric $u^* in cal(U)(cal(X))$ satisfying $u^* <= d$ pointwise, given simultaneously by:
     $ u^*(x, y) &= min_(P: x ~> y) max_(e in P) d(e) \
       &= max_(e in P_"MST" (x, y)) d(e) \
       &= inf {epsilon >= 0 : [x] = [y] " in " H_0 ("VR"_epsilon (cal(X)))}, $
     where $P$ ranges over all paths in the complete graph on $cal(X)$ and $P_"MST" (x, y)$ is the unique path between $x$ and $y$ in any Minimum Spanning Tree of $(cal(X), d)$. Hence Single Linkage (Alg. 56) is the *unique* sub-dominant, permutation-equivariant, scale-equivariant, excisive hierarchical functor @carlsson2010characterization.
  3. *(Mutual-Reachability Density Functor):* Fix $k >= 1$ and define the $k$-NN core distance $"core"_k (x) := d(x, x^((k)))$ and mutual-reachability dissimilarity
     $ d_("mreach", k) (x, y) := max("core"_k (x), "core"_k (y), d(x, y)). $
     Let $u_k^*$ be the maximal sub-dominant ultrametric of $(cal(X), d_("mreach", k))$. Then for every threshold $epsilon > 0$:
     (a) Every singleton ${x}$ with $"core"_k (x) > epsilon$ remains isolated in $Theta(u_k^*)(epsilon)$ (exact DBSCAN noise/border separation);
     (b) The non-singleton blocks of $Theta(u_k^*)(epsilon)$ coincide *identically* with the core-point connected components of DBSCAN (Alg. 31 @ester1996dbscan) at radius $epsilon$ and $"minPts" = k + 1$; and
     (c) Varying $epsilon$ across $(0, infinity)$ traces the exact HDBSCAN hierarchy (Alg. 32 @campello2013hdbscan), where the lifespan of each persistent $H_0$ bar corresponds to the Excess-of-Mass stability score $S(C) = sum_(x in C) (lambda_"max" (x, C) - lambda_"min" (C))$ with $lambda = 1 / epsilon$.
]

#proof[
  *(1) Dendrogram–Ultrametric Bijection:* Given $theta in cal(D)(cal(X))$, symmetry and $Upsilon(theta)(x, x) = 0$ are immediate. For any $x, y, z in cal(X)$, set $t_0 = max(Upsilon(theta)(x, y), Upsilon(theta)(y, z))$. By right-continuity and monotonicity of $theta$, both $x tilde.op_(theta(t_0)) y$ and $y tilde.op_(theta(t_0)) z$ hold; since $theta(t_0)$ is an equivalence relation, transitivity gives $x tilde.op_(theta(t_0)) z$, hence $Upsilon(theta)(x, z) <= t_0$, proving $Upsilon(theta) in cal(U)(cal(X))$. Conversely, given $u in cal(U)(cal(X))$, the strong triangle inequality directly guarantees that $u(x, y) <= t$ is transitive (hence an equivalence relation) for every $t >= 0$, monotonic in $t$, and right-continuous because $cal(X)$ is finite. Finally, $x tilde.op_(Theta(Upsilon(theta))(t)) y <=> inf {s : x tilde.op_(theta(s)) y} <= t <=> x tilde.op_(theta(t)) y$ by right-continuity, and $Upsilon(Theta(u))(x, y) = inf {t : u(x, y) <= t} = u(x, y)$, proving $Upsilon$ and $Theta$ are bijective inverses.

  *(2) Maximal Sub-Dominant Ultrametric & $H_0$ Equivalence:* Define $u_"path" (x, y) := min_(P: x ~> y) max_(e in P) d(e)$. Concatenating minimax paths from $x$ to $y$ and $y$ to $z$ yields a walk from $x$ to $z$ whose heaviest edge is at most $max(u_"path" (x, y), u_"path" (y, z))$, so $u_"path" in cal(U)(cal(X))$. Taking the single-edge path $(x, y)$ gives $u_"path" (x, y) <= d(x, y)$. Conversely, for *any* $u in cal(U)(cal(X))$ with $u <= d$ and any path $P = (v_0, v_1, dots, v_r)$ from $x = v_0$ to $y = v_r$, induction on the strong triangle inequality yields
  $ u(x, y) <= max_(0 <= s < r) u(v_s, v_(s+1)) <= max_(0 <= s < r) d(v_s, v_(s+1)). $
  Minimizing over all paths $P$ forces $u(x, y) <= u_"path" (x, y)$, proving $u^* = u_"path"$ is the unique pointwise maximum in ${u in cal(U)(cal(X)) : u <= d}$.

  Next, let $T$ be any Minimum Spanning Tree of $(cal(X), d)$ and let $e^*$ be the heaviest edge on $P_"MST" (x, y)$. Removing $e^*$ disconnects $T$ into two components $A$ and $B$ with $x in A, y in B$. By the MST cut property, every edge crossing $(A, B)$ has weight at least $d(e^*)$; since every path $P: x ~> y$ must cross $(A, B)$, we have $max_(e in P) d(e) >= d(e^*)$, while $P_"MST" (x, y)$ achieves equality. Thus $u^*(x, y) = max_(e in P_"MST" (x, y)) d(e)$, which by Kruskal's algorithm coincides identically with the Single Linkage merge height @gower1969minimum, @sibson1973slink. Finally, two $0$-cycles $[x], [y]$ become homologous in $H_0 ("VR"_epsilon (cal(X)))$ if and only if $x$ and $y$ lie in the same connected component of the $1$-skeleton $G_epsilon = (cal(X), {e : d(e) <= epsilon})$, which holds iff there exists a path $P: x ~> y$ with $max_(e in P) d(e) <= epsilon$, i.e., $u^*(x, y) <= epsilon$. Taking the infimum over $epsilon$ recovers $u^*(x, y)$ identically.

  *(3) Mutual-Reachability Density Functor:* Fix $epsilon > 0$. By definition, $d_("mreach", k) (x, y) <= epsilon$ holds if and only if *all three* inequalities $"core"_k (x) <= epsilon$, $"core"_k (y) <= epsilon$, and $d(x, y) <= epsilon$ hold simultaneously. Notice that $"core"_k (x) <= epsilon$ is equivalent to $|B_epsilon (x) inter cal(X)| >= k + 1$, which is the exact DBSCAN core-point predicate @ester1996dbscan! Therefore, an edge ${x, y}$ enters the $1$-skeleton of $"VR"_epsilon (cal(X), d_("mreach", k))$ if and only if *both* $x$ and $y$ are DBSCAN core points within distance $epsilon$. Consequently: (a) any point $x$ with $"core"_k (x) > epsilon$ has degree $0$ in this graph and remains an isolated singleton in $Theta(u_k^*)(epsilon)$; (b) the connected components of size $>= 2$ (plus any isolated core point with $k = 0$) are precisely the maximal density-connected core sets of DBSCAN; and (c) as $lambda = 1 / epsilon$ decreases, tracking connected components of $Theta(u_k^*)(1 / lambda)$ reproduces the exact HDBSCAN condensed cluster tree @campello2013hdbscan, @campello2015hierarchical, where each component's topological persistence integrates point-wise lifespan $lambda_"max" (x, C) - lambda_"min" (C)$.
]

#subsec("A", "Morse–Smale Gradient Flows, Lance–Williams Ultrametricity & Cubical Grids")

Theorem 3 covers *connectivity-based* density and linkage algorithms, but Lens 3 encompasses three further pillars spanning Categories II, III, VI, VII, and VIII:

1. *Morse–Smale Mode-Seeking Flows:* Given a smooth kernel density estimator $hat(f)_h (x)$ with radially symmetric profile $K(x) = c_(k, d) g(norm(x)_2^2)$, the Mean Shift algorithm (Alg. 34 @cheng1995mean, @comaniciu2002mean) and DENCLUE 2.0 (Alg. 43 @hinneburg2007denclue2) follow the adaptive-step gradient ascent flow $x^((t+1)) = x^((t)) + m_h (x^((t)))$, where the mean-shift vector $m_h (x)$ is positively proportional to $(nabla hat(f)_h (x)) / (hat(f)_h (x))$. Topologically, the basins of attraction of local maxima of $hat(f)_h$ form the *ascending Morse–Smale complex* of $hat(f)_h$ @ariascastro2016meanshift. Quickshift (Alg. 41 @vedaldi2008quickshift), Medoidshift (Alg. 42), and Density Peaks Clustering (DPC, Alg. 51 @rodriguez2014dpc) discretize this continuous gradient flow onto a directed parent forest linking each observation $x_i$ to its nearest higher-density neighbor $arg min_(j : hat(f)(x_j) > hat(f)(x_i)) d(x_i, x_j)$, while Watershed Clustering (Alg. 167 @vincent1991watersheds) constructs the exact dual descending Morse–Smale basins on image intensity topographies.

2. *Lance–Williams Ultrametricity & Dasgupta Cost:* Beyond Single Linkage, all classical agglomerative methods in Category III (SLINK Alg. 56, CLINK Alg. 57 @defays1977clink, UPGMA Alg. 58, Ward Alg. 59, WPGMA Alg. 60, UPGMC Alg. 61, WPGMC Alg. 62, Flexible Beta Alg. 74) obey the Lance–Williams recursion @lance1967general:
   $ d(A union B, C) &= alpha_A d(A, C) + alpha_B d(B, C) \
     &quad + beta d(A, B) + gamma abs(d(A, C) - d(B, C)). $
   A merge sequence induces a valid ultrametric (no *dendrogram inversions*, i.e., merge heights are non-decreasing) whenever the Bruynooghe reducibility condition $alpha_A + alpha_B + beta >= 1$ and $alpha_A, alpha_B >= 0, min(alpha_A, alpha_B) + gamma >= 0$ holds—satisfied strictly by Single, Complete, Average, Ward, and Minimax Linkage (Alg. 73 @bien2011minimax), but violated by Centroid (UPGMC) and Median (WPGMC) where $beta < 0$ allows merge heights to decrease. To evaluate dendrograms globally rather than greedily, Dasgupta @dasgupta2016cost introduced the tree cost $cal(L)_"Dasg" (T; W) = sum_(i < j) w_(i j) |"leaves"(T_(i or j))|$, for which UPGMA (Alg. 58) achieves a $1/3$-approximation on similarity graphs @moseley2017approximation while Gini-regularized Genie (Alg. 69 @gagolewski2016genie) prevents Single Linkage's single-bridge chaining catastrophe.

3. *Cubical Complexes & Topological Neural Graphs:* Grid and subspace density methods—STING (Alg. 133 @wang1997sting), CLIQUE (Alg. 134 @agrawal1998clique), WaveCluster (Alg. 135 @sheikholeslami1998wavecluster), SUBCLU (Alg. 137 @kailing2004subclu), and D-Stream (Alg. 197 @chen2007dstream)—replace Vietoris–Rips balls with *cubical complexes* whose dense cells satisfy the Apriori anti-monotonicity property across coordinate projections: any dense $s$-dimensional cell projects onto dense $(s - 1)$-dimensional faces. Likewise, Growing Neural Gas (Alg. 177 @fritzke1995gng) and Growing Cell Structures (Alg. 180) dynamically construct an induced Delaunay / simplicial nerve approximating the $1$-skeleton of the underlying data manifold!

#sec("V", "Lens 4: Spectral Dirichlet Energies, Cheeger Inequalities, and Random-Walk Diffusion")

When data arrive directly as a relational network $G = (V, E, W)$—or when nonlinear manifolds are captured via a local $k$-NN or Gaussian similarity kernel $W_(i j) = exp(-norm(x_i - x_j)_2^2 / (2 sigma^2))$—clustering becomes an operator-theoretic problem on the Hilbert space $ell^2 (V)$ of graph signals @chung1997spectral, @vonluxburg2007tutorial.

Let $W = W^T in bb(R)_(>= 0)^(n times n)$ be a symmetric non-negative affinity matrix with node degrees $d_i = sum_(j=1)^n W_(i j) > 0$, degree matrix $D = "diag"(d_1, dots, d_n)$, and total volume $"vol"(S) := sum_(i in S) d_i$. The three canonical *graph Laplacians* are:
$ L := D - W, quad L_"sym" := D^(-1/2) L D^(-1/2), quad L_"rw" := D^(-1) L = I_n - P, $
where $P = D^(-1) W$ is the Markov transition matrix of the natural random walk on $G$ with reversible stationary measure $pi_i = d_i / "vol"(V)$. Every Laplacian encodes the *discrete Dirichlet energy* of a function $f in bb(R)^n$:
$ cal(E)_G (f) := chevron.l f, L f chevron.r = 1 / 2 sum_(i, j = 1)^n W_(i j) (f_i - f_j)^2 >= 0. $
Because $L$ is positive semidefinite with $L bold(1)_n = bold(0)$, its eigenvalues satisfy $0 = lambda_1 (L_"sym") <= lambda_2 (L_"sym") <= dots <= lambda_n (L_"sym") <= 2$, and the multiplicity of $0$ equals the exact number of connected components of $G$.

The following foundational theorem unifies discrete graph cuts (Algs. 103–105), conductance bottlenecks (Algs. 107–118), and weighted Kernel $k$-Means (Alg. 16) via Courant–Fischer variational analysis and the discrete Cheeger inequality.

#theorem(title: "Dirichlet Cut Relaxation, Cheeger Inequality & Kernel-Cut Duality", origin: "Classical — Cited")[
  Let $G = (V, E, W)$ be a connected weighted undirected graph with $d_i > 0$ for all $i in V$. Then @chung1997spectral, @vonluxburg2007tutorial, @dhillon2004kernel:
  1. *(Multi-Way Cut Trace Representation):* For any $k$-partition $cal(C) = {C_1, dots, C_k}$, define the normalized indicator matrix $H in bb(R)^(n times k)$ with entries $H_(i j) = "vol"(C_j)^(-1/2) bb(1)[i in C_j]$. Then $H^T D H = I_k$ and the multi-way Normalized Cut @shi2000ncut satisfies the exact identity
     $ "Ncut"(cal(C)) := sum_(j=1)^k ("cut"(C_j, overline(C)_j)) / ("vol"(C_j)) = "Tr"(H^T L H). $
     Setting $Y = D^(1/2) H$ transforms the continuous relaxation $min_(Y^T Y = I_k) "Tr"(Y^T L_"sym" Y)$ into the bottom-$k$ eigenspace of $L_"sym"$ (NJW Alg. 103 @ng2002spectral), while replacing $D$ by $I_n$ gives the Ratio Cut identity $"Rcut"(cal(C)) = "Tr"(H_"rc"^T L H_"rc")$ (Alg. 105 @hagen1992ratiocut).
  2. *(Discrete Cheeger Conductance Inequality):* Define the Cheeger conductance constant
     $ h_G := min_(emptyset != S subset V : "vol"(S) <= "vol"(V) / 2) ("cut"(S, overline(S))) / ("vol"(S)). $
     Then the Fiedler eigenvalue $lambda_2 (L_"sym")$ satisfies the two-sided sandwich bound
     $ h_G^2 / 2 <= lambda_2 (L_"sym") <= 2 h_G. $
  3. *(Exact Weighted Kernel K-Means Isomorphism):* Let $K = sigma D + W$ with $sigma >= 0$ chosen so $K succ.eq 0$, and let $phi_K$ be the RKHS feature map satisfying $chevron.l phi_K (i), phi_K (j) chevron.r = (D^(-1) K D^(-1))_(i j)$. Then minimizing the weighted Kernel $k$-Means objective (Alg. 16) with point weights $w_i = d_i$ over all $k$-partitions $cal(C)$ equals minimizing $"Ncut"(cal(C))$ up to an additive partition-independent constant:
     $ J_"WKKM" (cal(C)) = "Tr"(D^(-1/2) K D^(-1/2)) - k sigma - k + "Ncut"(cal(C)). $
]

#proof[
  *(1) Multi-Way Trace Identity:* Let $h_j$ denote column $j$ of $H$. Since $C_1, dots, C_k$ are disjoint, $h_j^T D h_ell = delta_(j ell) "vol"(C_j)^(-1) sum_(i in C_j) d_i = delta_(j ell)$, so $H^T D H = I_k$. Evaluating the quadratic form column-by-column gives
  $ h_j^T L h_j &= 1 / 2 sum_(a, b = 1)^n W_(a b) (h_(a j) - h_(b j))^2 \
    &= 1 / 2 (sum_(a in C_j, b in.not C_j) + sum_(a in.not C_j, b in C_j)) W_(a b) / ("vol"(C_j)) \
    &= ("cut"(C_j, overline(C)_j)) / ("vol"(C_j)). $
  Summing over $j = 1, dots, k$ proves $"Ncut"(cal(C)) = "Tr"(H^T L H)$. Substituting $Y = D^(1/2) H$ gives $Y^T Y = I_k$ and $"Tr"(H^T L H) = "Tr"(Y^T L_"sym" Y)$; by the Ky Fan minimum principle, its continuous relaxation over all orthonormal $Y in bb(R)^(n times k)$ equals $sum_(ell=1)^k lambda_ell (L_"sym")$, attained by the $k$ smallest eigenvectors of $L_"sym"$. Replacing $D$ with $I_n$ and $"vol"(C_j)$ with $|C_j|$ yields the Ratio Cut identity verbatim.

  *(2) Cheeger Conductance Inequality:* By the Courant–Fischer variational characterization,
  $ lambda_2 (L_"sym") = min_(f != bold(0) : chevron.l f, D bold(1)_n chevron.r = 0) (chevron.l f, L f chevron.r) / (chevron.l f, D f chevron.r). $
  *Upper bound ($lambda_2 <= 2 h_G$):* Let $S subset V$ achieve $h_G$ with $"vol"(S) <= "vol"(V) / 2$, and set $overline(S) = V backslash S$. Construct the $D$-orthogonal test vector $f_i = "vol"(overline(S))$ for $i in S$ and $f_i = -"vol"(S)$ for $i in overline(S)$. Direct evaluation gives $chevron.l f, D bold(1)_n chevron.r = "vol"(S) "vol"(overline(S)) - "vol"(overline(S)) "vol"(S) = 0$, while
  $ chevron.l f, L f chevron.r &= "cut"(S, overline(S)) ("vol"(S) + "vol"(overline(S)))^2 \
    &= "cut"(S, overline(S)) "vol"(V)^2, \
    chevron.l f, D f chevron.r &= "vol"(S) "vol"(overline(S))^2 + "vol"(overline(S)) "vol"(S)^2 \
    &= "vol"(S) "vol"(overline(S)) "vol"(V). $
  Since $"vol"(overline(S)) >= "vol"(V) / 2$, the Rayleigh quotient satisfies
  $ lambda_2 (L_"sym") <= ("cut"(S, overline(S)) "vol"(V)) / ("vol"(S) "vol"(overline(S))) <= 2 ("cut"(S, overline(S))) / ("vol"(S)) = 2 h_G. $

  *Lower bound ($h_G^2 / 2 <= lambda_2$):* Let $f$ be a Fiedler eigenvector of $L_"rw"$ ($L f = lambda_2 D f$, $f perp D bold(1)_n$). Shift $f$ by its $D$-weighted median $c_0$ so that $g := f - c_0 bold(1)_n$ satisfies $"vol"({i : g_i > 0}) <= "vol"(V) / 2$ and $"vol"({i : g_i < 0}) <= "vol"(V) / 2$. Split $g = g_+ - g_-$ with $g_+ = max(g, 0)$ and $g_- = max(-g, 0)$. Without loss of generality assume $g_+ != bold(0)$ and its Rayleigh quotient satisfies $R(g_+) <= R(g_-)$. Because $L g = L f = lambda_2 D f$ and $g_(+, i) g_(-, i) = 0$, note that $chevron.l g_+, L g_+ chevron.r <= chevron.l g_+, L g chevron.r$ (dropping negative cross-terms $-W_(i j) g_(+, i) g_(-, j)$ strictly reduces the quadratic form) while $chevron.l g, D f chevron.r = chevron.l f, D f chevron.r <= chevron.l g, D g chevron.r$; in fact checking $(chevron.l g_+, L g_+ chevron.r) / (chevron.l g_+, D g_+ chevron.r) <= lambda_2$ directly via $L g_+ <= L g$ on $"supp"(g_+)$ and multiplying by $g_+ >= 0$:
  $ chevron.l g_+, L g_+ chevron.r &<= chevron.l g_+, L f chevron.r = lambda_2 chevron.l g_+, D f chevron.r \
    &= lambda_2 chevron.l g_+, D (g_+ + c_0 bold(1)_n) chevron.r \
    &= lambda_2 chevron.l g_+, D g_+ chevron.r + lambda_2 c_0 chevron.l g_+, D bold(1)_n chevron.r. $
  Choosing the sign of $f$ so that $c_0 <= 0$ (if $c_0 > 0$, negate $f$ first, which swaps $g_+$ and $g_-$ and flips $c_0$ to $-c_0 < 0$) ensures $c_0 chevron.l g_+, D bold(1)_n chevron.r <= 0$, hence $chevron.l g_+, L g_+ chevron.r <= lambda_2 chevron.l g_+, D g_+ chevron.r$.

  Now relabel vertices so $g_(+, 1) >= g_(+, 2) >= dots >= g_(+, n) >= 0$, normalize so $max_i g_(+, i) = 1$, and define threshold sets $S_t := {i : g_(+, i)^2 >= t}$ for $t in (0, 1]$. By construction, $"vol"(S_t) <= "vol"({i : g_i > 0}) <= "vol"(V) / 2$, so $"cut"(S_t, overline(S)_t) >= h_G "vol"(S_t)$ for every $t > 0$. Integrating over $t in [0, 1]$ via the discrete co-area (layer-cake) identity yields
  $ sum_(i < j) W_(i j) abs(g_(+, i)^2 - g_(+, j)^2) &= integral_0^1 "cut"(S_t, overline(S)_t) dif t \
    &>= h_G integral_0^1 "vol"(S_t) dif t = h_G sum_(i=1)^n d_i g_(+, i)^2 \
    &= h_G chevron.l g_+, D g_+ chevron.r. $
  Applying Cauchy–Schwarz to the left-hand side using $|a^2 - b^2| = |a - b| (a + b)$ and $(a + b)^2 <= 2(a^2 + b^2)$ gives:
  $ (sum_(i < j) W_(i j) abs(g_(+, i)^2 - g_(+, j)^2))^2 \
    <= (sum_(i < j) W_(i j) (g_(+, i) - g_(+, j))^2) (sum_(i < j) W_(i j) (g_(+, i) + g_(+, j))^2) \
    <= chevron.l g_+, L g_+ chevron.r dot 2 sum_(i < j) W_(i j) (g_(+, i)^2 + g_(+, j)^2) \
    = 2 chevron.l g_+, L g_+ chevron.r chevron.l g_+, D g_+ chevron.r. $
  Combining both bounds yields $h_G^2 chevron.l g_+, D g_+ chevron.r^2 <= 2 chevron.l g_+, L g_+ chevron.r chevron.l g_+, D g_+ chevron.r$, and dividing by $2 chevron.l g_+, D g_+ chevron.r^2$ gives $h_G^2 / 2 <= (chevron.l g_+, L g_+ chevron.r) / (chevron.l g_+, D g_+ chevron.r) <= lambda_2 (L_"sym")$.

  *(3) Weighted Kernel K-Means Isomorphism:* Let $tilde(K) = D^(-1) K D^(-1)$ be the Gram matrix of $phi_K$. The weighted centroid of cluster $C_j$ is $mu_j = "vol"(C_j)^(-1) sum_(i in C_j) d_i phi_K (i)$. Expanding the weighted within-cluster inertia gives @dhillon2004kernel:
  $ J_"WKKM" (cal(C)) &= sum_(j=1)^k sum_(i in C_j) d_i norm(phi_K (i) - mu_j)^2 \
    &= sum_(i=1)^n d_i tilde(K)_(i i) - sum_(j=1)^k 1 / ("vol"(C_j)) sum_(a, b in C_j) d_a d_b tilde(K)_(a b). $
  Using $d_a d_b tilde(K)_(a b) = K_(a b)$ and $H_(a j) = "vol"(C_j)^(-1/2) bb(1)[a in C_j]$ with $Y = D^(1/2) H$, the second term equals $"Tr"(H^T K H) = "Tr"(Y^T D^(-1/2) K D^(-1/2) Y)$. Substituting $K = sigma D + W = (sigma + 1) D - L$ gives
  $ "Tr"(H^T K H) &= (sigma + 1) "Tr"(H^T D H) - "Tr"(H^T L H) \
    &= k (sigma + 1) - "Ncut"(cal(C)). $
  Since $sum_i d_i tilde(K)_(i i) = "Tr"(D^(-1/2) K D^(-1/2))$, substituting back yields $J_"WKKM" (cal(C)) = "Tr"(D^(-1/2) K D^(-1/2)) - k sigma - k + "Ncut"(cal(C))$, establishing exact algebraic equivalence.
]

#subsec("A", "Diffusion Geometry, Modularity Spectra & Self-Expressive Subspaces")

Lens 4 reaches far beyond classical spectral cuts (Algs. 103–105), unifying four seemingly distinct families of algorithms:
- *Multi-Scale Random-Walk Diffusion:* The $t$-step transition probabilities $P^t$ induce the *diffusion distance* between vertices $i$ and $j$:
  $ D_t^2 (i, j) &:= sum_(r=1)^n (P_(i r)^t - P_(j r)^t)^2 / pi_r \
    &= sum_(ell=2)^n (1 - lambda_ell (L_"sym"))^(2 t) (psi_ell (i) - psi_ell (j))^2, $
  where $psi_ell$ are the eigenvectors of $L_"rw"$. Walktrap (Alg. 114 @pons2006walktrap) performs Ward agglomerative linkage directly on $D_t$, Power Iteration Clustering (PIC, Alg. 132 @lin2010pic) approximates a linear combination of slowly decaying diffusion modes in $O(t |E|)$ time without eigendecomposition, and Markov Clustering (MCL, Alg. 109 @vandongen2000mcl) alternates linear diffusion ($P -> P^2$) with entrywise Hadamard inflation ($P_(i j) -> P_(i j)^r / sum_a P_(a j)^r$) to drive the walk toward absorbing community attractors.
- *Modularity, Spinglass & Stochastic Block Models:* Shifting the adjacency matrix by the degree-null configuration model yields Newman's modularity matrix $B = W - (2 m)^(-1) d d^T$, where $2 m = "vol"(V)$. Maximizing modularity $Q(cal(C)) = (2 m)^(-1) "Tr"(H_0^T B H_0)$—via Louvain (Alg. 107 @blondel2008louvain), Leiden (Alg. 108 @traag2019leiden), Leading Eigenvector (Alg. 116 @newman2006leading), or Potts Spinglass Hamiltonians (Alg. 117 @reichardt2006spinglass)—is a spectral-quadratic optimization dual to Ncut whose resolution scale is controlled by the rank-1 null-model penalty (causing the Fortunato–Barthélemy resolution limit @fortunato2007resolution, and coinciding with the profile log-likelihood of a uniform-planting Degree-Corrected SBM, Alg. 131 @karrer2011dcsbm, @abbe2018sbmreview).
- *Sparse Subspace Clustering & Co-Regularized Views:* In high-dimensional union-of-subspaces data (Category VI), Sparse Subspace Clustering (SSC, Alg. 150 @elhamifar2013ssc, @soltanolkotabi2012geometric) solves $min_(c_i) norm(c_i)_1$ subject to $x_i = X c_i, c_(i i) = 0$. When subspaces are sufficiently separated, the *subspace-detection property* guarantees $c_(i j) = 0$ whenever $x_i$ and $x_j$ belong to different subspaces, making the symmetrized affinity $W = |C| + |C|^T$ block-diagonal so that the nullspace of $L_"sym"$ recovers the exact subspace partition! Likewise, Spectral Biclustering (Alg. 123 @kluger2003spectral) solves the bipartite SVD of $D_r^(-1/2) A D_c^(-1/2)$, while Co-regularized Multi-View Spectral Clustering (Alg. 175 @kumar2011coregularized) couples per-view Dirichlet energies via a cross-view subspace projection penalty $-lambda "Tr"(U^((v)) U^((v) T) U^((w)) U^((w) T))$.

#sec("VI", "Lens 5: Minimum Description Length, Rate-Distortion, and Information Bottleneck")

All four lenses above focus on geometric distortion, thermodynamic energy, topological connectivity, or spectral cuts—yet none intrinsically answers *how many clusters $k$ are statistically warranted*, nor how to cluster without a metric when data consist solely of symbolic co-occurrences, trajectories, or paired augmentations. Lens 5 resolves these challenges by casting clustering as *lossy or lossless information compression* governed by Shannon entropy, Kolmogorov/MDL coding, Rate-Distortion theory, and Mutual Information maximization.

Under the *two-part Minimum Description Length (MDL)* principle, the optimal cluster model $cal(M)$ minimizes the total bit-length required to transmit both the model parameters and the residual data encoded under that model:
$ L_"MDL" (cal(M), cal(X)) := L(cal(M)) + L(cal(X) | cal(M)). $
By Kraft's inequality, every prefix-free code corresponds to a probability distribution with ideal codelength $-log_2 p(dot)$, converting statistical model selection into exact bit-budget accounting. In parallel, *Shannon Rate-Distortion theory* quantifies the minimal mutual information (rate) $R(D)$ needed to compress a source $X$ into discrete cluster codewords $C in {1, dots, k}$ at expected distortion at most $D$:
$ R(D) := min_(q(c | x) : bb(E)[d(X, mu_C)] <= D) I_q (X; C). $
Finally, Tishby's *Information Bottleneck* and *self-supervised mutual-information clustering* maximize $I(C(X); C(X'))$ across stochastically perturbed views $(X, X')$ or consecutive random-walk states $(X_t, X_(t+1))$ while compressing nuisance variation.

The following unifying theorem establishes three exact analytical bridges tying together X-Means (Alg. 10), MCLUST (Alg. 93), SNOB (Alg. 102), LBG vector quantization (Alg. 30), Infomap (Alg. 113), ENCLUS (Alg. 143), Consensus Ensembles (Algs. 171–173), IIC (Alg. 188), and Contrastive Clustering (Alg. 192).

#theorem(title: "Unified MDL–BIC–Rate-Distortion–Map Equation Information Duality", origin: "Original Unified Proof")[
  The following three information-theoretic dualities hold across probabilistic, prototype, graph, and deep clustering regimes:
  1. *(Laplace–MDL–BIC Asymptotic Equivalence):* Let $cal(M)_k$ be a regular parametric mixture family with parameter vector $theta in bb(R)^(p_k)$, smooth strictly positive prior $pi(theta)$, log-likelihood $ell_n (theta) = sum_(i=1)^n log p(x_i | theta)$, unique interior MLE $hat(theta)_n$, and observed Fisher information $hat(J)_n := -n^(-1) nabla^2 ell_n (hat(theta)_n) succ 0$. Then the exact Bayesian / two-part MDL marginal codelength satisfies @schwarz1978bic, @wallace1968snob:
     $ -log p(cal(X) | cal(M)_k) &= -ell_n (hat(theta)_n) + (p_k / 2) log n \
       &quad + 1 / 2 log det hat(J)_n - log pi(hat(theta)_n) \
       &quad - (p_k / 2) log (2 pi) + o(1), $
     recovering the Bayesian Information Criterion $"BIC"(cal(M)_k) = -2 ell_n (hat(theta)_n) + p_k log n$ (X-Means Alg. 10 @pelleg2000xmeans, MCLUST Alg. 93 @fraley2002mclust) at leading $O(log n)$ order and Wallace–Boulton's Minimum Message Length (SNOB Alg. 102 @wallace1968snob) when retaining the Fisher-determinant curvature term $1/2 log det hat(J)_n$.
  2. *(Rate-Distortion = Bregman Helmholtz Variational Equivalence):* For any Bregman divergence $D_phi$, the Lagrangian of the empirical Shannon rate-distortion problem $min_(q, mu) I_q (X; C) + beta bb(E)_q [D_phi (X || mu_C)]$ admits the exact variational upper-bound representation
     $ min_(q, pi, mu) quad 1 / n sum_(i=1)^n ("KL"(q_i || pi) + beta sum_(j=1)^k q_(i j) D_phi (x_i || mu_j)), $
     whose alternating Blahut–Arimoto updates coincide *identically* with soft Bregman EM (Theorem 2 at inverse temperature $beta$). Every stationary point traces the Pareto frontier of the rate-distortion curve $R(D)$ with tangent slope $(dif R) / (dif D) = -beta$, converging to LBG / K-Means vector quantization (Algs. 1, 30) as $beta -> infinity$.
  3. *(Map Equation Compression–Mutual Information Identity):* Let $(X_t)$ be a stationary random walk on a strongly connected graph $G$ with stationary distribution $p_i$, partitioned into modules $cal(M) = {C_1, dots, C_k}$ with module exit flows $q_(j arrow.tr) = sum_(i in C_j, r in.not C_j) p_i P_(i r)$, total inter-module transition rate $q_arrow.tr = sum_(j=1)^k q_(j arrow.tr)$, and module codebook weights $p_j^circle.small = q_(j arrow.tr) + sum_(i in C_j) p_i$. Then the Infomap Map Equation codelength @rosvall2008infomap
     $ L(cal(M)) := q_arrow.tr H(cal(Q)) + sum_(j=1)^k p_j^circle.small H(cal(P)_j) $
     admits the exact decomposition (in bits, with $C_t := C(X_t)$):
     $ L(cal(M)) = H(X_t) - I(X_t; C_t) + q_arrow.tr H(C_(t+1) | C_t != C_(t+1)) \
       + sum_(j=1)^k p_j^circle.small H_2 (q_(j arrow.tr) / p_j^circle.small), $
     where $H_2 (p) = -p log_2 p - (1 - p) log_2 (1 - p)$ is binary entropy. Consequently, minimizing $L(cal(M))$ simultaneously *maximizes* state-module mutual information $I(X_t; C_t)$ while *minimizing* cross-module transition entropy $q_arrow.tr$, yielding a strictly positive compression gain over the unpartitioned code $H(X_t)$ whenever modular confinement outweighs the index-codebook switching cost.
]

#proof[
  *(1) Laplace–MDL–BIC Expansion:* The marginal likelihood under $cal(M)_k$ is $p(cal(X) | cal(M)_k) = integral exp(ell_n (theta)) pi(theta) dif theta$. Performing a second-order Taylor expansion of $ell_n (theta)$ about the interior MLE $hat(theta)_n$ (where $nabla ell_n (hat(theta)_n) = bold(0)$) gives
  $ ell_n (theta) = ell_n (hat(theta)_n) - n / 2 (theta - hat(theta)_n)^T hat(J)_n (theta - hat(theta)_n) + R_n (theta), $
  with remainder $R_n (theta) = O(n norm(theta - hat(theta)_n)_2^3)$. Under standard regularity, the integral concentrates on a $O(n^(-1/2))$-neighborhood of $hat(theta)_n$, yielding the multidimensional Gaussian Laplace integral:
  $ p(cal(X) | cal(M)_k) &= exp(ell_n (hat(theta)_n)) pi(hat(theta)_n) \
    &quad times (2 pi / n)^(p_k / 2) (det hat(J)_n)^(-1/2) (1 + o(1)). $
  Taking $-log$ of both sides produces the exact asymptotic codelength expansion in Claim 1. In SNOB (Alg. 102 @wallace1968snob), the Fisher volume $(det hat(J)_n)^(1/2)$ specifies the optimal parameter quantization cell width $delta_theta prop (det(n hat(J)_n))^(-1/2)$, making MML invariant to smooth parameter reparameterizations, whereas X-Means (Alg. 10) and MCLUST (Alg. 93) truncate at the $O(log n)$ BIC penalty $(p_k / 2) log n$.

  *(2) Rate-Distortion = Bregman Helmholtz Equivalence:* Let $q(c | x_i) = q_(i c)$ denote the stochastic encoder from the uniform empirical distribution $hat(p)(x_i) = 1 / n$ to cluster index $c in {1, dots, k}$, with induced marginal $overline(q)_c := n^(-1) sum_(i=1)^n q_(i c)$. By definition of mutual information,
  $ I_q (X; C) = 1 / n sum_(i=1)^n sum_(j=1)^k q_(i j) log (q_(i j) / overline(q)_j) = 1 / n sum_(i=1)^n "KL"(q_i || overline(q)). $
  For any auxiliary distribution $pi in Delta_k$, note the identity
  $ 1 / n sum_(i=1)^n "KL"(q_i || pi) - I_q (X; C) = sum_(j=1)^k overline(q)_j log (overline(q)_j / pi_j) = "KL"(overline(q) || pi) >= 0, $
  with equality *if and only if* $pi = overline(q)$. Therefore,
  $ I_q (X; C) = min_(pi in Delta_k) 1 / n sum_(i=1)^n "KL"(q_i || pi). $
  Adding $beta bb(E)_q [D_phi (X || mu_C)] = (beta / n) sum_(i,j) q_(i j) D_phi (x_i || mu_j)$ produces the exact Bregman Helmholtz free energy of Theorem 2 (divided by $n beta$). Minimizing over $pi$ sets $pi_j = overline(q)_j$; minimizing over $mu_j$ yields the unique arithmetic mean $mu_j = (sum_i q_(i j) x_i) / (sum_i q_(i j))$ by Theorem 1(1); and minimizing over $q_i$ yields $q_(i j) prop pi_j exp(-beta D_phi (x_i || mu_j))$ by Theorem 2(2)—precisely the Blahut–Arimoto rate-distortion recursion! Since $R(D)$ is convex and non-increasing, supporting hyperplanes have slope $(dif R) / (dif D) = -beta$, and letting $beta -> infinity$ drives $q_(i j)$ onto deterministic Voronoi vertices, recovering Lloyd / LBG quantization @linde1980algorithm.

  *(3) Map Equation–Mutual Information Identity:* Let $p_(j) := sum_(i in C_j) p_i$ be the stationary probability of module $C_j$, so $p_j^circle.small = p_j + q_(j arrow.tr)$. Expand the module codebook entropy $p_j^circle.small H(cal(P)_j)$ by separating the exit symbol (probability $q_(j arrow.tr) / p_j^circle.small$) from the within-module node symbols (probabilities $p_i / p_j^circle.small$ for $i in C_j$):
  $ p_j^circle.small H(cal(P)_j) &= -q_(j arrow.tr) log_2 (q_(j arrow.tr) / p_j^circle.small) - sum_(i in C_j) p_i log_2 (p_i / p_j^circle.small). $
  Insert $p_j / p_j$ inside the second logarithm:
  $ -sum_(i in C_j) p_i log_2 (p_i / p_j^circle.small) = -sum_(i in C_j) p_i log_2 (p_i / p_j) - p_j log_2 (p_j / p_j^circle.small). $
  Recognizing $-q_(j arrow.tr) log_2 (q_(j arrow.tr) / p_j^circle.small) - p_j log_2 (p_j / p_j^circle.small) = p_j^circle.small H_2 (q_(j arrow.tr) / p_j^circle.small)$ and $-sum_(i in C_j) p_i log_2 (p_i / p_j) = p_j H(X_t | C_t = j)$, summing over $j = 1, dots, k$ yields
  $ sum_(j=1)^k p_j^circle.small H(cal(P)_j) = H(X_t | C_t) + sum_(j=1)^k p_j^circle.small H_2 (q_(j arrow.tr) / p_j^circle.small). $
  By the fundamental mutual-information identity $H(X_t | C_t) = H(X_t) - I(X_t; C_t)$ (which equals $H(X_t) - H(C_t)$ for deterministic partitions). Meanwhile, in the index codebook $cal(Q) = {q_(j arrow.tr) / q_arrow.tr}_(j=1)^k$, by detailed balance (or flow conservation $q_(j arrow.tr) = q_(j arrow.br)$ in ergodic walks), $q_(j arrow.tr) / q_arrow.tr = bb(P)(C_(t+1) = j | C_(t+1) != C_t)$, so $q_arrow.tr H(cal(Q)) = q_arrow.tr H(C_(t+1) | C_t != C_(t+1))$. Combining both terms establishes the exact decomposition in Claim 3.
]

#subsec("A", "Self-Supervised Mutual Information, Consensus Ensembles & Conceptual Utility")

Theorem 5 illuminates why information-theoretic objectives appear identically across four otherwise distant categories:
- *Self-Supervised Deep Clustering (Category VIII):* Invariant Information Clustering (IIC, Alg. 188 @ji2019iic) maximizes $I(Phi(x); Phi(g(x)))$ directly over paired augmentations $(x, g(x))$, expanding as $H(Phi(x)) - H(Phi(x) | Phi(g(x)))$: the conditional entropy aligns augmented pairs into identical clusters while the marginal entropy $H(Phi(x))$ automatically resists degenerate collapse! Contrastive Clustering (Alg. 192 @li2021contrastive) maximizes the InfoNCE variational lower bound on mutual information @oord2018infonce simultaneously at the instance level and the cluster-column level.
- *Ensemble Consensus (Category VII):* Strehl and Ghosh @strehl2002cluster defined the *Supra-Consensus Clustering* problem as maximizing average Normalized Mutual Information $sum_(r=1)^M "NMI"(cal(C)^*, cal(C)^((r)))$ across $M$ base partitions—approximated by CSPA, HGPA, and MCLA (Algs. 171–173)—which is literally an information-theoretic Fréchet mean in the space of discrete partitions.
- *Subspace & Conceptual Entropy (Categories IV & VI):* ENCLUS (Alg. 143 @cheng1999enclus) exploits downward closure of Shannon subspace entropy $H(S) <= H(S')$ for $S subset.eq S'$ to prune uniform noise subspaces, while COBWEB (Alg. 85 @fisher1987cobweb) maximizes Category Utility, which is a quadratic (Gini-impurity) second-order Taylor approximation to the mutual information $I(C; A_1, dots, A_d)$ between cluster labels and attribute values.

#sec("VII", "Lens 6: Entropic Optimal Transport and Balanced Assignment")

A subtle yet fatal pathology plagues unconstrained variational clustering (Lenses 1–2) whenever feature representations $Phi_omega (x)$ are learned jointly with cluster assignments (Category VIII), or whenever severe cardinality imbalance creates empty or degenerate clusters. If one naively minimizes $sum_(i,j) u_(i j) norm(Phi_omega (x_i) - mu_j)_2^2$ over $(omega, U, mu)$, the encoder collapses all points onto a single constant vector $Phi_omega (x_i) equiv mu_1$, achieving zero distortion with zero semantic information (the notorious *representation collapse* of naive DeepCluster Alg. 191 @caron2018deepcluster and unregularized DEC Alg. 181 @xie2016dec).

Lens 6 eliminates this degeneracy geometrically by replacing the row-stochastic simplex $cal(M)_"simp"$ with the *Kantorovich transportation polytope*
$ Pi(r, c) := {P in bb(R)_(>= 0)^(n times k) : P bold(1)_k = r, quad P^T bold(1)_n = c}, $
where $r in Delta_n$ specifies sample weights (typically $r = n^(-1) bold(1)_n$) and $c in Delta_k$ prescribes target cluster proportions (typically equipartition $c = k^(-1) bold(1)_k$). Given a cost matrix $C in bb(R)^(n times k)$—for example, negative cosine similarity $C_(i j) = -chevron.l z_i, mu_j chevron.r$ in SwAV (Alg. 190 @caron2020swav) or squared Euclidean distance $C_(i j) = norm(x_i - mu_j)_2^2$—*Entropic Optimal Transport (EOT)* @cuturi2013sinkhorn solves the strictly convex program:
$ W_epsilon (r, c; C) := min_(P in Pi(r, c)) chevron.l P, C chevron.r - epsilon H(P), $
where $H(P) := -sum_(i=1)^n sum_(j=1)^k P_(i j) (log P_(i j) - 1)$ and $epsilon > 0$ is the entropic smoothing parameter.

The following theorem establishes the dual diagonal scaling structure, exact geometric convergence rate in the Birkhoff–Hilbert projective metric, and mathematical collapse immunity of EOT clustering.

#theorem(title: "Duality, Hilbert-Metric Contraction & Collapse Immunity of Entropic OT Clustering", origin: "Original Unified Proof")[
  Let $C in bb(R)^(n times k)$ be any finite cost matrix, let $epsilon > 0$, $r in bb(R)_(> 0)^n$ with $sum_i r_i = 1$, $c in bb(R)_(> 0)^k$ with $sum_j c_j = 1$, and define the positive Gibbs kernel $K_(i j) := exp(-C_(i j) / epsilon) > 0$. Then @cuturi2013sinkhorn, @caron2020swav:
  1. *(Unique Dual Scaling Solution):* The problem $W_epsilon (r, c; C)$ admits a *unique* strictly positive global minimizer $P^* in Pi(r, c)$ of the diagonal-scaling form
     $ P^* = "diag"(u^*) K "diag"(v^*), $
     where $(u^*, v^*) in bb(R)_(> 0)^n times bb(R)_(> 0)^k$ solve the coupled Sinkhorn–Knopp equations $u^* = r slash.big (K v^*)$ and $v^* = c slash.big (K^T u^*)$ (componentwise division), unique up to reciprocal rescaling $(alpha u^*, alpha^(-1) v^*)$ for $alpha > 0$.
  2. *(Exact Birkhoff–Hilbert Geometric Contraction):* Equip the positive orthant $bb(R)_(> 0)^n$ with the Hilbert projective metric
     $ d_H (u, u') := log max_(1 <= i, i' <= n) (u_i u'_(i')) / (u_(i') u'_i). $
     Starting from any $v^((0)) in bb(R)_(> 0)^k$, the Sinkhorn iterates $u^((t+1)) = r slash.big (K v^((t)))$, $v^((t+1)) = c slash.big (K^T u^((t+1)))$ contract geometrically to the unique projective fixed point $[u^*]$ at rate
     $ d_H (u^((t)), u^*) <= kappa(K)^(2 t) d_H (u^((0)), u^*), $
     where $kappa(K) := tanh(Delta(K) / 4) < 1$ and the projective kernel oscillation diameter is
     $ Delta(K) := max_(i, i', j, j') (C_(i j') + C_(i' j) - C_(i j) - C_(i' j')) / epsilon < infinity. $
  3. *(Exact Anti-Collapse & Equipartition Bound):* Setting $r = n^(-1) bold(1)_n$ and $c = k^(-1) bold(1)_k$ in online deep prototype clustering (SwAV Alg. 190 @caron2020swav, DEPICT Alg. 185 @ghasedi2017depict) with soft assignment matrix $Q^* = n P^*$ guarantees:
     (a) maximal cluster marginal entropy $H(n^(-1) (Q^*)^T bold(1)_n) = log k$ (zero empty-cluster degeneracy); and
     (b) if the representation network collapses to a constant $z_i equiv z_0$, then $P^* = r c^T = (n k)^(-1) bold(1)_n bold(1)_k^T$ becomes uniform, imposing a strictly positive cross-entropy penalty $cal(L)_"swap" >= log k > 0$ that prevents trivial collapse whenever prototype vectors span distinct directions.
]

#proof[
  *(1) Strict Convexity & Fenchel–Rockafellar Duality:* Since $x |-> x (log x - 1)$ has second derivative $1 / x > 0$ on $(0, infinity)$, the objective $J_epsilon (P) = chevron.l P, C chevron.r + epsilon sum_(i,j) P_(i j) (log P_(i j) - 1)$ is strictly convex on the compact non-empty polytope $Pi(r, c)$ (which contains the strictly interior point $r c^T > 0$). Moreover, the gradient $(partial J_epsilon) / (partial P_(i j)) = C_(i j) + epsilon log P_(i j) -> -infinity$ as $P_(i j) -> 0^+$, so the unique minimizer $P^*$ must lie strictly in the interior $P_(i j)^* > 0$. Introducing dual multipliers $f in bb(R)^n$ and $g in bb(R)^k$ for $P bold(1)_k = r$ and $P^T bold(1)_n = c$, the Lagrangian is
  $ cal(L)(P, f, g) = chevron.l P, C chevron.r + epsilon sum_(i,j) P_(i j) (log P_(i j) - 1) \
    - chevron.l f, P bold(1)_k - r chevron.r - chevron.l g, P^T bold(1)_n - c chevron.r. $
  First-order KKT stationarity $(partial cal(L)) / (partial P_(i j)) = 0$ gives
  $ C_(i j) + epsilon log P_(i j)^* - f_i - g_j = 0 \
    => P_(i j)^* = exp(f_i / epsilon) exp(-C_(i j) / epsilon) exp(g_j / epsilon). $
  Defining $u_i^* := exp(f_i / epsilon) > 0$ and $v_j^* := exp(g_j / epsilon) > 0$ produces $P^* = "diag"(u^*) K "diag"(v^*)$. Enforcing the primal feasibility constraints $P^* bold(1)_k = r$ and $(P^*)^T bold(1)_n = c$ yields $u^* dot.o (K v^*) = r$ and $v^* dot.o (K^T u^*) = c$. Because $P^*$ is unique and $K > 0$, any alternative pair $(tilde(u), tilde(v))$ satisfies $tilde(u)_i tilde(v)_j = u_i^* v_j^*$ for all $(i, j)$, forcing $tilde(u)_i / u_i^* = v_j^* / tilde(v)_j = alpha$ for a constant $alpha > 0$.

  *(2) Birkhoff–Hilbert Contraction Rate:* Observe first two elementary invariances of the Hilbert metric $d_H$ on $bb(R)_(> 0)^m$: (i) *positive diagonal invariance*, $d_H ("diag"(a) x, "diag"(a) y) = d_H (x, y)$ for every $a in bb(R)_(> 0)^m$; and (ii) *componentwise inversion invariance*, $d_H (x^(-1), y^(-1)) = d_H (x, y)$, since $((x_i^(-1) y_(i')^(-1)) / (x_(i')^(-1) y_i^(-1))) = (x_(i') y_i) / (x_i y_(i'))$ has identical maximum over $(i, i')$. Therefore, the Sinkhorn half-step $u^((t+1)) = r slash.big (K v^((t)))$ satisfies
  $ d_H (u^((t+1)), u^*) = d_H (r slash.big (K v^((t))), r slash.big (K v^*)) = d_H (K v^((t)), K v^*). $
  Because $K$ is strictly positive, Birkhoff's classical projective contraction theorem states that for all $x, y in bb(R)_(> 0)^k$,
  $ d_H (K x, K y) <= tanh(Delta(K) / 4) d_H (x, y), $
  where the projective diameter of $K$ is
  $ Delta(K) := sup_(x, y > 0) d_H (K x, K y) = max_(i, i', j, j') log ((K_(i j) K_(i' j')) / (K_(i j') K_(i' j))). $
  To verify Birkhoff's factor directly for any two rows $i, i'$, let $a_j = K_(i j) y_j$, $b_j = K_(i' j) y_j$, and $r_j = x_j / y_j$, with $M = max_j r_j$ and $m = min_j r_j$ so $d_H (x, y) = log(M / m)$. By quasi-convexity of the linear-fractional ratio $((sum_j a_j r_j) / (sum_j b_j r_j)) / ((sum_j a_j) / (sum_j b_j))$ over the box $[m, M]^k$, the worst-case extremal occurs at a vertex $r_j in {m, M}$, splitting indices into $S = {j : r_j = M}$ and $overline(S) = {j : r_j = m}$. Letting $A_1 = sum_(j in S) a_j, A_2 = sum_(j in overline(S)) a_j, B_1 = sum_(j in S) b_j, B_2 = sum_(j in overline(S)) b_j$, note that $(A_1 B_2) / (A_2 B_1) <= exp(Delta(K))$. Setting $rho = exp(Delta(K) / 2)$ and optimizing the two-state cross-ratio over $A_1 / A_2, B_1 / B_2$ subject to their ratio being bounded by $rho^2$ yields the sharp bound
  $ log ((M A_1 + m A_2) / (M B_1 + m B_2) dot (B_1 + B_2) / (A_1 + A_2)) \
    <= log ((rho sqrt(M / m) + 1) / (rho + sqrt(M / m))) <= tanh(Delta(K) / 4) log (M / m), $
  since $g(s) = log((rho e^(s/2) + 1) / (rho + e^(s/2)))$ is concave on $s >= 0$ with $g(0) = 0$ and $g'(0) = (rho - 1) / (2(rho + 1)) = tanh(Delta(K) / 4) / 2$ (and doubling for the reciprocal worst pair gives $tanh(Delta(K) / 4)$). Substituting $K_(i j) = exp(-C_(i j) / epsilon)$ yields $Delta(K) = max_(i, i', j, j') (C_(i j') + C_(i' j) - C_(i j) - C_(i' j')) / epsilon$. Applying the same bound to the adjoint half-step $v^((t+1)) = c slash.big (K^T u^((t+1)))$ (noting $Delta(K^T) = Delta(K)$) yields
  $ d_H (u^((t+1)), u^*) <= kappa(K) d_H (v^((t)), v^*) <= kappa(K)^2 d_H (u^((t)), u^*), $
  and induction over $t$ establishes Claim 2.

  *(3) Collapse Immunity:* With $r = n^(-1) bold(1)_n$ and $c = k^(-1) bold(1)_k$, the column sums of $Q^* = n P^*$ satisfy $(Q^*)^T bold(1)_n = n c = (n / k) bold(1)_k$, so the empirical cluster frequency vector $hat(pi) = n^(-1) (Q^*)^T bold(1)_n = k^(-1) bold(1)_k$ has maximal Shannon entropy $H(hat(pi)) = log k$. Suppose an encoder collapses all embeddings to $z_i equiv z_0$. Then $C_(i j) = -chevron.l z_0, mu_j chevron.r = c_j$ is independent of $i$, so $K_(i j) = exp(-c_j / epsilon)$ is rank-1: $K = bold(1)_n w^T$. Taking $u^* = n^(-1) bold(1)_n$ and $v^* = k^(-1) w^(-1)$ satisfies $P^* = "diag"(u^*) K "diag"(v^*) = (n k)^(-1) bold(1)_n bold(1)_k^T in Pi(r, c)$, which by Claim 1 is the *unique* optimal transport code. Evaluating SwAV's swapped cross-entropy prediction loss $cal(L)_"swap" = -n^(-1) sum_(i,j) Q_(i j)^* log p_(i j)$ against the collapsed constant softmax $p_(i j) = p_j$ gives
  $ cal(L)_"swap" = -1 / k sum_(j=1)^k log p_j = log k + "KL"(k^(-1) bold(1)_k || p) >= log k, $
  strictly bounded away from the zero loss achieved by well-separated, view-consistent representations ($Q_(i j)^* -> bb(1)[j = j^*(i)]$, $p_(i, j^*(i)) -> 1$, $cal(L)_"swap" -> 0$).
]

#subsec("A", "Semi-Discrete Wasserstein Quantization, Elastic Barycenters & Balanced Cuts")

Entropic Optimal Transport extends beyond SwAV (Alg. 190) and DEPICT (Alg. 185) to unify three structural regimes across our 200-algorithm catalog:
- *Semi-Discrete Wasserstein Quantization:* When the column marginal constraint $c$ is left *free* and $epsilon -> 0^+$, the Wasserstein-2 quantization cost $min_(nu_k = sum_j c_j delta_(mu_j)) W_2^2 (hat(mu)_n, nu_k)$ equals the Euclidean $k$-means objective $n^(-1) J(cal(C), bold(mu))$ identically! Constraining $c = k^(-1) bold(1)_k$ therefore elevates K-Means into *balanced optimal-transport partitioning*, directly generalizing cardinality-balanced graph partitioners such as Kernighan–Lin (Alg. 120), METIS (Alg. 121 @karypis1998metis), Chaco (Alg. 122 @hendrickson1995chaco), and Gini-bounded Genie (Alg. 69 @gagolewski2016genie).
- *Differentiable Elastic Alignment Transport (Category VII):* In time-series clustering, Dynamic Time Warping (DTW @sakoe1978dtw) minimizes alignment cost over monotonic warping paths $A in cal(A)_(L_1, L_2)$. Replacing the hard minimum with Cuturi and Blondel's log-sum-exp entropic relaxation @cuturi2017softdtw produces Soft-DTW (Alg. 164):
  $ "DTW"_gamma (x, y) := -gamma log sum_(A in cal(A)) exp(-(chevron.l A, Delta(x, y) chevron.r) / gamma), $
  which plays the exact structural role of Sinkhorn entropic smoothing on path polytopes—transforming the non-smooth combinatorial DBA barycenter (Alg. 163 @petitjean2011dba) into a differentiable Fréchet mean!
- *Message-Passing Transport Duals (Category V):* Affinity Propagation (Alg. 106 @frey2007affinity) executes max-sum belief propagation on a factor graph enforcing exemplar consistency ($u_(i j) <= u_(j j)$), which is a combinatorial dual-ascent relaxation of the uncapacitated facility-location / $k$-medoids transport program.

#wide_figure(caption: [*Table I:* Objective Geometry, Induced Partition Structure, and Master Mathematical Lenses Across All Eight Clustering Families (Algorithms 1–200).])[
  #table(
    columns: (1.25fr, 1.35fr, 1.25fr, 1.35fr, 1.05fr, 1.15fr),
    align: (left, left, left, left, left, left),
    [Family (Algs.)], [Representative Algorithms], [Native Ambient / Induced Geometry], [Objective Functional Class], [Output Structure], [Primary Lenses],
    [*Cat. I:* Centroid & Partitioning (1–30)], [K-Means, K-Means++, PAM/FasterPAM, K-Medians, Spherical, Kernel, X/G-Means, Trimmed], [Euclidean, $L_1$, spherical $bb(S)^(d-1)$, metric medoid, or RKHS Voronoi cells], [Additive Bregman / metric distortion $sum_i D(x_i || mu_(c_i))$], [Hard flat $k$-partition (or trimmed set)], [Lens 1 (Bregman), Lens 2 ($beta -> infinity$), Lens 5 (BIC/VQ)],
    [*Cat. II:* Density & Mode-Seeking (31–55)], [DBSCAN, HDBSCAN, OPTICS, Mean Shift, DENCLUE 1/2, SNN, Quickshift, DPC, DeBaCl], [Non-convex density superlevel sets $L_lambda (hat(f)_h)$ & Morse–Smale mode basins], [Level-set connectivity, KDE gradient flow, or decision-graph peaks], [Arbitrary-shape clusters + noise set $C_(-1)$], [Lens 3 (Level Sets / Morse), Lens 4 (Heat Drift)],
    [*Cat. III:* Hierarchical & Tree-Based (56–80)], [SLINK, CLINK, UPGMA, Ward, DIANA, BIRCH, CURE, CHAMELEON, Genie, Minimax, BICO, Finch], [Ultrametric tree space $(cal(X), u)$, MST minimax paths, or CF coreset balls], [Lance–Williams merge cost, Dasgupta tree cost, or coreset inertia], [Multi-resolution nested dendrogram $theta(t)$], [Lens 3 (Ultrametrics / $H_0$), Lens 1 (Ward / BICO)],
    [*Cat. IV:* Distribution & Topic (81–102)], [GMM, DPGMM/IGMM, LDA, pLSA, Student-$t$, vMF, NMF, MCLUST, HDP, CTM, BTM, SNOB], [Exponential-family statistical manifolds, ellipsoidal Gaussians, topic simplex], [Marginal log-likelihood, ELBO, KL-NMF, BIC/MML, or Bayes posterior], [Soft posteriors $U in cal(M)_"simp"$ or mixed membership], [Lens 1 (Exp-Family), Lens 2 (Helmholtz EM), Lens 5 (MDL/MML)],
    [*Cat. V:* Graph, Spectral & Network (103–132)], [NJW, Ncut, Ratio Cut, Louvain, Leiden, MCL, Infomap, Walktrap, SBM/DC-SBM, BigCLAM, PIC], [Weighted graph $G=(V,E,W)$, Laplacian eigenspace, random-walk diffusion metric], [Dirichlet energy $"Tr"(H^T L H)$, modularity $Q$, Map Equation, SBM MLE], [Disjoint or overlapping network communities], [Lens 4 (Spectral / Cheeger), Lens 5 (Infomap MDL)],
    [*Cat. VI:* Grid, Subspace & Projection (133–150)], [STING, CLIQUE, WaveCluster, PROCLUS, SUBCLU, ORCLUS, ENCLUS, COSA, P3C, SSC], [Cubical cell grids, axis-parallel / oblique affine subspaces $union_j cal(S)_j$], [Cell density, subspace entropy, local PCA inertia, $ell_1$ self-expression], [Projected / overlapping subspace clusters], [Lens 3 (Cubical $H_0$), Lens 4 (SSC Spectral), Lens 5 (ENCLUS)],
    [*Cat. VII:* Fuzzy, Biclust. & Domain (151–175)], [FCM, PCM, GK, GG, Cheng–Church, FABIA, k-Shape, DBA, Soft-DTW, SLIC, Watershed, CSPA/MCLA], [Power-simplex fuzzy hulls, two-mode submatrices, elastic warp/shape manifolds], [Fuzzy/possibilistic inertia, MSR, Soft-DTW, hypergraph NMI consensus], [Soft / possibilistic $U$, biclusters, superpixels], [Lens 2 (Tsallis / Fuzzy), Lens 5 (NMI), Lens 6 (Soft-DTW)],
    [*Cat. VIII:* Neural, Deep, Stream & Evo. (176–200)], [SOM, GNG, DEC/IDEC, DCN, VaDE, IIC, SCAN, SwAV, CC, CluStream, DenStream, StreamKM++, GKA], [Learned nonlinear neural latent manifold $Phi_omega (cal(X))$, fading stream coresets], [Amortized ELBO, InfoNCE / MI, Sinkhorn OT loss, coreset / stochastic fitness], [Latent clusters, evolving micro-clusters, global optima], [Lens 2 (VaDE/DCN), Lens 5 (IIC/CC), Lens 6 (SwAV OT)],
  )
]

#sec("VIII", "Cross-Lens Synthesis: The Parametric Simplex and Master Six-Tuple Architecture")

Sections II–VII established six powerful mathematical lenses, each unifying several families. We now complete the grand unification by proving *two cross-lens meta-results*: first, Proposition 1 demonstrates that Bregman quantization (Lens 1), kernel mode-seeking (Lens 3), spectral diffusion (Lens 4), and entropic optimal transport (Lens 6) are not isolated formalisms—they form the continuous corners of a single *metric–diffusion–transport simplex* parameterized by bandwidth $h$, diffusion time $t$, and temperature $epsilon$. Second, Theorem 7 embeds *every one of the eight families (Algorithms 1–200)* into a canonical six-tuple master variational architecture.

#proposition(title: "The Metric–Diffusion–Transport Parametric Simplex", origin: "Original Unified Proof")[
  Let $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$ be drawn i.i.d. from a smooth strictly positive density $p in C^3 (bb(R)^d)$ with bounded derivatives, and form the Gaussian affinity graph $W_(i j)^((h)) = exp(-norm(x_i - x_j)_2^2 / (2 h^2))$ with random-walk Laplacian $L_"rw"^((h)) = I_n - (D^((h)))^(-1) W^((h))$ and transition kernel $P^((h)) = I_n - L_"rw"^((h))$. Then the following three continuous cross-lens transitions hold:
  1. *(Spectral–Density Drift Bridge, Lens 4 $<=>$ Lens 3):* For any smooth test function $g in C^3 (bb(R)^d)$, as $n -> infinity$ and $h -> 0^+$ with $(n h^(d+2)) / (log n) -> infinity$, the normalized graph Laplacian operator converges pointwise in probability to the Fokker–Planck drift-diffusion generator:
     $ (1 / h^2) (L_"rw"^((h)) g)(x_i) ->_p -1 / 2 Delta g(x_i) - chevron.l nabla log p(x_i), nabla g(x_i) chevron.r, $
     whose first-order drift field $+nabla log p(x)$ is *identically* the limiting continuous Mean Shift / DENCLUE gradient-ascent velocity field (Category II).
  2. *(Multi-Scale Diffusion–Spectral–Metric Bridge, Lens 4 $<=>$ Lens 1/3):* For fixed $h > 0$, the diffusion distance $D_t^2 (i, j) = sum_(ell=2)^n (1 - lambda_ell)^(2 t) (psi_ell (i) - psi_ell (j))^2$ satisfies:
     - as $t -> infinity$, $(D_t^2 (i, j)) / (1 - lambda_2)^(2 t) -> sum_(ell : lambda_ell = lambda_2) (psi_ell (i) - psi_ell (j))^2$, collapsing the ambient metric onto the Fiedler spectral-cut embedding (Category V);
     - at $t = 1$, $D_1^2 (i, j) = norm(P_(i, dot) - P_(j, dot))_(D^(-1))^2$ measures $1$-step shared-neighborhood overlap, linking Walktrap (Alg. 114) directly to Shared Nearest Neighbor (SNN, Alg. 36) similarity.
  3. *(Entropic Transport–Voronoi Corner, Lens 6 $<=>$ Lens 1/2):* For any cost matrix $C_(i j) = D_phi (x_i || mu_j)$ and uniform data weights $r = n^(-1) bold(1)_n$:
     - relaxing the column constraint $c$ to be optimized freely over $Delta_k$ reduces $W_epsilon (r, c; C)$ identically to the Bregman Helmholtz free energy $n^(-1) cal(F)_(1/epsilon, 1)$ (Lens 2), whose $epsilon -> 0^+$ limit recovers hard Bregman $k$-means (Lens 1);
     - enforcing $c = k^(-1) bold(1)_k$ while sending $epsilon -> 0^+$ yields the exact linear-assignment *balanced $k$-means* partition (Lens 6).
]

#proof[
  *(1) Fokker–Planck Drift-Diffusion Limit:* Conditioned on $x_i = x$, by the law of large numbers (with concentration rate $O((n h^d)^(-1/2))$ uniformly over $h -> 0^+$ when $(n h^(d+2)) / (log n) -> infinity$), the random-walk operator $(P^((h)) g)(x) = (sum_j W_(i j)^((h)) g(x_j)) / (sum_j W_(i j)^((h)))$ converges in probability to the continuous ratio of expectations:
  $ cal(P)_h g(x) := (integral_(bb(R)^d) exp(-norm(y - x)_2^2 / (2 h^2)) g(y) p(y) dif y) / (integral_(bb(R)^d) exp(-norm(y - x)_2^2 / (2 h^2)) p(y) dif y). $
  Change variables to $y = x + h z$ ($dif y = h^d dif z$) and expand $p(x + h z)$ and $g(x + h z)$ in Taylor series to order $h^2$:
  $ p(x + h z) &= p(x) + h chevron.l nabla p(x), z chevron.r + h^2 / 2 z^T nabla^2 p(x) z + O(h^3), \
    g(x + h z) &= g(x) + h chevron.l nabla g(x), z chevron.r + h^2 / 2 z^T nabla^2 g(x) z + O(h^3). $
  Integrating against the isotropic standard Gaussian measure $(2 pi)^(-d/2) e^(-norm(z)_2^2 / 2) dif z$—for which odd moments vanish ($bb(E)[Z] = bold(0)$) and $bb(E)[Z Z^T] = I_d$—gives the numerator and denominator:
  $ "Num" &= g(x) p(x) + h^2 / 2 [p(x) Delta g(x) + 2 chevron.l nabla p(x), nabla g(x) chevron.r \
    &quad + g(x) Delta p(x)] + O(h^3), \
    "Den" &= p(x) + h^2 / 2 Delta p(x) + O(h^3). $
  Dividing $"Num" / "Den"$ and applying $(1 + alpha h^2)^(-1) = 1 - alpha h^2 + O(h^4)$ cancels the $g(x) Delta p(x)$ term identically:
  $ cal(P)_h g(x) = g(x) + h^2 (1 / 2 Delta g(x) + chevron.l (nabla p(x)) / p(x), nabla g(x) chevron.r) + O(h^3). $
  Rearranging $h^(-2) (g(x) - cal(P)_h g(x))$ and letting $h -> 0^+$ yields $-1 / 2 Delta g(x) - chevron.l nabla log p(x), nabla g(x) chevron.r$. Notice that setting $g(y) = y$ (coordinate functions, for which $Delta g = bold(0)$ and $nabla g = I_d$) yields $h^(-2) (cal(P)_h ("id")(x) - x) -> nabla log p(x) = (nabla p(x)) / (p(x))$, which is precisely the population Mean Shift / DENCLUE gradient-ascent step @cheng1995mean, @comaniciu2002mean!

  *(2) Multi-Scale Diffusion Limits:* For a connected non-bipartite graph, $0 = lambda_1 < lambda_2 <= dots <= lambda_n < 2$, so $|1 - lambda_ell| < 1$ for all $ell >= 2$. Factoring out $(1 - lambda_2)^(2 t)$ from $D_t^2 (i, j)$ leaves
  $ (D_t^2 (i, j)) / (1 - lambda_2)^(2 t) &= sum_(ell : lambda_ell = lambda_2) (psi_ell (i) - psi_ell (j))^2 \
    &quad + sum_(ell : lambda_ell > lambda_2) ((1 - lambda_ell) / (1 - lambda_2))^(2 t) (psi_ell (i) - psi_ell (j))^2. $
  Because $|(1 - lambda_ell) / (1 - lambda_2)| < 1$ (assuming the positive spectral gap dominates any negative eigenvalue, which always holds for the lazy walk $tilde(P) = (I + P) / 2$), all higher-frequency terms decay exponentially to $0$ as $t -> infinity$, isolating the Fiedler eigenspace. Conversely, at $t = 1$, $D_1^2 (i, j) = sum_r (P_(i r) - P_(j r))^2 / pi_r$ directly measures the degree-weighted Euclidean discrepancy between the $1$-hop transition rows of $i$ and $j$, vanishing when $i$ and $j$ share identical local neighborhoods (the exact mechanism of SNN, Alg. 36).

  *(3) Entropic OT–Voronoi Corners:* Write $P_(i j) = r_i u_(i j) = n^(-1) u_(i j)$ where $u_i in Delta_k$. Then
  $ chevron.l P, C chevron.r - epsilon H(P) &= 1 / n sum_(i,j) u_(i j) C_(i j) \
    &quad + epsilon / n sum_(i,j) u_(i j) log u_(i j) - epsilon (log n + 1). $
  Up to the constant $-epsilon (log n + 1 + log k)$, minimizing over $U in cal(M)_"simp"$ without a fixed column marginal equals minimizing $n^(-1) cal(F)_(1/epsilon, 1) (U, Theta)$ with uniform prior $pi = k^(-1) bold(1)_k$. By Theorem 2(3), sending $epsilon -> 0^+$ drives $U$ onto hard Voronoi vertices, recovering Bregman $k$-means. When the column constraint $P^T bold(1)_n = k^(-1) bold(1)_k$ is retained, $P^*(epsilon)$ lies on the compact Birkhoff-type polytope $Pi(r, c)$; as $epsilon -> 0^+$, the entropic barrier $-epsilon H(P) -> 0$ uniformly (bounded by $epsilon (log(n k) + 1)$), so every cluster point of $P^*(epsilon)$ solves the linear program $min_(P in Pi(r, c)) chevron.l P, C chevron.r$, whose extreme points yield exact cardinality-balanced partitions when $k$ divides $n$.
]

#wide_figure(caption: [*Table II:* Regularizers, Entropy/Complexity Penalties, Constraint Polytopes, and Optimization Dynamics Across All Eight Clustering Families.])[
  #table(
    columns: (1.18fr, 1.15fr, 1.32fr, 1.25fr, 1.22fr, 1.28fr),
    align: (left, left, left, left, left, left),
    [Family (Algs.)], [Assignment Polytope $cal(M)$], [Entropy / Fuzziness / Transport Regularizer], [Model-Order / Complexity Penalty], [Anti-Collapse / Anti-Chaining], [Optimization Dynamics],
    [*Cat. I:* Centroid & Partitioning (1–30)], [Hard vertices $cal(M)_"hard"$ (or trimmed subset $S$)], [None ($beta -> infinity$); annealed power mean in Alg. 27], [BIC in X-Means (Alg. 10), AD-test in G-Means (Alg. 11), feature $ell_p$ weights], [Empty-cluster reseeding; K-Means++ $D^2$ seeding; $alpha$-trimming], [Alternating block coordinate descent, Hartigan swaps, PAM local search],
    [*Cat. II:* Density & Mode-Seeking (31–55)], [Hard + noise $C_(-1)$ or mode-basin partition], [Kernel bandwidth $h$ / radius $epsilon$ spatial smoothing], [Implicit via $h, epsilon, "minPts"$; Excess-of-Mass persistence in HDBSCAN], [Core-point density floor; mutual reachability (Alg. 32); SNN rank], [Graph BFS/Union-Find, Prim MST, fixed-point KDE gradient ascent],
    [*Cat. III:* Hierarchical & Tree-Based (56–80)], [Ultrametric cone $cal(U)(cal(X))$ / dendrogram], [Space-dilating $beta$ in Flexible Beta (Alg. 74); Gini cap in Genie (Alg. 69)], [Tree cut height, Dasgupta cost, CF-tree branching / threshold $T$], [Complete/Ward variance penalty; Gini balance (Alg. 69); CURE shrinkage], [Greedy agglomerative merge ($O(n^2)$ NN-chain), divisive split, online tree],
    [*Cat. IV:* Distribution & Topic (81–102)], [Simplex $cal(M)_"simp"$ or Dirichlet topic simplex], [Gibbs KL entropy $"KL"(u_i || pi)$; Dirichlet / NIW conjugate priors], [BIC (MCLUST/Mixmod), MML Fisher volume (SNOB), DP concentration $alpha$], [Covariance shrinkage / parsimonious eigen-models; Bayesian priors], [EM majorization-minimization, coordinate CAVI, collapsed Gibbs MCMC],
    [*Cat. V:* Graph, Spectral & Network (103–132)], [Orthonormal $Y^T Y = I_k$, discrete $cal(M)_"hard"$, or $bb(R)_(>= 0)^(n times k)$], [Degree volume normalization $D$; MCL inflation $r$; Spinglass $gamma$], [Modularity null $(d d^T) / (2m)$, Infomap Map Equation MDL, SBM entropy], [Leiden refinement guarantees connected communities; degree correction], [Lanczos / power eigensolver, multilevel greedy moves, belief propagation],
    [*Cat. VI:* Grid, Subspace & Projection (133–150)], [Axis-parallel cells, projected $cal(M)_"hard"$, or sparse $C$], [Cell resolution $xi$; subspace $ell_1$ sparsity norm $norm(C)_1$ (Alg. 150)], [Entropy threshold (ENCLUS), Poisson p-value (P3C), non-redundancy (INSCY)], [Apriori downward-closure pruning; diagonal $c_(i i) = 0$ in SSC], [Bottom-up Apriori search, local PCA / medoid projection, ADMM / LASSO + spectral],
    [*Cat. VII:* Fuzzy, Biclust. & Domain (151–175)], [Simplex $cal(M)_"simp"$, hypercube $cal(M)_"poss"$, or warp paths], [Fuzzifier $m > 1$; PCM typicality $eta_j$; Soft-DTW $gamma$; spatial $W_"sp"$], [Det-constrained covariance (GK), sparse factor prior (FABIA), multi-view $lambda$], [PCM breaks sum-to-one noise trap; co-regularization prevents view drift], [Alternating Picard iteration, dynamic programming, SVD / hypergraph cuts],
    [*Cat. VIII:* Neural, Deep, Stream & Evo. (176–200)], [Simplex $cal(M)_"simp"$, Sinkhorn $Pi(r, c)$, or coreset], [Sinkhorn $-epsilon H(P)$ (SwAV); VAE KL (VaDE); InfoNCE / MI (IIC, CC)], [Autoencoder reconstruction $Omega(omega)$, ART vigilance $rho$, decay $2^(-lambda t)$], [Sinkhorn equipartition (Alg. 190); IDEC/DCN reconstruction; MI entropy], [Stochastic gradient descent + alternating OT/EM, single-pass CF merge, GA/PSO],
  )
]

We are now ready to state and prove the culminating structural theorem of Comparison Chapter I: the *Canonical Six-Tuple Master Representation Theorem*.

#theorem(title: "Canonical Six-Tuple Representation of the 200-Algorithm Universe", origin: "Original Unified Proof")[
  Define a *canonical clustering architecture* as a six-tuple
  $ cal(A) = (cal(Z), Phi_omega, cal(D), cal(R), cal(M), cal(O)), $
  consisting of (1) a latent representation space $cal(Z)$, (2) an embedding map $Phi_omega: Omega -> cal(Z)$ (identity for non-deep methods, RKHS map for kernel methods, or neural encoder for deep methods), (3) a discrepancy, likelihood, or affinity operator $cal(D)$, (4) a complexity, entropy, or topological regularizer $cal(R)$, (5) an assignment manifold $cal(M) in {cal(M)_"hard", cal(M)_"simp", Pi(r, c), cal(M)_"poss", cal(U)(cal(X))}$, and (6) an optimization or functorial dynamics operator $cal(O)$. Then every clustering family $c in {"I", dots, "VIII"}$ in Parts I–VIII embeds rigorously into the *Master Composite Functional*
  $ min_(omega, Theta, U in cal(M)) quad cal(J)_"master" (omega, Theta, U) := cal(D)(Phi_omega (cal(X)), U, Theta) + cal(R)(U, Theta, omega), $
  or into its exact functorial filtration $epsilon |-> arg min_(U in cal(M)_epsilon) cal(J)_"master"$ over a scale parameter $epsilon in [0, infinity)$, satisfying three universal structural dichotomies:
  1. *(Global Convexity vs. Combinatorial Non-Convexity):* Whenever $cal(M) = cal(M)_"hard"$ with $k >= 2$ or whenever prototypes $Theta$ and assignments $U$ are jointly optimized in Euclidean space ($d >= 2$), $cal(J)_"master"$ is non-convex and NP-hard globally; global convexity is achieved *only* when fixing prototypes/features to optimize $U$ over a convex polytope ($cal(M)_"simp"$ or $Pi(r, c)$ with $epsilon > 0$), when relaxing $cal(M)_"hard"$ to the Stiefel manifold $Y^T Y = I_k$ (Lens 4), or when restricting to the unique sub-dominant ultrametric functor $u^*$ (Lens 3).
  2. *(Implicit Metric–Prior Equivalence):* Every hard or soft partitioner with $cal(D) = sum_(i,j) u_(i j) d(Phi(x_i), mu_j)$ induces an implicit generative density $p_j (x) prop exp(-beta d(Phi(x), mu_j))$ (via Theorem 1–2), proving that *no clustering algorithm is assumption-free*: geometric dissimilarity choice is mathematically isomorphic to a distributional likelihood prior.
  3. *(Invariance–Expressivity Conservation):* Any operator $cal(A)$ attaining non-convex manifold recovery (Categories II, III, V, VIII) must sacrifice either global parametric summary compactness (requiring $Omega(n)$ graph/density support instead of $O(k d)$ prototypes) or convex identifiability (introducing multi-modal local optima in $omega$).
]

#proof[
  We construct the explicit six-tuple embedding for all eight families and verify the three universal dichotomies directly from Theorems 1–6:

  - *Family I (Algs. 1–30):* $Phi_omega = "id"$ (or RKHS $phi_K$ in Alg. 16), $cal(D) = sum_(i,j) u_(i j) D_phi (x_i || mu_j)$ (or metric $d(x_i, m_j)$ for medoids), $cal(R) = 0$ (or BIC penalty $(p_k / 2) log n$ in Alg. 10, feature $ell_p$ penalty in Algs. 28–29), $cal(M) = cal(M)_"hard"$ (with cardinality constraint $sum_i u_(i, 0) <= alpha n$ for trimming in Alg. 26), and $cal(O)$ = alternating Voronoi/centroid projection or PAM swap descent.
  - *Family II (Algs. 31–55):* Realized via the level-set filtration $epsilon |-> arg min$ on $cal(U)(cal(X))$ with $cal(D) = d_("mreach", k)$ (Theorem 3(3) for DBSCAN/OPTICS/HDBSCAN/DeBaCl) or continuous gradient-ascent dynamics $cal(O)$ maximizing the KDE potential $-cal(D)(x) = hat(f)_h (x)$ (Proposition 1(1) for Mean Shift / DENCLUE / Quickshift / DPC).
  - *Family III (Algs. 56–80):* $cal(M) = cal(U)(cal(X))$ equipped with sub-dominant ultrametric projection $cal(D)(u, d) = norm(d - u)_(L^1)$ subject to $u <= d$ (Theorem 3(2) for Single Linkage) or Lance–Williams greedy step-minimization of cluster merger inertia / Dasgupta tree cost $cal(L)_"Dasg" (T; W)$ (Algs. 57–80), with BIRCH/BICO (Algs. 64, 78) inserting a coreset compression map $Phi_omega$ preserving Bregman sufficient statistics (Remark after Theorem 1).
  - *Family IV (Algs. 81–102):* $cal(M) = cal(M)_"simp"$, $cal(D) = sum_(i,j) u_(i j) D_phi (x_i || mu_j)$ (via Theorem 1(2) exponential-family duality, or heavy-tailed Student-$t$ / skew-normal log-densities), $cal(R) = sum_i "KL"(u_i || pi) + Omega(Theta)$ (Bayesian conjugate prior, Dirichlet process stick-breaking in DPGMM/HDP Algs. 82, 91, or MML Fisher volume in SNOB Alg. 102 via Theorem 5(1)), and $cal(O)$ = EM / CAVI / collapsed Gibbs.
  - *Family V (Algs. 103–132):* Graph affinity operator $cal(D)(U) = "Tr"(U^T L U)$ (Ncut/Ratio Cut/NJW via Theorem 4(1)), $-"Tr"(U^T B U)$ (Modularity/Louvain/Leiden/Spinglass), negative SBM Poisson/Bernoulli log-likelihood (Algs. 130–131), or Map Equation $L(cal(M))$ (Infomap Alg. 113 via Theorem 5(3)); $cal(M)$ is discrete $cal(M)_"hard"$, Stiefel-relaxed $Y^T Y = I_k$, or non-negative overlapping $bb(R)_(>= 0)^(n times k)$ (BigCLAM Alg. 125).
  - *Family VI (Algs. 133–150):* $Phi_omega$ projects onto axis-parallel cubical cells (STING/CLIQUE/WaveCluster Algs. 133–135), cluster-specific orthonormal subspaces $P_j$ (PROCLUS/ORCLUS Algs. 136, 138 with $cal(D) = sum_(i,j) u_(i j) norm(P_j (x_i - mu_j))_2^2$), or self-expressive sparse affinity $cal(R)(C) = norm(C)_1$ subject to $X = X C, "diag"(C) = bold(0)$ feeding into Lens 4 spectral cuts (SSC Alg. 150).
  - *Family VII (Algs. 151–175):* $cal(M) = cal(M)_"simp"$ or $cal(M)_"poss"$ with power-coupled Tsallis/fuzzy functional (Theorem 2(4) for FCM/PCM/GK/GG/KFCM/sFCM Algs. 151–156), two-mode row-column discrepancy $cal(D)(R, C)$ for biclustering (Algs. 157–161), elastic warp distance $"DTW"_gamma$ for time-series (Algs. 162–164), hypergraph NMI consensus (Algs. 170–173 via Lens 5), or multi-view coupled Dirichlet energy (Algs. 174–175).
  - *Family VIII (Algs. 176–200):* Parametric neural encoder $Phi_omega: bb(R)^d -> cal(Z)$ jointly trained with Helmholtz amortized ELBO (VaDE Alg. 184, DCN Alg. 183, DEC/IDEC Algs. 181–182 via Theorem 2), mutual-information / InfoNCE maximization (IIC Alg. 188, CC Alg. 192 via Lens 5), or Sinkhorn-Knopp optimal transport $cal(M) = Pi(r, c)$ (SwAV Alg. 190, DEPICT Alg. 185 via Theorem 6); streaming algorithms (Algs. 193–197) couple online exponentially decayed coresets with offline Lenses 1–3, while evolutionary metaheuristics (Algs. 198–200) replace gradient/block $cal(O)$ with global Markov-chain stochastic search over $cal(M)_"hard"$.

  *Verification of the Three Dichotomies:*
  (1) NP-hardness of joint Euclidean $k$-means on $cal(M)_"hard"$ holds even at $k = 2$ or $d = 2$ (as proved in Part I, Alg. 1); conversely, fixing $Theta$ and $omega$ makes $U |-> cal(J)_"master"$ strictly convex on $cal(M)_"simp"$ and $Pi(r, c)$ (Theorems 2(1) and 6(1)), while Stiefel spectral relaxation reduces to a symmetric eigenvalue problem (Theorem 4(1)) and sub-dominant ultrametric projection reduces to polynomial-time Kruskal MST (Theorem 3(2)).
  (2) Directly by Theorem 1(2) and Theorem 2(2), every nonnegative divergence $D(Phi(x), mu_j)$ defines a Gibbs Boltzmann factor $p_j (x) = Z_j^(-1) exp(-beta D(Phi(x), mu_j))$ whose MAP decision boundary coincides identically with the hard partitioner's Voronoi cells.
  (3) To represent an arbitrary smooth codimension-$0$ non-convex submanifold $S subset bb(R)^d$ without convex-hull leakage, a fixed-dimension linear/quadratic prototype set of size $O(k d)$ has VC-dimension $O(k d)$ and cannot separate concentric spheres or intertwined spirals with $k = 2$; hence non-convex separation strictly requires either an instance-dependent nonparametric kernel/graph support of size $Omega(n)$ (Categories II, III, V, or RKHS Kernel K-Means) or a nonlinear parametric map $Phi_omega$ of high expressive capacity (Category VIII), whose composition with clustering loss introduces non-convexity in $omega$ and requires explicit anti-collapse regularization (Theorem 6(3)).
]

#wide_figure(caption: [*Table III:* Implicit Statistical Models, Consistency Guarantees, Cluster Geometry Assumptions, Outlier Robustness, and High-Dimensional Concentration Across All Eight Clustering Families.])[
  #table(
    columns: (1.16fr, 1.32fr, 1.24fr, 1.24fr, 1.24fr, 1.20fr),
    align: (left, left, left, left, left, left),
    [Family (Algs.)], [Implicit Statistical Model / Target], [Consistency & Asymptotic Guarantees], [Cluster Geometry & Scale Assumptions], [Outlier / Noise Handling], [High-Dimensional Behavior ($d >> 1$)],
    [*Cat. I:* Centroid & Partitioning (1–30)], [Equal-variance exponential-family mixture (or metric medoids)], [Pollard strong consistency @pollard1981strong; $O(log k)$ K-Means++ ratio @arthur2007kmeanspp], [Convex Voronoi polyhedra; spherical / isotropic unless weighted], [Breakdown $1/n$ for K-Means; $approx 1/(2k)$ for PAM; $alpha$ for Trimmed], [Euclidean collapses under distance concentration @beyer1999nearest; Spherical / Weighted survive],
    [*Cat. II:* Density & Mode-Seeking (31–55)], [Non-parametric density $f(x)$ superlevel sets & modes], [Hartigan / Chaudhuri cluster-tree consistency @chaudhuri2014consistent; Mean Shift mode consistency], [Arbitrary non-convex connected manifolds; multi-scale in HDBSCAN/OPTICS], [Native explicit noise set $C_(-1)$; breakdown up to $1 - "minPts"/n$], [KDE suffers exponential curse $O(n^(-4/(d+4)))$; SNN / SNN-DPC partially mitigate],
    [*Cat. III:* Hierarchical & Tree-Based (56–80)], [Latent ultrametric / phylogenetic or nested tree process], [SLINK consistent for single level-set branch; UPGMA $1/3$-approx Dasgupta], [SLINK: elongated chains; CLINK/Ward: compact convex; CURE: arbitrary], [SLINK vulnerable to noise bridges; Genie/HCA-DBSCAN/CURE resist noise], [Pairwise metric distances concentrate; BICO/BIRCH/Finch rely on metric quality],
    [*Cat. IV:* Distribution & Topic (81–102)], [Parametric / nonparametric generative mixture $sum_j pi_j p(x | theta_j)$], [MLE consistency (regular mixtures); DP posterior consistency; BIC selection], [Ellipsoidal (GMM/MCLUST), directional (vMF), simplex/counts (LDA/Poisson)], [Gaussian MLE breakdown $1/n$; Student-$t$ TMM has bounded influence @peel2000robust], [Full $O(k d^2)$ covariance overfits when $d approx n$; MCLUST parsimonious / topic models scale],
    [*Cat. V:* Graph, Spectral & Network (103–132)], [Stochastic Block Model (SBM / DC-SBM) or Markov diffusion graph], [Exact SBM recovery thresholds @abbe2018sbmreview; $L_"sym"$ consistency @vonluxburg2007tutorial], [Arbitrary topological connectivity; balanced conductance / block density], [Robust to isolated outliers via $L_"sym"$ (unlike unnormalized $L$); OSLOM filters noise], [Depends entirely on affinity graph construction; sparse $k$-NN graphs resist ambient noise],
    [*Cat. VI:* Grid, Subspace & Projection (133–150)], [Union of low-dim subspaces $union_j cal(S)_j$ + uniform background noise], [SSC subspace-detection recovery @soltanolkotabi2012geometric; Apriori completeness], [Axis-parallel hyper-rectangles or oblique low-dim affine subspaces], [Explicit background / sparse-error outlier term ($ell_1$ noise in SSC, P3C)], [Specifically engineered to defeat distance concentration via local $d_j << d$],
    [*Cat. VII:* Fuzzy, Biclust. & Domain (151–175)], [Soft graded mixtures, additive/factor matrices, warped temporal templates], [Bezdek / Hathaway local convergence @hathaway1986fcm; Soft-DTW smooth Fréchet mean], [Overlapping soft boundaries, row-column submatrices, phase/warp-invariant curves], [FCM sensitive to outliers (sum-to-one trap); PCM/FABIA/k-Shape resist noise], [High-$d$ fuzzifier $m$ drives $u_(i j) -> 1/k$; ROCKET / multi-view / biclustering handle high-$d$],
    [*Cat. VIII:* Neural, Deep, Stream & Evo. (176–200)], [Deep latent manifold GMM / contrastive views / non-stationary stream], [GKA Markov global convergence @krishna1999gka; coreset $(1 + epsilon)$ distortion @harpeled2004coresets], [Arbitrary nonlinear semantic manifolds; evolving drifting stream geometries], [DenStream outlier micro-clusters; VaDE/SwAV rely on view augmentation & OT], [Excels on structured high-$d$ perceptual modalities (vision/audio/text) via learned $Phi_omega$],
  )
]

#subsec("A", "Interpreting Tables I–III: Cross-Paradigm Diagnostic Principles")

Reading Tables I, II, and III horizontally reveals four overarching laws that govern algorithm selection across all 200 methods:

1. *The Geometry–Regularization Coupling Law:* Notice across Tables I and II that *no* relaxation of the assignment polytope $cal(M)$ changes cluster shape unless the discrepancy $cal(D)$ or representation $Phi_omega$ changes as well. Moving from K-Means (Cat. I) to Fuzzy C-Means (Cat. VII) or isotropic GMM (Cat. IV) softens cluster boundaries via Lens 2 entropy, yet all three still induce *linear or quadratic Voronoi bisectors* in $bb(R)^d$. Breaking convexity strictly requires either a local neighborhood/density operator (Lenses 3–4: Cat. II, Cat. III, Cat. V, Kernel K-Means) or a nonlinear encoder $Phi_omega$ (Cat. VIII).
2. *The Anti-Collapse / Anti-Chaining Duality:* Table II exposes a striking symmetry between the two extreme failure modes of unsupervised learning. *Prototype and deep representation methods* suffer from *collapse* (all points merging into one prototype or trivial latent constant), cured by Sinkhorn marginal constraints (Lens 6), mutual-information marginal entropy (Lens 5), or reconstruction penalties (Lens 2). Conversely, *connectivity and graph methods* suffer from *chaining and resolution collapse* (weak noise bridges fusing distinct clusters), cured by mutual-reachability core thresholds (Lens 3: HDBSCAN), Gini size bounds (Genie), degree-normalized conductance (Lens 4: $L_"sym"$), or Leiden connectivity refinement.
3. *The Breakdown–Sufficient-Statistic Trade-Off:* Table III confirms the sharp divide predicted by Theorem 1: algorithms admitting exact additive sufficient statistics (K-Means, GMM, Ward, BIRCH, CluStream) achieve unmatched streaming speed ($O(n)$ single-pass) at the exact mathematical price of a $1/n$ breakdown point. Lifting robustness requires either trimming/medoids (Cat. I: PAM, Trimmed K-Means), heavy-tailed scale-mixtures (Cat. IV: Student-$t$ TMM), simplex relaxation (Cat. VII: PCM), or explicit low-density noise excision (Cat. II: DBSCAN/HDBSCAN).
4. *High-Dimensional Survival Regimes:* As established in Table III (and proved formally via Beyer's concentration theorem @beyer1999nearest in Comparison Chapter XI), full-dimensional isotropic Euclidean distances degenerate as $d -> infinity$. Only four mechanisms survive asymptotically: (i) directional normalization onto $bb(S)^(d-1)$ (Spherical K-Means, vMF); (ii) rank-based neighborhood graphs (SNN, $k$-NN spectral); (iii) cluster-specific subspace projection or $ell_1$ self-expression (Cat. VI: PROCLUS, ORCLUS, SSC); and (iv) low-dimensional semantic representation learning $Phi_omega: bb(R)^d -> bb(R)^p$ with $p << d$ (Cat. VIII: VaDE, SwAV, IIC, Contrastive Clustering).

#subsec("B", "Bridge to Comparison Chapters II–XIV")

This grand unified framework serves as the theoretical backbone for the remainder of Part IX:
- *Comparison Chapter II* makes every dashed edge in Figure 1 algebraically explicit, delivering 15 self-contained formal theorems proving exact cross-family equivalences and asymptotic limits (including Weighted Kernel K-Means $<=>$ Ncut, GMM EM $->$ Lloyd as $sigma -> 0^+$, FCM $->$ K-Means as $m -> 1^+$, Single Linkage $<=>$ MST $<=>$ core DBSCAN, Mean Shift $<=>$ Gaussian KDE EM, KL-NMF $<=>$ pLSA, Modularity $<=>$ Potts Spinglass $<=>$ DC-SBM, and SwAV $<=>$ Balanced Entropic OT).
- *Comparison Chapters III–X* zoom inside each of the eight categories to provide head-to-head mathematical trade-offs, theorem-backed failure counterexamples, and complete intra-family comparative tables across all 200 algorithms.
- *Comparison Chapter XI* establishes the universal impossibility and robustness boundaries constraining our entire master functional—proving Kleinberg's Impossibility Theorem @kleinberg2002impossibility, Carlsson–Mémoli functoriality @carlsson2010characterization, high-dimensional measure concentration @beyer1999nearest, stability-symmetry limits, and the comparative breakdown-point hierarchy.
- *Comparison Chapters XII–XIV* operationalize this theory through exhaustive 200-row reference matrices, real 800-run multi-geometry benchmarks from the native Rust `ClusterOverview` engine, and a domain-to-algorithm engineering architecture spanning ten scientific and industrial verticals.

#remark(title: "Unification Takeaway for Theorists and Practitioners")[
  Viewing all 200 algorithms through the six-tuple $(cal(Z), Phi_omega, cal(D), cal(R), cal(M), cal(O))$ converts algorithm selection from trial-and-error folklore into structured mathematical diagnosis. Before invoking any clusterer on empirical data, a practitioner need only identify five observables: (1) *intrinsic geometry* (convex Bregman hulls vs. non-convex density/graph manifolds vs. low-dimensional subspaces); (2) *assignment semantics* (hard disjoint, soft probabilistic, overlapping, or multi-resolution ultrametric); (3) *contamination level* (requiring trimming, heavy tails, or explicit density noise excision); (4) *cardinality balance and collapse risk* (especially in learned neural embeddings, where Sinkhorn OT or MI entropy is mandatory); and (5) *computational and streaming budget* (dictating whether additive Bregman sufficient statistics or sparse local graph sketches are required). Each observable maps directly onto one of our six master lenses—and through them, to the provably appropriate family among Algorithms 1–200.
]
