#import "../template/ieee.typ": *

#alg_entry(
  id: 11,
  name: "G-Means (Gaussian-Means)",
  category: "Cat. I: Centroid & Partitioning",
  year: "2003",
  paradigm: "Recursive Normality-Tested Partitioning",
  complexity: $O(n k d T + n k log n)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

A persistent limitation of Lloyd's K-Means (Alg. 1) @lloyd1982least is the requirement that the cluster count $k$ be specified a priori. The *Gaussian-Means (G-Means)* algorithm, introduced by Hamerly and Elkan @hamerly2003gmeans, automates model selection by coupling recursive binary splitting with a non-parametric goodness-of-fit test for multivariate Gaussianity. Formally, let $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$ be a finite dataset. G-Means models $cal(X)$ as a finite union of clusters whose empirical distributions are each consistent with a single—possibly anisotropic—multivariate normal distribution $cal(N)(mu_j, Sigma_j)$, while retaining the linear-time Voronoi geometry of K-Means rather than fitting full-covariance mixture models via EM @mclachlan2000finite.

Testing multivariate normality directly in $bb(R)^d$ requires $O(d^2)$ covariance degrees of freedom and suffers severe loss of statistical power as $d$ increases. G-Means circumvents this barrier via *directional one-dimensional projection*: for any candidate cluster $C_j subset.eq {1, dots, n}$ of size $m = |C_j|$ with centroid $mu_j$, G-Means executes a local 2-means split of $C_j$ yielding child prototypes $(c_1, c_2)$ and defines the split axis $v = c_1 - c_2 in bb(R)^d$. Projecting $C_j$ onto $v$ produces scalar coordinates
$ y_i = chevron.l x_i, v chevron.r / norm(v)_2^2, quad i in C_j, $
which are standardized to zero sample mean and unit sample variance via $z_i = (y_i - overline(y)) / s_y$, where $overline(y) = m^(-1) sum_(i in C_j) y_i$ and $s_y^2 = (m - 1)^(-1) sum_(i in C_j) (y_i - overline(y))^2$. Let $z_((1)) <= dots <= z_((m))$ denote the sorted standardized projections and $F_i = Phi(z_((i)))$ their images under the standard normal cumulative distribution function $Phi$. G-Means evaluates the *Anderson–Darling (AD) statistic* @anderson1952asymptotic
$ A^2 (Z) = -m - 1 / m sum_(i=1)^m (2 i - 1) [ln F_i + ln(1 - F_(m + 1 - i))], $
together with Stephens's finite-sample correction for estimated mean and variance @hamerly2003gmeans:
$ A_*^2 (Z) = A^2 (Z) (1 + 4 / m - 25 / m^2). $
At significance level $alpha in (0, 1)$ with critical threshold $c_alpha$, the null hypothesis $H_0: {x_i}_(i in C_j) ~ cal(N)(mu_j, Sigma_j)$ is rejected—and the split $(c_1, c_2)$ accepted—whenever $A_*^2 (Z) > c_alpha$.

#definition(title: "Anderson–Darling Gaussian Splitting Equilibrium")[
  Given dataset $cal(X) subset.eq bb(R)^d$, significance level $alpha in (0, 1)$ with critical value $c_alpha > 0$, and minimum splittable cardinality $m_"min" >= 2$, a partition-prototype pair $(cal(C), bold(mu))$ with $k = |cal(C)|$ is a *G-Means equilibrium* if: (i) $(cal(C), bold(mu))$ is a discrete centroidal Voronoi tessellation (CVT) of the Euclidean WCSS functional $J(cal(C), bold(mu)) = sum_(j=1)^k sum_(i in C_j) norm(x_i - mu_j)_2^2$; and (ii) for every cluster $C_j$ with $|C_j| >= m_"min"$, the local 2-means split axis $v_j = c_(j,1) - c_(j,2)$ yields standardized projections $Z_j$ satisfying $A_*^2 (Z_j) <= c_alpha$.
]

To align candidate splits with the direction of maximum variance, Hamerly and Elkan initialize local 2-means along the leading eigenvector $u_1$ of the sample covariance $hat(Sigma)_(C_j)$ using the exact half-normal displacement $plus.minus sqrt(2 lambda_1 / pi) u_1$ (derived in Lemma 1 below). Algorithm 1 formalizes the full recursive procedure.

#algbox(
  title: "G-Means Recursive Normality-Tested Clustering",
  inputs: [Dataset $cal(X) subset.eq bb(R)^d$, initial centers $bold(mu)^((0))$ (default $k_0 = 1$), significance level $alpha$ (critical value $c_alpha$), $m_"min"$, $k_"max"$.],
  outputs: [Discovered cluster count $k$, partition $cal(C) = {C_1, dots, C_k}$, and centroids $bold(mu)$.],
)[
  1. Set $bold(mu) <- bold(mu)^((0))$ and run global Lloyd's K-Means (Alg. 1) on $(cal(X), bold(mu))$ to obtain $(cal(C), bold(mu))$.
  2. *Repeat:* set $"split" <- "false"$ and $bold(mu)_"next" <- emptyset$.
  3. #h(1em) *For each* cluster $C_j in cal(C)$ with centroid $mu_j$ *do*:
  4. #h(2em) *If* $|C_j| < m_"min"$ or $|bold(mu)_"next"| >= k_"max"$ *then* append $mu_j$ to $bold(mu)_"next"$ and *continue*.
  5. #h(2em) Compute leading eigenpair $(lambda_1, u_1)$ of $hat(Sigma)_(C_j)$ (via power iteration) and set $c_1^((0)), c_2^((0)) <- mu_j plus.minus sqrt(2 lambda_1 / pi) u_1$.
  6. #h(2em) Run local 2-means on $C_j$ from $(c_1^((0)), c_2^((0)))$ to obtain $(c_1, c_2)$ and split axis $v <- c_1 - c_2$.
  7. #h(2em) Project $y_i <- chevron.l x_i, v chevron.r / norm(v)_2^2$ for $i in C_j$, standardize $z_i <- (y_i - overline(y)) / s_y$, sort $z_((1)) <= dots <= z_((|C_j|))$, and evaluate $A_*^2 (Z)$.
  8. #h(2em) *If* $A_*^2 (Z) > c_alpha$ *then* append both $c_1, c_2$ to $bold(mu)_"next"$ and set $"split" <- "true"$; *else* append $mu_j$ to $bold(mu)_"next"$.
  9. #h(1em) *If* $"split"$ *then* run global K-Means on $(cal(X), bold(mu)_"next")$ to update $(cal(C), bold(mu))$; *else break*.
  10. *Return* $(|cal(C)|, cal(C), bold(mu))$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five foundational results governing G-Means: (1) exact population optimality of the $sqrt(2 lambda_1 / pi)$ principal-axis split and full scale-sign-isometry invariance of $A_*^2$; (2) an exact closed-form derivation of the discrete Anderson–Darling statistic from its weighted Cramér–von Mises integral; (3) strict global WCSS descent and finite-round termination; (4) post-selection and family-wise Type-I error inflation bounds; and (5) asymptotic split consistency on separated mixtures contrasted with a constructive proof of unbounded over-splitting on unimodal non-Gaussian clusters.

#subsec("A", "Principal-Axis Split Optimality and Projection Invariance")

#lemma(title: "Exact Half-Gaussian Split Displacement and AD Invariance", origin: "Original Derivation")[
  Let $X ~ cal(N)(mu, Sigma)$ on $bb(R)^d$ with largest eigenvalue $lambda_1 > 0$ and unit eigenvector $u_1$.
  1. Partitioning $bb(R)^d$ by the median hyperplane $H = {x : chevron.l x - mu, u_1 chevron.r = 0}$ yields exact conditional centroids $mu_(plus.minus) = mu plus.minus sqrt(2 lambda_1 / pi) u_1$, reducing population WCSS by exactly $(2 / pi) lambda_1$.
  2. For any finite cluster $C_j$, any non-zero scalar $gamma != 0$, and any rigid Euclidean isometry $T(x) = Q x + b$ ($Q in O(d)$), the Anderson–Darling statistic $A_*^2 (Z)$ is strictly invariant under rescaling/negating the split vector $v |-> gamma v$ and under transforming $x_i |-> T(x_i)$.
]

#proof[
  (1) Decompose $X = mu + W u_1 + X_perp$, where $W = chevron.l X - mu, u_1 chevron.r ~ cal(N)(0, lambda_1)$ and $X_perp$ is orthogonal to $u_1$. Because $u_1$ is an eigenvector of $Sigma$, $"Cov"(W, X_perp) = (I_d - u_1 u_1^T) Sigma u_1 = lambda_1 (u_1 - u_1) = bold(0)$; joint Gaussianity therefore implies independence of $W$ and $X_perp$, so $bb(E)[X_perp | W > 0] = bold(0)$. Writing $W = sqrt(lambda_1) G$ with $G ~ cal(N)(0, 1)$,
  $ bb(E)[W | W > 0] &= 2 sqrt(lambda_1) integral_0^infinity g / sqrt(2 pi) e^(-g^2 / 2) dif g \
                     &= sqrt(2 lambda_1 / pi) [-e^(-g^2 / 2)]_0^infinity = sqrt(2 lambda_1 / pi). $
  By symmetry, $bb(E)[W | W <= 0] = -sqrt(2 lambda_1 / pi)$, proving $mu_(plus.minus) = mu plus.minus sqrt(2 lambda_1 / pi) u_1$. By the parallel-axis identity, splitting at $H$ reduces expected squared distortion by $1/2 norm(mu_+ - mu)_2^2 + 1/2 norm(mu_- - mu)_2^2 = (2 / pi) lambda_1$, which is maximal among all hyperplane cuts since $lambda_1 = max_(norm(u)_2 = 1) u^T Sigma u$.

  (2) Scaling $v |-> gamma v$ maps $y_i |-> gamma^(-1) y_i$. For $gamma > 0$, both $(y_i - overline(y))$ and $s_y$ scale by $gamma^(-1)$, leaving $z_i$ identically unchanged. For $gamma < 0$, $z_i |-> -z_i$, reversing the order statistics to $tilde(z)_((i)) = -z_((m + 1 - i))$. By symmetry of the standard normal CDF, $tilde(F)_i = Phi(-z_((m + 1 - i))) = 1 - F_(m + 1 - i)$ and $1 - tilde(F)_(m + 1 - i) = F_i$. Consequently,
  $ ln tilde(F)_i + ln(1 - tilde(F)_(m + 1 - i)) = ln(1 - F_(m + 1 - i)) + ln F_i $
  for every index $i = 1, dots, m$, leaving $A^2 (Z)$ and $A_*^2 (Z)$ strictly unchanged. Finally, under $T(x) = Q x + b$, the covariance transforms as $Q hat(Sigma)_(C_j) Q^T$, the split axis rotates to $Q v$, and $chevron.l Q x_i + b, Q v chevron.r = chevron.l x_i, v chevron.r + chevron.l b, Q v chevron.r$ shifts all projections by the constant $chevron.l b, Q v chevron.r$, which cancels identically upon centering $y_i - overline(y)$.
]

#subsec("B", [Integral Representation and Asymptotic Null Law of $A^2$])

#theorem(title: "Exact Evaluation and Asymptotics of the Anderson–Darling Statistic", origin: "Classical — Cited")[
  Let $F_m (z) = m^(-1) sum_(i=1)^m bb(1)_(z_i <= z)$ be the empirical CDF of standardized scores $z_1, dots, z_m$ with $F_i = Phi(z_((i))) in (0, 1)$. Then the weighted Cramér–von Mises functional @anderson1952asymptotic
  $ A^2 = m integral_(-infinity)^infinity (F_m (z) - Phi(z))^2 / (Phi(z) (1 - Phi(z))) dif Phi(z) $
  evaluates identically to
  $ A^2 = -m - 1 / m sum_(i=1)^m (2 i - 1) [ln F_i + ln(1 - F_(m + 1 - i))]. $
  Moreover, under a simple Gaussian null, $A^2 ->_d sum_(r=1)^infinity Z_r^2 / (r (r + 1))$ with $Z_r ~^("i.i.d.") cal(N)(0, 1)$; when $(mu, sigma)$ are estimated from the sample, $A^2$ converges to a tighter composite Gaussian-bridge quadratic form approximated to $O(m^(-3))$ by $A_*^2$ @anderson1952asymptotic @hamerly2003gmeans.
]

#proof[
  Substitute $u = Phi(z) in [0, 1]$ with order statistics $0 = u_0 < u_1 <= dots <= u_m < u_(m+1) = 1$, where $u_i = F_i$. On each interval $(u_i, u_(i+1))$, $F_m (Phi^(-1)(u)) = i / m$ is constant, so
  $ A^2 = m sum_(i=0)^m integral_(u_i)^(u_(i+1)) (i / m - u)^2 / (u (1 - u)) dif u. $
  Verify the exact partial-fraction decomposition
  $ (i / m - u)^2 / (u (1 - u)) = (i / m)^2 / u + (1 - i / m)^2 / (1 - u) - 1, $
  since $(i / m)^2 (1 - u) + (1 - i / m)^2 u - u (1 - u) = (i / m)^2 - 2 (i / m) u + u^2 = (i / m - u)^2$. Integrating the constant $-1$ over $[0, 1]$ and multiplying by $m$ contributes $-m$. For the $1 / u$ term (vanishing at $i = 0$ and noting $ln u_(m+1) = ln 1 = 0$), summation by parts gives
  $ &m sum_(i=1)^m (i / m)^2 (ln u_(i+1) - ln u_i) \
    &= 1 / m sum_(i=1)^m [(i - 1)^2 - i^2] ln u_i = -1 / m sum_(i=1)^m (2 i - 1) ln u_i. $
  Similarly, for the $1 / (1 - u)$ term (vanishing at $i = m$ with $ln(1 - u_0) = 0$),
  $ &m sum_(i=0)^(m-1) (1 - i / m)^2 [ln(1 - u_i) - ln(1 - u_(i+1))] \
    &= 1 / m sum_(j=1)^m [(m - j)^2 - (m - j + 1)^2] ln(1 - u_j) \
    &= -1 / m sum_(j=1)^m (2 (m - j) + 1) ln(1 - u_j). $
  Re-indexing via $i = m + 1 - j$ maps $2(m - j) + 1 |-> 2 i - 1$ and $u_j |-> u_(m + 1 - i)$. Summing all three parts yields the exact discrete formula. As $m -> infinity$, Anderson and Darling @anderson1952asymptotic establish weak convergence of $sqrt(m)(F_m (Phi^(-1)(u)) - u)$ to a Brownian bridge $B(u)$ whose Legendre eigenfunction expansion yields $sum_(r=1)^infinity Z_r^2 / (r(r+1))$; estimating $(overline(y), s_y)$ projects out the first two Hermite score components, shrinking the asymptotic quantiles and motivating Stephens's $(1 + 4/m - 25/m^2)$ correction @hamerly2003gmeans.
]

#remark(title: "Tail Sensitivity vs. Kolmogorov–Smirnov")[
  The weight $[u (1 - u)]^(-1)$ diverges at $u -> 0^+$ and $u -> 1^-$, amplifying tail discrepancies relative to the unweighted Cramér–von Mises or Kolmogorov–Smirnov tests. This boosts power against heavy-tailed or bimodal mixtures, yet also amplifies sensitivity to isolated outliers.
]

#subsec("C", "Global Energy Descent, Termination, and Multiplicity Bias")

#theorem(title: "Monotone WCSS Descent and Finite-Round Termination", origin: "Original Proof")[
  Let $J_r = J(cal(C)^((r)), bold(mu)^((r)))$ denote the global WCSS after outer splitting round $r >= 0$ of Algorithm 1, and suppose every accepted split produces non-empty distinct child clusters $C_(j,1), C_(j,2)$ with $c_1 != c_2$. Then:
  1. $J_(r+1) <= J_r - sum_(j in cal(S)_r) (|C_(j,1)| norm(c_(j,1) - mu_j)_2^2 + |C_(j,2)| norm(c_(j,2) - mu_j)_2^2) < J_r$, where $cal(S)_r != emptyset$ is the set of clusters split in round $r$.
  2. Algorithm 1 terminates in at most $min(k_"max", n) - k_0 + 1$ outer rounds at a G-Means equilibrium (or upon reaching $k_"max"$).
]

#proof[
  (1) Prior to global K-Means refinement in Step 9, replace each split cluster $C_j$ ($j in cal(S)_r$) by its local 2-means partition ${C_(j,1), C_(j,2)}$ with centroids $(c_(j,1), c_(j,2))$, holding all unsplit clusters fixed. Applying the parallel-axis identity to $C_(j,1)$ and $C_(j,2)$ around the parent centroid $mu_j$ yields
  $ &sum_(i in C_j) norm(x_i - mu_j)_2^2 - sum_(b=1)^2 sum_(i in C_(j,b)) norm(x_i - c_(j,b))_2^2 \
    &= |C_(j,1)| norm(c_(j,1) - mu_j)_2^2 \
    &quad + |C_(j,2)| norm(c_(j,2) - mu_j)_2^2 > 0, $
  where strict positivity holds because $|C_(j,1)|, |C_(j,2)| >= 1$ and $c_(j,1) != c_(j,2)$ implies at least one child differs from $mu_j$. Summing over $j in cal(S)_r$ gives the pre-refinement energy drop. Running global Lloyd's K-Means in Step 9 monotonically non-increases WCSS at every inner iteration (Theorem 1 of Alg. 1), establishing Claim 1.

  (2) Whenever $"split" = "true"$, $|cal(S)_r| >= 1$ and the number of centroids satisfies $k_(r+1) = k_r + |cal(S)_r| >= k_r + 1$. Because $k_r <= min(k_"max", n)$ (each non-empty cluster requires at least one point, and Step 4 halts splitting at $k_"max"$), at most $min(k_"max", n) - k_0$ rounds can execute a split, followed by at most one final verification round with $"split" = "false"$.
]

#proposition(title: "Family-Wise Over-Splitting and Data-Snooping Inflation", origin: "Original Derivation")[
  Suppose $cal(X)$ is generated by $K^*$ well-separated true Gaussian clusters tested independently at nominal level $alpha$ in a single round.
  1. Even ignoring data-dependent axis selection, the family-wise probability of at least one spurious split satisfies
     $ "FWER"(K^*, alpha) = 1 - (1 - alpha)^(K^*) >= 1 - exp(-alpha K^*). $
  2. Because the projection vector $hat(v) = c_1 - c_2$ is chosen by 2-means on the *same* sample $C_j$ to maximize between-group variance, $A_*^2 (Z(hat(v)))$ is stochastically larger than $A_*^2 (Z(v_"fix"))$ for an independent direction $v_"fix"$, inflating finite-sample Type-I error as $d / |C_j|$ grows.
]

#proof[
  (1) Under $K^*$ independent valid level-$alpha$ tests, no cluster splits with probability $(1 - alpha)^(K^*)$. Using the elementary inequality $1 - alpha <= exp(-alpha)$ for all $alpha in (0, 1)$ gives $1 - (1 - alpha)^(K^*) >= 1 - exp(-alpha K^*)$, which exceeds $0.63$ once $K^* >= 1 / alpha$ (e.g., $K^* = 100$ at $alpha = 0.01$); controlling FWER at level $alpha_0$ requires Bonferroni adjustment $alpha_"adj" = alpha_0 / k$.

  (2) Let $C_j$ consist of $m$ i.i.d. draws from $cal(N)(bold(0), I_d)$. Any fixed unit direction $v_"fix"$ satisfies $chevron.l X, v_"fix" chevron.r ~ cal(N)(0, 1)$ exactly. In contrast, the 2-means axis $hat(v)$ maximizes the empirical between-group sum of squares $B(v) = (|C_1| |C_2| / m) (overline(y)_1 - overline(y)_2)^2$, which aligns $hat(v)$ with the random empirical direction exhibiting the largest sample kurtosis deficit (bimodal fluctuation) in $bb(R)^d$. In the proportional regime $d / m -> rho > 0$, empirical projections onto variance-maximizing directions exhibit systematic excess separation over the fixed-direction null, pushing the empirical CDF away from $Phi$ near $z = 0$ and raising the rejection rate above $alpha$.
]

#subsec("D", "Mixture Consistency vs. Unimodal Non-Gaussian Over-Splitting")

#theorem(title: "Asymptotic Split Power and Unbounded Uniform Over-Splitting", origin: "Original Counterexample")[
  1. *Consistent Separation of Gaussian Mixtures:* For any symmetric two-component mixture $P_delta = 1/2 cal(N)(-delta e_1, I_d) + 1/2 cal(N)(delta e_1, I_d)$ with $delta > 0$, the standardized marginal CDF $F_delta != Phi$ yields a strictly positive AD divergence $D_"AD" (F_delta || Phi) > 0$, so $A_*^2 (Z) / m -> D_"AD" (F_delta || Phi)$ almost surely along $e_1$ and G-Means splits $P_delta$ with probability tending to $1$ as $m -> infinity$.
  2. *Unbounded Over-Splitting of Unimodal Non-Gaussian Clusters:* Let $X_1, dots, X_m ~^("i.i.d.") "Unif"[-a, a] subset.eq bb(R)$ be a *single* connected uniform cluster. Then its standardized CDF $G$ satisfies the strict lower bound
     $ D_"AD" (G || Phi) >= (1 - Phi(sqrt(3)))^2 > 1.7 times 10^(-3). $
     Consequently, $A_*^2 (Z) = Omega(m)$ almost surely, causing G-Means to reject normality on the parent and recursively on every uniform sub-interval, driving the estimated cluster count $hat(k)_m -> infinity$ as $m -> infinity$.
  3. *Outlier Breakdown:* The finite-sample breakdown point of G-Means is $epsilon^* = 1 / n$: sending a single point $norm(x_1)_2 -> infinity$ both corrupts its centroid and drives $A_*^2 (Z) -> +infinity$, triggering spurious splitting.
]

#proof[
  (1) Standardizing $P_delta$ along $e_1$ gives variance $1 + delta^2$ and excess kurtosis $-2 delta^4 / (1 + delta^2)^2 < 0$ for all $delta > 0$, so $F_delta != Phi$ on an open interval. By the Glivenko–Cantelli theorem and continuous mapping on the weighted $L^2$ integral, $A^2 (Z) / m -> D_"AD" (F_delta || Phi) = integral_(-infinity)^infinity (F_delta (z) - Phi(z))^2 / (Phi(z)(1 - Phi(z))) dif Phi(z) > 0$ almost surely, whence $A_*^2 (Z) > c_alpha$ eventually a.s.

  (2) Standardizing $"Unif"[-a, a]$ yields $U ~ "Unif"[-sqrt(3), sqrt(3)]$ with CDF $G(z) = 0$ for $z < -sqrt(3)$ and $G(z) = 1$ for $z > sqrt(3)$. Restricting the positive AD integrand to the symmetric tails $|z| > sqrt(3)$ and substituting $u = Phi(z)$ with $u_0 = Phi(sqrt(3)) approx 0.9582$ gives
  $ D_"AD" (G || Phi) &>= 2 integral_(u_0)^1 (1 - u)^2 / (u (1 - u)) dif u \
                      &= 2 integral_(u_0)^1 (1 - u) / u dif u >= 2 integral_(u_0)^1 (1 - u) dif u \
                      &= (1 - u_0)^2 = (1 - Phi(sqrt(3)))^2 > 1.74 times 10^(-3). $
  Because $A_*^2 (Z) / m -> D_"AD" (G || Phi)$ a.s., for any fixed $c_alpha$ the parent interval is split with probability approaching $1$. By symmetry, 2-means bisects $[-a, a]$ into $[-a, 0]$ and $[0, a]$, each of which is again uniform with size $approx m / 2$. By induction on recursion depth $ell$, as long as $m / 2^ell >> c_alpha / D_"AD" (G || Phi)$, every child interval rejects $H_0$ with probability tending to $1$, yielding $hat(k)_m = Omega(m D_"AD" (G || Phi) / c_alpha) -> infinity$.

  (3) Fixing $x_2, dots, x_n$ and sending $norm(x_1)_2 -> infinity$ drives its cluster centroid to infinity (Proposition 1(3) of Alg. 1). Moreover, along the split axis separating $x_1$, the standardized score $z_((m)) -> sqrt(m - 1)$ while $z_((1)), dots, z_((m-1)) -> -(m - 1)^(-1/2)$, making $1 - F_m = 1 - Phi(sqrt(m - 1))$ exponentially small ($approx exp(-(m - 1)/2)$) and driving $-(2m - 1) / m ln(1 - F_m) = Omega(m) -> infinity$ (and even for bounded $m$, pushing $x_1$ off-axis produces a degenerate two-point separation whose $A_*^2$ exceeds $c_alpha$ for any $m >= 4$).
]

#sec("III", "Complexity, Geometry and Operational Profile")

Let $h <= k$ denote the number of outer splitting rounds until G-Means reaches equilibrium, and let $k_r$ be the active cluster count in round $r$. Within round $r$, global Lloyd refinement costs $O(n k_r d T_r)$ operations, whereas local 2-means executes over disjoint clusters satisfying $sum_(j=1)^(k_r) |C_j| = n$, totaling only $O(n d T_"loc")$. Computing the leading eigenvector $u_1$ of $hat(Sigma)_(C_j)$ via $T_"pow"$ power iterations avoids materializing the $d times d$ covariance matrix, costing $O(|C_j| d T_"pow")$ per cluster and $O(n d T_"pow")$ overall—dramatically cheaper than $O(n d^2 + d^3)$ full SVD. Sorting the 1D projections across all clusters costs $sum_j O(|C_j| log |C_j|) <= O(n log n)$ per round. Aggregating over $h <= k$ rounds yields total runtime $O(n k d T + n h log n) subset.eq O(n k d T + n k log n)$ and strictly linear memory $O((n + k) d)$.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for G-Means],
    [Time / Space Complexity], [$O(n k d T + n k log n)$ time; $O((n + k) d)$ memory using power-iteration PCA seeding.],
    [Geometric Inductive Bias], [Arbitrary-covariance Gaussian clusters (anisotropic ellipsoids tolerated by 1D projection test).],
    [Hyperparameter Burden], [Replaces integer $k$ with statistical significance level $alpha$ (typically $alpha in [10^(-4), 10^(-2)]$).],
    [Invariances], [Invariant under $Q in O(d)$, translation, and split-axis scaling; sensitive to anisotropic feature scaling.],
    [Robustness / Breakdown], [$epsilon^* = 1/n$; single outlier inflates tail AD weight and triggers spurious splits (Theorem 3).],
    [Recommended Upgrades], [Bonferroni $alpha / k$ for large $k$; Student-$t$ mixtures @peel2000robust or Trimmed K-Means @cuesta1997trimmed under outliers.],
  ),
  caption: [Summary of theoretical guarantees and operational profile of G-Means.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Telemetry and Sensor Regime Discovery:* Industrial IoT and data-center streams operate across unknown numbers of load or fault regimes whose fluctuations around setpoints are approximately Gaussian; G-Means resolves modes automatically in linear time.
2. *Flow Cytometry and Single-Cell Gating:* Log- or arcsinh-transformed fluorescence intensities form correlated elliptical Gaussian clouds; 1D projection preserves Gaussianity under arbitrary covariance correlation (Proposition 2), avoiding over-fragmentation.
3. *Astronomical and Spectroscopic Source Separation:* Color-magnitude diagrams and spectral catalogs comprise unknown numbers of stellar populations broadened by Gaussian photometric noise, naturally governed by significance level $alpha$.
4. *Adaptive Vector Quantizer Codebook Sizing:* In speech and image compression, G-Means grows codebook size $k$ dynamically until residual quantization cells reduce to unimodal Gaussian noise, allocating bits proportionally to multimodality.
5. *Exploratory Factor and Customer Segmentation:* On standardized behavioral embeddings without a fixed tier count, G-Means outputs both a statistically defensible $k$ and an interpretable hierarchical split lineage.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

G-Means fails predictably under four regimes: (a) *non-Gaussian unimodal clusters* (uniform slabs, skew-normal or heavy-tailed power-law blobs), where Theorem 3(2) guarantees asymptotic over-splitting as $n$ grows; (b) *non-convex manifolds* (curved arcs, concentric rings), which G-Means approximates by tiling many small local Gaussian pieces; (c) *high-dimensional small-sample data* ($d / |C_j|$ large), where post-selection optimization of $hat(v)$ inflates $A_*^2$ (Proposition 1); and (d) *outlier-contaminated datasets*, where extreme values exploit the $[u(1-u)]^(-1)$ tail weight to trigger spurious singleton splits.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Why G-Means Fits / Prerequisite], [Caution / Alternative],
    [IoT / Telemetry], [Unknown $k$; Gaussian sensor noise], [Winsorize spikes first],
    [Flow Cytometry], [Handles elongated Gaussian clouds], [Use GMM @fraley2002mclust for overlap],
    [Astronomy / Spectra], [Principled $alpha$-controlled splitting], [Use $t$-mixtures @peel2000robust for tails],
    [Uniform / Spatial Grid], [Poor fit: splits flat regions], [Use DBSCAN / HDBSCAN],
    [High-D Text / NLP], [Over-splits sparse non-Gaussian data], [Use Spherical (Alg. 12)],
  ),
  caption: [Application suitability and domain engineering guide for G-Means.],
)

#sec("V", "Structural Relations and Empirical Behavior")

G-Means occupies a pivotal bridge between purely geometric partitioning (Lloyd's K-Means, Alg. 1; Bisecting K-Means, Alg. 9) and penalized-likelihood model selection (X-Means, Alg. 10 @pelleg2000xmeans; BIC-selected GMMs @fraley2002mclust @schwarz1978bic). The following result formalizes its decisive theoretical advantage over X-Means on correlated clusters.

#proposition(title: "Anisotropic Ellipsoid Immunity: G-Means vs. X-Means BIC", origin: "Original Derivation")[
  Let $X_1, dots, X_m ~^("i.i.d.") cal(N)(bold(0), Sigma)$ with $Sigma = "diag"(lambda_1, lambda_2, dots, lambda_d)$ and aspect ratio $lambda_1 > lambda_2 = dots = lambda_d = 1$.
  1. For every non-zero split vector $v in bb(R)^d$, the linear projection $chevron.l X, v chevron.r$ is strictly univariate Gaussian $cal(N)(0, v^T Sigma v)$, regardless of the condition number $kappa(Sigma) = lambda_1$. Hence G-Means maintains nominal null rejection behavior under arbitrary linear anisotropy.
  2. In contrast, X-Means @pelleg2000xmeans models each cluster with an isotropic covariance $hat(sigma)^2 I_d$ and splits whenever the Schwarz BIC @schwarz1978bic increases. Along $e_1$, a 2-means split reduces pooled isotropic variance from $overline(lambda) = (lambda_1 + d - 1) / d$ to $overline(lambda)_2 = overline(lambda) - (2 / (pi d)) lambda_1$, increasing isotropic log-likelihood by $Theta(m)$ while incurring only an $O(d ln m)$ BIC penalty—guaranteeing asymptotic over-splitting with probability $1$ whenever $lambda_1$ is large enough that the BIC gain is positive.
]

#proof[
  (1) By the closure of multivariate normal distributions under linear transformations, $X ~ cal(N)(bold(0), Sigma)$ implies $Y = chevron.l X, v chevron.r / norm(v)_2^2 ~ cal(N)(0, v^T Sigma v / norm(v)_2^4)$ for every fixed $v != bold(0)$, so standardization recovers $cal(N)(0, 1)$ identically for all $kappa(Sigma) in [1, infinity)$.

  (2) Under X-Means's spherical Gaussian likelihood @pelleg2000xmeans @bishop2006pattern, the maximum log-likelihood of $m$ points with pooled isotropic variance $hat(sigma)^2$ is $ell = -(m d / 2) ln(2 pi e hat(sigma)^2)$. By Lemma 1(1), splitting $cal(N)(bold(0), Sigma)$ along $e_1$ asymptotically reduces total variance $sum_(r=1)^d "Var"(X_r)$ from $lambda_1 + d - 1$ to $(1 - 2 / pi) lambda_1 + d - 1$, yielding a per-sample isotropic log-likelihood gain
  $ (Delta ell) / m -> d / 2 ln ( (lambda_1 + d - 1) / ((1 - 2 / pi) lambda_1 + d - 1) ) - ln 2. $
  The first term is strictly increasing in $lambda_1$ and converges to $(d / 2) ln(pi / (pi - 2)) approx 0.506 d$ as $lambda_1 -> infinity$; hence for every $d >= 2$ (where $0.506 d > ln 2 approx 0.693$) and sufficiently large $lambda_1$, $Delta ell = Omega(m)$ strictly dominates the BIC complexity penalty $(d + 1) / 2 ln m$, forcing X-Means to fracture the single ellipsoid into multiple spherical balls.
]

In the ClusterOverview 4-dataset empirical benchmark suite (using Silhouette @rousseeuw1987silhouettes and Adjusted Rand Index @hubert1985comparing), G-Means's theoretical behavior is reflected directly in its measured metrics: on *Iris*, G-Means attains *Sil = 0.586* and *ARI = 0.567*—achieving a substantially higher Silhouette than fixed-$k$ Spherical K-Means (0.444) and Elkan's K-Means (0.447) because its Anderson–Darling projection test cleanly separates Setosa while treating the overlapping Versicolor–Virginica cloud as a single elongated Gaussian-like ellipsoid (matching Proposition 2). On *Customers*, G-Means reaches *Sil = 0.342*, outperforming Spherical (0.319) and Elkan (0.334). On *Multi-Shape*, recursive splitting decomposes non-Gaussian geometries into locally compact convex pieces, yielding a strong *Sil = 0.760* (ARI = N/A). On *High-D Subspace+Noise*, G-Means achieves *Sil = 0.309* (ARI = N/A), clearly improving over fixed-$k$ Euclidean Elkan (0.239) by adaptively partitioning multimodal structure. Its mean runtime of *0.338 ms* is roughly $5$–$7 times$ higher than Elkan (0.070 ms) and Spherical K-Means (0.048 ms), directly capturing the overhead of recursive 2-means trials, covariance-axis estimation, 1D sorting, and Anderson–Darling evaluations.
