#import "../template/ieee.typ": *

#alg_entry(
  id: 164,
  name: "Soft-DTW Clustering",
  category: "Cat. VII: Fuzzy, Biclustering & Domain",
  year: "2017",
  paradigm: "Differentiable Elastic Barycenter Optimization",
  complexity: $O(T (K + I) n m^2)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Although Dynamic Time Warping (DTW) @sakoe1978dtw and DBA (Alg. 163) @petitjean2011dba provide effective elastic alignment for time-series clustering, the hard minimum over admissible warping paths makes $"DTW"^2 (mu, x)$ non-differentiable at every alignment-switching boundary and traps coordinate-alternating updates in sharp combinatorial local minima. Inspired by entropic regularization of optimal transport @cuturi2013sinkhorn, Cuturi and Blondel @cuturi2017softdtw introduced *Soft-DTW*, replacing the hard minimum in Bellman's recursion with a log-sum-exp soft-minimum parameterized by smoothing temperature $gamma > 0$. This yields a $C^infinity$-smooth elastic loss whose exact value and gradient with respect to both sequences are computable in quadratic time, enabling quasi-Newton (L-BFGS) and smooth majorization-minimization barycenter optimization inside time-series $k$-Means.

Let $X = {x_1, dots, x_n}$ be univariate time series with $x_i in bb(R)^(L_i)$ (extending directly to $bb(R)^(d times L_i)$ via Euclidean frame cost $norm(mu_p - x_(i, q))_2^2$). Given centroid $mu in bb(R)^m$ and sequence $x in bb(R)^L$, define the pairwise cost matrix $Delta(mu, x) in bb(R)^(m times L)$ by $Delta(mu, x)_(p, q) = (mu_p - x_q)^2$, and let $cal(A)_(m, L) subset {0, 1}^(m times L)$ denote the set of admissible DTW warping paths from $(1, 1)$ to $(m, L)$ with steps in ${(1, 0), (0, 1), (1, 1)}$.

#definition(title: "Soft-Minimum, Soft-DTW, and Smooth Barycenter Clustering")[
  For smoothing temperature $gamma > 0$, define the *soft-minimum operator*
  $ min^gamma {a_1, dots, a_r} = -gamma log sum_(j = 1)^r exp(-a_j / gamma), $
  with convention $min^0 {a_1, dots, a_r} = min_j a_j$. The *Soft-DTW discrepancy* @cuturi2017softdtw between $mu in bb(R)^m$ and $x in bb(R)^L$ is
  $ "dtw"_gamma (mu, x) & = min^gamma_(A in cal(A)_(m, L)) chevron.l A, Delta(mu, x) chevron.r \
    & = -gamma log sum_(A in cal(A)_(m, L)) exp(- chevron.l A, Delta(mu, x) chevron.r / gamma). $
  Associated with $"dtw"_gamma$ is the *Gibbs alignment distribution* over $cal(A)_(m, L)$:
  $ P_gamma (A | mu, x) = (exp(- chevron.l A, Delta(mu, x) chevron.r / gamma)) / (sum_(B in cal(A)_(m, L)) exp(- chevron.l B, Delta(mu, x) chevron.r / gamma)), $
  and its *expected alignment matrix* $E = bb(E)_(A tilde.op P_gamma) [A] in [0, 1]^(m times L)$. *Soft-DTW Clustering* seeks a partition $cal(C) = {C_1, dots, C_K}$ and centroids ${mu_1, dots, mu_K} subset bb(R)^m$ minimizing
  $ J_gamma (cal(C), {mu_k}) = sum_(k = 1)^K sum_(i in C_k) "dtw"_gamma (mu_k, x_i). $
]

Evaluating $"dtw"_gamma (mu, x)$ naively over $|cal(A)_(m, L)|$ paths is exponentially expensive, but dynamic programming computes it in $O(m L)$ operations via the forward recursion
$ r_(p, q) = Delta_(p, q) + min^gamma {r_(p - 1, q), r_(p, q - 1), r_(p - 1, q - 1)}, $
initialized at $r_(0, 0) = 0$ and $r_(p, 0) = r_(0, q) = +infinity$ for $p, q >= 1$, outputting $"dtw"_gamma (mu, x) = r_(m, L)$. Reverse-mode differentiation propagates sensitivities backward from $E_(m, L) = 1$ using soft-transition weights:
$ E_(p, q) & = E_(p + 1, q) e_(p + 1, q)^((1)) + E_(p, q + 1) e_(p, q + 1)^((2)) \
  & quad + E_(p + 1, q + 1) e_(p + 1, q + 1)^((3)), $
where $e_(p + 1, q)^((1)) = exp((r_(p + 1, q) - r_(p, q) - Delta_(p + 1, q)) / gamma)$ and analogously for horizontal and diagonal transitions (with boundary $E_(m + 1, q) = E_(p, L + 1) = 0$).

#algbox(
  title: "Soft-DTW k-Means with Differentiable Barycenters",
  inputs: [Series $X = {x_1, dots, x_n}$, clusters $K$, smoothing $gamma > 0$, outer cap $T$, solver steps $I$],
  outputs: [Partition $cal(C) = {C_1, dots, C_K}$ and smooth barycenters ${mu_1, dots, mu_K}$],
)[
  1. *Initialize*: Select $K$ initial centroids $mu_1^((0)), dots, mu_K^((0)) in bb(R)^m$ from $X$ (e.g., via $k$-Means++ or medoids).
  2. *Assignment*: For each $i in {1, dots, n}$, compute $r_(m, L_i)$ for all $k$ via forward DP and assign $z_i^((t)) = arg min_(1 <= k <= K) "dtw"_gamma (mu_k^((t - 1)), x_i)$.
  3. *Smooth Barycenter Optimization*: For each non-empty $C_k^((t))$, minimize $F_(k, gamma) (mu) = sum_(i in C_k^((t))) "dtw"_gamma (mu, x_i)$ for at most $I$ iterations using L-BFGS (or Soft-MM updates):
     (a) Run forward DP for each $i in C_k^((t))$ to obtain $r_(p, q)$ and $F_(k, gamma) (mu)$.
     (b) Run backward DP for each $i in C_k^((t))$ to obtain $E_i = nabla_Delta "dtw"_gamma (mu, x_i)$.
     (c) Assemble exact gradient $nabla_mu F_(k, gamma) (mu) = 2 sum_(i in C_k^((t))) ("diag"(E_i bold(1)_(L_i)) mu - E_i x_i)$ and perform L-BFGS line search (or closed-form Soft-MM step).
  4. *Termination*: Stop when $z^((t)) = z^((t - 1))$ or $t = T$; else increment $t$ and return to Step 2.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five formal guarantees: exactness of forward-backward DP differentiation, variational Gibbs duality with uniform DTW approximation bounds, the exact second-order Hessian revealing the large-$gamma$ convexification threshold, global Soft-MM descent, and an elementary proof of raw Soft-DTW's negativity and shrinkage anomaly.

#subsec("A", "Exact Forward-Backward Dynamic Programming and Gradient")

#theorem(title: "Exactness of Forward DP and Reverse-Mode Gradient", origin: "Classical — Cited")[
  For every $gamma > 0$, $mu in bb(R)^m$, and $x in bb(R)^L$ @cuturi2017softdtw:
  (i) The forward recursion satisfies $r_(m, L) = "dtw"_gamma (mu, x)$;
  (ii) $"dtw"_gamma (mu, x)$ is $C^infinity$-smooth in $Delta$ and $mu$, with $nabla_Delta "dtw"_gamma (mu, x) = bb(E)_(P_gamma) [A] = E$ computed exactly by the backward recursion in $O(m L)$ time;
  (iii) Under squared Euclidean cost $Delta_(p, q) = (mu_p - x_q)^2$,
  $ nabla_mu "dtw"_gamma (mu, x) = 2 ("diag"(E bold(1)_L) mu - E x). $
]

#proof[
  *(i)* Define partial partition functions $Z_(p, q) = exp(-r_(p, q) / gamma)$ for $1 <= p <= m$, $1 <= q <= L$, with $Z_(0, 0) = 1$ and $Z_(p, 0) = Z_(0, q) = 0$ for $p, q >= 1$. Exponentiating $-r_(p, q) / gamma$ converts the forward recursion into
  $ Z_(p, q) & = exp(-Delta_(p, q) / gamma) \
    & quad times (Z_(p - 1, q) + Z_(p, q - 1) + Z_(p - 1, q - 1)). $
  By induction on $p + q$, because every admissible path in $cal(A)_(p, q)$ ending at $(p, q)$ enters uniquely from $(p - 1, q)$, $(p, q - 1)$, or $(p - 1, q - 1)$, we have $Z_(p, q) = sum_(A in cal(A)_(p, q)) exp(-chevron.l A, Delta chevron.r / gamma)$. Taking $-gamma log Z_(m, L)$ proves $r_(m, L) = "dtw"_gamma (mu, x)$.

  *(ii)* Smoothness follows because $Z_(m, L) > 0$ is a finite sum of smooth positive exponentials. Differentiating $-gamma log Z_(m, L)$ directly with respect to $Delta_(p, q)$ gives
  $ (partial "dtw"_gamma) / (partial Delta_(p, q)) = (sum_(A in cal(A)_(m, L)) A_(p, q) exp(-chevron.l A, Delta chevron.r / gamma)) / Z_(m, L) = bb(E)_(P_gamma) [A_(p, q)]. $
  In the computational DAG of the forward DP, $Delta_(p, q)$ influences $r_(m, L)$ exclusively through $r_(p, q)$, which feeds only into $r_(p + 1, q)$, $r_(p, q + 1)$, and $r_(p + 1, q + 1)$. Since $partial min^gamma {a, b, c} / partial a = exp((min^gamma {a, b, c} - a) / gamma)$, we obtain $partial r_(p + 1, q) / partial r_(p, q) = exp((r_(p + 1, q) - Delta_(p + 1, q) - r_(p, q)) / gamma) = e_(p + 1, q)^((1))$, and symmetrically for the other two outgoing edges. The multivariate chain rule yields $E_(p, q) = partial r_(m, L) / partial r_(p, q) = partial r_(m, L) / partial Delta_(p, q)$ via the exact backward recursion in $O(m L)$ operations.

  *(iii)* Applying the chain rule to $Delta_(p, q) = (mu_p - x_q)^2$ gives
  $ (partial "dtw"_gamma) / (partial mu_p) = sum_(q = 1)^L E_(p, q) dot.op 2 (mu_p - x_q) = 2 (E bold(1)_L)_p mu_p - 2 (E x)_p, $
  which in vector form is $2 ("diag"(E bold(1)_L) mu - E x)$.
]

#subsec("B", "Gibbs Variational Duality and Tight Approximation Bounds")

#theorem(title: "Entropic Variational Duality and Uniform DTW Envelope", origin: "Original Proof")[
  Let $Delta(cal(A)_(m, L))$ be the probability simplex over $cal(A)_(m, L)$ and $H(p) = -sum_(A in cal(A)_(m, L)) p(A) log p(A)$ be the Shannon entropy. For every $gamma > 0$, $mu in bb(R)^m$, and $x in bb(R)^L$:
  (i) $"dtw"_gamma (mu, x) = min_(p in Delta(cal(A)_(m, L))) { sum_A p(A) chevron.l A, Delta(mu, x) chevron.r - gamma H(p) }$, uniquely minimized at $p^ast = P_gamma (dot.op | mu, x)$;
  (ii) The uniform smoothing sandwich holds everywhere:
  $ "DTW"^2 (mu, x) - gamma log |cal(A)_(m, L)| <= "dtw"_gamma (mu, x) <= "DTW"^2 (mu, x); $
  (iii) For $m = L$, $log |cal(A)_(m, m)| <= 2 (m - 1) log(1 + sqrt(2)) < 1.763 (m - 1)$, whence $"dtw"_gamma -> "DTW"^2$ uniformly on $bb(R)^m times bb(R)^L$ as $gamma -> 0^+$.
]

#proof[
  *(i)* For any distribution $p in Delta(cal(A)_(m, L))$, expand the Kullback-Leibler divergence against the Gibbs law $P_gamma (A) = Z^(-1) exp(-chevron.l A, Delta chevron.r / gamma)$:
  $ "KL"(p || P_gamma) & = sum_A p(A) log (p(A)) / (P_gamma (A)) \
    & = -H(p) + 1 / gamma sum_A p(A) chevron.l A, Delta chevron.r + log Z. $
  Multiplying by $gamma > 0$ and rearranging gives
  $ sum_A p(A) chevron.l A, Delta chevron.r - gamma H(p) & = -gamma log Z + gamma "KL"(p || P_gamma) \
    & >= "dtw"_gamma (mu, x), $
  with equality if and only if $"KL"(p || P_gamma) = 0$, i.e., $p = P_gamma$.

  *(ii)* Let $c^ast = "DTW"^2 (mu, x) = min_(A in cal(A)_(m, L)) chevron.l A, Delta(mu, x) chevron.r$. Every term in $Z = sum_(A in cal(A)_(m, L)) exp(-chevron.l A, Delta chevron.r / gamma)$ is positive and bounded above by $exp(-c^ast / gamma)$, while at least one optimal path attains $exp(-c^ast / gamma)$. Thus
  $ exp(-c^ast / gamma) <= Z <= |cal(A)_(m, L)| exp(-c^ast / gamma). $
  Applying the strictly decreasing map $u |-> -gamma log u$ yields $c^ast - gamma log |cal(A)_(m, L)| <= "dtw"_gamma (mu, x) <= c^ast$.

  *(iii)* For $m = L$, $|cal(A)_(m, m)|$ equals the central Delannoy number $D(m - 1, m - 1) = sum_(j = 0)^(m - 1) binom(m - 1, j)^2 2^j$. Because $sum_j a_j^2 <= (sum_j a_j)^2$ for $a_j = binom(m - 1, j) 2^(j / 2) >= 0$, we obtain $D(m - 1, m - 1) <= (1 + sqrt(2))^(2(m - 1))$, so $log |cal(A)_(m, m)| <= 2(m - 1) log(1 + sqrt(2)) < 1.763 (m - 1)$. Since the bound is independent of $(mu, x)$, $"sup"_(mu, x) |"dtw"_gamma (mu, x) - "DTW"^2 (mu, x)| <= gamma log |cal(A)_(m, L)| -> 0$ as $gamma -> 0^+$.
]

#subsec("C", "Exact Hessian and Large-Temperature Convexification")

#theorem(title: "Exact Second-Order Hessian and Convexity Threshold", origin: "Original Derivation")[
  Fix $x in bb(R)^L$ and $gamma > 0$, and let $g_A (mu) = nabla_mu chevron.l A, Delta(mu, x) chevron.r = 2 ("diag"(A bold(1)_L) mu - A x)$ for each $A in cal(A)_(m, L)$. Then the exact Hessian of $"dtw"_gamma (dot.op, x)$ is
  $ nabla_mu^2 "dtw"_gamma (mu, x) = 2 "diag"(E bold(1)_L) - 1 / gamma "Cov"_(A tilde.op P_gamma) (g_A (mu)). $
  Consequently:
  (i) Globally, $nabla_mu^2 "dtw"_gamma (mu, x) prec.eq 2 L I_m$;
  (ii) On any bounded hypercube $norm(mu)_infinity <= R$ with $norm(x)_infinity <= R$, $"dtw"_gamma (dot.op, x)$ is *strictly convex* whenever $gamma > 8 m L^2 R^2$;
  (iii) Conversely, when $gamma$ is sufficiently small near a multi-path boundary where $"Var"_(P_gamma) (v^top g_A (mu)) > 2 gamma v^top "diag"(E bold(1)_L) v$ for some $v != 0$, $nabla_mu^2 "dtw"_gamma (mu, x)$ acquires a strictly negative eigenvalue.
]

#proof[
  Write $phi_A (mu) = chevron.l A, Delta(mu, x) chevron.r$ so $nabla_mu "dtw"_gamma (mu, x) = sum_A P_gamma (A) g_A (mu)$ with $g_A (mu) = nabla phi_A (mu)$ and $nabla^2 phi_A (mu) = 2 "diag"(A bold(1)_L)$. Differentiating $P_gamma (A) = Z^(-1) exp(-phi_A (mu) / gamma)$ gives
  $ nabla_mu P_gamma (A) = -1 / gamma P_gamma (A) (g_A (mu) - bb(E)_(P_gamma) [g_B (mu)]). $
  Applying the product rule to $sum_A P_gamma (A) g_A (mu)$ yields
  $ nabla_mu^2 "dtw"_gamma (mu, x) & = sum_A P_gamma (A) (2 "diag"(A bold(1)_L)) \
    & quad - 1 / gamma sum_A P_gamma (A) (g_A - bar(g)) g_A^top \
    & = 2 "diag"(E bold(1)_L) - 1 / gamma "Cov"_(P_gamma) (g_A (mu)), $
  using $sum_A P_gamma (A) (g_A - bar(g)) bar(g)^top = 0$ where $bar(g) = bb(E)_(P_gamma) [g_B]$.

  *(i)* Every admissible path $A in cal(A)_(m, L)$ visits each row $p$ at least once and at most $L$ times, so $1 <= (A bold(1)_L)_p <= L$ deterministically, hence $I_m prec.eq "diag"(E bold(1)_L) prec.eq L I_m$. Since $"Cov"_(P_gamma) (g_A) succ.eq 0$, we obtain $nabla_mu^2 "dtw"_gamma prec.eq 2 L I_m$.

  *(ii)* On $[-R, R]^m$ with $norm(x)_infinity <= R$, $|mu_p - x_q| <= 2 R$, so $|(g_A (mu))_p| <= 2 L (2 R) = 4 L R$ and $norm(g_A (mu))_2^2 <= 16 m L^2 R^2$ for every $A$. Thus $lambda_max ("Cov"_(P_gamma) (g_A)) <= bb(E)_(P_gamma) [norm(g_A)_2^2] <= 16 m L^2 R^2$. Combined with $2 "diag"(E bold(1)_L) succ.eq 2 I_m$,
  $ nabla_mu^2 "dtw"_gamma (mu, x) succ.eq (2 - (16 m L^2 R^2) / gamma) I_m succ 0 $
  for all $gamma > 8 m L^2 R^2$. *(iii)* Evaluating the quadratic form $v^top nabla_mu^2 "dtw"_gamma v = 2 v^top "diag"(E bold(1)_L) v - gamma^(-1) "Var"_(P_gamma) (v^top g_A)$ directly establishes negativity whenever the variance term dominates.
]

Theorem 3 rigorously explains the empirical advantage of Soft-DTW over DBA: increasing $gamma$ suppresses the negative path-covariance curvature, convexifying the optimization landscape (and enabling homotopy continuation from large $gamma$ down to small $gamma$).

#subsec("D", "Soft-MM Descent and Smooth Barycenter Stationarity")

#proposition(title: "Global Soft-MM Surrogate and Fixed-Point Equivalence", origin: "Original Proof")[
  For cluster $C_k$ and iterate $mu^((r))$ with expected alignments $E_i^((r)) = nabla_Delta "dtw"_gamma (mu^((r)), x_i)$, define the quadratic surrogate
  $ Q_gamma (mu | mu^((r))) & = F_(k, gamma) (mu^((r))) \
    & quad + sum_(i in C_k) chevron.l E_i^((r)), Delta(mu, x_i) - Delta(mu^((r)), x_i) chevron.r. $
  Then: (i) $F_(k, gamma) (mu) <= Q_gamma (mu | mu^((r)))$ globally for all $mu in bb(R)^m$; (ii) minimizing $Q_gamma (dot.op | mu^((r)))$ yields the closed-form *Soft-DBA* step
  $ mu^((r + 1)) = (sum_(i in C_k) "diag"(E_i^((r)) bold(1)_(L_i)))^(-1) (sum_(i in C_k) E_i^((r)) x_i), $
  with strict Lyapunov descent $F_(k, gamma) (mu^((r))) - F_(k, gamma) (mu^((r + 1))) >= |C_k| norm(mu^((r + 1)) - mu^((r)))_2^2$; (iii) $mu^ast$ is a fixed point of Soft-DBA if and only if $nabla_mu F_(k, gamma) (mu^ast) = 0$, and every such stationary point lies in $[min_(i, q) x_(i, q), max_(i, q) x_(i, q)]^m$.
]

#proof[
  *(i)* By Theorem 2(i), $Delta |-> "dtw"_gamma (Delta)$ is a pointwise infimum of affine functions of $Delta$, hence globally *concave* in $Delta$. By the first-order concavity inequality with $E_i^((r)) = nabla_Delta "dtw"_gamma (Delta(mu^((r)), x_i))$,
  $ "dtw"_gamma (Delta(mu, x_i)) & <= "dtw"_gamma (Delta(mu^((r)), x_i)) \
    & quad + chevron.l E_i^((r)), Delta(mu, x_i) - Delta(mu^((r)), x_i) chevron.r. $
  Summing over $i in C_k$ gives $F_(k, gamma) (mu) <= Q_gamma (mu | mu^((r)))$ with equality at $mu = mu^((r))$.

  *(ii)* Up to an additive constant, $Q_gamma (mu | mu^((r)))$ equals $sum_(p = 1)^m sum_(i in C_k) sum_(q = 1)^(L_i) (E_i^((r)))_(p, q) (mu_p - x_(i, q))^2$. Because $(E_i^((r)) bold(1)_(L_i))_p >= 1$ for every $i, p$, the weight $w_p^((r)) = sum_(i in C_k) (E_i^((r)) bold(1)_(L_i))_p >= |C_k| > 0$. Setting the derivative to zero gives the unique minimizer $mu^((r + 1))$, and expanding $Q_gamma$ around $mu^((r + 1))$ yields
  $ F_(k, gamma) (mu^((r))) - F_(k, gamma) (mu^((r + 1))) & >= sum_(p = 1)^m w_p^((r)) (mu_p^((r)) - mu_p^((r + 1)))^2 \
    & >= |C_k| norm(mu^((r + 1)) - mu^((r)))_2^2. $

  *(iii)* By Theorem 1(iii), $nabla_mu F_(k, gamma) (mu^ast) = 2 sum_(i in C_k) ("diag"(E_i^ast bold(1)_(L_i)) mu^ast - E_i^ast x_i) = 0$ holds if and only if $mu^ast = (sum_i "diag"(E_i^ast bold(1)_(L_i)))^(-1) (sum_i E_i^ast x_i)$. Since every entry of $E_i^ast$ is nonnegative, each coordinate $mu_p^ast$ is a convex combination of ${x_(i, q)}$, enforcing range confinement.
]

#subsec("E", "Negativity Anomaly and Debiased Soft-DTW Divergence")

#proposition(title: "Negativity and Non-Minimal Self-Discrepancy of Raw Soft-DTW", origin: "Original Proof")[
  For every $gamma > 0$ and $m >= 2$: (i) any constant series $x = c bold(1)_m$ satisfies $"dtw"_gamma (x, x) = -gamma log D(m - 1, m - 1) < 0$; (ii) raw $"dtw"_gamma$ is NOT minimized at $y = x$: specifically, for $gamma = 1$, $x = (0, 1)^top$, and $y = (1 / 4, 3 / 4)^top$, we have $"dtw"_1 (x, y) < "dtw"_1 (x, x) < 0$.
]

#proof[
  *(i)* For $x = c bold(1)_m$, $Delta(x, x) = 0$, so every path in $cal(A)_(m, m)$ has zero cost and $Z = |cal(A)_(m, m)| = D(m - 1, m - 1) >= 3$, whence $"dtw"_gamma (x, x) = -gamma log D(m - 1, m - 1) < 0$.

  *(ii)* For $m = 2$, $cal(A)_(2, 2)$ consists of three paths: diagonal $A_1 = {(1, 1), (2, 2)}$, $A_2 = {(1, 1), (1, 2), (2, 2)}$, and $A_3 = {(1, 1), (2, 1), (2, 2)}$. At $x = (0, 1)^top$ and $gamma = 1$, the path costs of $(x, x)$ are $0, 1, 1$, so $Z(x, x) = 1 + 2 e^(-1)$. For $y = (1 / 4, 3 / 4)^top$, $Delta(x, y)_(1, 1) = Delta(x, y)_(2, 2) = 1 / 16$ and $Delta(x, y)_(1, 2) = Delta(x, y)_(2, 1) = 9 / 16$, giving path costs $1 / 8$, $11 / 16$, and $11 / 16$, hence $Z(x, y) = e^(-1 / 8) + 2 e^(-11 / 16)$. Applying the elementary inequality $e^u >= 1 + u$ to $e^(-1 / 8) >= 7 / 8$ and $e^(-11 / 16) = e^(-1) e^(5 / 16) >= e^(-1) (21 / 16)$ yields
  $ Z(x, y) - Z(x, x) >= (7 / 8 - 1) + 2 (21 / 16 - 1) e^(-1) = (5 e^(-1) - 1) / 8 > 0, $
  because $e < 3 < 5$. Taking $-log$ yields $"dtw"_1 (x, y) < "dtw"_1 (x, x) < 0$.
]

Proposition 2 explains why raw Soft-DTW barycenters exhibit *amplitude shrinkage* toward smoother or flatter series (which increase alignment path entropy), motivating the debiased *Soft-DTW divergence* $S_gamma (mu, x) = "dtw"_gamma (mu, x) - 1 / 2 "dtw"_gamma (mu, mu) - 1 / 2 "dtw"_gamma (x, x)$ whenever amplitude preservation at moderate-to-large $gamma$ is critical.

#sec("III", "Complexity, Geometry and Operational Profile")

Each assignment step evaluates $K n$ forward Soft-DTW passes in $O(K n m^2)$ time. Each barycenter update runs $I$ iterations of L-BFGS (or Soft-MM, Proposition 1), where every function-and-gradient evaluation requires one forward DP pass ($O(m^2)$) and one backward DP pass ($O(m^2)$) per series in $C_k$, totaling $O(I n m^2)$ across clusters per outer step and $O(T (K + I) n m^2)$ overall. Although the asymptotic scaling matches DBA (Alg. 163), the constant factor is $5$--$10 times$ larger in practice because each inner cell executes three `exp` and one `log` operation (stabilized via max-shifting $min^gamma {a, b, c} = m_0 - gamma log(sum e^(-(a_j - m_0) / gamma))$ with $m_0 = min{a, b, c}$ to prevent underflow) and L-BFGS performs multiple line-search evaluations per step.

Memory consumption is $O(n m + K m + m^2)$ because the forward matrix $R$ ($m times m$) must be retained during the backward sweep for $E$; for $m$ up to several thousand, $R$ and $E$ fit in L2/L3 cache or GPU shared memory, where anti-diagonal wavefront parallelism ($p + q = "const"$) achieves $O(m)$ parallel span across $O(m)$ GPU threads @cuturi2017softdtw.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Distance Primitive], [Entropic log-sum-exp $"dtw"_gamma (mu, x)$ over $cal(A)_(m, L)$],
    [Centroid Operator], [L-BFGS gradient descent or Soft-MM expected alignment averaging],
    [Time Complexity], [$O(T (K + I) n m^2)$ with higher `exp`/`log` constant factor],
    [Space Complexity], [$O(n m + K m + m^2)$ (full $R$ matrix kept for reverse-mode DP)],
    [Differentiability], [$C^infinity$-smooth in $mu$ and $x$ for all $gamma > 0$ (Theorem 1)],
    [Curvature Regime], [Non-convex for small $gamma$; strictly convex for $gamma > 8 m L^2 R^2$],
    [Key Hyperparameter], [Temperature $gamma > 0$ controlling smoothing vs. shrinkage bias],
  ),
  caption: [Operational and geometric profile of Soft-DTW Clustering (Alg. 164).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *End-to-End Deep Time-Series Clustering and Representation Learning*: Because $"dtw"_gamma$ provides exact analytical gradients through reverse-mode DP (Theorem 1), it plugs directly as a differentiable clustering loss into neural autoencoders and temporal convolutional networks, jointly training latent trajectories and elastic cluster centroids via backpropagation @cuturi2017softdtw.
2. *Multivariate Biosignal and EEG/MEG Trial Averaging*: Multi-channel neural evoked potentials exhibit variable inter-trial latency jitter alongside low signal-to-noise ratios. Averaging over the Gibbs distribution $P_gamma$ instead of a single hard Viterbi path reduces alignment variance and yields physiologically smooth multi-channel centroids without jagged hard-switching artifacts.
3. *Robotics Trajectory Learning from Demonstration (LfD)*: Human kinematic demonstrations of manipulation tasks vary in velocity profiles and pause durations. Soft-DTW barycenters generate $C^infinity$-smooth, dynamically feasible reference trajectories whose differentiability supports direct integration with smooth trajectory optimization and MPC controllers.
4. *Financial Regime and Intraday Volatility Curve Clustering*: Intraday volume and volatility trajectories exhibit noisy timing shifts around macro announcements and market opens. Moderate-$gamma$ smoothing mitigates single-tick noise spikes that would cause hard DTW paths to jump discontinuously between clusters.
5. *Audio, Speech Prosody, and Acoustic Spectrogram Prototyping*: Pitch contours and mel-spectrogram sequences possess both time-warping variability and continuous harmonic structure; smooth Soft-DTW barycenters preserve spectral continuity better than hard DBA averages, particularly when combined with homotopy annealing on $gamma$.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Soft-DTW Clustering presents four specific failure regimes: (i) *amplitude shrinkage and negativity bias* at large $gamma$ (Proposition 2), where raw $"dtw"_gamma$ pulls centroids toward flatter trajectories and can assign a series closer to a flattened prototype than to itself (remedied by Soft-DTW divergence $S_gamma$ or smaller $gamma$); (ii) *numerical overflow/underflow* for very small $gamma << min Delta_(p, q)$, requiring log-sum-exp max-shifting; (iii) *high computational latency* on CPUs for large $n$ or long $m$ due to transcendental DP operations inside iterative L-BFGS loops; and (iv) *lack of scale/offset invariance* on unnormalized signals, requiring explicit preprocessing if amplitude invariance is desired (unlike k-Shape, Alg. 162).

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Deep Temporal Clustering], [Very High: exact $C^infinity$ gradients], [Tune $gamma$ to avoid shrinkage],
    [Noisy Multi-Channel Bio], [Very High: Gibbs variance reduction], [Use debiased $S_gamma$ at large $gamma$],
    [Robot LfD Trajectories], [High: smooth elastic reference paths], [GPU wavefront DP for long $m$],
    [Real-Time Edge / IoT], [Low: high `exp`/`log` CPU latency], [Use ROCKET (Alg. 165) or k-Shape (Alg. 162)],
    [Strict Scale-Free Shapes], [Moderate: requires manual $z$-norm], [Use k-Shape (Alg. 162)],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for Soft-DTW.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Soft-DTW Clustering establishes a continuous one-parameter homotopy connecting combinatorial DTW partitioning (Alg. 163) to convex path-averaged quadratic clustering.

#proposition(title: "Asymptotic Temperature Limits of Soft-DTW Clustering", origin: "Original Derivation")[
  For any finite dataset $X = {x_1, dots, x_n}$ and centroids ${mu_k}$:
  (i) As $gamma -> 0^+$, $J_gamma (cal(C), {mu_k}) -> J_"DBA" (cal(C), {mu_k})$ uniformly on $bb(R)^(K m)$, recovering DBA (Alg. 163);
  (ii) As $gamma -> +infinity$, the shifted objective satisfies
  $ lim_(gamma -> +infinity) & ("dtw"_gamma (mu, x) + gamma log |cal(A)_(m, L)|) \
    & = 1 / (|cal(A)_(m, L)|) sum_(A in cal(A)_(m, L)) chevron.l A, Delta(mu, x) chevron.r, $
  which is strictly convex in $mu$ with unique closed-form cluster minimizer $mu_k = (sum_(i in C_k) "diag"(bar(A)_i bold(1)_(L_i)))^(-1) (sum_(i in C_k) bar(A)_i x_i)$, where $bar(A)_i = |cal(A)_(m, L_i)|^(-1) sum_(A in cal(A)_(m, L_i)) A$.
]

#proof[
  *(i)* Immediate by summing Theorem 2(ii)--(iii) across all $n$ data points: $|J_gamma - J_"DBA"| <= gamma sum_(i = 1)^n log |cal(A)_(m, L_i)| -> 0$ uniformly as $gamma -> 0^+$.
  *(ii)* Expanding $exp(-chevron.l A, Delta chevron.r / gamma) = 1 - chevron.l A, Delta chevron.r / gamma + O(gamma^(-2))$ inside $Z = sum_A exp(-chevron.l A, Delta chevron.r / gamma)$ gives
  $ Z / (|cal(A)_(m, L)|) = 1 - 1 / (gamma |cal(A)_(m, L)|) sum_A chevron.l A, Delta chevron.r + O(gamma^(-2)). $
  Applying $-gamma log(1 - u / gamma + O(gamma^(-2))) = u + O(gamma^(-1))$ yields the stated limit. Since $bar(A)_i$ is strictly positive row-wise with $(bar(A)_i bold(1)_(L_i))_p >= 1$, the limiting quadratic is strictly convex with the stated closed-form minimizer.
]

Empirical results from the ClusterOverview four-dataset benchmark highlight both the clustering strength and computational trade-off of Soft-DTW Clustering @rousseeuw1987silhouettes @hubert1985comparing. Across the four suites, Soft-DTW achieves Iris $"Sil" = 0.442$ ($"ARI" = 0.661$), Customers $"Sil" = 0.322$, Multi-Shape $"Sil" = 0.603$ ($"ARI" = "N/A"$), and High-D $"Sil" = 0.275$ ($"ARI" = "N/A"$), at a mean runtime of $147.236$ ms. Notably, Soft-DTW attains the *highest Iris ARI* ($0.661$ vs. $0.608$ for DBA, $0.515$ for ROCKET-Clustering, and $0.397$ for k-Shape) and the *highest High-D Silhouette* ($0.275$) among all four time-series algorithms (Algs. 162--165), while jumping dramatically over DBA on Multi-Shape ($0.603$ vs. $0.459$). This empirical gain directly reflects Theorem 3: entropic smoothing over competing alignments eliminates brittle hard-path switching and smooths away suboptimal local minima. The price is a $5.6 times$ runtime overhead relative to DBA ($26.159$ ms), $10 times$ relative to k-Shape ($14.685$ ms), and $135 times$ relative to ROCKET-Clustering ($1.091$ ms).
