#import "../template/ieee.typ": *

#alg_entry(
  id: 38,
  name: "DBCLASD (Distribution-Based Clustering)",
  category: "Cat. II: Density & Mode-Seeking",
  year: "1998",
  paradigm: "Nonparametric Poisson NN-Distance Testing",
  complexity: $O(n log n)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

A central operational weakness of density-based clustering algorithms such as DBSCAN @ester1996dbscan is their acute sensitivity to the user-supplied neighborhood radius $epsilon$ and minimum density count $"MinPts"$, whose optimal values vary drastically across datasets and scales. *DBCLASD (Distribution-Based Clustering of Large Spatial Databases)* @xu1998dbclasd eliminates all user-tuned density parameters by defining clusters *statistically* rather than via fixed geometric thresholds.

Let $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$ be a finite spatial dataset (originally formulated for $d = 2$ and directly generalizable to $d >= 2$). DBCLASD posits that points inside a true cluster $C subset.eq cal(X)$ occupying an arbitrary connected region $R subset.eq bb(R)^d$ of Lebesgue volume $A = "Vol"(R) > 0$ are distributed *uniformly* over $R$—equivalently, as a finite realization of a *homogeneous spatial Poisson point process* of unknown intensity $lambda = |C| / A > 0$. For any subset $C subset.eq cal(X)$ with $|C| >= 2$ and any point $x in C$, define its intra-cluster nearest-neighbor distance:
$ "NNDist"_C (x) = min_(y in C backslash {x}) norm(x - y)_2, $
and let $"NNDistSet"(C) = {"NNDist"_C (x) : x in C}$ denote the empirical multiset of nearest-neighbor distances within $C$.

#definition(title: "Distribution-Based Cluster (DBCLASD)")[
  Let $F_(N, A, d) (r)$ be the theoretical cumulative distribution function of the nearest-neighbor distance for $N = |C|$ points uniformly distributed in a region of volume $A$ in $bb(R)^d$, and let $hat(A)(C)$ be an empirical grid-based estimate of $"Vol"(R)$. A non-empty subset $C subset.eq cal(X)$ is a *distribution-based cluster* at significance level $alpha in (0, 1)$ if it satisfies three conditions @xu1998dbclasd:
  1. *Expected-Distribution Condition:* The empirical multiset $"NNDistSet"(C)$ passes a Pearson $chi^2$ goodness-of-fit test against $F_(|C|, hat(A)(C), d)$ at level $alpha$, i.e., $chi^2 (C) <= chi^2_(m-1, 1 - alpha)$ across $m$ distance bins.
  2. *Maximality:* Adding any unclassified neighboring point $p in cal(X) backslash C$ within the adaptive retrieval horizon $r_"max" (C)$ violates Condition 1.
  3. *Connectivity:* For every cell $g_1, g_2$ in the grid covering $C$, there exists a sequence of face/vertex-adjacent occupied grid cells connecting $g_1$ to $g_2$.
]

To approximate $A = "Vol"(R)$ for arbitrary non-convex shapes without parametric assumptions, DBCLASD overlays a regular axis-aligned grid of side length $g$ on $C$ and estimates $hat(A)_g (C) = |cal(G)_g (C)| dot g^d$, where $cal(G)_g (C)$ is the set of cells containing at least one point of $C$. Candidate points are fetched around newly accepted cluster members via spatial range queries of adaptive radius $r_"max" (N, hat(A), d)$ computed directly from the null tail probability (Algorithm 1). To mitigate order dependence, rejected candidates are retained in an `Unsuccessful` queue and re-tested as the cluster grows, and points may be reassigned across clusters if statistical fit changes.

#algbox(
  title: "DBCLASD Incremental Region-Growing Algorithm",
  inputs: [Spatial dataset $cal(X) subset.eq bb(R)^d$ indexed by an R\*-tree; significance level $alpha$ (default $0.05$, fixed internally).],
  outputs: [Cluster set $cal(C) = {C_1, dots, C_K}$ and noise set $cal(N)$.],
)[
  1. Mark all $x in cal(X)$ as unclassified; set $cal(C) <- emptyset$.
  2. *For each* unclassified seed $x in cal(X)$ *do*:
  3. #h(1em) Initialize candidate cluster $C$ with $x$ and its $N_"init"$ (e.g., $30$) nearest neighbors to form a statistically testable seed; initialize queues $"Cand" <- emptyset$, $"Unsucc" <- emptyset$.
  4. #h(1em) Estimate volume $hat(A)(C)$ via occupied grid cells; compute adaptive search radius $r_"max" (|C|, hat(A)(C), d)$.
  5. #h(1em) Enqueue all points within $r_"max"$ of $C$ into $"Cand"$.
  6. #h(1em) *While* $"Cand" != emptyset$ *do*:
  7. #h(2em) Pop $p$ from $"Cand"$; update grid volume $hat(A)(C union {p})$ and incremental $"NNDistSet"(C union {p})$.
  8. #h(2em) Evaluate Pearson statistic $chi^2 (C union {p})$ against $F_(|C|+1, hat(A), d)$.
  9. #h(2em) *If* $chi^2 (C union {p}) <= tau_(m, alpha)$ *then*:
  10. #h(3em) Accept $C <- C union {p}$; update $r_"max"$; retrieve new $r_"max"$-neighbors of $p$ into $"Cand"$; re-enqueue $"Unsucc"$ into $"Cand"$ and clear $"Unsucc"$.
  11. #h(2em) *Else*: move $p$ to $"Unsucc"$.
  12. #h(1em) *If* $|C| >= N_"init"$ and $C$ passes $chi^2$ test *then* add $C$ to $cal(C)$.
  13. Assign remaining unclassified points to noise $cal(N)$; *return* $(cal(C), cal(N))$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five formal theoretical results for DBCLASD: (1) exact finite-sample and asymptotic Weibull nearest-neighbor distance laws in $bb(R)^d$; (2) closed-form derivation and asymptotic scaling of the parameter-free retrieval radius $r_"max"$; (3) geometric boundary bias and exponential high-dimensional breakdown of grid volume estimation; (4) an explicit counterexample proving asymptotic rejection of unimodal Gaussian clusters; and (5) classical goodness-of-fit consistency bounds.

#subsec("A", "Exact NN-Distance Law and Poisson–Weibull Limit")

#theorem(title: "Exact Finite-Sample NN Law and Uniform Weibull Convergence", origin: "Original Derivation")[
  Let $X_1, dots, X_N$ be i.i.d. uniformly distributed on a measurable region $R subset.eq bb(R)^d$ of volume $A = "Vol"(R) > 0$, and let $V_d = pi^(d/2) / Gamma(d/2 + 1)$ be the volume of the unit Euclidean ball in $bb(R)^d$. For an interior point $x in R$ and radius $r >= 0$ with $B_r (x) subset.eq R$:
  1. The exact nearest-neighbor distance CDF is
     $ F_(N, A, d) (r) = 1 - (1 - (V_d r^d) / A)^(N - 1). $
  2. As $N, A -> infinity$ with $N / A -> lambda > 0$, $F_(N, A, d) (r)$ converges uniformly on $[0, infinity)$ to the Weibull CDF $F_lambda (r) = 1 - exp(-lambda V_d r^d)$, with finite-sample error bounded by
     $ sup_(r >= 0) |F_(N, A, d) (r) - F_lambda (r)| <= (1 + 2 e^(-1)) / N = O(1 / N) $
     (where $F_(N, A, d) (r) := 1$ for $V_d r^d >= A$).
  3. Under $F_lambda$, the $k$-th moment equals $bb(E)[R_1^k] = (lambda V_d)^(-k / d) Gamma(1 + k / d)$; for $d = 2$, $R_1$ is Rayleigh with $bb(E)[R_1] = (2 sqrt(lambda))^(-1)$ and $"Var"(R_1) = (4 - pi) / (4 pi lambda)$.
]

#proof[
  (1) Conditioning on $X_1 = x$, the remaining $N - 1$ points are independent and uniform on $R$. Because $B_r (x) subset.eq R$, each $X_j$ ($j >= 2$) falls inside $B_r (x)$ with exact probability $"Vol"(B_r (x)) / "Vol"(R) = V_d r^d / A$. Hence the event ${R_1 > r} = inter.big_(j=2)^N {X_j in.not B_r (x)}$ has probability $(1 - V_d r^d / A)^(N - 1)$, giving Claim 1.

  (2) Set $u = lambda V_d r^d >= 0$ with $lambda = N / A$, so $V_d r^d / A = u / N$. For $u >= N$, $F_(N, A, d) = 1$ and $1 - F_lambda = e^(-u) <= e^(-N) < 1 / N$. For $u in [0, N]$, note that $(1 - u / N)^N <= e^(-u)$ and simultaneously $(1 - u / N)^(N - 1) = (1 - u / N)^(-1) (1 - u / N)^N$. Using $0 <= e^(-u) - (1 - u / N)^N <= (u^2 / N) e^(-u)$ (from $ln(1 - t) >= -t - t^2 / (2(1-t))$ and convexity), together with $sup_(u >= 0) u^2 e^(-u) = 4 e^(-2)$ and $(1 - u/N)^(N-1) - (1 - u/N)^N = (u / N)(1 - u/N)^(N-1) <= (u / N) e^(-u (N-1)/N)$, we bound
  $ |(1 - u / N)^(N - 1) - e^(-u)| &<= |e^(-u) - (1 - u / N)^N| \
    &quad + (u / N) (1 - u / N)^(N - 1). $
  For $N >= 2$, the first term is at most $4 e^(-2) / N$ and the second attains its maximum over $u in [0, N]$ at $u = 1$, giving $(1/N)(1 - 1/N)^(N-1) <= 1 / N$ (in fact $-> e^(-1) / N$). Summing yields uniform error $<= (1 + 4 e^(-2)) / N < (1 + 2 e^(-1)) / N$.

  (3) Substituting $w = lambda V_d r^d$ so $r = (w / (lambda V_d))^(1/d)$ and $dif F_lambda (r) = e^(-w) dif w$ yields
  $ bb(E)[R_1^k] = integral_0^infinity (w / (lambda V_d))^(k / d) e^(-w) dif w = (lambda V_d)^(-k / d) Gamma(1 + k / d). $
  Setting $d = 2$ ($V_2 = pi$), $k = 1$ gives $Gamma(3/2) / sqrt(pi lambda) = (sqrt(pi) / 2) / sqrt(pi lambda) = (2 sqrt(lambda))^(-1)$, and $k = 2$ gives $bb(E)[R_1^2] = 1 / (pi lambda)$, so $"Var"(R_1) = 1 / (pi lambda) - 1 / (4 lambda) = (4 - pi) / (4 pi lambda)$.
]

#subsec("B", "Exact Parameter-Free Candidate Retrieval Radius")

#lemma(title: "Closed-Form Retrieval Radius and Extreme-Tail Bound", origin: "Original Proof")[
  Define DBCLASD's adaptive retrieval radius $r_"max" (N, A, d)$ as the unique root of the expected extreme-tail count equation $N (1 - F_(N, A, d) (r)) = 1$ for $N >= 2$ @xu1998dbclasd. Then:
  1. $r_"max" (N, A, d) = ((A / V_d) (1 - N^(-1 / (N - 1))))^(1 / d)$ (recovering $sqrt((A / pi)(1 - N^(-1/(N-1))))$ for $d = 2$).
  2. For fixed density $lambda = N / A$ and $N -> infinity$,
     $ r_"max" (N, A, d) = ((ln N) / (lambda V_d))^(1 / d) (1 + O((ln N) / N)). $
  3. Under the homogeneous Poisson model $F_lambda$, for any multiplier $c > 1$, the probability that any interior point of $C$ has nearest neighbor exceeding $c dot ((ln N) / (lambda V_d))^(1/d)$ is at most $N^(1 - c^d) -> 0$.
]

#proof[
  (1) Solving $N (1 - V_d r^d / A)^(N - 1) = 1$ gives $(1 - V_d r^d / A)^(N - 1) = N^(-1)$, hence $1 - V_d r^d / A = N^(-1 / (N - 1))$, and solving for $r > 0$ immediately yields the closed-form identity.

  (2) Rewrite $1 - N^(-1 / (N - 1)) = 1 - exp(-(ln N) / (N - 1))$. Using Taylor's theorem $1 - e^(-t) = t - O(t^2)$ with $t = (ln N) / (N - 1) = (ln N / N)(1 + O(1 / N))$, we obtain
  $ (A / V_d) (1 - N^(-1 / (N - 1))) = (N / (lambda V_d)) ((ln N) / N) (1 + O((ln N) / N)) \ = ((ln N) / (lambda V_d)) (1 + O((ln N) / N)). $
  Raising both sides to the $1 / d$ power establishes Claim 2.

  (3) Setting $r_c = c ((ln N) / (lambda V_d))^(1/d)$, each interior point satisfies $bb(P)(R_1 > r_c) = exp(-lambda V_d r_c^d) = exp(-c^d ln N) = N^(-c^d)$. Applying Boole's union bound over all $N$ points gives $bb(P)(max_(i <= N) R_(1, i) > r_c) <= N dot N^(-c^d) = N^(1 - c^d) -> 0$ as $N -> infinity$ for every $c > 1$. Thus $r_"max"$ is asymptotically the exact percolation/connectivity radius needed to traverse a homogeneous cluster without user-supplied $epsilon$.
]

#subsec("C", "Edge-Effect Bias and High-Dimensional Grid Breakdown")

#theorem(title: "Boundary Bias and Exponential Dimensional Distortion", origin: "Original Proof")[
  1. *Boundary Stochastic Dominance:* Let $x in partial R$ lie on a smooth flat boundary of a homogeneous Poisson cluster $R subset.eq bb(R)^d$ of intensity $lambda$. For small $r > 0$, $R_1 (x)$ obeys $F_"bnd" (r) = 1 - exp(-(1/2) lambda V_d r^d) < F_"int" (r)$ for all $r > 0$, inflating expected NN distance by a factor $2^(1 / d)$ near $partial R$.
  2. *Grid Volume Explosion in High $d$:* Cover $R = B_rho (0) subset.eq bb(R)^d$ by a uniformly shifted regular grid $g bb(Z)^d + u$ ($u ~ "Uniform"([0, g)^d)$) of cell side $g > 0$. In the full-coverage regime where every intersecting cell is occupied, the expected grid volume estimator $hat(A)_(g, u) (R)$ satisfies
     $ (bb(E)_u [hat(A)_(g, u) (R)]) / ("Vol"(R)) &>= (1 + g / (V_d^(1/d) rho))^d \
       &>= (1 + (g sqrt(d)) / (rho sqrt(2 pi e)))^d. $
     Maintaining bounded relative volume error forces $g = O(rho / d^(3/2))$, at which scale expected points per cell scale as $O(N d^(-d)) -> 0$, collapsing $hat(A)_g (C) -> |C| g^d$ and yielding a vacuous constant density $hat(lambda) equiv g^(-d)$.
]

#proof[
  (1) Near a $C^2$ boundary point $x in partial R$, the tangent hyperplane bisects $B_r (x)$, so $"Vol"(B_r (x) inter R) = (1/2) V_d r^d + O(r^(d+1))$. Under a homogeneous Poisson process of intensity $lambda$, $bb(P)(R_1 (x) > r) = exp(-lambda "Vol"(B_r (x) inter R)) = exp(-(1/2) lambda V_d r^d + O(r^(d+1)))$, strictly exceeding $exp(-lambda V_d r^d)$. Thus boundary points systematically inflate the empirical right tail of $"NNDistSet"(C)$, elevating $chi^2 (C)$ on elongated or high-perimeter shapes.

  (2) For shift $u in [0, g)^d$, the occupied-cell volume is
  $ hat(A)_(g, u) (R) = g^d sum_(k in bb(Z)^d) bb(I)((u + g k + [0, g]^d) inter B_rho (0) != emptyset). $
  Averaging $u$ over $[0, g)^d$ tiles all of $bb(R)^d$ by $x = u + g k$, giving the exact integral-geometric identity
  $ bb(E)_u [hat(A)_(g, u) (R)] &= integral_(bb(R)^d) bb(I)((x + [0, g]^d) inter B_rho (0) != emptyset) dif x \
    &= "Vol"(B_rho (0) + [-g, 0]^d) = "Vol"(B_rho (0) + [0, g]^d). $
  Applying the Brunn–Minkowski inequality to the convex bodies $B_rho (0)$ and $[0, g]^d$ yields
  $ "Vol"(B_rho (0) + [0, g]^d)^(1/d) >= "Vol"(B_rho (0))^(1/d) + g = V_d^(1/d) rho + g. $
  Dividing by $V_d^(1/d) rho$ and raising to the $d$-th power yields the first inequality. By Stirling's bound $Gamma(d/2 + 1) >= sqrt(pi d) (d / (2 e))^(d/2)$, we have $V_d^(1/d) = (pi^(d/2) / Gamma(d/2 + 1))^(1/d) <= sqrt(2 pi e / d)$, establishing the second inequality. Consequently, bounding the multiplicative volume bias requires $d^(3/2) g / rho = O(1)$, i.e., $g = O(rho / d^(3/2))$. At that resolution, the ratio of cell volume to ball volume obeys
  $ g^d / (V_d rho^d) = O(d^(-3 d / 2) / (2 pi e / d)^(d/2)) = O(d^(-d)) -> 0, $
  so for any fixed sample size $N$, expected occupancy per cell $N g^d / "Vol"(R) -> 0$ super-exponentially. Almost every occupied cell then holds a single point, forcing $|cal(G)_g (C)| = |C|$ and $hat(A)_g (C) = |C| g^d$, which renders $hat(lambda) = |C| / hat(A)_g (C) equiv g^(-d)$ completely blind to true spatial geometry.
]

#subsec("D", "Failure Mode Under Unimodal Gaussian Densities")

#proposition(title: "Asymptotic Rejection of Unimodal Gaussian Clusters", origin: "Original Counterexample")[
  Let $X_1, dots, X_N$ be i.i.d. draws from a single isotropic 2D Gaussian cluster $cal(N)(0, sigma^2 I_2)$. Under the local homogeneous-Poisson intensity $lambda(x) = N f(x) = (N / (2 pi sigma^2)) exp(-norm(x)_2^2 / (2 sigma^2))$, the marginal nearest-neighbor distance $R_1$ obeys the exact heavy-tailed law
  $ bb(P)(R_1 > r) = (1 - e^(-a)) / a, quad a := (N r^2) / (2 sigma^2), $
  whose tail decays algebraically as $2 sigma^2 / (N r^2)$ with infinite second moment ($bb(E)[R_1^2] = infinity$). Consequently, Pearson's $chi^2$ test comparing $"NNDistSet"(C)$ against *any* single homogeneous Rayleigh law $F_lambda (r) = 1 - exp(-pi lambda r^2)$ diverges as $chi^2 (C) = Omega_p (N)$, guaranteeing asymptotic rejection ($bb(P)("reject") -> 1$) as $N -> infinity$.
]

#proof[
  Conditional on radial norm $rho = norm(X_1)_2$, the local Poisson survival probability is $bb(P)(R_1 > r | rho) = exp(-pi lambda(rho) r^2) = exp(-a e^(-rho^2 / (2 sigma^2)))$ with $a = N r^2 / (2 sigma^2)$. Because $X_1 ~ cal(N)(0, sigma^2 I_2)$, the squared radius $rho^2 / (2 sigma^2)$ follows $"Exp"(1)$, hence $U := exp(-rho^2 / (2 sigma^2))$ is strictly $"Uniform"(0, 1)$. Marginalizing over $U$ gives the exact closed form
  $ bb(P)(R_1 > r) = integral_0^1 exp(-a u) dif u = (1 - e^(-a)) / a. $
  As $r -> infinity$, $(1 - e^(-a)) / a ~ 2 sigma^2 / (N r^2)$, so $bb(E)[R_1^2] = integral_0^infinity bb(P)(R_1^2 > t) dif t = +infinity$, whereas every homogeneous Rayleigh null model $F_lambda$ has finite second moment $(pi lambda)^(-1)$ and super-exponential tail $exp(-pi lambda r^2)$. Partitioning rescaled distances $sqrt(N) R_1$ into $m >= 3$ fixed intervals yields a true bin vector $p^*$ with strictly positive $L_1$ distance $delta_1 := inf_(lambda > 0) norm(p^* - p(lambda))_1 > 0$ from the 1-parameter Rayleigh family ${p(lambda)}$. By Cauchy–Schwarz (since $sum_(j=1)^m p_j (hat(lambda)) = 1$),
  $ chi^2 (C) &= N sum_(j=1)^m (hat(p)_j - p_j (hat(lambda)))^2 / (p_j (hat(lambda))) \
    &>= N (sum_(j=1)^m |hat(p)_j - p_j (hat(lambda))|)^2 = N norm(hat(p) - p(hat(lambda)))_1^2. $
  By consistency of empirical bin frequencies $hat(p) -> p^*$ in probability, $chi^2 (C) >= N (delta_1 - o_p (1))^2 = Omega_p (N) -> infinity$, exceeding any fixed critical value $chi^2_(m-1, 1 - alpha)$ with probability approaching $1$. Thus DBCLASD inevitably fragments Gaussian blobs or discards their low-density fringes as noise.
]

#subsec("E", "Goodness-of-Fit Calibration and Statistical Power")

#theorem(title: "Asymptotic Calibration and Consistency of the Goodness-of-Fit Test", origin: "Classical — Cited")[
  Let $C$ be a realization of size $N$ from a spatial point process on $R subset.eq bb(R)^d$ partitioned into $m$ fixed quantile intervals of the reference distribution $F_(N, hat(A), d)$ @xu1998dbclasd @anderson1952asymptotic. Under independent distance sampling (or negligible boundary/pair-correlation dependence) and consistent volume estimation $hat(A) -> "Vol"(R)$:
  1. Under the homogeneous Poisson null $H_0$, Pearson's statistic converges in distribution to $chi^2_(m - 1)$ (or $chi^2_(m - 2)$ when intensity $lambda$ is estimated from the same sample), controlling false cluster splitting at level $alpha$.
  2. Under any smooth inhomogeneous alternative density $f$ whose induced NN-distance distribution differs from ${F_lambda}_(lambda > 0)$ on the $m$ bins, the test is consistent with power $bb(P)(chi^2 (C) > tau_(m, alpha)) -> 1$ as $N -> infinity$.
]

#proof(label: "Proof Synthesis (Cited — Xu et al. 1998; Anderson & Darling 1952)")[
  Classical goodness-of-fit theory @anderson1952asymptotic shows that for $N$ observations partitioned into $m$ bins with expected probabilities $p_j > 0$, the standardized count vector $Z_j = (O_j - N p_j) / sqrt(N p_j)$ converges by the multivariate CLT to a singular Gaussian with covariance $I_m - sqrt(p) sqrt(p)^T$, whose squared norm $sum_j Z_j^2$ follows $chi^2_(m-1)$ (losing one degree of freedom if $lambda$ is estimated via minimum-$chi^2$ or MLE). In DBCLASD @xu1998dbclasd, because mutual nearest-neighbor pairs introduce mild positive correlation between $"NNDist"_C (x)$ and $"NNDist"_C (y)$ (each pair shares the same distance whenever $x$ and $y$ are reciprocal 1-NNs), the empirical variance of $O_j$ is slightly inflated relative to pure i.i.d. multinomial counts; Xu et al. compensate practically by employing a conservative threshold (e.g., $alpha = 0.05$ with $N_"init" >= 30$) while preserving consistency against inhomogeneous alternatives (as explicitly demonstrated in Proposition 1).
]

#sec("III", "Complexity, Geometry and Operational Profile")

With an R\*-tree index in $d = 2$–$3$, each point triggers an $r_"max"$-range query costing $O(log n + k_"cand")$, plus an incremental update to the grid cell hash table and $"NNDistSet"(C)$ (updating only existing cluster points whose previous NN distance exceeds their distance to the newly added point $p$). Because unsuccessful candidates may be re-tested multiple times as a cluster expands, Xu et al. @xu1998dbclasd observe empirically that DBCLASD executes in roughly $3 times$ the runtime of DBSCAN—still $O(n log n)$ with spatial indexing, or $O(n^2)$ without an index—using $O(n d)$ memory.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for DBCLASD],
    [Time / Space Complexity], [$O(n log n)$ indexed ($approx 3 times$ DBSCAN), $O(n^2)$ unindexed; $O(n d)$ memory.],
    [Hyperparameter Burden], [Strictly *parameter-free* (no $epsilon$, $"MinPts"$, or $k$; $r_"max"$ derived from $F_(N, A, d)$).],
    [Geometric Inductive Bias], [Arbitrary non-convex shapes with *locally uniform* (homogeneous Poisson) point density.],
    [Determinism / Order], [Mitigated via `Unsuccessful` re-testing and reassignment, though minor seed sensitivity persists.],
    [Failure Regimes], [Non-uniform Gaussian/peaked clusters (Prop. 1), thin high-perimeter arcs & $d >= 4$ (Theorem 2).],
    [Recommended Upgrades], [HDBSCAN (Alg. 32) or OPTICS (Alg. 33) for non-uniform densities; GDBSCAN (Alg. 39) for weighted objects.],
  ),
  caption: [Summary of theoretical properties and operational profile of DBCLASD.],
)

#sec("IV", "Engineering and Scientific Applications")

DBCLASD shines in low-dimensional spatial exploration where clusters have distinct, unknown scales or densities—so no single global DBSCAN $epsilon$ works—yet points within each cluster follow approximately uniform spatial sampling.

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Exploratory Geographic Information Systems (GIS) & Land-Use Zoning:* Parcel centroids, agricultural plots, and uniformly sampled urban building footprints in 2D exhibit roughly constant spatial density inside each zone while varying dramatically between dense downtown cores and sparse suburbs. Because $r_"max"$ adapts automatically to each cluster's estimated intensity $hat(lambda) = N / hat(A)$, DBCLASD discovers multi-density zones without manual $epsilon$ tuning.
2. *2D Astronomical & Ecological Quadrat Surveys:* Uniformly distributed stellar associations, forest-canopy tree maps, and benthic species surveys naturally match the spatial Poisson null model (Complete Spatial Randomness, CSR) inside homogeneous habitats, allowing DBCLASD to separate habitat patches from sparse background clutter.
3. *Semiconductor Wafer Defect & Surface Metrology Mapping:* Uniformly clustered scratch bands, particle clouds, and etch-defect regions on 2D silicon wafers occur at unpredictable densities across fabrication lots; parameter-free statistical acceptance eliminates fragile per-lot threshold calibration.
4. *Retail & Customer Footprint Analysis (Low-Dimensional Segments):* When two-dimensional customer coordinates (e.g., normalized income vs. spending score) occupy roughly uniform polygonal regions of differing area, DBCLASD adapts smoothly across cluster scales without forcing equal-radius hyperspheres.
5. *Automated Spatial Database Pre-Indexing:* As a zero-configuration first-pass partitioner over legacy 2D R\*-tree tables, DBCLASD identifies coherent spatial blocks and isolates background noise before downstream parametric modeling.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

DBCLASD exhibits three sharp, mathematically characterized failure modes: (a) *unimodal peaked (Gaussian or Student-$t$) clusters* fail the homogeneous Poisson test and are fragmented or stripped of their tails (Proposition 1); (b) *narrow filaments or thin curved arcs* suffer severe edge-effect NN-distance inflation (Theorem 2(1)); and (c) *moderate-to-high dimensions ($d >= 4$)* break the grid-based volume estimator $hat(A)_g (C)$ exponentially (Theorem 2(2)).

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Why DBCLASD Fits / Prerequisite], [Caution / Alternative],
    [2D GIS / Urban Zones], [Zero tuning; adapts to multi-density areas], [Use HDBSCAN (Alg. 32) for peaked density],
    [Wafer / Ecology Maps], [Exact fit to homogeneous CSR patches], [Use OPTICS (Alg. 33) for thin filaments],
    [2D Customer Segments], [Adapts to varying uniform segment areas], [Use GMM (Alg. 81) if strongly Gaussian],
    [Thin Non-Convex Arcs], [Edge-effect bias distorts NN distribution], [Use DBSCAN (Alg. 31) or GDBSCAN (Alg. 39)],
    [High-D ($d >= 4$) Data], [Grid volume $hat(A)_g$ collapses (Theorem 2)], [Use SNN (Alg. 36) or RNN-DBSCAN (Alg. 45)],
  ),
  caption: [Application suitability and failure-mode matrix for DBCLASD.],
)

#sec("V", "Structural Relations and Empirical Behavior")

DBCLASD occupies a unique bridge between nonparametric density clustering and formal spatial hypothesis testing.

#proposition(title: "DBCLASD as Statistically-Gated Adaptive-Radius DBSCAN", origin: "Original Derivation")[
  At any expansion step of Algorithm 1 on a cluster $C$ of size $N$ and estimated intensity $hat(lambda) = N / hat(A)_g (C)$, define the dynamic neighborhood predicate $"NPred"_C (x, p) := [norm(x - p)_2 <= r_"max" (N, hat(A)_g (C), d)]$ and the non-additive acceptance functional $Phi_C (p) := bb(I)[chi^2 (C union {p}) <= tau_(m, alpha)]$. Then DBCLASD's region-growing transition coincides with an *adaptive-radius* generalization of DBSCAN (Alg. 31) @ester1996dbscan where the static local cardinality threshold $|N_epsilon (p)| >= "MinPts"$ is replaced by the global distributional predicate $Phi_C (p) = 1$, and $r_"max" approx ((ln N) / (hat(lambda) V_d))^(1/d)$ scales automatically with local inverse density $hat(lambda)^(-1/d)$—linking DBCLASD directly to VDBSCAN (Alg. 50) @liu2007vdbscan and GDBSCAN (Alg. 39) @sander1998gdbscan.
]

#proof[
  In DBSCAN, cluster $C$ expands from $x in C$ to candidate $p$ iff $d(x, p) <= epsilon$ and $p$ satisfies the density criterion. In Algorithm 1, candidates around $x in C$ are retrieved iff $norm(x - p)_2 <= r_"max" (N, hat(A)_g (C), d)$ and accepted into $C$ iff $Phi_C (p) = 1$. By Lemma 1(2), $r_"max" (N, hat(A)_g (C), d) = ((ln N) / (hat(lambda) V_d))^(1/d) (1 + o(1))$, so the search ball automatically contracts in dense regions ($hat(lambda)$ large) and expands in sparse regions ($hat(lambda)$ small), mirroring multi-$epsilon$ adaptation in VDBSCAN @liu2007vdbscan while replacing fixed local counts with a global distributional invariant. Unlike GDBSCAN's fixed symmetric predicate, however, both $"NPred"_C$ and $Phi_C$ depend dynamically on the evolving cluster $C$, explaining why DBCLASD requires deferred candidate re-evaluation (`Unsuccessful`) rather than a single static graph traversal.
]

Empirical results on the ClusterOverview 4-dataset benchmark suite provide striking quantitative validation of Theorems 1–2 and Proposition 1. On the 2D *Customers* dataset, DBCLASD attains a strong Silhouette of $0.414$—decisively outperforming fixed-radius PDBSCAN ($0.296$)—precisely because adaptive $r_"max"$ accommodates segments of differing spatial densities without a single global $epsilon$. Conversely, on 4D *Iris*, DBCLASD collapses to Silhouette $= 0.060$ and ARI $= 0.314$, directly confirming both Proposition 1 (the unimodal Gaussian-like Iris species fail the homogeneous Poisson test) and Theorem 2(2) (4D grid-volume distortion on small $N approx 50$ classes). On *Multi-Shape* (Silhouette $= 0.334$, ARI N/A), thin curved rings incur severe perimeter-to-area edge-effect bias (Theorem 2(1)), trailing PDBSCAN ($0.717$) and GDBSCAN ($0.577$). On *High-D Subspace+Noise* (Silhouette $= 0.317$, ARI N/A), grid-volume collapse limits separation, though DBCLASD achieves the fastest mean runtime of the three algorithms at $0.381$ ms thanks to lightweight candidate pruning on small sample sizes.
