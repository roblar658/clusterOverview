#import "../template/ieee.typ": *

#alg_entry(
  id: 25,
  name: "CLARA-MED (Median-Robust CLARA)",
  category: "Cat. I: Centroid & Partitioning",
  year: "1990 / Ext.",
  paradigm: "Subsampled Medoid Search with Median-Based Objective Scoring",
  complexity: $O(q (s^2 T + n k))$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Classical PAM (Alg. 4 @kaufman1990finding) and even $O(k)$-accelerated FastPAM/FasterPAM (Alg. 24 @schubert2019fasterpam @schubert2021fast) require $Omega(n^2)$ time and distance evaluations, preventing application to large datasets ($n >= 10^5$). *CLARA* (Clustering Large Applications, Kaufman and Rousseeuw @srikantha2011clara @kaufman1990finding) achieves $O(n)$ scalability by drawing $q$ uniform subsamples $S_1, dots, S_q subset.eq X$ of fixed size $s << n$, computing a candidate medoid set $M^((r)) subset.eq S_r$ ($|M^((r))| = k$) on each subsample, and selecting the candidate minimizing the full-dataset *mean* nearest-medoid dissimilarity:
$ overline(D)(M) := 1 / n sum_(i=1)^n d(x_i, M), quad d(x, M) := min_(m in M) d(x, m). $

Although each individual medoid is a real data object, CLARA's outer model-selection functional $overline(D)(M)$ is an *unbounded arithmetic mean*. Consequently, even a single gross outlier $x_"out"$ at extreme distance $R -> infinity$ dominates $overline(D)(M)$, biasing CLARA to select a corrupted subsample that sacrifices a cluster medoid to absorb $x_"out"$ rather than a clean medoid configuration. *CLARA-MED* replaces mean-dissimilarity evaluation with a *median-based robust selection objective*, immunizing outer candidate selection against up to nearly $50\%$ extreme contamination while preserving linear $O(n)$ scaling.

#definition(title: "Median-Dissimilarity Objective & CLARA-MED Estimator")[
  Let $(X, d)$ be a finite dissimilarity space of size $n$. For any $k$-medoid set $M subset.eq X$, let $d_((1)) (M) <= d_((2)) (M) <= dots <= d_((n)) (M)$ denote the sorted nearest-medoid distances ${d(x_i, M)}_(i=1)^n$. The *CLARA-MED objective* (Least Median of Distances) is the sample median
  $ J_"med" (M) := 1 / 2 (d_((floor((n + 1)/2))) (M) + d_((ceil((n + 1)/2))) (M)). $
  To break ties smoothly without forfeiting robustness, CLARA-MED pairs $J_"med" (M)$ with the secondary *lower-half trimmed mean* $J_(1/2) (M) := ceil(n / 2)^(-1) sum_(r=1)^(ceil(n / 2)) d_((r)) (M)$. Given $q$ subsampled candidate medoid sets ${M^((r))}_(r=1)^q$, CLARA-MED returns $hat(M) = arg min_(1 <= r <= q) (J_"med" (M^((r))), J_(1/2) (M^((r))))$ lexicographically.
]

#algbox(
  title: "CLARA-MED (Median-Robust Subsampled Medoid Clustering)",
  inputs: [Dataset $X = {x_1, dots, x_n}$ with dissimilarity $d$, cluster count $k >= 1$, subsample count $q >= 1$, subsample size $s in [k, n]$ (recommended $s = Theta(k log k)$).],
  outputs: [Robust $k$-medoid set $hat(M) subset.eq X$ and full-dataset partition $cal(C)$.],
)[
  1. Initialize $(J^*, T^*, hat(M)) <- (+infinity, +infinity, emptyset)$.
  2. *For* $r = 1, dots, q$ *do*:
  3. #h(1em) Draw a uniform random subsample $S_r subset.eq X$ of size $s$ without replacement.
  4. #h(1em) Run FastPAM1 / FasterPAM (Alg. 24) or PAM (Alg. 4) on $(S_r, d)$ to obtain $k$ candidate medoids $M^((r)) subset.eq S_r$.
  5. #h(1em) *Full-dataset evaluation:* For each $i in {1, dots, n}$, compute $d_i^((r)) <- min_(m in M^((r))) d(x_i, m)$ and nearest index $c_i^((r))$.
  6. #h(1em) Compute $J_r <- J_"med" (M^((r)))$ and $T_r <- J_(1/2) (M^((r)))$ in $O(n)$ time via linear-time order-statistic selection (Quickselect / introselect).
  7. #h(1em) *If* $J_r < J^*$ or ($J_r = J^*$ and $T_r < T^*$) *then*:
  8. #h(2em) Set $(J^*, T^*, hat(M)) <- (J_r, T_r, M^((r)))$ and cache assignments $hat(c) <- c^((r))$.
  9. *Return* $hat(M)$ and clusters $C_j = {i in {1, dots, n} : hat(c)_i = j}$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five formal properties governing CLARA-MED: (1) its exact finite-sample $50\%$ scoring breakdown point compared to CLARA's $1/n$ collapse; (2) exact without-replacement subsample coverage bounds proving that classical $s = 40 + 2 k$ fails for large $k$ and $s = Omega(k log k)$ is required; (3) an explicit counterexample proving minority-cluster masking under global median scoring; (4) a moment-free $O(sqrt((log q)/(n - s)))$ finite-sample concentration bound; and (5) strong asymptotic consistency of trimmed/median quantizers.

#subsec("A", "Exact 50% Breakdown Point vs. Classical CLARA")

#theorem(title: "Exact Scoring Breakdown Point of CLARA-MED", origin: "Original Proof")[
  Let $(X, d)$ be a metric space consisting of $n - b$ clean points $X_"clean"$ contained in a ball of diameter $D_0 < infinity$ and $b$ adversarial outliers $X_"out"$ moved to distance $R -> infinity$ from $X_"clean"$. Assume $k >= 2$ and let $M_"clean" subset.eq X_"clean"$ ($|M_"clean"| = k$) be any clean medoid configuration. Then:
  1. For classical CLARA, $b = 1$ outlier suffices to drive $overline(D)(M_"clean") -> infinity$ while an outlier-hijacked set $M_"hijack"$ containing $x_"out"$ maintains $overline(D)(M_"hijack") <= D_0$, giving breakdown point $epsilon_"CLARA"^* = 1 / n$.
  2. For CLARA-MED, $J_"med" (M_"clean") <= D_0$ uniformly for all $R > 0$ if and only if $b <= ceil(n / 2) - 1$; conversely, $b = ceil(n / 2)$ outliers drive $J_"med" (M_"clean") -> infinity$. Hence the exact finite-sample breakdown point of CLARA-MED's selection functional is
     $ epsilon_"CLARA-MED"^* = (ceil(n / 2)) / n -> 1 / 2. $
]

#proof[
  (1) With $b = 1$, $d(x_"out", M_"clean") >= R$, so $overline(D)(M_"clean") >= R / n -> infinity$ as $R -> infinity$. Conversely, choosing $M_"hijack" = {x_"out"} union M'$ with $k - 1$ points $M' subset.eq X_"clean"$ gives $d(x_"out", M_"hijack") = 0$ and $d(x, M_"hijack") <= D_0$ for all $x in X_"clean"$, hence $overline(D)(M_"hijack") <= ((n - 1)/n) D_0 < R / n$ for $R > n D_0$. Thus classical CLARA prefers $M_"hijack"$ over $M_"clean"$ whenever $M_"hijack"$ appears in the candidate pool, proving $epsilon_"CLARA"^* = 1 / n$.

  (2) Let $u := ceil((n + 1) / 2) = floor(n / 2) + 1$ be the upper median order-statistic index in Definition 1. Notice $n - u = n - floor(n / 2) - 1 = ceil(n / 2) - 1$.
  - *If $b <= ceil(n / 2) - 1$:* Then the number of clean points satisfies $n - b >= u$. Every $x in X_"clean"$ has $d(x, M_"clean") <= D_0$, so at least $u$ sorted distances lie in $[0, D_0]$. Consequently both $d_((floor((n+1)/2))) (M_"clean") <= d_((u)) (M_"clean") <= D_0$, proving $J_"med" (M_"clean") <= D_0$ independently of $R$.
  - *If $b >= ceil(n / 2)$:* Then $n - b <= u - 1$, so at most $u - 1$ points lie in $X_"clean"$. Hence the $u$-th order statistic must belong to $X_"out"$, giving $d_((u)) (M_"clean") >= R$ and
    $ J_"med" (M_"clean") >= 1 / 2 d_((u)) (M_"clean") >= R / 2 -> infinity, $
    as $R -> infinity$. Therefore the smallest contaminating fraction causing unbounded explosion is exactly $ceil(n / 2) / n$.
]

#subsec("B", "Subsample Cluster Coverage and Necessary Sample Scaling")

Because medoids can only be chosen from $S_r$, CLARA-MED requires at least one subsample to simultaneously cover all $k$ true clusters and avoid outlier majority inside $S_r$.

#theorem(title: "Subsample Coverage & Necessity of s = Omega(k log k)", origin: "Original Derivation")[
  Suppose $X$ consists of $k$ disjoint clean clusters $C_1^*, dots, C_k^*$ with minimum relative mass $p_"min" := min_j |C_j^*| / n > 0$, plus an outlier fraction $epsilon = b / n < 1 / 2$. Let $S_r subset.eq X$ be drawn uniformly without replacement of size $s$, and let $cal(G)_r$ be the event that $S_r$ hits every $C_j^*$ and contains fewer than $s / 2$ outliers. Then:
  1. $bb(P)[cal(G)_r^c] <= k exp(- s p_"min") + exp(- 2 s (1/2 - epsilon)^2)$. In particular, choosing $s >= max(p_"min"^(-1) ln(4 k), (ln 4) / (2(1/2 - epsilon)^2))$ and $q >= ceil(log_2 (1 / delta))$ guarantees $bb(P)[union.big_(r=1)^q cal(G)_r] >= 1 - delta$.
  2. *Failure of Classical $s = 40 + 2 k$:* When $k >= 3$ clean clusters have equal size $m = n / k$ (with $b = 0$) and $n >= k s$, the probability that a single subsample covers all $k$ clusters satisfies
     $ bb(P)[forall j, space S_r inter C_j^* != emptyset] <= (e^(s / (k - 2))) / k. $
     Under Kaufman–Rousseeuw's classical heuristic $s = 40 + 2 k$ @kaufman1990finding, this coverage probability decays as $O(1 / k) -> 0$ as $k -> infinity$, proving that $s = Omega(k log k)$ is necessary for constant per-draw coverage.
]

#proof[
  (1) For a single cluster $C_j^*$ of size $n_j = p_j n$, sampling without replacement misses $C_j^*$ with probability
  $ bb(P)[S_r inter C_j^* = emptyset] &= product_(t=0)^(s-1) (n - n_j - t) / (n - t) \
    &<= (1 - p_j)^s <= e^(- s p_"min"), $
  since $(n - n_j - t)/(n - t) <= (n - n_j)/n = 1 - p_j$ for all $t >= 0$. A union bound over $j = 1, dots, k$ bounds the cluster-miss probability by $k e^(- s p_"min")$. Next, $|S_r inter X_"out"|$ follows a hypergeometric distribution with mean $s epsilon$; because Hoeffding's moment-generating bound holds without degradation for sampling without replacement, $bb(P)[|S_r inter X_"out"| >= s / 2] <= exp(- 2 s (1/2 - epsilon)^2)$. Summing the two failure events proves the first bound; for $s$ above the stated threshold, $bb(P)[cal(G)_r^c] <= 1/4 + 1/4 = 1/2$, so independence across $q$ draws yields $bb(P)[inter.big_(r=1)^q cal(G)_r^c] <= 2^(-q) <= delta$.

  (2) Let $E_j = {S_r inter C_j^* = emptyset}$ for $k >= 3$ equal clusters of size $m = n / k$, and let $N_"miss" = sum_(j=1)^k bb(1)_(E_j)$. For any pair $i != j$, exact hypergeometric expansion gives
  $ (bb(P)[E_i inter E_j]) / (bb(P)[E_i] bb(P)[E_j]) &= product_(t=0)^(s-1) ((n - 2 m - t)(n - t)) / ((n - m - t)^2) \
    &= product_(t=0)^(s-1) (1 - (m^2) / ((n - m - t)^2)) < 1. $
  Hence $"Cov"(bb(1)_(E_i), bb(1)_(E_j)) < 0$ strictly for every finite $n$, implying
  $ "Var"(N_"miss") <= sum_(j=1)^k bb(P)[E_j] (1 - bb(P)[E_j]) <= bb(E)[N_"miss"]. $
  Applying Chebyshev's inequality yields
  $ bb(P)[N_"miss" = 0] <= ("Var"(N_"miss")) / ((bb(E)[N_"miss"])^2) <= 1 / (bb(E)[N_"miss"]). $
  Since $n >= k s$ implies $t < s <= m$ for all $t in {0, dots, s - 1}$, we have $n - t > (k - 1) m$, so each factor satisfies
  $ (n - m - t) / (n - t) = 1 - m / (n - t) > 1 - 1 / (k - 1) = (k - 2) / (k - 1). $
  Because $ln((k - 1)/(k - 2)) = ln(1 + 1/(k - 2)) <= 1/(k - 2)$, taking the $s$-fold product yields $bb(P)[E_j] >= e^(- s / (k - 2))$ and $bb(E)[N_"miss"] >= k e^(- s / (k - 2))$. Substituting into the Chebyshev bound proves $bb(P)[N_"miss" = 0] <= e^(s / (k - 2)) / k$. For $s = 40 + 2 k$ and $k >= 42$, $s / (k - 2) <= 3.1$, forcing per-draw full coverage below $e^(3.1) / k -> 0$.
]

#subsec("C", "Minority-Cluster Blindness of Global Median Scoring")

Maximizing breakdown resistance to $50\%$ carries an exact geometric trade-off: a pure global median cannot penalize an uncovered minority cluster whose mass lies strictly below $50\%$.

#theorem(title: "Minority-Cluster Masking Counterexample", origin: "Original Counterexample")[
  For every integer $m >= 2$ and every target ratio $B > 0$, there exists a one-dimensional two-cluster dataset $X subset.eq bb(R)$ of size $n = 4 m - 1$ ($k = 2$) such that a corrupted medoid pair $M_"bad"$ placing *both* medoids inside the majority cluster strictly beats every cluster-covering medoid pair $M_"opt"$ under $J_"med"$, i.e., $J_"med" (M_"bad") < J_"med" (M_"opt")$, while incurring unbounded total-deviation distortion $"TD"(M_"bad") / "TD"(M_"opt") > B$.
]

#proof[
  Fix $m >= 2$ and $L > 2$. Construct $X subset.eq bb(R)$ of size $n = 4 m - 1$ consisting of: (i) a majority cluster $C_1$ of size $3 m$ with $m$ points each at $-1$, $0$, and $+1$; and (ii) a separated minority cluster $C_2$ of size $m - 1$ all colocated at $L$. Because $n = 4 m - 1$ is odd, $J_"med" (M)$ is uniquely the $2 m$-th smallest distance $d_((2 m)) (M)$.

  Evaluate the corrupted pair $M_"bad" = {-1, +1} subset.eq C_1$: the $2 m$ points at $-1$ and $+1$ have distance $0$; the $m$ points at $0$ have distance $1$; and the $m - 1$ points at $L$ have distance $L - 1$. Thus the $2 m$-th smallest distance is $J_"med" (M_"bad") = 0$, whereas $"TD"(M_"bad") = m + (m - 1)(L - 1)$.

  Now consider *any* covering medoid set $M_"opt" = {c, L}$ with $c in {-1, 0, +1}$ and one medoid on $C_2$. The points at $c$ ($m$ points) and $L$ ($m - 1$ points) have distance $0$—totaling only $2 m - 1$ zeros—while the remaining $2 m$ points of $C_1$ lie at distance $>= 1$ from $c$. Hence $d_((2 m)) (M_"opt") = 1$, giving $J_"med" (M_"opt") = 1 > 0 = J_"med" (M_"bad")$. Yet for $M_"opt" = {0, L}$, $"TD"(M_"opt") = 2 m$, so
  $ ("TD"(M_"bad")) / ("TD"(M_"opt")) &= (m + (m - 1)(L - 1)) / (2 m) \
    &-> infinity quad (L -> infinity). $
  Choosing $L > 1 + 2 m B / (m - 1)$ exceeds $B$ while CLARA-MED strictly selects $M_"bad"$.
]

#subsec("D", "Moment-Free Concentration and Trimmed Consistency")

Unlike K-Means (requiring $bb(E)[norm(X)_2^2] < infinity$ @pollard1981strong) and classical CLARA (requiring $bb(E)[d(X, M)] < infinity$), CLARA-MED's median score concentrates exponentially fast *without any finite moment assumption*.

#theorem(title: "Moment-Free Out-of-Subsample Generalization Bound", origin: "Original Proof")[
  Let $X_1, dots, X_n$ be i.i.d. draws from an arbitrary distribution $P$ on $(cal(M), d)$ (possibly heavy-tailed with infinite first moment). For any candidate set $M^((r))$ computed on subsample $S_r$ ($|S_r| = s$) and evaluated on the independent holdout $Y_r = X backslash S_r$ of size $m = n - s$, let $F_r (t) = bb(P)_(Y tilde.op P)[d(Y, M^((r))) <= t]$ have unique population median $t_r^*$ and local density $f_r (t) >= underline(f) > 0$ on $[t_r^* - eta, t_r^* + eta]$. Then with probability at least $1 - delta$,
  $ max_(1 <= r <= q) |hat(J)_"med" (M^((r)); Y_r) - t_r^*| \
    <= 1 / underline(f) sqrt((ln(2 q / delta)) / (2 (n - s))), $
  whenever the right-hand side does not exceed $eta$.
]

#proof[
  Condition on $S_r$ so $M^((r))$ is fixed and $Y_r = {Y_1, dots, Y_m}$ is an independent i.i.d. sample of size $m = n - s$. Let $hat(F)_m (t) = m^(-1) sum_(i=1)^m bb(1){d(Y_i, M^((r))) <= t}$. Fix $epsilon in (0, eta]$. If $hat(J)_"med" (M^((r)); Y_r) > t_r^* + epsilon$, then the upper median satisfies $d_((floor(m/2) + 1)) > t_r^* + epsilon$, implying at most $floor(m/2)$ points satisfy $d(Y_i, M^((r))) <= t_r^* + epsilon$, so $hat(F)_m (t_r^* + epsilon) <= 1 / 2$. Meanwhile,
  $ F_r (t_r^* + epsilon) &= 1 / 2 + integral_(t_r^* )^(t_r^* + epsilon) f_r (u) d u \
    &>= 1 / 2 + underline(f) epsilon. $
  By Hoeffding's inequality for bounded Bernoulli indicators in $[0, 1]$,
  $ &bb(P)[hat(J)_"med" > t_r^* + epsilon | S_r] \
    &<= bb(P)[hat(F)_m (t_r^* + epsilon) - F_r (t_r^* + epsilon) <= - underline(f) epsilon | S_r] \
    &<= exp(- 2 m (underline(f) epsilon)^2). $
  Symmetrically, $hat(J)_"med" < t_r^* - epsilon$ implies $d_((floor((m+1)/2))) < t_r^* - epsilon$; since $floor((m+1)/2) >= m / 2$ for all $m >= 1$, $hat(F)_m (t_r^* - epsilon) >= 1 / 2$ while $F_r (t_r^* - epsilon) <= 1 / 2 - underline(f) epsilon$, yielding the identical lower tail $exp(- 2 m (underline(f) epsilon)^2)$. Unconditioning over $S_r$ and applying a union bound over $r = 1, dots, q$ with $epsilon = underline(f)^(-1) sqrt(ln(2 q / delta) / (2 m))$ completes the proof.
]

#remark(title: "Strong Consistency of Trimmed Quantizers — Cited")[
  For continuous trimmed/quantile quantizers with trimming level $alpha in (0, 1)$, Cuesta-Albertos, Gordaliza, and Matrán @cuesta1997trimmed prove strong almost-sure consistency of empirical minimizers under mild continuity of $P$ (requiring no higher moments once $alpha > 0$), contrasting directly with Pollard's second-moment requirement for untrimmed K-Means @pollard1981strong.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Each of the $q$ iterations of CLARA-MED runs an inner medoid solver on $S_r$ of size $s$—costing $O(k s^2 T)$ with classical PAM (Alg. 4) or $O(s^2 T)$ with FastPAM1/FasterPAM (Alg. 24 @schubert2021fast)—followed by a full-dataset assignment scan costing $O(n k)$ distance evaluations and an exact linear-time $O(n)$ order-statistic selection (Quickselect / introselect) to compute $J_"med"$ and $J_(1/2)$. Total runtime is therefore $O(q (s^2 T + n k))$ (strictly *linear* in $n$ for fixed $s, k, q$!), and peak memory is $O(s^2 + n)$—completely eliminating the $O(n^2)$ pairwise matrix barrier of PAM and FastPAM.

Because all $q$ draws are mutually independent, CLARA-MED is embarrassingly parallel across $q$ threads or distributed workers, requiring only $O(q k)$ final medoid summaries to be merged.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for CLARA-MED],
    [Time / Space Complexity], [$O(q (s^2 T + n k))$ time (linear in $n$); $O(s^2 + n)$ memory (no $n times n$ matrix).],
    [Breakdown Robustness], [Exact $epsilon^* = ceil(n / 2)/n -> 50\%$ outer scoring breakdown (Theorem 1), vs. $1/n$ for CLARA.],
    [Moment Requirements], [None: moment-free $O(sqrt((log q)/(n - s)))$ concentration even on Cauchy tails (Theorem 4).],
    [Subsample Scaling], [Requires $s = Omega(k log k)$ and $q = Omega(log(1/delta))$ for full cluster coverage (Theorem 2).],
    [Geometric Trade-Off], [Global median can mask minority clusters of mass $< 50\%$ (Theorem 3); tune quantile $tau > 1 - p_"min"$.],
    [Recommended Companions], [FastPAM/FasterPAM (Alg. 24) as inner solver; Trimmed K-Means (Alg. 26) for Euclidean continuous trimming.],
  ),
  caption: [Theoretical guarantees and operational profile of CLARA-MED.],
)

#sec("IV", "Engineering and Scientific Applications")

CLARA-MED targets the intersection of *large-scale $n$* ($10^5$–$10^7$), *non-Euclidean dissimilarities*, and *heavy-tailed outlier corruption*, where both full-matrix medoid solvers (Alg. 24) and mean-based CLARA break down.

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Financial Fraud & High-Volume Account Profiling:* Retail banking and card-payment datasets span millions of accounts with Pareto-distributed whale transactions and fraudulent spikes. Mean-based CLARA and K-Means centroids are distorted by extreme amounts; CLARA-MED's $50\%$ breakdown preserves clean archetypal customer medoids in $O(n)$ time.
2. *Industrial IoT & Sensor Telemetry Zoning:* High-frequency sensor fleets suffer intermittent hardware dropouts, saturation spikes, and communication glitches generating near-infinite feature distances. Median scoring (Theorem 4) concentrates without requiring finite moments or fragile pre-cleaning heuristics.
3. *Large-Scale E-Commerce & Behavioral Cohorting:* Session-duration, basket-value, and clickstream edit-distance matrices on $n > 10^6$ users exhibit severe right skew. Using FastPAM inside $s = Theta(k log k)$ subsamples plus median full-dataset scoring yields real representative user profiles at linear scale.
4. *Cybersecurity & Network Flow Signature Extraction:* NetFlow logs mix dominant benign traffic classes with volumetric DDoS floods and malformed packet headers; CLARA-MED isolates typical protocol exemplars without allowing burst anomalies to hijack candidate selection.
5. *Large-Scale GPS & Mobility Trajectory Mining:* Urban fleet telemetry ($n > 10^5$ trajectories under DTW or Hausdorff distance) contains GPS multipath teleportation artifacts that blow up mean-distance scoring while leaving $J_"med"$ unaffected.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

CLARA-MED fails when: (a) *small minority clusters* ($p_"min" << 0.5$) must be discovered alongside wide majority clusters (Theorem 3; mitigate by replacing the $50\%$ median with a $tau$-quantile $Q_tau$ or $tau$-trimmed mean $J_tau$ where $1 - p_"min" < tau < 1 - epsilon$); (b) *$k$ is large while $s$ remains fixed at $40 + 2 k$* (Theorem 2(2)); (c) *$n$ is small* ($n < 1000$), where repeated subsampling is slower than running FasterPAM directly; and (d) *non-convex density manifolds* violate metric-ball assumptions.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Why CLARA-MED Fits], [Caution / Alternative],
    [Heavy-Tailed Finance ($n > 10^5$)], [$O(n)$ scaling + $50\%$ outlier immunity], [Use $tau$-trimmed score if rare niches exist],
    [IoT / Sensor Streams], [Moment-free concentration under spikes], [Set $s = Theta(k log k)$ (Thm. 2)],
    [Large Trajectory DTW], [Avoids $O(n^2)$ DTW matrix; resists GPS glitches], [Pair with FasterPAM (Alg. 24) inside $S_r$],
    [Small Clean Data ($n < 10^3$)], [Subsampling overhead unnecessary], [Use FasterPAM (Alg. 24) directly],
    [Rare-Anomaly Discovery], [Global median masks clusters $< 50\%$ (Thm. 3)], [Use DBSCAN (Alg. 31) or $Q_tau$ ($tau = 0.9$)],
  ),
  caption: [Domain suitability matrix and failure-mode mitigations for CLARA-MED.],
)

#sec("V", "Structural Relations and Empirical Behavior")

CLARA-MED resides on a continuous quantile–trimming spectrum unifying classical CLARA @srikantha2011clara @kaufman1990finding, Trimmed K-Means (Alg. 26 @cuesta1997trimmed), and FastPAM/FasterPAM (Alg. 24).

#proposition(title: "Unified Quantile–Trimming Sandwich Bound", origin: "Original Unified Proof")[
  For any $k$-medoid set $M subset.eq X$ and order index $m in {1, dots, n - 1}$ with fraction $tau_n = m / n in (0, 1)$, define the $tau_n$-quantile $Q_(tau_n) (M) := d_((m)) (M)$ and the $tau_n$-trimmed mean $J_(tau_n) (M) := m^(-1) sum_(r=1)^m d_((r)) (M)$. Then:
  1. $J_1 (M) = overline(D)(M)$ recovers classical CLARA, $Q_(ceil(n/2)/n) (M)$ recovers CLARA-MED, and minimizing $J_(1 - alpha)$ over squared Euclidean centroids recovers Trimmed K-Means (Alg. 26 @cuesta1997trimmed).
  2. For every $M subset.eq X$ and every $tau_n = m / n in (0, 1)$,
     $ J_(tau_n) (M) <= Q_(tau_n) (M) &<= (overline(D)(M) - tau_n J_(tau_n) (M)) / (1 - tau_n) \
       &<= (overline(D)(M)) / (1 - tau_n). $
]

#proof[
  Claim 1 follows immediately by substitution into Definition 1. For Claim 2, since $d_((1)) (M) <= dots <= d_((m)) (M)$, averaging over the first $m$ order statistics yields
  $ J_(tau_n) (M) = 1 / m sum_(r=1)^m d_((r)) (M) <= d_((m)) (M) = Q_(tau_n) (M). $
  For the upper bound, split the full sum $n overline(D)(M)$ at index $m$ and note that $d_((r)) (M) >= d_((m)) (M)$ for all $r > m$:
  $ n overline(D)(M) &= sum_(r=1)^m d_((r)) (M) + sum_(r=m+1)^n d_((r)) (M) \
    &>= m J_(tau_n) (M) + (n - m) Q_(tau_n) (M). $
  Dividing by $n$ gives $overline(D)(M) >= tau_n J_(tau_n) (M) + (1 - tau_n) Q_(tau_n) (M)$. Rearranging and using $J_(tau_n) (M) >= 0$ yields both upper inequalities identically.
]

Proposition 1 formalizes the exact asymmetry between mean and median scoring: every low-CLARA partition guarantees a bounded CLARA-MED score ($Q_(1/2) <= 2 overline(D)$), whereas under heavy tails $overline(D)$ can diverge while $Q_(1/2)$ and $J_(1/2)$ remain tightly bounded.

*Empirical Behavior on ClusterOverview Benchmarks:* On the ClusterOverview 4-dataset benchmark suite, CLARA-MED matches FastPAM / FasterPAM (Alg. 24) to three decimal places on *Customers* (Sil = *0.327*), *Multi-Shape* (Sil = *0.713*, ARI = N/A), and *High-D Subspace+Noise* (Sil = *0.302*, ARI = N/A), confirming that median-robust subsampled selection fully preserves high-cohesion medoid structure on skewed and noisy distributions while decisively outperforming K-Means|| (Alg. 23: Sil = 0.301, 0.609, 0.249). On *Iris*, CLARA-MED achieves slightly higher Silhouette (*0.446* vs. FastPAM's 0.437) and competitive ARI (*0.649* vs. FastPAM's 0.663 and K-Means||'s 0.431), reflecting the subtle trade-off between median cohesion and full-sum label alignment near species boundaries. Its mean runtime of *0.988 ms* ($3.7 times$ FastPAM's 0.268 ms and $9.4 times$ K-Means||'s 0.105 ms) reflects the fixed $q$-subsample multiplier on small datasets ($n <= 240$); as proven in Section III, CLARA-MED's $O(n)$ scaling overtakes $O(n^2)$ FastPAM once $n >> q s$.
