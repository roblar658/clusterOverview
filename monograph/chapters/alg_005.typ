#import "../template/ieee.typ": *

#alg_entry(
  id: 5,
  name: "CLARA (Clustering Large Applications)",
  category: "Cat. I: Centroid & Partitioning",
  year: "1986 / 1990",
  paradigm: "Subsampled k-Medoids Risk Minimization",
  complexity: $O(R (k s^2 T + n k))$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Let $(cal(X), d)$ be a finite metric space of cardinality $|cal(X)| = n$ with metric (or symmetric nonnegative dissimilarity) $d : cal(X) times cal(X) -> [0, infinity)$ satisfying $d(x, x) = 0$. Given $k in {1, dots, n}$, the *metric $k$-medoids problem* seeks a subset of $k$ representative observations $M = {m_1, dots, m_k} subset.eq cal(X)$ minimizing the total assignment dissimilarity @kaufman1990finding:
$ L(M) = sum_(x in cal(X)) min_(m in M) d(x, m), quad overline(L)(M) = 1 / n L(M). $
Each medoid set $M$ induces a discrete Voronoi partition $cal(C)(M) = {C_1, dots, C_k}$ where $C_j = {x in cal(X) : c(x) = j}$ and $c(x) = arg min_(j in {1, dots, k}) d(x, m_j)$ under deterministic tie-breaking.

Partitioning Around Medoids (PAM, Alg. 4) solves this combinatorial objective via greedy BUILD initialization followed by steepest-descent SWAP refinement over all $k(n - k)$ candidate swaps, requiring $Theta(n^2)$ memory for the dissimilarity matrix and $O(k(n - k)^2)$ time per pass @kaufman1990finding. On datasets with $n >= 10^4$, materializing $n^2$ pairwise distances becomes intractable. *CLARA* (Clustering Large Applications), introduced by Kaufman and Rousseeuw @srikantha2011clara @kaufman1990finding, overcomes this quadratic barrier by decoupling medoid *discovery* from medoid *evaluation*: it applies PAM independently to $R$ uniform random subsamples of fixed size $s << n$ (classically $s = 40 + 2k$, $R = 5$) and scores each resulting candidate set $M_r subset.eq S_r$ against the *entire* dataset $cal(X)$ in $O(n k)$ distance evaluations without ever forming the full $n times n$ matrix.

#definition(title: "Subsampled Empirical Medoid Risk")[
  For any subsample $S subset.eq cal(X)$ of size $s >= k$, define the *empirical subsample risk*
  $ hat(L)_S (M) = 1 / s sum_(x in S) min_(m in M) d(x, m), quad M subset.eq S, space |M| = k. $
  Let $cal(A)_"PAM" (S)$ denote the $k$-medoid subset returned by running PAM (BUILD + SWAP) on $(S, d|_(S times S))$. Given $R$ i.i.d. uniform without-replacement draws $S_1, dots, S_R tilde "Unif"(binom(cal(X), s))$, the *CLARA estimator* is
  $ M_"CLARA" = arg min_(r in {1, dots, R}) overline(L)(M_r), quad M_r = cal(A)_"PAM" (S_r). $
]

Because candidate scoring uses the exact full-population functional $overline(L)$ rather than the subsample proxy $hat(L)_(S_r)$, CLARA protects against selecting an anomalously dense subsample whose internal cost misrepresents the global dataset (Algorithm 1).

#algbox(
  title: "CLARA (Clustering Large Applications)",
  inputs: [Dataset $(cal(X), d)$ of size $n$, integer $k >= 1$, subsample size $s in [k, n]$, number of draws $R >= 1$.],
  outputs: [Medoid set $M_"best" subset.eq cal(X)$ and partition $cal(C)(M_"best")$.],
)[
  1. Initialize $L_"best" <- +infinity$ and $M_"best" <- emptyset$.
  2. *For* $r = 1, dots, R$ *do*:
  3. #h(1em) Draw $S_r subset.eq cal(X)$ uniformly without replacement ($|S_r| = s$).
  4. #h(1em) Form the local $s times s$ dissimilarity matrix $D^((r))_(u v) = d(u, v)$ for $u, v in S_r$.
  5. #h(1em) Run PAM (Alg. 4: BUILD + SWAP) on $(S_r, D^((r)))$ to obtain candidate medoids $M_r subset.eq S_r$ ($|M_r| = k$).
  6. #h(1em) Evaluate exact full-dataset cost in $O(n k)$ distance queries:
     $ L(M_r) = sum_(i=1)^n min_(m in M_r) d(x_i, m). $
  7. #h(1em) *If* $L(M_r) < L_"best"$, set $L_"best" <- L(M_r)$ and $M_"best" <- M_r$.
  8. Assign each $x_i in cal(X)$ to $c_i = arg min_(m in M_"best") d(x_i, m)$.
  9. *Return* $(M_"best", cal(C)(M_"best"))$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five foundational results for CLARA: (1) tight without-replacement cluster-coverage bounds quantifying both sufficiency and necessity of the sample budget $s R$; (2) a finite-sample metric excess-risk inequality decomposing estimation error from subsample medoid restriction; (3) an explicit adversarial counterexample proving unbounded approximation ratio under rare clusters; (4) exact invariance, within-cluster breakdown point $1/2$, and CLARA's probabilistic subsampling shield against isolated outliers; and (5) the classical $O(k)$ acceleration via FastPAM swap decomposition.

#subsec("A", "Without-Replacement Cluster Coverage and Sample Complexity")

#theorem(title: "Exact Cluster-Hit Sufficiency and Necessity Bounds", origin: "Original Proof")[
  Let $cal(C)^* = {C_1^*, dots, C_k^*}$ be an optimal partition of $cal(X)$ with cluster shares $p_j = |C_j^*| / n$ and minimum share $p_"min" = min_j p_j > 0$. For $R$ independent without-replacement subsamples $S_1, dots, S_R$ of size $s$, let $E_r = {S_r inter C_j^* != emptyset space forall j}$ be the event that draw $r$ hits every cluster. Then:
  1. For each draw, $bb(P)(E_r^c) <= sum_(j=1)^k (1 - p_j)^s <= k e^(-s p_"min")$.
  2. Across $R$ draws, $bb(P)(union.big_(r=1)^R E_r) >= 1 - [min(1, k e^(-s p_"min"))]^R$.
  3. Conversely, if $s <= n / 2$ and $p_"min" <= 1 / 4$, then achieving failure probability $bb(P)(inter.big_(r=1)^R E_r^c) <= delta$ necessarily requires
     $ R s >= 1 / (4 p_"min") ln(1 / delta). $
]

#proof[
  (1) Fix cluster $C_j^*$ of size $n_j = p_j n$. Under uniform sampling without replacement of $s$ objects from $n$, the exact miss probability is
  $ bb(P)(S_r inter C_j^* = emptyset) &= binom(n - n_j, s) / binom(n, s) \
    &= product_(t=0)^(s-1) (n - n_j - t) / (n - t) = product_(t=0)^(s-1) (1 - n_j / (n - t)). $
  Because $n - t <= n$ for every $t >= 0$, each factor satisfies $1 - n_j / (n - t) <= 1 - n_j / n = 1 - p_j <= e^(-p_j)$. Multiplying over $t = 0, dots, s - 1$ gives $bb(P)(S_r inter C_j^* = emptyset) <= (1 - p_j)^s <= e^(-s p_j)$. Applying a union bound over $j = 1, dots, k$ proves Claim 1.

  (2) By independence of the $R$ draws, $bb(P)(inter.big_(r=1)^R E_r^c) = [bb(P)(E_1^c)]^R <= [min(1, k e^(-s p_"min"))]^R$, yielding Claim 2 immediately.

  (3) Let $j^*$ attain $p_"min"$, with $n_0 = p_"min" n$. Because $s <= n / 2$, we have $n - t >= n - s + 1 > n / 2$ for every $0 <= t <= s - 1$, hence $u_t := n_0 / (n - t) <= 2 p_"min" <= 1 / 2$. Consider $f(u) = e^(-2 u) - (1 - u)$ on $[0, 1/2]$: $f(0) = 0$, $f(1/2) = e^(-1) - 1/2 < 0$, and $f''(u) = 4 e^(-2 u) > 0$, so convexity forces $1 - u >= e^(-2 u)$ for all $u in [0, 1/2]$. Consequently,
  $ bb(P)(E_r^c) &>= bb(P)(S_r inter C_(j^*)^* = emptyset) \
    &= product_(t=0)^(s-1) (1 - u_t) >= e^(-4 s p_"min"). $
  Raising to the $R$-th power yields $delta >= bb(P)(inter.big_(r=1)^R E_r^c) >= e^(-4 R s p_"min")$; taking logarithms gives $4 R s p_"min" >= ln(1 / delta)$, completing the proof.
]

#subsec("B", "Uniform Concentration and Excess-Risk Decomposition")

#lemma(title: "Uniform Without-Replacement Medoid Concentration", origin: "Original Proof")[
  Suppose $(cal(X), d)$ has finite diameter $Delta = max_(x, y in cal(X)) d(x, y) < infinity$. For a uniform without-replacement subsample $S$ of size $s >= k$ and any $delta in (0, 1)$, with probability at least $1 - delta$:
  $ sup_(M subset.eq cal(X), space |M|=k) |hat(L)_S (M) - overline(L)(M)| \
    <= Delta sqrt((k ln n + ln(2 / delta)) / (2 s)). $
]

#proof[
  Fix any $M subset.eq cal(X)$ with $|M| = k$ and set $Z_i = min_(m in M) d(X_i, m) in [0, Delta]$ for $X_1, dots, X_s$ drawn without replacement from $cal(X)$. By Hoeffding's convex domination theorem for sampling without replacement, for every $lambda in bb(R)$, $bb(E)[exp(lambda sum_(i=1)^s (Z_i - overline(L)(M)))] <= bb(E)[exp(lambda sum_(i=1)^s (Y_i - overline(L)(M)))]$, where $Y_1, dots, Y_s$ are i.i.d. uniform with-replacement draws from $cal(X)$. Applying Hoeffding's lemma to each bounded variable $Y_i in [0, Delta]$ and Chernoff's bound gives
  $ bb(P)(|hat(L)_S (M) - overline(L)(M)| > epsilon) <= 2 exp(- (2 s epsilon^2) / Delta^2). $
  There are exactly $binom(n, k) <= n^k$ candidate $k$-medoid sets $M subset.eq cal(X)$. By the union bound,
  $ bb(P)(sup_(|M|=k) |hat(L)_S (M) - overline(L)(M)| > epsilon) \
    <= 2 n^k exp(- (2 s epsilon^2) / Delta^2). $
  Setting the right-hand side equal to $delta$ and solving for $epsilon$ yields the stated uniform deviation bound.
]

#theorem(title: "Metric Excess-Risk Bound for CLARA", origin: "Original Derivation")[
  Let $M^* = {m_1^*, dots, m_k^*}$ be a globally optimal $k$-medoid set on $(cal(X), d)$ with cluster shares $p_j = |C_j^*| / n$, and suppose $cal(A)_"PAM" (S_r)$ returns an $alpha$-approximate minimizer of $hat(L)_(S_r)$ over $S_r$ (e.g., $alpha = 5$ for metric 1-swap local optima @kanungo2004local, or $alpha = 1$ for an exact subsample solver). Then for any $delta in (0, 1)$, with probability at least $1 - delta$,
  $ overline(L)(M_"CLARA") &<= alpha overline(L)(M^*) + alpha min_(1 <= r <= R) sum_(j=1)^k p_j d(m_j^*, S_r) \
    &quad + (alpha + 1) Delta sqrt((k ln n + ln(2 R / delta)) / (2 s)), $
  where $d(m_j^*, S_r) = min_(u in S_r) d(m_j^*, u)$ is the nearest-subsample projection distance to $m_j^*$.
]

#proof[
  Apply Lemma 1 with confidence parameter $delta / R$ to each subsample $S_r$ ($r = 1, dots, R$) and take a union bound; with probability at least $1 - delta$, the uniform error satisfies $|hat(L)_(S_r)(M) - overline(L)(M)| <= epsilon_(s, R)$ simultaneously for all $r$ and all $|M| = k$, where $epsilon_(s, R) = Delta sqrt((k ln n + ln(2 R / delta)) / (2 s))$.

  For each draw $r$, project each optimal medoid onto $S_r$ by choosing $tilde(m)_(r, j) in arg min_(u in S_r) d(m_j^*, u)$, and form a $k$-element subset $tilde(M)_r subset.eq S_r$ containing ${tilde(m)_(r, 1), dots, tilde(m)_(r, k)}$ (padding with arbitrary distinct points of $S_r$ if collisions occur, which can only decrease point-to-set distances since $s >= k$). By the triangle inequality, every $x in C_j^*$ satisfies
  $ min_(m in tilde(M)_r) d(x, m) <= d(x, tilde(m)_(r, j)) <= d(x, m_j^*) + d(m_j^*, S_r). $
  Averaging over all $x in cal(X)$ gives the *projection inequality*:
  $ overline(L)(tilde(M)_r) <= overline(L)(M^*) + sum_(j=1)^k p_j d(m_j^*, S_r). $
  Since $M_r subset.eq S_r$ is an $alpha$-approximate minimizer of $hat(L)_(S_r)$ and $tilde(M)_r subset.eq S_r$ is feasible on $S_r$, on the high-probability event we obtain
  $ overline(L)(M_r) &<= hat(L)_(S_r)(M_r) + epsilon_(s, R) \
    &<= alpha hat(L)_(S_r)(tilde(M)_r) + epsilon_(s, R) \
    &<= alpha (overline(L)(tilde(M)_r) + epsilon_(s, R)) + epsilon_(s, R) \
    &<= alpha overline(L)(M^*) + alpha sum_(j=1)^k p_j d(m_j^*, S_r) + (alpha + 1) epsilon_(s, R). $
  Because CLARA computes $overline(L)(M_r)$ *exactly* on $cal(X)$ and outputs $M_"CLARA" = arg min_r overline(L)(M_r)$, taking the minimum over $r in {1, dots, R}$ establishes the theorem.
]

#subsec("C", "Unbounded Suboptimality Under Minority Clusters")

Theorem 2 reveals that CLARA's excess risk is governed by $min_r sum_j p_j d(m_j^*, S_r)$, which explodes if no subsample lands near a distant low-weight cluster.

#theorem(title: "Unbounded Approximation Ratio of Fixed-Sample CLARA", origin: "Original Counterexample")[
  Fix any sample size $s >= k >= 2$, any draw count $R >= 1$, any target factor $M > 0$, and any failure tolerance $eta in (0, 1)$. There exists a finite one-dimensional Euclidean dataset $cal(X) subset.eq bb(R)$ of size $n$ such that with probability at least $1 - eta$, CLARA returns a medoid set satisfying $L(M_"CLARA") / L(M^*) > M$.
]

#proof[
  Choose $n > 2 R s / eta$ with $n >= s + 2$. Construct $cal(X) = C_"bulk" union C_"rare" subset.eq bb(R)$ where $C_"bulk"$ contains $n - 2$ distinct points evenly spaced in $(0, 1]$, and $C_"rare" = {D, D + 1}$ for a separation $D > 1$ to be chosen.

  For the global optimum $M^*$, placing $k - 1$ medoids inside $C_"bulk"$ (whose diameter is at most $1$) and one medoid at $D in C_"rare"$ incurs cost at most $1$ per point in $C_"bulk"$ and $1$ for $D + 1$, giving
  $ L(M^*) <= (n - 2) dot.c 1 + 1 = n - 1 < n. $
  Now consider any single CLARA draw $S_r$ of size $s$. By symmetry and union bound over the two points of $C_"rare"$,
  $ bb(P)(S_r inter C_"rare" != emptyset) <= 2 s / n. $
  Applying a union bound across all $R$ draws gives
  $ bb(P)(exists r in {1, dots, R} : S_r inter C_"rare" != emptyset) <= (2 R s) / n < eta. $
  Hence with probability at least $1 - eta$, every candidate medoid in $union.big_(r=1)^R S_r$ lies inside $(0, 1]$. In that event, both points of $C_"rare"$ lie at distance at least $D - 1$ from $M_"CLARA"$, forcing $L(M_"CLARA") >= 2(D - 1)$. Consequently,
  $ L(M_"CLARA") / L(M^*) > (2(D - 1)) / n. $
  Selecting $D = n M / 2 + 2$ yields $L(M_"CLARA") / L(M^*) > M$ with probability at least $1 - eta$, independent of how accurately PAM solves each subsample.
]

#subsec("D", "Invariances, Breakdown Point, and Subsampling Outlier Shield")

#proposition(title: "Invariance, Local Breakdown 1/2, and Outlier Shield", origin: "Original Proof")[
  For any finite metric space $(cal(X), d)$:
  1. *Positive Homothety & Isometry Invariance:* For any bijection $phi$ satisfying $d'(phi(x), phi(y)) = c dot.c d(x, y)$ ($c > 0$), CLARA's PAM swaps and full-dataset selection rule are strictly invariant (conditioned on fixed sample indices).
  2. *Within-Cluster Medoid Breakdown Point:* For a fixed cluster $C_j$ of size $m$, the Donoho–Huber breakdown point of its 1-medoid is $epsilon_"loc"^* = ceil(m / 2) / m -> 1 / 2$ (contrasting with $1 / m$ for K-Means centroids).
  3. *Joint Outlier Hijacking and Subsampling Shield:* For $k >= 2$, full PAM's global medoid-set breakdown point is $1 / n$, yet CLARA excludes $b$ isolated extreme outliers from all $R$ candidate pools with probability at least $(1 - b / (n - s + 1))^(R s)$ (exact $(1 - s / n)^R$ for $b = 1$), while their relative impact on comparing any two uncorrupted medoid sets in a ball of radius $rho$ is uniformly bounded by $2 b rho$ regardless of outlier magnitude.
]

#proof[
  (1) Multiplying $d$ by $c > 0$ scales both $hat(L)_(S_r)$ and $overline(L)$ by $c$, leaving every comparison and $arg min$ unchanged.

  (2) Let $m - b$ uncorrupted points of $C_j$ lie in a ball $B(x_0, rho)$ while $b <= floor((m - 1)/2)$ corrupted points move to distance $R_c -> infinity$ from $x_0$. Comparing any uncorrupted candidate $u in B(x_0, rho)$ against a corrupted point $v$ at distance $R_c$:
  $ sum_(x in C_j) d(x, u) &<= (m - b)(2 rho) + b(R_c + rho), \
    sum_(x in C_j) d(x, v) &>= (m - b)(R_c - rho). $
  Subtracting gives $L_j (v) - L_j (u) >= (m - 2 b) R_c - (3 m - b) rho -> +infinity$ because $m - 2 b >= 1$. Conversely, placing $b = ceil(m / 2)$ corrupted points together at $v$ with $R_c -> infinity$ flips the inequality ($m - 2 b <= 0$), establishing $epsilon_"loc"^* = ceil(m / 2) / m$.

  (3) For full PAM with $k >= 2$, moving a single point $x_"out"$ to distance $R_c > L_(k-1)^* (cal(X) backslash {x_"out"})$ strictly lowers total cost by placing a singleton medoid at $x_"out"$ (cost $< R_c$ vs. $>= R_c - 2 rho$), so full PAM's global set breakdown point is $1 / n$. In CLARA, the probability that $b$ outliers miss all $R$ draws is $[binom(n-b, s) / binom(n, s)]^R >= (1 - b / (n - s + 1))^(R s)$ (and $(1 - s/n)^R$ for $b = 1$). On this event, every candidate set $M_1, M_2$ lies inside $B(x_0, rho)$, and by the triangle inequality $|d(x_"out", M_1) - d(x_"out", M_2)| <= 2 rho$ for each outlier, bounding total objective perturbation by $2 b rho$ independently of $R_c$.
]

#subsec("E", "FastCLARA Acceleration via Cached Swap Decomposition")

#theorem(title: "Exact O(k) Acceleration of Subsample SWAP in FastCLARA", origin: "Classical — Cited")[
  Replacing classic PAM's SWAP loop inside CLARA with the FastPAM cached nearest/second-nearest decomposition @schubert2019fasterpam @schubert2021fast reduces each subsample SWAP iteration from $O(k(s - k)^2)$ to $O((s - k)^2)$ while returning the *exact same* best swap at every step, reducing total CLARA runtime from $O(R(k(s - k)^2 T + n k))$ to $O(R((s - k)^2 T + n k))$.
]

#proof(label: "Proof Synthesis (Cited)")[
  On subsample $S_r$, cache for each $x_o in S_r$ its nearest medoid index $c_1 (o)$ with distance $d_1 (o)$ and second-nearest distance $d_2 (o)$. Schubert and Rousseeuw @schubert2019fasterpam @schubert2021fast decompose the cost change of swapping medoid $m_i$ for candidate $x_c in S_r backslash M$ into a candidate-independent removal loss $Delta L^(-m_i) = sum_(o : c_1 (o) = i) (d_2 (o) - d_1 (o))$ plus an additive gain evaluated simultaneously across all $i in {1, dots, k}$ during a single $O(s)$ loop over $o in S_r$: each $x_o$ adds $min(0, d(x_o, x_c) - d_1 (o))$ to a shared accumulator $Delta_"add"$ (or adjusts entry $i = c_1 (o)$ by $min(d(x_o, x_c), d_2 (o)) - d_2 (o)$). Eliminating the inner loop over $i$ saves an exact factor of $k$, permitting $sqrt(k)$-larger subsamples $s$ at constant time cost.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Classical CLARA executes $R$ subsample PAM runs costing $O(R k(s - k)^2 T)$ plus $R$ full-dataset evaluations costing $O(R n k d)$ in $bb(R)^d$ (or $R n k$ oracle calls in general metric spaces). Crucially, space complexity is only $O(s^2 + n d)$—strictly linear in $n$—and for classical sizes $s in [40 + 2k, 200]$, the $s times s$ distance matrix occupies $< 320$ KB, remaining entirely inside L2 CPU cache. Because the $R$ draws are mutually independent, CLARA is embarrassingly parallel across $R$ workers, requiring only a single scalar reduction of $L(M_r)$. Geometrically, CLARA inherits PAM's metric Voronoi cells $V_j (M) = {x : d(x, m_j) <= d(x, m_ell) space forall ell}$ while restricting vertices to observed points; in high dimensions, metric concentration @beyer1999nearest shrinks relative separation between $d_1$ and $d_2$, increasing variance across subsamples and warranting larger $s$ (e.g., $s = 80 + 4k$ @schubert2021fast).

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for CLARA],
    [Time / Space Complexity], [$O(R(k s^2 T + n k))$ classic; $O(R(s^2 T + n k))$ FastCLARA; $O(s^2 + n d)$ space (no $O(n^2)$ matrix).],
    [Geometric Inductive Bias], [Metric Voronoi cells centered on actual dataset exemplars; handles arbitrary dissimilarities.],
    [Determinism / Seeding], [Randomized over $R$ subsamples; deterministic PAM inside each draw; exact global scoring.],
    [Invariances], [Invariant under metric isometries and positive scaling $c dot.c d$; depends on chosen metric.],
    [Outlier Robustness], [Within-cluster breakdown $1/2$; subsampling shields against isolated outlier hijacking (Prop. 1).],
    [Recommended Companions], [PAM (Alg. 4) for $n <= 3000$; CLARANS (Alg. 6) for dynamic full-graph search; FastPAM (Alg. 24); CLARA-MED (Alg. 25); BICO (Alg. 78).],
  ),
  caption: [Summary of theoretical guarantees and operational profile of CLARA.],
)

#sec("IV", "Engineering and Scientific Applications")

CLARA excels whenever datasets exceed PAM's quadratic memory ceiling ($n >> 10^4$), arithmetic centroids are undefined or unrepresentative, and actual exemplar objects are required.

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Large-Scale Spatial Facility Location & Logistics:* Minimizing total $L_1$ or road-network travel distance $sum_i d(x_i, m_(c_i))$ over millions of demand addresses directly models warehouse and emergency-depot siting. Because medoids are restricted to observed sites, every prototype is an immediately feasible physical location.
2. *Enterprise Customer & Telecom Archetyping:* Customer tables combining skewed spend, tenure, and categorical flags via Gower or Manhattan dissimilarities violate Euclidean K-Means assumptions. CLARA scales linearly in $n$ while outputting real inspectable customer personas resilient to billing spikes (Proposition 1).
3. *Genomic & Protein Sequence Representative Selection:* When clustering $10^5$ biological sequences under edit or alignment dissimilarity, computing the full $n^2$ pairwise alignment matrix is prohibitive. CLARA computes only $R binom(s, 2) + R n k$ alignments—orders of magnitude fewer—while yielding authentic consensus sequences.
4. *Industrial Sensor & Telemetry Summarization:* High-frequency IoT streams suffer transient spike artifacts that corrupt K-Means centroids ($epsilon^* = 1/n$ @lloyd1982least). CLARA's Lipschitz-1 metric objective and subsampling shield prevent isolated sensor faults from distorting operational regime prototypes.
5. *Metric-Space Database Pivot Selection:* Metric index structures (VP-trees, M-trees, and BUBBLE Alg. 71 @ganti1999bubble) require well-dispersed routing pivots across non-Euclidean spaces; CLARA selects $k$ balanced medoid pivots in a single linear pass over disk-resident tables.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

CLARA degrades sharply under (a) *severely imbalanced or rare clusters* where $p_"min" << 1 / (R s)$ (Theorem 1(3) and Theorem 3); (b) *large cluster counts $k$* where $s = 40 + 2k$ provides barely two sampled points per cluster, inflating projection error $d(m_j^*, S_r)$ in Theorem 2; (c) *non-convex density manifolds* (better served by DBSCAN Alg. 31 or HDBSCAN Alg. 32); and (d) *heavily overlapping clusters* where restricting medoids to static subsamples misses optimal boundary pivots (improved by CLARANS Alg. 6).

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Why CLARA Fits / Prerequisite], [Caution / Alternative],
    [Facility Location], [Real-site medoids; $L_1$/network metric], [Increase $s$ for small regions],
    [Customer Personas], [Gower/Manhattan + outlier shield], [Use CLARANS (Alg. 6) if overlap],
    [Bio-Sequences], [Avoids $O(n^2)$ pairwise alignments], [Use FastCLARA (Alg. 24) for big $k$],
    [IoT Telemetry], [Bounded spike influence; $O(n k)$ scan], [Use Trimmed (Alg. 26) for bursts],
    [Rare Fraud / Micro], [Fails if $p_"min" < (R s)^(-1)$ (Thm. 3)], [Use Coreset / BICO (Alg. 78)],
  ),
  caption: [Application suitability and engineering matrix for CLARA.],
)

#sec("V", "Structural Relations and Empirical Behavior")

CLARA sits directly between exact quadratic medoid optimization (PAM, Alg. 4), sensitivity-weighted coreset compression @harpeled2004coresets, and dynamic randomized graph search (CLARANS, Alg. 6 @ng2002clarans).

#proposition(title: "Exact PAM Limit and Coreset Contrast", origin: "Original Derivation")[
  As $s -> n$, CLARA with $R = 1$ reduces identically to PAM (Alg. 4), driving both the concentration term $epsilon_(s, R)$ and projection bias $d(m_j^*, S_r)$ in Theorem 2 to zero. Conversely, for fixed $s << n$, CLARA performs *uniform unweighted* subsampling whose worst-case error is unbounded (Theorem 3), whereas sensitivity-weighted coreset sampling @harpeled2004coresets @fichtenberger2013bico assigns weight $w(x) = 1 / q(x)$ inversely proportional to sampling probability $q(x)$, preserving $(1 plus.minus epsilon)$ multiplicative cost uniformly across all $M$ even when $p_"min" -> 0$.
]

#proof[
  When $s = n$, $S_1 = cal(X)$ deterministically, so $hat(L)_(S_1) equiv overline(L)$, $d(m_j^*, S_1) = 0$, and Step 5 of Algorithm 1 executes full PAM on $cal(X)$. Under sensitivity sampling, each point $x in cal(X)$ is drawn with probability $q(x) >= s_0 (x) / S_"tot"$ proportional to its worst-case sensitivity $s_0 (x) = sup_M d(x, M) / L(M)$, ensuring that isolated remote clusters of tiny share $p_"min"$ have large sensitivity and are sampled with high probability, eliminating the uniform-sampling blind spot proved in Theorem 3.
]

CLARA also inspired CLARA-MED (Alg. 25) for enhanced median-based scoring and directly motivated CLARANS (Alg. 6 @ng2002clarans), which removes CLARA's static restriction $M_r subset.eq S_r$ by subsampling SWAP transitions across the full dataset instead of subsampling observations upfront.

In empirical evaluation on the ClusterOverview 4-dataset benchmark suite (validated via Silhouette @rousseeuw1987silhouettes and ARI), CLARA achieves strong convex separation on Iris (`Sil = 0.439`, `ARI = 0.605`), solid clustering on Multi-Shape (`Sil = 0.712`, `ARI = N/A`), moderate structure on Customers (`Sil = 0.327`), and expected metric-concentration attenuation on High-D Subspace+Noise (`Sil = 0.299`, `ARI = N/A`), with a mean runtime of `0.845 ms`. Notably, while CLARA matches CLARANS's Silhouette on Iris (`0.439`) and Multi-Shape (`0.712`), its Iris ARI (`0.605`) trails CLARANS (`0.676`)—directly reflecting Theorem 2's static subsample projection bias $d(m_j^*, S_r)$ on overlapping species boundaries—and its runtime (`0.845 ms`) reflects executing full PAM convergence over multiple draws before global scoring.
