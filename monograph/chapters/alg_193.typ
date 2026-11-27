#import "../template/ieee.typ": *

#alg_entry(
  id: 193,
  name: "CluStream",
  category: "Cat. VIII: Neural, Deep, Stream & Evolutionary",
  year: "2003",
  paradigm: "Two-Phase Online Micro-Cluster & Pyramidal Offline Stream",
  complexity: $O(n q d + I_"off" q k d)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

CluStream @aggarwal2003clustream addresses clustering over an unbounded, evolving data stream $cal(S) = {(x_1, T_1), (x_2, T_2), dots}$ with $x_i in RR^d$ arriving at nondecreasing timestamps $T_1 <= T_2 <= dots$. Because stream concepts drift over time, an analyst rarely wants a single clustering over the entire infinite past; instead, queries request a $k$-partition over an arbitrary user-specified *time horizon* $h > 0$, i.e., the substream $X_h (t_c) = {x_i : t_c - h < T_i <= t_c}$ ending at current clock $t_c$.

CluStream decouples stream summarization from macro-level analysis into an *online micro-clustering component* and an *offline macro-clustering component*. Extending the Clustering Feature (CF) vector of BIRCH @zhang1996birch with temporal sufficient statistics, CluStream summarizes any set $C = {(x_(i_1), T_(i_1)), dots, (x_(i_n), T_(i_n))}$ by a $(2d + 3)$-dimensional *micro-cluster tuple*
$
  "CF"(C) = (overline("CF2")^x, overline("CF1")^x, "CF2"^t, "CF1"^t, n) ,
$
where $overline("CF2")^x = sum_(r=1)^n x_(i_r)^(dot.o 2) in RR^d$ (entrywise squares), $overline("CF1")^x = sum_(r=1)^n x_(i_r) in RR^d$, $"CF2"^t = sum_(r=1)^n T_(i_r)^2$, $"CF1"^t = sum_(r=1)^n T_(i_r)$, and $n = abs(C)$. From $"CF"(C)$ alone, CluStream evaluates the spatial centroid $c = overline("CF1")^x / n$, the spatial root-mean-square (RMS) radius
$
  R(C) = sqrt(1 / n sum_(j=1)^d "CF2"_j^x - norm(c)_2^2) ,
$
and the temporal mean and standard deviation $mu_t = "CF1"^t / n$, $sigma_t = sqrt("CF2"^t / n - mu_t^2)$.

Online, CluStream maintains a fixed memory budget of $q >> k$ micro-clusters ${M_1, dots, M_q}$. Each incoming point $x$ either joins its nearest micro-cluster $M_p$ if $norm(x - c_p)_2 <= t_"fac" R(M_p)$ (with boundary factor $t_"fac" >= 1$; for a singleton $n_p = 1$, $R(M_p)$ is set to the distance to its nearest neighboring micro-cluster), or seeds a new singleton micro-cluster. To preserve the count $q$, CluStream checks whether the minimum *relevance stamp* $r_m (M)$—defined as the $m / (2n)$-th normal quantile $mu_t + Phi^(-1) (m / (2n)) sigma_t$ for $n >= 2m$ (or $mu_t$ for $n < 2m$) approximating the average arrival time of the last $m$ points—falls below recency threshold $delta$. If so, the expired micro-cluster is deleted; otherwise, the two closest micro-clusters are merged and their identifier lists concatenated. Micro-cluster snapshots are stored on disk according to a *Pyramidal Time Frame (PTF)* parameterized by integer base $alpha >= 2$ and capacity parameter $l >= 1$: snapshots of order $i >= 0$ are recorded whenever $t equiv 0 mod alpha^i$, retaining only the $alpha^l + 1$ most recent snapshots per order.

#definition(title: "Horizon-Restricted Offline Macro-Clustering Objective")[
  Given micro-cluster snapshots at $t_c$ and $t_s approx t_c - h$, let $cal(M)_h = {M_1, dots, M_q}$ denote the horizon-differenced micro-clusters with centroids $c_r$ and positive counts $n_r$. Offline macro-clustering minimizes the weighted $k$-Means objective @lloyd1982least
  $
    Phi_"off" (bold(mu)) = sum_(r=1)^q n_r min_(j in {1, dots, k}) norm(c_r - mu_j)_2^2
  $
  over macro-centers $bold(mu) = (mu_1, dots, mu_k) in (RR^d)^k$.
]

#algbox(
  title: "CluStream (Online Micro-Clustering + Offline Macro-Clustering)",
  inputs: [Stream $cal(S)$, micro-cluster budget $q$, boundary factor $t_"fac"$, recency threshold $delta$, PTF parameters $(alpha, l)$, macro-clusters $k$, horizon $h$.],
  outputs: [Macro-cluster centers $bold(mu) = (mu_1, dots, mu_k)$ over horizon $h$.],
)[
  1. *Init:* Collect initial buffer $X_"init"$, run $k$-Means with $q$ centers to form ${M_1, dots, M_q}$ with unique IDs.
  2. *Online Phase:* For each arriving $(x, t)$:
  3. #h(1em) Find nearest micro-cluster $p = arg min_(r in [q]) norm(x - c_r)_2$ and compute its maximum boundary $B_p = t_"fac" R(M_p)$.
  4. #h(1em) *If* $norm(x - c_p)_2 <= B_p$: absorb $(x, t)$ into $M_p$ via CF vector addition.
  5. #h(1em) *Else:* find $r^* = arg min_r r_m (M_r)$. If $r_m (M_(r^*)) < delta$, delete $M_(r^*)$; else merge closest pair $(u, v) = arg min_(r < s) norm(c_r - c_s)_2$ by CF addition and ID-list union. Insert singleton ${(x, t)}$.
  6. #h(1em) At clock $t$, store PTF snapshot at all orders $i$ with $alpha^i | t$, evicting any order-$i$ snapshot beyond the $alpha^l + 1$ most recent.
  7. *Offline Phase (query $(t_c, h)$):* Retrieve snapshots $cal(S)(t_c)$ and stored $cal(S)(t_s)$ nearest $t_c - h$; subtract matched ID-list CF tuples to obtain $cal(M)_h$.
  8. Run weighted $k$-Means ($k$ centers, initialized via weighted seeding) on ${(c_r, n_r) : M_r in cal(M)_h}$ until convergence; *return* $bold(mu)$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

#subsec("A", "Exact CF Algebra and Parallel-Axis Merge Identity")

Both online merging and offline horizon differencing rely on the exact linear closure of CF tuples @zhang1996birch @aggarwal2003clustream.

#lemma(title: "CF Additivity, Subtractivity, and Exact Merge Variance Formula", origin: "Classical — Cited")[
  Let $C_1, C_2$ be disjoint timestamped point multisets with centroids $c_1, c_2$ and sizes $n_1, n_2$. Then:
  1. $"CF"(C_1 union C_2) = "CF"(C_1) + "CF"(C_2)$, and for $C_2 subset.eq C_1$, $"CF"(C_1 without C_2) = "CF"(C_1) - "CF"(C_2)$.
  2. The within-cluster sum of squared errors is recovered exactly from $"CF"(C)$ via $"SSE"(C) = sum_(j=1)^d "CF2"_j^x - (1 / n) norm(overline("CF1")^x)_2^2$.
  3. Merging $C_1$ and $C_2$ increases total squared error by the exact Ward-type increment
  $
    "SSE"(C_1 union C_2) = "SSE"(C_1) + "SSE"(C_2) + (n_1 n_2) / (n_1 + n_2) norm(c_1 - c_2)_2^2 .
  $
]

#proof[
  Part 1 follows directly from linearity of finite summation over disjoint index sets: each entry of $"CF"(C)$ is an additive functional $sum_((x, T) in C) phi(x, T)$ with $phi(x, T) = (x^(dot.o 2), x, T^2, T, 1)$. For part 2, expanding the squared norm gives
  $
    "SSE"(C) &= sum_(x in C) (norm(x)_2^2 - 2 chevron.l x, c chevron.r + norm(c)_2^2) \
             &= sum_(j=1)^d sum_(x in C) x_j^2 - 2 chevron.l n c, c chevron.r + n norm(c)_2^2 \
             &= sum_(j=1)^d "CF2"_j^x - (1 / n) norm(overline("CF1")^x)_2^2 .
  $
  For part 3, let $c = (n_1 c_1 + n_2 c_2) / (n_1 + n_2)$ be the merged centroid. Then $c_1 - c = (n_2 / (n_1 + n_2))(c_1 - c_2)$ and $c_2 - c = (-n_1 / (n_1 + n_2))(c_1 - c_2)$. Writing $x - c = (x - c_r) + (c_r - c)$ for $x in C_r$ ($r in {1, 2}$) and noting $sum_(x in C_r) (x - c_r) = bold(0)$ yields
  $
    "SSE"(C_1 union C_2) &= sum_(r=1)^2 ("SSE"(C_r) + n_r norm(c_r - c)_2^2) \
                         &= "SSE"(C_1) + "SSE"(C_2) + (n_1 n_2^2 + n_2 n_1^2) / (n_1 + n_2)^2 norm(c_1 - c_2)_2^2 ,
  $
  which simplifies identically to the stated formula.
]

#subsec("B", "Pyramidal Time Frame Space and Horizon Approximation")

The Pyramidal Time Frame guarantees logarithmic storage while answering any historical horizon query within a sharp relative error @aggarwal2003clustream.

#theorem(title: "PTF Storage Complexity and Exact Horizon Relative Error Bound", origin: "Classical — Cited")[
  For integer base $alpha >= 2$, capacity parameter $l >= 1$, and stream clock $T >= 1$:
  1. The total number of snapshots stored in the Pyramidal Time Frame is at most $(alpha^l + 1)(floor(log_alpha T) + 1)$.
  2. For every requested horizon $h in [1, T - 1]$, there exists a stored snapshot at time $t_s$ satisfying the exact two-sided relative horizon error bound
  $
    abs((T - t_s) - h) / h <= 1 / (2 alpha^(l - 1)) ,
  $
  and a stored snapshot $t_s^- <= T - h$ satisfying $((T - t_s^-) - h) / h < 1 / alpha^(l - 1)$.
]

#proof[
  Only orders $i in {0, 1, dots, floor(log_alpha T)}$ satisfy $alpha^i <= T$ and can hold snapshots by time $T$. Since each order retains at most $alpha^l + 1$ snapshots, total stored snapshots never exceed $(alpha^l + 1)(floor(log_alpha T) + 1)$, proving part 1.

  For part 2, let $m_i = floor(T / alpha^i)$. For any order $i$ with $m_i >= alpha^l$, the retained order-$i$ snapshots are precisely the consecutive multiples $t_(i, r) = (m_i - r) alpha^i$ for $r = 0, 1, dots, alpha^l$, spanning from $t_"old" (i) = (m_i - alpha^l) alpha^i$ up to $m_i alpha^i$. Because $T - m_i alpha^i >= 0$, the lookback reach satisfies $T - t_"old" (i) >= alpha^(i + l)$.

  If $h <= alpha^l$, order $i = 0$ retains every integer in $[T - alpha^l, T]$, so $t_s = T - h$ is stored with zero error. Otherwise, assume $h > alpha^l$ and choose the unique integer $i >= 1$ such that
  $
    alpha^(i - 1 + l) < h <= alpha^(i + l) .
  $
  Then $t_"old" (i) <= T - alpha^(i + l) <= T - h$. Define the floor multiple $t_s^- = floor((T - h) / alpha^i) alpha^i$ and ceiling multiple $t_s^+ = ceil((T - h) / alpha^i) alpha^i$. Clearly $t_"old" (i) <= t_s^- <= T - h <= t_s^+$. Moreover, since $l >= 1$, we have $h > alpha^(i - 1 + l) >= alpha^i$, whence $T - h < T - alpha^i < m_i alpha^i$; as $m_i alpha^i$ is a multiple of $alpha^i$, this implies $t_s^+ <= m_i alpha^i$. Therefore *both* $t_s^-$ and $t_s^+$ are currently retained in order $i$. Since $0 <= (T - h) - t_s^- < alpha^i$, selecting $t_s^-$ gives $((T - t_s^-) - h) / h < alpha^i / alpha^(i - 1 + l) = 1 / alpha^(l - 1)$. Selecting $t_s in {t_s^-, t_s^+}$ closer to $T - h$ ensures $abs((T - h) - t_s) <= alpha^i / 2$, dividing by $h > alpha^(i - 1 + l)$ yields $1 / (2 alpha^(l - 1))$.
]

#subsec("C", "Exact Offline Distortion Decomposition and Error Transfer")

We now establish the exact algebraic coupling between CluStream's online micro-cluster quantization energy and offline macro-clustering optimality.

#theorem(title: "Exact Pythagorean Decomposition and Macro-Clustering Approximation Bound", origin: "Original Derivation")[
  Let $cal(M)_h = {M_1, dots, M_q}$ partition horizon $X_h$ with counts $n_r$, centroids $c_r$, and micro-cluster quantization energy $C_"micro" = sum_(r=1)^q n_r R(M_r)^2$. Let $"OPT"_k (X_h)$ be the unconstrained point-level $k$-Means optimum on $X_h$.
  1. For any macro-cluster partition $Pi = {P_1, dots, P_k}$ of $cal(M)_h$ with centers $bold(mu) = (mu_1, dots, mu_k)$, the exact point-level cost satisfies
  $
    sum_(j=1)^k sum_(M_r in P_j) sum_(x in M_r) norm(x - mu_j)_2^2 = C_"micro" + sum_(j=1)^k sum_(M_r in P_j) n_r norm(c_r - mu_j)_2^2 .
  $
  2. If an offline solver produces $hat(bold(mu))$ satisfying $Phi_"off" (hat(bold(mu))) <= gamma min_(bold(mu)) Phi_"off" (bold(mu))$ for approximation factor $gamma >= 1$, then the true point-level $k$-Means cost $"Cost"(X_h, hat(bold(mu)))$ obeys
  $
    "Cost"(X_h, hat(bold(mu))) <= 2 gamma "OPT"_k (X_h) + (2 gamma + 1) C_"micro" .
  $
]

#proof[
  For each micro-cluster $M_r$ and any $mu in RR^d$, expanding $x - mu = (x - c_r) + (c_r - mu)$ gives
  $
    sum_(x in M_r) norm(x - mu)_2^2 &= sum_(x in M_r) norm(x - c_r)_2^2 + 2 chevron.l sum_(x in M_r) (x - c_r), c_r - mu chevron.r \
                                    &+ n_r norm(c_r - mu)_2^2 = n_r R(M_r)^2 + n_r norm(c_r - mu)_2^2 ,
  $
  because $sum_(x in M_r) (x - c_r) = bold(0)$. Summing over $M_r in P_j$ and $j = 1, dots, k$ establishes part 1.

  For part 2, let $bold(nu)^*$ be an optimal unconstrained $k$-center set for $X_h$, and denote by $nu(y) = arg min_(nu in bold(nu)^*) norm(y - nu)_2$ the nearest center in $bold(nu)^*$ to any $y in RR^d$. For every $M_r$ and every $x in M_r$, optimality of $nu(c_r)$ and the inequality $norm(a + b)_2^2 <= 2 norm(a)_2^2 + 2 norm(b)_2^2$ give
  $
    norm(c_r - nu(c_r))_2^2 <= norm(c_r - nu(x))_2^2 <= 2 norm(c_r - x)_2^2 + 2 norm(x - nu(x))_2^2 .
  $
  Summing over $x in M_r$ yields $n_r norm(c_r - nu(c_r))_2^2 <= 2 n_r R(M_r)^2 + 2 sum_(x in M_r) norm(x - nu(x))_2^2$. Summing over $r = 1, dots, q$ bounds the optimal offline objective:
  $
    min_(bold(mu)) Phi_"off" (bold(mu)) <= Phi_"off" (bold(nu)^*) <= 2 C_"micro" + 2 "OPT"_k (X_h) .
  $
  Finally, assigning each point $x in M_r$ to its nearest center in $hat(bold(mu))$ costs no more than assigning all of $M_r$ to $hat(mu)(c_r)$, so by part 1,
  $
    "Cost"(X_h, hat(bold(mu))) &<= C_"micro" + Phi_"off" (hat(bold(mu))) \
                               &<= C_"micro" + gamma (2 C_"micro" + 2 "OPT"_k (X_h)) ,
  $
  yielding $2 gamma "OPT"_k (X_h) + (2 gamma + 1) C_"micro"$.
]

#subsec("D", "Worst-Case Unbounded Distortion of Greedy Closest-Pair Merging")

Theorem 3 shows that CluStream's accuracy hinges on keeping $C_"micro"$ small. Unlike importance-sampled coresets @ackermann2012streamkmpp, CluStream's unweighted closest-pair merge rule can suffer unbounded approximation error under sparse outliers.

#proposition(title: "Unbounded Approximation Ratio Counterexample for CluStream", origin: "Original Counterexample")[
  Fix any micro-cluster budget $q >= 3$, boundary factor $t_"fac" >= 1$, and $k = 2$ in $RR^1$. For every $M > 0$, there exists an active-horizon stream (with no recency expirations) on which CluStream's online micro-cluster merge forces *every* offline $2$-Means solution on $cal(M)$ to satisfy
  $
    "Cost"(X, hat(bold(mu))) / "OPT"_2 (X) > M .
  $
]

#proof[
  Pick a large integer $N_0$ and set $delta = min(1 / (4 t_"fac"), 1 / sqrt(N_0))$. Construct a stream $X subset RR^1$ arriving inside a single active recency window (so no micro-cluster expires):
  1. First arrive $N_0$ points alternating at $-delta$ and $+delta$ (forming micro-cluster $M_1$ with $n_1 = N_0, c_1 = 0, R_1 = delta$, boundary $t_"fac" delta <= 1/4$).
  2. Next arrive $N_0$ points alternating at $1 - delta$ and $1 + delta$. Because $1 - delta > 1/4$, the first point seeds $M_2$; its singleton boundary is $t_"fac" (1 - delta) > 2 delta$, so the second point $1 + delta$ is absorbed, giving $R_2 = delta$ and absorbing all remaining $N_0 - 2$ points ($n_2 = N_0, c_2 = 1$).
  3. Finally arrive $q - 1$ outlier points at $z_r = 1 + (t_"fac" + 2)^r$ for $r = 1, dots, q - 1$.

  By induction on $r$, $z_1 - 1 = t_"fac" + 2 > 1/4$ creates singleton $M_3$ whose nearest neighbor is $M_2$ at distance $D_1 = t_"fac" + 2$. For $r >= 2$, $z_r - z_(r-1) = (t_"fac" + 1)(t_"fac" + 2)^(r-1) > t_"fac" D_(r-1)$, so $z_r$ lies outside the singleton boundary of $M_(r+1)$ and spawns a fresh micro-cluster. At $r = q - 1$, $q$ micro-clusters already exist ($M_1, M_2$ plus $q - 2$ singletons), forcing a closest-pair merge. Because $abs(c_1 - c_2) = 1$ whereas every other centroid pair is separated by at least $t_"fac" + 2 >= 3$, CluStream uniquely merges $M_1$ and $M_2$ into $M_(12)$ of weight $2 N_0$, centroid $1/2$, and internal energy
  $
    "SSE"(M_(12)) = 2 N_0 (1 / 4 + delta^2) > N_0 / 2 .
  $
  By Theorem 3, any offline $2$-Means macro-clustering over $cal(M)$ incurs cost at least $"SSE"(M_(12)) > N_0 / 2$. Conversely, placing unconstrained centers at $0$ and $1$ yields
  $
    "OPT"_2 (X) <= 2 N_0 delta^2 + sum_(r=1)^(q-1) (z_r - 1)^2 <= 2 + C(q, t_"fac") ,
  $
  where $C(q, t_"fac")$ is independent of $N_0$. Choosing $N_0 > 2 M (2 + C(q, t_"fac"))$ gives $"Cost"(X, hat(bold(mu))) / "OPT"_2 (X) > M$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Online processing of each incoming stream point requires computing Euclidean distances to $q$ micro-cluster centroids in $O(q d)$ time, plus $O(1)$ CF addition if absorbed. Whenever a new micro-cluster forces a merge, scanning pairwise micro-cluster distances takes $O(q^2 d)$ (or $O(q d)$ using an incremental priority queue of nearest pairs). Over $n$ points, online time is $O(n q d)$ in the absorption-dominated regime, while working RAM holds only $O(q d)$ floating-point scalars for active CF tuples. Historical disk storage under the Pyramidal Time Frame requires $O(q d alpha^l log_alpha T)$ space (Theorem 2). Offline macro-clustering executes $I_"off"$ weighted Lloyd iterations over $q$ weighted centroids in $O(I_"off" q k d)$ time—completely independent of the horizon point count $N_h$.

Geometrically, CluStream imposes *spherical micro-quantization* online and *convex Voronoi partitioning* offline. Because all CF tuples are stored in a contiguous array of size $q$ (typically $q in [10k, 100k]$), inner-loop centroid scans enjoy high L1/L2 cache locality and trivial SIMD vectorization.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Online Time / Point], [$O(q d)$ typical absorption; $O(q^2 d)$ worst-case merge event.],
    [Offline Query Time], [$O(I_"off" q k d)$, independent of horizon size $N_h$.],
    [RAM / Disk Space], [$O(q d)$ RAM; $O(q d alpha^l log_alpha T)$ PTF snapshot disk space.],
    [Horizon Support], [Arbitrary historical windows $[t_c - h, t_c]$ within $1 / (2 alpha^(l-1))$ error.],
    [Inductive Bias], [Spherical micro-clusters; convex macro-clusters via weighted $k$-Means.],
    [Concept Drift], [Handled via relevance-stamp expiry and horizon subtraction.],
    [Outlier Sensitivity], [Moderate--High: unweighted closest-pair merges vulnerable (Prop. 4).],
  ),
  caption: [Operational and computational profile of CluStream.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Network Intrusion Detection and Packet Flow Monitoring:* High-speed router telemetry generates millions of connection records with evolving attack signatures @aggarwal2003clustream. CluStream's $O(q d)$ online step tracks live traffic at line rate while PTF horizon queries compare last-hour flow clusters against last-week baselines via exact CF subtraction (Lemma 1).
2. *Financial Tick-Stream Regime Analysis:* Intraday order-book feature vectors undergo intraday volatility shifts and structural regime breaks. Arbitrary-horizon offline queries let quantitative analysts inspect multi-scale regimes (5-minute vs. full-session horizons) without re-ingesting historical ticks.
3. *Clickstream and E-Commerce Customer Segmentation:* User behavioral vectors evolve across promotional campaigns and seasonal cycles. Maintaining $q$ compact micro-clusters captures fine-grained micro-segments online while offline weighted $k$-Means updates marketing personas on demand.
4. *IoT Sensor Fleet Telemetry and Predictive Maintenance:* Industrial vibration and thermal streams exhibit gradual drift prior to mechanical failure. Micro-cluster relevance stamps naturally retire obsolete operating states while preserving compressed summaries of emerging fault modes.
5. *Cloud Datacenter Workload Characterization:* Virtual-machine CPU, memory, and I/O utilization streams fluctuate on diurnal cycles. PTF snapshots allow capacity planners to isolate peak-hour clusters from off-peak baselines using logarithmic storage (Theorem 2).

#subsec("B", "Failure Modes and Domain Suitability Matrix")

CluStream fails predictably in three regimes. First, *non-convex manifold streams*: because the offline stage executes weighted $k$-Means, crescent or ring-shaped clusters are bisected along linear hyperplanes even when micro-clusters cover the manifold accurately. Second, *heavy-tailed outlier bursts*: as proved in Proposition 4, rare extreme outliers consume micro-cluster slots and force heavy core clusters to merge prematurely, inflating $C_"micro"$. Third, *high-dimensional streams*: Euclidean radius thresholds $t_"fac" R$ lose selectivity as $d$ grows due to pairwise distance concentration @beyer1999nearest.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Evolving Multi-Horizon Streams], [Strong: PTF + CF subtraction enable retroactive windows.], [Monitor outlier rate; compare ClusTree @kranen2011clustree.],
    [Noisy / Non-Convex Streams], [Weak: offline $k$-Means splits curves; outliers trigger merges.], [Use DenStream @cao2006denstream or D-Stream @chen2007dstream.],
    [Stationary Large-Scale Streams], [Moderate: fast $O(q d)$, lacks worst-case coreset bound.], [Prefer StreamKM++ @ackermann2012streamkmpp or BICO @fichtenberger2013bico.],
    [High-Dimensional Embeddings], [Moderate--Weak: spherical boundary degrades with $d$.], [Apply PCA / random projection or Contrastive Clustering @li2021contrastive.],
    [Real-Time Low-D Telemetry], [Strong: fixed $O(q d)$ RAM, predictable sub-ms latency.], [Tune $q >= 10 k$, $t_"fac" approx 2$, and recency $delta$.],
  ),
  caption: [Domain suitability and failure-mode matrix for CluStream.],
)

#sec("V", "Structural Relations and Empirical Behavior")

CluStream serves as the architectural bridge between static CF-tree compression in BIRCH @zhang1996birch and modern evolving-stream algorithms such as DenStream @cao2006denstream, D-Stream @chen2007dstream, ClusTree @kranen2011clustree, and BICO @fichtenberger2013bico.

#proposition(title: "Exact Static-Limit Reduction to Leaf-Level BIRCH and Weighted Merge Remedy", origin: "Original Derivation")[
  Under a full-stream horizon ($h = T$) with no recency deletion ($delta = -oo$), CluStream's online spatial sufficient statistics $(overline("CF2")^x, overline("CF1")^x, n)$ coincide identically with a flat-array specialization of BIRCH's leaf CF tuples @zhang1996birch. Moreover, replacing CluStream's unweighted closest-pair merge rule with the minimum-increment rule
  $
    (u^*, v^*) = arg min_(r < s) (n_r n_s) / (n_r + n_s) norm(c_r - c_s)_2^2
  $
  strictly minimizes the one-step increase in $C_"micro"$ and eliminates the unweighted-merge pathology of Proposition 4 on the constructed instance.
]

#proof[
  Dropping temporal fields $("CF2"^t, "CF1"^t)$ when $h = T$ and $delta = -oo$ leaves $(n, overline("CF1")^x, overline("CF2")^x)$, which is algebraically identical to BIRCH's CF triple $(N, overline("LS"), "SS")$ (summing $overline("CF2")^x$ coordinates yields $"SS"$). By Lemma 1, merging any pair $(M_r, M_s)$ increases total micro-cluster quantization energy $C_"micro" = sum_u "SSE"(M_u)$ by precisely
  $
    Delta C_"micro" (r, s) = (n_r n_s) / (n_r + n_s) norm(c_r - c_s)_2^2 .
  $
  Thus selecting $(u^*, v^*)$ minimizes the exact one-step growth of $C_"micro"$. On the counterexample of Proposition 4, merging $(M_1, M_2)$ incurs $Delta C_"micro" = (N_0 / 2) dot.op 1^2 = N_0 / 2$, whereas merging any two singleton outliers $(M_r, M_s)$ incurs $Delta C_"micro" = (1 / 2)(z_(s-2) - z_(r-2))^2 = O(1)$ independent of $N_0$. For all sufficiently large $N_0$, the Ward-weighted merge strictly merges outliers instead of $(M_1, M_2)$, preserving $"SSE"(M_1) + "SSE"(M_2) = 2 N_0 delta^2 <= 2$.
]

*Empirical Benchmark Analysis:* On the ClusterOverview four-dataset benchmark, CluStream attains Iris Sil = 0.445, ARI = 0.608; Customers Sil = 0.344; Multi-Shape Sil = 0.609 (ARI = N/A); and High-D Sil = 0.253 (ARI = N/A), with a fast mean runtime of 2.246 ms. These metrics confirm our analytical profile: CluStream outperforms Contrastive Clustering on Customers (0.344 vs. 0.299) and runs 4.1$times$ faster (2.246 ms vs. 9.230 ms), while closely tracking spherical $k$-Means-family scores on Iris (ARI = 0.608). However, its Multi-Shape Silhouette (0.609) lags far behind density-based DenStream (0.793) because offline weighted $k$-Means imposes convex macro-boundaries, and its High-D Silhouette (0.253) trails StreamKM++ (0.302) because unweighted spherical micro-cluster boundaries suffer higher quantization energy $C_"micro"$ in high dimensions (Theorem 3).
