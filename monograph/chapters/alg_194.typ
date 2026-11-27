#import "../template/ieee.typ": *

#alg_entry(
  id: 194,
  name: "DenStream",
  category: "Cat. VIII: Neural, Deep, Stream & Evolutionary",
  year: "2006",
  paradigm: "Damped-Window Density-Based Two-Phase Stream Clustering",
  complexity: $O(n (abs(P) + abs(O)) d + abs(P)^2 d)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

DenStream @cao2006denstream extends density-based spatial clustering (DBSCAN @ester1996dbscan) to evolving data streams corrupted by arbitrary noise, eliminating both CluStream's fixed spherical macro-clusters @aggarwal2003clustream and its vulnerability to outlier bursts. DenStream operates under the *damped window model*: each stream point $p_i in RR^d$ arriving at timestamp $T_i$ carries an exponentially decaying weight $f(t - T_i) = 2^(-lambda (t - T_i))$ at current time $t >= T_i$, where $lambda > 0$ governs the forgetting rate.

Rather than retaining raw points, DenStream summarizes local dense neighborhoods via *weighted micro-clusters* bounded in radius by $epsilon > 0$. For a multiset of timestamped points $C = {(p_1, T_1), dots, (p_m, T_m)}$, its damped Clustering Feature tuple at time $t$ is $("CF1", "CF2", w)$, where
$
  w(t) &= sum_(j=1)^m f(t - T_j) , \
  overline("CF1")(t) &= sum_(j=1)^m f(t - T_j) p_j , quad overline("CF2")(t) = sum_(j=1)^m f(t - T_j) p_j^(dot.o 2) .
$
The weighted centroid is $c = overline("CF1") / w$ and the CF-computable weighted RMS radius is
$
  r_"rms" (C) = sqrt(1 / w sum_(l=1)^d "CF2"_l - norm(c)_2^2) .
$
DenStream separates genuine evolving clusters from transient noise by maintaining two distinct micro-cluster buffers controlled by core weight threshold $mu > 0$ and outlier factor $beta in (0, 1]$ (with $beta mu > 1$):
1. *Potential core-micro-clusters* ($p$-micro-clusters, set $P$): tuples $("CF1", "CF2", w)$ satisfying $w >= beta mu$ and $r_"rms" <= epsilon$.
2. *Outlier micro-clusters* ($o$-micro-clusters, set $O$): tuples $("CF1", "CF2", w, t_0)$ satisfying $w < beta mu$ and $r_"rms" <= epsilon$, annotated with creation timestamp $t_0$.

Every $T_p = ceil((1 / lambda) log_2 ((beta mu) / (beta mu - 1)))$ time steps, DenStream executes a *pruning sweep*: any $p$-micro-cluster whose decayed weight falls below $beta mu$ is removed, and any $o$-micro-cluster whose weight falls below the dynamic lower-envelope threshold
$
  xi(t_c, t_0) = (2^(-lambda (t_c - t_0 + T_p)) - 1) / (2^(-lambda T_p) - 1)
$
is permanently discarded as noise. On demand, the offline phase applies a weighted variant of DBSCAN @ester1996dbscan to the centroids of $P$: a $p$-micro-cluster $c_p in P$ is a *core* micro-cluster if $w(c_p) >= mu$, and $c_q in P$ is *directly density-reachable* from $c_p$ if $w(c_p) >= mu$ and $norm(c_p - c_q)_2 <= 2 epsilon$.

#definition(title: "Damped Core-Micro-Cluster and Offline Density Connectivity")[
  At time $t$, a weighted micro-cluster $c_p = ("CF1", "CF2", w)$ is a *core-micro-cluster* ($c$-micro-cluster) if $w >= mu$ and $r_"rms" (c_p) <= epsilon$. Two $p$-micro-clusters $c_a, c_b in P$ are *density-connected* if there exists a core $c_o in P$ from which both $c_a$ and $c_b$ are reachable via chains of directly density-reachable $p$-micro-clusters (step length $<= 2 epsilon$). Maximal density-connected subsets of $P$ form the offline clusters.
]

#algbox(
  title: "DenStream (Online Merging/Pruning + Offline DBSCAN)",
  inputs: [Stream $cal(S)$, radius $epsilon$, core weight $mu$, outlier factor $beta in (0, 1]$ ($beta mu > 1$), decay rate $lambda > 0$.],
  outputs: [Arbitrary-shape density-connected clusters at query time $t_c$.],
)[
  1. *Init:* Collect initial points, run DBSCAN to seed $p$-micro-cluster buffer $P$; set $O <- emptyset$ and $T_p <- ceil((1 / lambda) log_2 ((beta mu) / (beta mu - 1)))$.
  2. *Online Merging:* For each arriving point $p$ at timestamp $t$:
  3. #h(1em) Find nearest $c_p in P$ (decayed to $t$); compute trial radius $r_p'$ after adding $p$. If $r_p' <= epsilon$, absorb $p$ into $c_p$ and continue.
  4. #h(1em) Else find nearest $c_o in O$ (decayed to $t$); compute trial radius $r_o'$. If $r_o' <= epsilon$, absorb $p$ into $c_o$; if updated $w(c_o) >= beta mu$, move $c_o$ from $O$ to $P$.
  5. #h(1em) Else create new $o$-micro-cluster $c_"new" = (p, p^(dot.o 2), 1, t)$ and insert into $O$.
  6. *Periodic Pruning:* If $t equiv 0 mod T_p$:
  7. #h(1em) Decay all buffers to $t$; delete every $c_p in P$ with $w(c_p) < beta mu$, and every $c_o in O$ with $w(c_o) < xi(t, t_0 (c_o))$.
  8. *Offline Query (at $t_c$):* Decay $P$ to $t_c$; run $2 epsilon$-DBSCAN with core threshold $w >= mu$ on centroids ${c_p : c_p in P}$; *return* connected components.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

#subsec("A", "Decay Invariance and Exact Rank-1 Absorption Dynamics")

Lazy timestamped updating is valid in DenStream because exponential damping preserves both centroids and radii identically across idle intervals.

#lemma(title: "Decay Invariance, Exact Incremental Radius Update, and Mean-Radius Bound", origin: "Original Derivation")[
  Let $C$ be a weighted micro-cluster with state $("CF1"(t_1), "CF2"(t_1), w(t_1))$, centroid $c$, and RMS radius $r_"rms"$. Over an idle interval $Delta t = t_2 - t_1 >= 0$:
  1. Scaling $("CF1", "CF2", w)$ by $2^(-lambda Delta t)$ yields exact damped statistics at $t_2$, while $c(t_2) = c(t_1)$ and $r_"rms" (t_2) = r_"rms" (t_1)$ remain strictly invariant.
  2. Absorbing a new point $p$ at $t_2$ produces the exact updated squared RMS radius
  $
    r_"rms"'^2 = tilde(w) / (tilde(w) + 1) r_"rms"^2 + tilde(w) / (tilde(w) + 1)^2 norm(p - c)_2^2 ,
  $
  where $tilde(w) = 2^(-lambda Delta t) w(t_1)$.
  3. The RMS radius upper-bounds the weighted mean deviation: $r_"mean" (C) = (1 / w) sum_j w_j norm(p_j - c)_2 <= r_"rms" (C)$.
]

#proof[
  For any point $(p_j, T_j) in C$, $f(t_2 - T_j) = 2^(-lambda (t_2 - t_1)) f(t_1 - T_j) = 2^(-lambda Delta t) f(t_1 - T_j)$. By linearity, all three CF components multiply by $a = 2^(-lambda Delta t) > 0$. Thus $c(t_2) = (a overline("CF1")) / (a w) = c(t_1)$ and
  $
    r_"rms"^2 (t_2) = (a sum_l "CF2"_l) / (a w) - norm(c(t_2))_2^2 = r_"rms"^2 (t_1) ,
  $
  proving part 1.

  For part 2, after decaying to $tilde(w) = a w(t_1)$ and adding $p$ with weight $1$, the updated weight is $tilde(w) + 1$ and the new centroid is $c' = (tilde(w) c + p) / (tilde(w) + 1)$. Hence $c - c' = (c - p) / (tilde(w) + 1)$ and $p - c' = tilde(w)(p - c) / (tilde(w) + 1)$. Writing $p_j - c' = (p_j - c) + (c - c')$ and using $sum_j tilde(w)_j (p_j - c) = bold(0)$ yields
  $
    (tilde(w) + 1) r_"rms"'^2 &= tilde(w) r_"rms"^2 + tilde(w) norm(c - c')_2^2 + norm(p - c')_2^2 \
                              &= tilde(w) r_"rms"^2 + (tilde(w) + tilde(w)^2) / (tilde(w) + 1)^2 norm(p - c)_2^2 \
                              &= tilde(w) r_"rms"^2 + tilde(w) / (tilde(w) + 1) norm(p - c)_2^2 .
  $
  Dividing by $tilde(w) + 1$ yields the exact update formula. For part 3, applying Cauchy--Schwarz to $sum_j sqrt(w_j) (sqrt(w_j) norm(p_j - c)_2)$ gives
  $
    (sum_j w_j norm(p_j - c)_2)^2 <= (sum_j w_j) (sum_j w_j norm(p_j - c)_2^2) = w^2 r_"rms"^2 (C) ,
  $
  and dividing by $w^2$ proves $r_"mean" (C) <= r_"rms" (C)$.
]

#subsec("B", "Deterministic Stream-Length-Independent Memory Bound")

Unlike CluStream, DenStream requires no fixed cluster budget $q$; exponential damping guarantees an exact finite memory ceiling on $P$ regardless of stream duration @cao2006denstream.

#theorem(title: "Global Damped Weight Cap and Potential Buffer Bound", origin: "Classical — Cited")[
  Suppose the stream arrives at constant rate $v >= 1$ points per unit time with decay rate $lambda > 0$. Then for every time $t >= 0$:
  1. The total decayed weight of all stream points is strictly bounded by $W(t) < v / (1 - 2^(-lambda))$.
  2. The number of active $p$-micro-clusters satisfies the deterministic, stream-length-independent upper bound
  $
    abs(P) <= floor(v / (beta mu (1 - 2^(-lambda)))) .
  $
]

#proof[
  At integer clock $t$, summing the geometric decay over all timestamps $tau in {0, 1, dots, t}$ gives
  $
    W(t) = sum_(tau=0)^t v dot.op 2^(-lambda (t - tau)) = v (1 - 2^(-lambda (t + 1))) / (1 - 2^(-lambda)) < v / (1 - 2^(-lambda)) .
  $
  Because every point belongs to at most one active $p$-micro-cluster and each $c_p in P$ maintains weight $w(c_p) >= beta mu$, summing over $P$ yields
  $
    beta mu abs(P) <= sum_(c_p in P) w(c_p) <= W(t) < v / (1 - 2^(-lambda)) .
  $
  Dividing by $beta mu$ and taking the integer floor completes the proof.
]

#subsec("C", "Exact Invariant Trajectory of the Outlier Pruning Rule")

We prove that the pruning threshold $xi(t_c, t_0)$ is not heuristic: it is the exact invariant solution of the periodic-arrival recurrence separating viable clusters from doomed noise @cao2006denstream.

#theorem(title: "Telescoping Invariance and Soundness of Outlier Pruning", origin: "Original Proof")[
  Let $beta mu > 1$, $lambda > 0$, and $T_p = ceil((1 / lambda) log_2 ((beta mu) / (beta mu - 1)))$. Set $a = 2^(-lambda T_p) in (0, 1)$. Then:
  1. For all $t_c >= t_0$ and every integer $m >= 0$, the threshold $xi$ satisfies the exact affine recurrence identity
  $
    a^m xi(t_c, t_0) + sum_(j=0)^(m-1) a^j = xi(t_c + m T_p, t_0) < 1 / (1 - a) <= beta mu .
  $
  2. Consequently, any $o$-micro-cluster with $w(t_c) < xi(t_c, t_0)$ that subsequently receives at most one point every $T_p$ time steps can *never* reach the promotion threshold $beta mu$ at any future time.
]

#proof[
  Write $b = 2^(-lambda (t_c - t_0)) in (0, 1]$. By definition of $xi$, multiplying numerator and denominator by $-1$ gives $xi(t_c, t_0) = (1 - a b) / (1 - a)$. Using the finite geometric sum $sum_(j=0)^(m-1) a^j = (1 - a^m) / (1 - a)$, we compute
  $
    a^m xi(t_c, t_0) + sum_(j=0)^(m-1) a^j &= (a^m (1 - a b) + 1 - a^m) / (1 - a) \
                                           &= (1 - a^(m + 1) b) / (1 - a) .
  $
  Since $a^(m + 1) b = 2^(-lambda (t_c + m T_p - t_0 + T_p))$, the right-hand side equals $xi(t_c + m T_p, t_0)$ identically. Moreover, $a^(m+1) b > 0$ implies $xi(t_c + m T_p, t_0) < 1 / (1 - a)$. By choice of $T_p >= (1 / lambda) log_2 ((beta mu) / (beta mu - 1))$, we have
  $
    a = 2^(-lambda T_p) <= (beta mu - 1) / (beta mu) = 1 - 1 / (beta mu) ,
  $
  which rearranges to $1 / (1 - a) <= beta mu$, proving part 1.

  For part 2, suppose $w(t_c) < xi(t_c, t_0)$ and future arrivals into $c_o$ occur at spacing at least $T_p$. The worst-case (maximal) weight at $t_c + m T_p$ occurs when arrivals occur right at times $t_c + T_p, dots, t_c + m T_p$, giving
  $
    w(t_c + m T_p) <= a^m w(t_c) + sum_(j=0)^(m-1) a^j < a^m xi(t_c, t_0) + sum_(j=0)^(m-1) a^j = xi(t_c + m T_p, t_0) < beta mu .
  $
  Between arrival epochs $w$ decays monotonically, so $w(t) < beta mu$ holds for all $t >= t_c$, validating permanent pruning.
]

#subsec("D", "Exact Geometric Separation Envelope for Offline DBSCAN")

Because DenStream constrains RMS radii rather than worst-case point diameters, one might fear false bridging across separated manifolds. The next theorem proves an exact $4 epsilon$ separation guarantee using only CF-computable RMS radii.

#theorem(title: "Witness Existence and Exact 4-Epsilon Non-Bridging Guarantee", origin: "Original Derivation")[
  Let $P$ be the active $p$-micro-clusters at time $t_c$, each satisfying $r_"rms" (c_p) <= epsilon$.
  1. Every $c_p in P$ contains at least one actual stream point $x_p in c_p$ with $norm(x_p - c_p)_2 <= epsilon$, and for any $gamma > 1$, strictly more than a fraction $1 - gamma^(-2)$ of the decayed mass of $c_p$ lies inside the open ball $B_(gamma epsilon) (c_p)$.
  2. Whenever $c_q$ is directly density-reachable from $c_p$ ($norm(c_p - c_q)_2 <= 2 epsilon$), there exist data points $x_p in c_p$ and $x_q in c_q$ with $norm(x_p - x_q)_2 <= 4 epsilon$.
  3. Contrapositively, if two ground-truth subsets $A, B subset RR^d$ satisfy $"dist"(A, B) = inf_(x in A, y in B) norm(x - y)_2 > 4 epsilon$, then no $p$-micro-cluster supported in $A$ can ever be directly density-reachable from any $p$-micro-cluster supported in $B$.
]

#proof[
  Let $p_j in c_p$ have normalized probabilities $alpha_j = w_j / w(c_p) > 0$ with $sum_j alpha_j = 1$. Then $r_"rms"^2 (c_p) = sum_j alpha_j norm(p_j - c_p)_2^2 <= epsilon^2$. Because the minimum of a finite set cannot exceed any convex combination, there exists $x_p in c_p$ with $norm(x_p - c_p)_2^2 <= epsilon^2$, i.e., $norm(x_p - c_p)_2 <= epsilon$. Furthermore, by Markov's inequality applied to the nonnegative random variable $Y = norm(X - c_p)_2^2$ under law $alpha$,
  $
    "Pr"(norm(X - c_p)_2 >= gamma epsilon) = "Pr"(Y >= gamma^2 epsilon^2) <= (EE[Y]) / (gamma^2 epsilon^2) <= 1 / gamma^2 ,
  $
  proving part 1. For part 2, pick the guaranteed witnesses $x_p in c_p$ and $x_q in c_q$ with $norm(x_p - c_p)_2 <= epsilon$ and $norm(x_q - c_q)_2 <= epsilon$. If $norm(c_p - c_q)_2 <= 2 epsilon$, the triangle inequality gives
  $
    norm(x_p - x_q)_2 <= norm(x_p - c_p)_2 + norm(c_p - c_q)_2 + norm(c_q - x_q)_2 <= 4 epsilon .
  $
  Part 3 follows immediately by contraposition: if $c_p subset.eq A$ and $c_q subset.eq B$ satisfied $norm(c_p - c_q)_2 <= 2 epsilon$, their witnesses $x_p in A, x_q in B$ would violate $"dist"(A, B) > 4 epsilon$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Online per-point cost in DenStream is $O((abs(P) + abs(O)) d)$ for linear nearest-micro-cluster search plus $O(d)$ for the exact rank-1 radius check (Lemma 1). While Theorem 2 bounds $abs(P) <= v / (beta mu (1 - 2^(-lambda)))$, the outlier buffer $abs(O)$ can grow up to $O(v T_p)$ between pruning sweeps under pure noise—making $T_p$ the principal memory governor. Periodic pruning scans $P union O$ in $O(abs(P) + abs(O))$ time every $T_p$ ticks. Offline clustering runs DBSCAN @ester1996dbscan on $abs(P)$ centroids in $O(abs(P)^2 d)$ worst-case (or $O(abs(P) log abs(P))$ with spatial indexing in low $d$).

Geometrically, DenStream excels at tracking arbitrary non-convex manifolds (spirals, crescents, filaments) and automatically determines the number of clusters $k$ while isolating background noise. However, it depends critically on a single global scale parameter $epsilon$ across all clusters and dimensions.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Online Time / Point], [$O((abs(P) + abs(O)) d)$; rank-1 radius test in $O(d)$.],
    [Offline Query Time], [$O(abs(P)^2 d)$ DBSCAN over $p$-micro-clusters (independent of $n$).],
    [Memory Complexity], [$O((abs(P) + abs(O)) d)$ with $abs(P) <= v / (beta mu (1 - 2^(-lambda)))$.],
    [Cluster Count $k$], [Automatic; adapts dynamically as manifolds appear, split, or fade.],
    [Inductive Bias], [Arbitrary non-convex density-connected shapes with explicit noise rejection.],
    [Key Hyperparameters], [Radius $epsilon$, core weight $mu$, outlier factor $beta$, decay rate $lambda$.],
    [High-D Vulnerability], [Severe: fixed $epsilon$ causes exponential merge collapse (Prop. 5).],
  ),
  caption: [Operational and computational profile of DenStream.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Real-Time Maritime, Aviation, and Mobility Trajectory Streams:* GPS pings from vehicles trace curved corridors, shipping lanes, and intersections contaminated by sensor glitches. DenStream's $2 epsilon$-density connectivity reconstructs non-convex traffic arteries while quarantining erroneous GPS jumps inside the $o$-micro-cluster buffer.
2. *Cybersecurity Anomaly and Botnet Detection:* Network traffic contains bursty port scans and DDoS outliers alongside shifting benign workloads @cao2006denstream. The dual $(P, O)$ buffer architecture with provable pruning (Theorem 3) cleanly flags transient anomalous bursts as unpromoted $o$-micro-clusters without corrupting baseline traffic clusters.
3. *Environmental and Geophysical Sensor Networks:* Seismic tremor catalogs and atmospheric aerosol monitors produce non-spherical spatiotemporal clusters whose activity waxes and wanes. Exponential damping ($2^(-lambda Delta t)$) tracks moving swarm fronts while fading dormant fault zones automatically.
4. *Industrial Robotics and LiDAR Obstacle Tracking:* Streaming 2D/3D point clouds of moving obstacles exhibit non-convex footprints and clutter reflections. Fixed-radius micro-clusters compress dense object surfaces in sub-millisecond online steps while offline DBSCAN segments individual obstacles without knowing object counts a priori.
5. *Clickstream Fraud and Micro-Trend Detection:* Emerging fraudulent campaigns start small and rapidly accumulate density; DenStream's $O -> P$ promotion threshold $beta mu$ naturally detects nascent trends as soon as arrival intensity exceeds $1 / T_p$, while pruning one-off noisy sessions.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

DenStream suffers from three structural failure modes. First, *high-dimensional collapse*: as proved in Proposition 5 below, fixed-$epsilon$ balls lose absorption probability exponentially in $d$, trapping all points in singleton $o$-micro-clusters that expire before promotion. Second, *multi-density heterogeneity*: like DBSCAN @ester1996dbscan @schubert2017dbscanrevisited, a single $(epsilon, mu)$ pair either fragments sparse clusters or bridges dense adjacent ones. Third, *unbounded outlier buffer spikes*: if $beta mu$ is close to $1$, $T_p = ceil((1/lambda) log_2 (beta mu / (beta mu - 1)))$ diverges, delaying pruning and inflating $abs(O)$.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [2D/3D Non-Convex Noisy Streams], [Excellent: arbitrary shapes + dynamic $k$ + outlier buffer.], [Tune $epsilon$ via $k$-distance plot; keep $beta mu >= 2$.],
    [High-Dimensional Feature Streams], [Poor: exponential absorption failure in $d$ (Prop. 5).], [Use StreamKM++ @ackermann2012streamkmpp or Contrastive Clustering @li2021contrastive.],
    [Multi-Scale Variable-Density Data], [Weak: single global $(epsilon, mu)$ cannot adapt locally.], [Use HDBSCAN @campello2013hdbscan or OPTICS @ankerst1999optics offline.],
    [Retroactive Multi-Horizon Queries], [Weak: damped window only queries current decay state.], [Use CluStream @aggarwal2003clustream with Pyramidal Time Frame.],
    [Stationary Spherical Big Data], [Moderate: fast, but lacks $(1 plus.minus epsilon)$ SSE guarantee.], [Prefer StreamKM++ @ackermann2012streamkmpp or BICO @fichtenberger2013bico.],
  ),
  caption: [Domain suitability and failure-mode matrix for DenStream.],
)

#sec("V", "Structural Relations and Empirical Behavior")

DenStream generalizes DBSCAN @ester1996dbscan to evolving streams (recovering static DBSCAN on micro-clusters as $lambda -> 0^+$) and directly inspired grid-based D-Stream @chen2007dstream and hierarchical ClusTree @kranen2011clustree. Crucially, its contrast against $k$-Means streaming methods (CluStream, StreamKM++) hinges on the geometry of fixed-radius balls in high dimensions @beyer1999nearest.

#proposition(title: "Exponential High-Dimensional Absorption Collapse under Fixed Radius", origin: "Original Derivation")[
  Let consecutive stream points be drawn i.i.d. from an isotropic Gaussian cluster $X_1, X_2 tilde.op cal(N)(bold(0), sigma^2 bold(I)_d)$ in $RR^d$. For any fixed DenStream radius $epsilon < sigma sqrt(d / 2)$, let $c = 2 epsilon^2 / (d sigma^2) in (0, 1)$. The probability that $X_2$ is absorbed into the singleton $o$-micro-cluster ${X_1}$ (requiring $r_"rms"' <= epsilon$) decays exponentially in dimension $d$:
  $
    "Pr"(r_"rms"' <= epsilon) <= (c e^(1 - c))^(d / 2) <= exp(-d / 4 (1 - c)^2) .
  $
]

#proof[
  At $tilde(w) = 1$ (zero inter-arrival decay, maximizing absorption probability) and $r_"rms" = 0$, Lemma 1 yields $r_"rms"' = (1 / 2) norm(X_2 - X_1)_2$, so absorption occurs iff $norm(X_2 - X_1)_2^2 <= 4 epsilon^2$. Since $X_2 - X_1 tilde.op cal(N)(bold(0), 2 sigma^2 bold(I)_d)$, the rescaled variable $Y = norm(X_2 - X_1)_2^2 / (2 sigma^2)$ follows a chi-squared distribution $chi_d^2$ with moment generating function $EE[e^(-s Y)] = (1 + 2 s)^(-d / 2)$ for $s > 0$, and the absorption event is $Y <= 2 epsilon^2 / sigma^2 = c d$.

  By Markov's inequality for any $s > 0$,
  $
    "Pr"(Y <= c d) = "Pr"(e^(-s Y) >= e^(-s c d)) <= e^(s c d) (1 + 2 s)^(-d / 2) .
  $
  Minimizing the exponent $s c d - (d / 2) log(1 + 2 s)$ sets $1 + 2 s^* = 1 / c$, i.e., $s^* = (1 - c) / (2 c) > 0$. Substituting $s^*$ gives
  $
    "Pr"(Y <= c d) <= exp((d / 2)(1 - c) + (d / 2) log c) = (c e^(1 - c))^(d / 2) .
  $
  Finally, for $u in (0, 1]$, $(1 - u)^2 >= 0$ implies $1 / u >= 2 - u$; integrating from $c$ to $1$ yields $-log c >= (1 - c) + (1 / 2)(1 - c)^2$, hence $1 - c + log c <= -(1 / 2)(1 - c)^2$, proving the second inequality.
]

*Empirical Benchmark Analysis:* The ClusterOverview four-dataset benchmark vividly confirms every theoretical property proved above. DenStream attains Iris Sil = 0.586 (highest among Algs. 192--195), ARI = 0.567; Customers Sil = 0.315; Multi-Shape Sil = 0.793 (ARI = N/A); High-D Sil = N/A (ARI = N/A); and mean runtime = 1.933 ms. On Multi-Shape, DenStream's Silhouette of 0.793 dramatically surpasses Contrastive Clustering (0.611), CluStream (0.609), and StreamKM++ (0.611), validating Theorem 4's non-convex density-connectivity fidelity. Simultaneously, the undefined High-D Silhouette (N/A) is the direct empirical signature of Proposition 5: under fixed $epsilon$ in high dimensions, pairwise distances concentrate above $2 epsilon$, starving $p$-micro-cluster promotion so fewer than two macro-clusters survive. Its 1.933 ms runtime remains 4.8$times$ faster than Contrastive Clustering and slightly faster than CluStream.
