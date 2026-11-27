#import "../template/ieee.typ": *

#alg_entry(
  id: 133,
  name: "STING (Statistical Information Grid)",
  category: "Cat. VI: Grid, Subspace & Projection",
  year: "1997",
  paradigm: "Hierarchical Grid-Based Statistical Query Clustering",
  complexity: $cal(O)(n + K)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Let $cal(X) = {bold(x)_1, dots.h, bold(x)_n} subset Omega$ be a finite dataset of $n$ spatial observations inside a bounded $d$-dimensional axis-aligned domain $Omega = product_(j=1)^d [a_j, b_j] subset RR^d$, where each point $bold(x)_i$ may additionally carry a scalar attribute value $y_i in RR$ (for pure spatial density clustering, $y_i equiv 1$ or spatial coordinates serve directly as attributes). Rather than evaluating pairwise distances or scanning raw observations per query, STING (STatistical INformation Grid) @wang1997sting constructs a multi-resolution hierarchical tree $cal(T)$ of disjoint rectangular cells partitioning $Omega$ across $L$ nested resolution levels $ell in {1, dots.h, L}$, storing compact parametric sufficient statistics in every cell.

At the root ($ell = 1$), a single cell $c_"root" = Omega$ covers the entire domain. Each non-leaf cell $c$ at level $ell < L$ is partitioned into $r$ disjoint rectangular children $"ch"(c)$ (classically $r = 4$ for 2D quadtree subdivisions or $r = 2^d$ for $d$-dimensional orthrees), inducing $K_ell = r^(ell - 1)$ cells at level $ell$ and $K = K_L = r^(L - 1)$ leaf cells at the finest resolution. Let $cal(X)(c) = {bold(x)_i in cal(X) : bold(x)_i in c}$ denote the subset of observations falling inside cell $c$.

#definition(title: "Statistical Information Grid and Query-Relevant Partition")[
  A *STING hierarchy* $(cal(T), theta)$ of depth $L$ and branching factor $r >= 2$ assigns to every cell $c in cal(T)$ a six-tuple of statistical parameters
  $ theta(c) = (n_c, m_c, s_c, "min"_c, "max"_c, cal(D)_c), $
  where $n_c = |cal(X)(c)|$ is the cell cardinality, $m_c$ and $s_c$ are the attribute mean and standard deviation over $cal(X)(c)$, $"min"_c$ and $"max"_c$ are the attribute extrema, and $cal(D)_c in {"Normal", "Uniform", "Exponential", "None"}$ labels the empirical distribution family. Given a spatial-attribute query predicate $Q$ with density threshold $tau > 0$ (or attribute condition $A subset.eq RR$ and selectivity threshold $p_0 in (0, 1]$) and confidence level $1 - alpha$, a cell $c$ is *relevant* if its upper confidence bound $U_c (alpha)$ for the fraction of query-satisfying points satisfies $U_c (alpha) >= p_0$ (equivalently, estimated satisfying count $hat(n)_c^Q >= tau_ell$). The *STING clustering* induced by $(Q, tau)$ is the family of maximal connected components of relevant leaf cells in the leaf adjacency graph $G_L = (cal(R)_L, E_L)$.
]

Leaf statistics $theta(c)$ for $c in cal(T)_L$ are accumulated in a single scan over $cal(X)$. All internal cell parameters ($ell = L - 1, dots.h, 1$) are computed strictly bottom-up from children $"ch"(c) = {c_1, dots.h, c_r}$ without touching the raw data:
$ n_c = sum_(i=1)^r n_(c_i), quad m_c = 1 / n_c sum_(i=1)^r n_(c_i) m_(c_i), $
$ s_c = sqrt(1 / n_c sum_(i=1)^r n_(c_i) (s_(c_i)^2 + m_(c_i)^2) - m_c^2), $
$ "min"_c = min_(1 <= i <= r) "min"_(c_i), quad "max"_c = max_(1 <= i <= r) "max"_(c_i). $
The parent distribution tag $cal(D)_c$ is determined by weighting child distributions by $n_(c_i) / n_c$; if children disagree or share a non-matching family, a $chi^2$ goodness-of-fit statistic compares the merged mixture against the candidate family, falling back to $"None"$ whenever the conflict exceeds a critical threshold @wang1997sting.

At query time, STING descends top-down from a starting layer $ell_0 in {1, dots.h, L}$ (chosen so $K_(ell_0)$ is small yet informative). For each candidate cell $c$, the algorithm estimates the probability $hat(p)_c = "Pr"(y in A | theta(c))$ from $cal(D)_c$ and constructs a $(1 - alpha)$ confidence interval $[L_c, U_c]$ around the satisfying count $n_c hat(p)_c$. Cells with $U_c < tau_ell$ are pruned immediately along with their entire subtrees; surviving cells expand their children until level $L$. Finally, a breadth-first search (BFS) links face-adjacent or Moore-adjacent relevant leaf cells into spatial clusters.

#algbox(
  title: "STING Hierarchical Grid Construction and Query Clustering",
  inputs: [Dataset $cal(X) subset Omega subset RR^d$, hierarchy depth $L$, branching factor $r$, query predicate $Q$, threshold $tau$, confidence $1 - alpha$.],
  outputs: [Connected spatial cluster regions $cal(C) = {C_1, dots.h, C_M}$ of relevant leaf cells.],
)[
  + *Leaf Accumulation:* Partition $Omega$ into $K = r^(L-1)$ leaf cells; in one pass over $cal(X)$, hash each $bold(x)_i$ to its unique leaf $c$ and update $(n_c, m_c, s_c, "min"_c, "max"_c, cal(D)_c)$.
  + *Bottom-Up Rollup:* For $ell = L - 1$ down to $1$ and each cell $c in cal(T)_ell$, aggregate $theta(c)$ from ${theta(c') : c' in "ch"(c)}$ via the exact moment recurrences and the $chi^2$ distribution-merging test.
  + *Top-Down Branch-and-Bound:* Initialize active frontier $cal(F)_(ell_0) <- cal(T)_(ell_0)$. For $ell = ell_0$ to $L - 1$:
    - For each $c in cal(F)_ell$, compute $hat(p)_c$ and confidence interval $[L_c, U_c]$ from $theta(c)$ under $Q$.
    - Mark $c$ relevant iff $U_c >= tau_ell$; set $cal(F)_(ell+1) <- union.big_(c "relevant") "ch"(c)$.
  + *Leaf Connected Components:* Retain $cal(R)_L = {c in cal(F)_L : U_c >= tau_L}$ and run BFS over leaf adjacency graph $G_L = (cal(R)_L, E_L)$ to merge adjacent relevant cells into regions $C_1, dots.h, C_M$.
  + *Return* cluster regions $cal(C) = {C_1, dots.h, C_M}$ and their constituent point sets.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five formal results governing STING's algebraic aggregation exactness, query-time decoupling from dataset cardinality, distribution-free pruning soundness, geometric non-invariance under rigid motions and high dimensions, and streaming order invariance.

#subsec("A", "Exact Bottom-Up Sufficient Statistics Invariance")

#theorem(title: "Zero-Error Hierarchical Moment Aggregation", origin: "Original Proof")[
  Let $cal(T)$ be any finite partition tree over $Omega$ whose non-leaf nodes satisfy $c = union.big_(i=1)^r c_i$ with pairwise disjoint children $c_i inter c_j = emptyset$ ($i != j$). Then for every internal node $c in cal(T)$ at every level $ell in {1, dots.h, L - 1}$, the bottom-up STING recurrences recover the exact raw-data population moments $(n_c, m_c, s_c)$ and extrema $("min"_c, "max"_c)$ of $cal(X)(c)$ with zero approximation error. Furthermore, if unbiased sample variances $tilde(s)_c^2 = (n_c - 1)^(-1) sum_(x in cal(X)(c)) (y_x - m_c)^2$ are desired for $n_c >= 2$, exact unbiasedness is preserved by replacing $n_(c_i) s_(c_i)^2$ with $(n_(c_i) - 1) tilde(s)_(c_i)^2$ and dividing by $n_c - 1$.
]

#proof[
  We proceed by backward induction on tree level $ell = L, L - 1, dots.h, 1$. For the base case $ell = L$, leaf parameters are computed directly from raw observations $cal(X)(c)$, so the claim holds identically.

  For the inductive step, fix an internal cell $c$ at level $ell < L$ with disjoint children $c_1, dots.h, c_r$ at level $ell + 1$, and assume every child tuple $(n_(c_i), m_(c_i), s_(c_i), "min"_(c_i), "max"_(c_i))$ is exact for $cal(X)(c_i)$. Disjointness implies $cal(X)(c) = union.big_(i=1)^r cal(X)(c_i)$ as a disjoint union, immediately yielding $n_c = sum_(i=1)^r n_(c_i)$, $"min"_c = min_i "min"_(c_i)$, and $"max"_c = max_i "max"_(c_i)$. For the first moment,
  $ n_c m_c = sum_(x in cal(X)(c)) y_x = sum_(i=1)^r sum_(x in cal(X)(c_i)) y_x = sum_(i=1)^r n_(c_i) m_(c_i). $
  For the second central moment, expand the raw sum of squared deviations around $m_c$ by inserting $plus.minus m_(c_i)$:
  $ n_c s_c^2 &= sum_(i=1)^r sum_(x in cal(X)(c_i)) ((y_x - m_(c_i)) + (m_(c_i) - m_c))^2 \
    &= sum_(i=1)^r sum_(x in cal(X)(c_i)) (y_x - m_(c_i))^2 + sum_(i=1)^r n_(c_i) (m_(c_i) - m_c)^2 \
    &quad + 2 sum_(i=1)^r (m_(c_i) - m_c) underbrace(sum_(x in cal(X)(c_i)) (y_x - m_(c_i)), = 0). $
  Since the cross term vanishes identically by definition of $m_(c_i)$, we obtain the discrete Law of Total Variance (parallel-axis identity):
  $ n_c s_c^2 = sum_(i=1)^r n_(c_i) s_(c_i)^2 + sum_(i=1)^r n_(c_i) (m_(c_i) - m_c)^2. $
  Expanding $(m_(c_i) - m_c)^2 = m_(c_i)^2 - 2 m_(c_i) m_c + m_c^2$ and using $sum_(i=1)^r n_(c_i) m_(c_i) = n_c m_c$ gives
  $ n_c s_c^2 = sum_(i=1)^r n_(c_i) (s_(c_i)^2 + m_(c_i)^2) - n_c m_c^2, $
  which divides by $n_c$ to recover STING's exact recurrence for $s_c$. Replacing $n_(c_i) s_(c_i)^2$ with $(n_(c_i) - 1) tilde(s)_(c_i)^2$ in the parallel-axis form preserves the exact numerator $sum_(x in cal(X)(c)) (y_x - m_c)^2$, so dividing by $n_c - 1$ preserves exact Bessel-corrected unbiasedness across arbitrarily imbalanced child counts.
]

#subsec("B", "Query-Time Decoupling and Tree-Work Bound")

#theorem(title: "Cardinality-Independent Query Complexity", origin: "Classical — Cited")[
  Following a single $cal(O)(n + K)$ offline construction pass, any top-down STING query over a branching-$r$ hierarchy ($r >= 2$) with $K = r^(L-1)$ leaf cells runs in $cal(O)(K)$ worst-case time, completely independent of dataset cardinality $n$ @wang1997sting. More tightly, if $M_ell <= r^(ell-1)$ cells are deemed relevant at level $ell$ and each leaf has $d_"adj"$ spatial neighbors, total query work obeys
  $ T_"query" <= c_0 (K_(ell_0) + r sum_(ell=ell_0)^(L-1) M_ell) + c_1 d_"adj" M_L <= c_0 (r K - 1) / (r - 1) + c_1 d_"adj" K. $
]

#proof[
  Let $c_0$ be the $cal(O)(1)$ cost of evaluating the parametric CDF/confidence test on a single cell tuple $theta(c)$. At starting level $ell_0$, STING evaluates all $K_(ell_0)$ cells. At every subsequent level $ell + 1 in {ell_0 + 1, dots.h, L}$, only the $r$ children of the $M_ell$ relevant parent cells at level $ell$ enter the frontier $cal(F)_(ell+1)$, so $|cal(F)_(ell+1)| = r M_ell$. Summing over all evaluated nodes yields
  $ N_"eval" = K_(ell_0) + r sum_(ell=ell_0)^(L-1) M_ell <= sum_(ell=1)^L r^(ell - 1) = (r^L - 1) / (r - 1) = (r K - 1) / (r - 1). $
  Since $r >= 2$, $(r K - 1) / (r - 1) < 2 K$. In Step 4, BFS visits each relevant leaf cell $c in cal(R)_L$ once and inspects at most $d_"adj"$ adjacent leaf cells ($2 d$ for face-adjacency, $3^d - 1$ for Moore adjacency), taking at most $c_1 d_"adj" M_L <= c_1 d_"adj" K$ operations. Neither term depends on $n$ because every statistical test reads only the constant-size summary $theta(c)$. For a localized query covering volume fraction $rho in (0, 1)$, $M_ell = cal(O)(rho r^(ell-1) + r^(ell(d-1)/d))$, reducing query work to $cal(O)(rho K + K^((d-1)/d))$.
]

#subsec("C", "Distribution-Free Pruning Soundness via Cantelli's Bound")

A critical subtlety of STING is that replacing a multimodal mixture of $r$ child distributions with a single parametric Normal approximation $cal(N)(m_c, s_c^2)$ at parent level $ell$ can underestimate tail mass and trigger *false pruning* of a dense leaf cell. We derive an exact distribution-free upper bound guaranteeing zero false negatives whenever $cal(D)_c = "None"$.

#proposition(title: "Cantelli Zero-False-Pruning Guarantee and Mixture Dilution", origin: "Original Derivation")[
  Consider a one-sided attribute query $Q: y >= u$ evaluated on a parent cell $c$ with exact moments $(n_c, m_c, s_c)$ and arbitrary unknown empirical distribution.
  + If $u > m_c$, the true fraction $p_c = n_c^(-1) |{x in cal(X)(c) : y_x >= u}|$ satisfies the sharp Cantelli inequality
    $ p_c <= U_c^"Cant" := s_c^2 / (s_c^2 + (u - m_c)^2), $
    and pruning $c$ only when $n_c U_c^"Cant" < tau_L$ guarantees *zero false pruning* of any descendant leaf cell $c^* subset c$ containing at least $tau_L$ satisfying observations.
  + Conversely, under parametric Normal pruning $hat(p)_c^"Norm" = 1 - Phi((u - m_c) / s_c)$, there exists a valid bimodal cell dataset with exact moments $(m_c, s_c)$ whose true tail proportion exceeds $hat(p)_c^"Norm"$ by a factor growing super-exponentially in $z = (u - m_c) / s_c$.
]

#proof[
  For part (1), let $Y$ denote the empirical random variable drawn uniformly from ${y_x : x in cal(X)(c)}$, so $EE[Y] = m_c$ and $"Var"(Y) = s_c^2$. Let $Z = Y - m_c$ and $t = u - m_c > 0$. For any $lambda > 0$,
  $ "Pr"(Y >= u) &= "Pr"(Z + lambda >= t + lambda) <= "Pr"((Z + lambda)^2 >= (t + lambda)^2) \
    &<= (EE[(Z + lambda)^2]) / (t + lambda)^2 = (s_c^2 + lambda^2) / (t + lambda)^2, $
  by Markov's inequality and $EE[Z] = 0$. Minimizing the right-hand side over $lambda > 0$ via calculus sets $dif / (dif lambda) [(s_c^2 + lambda^2) / (t + lambda)^2] = 0$, giving $lambda^* = s_c^2 / t$. Substituting $lambda^*$ yields
  $ p_c <= (s_c^2 + s_c^4 / t^2) / (t + s_c^2 / t)^2 = s_c^2 / (s_c^2 + t^2) = U_c^"Cant". $
  Because every descendant leaf $c^* subset c$ satisfies $|{x in cal(X)(c^*) : y_x >= u}| <= |{x in cal(X)(c) : y_x >= u}| = n_c p_c <= n_c U_c^"Cant"$, having $n_c U_c^"Cant" < tau_L$ strictly implies no leaf under $c$ can reach $tau_L$, proving zero false pruning.

  For part (2), construct a two-point empirical measure on $c$ placing mass $p = (1 + z^2)^(-1)$ at $y_+ = m_c + z s_c$ and mass $1 - p$ at $y_- = m_c - s_c / z$ for any $z > 0$. Direct calculation confirms $EE[Y] = m_c$, $"Var"(Y) = s_c^2$, and exact tail probability $p_c = "Pr"(Y >= m_c + z s_c) = (1 + z^2)^(-1)$ (achieving Cantelli's bound with equality). By contrast, Mill's ratio bounds the Gaussian tail by $hat(p)_c^"Norm" <= (z sqrt(2 pi))^(-1) exp(-z^2 / 2)$. The ratio $p_c / hat(p)_c^"Norm" >= (z sqrt(2 pi) / (1 + z^2)) exp(z^2 / 2)$ diverges super-exponentially as $z -> infinity$; thus any Normal-tail pruner drops $c$ while a leaf holding all $n_c / (1 + z^2)$ upper-mode points still exceeds $tau_L$.
]

#subsec("D", "Grid-Phase Aliasing, Diagonal Staircasing, and Dimensionality")

#theorem(title: "Non-Isotropy, Boundary-Splitting Erasure, and Curse of Dimensionality", origin: "Original Counterexample")[
  Let STING operate on a uniform leaf grid of side length $delta > 0$ in $RR^d$ with leaf density threshold $tau_L > 0$.
  + *Translation / Phase Erasure:* There exist two isometric datasets $cal(X), cal(X)' subset RR^2$ differing solely by a rigid translation of $delta / 2$ such that STING outputs a single intact cluster on $cal(X)$ and the empty clustering $emptyset$ on $cal(X)'$.
  + *Rotational Staircasing:* Any smooth cluster boundary segment of Euclidean length $L_b$ inclined at angle $pi / 4$ to the coordinate axes incurs Hausdorff approximation error at least $delta / (2 sqrt(2))$ and grid-boundary perimeter inflation by exact factor $sqrt(2)$.
  + *Dimensionality Barrier:* Maintaining $m >= 2$ intervals per axis across $d$ dimensions forces $K = m^d$ leaf cells, so average leaf occupancy $n / m^d -> 0$ exponentially in $d$, driving spatial connectivity to collapse when $d > log_m (n / tau_L)$.
]

#proof[
  For (1), fix $delta > 0$ and any integer threshold $tau_L >= 2$. Place $n_0 = tau_L$ observations tightly inside a ball $B((delta/2, delta/2), epsilon)$ with $0 < epsilon < delta / 4$, centered inside the leaf cell $c_(0,0) = [0, delta)^2$. Every point falls in $c_(0,0)$, yielding $n_(c_(0,0)) = tau_L >= tau_L$, so STING returns $C_1 = {c_(0,0)}$. Now translate $cal(X)$ rigidly by $(delta/2, 0)$ so the ball is centered at $(delta, delta/2)$ on the shared vertical boundary between $c_(0,0)$ and $c_(1,0) = [delta, 2 delta) times [0, delta)$, splitting the points evenly so $floor(tau_L / 2)$ lie left of $x_1 = delta$ and $ceil(tau_L / 2)$ lie right. For every integer $tau_L >= 2$, $ceil(tau_L / 2) <= (tau_L + 1) / 2 < tau_L$, so both leaf cells fail the threshold and STING returns $emptyset$, violating both translation invariance and Kleinberg's consistency intuition @kleinberg2002impossibility.

  For (2), consider the half-space boundary $x_1 + x_2 = beta$ of length $L_b$ crossing a square grid of mesh $delta$. Each diagonal step advancing $sqrt(2) delta$ along the line traverses one horizontal and one vertical grid edge of total Manhattan length $2 delta$, giving exact perimeter ratio $2 delta / (sqrt(2) delta) = sqrt(2)$, while the midpoint of any bisected grid cell lies at perpendicular distance $delta / (2 sqrt(2))$ from the nearest cell vertex on the opposite side of the line.

  For (3), an axis-parallel grid with $m$ subdivisions per dimension partitions $Omega$ into $K = m^d$ disjoint leaves. By the pigeonhole principle, the mean count per leaf is $bar(n)_L = n / m^d$. For any uniform or diffuse high-dimensional cluster supported on a constant fraction $gamma in (0, 1]$ of $Omega$, expected leaf counts fall below $tau_L$ as soon as $gamma n / m^d < tau_L$, i.e., $d > log_m (gamma n / tau_L)$, making dense-leaf percolation impossible without exponentially large $n$ @beyer1999nearest.
]

#subsec("E", "Permutation Invariance and Exact Streaming Maintenance")

#proposition(title: "Order-Invariance and Online Insertion Complexity", origin: "Original Proof")[
  For any fixed grid geometry $cal(T)$ over $Omega$, the statistical summary tree $theta(cal(T))$ and all downstream STING query outputs are strictly invariant under every permutation $pi in S_n$ of input arrival order. Moreover, inserting a new point $(bold(x)_(n+1), y_(n+1)) in Omega$ updates the entire hierarchy exactly in $cal(O)(L)$ time using Welford-Chan running recurrences.
]

#proof[
  Each observation $(bold(x)_i, y_i)$ maps deterministically to a unique leaf cell $c(bold(x)_i)$ via coordinate floor indexing $u_j = floor((x_(i,j) - a_j) / delta_j)$, independent of other points. Within each leaf, $n_c$, $sum y_x$, $sum y_x^2$, $min y_x$, and $max y_x$ are commutative associative monoid reductions, hence invariant under all $pi in S_n$. By Theorem II.A, all ancestor nodes depend solely on leaf tuples, establishing global permutation invariance (unlike order-dependent CF-trees in BIRCH @zhang1996birch).

  Upon arrival of $(bold(x)_(n+1), y_(n+1))$, the unique root-to-leaf path contains exactly $L$ cells. For each cell $c$ on this path, setting $n'_c = n_c + 1$, $delta_y = y_(n+1) - m_c$, $m'_c = m_c + delta_y / n'_c$, and $M'_(2,c) = M_(2,c) + delta_y (y_(n+1) - m'_c)$ (where $M_(2,c) = n_c s_c^2$) updates all moments in $cal(O)(1)$ per level without catastrophic cancellation, giving $cal(O)(L)$ total update time.
]

#sec("III", "Complexity, Geometry and Operational Profile")

STING separates offline spatial summarization from online query-driven clustering. Given $n$ points and $K = r^(L-1)$ leaf cells, computing leaf cell indices via integer arithmetic takes $cal(O)(n d)$ time and rolling up the $(K - 1) / (r - 1)$ internal nodes takes $cal(O)(K)$ time, yielding $cal(O)(n d + K)$ build cost (simply $cal(O)(n + K)$ for fixed low dimension $d$). Because the hierarchy is static and array-indexable via Morton (Z-order) space-filling curves, all tree levels reside in contiguous memory buffers with high L1/L2 cache locality and zero pointer-chasing overhead. Storing 6 scalars plus a 1-byte distribution tag per cell requires roughly $49$ bytes per node, so $K = 4^8 = 65536$ leaf cells (plus $21845$ internal nodes) occupies only $approx 4.3$ MB of RAM regardless of whether $n = 10^4$ or $n = 10^9$.

Geometrically, STING induces *isothetic (axis-parallel) orthogonal polytope* clusters formed by unions of congruent leaf hyper-rectangles. It handles arbitrary non-convex, elongated, and multi-connected shapes effortlessly in 2D/3D (subject to the $delta$-mesh staircasing of Theorem II.D) and inherently filters sparse noise points whose cells remain below $tau_L$. Its primary trade-off is that resolution $L$ must be fixed *a priori* and memory scales as $m^d$.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [*Property*], [*Characterization*],
    [Build Complexity], [$cal(O)(n d + K)$ single-pass leaf hashing + bottom-up rollup],
    [Query Complexity], [$cal(O)(K)$ worst-case; $cal(O)(rho K + K^((d-1)/d))$ localized; independent of $n$],
    [Space Complexity], [$cal(O)(K)$ total ($approx r K / (r - 1)$ nodes, constant 6-tuple per cell)],
    [Cluster Geometry], [Non-convex unions of axis-parallel leaf hyper-rectangles],
    [Order / Parallelism], [Strictly permutation-invariant; embarrassingly parallel map-reduce build],
    [Noise & Outliers], [Implicitly rejected via cell density/confidence pruning ($U_c < tau_ell$)],
    [Primary Hyperparams], [Hierarchy depth $L$, branching factor $r$, query threshold $tau$, confidence $1 - alpha$],
    [Key Limitations], [Exponential grid growth $m^d$; axis-aligned staircasing; grid-phase sensitivity],
  ),
  caption: [Operational and computational profile of STING @wang1997sting.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

*1. Geographic Information Systems (GIS) and Urban Spatial Databases:* Urban planning queries (e.g., locating contiguous residential blocks whose mean property value exceeds a threshold with 95% confidence) match STING's two-layer spatial-plus-attribute architecture directly @wang1997sting. Quadtree cells ($r = 4$) map natively onto map tiles, enabling millisecond interactive OLAP region rollups over tens of millions of land parcels without re-scanning raw records.

*2. Remote Sensing and Satellite Imagery Tiling:* Multispectral Earth-observation archives organize reflectance measurements over fixed latitude-longitude grids. Because STING's bottom-up moments are exact (Theorem II.A) and map-reduce additive, coarse-to-fine pruning rapidly isolates agricultural stress zones, wildfire thermal anomalies, or urban heat islands across multi-terabyte rasters in $cal(O)(K)$ query time.

*3. Astronomical Sky-Survey Catalog Indexing:* Wide-field sky surveys (such as SDSS or Gaia) record billions of celestial objects indexed by right ascension and declination alongside photometric attributes (luminosity, redshift). Hierarchical spatial grids with stored cell-level moments enable rapid identification of dense stellar clusters and galaxy overdensities while filtering diffuse sky background via confidence bounds.

*4. Meteorological and Environmental Sensor Telemetry:* Distributed weather stations, air-quality monitors, and oceanographic buoy arrays stream continuous scalar readings (temperature, $"PM"_(2.5)$, salinity) at fixed geographic coordinates. By Proposition II.E, STING accommodates incremental streaming updates in $cal(O)(L)$ time per reading and answers real-time regional hazard queries without reconstructing spatial indices.

*5. Large-Scale Spatial OLAP and Database Query Optimization:* Inside relational/columnar analytical engines, STING acts as a materialized multi-resolution summary index @xu2005survey. Analysts issuing dozens of exploratory clustering predicates with varying thresholds $tau$ amortize the $cal(O)(n)$ construction cost across all subsequent $cal(O)(K)$ queries, outperforming repeated runs of DBSCAN @ester1996dbscan by orders of magnitude when $n >> K$.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

STING fails predictably in three mathematical regimes: (i) *high-dimensional feature spaces* ($d >= 6$), where $K = m^d$ either exhausts memory or produces sparse occupancy $n / m^d << 1$ (Theorem II.D); (ii) *oblique or finely separated clusters* whose inter-cluster gap is smaller than $delta$, which merge across coarse cells or suffer boundary-splitting erasure; and (iii) *heavy-tailed multimodal parent cells* where parametric Normal pruning violates recall unless distribution-free Cantelli bounds (Proposition II.C) are substituted.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [*Domain*], [*Mathematical Fit*], [*Caution / Alternative*],
    [2D/3D GIS & Spatial OLAP], [Exact quadtree/octree moment rollup; $cal(O)(K)$ repeated queries], [Grid-phase shifts; refine edges via Grid-DBSCAN @mahran2008gridbscan],
    [Satellite Raster & Sky Surveys], [Single-pass $cal(O)(n)$ build; native multi-scale thresholding], [Oblique filaments staircase; use WaveCluster @sheikholeslami1998wavecluster],
    [Streaming Sensor Telemetry], [$cal(O)(L)$ exact Welford updates; permutation-invariant], [Temporal concept drift; use D-Stream @chen2007dstream],
    [High-D Tabular / Omics ($d >= 10$)], [Poor: $K = m^d$ sparsity destroys leaf connectivity], [Use CLIQUE @agrawal1998clique or PROCLUS @aggarwal1999proclus],
    [Variable-Density Point Clouds], [Single global $tau_L$ misses sparse or merges dense clusters], [Use OPTIGRID @hinneburg1999optigrid or HDBSCAN],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for STING.],
)

#sec("V", "Structural Relations and Empirical Behavior")

STING occupies the foundational vertex of grid-based clustering, connecting directly to density-based and wavelet-based paradigms through leaf-level density superlevel sets.

#proposition(title: "Equivalence to Grid-Quantized DBSCAN and Zero-Pass WaveCluster", origin: "Original Proof")[
  Consider pure spatial density clustering ($y_i equiv 1$) on a uniform leaf grid of mesh $delta$ in $RR^d$ with leaf threshold $tau_L = "MinPts"$ and Moore neighborhood adjacency ($L_infinity$-distance $1$ between cells).
  + Every STING leaf cluster $C subset cal(R)_L$ is a subgraph of a DBSCAN clustering @ester1996dbscan with radius $epsilon = 2 sqrt(d) delta$ and threshold $"MinPts" = tau_L$, and conversely every DBSCAN core-point component with radius $epsilon' = delta / 2$ induces a connected union of dense STING cells when shifted into cell interiors @mahran2008gridbscan.
  + Restricting STING's hierarchy to $r = 2^d$ dyadic cells with exact counts $n_c$ at level $L - s$ is algebraically identical to the $s$-level Haar Discrete Wavelet Transform approximation subband of WaveCluster @sheikholeslami1998wavecluster scaled by $2^(-s d / 2)$, with zero high-pass detail filtering.
]

#proof[
  For (1), if two leaf cells $c, c' in cal(R)_L$ are Moore-adjacent, the maximum Euclidean distance between any $bold(x) in c$ and $bold(x)' in c'$ is at most $sqrt(sum_(j=1)^d (2 delta)^2) = 2 sqrt(d) delta$. Because each cell holds at least $tau_L$ points, every point in $c union c'$ has at least $tau_L$ neighbors within radius $2 sqrt(d) delta$, placing $cal(X)(c) union cal(X)(c')$ in a single DBSCAN density-connected component. Conversely, any ball of radius $delta / 2$ centered at a cell midpoint lies strictly inside that leaf, so a DBSCAN core point with $epsilon' = delta / 2$ guarantees $n_c >= tau_L$.

  For (2), the 1D unitary Haar scaling filter $h = (1 / sqrt(2), 1 / sqrt(2))$ sums adjacent dyadic pairs and divides by $sqrt(2)$. Applying $h$ separably across $d$ axes for $s$ scales computes $A_s [bold(v)] = 2^(-s d / 2) sum_(bold(u) in "block"(bold(v))) f[bold(u)]$, which equals $2^(-s d / 2) n_c$ for the level-$(L - s)$ STING ancestor cell $c$ covering the $2^(s d)$ leaf block. Thus thresholding $A_s [bold(v)] >= 2^(-s d / 2) tau_(L-s)$ produces the exact same relevant cells and BFS components as STING at level $L - s$.
]

*Related Grid Architectures:* BANG @schikuta1998bang replaces STING's static uniform cells with an adaptive balanced grid-file tree to prevent empty-cell explosion in skewed regions; OPTIGRID @hinneburg1999optigrid chooses optimal non-uniform cutting hyperplanes along low-density valleys to handle moderate dimensions; D-Stream @chen2007dstream generalizes STING's static cell counts to exponentially decayed density grids for evolving data streams; and CLIQUE @agrawal1998clique (Alg. 134) restricts grid density search to lower-dimensional subspaces to overcome Theorem II.D(3) @aggarwal2013data @jain2010data.

*Empirical Benchmark Behavior:* On the ClusterOverview 4-dataset benchmark evaluated via Silhouette @rousseeuw1987silhouettes and Adjusted Rand Index (ARI) @hubert1985comparing, STING achieves *Iris* $"Sil" = 0.517$, $"ARI" = 0.527$; *Customers* $"Sil" = 0.449$; *Multi-Shape* $"Sil" = 0.593$ ($"ARI" = "N/A"$); and *High-D* $"Sil" = 0.257$ ($"ARI" = "N/A"$), with a mean runtime of just *$0.052$ ms* — the fastest across the grid/subspace quartet (Algs. 133–136). These measurements mirror our theoretical predictions:
- *Non-Convex Dominance & Speed:* STING records its highest quality on *Multi-Shape* ($"Sil" = 0.593$), where BFS over dense leaf cells faithfully traces curved non-convex geometry in $cal(O)(n + K)$ time without iterative distance evaluations.
- *Tabular Density Fit:* On *Customers* ($"Sil" = 0.449$), STING outperforms CLIQUE ($0.176$), WaveCluster ($0.353$), and PROCLUS ($0.273$), reflecting the clean preservation of sharp cell boundaries without wavelet smoothing blur.
- *High-Dimensional Degradation:* On *High-D*, Silhouette drops sharply to $0.257$, empirically confirming Theorem II.D(3)'s sparsity decay in full-space grids (lagging PROCLUS's $0.294$). Notably, STING matches WaveCluster's scores identically on *Iris* ($0.517 / 0.527$), *Multi-Shape* ($0.593$), and *High-D* ($0.257$), directly instantiating Proposition V(2)'s Haar-equivalence regime when wavelet filtering leaves the active leaf superlevel set intact.
