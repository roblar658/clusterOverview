#import "../template/ieee.typ": *

#alg_entry(
  id: 50,
  name: "VDBSCAN: Varied-Density Spatial Clustering with Adaptive Radii",
  category: "Cat. II: Density & Mode-Seeking",
  year: "2007",
  paradigm: "Multi-Scale k-Distance Level-Set Peeling",
  complexity: $cal(O)(n^2 + L n)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Classical DBSCAN @ester1996dbscan estimates connected components of a single upper density level set ${x in bb(R)^d : f(x) >= lambda}$ using a global neighborhood radius $epsilon > 0$ and cardinality threshold $"MinPts"$. While effective on uniform-density manifolds, a single global radius faces an irreconcilable trade-off whenever a dataset $cal(X) = {x_1, dots, x_n} subset bb(R)^d$ contains clusters spanning multiple density regimes $f_1 > f_2 > dots > f_L$: setting $epsilon$ small enough to prevent merging closely spaced dense clusters shatters sparse clusters into noise, whereas enlarging $epsilon$ to recover sparse clusters coalesces adjacent dense clusters into a single component @schubert2017dbscanrevisited.

Varied-Density DBSCAN (VDBSCAN) @liu2007vdbscan addresses this limitation by extracting a discrete, strictly increasing sequence of $L$ neighborhood radii $0 < epsilon_1 < epsilon_2 < dots < epsilon_L$ directly from the empirical $k$-nearest-neighbor distance distribution ($k$-dist plot), and then executing a sequential *density-peeling* (deflation) cascade. At each stage $ell = 1, dots, L$, DBSCAN runs with radius $epsilon_ell$ exclusively on the *unclustered* residual subset $U_ell subset.eq cal(X)$; points assigned to clusters at density stratum $ell$ are frozen and removed before coarser radii $epsilon_(ell+1), dots, epsilon_L$ are applied.

Let $d(dot, dot)$ be a metric on $bb(R)^d$. For a fixed integer $k >= 2$, let $x_((k)) (x) in cal(X) without {x}$ denote the $k$-th nearest neighbor of $x$ in $cal(X)$, and define the *$k$-distance* function
$
  d_k (x) = d(x, x_((k)) (x)), \
  s_1 <= s_2 <= dots <= s_n,
$
where $(s_1, dots, s_n)$ sorts ${d_k (x_i)}_(i=1)^n$ in non-decreasing order (equivalently, the order-reversed $k$-dist curve of @ester1996dbscan @liu2007vdbscan). Regions of approximately constant local density manifest as flat horizontal *plateaus* on the sorted $k$-dist curve, separated by abrupt upward steps (knees/elbows) marking transitions to sparser strata or background noise.

#definition(title: "Multi-Scale k-Dist Strata and Sequential Deflation Operator")[
  Given sorted $k$-distances $s_1 <= dots <= s_n$, an *adaptive radius schedule* $cal(E) = (epsilon_1, dots, epsilon_L)$ with $0 < epsilon_1 < dots < epsilon_L$ is obtained by locating $L$ prominent curvature knees $1 <= i_1 < dots < i_L <= n$ of the discrete $k$-dist profile (e.g., via local maxima of normalized second differences $Delta^2 s_i = s_(i+1) - 2 s_i + s_(i-1)$ or piecewise-linear regression change points) and setting $epsilon_ell = s_(i_ell)$. Initializing the unclustered pool $U_1 = cal(X)$ and threshold $"MinPts" = k + 1$, the *VDBSCAN peeling operator* iteratively computes, for $ell = 1, dots, L$:
  $
    (cal(C)_ell, N_ell) &= "DBSCAN"(U_ell, epsilon_ell, "MinPts"), \
    U_(ell+1) &= U_ell without (union.big_(C in cal(C)_ell) C),
  $
  where $cal(C)_ell$ is the collection of clusters discovered in $U_ell$ at scale $epsilon_ell$ and $N_ell = U_(ell+1)$ is the residual unclustered set. The final clustering is $cal(C) = union.big_(ell=1)^L cal(C)_ell$, with final noise set $cal(N) = U_(L+1)$.
]

Crucially, because $epsilon_1 < epsilon_2 < dots < epsilon_L$ proceeds from *highest* density to *lowest* density while deleting clustered points after every pass, dense clusters are isolated before the search radius ever expands enough to bridge them.

#algbox(
  title: "VDBSCAN: Varied-Density Spatial Clustering",
  inputs: [Dataset $cal(X) = {x_1, dots, x_n} subset bb(R)^d$, neighbor order $k >= 2$ (setting $"MinPts" = k + 1$), knee-sensitivity threshold $tau > 0$ (or stratum count $L$).],
  outputs: [Multi-density cluster collection $cal(C) = union.big_(ell=1)^L cal(C)_ell$ and noise set $cal(N)$.],
)[
  1. *Compute $k$-Distance Profile:* For every $x_i in cal(X)$, compute pairwise distances (or query a spatial index) to obtain its $k$-nearest distance $d_k (x_i)$ and sort into $s_1 <= s_2 <= dots <= s_n$.
  2. *Detect Multi-Density Radii:* Identify $L$ sharp upward slope transitions (knees) $i_1 < dots < i_L$ on the sorted $k$-dist curve exceeding curvature threshold $tau$, setting $epsilon_ell = s_(i_ell)$ for $ell = 1, dots, L$.
  3. *Initialize Active Pool:* Set $U_1 <- cal(X)$ and $cal(C) <- emptyset$.
  4. *Iterative Density Peeling:* For $ell = 1, dots, L$:
     (a) Run sequential DBSCAN on $U_ell$ with parameters $(epsilon_ell, "MinPts")$, yielding clusters $cal(C)_ell$.
     (b) Append each cluster $C in cal(C)_ell$ to $cal(C)$ with a fresh global cluster identifier.
     (c) Deflate active pool: $U_(ell+1) <- U_ell without (union.big_(C in cal(C)_ell) C)$.
     (d) If $|U_(ell+1)| < "MinPts"$, terminate early.
  5. *Return:* Output $cal(C)$ and label remaining points in $cal(N) = U_(L+1)$ as noise.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish four mathematical properties of VDBSCAN: finite-sample concentration of $k$-dist plateaus around density strata, a sufficient asymmetric separation theorem proving exact multi-density recovery and strict dominance over single-scale DBSCAN, an explicit counterexample exposing halo-erosion non-monotonicity, and exact similarity-group equivariance alongside high-dimensional plateau collapse.

#subsec("A", "Concentration of k-Distance Plateaus Around Density Strata")

#lemma(title: "Finite-Sample k-Distance Concentration on Constant-Density Strata", origin: "Original Derivation")[
  Let $X_1, dots, X_n$ be drawn i.i.d. from a density $f$ on $bb(R)^d$, and suppose a compact cluster support $Omega_ell subset bb(R)^d$ has constant interior density $f(y) = f_ell > 0$ on $B_R (x) subset Omega_ell$ for an interior point $x in Omega_ell$. Let $V_d = pi^(d/2) / Gamma(d/2 + 1)$ be the unit Euclidean ball volume, and define the characteristic $k$-NN radius
  $
    r_ell = ( k / ((n - 1) V_d f_ell) )^(1/d).
  $
  If $r_ell (1 + eta) <= R$ for some $eta in (0, 1)$, then with $eta_- = 1 - (1 - eta)^d in (0, 1)$ the $k$-distance $d_k (x)$ satisfies
  $
    "Pr"( |d_k (x) / r_ell - 1| > eta ) <= 2 exp( - (k eta_-^2) / 2 ).
  $
  Consequently, adjacent density plateaus $f_ell > f_(ell+1)$ exhibit asymptotic $k$-dist jump ratio $r_(ell+1) / r_ell = (f_ell / f_(ell+1))^(1/d)$.
]
#proof[
  Condition on $X_1 = x$. For any radius $r <= R$, the count $M(r) = sum_(j=2)^n upright(I) [||X_j - x||_2 <= r]$ follows $"Binomial"(n - 1, p(r))$ with $p(r) = V_d f_ell r^d$ because $B_r (x) subset Omega_ell$. By definition of $d_k (x)$,
  $
    {d_k (x) > r_+} = {M(r_+) < k}, \
    {d_k (x) < r_-} = {M(r_-) >= k}.
  $
  Set $r_+ = (1 + eta) r_ell$ and $r_- = (1 - eta) r_ell$. Using $(n - 1) V_d f_ell r_ell^d = k$, the expected ball counts are $mu_+ = k (1 + eta)^d$ and $mu_- = k (1 - eta)^d = k (1 - eta_-)$.

  For the upper tail of $d_k (x)$, note that $(1 - eta)(1 + eta) = 1 - eta^2 < 1$ implies $(1 + eta)^d > (1 - eta)^(-d) = 1 / (1 - eta_-)$, so $k < (1 - eta_-) mu_+$ and $mu_+ > k$. By the multiplicative Chernoff lower-tail bound for binomial sums,
  $
    "Pr"(M(r_+) < k) &<= "Pr"(M(r_+) <= (1 - eta_-) mu_+) \
    &<= exp( - (mu_+ eta_-^2) / 2 ) <= exp( - (k eta_-^2) / 2 ).
  $
  For the lower tail of $d_k (x)$, write $k = mu_- (1 + delta_-)$ with $delta_- = eta_- / (1 - eta_-) > 0$. Applying the unconditional Chernoff upper-tail bound $"Pr"(M >= (1 + delta) mu) <= exp(- mu delta^2 / (2 + delta))$ (valid for all $delta > 0$) gives
  $
    "Pr"(M(r_-) >= k) &<= exp( - (mu_- delta_-^2) / (2 + delta_-) ) \
    &= exp( - (k eta_-^2) / (2 - eta_-) ) <= exp( - (k eta_-^2) / 2 ),
  $
  using $mu_- delta_- = k eta_-$ and $(1 - eta_-)(2 + delta_-) = 2 - eta_- < 2$. A union bound over the two tails yields the stated concentration inequality. Taking ratios of $r_(ell+1)$ and $r_ell$ cancels $k / ((n - 1) V_d)$ and yields $(f_ell / f_(ell+1))^(1/d)$.
]

Lemma 1 explains both why VDBSCAN works in low dimensions (when $k = Omega(log n)$, points in stratum $Omega_ell$ concentrate tightly around $r_ell$, producing flat plateaus) and why contrast attenuates exponentially in $d$: a $100 times$ density gap produces a $10 times$ plateau step in $d = 2$ but only a $100^(1/20) approx 1.26 times$ step in $d = 20$.

#subsec("B", "Asymmetric Separation Guarantee and Dominance Over DBSCAN")

#theorem(title: "Exact Multi-Density Recovery Under Asymmetric Separation", origin: "Original Proof")[
  Let $cal(X) = C_1 union C_2 union dots union C_L$ be a disjoint union of $L$ finite subsets with schedule $0 < epsilon_1 < epsilon_2 < dots < epsilon_L$ and threshold $"MinPts" >= 2$. Suppose that for every stratum $ell in {1, dots, L}$:
  1. *(Intra-Stratum Core Connectivity)* Every point $x in C_ell$ is $(epsilon_ell, "MinPts")$-core within $C_ell$, and the core graph $G_(epsilon_ell) (C_ell)$ is connected.
  2. *(Strict Stratum Density Ordering)* For every coarser stratum $m > ell$ and every $y in C_m$, $|N_(epsilon_ell) (y) inter (union.big_(q >= ell) C_q)| < "MinPts"$.
  3. *(Asymmetric Inter-Stratum Separation)* For every $m > ell$, $"dist"(C_ell, C_m) = min_(x in C_ell, y in C_m) d(x, y) > epsilon_ell$.
  Then VDBSCAN with schedule $(epsilon_1, dots, epsilon_L)$ recovers the exact ground-truth partition ${C_1, dots, C_L}$ with zero noise.
]
#proof[
  We prove by induction on $ell = 1, dots, L$ that at the start of stage $ell$, the active pool equals $U_ell = union.big_(q=ell)^L C_q$, and stage $ell$ extracts $cal(C)_ell = {C_ell}$.

  The base case $ell = 1$ holds immediately since $U_1 = cal(X) = union.big_(q=1)^L C_q$. Assume $U_ell = union.big_(q=ell)^L C_q$ at the beginning of stage $ell$. Consider the $(epsilon_ell, "MinPts")$-DBSCAN execution on $U_ell$:
  - By condition (1), every $x in C_ell$ has $|N_(epsilon_ell) (x) inter U_ell| >= |N_(epsilon_ell) (x) inter C_ell| >= "MinPts"$, so all points of $C_ell$ are core points in $U_ell$ and form a single core-connected component.
  - By condition (3), for every $x in C_ell$ and $y in U_ell without C_ell = union.big_(m > ell) C_m$, we have $d(x, y) >= "dist"(C_ell, C_m) > epsilon_ell$. Thus no point of $U_ell without C_ell$ lies in $N_(epsilon_ell) (x)$ for any $x in C_ell$, so $C_ell$ neither merges with nor absorbs any border point from $union.big_(m > ell) C_m$.
  - By condition (2), no point $y in union.big_(m > ell) C_m$ attains $"MinPts"$ neighbors in $U_ell$ within radius $epsilon_ell$, so no point outside $C_ell$ is a core point at stage $ell$.

  Consequently, the sole core component at stage $ell$ is $C_ell$, and no outside point is density-reachable from $C_ell$. Hence $cal(C)_ell = {C_ell}$, and the deflation step yields $U_(ell+1) = U_ell without C_ell = union.big_(q=ell+1)^L C_q$, completing the induction. At termination, $U_(L+1) = emptyset$ and $cal(C) = {C_1, dots, C_L}$.
]

Notice the remarkable asymmetry in condition (3): separation only needs to exceed the *smaller* radius $epsilon_ell$, NOT the coarser radius $epsilon_m$! This produces a strict representability advantage over standard DBSCAN.

#corollary(title: "Strict Representability Dominance Over Single-Scale DBSCAN", origin: "Original Proof")[
  For every $"MinPts" = m >= 2$ and every $0 < epsilon_1 < Delta < epsilon_2$, there exists a finite dataset $cal(X) = C_1 union C_2 subset bb(R)^(m-1)$ that VDBSCAN recovers identically in two stages, whereas *no* single global radius $epsilon > 0$ allows standard DBSCAN to recover ${C_1, C_2}$.
]
#proof[
  In $bb(R)^(m-1)$, place $C_1$ as the $m$ vertices of a regular simplex of edge length $r_1 = epsilon_1 / 2$, and place $C_2$ as the $m$ vertices of a regular simplex of edge length $r_2 = (Delta + epsilon_2) / 2 in (Delta, epsilon_2)$, translated so that the closest inter-cluster distance equals $"dist"(C_1, C_2) = Delta$.

  All three hypotheses of Theorem 1 hold for schedule $(epsilon_1, epsilon_2)$: (1) every pairwise distance inside $C_1$ equals $r_1 < epsilon_1$ and inside $C_2$ equals $r_2 < epsilon_2$, making both simplices core-connected at their respective scales; (2) for every $y in C_2$, its nearest distinct point in $cal(X)$ lies at distance $min(Delta, r_2) = Delta > epsilon_1$, so $|N_(epsilon_1) (y)| = 1 < m$; and (3) $"dist"(C_1, C_2) = Delta > epsilon_1$. Hence VDBSCAN extracts $C_1$ at stage 1, removes $C_1$, and recovers $C_2$ at stage 2 identically.

  By contrast, consider standard DBSCAN with *any* single global radius $epsilon > 0$:
  - *Case 1 ($epsilon < r_2$):* No two distinct points of $C_2$ lie within distance $epsilon$, so $|N_epsilon (y) inter C_2| = 1$ for every $y in C_2$. Thus $C_2$ has no internal edges and either dissolves into noise (if $epsilon < Delta$) or partially attaches to $C_1$ as border points (if $Delta <= epsilon < r_2$).
  - *Case 2 ($epsilon >= r_2$):* Since $r_2 > Delta = "dist"(C_1, C_2)$, every point of $C_1$ and $C_2$ is globally core and at least one pair $(x, y) in C_1 times C_2$ satisfies $d(x, y) = Delta < r_2 <= epsilon$, forcing DBSCAN to merge $C_1$ and $C_2$ into a single cluster.
  Hence no single $epsilon > 0$ recovers ${C_1, C_2}$ in standard DBSCAN.
]

#subsec("C", "Halo Erosion and Violation of Level-Set Nesting")

Because VDBSCAN deletes clusters $cal(C)_ell$ from $U_ell$ after each pass, its behavior diverges fundamentally from hierarchical level-set trees @hartigan1981consistency @campello2013hdbscan.

#proposition(title: "Concentric Halo Erosion Under Sequential Peeling", origin: "Original Counterexample")[
  For any $"MinPts" = m >= 4$ and $0 < epsilon_1 < epsilon_2$, there exists a 1D dataset $cal(X) = C_1 union C_2$ consisting of a dense core $C_1$ and an adjacent sparser halo $C_2$ at distance $<= epsilon_1$ from $C_1$ such that stage 1 absorbs part of $C_2$ as border points of $C_1$, causing $C_2 without C_1^+$ to lose core status and collapse to noise at stage 2.
]
#proof[
  Fix $m >= 4$ and $0 < epsilon_1 < epsilon_2$ with $epsilon_2 > 2 epsilon_1$. Place $C_1$ as $m$ points at $0$, and place $C_2 = {y_1, dots, y_m}$ where $y_1 = epsilon_1$ and $y_2, dots, y_m$ lie at $epsilon_1 + epsilon_2 / 2$ (choosing $epsilon_2 / 2 > epsilon_1$). Notice that *in isolation*, $C_2$ has diameter $epsilon_2 / 2 < epsilon_2$ and cardinality $m$, so every point of $C_2$ is $(epsilon_2, m)$-core and forms a valid cluster of size $m$. However, at VDBSCAN stage 1 (radius $epsilon_1$), $C_1$ is a core clique and $d(0, y_1) = epsilon_1$, so $y_1$ is directly density-reachable from $C_1$ and is absorbed as a *border point* of $C_1$ ($C_1^+ = C_1 union {y_1}$). Deflation removes $y_1$, leaving $U_2 = {y_2, dots, y_m}$ of cardinality $m - 1 < m = "MinPts"$. Consequently, at stage 2 (radius $epsilon_2$), no point in $U_2$ can attain $m$ neighbors, and the entire sparse halo $C_2$ is discarded as noise.
]

#subsec("D", "Similarity-Group Equivariance and High-Dimensional Collapse")

#theorem(title: "Scale Equivariance and High-Dimensional k-Dist Flattening", origin: "Classical — Cited")[
  1. *(Similarity Equivariance)* Under any Euclidean similarity map $T(x) = alpha Q x + b$ ($alpha > 0$, orthogonal $Q in O(d)$, $b in bb(R)^d$), the sorted $k$-dist curve scales homogeneously as $s'_i = alpha s_i$. Any scale-invariant knee detector (operating on normalized ratios $s_(i+1) / s_i$ or log-differences $log s_(i+1) - log s_i$) yields $epsilon'_ell = alpha epsilon_ell$ and produces an identical cluster partition.
  2. *(High-Dimensional Contrast Collapse)* Under the Beyer-Goldstein-Ramakrishnan-Shaft concentration regime @beyer1999nearest where $"Var"(||X||_p) / (bb(E)[||X||_p])^2 -> 0$ as $d -> infinity$, the global $k$-dist dynamic range collapses in probability:
  $
    s_n / s_1 ->^p 1 quad (d -> infinity),
  $
  eliminating $k$-dist elbows and degenerating VDBSCAN to single-scale DBSCAN ($L -> 1$) or indiscriminate noise/giant-component collapse.
]
#proof[
  *(1)* Orthogonal transformations and translations preserve Euclidean distances, so $d(T(x), T(y)) = alpha d(x, y)$ for all $x, y in cal(X)$. Hence $k$-nearest neighbor rankings are strictly preserved, $d'_k (T(x)) = alpha d_k (x)$, and $s'_i = alpha s_i$. Any knee criterion invariant under positive scaling selects identical indices $i_1, dots, i_L$, giving $epsilon'_ell = alpha epsilon_ell$. For every stage $ell$, $d(T(x), T(y)) <= alpha epsilon_ell <=> d(x, y) <= epsilon_ell$, so every DBSCAN neighborhood graph and deflation step is isomorphic to the untransformed run.

  *(2)* By Theorem 1 of @beyer1999nearest, when the relative variance of distances vanishes as $d -> infinity$, the ratio of farthest to nearest pairwise distance converges in probability to $1$: $D_"max" / D_"min" ->^p 1$, where $D_"min" = min_(i != j) d(x_i, x_j)$ and $D_"max" = max_(i != j) d(x_i, x_j)$. Since every $k$-distance satisfies $D_"min" <= s_1 <= s_n <= D_"max"$, squeezing gives
  $
    1 <= s_n / s_1 <= D_"max" / D_"min" ->^p 1.
  $
  Consequently, all normalized slope jumps satisfy $(s_(i+1) - s_i) / s_i <= s_n / s_1 - 1 ->^p 0$, erasing the multi-density elbow structure required to separate $epsilon_1 < dots < epsilon_L$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Let $L$ be the number of detected density strata. Computing exact pairwise distances takes $cal(O)(n^2 d)$ time and $cal(O)(n^2)$ space; extracting $k$-distances takes $cal(O)(n^2)$ via linear selection (or $cal(O)(n log n)$ with an R\*-tree or $k$-d tree in low $d$). Sorting the $k$-dist profile and running curvature knee detection take $cal(O)(n log n)$ and $cal(O)(n)$ time, respectively.

In the peeling loop, naive re-indexing executes $L$ independent DBSCAN calls costing $sum_(ell=1)^L cal(O)(|U_ell|^2 d)$ (worst-case $cal(O)(L n^2 d)$). However, precomputing the sorted neighbor adjacency lists (or retaining the $n times n$ distance matrix once) allows each stage $ell$ to traverse only active edges with length $<= epsilon_ell$ in $cal(O)(|E_(epsilon_ell) (U_ell)|)$ time; since each clustered point is processed and removed once, amortized peeling overhead is $cal(O)(n^2 + L n)$ in general metric spaces and $cal(O)(L n log n)$ in low-dimensional indexed spaces.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$cal(O)(n^2 d + L n)$ with cached distances; $cal(O)(L n log n)$ with low-$d$ spatial index.],
    [Space Complexity], [$cal(O)(n^2)$ with distance cache, or $cal(O)(n (d + k))$ with spatial index.],
    [Inductive Bias], [Piecewise-constant multi-density level sets with arbitrary non-convex shapes and noise rejection.],
    [Hyperparameters], [Neighbor order $k$ ($"MinPts" = k + 1$), knee detector threshold $tau$ (or stratum count $L$).],
    [Invariances], [Exact Euclidean similarity equivariance (translation, rotation, uniform scaling); permutation-invariant on core points.],
    [Main Fragility], [Sensitive to $k$-dist knee detection on smooth gradients and halo erosion at touching strata (Prop. 3).],
  ),
  caption: [Operational and complexity profile of VDBSCAN.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Urban-Suburban Geospatial Point-of-Interest Clustering:* Spatial datasets combining dense downtown retail corridors with sparse suburban shopping nodes violate DBSCAN's single-$epsilon$ assumption @liu2007vdbscan. VDBSCAN peels tight metropolitan blocks at $epsilon_1$ before expanding to suburban spacings at $epsilon_2, epsilon_3$, recovering both hierarchies cleanly.
2. *Astronomical Star-Cluster and Galaxy-Group Identification:* Deep-field sky surveys contain compact globular clusters alongside diffuse stellar associations and filamentary superclusters spanning orders of magnitude in number density. Multi-scale $k$-dist plateaus directly encode these astrophysical density scales while preserving non-spherical morphology.
3. *Single-Cell Cytometry and Rare-Population Profiling:* Flow and mass cytometry panels exhibit abundant major immune lineages alongside sparse rare progenitor populations in low-to-moderate marker subspaces. Peeling abundant dense populations first prevents large radii from fusing major T-cell subsets while still detecting sparse rare clusters.
4. *Seismic Micro-Earthquake Fault Mapping:* Aftershock swarms exhibit extreme spatio-temporal density contrast between active mainshock rupture zones and quieter secondary fault segments. Sequential peeling isolates primary fracture cores without discarding sparse background fault lineaments.
5. *Customer Mobility & Multi-Tier Market Segmentation:* Retail footfall and transaction-behavior spaces frequently comprise tight high-frequency loyalty cohorts alongside dispersed occasional shoppers; VDBSCAN separates compact core segments prior to grouping broader peripheral cohorts.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Four concrete failure regimes limit VDBSCAN in practice: (i) *smooth continuous density gradients*, where the $k$-dist plot rises smoothly without sharp elbows, causing automated knee detectors to pick arbitrary thresholds $epsilon_ell$ that over-split unimodal clusters; (ii) *touching core-halo topologies* (Proposition 3), where earlier passes erode border fringes of adjacent sparser strata; (iii) *high-dimensional plateau collapse* (Theorem 4), where $s_n / s_1 -> 1$ erases density contrast; and (iv) *uniform-density overkill*, where forcing $L > 1$ on homoscedastic data degrades cluster cohesion relative to a single tuned $epsilon$.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Multi-Density 2D/3D Spatial Data], [Strong: sharp $k$-dist plateaus yield clean multi-scale recovery.], [Verify gap between touching strata (Prop. 3).],
    [Smoothly Varying Density Gradients], [Weak: lacks discrete $k$-dist knees; prone to over-splitting.], [Use HDBSCAN @campello2013hdbscan or OPTICS @ankerst1999optics.],
    [Uniform-Density Non-Convex Manifolds], [Moderate: extra strata can fragment boundary arcs.], [Prefer standard DBSCAN @ester1996dbscan or MR-DBSCAN (Alg. 49).],
    [High-Dimensional Feature Spaces ($d > 15$)], [Poor: $k$-dist contrast collapses (Theorem 4).], [Use SNN-DPC @liu2018snndpc or RNN-DBSCAN @bryant2018rnndbscan.],
    [Nested Core-Within-Halo Structures], [Fragile: Stage 1 border absorption erodes sparse halos.], [Use HDBSCAN @campello2015hierarchical or SNN @ertoz2003snn.],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for VDBSCAN.],
)

#sec("V", "Structural Relations and Empirical Behavior")

#proposition(title: "Structural Relation Among VDBSCAN, OPTICS, and HDBSCAN", origin: "Original Derivation")[
  Let $cal(T)_"OPT"$ denote the reachability-distance plot of OPTICS @ankerst1999optics and $cal(T)_"HDB"$ the mutual-reachability minimum spanning tree of HDBSCAN @campello2013hdbscan for fixed $"MinPts" = k + 1$. Then:
  1. For $L = 1$, VDBSCAN coincides identically with a single horizontal cut at height $epsilon_1$ of both $cal(T)_"OPT"$ and the core-distance-unweighted single-linkage hierarchy (standard DBSCAN).
  2. For $L >= 2$, VDBSCAN differs structurally from taking $L$ horizontal cuts of $cal(T)_"OPT"$: horizontal cuts of $cal(T)_"OPT"$ form a *nested* laminar family ($C subset C'$ or $C inter C' = emptyset$), whereas VDBSCAN's deflation operator $U_(ell+1) = U_ell without (union.big cal(C)_ell)$ produces *non-nested* disjoint clusters by excising dense subtrees prior to coarsening.
]
#proof[
  Claim (1) is immediate because $U_1 = cal(X)$ runs standard DBSCAN$(cal(X), epsilon_1, k + 1)$, which by Theorem 1 of @ankerst1999optics is equivalent to extracting contiguous valleys below $epsilon_1$ in the OPTICS reachability plot. For Claim (2), any two horizontal thresholds $epsilon_1 < epsilon_2$ on a single-linkage or OPTICS dendrogram yield equivalence relations satisfying $x ~_(epsilon_1) y => x ~_(epsilon_2) y$, so every $epsilon_1$-cluster is strictly contained inside an $epsilon_2$-cluster. By contrast, in VDBSCAN, every point clustered at $epsilon_1$ is removed from $U_2$, so $C_1 in cal(C)_1$ and $C_2 in cal(C)_2$ are disjoint even when $"dist"(C_1, C_2) < epsilon_2$ (as demonstrated constructively in Corollary 2).
]

This non-nested excision explains VDBSCAN's distinct position relative to OPTICS @ankerst1999optics, HDBSCAN @campello2013hdbscan @campello2015hierarchical, SNN @ertoz2003snn, and RNN-DBSCAN @bryant2018rnndbscan: VDBSCAN avoids building a full reachability dendrogram or mutual-reachability MST by executing a lightweight discrete schedule of deflated DBSCAN passes, trading full hierarchical stability for speed and simplicity.

*Empirical Benchmark Analysis:* On the ClusterOverview four-dataset benchmark (Silhouette @rousseeuw1987silhouettes and Adjusted Rand Index @hubert1985comparing), VDBSCAN obtains:
- *Iris:* Silhouette $= 0.400$, $"ARI" = 0.432$;
- *Customers:* Silhouette $= 0.589$;
- *Multi-Shape:* Silhouette $= 0.614$ ($"ARI" =$ N/A);
- *High-D:* Silhouette $= 0.318$ ($"ARI" =$ N/A);
- *Mean Runtime:* $4.156$ ms.

These results align closely with the analytical trade-offs established in Section II. First, VDBSCAN runs more than twice as fast as MR-DBSCAN ($4.156$ ms vs. $9.373$ ms) because it operates in memory with zero distributed partition-merge overhead and shrinks $|U_ell|$ after each stratum. Second, on *Customers* ($0.589$) and *Multi-Shape* ($0.614$), VDBSCAN achieves solid non-convex and varied-density separation, though it trails MR-DBSCAN's exact single-scale DBSCAN ($0.663$ and $0.806$) and DPC's *Multi-Shape* score ($0.718$). That gap reflects Proposition 3 and the smooth-gradient failure mode: on moderately uniform or smoothly graded synthetic shapes, automatic $k$-dist knee detection extracts secondary radii that erode cluster halos or over-fragment smooth boundaries, also depressing *Iris* ($"Sil" = 0.400$, $"ARI" = 0.432$) where Versicolor and Virginica lack a crisp $k$-dist plateau jump. Third, *High-D* drops to $0.318$, directly reflecting Theorem 4's $k$-dist contrast collapse ($s_n / s_1 -> 1$) in high dimensions.
