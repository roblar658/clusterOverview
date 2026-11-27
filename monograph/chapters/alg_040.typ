#import "../template/ieee.typ": *

#alg_entry(
  id: 40,
  name: "ST-DBSCAN (Spatiotemporal DBSCAN)",
  category: "Cat. II: Density & Mode-Seeking",
  year: "2007",
  paradigm: "Product-Space Density-Reachability",
  complexity: $O(n log n)$ + $O(n W_"max" d_S)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Spatiotemporal datasets record events or sensor measurements indexed jointly over geographic or physical space and time, often accompanied by non-spatial physical attributes (such as sea-surface temperature, pollutant concentration, or seismic magnitude). Standard density-based clustering (DBSCAN, Alg. 31 @ester1996dbscan) collapses heterogeneous coordinates into a single isotropic radius $epsilon$, conflating spatial meters with temporal seconds and failing whenever adjacent clusters drift gradually across non-spatial values. ST-DBSCAN, introduced by Birant and Kut @birant2007stdbscan, resolves these incompatibilities by decoupling spatial, temporal, and non-spatial proximity thresholds.

Formally, let $cal(X) = {p_1, dots, p_n} subset.eq cal(S) times cal(T) times cal(A)$ be a finite spatiotemporal dataset where each observation $p_i = (s_i, t_i, a_i)$ comprises: (i) a spatial location $s_i in (cal(S), d_S)$ (typically $bb(R)^(d_S)$ with Euclidean or geodesic distance), (ii) a timestamp $t_i in (cal(T), d_T)$ (typically $bb(R)$ with $d_T (t_i, t_j) = |t_i - t_j|$ or circular distance on $bb(S)^1$ for seasonal periodicity), and (iii) an optional non-spatial attribute vector $a_i in (cal(A), d_A)$ (commonly $a_i in bb(R)$ with $d_A (u, v) = |u - v|$). Given spatial radius $epsilon_1 > 0$, temporal radius $epsilon_2 > 0$, minimum cardinality $"MinPts" in bb(N)_(>= 2)$, and optional attribute-homogeneity tolerance $Delta epsilon in (0, infinity]$, define the *spatiotemporal cylinder neighborhood* of $p_i$ as
$ N_(epsilon_1, epsilon_2) (p_i) = { p_j in cal(X) : d_S (s_i, s_j) <= epsilon_1 space and space d_T (t_i, t_j) <= epsilon_2 }. $
When non-spatial attributes $a_i$ are absent or unconstrained, one sets $Delta epsilon = infinity$ (the *pure spatiotemporal* regime); when active ($Delta epsilon < infinity$), a candidate neighbor $p_j$ is admitted into a growing cluster $C$ only if its attribute deviation from the running cluster mean $overline(a)_C = |C|^(-1) sum_(p in C) a_p$ satisfies $d_A (a_j, overline(a)_C) <= Delta epsilon$.

#definition(title: "ST-Core Points, Reachability, and Clusters")[
  Fix parameters $(epsilon_1, epsilon_2, "MinPts")$ and attribute bound $Delta epsilon = infinity$. A point $p in cal(X)$ is an *ST-core point* ($p in cal(X)_"core"$) if $|N_(epsilon_1, epsilon_2) (p)| >= "MinPts"$. Point $q$ is *directly ST-density-reachable* from $p$ ($p ->_"ST" q$) if $p in cal(X)_"core"$ and $q in N_(epsilon_1, epsilon_2) (p)$. Point $q$ is *ST-density-reachable* from $p$ ($p ->_"ST"^* q$) if there exists a finite chain $p = r_0, r_1, dots, r_m = q$ with $r_(ell-1) ->_"ST" r_ell$ for all $ell = 1, dots, m$. Points $p, q$ are *ST-density-connected* if there exists $o in cal(X)_"core"$ such that $o ->_"ST"^* p$ and $o ->_"ST"^* q$. A non-empty subset $C subset.eq cal(X)$ is an *ST-cluster* if it satisfies maximality ($p in C, space p ->_"ST"^* q => q in C$) and pairwise ST-density-connectivity; unassigned points form the *noise set* $cal(N) = cal(X) backslash union.big_c C_c$.
]

To handle adjacent regions of differing density, Birant and Kut @birant2007stdbscan additionally define a cluster-level density factor and border-assignment rule executed via breadth-first expansion (Algorithm 1).

#algbox(
  title: "ST-DBSCAN Clustering Algorithm",
  inputs: [Spatiotemporal set $cal(X) = {(s_i, t_i, a_i)}_(i=1)^n$, spatial radius $epsilon_1 > 0$, temporal radius $epsilon_2 > 0$, integer $"MinPts" >= 2$, attribute tolerance $Delta epsilon in (0, infinity]$.],
  outputs: [Cluster labels $L in {-1, 1, dots, K}^n$ ($-1$ denotes noise).],
)[
  1. Initialize $L_i <- "UNVISITED"$ for all $i in {1, dots, n}$, cluster counter $c <- 0$.
  2. *For* $i = 1, dots, n$ *do*:
  3. #h(1em) *If* $L_i != "UNVISITED"$, *continue*.
  4. #h(1em) Retrieve cylinder neighborhood $N_i <- N_(epsilon_1, epsilon_2) (p_i)$.
  5. #h(1em) *If* $|N_i| < "MinPts"$, set $L_i <- -1$ (provisional noise) and *continue*.
  6. #h(1em) Increment $c <- c + 1$; set $C_c <- {i}$, $L_i <- c$, $overline(a)_c <- a_i$, and push $N_i backslash {i}$ onto FIFO queue $cal(Q)$.
  7. #h(1em) *While* $cal(Q) != emptyset$ *do*:
  8. #h(2em) Pop $u$ from $cal(Q)$.
  9. #h(2em) *If* $L_u in {"UNVISITED", -1}$ *and* $d_A (a_u, overline(a)_c) <= Delta epsilon$:
  10. #h(3em) Assign $L_u <- c$, update $C_c <- C_c union {u}$ and running mean $overline(a)_c$.
  11. #h(3em) Retrieve $N_u <- N_(epsilon_1, epsilon_2) (p_u)$; *if* $|N_u| >= "MinPts"$, enqueue every $v in N_u$ with $L_v in {"UNVISITED", -1}$ satisfying $d_A (a_v, overline(a)_c) <= Delta epsilon$.
  12. *Return* label vector $L$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five formal properties governing ST-DBSCAN: (1) isometric embedding of the cylinder predicate into a product $L_infinity$-metric space; (2) exact core-component order invariance when $Delta epsilon = infinity$ contrasted with a constructive counterexample proving core order-dependence when $Delta epsilon < infinity$; (3) decoupled space-time scaling invariance and a sharp kinematic tracking bound; (4) statistical consistency for spatiotemporal density level sets; and (5) exact chronological sort-sweep complexity alongside adversarial breakdown bounds.

#subsec("A", "Product-Metric Equivalence and Order-Dependence Dichotomy")

#lemma(title: "Product-Norm Metric Embedding and Cylinder Geometry", origin: "Original Proof")[
  Let $(cal(S), d_S)$ and $(cal(T), d_T)$ be metric spaces and $epsilon_1, epsilon_2 > 0$. Define the rescaled Chebyshev product function on $cal(S) times cal(T)$:
  $ d_(epsilon_1, epsilon_2) (p, q) := max( (d_S (s_p, s_q)) / epsilon_1, space (d_T (t_p, t_q)) / epsilon_2 ). $
  Then: (1) $d_(epsilon_1, epsilon_2)$ is a genuine metric on $cal(S) times cal(T)$; (2) for $Delta epsilon = infinity$, ST-DBSCAN$(epsilon_1, epsilon_2, "MinPts")$ is combinatorially identical to standard DBSCAN$(epsilon = 1, "MinPts")$ on $(cal(X), d_(epsilon_1, epsilon_2))$; and (3) for $cal(S) = bb(R)^(d_S)$ ($d_S >= 1$) and $cal(T) = bb(R)$, the cylinder neighborhood $N_(epsilon_1, epsilon_2) (p)$ has $(d_S + 1)$-dimensional Lebesgue volume
  $ V(epsilon_1, epsilon_2) = 2 omega_(d_S) epsilon_1^(d_S) epsilon_2, quad omega_(d_S) = pi^(d_S / 2) / Gamma(d_S / 2 + 1), $
  whose Euclidean diameter $2 sqrt(epsilon_1^2 + epsilon_2^2)$ exceeds the diameter of any inscribed isotropic ball by a factor $sqrt(1 + (max(epsilon_1, epsilon_2) / min(epsilon_1, epsilon_2))^2) >= sqrt(2)$.
]

#proof[
  (1) Non-negativity and symmetry follow directly from $d_S, d_T >= 0$ and $epsilon_1, epsilon_2 > 0$. Moreover, $d_(epsilon_1, epsilon_2) (p, q) = 0$ iff $d_S (s_p, s_q) = 0$ and $d_T (t_p, t_q) = 0$, i.e., $(s_p, t_p) = (s_q, t_q)$. For any three points $p, q, r in cal(S) times cal(T)$, the triangle inequalities for $d_S$ and $d_T$ yield
  $ d_(epsilon_1, epsilon_2) (p, r) &= max( (d_S (s_p, s_r)) / epsilon_1, space (d_T (t_p, t_r)) / epsilon_2 ) \
    &<= max( (d_S (s_p, s_q) + d_S (s_q, s_r)) / epsilon_1, space (d_T (t_p, t_q) + d_T (t_q, t_r)) / epsilon_2 ) \
    &<= d_(epsilon_1, epsilon_2) (p, q) + d_(epsilon_1, epsilon_2) (q, r), $
  using $max(a + b, c + d) <= max(a, c) + max(b, d)$.

  (2) By definition of the maximum, $d_(epsilon_1, epsilon_2) (p, q) <= 1$ holds if and only if $d_S (s_p, s_q) <= epsilon_1$ and $d_T (t_p, t_q) <= epsilon_2$. Hence $N_(epsilon_1, epsilon_2) (p)$ coincides identically with the unit closed ball $B_(d_(epsilon_1, epsilon_2)) (p, 1)$, making core sets and reachability chains identical to DBSCAN$(1, "MinPts")$ under $d_(epsilon_1, epsilon_2)$.

  (3) The cylinder $B_(d_S) (s, epsilon_1) times [t - epsilon_2, t + epsilon_2]$ is a Cartesian product of a $d_S$-ball of volume $omega_(d_S) epsilon_1^(d_S)$ and a temporal interval of length $2 epsilon_2$, giving $V(epsilon_1, epsilon_2) = 2 omega_(d_S) epsilon_1^(d_S) epsilon_2$ by Fubini's theorem. Its maximal chord joins opposite rim points $(s + epsilon_1 u, t + epsilon_2)$ and $(s - epsilon_1 u, t - epsilon_2)$ (unit $u in bb(R)^(d_S)$) with length $2 sqrt(epsilon_1^2 + epsilon_2^2)$, whereas the largest inscribed Euclidean ball has diameter $2 min(epsilon_1, epsilon_2)$, yielding ratio $sqrt(1 + (max(epsilon_1, epsilon_2) / min(epsilon_1, epsilon_2))^2) >= sqrt(2)$.
]

#theorem(title: "Core Order Invariance vs. Attribute-Guard Order Failure", origin: "Original Counterexample")[
  Let $G_"core" = (cal(X)_"core", E_"core")$ be the undirected graph on ST-core points with $(p, q) in E_"core" <=> d_(epsilon_1, epsilon_2) (p, q) <= 1$.
  1. When $Delta epsilon = infinity$, the restriction of Algorithm 1's output to $cal(X)_"core"$ is *strictly invariant* under all $n!$ input permutations and equals the connected-component partition of $G_"core"$; only non-core border points adjacent to $>= 2$ core components depend on visitation order.
  2. When $Delta epsilon < infinity$, *core-point clusterings are no longer permutation-invariant*: there exists a 5-point dataset in which every point is an ST-core point, yet two orderings of the identical dataset produce different numbers of clusters (1 vs. 2).
]

#proof[
  (1) By Lemma 1(1), $d_(epsilon_1, epsilon_2)$ is symmetric, so the relation $p ->_"ST" q$ restricted to $cal(X)_"core"$ is reflexive and symmetric; its transitive closure $->_"ST"^*$ is therefore an equivalence relation whose equivalence classes are precisely the connected components of $G_"core"$. When Algorithm 1 seeds a new cluster $C_c$ at any unvisited $p in cal(X)_"core"$, BFS enqueues every core neighbor iteratively until all vertices in $p$'s connected component of $G_"core"$ receive label $c$, and no core point from another component is ever reached. Since a border point $b in.not cal(X)_"core"$ cannot enqueue neighbors, only border points linked to multiple components can inherit whichever core component visits them first.

  (2) Construct five points in $bb(R) times bb(R) times bb(R)$ with identical timestamp $t_i = 0$, parameters $epsilon_1 = 1$, $epsilon_2 = 1$, $"MinPts" = 2$, and attribute tolerance $Delta epsilon = 1.55$:
  $ &p_1 = (0, 0, 0.0), quad p_2 = (0, 0, 0.0), quad p_3 = (1, 0, 1.0), \
    &p_4 = (2, 0, 2.0), quad p_5 = (2, 0, 2.0). $
  Check core degrees: $|N(p_1)| = |N(p_2)| = |{p_1, p_2, p_3}| = 3 >= 2$, $|N(p_3)| = 5 >= 2$, and $|N(p_4)| = |N(p_5)| = |{p_3, p_4, p_5}| = 3 >= 2$. Thus all five points belong to $cal(X)_"core"$.
  - *Ordering A* $(p_1, p_2, p_3, p_4, p_5)$: Cluster 1 starts at $p_1$ ($overline(a)_1 = 0$), absorbs $p_2$ ($|0 - 0| = 0 <= 1.55$, $overline(a)_1 = 0$), and absorbs $p_3$ ($|1.0 - 0| = 1.0 <= 1.55$, $overline(a)_1 = 1/3$). Expanding $N(p_3)$ tests $p_4$ and $p_5$, both with attribute $2.0$:
    $ |2.0 - overline(a)_1| = |2.0 - 1 / 3| = 5 / 3 approx 1.667 > 1.55. $
    Both are rejected, terminating Cluster 1 as ${p_1, p_2, p_3}$. Subsequently $p_4$ seeds Cluster 2 $= {p_4, p_5}$. Result: *two* core clusters.
  - *Ordering B* $(p_3, p_1, p_4, p_2, p_5)$: Cluster 1 starts at $p_3$ ($overline(a)_1 = 1.0$) and scans $N(p_3)$ in list order $(p_1, p_4, p_2, p_5)$:
    - $p_1$: $|0.0 - 1.0| = 1.0 <= 1.55 =>$ added, $overline(a)_1 = 0.50$;
    - $p_4$: $|2.0 - 0.50| = 1.50 <= 1.55 =>$ added, $overline(a)_1 = 1.00$;
    - $p_2$: $|0.0 - 1.00| = 1.00 <= 1.55 =>$ added, $overline(a)_1 = 0.75$;
    - $p_5$: $|2.0 - 0.75| = 1.25 <= 1.55 =>$ added, $overline(a)_1 = 1.00$.
    Result: *one* merged cluster ${p_1, p_2, p_3, p_4, p_5}$.
  Because the running mean $overline(a)_c$ introduces path-dependent hysteresis into edge acceptance, finite $Delta epsilon$ breaks core-graph permutation invariance.
]

#subsec("B", "Decoupled Scaling Invariance and Kinematic CFL Bound")

#theorem(title: "Group Equivariance and Kinematic Trajectory Tracking Bound", origin: "Original Proof")[
  Consider the spatiotemporal transformation group $T(s, t) = (alpha Q s + b, space beta t + tau_0)$ with $Q in O(d_S)$, $b in bb(R)^(d_S)$, $tau_0 in bb(R)$, and independent positive scales $alpha, beta > 0$.
  1. ST-DBSCAN with $Delta epsilon = infinity$ is *strictly equivariant* under $T$ via $(epsilon_1', epsilon_2') = (alpha epsilon_1, beta epsilon_2)$, and invariant under rigid spatial isometries and time translations ($alpha = beta = 1$).
  2. Let a coherent moving group with spatial velocity $v in bb(R)^(d_S)$ be sampled uniformly every $delta t > 0$ along $s(t_m) = s_0 + v (m delta t)$ with at least $"MinPts"$ simultaneous observations per epoch. Consecutive epochs remain directly ST-density-reachable if and only if the *Courant–Friedrichs–Lewy (CFL) kinematic condition* holds:
     $ delta t <= epsilon_2 quad and quad norm(v)_2 <= epsilon_1 / (delta t). $
]

#proof[
  (1) Because $Q^T Q = I_(d_S)$, we have $norm((alpha Q s_i + b) - (alpha Q s_j + b))_2 = alpha norm(s_i - s_j)_2$, and $|(beta t_i + tau_0) - (beta t_j + tau_0)| = beta |t_i - t_j|$. Dividing by $(alpha epsilon_1, beta epsilon_2)$ cancels $alpha$ and $beta$ identically, preserving $d_(epsilon_1, epsilon_2) (p_i, p_j)$ for all pairs $(i, j)$ and therefore preserving the entire neighborhood structure.

  (2) Within epoch $m$, all observations at $(s(t_m), t_m)$ have zero mutual distance and form a core clique of size $>= "MinPts"$. For epoch $m$ and epoch $m+1$ to connect directly, their separation must satisfy both $d_T (t_m, t_(m+1)) = delta t <= epsilon_2$ and $d_S (s(t_m), s(t_(m+1))) = norm(v)_2 delta t <= epsilon_1$. Rearranging yields $delta t <= epsilon_2$ and $norm(v)_2 <= epsilon_1 / delta t$. If either inequality fails, consecutive epochs have disjoint neighborhoods and the moving cluster fractures into disconnected temporal slices. Conversely, in isotropic joint Euclidean DBSCAN ($norm(s_i - s_j)_2^2 + (t_i - t_j)^2 <= epsilon^2$), rescaling temporal units by $beta -> 0^+$ collapses the time axis and merges spatially colocated events separated by arbitrary time intervals.
]

#subsec("C", "Statistical Consistency of Spatiotemporal Level Sets")

#theorem(title: "Product-Space Density Level-Set Consistency", origin: "Classical — Cited")[
  Let $(S_i, T_i)_(i=1)^n$ be i.i.d. draws from a Lipschitz-continuous joint density $f(s, t)$ on $bb(R)^(d_S) times bb(R)$ whose upper level set $L_lambda = {(s, t) : f(s, t) >= lambda}$ has finitely many compact, positively separated connected components. Define the empirical cylinder density $hat(f)_n (s, t) = |N_(epsilon_1, epsilon_2) (s, t)| / (n V(epsilon_1, epsilon_2))$ with $V(epsilon_1, epsilon_2) = 2 omega_(d_S) epsilon_1^(d_S) epsilon_2$. If $epsilon_1, epsilon_2 -> 0$ such that $n epsilon_1^(d_S) epsilon_2 / log n -> infinity$ and $"MinPts" = ceil(lambda n V(epsilon_1, epsilon_2))$, then as $n -> infinity$, the connected components of the ST-DBSCAN core graph $G_"core"$ converge consistently to the connected components of $L_lambda$ @hartigan1981consistency @chaudhuri2014consistent @schubert2017dbscanrevisited.
]

#proof(label: "Proof Synthesis (Cited)")[
  By Lemma 1(3), the family of spatiotemporal cylinders $cal(B) = {B_(d_S) (s, epsilon_1) times [t - epsilon_2, t + epsilon_2]}$ is the intersection of Euclidean balls in $bb(R)^(d_S)$ (VC dimension $d_S + 2$) and intervals in $bb(R)$ (VC dimension $2$), hence $cal(B)$ has finite VC dimension $v = O(d_S)$. Standard uniform relative-deviation inequalities for VC classes @chaudhuri2014consistent yield with probability at least $1 - delta$:
  $ sup_(s, t) |hat(f)_n (s, t) - bb(E)[hat(f)_n (s, t)]| = O( sqrt( (d_S log n + log(1 / delta)) / (n V(epsilon_1, epsilon_2)) ) ) -> 0. $
  Meanwhile, Lipschitz continuity $|f(s', t') - f(s, t)| <= L (norm(s' - s)_2 + |t' - t|)$ bounds the deterministic smoothing bias $|bb(E)[hat(f)_n (s, t)] - f(s, t)| <= L (epsilon_1 + epsilon_2) -> 0$. Consequently, $hat(f)_n -> f$ uniformly in probability. Because $epsilon_1, epsilon_2$ shrink below the inter-component separation of $L_lambda$ while remaining larger than the maximal intra-component nearest-neighbor spacing inside $L_(lambda + eta)$, every true component of $L_lambda$ forms a single connected component in $G_"core"$ while distinct components remain disconnected @hartigan1981consistency @chaudhuri2014consistent.
]

#subsec("D", "Chronological Sweep Complexity and Adversarial Robustness")

#proposition(title: "Sort-Sweep Complexity and Bridge Breakdown Bounds", origin: "Original Derivation")[
  Let $W_"max" = max_(1 <= i <= n) |{j in {1, dots, n} : |t_i - t_j| <= epsilon_2}|$ be the peak temporal-window occupancy.
  1. Sorting $cal(X)$ by timestamp $t_i$ and restricting spatial distance queries to the sliding temporal window $[t_i - epsilon_2, t_i + epsilon_2]$ executes exact ST-DBSCAN in $O(n log n + n W_"max" d_S)$ time and $O(n d_S)$ space, without requiring high-dimensional spatial trees.
  2. Any adversarial injection of $m < "MinPts"$ isolated spatiotemporal outliers is completely rejected as noise (zero corruption of existing clusters). However, two ST-clusters at identical timestamp separated by spatial gap $Delta > epsilon_1$ can be adversarially merged by injecting at most $(ceil(Delta / epsilon_1) - 1) dot "MinPts"$ bridging points.
]

#proof[
  (1) Chronological sorting requires $O(n log n)$ comparisons. Maintaining two monotonic pointers $(ell_i, r_i)$ satisfying $t_(ell_i) >= t_i - epsilon_2$ and $t_(r_i) <= t_i + epsilon_2$ advances each pointer at most $n$ times across all $i = 1, dots, n$, incurring $O(n)$ index overhead. For each $i$, at most $W_"max"$ spatial distances in $bb(R)^(d_S)$ are evaluated at $O(d_S)$ cost each, totaling $O(n log n + n W_"max" d_S)$ time. BFS visits each discovered edge at most twice, bounded by the same neighbor enumeration cost, and stores only coordinates and $O(n)$ labels/queue entries.

  (2) If $m < "MinPts"$ corrupted points lie at spatial distance $> epsilon_1$ or temporal distance $> epsilon_2$ from all inliers, each corrupted point has neighborhood size at most $m < "MinPts"$, so none becomes an ST-core point and none alters any inlier neighborhood. Conversely, let $p in C_1, q in C_2$ be core points with $t_p = t_q$ and $norm(s_p - s_q)_2 = Delta > epsilon_1$. Place $r = ceil(Delta / epsilon_1) - 1$ intermediate stations uniformly along the line segment $[s_p, s_q]$ at spacing $Delta / (r + 1) <= epsilon_1$, each populated with $"MinPts"$ coincident points at timestamp $t_p$. Every station forms an ST-core clique directly reachable from its adjacent stations, creating a core chain connecting $p$ to $q$ with $r dot "MinPts"$ adversarial points.
]

#sec("III", "Complexity, Geometry and Operational Profile")

When temporal data spans long horizons relative to $epsilon_2$, chronological sort-and-sweep (Proposition 1) outperforms generic $O(n^2)$ pairwise scans dramatically: for uniform event arrival rate $lambda_t$, $W_"max" = O(lambda_t epsilon_2)$, yielding linear-logarithmic $O(n log n)$ scaling whenever $lambda_t epsilon_2 = O(log n)$. When spatial regions are large relative to $epsilon_1$ within dense time slices ($W_"max" = Omega(n)$), Birant and Kut @birant2007stdbscan accelerate range queries using a spatial $R^*$-tree or $(d_S + 1)$-dimensional spatiotemporal R-tree, achieving $O(n log n)$ average complexity for low $d_S <= 3$ while degrading toward $O(n^2 d_S)$ in high spatial dimensions due to bounding-box overlap. Because the sliding temporal window touches contiguous memory in chronologically sorted arrays, sort-sweep exhibits superior L1/L2 cache locality compared to pointer-heavy tree traversals.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for ST-DBSCAN],
    [Time / Space Complexity], [$O(n log n + n W_"max" d_S)$ via time-sweep; $O(n log n)$ with low-$d_S$ R-tree; $O(n^2 d_S)$ worst case; $O(n d_S)$ memory.],
    [Geometric Inductive Bias], [Arbitrary-shape connected tubes/trajectories in product space $cal(S) times cal(T)$ via cylinder $B_(epsilon_1) times [-epsilon_2, epsilon_2]$.],
    [Determinism], [Core partition permutation-invariant iff $Delta epsilon = infinity$; finite $Delta epsilon$ induces order dependence (Theorem 1).],
    [Invariances], [Invariant under $Q in O(d_S)$, spatial/time shifts; equivariant under independent $(alpha epsilon_1, beta epsilon_2)$ scaling.],
    [Outlier Robustness], [Rejects isolated noise ($< "MinPts"$); vulnerable to $O((Delta / epsilon_1) "MinPts")$ density bridges.],
    [Recommended Companions], [OPTICS (Alg. 33) or HDBSCAN (Alg. 32) for varying density; TRACLUS (Alg. 169) for line-segment trajectories.],
  ),
  caption: [Summary of theoretical guarantees and operational properties of ST-DBSCAN.],
)

#sec("IV", "Engineering and Scientific Applications")

ST-DBSCAN excels whenever observations carry physically incommensurable spatial and temporal coordinates and clusters evolve as localized space-time bursts, propagating fronts, or recurring seasonal anomalies rather than static convex clouds.

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Oceanography and Climate Anomaly Detection:* Originally validated by Birant and Kut @birant2007stdbscan on sea-surface temperature (SST), sea-surface height, and wave-height buoy records, ST-DBSCAN uses geographic distance $epsilon_1$, temporal window $epsilon_2$, and thermal homogeneity guard $Delta epsilon$ to isolate persistent upwelling plumes and marine heatwaves without smearing summer and winter records at identical coordinates.
2. *Seismology and Aftershock Sequence Clustering:* Earthquake catalogs reside naturally in $("lat", "lon", "depth", t)$ space where mainshock-aftershock swarms obey Omori temporal decay within fault-zone spatial corridors. Decoupled $(epsilon_1, epsilon_2)$ thresholds separate distinct aftershock sequences along the same tectonic plate boundary while filtering background micro-seismicity as noise ($L_i = -1$).
3. *Urban Mobility and Ride-Hailing Hotspot Dynamics:* GPS pickup/drop-off streams exhibit sharp spatiotemporal localization (e.g., stadium egresses, rush-hour transit hubs). The product cylinder $B_(epsilon_1) times [-epsilon_2, epsilon_2]$ tracks shifting congestion cores provided sampling satisfies the CFL velocity bound (Theorem 2), while circular temporal distance $d_T$ on $bb(S)^1$ captures recurring weekly demand regimes.
4. *Epidemiological Outbreak and Public-Health Surveillance:* Spatiotemporal disease incidence logs require separating genuine localized transmission clusters from sporadic background cases. Non-parametric cylinder connectivity adapts to irregular road networks or river corridors, and explicit noise labeling prevents isolated cases from distorting outbreak boundaries.
5. *Maritime AIS Tracking and Wildlife Telemetry:* Automatic Identification System (AIS) vessel pings and GPS animal collars produce dense loitering bursts (fishing grounds, foraging stops) linked by fast transit legs. By Theorem 2, setting $epsilon_1 / delta t$ below transit speed automatically isolates stationary loitering clusters while classifying high-speed transit pings as inter-cluster noise.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

ST-DBSCAN degrades systematically under four regimes: (a) *multi-scale variable density*, where a single global pair $(epsilon_1, epsilon_2)$ cannot simultaneously resolve dense urban cores and sparse rural events (addressed by HDBSCAN Alg. 32 @campello2013hdbscan, VDBSCAN Alg. 50 @liu2007vdbscan, or RNN-DBSCAN Alg. 45 @bryant2018rnndbscan); (b) *gradual attribute drift with finite $Delta epsilon$*, where order-dependent running averages (Theorem 1(2)) yield non-reproducible splits; (c) *high-velocity aliasing*, where moving clusters violating $norm(v)_2 <= epsilon_1 / delta t$ fracture into disjoint slices; and (d) *high-dimensional non-spatial spaces*, where distance concentration @beyer1999nearest collapses the contrast between core and noise radii.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Why ST-DBSCAN Fits / Requirement], [Caution / Alternative],
    [Climate / SST Fields], [Separates spatial, seasonal & value scales], [Order-sensitive if $Delta epsilon < infinity$],
    [Seismic Catalogs], [Fault-zone shapes + Omori time window], [Use OPTICS (Alg. 33) for multi-density],
    [Urban GPS Streams], [Fast sort-sweep; isolates surge hotspots], [Use DenStream (Alg. 194) online],
    [Full Trajectory Paths], [Isolates stop/loiter bursts via CFL bound], [Use TRACLUS (Alg. 169) for paths],
    [High-D Static Tables], [Cylinder collapses when time absent], [Use SNN (Alg. 36) or HDBSCAN (Alg. 32)],
  ),
  caption: [Application suitability and engineering guidelines for ST-DBSCAN.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Within Category II, ST-DBSCAN sits directly inside the generalized density-reachability hierarchy linking DBSCAN (Alg. 31 @ester1996dbscan), GDBSCAN (Alg. 39 @sander1998gdbscan), and spatiotemporal streaming methods.

#proposition(title: "Exact Reductions to GDBSCAN and Classical DBSCAN", origin: "Original Derivation")[
  Let $Delta epsilon = infinity$. Then: (1) ST-DBSCAN is a strict specialization of GDBSCAN (Alg. 39 @sander1998gdbscan) obtained by setting the binary symmetric neighborhood predicate to $"NPred"(p, q) := (d_S (s_p, s_q) <= epsilon_1) and (d_T (t_p, t_q) <= epsilon_2)$ and uniform unit weights $w(p) equiv 1$ with cardinality threshold $"MinCard" = "MinPts"$; and (2) taking the infinite temporal-window limit $epsilon_2 -> infinity$ recovers classical spatial DBSCAN$(epsilon_1, "MinPts")$ (Alg. 31) identically, whereas taking $epsilon_1 -> infinity$ reduces ST-DBSCAN to one-dimensional temporal burst segmentation.
]

#proof[
  (1) GDBSCAN requires only a reflexive, symmetric neighborhood predicate $"NPred"$ and a monotone weight aggregation function $w_"Card" (S) = sum_(p in S) w(p)$ @sander1998gdbscan. By Lemma 1, $"NPred"(p, q) <=> d_(epsilon_1, epsilon_2) (p, q) <= 1$ is reflexive and symmetric, and with $w equiv 1$, $w_"Card" (N_"NPred" (p)) = |N_(epsilon_1, epsilon_2) (p)| >= "MinPts"$ reproduces Definition 1 verbatim. (2) For any finite dataset, letting $epsilon_2 >= max_(i, j) d_T (t_i, t_j)$ makes $d_T (t_i, t_j) <= epsilon_2$ tautologically true for all pairs, reducing $N_(epsilon_1, epsilon_2) (p_i)$ to the purely spatial $epsilon_1$-ball of DBSCAN; symmetrically, $epsilon_1 >= max_(i, j) d_S (s_i, s_j)$ eliminates spatial filtering and clusters solely contiguous temporal intervals.
]

Beyond GDBSCAN, ST-DBSCAN complements trajectory-segment algorithms such as TRACLUS (Alg. 169 @lee2007traclus), which cluster line segments rather than point events, as well as online density stream engines such as DenStream (Alg. 194 @cao2006denstream) and D-Stream (Alg. 197 @chen2007dstream), which replace hard temporal windows $[-epsilon_2, epsilon_2]$ with continuous exponential time-decay weights $2^(-lambda Delta t)$.

On the ClusterOverview 4-dataset empirical benchmark (evaluating Silhouette @rousseeuw1987silhouettes and Adjusted Rand Index @hubert1985comparing with normalized temporal window $epsilon_2 = 0.45$ and $k$-distance spatial threshold $epsilon_1$), ST-DBSCAN records *Iris* $"Sil" = 0.615$ and $"ARI" = 0.547$, *Customers* $"Sil" = 0.492$, *Multi-Shape* $"Sil" = 0.630$ ($"ARI" =$ N/A), and *High-D Subspace+Noise* $"Sil" = 0.359$ ($"ARI" =$ N/A), at a mean runtime of *1.028 ms*. These numbers mirror our theoretical profile directly: strong Silhouette on *Multi-Shape* ($0.630$) and *Iris* ($0.615$) confirms non-convex density-reachability and noise filtering; the moderate *Iris* $"ARI" = 0.547$ reflects the well-known merging of the closely touching _Versicolor_ and _Virginica_ density lobes under a single global $epsilon_1$; the drop on *High-D Subspace+Noise* ($0.359$) reflects high-dimensional distance concentration @beyer1999nearest; and the $1.028$ ms runtime reflects fast temporal-window filtering consistent with Proposition 1.
