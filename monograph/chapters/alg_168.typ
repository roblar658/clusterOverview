#import "../template/ieee.typ": *

#alg_entry(
  id: 168,
  name: "DBSCAN-3D / Euclidean Cluster Extraction",
  category: "Cat. VII: Fuzzy, Biclustering & Domain",
  year: "2011",
  paradigm: "3D Spatial-Index Density & Euclidean Connected-Component Extraction",
  complexity: $O(N log N)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Unorganized 3D point clouds from spinning LiDAR, RGB-D sensors, and photogrammetry consist of millions of sparse spatial coordinates lacking mesh topology. Isolating physical objects—vehicles, pedestrians, tabletop items, or buildings—requires fast geometric segmentation handling arbitrary non-convex 3D shapes without pre-specifying object counts. The Point Cloud Library (PCL) Euclidean Cluster Extraction (ECE) algorithm of Rusu and Cousins @rusu2011pcl and its core-density generalization DBSCAN-3D (adapting Ester et al. @ester1996dbscan and Sander et al. @sander1998gdbscan to 3D spatial indices) extract density-gated Euclidean connected components via balanced 3D $k$-d trees or voxel octrees.

Let $cal(P) = \{bold(p)_1, dots, bold(p)_N\} subset bb(R)^3$ be a 3D point cloud. For cluster tolerance $epsilon > 0$, define the closed Euclidean ball $B_epsilon (bold(p)) = \{bold(q) in cal(P) : norm(bold(p) - bold(q))_2 <= epsilon\}$ and the *3D $epsilon$-proximity graph* $G_epsilon (cal(P)) = (cal(P), E_epsilon)$ where $\{bold(p), bold(q)\} in E_epsilon <==> 0 < norm(bold(p) - bold(q))_2 <= epsilon$. PCL ECE applies cluster cardinality bounds $1 <= S_"min" <= S_"max" <= N$, whereas DBSCAN-3D adds a local core density threshold $m_"pts" in bb(N)$ (`"MinPts"`).

#definition(title: "DBSCAN-3D and PCL Euclidean Cluster Extraction")[
  Fix tolerance $epsilon > 0$, core threshold $m_"pts" >= 1$, and cardinality bounds $1 <= S_"min" <= S_"max"$.
  + Point $bold(p) in cal(P)$ is a *core point* if $|B_epsilon (bold(p))| >= m_"pts"$; denote the core set by $cal(P)_"core" (epsilon, m_"pts")$.
  + Point $bold(q)$ is *directly density-reachable* from $bold(p)$ if $bold(p) in cal(P)_"core"$ and $bold(q) in B_epsilon (bold(p))$; *density-reachable* via a finite chain of directly density-reachable steps; and *density-connected* to $bold(r)$ if both are density-reachable from a shared core point $bold(o)$.
  + A *DBSCAN-3D cluster* is a maximal non-empty subset $C subset.eq cal(P)$ closed under density-reachability and pairwise density-connectivity with $S_"min" <= |C| <= S_"max"$; unclustered points form noise $cal(N)$.
  + *PCL Euclidean Cluster Extraction (ECE)* @rusu2011pcl is the exact specialization $m_"pts" = 1$: every point is core, so clusters are precisely the connected components $C$ of $G_epsilon (cal(P))$ with $S_"min" <= |C| <= S_"max"$.
]

Production LiDAR pipelines typically precede clustering by ground-plane removal (RANSAC) and cubic voxel-grid downsampling of leaf edge $v > 0$, indexing the resulting centroids in a 3D $k$-d tree $cal(T)$ supporting `radiusSearch`($bold(p), epsilon$).

#algbox(
  title: "DBSCAN-3D / PCL Euclidean Cluster Extraction",
  inputs: [3D point cloud $cal(P) subset bb(R)^3$, tolerance $epsilon > 0$, core threshold $m_"pts"$ ($= 1$ for ECE), size bounds $[S_"min", S_"max"]$, optional voxel leaf $v >= 0$.],
  outputs: [Valid 3D clusters $\{C_1, dots, C_K\}$ and noise set $cal(N)$.],
)[
  + *(Optional Voxelization)* If $v > 0$, hash points into cubic voxels of edge $v$ and replace each occupied voxel by its centroid to form $cal(P)_v$.
  + Build a balanced 3D $k$-d tree $cal(T)$ over the active cloud in $O(N log N)$ time.
  + Set $"visited"(i) <- "false"$, $"lab"(i) <- 0$ for all $i$, and $cal(C) <- emptyset$.
  + *For each* unvisited point $bold(p)_i$:
    - Mark visited and query $N_i <- cal(T)."radiusSearch"(bold(p)_i, epsilon)$.
    - If $|N_i| < m_"pts"$, mark $bold(p)_i$ provisionally noise and continue.
    - Initialize $C <- \{bold(p)_i\}$ and FIFO queue $Q <- N_i without \{bold(p)_i\}$.
    - *While* $Q != emptyset$: pop $bold(q)$; if unvisited, mark visited, compute $N_q <- cal(T)."radiusSearch"(bold(q), epsilon)$, and if $|N_q| >= m_"pts"$, append unqueued points of $N_q$ to $Q$; if unassigned, add $bold(q)$ to $C$.
    - Retain $C$ in $cal(C)$ iff $S_"min" <= |C| <= S_"max"$; else label $C$ as noise.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five rigorous results: exact EMST equivalence of ECE and bridge-chaining separation from DBSCAN-3D; exact voxel-grid topological preservation and packing bounds; spinning-LiDAR $Theta(r^(-2))$ density decay and range-adaptive tolerance; $E(3)$ invariance and border-point order sensitivity; and 3D vs. high-dimensional $k$-d tree complexity.

#subsec("A", "Single-Linkage Equivalence and Bridge-Chaining Separation")

#theorem(title: "EMST-Cut Equivalence of ECE and Bridge-Suppression Separation of DBSCAN-3D", origin: "Original Proof")[
  Let $cal(P) subset bb(R)^3$ be finite and $T^* (cal(P))$ be any Euclidean Minimum Spanning Tree (EMST) of $cal(P)$.
  + *(ECE--EMST Identity)* Ignoring $S_"max"$, PCL ECE ($m_"pts" = 1$) outputs identically the connected components of size $>= S_"min"$ formed by deleting all edges of weight $> epsilon$ from $T^* (cal(P))$ (the $epsilon$-cut of Single-Linkage @sibson1973slink).
  + *(Bridge Separation)* There exists a 3D configuration of two dense objects $A, B$ separated by $> epsilon$ linked by a 1D noise bridge that ECE merges for all $S_"min" <= |A| + |B| + L$, whereas DBSCAN-3D with $m_"pts" = 4$ recovers $A$ and $B$ strictly separated.
]
#proof[
  *(Part 1)* Let $T_epsilon$ retain all edges of $T^* (cal(P))$ of weight $<= epsilon$. Clearly $E(T_epsilon) subset.eq E_epsilon$, so connectivity in $T_epsilon$ implies connectivity in $G_epsilon (cal(P))$. Conversely, if $bold(p), bold(q)$ connect in $G_epsilon (cal(P))$ via hops $bold(u)_(t-1) - bold(u)_t$ of length $<= epsilon$, the MST cycle property forces every edge on the unique path between $bold(u)_(t-1)$ and $bold(u)_t$ in $T^* (cal(P))$ to have weight $<= norm(bold(u)_(t-1) - bold(u)_t)_2 <= epsilon$ (otherwise swapping lowers MST weight). Hence every hop stays connected in $T_epsilon$, proving identical connected components prior to the $|C| >= S_"min"$ filter.

  *(Part 2)* Fix $epsilon = 1$ and $L >= 3$. Place tight clusters $A, B$ of $M >= 5$ points inside radius-$0.2$ balls centered at $(0, 0, 0)$ and $(L + 1, 0, 0)$ (with boundary points $(0.2, 0, 0) in A$, $(L + 0.8, 0, 0) in B$), joined by bridge points $bold(b)_t = (t, 0, 0)$ for $t = 1, dots, L$. Because all consecutive distances are $<= 1$, $G_1 (cal(P))$ is one component of size $2 M + L$, so ECE merges $A union B$. Under DBSCAN-3D with $m_"pts" = 4$, every point in $A$ and $B$ has $|B_1| >= M >= 5$ (core), whereas every interior bridge point $bold(b)_t$ ($2 <= t <= L - 1$) satisfies $B_1 (bold(b)_t) = \{bold(b)_(t-1), bold(b)_t, bold(b)_(t+1)\}$, so $|B_1 (bold(b)_t)| = 3 < 4$ (non-core). Since density-reachable chains cannot traverse non-core points, DBSCAN-3D discards the bridge and isolates $A$ and $B$.
]

#subsec("B", "Voxel-Grid Metric Distortion and Exact Separation Preservation")

#theorem(title: "Exact Voxel-Grid Topological Preservation and Packing Bound", origin: "Original Derivation")[
  Let $cal(P)_v$ replace every occupied cubic voxel of edge $v > 0$ in $cal(P) subset bb(R)^3$ by its centroid.
  + *(Distortion Bound)* For all $bold(p), bold(q) in cal(P)$ with centroids $bar(bold(p)), bar(bold(q)) in cal(P)_v$,
    $ |norm(bar(bold(p)) - bar(bold(q)))_2 - norm(bold(p) - bold(q))_2| <= 2 sqrt(3) v. $
  + *(Exact Preservation Margin)* If raw objects $C_1, dots, C_K$ are each $epsilon_0$-connected and pairwise separated by minimum Euclidean gap $delta$, then whenever $delta - epsilon_0 > 4 sqrt(3) v$, any $epsilon in [epsilon_0 + 2 sqrt(3) v, thick delta - 2 sqrt(3) v)$ guarantees zero false splits and zero false merges in ECE on $cal(P)_v$.
  + *(Uniform Degree Bound)* Every Euclidean $epsilon$-ball in $cal(P)_v$ contains at most
    $ |B_epsilon (bar(bold(p)))| <= (4 pi) / 3 (epsilon / v + sqrt(3))^3 $
    voxel centroids, certifying $O(1)$ degree whenever $epsilon \/ v = O(1)$.
]
#proof[
  *(Part 1)* Because a cubic voxel of edge $v$ has Euclidean diameter $sqrt(3) v$ and centroid $bar(bold(p))$ lies in the same cube, $norm(bold(p) - bar(bold(p)))_2 <= sqrt(3) v$. Two applications of the triangle inequality yield $|norm(bar(bold(p)) - bar(bold(q)))_2 - norm(bold(p) - bold(q))_2| <= norm(bar(bold(p)) - bold(p))_2 + norm(bar(bold(q)) - bold(q))_2 <= 2 sqrt(3) v$.

  *(Part 2)* For intra-object pairs with $norm(bold(p) - bold(q))_2 <= epsilon_0$, Part 1 gives $norm(bar(bold(p)) - bar(bold(q)))_2 <= epsilon_0 + 2 sqrt(3) v <= epsilon$ (preserving all connections). For inter-object pairs with $norm(bold(p) - bold(q))_2 >= delta$, Part 1 gives $norm(bar(bold(p)) - bar(bold(q)))_2 >= delta - 2 sqrt(3) v > epsilon$ (preventing inter-object edges and shared voxels since $delta > 4 sqrt(3) v > sqrt(3) v$).

  *(Part 3)* Each $bar(bold(q)) in B_epsilon (bar(bold(p)))$ corresponds to a distinct occupied cube $V_q$ of volume $v^3$. Every point $bold(x) in V_q$ satisfies $norm(bold(x) - bar(bold(q)))_2 <= sqrt(3) v$, hence $norm(bold(x) - bar(bold(p)))_2 <= epsilon + sqrt(3) v$. Thus the pairwise-disjoint cubes $\{V_q : bar(bold(q)) in B_epsilon (bar(bold(p)))\}$ pack inside a sphere of radius $epsilon + sqrt(3) v$ centered at $bar(bold(p))$. Comparing volumes gives $|B_epsilon (bar(bold(p)))| dot.c v^3 <= (4 pi \/ 3) (epsilon + sqrt(3) v)^3$; dividing by $v^3$ yields the stated bound.
]

#subsec("C", "Spinning-LiDAR Range Decay and Adaptive Tolerance")

#theorem(title: "LiDAR Inverse-Square Density Decay and Critical Fragmentation Range", origin: "Original Derivation")[
  Consider a spinning 3D LiDAR at the origin with azimuth step $Delta theta > 0$ and elevation step $Delta phi > 0$ scanning a locally planar surface at range $r > 0$ with normal incidence angle $psi in [0, pi \/ 2)$.
  + Local surface density decays as $rho_"surf" (r) = (cos psi) \/ (r^2 Delta theta Delta phi) = Theta(r^(-2))$, while adjacent beam spacings grow as $Delta s_h (r) approx r Delta theta \/ cos psi_h$ and $Delta s_v (r) approx r Delta phi \/ cos psi_v$.
  + Under fixed $epsilon > 0$, a rigid surface fragments into disconnected scan rings once range exceeds $r_"crit" (epsilon) = (epsilon cos psi) \/ max(Delta theta, Delta phi)$.
  + Setting range-adaptive tolerance $epsilon(r) = kappa r max(Delta theta, Delta phi)$ with $kappa > 1 \/ cos psi_"max"$ preserves connectivity and yields range-invariant expected neighborhood count $bb(E)[|B_(epsilon(r)) (bold(p))|] = Theta(kappa^2)$ for all $r > 0$.
]
#proof[
  Solid angle $Delta theta Delta phi$ subtends orthogonal area $r^2 Delta theta Delta phi$ at range $r$, projecting onto a surface tilted by $psi$ to area $(r^2 Delta theta Delta phi) \/ cos psi$, giving $rho_"surf" (r) = (cos psi) \/ (r^2 Delta theta Delta phi)$ and inter-ray spacings $r Delta phi \/ cos psi_v$, $r Delta theta \/ cos psi_h$ to first order. For $r > r_"crit" (epsilon)$, the coarser spacing exceeds $epsilon$, severing inter-ring edges in $G_epsilon (cal(P))$. Conversely, under $epsilon(r) = kappa r max(Delta theta, Delta phi)$ ($kappa > 1 \/ cos psi_"max"$), both spacings stay strictly below $epsilon(r)$ and the expected return count inside a surface disk of radius $epsilon(r)$ is
  $ bb(E)[|B_(epsilon(r)) (bold(p))|] approx pi epsilon(r)^2 rho_"surf" (r) = pi kappa^2 cos psi dot.c (max(Delta theta, Delta phi)^2) / (Delta theta Delta phi), $
  eliminating range dependence $r$ identically.
]

#subsec("D", "Rigid-Motion Invariance and Border-Point Determinism")

#proposition(title: "E(3) Invariance and Order Sensitivity of Core vs. Border Assignments", origin: "Original Proof")[
  Both PCL ECE and DBSCAN-3D core classes are exactly invariant under 3D rigid motions $T(bold(p)) = R bold(p) + bold(t)$ ($R in O(3)$) and equivariant under uniform scaling $bold(p) |-> alpha bold(p)$ iff $epsilon |-> alpha epsilon$. PCL ECE ($m_"pts" = 1$) is strictly permutation-invariant, whereas DBSCAN-3D ($m_"pts" >= 2$) is permutation-invariant on $cal(P)_"core"$ but order-dependent on border points $bold(b) in B_epsilon (bold(c)_1) inter B_epsilon (bold(c)_2)$ shared by distinct core clusters @schubert2017dbscanrevisited.
]
#proof[
  Rigid isometries preserve all Euclidean norms, leaving $B_epsilon$, $G_epsilon$, and $cal(P)_"core"$ unchanged; uniform scaling multiplies norms by $alpha$. For ECE, clusters are undirected graph components of $G_epsilon (cal(P))$, independent of BFS seed order. For DBSCAN-3D ($m_"pts" >= 2$), density-connectivity forms an equivalence relation on $cal(P)_"core"$ (order-invariant @schubert2017dbscanrevisited), but a shared non-core border point $bold(b)$ is permanently claimed by whichever core cluster's BFS reaches $bold(b)$ first.
]

#subsec("E", "3D Spatial-Index Query Complexity and High-Dimensional Collapse")

#theorem(title: "3D k-d Tree Complexity and Curse-of-Dimensionality Degradation", origin: "Classical — Cited")[
  Constructing a balanced $d$-dimensional $k$-d tree takes $O(d N log N)$ time and $O(d N)$ space. In $d = 3$, DBSCAN-3D / ECE runs in $O(N log N)$ expected time when $sum_(i=1)^N |B_epsilon (bold(p)_i)| = O(N)$ (guaranteed post-voxelization by Theorem 2(3)) and $O(N^2)$ worst case @rusu2011pcl @schubert2017dbscanrevisited. In general dimension $d$, worst-case radius queries degrade to $O(d N^(1 - 1 \/ d) + |B_epsilon|)$ while pairwise distance concentration collapses fixed-$epsilon$ separation @beyer1999nearest.
]
#proof[
  Median-split construction satisfies $T(N) = 2 T(N \/ 2) + O(d N) = O(d N log N)$. Each point issues at most one `radiusSearch`; under bounded 3D local crowding, each query visits $O(log N)$ internal nodes plus $O(|B_epsilon (bold(p)_i)|)$ reported points, totaling $O(N log N)$ (or $Theta(N^2)$ if all points lie within radius $epsilon$ @schubert2017dbscanrevisited). In $bb(R)^d$, an axis-aligned hypercube bounding an $epsilon$-ball intersects $O(N^(1 - 1 \/ d))$ $k$-d cells worst case (approaching linear $O(d N)$ scan for large $d$), while $(d_"max" - d_"min") \/ d_"min"$ converges in probability to $0$ under broad i.i.d. distributions @beyer1999nearest, destroying radius contrast.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Combining voxel-grid downsampling (leaf $v$) with a 3D $k$-d tree or octree keeps DBSCAN-3D / ECE in the $O(N log N)$ regime with $O(N)$ memory, because Theorem 2(3) caps every neighborhood at $O((epsilon \/ v)^3)$. Flat grid hashing (`GriDBSCAN` @mahran2008gridbscan) or Morton-code sorting further boosts cache locality to expected $O(N)$ lookup for fixed $epsilon$. Geometrically, DBSCAN-3D / ECE isolates arbitrary non-convex 3D shells and irregular objects while filtering outliers via $m_"pts"$ and $S_"min"$, though unremoved support planes (roads, tabletops) topologically short-circuit all resting objects unless segmented beforehand.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$O(N log N)$ in 3D with bounded density ($k$-d tree / octree); $O(N^2)$ worst case],
    [Space Complexity], [$O(N)$ for 3D $k$-d tree, visited flags, BFS queue, and cluster indices],
    [Inductive Bias], [3D Euclidean connected components (ECE) or core-gated density chains],
    [Cluster Shape], [Arbitrary non-convex 3D surfaces, shells, and physical objects],
    [Invariance], [Exact $E(3)$ rigid-motion invariance; ECE strictly order-invariant (Prop. 1)],
    [Key Hyperparameters], [Tolerance $epsilon$; core threshold $m_"pts"$ ($1$ for ECE); size filter $[S_"min", S_"max"]$],
    [LiDAR Range Scaling], [Fixed $epsilon$ fragments beyond $r_"crit" prop epsilon \/ Delta phi$; use adaptive $epsilon(r) prop r$],
    [Dimensional Scaling], [Optimal in $d = 2, 3$; $k$-d tree pruning and radius contrast collapse in high $d$],
  ),
  caption: [Operational and theoretical profile of DBSCAN-3D / ECE (Alg. 168).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

+ *Autonomous Driving LiDAR Perception:* After RANSAC ground-plane removal, vehicles, cyclists, and pedestrians float as well-separated 3D point clouds. Voxelized DBSCAN-3D with adaptive $epsilon(r)$ (Theorem 3) extracts class-agnostic obstacle proposals in real time, catching out-of-distribution hazards missed by supervised detectors @rusu2011pcl.
+ *Robotic Manipulation & Industrial Bin-Picking:* On RGB-D tabletop scans, removing the supporting plane lets PCL ECE isolate graspable parts of arbitrary pose/shape while $[S_"min", S_"max"]$ discards sensor speckles and bin walls @rusu2011pcl.
+ *Airborne LiDAR Urban & Forestry Mapping:* Stripping ground returns allows DBSCAN-3D ($m_"pts" >= 4$) to reject stray returns (Theorem 1(2)) while grouping buildings, pylons, and tree crowns across kilometer tiles.
+ *3D Industrial Metrology & Additive Inspection:* Structured-light scans use ECE directly in millimeter units to isolate the primary workpiece from spatter or detached supports.
+ *3D Microscopy & Astrophysical Structure Finding:* Cryo-ET, 3D single-molecule microscopy, and cosmological Friends-of-Friends halo finders (mathematically identical to ECE) extract non-spherical complexes and dark-matter filaments in $O(N log N)$ time.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Four core failure modes arise: (1) *Support-Plane Short-Circuiting* when unremoved ground merges all resting objects; (2) *Spinning-LiDAR Fragmentation* beyond $r_"crit"$ under fixed $epsilon$ (Theorem 3); (3) *Single-Point Bridge Chaining in ECE* through rain/dust unless $m_"pts" >= 4$ (Theorem 1(2)); and (4) *High-Dimensional Collapse* where $k$-d trees approach $O(d N^2)$ and distance concentration erodes density gaps @beyer1999nearest.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Autonomous LiDAR (No Ground)], [High: non-convex 3D objects + real-time $k$-d tree speed], [Use adaptive $epsilon(r)$ (Thm. 3) to stop far splits],
    [Tabletop / Bin-Picking RGB-D], [Excellent: exact physical gap isolation in PCL @rusu2011pcl], [Remove table plane first; set $m_"pts" >= 4$ for noise],
    [Rain / Fog LiDAR Clouds], [Moderate: ECE chains through airborne droplets], [Use DBSCAN-3D ($m_"pts" >= 5$) or outlier pre-filter],
    [Multi-Density Spatial Clouds], [Low--Moderate: global $epsilon$ misses sparse clusters], [Use OPTICS @ankerst1999optics or HDBSCAN @campello2013hdbscan],
    [High-D Generic Feature Data], [Poor: $k$-d tree degrades to $O(d N^2)$; distances concentrate], [Use HDBSCAN on UMAP/PCA or subspace methods],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for DBSCAN-3D / ECE.],
)

#sec("V", "Structural Relations and Empirical Behavior")

DBSCAN-3D / ECE anchors spatial density clustering, connecting planar DBSCAN @ester1996dbscan, OPTICS @ankerst1999optics, HDBSCAN @campello2013hdbscan, GriDBSCAN @mahran2008gridbscan, and Single-Linkage @sibson1973slink.

#proposition(title: "Inclusion Hierarchy Across ECE, DBSCAN-3D, and HDBSCAN", origin: "Original Proof")[
  For any finite $cal(P) subset bb(R)^3$ and $epsilon > 0$, let $cal(C)_"DB" (epsilon, m_"pts")$ be the DBSCAN-3D core partition and $cal(C)_"SL" (epsilon)$ the $epsilon$-cut of Single-Linkage (ECE without size filtering). Then: (i) $cal(C)_"DB" (epsilon, m_"pts")$ strictly refines the restriction of $cal(C)_"SL" (epsilon)$ to $cal(P)_"core"$, with universal equality iff $m_"pts" = 1$; and (ii) replacing Euclidean weights with mutual-reachability distances $d_"mreach" (bold(p), bold(q)) = max("core"_k (bold(p)), "core"_k (bold(q)), norm(bold(p) - bold(q))_2)$ lifts DBSCAN-3D into the HDBSCAN hierarchy @campello2013hdbscan.
]
#proof[
  Every directly density-reachable core pair satisfies $norm(bold(p) - bold(q))_2 <= epsilon$, so density-connected core clusters stay inside single components of $G_epsilon (cal(P))$ ($= cal(C)_"SL" (epsilon)$ by Theorem 1(1)), with strictness shown by Theorem 1(2) for $m_"pts" >= 4$ and equality for $m_"pts" = 1$. For (ii), $d_"mreach" (bold(p), bold(q)) <= epsilon$ holds iff both endpoints are core for $m_"pts" = k$ and within distance $epsilon$, so the $epsilon$-cut of the mutual-reachability MST reproduces $cal(C)_"DB" (epsilon, k)$ identically @campello2013hdbscan.
]

#corollary(title: "Voxel-Grid Pre-Downsampling Distortion Bound for 3D ECE", origin: "Original Proof")[
  Let $cal(P)_v$ be obtained from a 3D point cloud $cal(P) subset bb(R)^3$ by mapping all points inside each non-empty cubic voxel of edge length $v > 0$ to its geometric voxel center $hat(bold(p)) = bold(c)_p$ (or to its voxel centroid with factor $2 sqrt(3) v$), and let $epsilon_0 > sqrt(3) v$. Then every pair $(bold(p), bold(q)) in cal(P)^2$ satisfies
  $ norm(bold(p) - bold(q))_2 - sqrt(3) v <= norm(hat(bold(p)) - hat(bold(q)))_2 <= norm(bold(p) - bold(q))_2 + sqrt(3) v. $
  Consequently, ECE on $cal(P)_v$ with radius $epsilon_0$ sandwiches exact unquantized ECE on $cal(P)$: every component of $G_(epsilon_0 - sqrt(3) v) (cal(P))$ remains connected in $G_(epsilon_0) (cal(P)_v)$, and every component of $G_(epsilon_0) (cal(P)_v)$ is contained in a single component of $G_(epsilon_0 + sqrt(3) v) (cal(P))$.
]
#proof[
  Inside a cubic voxel of side $v$ in $bb(R)^3$, the maximum Euclidean distance from any point $bold(p)$ in the voxel to its geometric center $hat(bold(p)) = bold(c)_p$ equals half the space diagonal, i.e., $norm(bold(p) - hat(bold(p)))_2 <= (sqrt(3) / 2) v$. Applying the reverse and forward triangle inequalities to any pair $(bold(p), bold(q)) in cal(P)^2$ yields
  $ |norm(hat(bold(p)) - hat(bold(q)))_2 - norm(bold(p) - bold(q))_2| <= norm(bold(p) - hat(bold(p)))_2 + norm(bold(q) - hat(bold(q)))_2 <= sqrt(3) v. $
  Therefore, if $norm(bold(p) - bold(q))_2 <= epsilon_0 - sqrt(3) v$, then $norm(hat(bold(p)) - hat(bold(q)))_2 <= epsilon_0$, so every edge of $G_(epsilon_0 - sqrt(3) v) (cal(P))$ persists in $G_(epsilon_0) (cal(P)_v)$; symmetrically, every edge of $G_(epsilon_0) (cal(P)_v)$ implies $norm(bold(p) - bold(q))_2 <= epsilon_0 + sqrt(3) v$, establishing the exact two-sided component-sandwich guarantee while reducing sample cardinality by $O(rho v^3)$. (For voxel centroids instead of geometric centers, replacing the circumradius $(sqrt(3)/2) v$ by the full voxel diameter $sqrt(3) v$ yields the analogous $2 sqrt(3) v$ bound.)
]

On the *ClusterOverview* four-dataset benchmark, DBSCAN-3D / ECE achieves *Iris* Silhouette @rousseeuw1987silhouettes $0.426$ and ARI @hubert1985comparing $0.442$, *Customers* Silhouette $0.327$, *Multi-Shape* Silhouette $0.596$ (ARI N/A), and *High-D* Silhouette $0.183$ (ARI N/A), at a mean runtime of $1.495$ ms. These numbers directly confirm Theorems 1 and 4: strong *Multi-Shape* recovery ($0.596$) reflects non-convex density chaining, moderate *Iris* ARI ($0.442$) stems from global-$epsilon$ bridging across overlapping _versicolor_ / _virginica_ populations, and the sharp *High-D* drop to $0.183$ illustrates distance concentration and $k$-d tree degradation @beyer1999nearest far outside its native 3D point-cloud domain.
