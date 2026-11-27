#import "../template/ieee.typ": *

#alg_entry(
  id: 167,
  name: "Watershed Clustering",
  category: "Cat. VII: Fuzzy, Biclustering & Domain",
  year: "1991",
  paradigm: "Topographic Immersion & Geodesic Influence-Zone Partitioning",
  complexity: $O(N log N)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Watershed clustering originates in mathematical morphology as a non-parametric topographic partitioning method treating a scalar image or negated probability density as an elevation relief. In the seminal digital immersion formulation of Vincent and Soille @vincent1991watersheds, one imagines piercing every regional minimum of a surface and slowly immersing it in water: rising water floods outward from each minimum to form *catchment basins*, and wherever distinct basins are about to merge, a dam—the *watershed line*—separates them.

Let $G = (V, E)$ be a finite, connected, undirected neighborhood graph over $N = |V|$ elements. In image segmentation, $V$ is a 2D/3D pixel lattice (4-, 8-, or 26-connected) with elevation $f: V -> bb(R)$ given by gradient magnitude $norm(nabla I(p))_2$ or negated distance transform; in general vector clustering over $\{bold(x)_1, dots, bold(x)_N\} subset bb(R)^d$, $G$ is a $k$-NN or $epsilon$-graph and $f(p) = -hat(p)_h (bold(x)_p)$ is a negated kernel density estimate mapping density modes to valleys.

Enumerate the distinct values of $f(V)$ as $h_"min" = h_0 < h_1 < dots < h_L = h_"max"$. For any level $h in bb(R)$, define the *sublevel set* $T_h (f) = \{p in V : f(p) <= h\}$. A *flat zone* is a maximal connected constant-elevation component, and a *regional minimum* at level $h$ is a flat zone $M subset V$ whose external neighbors $q in N(M) without M$ all satisfy $f(q) > h$. For $A subset.eq V$, let $d_A (p, q)$ denote shortest-path geodesic hop distance in the induced subgraph $G[A]$ ($infinity$ if disconnected).

#definition(title: "Geodesic Influence Zones and Immersion Watershed Partition")[
  For disjoint connected subsets $B = union.big_(i=1)^k B_i subset.eq A subset.eq V$, the *geodesic influence zone* of $B_i$ in $A$ is
  $ "IZ"_A (B_i) = \{p in A : forall j != i, thick d_A (p, B_i) < d_A (p, B_j)\}, $
  with skeleton $"SKIZ"_A (B) = A without union.big_(i=1)^k "IZ"_A (B_i)$. Letting $M_h$ be the union of all regional minima at level $h$, the Vincent--Soille immersion basins evolve recursively over $r = 0, dots, L - 1$ as @vincent1991watersheds:
  $ X_(h_0) = T_(h_0) (f), quad X_(h_(r+1)) = M_(h_(r+1)) union "IZ"_(T_(h_(r+1)) (f)) (X_(h_r)). $
  The connected components of $X_(h_L)$ form the final *catchment basins* $\{B_1, dots, B_K\}$, and $W(f) = V without X_(h_L)$ is the *watershed ridge set*.
]

The cluster count $K$ is determined non-parametrically by the number of regional minima (or user-supplied markers in marker-controlled watershed).

#algbox(
  title: "Vincent–Soille Immersion Watershed Clustering",
  inputs: [Connected graph $G = (V, E)$, elevation $f: V -> bb(R)$, optional $h$-minima smoothing parameter $delta_h >= 0$.],
  outputs: [Labels $"lab"(p) in \{1, dots, K\} union \{"WSHED"\}$ for all $p in V$.],
)[
  + *(Optional Denoising)* If $delta_h > 0$, replace $f$ with morphological reconstruction by erosion $R_f^"er" (f + delta_h)$ to prune minima of dynamics $<= delta_h$.
  + Sort vertices in non-decreasing order of $f(p)$ (radix/bucket sort for $G$ discrete gray levels or $O(N log N)$ sort for continuous $f$).
  + Initialize $"lab"(p) <- "INIT"$, $"dist"(p) <- 0$ for all $p$, $K <- 0$, FIFO queue $Q <- emptyset$.
  + *For each* distinct elevation level $h in \{h_0, dots, h_L\}$:
    - *(Mask & Seed)* For each $p$ with $f(p) = h$: set $"lab"(p) <- "MASK"$; if any $q in N(p)$ has $"lab"(q) > 0$ or `"WSHED"`, set $"dist"(p) <- 1$ and push $p$ to $Q$.
    - *(Geodesic BFS Flooding)* Push sentinel into $Q$, set $d_"cur" <- 1$. Pop $p$ until $Q$ empties: on sentinel, increment $d_"cur"$ and re-push if $Q != emptyset$; else propagate labels from neighbors $q$ with $"dist"(q) < d_"cur"$ (marking basin ties `"WSHED"`), and enqueue `"MASK"` neighbors at distance $d_"cur" + 1$.
    - *(New Regional Minima)* For each $p$ at level $h$ still `"MASK"`: set $K <- K + 1$, $"lab"(p) <- K$, and BFS-label its connected `"MASK"` flat zone as $K$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five formal guarantees: strict ordinal invariance, exact equivalence to minima-rooted Minimum Spanning Forests (MSF), zero-threshold noise breakdown with exact $h$-minima persistence pruning, plateau parity/thickening counterexamples, and wavefront algorithmic complexity.

#subsec("A", "Ordinal Invariance Under Monotonic Elevation Rescaling")

#theorem(title: "Exact Invariance Under Strictly Increasing Elevation Transformations", origin: "Original Proof")[
  For any connected graph $G = (V, E)$, elevation $f: V -> bb(R)$, and strictly increasing map $phi: bb(R) -> bb(R)$, the transformed landscape $g = phi compose f$ induces identical catchment basins $X_(phi(h)) (g) = X_h (f)$ at every level $h$ and an identical watershed set $W(g) = W(f)$.
]
#proof[
  Strict monotonicity maps the ordered levels $h_0 < dots < h_L$ of $f$ identically onto $phi(h_0) < dots < phi(h_L)$ of $g$. Because $phi(f(p)) <= phi(h_r) <==> f(p) <= h_r$, every sublevel set coincides: $T_(phi(h_r)) (g) = T_(h_r) (f)$. Likewise, $g(p) = g(q) <==> f(p) = f(q)$ and $g(q) > g(p) <==> f(q) > f(p)$, preserving all flat zones and regional minima ($M_(phi(h_r)) (g) = M_(h_r) (f)$). By induction on $r$: the base $X_(phi(h_0)) (g) = T_(phi(h_0)) (g) = X_(h_0) (f)$ holds immediately; assuming $X_(phi(h_r)) (g) = X_(h_r) (f)$, equality of induced graphs $G[T_(phi(h_(r+1))) (g)] = G[T_(h_(r+1)) (f)]$ yields identical geodesic distances and influence zones, so $X_(phi(h_(r+1))) (g) = X_(h_(r+1)) (f)$. At $r = L$, taking complements gives $W(g) = W(f)$.
]

Theorem 1 guarantees that density watershed clustering on $-hat(p)_h$, $-log hat(p)_h$, or rank-normalized densities produces the exact same partition, and image watersheds are completely immune to monotonic contrast/gamma transforms.

#subsec("B", "Equivalence to Minima-Rooted Minimum Spanning Forests")

#theorem(title: "Equivalence to Minima-Rooted Minimum Spanning Forests on Injective Landscapes", origin: "Original Unified Proof")[
  Let $G = (V, E)$ be connected and $f: V -> bb(R)$ be injective with regional minima $M = \{m_1, dots, m_K\}$. Assign edge pass weights $w(u, v) = max(f(u), f(v))$ and define $N^-(p) = \{q in N(p) : f(q) < f(p)\}$. Then:
  + Every $p in V without M$ belongs to basin $B_i$ iff $N^-(p)$ intersects $B_i$ and no other $B_j$ ($j != i$); otherwise $p in W(f)$.
  + Each basin $B_i$ is spanned by a tree $T_i$ whose paths to $m_i$ are strictly decreasing in $f$.
  + Extending $union.big_(i=1)^K T_i$ by assigning each $p in W(f)$ to any $q in N^-(p)$ yields an $M$-rooted spanning forest $F^*$ achieving exact global minimum weight:
    $ w(F^*) = sum_(p in V without M) f(p) = min_(F "rooted at" M) sum_(e in E(F)) w(e). $
]
#proof[
  Injectivity makes every regional minimum a singleton $\{m_i\}$ and adds one vertex $p_r$ per level: $T_(h_r) (f) = T_(h_(r-1)) (f) union \{p_r\}$.

  *(Part 1)* Induct on $r$. If $p_r in M$, $N^-(p_r) = emptyset$, forming a new singleton basin. If $p_r in.not M$, every neighbor of $p_r$ in $G[T_(h_r) (f)]$ lies in $N^-(p_r) != emptyset$, so $d_(T_(h_r)) (p_r, B_j) = 1$ iff $N^-(p_r) inter B_j != emptyset$ (and $>= 2$ otherwise). Thus: if $N^-(p_r)$ hits a unique $B_i$, $d_(T_(h_r)) (p_r, B_i) = 1 < 2 <= d_(T_(h_r)) (p_r, B_j)$ ($j != i$), giving $p_r in B_i$; if $N^-(p_r)$ hits two distinct basins, distance $1$ ties place $p_r in W(f)$; if $N^-(p_r) subset.eq W(f)$, every shortest path from $p_r$ to any basin traverses some $w in N^-(p_r) subset.eq W(f)$ which by induction is equidistant to $>= 2$ basins inside $T_(f(w)) subset.eq T_(h_r)$, keeping $p_r$ tied in $W(f)$.

  *(Part 2)* By Part 1, every non-minimum $p in B_i$ has a parent $"par"(p) in N^-(p) inter B_i$ with $f("par"(p)) < f(p)$; choosing one parent per vertex yields an acyclic directed spanning tree $T_i$ rooted at $m_i$.

  *(Part 3)* Orienting any $M$-rooted spanning forest $F$ toward its roots assigns each $p in V without M$ a unique outgoing edge $e_p = (p, "par"_F (p))$ (a bijection $V without M <-> E(F)$), so
  $ w(F) = sum_(p in V without M) max(f(p), f("par"_F (p))) >= sum_(p in V without M) f(p). $
  Because $F^*$ links every $p in V without M$ strictly downward ($f("par"(p)) < f(p)$), $w(p, "par"(p)) = f(p)$ everywhere, attaining the exact global minimum @gower1969minimum.
]

#subsec("C", "Noise-Induced Over-Segmentation and h-Minima Persistence")

#theorem(title: "Zero-Threshold Noise Breakdown and Exact h-Minima Persistence Bound", origin: "Original Derivation")[
  + *(Catastrophic Over-Segmentation)* Let $f_0 equiv c$ on a $d$-regular interior region $V_0 subset.eq V$ ($|V_0| = N_0$) be perturbed by i.i.d. continuous noise $epsilon_p$ of arbitrarily small variance $sigma^2 > 0$. Then $bb(E)[K_"noise"] = N_0 \/ (d + 1)$ (i.e., $N_0 \/ 5$ on 4-connected grids and $N_0 \/ 9$ on 8-connected grids), independent of $sigma > 0$.
  + *(Persistence Filtering)* Define the topological dynamics of non-global minimum $M_i$ by
    $ "dyn"(M_i) = inf_(gamma : M_i -> q, thick f(q) < f(M_i)) (max_(v in gamma) f(v) - f(M_i)). $
    Erosion reconstruction $tilde(f) = R_f^"er" (f + delta_h)$ eliminates every minimum with $"dyn"(M_i) <= delta_h$ and preserves every minimum with $"dyn"(M_i) > delta_h$. Under bounded noise $norm(epsilon)_infinity < eta$, setting $delta_h = 2 eta$ removes all flat-region spurious minima while preserving every true basin of depth $Delta > 4 eta$.
]
#proof[
  *(Part 1)* On each interior vertex's closed neighborhood $overline(N)(p)$ of size $d + 1$, exchangeability of i.i.d. continuous perturbations (zero ties a.s.) implies $"Pr"[p "is a minimum"] = 1 \/ (d + 1)$. Linearity of expectation gives $bb(E)[K_"noise"] = N_0 \/ (d + 1)$ for every $sigma > 0$.

  *(Part 2)* Erosion reconstruction $tilde(f) = R_f^"er" (f + delta_h)$ equals the minimax pass elevation from the lifted marker $f + delta_h$:
  $ tilde(f)(p) = min_(q in V) max(f(q) + delta_h, thick min_(gamma : p -> q) max_(v in gamma) f(v)). $
  If $"dyn"(M_i) <= delta_h$, some path $gamma$ from $M_i$ to $q$ ($f(q) < f(M_i)$) has saddle $s = max_(v in gamma) f(v) <= f(M_i) + delta_h$; every vertex on the sub-path to the saddle attains $tilde(f) <= max(f(q) + delta_h, s) <= f(M_i) + delta_h$, filling $M_i$ into $q$'s basin so $M_i$ ceases to be a regional minimum. Conversely, if $"dyn"(M_i) > delta_h$, every escape path to a deeper point rises strictly above $f(M_i) + delta_h$, leaving $M_i$ a strict regional minimum at level $f(M_i) + delta_h$. Under $|epsilon_p| < eta$, any spurious flat-region minimum has dynamics $< (c + eta) - (c - eta) = 2 eta = delta_h$ (eliminated), while a true basin of depth $Delta > 4 eta$ has perturbed dynamics $>= (c_0 + Delta - eta) - (c_0 + eta) = Delta - 2 eta > 2 eta = delta_h$ (preserved).
]

#subsec("D", "Plateau Parity and Quadratic Ridge Thickening")

#proposition(title: "Plateau Parity Non-Separation and O(L^2) Watershed Thickening", origin: "Original Counterexample")[
  On discrete landscapes with flat plateaus, Vincent--Soille immersion $W(f)$ exhibits two geometric anomalies: (i) on an even-length 1D plateau of $2 r$ vertices between two basins, $W(f) = emptyset$ (adjacent vertices carry different basin labels without a watershed pixel); and (ii) on a 2D 4-connected grid $\{0, dots, L\}^2$, just two single-pixel minima can induce a thick 2D watershed block of exact cardinality $L^2$.
]
#proof[
  *(i)* On path $v_0 - v_1 - dots - v_(2 r + 1)$ with $f(v_0) = f(v_(2 r + 1)) = 0$ and $f(v_1) = dots = f(v_(2 r)) = 1$, at $h = 1$ we have $d_V (v_i, v_0) = i$ and $d_V (v_i, v_(2 r + 1)) = 2 r + 1 - i$. Because $2 r + 1$ is odd, $i != 2 r + 1 - i$ for all $i$, assigning $\{v_1, dots, v_r\} subset.eq B_1$ and $\{v_(r+1), dots, v_(2 r)\} subset.eq B_2$ with $W(f) = emptyset$.

  *(ii)* On the 4-connected grid $\{0, dots, L\}^2$, set two minima $f(0, 1) = f(1, 0) = 0$, block the corner $f(0, 0) = 2$, and set plateau height $f(x, y) = 1$ everywhere else. At $h = 1$, every vertex $(x, y) in \{1, dots, L\}^2$ has exact $ell_1$ geodesic distance $x + (y - 1) = x + y - 1$ to $(0, 1)$ and $(x - 1) + y = x + y - 1$ to $(1, 0)$ (avoiding the blocked $(0, 0)$). Hence all $L^2$ interior vertices of $\{1, dots, L\}^2$ tie simultaneously and belong to $W(f)$.
]

#subsec("E", "Correctness and Linear / Log-Linear Complexity")

#theorem(title: "Wavefront Invariant and O(N + G) / O(N log N) Complexity", origin: "Classical — Cited")[
  The sentinel-delimited FIFO queue in Vincent--Soille immersion @vincent1991watersheds preserves the exact geodesic breadth-first invariant at every level $h$, computing $X_(h_L) (f)$ in $O(N + |E| + G)$ time for $G$ discrete integer levels (radix sort) and $O(N log N + |E|)$ time for real elevations, with $O(N + |E|)$ space.
]
#proof[
  By induction on $d_"cur"$, vertices popped between sentinels $d_"cur" - 1$ and $d_"cur"$ lie at exact geodesic distance $d_"cur"$ from $X_(h^-) union W_(h^-)$ in $G[T_h (f)]$ and inspect all neighbors at distance $< d_"cur"$ before finalizing labels, realizing exact $"IZ"_(T_h (f)) (X_(h^-))$ flooding @vincent1991watersheds. Each vertex enters the level scan once, is queued at most twice (geodesic propagation and new-minimum BFS), and has incident edges inspected $O(1)$ times, yielding $O(N + |E|)$ flooding time plus $O(N + G)$ radix or $O(N log N)$ comparison sorting.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Watershed runs in $O(N log N + |E|)$ time for continuous densities on $k$-NN graphs ($|E| = O(k N)$) and linear $O(N + G)$ time on quantized 2D/3D lattices ($|E| = O(N)$), using $O(N + |E|)$ space ($O(N)$ on implicit grids). It makes *zero convexity assumptions*: catchment basins conform to curved manifolds, rings, and filaments bounded by elevation ridges (or density troughs). Rather than choosing $K$, practitioners control topological scale via $h$-minima filtering $delta_h$ (Theorem 3) or marker seeding. On high-dimensional tabular data, quality depends strongly on $k$-NN graph fidelity and kernel bandwidth.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$O(N + G)$ on quantized grids ($G$ levels); $O(N log N + |E|)$ on general graphs],
    [Space Complexity], [$O(N)$ on implicit pixel/voxel grids; $O(N + |E|)$ on explicit $k$-NN graphs],
    [Inductive Bias], [Topographic basins of attraction separated by pass ridges / density troughs],
    [Cluster Shape], [Arbitrary non-convex, curved, and elongated geometries],
    [Invariance], [Exact invariance to all strictly increasing transforms $phi compose f$ (Theorem 1)],
    [Noise Robustness], [Zero breakdown threshold raw ($N_0 \/ (d+1)$); stable after $h$-minima pruning],
    [Cluster Count $K$], [Non-parametric (equals regional minima count or user seed markers)],
    [Boundary Structure], [Vertex $W(f)$ may omit even-plateau cuts or thicken to $Theta(L^2)$ on 2D plateaus],
  ),
  caption: [Operational and theoretical profile of Watershed Clustering (Alg. 167).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

+ *Biomedical Cell & Nuclear Instance Segmentation:* Overlapping fluorescent nuclei share identical intensities, defeating color clustering. Negating the Euclidean distance transform $f(p) = -d(p, partial Omega_"fg")$ turns each cell core into a minimum and each contact neck into a saddle, letting marker-controlled watershed split clumps cleanly.
+ *Materials Science & 3D Micro-CT Pore Partitioning:* In porous rock, battery electrodes, and metallic foams, negated 3D distance-transform watershed segments interconnected void spaces precisely at minimal throat cross-sections matching pore-network flow physics.
+ *Terrain Hydrology & DEM Catchment Delineation:* On digital elevation models, Theorem 2 guarantees that immersion traces strictly descending drainage paths into sinks, delineating physical watersheds and river divides.
+ *Airborne LiDAR Individual Tree Crown Segmentation:* Inverting Canopy Height Models converts treetops into valleys and inter-crown gaps into ridges; $h$-minima pruning (Theorem 3) suppresses branch noise while isolating individual crowns.
+ *Non-Parametric Feature-Space Density Mode Clustering:* Evaluating $f(i) = -hat(p)_h (bold(x)_i)$ on a $k$-NN graph extracts arbitrary-shape density modes without iterative gradient ascent @comaniciu2002mean @hinneburg1998denclue and with exact log-density invariance (Theorem 1).

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Key failure modes include: (1) *Noise-Induced Over-Segmentation* generating $approx N \/ (d + 1)$ basins raw (Theorem 3); (2) *Basin Leakage* across a single weak boundary gap allowing background to flood foreground; (3) *High-D Graph Sparsity* causing noisy density spikes; and (4) *Plateau Thickening* leaving unassigned $W(f)$ blocks (Proposition 1).

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Touching Cells / Micro-CT Pores], [Excellent: distance-transform valleys isolate necks cleanly], [Apply $h$-minima or deep markers to curb over-splits],
    [DEM / LiDAR Canopy Mapping], [Exact: relief directly matches physical runoff / crowns], [Fill DEM micro-pits via erosion reconstruction],
    [Natural Textured Images], [Moderate: raw gradient yields severe over-segmentation], [Pre-smooth or use SLIC @achanta2012slic for regular tiles],
    [Non-Convex Low-D Feature Data], [High: negated KDE on $k$-NN graph captures curved shapes], [Compare HDBSCAN @campello2013hdbscan or DPC @rodriguez2014dpc],
    [High-D Sparse Tabular Data], [Low--Moderate: noisy density minima & graph hubs], [Reduce dimension first or use spectral / density trees],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for Watershed Clustering.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Watershed unifies mathematical morphology, hierarchical density clustering, and mode-seeking vector algorithms @hartigan1981consistency @comaniciu2002mean @rodriguez2014dpc @campello2013hdbscan.

#proposition(title: "Links to Mode Seeking, Density Peaks, and Single-Linkage Trees", origin: "Original Proof")[
  On a neighborhood graph $G$ with injective negated density $f(p) = -hat(p)_h (bold(x)_p)$: (i) every basin $B_i$ discretizes a Mean Shift @comaniciu2002mean / DENCLUE @hinneburg1998denclue basin of attraction via strictly density-increasing paths to mode $m_i$ (Theorem 2); (ii) the saddle-merge tree of catchment basins is isomorphic to both Hartigan's high-density cluster tree @hartigan1981consistency and Single-Linkage @sibson1973slink on $w(u, v) = max(f(u), f(v))$; and (iii) linking each $p$ to $arg min_(q in N^-(p)) f(q)$ yields an optimal forest mirroring Density Peaks parent pointers @rodriguez2014dpc.
]
#proof[
  Theorem 2(2) ensures every $p in B_i$ reaches mode $m_i$ along a strictly increasing $hat(p)_h$-path, discretizing Mean Shift gradient ascent. Because $w(u, v) <= h <==> f(u) <= h and f(v) <= h$, connected components of $G[T_h (f)]$ coincide identically with edges $\{e : w(e) <= h\}$, making the basin-merge tree isomorphic to Single-Linkage under $w$ and the superlevel connectivity tree $\{hat(p)_h >= -h\}$ @hartigan1981consistency. Selecting steepest neighbor $"par"(p) = arg min_(q in N^-(p)) f(q)$ achieves the optimal weight of Theorem 2(3), matching Density Peaks' higher-density pointer tree @rodriguez2014dpc.
]

#corollary(title: "Persistence Stability Under Uniform Relief Perturbations", origin: "Original Proof")[
  Let $tilde(f) : V -> bb(R)$ be a noisy relief satisfying $norm(tilde(f) - f)_infinity <= delta$. Suppose every true catchment basin $B_i$ of $f$ has topological persistence (saddle height minus minimum height) $pi_i := h_"sad" (B_i) - f(m_i) > 2 delta$, while every spurious noise basin has persistence $< pi_"noise"$. Then $h$-minima watershed pruning with threshold $h in (pi_"noise" + 2 delta, min_i pi_i - 2 delta)$ recovers the exact number of true basins and preserves every true marker up to $2 delta$-sublevel containment.
]
#proof[
  For any vertex $p in V$, $f(p) - delta <= tilde(f)(p) <= f(p) + delta$, so every sublevel set obeys the interleaving $T_(alpha - delta) (f) subset.eq T_alpha (tilde(f)) subset.eq T_(alpha + delta) (f)$. By the stability of 0-dimensional persistent homology (the Elder Rule on component merge trees), the birth height $tilde(b)_i$ and death (saddle) height $tilde(d)_i$ of every basin shift by at most $delta$: $|tilde(b)_i - f(m_i)| <= delta$ and $|tilde(d)_i - h_"sad" (B_i)| <= delta$. Consequently, the perturbed persistence $tilde(pi)_i = tilde(d)_i - tilde(b)_i$ satisfies $|tilde(pi)_i - pi_i| <= 2 delta$. Every true basin therefore retains $tilde(pi)_i >= min_j pi_j - 2 delta > h$, whereas every noise basin has $tilde(pi)_"noise" <= pi_"noise" + 2 delta < h$ and is annihilated during $h$-minima reconstruction. Thus the surviving marker count equals the exact ground-truth basin count with zero topological false merges or splits.
]

On the *ClusterOverview* benchmark, Watershed Clustering achieves *Iris* Silhouette @rousseeuw1987silhouettes $0.586$ and ARI @hubert1985comparing $0.567$, *Customers* Silhouette $0.247$, *Multi-Shape* Silhouette $0.712$ (ARI N/A), and *High-D* Silhouette $0.302$ (ARI N/A), with mean runtime $1.376$ ms. Notably, Watershed attains the *highest* Silhouette among Algs. 166--169 on *Iris* ($0.586$), *Multi-Shape* ($0.712$), and *High-D* ($0.302$). Its $0.712$ *Multi-Shape* score directly reflects Theorem 2's non-convex topographic basins tracing curved manifolds along density troughs, while its lower *Customers* score ($0.247$) reflects Theorem 3's sensitivity to flat weakly modal distributions lacking sharp ridges. Its $1.376$ ms runtime matches Theorem 4's $O(N log N + |E|)$ graph-flooding cost.
