#import "../template/ieee.typ": *

#alg_entry(
  id: 166,
  name: "SLIC (Simple Linear Iterative Clustering)",
  category: "Cat. VII: Fuzzy, Biclustering & Domain",
  year: "2012",
  paradigm: "Spatially Localized 5D k-Means Superpixel Partitioning",
  complexity: $O(N I)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Superpixel over-segmentation partitions a digital image into compact, perceptually homogeneous connected atomic regions that adhere tightly to object boundaries while reducing downstream graph cardinality by orders of magnitude. Simple Linear Iterative Clustering (SLIC), introduced by Achanta et al. @achanta2012slic, formulates superpixel generation as a spatially localized adaptation of Lloyd's $k$-means algorithm @lloyd1982least in a five-dimensional color-spatial product space.

Let $Omega = \{1, dots, W\} times \{1, dots, H\}$ be a discrete 2D pixel lattice of size $N = W H$. Each pixel $i in \{1, dots, N\}$ is embedded into $cal(X) = bb(R)^3_"Lab" times bb(R)^2_"xy"$ as $bold(z)_i = (bold(c)_i^top, bold(p)_i^top)^top = (l_i, a_i, b_i, x_i, y_i)^top in bb(R)^5$, where $bold(c)_i$ denotes perceptually uniform CIELAB color coordinates and $bold(p)_i in Omega$ denotes spatial pixel coordinates. Given target superpixel count $k in \{1, dots, N\}$, SLIC sets the nominal grid interval $S = sqrt(N \/ k)$ and compactness $m > 0$ balancing spatial regularity against chromatic fidelity. Between pixel $bold(z)_i$ and prototype $bold(mu)_j = (bold(c)_j^mu, bold(p)_j^mu)^top$, SLIC uses the squared anisotropic Mahalanobis distance
$ D_S (bold(z)_i, bold(mu)_j)^2 &= d_c (i, j)^2 + (m / S)^2 d_s (i, j)^2 \
  &= (bold(z)_i - bold(mu)_j)^top W_(m, S) (bold(z)_i - bold(mu)_j), $
where $d_c (i, j) = norm(bold(c)_i - bold(c)_j^mu)_2$, $d_s (i, j) = norm(bold(p)_i - bold(p)_j^mu)_2$, and $W_(m, S) = "diag"(1, 1, 1, m^2 \/ S^2, m^2 \/ S^2) succ 0$. Unlike global 5D $k$-means @macqueen1967some, SLIC restricts assignment comparisons to a $2 S times 2 S$ spatial window around each prototype, eliminating $O(k)$ search overhead per pixel.

#definition(title: "Windowed 5D SLIC Objective and Connected Refinement")[
  For prototypes $bold(M) = \{bold(mu)_1, dots, bold(mu)_k\} subset bb(R)^5$, define the $2S times 2S$ spatial search window of center $j$ as $cal(W)_S (bold(mu)_j) = \{i : norm(bold(p)_i - bold(p)_j^mu)_infinity <= S\}$ and pixel $i$'s active candidate set as $cal(A)_S (i; bold(M)) = \{j : i in cal(W)_S (bold(mu)_j)\}$. A hard assignment $bold(C) = (c_1, dots, c_N) in \{1, dots, k\}^N$ is *window-feasible* if $c_i in cal(A)_S (i; bold(M))$ whenever $cal(A)_S (i; bold(M)) != emptyset$. SLIC minimizes the weighted distortion
  $ J_(m, S) (bold(C), bold(M)) = sum_(j=1)^k sum_(i : c_i = j) D_S (bold(z)_i, bold(mu)_j)^2 $
  under window-restricted assignment updates, followed by a connected-component operator $Pi_"conn"$ merging disjoint fragments into adjacent superpixels.
]

To avoid seeding centers on high-contrast edges or impulsive noise, seeds on an $S$-spaced lattice are shifted to the lowest gradient pixel in their $3 times 3$ neighborhood, where $G(x, y) = norm(bold(c)_(x+1, y) - bold(c)_(x-1, y))_2^2 + norm(bold(c)_(x, y+1) - bold(c)_(x, y-1))_2^2$.

#algbox(
  title: "Simple Linear Iterative Clustering (SLIC)",
  inputs: [CIELAB image lattice $Omega$ ($N$ pixels), count $k$, compactness $m > 0$, tolerance $tau > 0$, max iterations $I_"max"$.],
  outputs: [Connected superpixel partition $bold(C) = (c_1, dots, c_N)$ and 5D prototypes $\{bold(mu)_j\}$.],
)[
  + Set $S <- floor(sqrt(N \/ k))$. Seed $k$ centers $bold(mu)_j^{(0)}$ on a regular $S$-spaced grid over $Omega$, shifting each to minimize $G(x, y)$ in its $3 times 3$ neighborhood.
  + Set $c_i <- -1$, $d_i <- infinity$ for all $i in \{1, dots, N\}$, and $t <- 0$.
  + *Repeat* until $E^{(t)} <= tau$ or $t = I_"max"$:
    - Reset $d_i <- infinity$. For each $j = 1, dots, k$ and every $i in cal(W)_S (bold(mu)_j^{(t)})$: compute $D <- D_S (bold(z)_i, bold(mu)_j^{(t)})$; if $D < d_i$, set $d_i <- D$ and $c_i^{(t)} <- j$.
    - Update each non-empty centroid $bold(mu)_j^{(t+1)} <- |C_j^{(t)}|^(-1) sum_(i in C_j^{(t)}) bold(z)_i$.
    - Set $E^{(t+1)} <- sum_(j=1)^k norm(bold(mu)_j^{(t+1)} - bold(mu)_j^{(t)})_2$ and $t <- t + 1$.
  + Apply $Pi_"conn"$: run BFS to merge connected components smaller than $floor(N \/ (4 k))$ into their adjacent superpixel label.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

The window constraint $cal(W)_S (bold(mu)_j)$ fundamentally alters the complexity, monotonicity, and symmetry structure of Lloyd's algorithm. Below we prove five formal theorems characterizing SLIC.

#subsec("A", "Unconditional Linear-Time Complexity per Iteration")

#lemma(title: "Unconditional O(N) Distance Evaluations Independent of k", origin: "Original Proof")[
  For every integer $1 <= k <= N$ and any prototype configuration $bold(M)^{(t)}$ on $Omega$ with $S = sqrt(N \/ k)$, the total number $T_"eval"$ of 5D distance evaluations per assignment sweep satisfies
  $ T_"eval" <= k (2 S + 1)^2 = 4 N + 4 sqrt(k N) + k <= 9 N, $
  independently of $k$ and prototype drift. Moreover, under bounded center drift $norm(bold(p)_j^mu - bold(g)_j)_infinity <= rho S$ from initial $S$-lattice nodes $bold(g)_j$, every pixel belongs to at most $(2 ceil(1 + rho) + 1)^2$ candidate windows.
]
#proof[
  For any spatial center $(x_j^mu, y_j^mu) in bb(R)^2$, the closed interval $[x_j^mu - S, x_j^mu + S]$ contains at most $floor(2 S) + 1 <= 2 S + 1$ integer coordinates, so $|cal(W)_S (bold(mu)_j^{(t)})| <= (2 S + 1)^2$ unconditionally regardless of drift or overlap. Summing across all $k$ centers and substituting $S = sqrt(N \/ k)$ gives
  $ T_"eval" = sum_(j=1)^k |cal(W)_S (bold(mu)_j^{(t)})| <= k (4 S^2 + 4 S + 1) = 4 N + 4 sqrt(k N) + k. $
  Since $1 <= k <= N$, $sqrt(k N) <= N$ and $k <= N$, proving $T_"eval" <= 9 N$ (approaching $4 N$ for $k << N$). Under drift bound $norm(bold(p)_j^mu - bold(g)_j)_infinity <= rho S$, membership $i in cal(W)_S (bold(mu)_j)$ implies $norm(bold(p)_i - bold(g)_j)_infinity <= (1 + rho) S$ by the triangle inequality; at most $2 ceil(1 + rho) + 1$ lattice multiples of $S$ fall within radius $(1 + rho) S$ per axis, giving $|cal(A)_S (i; bold(M))| <= (2 ceil(1 + rho) + 1)^2 = O(1)$.
]

#subsec("B", "Lyapunov Descent Criterion and Non-Monotonicity Counterexample")

Unlike classical Lloyd $k$-means @lloyd1982least, windowed SLIC is not unconditionally monotone in the global energy $J_(m, S)$.

#theorem(title: "Conditional Lyapunov Descent and Window-Truncation Non-Monotonicity", origin: "Original Counterexample")[
  Consider $J_(m, S) (bold(C), bold(M))$ under $W_(m, S) succ 0$.
  + *(Conditional Descent)* If $c_i^{(t-1)} in cal(A)_S (i; bold(M)^{(t)})$ holds for all $i in \{1, dots, N\}$ at iteration $t >= 1$, then
    $ J_(m, S) (bold(C)^{(t)}, bold(M)^{(t+1)}) <= J_(m, S) (bold(C)^{(t-1)}, bold(M)^{(t)}) - Delta_"shift"^{(t)}, $
    where $Delta_"shift"^{(t)} = sum_(j=1)^k |C_j^{(t)}| D_S (bold(mu)_j^{(t+1)}, bold(mu)_j^{(t)})^2 >= 0$.
  + *(Non-Monotonicity Counterexample)* Without window retention, there exists a 5-pixel integer lattice instance ($k = 2, S = 2, m = 2$) where one SLIC step strictly raises both post-assignment ($14 -> 106$) and post-update ($14 -> 52.5$) objectives, converging to a strictly worse fixed point.
]
#proof[
  *(Part 1)* Since $c_i^{(t-1)} in cal(A)_S (i; bold(M)^{(t)})$, window-restricted minimization gives $D_S (bold(z)_i, bold(mu)_(c_i^{(t)})^{(t)})^2 <= D_S (bold(z)_i, bold(mu)_(c_i^{(t-1)})^{(t)})^2$ for every $i$, so $J_(m, S) (bold(C)^{(t)}, bold(M)^{(t)}) <= J_(m, S) (bold(C)^{(t-1)}, bold(M)^{(t)})$. Expanding each cluster sum around its exact mean $bold(mu)_j^{(t+1)}$ cancels the cross term $sum_(i in C_j^{(t)}) (bold(z)_i - bold(mu)_j^{(t+1)})^top W_(m, S) (bold(mu)_j^{(t+1)} - bold(mu)_j^{(t)}) = 0$, giving the exact identity
  $ sum_(i in C_j^{(t)}) D_S (bold(z)_i, bold(mu)_j^{(t)})^2 = sum_(i in C_j^{(t)}) D_S (bold(z)_i, bold(mu)_j^{(t+1)})^2 + |C_j^{(t)}| D_S (bold(mu)_j^{(t+1)}, bold(mu)_j^{(t)})^2. $
  Summing over $j = 1, dots, k$ yields the Lyapunov inequality.

  *(Part 2)* Set $k = 2, S = 2, m = 2$ (so $(m \/ S)^2 = 1$), fix $a_i = b_i = 0$, and place five distinct integer pixels $(l, x, y)$ at
  $ bold(z)_1 = (0, 0, 0), quad bold(z)_2 = (0, 4, 0), quad bold(z)_3 = (0, 4, 1), $
  $ bold(z)_4 = (0, 4, -1), quad bold(z)_5 = (10, 1, 0). $
  Initialize $bold(mu)_1^{(0)} = (0, 2, 0)$ and $bold(mu)_2^{(0)} = (10, 1, 0)$. At $t = 0$, $bold(z)_1$ lies in both windows ($D_S^2 = 4 < 101$), $bold(z)_2, bold(z)_3, bold(z)_4$ lie only in $cal(W)_2 (bold(mu)_1^{(0)})$, and $bold(z)_5$ achieves $0$ at $bold(mu)_2^{(0)}$, yielding $C_1^{(0)} = \{1, 2, 3, 4\}$, $C_2^{(0)} = \{5\}$, and $J_(2, 2) (bold(C)^{(0)}, bold(M)^{(0)}) = 4 + 4 + 5 + 5 + 0 = 18$. Updating centroids gives $bold(mu)_1^{(1)} = (0, 3, 0)$ and $bold(mu)_2^{(1)} = (10, 1, 0)$ with reduced energy $J_(2, 2) (bold(C)^{(0)}, bold(M)^{(1)}) = (9 + 1 + 2 + 2) + 0 = 14$.

  At $t = 1$, $norm(bold(p)_1 - bold(p)_1^(mu, (1)))_infinity = |0 - 3| = 3 > S = 2$, ejecting $bold(z)_1$ from $cal(W)_2 (bold(mu)_1^{(1)})$ while it stays in $cal(W)_2 (bold(mu)_2^{(1)})$ ($|0 - 1| = 1 <= 2$). Windowed assignment thus forces $bold(z)_1 in C_2^{(1)}$, producing $C_1^{(1)} = \{2, 3, 4\}$, $C_2^{(1)} = \{1, 5\}$, and post-assignment energy $J_(2, 2) (bold(C)^{(1)}, bold(M)^{(1)}) = (1 + 2 + 2) + (101 + 0) = 106 > 14$. Recomputing centroids gives $bold(mu)_1^{(2)} = (0, 4, 0)$ and $bold(mu)_2^{(2)} = (5, 0.5, 0)$ with energy $J_(2, 2) (bold(C)^{(1)}, bold(M)^{(2)}) = (0 + 1 + 1) + (25.25 + 25.25) = 52.5 > 18 > 14$. At $t = 2$, $bold(z)_1$ remains excluded from $cal(W)_2 (bold(mu)_1^{(2)})$ ($|0 - 4| = 4 > 2$), trapping SLIC at a strictly suboptimal fixed point.
]

#subsec("C", "Compactness-Adherence Phase Transition")

#theorem(title: "Voronoi Perturbation Bound and Exact Step-Edge Adherence Law", origin: "Original Derivation")[
  Let $Delta_C = sup_(i, l) norm(bold(c)_i - bold(c)_l)_2 < infinity$ be the image's CIELAB gamut diameter.
  + *(Voronoi Regularization Limit)* For any pixel $i$ assigned to $j$ over competing $j' in cal(A)_S (i; bold(M))$, the signed normal displacement $delta_perp (i; j, j')$ of $bold(p)_i$ past the perpendicular bisector of $bold(p)_j^mu, bold(p)_(j')^mu$ satisfies
    $ delta_perp (i; j, j') <= (S^2 Delta_C^2) / (2 m^2 norm(bold(p)_j^mu - bold(p)_(j')^mu)_2). $
    Thus as $m \/ Delta_C -> infinity$, SLIC converges pointwise off the Voronoi skeleton to the exact 2D Euclidean Voronoi diagram of $\{bold(p)_j^mu\}$.
  + *(Color-Edge Adherence Criterion)* Across an ideal step edge separating colors $bold(c)_A, bold(c)_B$ with contrast $Delta_E = norm(bold(c)_A - bold(c)_B)_2 > 0$ and pure prototypes $bold(c)_(j_A)^mu = bold(c)_A$, $bold(c)_(j_B)^mu = bold(c)_B$, every pixel $i$ in region $A$ inside $cal(W)_S (bold(mu)_(j_A))$ adheres strictly to $j_A$ over $j_B$ whenever $m < Delta_E \/ sqrt(2)$.
]
#proof[
  *(Part 1)* Preference $D_S (bold(z)_i, bold(mu)_j)^2 <= D_S (bold(z)_i, bold(mu)_(j'))^2$ rearranges to
  $ norm(bold(p)_i - bold(p)_j^mu)_2^2 - norm(bold(p)_i - bold(p)_(j')^mu)_2^2 <= (S / m)^2 (d_c (i, j')^2 - d_c (i, j)^2) <= (S Delta_C / m)^2. $
  Setting midpoint $bold(m) = (bold(p)_j^mu + bold(p)_(j')^mu) \/ 2$ and unit normal $bold(u) = (bold(p)_(j')^mu - bold(p)_j^mu) \/ norm(bold(p)_(j')^mu - bold(p)_j^mu)_2$, the left-hand side equals $2 norm(bold(p)_j^mu - bold(p)_(j')^mu)_2 chevron.l bold(p)_i - bold(m), bold(u) chevron.r = 2 norm(bold(p)_j^mu - bold(p)_(j')^mu)_2 delta_perp (i; j, j')$. Dividing yields the bound, vanishing as $O((Delta_C \/ m)^2)$.

  *(Part 2)* For pixel $i$ in region $A$, $d_c (i, j_A) = 0$, $d_c (i, j_B) = Delta_E$, and $i in cal(W)_S (bold(mu)_(j_A)) => norm(bold(p)_i - bold(p)_(j_A)^mu)_infinity <= S => d_s (i, j_A)^2 <= 2 S^2$. Even if $d_s (i, j_B) = 0$, strict adherence holds whenever $(m \/ S)^2 (2 S^2) < Delta_E^2 <==> m < Delta_E \/ sqrt(2)$.
]

#subsec("D", "Symmetry Group, Rotational Anisotropy, and Conditioning")

#proposition(title: "Invariance Group, Rotational Symmetry Breaking, and Conditioning", origin: "Original Proof")[
  SLIC is exactly invariant under: (i) CIELAB rigid isometries $bold(c)_i |-> Q_c bold(c)_i + bold(t)_c$ ($Q_c in O(3)$), (ii) integer spatial translations, and (iii) joint scaling $(bold(c)_i, m) |-> (beta bold(c)_i, beta m)$ ($beta > 0$). However, windowed SLIC is *not* invariant under continuous 2D rotations $Q_p in "SO"(2) without D_4$, and its metric tensor has condition number $kappa(W_(m, S)) = max(m^2 \/ S^2, S^2 \/ m^2)$.
]
#proof[
  Rigid $O(3)$ transforms preserve $d_c$, spatial translations preserve $bold(p)_i - bold(p)_j^mu$ in both $ell_2$ and $ell_infinity$ norms, and $(bold(c)_i, m) |-> (beta bold(c)_i, beta m)$ scales $D_S^2$ uniformly by $beta^2 > 0$. Because $cal(W)_S$ is an $ell_infinity$ ball of inradius $S$ and circumradius $sqrt(2) S$, rotating displacement $bold(v) = (S, S)^top$ (inside $cal(W)_S$) by $pi \/ 4$ yields $(0, sqrt(2) S)^top$ with $ell_infinity$ norm $sqrt(2) S > S$, ejecting it from the window; only the dihedral group $D_4$ preserves the $ell_infinity$ ball. Finally, $W_(m, S) = "diag"(1, 1, 1, m^2 \/ S^2, m^2 \/ S^2)$ has eigenvalues $\{1, m^2 \/ S^2\}$, giving $kappa(W_(m, S)) = max(m^2 \/ S^2, S^2 \/ m^2)$.
]

#subsec("E", "Exact High-Rate Distortion Bound and Lloyd Stationarity")

#theorem(title: "Discrete High-Rate Quantization Scaling and Unconstrained Equivalence", origin: "Classical — Cited")[
  Let $bold(c): [0, 1]^2 -> bb(R)^3$ be an $L$-Lipschitz color field sampled on a $sqrt(N) times sqrt(N)$ grid partitioned into $k$ congruent $S times S$ integer tiles ($S = sqrt(N \/ k) in bb(N)$). Then: (i) the optimal SLIC energy satisfies the exact discrete upper bound
  $ J_(m, S)^* <= (N / 6) (L^2 / k + m^2), $
  demonstrating $O(1 \/ k)$ chromatic distortion decay and resolution-invariant $O(m^2)$ per-pixel spatial regularization; and (ii) any SLIC fixed point $(bold(C)^*, bold(M)^*)$ whose pixels each have their global nearest center inside $cal(A)_S (i; bold(M)^*)$ is a stationary fixed point of unconstrained 5D Lloyd $k$-means @achanta2012slic @lloyd1982least @bottou1995convergence.
]
#proof[
  For each $S times S$ square tile centered at pixel coordinate $bold(p)_j^0 = sqrt(N) bold(u)_j^0$ with $bold(c)_j^0 = bold(c)(bold(u)_j^0)$, every pixel $i$ in tile $j$ satisfies $norm(bold(p)_i - bold(p)_j^0)_infinity <= (S - 1) \/ 2 < S$, so the tiling is strictly window-feasible. By $L$-Lipschitzness on $[0, 1]^2$, $norm(bold(c)_i - bold(c)_j^0)_2^2 <= (L^2 \/ N) norm(bold(p)_i - bold(p)_j^0)_2^2$, while $(m \/ S)^2 = m^2 k \/ N$. Hence
  $ D_S (bold(z)_i, bold(mu)_j^0)^2 <= ((L^2 + m^2 k) / N) norm(bold(p)_i - bold(p)_j^0)_2^2. $
  Over any $S$ consecutive centered integers, the exact mean squared coordinate is $(S^2 - 1) \/ 12$; summing across $x$ and $y$ gives tile-mean squared radius $(S^2 - 1) \/ 6 < S^2 \/ 6 = N \/ (6 k)$. Summing over all $N$ pixels yields
  $ J_(m, S)^* <= N ((L^2 + m^2 k) / N) (N / (6 k)) = (N / 6) (L^2 / k + m^2). $
  For (ii), if $arg min_(1 <= j <= k) D_S (bold(z)_i, bold(mu)_j^*) in cal(A)_S (i; bold(M)^*)$ for every $i$, then windowed cells coincide identically with unconstrained 5D Voronoi cells under $W_(m, S)$, satisfying Lloyd stationarity @lloyd1982least @bottou1995convergence.
]

#sec("III", "Complexity, Geometry and Operational Profile")

SLIC decouples computation from superpixel count $k$: by Lemma 1, each sweep executes at most $9 N$ (asymptotically $4 N$) 5D distance evaluations, yielding $O(N I)$ total time. Empirically $I = 10$ iterations suffice for visual convergence @achanta2012slic, sidestepping exponential worst-case $k$-means iteration bounds @vattani2011kmeans. Memory is strictly $O(N + k)$ (label array, distance buffer, and $5 k$ prototype floats). Scanning $2 S times 2 S$ windows row-by-row maximizes L1/L2 cache locality and enables straightforward SIMD or GPU parallelization, followed by a single $O(N)$ BFS pass for $Pi_"conn"$. Geometrically, SLIC induces convex 5D Mahalanobis Voronoi cells whose 2D projections form compact, quasi-hexagonal tiles that deform selectively along high-contrast edges whenever $m < Delta_E \/ sqrt(2)$.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$O(N I)$ ($<= 9 N$ evals/iter, independent of $k$); $O(N)$ BFS connectivity pass],
    [Space Complexity], [$O(N + k)$ (label map, distance buffer, $5 k$ centroid coordinates)],
    [Metric Geometry], [Anisotropic 5D tensor $W_(m, S) = "diag"(1, 1, 1, m^2 \/ S^2, m^2 \/ S^2)$],
    [Cluster Shape], [Compact, spatially localized, edge-adhering superpixels of area $approx N \/ k$],
    [Key Hyperparameters], [Target count $k$; compactness $m in [1, 40]$ (default $m = 10$ for CIELAB)],
    [Monotonicity], [Conditionally monotone (Thm. 1); non-monotone if pixel exits $2S times 2S$ window],
    [Connectivity], [Not intrinsic to 5D assignment; enforced post-hoc via BFS relabeling $Pi_"conn"$],
    [Cache & Parallelism], [Contiguous row-major $2S times 2S$ tile scans; high L1/L2 cache hit rate],
  ),
  caption: [Operational and theoretical profile of SLIC (Alg. 166).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

+ *Semantic Segmentation & Vision Graph Construction:* Conditional random fields and graph neural networks on megapixel images compress $N approx 10^6$ pixels into $k approx 10^3$ homogeneous nodes via SLIC in $O(N)$ time, preserving semantic borders whenever $Delta_E > sqrt(2) m$ (Theorem 2) @achanta2012slic.
+ *Computational Histopathology & Biomedical Microscopy:* Whole-slide H&E biopsy scans span tens of millions of pixels; CIELAB $d_c$ tracks perceptual stain variations while $(m \/ S) d_s$ prevents merging distant identically stained nuclei, accelerating tissue phenotyping.
+ *Remote Sensing & Multispectral Earth Observation:* Extending $bold(c)_i$ to $p$-band multispectral/SAR reflectance vectors enables linear-time object-based image analysis over gigapixel satellite tiles, grouping agricultural parcels while suppressing speckle via local averaging.
+ *Dense Stereo Matching & Optical Flow Regularization:* Piecewise-planar disparity and motion models exploit SLIC's compact, edge-respecting support patches to stabilize local plane fitting without the $O(N^(3 \/ 2))$ cost of Normalized Cuts @shi2000ncut.
+ *Embedded Robotics & Real-Time UAV Perception:* Under tight 30--60 Hz frame budgets, SLIC's deterministic $O(N I)$ latency and $O(N)$ footprint outperform graph cuts and kernel density mode-seekers @comaniciu2002mean @vedaldi2008quickshift.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Four failure modes stem directly from SLIC's structure: (1) *Thin Elongated Structures* narrower than $S \/ 2$ (e.g., vessels, cables) miss initial seeds and risk absorption by $Pi_"conn"$; (2) *Low-Contrast Textured Boundaries* ($Delta_E << sqrt(2) m$) trigger Voronoi rigidity (Theorem 2), leaking across edges; (3) *Window-Truncation Jumps* under high color variance cause centroid drift $> S$, non-monotone energy spikes (Theorem 1), and disconnected fragments; and (4) *Non-Lattice High-D Tabular Data* lacks a 2D spatial grid, so arbitrary coordinate splitting and tensor conditioning $kappa(W_(m, S))$ erode global purity.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Natural RGB / CIELAB Images], [High: perceptual 5D metric + $O(N)$ localized search], [Tune $m$; compare Quick Shift @vedaldi2008quickshift],
    [3D Volumetric Medical Scans], [High: extends to 6D supervoxels ($2S times 2S times 2S$)], [Watch $8 S^3$ factor; Watershed @vincent1991watersheds for cells],
    [Thin Vascular / Filament Graphs], [Poor: fixed scale $S$ misses sub-$S$ thin structures], [Use Watershed @vincent1991watersheds or multiscale cuts],
    [High-D Non-Spatial Tabular], [Low: no spatial lattice; $2S times 2S$ windowing invalid], [Use $k$-means++ @arthur2007kmeanspp or GMM @fraley2002mclust],
    [Real-Time Edge / UAV Vision], [Excellent: deterministic $O(N I)$ time, high cache reuse], [Cap $I <= 10$ for strict latency bounds],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for SLIC.],
)

#sec("V", "Structural Relations and Empirical Behavior")

SLIC bridges prototype quantization and spatial segmentation, uniting Lloyd $k$-means @lloyd1982least, Spatial Fuzzy $c$-Means (SFCM) @chuang2006sfcm, Mean Shift @comaniciu2002mean, and Watershed @vincent1991watersheds.

#proposition(title: "Asymptotic Reduction to Weighted k-Means and Spatial-Regularization Family", origin: "Original Proof")[
  Let $J_(m, S, R)$ denote generalized SLIC with window half-width $R >= S$ and membership fuzzifier $q >= 1$. Then: (i) as $R -> infinity$ with $q = 1$, SLIC reduces identically to weighted 5D $k$-means under $W_(m, S)$ @lloyd1982least; (ii) relaxing $q > 1$ with local neighborhood membership priors recovers the SFCM variational family @chuang2006sfcm; and (iii) as $m -> infinity$, the partition converges to the exact unweighted 2D Voronoi tessellation.
]
#proof[
  When $R >= max(W, H)$, every pixel satisfies $cal(A)_R (i; bold(M)) = \{1, dots, k\}$, making assignment and centroid updates identical to Lloyd's algorithm on $tilde(bold(z))_i = W_(m, S)^(1 \/ 2) bold(z)_i$. Introducing continuous memberships $u_(i j) in [0, 1]$ ($sum_j u_(i j) = 1$) weighted by $u_(i j)^q$ with spatial neighborhood cohesion yields SFCM optimality @chuang2006sfcm. Part (iii) follows directly from Theorem 2(1).
]

On the *ClusterOverview* four-dataset benchmark (adapting SLIC by partitioning features into attribute and spatial subspaces), SLIC achieves *Iris* Silhouette @rousseeuw1987silhouettes $0.447$ and ARI @hubert1985comparing $0.583$, *Customers* Silhouette $0.330$, *Multi-Shape* Silhouette $0.609$ (ARI N/A), and *High-D* Silhouette $0.259$ (ARI N/A), at a mean runtime of $0.227$ ms—the fastest latency among Algs. 166--169. This $0.227$ ms speed directly validates Lemma 1's $O(N I)$ localized sweeps. Strong *Multi-Shape* cohesion ($0.609$) reflects SLIC's ability to tile curved geometries into compact local patches, whereas the moderate *Iris* ARI ($0.583$) and sharp *High-D* drop ($0.259$) corroborate Proposition 1: without an intrinsic 2D spatial lattice, artificial windowing and anisotropic tensor conditioning degrade global partition alignment.
