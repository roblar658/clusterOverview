#import "../template/ieee.typ": *

#alg_entry(
  id: 41,
  name: "Quickshift",
  category: "Cat. II: Density & Mode-Seeking",
  year: "2008",
  paradigm: "Directed-Tree Nonparametric Mode-Seeking",
  complexity: $O(n^2 d)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Nonparametric mode-seeking algorithms define clusters as basins of attraction of local maxima (modes) of an empirical probability density estimator, liberating clustering from convexity assumptions and pre-specified cluster counts $k$. Continuous Mean Shift (Alg. 34 @cheng1995mean @comaniciu2002mean) ascends the density surface by iterating weighted averages in $bb(R)^d$, requiring $O(n^2 d T)$ operations and ambient vector-space geometry; Medoidshift (Alg. 42 @sheikh2007medoidshift) eliminates vector averaging by shifting strictly between dataset exemplars, but incurs $O(n^3)$ matrix-product overhead and can stall prematurely at sparse points. Introduced by Vedaldi and Soatto @vedaldi2008quickshift, *Quickshift* achieves a single-pass, non-iterative discrete mode-seeking procedure valid in arbitrary finite metric spaces.

Let $cal(X) = {x_1, dots, x_n}$ be a finite dataset embedded in a metric (or symmetric dissimilarity) space $(cal(M), d)$, typically $bb(R)^d$ with $d(x_i, x_j) = norm(x_i - x_j)_2$. Given a positive kernel bandwidth $sigma > 0$, Quickshift first evaluates a Parzen–Rosenblatt kernel density estimate (KDE) at every sample point:
$ hat(p)_sigma (x_i) = 1 / n sum_(j=1)^n K_sigma (d(x_i, x_j)), $
where for the canonical isotropic Gaussian kernel on $bb(R)^d$,
$ K_sigma (r) = (2 pi sigma^2)^(-d / 2) exp( - r^2 / (2 sigma^2) ). $
(Any symmetric radially non-increasing kernel $K_sigma$, including uniform, Epanechnikov, or domain-specific positive-definite Mercer kernels, may be substituted.) Rather than computing iterative gradient steps, Quickshift defines for each index $i in {1, dots, n}$ the *strictly higher-density upper set*
$ cal(U)_i := { j in {1, dots, n} : hat(p)_sigma (x_j) > hat(p)_sigma (x_i) }, $
and links $i$ directly to its *nearest neighbor in $cal(U)_i$*:
$ pi(i) := arg min_(j in cal(U)_i) d(x_i, x_j), quad delta_i := d(x_i, x_(pi(i))), $
with the convention $pi(i) = i$ and $delta_i = +infinity$ when $cal(U)_i = emptyset$ (i.e., at a global density maximizer). To guarantee full determinism under exact density or distance ties, one breaks density ties using the strict lexicographic order $(hat(p)_sigma (x_j), -j) > (hat(p)_sigma (x_i), -i)$ and distance ties by smallest index $j$.

#definition(title: "Quickshift Directed Forest and Scale-Cut Modes")[
  Given a maximum link-length threshold $tau in (0, infinity]$, the *pruned Quickshift parent map* $pi_tau : {1, dots, n} -> {1, dots, n}$ is defined by
  $ pi_tau (i) := cases( pi(i) & "if" cal(U)_i != emptyset space and space delta_i <= tau, i & "otherwise." ) $
  The induced directed graph $G_tau = (V, E_tau)$ with $V = {1, dots, n}$ and $E_tau = {(i, pi_tau (i)) : pi_tau (i) != i}$ is the *Quickshift forest* at scale $tau$. Its root set $cal(R)_tau = {i in V : pi_tau (i) = i}$ constitutes the discovered *density modes*, and the resulting cluster partition $cal(C)(tau)$ consists of the weakly connected components (directed in-trees) of $G_tau$.
]

Varying $tau$ slices the single unpruned tree $G_infinity$ into a complete multi-scale hierarchy of mode-seeking partitions without re-evaluating densities or pairwise distances (Algorithm 1).

#algbox(
  title: "Quickshift Mode-Seeking Algorithm",
  inputs: [Dataset $cal(X) = {x_1, dots, x_n}$ with metric $d$, bandwidth $sigma > 0$, link-cut threshold $tau > 0$ (or target cluster count $k >= 2$).],
  outputs: [Mode roots $cal(R)_tau$ and cluster partition $cal(C)(tau)$.],
)[
  1. *Density estimation:* For each $i in {1, dots, n}$, compute pairwise distances $d_(i j) = d(x_i, x_j)$ and KDE score $hat(p)_sigma (x_i) = sum_(j=1)^n exp(-d_(i j)^2 / (2 sigma^2))$.
  2. *Nearest higher-density linking:* For each $i in {1, dots, n}$:
  3. #h(1em) Find $cal(U)_i = {j : (hat(p)_sigma (x_j), -j) > (hat(p)_sigma (x_i), -i)}$.
  4. #h(1em) *If* $cal(U)_i = emptyset$, set $pi(i) <- i$ and $delta_i <- +infinity$.
  5. #h(1em) *Else*, set $pi(i) <- arg min_(j in cal(U)_i) d_(i j)$ and $delta_i <- d_(i, pi(i))$.
  6. *Scale pruning:* For each $i$, set $pi_tau (i) <- i$ if $delta_i > tau$, else $pi_tau (i) <- pi(i)$. (In fixed-$k$ mode, instead sever the $k - 1$ finite edges with largest $delta_i$.)
  7. *Root extraction:* For each $i in {1, dots, n}$, follow $u <- pi_tau (u)$ until $pi_tau (u) = u$ (using path compression) and assign label $L_i <- u$.
  8. *Return* roots $cal(R)_tau = {i : pi_tau (i) = i}$ and partition induced by $L$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five foundational results for Quickshift: (1) strict acyclicity, spanning in-tree topology, and exact hierarchical nesting across all cut thresholds $tau$; (2) exact equivalence between Quickshift roots and $tau$-local density modes together with an asymptotic gradient-alignment bound; (3) monotone-density invariance and a 7-point constructive counterexample proving an inescapable single-threshold $tau$ trade-off; and (4) statistical consistency of Quickshift mode-tree recovery.

#subsec("A", "Strict Acyclicity and Multi-Scale Hierarchical Nesting")

#theorem(title: "Spanning In-Tree Topology and Monotone Refinement", origin: "Original Proof")[
  Under strict lexicographic density ordering $succ$, let $i^* = arg max_i (hat(p)_sigma (x_i), -i)$ be the unique global density maximizer, and let $G_tau = (V, E_tau)$ be the Quickshift graph at threshold $tau in (0, infinity]$. Then:
  1. $G_infinity$ is a *directed acyclic spanning in-tree* (arborescence) rooted at $i^*$, with $|E_infinity| = n - 1$ and every directed path having length at most $n - 1$ along strictly ascending density.
  2. For every $tau in (0, infinity]$, the number of discovered modes equals exactly
     $ |cal(R)_tau| = 1 + |{i in V backslash {i^*} : delta_i > tau}|. $
  3. For any $0 < tau_1 <= tau_2 <= infinity$, the partition $cal(C)(tau_1)$ is a *strict hierarchical refinement* of $cal(C)(tau_2)$: every cluster $C in cal(C)(tau_1)$ satisfies $C subset.eq C'$ for a unique $C' in cal(C)(tau_2)$, and $|cal(R)_(tau_1)| >= |cal(R)_(tau_2)|$.
]

#proof[
  (1) Every directed edge $(i, pi(i)) in E_infinity$ satisfies $pi(i) in cal(U)_i$, so $(hat(p)_sigma (x_(pi(i))), -pi(i)) succ (hat(p)_sigma (x_i), -i)$. Suppose for contradiction that $G_infinity$ contained a directed cycle $v_0 -> v_1 -> dots -> v_m = v_0$ ($m >= 1$). By transitivity of the strict total order $succ$, $(hat(p)_sigma (x_(v_0)), -v_0) succ (hat(p)_sigma (x_(v_0)), -v_0)$, which is impossible. Furthermore, for every $i != i^*$, $i^* in cal(U)_i != emptyset$, so $i$ has out-degree exactly $1$, while $i^*$ has out-degree $0$. Thus $|E_infinity| = n - 1$. In any finite acyclic digraph where $n - 1$ vertices have out-degree $1$ and one vertex $i^*$ has out-degree $0$, following outgoing edges from any vertex cannot repeat vertices and must terminate at the unique sink $i^*$ in at most $n - 1$ steps; hence $G_infinity$ is a spanning in-tree rooted at $i^*$.

  (2) Pruning removes precisely the edge set $D(tau) = {(i, pi(i)) : i != i^*, space delta_i > tau}$ from the spanning tree $G_infinity$. Deleting $m = |D(tau)|$ edges from an $n$-vertex tree partitions it into exactly $m + 1$ connected tree components, each with a unique vertex of out-degree $0$ (its root), giving $|cal(R)_tau| = 1 + |D(tau)|$.

  (3) Whenever $tau_1 <= tau_2$, any edge satisfying $delta_i <= tau_1$ also satisfies $delta_i <= tau_2$, hence $E_(tau_1) subset.eq E_(tau_2)$ and $D(tau_2) subset.eq D(tau_1)$. Consequently every path in $G_(tau_1)$ persists in $G_(tau_2)$, so each connected component of $G_(tau_1)$ is contained within a single connected component of $G_(tau_2)$, and $|cal(R)_(tau_1)| = 1 + |D(tau_1)| >= 1 + |D(tau_2)| = |cal(R)_(tau_2)|$.
]

#subsec("B", "Local-Mode Exactness and Gradient-Ascent Alignment")

#proposition(title: "Exact Mode Recovery and Half-Space Ascent Bound", origin: "Original Proof")[
  Assume distinct empirical densities $hat(p)_sigma (x_i)$. Call $x_m in cal(X)$ a *$tau$-local empirical mode* if $hat(p)_sigma (x_m) >= hat(p)_sigma (x_j)$ for all $x_j in cal(X)$ with $d(x_m, x_j) <= tau$.
  1. $m in cal(R)_tau$ if and only if $x_m$ is a $tau$-local empirical mode of $hat(p)_sigma$.
  2. In $bb(R)^d$ with $hat(p)_sigma in C^2 (bb(R)^d)$ and uniform Hessian operator bound $L_H = sup_x norm(nabla^2 hat(p)_sigma (x))_"op" < infinity$, every non-root step $u_i = x_(pi(i)) - x_i$ (length $delta_i = norm(u_i)_2$) from a point with $g_i = nabla hat(p)_sigma (x_i) != bold(0)$ satisfies the directional cosine bound
     $ chevron.l g_i, u_i chevron.r / (norm(g_i)_2 norm(u_i)_2) > - (L_H delta_i) / (2 norm(g_i)_2). $
     In particular, as the local spacing $delta_i -> 0$, every Quickshift link lies in the ascending gradient half-space $chevron.l g_i, u_i chevron.r >= 0$.
]

#proof[
  (1) ($=>$) Suppose $m in cal(R)_tau$. Then either $cal(U)_m = emptyset$ (so $hat(p)_sigma (x_m)$ is globally maximal and trivially $tau$-locally maximal) or $delta_m = min_(j in cal(U)_m) d(x_m, x_j) > tau$. In the latter case, no point $x_j$ with $d(x_m, x_j) <= tau$ can belong to $cal(U)_m$, so $hat(p)_sigma (x_j) <= hat(p)_sigma (x_m)$ for all such $x_j$. ($arrow.l.double$) Conversely, if $hat(p)_sigma (x_m) >= hat(p)_sigma (x_j)$ for all $d(x_m, x_j) <= tau$, then every $j in cal(U)_m$ necessarily has $d(x_m, x_j) > tau$; hence either $cal(U)_m = emptyset$ or $delta_m > tau$, both yielding $pi_tau (m) = m$ and $m in cal(R)_tau$.

  (2) By Taylor's theorem with Lagrange remainder at some $xi_i$ on the segment $[x_i, x_(pi(i))]$,
  $ 0 &< hat(p)_sigma (x_(pi(i))) - hat(p)_sigma (x_i) \
      &= chevron.l g_i, u_i chevron.r + 1 / 2 u_i^T nabla^2 hat(p)_sigma (xi_i) u_i \
      &<= chevron.l g_i, u_i chevron.r + 1 / 2 L_H delta_i^2. $
  Rearranging gives $chevron.l g_i, u_i chevron.r > - 1 / 2 L_H delta_i^2$; dividing both sides by $norm(g_i)_2 delta_i > 0$ yields the stated directional cosine lower bound, which converges to $0$ as $delta_i -> 0$ at any non-stationary point ($norm(g_i)_2 > 0$).
]

#subsec("C", "Invariances and Single-Threshold Scale Trade-Off")

#theorem(title: "Monotone-Density Invariance & No-Valid-Threshold Counterexample", origin: "Original Counterexample")[
  1. *Monotone Density & Metric Invariance:* Quickshift is strictly invariant under any strictly increasing transformation $phi(hat(p)_sigma)$ of the density values (e.g., dropping normalizing constants or taking $log hat(p)_sigma$) and under all isometries of $(cal(M), d)$; it operates directly on arbitrary dissimilarity matrices without vector arithmetic.
  2. *Impossibility of Universal Single-$tau$ Recovery:* There exists a one-dimensional dataset $cal(X) subset.eq bb(R)$ of $n = 7$ points partitioning cleanly into two well-separated clusters $C_1, C_2$ (with Gaussian KDE bandwidth $sigma = 0.5$) such that *no threshold $tau > 0$* recovers the two-cluster partition ${C_1, C_2}$: every $tau$ either fragments $C_1$'s periphery into a spurious singleton mode or merges $C_1$'s mode into $C_2$.
]

#proof[
  (1) For any strictly increasing $phi$, $phi(hat(p)_sigma (x_j)) > phi(hat(p)_sigma (x_i)) <=> hat(p)_sigma (x_j) > hat(p)_sigma (x_i)$, leaving $cal(U)_i$, $pi(i)$, and $delta_i$ unchanged. Because $K_sigma (d(x_i, x_j))$ and $arg min_(j in cal(U)_i) d(x_i, x_j)$ depend exclusively on pairwise distances, any metric isometry commutes identically with Quickshift.

  (2) Fix $sigma = 0.5$ (so $2 sigma^2 = 0.5$) and define the 7-point 1D dataset consisting of Cluster $C_1 = {-3.0, 0.0, 0.2}$ (a small cluster with one sparse tail point) and Cluster $C_2 = {1.8, 2.0, 2.1, 2.2}$ (a denser adjacent cluster). Evaluating unnormalized Gaussian KDE scores $tilde(p)(x_i) = sum_(j=1)^7 exp(-(x_i - x_j)^2 / 0.5)$ yields:
  - $tilde(p)(-3.0) = 1 + e^(-18) + e^(-20.48) + dots < 1.0001$;
  - $tilde(p)(0.0) = 1 + e^(-0.08) + e^(-6.48) + dots approx 1.9247$;
  - $tilde(p)(0.2) = 1 + e^(-0.08) + e^(-5.12) + dots approx 1.9309$;
  - $tilde(p)(1.8) = 1 + e^(-0.08) + e^(-0.18) + e^(-0.32) + dots approx 3.490$.
  Observe the strict ordering $tilde(p)(-3.0) < tilde(p)(0.0) < tilde(p)(0.2) < tilde(p)(1.8)$ (and $tilde(p)(2.0), tilde(p)(2.1), tilde(p)(2.2) > 3.49$). Now compute the critical link lengths $delta_i$:
  - For the tail point $x_1 = -3.0$, its nearest strictly denser neighbor is $x_2 = 0.0$, giving $delta_1 = 3.0$.
  - For $C_1$'s local mode $x_3 = 0.2$, its nearest strictly denser neighbor lies in $C_2$ at $x_4 = 1.8$, giving $delta_3 = 1.8 - 0.2 = 1.6$.
  Since $delta_3 = 1.6 < 3.0 = delta_1$:
  - If $tau < 1.6$, edge $(x_1, x_2)$ is severed, isolating $x_1 = -3.0$ as a spurious third singleton cluster (over-fragmentation, $|cal(R)_tau| >= 3$).
  - If $tau >= 1.6$, edge $(x_3, x_4)$ is retained, merging all of ${0.0, 0.2}$ into $C_2$ (mode collapse, destroying $C_1$).
  Hence no single scalar $tau > 0$ produces ${C_1, C_2}$.
]

#subsec("D", "Statistical Consistency of Mode-Tree Estimation")

#theorem(title: "Consistent Recovery of Density Modes and Basins", origin: "Classical — Cited")[
  Let $X_1, dots, X_n$ be i.i.d. from a compactly supported $C^2$ density $p$ on $bb(R)^d$ possessing finitely many strict local maxima $m_1, dots, m_K$ with negative-definite Hessians and non-degenerate saddle separation. Under the Vedaldi–Soatto KDE tree construction @vedaldi2008quickshift and cluster-tree consistency theory @chaudhuri2014consistent @ariascastro2016meanshift, if the bandwidth $sigma_n -> 0$ satisfies $n sigma_n^(d + 2) / log n -> infinity$ and the cut threshold $tau_n -> 0$ satisfies $tau_n = omega((log n / n)^(1 / d))$ while remaining below the minimum basin radius $r_0$, then with probability tending to $1$ as $n -> infinity$, Quickshift recovers exactly $K$ modes converging to ${m_1, dots, m_K}$ inside high-density regions.
]

#proof(label: "Proof Synthesis (Cited)")[
  Standard uniform kernel-density and gradient convergence bounds @ariascastro2016meanshift ensure that under $n sigma_n^(d+2) / log n -> infinity$, both $norm(hat(p)_(sigma_n) - p)_infinity -> 0$ and $norm(nabla hat(p)_(sigma_n) - nabla p)_infinity -> 0$ in probability. Within any upper level set $L_lambda$ containing the modes $m_1, dots, m_K$, a covering argument on balls of radius $r_n = C (log n / n)^(1 / d) << tau_n$ @chaudhuri2014consistent guarantees that every non-mode sample $x_i$ with $norm(nabla p(x_i))_2 >= eta > 0$ finds a sample $x_j$ in its ascending cone at distance $<= r_n <= tau_n$ with $hat(p)_(sigma_n) (x_j) > hat(p)_(sigma_n) (x_i)$, ruling out spurious internal roots by Proposition 1(1). Conversely, near each true local maximum $m_c$, the sample maximizing $hat(p)_(sigma_n)$ within $B(m_c, r_0)$ has no strictly denser sample within distance $tau_n < r_0$, guaranteeing it survives as a unique mode root converging to $m_c$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

In full generality on an unstructured $n$-point dataset in $bb(R)^d$, Step 1 computes $binom(n, 2)$ pairwise distances and kernel sums in $Theta(n^2 d)$ time, Step 2 scans all $n^2$ pairs to find $pi(i)$ in $Theta(n^2)$ time, and Step 4 extracts labels via path compression in $O(n)$ time. Total runtime is therefore $Theta(n^2 d)$—strictly avoiding both the $T$-iteration loop of Mean Shift ($O(n^2 d T)$) and the cubic matrix multiplication of Medoidshift ($O(n^3)$). Memory is $O(n d)$ if rows of distances are streamed on the fly, or $O(n^2)$ if cached. Crucially, on structured 2D/3D image grids (Vedaldi and Soatto's original computer vision setting @vedaldi2008quickshift), truncating kernel support and parent search to a local pixel window of side $2 R + 1$ slashes complexity to $O(n R^2 d)$—linear in image pixel count $n$ and highly SIMD/GPU parallelizable because every point's density and parent pointer are completely independent!

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for Quickshift],
    [Time / Space Complexity], [$O(n^2 d)$ general; $O(n R^2 d)$ on windowed image grids; $O(n d)$ streaming or $O(n^2)$ cached space.],
    [Geometric Inductive Bias], [Arbitrary non-convex density basins linked by single-pass nearest-higher-density arborescence.],
    [Determinism], [Strictly deterministic given fixed tie-breaking; non-iterative (zero convergence tolerance needed).],
    [Invariances], [Invariant under metric isometries and any strictly increasing density transform $phi(hat(p)_sigma)$; scale-sensitive via $(sigma, tau)$.],
    [Outlier Behavior], [Far outliers become singleton roots ($delta_i > tau$) rather than distorting inlier modes; no native noise label.],
    [Recommended Companions], [DPC (Alg. 51) to filter outlier roots via $rho_i delta_i$; SLIC (Alg. 166) for compact regular superpixels; HDBSCAN (Alg. 32) for noise.],
  ),
  caption: [Summary of theoretical guarantees and operational properties of Quickshift.],
)

#sec("IV", "Engineering and Scientific Applications")

Because Quickshift combines nonparametric arbitrary-shape mode discovery with single-pass deterministic execution and zero vector-averaging requirements, it has become a workhorse primitive across vision, manifold analysis, and interactive multi-scale data exploration.

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Superpixel Image Over-Segmentation:* In its canonical deployment @vedaldi2008quickshift, each image pixel is mapped to a 5D color-spatial vector $(L, a, b, alpha x, alpha y)$. Windowed Quickshift ($O(n R^2)$) climbs color-spatial density modes while snapping naturally to sharp luminance boundaries—since crossing an edge drops KDE density or exceeds $tau$—serving as a leading nonparametric competitor to SLIC (Alg. 166 @achanta2012slic) and Watershed (Alg. 167 @vincent1991watersheds).
2. *Interactive Multi-Scale Mode Exploration:* By Theorem 1(3), constructing $G_infinity$ once lets analysts sweep the threshold $tau$ continuously across all scales in $O(n)$ time per cut, visualizing mode births and merges without repeating $O(n^2 d)$ density evaluations.
3. *Non-Euclidean & Kernel-Space Clustering:* Because Theorem 2(1) eliminates vector centroids, Quickshift clusters strings, biological sequences, probability histograms (via Jensen–Shannon or Wasserstein metrics), and positive-definite kernel embeddings directly from pairwise dissimilarity matrices.
4. *Hyperspectral & Remote-Sensing Region Growing:* Satellite scenes combine high spectral dimensionality with spatially coherent land-cover polygons; Quickshift's product spatial-spectral KDE avoids forcing convex ellipsoidal shapes onto meandering rivers, roads, and agricultural strips.
5. *Biomedical Microscopy & Flow Cytometry Gating:* Automated cell-nuclei segmentation and high-throughput cytometry populations exhibit skewed, crescent-shaped, and multi-modal distributions where Quickshift isolates distinct cellular phenotypes without assuming Gaussianity or pre-specifying $k$.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Quickshift's chief failure modes follow directly from our theory: (a) *spurious outlier singleton modes*—isolated noise points automatically satisfy $delta_i > tau$ and become roots (Theorem 2(2)) unless post-filtered by minimum cluster size or density threshold; (b) *bandwidth–threshold coupling*, where dense cores and sparse peripheries cannot share one $tau$ (Theorem 2(2)); (c) *high-dimensional KDE collapse*, where $n sigma^d$ variance blows up and distances concentrate @beyer1999nearest; and (d) *flat-plateau sensitivity*, where uniform-density regions depend heavily on arbitrary tie-breaking.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Why Quickshift Fits / Prerequisite], [Caution / Alternative],
    [Vision Superpixels], [Windowed $O(n R^2)$ speed; edge-aware], [Use SLIC (Alg. 166) for uniform size],
    [Non-Convex 2D/3D Data], [Captures arcs/crescents in one pass], [Filter small outlier roots],
    [Non-Euclidean Metrics], [Zero centroids; arbitrary metric $d$], [Pairwise $O(n^2)$ limits $n > 10^5$],
    [Noisy Outlier Streams], [Modes stay uncorrupted by far outliers], [Use HDBSCAN (Alg. 32) or DPC (Alg. 51)],
    [High-D Embeddings], [KDE variance explodes in $d >> 10$], [Use SNN-DPC (Alg. 52) or UMAP+HDBSCAN],
  ),
  caption: [Application suitability and engineering guidelines for Quickshift.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Quickshift occupies a pivotal structural bridge uniting continuous mode-seeking (Mean Shift, Alg. 34), exemplar mode-seeking (Medoidshift, Alg. 42), minimum spanning trees (SLINK, Alg. 56 @sibson1973slink @gower1969minimum), and Density Peaks Clustering (DPC, Alg. 51 @rodriguez2014dpc).

#proposition(title: "Exact Structural Isomorphism with Density Peaks Clustering", origin: "Original Proof")[
  Let $cal(X) = {x_1, dots, x_n}$ be equipped with symmetric distances $d_(i j)$ and identical local density estimates $rho_i := hat(p)_sigma (x_i)$ (strictly ordered). Define for each point $i$ its nearest higher-density neighbor $n_i = arg min_(j : rho_j > rho_i) d_(i j)$ and separation $delta_i = d_(i, n_i)$ (with $delta_(i^*) = max_j d_(i^*, j)$ or $+infinity$ at the global maximum $i^*$). Then:
  1. The directed parent graph $G_infinity$ constructed by Quickshift (2008 @vedaldi2008quickshift) is *graph-isomorphic* (edge-for-edge identical) to the assignment pointer forest of Density Peaks Clustering (DPC, 2014 @rodriguez2014dpc).
  2. Quickshift and DPC differ *solely in their root-selection criterion* in the $(rho_i, delta_i)$ decision plane: Quickshift selects the half-plane $cal(R)_tau^"QS" = {i : delta_i > tau}$, whereas DPC selects the upper-right quadrant $cal(R)^"DPC" = {i : rho_i > rho_"min" space and space delta_i > delta_"min"}$ (or top-$k$ in $gamma_i = rho_i delta_i$), pruning low-density outlier roots that Quickshift retains.
]

#proof[
  (1) By direct comparison of Algorithm 1 with Rodriguez and Laio @rodriguez2014dpc, once the local density values $rho_i = hat(p)_sigma (x_i)$ coincide, the upper set ${j : rho_j > rho_i}$ equals $cal(U)_i$ identically for every $i$. Consequently $n_i = pi(i)$ and $delta_i = d(x_i, x_(pi(i)))$ coincide at every non-maximal vertex $i != i^*$. In DPC, every non-center point is assigned sequentially in descending order of $rho_i$ to the same cluster as its nearest higher-density neighbor $n_i$; inductively, this assigns every vertex $i$ to the exact root reached by following the directed pointer chain $i -> pi(i) -> pi(pi(i)) -> dots$ until hitting a designated center.

  (2) Removing outgoing edges from a chosen root set $cal(R) subset.eq V$ yields identical clusters in both algorithms. Quickshift's rule $cal(R)_tau^"QS" = {i : delta_i > tau}$ depends purely on $delta_i$, admitting any sparse outlier with low $rho_i$ and large $delta_i$ as a singleton mode (matching Theorem 2(2)). DPC's additional lower bound $rho_i > rho_"min"$ (or ranking by $gamma_i = rho_i delta_i$) excludes the low-$rho_i$, high-$delta_i$ region, reclassifying isolated outliers into the halo/noise set rather than distinct cluster modes.
]

Compared to Medoidshift (Alg. 42), Quickshift trades exact local Fréchet-medoid minimization for direct nearest-higher-density jumps—cutting general complexity from $O(n^3)$ to $O(n^2 d)$ and preventing discrete mode-locking inside sparse basins at the cost of introducing the explicit scale parameter $tau$.

On the ClusterOverview 4-dataset empirical benchmark (evaluated with adaptive median $k$-distance bandwidth $sigma$ and $tau = 2.2 sigma$), Quickshift achieves *Iris* $"Sil" = 0.586$ and $"ARI" = 0.567$, *Customers* $"Sil" = 0.312$, *Multi-Shape* $"Sil" = 0.787$ ($"ARI" =$ N/A), and *High-D Subspace+Noise* $"Sil" = 0.254$ ($"ARI" =$ N/A), with a fast mean runtime of *1.073 ms*. The standout *Multi-Shape* Silhouette of $0.787$—substantially outperforming both ST-DBSCAN ($0.630$) and Medoidshift ($0.262$)—validates Quickshift's ability to trace curved, non-convex density manifolds in a single pass without mode-locking. Conversely, its drop to $0.254$ on *High-D Subspace+Noise* confirms both KDE variance degradation in higher dimensions @beyer1999nearest and Theorem 2(2)'s tendency to spawn singleton outlier roots when noise points exceed $tau$, while its $1.073$ ms runtime confirms an order-of-magnitude empirical speedup over cubic Medoidshift ($11.398$ ms).
