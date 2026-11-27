#import "../template/ieee.typ": *

#alg_entry(
  id: 51,
  name: "DPC: Clustering by Fast Search and Find of Density Peaks",
  category: "Cat. II: Density & Mode-Seeking",
  year: "2014",
  paradigm: "Decision-Graph Mode Detection & Monotone Parent Forest",
  complexity: $cal(O)(n^2)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Mode-seeking clustering methods such as Mean Shift @cheng1995mean @comaniciu2002mean and DENCLUE @hinneburg1998denclue identify clusters with basins of attraction of local maxima of a continuous kernel density estimate, requiring iterative gradient ascent for every data point. Rodriguez and Laio's Density Peaks Clustering (DPC) @rodriguez2014dpc replaces iterative optimization with a non-iterative geometric characterization built directly from the pairwise distance matrix $D = (d_(i j))_(i,j=1)^n$ on a finite dataset $cal(X) = {x_1, dots, x_n}$.

DPC rests on two structural axioms about cluster centers (density peaks): (i) centers are surrounded by neighbors of lower local density, and (ii) centers lie at a relatively large distance from *any* point possessing strictly higher local density. Unlike centroid-based models ($k$-means) or single-threshold density connectivity (DBSCAN @ester1996dbscan), DPC captures arbitrarily shaped non-convex clusters in a single non-iterative pass while relying exclusively on pairwise distances rather than Euclidean coordinates.

#definition(title: "Local Density, Higher-Density Separation, and Decision Graph")[
  Given a symmetric distance matrix $D = (d_(i j))_(i,j=1)^n$ with $d_(i i) = 0$ and a cutoff scale $d_c > 0$, define the *local density* $rho_i$ of point $i in {1, dots, n}$ via either the step kernel or the smooth Gaussian kernel @rodriguez2014dpc:
  $
    rho_i^"cut" = sum_(j != i) upright(I) [d_(i j) < d_c], \
    rho_i^"gau" = sum_(j != i) exp( - (d_(i j) / d_c)^2 ).
  $
  Let $pi: {1, dots, n} -> {1, dots, n}$ be a strict total order refining local density ($rho_i > rho_j => pi(i) < pi(j)$, breaking ties deterministically by index), with unique global density maximum $i^* = pi^(-1)(1)$. For every $i != i^*$, define the *higher-density separation* $delta_i$ and *nearest higher-density parent* $p(i)$ by
  $
    delta_i = min_(j : pi(j) < pi(i)) d_(i j), \
    p(i) = arg min_(j : pi(j) < pi(i)) d_(i j),
  $
  and set $delta_(i^*) = max_(1 <= j <= n) d_(i^* j)$ (or $max_(j != i^*) delta_j$). The *decision graph* is the 2D scatter plot ${(rho_i, delta_i)}_(i=1)^n$, and the *composite peak score* is $gamma_i = rho_i delta_i$.
]

Cluster centers $cal(S) = {c_1, dots, c_K}$ are identified as upper-right outliers on the decision graph (simultaneously large $rho_i$ and $delta_i$, commonly operationalized by selecting the top-$K$ points in $gamma_i$ or thresholding $rho_i > rho_"min"$ and $delta_i > delta_"min"$). Once $cal(S)$ is fixed, every remaining point is assigned in a single descending-density sweep to the same cluster as its nearest higher-density parent $p(i)$.

Optionally, DPC separates reliable cluster *cores* from peripheral *halo* (noise) points: for each cluster $C_k$, define its border region $B_k = {i in C_k : exists j in.not C_k "with" d_(i j) < d_c}$ and border density threshold $rho_b (k) = max_(i in B_k) rho_i$ (with $rho_b (k) = 0$ if $B_k = emptyset$). Points in $C_k$ satisfying $rho_i > rho_b (k)$ form the cluster core, while those with $rho_i <= rho_b (k)$ form the halo.

#algbox(
  title: "DPC: Fast Search and Find of Density Peaks",
  inputs: [Pairwise distance matrix $D in bb(R)_(>= 0)^(n times n)$, cutoff percentile $p in (0, 100)$ (typically $1$–$2\%$), cluster count $K$ (or thresholds $rho_"min", delta_"min"$).],
  outputs: [Cluster centers $cal(S) = {c_1, dots, c_K}$, labels $ell: {1, dots, n} -> {1, dots, K}$, and optional core/halo flags.],
)[
  1. *Select Cutoff Scale:* Set $d_c$ to the $p$-th percentile of off-diagonal entries ${d_(i j) : i < j}$.
  2. *Compute Local Densities:* Evaluate $rho_i$ (Gaussian or step kernel) for all $i = 1, dots, n$, and sort indices into descending density order $pi^(-1)(1), dots, pi^(-1)(n)$.
  3. *Compute Higher-Density Separation & Parents:* Set $delta_(i^*) = max_j d_(i^* j)$ for $i^* = pi^(-1)(1)$. For rank $r = 2, dots, n$ with $i = pi^(-1)(r)$, find $p(i) = arg min_(j : pi(j) < r) d_(i j)$ and set $delta_i = d_(i, p(i))$.
  4. *Identify Density Peaks:* Form decision graph ${(rho_i, delta_i)}$ and select $K$ centers $cal(S) = {c_1, dots, c_K}$ maximizing $gamma_i = rho_i delta_i$ (ensuring $i^* in cal(S)$); assign $ell(c_k) = k$.
  5. *Single-Pass Parent Assignment:* Traverse $r = 2, dots, n$ with $i = pi^(-1)(r)$; if $i in.not cal(S)$, assign $ell(i) <- ell(p(i))$.
  6. *Optional Halo Filtering:* Compute $rho_b (k)$ over inter-cluster boundary pairs within $d_c$, and mark $i in C_k$ with $rho_i <= rho_b (k)$ as halo.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We present four formal results: exact directed-forest well-posedness of single-pass assignment, a variational optimality theorem characterizing DPC as a minimum-weight monotone spanning forest, a finite-sample mode-separation recovery guarantee, and constructive counterexamples proving domino error propagation and multi-density peak suppression.

#subsec("A", "Acyclicity and Spanning-Forest Structure")

#theorem(title: "Directed Acyclicity and Order-Independent Assignment", origin: "Original Proof")[
  Let $T_0 = (V, E_0)$ be the directed graph on $V = {1, dots, n}$ with edges $E_0 = {(i, p(i)) : i in V without {i^*}}$. For any center set $cal(S) subset.eq V$ containing $i^*$, let $T_cal(S) = (V, E_cal(S))$ be obtained by deleting all outgoing edges from $cal(S)$. Then:
  1. $T_0$ is a spanning in-arborescence rooted at $i^*$ (acyclic, with unique directed path from every $i in V$ to $i^*$).
  2. $T_cal(S)$ is a disjoint forest of $|cal(S)| = K$ in-arborescences rooted at $c_1, dots, c_K$, and the single-pass assignment assigns every $i in V$ to the unique root $r(i) in cal(S)$ reachable from $i$ in $T_cal(S)$.
  3. The labeling $ell$ is invariant across all topological orderings of $T_cal(S)$ and completes in exactly $n - K$ pointer evaluations.
]
#proof[
  *(1)* Every directed edge $(i, p(i)) in E_0$ satisfies $pi(p(i)) < pi(i)$ by construction. Suppose for contradiction that $T_0$ contained a directed cycle $v_0 -> v_1 -> dots -> v_m = v_0$ ($m >= 1$); then $pi(v_0) > pi(v_1) > dots > pi(v_m) = pi(v_0)$, a contradiction. Because $T_0$ is acyclic, every vertex $i != i^*$ has out-degree exactly $1$, and $i^*$ has out-degree $0$, the unique path starting from any $i in V$ strictly decreases $pi$ at every hop and must terminate at the unique sink $i^*$ in at most $pi(i) - 1 < n$ steps. Thus $T_0$ is a spanning in-arborescence rooted at $i^*$.

  *(2)* Deleting the outgoing edges of $cal(S) without {i^*}$ turns every $c_k in cal(S)$ into a sink of out-degree $0$, while every $i in.not cal(S)$ retains out-degree $1$ with strictly decreasing rank $pi$. Hence every $i in V$ has a unique directed path to a first encountered center $r(i) in cal(S)$, partitioning $V$ into $K$ disjoint rooted subtrees.

  *(3)* In any topological ordering of $T_cal(S)$ (in particular, sorting by ascending $pi$, i.e., descending $rho$), every parent $p(i)$ is processed strictly before $i$ since $pi(p(i)) < pi(i)$. By strong induction on $pi(i)$, $ell(i) = ell(p(i)) = ell(r(p(i))) = ell(r(i))$, independent of tie-breaking among incomparable vertices. Exactly one lookup is executed for each of the $n - K$ non-center points.
]

#subsec("B", "Variational Characterization: Minimum Monotone Spanning Tree")

#theorem(title: "Variational Optimality Over Density-Monotone Spanning Forests", origin: "Original Derivation")[
  Let $cal(F)_(pi, cal(S))$ denote the set of all directed spanning forests $F = (V, E_F)$ rooted at $cal(S)$ (with $i^* in cal(S)$) such that every non-root vertex $i in V without cal(S)$ has a unique outgoing edge $(i, q(i)) in E_F$ obeying density ascent $pi(q(i)) < pi(i)$. Define the total linkage cost $W(F) = sum_(i in V without cal(S)) d_(i, q(i))$. Then DPC's assignment forest $T_cal(S)$ is an exact global minimizer of $W$ over $cal(F)_(pi, cal(S))$:
  $
    W(T_cal(S)) = min_(F in cal(F)_(pi, cal(S))) W(F) = sum_(i in V without cal(S)) delta_i.
  $
]
#proof[
  First observe that *any* map $q: V without cal(S) -> V$ satisfying $pi(q(i)) < pi(i)$ for all $i in.not cal(S)$ automatically defines a cycle-free spanning forest rooted at $cal(S)$: no directed cycle can exist because $pi$ strictly decreases along every edge, every non-root has out-degree $1$, and every sink lies in $cal(S)$ (since $i^* in cal(S)$ is the only vertex with empty higher-rank set). Consequently, the global structural constraint $F in cal(F)_(pi, cal(S))$ factorizes completely into $n - K$ independent vertex constraints $q(i) in {j in V : pi(j) < pi(i)}$. Therefore,
  $
    min_(F in cal(F)_(pi, cal(S))) W(F) &= sum_(i in V without cal(S)) min_(j : pi(j) < pi(i)) d_(i j) \
    &= sum_(i in V without cal(S)) delta_i = W(T_cal(S)).
  $
  In particular, for $cal(S) = {i^*}$, the full parent tree $T_0$ is the minimum-weight density-monotone spanning arborescence, and choosing $cal(S)$ to maximize $sum_(c in cal(S) without {i^*}) delta_c$ is equivalent to cutting the $K - 1$ longest edges of $T_0$ to minimize residual tree weight $W(T_cal(S))$.
]

Theorem 2 uncovers the exact discrete optimization problem solved by DPC: unlike Kruskal/Prim's unconstrained minimum spanning tree (single-linkage clustering), DPC solves a *density-directed* minimum spanning forest problem where edge orientation is constrained by the empirical density hierarchy $pi$.

#subsec("C", "Finite-Sample Island Separation and Mode Consistency")

#theorem(title: "Exact Island Recovery and Decision-Graph Margin", origin: "Classical — Cited")[
  Let $cal(X) = C_1 union dots union C_K$ be partitioned into $K$ disjoint islands, and let $c_k = arg min_(i in C_k) pi(i)$ be the unique local density peak of $C_k$. Suppose there exist radii $0 < r_"in" < r_"out"$ such that:
  1. *(Intra-Island Ascent)* Every non-peak $i in C_k without {c_k}$ possesses at least one strictly denser point $j in C_k$ ($pi(j) < pi(i)$) with $d_(i j) <= r_"in"$.
  2. *(Inter-Island Separation)* For all $k != m$, $"dist"(C_k, C_m) >= r_"out"$.
  Then: (a) every non-peak $i in.not cal(S) = {c_1, dots, c_K}$ satisfies $delta_i <= r_"in"$ and $p(i) in C_k$; (b) every island peak $c_k$ satisfies $delta_(c_k) >= r_"out"$ (assuming $r_"out" <= max_j d_(i^* j)$); and (c) thresholding $delta_i in (r_"in", r_"out")$ identifies $cal(S)$ uniquely and recovers ${C_1, dots, C_K}$ with zero error @rodriguez2014dpc @chaudhuri2014consistent.
]
#proof[
  Fix any island $C_k$ and any non-peak point $i in C_k without {c_k}$. By (1), there exists $j in C_k$ with $pi(j) < pi(i)$ and $d_(i j) <= r_"in"$. Hence $delta_i = min_(j' : pi(j') < pi(i)) d_(i j') <= r_"in"$. Meanwhile, for every point $u in.not C_k$, condition (2) gives $d_(i u) >= r_"out" > r_"in" >= delta_i$, so the minimizer $p(i)$ must lie strictly inside $C_k$.

  Now consider any island peak $c_k$ other than the global maximum $i^*$. Since $c_k$ achieves minimum rank $pi$ within $C_k$, every point $j$ with $pi(j) < pi(c_k)$ must belong to a different island $C_m$ ($m != k$). By (2), $d_(c_k, j) >= r_"out"$ for all such $j$, whence $delta_(c_k) >= r_"out"$. For $i^* = c_(k^*)$, $delta_(i^*) = max_j d_(i^* j) >= "dist"(C_(k^*), C_m) >= r_"out"$ whenever $K >= 2$. Thus $cal(S) = {i : delta_i > tau}$ for any $tau in (r_"in", r_"out")$. Finally, because $p(i) in C_k$ for all $i in C_k without {c_k}$, every directed path in $T_cal(S)$ starting in $C_k$ remains inside $C_k$ and terminates at $c_k$, recovering ${C_1, dots, C_K}$ identically.
]

#subsec("D", "Domino Error Cascade and Multi-Density Peak Suppression")

#proposition(title: "Domino Propagation and Multi-Density Suppression of Gamma Selection", origin: "Original Counterexample")[
  1. *(Domino Cascade)* In DPC's assignment tree $T_cal(S)$, a single distance perturbation altering the parent pointer $p(u)$ of one non-center vertex $u$ misclassifies its entire descendant subtree $T_u$, which can have cardinality $Omega(n)$.
  2. *(Multi-Density Suppression)* Under the composite score $gamma_i = rho_i delta_i$ with fixed cutoff $d_c > 0$, there exists a well-separated two-cluster 1D mixture $cal(X) = C_1 union C_2$ (with $C_1$ dense and $C_2$ sparse) such that a secondary shoulder point $s in C_1$ attains $gamma_s > gamma_(c_2)$, causing top-$2$ selection on $gamma_i$ to split $C_1$ and miss $C_2$ completely.
]
#proof[
  *(1)* Consider a 1D chain of $m$ points $u_1, dots, u_m$ with strictly decreasing densities $rho_(u_1) > dots > rho_(u_m)$ spaced by $h$ so $p(u_(t+1)) = u_t$ for $t = 1, dots, m - 1$. All $m$ points lie in the subtree $T_(u_1)$ rooted at $u_1$. If $u_1$ lies near the saddle between centers $c_1$ and $c_2$ with $d_(u_1, c_1) = d_(u_1, c_2) + eta$ for arbitrarily small $eta > 0$, then $u_1$ attaches to $c_2$, and by Theorem 1(2) all $m = Omega(n)$ descendants $u_1, dots, u_m$ inherit label $c_2$ simultaneously—no local voting dampens the single-edge error.

  *(2)* Fix $d_c = 1$ and the step kernel $rho_i^"cut" = sum_(j != i) upright(I) [d_(i j) < 1]$. Construct a single connected dense cluster $C_1$ from three collinear sub-clumps (for integer $M >= 30$): $M$ points at $x = 0$ (left shoulder $A$), $M - 3$ bridge points at $x = 0.75$ (center $B$), and $M - 1$ points at $x = 1.5$ (right shoulder $S$). Because $d(0, 0.75) = d(0.75, 1.5) = 0.75 < 1 < 1.5 = d(0, 1.5)$, every point in $B$ is within $d_c$ of all points in $C_1$, yielding exact local densities
  $
    rho(B) &= 3 M - 5, \
    rho(A) &= 2 M - 4, \
    rho(S) &= 2 M - 5,
  $
  so $rho(B) > rho(A) > rho(S)$. Every shoulder point $a in A$ has nearest strictly denser neighbor in $B$ at distance $delta_a = 0.75$, giving $gamma_a = 0.75 (2 M - 4) = 1.5 M - 3$. Now place a well-separated sparse cluster $C_2$ of $5$ points within radius $0.2$ of $x = 10$ (with its peak $c_2$ at $x = 10$), yielding $rho_(c_2) = 4$ and $delta_(c_2) = 10 - 1.5 = 8.5$, hence $gamma_(c_2) = 4 dot 8.5 = 34$. Setting $M = 30$ yields
  $
    gamma_a = 42 > 34 = gamma_(c_2),
  $
  even though $C_2$ is separated by $8.5 d_c$ while $a in A$ lies inside $C_1$'s single connected density ball. Consequently, top-$2$ ranking on $gamma_i$ selects both $b in B$ and $a in A$ inside $C_1$, completely missing $C_2$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Given raw coordinates $cal(X) subset bb(R)^d$, forming the symmetric distance matrix $D$ requires $cal(O)(n^2 d)$ time and $cal(O)(n^2)$ memory. Computing $rho_i$ for all $n$ points via vectorized kernel evaluation takes $cal(O)(n^2)$ operations; sorting $rho$ requires $cal(O)(n log n)$; scanning higher-density neighbors to compute $(delta_i, p(i))$ requires $sum_(r=2)^n (r - 1) = n (n - 1) / 2 = cal(O)(n^2)$ comparisons; selecting $K$ peaks and executing single-pass assignment takes $cal(O)(n)$ time.

Critically, DPC performs *zero iterative epochs* (unlike Lloyd's $k$-means, EM, or Mean Shift @comaniciu2002mean) and requires only cache-friendly BLAS/vectorized distance sweeps, explaining its exceptionally low empirical latency on moderate datasets ($n <= 10^4$). Its primary scalability bottleneck is the $cal(O)(n^2)$ space cost of storing $D$, which exceeds $32$ GB at $n approx 6.5 times 10^4$ in double precision unless mitigated by $k$-NN graphs or spatial hashing.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$cal(O)(n^2 d)$ from coordinates; $cal(O)(n^2)$ from precomputed distance matrix $D$.],
    [Space Complexity], [$cal(O)(n^2)$ for full distance matrix $D$; $cal(O)(n)$ for $(rho, delta, p, ell)$ vectors.],
    [Inductive Bias], [Unimodal density basins (possibly non-convex) with separated peaks; non-iterative tree propagation.],
    [Hyperparameters], [Cutoff distance $d_c$ (typically $1$–$2\%$ distance quantile) and peak count $K$ (or $(rho_"min", delta_"min")$).],
    [Metric Generality], [Operates on arbitrary symmetric dissimilarity matrices (no vector-space mean required).],
    [Main Fragilities], [Multi-density peak suppression (Prop. 4), domino misassignment cascades, and $cal(O)(n^2)$ memory wall.],
  ),
  caption: [Operational and complexity profile of DPC.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Molecular Dynamics Conformational State Discovery:* Protein-folding trajectories inhabit rugged free-energy landscapes whose metastable states are high-density basins separated by low-density transition barriers @rodriguez2014dpc. Because DPC operates directly on RMSD or dihedral-angle dissimilarity matrices without Euclidean averaging (which yields unphysical conformations), it identifies true representative medoid conformations and separates transition-state halos.
2. *Face Recognition and Image-Patch Categorization:* On the Olivetti face database benchmark @rodriguez2014dpc, pairwise visual similarity matrices are non-Euclidean and classes vary strongly in pose manifolds. DPC selects canonical exemplar faces directly from the decision graph and assigns variations along geodesic density paths, avoiding centroid blurring.
3. *Single-Cell Transcriptomics & Flow Cytometry Phenotyping:* After nonlinear embedding (e.g., UMAP/t-SNE or diffusion maps), distinct cell phenotypes appear as compact density peaks with curved differentiation trajectories; DPC rapidly extracts major lineages and flags ambiguous intermediate states via halo filtering.
4. *Wireless Sensor Network & Spatial Event Localization:* Fast single-pass latency ($< 1$ ms on $n approx 10^3$) allows DPC to detect localized acoustic, thermal, or RF emitter hotspots directly from pairwise signal-strength dissimilarities without iterative message-passing convergence.
5. *Power-Grid Load-Curve & Financial Regime Prototypes:* Time-series clusters under dynamic time warping or correlation distance lack closed-form means; DPC directly takes the precomputed dissimilarity matrix and selects actual historical regime days as cluster peaks.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Five failure mechanisms characterize DPC: (i) *multi-density suppression* (Proposition 4), where global $d_c$ inflates dense-cluster $rho_i$ so shoulders outrank sparse peaks; (ii) *multi-modal elongated manifolds*, where a flat ring or wavy ribbon contains multiple local density ripples, causing DPC to split one manifold across several peaks; (iii) *domino error cascades* (Proposition 4), where one erroneous edge across a saddle redirects an entire subtree; (iv) *high-dimensional distance concentration* @beyer1999nearest, flattening both $rho_i$ and $delta_i$; and (v) *quadratic memory exhaustion* above $n approx 5 times 10^4$.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Conformational MD / Non-Euclidean Dissimilarities], [Ideal: metric-free exemplar peaks plus transition halo filtering.], [Use Gaussian kernel for small $n$ to avoid ties.],
    [Well-Separated Unimodal Non-Convex Clusters], [Strong: fast single-pass recovery (Theorem 3).], [Verify single peak per manifold on decision graph.],
    [Multi-Density Mixtures], [Poor: global $d_c$ suppresses sparse peaks (Prop. 4).], [Use SNN-DPC (Alg. 52, @liu2018snndpc) or HDBSCAN.],
    [Large-Scale Data ($n > 5 times 10^4$)], [Blocked by $cal(O)(n^2)$ distance matrix RAM.], [Use MR-DBSCAN (Alg. 49) or approximate $k$-NN DPC.],
    [High-Dimensional Noisy Spaces], [Weak: $d_c$ kernel and Euclidean $delta_i$ lose contrast.], [Use SNN-DPC @liu2018snndpc or spectral embeddings.],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for DPC.],
)

#sec("V", "Structural Relations and Empirical Behavior")

#proposition(title: "Exact Algorithmic Equivalence Between DPC and Unbounded QuickShift", origin: "Original Derivation")[
  Let QuickShift @vedaldi2008quickshift be executed on $cal(X)$ using the exact Parzen density estimate $hat(f)(x_i) prop rho_i$ and unbounded mode-jump threshold $tau = infinity$, linking each non-maximal point $i$ to $p_"QS" (i) = arg min_(j : hat(f)(x_j) > hat(f)(x_i)) d_(i j)$. Then QuickShift constructs the exact spanning arborescence $T_0$ of DPC. Consequently, DPC with center set $cal(S)$ is identical to QuickShift with an adaptive edge-cut set consisting of the outgoing edges of $cal(S) without {i^*}$ rather than a uniform length threshold $tau$ (which recovers DPC when $cal(S) = {i : delta_i > tau}$).
]
#proof[
  Since $hat(f)(x_i) = c dot rho_i$ for a positive normalization constant $c > 0$, the strict density comparison $hat(f)(x_j) > hat(f)(x_i)$ coincides identically with $rho_j > rho_i$. With $tau = infinity$, QuickShift imposes no finite cutoff on $d_(i j)$ and sets every non-maximal vertex's parent to $arg min_(j : rho_j > rho_i) d_(i j) = p(i)$, yielding the identical directed tree $T_0$ (under identical tie-breaking $pi$). In standard QuickShift, pruning edges of length $> tau$ produces roots ${i : delta_i > tau}$; DPC generalizes this threshold rule to joint $(rho_i, delta_i)$ selection on the decision graph while keeping the underlying parent forest $T_0$ unchanged.
]

DPC also contrasts cleanly with Mean Shift @comaniciu2002mean @ariascastro2016meanshift and Medoid Shift @sheikh2007medoidshift (which iterate continuous or discrete mode shifts until convergence, whereas DPC takes discrete nearest-higher-density jumps in one pass), Affinity Propagation @frey2007affinity (which selects exemplars via message passing optimizing sum of similarities rather than density peaks), and SNN-DPC (Alg. 52, @liu2018snndpc), which replaces DPC's global $d_c$ and single-pass tree with shared-nearest-neighbor density and two-step assignment.

*Empirical Benchmark Analysis:* On the ClusterOverview four-dataset benchmark (Silhouette @rousseeuw1987silhouettes and Adjusted Rand Index @hubert1985comparing), DPC achieves:
- *Iris:* Silhouette $= 0.278$, $"ARI" = 0.569$;
- *Customers:* Silhouette $= 0.247$;
- *Multi-Shape:* Silhouette $= 0.718$ ($"ARI" =$ N/A);
- *High-D:* Silhouette $= 0.296$ ($"ARI" =$ N/A);
- *Mean Runtime:* $0.888$ ms.

These figures illuminate every theoretical strength and weakness proved in Section II. First, DPC is by far the fastest algorithm of Algs. 49–52 ($0.888$ ms mean runtime—$2.4 times$ faster than SNN-DPC, $4.7 times$ faster than VDBSCAN, and $10.6 times$ faster than MR-DBSCAN), directly reflecting its non-iterative vectorized single-pass design. Second, on *Multi-Shape*, DPC attains a strong Silhouette of $0.718$, confirming Theorem 3: monotone parent-tree assignment follows non-convex manifolds effectively whenever each shape has a clean peak. Third, *Iris* displays an instructive divergence: DPC achieves the *highest* $"ARI"$ ($0.569$) among all four algorithms despite a low Silhouette ($0.278$), because parent-pointer assignment separates ground-truth species along density basins even where Versicolor and Virginica overlap closely in Euclidean margin. Conversely, *Customers* ($0.247$) and *High-D* ($0.296$) expose Proposition 4's multi-density suppression and Euclidean distance concentration under a single global $d_c$, directly motivating SNN-DPC (Alg. 52).
