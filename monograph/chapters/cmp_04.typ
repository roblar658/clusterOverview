#import "../template/ieee.typ": *

#cmp_entry(
  roman: "IV",
  title: "Intra-Category Comparative Analysis: Density and Mode-Seeking (Algs. 31-55)",
  scope_label: "Head-to-Head Theoretical & Practical Comparison of Algorithms 31-55",
  provenance: "Original Unifying Proofs + Cited Classical Theorems",
)

The twenty-five algorithms of Category II (Algs. 31--55) abandon the parametric distortion functionals and convex Voronoi partitions of Part I in favor of a nonparametric topological premise: clusters correspond either to maximal connected components of empirical upper density level sets $L_f (lambda) = {x in bb(R)^d : f(x) >= lambda}$ or to basins of attraction of local modes of a smooth or discrete density field $hat(f)$. Freed from pre-specified cluster counts $k$ and spherical symmetry assumptions, density-based and mode-seeking methods recover non-convex manifolds, nested filaments, and arbitrary noise backgrounds @ester1996dbscan @comaniciu2002mean @jain2010data @aggarwal2013data. Yet this expressive power introduces four fundamental mathematical fault lines that partition Category II into distinct algorithmic families:
1. *Single-Scale vs. Multi-Scale Density Connectivity:* How an algorithm responds when true clusters span distinct local intensities $lambda_1 >> lambda_2$ or reside in high-dimensional spaces where $L_2$ distances concentrate @beyer1999nearest @schubert2017dbscanrevisited.
2. *Continuous Gradient Flows vs. Discrete Mode-Seeking:* Whether mode discovery proceeds via fixed-step ascent, self-normalizing Expectation--Maximization (EM) fixed-point maps, dataset-contracting blurring operators, Cartesian grid quantization, or non-iterative discrete graphs on $cal(X)$.
3. *Level-Set Tree Consistency vs. Decision-Graph Heuristics:* Whether modal hierarchy and outlier scoring emerge from provably consistent filtrations of Hartigan's cluster tree @hartigan1981consistency @chaudhuri2014consistent or from geometric parent-pointer decision graphs @rodriguez2014dpc @liu2018snndpc.
4. *Exact Distributed Preservation vs. Subquadratic Approximation:* How quadratic pairwise query complexity $Theta(n^2 d)$ is broken across parallel clusters, axis-aligned grids, Voronoi decompositions, or subsampled core covers without destroying global topological connectivity.

This comparative chapter provides a unified mathematical dissection of all twenty-five methods in Category II across four structured analytical sections—(I) density-reachability and variable-density extensions (Algs. 31--33, 36, 38--40, 45, 50); (II) kernel mode-seeking and gradient-ascent flows (Algs. 34--35, 41--43, 54--55); (III) density peaks, level-set trees, and local outlier scoring (Algs. 46, 51--53); and (IV) parallel, grid, and subsampled scalability architectures (Algs. 37, 44, 47--49)—culminating in (V) a complete 25-algorithm head-to-head comparative matrix and mathematically grounded engineering decision guide.

// =============================================================================
// SECTION I
// =============================================================================
#sec("I", "Density-Reachability and Variable-Density Families")

The foundational lineage of Category II stems from DBSCAN (Alg. 31, @ester1996dbscan), which operationalizes level-set connectivity at a fixed metric radius $epsilon > 0$ and cardinality threshold $m = "MinPts" in bb(N)_(>= 2)$. Eight descendants in Part II attack orthogonal limitations of this single-scale paradigm: HDBSCAN (Alg. 32, @campello2013hdbscan @campello2015hierarchical), OPTICS (Alg. 33, @ankerst1999optics), and VDBSCAN (Alg. 50, @liu2007vdbscan) generalize the *density filtration*; SNN (Alg. 36, @ertoz2003snn) and RNN-DBSCAN (Alg. 45, @bryant2018rnndbscan) replace metric radii with *rank-adaptive neighborhood topology*; GDBSCAN (Alg. 39, @sander1998gdbscan) and ST-DBSCAN (Alg. 40, @birant2007stdbscan) generalize *predicates, weights, and space--time domains*; and DBCLASD (Alg. 38, @xu1998dbclasd) replaces geometric thresholds with *nonparametric Poisson goodness-of-fit testing*.

#subsec("A", "Impossibility of Fixed-Scale Recovery Under Multi-Density Mixtures")

A recurring empirical critique of DBSCAN @schubert2017dbscanrevisited is its inability to coexist with clusters of disparate densities. We begin by proving a sharp geometric impossibility boundary that rigorously explains *why* every single-radius level-set algorithm fails simultaneously on dense adjacent clusters and sparse distant clusters.

#definition(title: "Intra-Cluster Core Scale and Inter-Cluster Gap")[
  For any finite subset $C subset.eq cal(X)$ in a metric space $(cal(M), d)$ with $|C| >= m >= 2$, let $d_m^((C)) (x)$ denote the distance from $x in C$ to its $m$-th nearest neighbor inside $C$ (counting $x$ itself). Define the *upper* and *lower* internal $m$-core thresholds of $C$ by
  $
    overline(tau)_m (C) := max_(x in C) d_m^((C)) (x), \
    underline(tau)_m (C) := min_(x in C) d_m^((C)) (x),
  $
  and for disjoint $A, B subset.eq cal(X)$ denote the inter-cluster gap by $"dist"(A, B) := min_(x in A, y in B) d(x, y)$.
]

#theorem(title: "Multi-Scale Separation Impossibility for Fixed-Radius DBSCAN", origin: "Original Proof")[
  Let $cal(X) = C_1 union.sq C_2 union.sq C_3$ consist of three disjoint clusters in $(cal(M), d)$ with $|C_r| >= m >= 2$ for $r in {1, 2, 3}$. Write $g_(12) := "dist"(C_1, C_2)$ and assume that $C_3$ is isolated from $C_1 union C_2$ by $"dist"(C_3, C_1 union C_2) >= underline(tau)_m (C_3)$. If the dense-cluster core scales, inter-cluster gap, and sparse-cluster core scale satisfy
  $
    max { overline(tau)_m (C_1), space overline(tau)_m (C_2) } < g_(12) < underline(tau)_m (C_3),
  $
  then there exists *no* radius $epsilon > 0$ for which DBSCAN$(cal(X), epsilon, m)$ recovers all three clusters as distinct non-empty components: specifically,
  1. for every $epsilon < underline(tau)_m (C_3)$, the entire sparse cluster collapses to noise ($C_3 subset.eq cal(X)_"noise"$);
  2. for every $epsilon >= underline(tau)_m (C_3)$, clusters $C_1$ and $C_2$ are deterministically fused into a single density-connected component.
]

#proof[
  Consider the exhaustive dichotomy on $epsilon > 0$.

  *Case 1 ($epsilon < underline(tau)_m (C_3)$):* Take an arbitrary $z in C_3$. Because $"dist"(C_3, C_1 union C_2) >= underline(tau)_m (C_3) > epsilon$, the closed ball $B_epsilon (z)$ contains zero points from $C_1 union C_2$, so $N_epsilon (z) = B_epsilon (z) inter C_3$. By definition of $underline(tau)_m (C_3) <= d_m^((C_3)) (z)$, strictly fewer than $m$ points of $C_3$ lie within distance $epsilon$ of $z$, yielding $|N_epsilon (z)| < m$. Hence $C_3 inter cal(X)_"core" = emptyset$. Moreover, since $d(z, x) > epsilon$ for every $x in C_1 union C_2$, no point of $C_3$ can be a border point of any core in $C_1 union C_2$. Thus $C_3 subset.eq cal(X)_"noise"$.

  *Case 2 ($epsilon >= underline(tau)_m (C_3)$):* By hypothesis, $epsilon >= underline(tau)_m (C_3) > g_(12) > max{overline(tau)_m (C_1), overline(tau)_m (C_2)}$. Consequently, for *every* $x in C_1$, $|N_epsilon (x)| >= |B_epsilon (x) inter C_1| >= m$ because $epsilon > overline(tau)_m (C_1) >= d_m^((C_1)) (x)$; symmetrically, every $y in C_2$ satisfies $|N_epsilon (y)| >= m$. Thus $C_1 union C_2 subset.eq cal(X)_"core"$. Let $(x_1^*, x_2^*) in C_1 times C_2$ realize $d(x_1^*, x_2^*) = g_(12)$. Because $x_1^*, x_2^* in cal(X)_"core"$ and $d(x_1^*, x_2^*) = g_(12) < epsilon$, $x_1^*$ and $x_2^*$ are mutually directly density-reachable. By maximality of DBSCAN clusters, every component intersecting $C_1$ containing $x_1^*$ merges with the component intersecting $C_2$ containing $x_2^*$, preventing $C_1$ and $C_2$ from remaining disjoint.
]

In Euclidean space $bb(R)^d$, if $C_1, C_2$ are drawn from homogeneous Poisson intensities $lambda_1$ and $C_3$ from $lambda_3$, then $d_m^((C_r)) asymp (m / (lambda_r V_d))^(1/d)$, where $V_d = pi^(d/2) / Gamma(d/2 + 1)$. Theorem 1 therefore triggers whenever the density contrast satisfies
$
  lambda_1 / lambda_3 > (g_(12) / (overline(tau)_m (C_1)))^d,
$
showing that even modest geometric gaps force fixed-scale DBSCAN into a fatal choice between *over-merging* dense clusters and *erasing* sparse clusters.

#subsec("B", "Nested Level-Set Recovery: HDBSCAN, OPTICS, and VDBSCAN")

HDBSCAN (Alg. 32), OPTICS (Alg. 33), and VDBSCAN (Alg. 50) each escape Theorem 1 by evaluating density connectivity across *multiple* radii, yet they induce mathematically distinct topological objects. Let $d_"core" (x) := d_m (x)$ on $cal(X)$, let $G_("DB"^*) (epsilon) := (V_"core" (epsilon), E_"core" (epsilon))$ denote the DBSCAN\* core graph (omitting border points, @campello2015hierarchical) with $V_"core" (epsilon) = {x : d_"core" (x) <= epsilon}$ and ${u, v} in E_"core" (epsilon) <=> d(u, v) <= epsilon$, and recall HDBSCAN's symmetric mutual-reachability metric
$
  d_"mreach" (u, v) = max { d_"core" (u), space d_"core" (v), space d(u, v) }, quad u != v,
$
contrasted with OPTICS's *asymmetric* reachability distance $d_"reach" (q | p) = max{d_"core" (p), d(p, q)}$.

#theorem(title: "Unified Filtration Equivalence of DBSCAN*, HDBSCAN, and OPTICS", origin: "Original Unified Proof")[
  Fix a finite metric space $(cal(X), d)$ and integer $m >= 2$. Let $cal(T)_"MST"$ be any Minimum Spanning Tree of the complete graph $(cal(X), d_"mreach")$, and let $(pi, c, r)$ be any OPTICS ordering with $epsilon_"max" = infinity$.
  1. *(Exact HDBSCAN--DBSCAN\* Equivalence):* For every scale $epsilon > 0$, the connected components of $G_("DB"^*) (epsilon)$ coincide identically with the non-singleton connected components (and isolated vertices with $d_"core" (x) <= epsilon$) of the forest $cal(T)_"MST" (epsilon)$ obtained by deleting all edges of $cal(T)_"MST"$ with weight strictly exceeding $epsilon$.
  2. *(Contiguous Core Embedding & Border Asymmetry in OPTICS):* Every connected component $C$ of $G_("DB"^*) (epsilon)$ is visited as a contiguous interval modulo border points in $pi$: once the first vertex $o_(t_0) in C$ is processed, every remaining vertex of $C$ is output with reachability $r_t <= epsilon$ before any vertex $z in.not C$ having $min_(u in C) d_"reach" (z | u) > epsilon$ is popped from the seed queue. However, OPTICS admits non-core border points $b$ ($d_"core" (b) > epsilon$) with $r(b) <= epsilon$, whereas $cal(T)_"MST" (epsilon)$ strictly classifies every such $b$ as isolated noise.
  3. *(Non-Horizontal Recovery of Theorem 1):* On the three-cluster configuration of Theorem 1 with $D_3 := "dist"(C_3, C_1 union C_2) > overline(tau)_m (C_3)$ and minimum cluster size $m_"cl" <= min_r |C_r|$, HDBSCAN's Excess-of-Mass (EoM) optimization recovers the exact three-cluster partition ${C_1, C_2, C_3}$ whenever
     $
       1 / (overline(tau)_m (C_r)) > 2 / g_(12) - 1 / D_3 quad (r in {1, 2}),
     $
     a condition automatically satisfied whenever $g_(12) > 2 max{overline(tau)_m (C_1), overline(tau)_m (C_2)}$.
]

#proof[
  *Part 1:* By definition of mutual reachability, an edge ${u, v}$ with $u != v$ satisfies $d_"mreach" (u, v) <= epsilon$ if and only if $d_"core" (u) <= epsilon$, $d_"core" (v) <= epsilon$, and $d(u, v) <= epsilon$—precisely the defining condition of $E_"core" (epsilon)$. Thus the $epsilon$-threshold graph $G_"mreach" (epsilon)$ restricted to $V_"core" (epsilon)$ equals $G_("DB"^*) (epsilon)$, while every vertex with $d_"core" (w) > epsilon$ has all incident $d_"mreach"$-edges strictly greater than $epsilon$ and is therefore isolated. By the classical cycle property of MSTs, for any edge-weighted graph, the connected components of its $epsilon$-threshold subgraph coincide with those of $cal(T)_"MST" (epsilon)$ because any edge $e = {u, v}$ of weight $<= epsilon$ absent from $cal(T)_"MST"$ closes a cycle whose every tree edge has weight $<= w(e) <= epsilon$.

  *Part 2:* Suppose $o_(t_0)$ is the earliest vertex of a DBSCAN\*($epsilon$) core component $C$ popped into $pi$, and assume inductively that the processed subset $S inter C$ is non-empty while $C backslash S != emptyset$. Because $C$ is connected in $G_("DB"^*) (epsilon)$, there exist $u in S inter C$ and $v in C backslash S$ with $d(u, v) <= epsilon$ and $d_"core" (u) <= epsilon$. Upon processing $u$, OPTICS updates $r[v] <= d_"reach" (v | u) = max{d_"core" (u), d(u, v)} <= epsilon$. Because OPTICS always extracts the unprocessed vertex minimizing $r[dot]$, no vertex $z$ whose current reachability exceeds $epsilon$ can be popped while $C backslash S != emptyset$; hence every remaining $v in C$ enters $pi$ with $r[v] <= epsilon$ prior to exiting the valley. For any non-core border point $b$ with $d_"core" (b) > epsilon$ lying within $epsilon$ of some $u in C$, OPTICS computes $d_"reach" (b | u) = max{d_"core" (u), d(u, b)} <= epsilon$, allowing $b$ to enter the valley before the valley terminates; by contrast, $d_"mreach" (b, u) >= d_"core" (b) > epsilon$ excludes $b$ from every edge of $cal(T)_"MST" (epsilon)$.

  *Part 3:* Indexing HDBSCAN's condensed tree $cal(T)_"cond"$ by $lambda = 1 / epsilon$, the root splits at $lambda_"root" = 1 / D_3$ into two valid branches ($|dot| >= m_"cl"$): the leaf $C_3$ (persisting up to at least $lambda_3 >= 1 / (overline(tau)_m (C_3))$) and the composite branch $B_(12) = C_1 union C_2$. At $lambda_"gap" = 1 / g_(12) > lambda_"root"$, $B_(12)$ bifurcates into valid leaves $C_1$ and $C_2$, whose points remain connected until at least $lambda_r := 1 / (overline(tau)_m (C_r))$ ($r = 1, 2$). By definition of HDBSCAN stability $S(C) = sum_(x in C) (lambda_"max" (x, C) - lambda_"min" (C))$ @campello2015hierarchical,
  $
    S(B_(12)) &= (|C_1| + |C_2|)(lambda_"gap" - lambda_"root"), \
    S(C_1) + S(C_2) &>= |C_1|(lambda_1 - lambda_"gap") + |C_2|(lambda_2 - lambda_"gap").
  $
  The hypothesis $lambda_r > 2 lambda_"gap" - lambda_"root"$ implies $lambda_r - lambda_"gap" > lambda_"gap" - lambda_"root"$ for both $r = 1, 2$, yielding $S(C_1) + S(C_2) > S(B_(12))$ strictly. Consequently, HDBSCAN's bottom-up dynamic program chooses ${C_1, C_2}$ over $B_(12)$ while simultaneously retaining leaf $C_3$ (as $S(C_3) > 0$), selecting a *non-horizontal* antichain cut in $cal(T)_"cond"$ that recovers ${C_1, C_2, C_3}$ exactly.
]

*Comparison with VDBSCAN (Alg. 50):* VDBSCAN @liu2007vdbscan bypasses MST and heap construction by sorting empirical $k$-distances $s_1 <= dots <= s_n$, extracting $L$ knee radii $epsilon_1 < dots < epsilon_L$, and running sequential *density peeling*: $(cal(C)_ell, U_(ell+1)) = "DBSCAN"(U_ell, epsilon_ell, k + 1)$. On Theorem 1's configuration, if the $k$-dist curve exhibits clean plateaus yielding $max{overline(tau)_m (C_1), overline(tau)_m (C_2)} <= epsilon_1 < g_(12)$ and $epsilon_2 in [overline(tau)_m (C_3), D_3)$, pass $ell = 1$ isolates and deletes $C_1, C_2$ completely before pass $ell = 2$ expands to $epsilon_2$ and recovers $C_3$. However, VDBSCAN suffers two structural weaknesses absent in HDBSCAN: (i) *plateau-masking under continuous density gradients* (e.g., Gaussian mixtures), where $s_i$ rises smoothly without discrete knees; and (ii) *irreversible peeling leakage*: if pass $ell = 1$ leaves unclustered low-density fringe points of $C_1$ and $C_2$ inside $U_2$, pass $ell = 2$ at coarser radius $epsilon_2 > g_(12)$ can stitch those residual fringes into a spurious "halo ring" cluster wrapping around $C_1 union C_2$! HDBSCAN avoids residual fringes completely because condensation tracks point departure monotonically along a single unified tree.

#subsec("C", "Rank-Adaptive Scale Invariance: SNN vs. RNN-DBSCAN")

Instead of exploring a hierarchy of metric radii, SNN (Alg. 36, @ertoz2003snn) and RNN-DBSCAN (Alg. 45, @bryant2018rnndbscan) abolish metric scale altogether by operating purely on the directed $k$-nearest-neighbor ($k$-NN) graph $G_k = (cal(X), E_k)$ where $(u, v) in E_k <=> v in N_k (u)$. SNN defines edge weights by mutual neighborhood overlap $s_"SNN" (i, j) = bold(1)[j in N_k (i) and i in N_k (j)] |N_k (i) inter N_k (j)|$ and requires three integers $(k, "Eps", "MinPts")$; RNN-DBSCAN compresses the entire specification to a *single* integer $k$, declaring $x$ a core point iff its reverse-$k$-NN in-degree satisfies $|R_k (x)| >= k$, connecting core components along $G_k$ edges (symmetrized over core vertices), and attaching border points within each cluster's adaptive internal radius $rho(C)$.

#theorem(title: "Local Homothety Invariance, Non-Empty Core Conservation, and Hubness Divergence", origin: "Original Proof")[
  Let $cal(X) = union.sq.big_(r=1)^K C_r subset.eq bb(R)^d$ satisfy $|C_r| > k$ and the isolation margin $"dist"(C_r, cal(X) backslash C_r) > overline(tau)_(k+1) (C_r)$ for every $r = 1, dots, K$.
  1. *(Multi-Scale Homothety Invariance):* Apply arbitrary independent cluster-specific scalings $tilde(C)_r = {c_r + alpha_r (x - c_r) : x in C_r}$ with $alpha_r > 0$ that preserve $"dist"(tilde(C)_r, tilde(cal(X)) backslash tilde(C)_r) > alpha_r overline(tau)_(k+1) (C_r)$ for all $r$. Then the directed $k$-NN graph $G_k$, the SNN similarity matrix $s_"SNN"$, the reverse-$k$-NN counts $|R_k (x)|$, and the core-cluster partitions of both SNN (Alg. 36) and RNN-DBSCAN (Alg. 45) are *identically invariant* to $(alpha_1, dots, alpha_K)$.
  2. *(Exact In-Degree Conservation & Guaranteed Non-Empty Core):* On *any* finite dataset $cal(X)$ of size $n > k >= 1$, RNN-DBSCAN satisfies the exact global conservation law
     $
       1 / n sum_(x in cal(X)) |R_k (x)| = k,
     $
     which implies $cal(C)_k = {x : |R_k (x)| >= k} != emptyset$ unconditionally (unlike DBSCAN and SNN, whose core sets vanish for overly strict $epsilon$ or $"Eps"$).
  3. *(High-Dimensional Hubness Asymmetry):* In $G_k$, every vertex has SNN-graph degree strictly bounded by $rho_"SNN" (i) <= k$, whereas in RNN-DBSCAN a single hub $h$ can attain $|R_k (h)| = n - 1$, enabling one-way $k$-NN arcs from distant clusters to bridge distinct core components in $H_k$.
]

#proof[
  *Part 1:* Because $"dist"(tilde(C)_r, tilde(cal(X)) backslash tilde(C)_r) > alpha_r overline(tau)_(k+1) (C_r)$, for every $x in C_r$, all $k$ nearest neighbors of its scaled image $tilde(x) in tilde(C)_r$ remain strictly inside $tilde(C)_r$. Within $tilde(C)_r$, all pairwise distances scale uniformly by the positive constant $alpha_r$: $norm(tilde(x) - tilde(y))_2 = alpha_r norm(x - y)_2$. Because positive scaling preserves strict distance ordering, $N_k (tilde(x)) = {tilde(y) : y in N_k (x)}$ with identical neighbor ranks $r_i (u)$. Since both SNN ($s_"SNN"$, $s_"SNN"^w$, $rho_"SNN"$) and RNN-DBSCAN ($R_k$, core set $cal(C)_k$, core graph $H_k$, and cluster-normalized border test $norm(tilde(z) - tilde(y))_2 <= rho(tilde(C)_r) = alpha_r rho(C_r)$) depend strictly on $G_k$ and intra-cluster distance ratios, both partitions are identically invariant to $(alpha_1, dots, alpha_K)$.

  *Part 2:* Count the directed arcs of $G_k$ in two ways (double-counting over out-degrees and in-degrees):
  $
    sum_(x in cal(X)) |R_k (x)| &= sum_(x in cal(X)) sum_(y in cal(X) backslash {x}) bold(1)[x in N_k (y)] \
    &= sum_(y in cal(X)) |N_k (y)| = n k.
  $
  Dividing by $n$ gives the exact mean $k$. By the pigeonhole principle, $max_(x in cal(X)) |R_k (x)| >= k$, guaranteeing $cal(C)_k != emptyset$ for every dataset and every $k in {1, dots, n - 1}$.

  *Part 3:* In SNN, $s_"SNN" (i, j) > 0$ requires the mutual-neighborhood gate $j in N_k (i) and i in N_k (j)$; hence $N_"SNN" (i) subset.eq N_k (i)$ and $rho_"SNN" (i) = |N_"SNN" (i)| <= |N_k (i)| = k$ identically. Conversely, consider a star-like high-dimensional hub configuration where $h$ lies near the global centroid of $K$ orthogonal spheres in $bb(R)^d$ ($d >> 1$) so that points across multiple clusters include $h$ in their $k$-NN lists while $h$ only selects $k$ nearest neighbors; then $|R_k (h)|$ grows up to $n - 1$. Because RNN-DBSCAN connects core vertices $u, v in cal(C)_k$ whenever $v in N_k (u) or u in N_k (v)$ (one-sided symmetrization), any core point $u$ pointing to $h$ merges with $h$'s component even when $h$ does not reciprocate $u in N_k (h)$.
]

Theorem 3 clarifies the exact division of labor between SNN and RNN-DBSCAN: in low-to-moderate dimensions ($d <= 20$) with sharp multi-density variation, RNN-DBSCAN achieves near-optimal parsimony (1 integer hyperparameter $k$, guaranteed non-empty cores, and adaptive cluster radii $rho(C)$); in high-dimensional sparse or text/embedding regimes ($d >= 50$) where distance concentration induces severe hubness @beyer1999nearest, SNN's mutual-overlap gate and rank-weighted similarity $s_"SNN"^w$ act as a structural firewall against hub-induced cluster collapse, at the cost of tuning three coupled thresholds $(k, "Eps", "MinPts")$.

#subsec("D", "Axiomatic, Spatiotemporal, and Distributional Generalizations: GDBSCAN, ST-DBSCAN, and DBCLASD")

Three additional extensions reshape the *domain semantics* and *statistical null model* of density connectivity.

#proposition(title: "Predicate Hierarchy, Order-Invariance Bifurcation, and Poisson Over-Dispersion Rejection", origin: "Original Unified Proof")[
  Let $cal(D)$ be a finite object database.
  1. *(Strict Axiomatic Inclusion & Order Bifurcation):* Standard DBSCAN (Alg. 31) and pure ST-DBSCAN (Alg. 40 with attribute bound $Delta epsilon = infinity$) are exact special cases of GDBSCAN (Alg. 39) with reflexive, symmetric predicates $"NPred"$ and additive cardinality $"wCard"(S) = |S|$; consequently, all three enjoy *strict core-component permutation invariance*. In contrast, full ST-DBSCAN with finite $Delta epsilon < infinity$ violates predicate symmetry and *destroys core-component order invariance*.
  2. *(Over-Dispersion Breakdown of DBCLASD on Non-Uniform Modes):* Let $X_1, dots, X_N$ be i.i.d. draws from a smooth non-uniform density $f$ on $bb(R)^d$ with local intensity $Lambda = N f(X)$ satisfying $"Var"(Lambda^(-1)) > 0$, and let $U = V_d R_1^d$ denote the rescaled nearest-neighbor ball volume. Under DBCLASD's homogeneous Poisson null, $"CV"^2 (U) = "Var"(U) / (bb(E)[U])^2 = 1$; under the inhomogeneous density $f$, as $N -> infinity$,
     $
       lim_(N -> infinity) "CV"^2 (N V_d R_1^d) = 1 + 2 ("Var"(f(X)^(-1))) / (bb(E)[f(X)^(-1)])^2 > 1,
     $
     causing DBCLASD's Pearson $chi^2$ goodness-of-fit test against $F_(N, hat(A), d)$ to reject any non-uniform unimodal cluster (e.g., an isotropic Gaussian) with probability approaching $1$ as $N -> infinity$.
]

#proof[
  *Part 1:* For pure ST-DBSCAN ($Delta epsilon = infinity$), define on $cal(S) times cal(T)$ the normalized $L_infinity$ product metric
  $
    d_"ST" ((s_i, t_i), (s_j, t_j)) := max { (d_S (s_i, s_j)) / epsilon_1, space (d_T (t_i, t_j)) / epsilon_2 }.
  $
  Then $p_j in N_(epsilon_1, epsilon_2) (p_i) <=> d_"ST" (p_i, p_j) <= 1$, which is reflexive and symmetric, placing pure ST-DBSCAN identically inside GDBSCAN (and inside metric DBSCAN at radius $1$). For any reflexive, symmetric $"NPred"$ with fixed $"wCard"$, the core set $cal(D)_"core" = {o : "wCard"(N_"NPred" (o)) >= "MinCard"}$ and undirected core graph $(cal(D)_"core", {{u, v} : "NPred"(u, v)})$ are independent of object visitation order, so BFS extracts identical core components under any permutation of $cal(D)$.

  To prove order dependence of ST-DBSCAN when $Delta epsilon < infinity$, consider three collinear spatial--temporal core points $p_1, p_2, p_3$ all pairwise within $(epsilon_1, epsilon_2)$ (padded with coincident duplicates at each attribute value so $|N_(epsilon_1, epsilon_2)| >= "MinPts"$), carrying scalar attributes $a_1 = 0$, $a_2 = 1$, $a_3 = 1.9$, and set $Delta epsilon = 1.2$. If BFS seeds at $p_1$: initially $overline(a) = 0$; $p_2$ is admitted ($|1 - 0| = 1 <= 1.2$), shifting $overline(a)$ to $0.5$; then $p_3$ is *rejected* because $|1.9 - 0.5| = 1.4 > 1.2$, splitting $p_3$ into a separate cluster. If BFS instead seeds at $p_2$: initially $overline(a) = 1$; both $p_1$ ($|0 - 1| = 1 <= 1.2$) and $p_3$ ($|1.9 - 1| = 0.9 <= 1.2$) are admitted into a *single* cluster! Thus running-mean gating breaks core order invariance.

  *Part 2:* Under a homogeneous spatial Poisson process of constant intensity $lambda > 0$ in $bb(R)^d$, the probability that no other point falls inside a ball of volume $u = V_d r^d$ is $exp(-lambda u)$, so $U = V_d R_1^d tilde.op "Exp"(lambda)$ with $bb(E)[U] = lambda^(-1)$, $"Var"(U) = lambda^(-2)$, and $"CV"^2 (U) = 1$. Now let $X_1, dots, X_N$ be i.i.d. from a smooth continuous density $f$ with $bb(E)[f(X)^(-2)] < infinity$. Conditioned on $X_1 = x$ (where $f(x) > 0$), the rescaled nearest-neighbor volume $W_N := N V_d R_1^d$ satisfies
  $
    bb(P)(W_N > w | X_1 = x) &= (1 - integral_(B(x, (w/(N V_d))^(1/d))) f(y) dif y)^(N-1) \
    &-> exp(-f(x) w),
  $
  so asymptotically $W | X tilde.op "Exp"(f(X))$, giving $bb(E)[W | X] = f(X)^(-1)$ and $"Var"(W | X) = f(X)^(-2)$. Applying the Law of Total Variance yields
  $
    "Var"(W) &= bb(E)["Var"(W | X)] + "Var"(bb(E)[W | X]) \
             &= bb(E)[f(X)^(-2)] + "Var"(f(X)^(-1)) \
             &= (bb(E)[f(X)^(-1)])^2 + 2 "Var"(f(X)^(-1)).
  $
  Dividing by $(bb(E)[W])^2 = (bb(E)[f(X)^(-1)])^2$ yields the exact identity $"CV"^2 (W) = 1 + 2 "Var"(f(X)^(-1)) / (bb(E)[f(X)^(-1)])^2 > 1$. Because DBCLASD fits a single homogeneous parameter $hat(lambda) = N / hat(A)$ (forcing null $"CV"^2 = 1$), the empirical NN histogram diverges by a strictly positive total-variation distance from $F_(N, hat(A), d)$, so the Pearson statistic satisfies $chi^2 (C) = Omega(N) -> infinity$, guaranteeing rejection of any non-uniform mode at fixed level $alpha$.
]

Proposition 1 exposes a critical practical insight often glossed over in surveys: DBCLASD's zero-parameter appeal @xu1998dbclasd holds *only* for spatially uniform regions (e.g., uniform urban building footprints or planar mineral deposits); on Gaussian or peaked statistical mixtures, its own statistical test actively fragments true clusters into concentric density rings!

#wide_figure(
  caption: [*Table I:* Head-to-Head Structural Comparison of Density-Reachability and Variable-Density Algorithms (Section I).],
  table(
    columns: (1.08fr, 0.95fr, 1.15fr, 0.85fr, 0.85fr, 1.12fr),
    align: (left, left, left, center, center, left),
    [Algorithm], [Density Primitive], [Key Parameters], [Multi-Density], [Core Order-Inv.], [Primary Failure Mode],
    [DBSCAN (\#31) @ester1996dbscan], [$epsilon$-ball count], [$epsilon$, $"MinPts"$], [No (Thm. 1)], [Yes], [Multi-scale merge/erasure; high-$d$],
    [HDBSCAN (\#32) @campello2015hierarchical], [Sym. $d_"mreach"$ MST], [$m_"pts"$, $m_"cl"$], [Yes (EoM)], [Yes], [Slow $O(n^2)$ in high $d$; flat-cut bias],
    [OPTICS (\#33) @ankerst1999optics], [Asym. $d_"reach"$ heap], [$epsilon_"max"$, $"MinPts"$, $xi$], [Yes (plot)], [Core only], [Sensitive $xi$ tuning; border asymmetry],
    [SNN (\#36) @ertoz2003snn], [Mutual $k$-NN overlap], [$k$, $"Eps"$, $"MinPts"$], [Yes (Thm. 3)], [Yes], [3 coupled params; $O(n^2 d + n k^2)$],
    [DBCLASD (\#38) @xu1998dbclasd], [Poisson NN $chi^2$ test], [None (internal $alpha$)], [Partial], [No], [Over-dispersion splits Gaussians (Prop. 1)],
    [GDBSCAN (\#39) @sander1998gdbscan], [Custom $"NPred", "wCard"$], [$"NPred"$, $"MinCard"$], [Via $"NPred"$], [Yes (sym.)], [Index availability for arbitrary $"NPred"$],
    [ST-DBSCAN (\#40) @birant2007stdbscan], [Cylinder + attr. mean], [$epsilon_1, epsilon_2, "MinPts", Delta epsilon$], [No], [Only $Delta epsilon = infinity$], [Order-dependent drift if $Delta epsilon < infinity$],
    [RNN-DBSCAN (\#45) @bryant2018rnndbscan], [Reverse $k$-NN $|R_k|$], [$k$ (single param)], [Yes (Thm. 3)], [Yes], [High-$d$ hub bridging (Thm. 3.3)],
    [VDBSCAN (\#50) @liu2007vdbscan], [Peeling $k$-dist strata], [$k$, knee $tau$ / $L$], [Discrete $L$], [Core/pass], [Smooth gradients; halo-ring leakage],
  ),
)

// =============================================================================
// SECTION II
// =============================================================================
#sec("II", "Kernel Mode-Seeking and Gradient Ascent Families")

The second major pillar of Category II defines clusters not via connected superlevel components, but via *basins of attraction* of local maxima of a nonparametric kernel density estimator $hat(f)_h (x) = (n h^d)^(-1) sum_(i=1)^n K((x - x_i) / h)$. Seven algorithms span this family: classical Mean Shift (Alg. 34, @cheng1995mean @comaniciu2002mean), DENCLUE 1.0 (Alg. 35, @hinneburg1998denclue), DENCLUE 2.0 (Alg. 43, @hinneburg2007denclue2), Blurring Mean Shift (Alg. 54, @cheng1995mean @carreira2006fast), Mean Shift++ (Alg. 55, @jang2021meanshiftpp), Medoidshift (Alg. 42, @sheikh2007medoidshift), and Quickshift (Alg. 41, @vedaldi2008quickshift). Comparing them rigorously exposes three foundational tradeoffs: (A) fixed-step normalized ascent versus adaptive EM fixed-point dynamics; (B) frozen-KDE linear convergence versus blurring cubic collapse and Cartesian grid quantization; and (C) continuous Euclidean shifting versus discrete exemplar medoids and acyclic parent trees.

#subsec("A", "Fixed-Step vs. Adaptive EM Fixed-Point Ascent: DENCLUE 1.0, DENCLUE 2.0, and Mean Shift")

Although DENCLUE 1.0 @hinneburg1998denclue, Mean Shift @cheng1995mean @comaniciu2002mean, and DENCLUE 2.0 @hinneburg2007denclue2 emerged from distinct research communities (spatial databases vs. computer vision), their mathematical relationship is remarkably tight under an isotropic Gaussian kernel $K_h (u) = (2 pi h^2)^(-d/2) exp(-norm(u)_2^2 / (2 h^2))$.

#theorem(title: "Exact Equivalence of Gaussian Mean Shift and DENCLUE 2.0 vs. Limit-Cycle Oscillation of DENCLUE 1.0", origin: "Original Unified Proof")[
  Let $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$ and let $hat(f)_h (x) = n^(-1) sum_(i=1)^n K_h (x - x_i)$ be the isotropic Gaussian KDE with bandwidth $h > 0$.
  1. *(Exact Operator Identity & Logarithmic EM Ascent):* The DENCLUE 2.0 adaptive update $T(x) = sum_(i=1)^n w_i (x) x_i$ with softmax weights $w_i (x) = (K_h (x - x_i)) / (sum_j K_h (x - x_j))$ is identically equal to the Gaussian Mean Shift update $x + m_(h, G) (x) = x + (h^2 / (hat(f)_h (x))) nabla hat(f)_h (x)$. Moreover, every iteration satisfies the strict quantitative ascent bound
     $
       log (hat(f)_h (T(x))) / (hat(f)_h (x)) >= (norm(T(x) - x)_2^2) / (2 h^2),
     $
     implying $sum_(t=0)^infinity norm(x^((t+1)) - x^((t)))_2^2 <= 2 h^2 log (hat(f)_"max" / (hat(f)_h (x^((0))))) < infinity$ and $lim_(t -> infinity) norm(nabla hat(f)_h (x^((t))))_2 = 0$.
  2. *(Irreducible Localization Barrier of DENCLUE 1.0):* Under DENCLUE 1.0's fixed-step normalized ascent $z^((t+1)) = z^((t)) + delta (nabla hat(f)_h (z^((t)))) / norm(nabla hat(f)_h (z^((t))))_2$ ($delta > 0$), for *any* isolated strict local maximum $x^*$ of $hat(f)_h$, every trajectory that does not hit a stationary point in finitely many steps satisfies
     $
       limsup_(t -> infinity) norm(z^((t)) - x^*)_2 >= delta / 2,
     $
     precluding asymptotic convergence to $x^*$ without auxiliary step-size decay.
]

#proof[
  *Part 1:* Differentiating $K_h (x - x_i)$ gives $nabla K_h (x - x_i) = -h^(-2) K_h (x - x_i)(x - x_i)$, hence
  $
    nabla hat(f)_h (x) = 1 / (n h^2) sum_(i=1)^n K_h (x - x_i)(x_i - x) = (hat(f)_h (x)) / h^2 (T(x) - x).
  $
  Solving for $T(x)$ immediately yields $T(x) = x + (h^2 / (hat(f)_h (x))) nabla hat(f)_h (x)$, which matches the Mean Shift vector $m_(h, G) (x) = T(x) - x$ for the self-shadow Gaussian profile $g(r) = -k'(r) = (1/2) k(r)$. For any candidate $x' in bb(R)^d$, rewrite the log-density ratio using $sum_i w_i (x) = 1$ and apply Jensen's inequality to the strictly concave logarithm:
  $
    log (hat(f)_h (x')) / (hat(f)_h (x)) &= log sum_(i=1)^n w_i (x) (K_h (x' - x_i)) / (K_h (x - x_i)) \
    &>= sum_(i=1)^n w_i (x) log (K_h (x' - x_i)) / (K_h (x - x_i)) \
    &= 1 / (2 h^2) sum_(i=1)^n w_i (x) (norm(x - x_i)_2^2 - norm(x' - x_i)_2^2).
  $
  Setting $x' = T(x) = sum_i w_i (x) x_i$, the weighted parallel-axis identity gives $sum_i w_i (x) norm(x - x_i)_2^2 = sum_i w_i (x) norm(x' - x_i)_2^2 + norm(x' - x)_2^2$, reducing the lower bound to $norm(T(x) - x)_2^2 / (2 h^2) >= 0$ with strict positivity whenever $nabla hat(f)_h (x) != bold(0)$. Telescoping over $t = 0, dots, N - 1$ yields
  $
    sum_(t=0)^(N-1) norm(x^((t+1)) - x^((t)))_2^2 <= 2 h^2 log (hat(f)_h (x^((N)))) / (hat(f)_h (x^((0)))) <= 2 h^2 log (hat(f)_"max") / (hat(f)_h (x^((0)))),
  $
  since $hat(f)_h (x) <= hat(f)_"max" := (2 pi h^2)^(-d/2)$. Letting $N -> infinity$ gives $norm(x^((t+1)) - x^((t)))_2 -> 0$, and because $hat(f)_h (x^((t))) >= hat(f)_h (x^((0))) > 0$, $norm(nabla hat(f)_h (x^((t))))_2 = ((hat(f)_h (x^((t)))) / h^2) norm(x^((t+1)) - x^((t)))_2 -> 0$.

  *Part 2:* Along any non-terminating DENCLUE 1.0 orbit, $nabla hat(f)_h (z^((t))) != bold(0)$ for all $t >= 0$, so every step has exact Euclidean length $norm(z^((t+1)) - z^((t)))_2 = delta$. By the triangle inequality around $x^*$,
  $
    delta = norm(z^((t+1)) - z^((t)))_2 <= norm(z^((t+1)) - x^*)_2 + norm(z^((t)) - x^*)_2,
  $
  whence $max{norm(z^((t)) - x^*)_2, norm(z^((t+1)) - x^*)_2} >= delta / 2$ for every $t >= 0$. Taking the limit superior yields $limsup_(t -> infinity) norm(z^((t)) - x^*)_2 >= delta / 2$, proving persistent limit-cycle oscillation around $x^*$.
]

*Architectural Distinction Beyond the Update Rule:* While Gaussian Mean Shift and DENCLUE 2.0 share the identical hill-climbing map $T$, they diverge operationally in *spatial pruning* and *cluster topology*: (i) DENCLUE 1.0/2.0 @hinneburg1998denclue @hinneburg2007denclue2 pre-index occupied $2 h$-hypercubes and discard cubes below density threshold $xi$ (yielding explicit noise isolation and subquadratic local sums in low $d$, whereas standard Mean Shift evaluates all $n$ kernels and partitions 100% of $cal(X)$ without noise pruning); and (ii) DENCLUE explicitly supports *arbitrary-shape clusters* by merging basins whose attractors $x_r^*, x_s^*$ are linked by a superlevel path $hat(f)_h (gamma(s)) >= xi$—bridging mode-seeking directly back to Section I's level-set connectivity! Moreover, Mean Shift extends naturally to compactly supported non-Gaussian kernels (Epanechnikov, uniform) that attain *finite-step convergence* @comaniciu2002mean @ariascastro2016meanshift, whereas DENCLUE 2.0's exact EM closed form relies specifically on the Gaussian (or log-concave exponential-family) kernel.

#subsec("B", "Linear vs. Cubic Contraction and Grid Quantization: Mean Shift, Blurring Mean Shift, and Mean Shift++")

Standard Mean Shift (Alg. 34) freezes the reference point cloud $cal(X)^((0))$ while moving query trajectories $y_i^((t))$. Blurring Mean Shift (BMS, Alg. 54, @cheng1995mean @carreira2006fast) synchronously replaces $cal(X)^((t))$ by the shifted points $bold(X)^((t+1)) = (bold(D)^((t)))^(-1) bold(W)^((t)) bold(X)^((t))$ at every step, dramatically accelerating intra-cluster collapse at the price of a subtle asymptotic pathology. Mean Shift++ (Alg. 55, @jang2021meanshiftpp) in turn approximates blurring shifts over a hashed Cartesian grid of width $h$ to achieve linear $O(n)$ runtime in fixed dimension $d$. We prove both the exact contraction gap and the global consensus collapse of BMS, followed by the exact quantization sandwich and crossover dimension of Mean Shift++.

#theorem(title: "Cubic vs. Linear Local Contraction and Inevitable Global Collapse of Blurring Mean Shift", origin: "Original Proof")[
  Let $h > 0$ be the isotropic Gaussian bandwidth.
  1. *(Population Rate Gap — Linear vs. Cubic):* For a Gaussian population $X^((0)) tilde.op cal(N)(mu, Sigma_0)$ with covariance eigenvalues $sigma_(0, 1)^2, dots, sigma_(0, d)^2 > 0$, standard frozen Mean Shift contracts the $j$-th principal standard deviation *linearly*:
     $
       sigma_(t+1, j)^"MS" = (sigma_(0, j)^2) / (sigma_(0, j)^2 + h^2) sigma_(t, j)^"MS",
     $
     whereas Gaussian Blurring Mean Shift contracts it *cubically*:
     $
       sigma_(t+1, j)^"BMS" = ((sigma_(t, j)^"BMS")^3) / ((sigma_(t, j)^"BMS")^2 + h^2) <= ((sigma_(t, j)^"BMS")^3) / h^2.
     $
  2. *(Finite-Sample Global Consensus Collapse):* For *any* finite dataset $cal(X)^((0)) subset.eq bb(R)^d$ of diameter $D_0 := "diam"(cal(X)^((0))) < infinity$, Gaussian Blurring Mean Shift satisfies the nested convex-hull inclusion $cal(X)^((t+1)) subset.eq "conv"(cal(X)^((t)))$ and strict geometric diameter contraction
     $
       "diam"(cal(X)^((t+1))) <= (1 - exp(- D_0^2 / (2 h^2))) "diam"(cal(X)^((t)))
     $
     for all $t >= 0$. Hence there exists a single consensus point $c^* in "conv"(cal(X)^((0)))$ such that $lim_(t -> infinity) x_i^((t)) = c^*$ for *every* $i in {1, dots, n}$, proving mathematically that Gaussian BMS *must* employ an early-stopping criterion to avoid destroying multi-cluster structure.
]

#proof[
  *Part 1:* Without loss of generality set $mu = bold(0)$ and rotate coordinates to diagonalize $Sigma_t = "diag"(sigma_(t, 1)^2, dots, sigma_(t, d)^2)$. Because both the Gaussian density and isotropic Gaussian kernel factor across coordinates, each axis decouples independently. For a 1D distribution $X tilde.op cal(N)(0, sigma^2)$ evaluated at query $y$, completing the square in the exponent gives
  $
    (x^2) / (2 sigma^2) + ((y - x)^2) / (2 h^2) = ((sigma^2 + h^2) / (2 sigma^2 h^2)) (x - (sigma^2) / (sigma^2 + h^2) y)^2 + (y^2) / (2 (sigma^2 + h^2)).
  $
  Thus the population Gaussian Mean Shift operator is the exact linear map $T_sigma (y) = (sigma^2 / (sigma^2 + h^2)) y$. In standard Mean Shift, the reference variance remains frozen at $sigma = sigma_(0, j)$, so $Y^((t+1)) = (sigma_(0, j)^2 / (sigma_(0, j)^2 + h^2)) Y^((t))$, yielding $sigma_(t+1, j)^"MS" = (sigma_(0, j)^2 / (sigma_(0, j)^2 + h^2)) sigma_(t, j)^"MS"$. In Blurring Mean Shift, the reference distribution updates to $X^((t)) tilde.op cal(N)(0, sigma_(t, j)^2)$, so $X^((t+1)) = (sigma_(t, j)^2 / (sigma_(t, j)^2 + h^2)) X^((t))$, multiplying the current standard deviation $sigma_(t, j)$ by $sigma_(t, j)^2 / (sigma_(t, j)^2 + h^2)$ to yield the exact cubic law $sigma_(t+1, j)^"BMS" = (sigma_(t, j)^"BMS")^3 / ((sigma_(t, j)^"BMS")^2 + h^2)$ @carreira2006fast.

  *Part 2:* At step $t$, each row of $bold(P)^((t)) = (bold(D)^((t)))^(-1) bold(W)^((t))$ defines convex combination weights over $cal(X)^((t))$, so $x_i^((t+1)) in "conv"(cal(X)^((t)))$ and hence $"diam"(cal(X)^((t))) <= D_0$ inductively for all $t >= 0$. Consequently, every Gaussian weight satisfies $exp(-D_0^2 / (2 h^2)) <= W_(i j)^((t)) <= 1$, which bounds every row-stochastic entry below by
  $
    P_(i j)^((t)) = (W_(i j)^((t))) / (sum_(ell=1)^n W_(i ell)^((t))) >= 1 / n exp(- D_0^2 / (2 h^2)) =: eta_0 > 0.
  $
  Fix any pair $i, k in {1, dots, n}$ and define the overlap mass $alpha_(i k) := sum_(j=1)^n min(P_(i j)^((t)), P_(k j)^((t))) >= n eta_0 = exp(-D_0^2 / (2 h^2))$. If $alpha_(i k) = 1$, then $x_i^((t+1)) = x_k^((t+1))$. Otherwise set $p_j^+ := (P_(i j)^((t)) - P_(k j)^((t)))_+$ and $p_j^- := (P_(k j)^((t)) - P_(i j)^((t)))_+$; both non-negative vectors sum to $1 - alpha_(i k) > 0$. Defining $u^+ := (1 - alpha_(i k))^(-1) sum_j p_j^+ x_j^((t))$ and $u^- := (1 - alpha_(i k))^(-1) sum_j p_j^- x_j^((t))$, both $u^+$ and $u^-$ lie inside $"conv"(cal(X)^((t)))$, and
  $
    norm(x_i^((t+1)) - x_k^((t+1)))_2 &= (1 - alpha_(i k)) norm(u^+ - u^-)_2 \
    &<= (1 - exp(-D_0^2 / (2 h^2))) "diam"(cal(X)^((t))),
  $
  since $"diam"("conv"(cal(X)^((t)))) = "diam"(cal(X)^((t)))$. Taking the maximum over $(i, k)$ yields the stated contraction. Because the nested compact convex hulls $"conv"(cal(X)^((t)))$ shrink to diameter zero, Cantor's intersection theorem yields a unique consensus limit $c^* = lim_(t -> infinity) x_i^((t))$ for all $i$.
]

Theorem 5 illuminates the stark practical contrast between Algs. 34 and 54: on well-separated clusters ($exp(-D_0^2/(2h^2)) << 1$), BMS exhibits a *two-timescale metastability*: intra-cluster modes collapse to machine precision in $3$--$6$ cubic steps while inter-cluster drift remains exponentially frozen at $O(exp(-D_"gap"^2 / (2 h^2)))$ per step, enabling dramatic iteration savings @carreira2006fast; yet running BMS blindly to convergence collapses all clusters into a single point! Compactly supported kernels (or Mean Shift++'s finite grid neighborhood) break this global coupling whenever inter-cluster gaps exceed $2 h$, turning clusters into strictly invariant disconnected components.

#proposition(title: "Grid-Window Sandwich and Exact Dimensional Crossover of Mean Shift++", origin: "Original Proof")[
  Let $q_h (x) = (floor(x_1 / h), dots, floor(x_d / h))$ be the Mean Shift++ cell map (Alg. 55, @jang2021meanshiftpp), with hypercube $Q_h (c) = product_(r=1)^d [c_r h, (c_r + 1) h)$ and $3^d$-cell Moore neighborhood $M_h (c) = {c + Delta : Delta in {-1, 0, 1}^d}$.
  1. *(Continuous Box-Kernel Sandwich):* For every $x in bb(R)^d$, the union of Moore cells $U_h (x) := union.big_(c' in M_h (q_h (x))) Q_h (c')$ satisfies the tight $L_infinity$-ball containment
     $
       B_infinity^circle.tiny (x, h) subset.eq U_h (x) subset.eq B_infinity^circle.tiny (x, 2 h).
     $
     Hence each Mean Shift++ step acts as a cell-quantized flat-kernel blurring update with effective axis-aligned bandwidth $in [h, 2 h)$, collapsing all points within the same occupied cell $Q_h (c)$ identically onto a single vector in one iteration.
  2. *(Asymptotic Dimensional Crossover):* Comparing per-iteration arithmetic $O(n d 3^d)$ of Mean Shift++ against $O(n^2 d)$ of Mean Shift / BMS yields the strict speedup regime
     $
       3^d < n quad <=> quad d < log_3 n,
     $
     so Mean Shift++ dominates in low dimensions (e.g., $d <= 6$ for images/LiDAR, where $3^d <= 729 << n$) but becomes asymptotically *slower* than naive quadratic Mean Shift whenever $d > log_3 n$ (e.g., $d >= 10$ at $n = 5 times 10^4$).
]

#proof[
  *Part 1:* Let $c = q_h (x)$, so $c_r h <= x_r < (c_r + 1) h$ for every coordinate $r in {1, dots, d}$. If $y in B_infinity^circle.tiny (x, h)$, then $|y_r - x_r| < h$, implying $(c_r - 1) h < y_r < (c_r + 2) h$, whence $floor(y_r / h) in {c_r - 1, c_r, c_r + 1}$ for all $r$; thus $q_h (y) in M_h (c)$ and $y in U_h (x)$. Conversely, if $y in U_h (x)$, then $y_r in [(c_r - 1) h, (c_r + 2) h)$ while $x_r in [c_r h, (c_r + 1) h)$, whence $|y_r - x_r| < 2 h$ for all $r$, giving $y in B_infinity^circle.tiny (x, 2 h)$. Because $U_h (x)$ depends on $x$ solely through $q_h (x)$, every point sharing cell $Q_h (c)$ receives the identical unweighted mean over $cal(X)^((t)) inter U_h (x)$.

  *Part 2:* Hashing $n$ points and accumulating cell sums costs $O(n d)$, and querying $3^d$ neighbors per occupied cell (at most $n$ occupied cells) costs $O(n d 3^d)$ per iteration, compared to $Theta(n^2 d)$ pairwise distances for Mean Shift / BMS. The ratio $(n d 3^d) / (n^2 d) = 3^d / n$ is strictly below $1$ iff $d < log_3 n$.
]

#subsec("C", "Discrete Exemplar vs. Directed-Tree Mode-Seeking: Medoidshift vs. Quickshift")

Both Medoidshift (Alg. 42, @sheikh2007medoidshift) and Quickshift (Alg. 41, @vedaldi2008quickshift) eliminate the Euclidean vector-space requirement of Mean Shift, operating on *arbitrary* symmetric dissimilarity matrices $D = (d_(i j))$ while eliminating iterative kernel re-evaluations and stopping tolerances. Despite this shared motivation, their mechanics induce radically different complexity and topological fidelity.

#theorem(title: "Euclidean Projection Identity, Locking-Ball Fragmentation of Medoidshift, and Quickshift Acyclicity", origin: "Original Proof")[
  Let $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$ with Euclidean metric $d(x_i, x_j) = norm(x_i - x_j)_2$, shadow kernel weights $W_(i j) = g(norm(x_i - x_j)_2^2) >= 0$ with $Z_i := sum_(j=1)^n W_(i j) > 0$, and continuous Mean Shift target $mu_i := Z_i^(-1) sum_(j=1)^n W_(i j) x_j$.
  1. *(Exact Nearest-Sample Mean Shift Identity):* Medoidshift's weighted-medoid objective satisfies
     $
       Q(k | i) := sum_(j=1)^n W_(i j) norm(x_k - x_j)_2^2 = Z_i norm(x_k - mu_i)_2^2 + V_i,
     $
     where $V_i := sum_(j=1)^n W_(i j) norm(x_j - mu_i)_2^2$ is independent of $k$. Hence $M(i) = arg min_k norm(x_k - mu_i)_2$, proving Euclidean Medoidshift equals nearest-sample projected Mean Shift.
  2. *(Locking-Ball Fragmentation Counterexample):* A sample $x_i$ locks as a spurious discrete mode ($M(i) = i$) if and only if the open ball $B^circle.tiny (mu_i, norm(x_i - mu_i)_2)$ contains no other point of $cal(X)$. In particular, on a unimodal 1D cluster containing $n - 1 >= 2$ core points at $x = 1$ and an outer flank point $x_0 = 0$, $x_0$ locks as a spurious singleton mode whenever $mu_0 <= 1 / 2$, i.e., for every Gaussian bandwidth $sigma <= (2 log(n - 1))^(-1/2)$.
  3. *(Quickshift Non-Locking & Complexity Domination):* Quickshift (Alg. 41) executes in $O(n^2 d)$ time versus Medoidshift's $O(n^3)$ matrix product $S = W D^((2))$, produces an unconditionally acyclic directed forest $G_tau$ for every $tau in (0, infinity]$, and never locks on the flank counterexample whenever $hat(p)_sigma (x_1) > hat(p)_sigma (x_0)$ and $tau >= 1$.
]

#proof[
  *Part 1:* Expand $norm(x_k - x_j)_2^2 = norm(x_k - mu_i)_2^2 + 2 chevron.l x_k - mu_i, mu_i - x_j chevron.r + norm(mu_i - x_j)_2^2$ and sum against $W_(i j)$. Because $sum_(j=1)^n W_(i j)(mu_i - x_j) = Z_i mu_i - Z_i mu_i = bold(0)$, the cross-term vanishes identically, leaving $Q(k | i) = Z_i norm(x_k - mu_i)_2^2 + V_i$. Since $Z_i > 0$ and $V_i$ is constant in $k$, minimizing $Q(k | i)$ over $k$ is equivalent to minimizing $norm(x_k - mu_i)_2$.

  *Part 2:* Under Medoidshift's self-preferring tie-breaking, $M(i) = i$ iff $norm(x_i - mu_i)_2 <= norm(x_k - mu_i)_2$ for all $k != i$, equivalently $cal(X) inter B^circle.tiny (mu_i, norm(x_i - mu_i)_2) = emptyset$. For $x_0 = 0$ and $n - 1$ points located at $x_1 = dots = x_(n-1) = 1$, the Gaussian shadow weights satisfy $W_(0 0) = 1$ and $W_(0 j) = e^(-1 / (2 sigma^2))$ ($j >= 1$), yielding
  $
    mu_0 = ((n - 1) e^(-1 / (2 sigma^2))) / (1 + (n - 1) e^(-1 / (2 sigma^2))).
  $
  Thus $mu_0 <= 1 / 2 <=> (n - 1) e^(-1 / (2 sigma^2)) <= 1 <=> sigma <= (2 log(n - 1))^(-1/2)$, in which case $B^circle.tiny (mu_0, mu_0) = (0, 2 mu_0) subset.eq (0, 1)$ contains zero sample points and $M(0) = 0$, fragmenting $x_0$ into a spurious mode!

  *Part 3:* Quickshift computes pairwise distances and KDE scores in $O(n^2 d)$ and scans for each $i$ its nearest strictly higher-density point $pi(i) = arg min_(j in cal(U)_i) d_(i j)$ in $O(n^2)$ total time, completely avoiding Medoidshift's $O(n^3)$ matrix multiplication $W D^((2))$. Because every directed edge $(i, pi_tau (i))$ strictly increases the lexicographic order $(hat(p)_sigma (x_j), -j)$, no directed cycle can exist, so $G_tau$ is a forest of directed in-trees. On the 1D flank example, $hat(p)_sigma (x_1) = (n - 1 + e^(-1/(2 sigma^2))) / n > (1 + (n - 1) e^(-1/(2 sigma^2))) / n = hat(p)_sigma (x_0)$ for all $n >= 3$ and every $sigma > 0$, so $pi(0) = 1$ and $x_0$ merges directly into the cluster mode whenever $tau >= 1$.
]

*Tradeoff Summary:* Medoidshift retains one narrow niche over Quickshift: its Minorization--Maximization (MM) descent on the Fréchet energy $Q(k | i)$ guarantees genuine local weighted-medoid optimality without requiring a manual link-cut radius $tau$ (which in Quickshift forces a delicate balance between under-splitting close modes with separation $< tau$ and over-splitting sparse fringes with spacing $> tau$). Yet Medoidshift's $O(n^3)$ wall and locking-ball fragmentation render Quickshift and Mean Shift++ decisively superior for large-scale mode-seeking.

#wide_figure(
  caption: [*Table II:* Head-to-Head Theoretical Comparison of Kernel Mode-Seeking and Gradient Ascent Algorithms (Section II).],
  table(
    columns: (1.12fr, 1.05fr, 0.92fr, 0.90fr, 0.88fr, 1.13fr),
    align: (left, left, left, center, center, left),
    [Algorithm], [Ascent / Shift Mechanism], [Complexity], [Convergence], [Space / Noise], [Dominant Tradeoff],
    [Mean Shift (\#34) @comaniciu2002mean], [Frozen-KDE EM shift], [$O(n^2 d T)$], [Linear (Gauss.)], [Euclid / No], [Exact stable modes; slow quadratic],
    [DENCLUE 1.0 (\#35) @hinneburg1998denclue], [Fixed-$delta$ norm. grad.], [$O(n log n + n d T)$], [Oscillates (Thm. 4)], [Euclid / Yes ($xi$)], [Fast $2 h$ cubes; $delta/2$ limit-cycle error],
    [DENCLUE 2.0 (\#43) @hinneburg2007denclue2], [Adaptive EM (equiv. MS)], [$O(n^2 d T_"avg")$], [Linear + monotonic], [Euclid / Yes ($xi$)], [Exact step-free EM + path merging],
    [Blurring MS (\#54) @carreira2006fast], [Dataset-updating $P^((t)) X^((t))$], [$O(T n^2 d)$], [Cubic (Thm. 5)], [Euclid / No], [3--6 iter.; collapses globally if unstopped],
    [Mean Shift++ (\#55) @jang2021meanshiftpp], [$3^d$ Moore grid blur], [$O(T n d 3^d)$], [Fast discrete], [Euclid / No], [Linear in $n$ for $d < log_3 n$; $3^d$ barrier],
    [Quickshift (\#41) @vedaldi2008quickshift], [Nearest higher-$hat(p)$ tree], [$O(n^2 d)$], [Non-iterative], [Metric / Partial], [No locking; sensitive to single cut $tau$],
    [Medoidshift (\#42) @sheikh2007medoidshift], [Discrete medoid $W D^((2))$], [$O(n^2 (d + n))$], [Finite MM], [Dissim. / No], [Zero $tau$ tuning; $O(n^3)$ + locking balls],
  ),
)

// =============================================================================
// SECTION III
// =============================================================================
#sec("III", "Density Peaks, Level-Set Trees, and Local Outlier Scoring")

Section III examines four modern Category II methods—DPC (Alg. 51, @rodriguez2014dpc), SNN-DPC (Alg. 52, @liu2018snndpc), DeBaCl (Alg. 46, @kent2013debacl), and OPTICS-OF (Alg. 53, @breunig1999opticsof)—that bridge Sections I and II by combining density estimation with graph-theoretic tree extraction, 2D decision-graph mode selection, and continuous local anomaly quantification.

#subsec("A", "Decision-Graph Geometry, Quickshift Duality, and SNN-DPC Scale Compensation")

A remarkable structural unification connects DPC (Alg. 51) directly to Quickshift (Alg. 41), while illuminating precisely why vanilla DPC fails on multi-density mixtures and how SNN-DPC (Alg. 52) eliminates those pathologies algebraically.

#theorem(title: "Exact Forest Isomorphism of DPC and Quickshift, Peak Masking, and Domino Propagation", origin: "Original Proof")[
  Let $(cal(X), d)$ be a finite metric space with pairwise distances $d_(i j)$ and Gaussian kernel density $rho_i^"gau" = sum_(j != i) exp(-(d_(i j) / d_c)^2)$.
  1. *(Exact Parent-Tree Isomorphism):* Setting Quickshift's bandwidth to $sigma = d_c / sqrt(2)$ (with identical tie-breaking order $pi$), the unpruned Quickshift directed tree $G_infinity$ (Alg. 41) and the DPC parent-pointer tree $T_"DPC" = (V, {(i, p(i)) : i != i^*})$ (Alg. 51) are *graph-theoretically identical*. Consequently, DPC and Quickshift differ exclusively in their *tree-cutting functional*: Quickshift cuts edges by a 1D threshold $delta_i > tau$, whereas DPC cuts at the top-$K$ vertices maximizing the 2D decision score $gamma_i = rho_i delta_i$ (or region $rho_i > rho_"min", delta_i > delta_"min"$).
  2. *(Multi-Density Peak Masking):* Suppose $cal(X) = C_1 union.sq C_2$ consists of a dense cluster $C_1$ (scale $alpha_1$, peak density $rho_1$, shoulder point $s in C_1$ with density $rho_s$ and intra-cluster separation $delta_s$) and a sparse cluster $C_2$ (scale $alpha_2 > alpha_1$, peak $c_2$ with density $rho_2$ and inter-cluster separation $delta_2 = D_(12)$). Whenever the density ratio satisfies
     $
       rho_s / rho_2 > D_(12) / delta_s,
     $
     the decision score ranks the non-peak shoulder $s in C_1$ strictly above the true cluster center $c_2$ ($gamma_s > gamma_(c_2)$), causing top-$K$ DPC to split $C_1$ while absorbing $C_2$ into $C_1$.
  3. *(Domino Subtree Error):* Under DPC's single-hop assignment $ell(i) = ell(p(i))$, reversing the parent pointer of a single saddle vertex $u$ misclassifies the entire descendant subtree $T_"DPC" (u)$, inflicting worst-case classification error $|T_"DPC" (u)| = Omega(n)$.
]

#proof[
  *Part 1:* Because $hat(p)_sigma (x_i) = (n (2 pi sigma^2)^(d/2))^(-1) (1 + rho_i^"gau")$ when $d_c = sqrt(2) sigma$, the mapping $rho_i^"gau" |-> hat(p)_sigma (x_i)$ is strictly affine-increasing, preserving the exact total order $pi$ on ${1, dots, n}$. Thus for every $i != i^*$, the upper set $cal(U)_i = {j : pi(j) < pi(i)}$ is identical, so $p(i) = arg min_(j in cal(U)_i) d_(i j) = pi_"QS" (i)$ and $delta_i = d_(i, p(i))$ coincide vertex-by-vertex. Deleting the outgoing arc of each selected center $c_m in cal(S)$ partitions $T_"DPC"$ into $K$ directed in-trees whose roots are $cal(S)$; Quickshift performs the identical edge-deletion operation with root set ${i : delta_i > tau}$.

  *Part 2:* Direct multiplication yields $gamma_s = rho_s delta_s > rho_2 D_(12) = gamma_(c_2)$ iff $rho_s / rho_2 > D_(12) / delta_s$. In Gaussian DPC with $d_c = alpha_1$ tuned to resolve $C_1$, a dense cluster of $m_1$ points at spacing $alpha_1$ yields $rho_s = Omega(m_1)$ and $delta_s = Theta(alpha_1)$, whereas a sparse cluster of spacing $alpha_2 >> alpha_1$ has $rho_2 <= m_2 exp(-(alpha_2 / alpha_1)^2)$. Thus $rho_s / rho_2 = Omega((m_1 / m_2) exp((alpha_2 / alpha_1)^2))$, which exponentially dwarfs $D_(12) / alpha_1$ even for moderate scale contrast $alpha_2 / alpha_1 >= 3$, masking $c_2$ completely.

  *Part 3:* By induction on depth in $T_"DPC"$, if $u in.not cal(S)$ and $p(u)$ is assigned label $ell_1$ instead of true label $ell_2$, then every vertex $v$ with a directed path $v -> dots -> u$ containing no center in $cal(S)$ inherits $ell(v) = ell(u) = ell_1$. On a crescent or density-bridge manifold whose low-density tail funnels through a single bottleneck $u$, $|T_"DPC" (u)| = Theta(n)$ points fail simultaneously.
]

SNN-DPC (Alg. 52, @liu2018snndpc) redesigns all three components of DPC: (i) local density $rho_i = sum_(j in L_k (i)) "Sim"(i, j)$ with $"Sim"(i, j) = (|"SNN"(i, j)|^2) / (sum_(p in "SNN"(i, j)) (d_(i p) + d_(j p)))$ on mutual $k$-NNs; (ii) *scale-compensated separation* $delta_i = min_(j : rho_j > rho_i) d_(i j) (S_i + S_j)$ where $S_i = sum_(p in Gamma(i)) d_(i p)$; and (iii) two-tier assignment via strong-overlap BFS ($|"SNN"(p, q)| >= k / 2$) followed by consensus $k$-NN voting. We now prove algebraically *why* SNN-DPC's seemingly heuristic factor $(S_i + S_j)$ achieves exact multi-scale invariance for sparse peaks.

#theorem(title: "Exact Scale-Compensation Law and Domino-Barrier Guarantee of SNN-DPC", origin: "Original Proof")[
  Let $cal(X) = C_1 union.sq C_2 subset.eq bb(R)^d$ consist of two isolated clusters ($|C_r| > k$), each generated as a homothetic copy $C_r = c_r + alpha_r T_0$ of a bounded base template $T_0$ with scales $0 < alpha_1 < alpha_2$ and inter-cluster separation $D_(12) := "dist"(C_1, C_2) > alpha_2 "diam"(T_0)$. Let $rho^((0))_i$ and $S^((0))_i$ denote the SNN-DPC densities and $k$-NN distance sums on $T_0$.
  1. *(Scale-Invariant Peak Dominance):* Every non-peak interior shoulder $s in C_1$ satisfies
     $
       gamma_s = alpha_1 dot rho^((0))_s d^((0))_(s, j(s)) (S^((0))_s + S^((0))_(j(s))) = Theta(alpha_1),
     $
     whereas the sparse cluster's true peak $c_2 in C_2$ satisfies
     $
       gamma_(c_2) >= rho^((0))_(c_2) S^((0))_(c_2) D_(12) = Omega(D_(12)),
     $
     *completely independent* of the sparse scale $alpha_2$. Consequently, the peak-to-shoulder ratio satisfies $gamma_(c_2) / gamma_s = Omega(D_(12) / alpha_1) -> infinity$ as $alpha_1 -> 0$, strictly defeating the masking pathology of Theorem 7.
  2. *(Inevitably Subordinate Domino Barrier):* Across any two clusters whose closed $k$-NN sets are disjoint ($Gamma(p) inter Gamma(q) = emptyset$ for $p in C_1, q in C_2$), every cross-cluster pair satisfies $|"SNN"(p, q)| = 0 < k / 2$; hence SNN-DPC's Step 1 BFS can *never* leak labels across the inter-cluster gap, confining any residual ambiguity strictly to consensus-voted Step 2 boundary points.
]

#proof[
  *Part 1:* Because $D_(12) > alpha_2 "diam"(T_0) > alpha_1 "diam"(T_0)$, every point $i in C_r$ has all $k$ nearest neighbors strictly inside $C_r$. Thus the combinatorial overlap $|"SNN"(i, j)|$ is invariant to $alpha_r$, while all intra-cluster distances scale linearly by $alpha_r$. Substituting into Definition 1 of Alg. 52 gives
  $
    "Sim"(i, j) &= (|"SNN"(i, j)|^2) / (alpha_r sum_(p in "SNN"(i, j)) (d^((0))_(i p) + d^((0))_(j p))) = alpha_r^(-1) "Sim"_0 (i, j), \
    rho_i &= alpha_r^(-1) rho^((0))_i, \
    S_i &= sum_(p in Gamma(i)) alpha_r d^((0))_(i p) = alpha_r S^((0))_i.
  $
  Consider any non-peak shoulder $s in C_1$. Because $C_1$ contains a strictly denser point $j in C_1$ at distance $d_(s j) <= alpha_1 "diam"(T_0) < D_(12)$ while any higher-density point in $C_2$ (if one exists) lies at distance $>= D_(12)$ with $S_j > S_s$, the minimizer $j(s) = arg min_(j : rho_j > rho_s) d_(s j)(S_s + S_j)$ lies in $C_1$ (or yields an even smaller product than an intra-$C_1$ minimizer). Hence
  $
    delta_s <= d_(s, j_1)(S_s + S_(j_1)) = alpha_1^2 d^((0))_(s, j_1)(S^((0))_s + S^((0))_(j_1)),
  $
  where $j_1 in C_1$ has $rho_(j_1) > rho_s$. Multiplying by $rho_s = alpha_1^(-1) rho^((0))_s$ yields
  $
    gamma_s = rho_s delta_s <= alpha_1 dot rho^((0))_s d^((0))_(s, j_1)(S^((0))_s + S^((0))_(j_1)) = O(alpha_1).
  $
  Now consider the unique peak $c_2 = arg max_(i in C_2) rho_i$. Since $alpha_1 < alpha_2$ and $max_(T_0) rho^((0))$ is achieved on $C_1$, $rho_(c_1) = alpha_1^(-1) rho^((0))_"max" > alpha_2^(-1) rho^((0))_"max" = rho_(c_2)$, so every $j$ with $rho_j > rho_(c_2)$ lies strictly in $C_1$ at distance $d_(c_2, j) >= D_(12)$. Therefore,
  $
    delta_(c_2) &= min_(j in C_1 : rho_j > rho_(c_2)) d_(c_2, j)(S_(c_2) + S_j) \
    &>= D_(12) (alpha_2 S^((0))_(c_2) + alpha_1 S^((0))_"min").
  $
  Multiplying by $rho_(c_2) = alpha_2^(-1) rho^((0))_(c_2)$ cancels $alpha_2$ identically:
  $
    gamma_(c_2) >= rho^((0))_(c_2) D_(12) (S^((0))_(c_2) + (alpha_1 / alpha_2) S^((0))_"min") >= rho^((0))_(c_2) S^((0))_(c_2) D_(12).
  $
  Dividing the lower bound on $gamma_(c_2)$ by the upper bound on $gamma_s$ gives $gamma_(c_2) / gamma_s = Omega(D_(12) / alpha_1)$, independent of $alpha_2$.

  *Part 2:* Step 1 BFS expands from $p$ to $q in Gamma(p)$ exclusively if $|"SNN"(p, q)| >= k / 2 >= 1$. For any $p in C_1$ and $q in C_2$ with disjoint closed $k$-NN sets, $|"SNN"(p, q)| = |Gamma(p) inter Gamma(q)| = 0 < k / 2$, so no cross-cluster edge can ever be traversed in Step 1.
]

Theorem 8 explains the remarkable empirical superiority of SNN-DPC over DPC across multi-density and manifold benchmarks @liu2018snndpc, while highlighting its exact price: $O(n^2 + n k^2)$ cost and explicit reliance on the user specifying the target peak count $K$ (or reading the decision graph manually).

#subsec("B", "Statistical Level-Set Trees vs. Decision Graphs: DeBaCl vs. HDBSCAN vs. DPC")

Whereas DPC and SNN-DPC construct a single flat mode forest requiring a manual choice of $K$, DeBaCl (Alg. 46, @kent2013debacl) directly estimates Hartigan's *population level-set cluster tree* $cal(T)_f = union.big_(lambda >= 0) pi_0 (L_f (lambda))$ @hartigan1981consistency. A crisp algebraic comparison reveals how DeBaCl's filtration differs structurally from both HDBSCAN (Alg. 32) and DPC (Alg. 51).

#theorem(title: "Vertex- vs. Edge-Filtration Duality and Finite-Sample Tree Consistency of DeBaCl", origin: "Classical — Cited")[
  Let $cal(X) = {x_1, dots, x_n}$ be i.i.d. draws from a continuous density $f$ on $bb(R)^d$ with compact support and well-separated modal components @chaudhuri2014consistent @kent2013debacl.
  1. *(Filtration Contrast with HDBSCAN):* DeBaCl estimates $cal(T)_f$ via a *vertex-weighted superlevel filtration* $G_lambda = G_k [{i : hat(f)_k (x_i) >= lambda}]$ on a *fixed* sparse graph $G_k$ (symmetric $k$-NN graph), whereas HDBSCAN estimates $cal(T)_f$ via an *edge-weighted sublevel filtration* $G_"mreach" (epsilon)$ on the complete graph where both vertex inclusion ($d_"core" (x_i) <= epsilon$) and metric edge radius ($d(x_i, x_j) <= epsilon$) co-vary with $epsilon = lambda^(-1/d)$.
  2. *(Hartigan Consistency Guarantee, Chaudhuri et al. @chaudhuri2014consistent / Kent et al. @kent2013debacl):* If $k = k_n$ satisfies $k_n / (d log n) -> infinity$ and $k_n / n -> 0$ (e.g., $k_n = C d log n$ for a sufficiently large constant $C$) and mass-pruning threshold $gamma_n = o(n)$ exceeds finite-sample ripple size, then with probability $1 - O(1 / n)$, the pruned DeBaCl tree $cal(P)_(gamma_n) (hat(cal(T)))$ is *Hartigan consistent*: for any two disjoint population modes $A, A' in cal(T)_f$ separated by a valley of depth $Delta f > 0$, the empirical tree recovers disjoint branches containing $A inter cal(X)$ and $A' inter cal(X)$ that merge strictly below the saddle level.
  3. *(Complexity Advantage):* With a spatial index for $k$-NN search, DeBaCl builds the complete level-set tree in $O(n k log n + n k alpha(n))$ time using Union--Find over $O(n k)$ edges, avoiding the $O(n^2)$ complete-graph computations of DPC, SNN-DPC, and brute-force HDBSCAN.
]

#proof[
  *Part 1:* In DeBaCl, the adjacency list $E(G_k)$ is frozen prior to the sweep; decreasing $lambda$ activates vertices in descending order of $hat(f)_k (x_i) = k / (n V_d r_k (x_i)^d)$ (equivalently ascending $r_k (x_i)$), uniting adjacent active vertices via Disjoint-Set Union whenever ${u, v} in E(G_k)$. Thus two active vertices $u, v$ with $r_k (u), r_k (v) <= epsilon$ connect iff $v in N_k (u) or u in N_k (v)$—which implies $d(u, v) <= max{r_k (u), r_k (v)} <= epsilon$! By contrast, HDBSCAN connects $u, v$ whenever $d(u, v) <= epsilon$ *even if* neither lies in the other's $k$-NN set (which can occur when $epsilon$ exceeds both $r_k (u)$ and $r_k (v)$, i.e., across a saddle at lower density). Hence DeBaCl's $k$-NN graph enforces a locally *tighter* topological skeleton at low density levels while remaining decoupled from the choice of density estimator $hat(f)$ (accepting KDE, $k$-NN, or adaptive estimators).

  *Part 2 (Sketch of @chaudhuri2014consistent):* By uniform Vapnik--Chervonenkis concentration over Euclidean balls of $bb(R)^d$, with probability $1 - delta$, every ball $B$ of true probability mass $p(B) >= C_0 (d log n + log(1/delta)) / n$ has empirical count within $O(sqrt(n p(B) d log n))$ of $n p(B)$. Setting $k_n >= C d log n$ guarantees that $hat(f)_k (x_i)$ concentrates uniformly around $f(x_i)$ within error $o(Delta f)$ on every level set above the saddle. Consequently, at any threshold $lambda$ inside the valley gap separating $A$ and $A'$, all sample points in the valley satisfy $hat(f)_k < lambda$ (preventing premature merges in $G_lambda$), while every pair of neighboring points inside $A$ lies within $r_(k_n)$ of one another (ensuring $G_lambda [A inter cal(X)]$ is connected). Spurious finite-sample bumps of cardinality $< gamma_n$ are contracted by iterative leaf pruning without altering the macro-topology of persistent modes.

  *Part 3:* Computing $N_k (i)$ via $k$-d / cover trees takes $O(n k log n)$ in low intrinsic dimension; sorting $n$ density values takes $O(n log n)$; and sweeping $|E(G_k)| <= n k$ edges through Union--Find with path compression and union by rank takes $O(n k alpha(n))$ time, where $alpha(dot)$ is the inverse Ackermann function.
]

*Head-to-Head Verdict (DeBaCl vs. HDBSCAN vs. DPC):* DeBaCl provides the strongest statistical tree-recovery guarantees and decouples density estimation from connectivity graph construction; HDBSCAN adds an automated, parameter-free *non-horizontal* extraction objective (EoM stability) whereas DeBaCl defaults to leaf extraction or user-selected cuts; DPC and SNN-DPC excel when a human analyst inspects the 2D decision graph to pick $K$ visually on moderate-sized metric datasets ($n <= 10^4$), with SNN-DPC strictly dominating vanilla DPC under heterogeneous densities.

#subsec("C", "Continuous Local Outlier Scoring vs. Binary Noise: OPTICS-OF")

Every density clusterer in Section I except OPTICS-OF (Alg. 53, @breunig1999opticsof) treats anomalies as an undifferentiated binary class ($cal(X)_"noise"$). OPTICS-OF augments the exact OPTICS traversal with the local reachability distance $"lrd"_(m, epsilon) (p) = |N_epsilon (p)|^(-1) sum_(o in N_epsilon (p)) d_"reach" (p, o)$ and scale-relative outlier factor $"OF"_(m, epsilon) (p) = |N_epsilon (p)|^(-1) sum_(o in N_epsilon (p)) ("lrd"_(m, epsilon) (p)) / ("lrd"_(m, epsilon) (o))$.

#theorem(title: "Scale-Relative Outlier Separation of OPTICS-OF vs. Absolute-Reachability Failure", origin: "Original Proof")[
  Let $cal(D) = C_1 union.sq C_2 union.sq {q} subset.eq bb(R)^d$ with $|C_1|, |C_2| >= m >= 2$, where $C_1$ is a dense homogeneous cluster of diameter $D_1$ whose points have intra-cluster $k$-NN reachabilities in $[r_1, (1 + eta) r_1]$ ($k = m - 1, eta in (0, 1/4)$), $C_2$ is a sparse homogeneous cluster with reachabilities in $[r_2, (1 + eta) r_2]$, and $q$ is a local anomaly placed at distance $Delta$ from $C_1$ such that
  $
    (1 + eta)(1 + 2 eta) r_1 < Delta < r_2 < "dist"(C_2, C_1 union {q}) - r_2.
  $
  Choose generating radius $epsilon >= max{r_2 (1 + eta), Delta + D_1}$ so every point has $|N_epsilon (p)| >= k$, while $N_epsilon (q) = C_1$, $N_epsilon (x) subset.eq C_2$ for all $x in C_2$, and $|C_1| >= (Delta + D_1) / (eta r_1)$. Then:
  1. *(Absolute-Reachability Failure of OPTICS & DBSCAN):* The OPTICS reachability of the anomaly $q$ satisfies $r(q) <= Delta < r_2 <= r(x)$ for every $x in C_2$; hence *no* global reachability threshold (and no single DBSCAN radius) can flag $q$ as an outlier without destroying the entire sparse cluster $C_2$.
  2. *(Exact Local Scale Separation of OPTICS-OF):* Every interior point $x in C_2$ satisfies $(1 + eta)^(-1) <= "OF"_(m, epsilon) (x) <= 1 + eta$, whereas the anomaly $q$ satisfies
     $
       "OF"_(m, epsilon) (q) >= Delta / ((1 + 2 eta) r_1) > 1 + eta >= max_(x in C_2) "OF"_(m, epsilon) (x),
     $
     strictly separating $q$ from *both* $C_1$ and $C_2$ via a single dimensionless threshold $tau in (1 + eta, Delta / ((1 + 2 eta) r_1))$.
]

#proof[
  *Part 1:* Because $q$ lies at distance $Delta$ from its nearest neighbor $o^* in C_1$ and $d_"core" (o^*) <= (1 + eta) r_1 < Delta$, once $o^*$ is processed, OPTICS sets $r(q) <= d_"reach" (q, o^*) = max{d_"core" (o^*), d(q, o^*)} = Delta$. Conversely, every point $x in C_2$ (after the first seed of $C_2$, whose reachability exceeds $r_2$) has all neighbors in $C_2$ at reachability at least $r_2$, so $r(x) >= r_2 > Delta >= r(q)$. Any cutoff rejecting $q$ must therefore lie below $Delta < r_2$, discarding all of $C_2$.

  *Part 2:* For every $x in C_2$, every $o in N_epsilon (x) subset.eq C_2$ has $d_"reach" (x, o) in [r_2, (1 + eta) r_2]$, so $"lrd"_(m, epsilon) (x) in [r_2, (1 + eta) r_2]$. Taking the ratio across $o in N_epsilon (x)$ yields
  $
    1 / (1 + eta) = r_2 / ((1 + eta) r_2) <= ("lrd"_(m, epsilon) (x)) / ("lrd"_(m, epsilon) (o)) <= ((1 + eta) r_2) / r_2 = 1 + eta,
  $
  and averaging over $N_epsilon (x)$ preserves $(1 + eta)^(-1) <= "OF"_(m, epsilon) (x) <= 1 + eta$. For the anomaly $q$, every neighbor $o in N_epsilon (q) = C_1$ lies at distance $d(q, o) >= Delta$, so $"lrd"_(m, epsilon) (q) >= Delta$. Meanwhile, each $o in C_1$ has $|C_1| - 1$ intra-cluster neighbors at reachability $<= (1 + eta) r_1$ plus the single anomaly $q$ with $d_"reach" (o, q) = max{d_"core" (q), d(o, q)} <= Delta + D_1$. Averaging over $|N_epsilon (o)| = |C_1|$ neighbors and using $|C_1| >= (Delta + D_1) / (eta r_1)$ yields
  $
    "lrd"_(m, epsilon) (o) <= ((|C_1| - 1)(1 + eta) r_1 + Delta + D_1) / (|C_1|) <= (1 + 2 eta) r_1.
  $
  Substituting both bounds into the OPTICS-OF score gives
  $
    "OF"_(m, epsilon) (q) = 1 / (|C_1|) sum_(o in C_1) ("lrd"_(m, epsilon) (q)) / ("lrd"_(m, epsilon) (o)) >= Delta / ((1 + 2 eta) r_1) > 1 + eta,
  $
  strictly above $max_(x in C_2) "OF"_(m, epsilon) (x)$.
]

Theorem 10 identifies both the triumph and the limitation of OPTICS-OF @breunig1999opticsof: it delivers joint multi-scale clustering *and* continuous scale-relative outlier ranking in a single priority-queue pass, yet requires the generating radius $epsilon$ to be set large enough that *every* candidate anomaly still satisfies $|N_epsilon (q)| >= k$ (otherwise $"lrd"$ becomes undefined or heavily truncated)—a boundary condition later eliminated by pure $k$-NN Local Outlier Factor (LOF).

#wide_figure(
  caption: [*Table III:* Head-to-Head Comparison of Density Peaks, Level-Set Tree, and Local Outlier Algorithms (Section III).],
  table(
    columns: (1.05fr, 1.02fr, 0.90fr, 0.90fr, 0.95fr, 1.18fr),
    align: (left, left, left, center, center, left),
    [Algorithm], [Core Representation], [Complexity], [Auto-$K$ / Hierarchy], [Outlier Model], [Key Theoretical Distinction],
    [DeBaCl (\#46) @kent2013debacl], [Vertex-filtered $G_k$ tree], [$O(n k log n)$], [Yes (full tree)], [Background / leaves], [Provably Hartigan consistent (Thm. 9)],
    [DPC (\#51) @rodriguez2014dpc], [$(rho_i, delta_i)$ parent forest], [$O(n^2)$], [No (manual $K$)], [Border halo $rho_b (k)$], [Isomorphic to Quickshift; masks sparse peaks],
    [SNN-DPC (\#52) @liu2018snndpc], [SNN $rho_i$ + scaled $delta_i$], [$O(n^2 + n k^2)$], [No (manual $K$)], [None (full partition)], [Exact $alpha_2$-invariance + domino barrier (Thm. 8)],
    [OPTICS-OF (\#53) @breunig1999opticsof], [Reachability + $"OF"$ score], [$O(n log n)$--$O(n^2)$], [Yes (valley plot)], [Continuous $"OF" > 1$], [Separates local anomalies in multi-density (Thm. 10)],
  ),
)

// =============================================================================
// SECTION IV
// =============================================================================
#sec("IV", "Parallel, Grid, and Subsampled Scalability Architectures")

Sequential DBSCAN (Alg. 31) hits two hard computational ceilings: (i) single-machine memory and sequential BFS traversal once $n >= 10^7$, and (ii) spatial-index degeneration to $Theta(n^2 d)$ once $d >= 15$ @schubert2017dbscanrevisited. Five algorithms in Category II attack these bottlenecks along two orthogonal axes:
- *Exact Global-DBSCAN Preservation:* PDBSCAN (Alg. 37, @xu1999pdbscan), GriDBSCAN (Alg. 48, @mahran2008gridbscan), and MR-DBSCAN (Alg. 49, @he2011mrdbscan) guarantee *zero approximation error* on global core clusters while parallelizing local range queries across $P$ processors or $M$ grid cells.
- *Subquadratic Approximate / Hybrid Acceleration:* DBSCAN++ (Alg. 47, @jang2019dbscanpp) and K-DBSCAN (Alg. 44) break the $Theta(n^2 d)$ barrier even in moderate-to-high dimensions by sacrificing exact global DBSCAN equivalence in exchange for provable covering-radius or Voronoi-partitioned speedup.

#subsec("A", "Exact Core-Graph Equivalence and Halo-Replication Tradeoffs: PDBSCAN, GriDBSCAN, and MR-DBSCAN")

Despite originating across three eras of parallel computing (shared-nothing distributed databases in 1999, multi-core spatial grids in 2008, and Hadoop/MapReduce clusters in 2011), PDBSCAN, GriDBSCAN, and MR-DBSCAN all instantiate a single universal theorem of *halo-augmented core-bridge merging*.

#theorem(title: "Unified Exact Core-Graph Equivalence, Core-Gating Necessity, and Optimal Grid Halo Barrier", origin: "Original Unified Proof")[
  Let $cal(X) subset.eq Omega = product_(j=1)^d [a_j, b_j]$ be partitioned into $P$ interior-disjoint regions $cal(R) = {R_1, dots, R_P}$ with unique home ownership $h(x) in {1, dots, P}$, and define for each $r$ the *$epsilon$-extended region* $R_r^(+epsilon) := {y in Omega : "dist"(y, R_r) <= epsilon}$ and augmented dataset $tilde(cal(X))_r := cal(X) inter R_r^(+epsilon)$. Execute exact sequential DBSCAN$(epsilon, m)$ independently on every $tilde(cal(X))_r$, obtaining local core indicators $"core"_r (x)$ and local clusters ${C_(r, c)}$.
  1. *(Universal Exact Core-Cluster Recovery):* Define the undirected merge graph $cal(H)_"core"$ on local clusters ${(r, c)}$ with an edge between $(r, c)$ and $(r', c')$ iff there exists a shared point $z in C_(r, c) inter C_(r', c')$ (or a witness pair $(u, v) in C_(r, c) times C_(r', c')$ with $d(u, v) <= epsilon$) such that $z$ (resp. both $u, v$) is a core point in its home partition. Then the connected components of $cal(H)_"core"$ induce the *exact same* global core set $V_"core"$ and global core clusters as sequential DBSCAN$(cal(X), epsilon, m)$.
  2. *(Necessity of Core-Gating & Multi-Label Halos):* Merging $(r, c)$ and $(r', c')$ whenever they share an arbitrary *non-core border point* $b$ causes false fusion of globally disjoint DBSCAN clusters; conversely, storing only a *single* local cluster ID per halo point $y in tilde(cal(X))_r backslash cal(X)_r$ (as in naive GriDBSCAN) fails to merge globally connected core components separated by a narrow cell boundary.
  3. *(Optimal GriDBSCAN Width & Exponential Halo Barrier):* For $n$ points uniformly distributed in $[0, L]^d$ (density $lambda = n / L^d$) partitioned by GriDBSCAN into uniform cubes of side $w >= epsilon$, the expected total brute-force pairwise work $W(w) = (L / w)^d dot lambda^2 (w + 2 epsilon)^(2 d)$ is *uniquely minimized* at $w^* = 2 epsilon$, yielding total work $O(n lambda (8 epsilon)^d)$ and an irreducible *dimensional halo replication factor*
     $
       "Rep"(w^*) = ((w^* + 2 epsilon) / w^*)^d = 2^d,
     $
     which explodes exponentially with $d$ (and degrades to $Omega(n^2)$ local work under skewed data unless replaced by MR-DBSCAN's cost-balanced PRBP cuts).
]

#proof[
  *Part 1:* Take any interior point $x in cal(X) inter R_r$. For every $y in cal(X)$ with $d(x, y) <= epsilon$, $"dist"(y, R_r) <= d(x, y) <= epsilon$, so $y in tilde(cal(X))_r$. Hence $N_epsilon^(r) (x) = N_epsilon (x)$, implying $"core"_r (x) = "core"(x)$ identically for all $x in cal(X) inter R_r$, so $union.big_(r=1)^P {x in R_r : "core"_r (x) = 1} = V_"core"$. Next, consider any global core edge ${u, v} in E_"core"$ ($u, v in V_"core"$, $d(u, v) <= epsilon$). Let $r = h(u)$ and $r' = h(v)$. Because $d(u, v) <= epsilon$, we have $v in tilde(cal(X))_r$ and $u in tilde(cal(X))_(r')$. Inside reducer $r$, $u$ is an interior core point ($"core"_r (u) = 1$) adjacent to $v$, so $v$ belongs to $u$'s local cluster $C_(r, c_u)$ (even if $v$ appears non-core locally in $r$ due to missing neighbors beyond $R_r^(+epsilon)$). Meanwhile in $r'$, $v$ is an interior core point ($"core"_(r') (v) = 1$) in $C_(r', c_v)$. Thus $v in C_(r, c_u) inter C_(r', c_v)$ with $"core"_(r') (v) = 1$, inserting the exact bridge edge $((r, c_u), (r', c_v)) in cal(H)_"core"$. Conversely, every edge of $cal(H)_"core"$ is witnessed by a global core point $z$ (since $"core"_r (z) = 1 => |N_epsilon (z)| >= |N_epsilon^(r) (z)| >= m => z in V_"core"$) lying within $epsilon$ of core members of both $C_(r, c)$ and $C_(r', c')$, certifying a valid global core path. Thus the connected components of $cal(H)_"core"$ match $pi_0 (G_"core")$ identically.

  *Part 2:* For false fusion without core-gating, place two disjoint global core clusters $K_1 subset.eq R_1$ and $K_2 subset.eq R_2$ separated by $1.6 epsilon$, and place a single non-core border point $b$ on the interface $partial R_1 inter partial R_2$ at distance $0.8 epsilon$ from both $K_1$ and $K_2$. Then $b$ is absorbed as a border point by $K_1$ in $tilde(cal(X))_1$ and by $K_2$ in $tilde(cal(X))_2$; merging on shared non-core $b$ falsely fuses $K_1$ and $K_2$. For single-label halo failure, place two locally disconnected pieces $K_a, K_b subset.eq R_1$ whose only global bridge is a core point $y in R_2$ lying within $epsilon$ of both $K_a$ and $K_b$ (so $y$ is in the halo $H_1$ of $R_1$); because $y$ may lack sufficient neighbors inside $tilde(cal(X))_1$ to be locally core in $R_1$, local DBSCAN on $tilde(cal(X))_1$ assigns $y$ as a *border point* to whichever of $K_a, K_b$ visits $y$ first, dropping the second cluster ID unless the full set $cal(R)_1 (y)$ of touching local core clusters is recorded!

  *Part 3:* With $M = (L / w)^d$ grid cells, each $epsilon$-haloed box $tilde(G)_g$ has side $w + 2 epsilon$ and expected cardinality $lambda (w + 2 epsilon)^d$. Approximating second moment by squared mean gives
  $
    W(w) = M dot lambda^2 (w + 2 epsilon)^(2 d) = lambda^2 L^d ( ((w + 2 epsilon)^2) / w )^d.
  $
  Let $phi(w) := (w + 2 epsilon)^2 / w = w + 4 epsilon + 4 epsilon^2 / w$ for $w >= epsilon$. Differentiating yields $phi'(w) = 1 - 4 epsilon^2 / w^2 = 0 <=> w^* = 2 epsilon$, with $phi''(w^*) = 1 / epsilon > 0$, proving $w^* = 2 epsilon$ is the unique global minimizer with $phi(w^*) = 8 epsilon$ and $W(w^*) = n lambda (8 epsilon)^d$. The expected replication factor per point is
  $
    "Rep"(w) = (M dot lambda (w + 2 epsilon)^d) / n = ((w + 2 epsilon) / w)^d,
  $
  which evaluates to $2^d$ at $w^* = 2 epsilon$ (and $3^d$ at $w = epsilon$).
]

*Architectural Comparison of the Three Exact Parallel Engines:*
- *PDBSCAN (Alg. 37, @xu1999pdbscan)* uses a shared-nothing *distributed R\*-tree (dR\*-tree)* over Hilbert-ordered shards, avoiding upfront replication of entire margins by issuing remote RPCs only when a core query $B_epsilon (x)$ crosses $partial R_r$. This saves memory when boundaries are sparse, achieving $O((n / P) log n)$ local query time, but couples worker execution via fine-grained synchronous network RPCs that stall on stragglers and fail under commodity node faults.
- *GriDBSCAN (Alg. 48, @mahran2008gridbscan)* replaces pointer-chasing dR\*-trees with embarrassingly parallel uniform $2 epsilon$-cells (Theorem 11.3), achieving ultra-fast single-node multi-core execution in $d in {2, 3, 4}$ on uniform spatial point clouds, yet collapses both under high dimension ($2^d$ halo blowup) and under severe spatial skew (one dense urban cell receiving $Theta(n)$ points).
- *MR-DBSCAN (Alg. 49, @he2011mrdbscan)* resolves both PDBSCAN's RPC coupling and GriDBSCAN's skew bottleneck within three fault-tolerant MapReduce jobs: PRBP bisects a coarse grid census to minimize $|Phi(R_L) - Phi(R_R)| + lambda B(R; j, t)$, balancing reducer compute work $Phi(R_r) approx Phi(Omega) / P$ while explicitly penalizing cross-cut margin replication $B(R; j, t)$, achieving $O((n / P) log(n / P) + m_"ov" alpha(m_"ov"))$ wall-clock time at web scale (still bounded by $d <= 10$ due to local R\*-tree and $epsilon$-margin geometry).

#subsec("B", "Subquadratic Approximation Tradeoffs: DBSCAN++ vs. K-DBSCAN")

When dimension $d$ exceeds $10$--$15$, spatial indices and $epsilon$-grid halos all suffer exponential degradation ($2^d$ replication or near-linear tree scans). DBSCAN++ (Alg. 47, @jang2019dbscanpp) and K-DBSCAN (Alg. 44) break $Theta(n^2 d)$ arithmetic directly in *arbitrary Euclidean dimension* $d$, trading exact DBSCAN equivalence for sharply quantifiable geometric guarantees.

#theorem(title: "DBSCAN++ Covering-Radius Preservation vs. K-DBSCAN Speedup-Fracture Tradeoff", origin: "Original Proof")[
  Let $cal(X) subset.eq bb(R)^d$ ($|cal(X)| = n$) with Euclidean metric.
  1. *(DBSCAN++ Deterministic Core-Preservation Bound):* Let $S subset.eq cal(X)$ of size $m = p n << n$ be a $delta$-cover of the global DBSCAN core set $V_"core" (epsilon)$ (attained within factor $2$ of optimal covering radius by $O(m n d)$ greedy $k$-center seeding). Running DBSCAN++ (Alg. 47) with density radius $epsilon_d = epsilon + delta$ and core-connectivity radius $epsilon_c = epsilon + 2 delta$ guarantees:
     - *(Zero False Splits):* Every connected component $K_a$ of $G_"core" (epsilon)$ maps into a *single* connected component of the subsampled core graph $G_C (epsilon_c)$, and every point $x in K_a$ is assigned to that component.
     - *(Zero False Merges):* If all global core components satisfy pairwise separation $"dist"(K_a, K_b) > epsilon + 4 delta$ and no spurious points outside $union.big_a B(K_a, delta)$ satisfy $|B(s, epsilon + delta) inter cal(X)| >= m_"pts"$, then distinct global components remain strictly disjoint in $G_C (epsilon_c)$.
     Conversely, under *uniform random subsampling* of size $m = p n$, a minority core cluster of size $n_"min"$ is missed completely with probability $(1 - n_"min" / n)^m approx exp(-p n_"min")$.
  2. *(K-DBSCAN Optimal Speedup & Voronoi Bisector Fracture):* Partitioning $cal(X)$ into $k$ balanced Voronoi cells via $T$ Lloyd iterations followed by intra-cell brute-force DBSCAN reduces arithmetic complexity from $Theta(n^2 d)$ to $O(n k d T + (n^2 / k) d)$, uniquely minimized at
     $
       k^* = Theta(sqrt(n / T)) quad => quad T_"K-DBSCAN" = O(n^(3/2) d sqrt(T)).
     $
     However, without homogeneous $epsilon$-collar reconciliation, *every* connected core manifold intersecting a Voronoi bisector hyperplane $H_(r s)$ of width exceeding $epsilon_r$ is deterministically fractured into disconnected subclusters; conversely, adding $epsilon$-collar merging restores topological continuity only by forfeiting K-DBSCAN's heterogeneous local-density radii ${epsilon_r}_(r=1)^k$.
]

#proof[
  *Part 1:* Because $S$ is a $delta$-cover of $V_"core" (epsilon)$, for every $u in V_"core" (epsilon)$ there exists $s(u) in S$ with $norm(u - s(u))_2 <= delta$. By the triangle inequality, $B(u, epsilon) subset.eq B(s(u), epsilon + delta)$, whence
  $
    |B(s(u), epsilon_d) inter cal(X)| >= |B(u, epsilon) inter cal(X)| >= "MinPts",
  $
  certifying $s(u) in C_S$. Now let ${u, v} in E_"core" (epsilon)$ be any edge inside a global core component $K_a$ (so $norm(u - v)_2 <= epsilon$). By the triangle inequality,
  $
    norm(s(u) - s(v))_2 <= norm(s(u) - u)_2 + norm(u - v)_2 + norm(v - s(v))_2 <= epsilon + 2 delta = epsilon_c,
  $
  so ${s(u), s(v)} in E_C$. Inducting over paths in $K_a$ shows that ${s(u) : u in K_a}$ lies inside a single connected component of $G_C (epsilon_c)$; moreover, every $u in K_a$ satisfies $norm(u - s(u))_2 <= delta <= epsilon_c$, so nearest-core assignment attaches $u$ to that component. For separation, any $s_a in B(K_a, delta)$ and $s_b in B(K_b, delta)$ satisfy $norm(s_a - s_b)_2 >= "dist"(K_a, K_b) - 2 delta > epsilon + 2 delta = epsilon_c$, preventing any cross-component edge in $G_C (epsilon_c)$. Finally, under uniform sampling without replacement, the exact miss probability for a minority subset of size $n_"min"$ is
  $
    product_(j=0)^(m-1) (n - n_"min" - j) / (n - j) <= (1 - n_"min" / n)^m <= exp(- p n_"min").
  $

  *Part 2:* Stage 1 Lloyd clustering costs $O(n k d T)$. In Stage 2, if cell sizes are balanced ($n_r = Theta(n / k)$), intra-cell distance computation across $k$ cells costs $sum_(r=1)^k O(n_r^2 d) = O((n^2 / k) d)$. Minimizing $n k d T + (n^2 / k) d$ with respect to $k$ sets $n d T = (n^2 / k^2) d <=> k^* = sqrt(n / T)$, yielding $O(n^(3/2) d sqrt(T))$. For the fracture statement, consider two concentric circular rings of radii $R_1 < R_2$ centered at the origin with $k >= 2$: every Lloyd centroid partition $V_1, dots, V_k$ cuts the outer ring across at least two radial Voronoi bisector rays $H_(r s)$. Because Stage 2 executes DBSCAN strictly inside disjoint $P_r = cal(X) inter V_r$, zero edges cross $H_(r s)$, deterministically shattering the connected ring into at least $k$ arcs. Reconciling across $epsilon$-collars $tilde(P)_r inter tilde(P)_s$ requires a *single global* radius $epsilon$ so core predicates agree on both sides of $H_(r s)$ (Theorem 11.1), negating cell-specific local radii $epsilon_r$.
]

*Head-to-Head Scalability Synthesis:* Theorem 12 reveals why DBSCAN++ with $k$-center seeding @jang2019dbscanpp is generally preferable to disjoint K-DBSCAN on *non-convex* high-dimensional manifolds: DBSCAN++ evaluates density against the *entire* dataset $cal(X)$ (avoiding Voronoi boundary truncation!) at $O(m n d)$ cost—matching minimax-optimal Hausdorff level-set recovery at sublinear $m = O(n^(d/(d+4)))$ @jang2019dbscanpp—whereas K-DBSCAN shines specifically when clusters are *spatially well-separated macro-regions* of sharply different local densities where Stage 1 Voronoi cells naturally isolate distinct density regimes so each cell's quantile radius $epsilon_r$ adapts without cutting across manifolds!

#wide_figure(
  caption: [*Table IV:* Scalability, Exactness, Communication, and Dimensional Scaling Comparison (Section IV + Baselines).],
  table(
    columns: (1.10fr, 1.05fr, 0.90fr, 0.85fr, 0.95fr, 1.15fr),
    align: (left, left, center, center, center, left),
    [Algorithm], [Parallel / Accel. Mechanism], [Exact DBSCAN?], [Arithmetic Cost], [High-$d$ ($d >= 20$)], [Communication / Approximation Bottleneck],
    [DBSCAN (\#31) @ester1996dbscan], [Sequential spatial index], [Baseline], [$O(n log n)$ / $O(n^2)$], [Degrades $O(n^2)$], [Single-node sequential BFS],
    [PDBSCAN (\#37) @xu1999pdbscan], [Shared-nothing dR\*-tree], [Yes (cores)], [$O((n / P) log n)$], [Degrades], [Fine-grained cross-shard RPC latency],
    [GriDBSCAN (\#48) @mahran2008gridbscan], [Uniform $2 epsilon$ halo grid], [Yes (multi-lbl)], [$O(n lambda (8 epsilon)^d)$], [Fails ($2^d$ halo)], [$2^d$ replication + skew load imbalance],
    [MR-DBSCAN (\#49) @he2011mrdbscan], [PRBP cost cuts + 3 MR jobs], [Yes (cores)], [$O((n/P) log(n/P))$], [Degrades ($d > 10$)], [Shuffle volume $O(n + sum |M_r|)$],
    [DBSCAN++ (\#47) @jang2019dbscanpp], [$k$-center / unif. subsample], [Approx. (Thm. 12)], [$O(m n d)$], [Strong ($m << n$)], [$delta$-cover inflation; uniform misses small],
    [K-DBSCAN (\#44)], [Voronoi Lloyd + local DBSCAN], [Approx. (Thm. 12)], [$O(n^(3/2) d sqrt(T))$], [Strong], [Bisector manifold fracture (no collar)],
    [Mean Shift++ (\#55) @jang2021meanshiftpp], [$3^d$ hashed grid blur], [Approx. MS], [$O(T n d 3^d)$], [Fails ($d > log_3 n$)], [$3^d$ neighbor lookups exceed $n$],
  ),
)

// =============================================================================
// SECTION V
// =============================================================================
#sec("V", "Complete 25-Algorithm Comparison Matrix and Decision Guide")

To synthesize all theoretical derivations from Sections I--IV into an actionable reference, Table V presents a complete head-to-head matrix across all twenty-five Category II algorithms (Algs. 31--55), followed by a mathematically grounded selection architecture.

#wide_figure(
  caption: [*Table V:* Master Head-to-Head Theoretical and Practical Comparison Matrix of All 25 Density and Mode-Seeking Algorithms (Category II, Algorithms 31--55). Abbreviations: *FD* = Fixed-Density Reachability, *VD* = Variable-Density / Multi-Scale, *KM* = Kernel Mode-Seeking, *DP/LT* = Density Peaks / Level-Set Trees, *SC* = Scalable / Distributed / Hybrid; *Bin.* = Binary Noise, *Cont.* = Continuous Outlier Score, *Halo* = Border Halo Filtering.],
  {
    set text(size: 6.9pt)
    set table(inset: (x: 3.2pt, y: 2.3pt))
    table(
      columns: (1.18fr, 0.44fr, 1.18fr, 0.94fr, 0.52fr, 0.56fr, 0.66fr, 0.60fr),
      align: (left, center, left, left, center, center, center, center),
      [Algorithm], [Grp.], [Key Parameters], [Time Complexity], [Multi-Dens.], [Noise / Out.], [High-$d$], [Space],
      [\#31 DBSCAN @ester1996dbscan], [FD], [$epsilon, "MinPts"$], [$O(n log n) "/" O(n^2)$], [No], [Bin.], [Low], [Metric],
      [\#32 HDBSCAN @campello2015hierarchical], [VD], [$m_"pts", m_"cl"$], [$O(n^2) "/" O(n log n)$], [Yes], [Bin. + prob], [Med], [Metric],
      [\#33 OPTICS @ankerst1999optics], [VD], [$epsilon_"max", "MinPts", xi$], [$O(n log n) "/" O(n^2)$], [Yes], [Bin.], [Low], [Metric],
      [\#34 Mean Shift @comaniciu2002mean], [KM], [Bandwidth $h$, kernel], [$O(n^2 d T)$], [No], [None], [Low], [Euclid.],
      [\#35 DENCLUE 1.0 @hinneburg1998denclue], [KM], [$sigma, xi, delta$], [$O(n log n + n d T)$], [No], [Bin. ($xi$)], [Low], [Euclid.],
      [\#36 SNN @ertoz2003snn], [VD], [$k, "Eps", "MinPts"$], [$O(n^2 d + n k^2)$], [Yes], [Bin.], [High], [Metric],
      [\#37 PDBSCAN @xu1999pdbscan], [SC], [$epsilon, "MinPts", P$], [$O((n / P) log n)$], [No], [Bin.], [Low], [Euclid.],
      [\#38 DBCLASD @xu1998dbclasd], [VD], [None (internal $alpha$)], [$O(n log n)$], [Partial], [Bin.], [Low], [Euclid.],
      [\#39 GDBSCAN @sander1998gdbscan], [FD], [$"NPred", "wCard", "MinCard"$], [$O(n T_"pred")$], [Custom], [Bin.], [Custom], [General],
      [\#40 ST-DBSCAN @birant2007stdbscan], [FD], [$epsilon_1, epsilon_2, "MinPts", Delta epsilon$], [$O(n log n + n W d_S)$], [No], [Bin.], [Low], [Space--Time],
      [\#41 Quickshift @vedaldi2008quickshift], [KM], [Bandwidth $sigma$, cut $tau$], [$O(n^2 d)$], [Partial], [Partial], [Med], [Metric],
      [\#42 Medoidshift @sheikh2007medoidshift], [KM], [Bandwidth $sigma$], [$O(n^2 (d + n))$], [No], [None], [Low], [Dissim.],
      [\#43 DENCLUE 2.0 @hinneburg2007denclue2], [KM], [$h, xi, tau_m$], [$O(n^2 d T_"avg")$], [No], [Bin. ($xi$)], [Low], [Euclid.],
      [\#44 K-DBSCAN], [SC], [$k, q, "minPts"$], [$O(n k d T + sum n_r^2 d)$], [Regional], [Bin.], [Med], [Euclid.],
      [\#45 RNN-DBSCAN @bryant2018rnndbscan], [VD], [$k$ (single integer)], [$O(n^2 d)$], [Yes], [Bin.], [Med], [Metric],
      [\#46 DeBaCl @kent2013debacl], [DP/LT], [$k$, pruning mass $gamma$], [$O(n k log n)$], [Yes], [Bg / leaves], [Med], [Euclid.],
      [\#47 DBSCAN++ @jang2019dbscanpp], [SC], [$m, epsilon_d, epsilon_c, "minPts"$], [$O(m n d)$], [No], [Bin.], [High], [Euclid.],
      [\#48 GriDBSCAN @mahran2008gridbscan], [SC], [$epsilon, "minPts", w_r$], [$O(n d + sum |tilde(P)_g|^2 d)$], [No], [Bin.], [Low], [Euclid.],
      [\#49 MR-DBSCAN @he2011mrdbscan], [SC], [$epsilon, "MinPts", P, lambda$], [$O((n/P) log(n/P))$], [No], [Bin.], [Low], [Euclid.],
      [\#50 VDBSCAN @liu2007vdbscan], [VD], [$k$, knee $tau$ (or $L$)], [$O(n^2 + L n)$], [Discrete], [Bin.], [Low], [Metric],
      [\#51 DPC @rodriguez2014dpc], [DP/LT], [Cutoff $d_c$, count $K$], [$O(n^2)$], [No], [Halo], [Med], [Dissim.],
      [\#52 SNN-DPC @liu2018snndpc], [DP/LT], [Neighbor $k$, count $K$], [$O(n^2 + n k^2)$], [Yes], [None], [High], [Metric],
      [\#53 OPTICS-OF @breunig1999opticsof], [DP/LT], [$epsilon, "MinPts", tau$], [$O(n log n) "/" O(n^2)$], [Yes], [Cont. ($"OF"$)], [Low], [Metric],
      [\#54 Blurring MS @carreira2006fast], [KM], [$h$, stopping tol.], [$O(T n^2 d)$], [No], [None], [Low], [Euclid.],
      [\#55 Mean Shift++ @jang2021meanshiftpp], [KM], [Grid width $h$], [$O(T n d 3^d)$], [No], [None], [Low ($d <= 8$)], [Euclid.],
    )
  },
)

#subsec("A", "Mathematically Grounded Practitioner's Decision Guide")

Selecting among Algorithms 31--55 should follow the provable theoretical invariants established in Theorems 1--12 rather than ad-hoc trial-and-error. We organize engineering selection into seven mutually informative operational regimes:

1. *Regime 1 — Homogeneous-Density Spatial Manifolds ($d <= 6$, $n <= 10^6$):* Use *DBSCAN (\#31)* with an R\*-tree or $k$-d tree ($O(n log n)$) whenever clusters share a common density scale ($overline(tau)_m (C_r) < "dist"(C_r, C_s)$ for all $r != s$). If exact multi-core speed is needed on uniform spatial data, upgrade to *GriDBSCAN (\#48)* at its provably optimal cell width $w^* = 2 epsilon$ (Theorem 11.3); if data exhibit severe urban/astronomical skew at cluster scale ($n >= 10^7$), deploy *MR-DBSCAN (\#49)* over *PDBSCAN (\#37)* to eliminate RPC stragglers while preserving exact core clusters.
2. *Regime 2 — Heterogeneous Multi-Density Clusters Without Known $K$ ($d <= 20$):* Fixed-$epsilon$ DBSCAN is provably impossible here (Theorem 1). Choose *HDBSCAN (\#32)* as the default workhorse: it extracts non-horizontal optimal cuts via Excess-of-Mass stability (Theorem 2.3), requires only two intuitive size integers $(m_"pts", m_"cl")$, and produces cluster membership probabilities. Choose *DeBaCl (\#46)* when formal statistical estimation of Hartigan's cluster tree or custom KDE inputs are required (Theorem 9), *OPTICS (\#33)* when visual interactive inspection of the 1D reachability plot across all scales is desired, or *RNN-DBSCAN (\#45)* when tuning budget allows only a single scale-free integer $k$ (Theorem 3). Reserve *VDBSCAN (\#50)* strictly for datasets with sharp, discrete $k$-dist density plateaus.
3. *Regime 3 — High-Dimensional Embeddings, Text, and Genomics ($d >= 30$):* Metric balls and grids fail simultaneously via distance concentration @beyer1999nearest and exponential halo explosion ($2^d$ in GriDBSCAN, $3^d$ in Mean Shift++). If cluster count $K$ is known or visually selectable from a decision graph, select *SNN-DPC (\#52)* for exact multi-scale peak invariance and domino protection (Theorem 8); if $K$ is unknown and noise pruning is required, select *SNN (\#36)*, whose mutual-overlap gate caps hub degree at $k$ (Theorem 3.3); if uniform-density scalability at large $n$ is the sole bottleneck, use *DBSCAN++ (\#47)* with greedy $k$-center seeding ($O(m n d)$, Theorem 12.1) rather than uniform subsampling.
4. *Regime 4 — Continuous Mode-Seeking, Image Segmentation, and Signal Filtering:* For low-dimensional ($d <= 6$) megapixel images or 3D LiDAR point clouds where $3^d << n$, *Mean Shift++ (\#55)* provides orders-of-magnitude linear speedup over quadratic mode-seekers (Proposition 2.2). For moderate $n$ ($<= 2 times 10^4$) where exact continuous KDE modes are needed, use *DENCLUE 2.0 (\#43)* (with noise threshold $xi$ and arbitrary-shape path merging) or *Mean Shift (\#34)* (with compact Epanechnikov kernel for finite termination); *never* deploy fixed-step *DENCLUE 1.0 (\#35)* without step decay (Theorem 4.2). Use *Blurring Mean Shift (\#54)* strictly with a validated early-stopping plateau criterion (Theorem 5.2).
5. *Regime 5 — Non-Euclidean Dissimilarity Spaces & Structured Domains:* When data live in arbitrary metric or non-metric dissimilarity spaces (sequences, graphs, shapes) lacking vector means, use *Quickshift (\#41)* for fast non-iterative hierarchical mode-seeking ($O(n^2 d)$, Theorem 6.3), *SNN-DPC (\#52)* or *DPC (\#51)* for exemplar peak detection, or *HDBSCAN (\#32)* / *SNN (\#36)* for arbitrary-shape connectivity. Avoid *Medoidshift (\#42)* unless $n <= 2000$ and a tuning-free MM medoid fixed point is specifically required, due to its $O(n^3)$ wall and locking-ball fragmentation (Theorem 6.2). For extended spatial objects (polygons, trajectories, weighted sensors), invoke *GDBSCAN (\#39)*; for space--time event streams, invoke *ST-DBSCAN (\#40)*—keeping $Delta epsilon = infinity$ (or pre-filtering attributes symmetrically) whenever strict permutation invariance is mandatory (Proposition 1.1).
6. *Regime 6 — Joint Multi-Density Clustering and Continuous Anomaly Detection:* Whenever local outliers near dense clusters have smaller absolute distances than interior points of sparse clusters, standard binary-noise DBSCAN and reachability-threshold OPTICS provably fail (Theorem 10.1). Deploy *OPTICS-OF (\#53)* with $epsilon >= max_p d_k (p)$ to obtain both the multi-scale cluster hierarchy and a scale-invariant continuous outlier factor $"OF"_(m, epsilon) (p)$ (Theorem 10.2).
7. *Regime 7 — Explicit Elimination Rules (When NOT to Deploy):* Never deploy *DBCLASD (\#38)* on Gaussian or smoothly peaked distributions (Proposition 1.2 proves asymptotic $chi^2$ self-fragmentation); never deploy *K-DBSCAN (\#44)* without collar merging on non-convex intertwined manifolds (Theorem 12.2); never deploy *DPC (\#51)* unmodified on datasets whose cluster density contrast exceeds $D_(12) / delta_s$ (Theorem 7.2); and never deploy *Mean Shift++ (\#55)* or *GriDBSCAN (\#48)* once $d > log_3 n$ or $2^d > n / P$.

#subsec("B", "Synthesis of Open Theoretical Frontiers Across Category II")

Unifying Algorithms 31--55 also crystallizes four deep mathematical problems that remain open across contemporary unsupervised learning research:
1. *Minimax Multi-Scale Level-Set Recovery on Low-Intrinsic-Dimension Manifolds:* While DeBaCl @chaudhuri2014consistent and #box[DBSCAN++] @jang2019dbscanpp attain minimax rates in ambient dimension $d$, proving sharp finite-sample tree-recovery rates for HDBSCAN and SNN-DPC that depend exclusively on intrinsic manifold dimension $d_0 << d$—without requiring prior knowledge of $d_0$—remains largely open.
2. *Provable Stopping-Time Criteria for Blurring Dynamics:* Theorem 5 proved that Gaussian Blurring Mean Shift achieves cubic local collapse yet inevitable global consensus; constructing a universally consistent, data-driven stopping rule $T^* (cal(X), h)$ that halts strictly inside the metastable multi-cluster plateau (preserving the exact modal count of frozen Mean Shift) lacks a general finite-sample guarantee.
3. *Subquadratic High-Dimensional Mutual-Reachability & SNN Graphs:* HDBSCAN, RNN-DBSCAN, SNN, and SNN-DPC all deliver strong multi-density robustness yet inherit $Theta(n^2 d)$ worst-case graph construction once $d >= 25$. Bridging approximate nearest-neighbor (ANN) graphs with provable preservation of mutual-reachability MST connectivity and SNN peak ranks—analogous to Theorem 12's $delta$-cover guarantee for DBSCAN++—is a high-impact algorithmic frontier.
4. *Cross-Paradigm Density--Hierarchical--Spectral Unification:* The mutual-reachability MST of HDBSCAN (Alg. 32) directly links Category II to Single Linkage in Category III (Alg. 56, @carlsson2010characterization), while the row-stochastic Gaussian affinity matrix $bold(P)^((t)) = (bold(D)^((t)))^(-1) bold(W)^((t))$ driving Blurring Mean Shift (Alg. 54) is identically the random-walk transition operator $I - L_"rw"$ governing Spectral Clustering and Diffusion Maps in Category V (Algs. 103--105, @vonluxburg2007tutorial). These cross-category equivalences are established rigorously in Comparison Chapters II, V, and VII of Part IX, alongside empirical multi-geometry validation across all density families in Comparison Chapter XIII.
