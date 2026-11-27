#import "../template/ieee.typ": *

#alg_entry(
  id: 32,
  name: "HDBSCAN (Hierarchical Density-Based Spatial Clustering of Applications with Noise)",
  category: "Cat. II: Density & Mode-Seeking",
  year: "2013 / 2015",
  paradigm: "Mutual-Reachability MST & Excess-of-Mass Optimization",
  complexity: $O(n^2) "/" O(n log n)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Fixed-scale DBSCAN (Alg. 31) extracts a single horizontal cut of the density level-set hierarchy at radius $epsilon$, failing whenever genuine clusters reside at distinct local densities. HDBSCAN @campello2013hdbscan @campello2015hierarchical eliminates the scale hyperparameter $epsilon$ entirely by constructing the full continuous hierarchy of DBSCAN\* level sets over all $epsilon in (0, infinity)$ and solving a discrete excess-of-mass optimization problem to select an optimal non-horizontal cut.

Let $(cal(X), d)$ be a finite metric space of size $n$, parameterized by a smoothing parameter $m_"pts" in bb(N)_(>= 2)$ and a minimum cluster size $m_"cl" in bb(N)_(>= 2)$ (often set equal: $m_"pts" = m_"cl"$). For each $x_i in cal(X)$, define its *core distance* $d_"core" (x_i)$ as the distance from $x_i$ to its $m_"pts"$-th nearest neighbor in $cal(X)$ (including $x_i$ itself). HDBSCAN transforms the ambient space via the *mutual reachability distance*:
$ d_"mreach" (x_i, x_j) = cases(0 & "if" i = j, max { d_"core" (x_i), space d_"core" (x_j), space d(x_i, x_j) } & "if" i != j.) $
This transformation leaves dense interiors unchanged while pushing sparse boundary points away from every neighbor to at least their local $m_"pts"$-NN scale. Next, HDBSCAN builds a Minimum Spanning Tree (MST) $cal(T)_"MST"$ of the complete weighted graph $G_"mreach" = (cal(X), E, d_"mreach")$ and converts $cal(T)_"MST"$ into a complete single-linkage dendrogram $cal(H)$ indexed by density level $lambda = 1 / epsilon$ (with $lambda = 0$ at the root).

Because most dendrogram splits merely shed isolated noise points rather than bifurcating into genuine modes, HDBSCAN condenses $cal(H)$ top-down using $m_"cl"$ into a *condensed cluster tree* $cal(T)_"cond"$. When a parent cluster $C$ splits at density $lambda = 1 / epsilon$ into children $C_L, C_R$: (i) if $|C_L| >= m_"cl"$ and $|C_R| >= m_"cl"$, both become child nodes of $C$ in $cal(T)_"cond"$ with birth density $lambda_"min" (C_L) = lambda_"min" (C_R) = lambda$; (ii) if only one child satisfies size $>= m_"cl"$, the smaller subset "falls out" of $C$ as noise at departure density $lambda_"max" (x, C) = lambda$ while $C$ continues intact; (iii) if both have size $< m_"cl"$, $C$ terminates as a leaf and all remaining points depart at $lambda$.

#definition(title: "Excess-of-Mass Stability and Flat Cut Program")[
  For each node $C_i in cal(V)^* = cal(V)(cal(T)_"cond") backslash {C_"root"}$, define its *relative excess-of-mass stability* @campello2013hdbscan @campello2015hierarchical by
  $ S(C_i) = sum_(x in C_i) ( lambda_"max" (x, C_i) - lambda_"min" (C_i) ) >= 0. $
  Let $cal(L)$ be the leaves of $cal(T)_"cond"$ and $P(ell) subset.eq cal(V)^*$ be the non-root nodes on the unique path from $C_"root"$ to leaf $ell$. HDBSCAN's *Excess of Mass (EOM)* flat clustering solves the integer program:
  $ max_(delta in {0, 1}^(|cal(V)^*|)) &sum_(i in cal(V)^*) delta_i S(C_i) \
    "s.t." &sum_(j in P(ell)) delta_j <= 1 quad forall ell in cal(L). $
]

Every point belonging to a selected node ($delta_i = 1$) receives label $i$; points that fell out above all selected clusters are labeled noise.

#algbox(
  title: "HDBSCAN\* Algorithm (Campello et al., 2013/2015)",
  inputs: [Metric dataset $(cal(X), d)$, smoothing $m_"pts" >= 2$, min cluster size $m_"cl" >= 2$.],
  outputs: [Flat labels $L(x) in {1, dots, k} union {"noise"}$, condensed tree $cal(T)_"cond"$, GLOSH scores.],
)[
  1. *Core distances:* For each $x_i in cal(X)$, compute $d_"core" (x_i)$ via $m_"pts"$-NN search.
  2. *Mutual reachability MST:* Construct $cal(T)_"MST"$ on $(cal(X), d_"mreach")$ via Prim's algorithm or Dual-Tree Borůvka.
  3. *Hierarchy:* Sort the $n - 1$ edges of $cal(T)_"MST"$ in ascending weight order and build the single-linkage hierarchy $cal(H)$ via Union-Find.
  4. *Condensation:* Traverse $cal(H)$ top-down with threshold $m_"cl"$ to form $cal(T)_"cond"$, recording $lambda_"min" (C_i)$ and $lambda_"max" (x, C_i)$.
  5. *Stability:* Compute $S(C_i) = sum_(x in C_i) (lambda_"max" (x, C_i) - lambda_"min" (C_i))$ for all $i in cal(V)^*$.
  6. *Bottom-up EOM DP:* Initialize $hat(S)(C_i) <- S(C_i)$ and $delta_i <- 1$ for leaves $C_i in cal(L)$. Processing internal nodes $i in cal(V)^*$ bottom-up:
     - Let $Sigma_i <- sum_(j in "ch"(i)) hat(S)(C_j)$.
     - *If* $S(C_i) >= Sigma_i$: set $hat(S)(C_i) <- S(C_i)$, $delta_i <- 1$, and $delta_u <- 0$ for all proper descendants $u$ of $i$.
     - *Else:* set $hat(S)(C_i) <- Sigma_i$ and $delta_i <- 0$.
  7. *Return* selected clusters ${C_i : delta_i = 1}$, noise labels, and GLOSH scores.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five formal properties of HDBSCAN: (1) metric validity of mutual reachability and subdominant ultrametric sufficiency of the MST; (2) exact global optimality and linear complexity of bottom-up EOM dynamic programming; (3) the layer-cake integral representation of stability and constructive resolution of DBSCAN's multi-density impossibility; (4) full Euclidean similarity invariance alongside an exact characterization of EOM's large-cluster bias; and (5) asymptotic cluster-tree consistency.

#subsec("A", "Mutual-Reachability Metric Axioms and MST Ultrametric Sufficiency")

#lemma(title: "Metric Space Preservation of Mutual Reachability", origin: "Original Proof")[
  If $(cal(X), d)$ is a metric space and $m_"pts" >= 1$, then $d_"mreach"$ defines a valid metric on $cal(X)$.
]

#proof[
  Non-negativity ($d_"mreach"(x, y) >= 0$) and symmetry ($d_"mreach"(x, y) = d_"mreach"(y, x)$) are immediate from symmetry of $d$ and the maximum operator. If $x != y$, then $d(x, y) > 0$, whence $d_"mreach"(x, y) >= d(x, y) > 0$; combined with $d_"mreach"(x, x) = 0$, identity of indiscernibles holds.

  For the triangle inequality $d_"mreach"(x, z) <= d_"mreach"(x, y) + d_"mreach"(y, z)$, assume $x, y, z$ are pairwise distinct (otherwise equality or triviality holds). Because $d_"mreach"$ is non-negative and dominates each of its arguments, all three terms inside $d_"mreach"(x, z)$ satisfy:
  $ d_"core" (x) &<= d_"mreach" (x, y) <= d_"mreach" (x, y) + d_"mreach" (y, z), \
    d_"core" (z) &<= d_"mreach" (y, z) <= d_"mreach" (x, y) + d_"mreach" (y, z), \
    d(x, z) &<= d(x, y) + d(y, z) <= d_"mreach" (x, y) + d_"mreach" (y, z), $
  where the third line invokes the triangle inequality of $d$. Taking the maximum of the three left-hand sides establishes $d_"mreach"(x, z) <= d_"mreach"(x, y) + d_"mreach"(y, z)$.
]

Crucially, Lemma 1 justifies metric-tree spatial pruning (Ball-trees and Dual-Tree Borůvka) directly in mutual-reachability space.

#proposition(title: "Subdominant Ultrametric Sufficiency of the Mutual-Reachability MST", origin: "Original Derivation")[
  Let $u(x, y) = min_(pi: x -> y) max_(e in pi) d_"mreach" (e)$ denote the minimax path distance in $G_"mreach"$. Then:
  1. $u$ is the unique *subdominant ultrametric* of $d_"mreach"$ (the pointwise largest ultrametric satisfying $u <= d_"mreach"$).
  2. For every $x, y in cal(X)$, $u(x, y)$ equals the maximum edge weight on the unique path $P_"MST" (x, y)$ in any MST $cal(T)_"MST"$ of $G_"mreach"$.
  Hence $cal(T)_"MST"$ losslessly encodes the entire continuum of DBSCAN\* hierarchies across all $epsilon > 0$.
]

#proof[
  (1) For any $x, y, z$, concatenating minimizing paths $pi_1: x -> y$ and $pi_2: y -> z$ produces an $x$-to-$z$ walk of maximum edge weight $max{u(x, y), u(y, z)}$, proving the ultrametric inequality $u(x, z) <= max{u(x, y), u(y, z)}$. For maximality, let $v <= d_"mreach"$ be any ultrametric and let $pi = (v_0, dots, v_r)$ be any path from $x$ to $y$. Induction on $r$ gives $v(x, y) <= max_(1 <= j <= r) v(v_(j-1), v_j) <= max_(1 <= j <= r) d_"mreach"(v_(j-1), v_j)$. Minimizing over $pi$ yields $v(x, y) <= u(x, y)$.

  (2) Let $w(e) := d_"mreach"(e)$. By the MST cycle property, every non-tree edge $e = (a, b)$ satisfies $w(e') <= w(e)$ for every edge $e'$ on the unique tree path $P_"MST"(a, b)$ (otherwise replacing $e'$ by $e$ strictly lowers total tree weight). Given any $x$-to-$y$ path $pi$ in $G_"mreach"$, replacing each non-tree edge $e in pi$ with $P_"MST"(e)$ creates an $x$-to-$y$ walk inside $cal(T)_"MST"$ whose maximum edge weight is at most $max_(e in pi) w(e)$. Deleting backtracking loops yields the unique simple tree path $P_"MST"(x, y)$, whence $max_(e in P_"MST"(x, y)) w(e) <= max_(e in pi) w(e)$. Since $P_"MST"(x, y)$ is itself a valid path, equality holds.
]

#subsec("B", "Exact Global Optimality of EOM Dynamic Programming")

#theorem(title: "Global Optimality and O(n / m_cl) Exactness of Bottom-Up EOM", origin: "Original Proof")[
  For any condensed cluster tree $cal(T)_"cond"$ with non-negative stabilities $S(C_i) >= 0$, Step 6 of Algorithm 1 computes an *exact global maximizer* of the EOM integer linear program in Definition 1 in $O(|cal(V)(cal(T)_"cond")|) = O(n / m_"cl")$ operations.
]

#proof[
  For each node $i in cal(V)^*$, let $T_i$ be the subtree rooted at $i$ and let $Phi^*(i)$ denote the maximum of $sum_(u in T_i) delta_u S(C_u)$ subject to the path-exclusivity constraint $sum_(j in P(ell) inter T_i) delta_j <= 1$ for all leaves $ell in T_i$. We prove $hat(S)(C_i) = Phi^*(i)$ by structural induction on the height $h(i)$ of $T_i$.

  *Base case ($h(i) = 0$, $C_i in cal(L)$):* Because $S(C_i) >= 0$, the unique optimal assignment on the singleton subtree $T_i$ is $delta_i = 1$, attaining $Phi^*(i) = S(C_i) = hat(S)(C_i)$.

  *Inductive step ($h(i) >= 1$):* Partition feasible assignments on $T_i$ into two exhaustive cases:
  - *Case 1 ($delta_i = 1$):* Because $i$ lies on every root-to-leaf path through $T_i$, feasibility forces $delta_u = 0$ for every proper descendant $u in T_i backslash {i}$, yielding objective exactly $S(C_i)$.
  - *Case 2 ($delta_i = 0$):* The child subtrees ${T_j : j in "ch"(i)}$ are pairwise vertex-disjoint and share no path-exclusivity constraints below $i$. Consequently, the feasible set factors as a Cartesian product across $j in "ch"(i)$, and the optimal value over Case 2 is
    $ sum_(j in "ch"(i)) Phi^*(j) = sum_(j in "ch"(i)) hat(S)(C_j) = Sigma_i $
    by the induction hypothesis.

  Taking the maximum of Cases 1 and 2 yields $Phi^*(i) = max{S(C_i), Sigma_i} = hat(S)(C_i)$, and the indicator update in Step 6 maintains a realizing feasible solution on $T_i$. Finally, because $C_"root"$ is excluded from $cal(V)^*$, the full program decouples across children of $C_"root"$, achieving global optimum $sum_(j in "ch"(C_"root")) hat(S)(C_j)$.

  For runtime, every leaf of $cal(T)_"cond"$ contains at least $m_"cl"$ disjoint points from $cal(X)$, so $|cal(L)| <= floor(n / m_"cl")$. Because condensation suppresses unary (non-splitting) edges, every internal node has degree $>= 2$, giving $|cal(V)(cal(T)_"cond")| <= 2 |cal(L)| - 1 <= 2 floor(n / m_"cl")$. A single post-order pass (with lazy top-down marker propagation instead of eager descendant clearing) executes in $O(n / m_"cl")$ time.
]

#subsec("C", "Excess-of-Mass Layer-Cake Integral and Multi-Density Resolution")

#theorem(title: "Layer-Cake Integral Representation of Cluster Stability", origin: "Classical — Cited")[
  For every node $C_i in cal(T)_"cond"$, let $C_i (lambda) = {x in C_i : lambda_"max" (x, C_i) >= lambda}$ denote the surviving subset of $C_i$ at density level $lambda >= lambda_"min" (C_i)$. Then HDBSCAN's stability admits the exact integral identity @campello2015hierarchical:
  $ S(C_i) = integral_(lambda_"min" (C_i))^infinity |C_i (lambda)| d lambda. $
  Under the volume-rescaled density coordinate $tilde(lambda)(epsilon) = m_"pts" / (n V_d epsilon^d)$, $n^(-1) S(C_i)$ is the empirical Monte Carlo estimator of Hartigan's population excess of mass $integral_(C_i) (f(x) - lambda_"min"(C_i)) d x$.
]

#proof[
  For each point $x in C_i$, since $lambda_"max"(x, C_i) >= lambda_"min"(C_i)$, the fundamental theorem of calculus writes its persistence length as
  $ lambda_"max" (x, C_i) - lambda_"min" (C_i) = integral_(lambda_"min" (C_i))^infinity bb(I)[x in C_i (lambda)] d lambda. $
  Summing over all $x in C_i$ and interchanging the finite sum with the non-negative integral (Tonelli's theorem) yields
  $ S(C_i) = integral_(lambda_"min" (C_i))^infinity sum_(x in C_i) bb(I)[x in C_i (lambda)] d lambda = integral_(lambda_"min" (C_i))^infinity |C_i (lambda)| d lambda. $
  Substituting $tilde(lambda) = hat(f)(x)$ replaces the default $1/epsilon$ scale (corresponding to $f^(1/d)$) with the $d$-dimensional density level, converging to the population excess of mass as $n -> infinity$ @campello2015hierarchical.
]

#corollary(title: "Multi-Density Recovery Across Heterogeneous Scales", origin: "Original Proof")[
  Consider three uniform 1D clusters $C_1, C_2, C_3$ of size $N >= max(m_"cl", 2 m_"pts")$: dense clusters $C_1, C_2$ (spacing $delta$) separated by gap $g$, and sparse cluster $C_3$ (spacing $Delta > g$) separated from $C_1 union C_2$ by $R >> Delta$. Whereas no fixed-$epsilon$ DBSCAN cut can recover ${C_1, C_2, C_3}$, HDBSCAN with EOM selection recovers ${C_1, C_2, C_3}$ identically whenever $g > 2 (m_"pts" - 1) delta$.
]

#proof[
  Inside $C_1$ and $C_2$, every point has at least $m_"pts" - 1$ neighbors within distance $(m_"pts" - 1) delta$, so $d_"core"(x) <= (m_"pts" - 1) delta < g$. Inside $C_3$, $d_"core"(x) <= (m_"pts" - 1) Delta << R$. Thus the two heaviest MST edges connect $C_(12) = C_1 union C_2$ to $C_3$ at weight $R$ and $C_1$ to $C_2$ at weight $g$. Because $|C_1| = |C_2| = |C_3| = N >= m_"cl"$, $cal(T)_"cond"$ has root children $C_(12)$ and $C_3$ born at $lambda_0 = 1/R$; $C_(12)$ splits into leaves $C_1, C_2$ at $lambda_1 = 1/g$. Their stabilities satisfy:
  $ S(C_(12)) &= 2 N (1 / g - 1 / R), \
    S(C_1) + S(C_2) &>= 2 N ( 1 / ((m_"pts" - 1) delta) - 1 / g ). $
  Whenever $g > 2 (m_"pts" - 1) delta$, we have $1 / ((m_"pts" - 1) delta) - 1 / g > 1 / g > 1 / g - 1 / R$, hence $S(C_1) + S(C_2) > S(C_(12))$. Because $S(C_3) > 0$, Theorem 1 selects the exact non-horizontal cut ${C_1, C_2, C_3}$.
]

#subsec("D", "Similarity Invariance, Root-Dominance Bias, and Consistency")

#proposition(title: "Similarity Invariance and EOM Large-Cluster Bias", origin: "Original Proof")[
  1. *Full Euclidean Similarity Invariance:* For any similarity map $T(x) = alpha Q x + b$ ($alpha > 0$, $Q in O(d)$, $b in bb(R)^d$), HDBSCAN's condensed tree topology and EOM flat partition are strictly invariant (modulo equal-weight MST ties).
  2. *EOM Plateau Merging Threshold:* Suppose a condensed parent node $C_0$ (born at $lambda_0$) splits at $lambda_1$ without prior point fallout into $k >= 2$ equal-mass leaf modes $C_1, dots, C_k$ dying at $lambda_2$. Then EOM merges all $k$ modes back into $C_0$ if and only if $lambda_2 < 2 lambda_1 - lambda_0$, independently of $k$.
]

#proof[
  (1) Because $norm(T(x) - T(y))_2 = alpha norm(x - y)_2$, every core distance and mutual reachability distance scales uniformly by $alpha$. Hence the MST edge ordering is identical, every density level $lambda = 1 / epsilon$ scales by $alpha^(-1)$, and every node stability scales uniformly: $S_T (C_i) = alpha^(-1) S(C_i)$. Multiplying the linear objective in Definition 1 by $alpha^(-1) > 0$ preserves its exact maximizers.

  (2) Under the hypothesis, $|C_0| = k N$ and $|C_j| = N$, giving $S(C_0) = k N (lambda_1 - lambda_0)$ and $sum_(j=1)^k S(C_j) = k N (lambda_2 - lambda_1)$. Canceling $k N > 0$ gives $S(C_0) > sum_(j=1)^k S(C_j) <=> lambda_1 - lambda_0 > lambda_2 - lambda_1 <=> lambda_2 < 2 lambda_1 - lambda_0$.
]

Under smooth densities on $bb(R)^d$, Chaudhuri et al. @chaudhuri2014consistent establish that pruned single-linkage trees built on $k$-NN density graphs with $m_"pts" = Theta(log n)$ (or $m_"pts" / log n -> infinity$, $m_"pts" / n -> 0$) consistently recover the Hartigan cluster tree @hartigan1981consistency in the minimax-optimal sense.

#sec("III", "Complexity, Geometry and Operational Profile")

HDBSCAN decomposes into four clean computational stages: (1) $m_"pts"$-NN core queries ($O(n log n)$ via KD/Ball-trees in low $d$, $O(n^2 d)$ in high $d$); (2) MST construction on the implicit metric graph $(cal(X), d_"mreach")$, taking $O(n^2)$ time and strictly $O(n)$ space via Prim's algorithm (or empirical $O(n log n)$ via Dual-Tree Borůvka leveraging Lemma 1's metric triangle inequality); (3) sorting $n - 1$ MST edges and Union-Find dendrogram assembly in $O(n log n)$; and (4) condensation and EOM dynamic programming in $O(n)$ time (Theorem 1). Critically, HDBSCAN never materializes the full $O(n^2)$ pairwise matrix when using Prim or Borůvka.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for HDBSCAN],
    [Time / Space Complexity], [$O(n log n)$ empirical (Dual-Tree Borůvka, low $d$); $O(n^2)$ worst-case Prim; $O(n)$ space.],
    [Geometric Inductive Bias], [Arbitrary non-convex manifolds spanning *variable* local densities; isolates sparse background.],
    [Determinism / Parameters], [Deterministic (up to equal-weight ties); only intuitive size parameters $(m_"cl", m_"pts")$, no $epsilon$.],
    [Invariances], [Strictly invariant under all Euclidean similarities $alpha Q x + b$ (Prop. 2(1)).],
    [Outlier Scoring], [Explicit noise rejection + continuous local-global anomaly score GLOSH @campello2015hierarchical.],
    [Recommended Upgrades], [`leaf` selection when EOM over-merges (Prop. 2(2)); UMAP reduction before high-$d$; OPTICS (Alg. 33) for plots.],
  ),
  caption: [Summary of theoretical guarantees and operational profile of HDBSCAN.],
)

#sec("IV", "Engineering and Scientific Applications")

By removing DBSCAN's brittle $epsilon$ parameter while retaining non-convex topology and noise rejection, HDBSCAN has become the premier density-clustering workhorse in modern data science.

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Single-Cell RNA-seq and CyTOF Population Discovery:* Post-UMAP cellular landscapes combine ultra-dense abundant lineages with sparse rare progenitor niches and technical doublets. Corollary 1 enables joint extraction across density scales while assigning doublets and debris to noise.
2. *Dense Neural Embedding & Topic Discovery (BERTopic):* Reduced sentence-transformer embeddings form variable-density semantic islands amid off-topic noise documents. Similarity invariance (Proposition 2) and noise filtering prevent forcing unrelated documents into coherent topics.
3. *Gaia Astrometry & Galactic Star-Cluster Mapping:* Open star clusters at varying heliocentric distances exhibit orders-of-magnitude variation in apparent phase-space density; HDBSCAN extracts both compact nearby cores and diffuse distant associations in one pass.
4. *Cybersecurity & Financial Anomaly Triage via GLOSH:* Because HDBSCAN computes the entire condensed hierarchy, every point receives a continuous Global-Local Outlier Score from Hierarchies (GLOSH @campello2015hierarchical) calibrated to its parent cluster's peak density.
5. *Geospatial Mobility Hubs & Ecological Patch Analysis:* Dense metropolitan transit nodes and sparse rural waypoints are resolved jointly, with stability $S(C_i)$ quantifying statistical persistence across spatial scales.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

HDBSCAN encounters limitations under: (a) *high ambient dimension* ($d > 15$), where $k$-NN core distances suffer concentration @beyer1999nearest and Dual-Tree Borůvka degrades to $O(n^2)$; (b) *EOM root dominance* on broad flat plateaus (Proposition 2(2)), mitigated by `leaf` extraction or setting `cluster_selection_epsilon`; (c) *micro-clusters smaller than $m_"cl"$*, which are systematically absorbed or discarded as noise; and (d) *single-linkage bridging* along high-density filaments connecting otherwise distinct populations.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Why HDBSCAN Fits / Prerequisite], [Caution / Alternative],
    [NLP / BERTopic], [Variable topic density; noise isolation], [Use UMAP first; `leaf` if over-merged],
    [Single-Cell Omics], [Rare + abundant cell niches simultaneously], [Tune $m_"cl"$ below rarest population],
    [Astrophysics / GIS], [Full similarity invariance; multi-scale], [Monitor chaining; compare DeBaCl (Alg. 46)],
    [Fraud / Outlier Triage], [Native GLOSH local-global scores], [Check small-ring size vs. $m_"cl"$],
    [Raw High-D Tabular], [Core distances concentrate in $d > 20$], [Use SNN (Alg. 36) or subspace methods],
  ),
  caption: [Application suitability and engineering matrix for HDBSCAN.],
)

#sec("V", "Structural Relations and Empirical Behavior")

HDBSCAN unifies fixed-radius density clustering, agglomerative single linkage, and hierarchical outlier detection inside a single algebraic framework.

#proposition(title: "Exact Reductions to SLINK, DBSCAN\*, and Similarity-Invariant GLOSH", origin: "Original Unified Proof")[
  1. Setting $m_"pts" = 1$ forces $d_"core" (x) = 0$ identically, so $d_"mreach" = d$ and HDBSCAN's uncondensed hierarchy $cal(H)$ coincides identically with classical Single-Linkage (SLINK, Alg. 56 @sibson1973slink @gower1969minimum).
  2. Taking a horizontal cut of $cal(H)$ at level $epsilon$ (instead of EOM selection) recovers DBSCAN\*($epsilon, m_"pts"$) (Alg. 31) identically.
  3. Define the *GLOSH* outlier score @campello2015hierarchical of $x in cal(X)$ by
     $ "GLOSH"(x) = 1 - lambda_"max" (x) / lambda_"max" (C_x), $
     where $lambda_"max"(x)$ is the density at which $x$ departs its lowest-level surviving cluster $C_x$ in $cal(T)_"cond"$ and $lambda_"max"(C_x) = sup_(y in C_x) lambda_"max"(y)$. Then $"GLOSH"(x) in [0, 1)$, is strictly invariant under all Euclidean similarities $T(x) = alpha Q x + b$, and equals $0$ iff $x$ attains the peak local density of $C_x$.
]

#proof[
  (1) For $m_"pts" = 1$, every point is its own 1st nearest neighbor at distance $0$, so $d_"mreach"(x, y) = max{0, 0, d(x, y)} = d(x, y)$. By Proposition 1, the MST of $(cal(X), d)$ generates the exact Single-Linkage dendrogram @gower1969minimum. (2) Follows directly from the mutual-reachability cut equivalence in Alg. 31. (3) By definition $0 < lambda_"max"(x) <= lambda_"max"(C_x)$, so $lambda_"max"(x) / lambda_"max"(C_x) in (0, 1]$ and $"GLOSH"(x) in [0, 1)$, vanishing iff $lambda_"max"(x) = lambda_"max"(C_x)$. Under any similarity $T(x) = alpha Q x + b$, both numerator and denominator scale by $alpha^(-1)$, canceling identically.
]

HDBSCAN further connects to OPTICS (Alg. 33, which traverses an asymmetric reachability ordering), DeBaCl (Alg. 46 @kent2013debacl, level-set trees via $k$-NN graphs), Genie (Alg. 69 @gagolewski2016genie, Gini-regularized single linkage), and HCA-DBSCAN (Alg. 70).

*Empirical Behavior on the ClusterOverview 4-Dataset Suite:* HDBSCAN delivers standout performance on *Customers* ($"Sil" = 0.777$), outperforming both DBSCAN ($0.450$) and OPTICS ($0.328$) by over $+0.32$ Silhouette—vividly corroborating Corollary 1's multi-density recovery across heterogeneous customer spending regimes. On *Multi-Shape* ($"Sil" = 0.795$, $"ARI" = "N/A"$), it matches DBSCAN ($0.803$) on non-convex arcs without requiring any manual $epsilon$ tuning. On *Iris* ($"Sil" = 0.608$, $"ARI" = 0.542$), it slightly edges DBSCAN's ARI ($0.534$) while retaining high cluster separation. On *High-D Subspace+Noise* ($"Sil" = 0.350$, $"ARI" = "N/A"$), it shares DBSCAN's high-dimensional $k$-NN concentration floor ($0.352$). Notably, its mean runtime of $2.469$ ms is $5.8 times$ faster than DBSCAN ($14.321$ ms), reflecting the speed of single-pass MST construction and $O(n / m_"cl")$ EOM dynamic programming.
