#import "../template/ieee.typ": *

#alg_entry(
  id: 46,
  name: "DeBaCl (Density-Based Clustering via Level Set Trees)",
  category: "Cat. II: Density & Mode-Seeking",
  year: "2013",
  paradigm: "Nonparametric Level-Set Cluster Tree Estimation",
  complexity: $O(n k log n + n k alpha(n))$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Classical flat density algorithms such as DBSCAN (Alg. 31) @ester1996dbscan partition data at a single density cutoff $lambda$, discarding multi-scale modal hierarchy and forcing a rigid trade-off between separating close modes and merging sparse clusters. Introduced by Kent, Rinaldo, and Verstynen @kent2013debacl, *DeBaCl* operationalizes Hartigan's @hartigan1981consistency nonparametric *density cluster tree* by estimating the full nested hierarchy of connected components across all upper level sets of an unknown probability density $f: bb(R)^d -> bb(R)_+$.

Let $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$ be an i.i.d. sample drawn from a continuous density $f$ with respect to Lebesgue measure on $bb(R)^d$. For any density threshold $lambda >= 0$, define the *population upper level set*
$ L_f (lambda) := { x in bb(R)^d : f(x) >= lambda }. $
Denote by $pi_0 (A)$ the set of path-connected components of a topological subspace $A subset.eq bb(R)^d$. Hartigan's *population level-set cluster tree* is the laminar collection
$ cal(T)_f := union.big_(lambda >= 0) pi_0 (L_f (lambda)), $
partially ordered by set inclusion: as $lambda$ increases from $0$ toward $norm(f)_infinity$, connected components either shrink continuously, split into child branches when $lambda$ crosses a saddle separating two local modes, or vanish once $lambda$ exceeds a mode's peak density.

To estimate $cal(T)_f$ from finite samples without parametric shape assumptions, DeBaCl decouples two statistical primitives: (i) a pointwise nonparametric *density estimator* $hat(f): cal(X) -> bb(R)_+$, canonically the $k$-nearest-neighbor ($k$-NN) estimator
$ hat(f)_k (x_i) = k / (n v_d r_k (x_i)^d), quad v_d = pi^(d / 2) / Gamma(d / 2 + 1), $
where $r_k (x_i)$ is the Euclidean distance from $x_i$ to its $k$-th nearest neighbor in $cal(X) backslash {x_i}$ (a kernel density estimator $hat(f)_h (x_i)$ may also be substituted); and (ii) a sparse *geometric similarity graph* $G = (V, E)$ on $V = {1, dots, n}$, canonically the symmetric $k$-NN graph ($E = {{i, j} : j in N_k (i) or i in N_k (j)}$) or the $epsilon$-neighborhood graph $G_epsilon = (V, {{i, j} : norm(x_i - x_j)_2 <= epsilon})$.

#definition(title: "Empirical Level-Set Tree & Mass-Pruned Tree")[
  For a finite graph $G = (V, E)$ and vertex density scores $hat(f): V -> bb(R)_+$, let $V_lambda := {i in V : hat(f)(x_i) >= lambda}$ and let $G_lambda := G[V_lambda]$ denote the subgraph induced by $V_lambda$. The *empirical level-set tree* is $hat(cal(T))(G, hat(f)) := union.big_(lambda >= 0) pi_0 (G_lambda)$, indexed equivalently by density level $lambda$ or upper-mass level $alpha = 1 - |V_lambda| / n in [0, 1]$. Given a *pruning threshold* $gamma in {1, dots, n}$, the *mass-pruned tree* $cal(P)_gamma (hat(cal(T)))$ is obtained by iteratively contracting every leaf branch $B$ with $|B| < gamma$ into its parent branch until every surviving leaf satisfies $|B| >= gamma$.
]

DeBaCl constructs $hat(cal(T))$ via a single descending-density sweep backed by Disjoint-Set Union (DSU), prunes finite-sample ripples below mass $gamma$, extracts high-density cluster cores (canonically the leaves of $cal(P)_gamma (hat(cal(T)))$, or a cut at mass level $alpha$ / first $K$ branches), and optionally assigns lower-density background points via nearest-neighbor or $k$-NN label propagation (Algorithm 1).

#algbox(
  title: "DeBaCl Level-Set Cluster Tree Construction & Extraction",
  inputs: [Dataset $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$, neighborhood integer $k >= 1$, pruning mass $gamma >= 1$, extraction mode (`leaf`, `upper`, or `first-k`), background assignment rule.],
  outputs: [Pruned level-set tree $cal(P)_gamma (hat(cal(T)))$ and cluster partition $cal(C)$.],
)[
  1. *Density & Similarity Graph:* Compute $k$-NN sets $N_k (i)$, radii $r_k (x_i)$, densities $hat(f)_k (x_i) = k / (n v_d r_k (x_i)^d)$, and symmetric graph $G = (V, E)$ with $|E| <= n k$.
  2. *Descending Density Sort:* Permute indices as $pi(1), dots, pi(n)$ such that $hat(f)_k (x_(pi(1))) >= hat(f)_k (x_(pi(2))) >= dots >= hat(f)_k (x_(pi(n)))$.
  3. *Union-Find Superlevel Sweep:* Initialize empty active set $A <- emptyset$, empty DSU forest, and empty tree $hat(cal(T))$.
  4. *For* $t = 1, dots, n$ *do*:
  5. #h(1em) Let $u = pi(t)$, $lambda_t = hat(f)_k (x_u)$, and $cal(M)_t = {"find"(v) : {u, v} in E, space v in A}$.
  6. #h(1em) Insert $u$ into $A$ as a singleton DSU set.
  7. #h(1em) *If* $cal(M)_t = emptyset$: birth a new leaf node $b_u$ in $hat(cal(T))$ at birth density $lambda_t$.
  8. #h(1em) *Else if* $|cal(M)_t| = 1$: attach $u$ to the unique active branch $b in cal(M)_t$ via $"union"(u, b)$.
  9. #h(1em) *Else* ($|cal(M)_t| >= 2$): close all branches $b in cal(M)_t$ at merge level $lambda_t$, create a new parent branch $b_"new"$ with children $cal(M)_t$, and union $u$ with all $b in cal(M)_t$.
  10. *Mass Pruning:* Iteratively collapse any leaf branch $B$ with $|B| < gamma$ into its parent (merging single-child internal chains) to obtain $cal(P)_gamma (hat(cal(T)))$.
  11. *Cluster Extraction:* Extract upper-level core sets $C_1^"core", dots, C_m^"core"$ from $cal(P)_gamma (hat(cal(T)))$ (e.g., all surviving leaves); assign each background point $x in cal(X) backslash union.big_j C_j^"core"$ to its nearest core label (or retain as background noise).
  12. *Return* $(cal(P)_gamma (hat(cal(T))), cal(C))$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish four foundational mathematical properties of DeBaCl: (A) exact laminar ultrametricity and algorithmic equivalence of the Union-Find sweep to a vertex-min maximum spanning forest; (B) classical Hartigan consistency and finite-sample cluster-tree recovery rates; (C) exact confluence, cardinality upper bounds, and resolution thresholds of $gamma$-mass pruning; and (D) full similarity/rank-order invariance alongside an exact finite-sample outlier breakdown bound for surviving upper-level cores.

#subsec("A", "Laminar Ultrametric Structure & Exactness of DSU Sweep")

#lemma(title: "Laminar Hierarchy & Superlevel Bottleneck Ultrametric", origin: "Original Proof")[
  For any finite graph $G = (V, E)$ and arbitrary vertex weights $hat(f): V -> bb(R)_+$, the empirical level-set collection $hat(cal(T))(G, hat(f)) = union.big_(lambda >= 0) pi_0 (G_lambda)$ is strictly *laminar*: for all $A, B in hat(cal(T))$,
  $ A inter B != emptyset ==> A subset.eq B space or space B subset.eq A. $
  Moreover, defining the *bottleneck-density similarity* between $u, v in V$ (with $s_G (u, v) = 0$ if disconnected in $G$) by
  $ s_G (u, v) := max_(P : u ~> v) min_(w in V(P)) hat(f)(w), $
  $s_G$ satisfies the strong ultrametric inequality
  $ s_G (u, v) >= min{s_G (u, z), space s_G (z, v)} quad forall u, v, z in V, $
  and $u, v$ belong to the same connected component of $G_lambda$ if and only if $s_G (u, v) >= lambda$.
]

#proof[
  Let $A in pi_0 (G_(lambda_1))$ and $B in pi_0 (G_(lambda_2))$ with $A inter B != emptyset$. Without loss of generality assume $lambda_1 <= lambda_2$. Then $V_(lambda_2) subset.eq V_(lambda_1)$, so $G_(lambda_2)$ is an induced subgraph of $G_(lambda_1)$. Since $B$ is connected in $G_(lambda_2)$, it remains connected in $G_(lambda_1)$. Because $A$ is a *maximal* connected subgraph of $G_(lambda_1)$ and shares at least one vertex with $B$, maximality forces $B subset.eq A$, proving laminarity.

  Next, let $u, v, z in V$. If $min{s_G (u, z), s_G (z, v)} = 0$, the ultrametric inequality holds trivially since $hat(f) >= 0$. Otherwise, let $P_1$ attain $s_G (u, z)$ from $u$ to $z$ and let $P_2$ attain $s_G (z, v)$ from $z$ to $v$. Concatenating $P_1$ and $P_2$ yields a walk from $u$ to $v$, which contains a simple path $P^*$ with $V(P^*) subset.eq V(P_1) union V(P_2)$. Therefore,
  $ s_G (u, v) &>= min_(w in V(P^*)) hat(f)(w) >= min_(w in V(P_1) union V(P_2)) hat(f)(w) \
               &= min{s_G (u, z), space s_G (z, v)}. $
  Finally, $u, v$ lie in the same component of $G_lambda$ iff there exists a path $P: u ~> v$ in $G$ whose vertices all satisfy $hat(f)(w) >= lambda$, which holds iff $max_(P: u ~> v) min_(w in V(P)) hat(f)(w) >= lambda$, i.e., $s_G (u, v) >= lambda$.
]

#theorem(title: "Exactness & Near-Linear Complexity of the DSU Sweep", origin: "Original Proof")[
  Define the symmetric edge-weighted graph $(G, w)$ with $w(u, v) := min{hat(f)(u), hat(f)(v)}$ for every ${u, v} in E$. Then:
  1. For every $t in {1, dots, n}$, the connected components maintained by Algorithm 1 (Steps 4–9) after processing all vertices of density $>= lambda_t$ coincide *identically* with $pi_0 (G_(lambda_t))$.
  2. The merge hierarchy $hat(cal(T))$ constructed by Algorithm 1 is isomorphic to the single-linkage hierarchy of any *maximum spanning forest* $F^*$ of $(G, w)$.
  3. Given $G = (V, E)$, Steps 2–9 execute in $O(n log n + (n + |E|) alpha(n))$ time and $O(n + |E|)$ memory, where $alpha(n)$ is the inverse Ackermann function.
]

#proof[
  (1) Proceed by induction on the sweep step $t$, assuming distinct densities (ties are resolved identically by batching equal-density vertices before recording merges). At $t = 1$, $A_1 = {pi(1)}$ is a singleton component matching $pi_0 (G_(lambda_1))$. Assume after step $t - 1$ the DSU partition equals $pi_0 (G[A_(t-1)])$. Inserting $u = pi(t)$ activates precisely those edges ${u, v} in E$ with $v in A_(t-1)$ (equivalently $hat(f)(v) >= hat(f)(u)$). Every other edge in $G[A_t]$ already lay inside $G[A_(t-1)]$. Merging $u$ with the DSU roots in $cal(M)_t$ therefore forms the exact transitive closure of $G[A_t]$, preserving the invariant $pi_0 (G[A_t])$.

  (2) For any edge $e = {u, v} in E$ with $hat(f)(v) >= hat(f)(u)$, $w(e) = hat(f)(u)$, so $e$ is first inspected precisely when $u$ enters $A$ in descending order of $w(e)$. Algorithm 1 unites two distinct DSU components via $e$ if and only if $e$ connects two components not already joined by heavier edges—which is literally *Kruskal's algorithm* for a maximum spanning forest $F^*$ of $(G, w)$ @gower1969minimum. By the cycle property of spanning forests, for every $u, v in V$, the unique path in $F^*$ maximizes the minimum edge weight $min_(e in P) w(e) = min_(w in V(P)) hat(f)(w) = s_G (u, v)$, yielding an isomorphic merge tree.

  (3) Sorting $n$ scalar densities takes $O(n log n)$ time. Each edge ${u, v} in E$ is examined at most twice (once from each endpoint), invoking at most $2 |E|$ `find` operations and at most $n - 1$ `union` operations. With union-by-rank and path compression, Tarjan's bound yields $O((n + |E|) alpha(n))$ total DSU time, and maintaining the tree nodes (at most $2 n - 1$ nodes) requires $O(n)$ space alongside $O(n + |E|)$ for $G$.
]

#subsec("B", "Statistical Consistency of Level-Set Cluster Tree Recovery")

Hartigan @hartigan1981consistency proved that classical Euclidean single-linkage clustering (Alg. 56) is *not* consistent for recovering $cal(T)_f$ in dimensions $d >= 2$, because sparse outlier chains create premature bridges across low-density valleys. By thresholding vertex density before connecting edges, DeBaCl restores full statistical consistency.

#theorem(title: "Consistent Level-Set Tree Recovery", origin: "Classical — Cited")[
  Let $X_1, dots, X_n$ be i.i.d. draws from a continuous density $f$ on $bb(R)^d$, and consider a $k$-NN level-set cluster tree estimator (with connectivity radius scaled proportionally to $r_k (x)$ or $k$-NN graph) @chaudhuri2014consistent @kent2013debacl. Suppose two disjoint connected components $A, A' in pi_0 (L_f (lambda))$ are separated by a *$(sigma, epsilon)$-valley*—that is, every continuous path from $A$ to $A'$ traverses a region $S$ of width $sigma > 0$ on which $f(x) <= (1 - epsilon) lambda$. Then there exist universal constants $c_1, c_2 > 0$ such that whenever
  $ k >= c_1 (d log n + log(1 / delta)) / epsilon^2 quad "and" quad k / n <= c_2 lambda v_d sigma^d, $
  with probability at least $1 - delta$ the empirical level-set tree $hat(cal(T))$ contains disjoint branches $hat(A), hat(A')$ satisfying $A inter cal(X) subset.eq hat(A)$ and $A' inter cal(X) subset.eq hat(A')$ simultaneously for all such separated pairs $(A, A')$.
]

#proof(label: "Proof Synthesis (Cited — Chaudhuri et al. 2014)")[
  Because the family of closed Euclidean balls in $bb(R)^d$ has Vapnik–Chervonenkis (VC) dimension $d + 2$, uniform relative deviation bounds @chaudhuri2014consistent ensure that with probability $>= 1 - delta$, every ball $B(x, r)$ with true mass $P(B(x, r)) >= k / n$ obeys $|P_n (B(x, r)) - P(B(x, r))| <= (epsilon / 4) P(B(x, r))$ once $k >= c_1 (d log n + log(1 / delta)) / epsilon^2$. On the high-density components $A$ and $A'$, where $f >= lambda$, the condition $k / n <= c_2 lambda v_d sigma^d$ forces $r_k (x) < sigma / 3$, guaranteeing both that $hat(f)_k (x) >= (1 - epsilon / 2) lambda$ throughout $A union A'$ and that overlapping $r_k$-balls cover $A$ and $A'$ into connected subgraphs. Conversely, inside the separating valley $S$ where $f <= (1 - epsilon) lambda$, uniform concentration implies every sample point $z in S$ has $hat(f)_k (z) < (1 - epsilon / 2) lambda$, while the valley thickness $sigma > 2 max_(x in A union A') r_k (x)$ prevents direct edges jumping across $S$. Hence $A inter cal(X)$ and $A' inter cal(X)$ enter $G_(hat(lambda))$ strictly before any bridging vertex in $S$, producing disjoint branches $hat(A), hat(A')$ in $hat(cal(T))$ @chaudhuri2014consistent @hartigan1981consistency.
]

#subsec("C", "Mass-Pruning Confluence, Leaf Bounds & Resolution Limit")

Finite-sample noise induces numerous small spurious leaves near flat density plateaus. DeBaCl addresses this via bottom-up $gamma$-mass pruning $cal(P)_gamma (hat(cal(T)))$.

#proposition(title: "Confluence, Leaf Bound & Minority-Mode Erasure of Mass Pruning", origin: "Original Derivation")[
  Let $hat(cal(T))$ be any finite level-set tree on $n$ vertices and let $gamma in {1, dots, n}$. Under DeBaCl's $gamma$-mass pruning operator $cal(P)_gamma$:
  1. *Order-Invariance (Confluence):* The pruned tree $cal(P)_gamma (hat(cal(T)))$ is uniquely determined and independent of the order in which sub-threshold leaves ($|B| < gamma$) are contracted.
  2. *Exact Leaf Cardinality Bound:* The number of surviving leaves $L(gamma)$ is monotonically non-increasing in $gamma$ and satisfies $L(gamma) <= floor(n / gamma)$ (for $n >= gamma$).
  3. *Hard Resolution Floor:* Any true cluster mode whose upper-level component above its saddle merge level contains strictly fewer than $gamma$ sample points is deterministically erased from $cal(P)_gamma (hat(cal(T)))$.
]

#proof[
  (1) For every node $B in hat(cal(T))$, define its total subtree mass $m(B) := |B|$ (the total number of points contained in $B$ and its descendants). Observe that contracting a leaf $B$ with $m(B) < gamma$ into its parent $p(B)$ transfers $B$'s points to $p(B)$, leaving $m(p(B))$ *strictly unchanged*! Hence the predicate "node $B$ has $m(B) >= gamma$ and possesses at least two children $C_1, C_2$ with $m(C_1) >= gamma$ and $m(C_2) >= gamma$" is an invariant structural property of $hat(cal(T))$ itself, independent of contraction order. Recursively, a branch survives in $cal(P)_gamma (hat(cal(T)))$ if and only if its subtree mass is at least $gamma$ and it is either the root or a child of a surviving split node; all nodes with $m(B) < gamma$ are absorbed into their lowest ancestor with $m >= gamma$, proving strict confluence.

  (2) By Lemma 1, distinct surviving leaves $B_1, dots, B_(L(gamma))$ of $cal(P)_gamma (hat(cal(T)))$ are pairwise disjoint subsets of $V = {1, dots, n}$ satisfying $|B_j| >= gamma$ for each $j$. Summing cardinalities yields
  $ n >= abs(union.big_(j=1)^(L(gamma)) B_j) = sum_(j=1)^(L(gamma)) |B_j| >= L(gamma) dot gamma ==> L(gamma) <= floor(n / gamma). $
  Furthermore, if $gamma_1 <= gamma_2$, every surviving split node at $gamma_2$ has at least two children of mass $>= gamma_2 >= gamma_1$, so every leaf of $cal(P)_(gamma_2)$ contains at least one leaf of $cal(P)_(gamma_1)$, implying $L(gamma_2) <= L(gamma_1)$.

  (3) Consider a bimodal 1D distribution with a major mode of $n - m$ points in $[0, 1]$ and a well-separated minority mode $C_"min"$ of $m < gamma$ points in $[10, 11]$ (with $m > k$ and separation $9 >> r_k$). Although $C_"min"$ forms a genuine disconnected leaf branch $B_"min"$ in the unpruned tree $hat(cal(T))$, its subtree mass is $m(B_"min") = m < gamma$. By the invariant characterization in (1), $B_"min"$ is deterministically contracted into the root, erasing $C_"min"$ completely.
]

#subsec("D", "Similarity Invariance & Finite-Sample Outlier Breakdown")

#proposition(title: "Similarity/Rank Invariance & Core Breakdown Bound", origin: "Original Proof")[
  Let DeBaCl use the symmetric $k$-NN graph $G$ and $k$-NN density estimator $hat(f)_k$ on $cal(X)$ with $n >= k + 1$ and pruning threshold $gamma >= 2$. Then:
  1. *Similarity & Monotone-Rank Invariance:* Under any Euclidean similarity $T(x) = c Q x + b$ ($c > 0$, $Q in O(d)$, $b in bb(R)^d$), or any strictly increasing reparameterization $phi(hat(f)_k (x))$, the pruned tree topology $cal(P)_gamma (hat(cal(T)))$ and mass-indexed level cuts are *strictly invariant*.
  2. *Anisotropic Scale Sensitivity:* Diagonal feature scaling $S(x) = "diag"(s_1, dots, s_d) x$ ($s_r != s_q > 0$) can alter both $G$ and the leaf structure of $cal(P)_gamma (hat(cal(T)))$.
  3. *Adversarial Outlier Breakdown Bound:* Let $D := "diam"(cal(X)) > 0$. Appending any set $Z = {z_1, dots, z_m}$ of $m < gamma$ arbitrary (even mutually collocated) outliers satisfying $min_(z in Z, x in cal(X)) norm(z - x)_2 > D$ leaves the set of surviving leaves above $lambda_"min" (cal(X)) := min_(x in cal(X)) hat(f)_k (x)$ and their branching topology on $cal(X)$ *identically unchanged*. Thus the upper-level core breakdown point satisfies $epsilon^* >= gamma / (n + gamma)$.
]

#proof[
  (1) Since $norm(T(x_i) - T(x_j))_2 = c norm(x_i - x_j)_2$ for $c > 0$, all pairwise distance orderings are preserved identically, leaving $N_k (i)$ and $G$ unchanged while scaling $r_k (T(x_i)) = c r_k (x_i)$ and $hat(f)_k (T(x_i)) = c^(-d) hat(f)_k (x_i)$. Because $c^(-d) > 0$ and any strictly increasing $phi$ preserve the descending vertex ordering $pi$, Algorithm 1 executes the exact same sequence of DSU merges and mass prunings.

  (2) Take two dense horizontal rows in $bb(R)^2$ of $m >= gamma$ points each along $y = 0$ and $y = 1$ with horizontal spacing $2$ and $k = 2$. Under identity scaling, vertical neighbors (distance $1$) connect the rows into one component. Applying $S = "diag"(1, 10)$ stretches vertical separation to $10 > 2$, disconnecting the $k$-NN graph into two distinct leaf clusters.

  (3) For every original point $x in cal(X)$, there exist $n - 1 >= k$ points in $cal(X) backslash {x}$ at distance $<= D < "dist"(x, Z)$. Hence neither $N_k (x)$ nor $r_k (x)$ changes for any $x in cal(X)$, preserving $hat(f)_k (x)$ and every edge within $cal(X)$ identically. Now partition $Z$ into $Z_"high" := {z in Z : hat(f)_k (z) >= lambda_"min" (cal(X))}$ and $Z_"low" := Z backslash Z_"high"$. Suppose some $z in Z$ acquires an edge to $x in cal(X)$ in the symmetric $k$-NN graph. Since $z in.not N_k (x)$, we must have $x in N_k (z)$, forcing
  $ r_k (z) >= norm(z - x)_2 > D >= max_(x' in cal(X)) r_k (x') ==> hat(f)_k (z) < lambda_"min" (cal(X)), $
  so $z in Z_"low"$. Therefore, *no edge* connects $Z_"high"$ to $cal(X)$! Consequently, during the sweep down to $lambda_"min" (cal(X))$, every component intersecting $Z_"high"$ is disjoint from $cal(X)$ and has size $<= |Z_"high"| <= m < gamma$, so $cal(P)_gamma$ prunes it completely. Vertices in $Z_"low"$ enter strictly below $lambda_"min" (cal(X))$, after all splits and merges of $hat(cal(T))(cal(X))$ are finalized. Thus every surviving upper-level core on $cal(X)$ is strictly preserved for all $m < gamma$, establishing $epsilon^* >= gamma / (n + gamma)$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

DeBaCl's computational cost splits cleanly into two phases: (i) $k$-NN graph and density construction, requiring $O(n k log n)$ expected time via $k$-d/ball trees in low ambient dimension $d$ and $O(n^2 d)$ exact time in high dimensions where spatial trees degenerate @beyer1999nearest; and (ii) level-set tree extraction and pruning, which runs in $O(n log n + n k alpha(n))$ time by Theorem 1. Crucially, once $hat(cal(T))$ is built, re-pruning at arbitrary thresholds $gamma$ or slicing at new mass levels $alpha$ requires only $O(n)$ tree traversal without recomputing distances—enabling rapid interactive exploratory analysis @kent2013debacl. Memory usage is $O(n (d + k))$ for the dataset and sparse graph, avoiding dense $O(n^2)$ pairwise matrices.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for DeBaCl],
    [Time / Space Complexity], [$O(n k log n + n k alpha(n))$ low-$d$ ($O(n^2 d)$ high-$d$); $O(n(d+k))$ space; $O(n)$ re-pruning.],
    [Geometric Inductive Bias], [Arbitrary non-convex, multi-scale superlevel connected components; no spherical assumption.],
    [Determinism / Seeding], [Fully deterministic (up to exact $k$-NN density ties); no random initialization.],
    [Invariances], [Invariant under rigid isometries, uniform scaling, and monotone density transforms; scale-sensitive.],
    [Outlier Breakdown Point], [Upper-level core breakdown $epsilon^* >= gamma / (n + gamma)$ (Proposition 2); background points assignable or isolated.],
    [Recommended Companions], [HDBSCAN (Alg. 32) for excess-of-mass extraction; OPTICS (Alg. 33) for reachability plots; DBSCAN++ (Alg. 47) for $n > 10^6$.],
  ),
  caption: [Summary of theoretical guarantees and operational properties of DeBaCl.],
)

#sec("IV", "Engineering and Scientific Applications")

DeBaCl excels when domain scientists must inspect *multi-scale modal structure*—discovering both coarse macro-clusters and nested high-density sub-modes—without committing a priori to a fixed cluster count $k$ or global density threshold $lambda$.

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Diffusion MRI White-Matter Fiber-Tract Segmentation:* Originally motivated by neuroimaging tractography @kent2013debacl, axonal streamline bundles form curved, branching, non-convex manifolds with strongly varying local fiber densities. DeBaCl's level-set hierarchy separates major fasciculi at low density levels while resolving finer sub-bundles at higher $lambda$ without breaking curved trajectories across Voronoi hyperplanes.
2. *Single-Cell Transcriptomics (scRNA-seq) & Flow Cytometry Gating:* After UMAP/PCA dimensionality reduction, immune and progenitor cell populations form nested density hierarchies (major lineages splitting into rare activation states). Mass pruning $gamma$ controls minimum detectable subpopulation size (Proposition 1) while preserving non-convex differentiation trajectories.
3. *Astrophysical Halo & Cosmic-Web Substructure Discovery:* Dark-matter simulations and galaxy redshift surveys exhibit hierarchical filaments, halos, and sub-halos across orders of magnitude in density. Because DeBaCl is invariant under monotone density transformations (Proposition 2), log-density and $k$-NN density estimates yield identical sub-halo merge trees.
4. *Geospatial Event & Mobility Hotspot Hierarchies:* Urban GPS traces, seismic epicenters, and crime logs contain high-density downtown peaks nested inside broader metropolitan corridors. DeBaCl's $O(n)$ interactive re-slicing lets analysts inspect multi-resolution spatial footprints without re-running clustering.
5. *Exploratory Anomaly & Background Separation in Industrial Telemetry:* By separating high-density operating regimes (surviving leaves of $cal(P)_gamma$) from low-density transitional states and isolated spikes (pruned by Proposition 2(3)), DeBaCl isolates healthy multimodal machine states from transient fault excursions.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

DeBaCl degrades in four identifiable regimes: (a) *high ambient dimensionality* ($d > 15$), where $k$-NN volume $r_k (x)^d$ suffers distance concentration @beyer1999nearest and noise dimensions flatten density contrast; (b) *rare minority clusters* with fewer than $gamma$ points, deterministically erased by Proposition 1(3); (c) *flat unimodal plateaus or heavily overlapping Gaussians without a valley*, which merge into a single leaf; and (d) *full-assignment noise inflation*, where assigning low-density background points to nearest leaves drags down geometric compactness (Silhouette) when true outliers should be left unassigned.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Why DeBaCl Fits / Prerequisite], [Caution / Alternative],
    [Neuroimaging / Tracts], [Curved branching manifolds; multi-scale leaves], [Tune $gamma$ for minor tracts],
    [Single-Cell / Cytometry], [Nested cell subtypes after PCA/UMAP], [Use HDBSCAN (Alg. 32) for auto-cut],
    [Astrophysics / Halos], [Scale-free rank invariance over density], [Use KD-trees in $d <= 6$],
    [High-D Tabular ($d > 20$)], [Requires prior subspace/PCA projection], [Use SNN (Alg. 36) or SUBCLU (Alg. 137)],
    [Streaming / Massive $n$], [Exact $k$-NN bottleneck on $n > 10^6$], [Use DBSCAN++ (Alg. 47) or BIRCH (Alg. 64)],
  ),
  caption: [Application suitability and engineering trade-off matrix for DeBaCl.],
)

#sec("V", "Structural Relations and Empirical Behavior")

DeBaCl sits at the exact mathematical intersection of density-based clustering (DBSCAN, HDBSCAN, OPTICS) and graph-based single-linkage hierarchies (SLINK, Alg. 56) @sibson1973slink @campello2015hierarchical.

#proposition(title: "Exact Diagonal Coupling Equivalence to HDBSCAN & DBSCAN", origin: "Original Unified Proof")[
  Let $r_k (x_i)$ be the $k$-NN radius on $cal(X)$ and define HDBSCAN's mutual-reachability metric (Alg. 32) @campello2013hdbscan by
  $ d_"mreach" (x_i, x_j) := max{r_k (x_i), space r_k (x_j), space norm(x_i - x_j)_2}. $
  For any scale $r > 0$, set $lambda(r) := k / (n v_d r^d)$ and let $G_r$ be the Euclidean $r$-ball graph on $cal(X)$. Then:
  1. The connected components of the HDBSCAN single-linkage hierarchy at threshold $r$ coincide *identically* with DeBaCl's superlevel components $pi_0 (G_r [V_(lambda(r))])$ evaluated along the coupled diagonal $(epsilon, lambda) = (r, lambda(r))$.
  2. Evaluating a single fixed level $lambda_0 = lambda(epsilon)$ on $G_epsilon$ yields the exact core-point clusters of DBSCAN$(epsilon, k + 1)$ (Alg. 31) @ester1996dbscan.
  3. Setting $hat(f) equiv 1$ (constant density) and sweeping $epsilon$ reduces the superlevel tree to Euclidean Single Linkage (Alg. 56) @sibson1973slink.
]

#proof[
  For any pair $(i, j)$, $d_"mreach" (x_i, x_j) <= r$ holds if and only if three inequalities hold simultaneously: (a) $r_k (x_i) <= r$, (b) $r_k (x_j) <= r$, and (c) $norm(x_i - x_j)_2 <= r$. Since $r |-> k / (n v_d r^d)$ is strictly decreasing on $(0, infinity)$, conditions (a)–(b) are equivalent to $hat(f)_k (x_i) >= lambda(r)$ and $hat(f)_k (x_j) >= lambda(r)$, i.e., $i, j in V_(lambda(r))$, while (c) is equivalent to ${i, j} in E(G_r)$. Thus the threshold graph ${d_"mreach" <= r}$ equals $G_r [V_(lambda(r))]$ vertex-for-vertex and edge-for-edge, proving (1). Setting $r = epsilon$ fixes $V_(lambda(epsilon)) = {i : |B(x_i, epsilon) inter cal(X)| >= k + 1}$, which is identically DBSCAN's core set, proving (2); taking constant $hat(f) equiv 1$ removes the vertex density filter so $d_"mreach"$ collapses to $norm(x_i - x_j)_2$, yielding Single Linkage (3).
]

This coupling reveals the key conceptual distinction: HDBSCAN @campello2015hierarchical sweeps along the 1D diagonal $epsilon = r(lambda)$ and selects clusters maximizing excess-of-mass stability, whereas DeBaCl @kent2013debacl decouples graph connectivity ($k$-NN or fixed $epsilon$) from the density sweep $lambda$, allowing independent smoothing control and explicit mass pruning $gamma$. Compared to mode-seeking gradient ascent (Mean Shift Alg. 34 @comaniciu2002mean; DPC Alg. 51 @rodriguez2014dpc), DeBaCl requires no step-size tuning and directly exposes saddle merge densities.

*Empirical Behavior on ClusterOverview Benchmarks:* Across the 4-dataset ClusterOverview benchmark suite (measuring Silhouette @rousseeuw1987silhouettes and Adjusted Rand Index @hubert1985comparing), DeBaCl achieves strong overall structure recovery with fast mean execution time of *1.836 ms*. On *Iris*, DeBaCl attains *Sil = 0.586* and *ARI = 0.567*, cleanly separating the isolated setosa mode while merging the overlapping versicolor–virginica ridge (which lacks a deep finite-sample valley) into a unified superlevel branch. On *Multi-Shape*, DeBaCl reaches *Sil = 0.552* (ARI = N/A), tracing non-convex arcs effectively while background-to-leaf assignment absorbs peripheral points rather than discarding them. On *Customers* (*Sil = 0.346*) and especially *High-D Subspace+Noise* (*Sil = 0.217*, ARI = N/A), performance drops predictably in accordance with Section IV-B: in $d = 12$ with ambient noise dimensions, $k$-NN volume $r_k (x)^d$ concentrates and full-point leaf assignment pulls noisy background observations into the nearest density branch, confirming the necessity of subspace selection or explicit noise retention in high-dimensional noisy regimes.
