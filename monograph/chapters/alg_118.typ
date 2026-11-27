#import "../template/ieee.typ": *

#alg_entry(
  id: 118,
  name: "HCS (Highly Connected Subgraphs)",
  category: "Cat. V: Graph, Spectral & Network",
  year: "2000",
  paradigm: "Recursive Edge-Connectivity Min-Cut Partitioning",
  complexity: $cal(O)(N_("cut") (m n + n^2 log n))$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Unlike modularity-based methods (Algs. 115--117) that optimize a global null-model contrast, the *Highly Connected Subgraphs (HCS)* algorithm introduced by Hartuv and Shamir @hartuv2000hcs defines clusters purely through an intrinsic combinatorial edge-connectivity threshold. Let $G = (V, E)$ be a finite simple undirected graph with $n = |V|$ vertices and $m = |E|$ edges (typically constructed from pairwise similarity scores via an unweighted threshold graph $A_(i j) = bold(1)[s(x_i, x_j) >= tau]$ or symmetric $k$NN graph). For any non-empty proper subset $S subset V$, let $E(S, bar(S))$ denote the cut set of edges with one endpoint in $S$ and one in $bar(S) = V without S$.

The *edge connectivity* of a graph $H = (V_H, E_H)$ ($|V_H| >= 2$) is the minimum cardinality of an edge cut disconnecting $H$:
$ kappa'(H) = min_(emptyset != S subset V_H) |E_H (S, V_H without S)|. $
By convention, $kappa'(H) = 0$ if $H$ is disconnected, and every vertex's degree bounds edge connectivity from above: $kappa'(H) <= delta(H) = min_(v in V_H) "deg"_H (v)$.

#definition(title: "Highly Connected Graph and HCS Clustering")[
  A graph $H = (V_H, E_H)$ with $n_H = |V_H|$ vertices is *highly connected* if its edge connectivity strictly exceeds half its vertex count:
  $ kappa'(H) > n_H / 2, \
    "equivalently" quad kappa'(H) >= floor(n_H / 2) + 1. $
  A *highly connected subgraph* of $G$ is an induced subgraph $G[C]$ ($C subset.eq V$) that is highly connected. By inspection, singletons ($n_H = 1$) are unclustered outliers and isolated edges ($K_2$, $n_H = 2, kappa'(K_2) = 1$) are not highly connected; the smallest non-trivial highly connected graph is the triangle $K_3$ ($kappa'(K_3) = 2 > 3 / 2$).
]

The core HCS paradigm recursively computes a global minimum cut $(S, bar(S))$ of the current subgraph $H$: if $|E_H (S, bar(S))| > |V_H| / 2$, then $H$ is certified highly connected and emitted as a cluster; otherwise, the cut edges are removed and HCS recurses independently on $H[S]$ and $H[bar(S)]$. Hartuv and Shamir @hartuv2000hcs also proposed three practical enhancements—*low-degree pruning*, *iterated HCS*, and *singleton adoption* (reattaching removed singletons to the cluster maximizing incident edges if at least $tau_("adopt")$ links exist)—to counter sensitivity to peripheral vertices.

#algbox(
  title: "Hartuv-Shamir HCS Recursive Min-Cut Algorithm",
  inputs: [Simple undirected graph $G = (V, E)$; optional singleton adoption threshold $tau_("adopt")$ (default disabled for pure HCS).],
  outputs: [Disjoint clusters $cal(C) = {C_1, dots, C_K}$ where each $G[C_r]$ is highly connected, plus singleton set $cal(S)$.],
)[
  1. *Procedure* $"HCS"(H = (V_H, E_H))$:
     - If $|V_H| <= 1$, add $V_H$ to singleton set $cal(S)$ and return.
     - If $H$ is disconnected, recurse on each connected component of $H$ and return.
     - Compute a global minimum edge cut $(S, bar(S))$ of $H$ (e.g., via Stoer-Wagner or Nagamochi-Ibaraki) with cut size $c = |E_H (S, bar(S))| = kappa'(H)$.
     - *Connectivity certificate:* If $c > |V_H| / 2$, emit $V_H$ as a valid cluster in $cal(C)$ and return.
     - *Recursive split:* Delete cut edges $E_H (S, bar(S))$ and invoke $"HCS"(H[S])$ and $"HCS"(H[bar(S)])$.
  2. *Optional post-processing (@hartuv2000hcs):* Iteratively adopt any singleton $v in cal(S)$ into cluster $C^* = arg max_(C in cal(C)) |N_G (v) inter C|$ whenever $|N_G (v) inter C^*| >= tau_("adopt")$.
  3. *Return:* $(cal(C), cal(S))$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

#subsec("A", "Diameter-2, Common-Neighbor, and Exact Degree Equivalence")

Every non-trivial cluster returned by pure HCS obeys tight deterministic guarantees on diameter, overlap, edge density, and Hamiltonicity, alongside a remarkable exact degree characterization.

#theorem(title: "Structural Cohesion Bounds and Exact Minimum-Degree Characterization", origin: "Classical — Cited")[
  Let $H = (V_H, E_H)$ be any simple graph with $n_H = |V_H| >= 2$ @hartuv2000hcs. Then:
  (i) *Exact degree equivalence:* $H$ is highly connected ($kappa'(H) > n_H / 2$) if and only if $delta(H) > n_H / 2$, in which case $kappa'(H) = delta(H) >= floor(n_H / 2) + 1$;
  (ii) *Diameter and triangle redundancy:* Whenever $n_H >= 3$ and $H$ is highly connected, $"diam"(H) <= 2$, and every pair of non-adjacent distinct vertices $u, v in V_H$ shares at least $3$ common neighbors (at least $4$ when $n_H$ is even);
  (iii) *Quadratic edge floor and density:* $|E_H| >= n_H (floor(n_H / 2) + 1) / 2 > n_H^2 / 4$, implying edge density $rho(H) = |E_H| / binom(n_H, 2) > n_H / (2(n_H - 1)) > 1 / 2$, and $H$ is Hamiltonian.
]

#proof[
  For (i), disconnecting a minimum-degree vertex yields $delta(H) >= kappa'(H)$, so $kappa'(H) > n_H / 2 => delta(H) > n_H / 2$. Conversely, suppose $delta(H) >= floor(n_H / 2) + 1$ (which requires $n_H >= 3$). Take any cut $(S, V_H without S)$ with $1 <= s = |S| <= floor(n_H / 2)$. If $s = 1$, $|E_H (S, V_H without S)| >= delta(H)$. If $2 <= s <= floor(n_H / 2)$, each $u in S$ has at most $s - 1$ internal neighbors and thus at least $delta(H) - s + 1$ crossing edges, yielding
  $ |E_H (S, V_H without S)| - delta(H) &>= s (delta(H) - s + 1) - delta(H) \
                                        &= (s - 1)(delta(H) - s) > 0, $
  because $s >= 2$ and $delta(H) >= floor(n_H / 2) + 1 > s$. Hence every minimum cut isolates a minimum-degree vertex, giving $kappa'(H) = delta(H) > n_H / 2$. (Note that while *testing* high connectivity reduces to an $cal(O)(n_H)$ degree check, *splitting* a non-HC subgraph still requires a global min-cut when $delta(H) <= n_H / 2$.)

  For (ii), if $n_H in {3, 4}$, (i) forces $delta(H) = n_H - 1$, so $H tilde.op K_(n_H)$ and $"diam"(H) = 1$. For $n_H >= 5$, any distinct non-adjacent $u, v in V_H$ have neighborhoods $N_H (u), N_H (v) subset.eq V_H without {u, v}$ (cardinality $n_H - 2$). By inclusion-exclusion and (i),
  $ |N_H (u) inter N_H (v)| &= |N_H (u)| + |N_H (v)| - |N_H (u) union N_H (v)| \
                          &>= 2 (floor(n_H / 2) + 1) - (n_H - 2), $
  which equals $4$ for even $n_H$ and $3$ for odd $n_H$, proving $"diam"(H) <= 2$.

  For (iii), the handshaking lemma and (i) yield $|E_H| = (1 / 2) sum_(v in V_H) "deg"_H (v) >= n_H (floor(n_H / 2) + 1) / 2 > n_H^2 / 4$, dividing by $binom(n_H, 2)$ gives $rho(H) > n_H / (2(n_H - 1)) > 1 / 2$, and Dirac's theorem ($delta(H) > n_H / 2$ for $n_H >= 3$) ensures Hamiltonicity.
]

#subsec("B", "Cut-Budget Bound, Tight Edge Expansion, and Spectral Floor")

We next prove exact bounds on both the global recursion tree (total calls and total removed edges) and the tight isoperimetric / spectral expansion of every HCS cluster.

#theorem(title: "Exact Cut-Budget Bound and Tight Isoperimetric / Laplacian Floor", origin: "Original Proof")[
  (i) *Recursion and edge-cut budget:* On any $n$-vertex graph $G$ ($n >= 2$), pure HCS terminates after at most $n - 1$ splitting min-cuts and at most $2 n - 1$ total min-cut calls, removing at most $floor(n^2 / 4)$ edges in aggregate across the entire recursion tree;
  (ii) *Tight edge expansion:* Every non-trivial HCS cluster $H$ ($n_H >= 3$) has isoperimetric edge-expansion number
  $ h(H) = min_(1 <= |S| <= floor(n_H / 2)) (|E_H (S, V_H without S)|) / |S| >= 2, $
  and this lower bound $2$ is attained with equality on an infinite family of highly connected graphs;
  (iii) *Algebraic connectivity floor:* The second-smallest eigenvalue $lambda_2 (L_H)$ of the unnormalized Laplacian $L_H = D_H - A_H$ satisfies
  $ lambda_2 (L_H) >= (h(H)^2) / (2 Delta(H)) >= 2 / (n_H - 1), $
  where $Delta(H) = max_v "deg"_H (v) <= n_H - 1$.
]

#proof[
  For (i), let $cal(T)$ be the binary recursion tree whose leaves are the $L <= n$ final emitted clusters and singletons. Every non-trivial split creates two strictly smaller non-empty subgraphs, so $cal(T)$ is a full binary tree with $L$ leaves, exactly $L - 1 <= n - 1$ internal split nodes, and $2 L - 1 <= 2 n - 1$ total nodes. At any internal node $H$ of size $s = |V_H| >= 2$, splitting requires $|E_H (S, bar(S))| <= floor(s / 2)$. Let $F(s)$ denote the maximum total cut edges removed over any valid binary recursion tree on $s$ vertices: $F(1) = 0$ and
  $ F(s) = max_(1 <= a <= s - 1) (floor(s / 2) + F(a) + F(s - a)). $
  We prove $F(s) <= floor(s^2 / 4)$ by strong induction on $s$. Base $s = 1$: $0 <= 0$. Assuming $F(a) <= floor(a^2 / 4)$ for $a < s$, observe that $floor(a^2 / 4) + floor((s - a)^2 / 4)$ is convex in $a$ and maximized at the boundary $a = 1, s - a = s - 1$, where
  $ floor(s / 2) + 0 + floor((s - 1)^2 / 4) = floor(s^2 / 4) $
  (verifiable directly for $s = 2 r$: $r + r(r - 1) = r^2$; and for $s = 2 r + 1$: $r + r^2 = floor((2 r + 1)^2 / 4)$). Thus total cut edges never exceed $floor(n^2 / 4)$.

  For (ii), take any subset $S subset V_H$ with $1 <= s = |S| <= floor(n_H / 2)$. By Theorem II.A(i), every vertex $u in S$ has at least $floor(n_H / 2) + 1$ incident edges in $H$, of which at most $s - 1$ can lie inside $S$. Hence each $u in S$ contributes at least $floor(n_H / 2) + 1 - (s - 1) = floor(n_H / 2) + 2 - s$ edges to $E_H (S, V_H without S)$. Summing over $u in S$ and dividing by $s$ yields
  $ (|E_H (S, V_H without S)|) / s >= floor(n_H / 2) + 2 - s >= 2, $
  since $s <= floor(n_H / 2)$. To show tightness, for any $r >= 3$ construct $H_(2 r)$ ($n_H = 2 r$) from two cliques $C_1, C_2 tilde.op K_r$ linked by a 2-regular bipartite cycle of $2 r$ edges between $C_1$ and $C_2$. Every vertex has degree $(r - 1) + 2 = r + 1 > n_H / 2$, so by Theorem II.A(i), $kappa'(H_(2 r)) = r + 1 > n_H / 2$ (highly connected), while $S = C_1$ achieves $|E(C_1, C_2)| / |C_1| = 2 r / r = 2$ exactly.

  For (iii), applying the discrete Cheeger inequality @chung1997spectral to the unnormalized Laplacian $L_H$ gives $h(H) <= sqrt(2 Delta(H) lambda_2 (L_H))$, so $lambda_2 (L_H) >= h(H)^2 / (2 Delta(H)) >= 4 / (2(n_H - 1)) = 2 / (n_H - 1)$.
]

#subsec("C", "Non-Uniqueness Under Minimum-Cut Tie-Breaking")

Because a graph can admit multiple distinct minimum cuts of identical minimum cardinality, pure HCS is not invariant to which minimum cut an arbitrary min-cut solver returns. We exhibit a minimal 7-vertex graph proving that tie-breaking alone can destroy an entire true cluster.

#theorem(title: "Minimum-Cut Non-Uniqueness and Cluster Shattering Counterexample", origin: "Original Counterexample")[
  Let $G_7$ be the 7-vertex, 11-edge simple graph formed by a clique $A = {1, 2, 3, 4} tilde.op K_4$, a triangle $B = {5, 6, 7} tilde.op K_3$, and two bridge edges $(4, 5)$ and $(4, 6)$. Then $kappa'(G_7) = 2 <= 7 / 2$, and two valid minimum cuts of size $2$ yield strictly incompatible HCS outputs:
  (i) Cutting $(A, B)$ produces the two-cluster partition $cal(C)_1 = {{1, 2, 3, 4}, {5, 6, 7}}$ with zero singletons;
  (ii) Cutting $({7}, V without {7})$ forces pure HCS—under every subsequent valid min-cut choice—to output $cal(C)_2 = {{1, 2, 3, 4}}$ and three singletons ${5}, {6}, {7}$, destroying $B$ completely.
]

#proof[
  Degrees in $G_7$ are $"deg"(1) = "deg"(2) = "deg"(3) = 3$, $"deg"(4) = 5$, $"deg"(5) = "deg"(6) = 3$, and $"deg"(7) = 2$. Since $A tilde.op K_4$ is 3-edge-connected, $B tilde.op K_3$ is 2-edge-connected, and two edges bridge $A$ and $B$, no single edge disconnects $G_7$; thus $kappa'(G_7) = delta(G_7) = 2 <= 7 / 2$.

  *Branch 1:* The cut $E(A, B) = {(4, 5), (4, 6)}$ has size $2 = kappa'(G_7)$. Its two induced sides are $G_7[A] tilde.op K_4$ ($kappa' = 3 > 4 / 2$) and $G_7[B] tilde.op K_3$ ($kappa' = 2 > 3 / 2$). Both pass the HCS certificate immediately, outputting $cal(C)_1 = {A, B}$ and $cal(S)_1 = emptyset$.

  *Branch 2:* The vertex-isolation cut $E({7}, V without {7}) = {(5, 7), (6, 7)}$ also has size $2 = kappa'(G_7)$. Vertex $7$ becomes a singleton. On the remaining 6-vertex graph $H = G_7[{1, 2, 3, 4, 5, 6}]$, vertices $5$ and $6$ now have degree $2$, so $kappa'(H) = 2 <= 6 / 2$ (not highly connected). Every minimum cut of $H$ has size $2$, and there are exactly three:
  - $E({5, 6}, A) = {(4, 5), (4, 6)}$: leaves $A tilde.op K_4$ (highly connected) and $K_2 = {5, 6}$ ($kappa'(K_2) = 1 <= 2 / 2$, which splits into singletons ${5}, {6}$);
  - $E({5}, V_H without {5})$: peels singleton ${5}$, leaving vertex $6$ with degree $1$ (only edge $(4, 6)$), so the unique next min-cut peels singleton ${6}$ and emits $A tilde.op K_4$;
  - $E({6}, V_H without {6})$: symmetric, yielding $A tilde.op K_4$ and singletons ${5}, {6}$.
  In every sub-branch, Branch 2 terminates at $cal(C)_2 = {A}$ and $cal(S)_2 = {{5}, {6}, {7}}$, shattering $42.9\%$ of vertices into singletons.
]

#subsec("D", "Degree-Bounded Shattering on Sparse Metric kNN Graphs")

#proposition(title: "Upper Bound on HCS Cluster Size Under Bounded Graph Degree", origin: "Original Derivation")[
  Let $G$ be any simple graph with maximum degree $Delta(G) = d_("max")$ (e.g., an unweighted symmetric $k$NN graph on a point cloud with $d_("max") = cal(O)(k)$). Then every cluster $C$ output by pure HCS satisfies the strict cardinality cap
  $ |C| <= 2 d_("max") - 1. $
  Consequently, any true ground-truth cluster of size $n_c > 2 d_("max") - 1$—such as an elongated manifold or large dense cloud in a fixed-$k$ nearest-neighbor graph—must shatter into at least $ceil(n_c / (2 d_("max") - 1))$ sub-clusters or singletons.
]

#proof[
  Let $H = G[C]$ be any non-trivial HCS cluster of size $n_H = |C|$. By Theorem II.A(i) and monotonicity of degree under taking induced subgraphs,
  $ floor(n_H / 2) + 1 <= kappa'(H) <= delta(H) <= Delta(H) <= Delta(G) = d_("max"). $
  Therefore $floor(n_H / 2) <= d_("max") - 1$, which implies $n_H <= 2(d_("max") - 1) + 1 = 2 d_("max") - 1$. In particular, on a sparse fixed-$k$ graph with $d_("max") = cal(O)(k)$ independent of $n$, cluster sizes are uniformly bounded by $cal(O)(k)$ regardless of the true cluster cardinality $n_c$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

*Time and Space Complexity.* Let $f(n, m)$ denote the runtime of a single global minimum edge-cut subroutine. Using the Stoer-Wagner algorithm with Fibonacci heaps, $f(n, m) = cal(O)(m n + n^2 log n)$ (or $cal(O)(m n)$ for unweighted graphs via Nagamochi-Ibaraki / Matula maximum-adjacency orderings @hartuv2000hcs). Note that by Theorem II.A(i), certifying a leaf that is *already* highly connected takes only $cal(O)(n_H)$ via $delta(H) > n_H / 2$, saving leaf min-cut calls; however, every internal split still invokes a global min-cut ($N_("cut") <= n - 1$ splits by Theorem II.B), giving worst-case runtime $cal(O)(N_("cut") (m n + n^2 log n))$ and $cal(O)(m + n)$ space. Because exact min-cut subroutines execute $n_H$ contraction phases per split, HCS is substantially slower than sparse greedy or Lanczos modularity methods on large graphs.

*Geometric Inductive Bias.* HCS enforces an *ultradense, compact, spherical-in-graph-metric* inductive bias: every cluster has graph diameter $<= 2$ and density $> 1 / 2$ (Theorem II.A). On similarity threshold graphs (where edges connect pairs within radius $tau$), HCS recovers tight, well-separated balls of diameter $<= 2 tau$ while cleanly shedding low-degree boundary noise as singletons. Unlike modularity methods, HCS's stopping criterion $kappa'(H) > |V_H| / 2$ is strictly local to $H$ and completely *immune to the global resolution limit* @fortunato2007resolution. However, by Proposition II.D, it requires dense intra-cluster connectivity ($d_("max") > n_c / 2$) and breaks down on sparse fixed-$k$ graphs or elongated manifolds of graph diameter $> 2$.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$cal(O)(N_("cut") (m n + n^2 log n))$ via Stoer-Wagner; $N_("cut") <= n - 1$ splits],
    [Space Complexity], [$cal(O)(m + n)$ for graph and binary recursion stack],
    [Cluster Count $K$], [Automatic from exact threshold $kappa'(H) > |V_H| / 2$ ($<=> delta(H) > |V_H| / 2$)],
    [Determinism], [Depends on min-cut tie-breaking (Thm. II.C)],
    [Connectedness], [Stronger: diameter $<= 2$, density $> 1 / 2$, Hamiltonian (Thm. II.A)],
    [Invariances], [Permutation equivariant; strictly local (no resolution limit)],
    [Key Limitation], [High min-cut cost; shatters sparse $k$NN graphs ($|C| <= 2 d_("max") - 1$)],
  ),
  caption: [Operational and complexity profile of HCS (Highly Connected Subgraphs).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *cDNA Microarray & Gene Expression Clustering:* Hartuv and Shamir's original application @hartuv2000hcs constructed similarity-threshold graphs over gene expression profiles ($n tilde.op 10^2$--$10^3$). Because co-regulated gene families form dense, nearly complete cliques with diameter $<= 2$ embedded in noisy backgrounds, HCS isolates tight co-expression modules while rejecting uninformative noise genes as singletons.
2. *Protein Sequence Homology & Ortholog Family Grouping:* Thresholding pairwise BLAST similarity scores produces dense homology clusters mixed with multi-domain "bridge" proteins. Min-cut splitting cleanly severs multi-domain bridges while guaranteeing every protein family obeys the diameter-2, $> 50\%$ density floor (Theorem II.A).
3. *High-Precision Entity Resolution and Record Deduplication:* In database deduplication graphs, true duplicate records form near-cliques whereas false-positive string matches form sparse bridges. HCS guarantees that merged entity blocks never suffer transitive "chaining" (unlike single-linkage clustering @sibson1973slink), because every pair in a cluster shares at least $3$ common matches.
4. *Financial Asset Core-Group Identification:* Thresholding correlation matrices above a high correlation cutoff isolates tightly synchronized asset baskets (density $> 0.5$) without merging small niche sectors into giant indices, exploiting HCS's total immunity to the modularity resolution limit.
5. *Small-to-Medium Compact Metric Clouds:* On well-separated low-dimensional datasets ($n <= 10^3$) converted to dense epsilon-neighborhood or generous $k$NN graphs ($k > n_c / 2$), HCS acts as a non-parametric mode isolator that automatically rejects boundary overlap points before optional singleton adoption.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

HCS degrades in four distinct regimes: (i) *sparse fixed-$k$ nearest-neighbor graphs* where $|C| <= 2 d_("max") - 1$ forces large clusters to shatter (Proposition II.D); (ii) *elongated manifolds, chains, and rings* whose intrinsic graph diameter vastly exceeds $2$; (iii) *graphs with peripheral low-degree attachments* where arbitrary minimum-cut tie-breaking peels singletons and shatters adjacent clusters unless low-degree pruning and singleton adoption are applied (Theorem II.C); and (iv) *large-scale networks ($n > 10^4$)* where repeated $cal(O)(m n)$ global min-cuts become computationally prohibitive.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Dense Similarity / Gene Graphs], [Very high: diam $<= 2$, density $> 1 / 2$], [Use singleton adoption @hartuv2000hcs],
    [Entity Deduplication], [Very high: prevents single-link chaining], [Resolve min-cut ties by balance],
    [Sparse Fixed-$k$ Manifolds], [Low: strict cap $|C| <= 2 d_("max") - 1$], [Use CNM (Alg. 115) or CHAMELEON @karypis1999chameleon],
    [High-Dim Noisy Feature Data], [Low: distance concentration blurs cuts], [Use subspace / spectral methods @ng2002spectral],
    [Large Networks ($n > 10^4$)], [Low: $cal(O)(N_("cut") m n)$ bottleneck], [Use Leiden @traag2019leiden or METIS @karypis1998metis],
  ),
  caption: [Domain suitability matrix for HCS (Highly Connected Subgraphs).],
)

#sec("V", "Structural Relations and Empirical Behavior")

HCS occupies a unique combinatorial corner in Category V: whereas spectral partitioning (Ratio Cut @hagen1992ratiocut, Normalized Cut @shi2000ncut, Leading Eigenvector Alg. 116) optimizes relative cut ratios and modularity algorithms (CNM Alg. 115, Spinglass Alg. 117) compare against a global random null, HCS enforces an *absolute, scale-invariant local connectivity certificate*.

#proposition(title: "Resolution-Limit Immunity and Separation from Clique Percolation", origin: "Original Proof")[
  (i) *Strict locality (no resolution limit):* For any graph $H$ and any auxiliary graph $G_("ext")$, running pure HCS on the disjoint union $H union.sq G_("ext")$ produces the exact same partition of $V_H$ as running pure HCS on $H$ alone;
  (ii) *Strict separation from Clique Percolation (CPM):* Neither HCS nor $k$-clique percolation @palla2005cpm subsumes the other: a chain of $r >= 4$ edge-overlapping triangles forms a single 3-clique CPM community yet is split by HCS ($"diam" >= 3$, $kappa' = 2 <= n / 2$), whereas the 6-vertex highly connected graph $H_6$ of Theorem II.B ($kappa'(H_6) = 4 > 3$) is accepted intact by HCS yet contains no 4-clique (failing $k$-CPM for all $k >= 4$).
]

#proof[
  For (i), Step 1 of HCS immediately splits any disconnected graph into connected components without invoking cross-component quantities; on component $H$, the stopping predicate $kappa'(H') > |V_(H')| / 2$ depends exclusively on $H'$, guaranteeing complete invariance to $G_("ext")$ and total immunity to the global resolution limit @fortunato2007resolution.

  For (ii), consider a chain of $r >= 4$ triangles sharing consecutive edges ($n = r + 2 >= 6$ vertices). Every consecutive pair of triangles shares a $K_2$, so the entire chain forms a single 3-clique community in CPM @palla2005cpm. However, its diameter is at least $3 > 2$ and its edge connectivity is $kappa' = 2 <= n / 2$, violating Theorem II.A so HCS splits it. Conversely, on the 6-vertex highly connected graph $H_6$ from Theorem II.B (two triangles $C_1, C_2$ joined by a 6-cycle matching, $delta(H_6) = kappa'(H_6) = 4 > 3$), $H_6$ is accepted intact by HCS, whereas no 4-clique exists anywhere in $H_6$ (every vertex's neighborhood induces the triangle-free 4-cycle $C_4$), separating HCS strictly from higher-order clique percolation as well as from edge-betweenness division @girvan2002community.
]

*Empirical Behavior on ClusterOverview Benchmarks.* Across the ClusterOverview 4-dataset suite, HCS obtains *Iris* $"Sil" = 0.477$, $"ARI" = 0.564$; *Customers* $"Sil" = 0.319$; *Multi-Shape* $"Sil" = 0.543$; and *High-D* $"Sil" = 0.254$, with a mean runtime of $17.121$ ms. These measurements crystallize every analytical property established above: (1) HCS achieves the *highest Iris Silhouette* ($0.477$) of all four Category V algorithms (ahead of Spinglass $0.411$, CNM $0.323$, and Leading Eigenvector $0.309$) and strong ARI ($0.564$), because compact 50-point biological clusters naturally meet the diameter-2, density-$> 1 / 2$ certificate once boundary overlap edges are severed by min-cuts; (2) on *Multi-Shape* ($0.543$), HCS cleanly isolates dense convex blobs—far outperforming Leading Eigenvector ($0.267$) and Spinglass ($0.088$)—while lagging CNM ($0.782$) due to partial shattering of elongated shapes under Theorem II.A ($"diam" <= 2$) and Proposition II.D ($|C| <= 2 d_("max") - 1$); (3) on *High-D* ($0.254$), high-dimensional distance concentration equalizes pairwise similarities and degrades min-cut separation; and (4) repeated exact Stoer-Wagner min-cuts make HCS $5 times$ slower than CNM ($3.417$ ms) and $31 times$ slower than Spinglass ($0.546$ ms), positioning HCS as a high-precision structural extractor for compact, moderate-size similarity graphs.
