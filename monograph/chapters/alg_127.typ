#import "../template/ieee.typ": *

#alg_entry(
  id: 127,
  name: "CPM (Clique Percolation Method)",
  category: "Cat. V: Graph, Spectral & Network",
  year: "2005",
  paradigm: "Overlapping k-Clique Graph Percolation",
  complexity: $cal(O)(d n 3^(d / 3) + c_k^2 k)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Most partition-based network clustering methods assign each vertex to a single
module, forcing bridge vertices in biological, social, and information networks
into artificial hard cuts. The Clique Percolation Method (CPM), introduced by
Palla et al. @palla2005cpm, formalizes overlapping community discovery through
local topological cohesion: communities are maximal regions reachable by
"rolling" a complete subgraph of order $k$ across shared $(k - 1)$-cliques.

Let $G = (V, E)$ be a simple undirected graph with $abs(V) = n$ vertices,
$abs(E) = m$ edges, and fixed integer clique order $k >= 3$. A *$k$-clique* is
a subset $K subset.eq V$ with $abs(K) = k$ inducing the complete graph $K_k$.
Let $cal(K)_k (G)$ denote the set of all $k$-cliques in $G$, and let
$cal(M)_k (G) = {M_1, dots, M_r}$ denote the family of all *maximal* cliques of
cardinality $abs(M_i) >= k$.

#definition(title: "k-Clique Adjacency, Percolation Cover, and Overlap Matrix")[
  Two $k$-cliques $K_1, K_2 in cal(K)_k (G)$ are *$k$-clique adjacent*, denoted
  $K_1 tilde.op_k K_2$, if and only if $abs(K_1 inter K_2) = k - 1$. The
  *$k$-clique adjacency graph* is $H_k (G) = (cal(K)_k (G), E_k)$ with
  ${K_1, K_2} in E_k <=> K_1 tilde.op_k K_2$. A *$k$-clique percolation
  community* is the vertex union
  $ C_ell = union.big_(K in cal(C)_ell) K subset.eq V, $
  where $cal(C)_ell$ is a connected component of $H_k (G)$. The *clique--clique
  overlap matrix* $O in bb(Z)_(>= 0)^(r times r)$ on $cal(M)_k (G)$ is defined
  by $O_(i j) = abs(M_i inter M_j)$, and its $(k - 1)$-thresholded adjacency
  matrix is $B^((k))_(i j) = bb(I)(O_(i j) >= k - 1)$ for $i != j$ (with zero
  diagonal).
]

Rather than optimizing a global scalar objective such as modularity
@newman2004finding @brandes2008modularity, CPM defines communities as a pure
simplicial connectivity invariant of $G$. Vertices outside every $k$-clique
remain unassigned (structural background), whereas a vertex lying in multiple
maximal cliques whose pairwise intersections have cardinality at most $k - 2$
belongs simultaneously to multiple overlapping communities. When $G$ is
weighted, edges are first filtered at a weight threshold $w_0$, or $k$-cliques
are retained only when their geometric mean edge weight exceeds an intensity
threshold $I_0$ @palla2005cpm.

#algbox(
  title: "Clique Percolation Method (CPM) via Maximal-Clique Overlap",
  inputs: [Simple graph $G = (V, E)$, integer clique size $k >= 3$],
  outputs: [Overlapping cover $cal(C)^((k)) = {C_1, dots, C_L}$ and noise set $V_"bg"$],
)[
  + Compute a degeneracy ordering of $V$ and enumerate all maximal cliques
    $cal(M)(G)$ via pivoted Bron--Kerbosch; retain
    $cal(M)_k (G) = {M_1, dots, M_r}$ satisfying $abs(M_i) >= k$.
  + Build inverted posting lists $P(v) = {i : v in M_i}$ for each $v in V$.
  + Accumulate non-zero off-diagonal overlaps $O_(i j) = abs(M_i inter M_j)$ by
    iterating over pairs $(i, j)$ co-occurring in $P(v)$ across all $v in V$.
  + Form the thresholded clique graph $B^((k))$ with edge $(i, j)$ whenever
    $i != j$ and $O_(i j) >= k - 1$.
  + Extract connected components ${cal(U)_1, dots, cal(U)_L}$ of $B^((k))$ via
    breadth-first search or union--find.
  + Project each component to $C_ell = union.big_(M_i in cal(U)_ell) M_i$ and
    set $V_"bg" = V without union.big_(ell=1)^L C_ell$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish four foundational results governing CPM: exact algebraic
equivalence between explicit $k$-clique rolling and maximal-clique overlap
thresholding, simplicial hierarchy and overlap structure, the classical random
graph percolation threshold, and exact locality paired with constructive
sparse-cut failure modes.

#subsec("A", "Equivalence of k-Clique Rolling and Maximal-Clique Thresholding")

Enumerating all $binom(s, k)$ $k$-subcliques inside a large clique of size $s$
is combinatorially prohibitive. The practical feasibility of CPM rests on the
following structural reduction.

#theorem(
  title: "Equivalence of H_k(G) and Maximal-Clique Overlap Connectivity",
  origin: "Original Proof",
)[
  For every simple graph $G = (V, E)$ and integer $k >= 3$, two $k$-cliques
  $K, K' in cal(K)_k (G)$ lie in the same connected component of $H_k (G)$ if
  and only if every pair of maximal cliques $M, M' in cal(M)_k (G)$ satisfying
  $K subset.eq M$ and $K' subset.eq M'$ lies in the same connected component of
  $B^((k))$. Consequently, projecting the connected components of $B^((k))$ to
  $V$ yields the exact family of $k$-clique percolation communities
  $cal(C)^((k))$.
]

#proof[
  We first establish an internal rolling lemma: if $K_a, K_b subset.eq M$ are
  any two $k$-cliques contained in the same maximal clique $M in cal(M)_k (G)$,
  then $K_a$ and $K_b$ are connected in $H_k (G)$. Write
  $K_a without K_b = {u_1, dots, u_t}$ and
  $K_b without K_a = {w_1, dots, w_t}$ where $t = k - abs(K_a inter K_b)$.
  Define $J_0 = K_a$ and $J_s = (J_(s-1) without {u_s}) union {w_s}$ for
  $s = 1, dots, t$. Since $J_s subset.eq M$ and $M$ is complete, every $J_s$ is
  a $k$-clique, $abs(J_(s-1) inter J_s) = k - 1$, and $J_t = K_b$, giving a
  length-$t$ path in $H_k (G)$.

  *($=>$)* Suppose $K = K_0, K_1, dots, K_p = K'$ is a path in $H_k (G)$. For
  each $s in {0, dots, p}$, choose a maximal clique $M_s in cal(M)_k (G)$ with
  $K_s subset.eq M_s$. For consecutive cliques,
  $ abs(M_s inter M_(s+1)) >= abs(K_s inter K_(s+1)) = k - 1, $
  so either $M_s = M_(s+1)$ or $(M_s, M_(s+1))$ is an edge of $B^((k))$.
  Moreover, for any $M supset.eq K$, $abs(M inter M_0) >= abs(K) = k > k - 1$,
  so $M$ and $M_0$ are equal or adjacent in $B^((k))$; symmetrically $M_p$ and
  $M'$ are equal or adjacent. Concatenating yields a walk from $M$ to $M'$ in
  $B^((k))$.

  *($<=$)* Conversely, suppose $M = M_0, M_1, dots, M_q = M'$ is a path in
  $B^((k))$. For each step $r in {0, dots, q - 1}$,
  $abs(M_r inter M_(r+1)) >= k - 1$, so we may pick a $(k - 1)$-subset
  $S_r subset.eq M_r inter M_(r+1)$. Since $abs(M_r), abs(M_(r+1)) >= k$, there
  exist vertices $x_r in M_r without S_r$ and $y_r in M_(r+1) without S_r$.
  Then $K_r^- = S_r union {x_r} subset.eq M_r$ and
  $K_r^+ = S_r union {y_r} subset.eq M_(r+1)$ are $k$-cliques satisfying
  $abs(K_r^- inter K_r^+) >= k - 1$, hence equal or adjacent in $H_k (G)$.
  Joining $K$ to $K_0^-$ inside $M_0$, each $K_(r-1)^+$ to $K_r^-$ inside
  $M_r$, and $K_(q-1)^+$ to $K'$ inside $M_q$ via the internal rolling lemma
  produces a path in $H_k (G)$ from $K$ to $K'$. Finally, every vertex $v in M$
  with $abs(M) >= k$ belongs to at least one $k$-subset of $M$, so the vertex
  unions of corresponding components of $B^((k))$ and $H_k (G)$ coincide
  identically.
]

#subsec("B", "Nested Hierarchy and Disjointness of (k - 1)-Skeletons")

Varying $k$ induces a strict hierarchical refinement across resolutions while
enforcing a sharp combinatorial barrier between overlapping communities.

#proposition(
  title: "Monotonic Nesting in k and (k - 1)-Face Partitioning",
  origin: "Original Proof",
)[
  Let $cal(C)^((k))$ be the family of $k$-clique communities of $G$ for
  $k >= 3$, and for each component $cal(C)_ell$ of $H_k (G)$ define its
  $(k - 1)$-skeleton
  $ Delta_(k-1) (cal(C)_ell) = { S subset.eq V : abs(S) = k - 1, exists K in cal(C)_ell, S subset K }. $
  Then:
  + *(Nested refinement)* Every community $C in cal(C)^((k+1))$ is contained in
    a unique community $C' in cal(C)^((k))$.
  + *(Face disjointness)* For distinct components $cal(C)_ell != cal(C)_(ell')$
    of $H_k (G)$, $Delta_(k-1) (cal(C)_ell) inter Delta_(k-1) (cal(C)_(ell')) = emptyset$.
    In particular, any two maximal cliques $M_1, M_2 in cal(M)_k (G)$ assigned
    to distinct communities satisfy $abs(M_1 inter M_2) <= k - 2$, and for
    $k = 3$ the triangle-bearing edges of $G$ are strictly partitioned.
]

#proof[
  *(1)* Let $cal(D)$ be a connected component of $H_(k+1) (G)$ with vertex
  union $C$. If $K_1, K_2 in cal(D)$ are adjacent $(k + 1)$-cliques, then
  $S = K_1 inter K_2$ is a $k$-clique. Any $k$-subset $J subset K_1$ differs
  from $S$ by at most one vertex, so $J = S$ or $abs(J inter S) = k - 1$,
  placing $J$ in the same component of $H_k (G)$ as $S$; the same holds for
  every $k$-subset of $K_2$. By induction along paths in $cal(D)$, all
  $k$-subsets of all $(k + 1)$-cliques in $cal(D)$ belong to a single unique
  connected component $cal(C)'$ of $H_k (G)$ with vertex union $C'$. Since
  every vertex of a $(k + 1)$-clique belongs to $k$ of its $k$-subsets,
  $C subset.eq C'$. Uniqueness of $C'$ holds because $cal(D)$ contains at least
  one $k$-clique, which uniquely determines $cal(C)'$.

  *(2)* Suppose for contradiction that
  $S in Delta_(k-1) (cal(C)_ell) inter Delta_(k-1) (cal(C)_(ell'))$. Then there
  exist $K_1 in cal(C)_ell$ and $K_2 in cal(C)_(ell')$ such that
  $S subset K_1$ and $S subset K_2$. Hence
  $abs(K_1 inter K_2) >= abs(S) = k - 1$, implying either $K_1 = K_2$ or
  ${K_1, K_2} in E_k$, which contradicts $cal(C)_ell != cal(C)_(ell')$. If
  maximal cliques $M_1, M_2$ in distinct components satisfied
  $abs(M_1 inter M_2) >= k - 1$, any $(k - 1)$-subset of $M_1 inter M_2$ would
  lie in both skeletons, again a contradiction. Setting $k = 3$, the
  $2$-skeleton $Delta_2 (cal(C)_ell)$ is precisely the set of edges
  participating in at least one triangle of $cal(C)_ell$, proving an exact
  partition of all triangle-bearing edges.
]

Note the subtlety in Part (2): two distinct $k$-clique communities $C_1$ and
$C_2$ *can* share an arbitrarily large set of bridge vertices (and even chord
edges not contained in triangles of both communities), provided no
$(k - 1)$-clique is shared by $k$-cliques from both components.

#subsec("C", "Critical Percolation Threshold on Random Graphs")

Choosing $k$ requires balancing fragmentation (large $k$) against giant-cluster
percolation (small $k$). Palla et al. @palla2005cpm established the sharp
Erdos--Renyi phase transition separating local clusters from a spanning giant
component.

#theorem(
  title: "Critical k-Clique Percolation Threshold in G(n, p)",
  origin: "Classical — Cited",
)[
  Let $G tilde.op G(n, p)$ be an Erdos--Renyi random graph on $n$ vertices with
  fixed clique order $k >= 3$ @palla2005cpm. In the sparse scaling
  $p = c n^(-1 / (k - 1))$, the critical edge probability for the emergence of
  a giant $k$-clique percolation cluster is
  $ p_c (k) = ((k - 1) n)^(-1 / (k - 1)) (1 + o(1)). $
]

#proof[
  Consider exploring a connected component of $H_k (G)$ via a branching
  process starting from an active $k$-clique $K$ entered across one of its
  $(k - 1)$-faces $S_"in"$. There remain $k - 1$ distinct unused
  $(k - 1)$-faces $S' subset K$. For each unused face $S'$ and each external
  vertex $v in V without K$ (of which there are $n - k$), the set
  $S' union {v}$ forms a newly discovered adjacent $k$-clique if and only if
  all $k - 1$ edges between $v$ and $S'$ are present, an event of probability
  $p^(k - 1)$.

  Could a single external vertex $v$ simultaneously attach to two distinct
  faces of $K$ (creating a $(k + 1)$-clique and double-counting offspring)?
  Such a vertex must be adjacent to all $k$ vertices of $K$, which occurs with
  probability $p^k$. Across all $n - k$ candidates, the expected number of such
  over-connected vertices is
  $ (n - k) p^k = c^k n^(1 - k / (k - 1)) (1 + o(1)) = cal(O)(n^(-1 / (k - 1))) = o(1). $
  Thus local neighborhoods in $H_k (G)$ are asymptotically tree-like. By
  linearity of expectation, the mean offspring number (expected adjacent
  $k$-cliques reachable from $K$ via its $k - 1$ unused faces) is
  $ mu_k (n, p) = (k - 1)(n - k) p^(k - 1). $
  The Galton--Watson criticality threshold $mu_k (n, p) = 1$ yields
  $ p^(k - 1) = 1 / ((k - 1)(n - k)) => p_c (k) = ((k - 1) n)^(-1 / (k - 1)) (1 + cal(O)(k / n)). $
  For $p < (1 - epsilon) p_c (k)$, the branching process is subcritical and
  components remain $cal(O)(1)$ in expectation; for $p > (1 + epsilon) p_c (k)$,
  a strictly positive survival probability produces a giant percolating
  component.
]

For practical calibration @palla2005cpm, selecting $k$ so that background edge
density satisfies $p_"out" < p_c (k)$ while intra-community density exceeds
$p_"in" > p_c (k)$ suppresses false percolation while preserving community
cohesion.

#subsec("D", "Exact Locality vs. Sparse-Cut Blindness and Chaining")

Unlike modularity maximization, which suffers from a global scale-dependent
resolution limit @fortunato2007resolution, CPM is completely local; however,
this locality carries two sharp structural pathologies.

#theorem(
  title: "Resolution-Limit Freedom and Constructive Failure Modes",
  origin: "Original Counterexample",
)[
  For every integer $k >= 3$:
  + *(Exact locality)* For any graph $G$ and any disjoint graph $G'$,
    $cal(C)^((k)) (G union.sq G') = cal(C)^((k)) (G) union.sq cal(C)^((k)) (G')$.
    More generally, attaching $G'$ to $G$ through edges belonging to no
    $k$-clique leaves $cal(C)^((k)) (G)$ unchanged.
  + *(Dense triangle-free blindness)* For every integer $s >= 2$, the
    disconnected graph $G_s = K_(s, s) union.sq K_(s, s)$ on $n = 4 s$ vertices
    has two zero-cut communities of internal density $s / (2 s - 1) > 1 / 2$
    and optimal modularity $Q = 1 / 2$, yet $cal(C)^((k)) (G_s) = emptyset$
    (0\% vertex coverage).
  + *(Single-chain collapse)* Linking two cliques $K_s, K_s$ ($s > k$) by a
    minimal $(k - 1)$-clique bridge of $k - 1$ shared vertices or a single
    rolling chain of $k$-cliques merges them into one community even as the
    inter-community conductance vanishes at rate $cal(O)(k^2 / s^2) -> 0$.
]

#proof[
  *(1)* A $k$-clique $K$ is complete and connected, so $K$ cannot contain
  vertices from both $G$ and $G'$ unless every cross-pair is an edge; when no
  bridge edge lies in a $k$-clique, $cal(K)_k (G union_B G') = cal(K)_k (G)
  union.sq cal(K)_k (G')$ with no cross-edges in $H_k$. Thus connected
  components of $H_k$ decouple identically, immunizing CPM against global
  network size $n$ and total edge count $m$.

  *(2)* Each complete bipartite component $K_(s, s)$ contains $s^2$ edges among
  $binom(2 s, 2) = s(2 s - 1)$ pairs, giving internal density
  $s / (2 s - 1) > 1 / 2$ and zero inter-component cut; the modularity of the
  exact 2-partition is $Q = 2 (1 / 2 - (1 / 2)^2) = 1 / 2$. Because $K_(s, s)$
  is bipartite, it contains no odd cycles and therefore no triangles, so its
  clique number satisfies $omega(G_s) = 2 < k$. Hence $cal(K)_k (G_s) = emptyset$
  and $cal(C)^((k)) (G_s) = emptyset$.

  *(3)* Take two disjoint cliques $K_s^((1)), K_s^((2))$ ($s > k$), select
  $(k - 1)$-subsets $S_1 = {u_1, dots, u_(k-1)} subset K_s^((1))$ and
  $S_2 = {w_1, dots, w_(k-1)} subset K_s^((2))$, and add the
  $binom(k, 2) = k(k - 1) / 2$ cross-edges ${{u_i, w_j} : 1 <= j <= i <= k - 1}$.
  For each $t = 1, dots, k - 1$, the set
  $K_t = {u_t, dots, u_(k-1), w_1, dots, w_t}$ is a $k$-clique: intra-$S_1$
  and intra-$S_2$ pairs are edges inside $K_s^((1))$ and $K_s^((2))$, while
  every cross-pair satisfies $j <= t <= i$. Moreover, $K_1$ shares $S_1$ with
  $K_s^((1))$, $K_t$ shares $k - 1$ vertices with $K_(t+1)$, and $K_(k-1)$
  shares $S_2$ with $K_s^((2))$. Thus $K_s^((1))$ and $K_s^((2))$ collapse
  into a single $k$-clique community even though the cut conductance obeys
  $ phi(K_s^((1))) = (k(k - 1)) / (2 s(s - 1) + k(k - 1)) = cal(O)(k^2 / s^2) -> 0 $
  as $s -> oo$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

The computational bottleneck of CPM is maximal-clique enumeration followed by
clique-overlap construction. On arbitrary graphs, Moon--Moser graphs achieve
$Theta(3^(n / 3))$ maximal cliques, making worst-case enumeration exponential in
$n$. However, real sparse networks exhibit low *degeneracy* $d$ (the smallest
integer such that every induced subgraph has a vertex of degree at most $d$).
Using degeneracy ordering with Tomita pivoting inside Bron--Kerbosch, maximal
clique enumeration runs in $cal(O)(d n 3^(d / 3))$ time, which is fixed-parameter
tractable in $d$ and linear in $n$ for bounded-degeneracy graphs.

Let $c_k = abs(cal(M)_k (G))$ be the number of maximal cliques of size at least
$k$, and $s_"max" = omega(G)$ the clique number. Building inverted posting
lists takes $cal(O)(sum_(i=1)^(c_k) abs(M_i)) = cal(O)(c_k s_"max")$ time and
memory. Accumulating the sparse overlap matrix $O$ via posting lists requires
$cal(O)(sum_(v in V) abs(P(v))^2)$ time, bounded in the worst case by
$cal(O)(c_k^2 k)$ (or $cal(O)(c_k^2 s_"max")$ if dense), and far smaller when
vertex clique-membership counts $abs(P(v))$ are bounded. Connected components
of $B^((k))$ take $cal(O)(c_k + abs(E(B^((k)))))$ time via union--find.

Geometrically, CPM is invariant under graph isomorphism and edge-order
permutations, requires no random seed (deterministic output), and does not fix
the number of communities $L$ a priori. Its sole discrete hyperparameter is $k$
(plus an optional edge-weight threshold $w_0$ on weighted graphs). Memory is
dominated by storing $cal(M)_k (G)$ and the sparse overlap graph, scaling as
$cal(O)(n + m + c_k s_"max" + abs(E(B^((k)))))$.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Worst-case time], [$cal(O)(3^(n / 3) + c_k^2 k)$ (Moon--Moser graphs)],
    [Parameterized time], [$cal(O)(d n 3^(d / 3) + sum_v abs(P(v))^2)$ in degeneracy $d$],
    [Space complexity], [$cal(O)(n + m + c_k s_"max" + abs(E(B^((k)))))$],
    [Partition type], [Soft overlapping cover with unassigned background],
    [Determinism], [Fully deterministic; invariant to vertex labeling],
    [Resolution limit], [Absent (strictly local simplicial invariant)],
    [Hyperparameters], [Clique order $k >= 3$ (and optional weight cut $w_0$)],
    [Parallelism], [Independent Bron--Kerbosch subtrees & posting lists],
  ),
  caption: [Operational and complexity profile of CPM (Alg. 127).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

CPM excels whenever communities are defined by tight multi-way cohesion and
shared boundary entities play dual functional roles across modules:

+ *Protein--Protein Interaction (PPI) Complexes @palla2005cpm:* Multiprotein
  molecular machines form densely interacting physical complexes (near-cliques
  of order $k in {3, 4, 5}$), while pleiotropic "moonlighting" proteins bind
  multiple complexes sequentially. CPM's overlapping cover recovers multi-complex
  proteins without forcing artificial hard splits, while leaving noisy
  high-throughput yeast two-hybrid false positives in $V_"bg"$.
+ *Co-Authorship and Scientific Collaboration Networks @palla2005cpm:* Multi-author
  papers induce literal cliques in collaboration graphs. A researcher bridging
  two disciplines participates in cliques from both fields whose pairwise
  intersection is only ${v}$ (size $1 <= k - 2$), matching the overlap
  geometry of Section II-B exactly.
+ *Semantic Word-Association and Concept Graphs @palla2005cpm:* Polysemous
  terms (e.g., "bright" linking intelligence, optics, and astronomy) sit at the
  intersection of dense semantic cliques. Rolling $k$-cliques cleanly separates
  distinct word senses while keeping the polysemous hub in every relevant sense
  cluster.
+ *Financial Contagion and Syndicated Lending Cartels:* Banks co-investing in
  syndicated loans form dense weighted cliques; filtering by loan volume $w_0$
  and running $k$-clique percolation isolates overlapping institutional risk
  blocs and systemic bridge banks.
+ *Functional Brain Connectivity (fMRI Regions):* Thresholded correlation
  networks exhibit strongly clustered cortical sub-networks (e.g., default-mode
  and frontoparietal control modules) linked by polymodal hub regions; CPM
  identifies overlapping functional assemblies without imposing modularity's
  global scale bias @fortunato2007resolution.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Three operational failure regimes require caution: (i) *locally tree-like or
bipartite graphs* (e.g., power grids, road networks, user--item bipartite
graphs), where triangle deficiency causes near-zero coverage (Section II-D);
(ii) *supercritical dense cores* with edge density $p > p_c (k)$, where
percolation merges distinct semantic modules into one giant hairball; and
(iii) *high-degeneracy graphs* containing massive near-cliques (e.g., dense
correlation matrices without sparsification), where $c_k$ explodes
combinatorially.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [PPI complexes], [High: dense cliques + moonlighting overlaps], [Filter noisy edges; try OSLOM @lancichinetti2011oslom],
    [Co-authorship], [High: papers generate exact cliques], [Tune $k$ above $p_c (k)$; try BIGCLAM @yang2013bigclam],
    [Word-sense graphs], [High: polysemous hubs at clique joins], [Large $k$ drops rare senses; try SLPA @xie2011slpa],
    [Bipartite / grid nets], [Poor: $omega(G) = 2 < k$ gives empty cover], [Use SBM (Alg. 130) or Louvain (Alg. 129)],
    [Dense $k$-NN graphs], [Moderate: local cliques but chaining risk], [High $d$ slows BK; use Leiden @traag2019leiden],
  ),
  caption: [Domain suitability and failure-mode matrix for CPM (Alg. 127).],
)

#sec("V", "Structural Relations and Empirical Behavior")

CPM occupies a clean topological niche connecting hierarchical single-linkage
clustering @sibson1973slink @gower1969minimum to modern overlapping community
detectors such as DEMON (Alg. 128) @coscia2012demon, OSLOM
@lancichinetti2011oslom, SLPA @xie2011slpa, and BIGCLAM @yang2013bigclam.

#proposition(
  title: "CPM as Single-Linkage Clustering on the k-Clique Complex",
  origin: "Original Derivation",
)[
  Define the normalized overlap dissimilarity on $cal(K)_k (G)$ by
  $ d_k (K_1, K_2) = 1 - abs(K_1 inter K_2) / k in {0, 1 / k, dots, 1}. $
  Then the connected components of $H_k (G)$ are precisely the single-linkage
  clusters of $(cal(K)_k (G), d_k)$ cut at threshold $tau = 1 / k$.
  Consequently, CPM inherits single-linkage's exact chaining behavior on
  $cal(K)_k (G)$ while discarding all vertices outside $union.big cal(K)_k (G)$.
]

#proof[
  For distinct $K_1, K_2 in cal(K)_k (G)$, $abs(K_1 inter K_2) <= k - 1$, so
  $d_k (K_1, K_2) >= 1 / k$ with equality if and only if
  $abs(K_1 inter K_2) = k - 1$, i.e., ${K_1, K_2} in E_k$. Single-linkage at
  height $tau = 1 / k$ forms connected components of the threshold graph
  ${{K_1, K_2} : d_k (K_1, K_2) <= 1 / k}$, which is identically $H_k (G)$.
]

Relative to Fast Unfolding (Alg. 129) @blondel2008louvain, CPM avoids the
global resolution limit and supports vertex overlap, at the expense of partial
coverage and sensitivity to clique chaining. Relative to DEMON (Alg. 128)
@coscia2012demon, both operate on dense triangle-rich local neighborhoods, but
CPM is deterministic and requires rigid $k$-cliques, whereas DEMON relaxes
cliques to majority-voted ego-network clusters.

*Empirical Benchmark Analysis (ClusterOverview 4-Dataset Suite).* Evaluated on
$k$-NN proximity graphs constructed from the four ClusterOverview benchmark
datasets (measuring Silhouette @rousseeuw1987silhouettes and Adjusted Rand
Index @hubert1985comparing):
- *Iris (`Sil = 0.586`, `ARI = 0.567`):* CPM ties Fast Unfolding (Alg. 129) in
  both Silhouette and ARI, cleanly isolating the dense Setosa $k$-NN clique
  component while experiencing partial bridging across the overlapping
  Versicolor--Virginica boundary (an instance of single-linkage clique chaining).
- *Customers (`Sil = 0.262`):* CPM outperforms both Fast Unfolding (`0.190`)
  and DEMON (`0.250`) by carving out dense transactional cores rather than
  merging small niches under a global modularity null model, though trailing
  SBM (`0.356`), which models heterogeneous core--periphery structure directly.
- *Multi-Shape (`Sil = 0.554`, `ARI = N/A`):* Local $k$-clique rolling tracks
  curved non-convex density ridges substantially better than DEMON (`0.396`),
  yet falls behind Fast Unfolding and SBM (`0.713`) where narrow manifold necks
  either pinch below $k$-clique width or chain adjacent shapes.
- *High-D (`Sil = 0.188`, `ARI = N/A`):* Concentration of Euclidean distances
  @beyer1999nearest degrades $k$-NN clique regularity, matching Fast Unfolding
  (`0.188`) and trailing SBM (`0.303`).
- *Runtime (`Mean = 2.800 ms`):* On sparse proximity graphs of modest
  degeneracy, pivoted Bron--Kerbosch plus inverted posting-list overlap runs
  $1.6 times$ faster than DEMON (`4.424 ms`) and $46.6 times$ faster than SBM
  (`130.533 ms`), while remaining $4.1 times$ slower than Fast Unfolding
  (`0.681 ms`) due to clique enumeration overhead.
