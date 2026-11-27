#import "../template/ieee.typ": *

#alg_entry(
  id: 107,
  name: "Louvain Method",
  category: "Cat. V: Graph, Spectral & Network",
  year: "2008",
  paradigm: "Hierarchical Greedy Modularity Optimization",
  complexity: $O(L m)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Community detection in large weighted networks seeks a hard vertex
partition whose internal link density significantly exceeds that of a
degree-preserving random null graph. Let $G = (V, E, w)$ be an
undirected weighted graph on $n = |V|$ vertices and $m = |E|$ edges,
encoded by a symmetric nonnegative adjacency matrix
$A = (A_(i j)) in RR_(>= 0)^(n times n)$ where diagonal entries
$A_(i i) >= 0$ store total self-loop weight (twice the undirected
loop multiplicity, matching the two directed endpoints of every
non-loop edge). Define vertex strength $k_i = sum_(j = 1)^n A_(i j)$
and total edge weight $2 m_w = sum_(i = 1)^n k_i = sum_(i, j) A_(i j)$,
assumed strictly positive.

The *Louvain method* @blondel2008louvain maximizes the
resolution-parameterized Newman--Girvan modularity functional
@newman2004finding @reichardt2006spinglass through a two-phase
multilevel ascent that alternates greedy single-vertex moves with
exact graph coarsening.

#definition(title: "Resolution-Parameterized Modularity Functional")[
  Let $cal(P) = {C_1, dots, C_K}$ be a partition of $V$ with cluster
  indicator matrix $H in {0, 1}^(n times K)$ ($H_(i c) = 1$ iff
  $i in C_c$, denoted $c_i = c$), and let $gamma > 0$ be a resolution
  parameter. Define the *modularity matrix*
  $B_gamma = A - (gamma / (2 m_w)) bold(k) bold(k)^top$, where
  $bold(k) = (k_1, dots, k_n)^top$. The modularity of $cal(P)$ is
  $ Q_gamma (cal(P)) &= 1 / (2 m_w) sum_(i, j in V)
    (A_(i j) - gamma (k_i k_j) / (2 m_w)) delta(c_i, c_j) \
    &= 1 / (2 m_w) "Tr"(H^top B_gamma H) \
    &= sum_(c = 1)^K [ (Sigma_("in") (c)) / (2 m_w)
    - gamma ((Sigma_("tot") (c)) / (2 m_w))^2 ], $
  where $Sigma_("in") (c) = sum_(i, j in C_c) A_(i j)$ is the total
  internal weight of $C_c$ and
  $Sigma_("tot") (c) = sum_(i in C_c) k_i$ is its total incident
  strength. Setting $gamma = 1$ recovers classical modularity
  @newman2004finding.
]

Because $Q_gamma (cal(P))$ is invariant under coarsening communities
into weighted super-nodes with self-loops (Theorem II.2), Louvain
escapes single-vertex local traps by lifting each converged partition
into a compressed super-graph and repeating greedy local ascent until
no further improvement exists.

#algbox(
  title: "Louvain Multilevel Modularity Optimization",
  inputs: [Graph $G = (V, E, w)$, resolution $gamma > 0$, tolerance $epsilon > 0$.],
  outputs: [Hierarchical community partitions ${cal(P)^((h))}$ and final partition $cal(P)^*$.],
)[
  1. Set $G^((0)) <- G$, level $h <- 0$, and mapping $pi^((0)) (i) = i$.
  2. *Phase 1 (Greedy Local Moving)* on $G^((h)) = (V^((h)), A^((h)))$:
     initialize singletons $c_i <- i$, $Sigma_("tot") (i) <- k_i^((h))$.
  3. Repeat sweeps over $V^((h))$ (in cyclic or random order): for each
     $i$, compute $Delta Q_gamma (i : c_i -> c)$ for all communities
     $c$ adjacent to $i$ plus an empty community $emptyset$; set
     $c^* <- arg max_c Delta Q_gamma (i : c_i -> c)$ and move $i$ to
     $c^*$ if $Delta Q_gamma (i : c_i -> c^*) > epsilon$, updating
     $Sigma_("tot")$ in $O(1)$ time.
  4. Terminate Phase 1 when a complete sweep yields no move. If every
     community has size $1$, stop and return the lifted partition
     $cal(P)^*$ on $V$.
  5. *Phase 2 (Exact Graph Aggregation)*: let $cal(P)^((h)) = {C_1, dots, C_K}$
     be the Phase-1 partition. Form $G^((h + 1))$ on $K$ super-nodes
     with edge weights
     $A_(c d)^((h + 1)) = sum_(i in C_c) sum_(j in C_d) A_(i j)^((h))$,
     including self-loops $A_(c c)^((h + 1)) = Sigma_("in")^((h)) (c)$.
  6. Update $pi^((h + 1)) <- pi_("agg") compose pi^((h))$, set
     $h <- h + 1$, and return to Step 2.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish exact incremental gain identities, coarsening invariance,
spectral upper bounds, resolution-limit thresholds, and an explicit
constructive counterexample demonstrating Louvain's disconnected
community pathology.

#subsec("A", "Exact Incremental Gain and Coarsening Invariance")

#lemma(
  title: "Exact O(deg(i)) Modularity Gain With Self-Loops",
  origin: "Original Derivation",
)[
  Let $i in C_a$ be reassigned to a distinct community $C_b$ ($b != a$,
  possibly empty). Define $k_(i, c) = sum_(j in C_c without {i}) A_(i j)$
  and $Sigma_("tot")^((-i)) (a) = Sigma_("tot") (a) - k_i$. Even in the
  presence of arbitrary self-loops $A_(i i) >= 0$, the exact change in
  $Q_gamma$ is
  $ Delta Q_gamma (i : C_a -> C_b)
    &= (k_(i, b) - k_(i, a)) / m_w \
    &+ (gamma k_i) / (2 m_w^2)
    (Sigma_("tot")^((-i)) (a) - Sigma_("tot") (b)). $
]

#proof[
  Only communities $C_a$ and $C_b$ change their contributions to
  $Q_gamma (cal(P))$. Using symmetry $A_(i j) = A_(j i)$, removing $i$
  from $C_a$ changes its internal weight from
  $Sigma_("in") (a) = Sigma_("in")^((-i)) (a) + 2 k_(i, a) + A_(i i)$
  to $Sigma_("in")^((-i)) (a)$, while inserting $i$ into $C_b$
  increases its internal weight by $2 k_(i, b) + A_(i i)$. The
  self-loop terms $+A_(i i)$ and $-A_(i i)$ cancel identically, yielding
  the linear edge-weight change
  $ Delta Q_("lin") = (2 k_(i, b) - 2 k_(i, a)) / (2 m_w)
    = (k_(i, b) - k_(i, a)) / m_w. $
  For the null-model term, removing $i$ from $C_a$ reduces
  $(Sigma_("tot") (a))^2$ from $(Sigma_("tot")^((-i)) (a) + k_i)^2$
  to $(Sigma_("tot")^((-i)) (a))^2$, a change of
  $-(2 k_i Sigma_("tot")^((-i)) (a) + k_i^2)$. Inserting $i$ into
  $C_b$ raises $(Sigma_("tot") (b))^2$ by
  $+(2 k_i Sigma_("tot") (b) + k_i^2)$. Again $k_i^2$ cancels
  identically, giving
  $ Delta Q_("null") = - gamma / (4 m_w^2)
    [2 k_i Sigma_("tot") (b) - 2 k_i Sigma_("tot")^((-i)) (a)]. $
  Summing $Delta Q_("lin") + Delta Q_("null")$ yields the claimed
  formula. Computing $k_(i, c)$ for all adjacent communities requires
  visiting each incident edge of $i$ once, costing $O("deg"(i))$ time.
]

#theorem(
  title: "Exact Modularity Preservation Under Phase-2 Aggregation",
  origin: "Original Proof",
)[
  Let $cal(P)_1 = {C_1, dots, C_K}$ be any partition of $G$ with
  indicator matrix $H_1 in {0, 1}^(n times K)$, and let $G'$ be the
  Phase-2 coarsened graph with adjacency matrix $A' = H_1^top A H_1$.
  For every partition $cal(P)'$ of the $K$ super-nodes with indicator
  $H' in {0, 1}^(K times M)$, let $cal(P) = pi^(-1) (cal(P)')$ be the
  lifted partition on $V$ (indicator $H = H_1 H'$). Then
  $ Q_gamma^(G') (cal(P)') = Q_gamma^(G) (cal(P)). $
  Consequently, Louvain executes strict monotone ascent on $Q_gamma$
  across all levels and terminates in finitely many moves.
]

#proof[
  Because $H_1 bold(1)_K = bold(1)_n$, the super-node strength vector
  satisfies
  $ bold(k)' = A' bold(1)_K = H_1^top A H_1 bold(1)_K
    = H_1^top A bold(1)_n = H_1^top bold(k), $
  so $2 m'_w = bold(1)_K^top bold(k)' = bold(1)_n^top bold(k) = 2 m_w$
  is strictly preserved. Substituting $A' = H_1^top A H_1$ and
  $bold(k)' = H_1^top bold(k)$ into the super-graph modularity matrix
  gives
  $ B'_gamma = H_1^top A H_1 - gamma / (2 m_w) H_1^top bold(k)
    bold(k)^top H_1 = H_1^top B_gamma H_1. $
  Evaluating the trace identity yields
  $ Q_gamma^(G') (cal(P)') &= 1 / (2 m'_w) "Tr"((H')^top B'_gamma H') \
    &= 1 / (2 m_w) "Tr"((H_1 H')^top B_gamma (H_1 H'))
    = Q_gamma^(G) (cal(P)). $
  Every accepted Phase-1 move strictly increases $Q_gamma$ by more than
  $epsilon > 0$, while Phase-2 coarsening preserves $Q_gamma$ exactly.
  Because the lattice of partitions of $V$ is finite (cardinality equal
  to the $n$-th Bell number $B_n$), no partition can repeat and
  Louvain halts in finitely many steps.
]

#subsec("B", "Invariances and Spectral Upper Bound")

#proposition(
  title: "Scale Invariance, Equivariance, and Spectral Bound",
  origin: "Original Proof",
)[
  (i) For every scalar $alpha > 0$ and permutation matrix $P$,
  $Q_gamma$ is invariant under $A |-> alpha A$ and equivariant under
  $A |-> P A P^top$. (ii) For $gamma = 1$, let $lambda_("max") (B_1)$
  be the largest eigenvalue of $B_1 = A - (1 / (2 m_w)) bold(k) bold(k)^top$.
  Every $K$-community partition $cal(P)$ satisfies
  $ Q_1 (cal(P)) <= n / (2 m_w) (1 - 1 / K) lambda_("max") (B_1), $
  which is tight at $K = 1$ and matches Newman's spectral relaxation
  @newman2006leading at $K = 2$.
]

#proof[
  (i) Scaling $A |-> alpha A$ multiplies $bold(k)$ and $m_w$ by
  $alpha$, so $(alpha k_i alpha k_j) / (2 alpha m_w) = alpha (k_i k_j) / (2 m_w)$
  and $B_gamma |-> alpha B_gamma$, canceling $1 / (2 alpha m_w)$ in
  $Q_gamma$. Relabeling $A |-> P A P^top$ maps $H |-> P H$, leaving
  $"Tr"(H^top B_gamma H)$ invariant.
  (ii) Observe that
  $(B_1 bold(1))_i = k_i - (k_i / (2 m_w)) sum_j k_j = 0$, so
  $B_1 bold(1) = bold(0)$. For each cluster $C_c$ with indicator
  $bold(h)_c in {0, 1}^n$, set the centered vector
  $bold(v)_c = bold(h)_c - (|C_c| / n) bold(1)$. By symmetry of $B_1$
  and $B_1 bold(1) = bold(0)$,
  $ bold(h)_c^top B_1 bold(h)_c = bold(v)_c^top B_1 bold(v)_c
    <= lambda_("max") (B_1) norm(bold(v)_c)_2^2. $
  Expanding the squared norm gives
  $norm(bold(v)_c)_2^2 = |C_c| - |C_c|^2 / n$. Summing over
  $c = 1, dots, K$ and applying Cauchy--Schwarz
  ($sum_(c = 1)^K |C_c|^2 >= (sum_c |C_c|)^2 / K = n^2 / K$) yields
  $ sum_(c = 1)^K norm(bold(v)_c)_2^2
    = n - 1 / n sum_(c = 1)^K |C_c|^2 <= n (1 - 1 / K). $
  Since $B_1 bold(1) = bold(0)$ ensures $lambda_("max") (B_1) >= 0$,
  dividing by $2 m_w$ establishes the bound.
]

#subsec("C", "NP-Hardness and the Modularity Resolution Limit")

#theorem(
  title: "NP-Completeness and Fortunato--Barthelemy Resolution Limit",
  origin: "Classical — Cited",
)[
  (i) Deciding whether a graph admits a partition with $Q_1 (cal(P)) >= tau$
  is NP-complete @brandes2008modularity. (ii) Merging two disjoint
  communities $C_1, C_2$ connected by undirected edge weight
  $e_(1 2) = sum_(i in C_1, j in C_2) A_(i j)$ strictly increases
  $Q_gamma$ iff $e_(1 2) > gamma Sigma_("tot") (C_1) Sigma_("tot") (C_2) / (2 m_w)$.
  In particular, on a ring of $q$ cliques $K_r$ ($r >= 3$) joined by
  single bridge edges, merging adjacent clique pairs strictly increases
  $Q_gamma$ whenever $gamma (r (r - 1) + 2) < q$ @fortunato2007resolution.
]

#proof[
  (i) Reduction from 3-Partition onto clique-star gadgets is proved in
  @brandes2008modularity. (ii) Merging $C_1$ and $C_2$ raises
  $Sigma_("in")$ by $2 e_(1 2)$ and increases the sum of squared
  strengths by $(Sigma_("tot") (C_1) + Sigma_("tot") (C_2))^2
  - Sigma_("tot") (C_1)^2 - Sigma_("tot") (C_2)^2
  = 2 Sigma_("tot") (C_1) Sigma_("tot") (C_2)$. Hence
  $ Delta Q_gamma = 1 / m_w
    [ e_(1 2) - gamma (Sigma_("tot") (C_1) Sigma_("tot") (C_2)) / (2 m_w) ], $
  which is strictly positive iff the stated inequality holds. On the
  clique ring, each clique has strength $d_C = r (r - 1) + 2$,
  $2 m_w = q d_C$, and $e_(1 2) = 1$. Thus $Delta Q_gamma > 0$ iff
  $1 > gamma d_C^2 / (q d_C) = gamma d_C / q$, i.e.,
  $gamma (r (r - 1) + 2) < q$.
]

#subsec("D", "Constructive Disconnected-Community Counterexample")

A fundamental structural defect of Louvain, identified empirically and
theoretically by Traag, Waltman, and van Eck @traag2019leiden, is that
bridge vertices may depart a community during Phase 1 after holding it
together, leaving permanently disconnected communities. We exhibit a
minimal explicit 26-vertex unweighted graph proving this pathology from
singleton initialization.

#proposition(
  title: "Constructive Counterexample: Disconnected Louvain Output",
  origin: "Original Counterexample",
)[
  There exists a 26-vertex unweighted graph $G^*$ ($m_w = 90$) on which
  Louvain at $gamma = 1$ initialized from singletons terminates
  permanently (across both Phase 1 and Phase 2) with an internally
  disconnected community $C = A union B$ having zero edges between $A$
  and $B$, and splitting $C$ into $A$ and $B$ strictly increases $Q_1$
  by $+256 / 16200 approx +0.0158$.
]

#proof[
  Construct $G^*$ from four cliques $A = {a_1, dots, a_4} tilde.op K_4$,
  $B = {b_1, dots, b_4} tilde.op K_4$, $D = {d_1, dots, d_7} tilde.op K_7$,
  $E tilde.op K_9$, and a bridge pair $U = {u_1, u_2} tilde.op K_2$.
  Add edges from both $u_1$ and $u_2$ to ${a_1, a_2}$, to ${b_1, b_2}$,
  and to ${d_1, dots, d_6}$, leaving $E$ disjoint. Edge count is
  $m_w = 6 + 6 + 1 + 21 + 4 + 4 + 12 + 36 = 90$, so $2 m_w = 180$.
  Degrees are $k(a_1) = k(a_2) = k(b_1) = k(b_2) = 5$,
  $k(a_3) = k(a_4) = k(b_3) = k(b_4) = 3$, $k(u_1) = k(u_2) = 11$,
  $k(d_1) = dots = k(d_6) = 8$, $k(d_7) = 6$, and $k(e) = 8$ for
  $e in E$; strengths sum to
  $Sigma_("tot") (A) = 16$, $Sigma_("tot") (B) = 16$,
  $Sigma_("tot") (U) = 22$, $Sigma_("tot") (D) = 54$,
  $Sigma_("tot") (E) = 72$ (total $180$). By Lemma II.1, the scaled
  gain $S(i : a -> b) = 2 m_w^2 Delta Q_1$ equals
  $180 (k_(i, b) - k_(i, a)) + k_i (Sigma_("tot")^((-i)) (a) - Sigma_("tot") (b))$.

  *Pass 1* (order $u_1, u_2, a_2, a_3, a_4, b_1, b_2, b_3, b_4$, then
  $d_7, dots, d_1$, then $E$):
  (1) $u_1$ sees singletons of smallest strength $5$ ($a_1, a_2, b_1, b_2$)
  and joins ${a_1}$ ($S = 180 - 55 = 125$), forming $C = {a_1, u_1}$
  ($Sigma_("tot") = 16$).
  (2) $u_2$ has $k_(u_2, C) = 2$, scoring $180(2) - 11(16) = 184$,
  strictly beating every singleton ($<= 125$), so $u_2 -> C$
  ($Sigma_("tot") = 27$).
  (3) $a_2, a_3, a_4$ sequentially join $C$ with scores $405, 264, 435$
  (all exceeding any singleton move $<= 171$), giving $C = A union U$
  ($Sigma_("tot") = 38$).
  (4) $b_1$ has $k_(b_1, C) = 2$, scoring $180(2) - 5(38) = 170$,
  strictly beating singleton ${b_3}$ ($180 - 15 = 165$), so $b_1 -> C$
  ($Sigma_("tot") = 43$). Then $b_2, b_3, b_4$ join $C$ with scores
  $325, 216, 387$, yielding $C = A union B union U$ ($Sigma_("tot") = 54$).
  (5) Next $d_7$ joins ${d_6}$ ($C_D$, $Sigma_("tot") = 14$); every
  $d_j$ ($j <= 5$) scores negative for $C$ ($360 - 8(54) = -72$) and
  joins $C_D$, yielding $C_D = D$ ($Sigma_("tot") = 54$), while $E$
  forms $C_E = E$ ($Sigma_("tot") = 72$).

  *Pass 2*: Revisit $u_1$ ($k_(u_1, C) = 5, Sigma^((-u_1)) (C) = 43$
  vs $k_(u_1, C_D) = 6, Sigma(C_D) = 54$):
  $S(u_1 : C -> C_D) = 180(1) + 11(43 - 54) = +59 > 0$, so $u_1$ moves
  to $C_D$. Revisit $u_2$ ($k_(u_2, C) = 4, Sigma^((-u_2)) (C) = 32$
  vs $k_(u_2, C_D) = 7, Sigma(C_D) = 65$):
  $S(u_2 : C -> C_D) = 180(3) + 11(32 - 65) = +177 > 0$, so $u_2$ moves
  to $C_D$. Now $C = A union B$ ($Sigma = 32$) is internally
  disconnected, and $C_D = D union U$ ($Sigma = 76$).

  *Pass 3 (Local optimality & Phase 2)*: For $a_1$ (and $a_2, b_1, b_2$),
  moving to $C_D$ yields $S = 180(2 - 3) + 5(27 - 76) = -425 < 0$, and
  moving to $emptyset$ yields $180(-3) + 5(27) = -405 < 0$. All other
  vertices in $C, C_D, C_E$ have all neighbors inside their current
  community or strictly negative exit scores ($u_1 -> C$ gives $-177$).
  Phase 1 halts. In Phase 2, $e(C, C_D) = 8 < (32 dot 76) / 180 = 13.51$,
  so by Theorem II.4 super-nodes $C$ and $C_D$ do not merge. Louvain
  terminates with disconnected $C = A union B$. Splitting $C$ into $A$
  and $B$ changes $Q_1$ by
  $+Sigma_("tot") (A) Sigma_("tot") (B) / (2 m_w^2) = 256 / 16200 > 0$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Using compressed sparse row (CSR) adjacency storage and a scratch array
of length $n$ whose touched community slots are cleared in
$O("deg"(i))$ time per vertex, each Phase-1 sweep evaluates all exact
gains in $O(m)$ arithmetic operations without hash-table overhead.
Phase-2 aggregation scans every edge once in $O(m)$ time, typically
reducing vertex count by a geometric factor $rho < 1$ across levels.
If $L$ denotes the total number of sweeps summed across all hierarchy
levels, total work is $O(L m)$; empirically $L = O(1)$ to $O(log n)$
on sparse real-world networks, scaling smoothly to $n > 10^8$ vertices
@blondel2008louvain. Memory footprint is strictly $O(n + m)$.

Geometrically, modularity optimization is non-metric and relies solely
on topological edge contrast against the configuration-model
expectation $gamma k_i k_j / (2 m_w)$. Because the objective landscape
of $Q_gamma$ exhibits an exponential degeneracy of near-optimal
partitions with similar scores, node visitation order induces moderate
stochastic variance across random shuffles.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time complexity], [$O(L m)$ empirical ($O(m)$ per sweep); $O(n^2)$ dense.],
    [Space complexity], [$O(n + m)$ CSR graph plus $O(n)$ community statistics.],
    [Objective function], [Non-convex discrete quadratic $Q_gamma (cal(P))$ (NP-hard).],
    [Hyperparameters], [Resolution $gamma > 0$ (default $1$); tolerance $epsilon$.],
    [Determinism], [Order-dependent; random permutations yield degenerate optima.],
    [Invariances], [Positive weight scale-invariant; vertex permutation-equivariant.],
    [Structural defects], [Resolution limit (Thm. II.4); disconnected clusters (Prop. II.5).],
    [Hierarchy output], [Native multilevel dendrogram across aggregation phases.],
  ),
  caption: [Operational and theoretical profile of the Louvain method.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Web-Scale Social and Citation Networks:* Social graphs and
   bibliometric citation networks exhibit power-law degree
   distributions and multi-scale homophily. The configuration-model
   null term $k_i k_j / (2 m_w)$ automatically discounts high-degree
   hubs while Phase-2 coarsening exposes nested macro-communities at
   near-linear cost @blondel2008louvain.
2. *Single-Cell Transcriptomic Graph Clustering:* In scRNA-seq
   pipelines, cells are embedded via PCA and linked into a shared
   nearest-neighbor ("SNN") graph. Louvain's topological objective
   recovers non-convex cellular manifolds and rare transcriptional
   lineages without requiring pre-specified cluster counts $K$.
3. *Functional Brain Connectomics (fMRI):* Thresholded resting-state
   correlation networks across cortical regions display small-world
   modular organization. Louvain extracts hierarchical resting-state
   subnetworks (default-mode, salience, motor) while remaining
   invariant to global correlation rescaling (Proposition II.3).
4. *Financial Fraud-Ring and Money-Laundering Detection:* Transaction
   graphs among accounts contain tightly knit collusive cycles embedded
   in sparse commerce flows. Weighted modularity flags anomalously
   dense subgraphs whose internal transfer volume exceeds random
   degree-matched expectations.
5. *Software Architecture and Microservice Decomposition:* Extracting
   cohesive modules from static call graphs and dependency networks
   maps directly onto maximizing internal call density minus random
   coupling, using Louvain's multilevel hierarchy to guide package and
   service boundaries.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Louvain fails systematically in four regimes: (a) *heterogeneous-scale
networks* where the Fortunato--Barthélemy resolution limit merges
small ground-truth cliques below $sqrt(2 m_w / gamma)$ unless $gamma$
is tuned; (b) *sparse bridge-rich topologies* where Phase-1 bridge
migration creates internally disconnected communities (Proposition II.5),
motivating the Leiden refinement fix (Alg. 108); (c) *disassortative or
bipartite graphs* (e.g., ecological food webs or bipartite recommender
graphs) where assortative modularity misaligns with cross-group block
structure, favoring Stochastic Block Models @karrer2011dcsbm
@abbe2018sbmreview; and (d) *high-dimensional Euclidean data with
severe hubness*, where naive k-NN graphs distort degree strengths.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Social / citation graphs], [Strong: hub-corrected null model + $O(L m)$ scale], [Use Leiden (Alg. 108) for connectivity],
    [scRNA-seq k-NN graphs], [High: non-convex manifolds, unknown $K$], [Disconnected clusters; prefer Leiden / CPM],
    [Brain connectomes], [High: multi-scale modular hierarchy], [Sweep $gamma$ for resolution limit],
    [Disassortative / bipartite], [Poor: assumes assortative density], [Use DC-SBM @karrer2011dcsbm or Spectral],
    [Streaming / dynamic graphs], [Moderate: batch coarsening needed], [Incremental LPA (Alg. 110) or streaming],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for Louvain.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Louvain stands at the nexus of greedy agglomerative modularity
@clauset2004cnm, spectral graph partitioning @shi2000ncut
@vonluxburg2007tutorial, and modern refinement-based community
detection @traag2019leiden. The following identity makes its
exact geometric bridge to minimum-cut clustering transparent.

#proposition(
  title: "Exact Cut-Balance Equivalence on Regular Graphs",
  origin: "Original Derivation",
)[
  On any $d$-regular weighted graph ($k_i = d$ for all $i in V$, so
  $2 m_w = n d$), modularity maximization is algebraically identical to
  minimizing a quadratically size-regularized graph cut:
  $ Q_gamma (cal(P)) = 1 - 1 / (n d) sum_(c = 1)^K
    "cut"(C_c, overline(C)_c) - gamma / n^2 sum_(c = 1)^K |C_c|^2, $
  where $"cut"(C_c, overline(C)_c) = sum_(i in C_c, j in.not C_c) A_(i j)$.
]

#proof[
  For every community $C_c$, conservation of incident strength gives
  $Sigma_("in") (c) = sum_(i in C_c) k_i - "cut"(C_c, overline(C)_c)
  = d |C_c| - "cut"(C_c, overline(C)_c)$ and $Sigma_("tot") (c) = d |C_c|$.
  Substituting $2 m_w = n d$ into Definition I.1 yields
  $ Q_gamma (cal(P)) &= sum_(c = 1)^K
    (d |C_c| - "cut"(C_c, overline(C)_c)) / (n d)
    - gamma sum_(c = 1)^K ((d |C_c|) / (n d))^2 \
    &= 1 - 1 / (n d) sum_(c = 1)^K "cut"(C_c, overline(C)_c)
    - gamma / n^2 sum_(c = 1)^K |C_c|^2, $
  using $sum_(c = 1)^K |C_c| = n$.
]

Thus on regular k-NN graphs, Louvain trades off cut weight against the
Simpson concentration $sum_c (|C_c| / n)^2$, closely paralleling
Normalized Cut @shi2000ncut while bypassing $O(n^3)$ eigendecomposition;
it directly accelerates Clauset--Newman--Moore @clauset2004cnm from
$O(m d_("tree") log n)$ to $O(L m)$, generalizes Label Propagation
(Alg. 110, which corresponds to the $gamma -> 0^+$ limit of Phase 1),
and is strictly improved in connectivity guarantees by the Leiden
algorithm (Alg. 108, @traag2019leiden) and flow-compression methods
such as Infomap @rosvall2008infomap and MCL (Alg. 109).

*Empirical Behavior (ClusterOverview 4-Dataset Benchmark).* Evaluated
on k-NN affinity graphs across four standardized datasets using the
Silhouette coefficient @rousseeuw1987silhouettes and Adjusted Rand
Index @hubert1985comparing, Louvain records *Iris* Sil = $0.354$,
ARI = $0.333$; *Customers* Sil = $0.430$; *Multi-Shape* Sil = $0.526$
(ARI = N/A); and *High-D* Sil = $0.232$ (ARI = N/A), with the fastest
mean runtime in this graph family at *$0.437$ ms*. This empirical
profile mirrors the theory directly: on *Multi-Shape* ($0.526$) and
*Customers* ($0.430$), Louvain's topological cut-balance objective
captures non-convex density-connected manifolds in sub-millisecond
time. On *Iris*, default $gamma = 1$ splits or merges across the
overlapping Versicolor--Virginica k-NN bridge (limiting ARI to $0.333$
despite respectable geometric separation Sil = $0.354$). Most tellingly,
on *High-D*, Louvain drops to Sil = $0.232$ — substantially below
Leiden ($0.303$), LPA ($0.303$), and MCL ($0.321$) — reflecting both
degree-hub distortion of the configuration null model and Phase-1
bridge-migration fragmentation (Proposition II.5) in hub-dominated
high-dimensional k-NN graphs.
