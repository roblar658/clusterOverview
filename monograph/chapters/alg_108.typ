#import "../template/ieee.typ": *

#alg_entry(
  id: 108,
  name: "Leiden Algorithm",
  category: "Cat. V: Graph, Spectral & Network",
  year: "2019",
  paradigm: "Refinement-Guaranteed Hierarchical Community Detection",
  complexity: $O(L m)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

The *Leiden algorithm* @traag2019leiden resolves a foundational
structural flaw of the Louvain method (Alg. 107, @blondel2008louvain):
greedy local modularity ascent can output arbitrarily poorly connected
or even internally disconnected communities (up to $25\%$ of clusters
on real networks). Leiden introduces an intermediate *partition
refinement* phase between local moving and graph coarsening, paired
with an active-node queue that accelerates convergence while unifying
Modularity and the resolution-limit-free Constant Potts Model (CPM)
@reichardt2006spinglass.

Let $G = (V, E, w)$ be an undirected weighted graph with symmetric
adjacency $A in RR_(>= 0)^(n times n)$, vertex degrees
$k_i = sum_j A_(i j)$, total weight $2 m_w = sum_i k_i$, and positive
node weights $w_i > 0$. For a community $C subset.eq V$, write
$Sigma_("in") (C) = sum_(i, j in C) A_(i j)$,
$e(C) = Sigma_("in") (C) / 2$, $W_C = sum_(i in C) w_i$, and
$e(S, T) = sum_(i in S, j in T) A_(i j)$ for disjoint $S, T$. Both
standard quality objectives unify into a single quadratic functional:
$ H_gamma (cal(P)) = 1 / (2 Z) sum_(c = 1)^K
  (Sigma_("in") (C_c) - gamma / nu W_(C_c)^2) + "const", $
where $(Z, nu, w_i, "const") = (m_w, 2 m_w, k_i, 0)$ yields
resolution-parameterized Modularity $Q_gamma$, while
$(Z, nu, w_i, "const") = (1, 1, 1, gamma n / 2)$ yields the
Constant Potts Model
$H_("CPM") (cal(P)) = sum_(c = 1)^K (e(C_c) - gamma binom(|C_c|, 2))$
(using $sum_c |C_c| = n$).

#definition(title: "Gamma-Separated, Gamma-Connected, and Gamma-Dense Partitions")[
  Relative to resolution $gamma > 0$ and normalization $nu > 0$, a
  partition $cal(P) = {C_1, dots, C_K}$ of $G$ is:
  (i) *$gamma$-separated* if every distinct pair satisfies
  $e(C_a, C_b) <= (gamma / nu) W_(C_a) W_(C_b)$;
  (ii) *$gamma$-connected* (well-connected) if every community $C in cal(P)$
  and every nontrivial subset $emptyset != S subset C$ satisfy
  $e(S, C without S) >= (gamma / nu) W_S W_(C without S)$;
  (iii) *subset optimal* if no nonempty subset $S subset.eq C$ can be
  reassigned to another community (or $emptyset$) to strictly increase
  $H_gamma$;
  (iv) *$gamma$-dense* (under CPM) if every $C in cal(P)$ satisfies
  $e(C) >= gamma binom(|C|, 2)$.
]

Crucially, Leiden coarsens the graph using a *refined* sub-partition
$cal(P)_("ref")$ rather than $cal(P)$ itself, while initializing each
super-node's cluster assignment to its parent community in $cal(P)$.
This decoupling gives the next level freedom to split poorly connected
parent communities without losing multilevel acceleration.

#algbox(
  title: "Leiden Algorithm (Fast Local Move + Refinement + Aggregation)",
  inputs: [Graph $G = (V, E, w)$, quality $H_gamma$ ($gamma > 0$), randomness $theta > 0$.],
  outputs: [Connected community partition $cal(P)^*$ of $V$.],
)[
  1. Initialize $G^((0)) <- G$, singleton partition $cal(P) <- {{i} : i in V}$, $h <- 0$.
  2. *Phase 1 (Fast Local Move)*: push all $v in V^((h))$ onto FIFO
     queue $Q$ (random order). While $Q != emptyset$, pop $v$; find
     $c^* = arg max_(c in cal(N)_C (v) union {emptyset}) Delta H_gamma (v : c_v -> c)$.
     If $Delta H_gamma (v : c_v -> c^*) > 0$, move $v$ to $c^*$ and
     enqueue every $u in N(v)$ with $c_u != c^*$ and $u in.not Q$.
  3. *Termination check*: if $|cal(P)| = |V^((h))|$, stop and return
     the lifted partition $cal(P)^*$ on $V$.
  4. *Phase 2 (Partition Refinement)*: set $cal(P)_("ref") <- {{v} : v in V^((h))}$.
     For each $C in cal(P)$, compute well-connected vertices
     $R(C) = {v in C : e(v, C without {v}) >= (gamma / nu) w_v (W_C - w_v)}$.
     For each $v in R(C)$ still a singleton in $cal(P)_("ref")$,
     collect well-connected target sub-communities
     $cal(T)_v = {S in cal(P)_("ref") : S subset C, e(S, C without S) >= (gamma / nu) W_S (W_C - W_S), Delta H_gamma (v -> S) >= 0}$.
     Sample $S^* in cal(T)_v$ with probability
     $"Pr"(S^*) prop exp(Delta H_gamma (v -> S^*) / theta)$ and merge
     $v$ into $S^*$.
  5. *Phase 3 (Refined Aggregation)*: coarsen $G^((h))$ along
     $cal(P)_("ref")$ to obtain $G^((h + 1))$ (summing edge and node
     weights $w'_S = W_S$); initialize each super-node $S$ in the
     parent community $C supset.eq S$ from $cal(P)$. Set $h <- h + 1$
     and go to Step 2.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish Leiden's internal-connectivity invariant, asymptotic
subset-optimality hierarchy, exact resolution-limit immunity of CPM,
and an original residual-gain bound characterizing the active-node
queue in Fast Local Move.

#subsec("A", "Guaranteed Internal Connectivity Invariant")

#theorem(
  title: "Guaranteed Internal Connectivity of Leiden Communities",
  origin: "Original Proof",
)[
  Let $G$ have nonnegative edge weights, strictly positive node weights
  $w_i > 0$, and resolution $gamma > 0$. Then:
  (i) Every refined sub-community $S in cal(P)_("ref")$ produced by
  Phase 2 induces a connected subgraph $G[S]$.
  (ii) Whenever a Leiden recursion terminates (Step 3), every output
  community $C in cal(P)^*$ induces a connected subgraph $G[C]$.
]

#proof[
  (i) We proceed by induction on the number of merges in Phase 2 at
  level $h$ (and inductively on level $h$, assuming every super-node
  of $G^((h))$ represents a connected vertex subset of $G$; base case
  $h = 0$ holds trivially for singletons). Initially every
  $S in cal(P)_("ref")$ is a single super-node, hence connected in $G$.
  Suppose singleton $v in R(C)$ merges into $S subset C$ with
  $G[S]$ connected. By construction of $cal(T)_v$, the merge satisfies
  $ Delta H_gamma (v -> S)
    = 1 / Z (e(v, S) - gamma / nu w_v W_S) >= 0. $
  Because $gamma > 0$, $w_v > 0$, and $W_S > 0$, this forces
  $ e(v, S) >= gamma / nu w_v W_S > 0. $
  Hence $v$ shares at least one positive-weight edge with $S$, so the
  lifted subgraph $G[S union {v}]$ is connected.

  (ii) Suppose for contradiction that at termination (Step 3 at final
  level $h^*$), some community $C in cal(P)^*$ is disconnected in $G$,
  partitioning as $C = C_1 union C_2$ with $e(C_1, C_2) = 0$ and
  $C_1, C_2 != emptyset$. Because Step 3 halts only when
  $|cal(P)| = |V^((h^*))|$, every output community $C$ corresponds to
  a single super-node of $G^((h^*))$, which by Step 5 was produced as
  one refined sub-community $S in cal(P)_("ref")$ at level $h^* - 1$.
  By part (i), $G[S]$ is connected, contradicting $e(C_1, C_2) = 0$.
  Moreover, if a disconnected community $C = C_1 union C_2$ ever arises
  after Phase 1 at any level $h$, Phase 2 keeps $C_1$ and $C_2$ in
  separate refined super-nodes (since $e(v, S) = 0 < (gamma / nu) w_v W_S$
  prevents cross-component merges). In $G^((h + 1))$, moving any
  disconnected super-node $S_1 subset.eq C_1$ from $C$ to $emptyset$
  yields
  $ Delta H_gamma (S_1 : C -> emptyset)
    = 1 / Z (-e(S_1, C without S_1)
    + gamma / nu W_(S_1) W_(C without S_1)) > 0 $
  whenever $S_1 = C_1$ (where $e(S_1, C without S_1) = 0$), strictly
  forcing Phase 1 to split disconnected components rather than stall.
]

#subsec("B", "Asymptotic Subset Optimality and Structural Bounds")

#theorem(
  title: "Asymptotic Subset Optimality and Gamma-Density",
  origin: "Classical — Cited",
)[
  Suppose Leiden is iterated repeatedly (re-initializing Step 1 from
  the previous iteration's output partition) with finite temperature
  $theta > 0$ @traag2019leiden. Then:
  (i) After every finite iteration, every community is connected and
  every pair is $gamma$-separated when no single-super-node merge
  improves $H_gamma$.
  (ii) With probability $1$ as iterations $t -> infinity$, the partition
  converges to a *subset-optimal* partition $cal(P)^*$, which is
  simultaneously $gamma$-connected, $gamma$-separated, and (under CPM)
  $gamma$-dense: $e(C) >= gamma binom(|C|, 2)$ for all $C in cal(P)^*$.
]

#proof[
  Because every accepted Phase-1 move strictly increases $H_gamma$ and
  Phase 3 preserves $H_gamma$ identically on refined partitions (by the
  exact block-indicator factorization of Theorem 107.II.2, with
  $W'_S = W_S$ invariant), $H_gamma (cal(P)^((t)))$ is non-decreasing
  over a finite partition lattice and stabilizes after finitely many
  iterations. Under $theta > 0$, every admissible merge sequence in
  Phase 2 has strictly positive probability, so any subset
  $S subset C$ whose detachment or transfer to another community $C'$
  strictly increases $H_gamma$ is eventually isolated as a refined
  super-node and moved in Phase 1 @traag2019leiden. At subset
  optimality, detaching any $emptyset != S subset C$ cannot improve
  $H_gamma$:
  $ Delta H_gamma (S : C -> emptyset)
    = 1 / Z (-e(S, C without S) + gamma / nu W_S W_(C without S)) <= 0, $
  proving $gamma$-connectedness. Merging any two communities $C_1, C_2$
  gives $Delta H_gamma = (1 / Z)(e(C_1, C_2) - (gamma / nu) W_(C_1) W_(C_2)) <= 0$,
  proving $gamma$-separation. Finally, summing $gamma$-connectedness
  over singletons $S = {v}$ for $v in C$ gives
  $ e(C) = 1 / 2 sum_(v in C) e(v, C without {v})
    >= gamma / (2 nu) (W_C^2 - sum_(v in C) w_v^2). $
  Under CPM ($w_v = 1, nu = 1$), $(|C|^2 - |C|) / 2 = binom(|C|, 2)$,
  yielding $e(C) >= gamma binom(|C|, 2)$.
]

#subsec("C", "Exact Subgraph Locality and Resolution-Limit Immunity of CPM")

Unlike modularity, whose global denominator $2 m_w$ couples distant
regions of the graph (Theorem 107.II.4), Leiden paired with CPM
exhibits strict subgraph locality.

#theorem(
  title: "Subgraph Locality and Resolution-Limit Freedom of CPM",
  origin: "Original Proof",
)[
  Let $cal(P)^*$ be a globally optimal CPM partition of $G$ at
  resolution $gamma > 0$, and let $U = union.big_(c in I) C_c^*$ be
  any union of communities from $cal(P)^*$. Then the restricted
  partition $cal(P)^*|_U = {C_c^* : c in I}$ is globally optimal for
  CPM on the induced subgraph $G[U]$ at the identical resolution
  $gamma$. Consequently, CPM is strictly immune to the
  Fortunato--Barthelemy resolution limit @fortunato2007resolution
  @traag2019leiden.
]

#proof[
  Suppose to the contrary that there exists a partition $cal(Q)$ of $U$
  with $H_("CPM")^(G[U]) (cal(Q)) > H_("CPM")^(G[U]) (cal(P)^*|_U)$.
  Form the hybrid partition $cal(P)' = cal(Q) union {C_c^* : c in.not I}$
  of $V$. Because $H_("CPM")$ sums community terms
  $e(C) - gamma binom(|C|, 2)$ that depend solely on internal edges
  and internal cardinalities — with no global graph parameter such as
  $m_w$ — every community outside $U$ contributes identically to
  $H_("CPM")^G (cal(P)')$ and $H_("CPM")^G (cal(P)^*)$, while edges
  between $U$ and $V without U$ are inter-community under both
  partitions and contribute zero. Therefore
  $ H_("CPM")^G (cal(P)') - H_("CPM")^G (cal(P)^*)
    = H_("CPM")^(G[U]) (cal(Q)) - H_("CPM")^(G[U]) (cal(P)^*|_U) > 0, $
  contradicting global optimality of $cal(P)^*$. In particular,
  embedding a clique ring inside an arbitrarily large disconnected or
  weakly linked network leaves the optimal CPM partition of the ring
  completely invariant.
]

#subsec("D", "Residual-Gain Bound for Queue-Based Fast Local Move")

Louvain repeatedly sweeps all $n$ vertices even when most neighborhoods
are static. Leiden's Fast Local Move enqueues only neighbors
$u in N(v) without C_(c^*)$ upon moving $v -> c^*$. Because moving $v$
also shifts community weights $W_(c_v)$ and $W_(c^*)$, non-adjacent
vertices experience a small null-model perturbation; we prove an exact
quantitative bound on any residual single-node gain when $Q$ empties.

#proposition(
  title: "Exact Residual-Gain Bound of the Active-Node Queue",
  origin: "Original Proof",
)[
  Suppose Fast Local Move terminates at step $T$ with empty queue $Q$.
  For any vertex $u in V$ in final community $a = c_u^((T))$, let
  $t_0 <= T$ be the last step at which $u$ was popped from $Q$ (so
  $c_u^((t_0)) = a$). Then for every alternative community $b != a$,
  $ Delta H_gamma^((T)) (u : a -> b)
    <= (gamma w_u) / (Z nu)
    [(W_a^((T)) - W_a^((t_0)))_+ + (W_b^((t_0)) - W_b^((T)))_+], $
  where $(x)_+ = max(0, x)$. Moreover, edge-count terms satisfy
  $e^((T)) (u, b) - e^((T)) (u, a) <= e^((t_0)) (u, b) - e^((t_0)) (u, a)$
  identically.
]

#proof[
  At step $t_0$, after $u$ made its greedy choice (staying in $a$ or
  moving into $a$), $u$ had no strictly improving move, so for every
  $b != a$:
  $ e^((t_0)) (u, b) - e^((t_0)) (u, a)
    + gamma / nu w_u (W_a^((-u), (t_0)) - W_b^((t_0))) <= 0. $
  Between $t_0$ and $T$, $u$ was never re-enqueued. By Step 2 of
  Leiden, if any neighbor $j in N(u)$ moved into a community $c != a$,
  then $u$ (having $c_u = a != c$) would have been enqueued — a
  contradiction. Thus no neighbor of $u$ left $a$ (so
  $e^((T)) (u, a) >= e^((t_0)) (u, a)$) and no neighbor of $u$ entered
  $b != a$ (so $e^((T)) (u, b) <= e^((t_0)) (u, b)$). Subtracting the
  step-$t_0$ optimality inequality from $Z Delta H_gamma^((T)) (u : a -> b)$
  eliminates the non-positive edge-difference increment and yields
  $ Z Delta H_gamma^((T)) (u : a -> b)
    &<= gamma / nu w_u
    [(W_a^((T)) - W_a^((t_0))) - (W_b^((T)) - W_b^((t_0)))] \
    &<= gamma / nu w_u
    [(W_a^((T)) - W_a^((t_0)))_+ + (W_b^((t_0)) - W_b^((T)))_+]. $
  Dividing by $Z$ completes the proof. On sparse graphs with
  $nu = 2 m_w$, the residual drift is $O(w_u Delta W / m_w^2)$ and is
  eliminated completely on subsequent Leiden iterations or a single
  $O(m)$ verification sweep.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Fast Local Move restricts attention to active boundary vertices, cutting
node visits by $2 times$ to $5 times$ over Louvain on large graphs once
dense cores stabilize @traag2019leiden. Phase-2 refinement touches each
vertex $v$ at most once when checking $e(v, C without {v})$ and once
when evaluating candidate sub-communities $S subset C$, visiting each
incident edge $O(1)$ times for $O(m)$ work per level. Phase-3
aggregation also runs in $O(m)$ time, giving overall complexity
$O(L m)$ and $O(n + m)$ space in CSR representation. On very small
graphs ($n < 10^3$), the extra refinement bookkeeping adds a small
constant factor relative to Louvain, whereas on massive networks
($n > 10^5$) Fast Local Move makes Leiden strictly faster in wall time.

Geometrically, under CPM at resolution $gamma$, every community has
internal edge density at least $gamma$ while inter-community density is
at most $gamma$ (Theorem II.2), carving the graph into cohesive density
level sets rather than degree-normalized ratios.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time complexity], [$O(L m)$ ($O(m)$ per refinement and aggregation level).],
    [Space complexity], [$O(n + m)$ CSR storage plus $O(n)$ queue / weight arrays.],
    [Supported objectives], [Modularity $Q_gamma$ and Constant Potts Model $H_("CPM")$.],
    [Connectivity guarantee], [100% internally connected communities at convergence (Thm. II.1).],
    [Asymptotic guarantee], [Subset-optimal, $gamma$-connected, $gamma$-separated (Thm. II.2).],
    [Hyperparameters], [Resolution $gamma > 0$, Gibbs temperature $theta > 0$ (default $0.01$).],
    [Resolution limit], [Present under $Q_gamma$; strictly absent under CPM (Thm. II.3).],
    [Determinism], [Stochastic (random queue order + Gibbs sampling in refinement).],
  ),
  caption: [Operational and theoretical profile of the Leiden algorithm.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Single-Cell RNA-seq and Spatial Transcriptomics:* Leiden is the
   modern default in Scanpy and Seurat, superseding Louvain because
   disconnected Louvain clusters erroneously merge biologically distinct
   rare cell states across k-NN manifold bridges. Theorem II.1 enforces
   connected transcriptional states, while CPM prevents resolution-limit
   loss of rare immune or stem subpopulations.
2. *Large-Scale Scientometric and Citation Mapping:* Developed at CWTS
   Leiden specifically for tens-of-millions-paper citation graphs
   @traag2019leiden, Leiden with CPM resolves fine-grained research
   specialties consistently across sub-disciplines without global
   network-size distortion (Theorem II.3).
3. *Protein--Protein Interaction (PPI) Complex Discovery:* Functional
   protein complexes are dense, internally connected subgraphs whose
   biological validity requires physical path connectivity. Leiden's
   $gamma$-density and $gamma$-connectedness guarantees prevent spurious
   fusion of disjoint complexes linked through promiscuous chaperone hubs.
4. *Urban Mobility and Commuting-Zone Delineation:* Spatial origin--
   destination flow networks require contiguous geographic regions.
   Refinement eliminates spatially fragmented enclaves that arise when
   high-traffic transit hubs switch communities during Louvain sweeps.
5. *Hierarchical GraphRAG and Knowledge-Graph Summarization:* Modern
   retrieval-augmented generation builds multilevel community summaries
   over entity-relation graphs; Leiden's connected multilevel hierarchy
   ensures every entity summary corresponds to a coherent relational
   component rather than unrelated fragments.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Despite dominating Louvain structurally, Leiden retains four operational
caveats: (a) running Leiden with Modularity $Q_gamma$ still inherits
the Fortunato--Barthélemy resolution limit (Theorem 107.II.4);
(b) running Leiden with CPM avoids the resolution limit but requires
tuning $gamma$ carefully on graphs with strongly heterogeneous cluster
densities, since a single global density threshold $gamma$ can
over-fragment sparse communities or merge dense cores; (c) bipartite or
disassortative block structures violate the assortative objective; and
(d) Gibbs sampling ($theta > 0$) introduces run-to-run variability
requiring consensus clustering @monti2003consensus for critical
scientific reporting.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [scRNA-seq / spatial omics], [Ideal: connected rare-cell clusters via CPM], [Tune $gamma$ across density scales],
    [Web / citation graphs], [Ideal: $O(L m)$ + subgraph locality (Thm. II.3)], [Fix random seed or consensus-aggregate],
    [PPI functional modules], [High: $gamma$-dense connected complexes], [Overlapping complexes: OSLOM / BIGCLAM],
    [Heterogeneous-density nets], [Moderate under single-$gamma$ CPM], [Multi-scale sweep or Infomap @rosvall2008infomap],
    [Bipartite / core-periphery], [Poor: assortative objective mismatch], [Use DC-SBM @karrer2011dcsbm or Spectral],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for Leiden.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Leiden strictly subsumes Louvain (Alg. 107) in structural quality while
connecting directly to Potts spin-glass clustering
@reichardt2006spinglass, Infomap @rosvall2008infomap, and spectral cut
relaxations @shi2000ncut @vonluxburg2007tutorial. We prove rigorously
that Leiden resolves Louvain's exact counterexample from Alg. 107.

#proposition(
  title: "Strict Dominance of Leiden Over Louvain on Disconnected Graphs",
  origin: "Original Proof",
)[
  On the 26-vertex counterexample graph $G^*$ of Proposition 107.II.5
  at $gamma = 1$, initialized from Louvain's converged disconnected
  partition $cal(P)_("Louv") = {A union B, D union U, E}$, one Leiden
  refinement and aggregation step strictly increases modularity by
  $Delta Q_1 = +256 / 16200 approx +0.0158$ and terminates at the
  internally connected partition $cal(P)^* = {A, B, D union U, E}$.
]

#proof[
  In $cal(P)_("Louv")$, the disconnected community $C = A union B$ has
  $Sigma_("tot") (C) = 32$, $e(A, B) = 0$, and $2 m_w = 180$. In
  Leiden's Phase 2 refinement inside $C$, every vertex $a in A$ has
  internal degree $e(a, C without {a}) = 3$. Check the well-connected
  condition $3 >= (k_a (32 - k_a)) / 180$: for $k_a = 5$,
  $5 dot 27 / 180 = 0.75 <= 3$; for $k_a = 3$, $3 dot 29 / 180 = 0.48 <= 3$.
  Hence $R(C) = A union B$. Because $e(a, S) = 0 < (k_a W_S) / 180$
  for any $a in A$ and $S subset.eq B$ (and symmetrically for $b in B$
  and $S subset.eq A$), no cross-clique merge has $Delta Q_1 >= 0$.
  Meanwhile every intra-clique merge inside $A$ (or $B$) has
  $e(v, S) >= 1 > 5 dot 16 / 180 = 0.44$, so Phase 2 produces two
  separate connected refined super-nodes $S_A = A$ and $S_B = B$
  (alongside connected super-nodes for $D union U$ and $E$).

  In Phase 3, $S_A$ and $S_B$ are initialized inside parent community
  $C$ with $e(S_A, S_B) = 0$ and $W_(S_A) = W_(S_B) = 16$. During
  Phase 1 on $G'$, moving $S_A$ from $C$ to $emptyset$ yields
  $ Delta Q_1 (S_A : C -> emptyset)
    = (16 dot 16) / (2 dot 90^2) = 256 / 16200 > 0. $
  Moving $S_A$ from $C$ into $C_D = D union U$ ($e(S_A, C_D) = 4, W_(C_D) = 76$)
  gives $Delta Q_1 = (1 / 90)(4 + (16 dot (16 - 76)) / 180) = -240 / 16200 < 0$.
  Hence $S_A$ moves strictly to $emptyset$, separating $A$ and $B$ into
  distinct connected communities with $Q_1 (cal(P)^*) = Q_1 (cal(P)_("Louv")) + 256 / 16200$.
]

*Empirical Behavior (ClusterOverview 4-Dataset Benchmark).* Across the
four ClusterOverview benchmarks (Silhouette @rousseeuw1987silhouettes,
ARI @hubert1985comparing), Leiden achieves *Iris* Sil = $0.341$,
ARI = $0.340$; *Customers* Sil = $0.430$; *Multi-Shape* Sil = $0.526$
(ARI = N/A); and *High-D* Sil = $0.303$ (ARI = N/A), at a mean runtime
of *$0.693$ ms*. Comparing directly against Louvain (Alg. 107) validates
Theorem II.1 and Proposition V.1 empirically: on *Customers* ($0.430$)
and *Multi-Shape* ($0.526$), well-behaved low-dimensional k-NN graphs
already yield connected Louvain partitions so both methods match
identically, whereas on *High-D* — where high-dimensional k-NN hubness
triggers bridge-migration disconnects — Leiden surges from Louvain's
$0.232$ to *$0.303$* ($+30.6\%$ relative Silhouette gain) while
improving *Iris* ARI ($0.340$ vs $0.333$). The modest $0.256$ ms
runtime overhead over Louvain ($0.693$ ms vs $0.437$ ms) reflects
Phase-2 refinement bookkeeping on small benchmark graphs, an overhead
that reverses in Leiden's favor on large-scale graphs via Fast Local
Move.
