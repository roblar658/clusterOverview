#import "../template/ieee.typ": *

#alg_entry(
  id: 129,
  name: "Fast Unfolding (Blondel)",
  category: "Cat. V: Graph, Spectral & Network",
  year: "2008",
  paradigm: "Multilevel Greedy Modularity Maximization",
  complexity: $cal(O)(L T m)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Early modularity maximization heuristics such as greedy agglomerative merging
@clauset2004cnm and spectral bisection @newman2006leading either scale poorly
on million-node graphs or become trapped in shallow local optima caused by
irreversible pairwise merges. Fast Unfolding, introduced by Blondel, Guillaume,
Lambiotte, and Lefebvre @blondel2008louvain—universally known today as the
*Louvain method*—resolves both bottlenecks through a two-phase multilevel
hierarchy alternating local single-vertex moves with exact super-graph
aggregation.

Let $G = (V, E, w)$ be a weighted undirected graph with symmetric non-negative
adjacency matrix $A in bb(R)_(>= 0)^(n times n)$, permitting self-loops
$A_(i i) >= 0$ that arise naturally at coarse hierarchy levels. Define the
weighted vertex degree $k_i = sum_(j=1)^n A_(i j)$ and total matrix weight
$2 m = sum_(i=1)^n k_i = sum_(i, j = 1)^n A_(i j)$.

#definition(title: "Newman--Girvan Modularity and Multilevel Super-Graph")[
  For a hard partition $cal(P) = {C_1, dots, C_r}$ of $V$ encoded by community
  labels $c : V -> {1, dots, r}$, the *modularity* @newman2004finding
  @blondel2008louvain is
  $ Q(cal(P); G) &= 1 / (2 m) sum_(i, j = 1)^n (A_(i j) - (k_i k_j) / (2 m)) delta(c_i, c_j) \
                 &= sum_(c=1)^r [ (Sigma_"in" (c)) / (2 m) - ((Sigma_"tot" (c)) / (2 m))^2 ], $
  where $Sigma_"in" (c) = sum_(i, j in C_c) A_(i j)$ is the sum of internal
  matrix entries (twice internal undirected edge weights plus self-loops) and
  $Sigma_"tot" (c) = sum_(i in C_c) k_i$ is the incident degree volume of
  community $c$. The *aggregated super-graph* $G' = (cal(P), E', A')$ induced
  by $cal(P)$ has one super-node per community $c in {1, dots, r}$ and
  symmetric weights
  $ A'_(c d) = sum_(i in C_c) sum_(j in C_d) A_(i j), quad A'_(c c) = Sigma_"in" (c). $
]

Because $Q(cal(P); G) in [-1 / 2, 1)$ compares observed intra-community edge
weight against the configuration-model null expectation $k_i k_j / (2 m)$,
maximizing $Q$ seeks dense, assortative modules without fixing the number of
communities $r$ a priori. Crucially, collapsing each Phase-1 community into a
weighted self-looped super-node enables subsequent Phase-1 sweeps to execute
collective multi-vertex moves while keeping the objective function exact.

#algbox(
  title: "Fast Unfolding (Blondel Multilevel Modularity Algorithm)",
  inputs: [Weighted graph $G = (V, E, A)$, gain tolerance $tau_"min" >= 0$],
  outputs: [Multilevel hierarchy ${cal(P)^((0)), dots, cal(P)^((L))}$ and flat partition],
)[
  + Initialize level $ell <- 0$, $G^((0)) <- G$, and singleton assignment
    $c_i <- i$ for all $i in V^((0))$.
  + *Phase 1 (Greedy Local Moving):* Maintain cached community volumes
    $Sigma_"tot" (c)$ for all active $c$. Repeatedly sweep through
    $i in V^((ell))$ (in fixed or shuffled order):
    - Temporarily detach $i$ from $a = c_i$ and compute incident weights
      $k_(i, c) = sum_(j in C_c without {i}) A^((ell))_(i j)$ to each
      neighboring community $c$.
    - Evaluate $Delta Q(i : a -> b)$ via Lemma 1 for all candidate
      communities $b$ adjacent to $i$ (plus $a$).
    - Move $i$ to $b^* = arg max_b Delta Q(i : a -> b)$ iff
      $Delta Q(i : a -> b^*) > tau_"min"$, updating $Sigma_"tot" (a)$ and
      $Sigma_"tot" (b^*)$ in $cal(O)(1)$ time.
    - Stop Phase 1 when a full sweep yields no move.
  + *Termination Check:* If Phase 1 performed zero moves on $G^((ell))$,
    terminate and return the hierarchy.
  + *Phase 2 (Network Aggregation):* Build super-graph $G^((ell+1))$ with
    super-weights $A^((ell+1))_(c d) = sum_(i in C_c, j in C_d) A^((ell))_(i j)$,
    set $ell <- ell + 1$, reset $c_u <- u$ on $V^((ell))$, and go to Step 2.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We present four rigorous results: the exact self-loop-aware $cal(O)(1)$
modularity gain formula, exact multilevel objective invariance and finite
convergence, the classical resolution limit and NP-completeness barrier, and a
fully verified unweighted counterexample proving that Fast Unfolding can
produce internally disconnected communities at a strict Phase-1 local optimum.

#subsec("A", "Exact O(1) Modularity Gain Formula With Self-Loops")

At hierarchy levels $ell >= 1$, super-nodes carry large self-loop weights
$A_(i i) = Sigma_"in" (i) > 0$. Deriving the exact move gain directly from the
matrix definition ensures zero discrepancy across all levels.

#lemma(
  title: "Exact Modularity Change for Single-Vertex Reassignment",
  origin: "Original Derivation",
)[
  Let $i in V$ have weighted degree $k_i = sum_j A_(i j)$ and self-loop weight
  $A_(i i) >= 0$, and let $k_(i, c) = sum_(j in C_c without {i}) A_(i j)$
  denote the weight of non-loop edges linking $i$ to $C_c without {i}$.
  + Inserting an isolated singleton ${i}$ into a community $C_c$ ($i in.not C_c$)
    changes modularity by
    $ Delta Q_"in" (i -> c) = k_(i, c) / m - (Sigma_"tot" (c) k_i) / (2 m^2). $
  + Moving $i$ directly from its current community $a = c_i$ to a distinct
    community $b != a$ yields the exact self-loop-independent net change
    $ Delta Q(i : a -> b) = (k_(i, b) - k_(i, a)) / m + (k_i (Sigma_"tot" (a) - k_i - Sigma_"tot" (b))) / (2 m^2). $
]

#proof[
  *(1)* Before inserting isolated vertex $i$ into $C_c$, communities ${i}$ and
  $C_c$ contribute
  $ Q_"before" = [ A_(i i) / (2 m) - (k_i / (2 m))^2 ] + [ (Sigma_"in" (c)) / (2 m) - ((Sigma_"tot" (c)) / (2 m))^2 ] $
  to $Q$, while all other communities are unaffected. After forming
  $C'_c = C_c union {i}$, the symmetry $A_(i j) = A_(j i)$ gives
  $ Sigma_"in" (c') = sum_(u, v in C_c union {i}) A_(u v) = Sigma_"in" (c) + A_(i i) + 2 k_(i, c), $
  and $Sigma_"tot" (c') = Sigma_"tot" (c) + k_i$. Subtracting $Q_"before"$ from
  $Q_"after" = Sigma_"in" (c') / (2 m) - (Sigma_"tot" (c') / (2 m))^2$, the
  terms $Sigma_"in" (c) / (2 m)$ and $A_(i i) / (2 m)$ cancel identically,
  leaving
  $ Delta Q_"in" (i -> c) &= (2 k_(i, c)) / (2 m) - (2 Sigma_"tot" (c) k_i) / (4 m^2) \
                          &= k_(i, c) / m - (Sigma_"tot" (c) k_i) / (2 m^2). $

  *(2)* Decompose moving $i$ from $a$ to $b != a$ into detaching $i$ from
  $C_a$ to form $C_a^- = C_a without {i}$ (whose volume is
  $Sigma_"tot" (a^-) = Sigma_"tot" (a) - k_i$) followed by inserting $i$ into
  $C_b$. The detachment changes $Q$ by $-Delta Q_"in" (i -> a^-)$, and the
  insertion changes $Q$ by $+Delta Q_"in" (i -> b)$. Summing yields
  $ Delta Q(i : a -> b) = (k_(i, b) - k_(i, a)) / m + (k_i (Sigma_"tot" (a) - k_i - Sigma_"tot" (b))) / (2 m^2), $
  which requires only $cal(O)(1)$ arithmetic operations once $k_(i, c)$ is
  accumulated over $i$'s incident edges.
]

#subsec("B", "Exact Multilevel Invariance and Finite Convergence")

The defining algorithmic guarantee of Fast Unfolding is that optimizing $Q$ on
successive aggregated super-graphs optimizes the *exact* original-graph
modularity without approximation error.

#theorem(
  title: "Exact Aggregation Invariance and Finite Monotonic Termination",
  origin: "Original Proof",
)[
  Let $G^((ell+1))$ be the Phase-2 super-graph constructed from partition
  $cal(P)^((ell))$ of $G^((ell))$.
  + *(Exact multilevel invariance)* For every partition $cal(Q)$ of
    $V^((ell+1))$ and its lifted coarsened partition $tilde(cal(Q))$ on the
    original vertex set $V^((0))$,
    $ Q(cal(Q); G^((ell+1))) = Q(tilde(cal(Q)); G^((0))). $
  + *(Finite monotonic convergence)* With $tau_"min" = 0$, every accepted move
    across all Phase-1 sweeps and all hierarchy levels strictly increases
    $Q(dot; G^((0)))$. Consequently, Fast Unfolding cannot cycle and terminates
    in finitely many steps at a partition locally optimal with respect to both
    individual vertex moves and super-node merges.
]

#proof[
  *(1)* By induction on $ell$, it suffices to verify one aggregation step from
  $G$ (partitioned into $cal(P) = {C_1, dots, C_r}$) to $G'$. The degree of
  super-node $c$ in $G'$ is
  $ k'_c = sum_(d=1)^r A'_(c d) = sum_(d=1)^r sum_(i in C_c) sum_(j in C_d) A_(i j) = sum_(i in C_c) k_i = Sigma_"tot" (c). $
  Hence total weight is strictly conserved:
  $2 m' = sum_(c=1)^r k'_c = sum_(i in V) k_i = 2 m$. For any super-community
  $S in cal(Q)$ with lifted vertex set $tilde(S) = union.big_(c in S) C_c$,
  $ Sigma'_"in" (S) &= sum_(c, d in S) A'_(c d) = sum_(i, j in tilde(S)) A_(i j) = Sigma_"in" (tilde(S)), \
    Sigma'_"tot" (S) &= sum_(c in S) k'_c = sum_(i in tilde(S)) k_i = Sigma_"tot" (tilde(S)). $
  Substituting into Definition 1 gives term-by-term identity
  $Q(cal(Q); G') = Q(tilde(cal(Q)); G)$.

  *(2)* Initializing $G^((ell+1))$ with singleton super-nodes corresponds to
  the lifted partition $cal(P)^((ell))$, so Phase 2 preserves $Q$ with zero
  jump. Each accepted Phase-1 move satisfies $Delta Q > 0$, producing a
  strictly increasing real sequence $Q(tilde(cal(P))_0; G^((0))) < Q(tilde(cal(P))_1; G^((0))) < dots$
  over the finite set of partitions of $V$ (of cardinality $B_n$, the $n$-th
  Bell number). Hence no partition can repeat, termination occurs in finite
  time, and stopping only when a full Phase-1 pass makes zero moves guarantees
  local optimality over single-vertex and super-node reassignments.
]

#subsec("C", "NP-Completeness and the Modularity Resolution Limit")

Because Fast Unfolding optimizes Newman--Girvan modularity, it inherits the
exact computational hardness and scale-dependent bias of the objective.

#theorem(
  title: "NP-Completeness and Fortunato--Barthelemy Resolution Limit",
  origin: "Classical — Cited",
)[
  + Globally maximizing $Q(cal(P); G)$ over all partitions of an undirected
    graph is NP-complete in the strong sense @brandes2008modularity.
  + Consider a ring of $p$ identical cliques $K_s$ ($s >= 3$), each linked to
    its two adjacent cliques by a single bridge edge, giving total edge count
    $m = p (binom(s, 2) + 1)$ @fortunato2007resolution. Merging two adjacent
    cliques strictly increases modularity if and only if
    $ p > s(s - 1) + 2 <==> Sigma_"tot" (K_s) < sqrt(2 m). $
    Thus modularity maximization fails to resolve well-defined cliques smaller
    than scale $Theta(sqrt(2 m))$.
]

#proof[
  NP-completeness is established by Brandes et al. @brandes2008modularity via
  reduction from 3-PARTITION. For the resolution limit
  @fortunato2007resolution, each clique $C_1, C_2 tilde.op K_s$ has
  $ell = binom(s, 2)$ internal edges and $2$ incident bridge edges, so
  $Sigma_"in" (C_1) = 2 ell$ and
  $d_C := Sigma_"tot" (C_1) = Sigma_"tot" (C_2) = 2 ell + 2 = s(s - 1) + 2$.
  Merging adjacent cliques $C_1$ and $C_2$ internalizes their shared bridge
  edge, raising $Sigma_"in"$ by $+2$ while combining their volumes into
  $2 d_C$. The exact modularity change is
  $ Delta Q_"merge" &= 2 / (2 m) - ((2 d_C) / (2 m))^2 + 2 (d_C / (2 m))^2 \
                    &= 1 / m - d_C^2 / (2 m^2) = (2 m - d_C^2) / (2 m^2). $
  Therefore $Delta Q_"merge" > 0 <=> d_C < sqrt(2 m)$. Substituting
  $2 m = p (2 ell + 2) = p d_C$ yields $p d_C > d_C^2 <=> p > d_C = s(s - 1) + 2$.
  For instance, on a ring of $K_5$ cliques ($d_C = 22$), any ring with
  $p >= 23$ cliques strictly favors pairwise merged pairs over the true
  $K_5$ partition.
]

#subsec("D", "Constructive Proof of Internally Disconnected Communities")

Traag et al. @traag2019leiden observed empirically that Louvain can output
badly connected or even *disconnected* communities when bridge vertices migrate
out of a community during Phase 1. Below we give a self-contained, arithmetic-verified
unweighted counterexample achieving this pathology at a strict local optimum.

#theorem(
  title: "Strict Phase-1 Local Optimum Containing a Disconnected Community",
  origin: "Original Counterexample",
)[
  There exists a connected unweighted simple graph $G$ on $n = 122$ vertices
  and a valid Phase-1 visitation order under which Fast Unfolding converges to
  a strict local optimum containing a community $C_1$ whose induced subgraph
  $G[C_1] = K_4 union.sq K_4$ is disconnected (zero edges between its two
  components), which Phase 2 then permanently locks into a single super-node.
]

#proof[
  Construct the connected unweighted graph $G$ from four cliques
  $A = {a_1, dots, a_4} tilde.op K_4$, $B = {b_1, dots, b_4} tilde.op K_4$,
  $D = {d_0, d_1, dots, d_(11)} tilde.op K_(12)$, and $H tilde.op K_(100)$,
  together with two adjacent bridge vertices $u, w$ (${u, w} in E$):
  - Connect both $u$ and $w$ to all $4$ vertices of $A$, all $4$ vertices of
    $B$, and the $10$ vertices ${d_1, dots, d_(10)} subset D$ (leaving $d_0$
    and $d_(11)$ non-adjacent to $u, w$).
  - Connect $d_0$ to a single vertex $h_1 in H$ by one bridge edge.

  Direct counting yields vertex degrees $k_(a_i) = k_(b_i) = 5$,
  $k_u = k_w = 19$, $k_(d_0) = 12$, $k_(d_(11)) = 11$,
  $k_(d_1) = dots = k_(d_(10)) = 13$, and total edges
  $m = 6 + 6 + 66 + 1 + 8 + 8 + 20 + 1 + 4950 = 5066$ ($2 m = 10132$).
  Every vertex outside $H$ has $k_i <= 19$ and combined volume
  $Sigma_"tot" (A union B union {u, w} union D) = 231$, so the null-model term
  in Lemma 1 obeys $k_i Sigma_"tot" / (2 m^2) < 4389 / (10132 m) < 0.44 / m$,
  ensuring every unit edge difference ($k_(i, b) - k_(i, a) >= 1$) strictly
  governs move direction.

  Trace Phase 1 from singletons in Sweep 1:
  + Visiting $u$ moves $u$ into minimal-degree neighbor ${a_1}$
    ($k_(a_1) = 5$), forming $C_1 = {u, a_1}$. Next $w$ has $k_(w, C_1) = 2$
    vs. $1$ to singletons, joining $C_1$. Visiting $a_2, a_3, a_4$ yields
    $k_(a_j, C_1) >= 3 > 1$, so $C_1 = {u, w} union A$.
  + Visiting $b_1$ gives $k_(b_1, C_1) = 2$ (to $u, w$) vs. $1$ to singletons
    ${b_2}, {b_3}, {b_4}$, so $b_1$ strictly joins $C_1$. Subsequently
    $b_2, b_3, b_4$ each have $k_(b_j, C_1) >= 3$, yielding
    $C_1 = {u, w} union A union B$ with $Sigma_"tot" (C_1) = 78$.
  + Visiting $d_0$ (non-adjacent to $u, w$) moves $d_0$ into ${d_(11)}$,
    forming $C_D = {d_0, d_(11)}$ ($Sigma_"tot" (C_D) = 23$). Visiting $d_1$
    yields $k_(d_1, C_D) = k_(d_1, C_1) = 2$; since
    $Sigma_"tot" (C_D) = 23 < 78 = Sigma_"tot" (C_1)$, Lemma 1 strictly favors
    $C_D$. Then $d_2, dots, d_(10)$ each have $k_(d_j, C_D) >= 3 > 2$, giving
    $C_D = D$ ($Sigma_"tot" (C_D) = 153$). Visiting $H$ forms $C_H = H$.

  Now trace Sweep 2:
  + Visiting $u in C_1$ ($k_(u, C_1^-) = 9$, $Sigma_"tot" (C_1^-) = 59$,
    $k_(u, C_D) = 10$):
    $ Delta Q(u : C_1 -> C_D) = 1 / m + (19 (59 - 153)) / (10132 m) = 0.8237 / m > 0, $
    so $u$ strictly migrates to $C_D$.
  + Visiting $w in C_1$ ($k_(w, C_1^-) = 8$, $Sigma_"tot" (C_1^-) = 40$,
    $k_(w, C_D) = 11$, $Sigma_"tot" (C_D) = 172$):
    $ Delta Q(w : C_1 -> C_D) = 3 / m + (19 (40 - 172)) / (10132 m) = 2.7524 / m > 0, $
    so $w$ strictly migrates to $C_D$.

  Verify strict local optimality of the resulting partition
  ${C_1 = A union B, C_D = D union {u, w}, C_H = H}$: every $a_i in A$ has
  $k_(a_i, C_1^-) = 3$ internal edges versus $k_(a_i, C_D) = 2$ edges to
  $C_D$ and $Sigma_"tot" (C_1^-) = 35 < 191 = Sigma_"tot" (C_D)$, so moving
  $a_i$ to $C_D$ (or a singleton) strictly decreases $Q$. By symmetry, no
  $b_i in B$ can move. Vertices in $C_D$ have $>= 11$ internal vs. $<= 8$
  external edges, and vertices in $H$ have $99$ internal vs. $<= 1$ external
  edge. Thus Phase 1 terminates cleanly with the disconnected community
  $C_1 = A union B$ ($G[C_1] = K_4 union.sq K_4$).
]

#sec("III", "Complexity, Geometry and Operational Profile")

Fast Unfolding owes its extraordinary speed to three structural efficiencies:
(i) each Phase-1 sweep over $G^((ell))$ visits every vertex $i$ once and scans
its adjacency list twice (once to accumulate $k_(i, c)$ into a sparse scratch
array, once to clear touched entries), taking exact $cal(O)(m_ell + n_ell)$
time per sweep; (ii) Phase-2 aggregation constructs $G^((ell+1))$ in a single
pass over the edges of $G^((ell))$ in $cal(O)(m_ell + n_ell)$ time via hash or
bucket accumulation; and (iii) empirical community sizes grow geometrically
across levels ($n_(ell+1) << n_ell$), so the first level $ell = 0$ dominates
total runtime and the number of hierarchy levels typically satisfies
$L = cal(O)(log n)$ @blondel2008louvain.

Across $L$ levels with at most $T$ Phase-1 sweeps per level, total runtime is
$cal(O)(L T m)$ (empirically near-linear $cal(O)(m)$ or $cal(O)(n log n)$ on
sparse networks, scaling to $10^8$ vertices and $10^9$ edges in minutes).
Space complexity is linear, $cal(O)(n + m)$, requiring only CSR adjacency
arrays, current community labels, cached volumes $Sigma_"tot" (c)$, and an
$cal(O)(n)$ scratch buffer for $k_(i, c)$. Because Phase 1 visits vertices in
order, the output partition is *order-dependent* (randomizing node order per
run is standard practice to probe the degenerate modularity landscape).
Moreover, the method is scale-invariant under positive global scaling
$A |-> alpha A$ ($alpha > 0$) since $Q$ depends only on normalized ratios
$A_(i j) / (2 m)$.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time complexity], [$cal(O)(L T m)$; empirically $cal(O)(m)$ / $cal(O)(n log n)$],
    [Space complexity], [$cal(O)(n + m)$ (CSR + volume array + scratch map)],
    [Partition type], [Hard multilevel hierarchy of nested flat partitions],
    [Determinism], [Order-dependent; greedy local optimum varies by seed],
    [Resolution limit], [Present: merges modules with $Sigma_"tot" < sqrt(2 m)$],
    [Connectivity], [Not guaranteed (Section II-D disconnected communities)],
    [Hyperparameters], [None (or resolution $gamma > 0$ in Potts generalization)],
    [Scale invariance], [Exact under uniform weight scaling $A |-> alpha A$],
  ),
  caption: [Operational and complexity profile of Fast Unfolding (Alg. 129).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

Fast Unfolding has become the de facto baseline across large-scale network
science whenever assortative, hard-partitioned mesoscale structure must be
extracted at linear cost:

+ *Telecommunications Call Graphs @blondel2008louvain:* Blondel et al.
  originally benchmarked Fast Unfolding on a 2.6-million-user Belgian mobile
  phone network, cleanly unfolding the Walloon / Flemish linguistic boundary
  and nested regional hierarchies in seconds because weighted call volumes
  strongly satisfy the assortative configuration null model.
+ *Single-Cell RNA Sequencing (scRNA-seq) Phenotyping:* In bioinformatics
  pipelines (e.g., Seurat / Scanpy), cells are embedded as a sparse $k$-NN or
  shared-nearest-neighbor (SNN) graph; Fast Unfolding rapidly partitions
  millions of transcriptomes into discrete cell types without prespecifying
  cluster count $r$.
+ *Web Graph and Citation Macro-Mapping @blondel2008louvain:* On web link
  graphs exceeding $10^8$ pages and $10^9$ hyperlinks, the $cal(O)(L T m)$
  streaming sweep and multilevel coarsening extract topical web domains and
  scientific disciplines where spectral methods @newman2006leading exhaust RAM.
+ *Knowledge-Graph and Entity-Resolution Pre-Partitioning:* Before running
  expensive pairwise entity matching or graph neural network mini-batching,
  Fast Unfolding shards massive graphs into low-cut balanced partitions,
  minimizing cross-machine communication volume.
+ *Brain Structural and Functional Connectomics:* Multilevel modularity
  unfolding maps nested cortical hierarchies from voxel-level fMRI correlation
  graphs into resting-state macro-systems while exposing intermediate
  resolutions $cal(P)^((0)), dots, cal(P)^((L))$.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Four structural failure modes dictate when Fast Unfolding should be modified
or replaced: (i) *multi-scale small communities* below the $sqrt(2 m)$
resolution limit (Section II-C), mitigated by tuning resolution $gamma > 1$
@reichardt2006spinglass or switching to Infomap @rosvall2008infomap;
(ii) *disconnected or poorly connected communities* (Section II-D), remedied by
the Leiden algorithm @traag2019leiden; (iii) *disassortative, bipartite, or
core--periphery networks*, where maximizing internal edge surplus fails
completely and SBM (Alg. 130) is required; and (iv) *overlapping communities*,
where CPM (Alg. 127) or DEMON (Alg. 128) is preferable.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Massive web / telecom], [High: $cal(O)(m)$ speed + nested hierarchy], [Resolution limit; use Leiden @traag2019leiden],
    [scRNA-seq $k$-NN], [High: fast auto-$r$ cell-type discovery], [Tune resolution $gamma$; check connectivity],
    [Multi-scale small nets], [Moderate: merges small cliques $< sqrt(2 m)$], [Use Infomap @rosvall2008infomap or Walktrap @pons2006walktrap],
    [Bipartite / disassort.], [Poor: assortative objective $Q$ fails], [Use SBM (Alg. 130) @karrer2011dcsbm],
    [Overlapping hubs], [Poor: hard partition forces bridge splits], [Use CPM (Alg. 127) or DEMON (Alg. 128)],
  ),
  caption: [Domain suitability and failure-mode matrix for Fast Unfolding (Alg. 129).],
)

#sec("V", "Structural Relations and Empirical Behavior")

Fast Unfolding directly links classical modularity heuristics
@clauset2004cnm @newman2006leading, statistical physics Potts models
@reichardt2006spinglass, and modern guaranteed-connected multilevel refinement
@traag2019leiden.

#proposition(
  title: "Equivalence to Potts Spin-Glass Hamiltonian and Leiden Refinement",
  origin: "Original Derivation",
)[
  Introduce a positive resolution parameter $gamma > 0$ into modularity:
  $ Q_gamma (cal(P); G) = 1 / (2 m) sum_(i, j = 1)^n (A_(i j) - gamma (k_i k_j) / (2 m)) delta(c_i, c_j). $
  Then:
  + Minimizing the Reichardt--Bornholdt Potts spin-glass Hamiltonian
    $cal(H)_gamma (cal(P)) = -sum_(i, j) (A_(i j) - gamma k_i k_j / (2 m)) delta(c_i, c_j)$
    @reichardt2006spinglass is identically equivalent to maximizing
    $Q_gamma (cal(P); G)$, with Lemma 1 generalizing simply by multiplying the
    quadratic volume term by $gamma$.
  + Inserting a *constrained refinement phase* between Phase 1 and Phase 2
    that re-partitions each $C in cal(P)^((ell))$ into well-connected
    sub-communities before aggregation (the Leiden algorithm @traag2019leiden)
    splits the disconnected community $C_1 = A union B$ of Section II-D into
    distinct super-nodes $A$ and $B$ and guarantees every output community is
    internally connected.
]

#proof[
  *(1)* Dividing $cal(H)_gamma (cal(P))$ by $-2 m < 0$ yields
  $Q_gamma (cal(P); G)$ identically. Replacing $k_i k_j / (2 m)$ by
  $gamma k_i k_j / (2 m)$ scales the quadratic volume penalty in Lemma 1 by
  $gamma$, while Section II-B's aggregation invariance remains exact because
  $Sigma'_"tot" (S) = Sigma_"tot" (tilde(S))$ is independent of $gamma$.

  *(2)* In Leiden @traag2019leiden, refinement initializes singletons within
  each Phase-1 community $C$ and permits merges only along intra-$C$ edges with
  positive gain. On Section II-D's $C_1 = A union B$, no edge connects $A$ to $B$,
  so refinement partitions $C_1$ into connected sub-communities $A$ and $B$
  prior to Phase 2, strictly increasing $Q$ by
  $2 (Sigma_"tot" (A) Sigma_"tot" (B)) / (4 m^2) > 0$ and restoring
  connectedness.
]

*Empirical Benchmark Analysis (ClusterOverview 4-Dataset Suite).* Across the
four ClusterOverview $k$-NN proximity benchmarks (Silhouette
@rousseeuw1987silhouettes and ARI @hubert1985comparing):
- *Iris (`Sil = 0.586`, `ARI = 0.567`):* Ties CPM (Alg. 127) at the top
  Silhouette and ARI among the four graph algorithms, cleanly unfolding
  Setosa while partially merging the contiguous Versicolor--Virginica region.
- *Customers (`Sil = 0.190`):* Lowest Silhouette among the four methods (below
  DEMON `0.250`, CPM `0.262`, and SBM `0.356`), directly manifesting Section II-C's
  resolution limit as small customer micro-segments coalesce into coarse
  modularity blocks.
- *Multi-Shape (`Sil = 0.713`, `ARI = N/A`):* Ties SBM (Alg. 130) for the best
  Silhouette, outperforming CPM (`0.554`) and DEMON (`0.396`) by smoothly
  aggregating curved high-conductance density bands across multilevel coarsening.
- *High-D (`Sil = 0.188`, `ARI = N/A`):* Ties CPM (`0.188`) and trails SBM
  (`0.303`), since high-dimensional hub edges equalize intra- and inter-cluster
  modularity gradients.
- *Runtime (`Mean = 0.681 ms`):* The fastest of all four network algorithms by
  a wide margin—$4.1 times$ faster than CPM (`2.800 ms`), $6.5 times$ faster
  than DEMON (`4.424 ms`), and $191.7 times$ faster than SBM (`130.533 ms`)—confirming
  its near-linear $cal(O)(L T m)$ computational superiority.
