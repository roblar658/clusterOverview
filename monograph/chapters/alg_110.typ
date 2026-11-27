#import "../template/ieee.typ": *

#alg_entry(
  id: 110,
  name: "Label Propagation (LPA)",
  category: "Cat. V: Graph, Spectral & Network",
  year: "2007",
  paradigm: "Asynchronous Local Majority-Vote Dynamics",
  complexity: $O(T m)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

The *Label Propagation Algorithm (LPA)* introduced by Raghavan, Albert,
and Kumara @raghavan2007lpa detects communities in near-linear time
using neither a pre-specified cluster count $K$ nor any global
objective hyperparameters. Instead, each vertex initializes with a
unique identifier and iteratively adopts the label carrying maximum
total edge weight among its immediate neighbors, driving the graph
toward a local consensus equilibrium.

Let $G = (V, E, w)$ be an undirected weighted graph on $n = |V|$
vertices and $m = |E|$ edges with symmetric nonnegative adjacency
matrix $A = (A_(i j)) in RR_(>= 0)^(n times n)$, zero diagonal
($A_(i i) = 0$), and total edge weight
$m_w = (1 / 2) sum_(i, j in V) A_(i j)$. Let
$ell : V -> {1, dots, n}$ denote a discrete vertex labeling, and for
each vertex $i in V$ and candidate label $c in {1, dots, n}$, define
the incident label weight
$ s_i (c; ell) = sum_(j in N(i)) A_(i j) delta(c, ell_j),
  quad M_i (ell) = arg max_(c in {1, dots, n}) s_i (c; ell). $

#definition(title: "Ferromagnetic Potts Energy and Local Majority Equilibrium")[
  For any vertex labeling $ell in {1, dots, n}^V$, define the
  *monochromatic-edge (ferromagnetic Potts) energy*
  $ Phi(ell) = 1 / 2 sum_(i, j in V) A_(i j) delta(ell_i, ell_j)
    = m_w - "cut"(ell), $
  where $"cut"(ell) = (1 / 2) sum_(i, j) A_(i j) (1 - delta(ell_i, ell_j))$
  is the total weight of inter-label edges. A labeling $ell^*$ is a
  *local majority equilibrium* (fixed point of LPA) if
  $ell_i^* in M_i (ell^*)$ for every $i in V$, i.e.,
  $ sum_(j in N(i)) A_(i j) delta(ell_i^*, ell_j^*)
    >= sum_(j in N(i)) A_(i j) delta(c, ell_j^*)
    quad forall i in V, forall c. $
]

Unlike modularity optimization (Algs. 107--108), $Phi(ell)$ contains
no antiferromagnetic community-size penalty; as we prove in Section II,
this eliminates hyperparameters and enables strict $O(m)$-step
asynchronous Lyapunov ascent, while creating two characteristic
phenomena: synchronous bipartite $2$-cycles and monster-community
collapse on sparse weakly modular graphs.

#algbox(
  title: "Asynchronous Label Propagation Algorithm (RAK-LPA)",
  inputs: [Weighted undirected graph $G = (V, E, w)$, max sweeps $T_("max")$.],
  outputs: [Community partition $cal(P)$ of $V$.],
)[
  1. Initialize unique labels $ell_i <- i$ for all $i = 1, dots, n$; set
     sweep counter $t <- 0$.
  2. While $t < T_("max")$: draw a uniformly random permutation
     $pi^((t))$ of $V$ and set $"changed" <- "false"$.
  3. For $r = 1, dots, n$, let $i = pi^((t)) (r)$:
     a. Accumulate $s_i (c) <- sum_(j in N(i)) A_(i j) delta(c, ell_j)$
        over all labels $c$ present in $N(i)$.
     b. Compute $M_i = arg max_c s_i (c)$. Under conservative
        tie-breaking, if $ell_i in.not M_i$, draw $c^* in M_i$ uniformly
        at random, set $ell_i <- c^*$, and set $"changed" <- "true"$.
  4. If $"changed" = "false"$ (every $i in V$ satisfies $ell_i in M_i$),
     break; else set $t <- t + 1$.
  5. *Connected-Component Post-Processing*: partition each non-empty
     label class $L_c = {i in V : ell_i = c}$ into connected components
     of $G[L_c]$ via BFS / DFS in $O(n + m)$ time and return $cal(P)$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We prove finite-time Lyapunov convergence with a tight $m$-update
bound on unweighted graphs, establish an exact period-$<= 2$ barrier
for synchronous LPA, characterize strong-community forward invariance
versus monster-community collapse, and exhibit a minimal 10-vertex
counterexample proving the necessity of Step 5 post-processing.

#subsec("A", "Lyapunov Monotonicity and Finite Convergence of Asynchronous LPA")

#theorem(
  title: "Monotone Potts Ascent and Exact m-Update Bound",
  origin: "Original Proof",
)[
  Under asynchronous LPA with conservative tie-breaking ($i$ changes
  label only when $ell_i in.not M_i (ell)$):
  (i) Every non-trivial update $ell_i -> c^*$ strictly increases the
  Potts energy by $Delta Phi = s_i (c^*; ell) - s_i (ell_i; ell) > 0$.
  (ii) On any finite weighted graph, asynchronous LPA terminates in
  finitely many updates at a local majority equilibrium.
  (iii) On any unweighted simple graph ($A_(i j) in {0, 1}$), the
  *total* number of label changes across the entire execution is
  strictly bounded by $m = |E|$.
]

#proof[
  (i) When vertex $i$ changes from $ell_i$ to $c^* != ell_i$ while all
  other labels remain fixed (asynchrony), only edges incident to $i$
  change their monochromatic status. Using $A_(i i) = 0$ and symmetry
  $A_(i j) = A_(j i)$, incident edges matching $ell_i$ (total weight
  $s_i (ell_i; ell)$) become bichromatic, while incident edges matching
  $c^*$ (total weight $s_i (c^*; ell)$) become monochromatic. Hence
  $ Phi(ell') - Phi(ell) = s_i (c^*; ell) - s_i (ell_i; ell). $
  Because conservative tie-breaking updates $i$ only if
  $ell_i in.not M_i (ell)$, we have
  $s_i (c^*; ell) > s_i (ell_i; ell)$, proving $Delta Phi > 0$.

  (ii) The set of attainable values ${Phi(ell) : ell in {1, dots, n}^V}$
  is finite and bounded above by $m_w$. Since $Phi$ strictly increases
  at every update, no labeling can recur and the algorithm halts in
  finitely many updates at a state where $ell_i in M_i (ell)$ for all
  $i in V$.

  (iii) On an unweighted graph, $Phi(ell) in {0, 1, dots, m}$ is
  integer-valued, initialized at $Phi(ell^((0))) = 0$ (all labels
  distinct). Every update raises $Phi$ by an integer $Delta Phi >= 1$.
  Therefore the cumulative count $U$ of label updates satisfies
  $U <= Phi(ell^((T))) - Phi(ell^((0))) <= m$.
]

#subsec("B", "Synchronous Oscillation and Exact Period-2 Barrier")

Raghavan et al. @raghavan2007lpa observed empirically that updating all
vertices simultaneously (*synchronous* LPA) can enter persistent
oscillations on bipartite subgraphs. We prove both a minimal
singleton-initialized $2$-cycle and a universal upper bound showing
that periods $>= 3$ are impossible on undirected graphs.

#theorem(
  title: "Bipartite 2-Cycle and Universal Period-2 Barrier for Synchronous LPA",
  origin: "Original Proof",
)[
  Consider synchronous LPA, where all vertices update simultaneously via
  $ell_i^((t + 1)) in M_i (ell^((t)))$.
  (i) *Constructive $2$-cycle from singletons*: on the 4-cycle
  $C_4 = (V_1 union V_2, E)$ with $V_1 = {1, 3}$, $V_2 = {2, 4}$,
  standard initialization $ell_i^((0)) = i$ and deterministic
  smallest-label tie-breaking enters a permanent period-$2$ cycle at
  $t = 1$.
  (ii) *Universal period-$<= 2$ theorem*: on *any* finite undirected
  weighted graph under conservative two-step tie-breaking (if
  $ell_i^((t - 1)) in M_i (ell^((t)))$, set $ell_i^((t + 1)) = ell_i^((t - 1))$),
  synchronous LPA converges in finitely many steps to either a fixed
  point (period $1$) or a period-$2$ limit cycle; cycles of period
  $>= 3$ cannot occur.
]

#proof[
  (i) At $t = 0$, $V_1$ holds $(1, 3)$ and $V_2$ holds $(2, 4)$. At
  $t = 1$, every $v in V_1$ sees neighbor labels ${2, 4}$ and chooses
  $min = 2$; every $u in V_2$ sees ${1, 3}$ and chooses $min = 1$.
  Thus $(ell_(V_1)^((1)), ell_(V_2)^((1))) = ((2, 2), (1, 1))$. At
  $t = 2$, $V_1$ sees ${1, 1} |-> 1$ and $V_2$ sees ${2, 2} |-> 2$,
  yielding $((1, 1), (2, 2))$. At $t = 3$, the state returns to
  $((2, 2), (1, 1))$, oscillating forever with period $2$.

  (ii) Define the cross-time Lyapunov functional
  $ Psi(ell^((t)), ell^((t - 1)))
    = sum_(i, j in V) A_(i j) delta(ell_i^((t)), ell_j^((t - 1))). $
  By symmetry $A_(i j) = A_(j i)$, rewrite
  $Psi(ell^((t)), ell^((t - 1))) = sum_(i in V) sum_(j in N(i)) A_(i j) delta(ell_i^((t - 1)), ell_j^((t)))$.
  Computing the one-step increment gives
  $ &Psi(ell^((t + 1)), ell^((t))) - Psi(ell^((t)), ell^((t - 1))) \
    &= sum_(i in V) [s_i (ell_i^((t + 1)); ell^((t)))
    - s_i (ell_i^((t - 1)); ell^((t)))]. $
  Because $ell_i^((t + 1)) in M_i (ell^((t)))$ maximizes
  $s_i (dot; ell^((t)))$, every summand is nonnegative. Under
  conservative two-step tie-breaking, if $ell^((t + 1)) != ell^((t - 1))$,
  at least one vertex $i$ has $ell_i^((t - 1)) in.not M_i (ell^((t)))$,
  making its summand strictly positive and hence
  $Psi(ell^((t + 1)), ell^((t))) > Psi(ell^((t)), ell^((t - 1)))$.
  Since $Psi <= 2 m_w$ assumes finitely many values, there exists $T_0$
  such that $ell^((t + 1)) = ell^((t - 1))$ for all $t >= T_0$, proving
  that every trajectory settles onto period $1$ or period $2$.
]

#subsec("C", "Global Triviality and Strong-Community Forward Invariance")

#proposition(
  title: "Trivial Global Optimum and Radicchi Strong-Community Stability",
  origin: "Original Proof",
)[
  (i) Every constant labeling $ell equiv c_0$ achieves the global
  maximum $Phi(ell) = m_w$ (zero cut); non-trivial LPA communities
  exist exclusively as metastable local maxima of $Phi$.
  (ii) Let $C subset.eq V$ be a *strong community* satisfying
  $ k_(i, C) = sum_(j in N(i) inter C) A_(i j)
    > sum_(j in N(i) without C) A_(i j) = k_(i, V without C)
    quad forall i in C. $
  Once every vertex in $C$ holds a common label $c_0$, the consensus
  on $C$ is *strictly forward-invariant* under both asynchronous and
  synchronous LPA: no vertex in $C$ can ever abandon $c_0$, regardless
  of external label assignments or update order.
]

#proof[
  (i) Because $delta(ell_i, ell_j) <= 1$ and $A_(i j) >= 0$,
  $Phi(ell) <= (1 / 2) sum_(i, j) A_(i j) = m_w$, attained with
  equality whenever $ell_i = ell_j$ on every connected component.
  (ii) Suppose $ell_j = c_0$ for all $j in C$. Take any $i in C$ and
  any competing label $c != c_0$. Every neighbor of $i$ holding label
  $c$ must lie outside $C$, whereas every neighbor in $C$ holds $c_0$.
  Hence
  $ s_i (c; ell) <= sum_(j in N(i) without C) A_(i j) = k_(i, V without C)
    < k_(i, C) <= s_i (c_0; ell). $
  Therefore $M_i (ell) = {c_0}$ is a strict singleton for every
  $i in C$, ensuring $ell_i = c_0$ persists across all future steps.
]

Proposition II.3 explains both the strength and Achilles' heel of LPA:
planted strong communities (e.g., dense stochastic blocks
@holland1983sbm @abbe2018sbmreview) lock permanently once locally
synchronized, whereas in sparse weakly modular graphs where
$k_(i, C) <= k_(i, V without C)$, a single early label can percolate
across bridges and collapse the network into one *monster community*
@raghavan2007lpa.

#subsec("D", "Disconnected Label-Class Counterexample")

Just as Louvain can emit disconnected communities, raw LPA label
classes $L_c = {i : ell_i = c}$ are not guaranteed to be connected,
even at a *strict* local majority equilibrium where every vertex has a
unique winning label.

#proposition(
  title: "Minimal 10-Vertex Disconnected Strict Equilibrium of LPA",
  origin: "Original Counterexample",
)[
  There exists a connected 10-vertex unweighted graph $G_(10)$ on which
  asynchronous LPA initialized from singletons converges to a *strict*
  local majority equilibrium ($|M_i (ell^*)| = 1$ for all $i in V$)
  whose label class $L_u = A union B$ consists of two disconnected
  components separated by a vertex $u$ holding a different label.
]

#proof[
  Let $G_(10)$ comprise three triangles $A = {a_1, a_2, a_3} tilde.op K_3$,
  $B = {b_1, b_2, b_3} tilde.op K_3$, $D = {d_1, d_2, d_3} tilde.op K_3$,
  and a single hub $u$ adjacent to $a_1, b_1, d_1, d_2, d_3$ ($14$ edges
  total). Starting from singleton labels $ell_v^((0)) = v$, execute
  asynchronous updates in the admissible sweep order
  $(a_1, a_2, a_3, b_1, b_2, b_3, d_1, d_2, d_3, u)$:
  (1) $a_1$ ties over ${u, a_2, a_3}$ and adopts $u$; $a_2$ ties over
  ${u, a_3}$ and adopts $u$; $a_3$ has strict majority $s_(a_3) (u) = 2$
  and adopts $u$. Now all of $A$ holds label $u$.
  (2) Symmetrically, $b_1, b_2, b_3$ adopt $u$, so all of $B$ holds $u$.
  (3) $d_1$ ties over ${u, d_2, d_3}$ and adopts $d_2$; $d_2$ has
  $s_(d_2) (d_2) = 1$ (from $d_1$), $s_(d_2) (d_3) = 1$, $s_(d_2) (u) = 1$,
  retaining $d_2$; $d_3$ has strict majority $s_(d_3) (d_2) = 2 > 1$
  and adopts $d_2$. Now all of $D$ holds $d_2$.
  (4) Finally, $u$ sees neighbor labels ${u, u, d_2, d_2, d_2}$, so
  $s_u (d_2) = 3 > 2 = s_u (u)$, and $u$ strictly adopts $d_2$.

  Now verify strict local majority equilibrium at every vertex:
  $a_1$ and $b_1$ each have $2$ neighbors with $u$ vs $1$ with $d_2$;
  $a_2, a_3, b_2, b_3$ each have $2$ neighbors with $u$; $u$ has $3$
  neighbors with $d_2$ vs $2$ with $u$; and $d_1, d_2, d_3$ each have
  $3$ neighbors with $d_2$. Thus $|M_i (ell^*)| = 1$ for all $10$
  vertices, yet $L_u = A union B$ has zero edges between $A$ and $B$.
  Step 5 BFS post-processing is therefore mathematically necessary to
  guarantee connected communities.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Each asynchronous sweep visits every edge twice ($2 m$ neighbor reads),
accumulating neighbor label weights into an $O(n)$ scratch buffer
cleared in $O("deg"(i))$ time per vertex, achieving strictly $O(m)$
time per sweep and $O(n + m)$ space. Empirically, $95\%$ of vertices
stabilize within $T approx 5$ sweeps @raghavan2007lpa, yielding
near-linear total runtime $O(T m)$ that scales seamlessly to
billion-edge networks. Moreover, graph coloring into $chi$ independent
sets permits lock-free parallel updates inside each color class while
preserving Theorem II.1's exact asynchronous Lyapunov monotonicity
(since no two simultaneously updated vertices share an edge).

Geometrically, LPA is a zero-temperature ferromagnetic Potts dynamics
with zero parameters: it imposes no spherical, balanced-size, or
global-scale prior, allowing arbitrary non-convex manifold clusters
to emerge whenever internal k-NN connectivity locally dominates
boundary degree (Proposition II.3). Its primary operational trade-off
is high stochastic variance across random update permutations and
tie-breaking choices.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time complexity], [$O(T m)$ ($O(m)$ per sweep; $<= m$ total unweighted changes).],
    [Space complexity], [$O(n + m)$ CSR graph plus $O(n)$ label / vote scratch arrays.],
    [Objective function], [Local maxima of Potts energy $Phi(ell) = m_w - "cut"(ell)$.],
    [Hyperparameters], [None (parameter-free; optional $T_("max")$ cap).],
    [Determinism], [Highly stochastic (depends on random node order and ties).],
    [Invariances], [Positive weight scale-invariant; permutation-equivariant.],
    [Parallelism], [Independent-set coloring preserves exact Lyapunov ascent.],
    [Structural caveats], [Monster-community collapse; Step 5 BFS needed (Prop. II.4).],
  ),
  caption: [Operational and theoretical profile of Label Propagation (LPA).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Billion-Edge Social Network Community Detection:* On distributed
   graph engines (Pregel, Apache Spark GraphX), LPA's purely local
   message-passing stencil avoids global synchronization of $m_w$ or
   multilevel coarsening tables, enabling rapid community discovery on
   graphs with $10^9+$ edges.
2. *Real-Time Semi-Supervised Fraud and Spam Propagation:* Pinning
   known fraudulent or benign seed accounts as fixed labels turns LPA
   into harmonic / discrete semi-supervised label diffusion, propagating
   risk labels across dense transactional neighborhoods in milliseconds.
3. *Non-Convex k-NN Manifold Segmentation:* When feature data lie on
   intertwined non-convex manifolds (moons, spirals, concentric shells),
   local majority voting propagates cleanly along high-density k-NN
   sheets without suffering modularity's resolution limit or quadratic
   size penalties.
4. *Fast Pre-Clustering of Massive Biological Interaction Graphs:*
   Before running expensive $O(n k^2)$ flow methods (MCL, Alg. 109) or
   stochastic block model inference @karrer2011dcsbm, LPA partitions
   massive interactomes into coarse strongly connected blocks in
   $O(m)$ time.
5. *Enterprise Entity Resolution and Deduplication:* Record-linkage
   similarity graphs consist of thousands of tiny, tightly connected
   duplicate cliques. Because each small clique is a strong community
   (Proposition II.3), LPA isolates tiny duplicate clusters accurately
   where modularity's resolution limit would erroneously merge them.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

LPA exhibits four principal failure modes: (a) *monster-community
collapse* on sparse, noisy, or Erdős--Rényi-like graphs where weak
boundaries fail the Radicchi strong-community inequality
(Proposition II.3), collapsing $> 90\%$ of vertices into one label;
(b) *high stochastic instability*, where different random seeds yield
markedly different local maxima of $Phi$ (requiring ensemble consensus
@monti2003consensus or deterministic tie-breaking); (c) *absence of
resolution control*, making it impossible to sweep granularity without
modifying the objective; and (d) *synchronous $2$-cycles* on bipartite
subgraphs if asynchronous or colored updates are omitted (Theorem II.2).

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Billion-edge social nets], [Ideal: $O(T m)$ local message passing], [Guard against monster collapse; use Leiden],
    [Non-convex k-NN manifolds], [Strong: local majority follows density sheets], [Sensitive to noisy cross-manifold bridges],
    [Entity deduplication], [High: preserves tiny strong cliques], [Run Step 5 BFS split (Prop. II.4)],
    [Weakly modular / noisy nets], [Poor: collapses into giant trivial label], [Use Leiden (Alg. 108) or MCL (Alg. 109)],
    [Multi-scale hierarchy needs], [Poor: flat output, no $gamma$ knob], [Use Louvain / Leiden or Walktrap],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for LPA.],
)

#sec("V", "Structural Relations and Empirical Behavior")

LPA sits at the zero-penalty boundary of Potts-based graph clustering,
unifying directly with Chinese Whispers @biemann2006chinese,
Speaker-Listener LPA (SLPA) for overlapping communities @xie2011slpa,
Fluid Communities @pares2018fluidc, and the Louvain / Leiden family
(Algs. 107--108). The following exact reduction formalizes its place
in the community-detection hierarchy.

#proposition(
  title: "Exact Equivalence of LPA to Resolution-Zero Louvain Phase 1",
  origin: "Original Derivation",
)[
  Let $Delta Q_gamma (i : C_a -> C_b)$ be Louvain's exact single-node
  modularity gain (Lemma 107.II.1) on a loopless weighted graph. In the
  zero-resolution limit $gamma -> 0^+$, every vertex move satisfies
  $ lim_(gamma -> 0^+) m_w Delta Q_gamma (i : C_a -> C_b)
    = s_i (b; ell) - s_i (a; ell) = Delta Phi(i : a -> b). $
  Consequently, asynchronous LPA with conservative tie-breaking is
  algebraically identical to Phase 1 of Louvain (and zero-temperature
  Potts Glauber dynamics @reichardt2006spinglass) at $gamma = 0$.
]

#proof[
  By Lemma 107.II.1, for any $i in C_a$ and $b != a$,
  $ m_w Delta Q_gamma (i : C_a -> C_b)
    = (k_(i, b) - k_(i, a))
    + (gamma k_i) / (2 m_w) (Sigma_("tot")^((-i)) (a) - Sigma_("tot") (b)). $
  Because $A_(i i) = 0$, the incident cluster weight $k_(i, c)$ equals
  $sum_(j in N(i)) A_(i j) delta(c, ell_j) = s_i (c; ell)$. Taking
  $gamma -> 0^+$ annihilates the degree-product null-model term,
  leaving $k_(i, b) - k_(i, a) = s_i (b; ell) - s_i (a; ell) = Delta Phi(i : a -> b)$.
  Hence $arg max_c Delta Q_0 (i : a -> c) = M_i (ell)$, proving exact
  equivalence of greedy update trajectories up to tie-breaking.
]

Viewed through Proposition V.1, the $gamma > 0$ configuration-model
penalty of Louvain (Alg. 107) and Leiden (Alg. 108) is precisely the
antiferromagnetic repulsive field added to LPA's ferromagnetic Potts
Hamiltonian $-Phi(ell)$ to prevent monster-community collapse, at the
cost of introducing a resolution parameter and (under modularity) the
resolution limit. Conversely, Chinese Whispers @biemann2006chinese is
essentially unweighted asynchronous LPA applied to lexical graphs,
while SLPA @xie2011slpa equips each vertex with a finite memory of
received labels to recover overlapping communities.

*Empirical Behavior (ClusterOverview 4-Dataset Benchmark).* Evaluated
on the four ClusterOverview benchmark datasets (Silhouette
@rousseeuw1987silhouettes, ARI @hubert1985comparing), LPA delivers
*Iris* Sil = $0.283$, ARI = *$0.481$*; *Customers* Sil = $0.430$;
*Multi-Shape* Sil = *$0.594$* (ARI = N/A); and *High-D* Sil = $0.303$
(ARI = N/A), at a mean runtime of *$0.471$ ms*. Remarkably, this
zero-hyperparameter algorithm achieves both the *highest Iris ARI*
($0.481$, outperforming MCL $0.448$, Leiden $0.340$, and Louvain
$0.333$) and the *highest Multi-Shape Silhouette* ($0.594$, well ahead
of MCL $0.534$ and Louvain / Leiden $0.526$), while matching Leiden on
*High-D* ($0.303$) and all three peers on *Customers* ($0.430$) in
under half a millisecond. This strong performance reflects Proposition
II.3 and Proposition V.1 directly: on clean k-NN graphs of distinct
non-convex shapes (*Multi-Shape*) and compact biological species
(*Iris*), the absence of a global modularity null penalty prevents
resolution-limit merges and allows local majority voting to lock onto
true strong-community manifolds — while the lower *Iris* Silhouette
($0.283$) reflects geometrically non-spherical boundary assignments
between overlapping Versicolor and Virginica samples.
