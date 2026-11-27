#import "../template/ieee.typ": *

#alg_entry(
  id: 111,
  name: "Chinese Whispers",
  category: "Cat. V: Graph, Spectral & Network",
  year: "2006",
  paradigm: "Asynchronous Local Majority Label Propagation",
  complexity: $cal(O)(T m)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Chinese Whispers (CW), introduced by Biemann @biemann2006chinese for
large-scale natural language processing graphs, is a parameter-free,
randomized, asynchronous bottom-up graph clustering algorithm. Unlike
spectral or modularity-based methods that require global matrix
factorizations or explicit null-model tuning, CW discovers community
structure purely through repeated local weighted-majority voting over
vertex neighborhoods.

Let $G = (V, E, w)$ be a finite, simple, undirected weighted graph with
$n = abs(V)$ vertices, $m = abs(E)$ edges, and symmetric non-negative
weights $w_(u v) = w_(v u) > 0$ for ${u, v} in E$ ($w_(u v) = 0$
otherwise, with no self-loops: $w_(v v) = 0$). Denote the open
neighborhood of $v in V$ by $N(v) = {u in V : {u, v} in E}$, its
weighted degree by $d_v = sum_(u in N(v)) w_(v u) > 0$, its unweighted
degree by $"deg"_0 (v) = abs(N(v))$, and the total edge weight by
$W = sum_({u, v} in E) w_(u v) = 1/2 sum_(v in V) d_v$. A hard
vertex coloring (cluster assignment) is a map
$c : V -> {1, dots, n}$, inducing the partition
$cal(C)(c) = {C_ell != emptyset : C_ell = c^(-1)(ell)}$.

To suppress hub dominance in scale-free graphs, Biemann
@biemann2006chinese defines a degree-downweighting function
$g : (0, infinity) -> (0, infinity)$ yielding three canonical
weighting schemes:
$
  g_"top" (d) = 1, quad
  g_"lin" (d) = d, quad
  g_"log" (d) = log_2 (d + 1).
$
Given a current coloring $c$, the *class score* of candidate label
$ell in {1, dots, n}$ at vertex $v in V$ is
$
  S_v (ell; c) = sum_(u in N(v)) (w_(v u)) / (g(d_u))
    bb(1)[c_u = ell],
$
where $bb(1)[dot.op]$ is the indicator function. Starting from the
singleton coloring $c_v^((0)) = v$, each epoch $t = 1, dots, T$ draws
a uniformly random permutation $sigma_t in S_n$ and updates vertices
asynchronously in order $v = sigma_t (1), dots, sigma_t (n)$ by
selecting $c_v <- ell^* in arg max_ell S_v (ell; c)$.

#definition(title: "Weighted Majority Equilibrium and Cut Potential")[
  Let $g : (0, infinity) -> (0, infinity)$ be positive. A coloring
  $c : V -> {1, dots, n}$ is a *weak locally stable equilibrium* of
  $g$-weighted Chinese Whispers if, for every $v in V$,
  $
    S_v (c_v; c) >= max_(ell in {1, dots, n}) S_v (ell; c),
  $
  and a *strict locally stable equilibrium* if $c_v$ is the unique
  maximizer of $S_v (dot.op; c)$ for all $v in V$. The *symmetrized
  internal-weight potential* of $c$ is
  $
    Phi_g (c) = sum_({u, v} in E) (w_(u v)) / (g(d_u) g(d_v))
      bb(1)[c_u = c_v].
  $
]

#algbox(
  title: "Chinese Whispers (Asynchronous Weighted Majority)",
  inputs: [Undirected graph $G = (V, E, w)$, weighting function $g$,
    max epochs $T$, tie-breaking mode (lazy or uniform random).],
  outputs: [Partition $cal(C)(c)$ induced by final coloring $c$.],
)[
  1. *Initialize:* Set $c_v <- v$ for all $v in V$; precompute
     $d_v = sum_(u in N(v)) w_(v u)$ and $gamma_v = 1 / g(d_v)$.
  2. *Epoch Loop:* For $t = 1, dots, T$:
     (a) Draw a uniformly random permutation $sigma_t$ of $V$; set
     $"changed" <- "false"$.
     (b) For $r = 1, dots, n$, let $v = sigma_t (r)$:
     compute $S_v (ell) = sum_(u in N(v)) w_(v u) gamma_u
     bb(1)[c_u = ell]$ over active neighbor labels $ell$; let
     $M_v = arg max_ell S_v (ell)$.
     (c) If lazy tie-breaking and $c_v in M_v$, keep $c_v$; else draw
     $ell^*$ uniformly from $M_v$ and set $c_v <- ell^*$. If $c_v$
     changed, set $"changed" <- "true"$.
     (d) If $"changed" = "false"$, terminate early.
  3. *Return:* Non-empty label classes $cal(C)(c)$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

Despite its heuristic origin, Chinese Whispers possesses an exact
game-theoretic and statistical-mechanics structure. We establish five
formal guarantees: a unified Lyapunov convergence theorem covering
all Biemann degree-normalization schemes, constructive counterexamples
for synchronous and directed updates, finite-sample stochastic block
model stability bounds, global-optimum degeneracy, and deterministic
linear-time complexity.

#subsec("A", "Unified Lyapunov Potential Across Weighting Schemes")

Although the `lin` and `log` weight matrices
$tilde(W)_(v u) = w_(v u) / g(d_u)$ are asymmetric on non-regular
graphs, row-scaling reveals an exact symmetric potential game.

#theorem(
  title: "Exact Unified Lyapunov Monotonicity and Finite Convergence",
  origin: "Original Unified Proof",
)[
  Let $G = (V, E, w)$ be a loopless undirected weighted graph and
  $g : (0, infinity) -> (0, infinity)$ be any positive function
  (including `top`, `lin`, and `log` modes). Then:
  (i) Every single asynchronous vertex update $c -> c'$ satisfies
  $Phi_g (c') >= Phi_g (c)$.
  (ii) Under lazy tie-breaking, every label change strictly increases
  $Phi_g$ by at least $delta_"min" > 0$, guaranteeing deterministic
  termination at a weak locally stable equilibrium after at most
  $floor(Phi_g^"max" / delta_"min")$ label flips. In particular, for
  unweighted `top` mode ($w_(u v) = 1$, $g equiv 1$), CW terminates
  after at most $m$ total vertex label changes.
  (iii) Under uniform random tie-breaking, $Phi_g (c^((t)))$
  stabilizes almost surely in finitely many steps, after which all
  subsequent updates preserve $Phi_g$ identically.
]

#proof[
  Suppose vertex $v in V$ updates its label from $a = c_v$ to
  $b = c'_v$ while $c'_u = c_u$ for all $u != v$. Because $G$ is
  simple and loopless ($w_(v v) = 0$), the only edges whose
  same-label indicator $bb(1)[c_u = c_w]$ can change are those
  incident to $v$. Using symmetry $w_(v u) = w_(u v)$,
  $
    Phi_g (c') - Phi_g (c) &= sum_(u in N(v)) (w_(v u)) /
      (g(d_v) g(d_u)) (bb(1)[c_u = b] - bb(1)[c_u = a]) \
      &= (1) / (g(d_v)) (S_v (b; c) - S_v (a; c)).
  $
  Since $b in arg max_ell S_v (ell; c)$, we have
  $S_v (b; c) >= S_v (a; c)$. Because $g(d_v) > 0$, it follows that
  $Phi_g (c') - Phi_g (c) >= 0$, establishing (i).

  For (ii), lazy tie-breaking modifies $c_v$ only when
  $S_v (b; c) > S_v (a; c)$. Because $G$ is finite, the finite set of
  achievable potential values
  $cal(V)_g = {Phi_g (c) : c in {1, dots, n}^V}$ has a strictly
  positive minimum nonzero gap
  $delta_"min" = min{x - y : x, y in cal(V)_g, x > y} > 0$. Every
  label flip therefore obeys $Phi_g (c') - Phi_g (c) >= delta_"min"$.
  Since $0 <= Phi_g (c) <= Phi_g^"max" = sum_({u, v} in E)
  w_(u v) / (g(d_u) g(d_v))$, at most
  $floor(Phi_g^"max" / delta_"min")$ label flips can occur before
  no vertex can strictly improve, at which point $c$ is a weak
  locally stable equilibrium. For unweighted `top` mode,
  $Phi_1 (c) in {0, 1, dots, m}$, so $delta_"min" = 1$ and
  $Phi_1^"max" = m$, giving at most $m$ total label flips.

  For (iii), the sequence $Phi_g (c^((k)))$ along individual updates
  is non-decreasing in the finite set $cal(V)_g$, hence must attain a
  constant value $Phi^*$ after some finite step $k_0$; all steps
  $k >= k_0$ satisfy $S_v (b; c) = S_v (a; c)$ and move only along
  zero-energy plateaus of $Phi_g$.
]

Theorem 1 reveals why `lin` mode empirically prevents mega-clusters:
it replaces unnormalized edge weights $w_(u v)$ by the
conductance-scaled symmetric kernel $w_(u v) / (d_u d_v)$, penalizing
edges joining high-degree hubs quadratically in degree while
preserving exact Lyapunov descent.

#subsec("B", "Failure under Synchronous and Directed Updates")

Both asynchrony and undirected symmetry are strictly necessary for
Lyapunov stability, as shown constructively below.

#proposition(
  title: "Limit Cycles under Synchronous and Directed CW",
  origin: "Original Counterexample",
)[
  (a) On the unweighted bipartite 4-cycle $K_(2, 2)$ with bipartition
  ${a_1, a_2} union {b_1, b_2}$, synchronous CW (`top` mode) enters
  a deterministic period-2 limit cycle with constant $Phi_1 = 0$ from
  the coloring $c(a_1) = c(a_2) = 1$, $c(b_1) = c(b_2) = 2$, whereas
  asynchronous CW reaches the monochromatic fixed point ($Phi_1 = 4$)
  within two vertex updates.
  (b) On the unweighted directed 3-cycle $1 -> 2 -> 3 -> 1$ (where
  each vertex aggregates incoming neighbors $N^-(v)$), asynchronous
  CW under fixed reverse-topological sweep order $(v_3, v_2, v_1)$
  admits a strict 6-step periodic orbit of unique best responses with
  no ties, proving that no strict Lyapunov potential exists on
  general directed graphs.
]

#proof[
  For (a), in $K_(2, 2)$ every $a_i$ has neighborhood ${b_1, b_2}$
  and every $b_j$ has neighborhood ${a_1, a_2}$. Under synchronous
  updates from $(c(a_1), c(a_2), c(b_1), c(b_2)) = (1, 1, 2, 2)$,
  each $a_i$ observes score $2$ for label $2$ (unique maximum) and
  switches to $2$, while simultaneously each $b_j$ observes score $2$
  for label $1$ and switches to $1$, yielding $(2, 2, 1, 1)$ and then
  $(1, 1, 2, 2)$ ad infinitum with $Phi_1 = 0$ (no ties occur). Under
  asynchronous updates, updating $a_1$ first yields $(2, 1, 2, 2)$
  ($Phi_1 = 2$); updating $a_2$ next yields $(2, 2, 2, 2)$
  ($Phi_1 = 4$), an absorbing equilibrium.

  For (b), on $1 -> 2 -> 3 -> 1$ the in-neighborhoods are
  $N^-(1) = {3}$, $N^-(2) = {1}$, $N^-(3) = {2}$, so every vertex
  deterministically copies its unique predecessor with zero ties.
  Starting from $c = (a, b, b)$ with $a != b$ and applying the
  asynchronous update order $(v_2, v_1, v_3)$ repeatedly produces
  $
    (a, b, b) &-> (a, a, b) -> (b, a, b) -> (b, a, a) \
              &-> (b, b, a) -> (a, b, a) -> (a, b, b).
  $
  Every transition is a strict unique best response. Any strict
  Lyapunov function $Psi$ would require
  $Psi(a, b, b) < Psi(a, a, b) < dots < Psi(a, b, b)$, a
  contradiction.
]

#subsec("C", "Planted-Community Stability in Stochastic Block Models")

While CW is not guaranteed to reach any specific local optimum from
arbitrary initializations, we can certify exactly when a ground-truth
partition is an absorbing fixed point.

#theorem(
  title: "Exact Equilibrium Criterion and SBM Recovery Threshold",
  origin: "Original Derivation",
)[
  (a) A partition $cal(C) = {C_1, dots, C_k}$ of $G = (V, E, w)$ is a
  strict locally stable equilibrium of `top`-mode CW if and only if
  every vertex $v in C_r$ satisfies the strong local dominance
  inequality
  $
    d_v^"in" := sum_(u in C_r) w_(v u)
      > max_(s != r) sum_(u in C_s) w_(v u).
  $
  (b) Let $G tilde.op "SBM"(n, k, p, q)$ be a symmetric $k$-block
  stochastic block model @holland1983sbm @abbe2018sbmreview with $k$
  equal blocks of size $s = n / k$, intra-block edge probability $p$,
  and inter-block probability $q < p$. If $s(p - q) > p$, then the
  planted partition fails to be a strict equilibrium of `top` CW
  with probability at most
  $
    "Pr"(cal(E)_"fail") <= n (k - 1)
      exp(- (k (s(p - q) - p)^2) / (n)).
  $
  In particular, whenever
  $p - q >= sqrt((1 + epsilon) k ln(n k) / n) + k p / n$ for any
  $epsilon > 0$, the planted partition is a strict CW fixed point
  with probability $1 - cal(O)((n k)^(-epsilon))$.
]

#proof[
  Part (a) follows immediately from Definition 1 with $g equiv 1$:
  label $r$ uniquely maximizes $S_v (ell; c)$ at $v in C_r$ if and
  only if $S_v (r; c) = d_v^"in" > S_v (s; c) = sum_(u in C_s) w_(v u)$
  for every $s != r$.

  For part (b), fix $v in C_r$ and a competing block $s' != r$. The
  random variables $D_v^"in" = sum_(u in C_r without {v}) A_(v u)$
  and $D_v^((s')) = sum_(u in C_(s')) A_(v u)$ are independent, with
  $D_v^"in" tilde.op "Bin"(s - 1, p)$ and
  $D_v^((s')) tilde.op "Bin"(s, q)$. Write their difference as
  $
    Z_(v, s') = D_v^"in" - D_v^((s'))
      = sum_(i = 1)^(s - 1) X_i + sum_(j = 1)^s (-Y_j),
  $
  a sum of $N = 2 s - 1$ independent random variables, each supported
  on an interval of length $1$ ($X_i in [0, 1]$ and
  $-Y_j in [-1, 0]$). Its expectation is
  $mu_Z = bb(E)[Z_(v, s')] = (s - 1) p - s q = s(p - q) - p > 0$.
  Applying Hoeffding's inequality to $Z_(v, s') - mu_Z <= -mu_Z$,
  $
    "Pr"(D_v^"in" <= D_v^((s')))
      &= "Pr"(Z_(v, s') - mu_Z <= -mu_Z) \
      &<= exp(- (2 mu_Z^2) / (2 s - 1))
      <= exp(- (mu_Z^2) / (s)),
  $
  since $2 s - 1 < 2 s = 2 n / k$. Substituting $s = n / k$ and
  taking a union bound over all $n$ vertices $v in V$ and all $k - 1$
  competing blocks $s' != r$ gives
  $
    "Pr"(cal(E)_"fail") <= n (k - 1)
      exp(- (k (s(p - q) - p)^2) / (n)).
  $
  Setting $s(p - q) - p >= sqrt((1 + epsilon) s ln(n k))$ yields
  $exp(-mu_Z^2 / s) <= (n k)^(-(1 + epsilon))$, so
  $"Pr"(cal(E)_"fail") <= (n k)^(-epsilon) -> 0$.
]

#subsec("D", "Global-Optimum Degeneracy and Equilibria Multiplicity")

Because $Phi_g$ lacks a null-model subtraction term (unlike
modularity @brandes2008modularity), its global optimum is trivial,
and non-trivial clustering relies entirely on early local trapping.

#proposition(
  title: "Monochromatic Global Optimum and Exponential Equilibria",
  origin: "Original Proof",
)[
  (a) On every connected weighted graph $G = (V, E, w)$ with
  $w_(u v) > 0$ on $E$, the unique global maximizer of $Phi_g (c)$
  (up to label relabeling) is the trivial monochromatic coloring
  $c_v equiv ell_0$, attaining $Phi_g^"max"$.
  (b) On the unweighted path graph $P_n$ ($n >= 2$), `top`-mode CW
  admits at least $F_(n - 1) = Theta(phi^n)$ distinct weak locally
  stable partitions (where $F_k$ is the $k$-th Fibonacci number and
  $phi = (1 + sqrt(5)) / 2$). Moreover, on a bridge-connected
  chain of $r$ cliques $K_4$ ($n = 4 r$ vertices), CW admits at least
  $2^(r - 1) = 2^(n / 4 - 1)$ *strict* locally stable partitions.
]

#proof[
  For (a), every term $w_(u v) / (g(d_u) g(d_v))$ in $Phi_g (c)$ is
  strictly positive, so $Phi_g (c) <= Phi_g^"max"$ with equality if
  and only if $c_u = c_v$ across every edge ${u, v} in E$. Because
  $G$ is connected, any two vertices are joined by a path, forcing
  $c_u = c_v$ for all $u, v in V$.

  For (b), partition the path $P_n$ ($1 - 2 - dots - n$) into
  contiguous blocks of lengths $b_1, dots, b_k >= 2$ summing to $n$.
  Every vertex has at least one neighbor in its own block and at most
  one neighbor in a distinct block, so its own label achieves score
  $>= 1 = max_ell S_v (ell; c)$ (weakly stable). The number $a_n$ of
  compositions of $n$ into parts $>= 2$ satisfies $a_2 = 1 = F_1$,
  $a_3 = 1 = F_2$, and $a_n = a_(n - 1) + a_(n - 2)$ (extend the
  final block by $1$, or append a new block of length $2$), hence
  $a_n = F_(n - 1) = Theta(phi^n)$. For the chain of $r$ cliques
  $K_4$ linked by single bridges between distinct portal vertices,
  merging any subset of the $r - 1$ bridge adjacencies yields
  $2^(r - 1)$ partitions in which every vertex has at least $3$
  same-label neighbors and at most $1$ external neighbor ($3 > 1$),
  making all $2^(r - 1)$ partitions strictly locally stable by
  Theorem 2(a).
]

#subsec("E", "Linear-Time Per-Sweep Complexity and Invariances")

#lemma(
  title: "Deterministic Linear Sweep Complexity and Invariances",
  origin: "Classical — Cited",
)[
  Following @biemann2006chinese and @raghavan2007lpa, using an
  epoch-stamped accumulator array of length $n$, each asynchronous CW
  sweep executes in deterministic worst-case time $cal(O)(n + m)$ and
  space $cal(O)(n + m)$. Furthermore, for `top` and `lin` modes,
  the transition law of CW is strictly invariant under uniform
  positive edge-weight scaling $w |-> alpha w$ ($alpha > 0$) and
  graph isomorphism.
]

#proof[
  Maintain two length-$n$ arrays $"acc"[1..n]$ and $"stamp"[1..n]$
  alongside an integer counter $kappa$ incremented at each vertex
  update $v = sigma_t (r)$. When scanning neighbor $u in N(v)$ with
  label $ell = c_u$, if $"stamp"[ell] != kappa$, set
  $"stamp"[ell] <- kappa$ and $"acc"[ell] <- 0$; then increment
  $"acc"[ell]$ by $w_(v u) gamma_u$ and track the running maximum.
  Clearing between vertices requires $cal(O)(1)$ operations, so
  updating $v$ takes $cal(O)(1 + "deg"_0 (v))$ worst-case time.
  Summing over $v in V$ gives
  $sum_(v in V) cal(O)(1 + "deg"_0 (v)) = cal(O)(n + 2 m) = cal(O)(n + m)$
  per sweep, and $cal(O)(T(n + m)) = cal(O)(T m)$ for connected
  graphs ($m >= n - 1$). Storing adjacency lists, labels, degrees,
  and the two accumulator arrays consumes $cal(O)(n + m)$ memory.
  Under $w |-> alpha w$ ($alpha > 0$), $d_v |-> alpha d_v$; in
  `top` mode $S_v (ell)$ scales uniformly by $alpha$, and in
  `lin` mode $w_(v u) / d_u$ is exactly unchanged, preserving
  $arg max_ell S_v (ell)$ identically (whereas `log` mode is not
  scale-invariant because $log_2 (alpha d_u + 1)$ is nonlinear).
  Equivariance under vertex relabeling follows from uniform sampling
  of $sigma_t in S_n$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Chinese Whispers is among the fastest graph clustering algorithms in
existence. Because each sweep touches every edge exactly twice with
no floating-point transcendental calls (in `top` and `lin` modes)
and no priority queue or eigensolver, wall-clock runtime is bounded
by L2 / L3 cache bandwidth during irregular adjacency-list fetches.
In sparse networks ($m = cal(O)(n)$), a typical budget of
$T in [10, 25]$ sweeps completes in linear $cal(O)(n)$ time and can
process tens of millions of edges on a single CPU core.

Geometrically, when applied to Euclidean data via a $k$-nearest
neighbor (kNN) similarity graph, CW acts as an isotropic
zero-temperature diffusion-majority filter on the data manifold. It
does not impose spherical or convex cluster shapes, allowing labels
to propagate around rings, crescents, and curved submanifolds so long
as inter-cluster bridges remain locally outnumbered at boundary
vertices (Theorem 2). However, because the global optimum is
monochromatic (Proposition 2), CW exhibits high run-to-run variance
and sensitivity to early permutation draws: a label that wins early
near a high-degree hub can avalanche across weak bottlenecks before
smaller peripheral communities consolidate.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time complexity],
    [$cal(O)(T(n + m))$; strictly linear per sweep via stamped
      arrays (Lemma 1); typically $T <= 25$ epochs.],
    [Space complexity],
    [$cal(O)(n + m)$ for CSR adjacency lists, label array, degree
      weights, and two $cal(O)(n)$ scratch buffers.],
    [Hyperparameters],
    [None required ($k$ automatic); optional weighting mode
      (`top`, `lin`, `log`) and epoch cap $T$.],
    [Determinism],
    [Stochastic (depends on random permutation $sigma_t$ and tie
      breaking); consensus or multi-restart advised.],
    [Cluster geometry],
    [Arbitrary non-convex manifold shapes via local graph
      connectivity; no Euclidean centroid assumption.],
    [Scale invariance],
    [Exact for `top` and `lin` under $w |-> alpha w$; broken
      for `log` mode.],
    [Resolution profile],
    [No global modularity resolution limit, but prone to giant
      monochromatic collapse across weak bridges.],
  ),
  caption: [Operational and theoretical profile of Chinese Whispers.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

*1. Unsupervised Word-Sense Induction and Lexical Semantics:*
Biemann's original domain @biemann2006chinese constructs word
co-occurrence graphs from billion-token corpora where polysemous
words act as high-degree bridges between topical cliques. Applying CW
to the ego-network of a target word separates senses cleanly because
local triangle density inside each semantic context satisfies
Theorem 2(a), while `lin` or `log` downweighting prevents
stop-words and general hubs from merging distinct senses.

*2. Large-Scale Face and Entity Identity Clustering:*
In computer vision pipelines (e.g., `dlib` face clustering over
millions of embeddings), identities form dense cosine-similarity
kNN cliques with unknown cardinality $K$ (ranging from $10^2$ to
$10^6$). CW's parameter-free automatic $K$ discovery and
$cal(O)(T m)$ complexity allow single-machine clustering of graphs
where hierarchical or spectral algorithms ($cal(O)(n^2)$ to
$cal(O)(n^3)$) are infeasible.

*3. Protein Homology and Sequence Similarity Networks:*
All-versus-all sequence alignment graphs (BLAST / MMseqs2 bit-score
networks) exhibit heavy-tailed degree distributions and wildly
varying family sizes (singletons alongside superfamilies of $10^4$
sequences). Because CW requires no global cluster-size balance prior
and operates strictly via local weighted majority, small isolated
protein families remain intact without suffering modularity's
resolution limit @fortunato2007resolution.

*4. Real-Time Social Network and Fraud-Ring Triage:*
In streaming transaction or interaction graphs with $10^7+$ edges,
`lin`-weighted CW serves as a low-latency first-pass community
extractor. The symmetrized potential $w_(u v) / (d_u d_v)$ naturally
discounts ubiquitous merchant or celebrity hubs while preserving
tightly knit collusive rings whose members share low-degree
counterparties.

*5. Graph Coarsening and Multilevel Preconditioning:*
Because a few sweeps ($T in {2, 3, 5}$) of asynchronous CW rapidly
contract dense local neighborhoods into balanced super-nodes in
strictly linear time, CW is frequently used as an ultra-fast
coarsening heuristic ahead of heavier partitioning or layout solvers
(such as multilevel spectral or modularity optimizers
@blondel2008louvain @traag2019leiden).

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Three mathematical failure regimes govern practical deployment:
(i) *Monster-cluster collapse* (Proposition 2a): in graphs with weak
community contrast (near or below the SBM stability threshold of
Theorem 2b) or dense hubs under `top` mode, a single early label
absorbs a macroscopic fraction of $V$; switching to `lin` mode or
capping $T$ mitigates but does not eliminate collapse.
(ii) *High stochastic variance* (Proposition 2b): on sparse,
low-clustering, or chain-like graphs (grids, road networks, chains),
exponentially many local equilibria yield erratic partitions across
random seeds.
(iii) *Directed and bipartite oscillation* (Proposition 1): omitting
random vertex shuffling or running synchronous updates on bipartite
or directed cyclic graphs triggers non-convergent limit cycles.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [NLP co-occurrence / lexical graphs],
    [High: strong local cliques, unknown $K$, `lin` tames hubs.],
    [Use multi-seed consensus; compare MCL @vandongen2000mcl.],
    [Embedding kNN identity graphs],
    [High: scales to $10^7$ edges with automatic $K$.],
    [Prune weak bridges; Leiden @traag2019leiden if over-merged.],
    [Sparse trees, meshes, road grids],
    [Poor: $Theta(phi^n)$ flat local equilibria (Prop. 2b).],
    [Use spectral @vonluxburg2007tutorial or Walktrap.],
    [Weak-contrast / noisy SBM graphs],
    [Poor: collapses to single giant component.],
    [Use Infomap @rosvall2008infomap or Leiden @traag2019leiden.],
    [Directed flow / citation graphs],
    [Poor: lacks Lyapunov potential (Prop. 1b).],
    [Use Infomap @rosvall2008infomap for directed flow.],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for
    Chinese Whispers.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Chinese Whispers sits directly at the foundation of label-propagation
and zero-temperature spin-model community detection, closely linked
to the Label Propagation Algorithm (LPA) @raghavan2007lpa,
Speaker-Listener LPA @xie2011slpa, Fluid Communities
@pares2018fluidc, and Potts spin-glass clustering
@reichardt2006spinglass.

#proposition(
  title: "Exact Equivalence to Symmetrized LPA and Zero-Temp Potts",
  origin: "Original Derivation",
)[
  Let $G = (V, E, w)$ be a loopless undirected weighted graph and $g$
  a positive degree-weighting function. Then:
  (i) Asynchronous `top`-mode CW ($g equiv 1$) is identically equal
  to asynchronous weighted LPA @raghavan2007lpa initialized with
  unique vertex labels $c_v^((0)) = v$.
  (ii) Asynchronous $g$-weighted CW on $G$ is algebraically
  equivalent to asynchronous `top`-mode LPA executed on the
  degree-symmetrized graph $hat(G)_g = (V, E, hat(w))$ with weights
  $hat(w)_(u v) = w_(u v) / (g(d_u) g(d_v))$.
  (iii) Every sweep of CW performs zero-temperature Glauber dynamics
  (greedy coordinate ascent) on the $n$-state ferromagnetic Potts
  Hamiltonian $cal(H)(c) = -Phi_g (c)$ with resolution parameter
  $gamma = 0$ in the Reichardt–Bornholdt family
  @reichardt2006spinglass.
]

#proof[
  Part (i) holds directly by comparing update rules: both algorithms
  initialize $c_v^((0)) = v$, traverse vertices in a uniformly random
  permutation $sigma_t$, and assign
  $c_v <- arg max_ell sum_(u in N(v)) w_(v u) bb(1)[c_u = ell]$.
  For (ii), for any vertex $v$ and candidate label $ell$,
  $
    hat(S)_v (ell; c) := sum_(u in N(v)) hat(w)_(v u) bb(1)[c_u = ell]
      = (1) / (g(d_v)) S_v (ell; c).
  $
  Because $1 / g(d_v) > 0$ is independent of $ell$, the maximizer set
  $arg max_ell hat(S)_v (ell; c) = arg max_ell S_v (ell; c)$ is
  identical at every update.
  For (iii), the Reichardt–Bornholdt Hamiltonian
  @reichardt2006spinglass on $hat(G)_g$ is
  $cal(H)_gamma (c) = -sum_({u, v} in E) (hat(w)_(u v) - gamma p_(u v))
  bb(1)[c_u = c_v]$; setting $gamma = 0$ yields
  $cal(H)_0 (c) = -Phi_g (c)$, whose single-spin zero-temperature
  Glauber update deterministically selects
  $c_v in arg min_ell cal(H)_0 (c_(-v), ell) = arg max_ell S_v (ell; c)$.
]

*Empirical Benchmark Behavior:* On the ClusterOverview 4-dataset
benchmark suite (constructed via Euclidean kNN similarity graphs
and evaluated with Silhouette @rousseeuw1987silhouettes and Adjusted
Rand Index @hubert1985comparing), Chinese Whispers achieves a mean
runtime of *0.670 ms*—by far the fastest among Algorithms 111–114
(over $11 times$ faster than Walktrap at 7.944 ms, $20 times$ faster
than Girvan–Newman at 13.872 ms, and $142 times$ faster than Infomap
at 95.352 ms), directly confirming the deterministic $cal(O)(T m)$
sweep complexity of Lemma 1. On *Multi-Shape*, CW attains a strong
Silhouette of *0.527* (ARI = N/A), reflecting its ability to
propagate labels smoothly along non-convex manifold neighborhoods
without centroid bias. On *Customers* (Sil = *0.430*) it ties
Infomap and substantially outperforms Girvan–Newman (0.260) and
Walktrap (0.356), while on *High-D* (Sil = *0.303*, ARI =
N/A) it matches the other graph methods closely. Conversely, on
*Iris*, CW obtains a modest Sil = *0.319* and ARI = *0.415*
(trailing Girvan–Newman's ARI = 0.657 and Walktrap's Sil =
0.476): because the Versicolor and Virginica classes overlap in the
kNN graph without a sharp bridge cut, the $gamma = 0$
monochromatic bias of Proposition 2 causes local majority voting to
merge borderline neighborhoods prematurely.
