#import "../template/ieee.typ": *

#alg_entry(
  id: 128,
  name: "DEMON",
  category: "Cat. V: Graph, Spectral & Network",
  year: "2012",
  paradigm: "Local Ego-Network Label Propagation & Merging",
  complexity: $cal(O)(T tau(G) + R M bar(s))$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Global community detection methods struggle on massive, evolving networks
because global objectives such as modularity couple distant regions of the
graph @fortunato2007resolution and ignore the local, egocentric perspective of
individual actors. DEMON (*Democratic Estimate of the Modular Organization of
a Network*), introduced by Coscia et al. @coscia2012demon, adopts a bottom-up
"local-first" paradigm: every vertex independently partitions its immediate
neighborhood into social circles via Label Propagation @raghavan2007lpa, and a
democratic merging step aggregates these local votes into a global overlapping
community cover.

Let $G = (V, E)$ be an undirected simple graph with $abs(V) = n$, $abs(E) = m$,
open neighborhood $N(v) = {u in V : {u, v} in E}$, and degree
$d_v = abs(N(v))$. DEMON takes a merging tolerance $epsilon in [0, 1]$ and a
minimum community size $s_"min" >= 1$.

#definition(title: "Ego-Minus-Ego Subgraph and Overlap-Tolerance Merging")[
  For each vertex $v in V$, the *ego-minus-ego network* is the induced subgraph
  $ G_v^- = G[N(v)] = (N(v), E_v^-), $
  where $E_v^- = {{u, w} in E : u, w in N(v)}$. Removing $v$ eliminates its
  central star edges so that $v$'s distinct local circles decouple cleanly.
  Executing asynchronous Label Propagation (LPA) @raghavan2007lpa on $G_v^-$
  partitions $N(v)$ into clusters ${L_(v, 1), dots, L_(v, r_v)}$; re-inserting
  ego $v$ produces $v$'s *local community proposals*
  $C_(v, j) = L_(v, j) union {v}$. Given two communities $C_1, C_2 subset.eq V$,
  define the *Szymkiewicz--Simpson overlap coefficient*
  $ "ov"(C_1, C_2) = abs(C_1 inter C_2) / min(abs(C_1), abs(C_2)) in [0, 1]. $
  Under Coscia et al.'s discrepancy tolerance $epsilon in [0, 1]$
  @coscia2012demon, $C_1$ and $C_2$ satisfy the *$epsilon$-merge criterion* iff
  at most an $epsilon$-fraction of the smaller set lies outside the larger:
  $ 1 - "ov"(C_1, C_2) <= epsilon <==> "ov"(C_1, C_2) >= 1 - epsilon. $
]

At $epsilon = 0$, two communities merge only under exact containment
($C_1 subset.eq C_2$ or $C_2 subset.eq C_1$); at $epsilon = 1$, every pair
merges into a single component. Intermediate values (typically
$epsilon in [0.2, 0.35]$) permit overlapping boundaries while consolidating
redundant ego proposals voted by multiple members of the same group.

#algbox(
  title: "DEMON (Local-First Democratic Community Discovery)",
  inputs: [Graph $G = (V, E)$, tolerance $epsilon in [0, 1]$, minimum size $s_"min"$, max LPA sweeps $T$],
  outputs: [Overlapping community cover $cal(C) = {C_1, dots, C_M}$],
)[
  + Initialize global community pool $cal(C) <- emptyset$.
  + *Local Ego-Voting Phase* (embarrassingly parallel over $v in V$):
    - Extract ego-minus-ego subgraph $G_v^- = (N(v), E_v^-)$.
    - Assign unique initial labels $ell(u) <- u$ for all $u in N(v)$.
    - For $t = 1, dots, T$, iterate over $u in N(v)$ in random order and update
      $ell(u) <- arg max_c sum_(w in N_(G_v^-) (u)) bb(I)(ell(w) = c)$
      (breaking ties uniformly at random), stopping early upon convergence.
    - For each label class $L_(v, j)$ on $G_v^-$, form proposal
      $C_(v, j) = L_(v, j) union {v}$.
  + *Democratic Merging Phase:* For each proposal $C$ emitted across $v in V$:
    - Search $cal(C)$ for a community $C'$ satisfying
      $"ov"(C, C') >= 1 - epsilon$.
    - If found, replace $C'$ with $C' union C$ (optionally propagating merges
      against remaining entries of $cal(C)$); otherwise insert $C$ into
      $cal(C)$.
  + Filter $cal(C) <- {C in cal(C) : abs(C) >= s_"min"}$ and return $cal(C)$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish four mathematical properties of DEMON: an exact global triangle
identity governing local LPA work and Lyapunov convergence, exact hub
disentanglement under ego removal, an original constructive proof of order
sensitivity and snowball drift in sequential $epsilon$-merging, and Coscia et
al.'s power-law complexity bound.

#subsec("A", "Exact Triangle Accounting and Ego-LPA Potential Monotonicity")

Unlike global LPA, whose per-sweep cost scales with $m$, the aggregate edge
volume across all $n$ ego-minus-ego subgraphs is governed strictly by the
global triangle count $tau(G)$.

#theorem(
  title: "Exact Triangle-Volume Identity and Asynchronous Ego-LPA Convergence",
  origin: "Original Proof",
)[
  Let $tau(G)$ denote the number of triangles ($3$-cliques) in $G = (V, E)$.
  + *(Exact ego-network volume)* The aggregate vertex and edge counts across
    all $n$ ego-minus-ego subgraphs satisfy
    $ sum_(v in V) abs(V(G_v^-)) = 2 m, quad sum_(v in V) abs(E(G_v^-)) = 3 tau(G). $
    Hence one asynchronous LPA sweep across all $n$ ego-networks inspects
    exactly $6 tau(G) + 2 m$ endpoint--edge incidences.
  + *(Exact potential game)* On each $G_v^-$, define the monochromatic edge
    potential $Phi_v (ell) = sum_({u, w} in E_v^-) bb(I)(ell(u) = ell(w))$.
    Under asynchronous single-vertex majority updates (retaining the current
    label whenever it attains the maximum), $Phi_v (ell)$ is monotonically
    non-decreasing and reaches a Nash fixed point in at most $abs(E_v^-)$
    strictly improving label moves.
]

#proof[
  *(1)* By the handshaking lemma,
  $sum_(v in V) abs(V(G_v^-)) = sum_(v in V) d_v = 2 m$. For edges, every
  unordered triangle ${a, b, c}$ in $G$ contributes the edge ${b, c}$ to
  $E_a^-$, ${a, c}$ to $E_b^-$, and ${a, b}$ to $E_c^-$, and contributes to no
  other ego-minus-ego graph. Conversely, every edge ${u, w} in E_v^-$ satisfies
  ${v, u}, {v, w}, {u, w} in E$ and thus uniquely identifies a triangle
  ${v, u, w}$ rooted at $v$. Double counting gives
  $sum_(v in V) abs(E_v^-) = 3 tau(G)$. In one sweep of asynchronous LPA on
  $G_v^-$, visiting each vertex $u in N(v)$ scans its $"deg"_(G_v^-) (u)$
  neighbors, summing to $2 abs(E_v^-) + abs(N(v))$; summing over $v in V$
  yields $6 tau(G) + 2 m$.

  *(2)* Fix $v in V$ and consider updating $u in N(v)$ from $c = ell(u)$ to
  $c' in arg max_gamma n_u (gamma)$, where
  $n_u (gamma) = abs({w in N_(G_v^-) (u) : ell(w) = gamma})$. Because only
  edges incident to $u$ can change monochromatic status,
  $ Phi_v (ell') - Phi_v (ell) = n_u (c') - n_u (c) >= 0. $
  With tie-retention ($c' = c$ whenever $n_u (c) = max_gamma n_u (gamma)$),
  every actual label change increases the integer-valued potential $Phi_v$ by
  at least $1$. Since $0 <= Phi_v (ell) <= abs(E_v^-)$, at most $abs(E_v^-)$
  label changes occur before reaching a fixed point where every vertex holds a
  majority label among its neighbors in $G_v^-$.
]

An immediate corollary of Part (1) is that on any triangle-free graph
($tau(G) = 0$), every ego-minus-ego graph is edgeless ($E_v^- = emptyset$), so
LPA produces only trivial size-$2$ proposals ${u, v}$ and DEMON returns
$emptyset$ whenever $s_"min" >= 3$.

#subsec("B", "Exact Hub-Overlap Disentanglement via Ego Removal")

Removing the ego $v$ is mathematically essential: it severs all length-$2$
cross-circle paths routed through $v$ and forces clique circles to converge
deterministically to distinct labels.

#proposition(
  title: "Deterministic Recovery of Clique Circles Around a Bridge Hub",
  origin: "Original Proof",
)[
  Let $v in V$ belong to $q >= 2$ complete subgraphs $K^((1)), dots, K^((q))$
  of orders $s_j = abs(K^((j))) >= 3$ intersecting pairwise only at
  $K^((i)) inter K^((j)) = {v}$, with no edges between
  $K^((i)) without {v}$ and $K^((j)) without {v}$ for $i != j$, and
  $N(v) = union.big_(j=1)^q (K^((j)) without {v})$. Then:
  + Every fixed point of asynchronous LPA on a clique $K_r$ ($r >= 2$) is
    monochromatic.
  + Asynchronous LPA on $G_v^-$ terminates almost surely at the exact $q$-way
    partition ${K^((1)) without {v}, dots, K^((q)) without {v}}$, so $v$'s
    proposals recover all $q$ overlapping circles ${K^((1)), dots, K^((q))}$
    exactly.
]

#proof[
  *(1)* Let $ell$ be any non-monochromatic labeling of a complete graph $K_r$
  ($r >= 2$). Let $c^*$ be a label achieving maximum multiplicity
  $m_1 = max_gamma abs({x in K_r : ell(x) = gamma}) < r$. Because $ell$ is not
  monochromatic, there exists a vertex $u in K_r$ with label $c = ell(u) != c^*$
  and class size $m_c <= m_1$. In $K_r$, $u$ is adjacent to all other $r - 1$
  vertices; hence $u$ observes $n_u (c^*) = m_1$ neighbors labeled $c^*$,
  whereas $u$ observes only $n_u (c) = m_c - 1 <= m_1 - 1 < m_1$ neighbors
  labeled $c$. Therefore $u$ is not at a majority label and strictly increases
  $Phi$ upon updating. Thus only monochromatic labelings are fixed points on
  $K_r$.

  *(2)* Under the stated hypotheses, removing $v$ produces the disjoint union
  of cliques $G_v^- = union.sq_(j=1)^q K_(s_j - 1)$ with $s_j - 1 >= 2$. Since
  LPA propagates labels only along edges, distinct connected components never
  share labels (each starts from disjoint vertex IDs). Within each component
  $K^((j)) without {v}$, random asynchronous updates break ties with positive
  probability and reach a unique monochromatic state by Part (1) and Theorem 1.
  Re-inserting $v$ yields the exact proposal family
  ${K^((1)), dots, K^((q))}$.
]

Had $v$ been retained in the closed ego graph $G[N[v]]$, a dominant circle
$K^((1))$ of large order $s_1$ would capture $v$'s label and bleed across $v$
into smaller circles; ego removal guarantees exact overlap preservation.

#subsec("C", "Order Sensitivity and Snowball Drift of Sequential Merging")

While local ego-voting is embarrassingly parallel and permutation-covariant,
the global democratic merging phase breaks associativity and confluence.

#theorem(
  title: "Non-Confluence and Snowball Drift of Greedy epsilon-Merging",
  origin: "Original Counterexample",
)[
  For every merging tolerance $epsilon in (0, 1)$:
  + *(Non-transitivity and order dependence)* There exists a concrete triple of
    local proposals ${C_1, C_2, C_3}$ such that even under exhaustive
    re-merging until no mergeable pair remains, processing in order
    $(C_2, C_1, C_3)$ yields a strictly different stable cover than processing
    in order $(C_3, C_1, C_2)$.
  + *(Snowball percolation drift)* For every integer $L >= 2$, there exists a
    sequence of equal-sized proposals $C_0, C_1, dots, C_(L-1)$ with
    $C_0 inter C_(L-1) = emptyset$ (zero direct overlap) that greedy
    $epsilon$-merging collapses into a single community of unbounded diameter
    $Theta(L)$.
]

#proof[
  *(1)* Fix $epsilon in (0, 1)$. Choose integers
  $s = ceil(1 / epsilon)$ and $L = ceil((s + 1) / (1 - epsilon)) + 1$, so that
  $ s / (s + 1) = 1 - 1 / (s + 1) > 1 - epsilon, quad s / L < 1 - epsilon. $
  Construct two sets $C_2, C_3$ of cardinality $L$ sharing an intersection
  $S = C_2 inter C_3$ of size $abs(S) = s$, and pick a fresh vertex
  $x in.not C_2 union C_3$. Set $C_1 = S union {x}$, so
  $abs(C_1) = s + 1 < L$. Checking the pairwise overlap coefficients gives:
  $ "ov"(C_1, C_2) = "ov"(C_1, C_3) = s / (s + 1) > 1 - epsilon, $
  so $C_1$ is mergeable with both $C_2$ and $C_3$, whereas
  $ "ov"(C_2 union {x}, C_3) = "ov"(C_3 union {x}, C_2) = s / L < 1 - epsilon. $
  Under order $(C_2, C_1, C_3)$, $C_1$ merges into $C_2$ forming
  $C_2 union {x}$, which cannot merge with $C_3$, terminating at the stable
  cover ${C_2 union {x}, C_3}$. Under order $(C_3, C_1, C_2)$, $C_1$ merges into
  $C_3$, terminating at the distinct stable cover ${C_2, C_3 union {x}}$ (the
  vertex $x$ flips membership solely due to arrival order).

  *(2)* Let $s = ceil(1 / epsilon)$ and step $r = floor(epsilon s) >= 1$ (note
  $s - r >= (1 - epsilon) s$). For $j = 0, 1, dots, L - 1$, define integer
  intervals $C_j = {j r + 1, dots, j r + s}$ of size $s$. At step $j >= 1$,
  let $U_(j-1) = union.big_(i=0)^(j-1) C_i$. Then
  $ abs(U_(j-1) inter C_j) >= abs(C_(j-1) inter C_j) = s - r, $
  while $min(abs(U_(j-1)), abs(C_j)) = s$, yielding
  $"ov"(U_(j-1), C_j) >= (s - r) / s >= 1 - epsilon$. By induction, all $L$
  proposals collapse into $U_(L-1)$ of size $s + (L - 1) r = Theta(L)$, even
  though $C_0 inter C_(L-1) = emptyset$ whenever $(L - 1) r >= s$.
]

#subsec("D", "Sub-Quadratic Scaling on Scale-Free Networks")

In realistic heavy-tailed graphs, local ego-voting remains tractable despite
high-degree hubs, as established by Coscia et al. @coscia2012demon.

#theorem(
  title: "Power-Law Ego-LPA Complexity Bound",
  origin: "Classical — Cited",
)[
  Suppose $G$ has $n$ vertices with empirical power-law degree distribution
  $P(d) = c d^(-alpha)$ supported on $d in {1, dots, K}$ with maximum degree
  $K = d_"max"$ and exponent $alpha in (2, 3)$ @coscia2012demon. Assuming
  $T = cal(O)(1)$ LPA sweeps per ego-network, the total computational work of
  the ego-network extraction and local LPA stage satisfies
  $ W_"ego" (n, K) = cal(O)(n K^(3 - alpha)). $
  Under the natural scale-free hub cutoff $K = Theta(n^(1 / (alpha - 1)))$,
  this simplifies to $W_"ego" = cal(O)(n^(2 / (alpha - 1))) = o(n^2)$, and for
  $alpha > 3$ it reduces to linear time $cal(O)(n)$.
]

#proof[
  For any vertex $v$ of degree $d$, its ego-minus-ego subgraph $G_v^-$ has at
  most $binom(d, 2) <= d^2 / 2$ edges. Constructing $G_v^-$ via hash-indexed
  adjacency membership and running $T$ sweeps of LPA costs at most
  $c_0 T d^2$ operations for an absolute constant $c_0$. Summing over the
  degree distribution gives
  $ W_"ego" (n, K) <= c_0 T n sum_(d=1)^K c d^(-alpha) d^2 = c_1 T n sum_(d=1)^K d^(2 - alpha). $
  Since $2 - alpha > -1$ for $alpha in (2, 3)$, bounding the sum by its integral
  yields
  $ sum_(d=1)^K d^(2 - alpha) <= 1 + integral_1^K x^(2 - alpha) dif x = 1 + (K^(3 - alpha) - 1) / (3 - alpha) = cal(O)(K^(3 - alpha)). $
  Substituting the natural cutoff $K = Theta(n^(1 / (alpha - 1)))$ gives
  exponent $1 + (3 - alpha) / (alpha - 1) = 2 / (alpha - 1) in (1, 2)$,
  strictly sub-quadratic in $n$. When $alpha > 3$, the series converges to the
  Riemann zeta constant $zeta(alpha - 2) = cal(O)(1)$, yielding
  $W_"ego" = cal(O)(n)$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

DEMON separates into two sharply contrasting computational stages. Stage 1
(ego-minus-ego extraction plus asynchronous LPA) requires
$cal(O)(T tau(G) + m)$ edge inspections by Theorem 1 (worst-case
$cal(O)(n d_"max"^2)$ on dense graphs, or $cal(O)(n K^(3 - alpha))$ on
scale-free networks by Section II-D). Because each ego-network $G_v^-$ is
independent, Stage 1 is *embarrassingly parallel* and natively supports
incremental streaming updates: inserting or deleting an edge ${u, v}$ modifies
only the ego-networks of ${u, v} union (N(u) inter N(v))$ @coscia2012demon.

Stage 2 (democratic merging) often dominates wall-clock runtime on graphs with
many raw proposals. Let $R = sum_(v in V) r_v <= 2 m$ be the number of raw
local proposals, $M$ the final community count, and $bar(s)$ the average
community size. Naive pairwise merge checking against the growing pool costs
$cal(O)(R M bar(s))$; inverted vertex-to-community indices restrict overlap
evaluations to pool communities sharing at least one vertex with $C$,
substantially accelerating sparse regimes. Space complexity is
$cal(O)(n + m + M bar(s))$ since ego-networks can be processed and discarded on
the fly without materializing all $R$ proposals simultaneously.

Geometrically, every raw local proposal has graph diameter at most $2$ (a star
centered at ego $v$), imparting a strong local-cohesion bias before merging
stitches adjacent neighborhoods into mesoscale structures. Due to random LPA
tie-breaking and order-sensitive merging (Section II-C), DEMON is *stochastic*
unless vertex iteration and merge orders are canonically seeded.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Stage 1 (Ego-LPA)], [$cal(O)(T tau(G) + m)$; $cal(O)(n K^(3 - alpha))$ scale-free],
    [Stage 2 (Merging)], [$cal(O)(R M bar(s))$ naive; accelerated via inverted index],
    [Space complexity], [$cal(O)(n + m + M bar(s))$ (streaming ego graphs)],
    [Partition type], [Soft overlapping cover; drops small groups $< s_"min"$],
    [Determinism], [Stochastic (LPA tie-breaking + merge order sensitivity)],
    [Dynamic updates], [Local recomputation restricted to $N(u) inter N(v)$],
    [Hyperparameters], [Merge tolerance $epsilon in [0, 1]$, min size $s_"min"$],
    [Parallelism], [Stage 1 fully parallel; Stage 2 sequential / sharded],
  ),
  caption: [Operational and complexity profile of DEMON (Alg. 128).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

DEMON's egocentric formulation aligns naturally with networks driven by
overlapping local circles and frequent edge mutations:

+ *Online Social Networks and Multi-Circle Profiling @coscia2012demon:* Users
  simultaneously belong to family, workplace, and hobby groups whose internal
  density is high but whose cross-circle connections pass primarily through the
  user. Ego removal (Section II-B) exposes each circle cleanly without a
  global modularity sweep over billions of edges.
+ *Dynamic E-Commerce Co-Purchase Streams @coscia2012demon:* In Amazon-style
  product co-purchasing graphs, multi-purpose items bridge distinct consumer
  niches. Because updating an edge only perturbs triangles incident to that
  edge, DEMON incrementally refreshes affected niches in real time.
+ *Actor Collaboration and Entertainment Knowledge Graphs @coscia2012demon:*
  Performers working across genres (e.g., IMDb collaboration networks) form
  overlapping, triangle-rich ego neighborhoods where strict cliques (as
  required by CPM) are overly rigid due to incomplete cast pairings, whereas
  local LPA accommodates missing intra-group edges.
+ *Decentralized Peer-to-Peer and Federated Graph Analytics:* Because each node
  requires only its $2$-hop neighborhood $G[N(v)]$ to cast local community
  votes, Stage 1 executes natively in distributed vertex-centric frameworks
  without centralizing the global adjacency matrix.
+ *Fraud Ring Detection in Payment and Telecom Networks:* Collusive accounts
  frequently reuse mule bridges across multiple independent scam cells. Local
  ego-voting isolates each fraud cell around shared mule accounts while
  avoiding the global resolution-limit washout of small rings inside massive
  transaction ledgers @fortunato2007resolution.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Four failure regimes limit DEMON: (i) *triangle-sparse or bipartite graphs*
($tau(G) approx 0$), where Theorem 1 implies empty ego-minus-ego edge sets;
(ii) *small $epsilon$ vs. large $epsilon$ phase collapse*, where excessive
tolerance triggers Section II-C snowball drift into giant amorphous blobs, while
overly strict $epsilon -> 0$ leaves thousands of near-duplicate ego fragments;
(iii) *high-dimensional $k$-NN hubness* @beyer1999nearest, where hub vertices
acquire huge heterogeneous ego-networks that mix unrelated manifolds; and
(iv) *reproducibility-critical pipelines*, where merge-order sensitivity causes
run-to-run variance unless canonical sorting and consensus runs are applied.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Social networks], [High: multi-circle hubs + high clustering], [Tune $epsilon$ against snowballing; try SLPA @xie2011slpa],
    [Dynamic co-purchase], [High: localized triangle-only updates], [Use inverted merge index; try BIGCLAM @yang2013bigclam],
    [Collaboration graphs], [High: relaxes rigid $k$-cliques via LPA], [Seed canonical merge order; compare CPM (Alg. 127)],
    [Bipartite / grid nets], [Poor: $tau(G) = 0 => E_v^- = emptyset$], [Use SBM (Alg. 130) or Louvain (Alg. 129)],
    [Dense $k$-NN manifolds], [Low--Moderate: radius-$1$ ego fragments], [Use Fast Unfolding (Alg. 129) or SBM (Alg. 130)],
  ),
  caption: [Domain suitability and failure-mode matrix for DEMON (Alg. 128).],
)

#sec("V", "Structural Relations and Empirical Behavior")

DEMON bridges triangle-based overlapping detection (CPM, Alg. 127
@palla2005cpm) and scalable label-propagation dynamics (global LPA
@raghavan2007lpa and SLPA @xie2011slpa).

#proposition(
  title: "Triangle-Subgraph Support Equivalence to 3-CPM and LPA Contrast",
  origin: "Original Derivation",
)[
  Let $G_Delta = (V_Delta, E_Delta)$ denote the subgraph of $G$ formed by all
  edges participating in at least one triangle. For $s_"min" >= 3$:
  + Every non-trivial local proposal $C_(v, j)$ retained by DEMON before
    merging is contained in $V_Delta$ and induces a subgraph of diameter at
    most $2$ in $G_Delta$. Thus DEMON shares the exact vertex support
    $V_Delta$ of $3$-CPM (Alg. 127), while strictly relaxing $3$-clique
    adjacency to local majority connectivity on $G_v^-$.
  + Unlike global LPA @raghavan2007lpa, which outputs a hard partition prone to
    monster-community collapse, localizing LPA to ${G_v^-}$ caps raw cluster
    radius at $1$ from ego $v$ and permits vertex overlap multiplicity up to
    $floor(d_v / 2)$ non-trivial circles.
]

#proof[
  *(1)* Let $C_(v, j) = L_(v, j) union {v}$ satisfy $abs(C_(v, j)) >= 3$, so
  $abs(L_(v, j)) >= 2$. Since asynchronous LPA initializes distinct labels on
  $N(v)$ and only copies labels across edges of $E_v^-$, every vertex
  $u in L_(v, j)$ must have at least one neighbor $w in L_(v, j)$ with
  ${u, w} in E_v^-$. Then ${v, u, w}$ is a triangle in $G$, placing $v, u, w$
  and all edges ${v, u}, {v, w}, {u, w}$ inside $G_Delta$. Because $v$ is
  adjacent in $G_Delta$ to every $u in L_(v, j)$, the induced diameter of
  $C_(v, j)$ is at most $2$. Since $3$-CPM covers precisely the vertices
  incident to at least one triangle, both algorithms operate on the identical
  support $V_Delta$.

  *(2)* Every non-trivial local cluster $L_(v, j)$ has cardinality at least
  $2$ and partitions $N(v)$ disjointly, so at most $floor(d_v / 2)$ non-trivial
  proposals can contain $v$ from its own ego vote, yielding overlapping
  membership bounded away from global monster collapse prior to merging.
]

Unlike modularity optimization (Fast Unfolding, Alg. 129 @blondel2008louvain)
or Girvan--Newman edge betweenness @girvan2002community, DEMON optimizes no
global partition score and instead trades global coherence for local
overlapping sensitivity and dynamic updatability (see also OSLOM
@lancichinetti2011oslom).

*Empirical Benchmark Analysis (ClusterOverview 4-Dataset Suite).* Evaluating
DEMON's induced clustering on $k$-NN proximity graphs across the four
ClusterOverview datasets (Silhouette @rousseeuw1987silhouettes and ARI
@hubert1985comparing):
- *Iris (`Sil = 0.157`, `ARI = 0.498`):* DEMON exhibits substantially lower
  Silhouette than CPM and Fast Unfolding (`0.586`) and SBM (`0.445`), with
  moderate ARI (`0.498`). On compact $k$-NN graphs of continuous Euclidean data,
  democratic $epsilon$-merging suffers from Section II-C snowball drift across
  the Versicolor--Virginica overlap zone, degrading hard-cluster margin
  separation even while partially preserving species grouping.
- *Customers (`Sil = 0.250`):* DEMON outperforms Fast Unfolding (`0.190`) and
  nearly matches CPM (`0.262`), validating that local ego-voting captures fine
  retail sub-cohorts without modularity's resolution-limit over-aggregation,
  though lagging behind SBM (`0.356`).
- *Multi-Shape (`Sil = 0.396`, `ARI = N/A`):* Lowest among the four network
  algorithms (trailing CPM `0.554` and Louvain / SBM `0.713`), reflecting the
  geometric mismatch between radius-$1$ ego-star proposals and elongated
  non-convex geometric manifolds.
- *High-D (`Sil = 0.163`, `ARI = N/A`):* Weakest Silhouette in the group
  (below CPM / Louvain `0.188` and SBM `0.303`), as high-dimensional distance
  concentration @beyer1999nearest inflates hub ego-networks while suppressing
  local triangle transitivity $tau(G)$.
- *Runtime (`Mean = 4.424 ms`):* Runs $29.5 times$ faster than SBM
  (`130.533 ms`), $1.6 times$ slower than CPM (`2.800 ms`), and $6.5 times$
  slower than Fast Unfolding (`0.681 ms`) due to $n$ repeated ego-LPA passes
  and sequential overlap checks.
