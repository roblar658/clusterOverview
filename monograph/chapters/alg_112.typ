#import "../template/ieee.typ": *

#alg_entry(
  id: 112,
  name: "Girvan-Newman",
  category: "Cat. V: Graph, Spectral & Network",
  year: "2002",
  paradigm: "Divisive Edge-Betweenness Hierarchical Partitioning",
  complexity: $cal(O)(m^2 n)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

The Girvan–Newman (GN) algorithm @girvan2002community
@newman2004finding inaugurated modern network community detection by
inverting classical agglomerative clustering: rather than merging
densely connected vertex pairs, GN operates *divisively*, iteratively
identifying and deleting inter-community bottleneck edges until the
graph fragments into a complete hierarchy (dendrogram) of nested
connected components.

Let $G_0 = (V, E_0)$ be a finite, simple, undirected unweighted graph
with $n = abs(V)$ vertices, $m = abs(E_0)$ edges, adjacency matrix
$A in {0, 1}^(n times n)$, and vertex degrees
$d_v = sum_(u in V) A_(v u)$. For any intermediate subgraph
$G = (V, E)$ ($E subset.eq E_0$) and any pair of vertices
$s, t in V$, let $"dist"_G (s, t)$ denote the shortest-path hop
distance (with $"dist"_G (s, t) = infinity$ if $s, t$ lie in distinct
connected components). Let $sigma_(s t)$ be the number of distinct
shortest paths between $s$ and $t$ in $G$ (with $sigma_(s s) = 1$ and
$sigma_(s t) = 0$ if disconnected), and for any edge $e in E$, let
$sigma_(s t) (e)$ be the number of shortest $s$–$t$ paths traversing
$e$.

#definition(title: "Edge Betweenness Centrality and Modularity")[
  The *shortest-path edge betweenness centrality* of edge $e in E$ in
  subgraph $G = (V, E)$ is the aggregate fraction of shortest paths
  passing through $e$ across all unordered connected vertex pairs:
  $
    C_B (e; G) = sum_({s, t} subset.eq V : sigma_(s t) > 0)
      (sigma_(s t) (e)) / (sigma_(s t))
      = 1/2 sum_(s in V) sum_(t in V without {s})
      (sigma_(s t) (e)) / (sigma_(s t)).
  $
  Given a candidate partition $cal(C) = {C_1, dots, C_k}$ of $V$
  produced during edge deletion, its *Newman–Girvan modularity*
  @newman2004finding evaluated on the original graph $G_0$ is
  $
    Q(cal(C)) = (1) / (2 m) sum_(u, v in V)
      (A_(u v) - (d_u d_v) / (2 m)) bb(1)[c_u = c_v]
      = sum_(r = 1)^k ( (e_(r r)) / (m) - ((a_r) / (2 m))^2 ),
  $
  where $e_(r r) = abs({{u, v} in E_0 : u, v in C_r})$ is the number
  of original edges inside $C_r$ and $a_r = sum_(v in C_r) d_v$ is
  its total original degree.
]

The core algorithmic principle of GN is *dynamic recomputation*:
after removing an edge $e^* in arg max_(e in E) C_B (e; G)$,
shortest paths reroute across surviving bridges, so $C_B$ must be
recomputed on the affected connected component before selecting the
next edge. Whenever edge deletion splits a component, a new level is
appended to the divisive dendrogram; the final flat partition is
selected by maximizing $Q(cal(C))$ over all dendrogram cuts.

#algbox(
  title: "Girvan-Newman Divisive Community Detection",
  inputs: [Undirected graph $G_0 = (V, E_0)$ with $n = abs(V)$ and
    $m = abs(E_0)$.],
  outputs: [Divisive dendrogram $cal(T)$ and modularity-optimal
    partition $cal(C)^*$.],
)[
  1. *Initialize:* Set $G <- G_0$, record initial connected-component
     partition $cal(C)_0$, and set $Q^* <- Q(cal(C)_0)$,
     $cal(C)^* <- cal(C)_0$.
  2. *Edge-Removal Loop:* While $E != emptyset$:
     (a) For every source $s in V$ (within each component altered in
     the prior step), execute BFS from $s$ to obtain distances
     $ell_s (v)$, path counts $sigma_(s v)$, and predecessor lists
     $P_s (w) = {v in N(w) : ell_s (v) + 1 = ell_s (w)}$.
     (b) Traverse vertices $w$ in non-increasing order of $ell_s (w)$
     and accumulate directed edge dependencies
     $delta_(s dot.op) (v, w) = (sigma_(s v) / sigma_(s w))
     (1 + sum_(u : w in P_s (u)) delta_(s dot.op) (w, u))$.
     (c) Sum $C_B (e) <- 1/2 sum_(s in V) delta_(s dot.op) (e)$ over
     all $e in E$.
     (d) Remove $e^* in arg max_(e in E) C_B (e)$ from $E$ (breaking
     ties deterministically or uniformly).
     (e) If removing $e^*$ increases the number of connected
     components, record the new partition $cal(C)_t$ in $cal(T)$ and
     update $(Q^*, cal(C)^*)$ if $Q(cal(C)_t) > Q^*$.
  3. *Return:* Dendrogram $cal(T)$ and optimal partition $cal(C)^*$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five formal results characterizing Girvan–Newman: exact
correctness of the $cal(O)(n m)$ dependency recurrence, an exact
conservation law and barbell bridge recovery theorem, a constructive
double-windmill proof that static betweenness fails without
recomputation, the exact modularity resolution limit on the GN
dendrogram, and tight best-case versus worst-case complexity bounds.

#subsec("A", "Exactness of the Shortest-Path Dependency Recurrence")

#theorem(
  title: "Exactness and Complexity of Edge-Betweenness Accumulation",
  origin: "Classical — Cited",
)[
  Following @girvan2002community and @newman2004finding, let
  $G = (V, E)$ be an unweighted undirected graph. For a fixed source
  $s in V$, define the single-source edge dependency of a shortest-path
  DAG edge $(v, w)$ (where $v in P_s (w)$) by
  $delta_(s dot.op) (v, w) = sum_(t in V) sigma_(s t) (v, w) / sigma_(s t)$.
  Then $delta_(s dot.op)$ satisfies the exact reverse-BFS recurrence
  $
    delta_(s dot.op) (v, w) = (sigma_(s v)) / (sigma_(s w))
      ( 1 + sum_(u : w in P_s (u)) delta_(s dot.op) (w, u) ),
  $
  and $C_B (e; G) = 1/2 sum_(s in V) delta_(s dot.op) (e)$ can be
  computed for all $e in E$ simultaneously in $cal(O)(n(n + m))$
  time ($cal(O)(n m)$ when $G$ is connected).
]

#proof[
  Fix $s in V$ and a directed edge $(v, w)$ on the shortest-path DAG
  rooted at $s$, so $ell_s (v) + 1 = ell_s (w)$. Partition the sum
  defining $delta_(s dot.op) (v, w)$ into $t = w$ and $t != w$. For
  $t = w$, every shortest $s$–$w$ path traversing $(v, w)$ consists
  of an arbitrary shortest $s$–$v$ path followed by $(v, w)$, so
  $sigma_(s w) (v, w) = sigma_(s v)$, contributing
  $sigma_(s v) / sigma_(s w)$.

  For $t != w$, any shortest $s$–$t$ path passing through $(v, w)$
  must leave $w$ along a unique outgoing DAG edge $(w, u)$ with
  $w in P_s (u)$ (since $ell_s$ increases strictly by $1$ at every
  hop). Denote by $sigma_(s t) (v, w, u)$ the number of shortest
  $s$–$t$ paths using both $(v, w)$ and $(w, u)$, and by
  $tau_(u t)$ the number of shortest $u$–$t$ suffixes of length
  $ell_s (t) - ell_s (u)$. Because prefixes and suffixes concatenate
  independently along shortest paths,
  $
    sigma_(s t) (v, w, u) = sigma_(s v) dot.op tau_(u t)
      = (sigma_(s v)) / (sigma_(s w)) (sigma_(s w) dot.op tau_(u t))
      = (sigma_(s v)) / (sigma_(s w)) sigma_(s t) (w, u).
  $
  Dividing by $sigma_(s t)$ and summing over all $t != w$ and all
  out-neighbors $u$ with $w in P_s (u)$ yields
  $
    delta_(s dot.op) (v, w) &= (sigma_(s v)) / (sigma_(s w))
      + sum_(u : w in P_s (u)) (sigma_(s v)) / (sigma_(s w))
      sum_(t in V) (sigma_(s t) (w, u)) / (sigma_(s t)) \
      &= (sigma_(s v)) / (sigma_(s w))
      ( 1 + sum_(u : w in P_s (u)) delta_(s dot.op) (w, u) ).
  $
  Since $ell_s (u) = ell_s (w) + 1 > ell_s (w)$, processing vertices
  in non-increasing order of $ell_s$ ensures all child dependencies
  $delta_(s dot.op) (w, u)$ are finalized before evaluating
  $delta_(s dot.op) (v, w)$. Each BFS visits each vertex once and
  each undirected edge twice ($cal(O)(n + m)$), and the reverse pass
  traverses each DAG edge once ($cal(O)(m)$). Repeating over all
  $s in V$ and dividing by $2$ (each unordered pair ${s, t}$ is
  counted from both endpoints) computes all $C_B (e; G)$ in
  $cal(O)(n(n + m))$ time.
]

#subsec("B", "Conservation Identity and Exact Barbell Bridge Recovery")

#theorem(
  title: "Betweenness Conservation Law and Exact Barbell Separation",
  origin: "Original Proof",
)[
  (a) For any connected unweighted graph $G = (V, E)$, the sum of
  edge betweenness centralities equals the Wiener index (total
  pairwise shortest-path distance):
  $
    sum_(e in E) C_B (e; G)
      = sum_({s, t} subset.eq V) "dist"_G (s, t)
      = binom(n, 2) bar(ell)(G).
  $
  (b) Let $G$ be the barbell graph consisting of two disjoint cliques
  $A = K_(n_1)$ and $B = K_(n_2)$ ($n_1, n_2 >= 2$) joined by a
  single bridge $e^* = {a, b}$ ($a in A$, $b in B$). Then the exact
  edge betweenness values are $C_B (e^*) = n_1 n_2$,
  $C_B ({a, u}) = 1 + n_2$ ($u in A without {a}$),
  $C_B ({b, v}) = 1 + n_1$ ($v in B without {b}$), and
  $C_B (e_"int") = 1$ for every interior clique edge. Consequently,
  $e^*$ is the strict unique maximizer of $C_B$, and GN recovers
  ${A, B}$ on Iteration 1.
]

#proof[
  For (a), fix any unordered pair ${s, t} subset.eq V$. Every
  shortest $s$–$t$ path $P$ traverses exactly $"dist"_G (s, t)$
  edges, so summing the indicator $bb(1)[e in P]$ over all $e in E$
  and all $sigma_(s t)$ shortest paths gives
  $sum_(e in E) sigma_(s t) (e) = "dist"_G (s, t) dot.op sigma_(s t)$.
  Dividing by $sigma_(s t)$ and summing over all ${s, t} subset.eq V$
  yields
  $
    sum_(e in E) C_B (e; G)
      = sum_({s, t} subset.eq V) sum_(e in E)
      (sigma_(s t) (e)) / (sigma_(s t))
      = sum_({s, t} subset.eq V) "dist"_G (s, t).
  $

  For (b), every cross-pair ${x, y} in A times B$ has a unique
  shortest path passing through $e^*$, while no intra-clique pair
  uses $e^*$ (cliques have diameter $1$), hence $C_B (e^*) = n_1 n_2$.
  For a portal edge ${a, u}$ ($u in A without {a}$), the shortest
  path from $u$ to $a$ is ${a, u}$ ($1$ pair) and the unique shortest
  path from $u$ to any $v in B$ is $u - a - b$ (for $v = b$) or
  $u - a - b - v$ (for $v in B without {b}$), contributing $n_2$
  cross-pairs; no pair not containing $u$ uses ${a, u}$ because $a$
  is directly adjacent to every $u' in A$. Thus
  $C_B ({a, u}) = 1 + n_2$, and by symmetry $C_B ({b, v}) = 1 + n_1$.
  For an interior edge ${u, u'}$ ($u, u' in A without {a}$), only the
  pair ${u, u'}$ itself uses it, giving $C_B ({u, u'}) = 1$. Since
  $n_1, n_2 >= 2$,
  $
    C_B (e^*) - C_B ({a, u})
      = n_1 n_2 - (1 + n_2) = (n_1 - 1) n_2 - 1 >= 1 > 0,
  $
  and symmetrically $C_B (e^*) > C_B ({b, v}) > 1$. Hence $e^*$ is
  uniquely selected and deleted at step 1, separating $A$ and $B$.
]

#subsec("C", "Necessity of Dynamic Recomputation (Static Failure)")

A classic temptation is to compute $C_B (e; G_0)$ once and delete
edges in static descending order in $cal(O)(n m)$ total time. The
next construction proves that static betweenness fails catastrophically
on elementary two-community graphs.

#proposition(
  title: "Double-Windmill Counterexample to Static Betweenness",
  origin: "Original Counterexample",
)[
  For any integer $r >= 2$, let $A$ and $B$ be identical windmill
  communities of $r$ triangles sharing hubs $h_A$ and $h_B$
  ($abs(A) = abs(B) = 2 r + 1$), joined by two inter-community
  bridges: hub-to-hub $e_1 = {h_A, h_B}$ and leaf-to-leaf
  $e_2 = {u_1, v_1}$ ($u_1 in A without {h_A}$,
  $v_1 in B without {h_B}$). Then:
  (i) Initially in $G_0$,
  $C_B (e_1) = 4 r^2 + 2 r - 3/2$, every spoke edge
  $e_"sp" = {h_A, u_2}$ outside $u_1$'s triangle has
  $C_B (e_"sp") = 4 r$, and the second bridge has
  $C_B (e_2) = 2 r + 5/2 < 4 r$. Thus static betweenness ranks all
  $4 r - 4$ internal spoke edges strictly above $e_2$ and shatters
  $A$ and $B$ before cutting $e_2$.
  (ii) After GN deletes $e_1$ and recomputes on $G_1 = G_0 without {e_1}$,
  $C_B (e_2; G_1) = (2 r + 1)^2 = 4 r^2 + 4 r + 1$ becomes the strict
  unique global maximizer, and GN recovers ${A, B}$ on Iteration 2.
]

#proof[
  Denote by $u_1'$ and $v_1'$ the triangle partners of $u_1$ in $A$
  and $v_1$ in $B$. In $G_0$, no intra-community pair uses $e_1$ or
  $e_2$ because both windmills have diameter $2$ whereas any detour
  between $A$ and $B$ has length $>= 3$. For cross-pairs
  $(x, y) in A times B$, the path through $e_1$ has length at most
  $3$. Enumerating every cross-pair where a shortest path uses $e_2$:
  (1) $(u_1, v_1)$ (length $1$ vs $3$: weight $1$);
  (2) $(u_1, v_1')$ and $(u_1', v_1)$ (length $2$ vs $3$: weight
  $1 + 1 = 2$);
  (3) $(h_A, v_1)$, $(u_1, h_B)$, and $(u_1', v_1')$ (tied at lengths
  $2, 2, 3$ respectively: weight $1/2 + 1/2 + 1/2 = 3/2$);
  (4) $(u_1, y)$ for the $2 r - 2$ vertices
  $y in B without {h_B, v_1, v_1'}$ (tied at length $3$ via
  $u_1 - v_1 - h_B - y$ and $u_1 - h_A - h_B - y$: weight $r - 1$);
  (5) symmetrically, $(x, v_1)$ for the $2 r - 2$ vertices
  $x in A without {h_A, u_1, u_1'}$ (weight $r - 1$).
  Every remaining cross-pair has length $>= 4$ via $e_2$ versus
  $<= 3$ via $e_1$, contributing $0$ to $e_2$. Summing gives
  $
    C_B (e_2; G_0) = 1 + 2 + 3/2 + 2(r - 1) = 2 r + 5/2,
  $
  and since all $(2 r + 1)^2$ cross-pairs use either $e_1$ or $e_2$,
  $C_B (e_1; G_0) = (2 r + 1)^2 - C_B (e_2; G_0) = 4 r^2 + 2 r - 3/2$.

  Now consider any internal spoke $e_"sp" = {h_A, u_2}$ where $u_2$
  lies outside ${u_1, u_1'}$ (with triangle partner $u_2'$). Because
  $u_2'$ is adjacent to $h_A$, no shortest path not ending at $u_2$
  traverses $e_"sp"$. Among pairs involving $u_2$: $u_2$ reaches
  $u_2'$ directly ($0$), reaches $h_A$ and the $2 r - 2$ other
  vertices of $A$ uniquely via $e_"sp"$ ($2 r - 1$), reaches every
  $y in B without {v_1}$ uniquely via $u_2 - h_A - h_B - y$ ($2 r$),
  and reaches $v_1$ via two length-$3$ paths
  ($u_2 - h_A - h_B - v_1$ and $u_2 - h_A - u_1 - v_1$) that both
  traverse $e_"sp"$ ($1$). Hence $C_B (e_"sp"; G_0) = 4 r > 2 r + 5/2$
  for all $r >= 2$, proving (i).

  After deleting $e_1$, edge $e_2$ is the sole bridge between $A$ and
  $B$, so $C_B (e_2; G_1) = (2 r + 1)^2 = 4 r^2 + 4 r + 1$. The
  highest-loaded internal edge in $G_1$ is ${h_A, u_1}$ (and
  symmetrically ${h_B, v_1}$), which carries the $2 r - 1$ vertices
  $A without {u_1, u_1'}$ to all $2 r + 1$ vertices of $B$ plus the
  $2 r - 1$ intra-$A$ pairs from $u_1$ to $A without {u_1, u_1'}$,
  totaling
  $
    C_B ({h_A, u_1}; G_1) = (2 r - 1)(2 r + 1) + (2 r - 1)
      = 4 r^2 + 2 r - 2.
  $
  Thus $C_B (e_2; G_1) - C_B ({h_A, u_1}; G_1) = 2 r + 3 > 0$,
  establishing (ii).
]

#subsec("D", "Modularity Resolution Limit on the GN Dendrogram")

Even when GN constructs a dendrogram whose exact split levels isolate
ground-truth communities, choosing the flat partition by maximizing
$Q(cal(C))$ inherits the Fortunato–Barthélemy resolution limit
@fortunato2007resolution.

#theorem(
  title: "Exact Ring-of-Cliques Resolution Limit of GN Cut Selection",
  origin: "Classical — Cited",
)[
  Following @fortunato2007resolution and @newman2004finding, let
  $G_(r, c)$ be a ring of $r >= 4$ cliques $K_c$ ($c >= 3$) linked
  cyclically by single bridge edges between portal vertices
  $p_i in C_i$, so each clique contributes $ell = binom(c, 2) + 1$
  edges and $m = r ell$. Then:
  (i) Girvan–Newman removes all $r$ bridge edges before any internal
  clique edge, so the true $r$-clique partition $cal(C)_r$ appears on
  the GN dendrogram.
  (ii) Merging two adjacent cliques $C_1, C_2$ along their single
  connecting bridge strictly increases Newman–Girvan modularity
  ($Delta Q > 0$) if and only if
  $
    r > 2 ell = c(c - 1) + 2 quad <==> quad ell < sqrt(m / 2).
  $
]

#proof[
  For (i), no inter-clique shortest path traverses any internal
  clique edge (all inter-clique traffic routes through the portals
  $p_i$). In the initial ring, summing bridge usage over all
  clique pairs via Theorem 2(a) and dividing by $r$ gives
  $C_B (e_"br") >= ((r^2 - 1) / 8) c^2$, whereas each portal spoke
  ${p_i, u}$ ($u in C_i without {p_i}$) carries exactly
  $1 + (r - 1) c$ pairs and non-portal edges carry $1$. Since
  $((r^2 - 1) / 8) c^2 - (r - 1) c - 1 > 0$ for all $r >= 4, c >= 3$,
  GN deletes a bridge first, opening the ring into a chain of $r$
  cliques. On any chain of $q >= 2$ cliques, the central bridge has
  $C_B = floor(q / 2) ceil(q / 2) c^2 >= (q - 1) c^2$, strictly
  exceeding every internal edge's load $<= 1 + (q - 1) c$ for
  $c >= 3$. By induction on $q$, GN removes all $r$ bridges before
  any internal edge, placing $cal(C)_r$ on $cal(T)$.

  For (ii), each clique $C_j$ has $e_(j j) = binom(c, 2)$ internal
  edges and original degree sum $a_j = 2 binom(c, 2) + 2 = 2 ell$,
  and adjacent cliques $C_1, C_2$ share $e_(1 2) = 1$ bridge edge.
  Merging $C_1$ and $C_2$ changes modularity by
  $
    Delta Q &= (e_(1 2)) / (m) - (2 a_1 a_2) / ((2 m)^2)
      = (1) / (r ell) - (2 (2 ell)^2) / (4 r^2 ell^2) \
      &= (1) / (r ell) - (2) / (r^2)
      = (r - 2 ell) / (r^2 ell).
  $
  Thus $Delta Q > 0$ if and only if $r > 2 ell = c(c - 1) + 2$;
  multiplying both sides by $ell$ gives $m = r ell > 2 ell^2$, i.e.,
  $ell < sqrt(m / 2)$. Whenever this inequality holds, maximizing $Q$
  along the GN dendrogram strictly rejects the true partition
  $cal(C)_r$ in favor of merged clique pairs.
]

#subsec("E", "Tight Worst-Case vs Balanced-Tree Complexity Bounds")

#proposition(
  title: "Tight Theta(m^2 n) Worst-Case and Theta(n^2) Balanced Bounds",
  origin: "Original Derivation",
)[
  Suppose GN recomputes $C_B$ exclusively on the connected component
  containing each deleted edge $e^*$. Then:
  (i) On every connected graph with $n$ vertices and $m$ edges, total
  time is $cal(O)(m^2 n)$ and space is $cal(O)(n + m)$.
  (ii) This bound is asymptotically tight even on trees ($m = n - 1$):
  on the star graph $K_(1, n - 1)$, GN executes in $Theta(n^3) =
  Theta(m^2 n)$ time, whereas on the unweighted path $P_n$, local
  component recomputation reduces total runtime to $Theta(n^2)$.
]

#proof[
  At deletion step $t in {1, dots, m}$, the active component has at
  most $n$ vertices and at most $m - t + 1$ edges, so Theorem 1 costs
  $cal(O)(n (m - t + 1))$; summing over $t = 1, dots, m$ gives
  $cal(O)(n m^2)$, while BFS queues, predecessor lists, and edge
  scores require $cal(O)(n + m)$ space.

  On the star $K_(1, n - 1)$ ($m = n - 1$), every edge is a bridge
  separating one leaf from the remaining $n - 1$ vertices, so all
  edges tie at $C_B = n - 1$. Deleting any edge detaches a singleton
  and leaves a star $K_(1, n - 2)$ of size $n - 1$. Recomputing on
  the surviving component at step $t$ costs
  $Theta((n - t)^2)$ (Theorem 1 on $n - t$ vertices and $n - t - 1$
  edges), giving total work
  $sum_(t = 0)^(n - 2) Theta((n - t)^2) = Theta(n^3) = Theta(m^2 n)$.

  Conversely, on the path $P_n$ ($1 - 2 - dots - n$), the $j$-th edge
  has $C_B (e_j) = j (n - j)$, uniquely maximized at the midpoint
  $j = floor(n / 2)$. Removing $e_j$ splits $P_n$ into two paths of
  size $floor(n / 2)$ and $ceil(n / 2)$, yielding the recursion
  $T(n) = c n^2 + T(floor(n / 2)) + T(ceil(n / 2))$. By the Master
  Theorem ($a = 2, b = 2, f(n) = c n^2 = Omega(n^(1 + 1))$),
  $T(n) = Theta(n^2)$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

The computational bottleneck of Girvan–Newman is its cubic-to-quintic
time scaling: $cal(O)(m^2 n)$ equals $cal(O)(n^3)$ on sparse graphs
($m = Theta(n)$) and degrades to $cal(O)(n^5)$ on dense graphs
($m = Theta(n^2)$). Consequently, exact GN is practical only up to
$n approx 10^3$–$5 times 10^3$ vertices (or $m <= 10^4$ edges) unless
truncated after the first $k$ splits or approximated by sampled
source BFS. Memory usage, however, is remarkably light at
$cal(O)(n + m)$, since the single-source dependency pass processes one
BFS tree at a time without materializing an $n times n$ distance
matrix.

Geometrically, GN is driven by *global geodesic bottlenecks* rather
than local density or Euclidean centroids. On geometric kNN
graphs, bridges connecting distinct non-convex manifolds (such as
concentric rings or two moons) concentrate shortest paths between
entire components, yielding sharp betweenness spikes (Theorem 2).
Unlike centroidal partitioners, GN can trace arbitrarily elongated,
curved, or unequal-variance shapes cleanly, provided inter-manifold
bridges are sparse relative to intra-manifold connectivity. Its
primary geometric weakness arises on overlapping spherical clusters
with dense boundary zones, where shortest paths disperse across many
parallel edges and peripheral spokes can outscore diffuse cuts.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time complexity],
    [$cal(O)(m^2 n)$ worst-case ($cal(O)(n^3)$ sparse, $cal(O)(n^5)$
      dense); $Theta(n^2)$ best-case on balanced trees (Prop. 2).],
    [Space complexity],
    [$cal(O)(n + m)$ using streaming single-source Brandes / NG
      BFS accumulation.],
    [Hyperparameters],
    [None required if selecting cut via max $Q$; optional target
      cluster count $k$ for early stopping.],
    [Determinism],
    [Deterministic up to exact $C_B$ ties (e.g., symmetric rings or
      stars).],
    [Output structure],
    [Complete divisive dendrogram $cal(T)$ plus optimal flat cut
      maximizing $Q(cal(C))$.],
    [Weighted extension],
    [Supported via Dijkstra ($cal(O)(m^2 n + m n^2 log n)$) using
      inverse similarity as edge length.],
    [Resolution profile],
    [Dendrogram resolves small modules cleanly, but max-$Q$ cut
      suffers $ell < sqrt(m / 2)$ limit (Thm. 3).],
  ),
  caption: [Operational and theoretical profile of Girvan–Newman.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

*1. Small-to-Medium Social and Organizational Networks:*
In seminal applications such as Zachary's karate club, college
football conferences, and scientific collaboration networks
@girvan2002community @newman2004finding ($n <= 10^3$), inter-group
liaisons carry disproportionate geodesic traffic. Edge betweenness
pinpoints those structural bridges directly, while the full divisive
dendrogram exposes nested sub-cliques (departments within faculties).

*2. Metabolic Pathways and Protein–Protein Interaction Modules:*
Biochemical reaction graphs and curated PPI subnetworks (e.g., yeast
complexes of $10^2$–$10^3$ proteins) exhibit modular functional
blocks linked by a handful of cross-talk signaling proteins. GN's
geodesic bottleneck criterion naturally isolates signaling bridges
without assuming equal module sizes or uniform node degrees.

*3. Power-Grid Islanding and Infrastructure Vulnerability Analysis:*
High-voltage transmission networks and transportation backbones are
sparse planar-like graphs ($m approx 1.5 n$) where cascade failures
propagate along high-betweenness tie-lines. GN simultaneously
discovers electrically cohesive control zones and ranks critical
cut-edges whose removal induces controlled islanding.

*4. Political Polarization and Legislative Co-Sponsorship Graphs:*
In roll-call voting and legislative co-sponsorship networks
($n approx 100$–$535$), bipartisan brokers span sparse ideological
bridges. Divisive betweenness cutting isolates partisan coalitions
and reveals sub-caucus fracture hierarchies as successive bridges are
severed.

*5. Non-Convex Manifold Separation on Low-Sample Geometric Graphs:*
When moderate-sized point clouds ($n <= 10^3$) form intertwined
ribbons, spirals, or moons, sparse mutual-kNN graphs exhibit
isolated spurious bridges between manifolds. Because those bridges
act as global geodesic bottlenecks (Theorem 2), GN severs them first,
often outperforming density or centroid methods once edge lengths are
suitably constructed.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Four failure regimes restrict Girvan–Newman:
(i) *Severe computational wall*: $cal(O)(m^2 n)$ renders exact GN
infeasible beyond a few thousand nodes; use Louvain
@blondel2008louvain, Leiden @traag2019leiden, or Infomap
@rosvall2008infomap on larger networks.
(ii) *Resolution-limit cut selection* (Theorem 3): selecting the
dendrogram cut via modularity merges small cliques when
$ell < sqrt(m / 2)$; inspecting multi-level dendrogram stability or
stopping at known $k$ avoids this artifact.
(iii) *Dense diffuse boundaries*: when two communities share many
parallel inter-cluster edges (high conductance, low modularity),
shortest-path load spreads across bridges so peripheral articulation
edges (leaves, dangling chains) can be peeled off first.
(iv) *Noise-bridge vulnerability on sparse graphs*: adding a single
long-range random shortcut can temporarily divert global shortest
paths (as in Proposition 1), though dynamic recomputation recovers
once the shortcut is excised.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Social / organizational graphs ($n < 2000$)],
    [Excellent: clear liaison bridges, full hierarchy desired.],
    [Check modularity resolution limit; compare Leiden.],
    [Power grids / sparse planar networks],
    [High: tie-lines carry peak geodesic load (Thm. 2).],
    [Use weighted Dijkstra for impedance distances.],
    [Non-convex kNN manifolds ($n < 10^3$)],
    [High: cuts inter-manifold bridges cleanly.],
    [Avoid overly dense kNN ($k$ too large).],
    [Web / social graphs ($n > 10^4$)],
    [Infeasible: $cal(O)(m^2 n)$ runtime explosion.],
    [Use Leiden @traag2019leiden or Infomap @rosvall2008infomap.],
    [Dense overlapping clusters],
    [Moderate / Poor: diffuse cuts peel leaves before split.],
    [Use Walktrap @pons2006walktrap or Spectral @shi2000ncut.],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for
    Girvan–Newman.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Girvan–Newman bridges classical graph cut optimization—Ratio Cut
@hagen1992ratiocut, Normalized Cut @shi2000ncut, and HCS min-cut
divisive clustering @hartuv2000hcs—with bottom-up modularity
optimizers such as Clauset–Newman–Moore @clauset2004cnm, Louvain
@blondel2008louvain, and Leiden @traag2019leiden. The following exact
inequality links GN's greedy step directly to Ratio Cut sparsity.

#proposition(
  title: "Exact Ratio-Cut Lower Bound on Maximum Edge Betweenness",
  origin: "Original Proof",
)[
  Let $G = (V, E)$ be any connected unweighted graph on $n$ vertices.
  For every non-trivial bipartition $(S, bar(S))$ ($emptyset != S
  subset.neq V$, $bar(S) = V without S$), define its Ratio Cut
  @hagen1992ratiocut by
  $
    "Rcut"(S, bar(S)) = abs(E(S, bar(S)))
      ( (1) / (abs(S)) + (1) / (abs(bar(S))) )
      = (n dot.op abs(E(S, bar(S)))) / (abs(S) dot.op abs(bar(S))).
  $
  Then the maximum edge betweenness in $G$ satisfies the universal
  sparsest-cut lower bound
  $
    max_(e in E) C_B (e; G)
      >= max_(emptyset != S subset.neq V)
      (n) / ("Rcut"(S, bar(S))),
  $
  with equality attained when $E(S, bar(S))$ is a single bridge
  separating two cliques (Theorem 2b).
]

#proof[
  Fix any non-trivial cut $(S, bar(S))$ with cut-edge set
  $F = E(S, bar(S))$ ($abs(F) >= 1$ since $G$ is connected). Every
  shortest path connecting a vertex $s in S$ to a vertex $t in bar(S)$
  must cross $F$ at least once, so for every $(s, t) in S times bar(S)$,
  $
    sum_(e in F) (sigma_(s t) (e)) / (sigma_(s t)) >= 1.
  $
  Summing this inequality over all $abs(S) dot.op abs(bar(S))$
  cross-cut pairs (and dropping non-negative contributions from
  intra-$S$ and intra-$bar(S)$ pairs) yields
  $
    sum_(e in F) C_B (e; G)
      >= sum_(s in S) sum_(t in bar(S)) sum_(e in F)
      (sigma_(s t) (e)) / (sigma_(s t))
      >= abs(S) dot.op abs(bar(S)).
  $
  By the pigeonhole principle (averaging over the $abs(F)$ cut edges),
  $
    max_(e in E) C_B (e; G)
      >= max_(e in F) C_B (e; G)
      >= (abs(S) dot.op abs(bar(S))) / (abs(F))
      = (n) / ("Rcut"(S, bar(S))).
  $
  Taking the maximum over all non-trivial cuts $(S, bar(S))$ completes
  the proof. When $F = {e^*}$ is a single bridge between two cliques
  $K_(n_1), K_(n_2)$, Theorem 2(b) gives
  $max_e C_B (e) = n_1 n_2 = n / "Rcut"(A, B)$, showing tightness.
]

Proposition 3 explains why Girvan–Newman acts as a greedy geodesic
relaxation of minimum Ratio Cut: any sufficiently sparse bottleneck
forces at least one cut edge to exhibit proportionally elevated
betweenness, guiding divisive splits without solving NP-hard
combinatorial cuts or computing eigenvectors @newman2006leading.

*Empirical Benchmark Behavior:* On the ClusterOverview 4-dataset
benchmark suite (Silhouette @rousseeuw1987silhouettes and Adjusted
Rand Index @hubert1985comparing), Girvan–Newman delivers the
strongest non-convex and ground-truth recovery scores in this four-
algorithm cohort: it achieves *Multi-Shape* Sil = *0.782*
(ARI = N/A)—substantially surpassing Chinese Whispers (0.527),
Infomap (0.478), and Walktrap (0.423)—and leads on *Iris* with
ARI = *0.657* (Sil = *0.451*). This empirical advantage
directly reflects Theorem 2 and Proposition 3: in moderate-size
kNN graphs, sparse geodesic bridges between distinct manifolds or
species clusters carry dominant shortest-path traffic and are severed
cleanly before internal edges. On *High-D*, GN obtains Sil =
*0.300* (ARI = N/A), comparable to the other graph methods.
Conversely, on *Customers* (Sil = *0.260*), diffuse overlapping
customer clouds lack sharp bridges, causing early peeling of
peripheral nodes and degrading silhouette separation relative to CW
and Infomap (both 0.430). Its mean runtime of *13.872 ms*—over
$20 times$ slower than Chinese Whispers (0.670 ms) and nearly
$2 times$ slower than Walktrap (7.944 ms)—quantifies the real cost of
iterative $cal(O)(m^2 n)$ dynamic BFS recomputation (Proposition 2).
