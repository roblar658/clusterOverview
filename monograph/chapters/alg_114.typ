#import "../template/ieee.typ": *

#alg_entry(
  id: 114,
  name: "Walktrap",
  category: "Cat. V: Graph, Spectral & Network",
  year: "2006",
  paradigm: "Random-Walk Spectral Distance Agglomerative Clustering",
  complexity: $cal(O)(m n H)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Walktrap, introduced by Pons and Latapy @pons2006walktrap, is a
bottom-up hierarchical community detection algorithm grounded in the
dynamics of short random walks. Its core hypothesis is geometric and
spectral: short random walks of fixed length $t$ (typically
$t in {3, 4, 5}$) tend to remain trapped within the same community,
so two vertices belonging to the same community induce nearly
identical $t$-step visitation distributions across the entire graph.

Let $G = (V, E, w)$ be a connected, simple, undirected weighted graph
with $n = abs(V)$ vertices, $m = abs(E)$ edges, symmetric adjacency
matrix $A in bb(R)_(>= 0)^(n times n)$ ($A_(u v) = w_(u v)$,
optionally augmented with positive self-loops to damp periodicity),
and diagonal degree matrix $D = "diag"(d_1, dots, d_n)$ with
$d_i = sum_(j = 1)^n A_(i j) > 0$ and total volume
$"vol"(V) = sum_(i = 1)^n d_i$. The one-step random-walk transition
matrix is $P = D^(-1) A$, with $t$-step transition probabilities
$P_(i j)^t = (P^t)_(i j)$ and stationary distribution
$pi_j = d_j / "vol"(V)$. Denote the $i$-th row of $P^t$ by the
column vector $P_(i dot.op)^t in bb(R)^n$.

#definition(title: "Pons-Latapy Random-Walk Distance and Ward Merge Cost")[
  For any walk length $t in bb(N)$, the *$t$-step random-walk
  distance* @pons2006walktrap between vertices $i, j in V$ is the
  degree-weighted $L_2$ metric between their $t$-step transition
  profiles:
  $
    r_(i j) (t) = sqrt(sum_(k = 1)^n
      ((P_(i k)^t - P_(j k)^t)^2) / (d_k))
      = norm(D^(-1/2) (P_(i dot.op)^t - P_(j dot.op)^t))_2.
  $
  For any non-empty community $C subset.eq V$, define its mean
  $t$-step transition vector
  $P_(C dot.op)^t = (1 / abs(C)) sum_(i in C) P_(i dot.op)^t$, and
  for disjoint communities $C_1, C_2 subset.eq V$, define
  $r_(C_1 C_2) (t) = norm(D^(-1/2) (P_(C_1 dot.op)^t
  - P_(C_2 dot.op)^t))_2$. Walktrap greedily merges the pair of
  *graph-adjacent* communities $(C_1, C_2)$ minimizing the Ward-type
  dispersion increase
  $
    Delta sigma(C_1, C_2) = (1) / (n)
      (abs(C_1) abs(C_2)) / (abs(C_1) + abs(C_2))
      r_(C_1 C_2)^2 (t).
  $
]

At each step $s = 1, dots, n - 1$, Walktrap merges the adjacent pair
$(C_1, C_2)$ minimizing $Delta sigma(C_1, C_2)$, updates the active
community graph and incident merge costs, and yields a full
agglomerative dendrogram $cal(T)$ from which the flat partition
maximizing Newman–Girvan modularity $Q$ @newman2004finding is
selected.

#algbox(
  title: "Walktrap Agglomerative Random-Walk Clustering",
  inputs: [Connected undirected graph $G = (V, E, w)$, walk length
    $t in bb(N)$ (default $t = 4$).],
  outputs: [Agglomerative dendrogram $cal(T)$ and modularity-optimal
    partition $cal(P)^*$.],
)[
  1. *Transition Profiles:* Compute $P = D^(-1) A$ and obtain rows
     $P_(i dot.op)^t$ via $t$ sparse vector-matrix multiplications
     from $e_i$ for all $i in V$.
  2. *Initialize:* Set singleton partition
     $cal(P)_1 = {{1}, dots, {n}}$. For every edge ${i, j} in E$,
     compute $r_(i j)^2 (t)$ and push $Delta sigma({i}, {j})$ onto a
     min-heap $cal(H)$.
  3. *Agglomerative Loop:* For $s = 1, dots, n - 1$:
     (a) Pop the adjacent pair $(C_1, C_2)$ with minimum
     $Delta sigma(C_1, C_2)$ from $cal(H)$ and form
     $C_3 = C_1 union C_2$.
     (b) Compute $P_(C_3 dot.op)^t = (abs(C_1) P_(C_1 dot.op)^t +
     abs(C_2) P_(C_2 dot.op)^t) / (abs(C_1) + abs(C_2))$.
     (c) For every active community $C$ adjacent to $C_3$: if $C$
     was adjacent to both $C_1$ and $C_2$, update $Delta sigma(C_3, C)$
     in $cal(O)(1)$ via the Lance–Williams recurrence (Theorem 2);
     otherwise evaluate $r_(C_3 C)^2 (t)$ directly in $cal(O)(n)$.
     (d) Record $cal(P)_(s + 1)$ in $cal(T)$ and track maximum
     modularity $Q(cal(P)_(s + 1))$.
  4. *Return:* Dendrogram $cal(T)$ and optimal partition $cal(P)^*$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five formal results for Walktrap: the exact spectral
embedding identity, a unified proof of Ward dispersion equivalence
and the Lance–Williams recurrence, an exact mixing-time collapse and
bipartite oscillation theorem, a 3-vertex constructive proof that
adjacency constraints induce dendrogram reversals, and the exact
$cal(O)(m n H)$ height-parameterized complexity bound.

#subsec("A", "Exact Spectral Decomposition and Euclidean Embedding")

#theorem(
  title: "Exact Spectral Representation of Random-Walk Distance",
  origin: "Classical — Cited",
)[
  Following @pons2006walktrap and @chung1997spectral, let
  $S = D^(1/2) P D^(-1/2) = D^(-1/2) A D^(-1/2)$ be the symmetric
  normalized adjacency matrix of connected $G$, with real eigenvalues
  $1 = lambda_1 > lambda_2 >= dots >= lambda_n >= -1$ and orthonormal
  eigenvectors $v_1, dots, v_n in bb(R)^n$ (where
  $v_1 (i) = sqrt(d_i / "vol"(V))$). Then for all $i, j in V$ and
  every $t in bb(N)$,
  $
    r_(i j)^2 (t) = sum_(alpha = 2)^n lambda_alpha^(2 t)
      ( (v_alpha (i)) / (sqrt(d_i))
      - (v_alpha (j)) / (sqrt(d_j)) )^2.
  $
  Hence $r_(i j) (t) = norm(Psi_t (i) - Psi_t (j))_2$ is an exact
  Euclidean pseudometric on $V$ under the $(n - 1)$-dimensional
  diffusion embedding
  $Psi_t (i) = (lambda_alpha^t v_alpha (i) / sqrt(d_i))_(alpha = 2)^n$.
]

#proof[
  Because $S = sum_(alpha = 1)^n lambda_alpha v_alpha v_alpha^top$ is
  real symmetric and $P = D^(-1/2) S D^(1/2)$, induction on $t$ gives
  $P^t = D^(-1/2) S^t D^(1/2) = sum_(alpha = 1)^n lambda_alpha^t
  (D^(-1/2) v_alpha) (D^(1/2) v_alpha)^top$. Reading entry $(i, k)$
  and dividing by $sqrt(d_k)$ yields
  $
    (P_(i k)^t) / (sqrt(d_k))
      = sum_(alpha = 1)^n lambda_alpha^t (v_alpha (i)) / (sqrt(d_i))
      v_alpha (k).
  $
  For $alpha = 1$, $v_1 (i) / sqrt(d_i) = 1 / sqrt("vol"(V))$ is
  constant across all $i in V$, so subtracting rows $i$ and $j$
  eliminates the $alpha = 1$ term identically:
  $
    (P_(i k)^t - P_(j k)^t) / (sqrt(d_k))
      = sum_(alpha = 2)^n lambda_alpha^t
      ( (v_alpha (i)) / (sqrt(d_i))
      - (v_alpha (j)) / (sqrt(d_j)) ) v_alpha (k).
  $
  Squaring and summing over $k = 1, dots, n$, then invoking
  orthonormality $sum_(k = 1)^n v_alpha (k) v_beta (k) =
  delta_(alpha beta)$, cancels all cross-terms $alpha != beta$ and
  produces the exact spectral formula. Since the right-hand side is
  the squared Euclidean distance between $Psi_t (i)$ and $Psi_t (j)$
  in $bb(R)^(n - 1)$, $r_(i j) (t)$ satisfies non-negativity, symmetry,
  and the triangle inequality (and is a strict metric whenever no two
  vertices share an identical row in $P^t$).
]

#subsec("B", "Ward Dispersion Equivalence and Lance-Williams Update")

#theorem(
  title: "Exact Ward Objective Equivalence and Lance-Williams Formula",
  origin: "Original Unified Proof",
)[
  Define the within-community spectral dispersion of partition
  $cal(P) = {C_1, dots, C_k}$ by
  $
    sigma(cal(P)) = (1) / (n) sum_(r = 1)^k sum_(i in C_r)
      r_(i C_r)^2 (t)
      = (1) / (n) sum_(r = 1)^k sum_(i in C_r)
      norm(Psi_t (i) - bar(Psi)_t (C_r))_2^2,
  $
  where $bar(Psi)_t (C) = (1 / abs(C)) sum_(i in C) Psi_t (i)$. Then:
  (i) Merging disjoint communities $C_1, C_2$ into
  $C_3 = C_1 union C_2$ increases $sigma$ by identically
  $sigma(cal(P)') - sigma(cal(P)) = Delta sigma(C_1, C_2)$
  @ward1963hierarchical.
  (ii) For any third disjoint community $C$ ($n_1 = abs(C_1)$,
  $n_2 = abs(C_2)$, $n_c = abs(C)$), the updated merge cost obeys the
  exact $cal(O)(1)$ Lance–Williams recurrence @lance1967general:
  $
    Delta sigma(C_3, C) =
      ((n_1 + n_c) Delta sigma(C_1, C)
      + (n_2 + n_c) Delta sigma(C_2, C)
      - n_c Delta sigma(C_1, C_2)) / (n_1 + n_2 + n_c).
  $
]

#proof[
  Let $mu_r = bar(Psi)_t (C_r)$ for $r in {1, 2, 3, c}$, so
  $mu_3 = (n_1 mu_1 + n_2 mu_2) / n_3$ with $n_3 = n_1 + n_2$. By
  Theorem 1, $r_(C_1 C_2) (t) = norm(mu_1 - mu_2)_2$. For any
  $i in C_1$, expanding
  $norm(Psi_t (i) - mu_3)_2^2 = norm((Psi_t (i) - mu_1)
  + (mu_1 - mu_3))_2^2$ and summing over $i in C_1$ cancels the
  cross-term because $sum_(i in C_1) (Psi_t (i) - mu_1) = 0$:
  $
    sum_(i in C_1) norm(Psi_t (i) - mu_3)_2^2
      = sum_(i in C_1) norm(Psi_t (i) - mu_1)_2^2
      + n_1 norm(mu_1 - mu_3)_2^2.
  $
  Substituting $mu_1 - mu_3 = (n_2 / n_3)(mu_1 - mu_2)$ and adding
  the symmetric identity for $C_2$ gives
  $
    n (sigma(cal(P)') - sigma(cal(P)))
      &= ( n_1 (n_2 / n_3)^2 + n_2 (n_1 / n_3)^2 )
      norm(mu_1 - mu_2)_2^2 \
      &= (n_1 n_2) / (n_1 + n_2) r_(C_1 C_2)^2 (t),
  $
  dividing by $n$ proves (i).

  For (ii), set $a = mu_1 - mu_c$, $b = mu_2 - mu_c$ (so
  $a - b = mu_1 - mu_2$), $alpha = n_1 / n_3$, and
  $beta = n_2 / n_3 = 1 - alpha$. Direct expansion verifies the
  Euclidean Stewart identity
  $
    norm(alpha a + beta b)_2^2
      = alpha norm(a)_2^2 + beta norm(b)_2^2
      - alpha beta norm(a - b)_2^2.
  $
  Since $mu_3 - mu_c = alpha a + beta b$, multiplying both sides by
  $n_3 n_c / (n (n_3 + n_c))$ yields
  $
    Delta sigma(C_3, C) &= (n_1 n_c) / (n (n_3 + n_c)) norm(a)_2^2
      + (n_2 n_c) / (n (n_3 + n_c)) norm(b)_2^2 \
      &quad - (n_1 n_2 n_c) / (n n_3 (n_3 + n_c))
      norm(mu_1 - mu_2)_2^2.
  $
  Substituting $(n_r n_c / n) norm(mu_r - mu_c)_2^2 =
  (n_r + n_c) Delta sigma(C_r, C)$ ($r in {1, 2}$) and
  $(n_1 n_2 / (n n_3)) norm(mu_1 - mu_2)_2^2 = Delta sigma(C_1, C_2)$
  yields the exact Lance–Williams identity.
]

#subsec("C", "Spectral Mixing Collapse and Bipartite Persistence")

Choosing the walk length $t$ governs a fundamental signal-to-noise
trade-off captured precisely by the spectrum of $S$.

#theorem(
  title: "Mixing-Time Collapse, Signal-Noise Window, and Bipartite Trap",
  origin: "Original Proof",
)[
  Let $lambda_* = max(abs(lambda_2), abs(lambda_n))$ be the
  second-largest eigenvalue modulus of $S = D^(-1/2) A D^(-1/2)$.
  (i) *Global collapse:* For all $i != j$ and all $t >= 1$,
  $
    r_(i j) (t) <= lambda_*^t sqrt((1) / (d_i) + (1) / (d_j)).
  $
  (ii) *Two-timescale filtering:* Suppose $S$ exhibits a spectral gap
  $1 = lambda_1 > lambda_2 >= dots >= lambda_k > lambda_"bulk" :=
  max_(alpha > k) abs(lambda_alpha) >= 0$. Decomposing
  $r_(i j)^2 (t) = r_(i j, "sig")^2 (t) + r_(i j, "noise")^2 (t)$
  into modes $alpha in {2, dots, k}$ versus $alpha > k$, the
  noise-to-signal ratio shrinks exponentially in $t$:
  $
    (r_(i j, "noise")^2 (t)) / (r_(i j, "sig")^2 (t))
      <= ( (lambda_"bulk") / (lambda_k) )^(2 t)
      (r_(i j, "noise")^2 (0)) / (r_(i j, "sig")^2 (0)).
  $
  (iii) *Bipartite failure without self-loops:* If $G$ is connected
  and bipartite ($V = V_+ union V_-$) with no self-loops, then
  $lambda_n = -1$ persists with unit modulus for every $t >= 1$,
  forcing $r_(i j)^2 (t) >= 4 / "vol"(V)$ for every cross-part pair
  $i in V_+, j in V_-$, whereas the lazy walk
  $P_"lazy" = 1/2 (I + P)$ eliminates negative eigenvalues completely.
]

#proof[
  For (i), Theorem 1 and $abs(lambda_alpha) <= lambda_*$ for all
  $alpha >= 2$ give
  $
    r_(i j)^2 (t) <= lambda_*^(2 t) sum_(alpha = 2)^n
      ( (v_alpha (i)) / (sqrt(d_i))
      - (v_alpha (j)) / (sqrt(d_j)) )^2.
  $
  Because $v_1 (i) / sqrt(d_i) - v_1 (j) / sqrt(d_j) = 0$, the sum
  can be extended to $alpha = 1, dots, n$. The orthogonal matrix
  $V = [v_1, dots, v_n]$ satisfies $V V^top = I$, so
  $sum_(alpha = 1)^n v_alpha (i)^2 = 1$ and
  $sum_(alpha = 1)^n v_alpha (i) v_alpha (j) = 0$ for $i != j$,
  yielding
  $
    sum_(alpha = 1)^n ( (v_alpha (i)) / (sqrt(d_i))
      - (v_alpha (j)) / (sqrt(d_j)) )^2
      = (1) / (d_i) + (1) / (d_j).
  $
  Taking square roots proves (i).

  For (ii), every mode $alpha in {2, dots, k}$ obeys
  $lambda_alpha^(2 t) >= lambda_k^(2 t)$, so
  $r_(i j, "sig")^2 (t) >= lambda_k^(2 t) r_(i j, "sig")^2 (0)$,
  while every mode $alpha > k$ obeys
  $lambda_alpha^(2 t) <= lambda_"bulk"^(2 t)$, so
  $r_(i j, "noise")^2 (t) <= lambda_"bulk"^(2 t) r_(i j, "noise")^2 (0)$.
  Dividing the two inequalities immediately gives (ii).

  For (iii), on any connected loopless bipartite graph $V_+ union V_-$,
  the vector $v_n (i) = +sqrt(d_i / "vol"(V))$ on $V_+$ and
  $-sqrt(d_i / "vol"(V))$ on $V_-$ is unit-norm and satisfies
  $S v_n = -v_n$, so $lambda_n = -1$ and $lambda_n^(2 t) = 1$ for all
  $t >= 1$. For $i in V_+, j in V_-$, the $alpha = n$ term in
  Theorem 1 alone equals
  $1 dot.op (2 / sqrt("vol"(V)))^2 = 4 / "vol"(V)$, proving
  $r_(i j)^2 (t) >= 4 / "vol"(V)$ permanently. Under
  $P_"lazy" = 1/2 (I + P)$, the symmetrized matrix is
  $S_"lazy" = 1/2 (I + S)$ with spectrum
  $mu_alpha = (1 + lambda_alpha) / 2 in [0, 1]$, mapping
  $lambda_n = -1 |-> 0$ and removing the bipartite mode entirely.
]

#subsec("D", "Adjacency-Constrained Dendrogram Reversals")

Classical unconstrained Euclidean Ward clustering @ward1963hierarchical
always produces monotonic (reversal-free) merge heights, yet
Walktrap's graph-adjacency constraint breaks monotonicity even on
three vertices.

#proposition(
  title: "Monotonicity Break under Graph-Adjacency Constraints",
  origin: "Original Counterexample",
)[
  (i) In *unconstrained* Ward clustering, if $(C_1, C_2)$ globally
  minimizes $Delta sigma$ over all pairs, then every subsequent merge
  cost with $C_3 = C_1 union C_2$ satisfies
  $Delta sigma(C_3, C) >= Delta sigma(C_1, C_2)$ (zero reversals).
  (ii) Under Walktrap's *graph-adjacency* constraint, merge heights
  can strictly decrease: on the 3-vertex path $P_3$ ($2 - 1 - 3$)
  with $P = D^(-1) A$, for *every* walk length $t >= 1$, Walktrap
  executes Step 1 at height $Delta sigma^((1)) = 1 / 6$ and Step 2 at
  strictly lower height $Delta sigma^((2)) = 1 / 18 =
  (1 / 3) Delta sigma^((1))$.
]

#proof[
  For (i), if $(C_1, C_2)$ is a global minimizer, then
  $Delta sigma(C_1, C) >= Delta sigma(C_1, C_2)$ and
  $Delta sigma(C_2, C) >= Delta sigma(C_1, C_2)$ for every $C$.
  Substituting into Theorem 2(ii) gives
  $
    Delta sigma(C_3, C)
      >= ((n_1 + n_c) + (n_2 + n_c) - n_c) / (n_1 + n_2 + n_c)
      Delta sigma(C_1, C_2)
      = Delta sigma(C_1, C_2).
  $

  For (ii), on $P_3$ ($2 - 1 - 3$, $d_1 = 2, d_2 = d_3 = 1$), the
  transition matrix is
  $P_(1 dot.op)^1 = (0, 1/2, 1/2)$ and
  $P_(2 dot.op)^1 = P_(3 dot.op)^1 = (1, 0, 0)$. Because $2$ and $3$
  are non-adjacent structural twins, induction yields
  $P_(2 dot.op)^t = P_(3 dot.op)^t$ for all $t >= 1$ (alternating
  between $(1, 0, 0)$ for odd $t$ and $(0, 1/2, 1/2)$ for even $t$,
  opposite to $P_(1 dot.op)^t$). Thus for every $t >= 1$,
  $r_(2 3)^2 (t) = 0$ (ineligible because ${2, 3} in.not E$) and
  $
    r_(1 2)^2 (t) = r_(1 3)^2 (t)
      = ((1 - 0)^2) / (2) + ((0 - 1/2)^2) / (1)
      + ((0 - 1/2)^2) / (1) = 1.
  $
  Step 1 must merge an adjacent pair, say $C_1 = {1}$ and
  $C_2 = {2}$, at height
  $Delta sigma^((1)) = (1 / 3)(1 dot.op 1 / 2) dot.op 1 = 1 / 6$.
  Applying Theorem 2(ii) with $C = {3}$ ($Delta sigma(C_1, C) = 1/6$,
  $Delta sigma(C_2, C) = 0$, $Delta sigma(C_1, C_2) = 1/6$) gives
  $
    Delta sigma^((2)) = Delta sigma(C_3, {3})
      = (2 (1 / 6) + 2 (0) - 1 (1 / 6)) / (3)
      = (1) / (18) < (1) / (6).
  $
]

#subsec("E", "Exact Dendrogram-Height Complexity Analysis")

#lemma(
  title: "Height-Parameterized Time and Quadratic Space Bounds",
  origin: "Classical — Cited",
)[
  Following @pons2006walktrap, let $H$ ($ceil(log_2 n) <= H <= n - 1$)
  be the height of Walktrap's merge dendrogram. Computing all $t$-step
  profiles and constructing the full hierarchy requires
  $cal(O)(t m n + m n H) = cal(O)(m n H)$ time (for $t = cal(O)(1)$)
  and $cal(O)(n^2)$ space. On sparse graphs ($m = cal(O)(n)$) with
  balanced dendrograms ($H = cal(O)(log n)$), runtime is
  $cal(O)(n^2 log n)$, while the worst case is $cal(O)(m n^2)$.
]

#proof[
  Computing one row $P_(i dot.op)^t$ from $e_i$ via $t$ sparse
  vector-matrix products traverses each edge at most twice per step
  ($cal(O)(t m)$), so all $n$ rows require $cal(O)(t m n)$ time and
  $cal(O)(n^2)$ memory; initializing $r_(i j)^2 (t)$ along all $m$
  edges costs $cal(O)(m n)$. During agglomeration, each of the
  $n - 1$ merges updates one length-$n$ vector $P_(C_3 dot.op)^t$
  ($cal(O)(n^2)$ total). For a community $C$ adjacent to $C_3$,
  updating $Delta sigma(C_3, C)$ via Theorem 2(ii) costs $cal(O)(1)$
  if $C$ was adjacent to both $C_1, C_2$, and $cal(O)(n)$ otherwise.
  Every such $cal(O)(n)$ evaluation is chargeable to an original edge
  ${u, v} in E$ joining $C$ to $C_3$ whose endpoint community just
  merged. Since each vertex ascends at most $H$ dendrogram levels,
  each of the $m$ edges triggers at most $2 H$ recomputations,
  yielding $cal(O)(m n H)$ distance-update time (plus
  $cal(O)(m H log n)$ heap operations, dominated by $cal(O)(m n H)$).
]

#sec("III", "Complexity, Geometry and Operational Profile")

Walktrap occupies a practical computational middle ground between
fast local heuristics (Chinese Whispers, Louvain) and costly global
geodesic solvers (Girvan–Newman): its $cal(O)(n^2 log n)$ typical
sparse complexity comfortably handles networks up to $n approx 10^4$
–$5 times 10^4$ vertices on standard hardware, though its
$cal(O)(n^2)$ dense storage of active $P_(C dot.op)^t$ vectors
eventually becomes a memory barrier on million-node graphs (where
sparse truncation or Leiden @traag2019leiden is required).

Geometrically, Theorem 1 and Theorem 2 prove that Walktrap performs
exact *Ward minimum-variance clustering in $t$-step diffusion-map
space* $Psi_t (V) subset.eq bb(R)^(n - 1)$, restricted to connected
graph merges. Consequently, Walktrap inherits Ward's strong inductive
bias toward *compact, roughly spherical, balanced-variance diffusion
basins*. In noisy Euclidean kNN graphs, averaging over $t = 4$
random-walk hops smooths out isolated spurious bridge edges far more
robustly than 1-hop label propagation, producing clean, high-silhouette
partitions on globular or mildly overlapping clusters. Conversely,
like Ward's method @ward1963hierarchical @mullner2011modern, Walktrap
penalizes long, high-diameter ribbons or concentric rings because
splitting an elongated manifold drastically reduces squared Euclidean
dispersion in $Psi_t$-space.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time complexity],
    [$cal(O)(m n H)$ ($cal(O)(n^2 log n)$ typical sparse balanced
      tree; $cal(O)(m n^2)$ worst case, Lemma 1).],
    [Space complexity],
    [$cal(O)(n^2)$ to store active community $t$-step vectors
      $P_(C dot.op)^t$ plus $cal(O)(m)$ heap.],
    [Hyperparameters],
    [Walk length $t in {3, 4, 5}$ (default $t = 4$); optional
      self-loop / lazy damping for bipartite graphs (Thm. 3iii).],
    [Determinism],
    [Fully deterministic (up to exact floating-point ties in
      $Delta sigma$).],
    [Output structure],
    [Complete agglomerative dendrogram $cal(T)$ (may contain
      adjacency reversals, Prop. 1ii) plus max-$Q$ cut.],
    [Scale invariance],
    [Strictly invariant under $w |-> alpha w$ ($alpha > 0$) for
      $P$, scaling $r_(i j)^2 (t)$ uniformly by $alpha^(-1)$.],
    [Resolution profile],
    [Dendrogram captures multi-scale structure; flat max-$Q$ cut
      inherits modularity resolution limit @fortunato2007resolution.],
  ),
  caption: [Operational and theoretical profile of Walktrap.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

*1. Single-Cell Transcriptomics and Cell-Type Hierarchies:*
In scRNA-seq pipelines (`igraph` / Seurat / Bioconductor), cells are
embedded as sparse kNN or shared-nearest-neighbor graphs of
$10^3$–$3 times 10^4$ nodes. Walktrap's deterministic $t$-step
diffusion averages out dropout noise on individual edges (Theorem 3ii)
while returning a full biological lineage dendrogram from broad cell
classes down to fine transcriptional subtypes.

*2. Gene Co-Expression and Functional Protein Modules:*
Weighted gene co-expression networks (WGCNA) exhibit dense noisy
correlation backgrounds where single-edge bridges are unreliable.
By integrating all length-$t$ weighted paths simultaneously through
$P^t$, Walktrap suppresses spurious pairwise correlations and isolates
compact co-regulated modules via Ward dispersion minimization.

*3. Medium-Scale Social, Citation, and Collaboration Communities:*
On benchmark networks of $10^3$–$5 times 10^4$ vertices (co-authorship,
co-purchasing, communication graphs) where Girvan–Newman's
$cal(O)(m^2 n)$ runtime is prohibitive, Walktrap preserves a
deterministic, reproducible community hierarchy without the seed
instability of Chinese Whispers or Louvain @blondel2008louvain.

*4. Image Superpixel and Region Merging on Affinity Graphs:*
On pixel or patch affinity graphs, short random walks capture
multi-path perceptual continuity across textured regions better than
single-edge cuts. Graph-constrained Ward merging enforces spatial
contiguity natively (every merged cluster remains a connected
subgraph) while minimizing spectral variance within each segment.

*5. Financial Asset Correlation Networks and Portfolio Sectoring:*
Stock-return correlation graphs exhibit nested hierarchical sector
and sub-industry structure with heavy market-mode noise. Because the
trivial eigenvector $v_1$ cancels identically in Theorem 1 and bulk
noise eigenvalues decay as $lambda_"bulk"^(2 t)$ (Theorem 3ii),
Walktrap cleanly extracts coherent industry sectors without requiring
explicit eigendecomposition or pre-specifying the sector count $k$.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Four mathematical failure modes govern Walktrap:
(i) *Quadratic memory bottleneck*: storing dense rows $P_(C dot.op)^t$
in $cal(O)(n^2)$ space exhausts RAM beyond $n approx 5 times 10^4$
vertices; use Louvain / Leiden or Infomap at web scale.
(ii) *Walk-length sensitivity* (Theorem 3): choosing $t$ too small
($t = 1, 2$) leaves intra-community noise unsmoothed, whereas $t$
larger than the mixing time collapses all distances exponentially
toward zero ($lambda_*^t -> 0$), causing numerical tie degeneracy.
(iii) *Bipartite and periodic graphs* (Theorem 3iii): without
self-loops or lazy transitions, $lambda_n = -1$ never decays and
splits bipartite graphs along their bipartition rather than true
communities.
(iv) *Elongated non-convex manifolds*: Ward's variance criterion
favors compact diffusion balls over long thin ribbons, causing
artificial bisection of curved shapes (e.g., long moons or rings).

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Single-cell / gene expression kNN graphs],
    [Excellent: deterministic hierarchy, $t$-step noise damping.],
    [Switch to Leiden @traag2019leiden above $5 times 10^4$ cells.],
    [Medium social / finance networks],
    [High: fast $cal(O)(n^2 log n)$, reproducible dendrogram.],
    [Check max-$Q$ cut against resolution limit.],
    [Bipartite / near-bipartite networks],
    [Poor without self-loops ($lambda_n = -1$, Thm. 3iii).],
    [Add self-loops / lazy walk $P_"lazy"$.],
    [Elongated non-convex manifolds],
    [Moderate / Poor: Ward bias splits long ribbons.],
    [Use Girvan–Newman or Spectral @shi2000ncut.],
    [Million-node web / directed graphs],
    [Infeasible: $cal(O)(n^2)$ RAM and undirected assumption.],
    [Use Infomap @rosvall2008infomap or Leiden @traag2019leiden.],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for
    Walktrap.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Walktrap forms a clean mathematical bridge between Normalized
Spectral Clustering @ng2002spectral @vonluxburg2007tutorial, Ward's
agglomerative method @ward1963hierarchical, Markov Clustering (MCL)
@vandongen2000mcl, and Infomap @rosvall2008infomap.

#proposition(
  title: "Exact Spectral-Filter Duality with Normalized Spectral Clustering",
  origin: "Original Derivation",
)[
  Let $L_"sym" = I - S$ be the symmetric normalized graph Laplacian
  with eigenvalues $0 = nu_1 < nu_2 <= dots <= nu_n <= 2$
  ($nu_alpha = 1 - lambda_alpha$) and eigenvectors $v_1, dots, v_n$.
  Define the spectral-filter kernel distance
  $
    d_f^2 (i, j) = sum_(alpha = 2)^n f(nu_alpha)^2
      ( (v_alpha (i)) / (sqrt(d_i))
      - (v_alpha (j)) / (sqrt(d_j)) )^2.
  $
  Then:
  (i) Unnormalized random-walk spectral clustering
  @shi2000ncut @vonluxburg2007tutorial embeds vertices via the *hard
  band-pass step filter* $f_"spec" (nu_alpha) = bb(1)[alpha <= k]$
  (requiring eigendecomposition and pre-specified $k$) and applies
  Euclidean $k$-means.
  (ii) Walktrap replaces the step filter with the *smooth polynomial
  filter* $f_"WT" (nu) = (1 - nu)^t$ (evaluated implicitly in
  $cal(O)(t m n)$ without eigendecomposition) and replaces $k$-means
  with graph-constrained Ward agglomeration.
  (iii) Replacing $P$ by the lazy walk $P_"lazy" = 1/2 (I + P)$
  strictly yields the monotone low-pass filter
  $f_"lazy" (nu) = (1 - nu / 2)^t in [0, 1]$, which decreases strictly
  with Laplacian frequency $nu in [0, 2]$.
]

#proof[
  In random-walk spectral clustering @vonluxburg2007tutorial, the
  eigenvectors of $L_"rw" = I - P = D^(-1/2) L_"sym" D^(1/2)$ are
  $u_alpha = D^(-1/2) v_alpha$ with coordinates
  $u_alpha (i) = v_alpha (i) / sqrt(d_i)$ and trivial constant mode
  $u_1 prop bold(1)$. Embedding vertex $i$ by
  $(u_2 (i), dots, u_k (i))$ yields squared Euclidean distance
  $sum_(alpha = 2)^k (v_alpha (i) / sqrt(d_i) -
  v_alpha (j) / sqrt(d_j))^2$, which matches $d_f^2 (i, j)$ with
  $f_"spec" (nu_alpha) = bb(1)[alpha <= k]$, proving (i).
  Substituting $lambda_alpha = 1 - nu_alpha$ into Theorem 1 directly
  yields $r_(i j)^2 (t) = d_f^2 (i, j)$ with
  $f_"WT" (nu) = (1 - nu)^t$, and Theorem 2(i) identifies Walktrap's
  merge step with graph-constrained Ward minimization in this space,
  proving (ii). Under $P_"lazy" = 1/2 (I + P)$,
  $mu_alpha = (1 + lambda_alpha) / 2 = 1 - nu_alpha / 2$, so
  $f_"lazy" (nu) = (1 - nu / 2)^t$, which is strictly decreasing on
  $[0, 2]$ with $f_"lazy" (0) = 1$ and $f_"lazy" (2) = 0$, proving
  (iii).
]

Proposition 2 clarifies why Walktrap often outperforms hard spectral
clustering when the true cluster count $k$ is unknown or the
eigengap is blurred: the smooth filter $(1 - nu)^t$ weights modes
continuously by persistence rather than imposing a brittle cutoff at
a guessed $k$.

*Empirical Benchmark Behavior:* On the ClusterOverview 4-dataset
benchmark suite (Silhouette @rousseeuw1987silhouettes and Adjusted
Rand Index @hubert1985comparing), Walktrap achieves the highest
*Iris* Silhouette (*0.476*, with ARI = *0.450*) among all four
chapters (surpassing Girvan–Newman's 0.451, Infomap's 0.339, and
Chinese Whispers' 0.319). This mirrors Theorem 3(ii) and Proposition
2: $t$-step diffusion damps noisy cross-class kNN edges while
Ward's minimum-variance criterion carves compact, well-separated
species clusters. On *High-D* (Sil = *0.302*, ARI = N/A) and
*Customers* (Sil = *0.356*), Walktrap performs solidly in the
middle of the cohort, avoiding Girvan–Newman's peripheral leaf
peeling on Customers (0.260) while trailing CW and Infomap (0.430).
Conversely, on *Multi-Shape* (Sil = *0.423*, ARI = N/A),
Walktrap scores lowest of the four methods—directly reflecting Ward's
dispersion bias (Theorem 2i), which bisects elongated non-convex
moons and rings into compact sub-segments rather than cutting sparse
inter-manifold bridges like Girvan–Newman (0.782). Its mean runtime
of *7.944 ms* makes Walktrap nearly $2 times$ faster than
Girvan–Newman (13.872 ms) and $12 times$ faster than Infomap
(95.352 ms), confirming the practical efficiency of Lance–Williams
updates (Theorem 2ii and Lemma 1).
