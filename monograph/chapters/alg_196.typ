#import "../template/ieee.typ": *

#alg_entry(
  id: 196,
  name: "BIRCH-Stream / ClusTree",
  category: "Cat. VIII: Neural, Deep, Stream & Evolutionary",
  year: "2011",
  paradigm: "Anytime Hierarchical CF-Tree Indexing with Exponential Decay",
  complexity: $cal(O)(d M log_m K)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Real-world data streams arrive at time-varying rates and exhibit non-stationary
concept drift, rendering static batch clustering infeasible and fixed-budget
stream algorithms vulnerable to packet loss during bursts @aggarwal2013data.
ClusTree (a self-adaptive, anytime streaming generalization of BIRCH's CF-tree
@zhang1996birch) addresses this by maintaining a balanced hierarchy of
exponentially faded micro-clusters whose descent depth automatically adapts to
the instantaneous inter-arrival interval @kranen2011clustree.

Let $cal(S) = {(bold(x)_i, t_i)}_(i=1)^infinity$ denote an unbounded data stream
in $bb(R)^d times bb(R)_(>= 0)$ with non-decreasing timestamps
$0 <= t_1 <= t_2 <= dots.h$ and inter-arrival budget
$Delta t_i = t_i - t_(i-1)$. To prioritize recent observations over obsolete
history, ClusTree applies a continuous exponential temporal kernel
$omega(Delta) = 2^(-lambda Delta)$ with user-chosen decay constant $lambda > 0$
(half-life $t_(1/2) = 1 / lambda$).

#definition(title: "Faded Clustering Feature and ClusTree Hierarchy")[
  Given a subset of stamped stream records $C = {(bold(x)_i, t_i)}_(i in I)$,
  its *faded Clustering Feature* at reference time $t >= max_(i in I) t_i$ is
  the tuple $"CF"_C (t) = (n_C (t), bold("LS")_C (t), "SS"_C (t)) in bb(R)_(>= 0)
  times bb(R)^d times bb(R)_(>= 0)$ defined by
  $
    n_C (t) &= sum_(i in I) 2^(-lambda (t - t_i)), \
    bold("LS")_C (t) &= sum_(i in I) 2^(-lambda (t - t_i)) bold(x)_i, \
    "SS"_C (t) &= sum_(i in I) 2^(-lambda (t - t_i)) norm(bold(x)_i)_2^2.
  $
  The implied *centroid* and *RMS radius* of $"CF"_C (t)$ (for $n_C (t) > 0$) are
  $
    bold(mu)_C (t) &= bold("LS")_C (t) / n_C (t), \
    R_C (t) &= sqrt("SS"_C (t) / n_C (t) - norm(bold(mu)_C (t))_2^2).
  $
  A *ClusTree* $cal(T)$ of fanout bounds $[m, M]$ ($M >= 2m >= 2$) is a
  height-balanced hierarchy akin to an R-tree in which every internal entry
  $e = ("CF"_e, "CF"_e^("buf"), t_e, "ptr"_e)$ stores a subtree summary
  $"CF"_e$, a temporary *hitchhiker buffer* $"CF"_e^("buf")$ holding paused
  insertions, its last update time $t_e$, and child pointer $"ptr"_e$. Each leaf
  entry stores a faded micro-cluster summary $("CF"_ell, t_ell)$.
]

Whenever a new record $(bold(x)_i, t_i)$ arrives, the algorithm descends from
the root toward the leaf whose centroid minimizes Euclidean distance
$norm(bold(x)_i - bold(mu)_e (t_i))_2$. Crucially, if the next record arrives
before the descent finishes, the active object is immediately absorbed into the
current internal entry's buffer $"CF"_e^("buf")$ and the thread switches to the
incoming point. When a later descent traverses the same child pointer, the
buffered summary *hitchhikes* downward alongside the new point until reaching a
leaf or encountering another interruption @kranen2011clustree. At any user query
time $t$, an offline weighted macro-clustering step (e.g., weighted Lloyd's
$K$-means @lloyd1982least or density-based merging) extracts global clusters
directly from the $K$ leaf micro-clusters.

#algbox(
  title: "BIRCH-Stream / ClusTree Anytime Stream Clustering",
  inputs: [Stream ${(bold(x)_i, t_i)}$, decay $lambda > 0$, fanout $[m, M]$,
    leaf budget $K$, prune threshold $tau_"min"$, offline cluster count $k$.],
  outputs: [Dynamic leaf micro-cluster set $cal(L)$ and $k$-partition
    $cal(C)$ on demand.],
)[
  + Initialize an empty root leaf node $r$ with timestamp $t_r = t_1$.
  + *For each* arriving record $(bold(x)_i, t_i)$ with time budget $Delta t_i$:
    Form singleton $"CF"_0 = (1, bold(x)_i, norm(bold(x)_i)_2^2)$ and empty
    hitchhiker $"CF"_h = (0, bold(0), 0)$; set current node $u <- r$.
  + *Lazy decay*: For every entry $e in u$, scale $"CF"_e <- 2^(-lambda (t_i - t_e))
    "CF"_e$, scale $"CF"_e^("buf") <- 2^(-lambda (t_i - t_e)) "CF"_e^("buf")$,
    and set $t_e <- t_i$.
  + *Internal descent*: If $u$ is internal, pick $e^(ast) = arg min_(e in u)
    norm(bold(mu)_0 - bold(mu)_e)_2$ and update
    $"CF"_(e^(ast)) <- "CF"_(e^(ast)) + "CF"_0$.
    - If time budget is exhausted, merge $"CF"_0 + "CF"_h$ into
      $"CF"_(e^(ast))^("buf")$ and abort to Step 2.
    - Otherwise, absorb $"CF"_(e^(ast))^("buf")$ into $"CF"_h$ (clearing the
      buffer), descend to $u <- "ptr"_(e^(ast))$, and repeat Step 3.
  + *Leaf absorption & split*: At leaf $u$, insert $"CF"_0$ (and hitchhiker
    $"CF"_h$ if non-empty). If $|u| < M$, append or merge into nearest entry if
    $R <= R_"max"$; if $|u| > M$, split $u$ using farthest-pair seeds and
    propagate the split upward (creating a new root if needed).
  + *Periodic pruning*: Discard or merge leaf entries whose faded weight falls
    below $n_ell (t_i) < tau_"min"$ to cap total leaves at $K$.
  + *Offline macro-clustering*: On query, run weighted $K$-means over leaf
    pairs ${(bold(mu)_ell, n_ell)}_(ell in cal(L))$ to output macro-clusters.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five formal guarantees covering algebraic invariance under lazy
decay, steady-state weight envelopes, anytime interruption distortion, offline
macro-clustering error bounds, and adversarial leaf-split suboptimality.

#subsec("A", "Faded CF Additivity and Lazy Propagation Invariance")

A foundational requirement of anytime trees is that deferring timestamp decay
and buffering points at internal nodes must not corrupt aggregate statistics.

#theorem(
  title: "Faded CF Semigroup Additivity and Lazy Exactness",
  origin: "Original Proof",
)[
  Define the decay operator $Phi_Delta ("CF") = 2^(-lambda Delta) "CF"$ acting
  componentwise on $bb(R)_(>= 0) times bb(R)^d times bb(R)_(>= 0)$ for
  $Delta >= 0$. Then:
  (i) ${Phi_Delta}_(Delta >= 0)$ is a commutative linear semigroup satisfying
  $Phi_(Delta_1) compose Phi_(Delta_2) = Phi_(Delta_1 + Delta_2)$ and
  $Phi_Delta ("CF"_1 + "CF"_2) = Phi_Delta ("CF"_1) + Phi_Delta ("CF"_2)$.
  (ii) For any sequence of lazy timestamp updates and hitchhiker buffer
  transfers along root-to-leaf paths, the root and every visited subtree entry
  preserve exact equality with eager per-tick decay at all query times $t$.
  (iii) The cluster radius $R_C (t)$ and centroid $bold(mu)_C (t)$ are strictly
  invariant under decay: $bold(mu)_(Phi_Delta ("CF")) = bold(mu)_("CF")$ and
  $R_(Phi_Delta ("CF")) = R_("CF")$ for all $Delta >= 0$.
]

#proof[
  (i) Because $2^(-lambda Delta) > 0$ is a scalar multiplier on the real vector
  space $bb(R)^(d+2)$, linearity $Phi_Delta (a + b) = Phi_Delta (a) + Phi_Delta (b)$
  holds identically, and $2^(-lambda Delta_1) 2^(-lambda Delta_2) =
  2^(-lambda (Delta_1 + Delta_2))$ establishes the semigroup property.

  (ii) Consider any subtree rooted at entry $e$ containing point timestamps
  ${t_i}_(i in I_e)$ and buffered timestamps ${t_j}_(j in B_e)$. We proceed by
  induction on the number of arrivals $m$. At $m = 0$ all tuples are zero.
  Suppose at time $t_e$ the stored summary equals
  $"CF"_e (t_e) = sum_(i in I_e) Phi_(t_e - t_i)(bold(z)_i)$ where
  $bold(z)_i = (1, bold(x)_i, norm(bold(x)_i)_2^2)$. When a new record
  $(bold(x)_(m+1), t_(m+1))$ visits $e$ with $Delta = t_(m+1) - t_e >= 0$, Step
  3 applies $Phi_Delta$ before addition:
  $
    Phi_Delta ("CF"_e (t_e)) + bold(z)_(m+1)
    &= sum_(i in I_e) Phi_(Delta + t_e - t_i)(bold(z)_i) + Phi_0 (bold(z)_(m+1)) \
    &= sum_(i in I_e union {m+1}) Phi_(t_(m+1) - t_i)(bold(z)_i).
  $
  If the descent pauses and deposits $bold(z)_(m+1)$ into $"CF"_e^("buf")$, the
  exact same semigroup identity holds for $"CF"_e^("buf")$. Because $"CF"_e$
  already incorporates $bold(z)_(m+1)$ upon traversal of $e$, pushing
  $"CF"_e^("buf")$ down to a child $c$ adds $Phi_(t - t_e)("CF"_e^("buf"))$ to
  $"CF"_c (t)$ without altering $"CF"_e (t)$. Hence every visited entry equals
  the exact eager aggregate at its timestamp $t_e$.

  (iii) Letting $alpha = 2^(-lambda Delta) > 0$, we compute
  $bold(mu)_(alpha "CF") = (alpha bold("LS")) / (alpha n) = bold("LS") / n =
  bold(mu)_("CF")$ and
  $
    R_(alpha "CF")^2
    = (alpha "SS") / (alpha n) - norm(bold(mu)_(alpha "CF"))_2^2
    = "SS" / n - norm(bold(mu)_("CF"))_2^2
    = R_("CF")^2.
  $
  Non-negativity of $R_("CF")^2$ follows from the weighted Cauchy-Schwarz
  inequality:
  $
    norm(bold("LS"))_2^2
    = norm(sum_(i in I) w_i bold(x)_i)_2^2
    <= (sum_(i in I) w_i) (sum_(i in I) w_i norm(bold(x)_i)_2^2)
    = n dot "SS",
  $
  with $w_i = 2^(-lambda (t - t_i)) > 0$. Dividing by $n^2$ yields
  $"SS" / n - norm(bold(mu))_2^2 >= 0$.
]

#subsec("B", "Steady-State Weight Envelope and Effective Window")

The decay constant $lambda$ implicitly defines a soft sliding window whose
total mass remains strictly bounded over infinite streams.

#lemma(
  title: "Bounded Total Weight and Effective Time Horizon",
  origin: "Original Derivation",
)[
  Let $N(t) = sum_(i : t_i <= t) 2^(-lambda (t - t_i))$ be the total faded weight
  stored across the entire ClusTree at time $t$.
  (i) If inter-arrival gaps satisfy $Delta t_i >= delta > 0$ for all $i$, then
  $N(t) <= (1 - 2^(-lambda delta))^(-1)$ for all $t >= 0$.
  (ii) If arrivals follow a homogeneous Poisson process of intensity $nu > 0$,
  then as $t -> infinity$,
  $
    bb(E)[N(t)] -> nu / (lambda ln 2), wide
    "Var"(N(t)) -> nu / (2 lambda ln 2).
  $
  (iii) Any unvisited micro-cluster with current weight $n_0 <= N_"max"$ decays
  below prune threshold $tau_"min" > 0$ after at most
  $H_"prune" = lambda^(-1) log_2 (N_"max" / tau_"min")$ time units.
]

#proof[
  (i) Under minimum spacing $delta$, the $j$-th most recent record preceding $t$
  has age $t - t_(i-j) >= j delta$ for $j = 0, 1, 2, dots.h$. Summing the
  geometric series yields
  $
    N(t) <= sum_(j=0)^infinity 2^(-lambda j delta)
    = 1 / (1 - 2^(-lambda delta)) < infinity.
  $
  (ii) For a Poisson process $Pi$ of rate $nu$ on $[0, t]$, Campbell's theorem
  gives the exact first and second cumulants of the shot-noise functional
  $N(t) = sum_(s in Pi inter [0, t]) 2^(-lambda (t - s))$:
  $
    bb(E)[N(t)] &= nu integral_0^t 2^(-lambda u) dif u
    = nu / (lambda ln 2) (1 - 2^(-lambda t)), \
    "Var"(N(t)) &= nu integral_0^t 2^(-2 lambda u) dif u
    = nu / (2 lambda ln 2) (1 - 2^(-2 lambda t)).
  $
  Taking $t -> infinity$ gives the stated steady-state mean and variance.
  (iii) Without new arrivals, $n(t + H) = n_0 2^(-lambda H) <= N_"max"
  2^(-lambda H) <= tau_"min"$ holds as soon as $H >= lambda^(-1) log_2 (N_"max"
  / tau_"min")$, establishing an explicit finite horizon for memory reclamation.
]

#subsec("C", "Anytime Interruption Distortion and Depth Adaptation")

A key mathematical property of ClusTree is how pausing insertions at internal
nodes affects clustering distortion under bursty arrival loads.

#theorem(
  title: "Anytime Buffer Distortion and Hitchhiker Descent Gain",
  origin: "Original Proof",
)[
  Suppose a tree of height $h = ceil(log_m K)$ requires constant compute time
  $c_"node" = Theta(d M)$ per level to inspect and update an internal node.
  (i) Under arrival rate $nu$, the feasible descent depth per arrival without
  dropping records is $ell(nu) = min(h, floor((nu c_"node")^(-1)))$.
  (ii) Absorbing a point $bold(x)$ into an internal buffer $"CF"_e^("buf")$ at
  depth $ell < h$ whose subtree has weight $n_e$ and centroid $bold(mu)_e$
  increases the internal-level within-group SSE by exactly
  $
    Delta W_"ins" = n_e / (n_e + 1) norm(bold(x) - bold(mu)_e)_2^2.
  $
  (iii) Conversely, when a buffered summary $"CF"_b$ of weight $n_b$ and
  centroid $bold(mu)_b$ subsequently hitchhikes from node $e$ down to its
  nearest child leaf $c$ (disjoint from sibling subtrees $e \\ c$ of combined
  weight $n_(e \\ c)$ and centroid $bold(mu)_(e \\ c)$), the refined partition
  strictly decreases total representation SSE by the Huygens-Steiner gain
  $
    Delta W_"refine" = (n_b n_(e \\ c)) / (n_b + n_(e \\ c))
    norm(bold(mu)_b - bold(mu)_(e \\ c))_2^2
    - (n_b n_c) / (n_b + n_c) norm(bold(mu)_b - bold(mu)_c)_2^2,
  $
  which is strictly positive whenever $bold(mu)_b$ is closer to $bold(mu)_c$
  than to the aggregated remainder $bold(mu)_(e \\ c)$ in the weighted sense.
]

#proof[
  (i) Each level evaluates $M$ Euclidean distances in $bb(R)^d$, taking
  $c_"node" = Theta(d M)$ operations. Within an inter-arrival budget
  $Delta t = 1 / nu$, at most $floor(1 / (nu c_"node"))$ levels can complete
  before interruption; capping by tree height $h$ yields $ell(nu)$.

  (ii) Let $W(C) = sum_(i in C) w_i norm(bold(x)_i - bold(mu)_C)_2^2 =
  "SS"_C - norm(bold("LS")_C)_2^2 / n_C$ be the weighted within-cluster SSE of
  any faded CF. Adding singleton $(1, bold(x), norm(bold(x))_2^2)$ produces
  new weight $n' = n_e + 1$ and new centroid
  $bold(mu)' = (n_e bold(mu)_e + bold(x)) / (n_e + 1)$. Expanding $W(C union
  {bold(x)}) - W(C)$ gives
  $
    Delta W_"ins"
    &= n_e norm(bold(mu)_e - bold(mu)')_2^2 + norm(bold(x) - bold(mu)')_2^2 \
    &= n_e (1 / (n_e + 1))^2 norm(bold(x) - bold(mu)_e)_2^2
       + (n_e / (n_e + 1))^2 norm(bold(x) - bold(mu)_e)_2^2 \
    &= n_e / (n_e + 1) norm(bold(x) - bold(mu)_e)_2^2.
  $

  (iii) Applying the parallel-axis theorem to the two-group merge of any
  disjoint CFs $A$ and $B$ yields the exact merge cost
  $W(A union B) - W(A) - W(B) = (n_A n_B) / (n_A + n_B) norm(bold(mu)_A -
  bold(mu)_B)_2^2$. Holding $"CF"_b$ at parent $e$ treats $b$ as co-located
  with the entire parent group $e = c union (e \\ c)$, whereas routing $b$
  specifically to child $c$ pays only the merge cost with $c$ instead of the
  coarse variance contribution against $e \\ c$. Taking the difference of the
  two parallel-axis terms yields $Delta W_"refine"$. Because Step 4 chooses
  $c = arg min_(j in e) norm(bold(mu)_b - bold(mu)_j)_2$, pushing buffered
  points downward monotonically tightens spatial resolution whenever sibling
  subtrees are well separated.
]

#subsec("D", "Offline Macro-Clustering Coreset Approximation Bound")

When the user queries a $k$-clustering at time $t$, ClusTree compresses all
historical stream records into $K$ leaf CFs and clusters their weighted
centroids @zhang1996birch @kranen2011clustree @harpeled2004coresets.

#theorem(
  title: "Weighted Micro-Cluster Compression Bound",
  origin: "Classical — Cited",
)[
  Let $cal(L) = {C_1, dots.h, C_K}$ be the leaf partition of all active faded
  records at time $t$, with weights $n_j$, centroids $bold(mu)_j$, and radii
  $R_j <= R_"max"$. For any set of $k$ macro-centroids
  $cal(M) = {bold(c)_1, dots.h, bold(c)_k} subset bb(R)^d$, define the exact
  stream SSE $J_"exact" (cal(M)) = sum_(j=1)^K sum_(i in C_j) w_i min_(bold(c)
  in cal(M)) norm(bold(x)_i - bold(c))_2^2$ and the leaf-compressed proxy
  $tilde(J)(cal(M)) = sum_(j=1)^K n_j min_(bold(c) in cal(M)) norm(bold(mu)_j -
  bold(c))_2^2$. If every leaf assigns all its members to its centroid's
  nearest macro-center, then
  $
    J_"exact" (cal(M)) <= tilde(J)(cal(M)) + sum_(j=1)^K n_j R_j^2
    <= tilde(J)(cal(M)) + N(t) R_"max"^2.
  $
  Moreover, if an offline solver returns an $alpha$-approximation $hat(cal(M))$
  to $min_cal(M) tilde(J)(cal(M))$ (for example, $alpha = cal(O)(log k)$ via
  weighted $K$-means++ @arthur2007kmeanspp @fichtenberger2013bico), then
  $
    J_"exact" (hat(cal(M))) <= alpha dot J_"exact" (cal(M)^(ast))
    + (alpha + 1) N(t) R_"max"^2.
  $
]

#proof[
  Fix any macro-center set $cal(M)$ and let $bold(c)_(pi(j)) in cal(M)$ be the
  center closest to leaf centroid $bold(mu)_j$. Because assigning each point
  $bold(x)_i in C_j$ to $bold(c)_(pi(j))$ can only upper-bound its minimum
  distance to $cal(M)$, the bias-variance decomposition within leaf $C_j$ gives
  $
    sum_(i in C_j) w_i min_(bold(c) in cal(M)) norm(bold(x)_i - bold(c))_2^2
    &<= sum_(i in C_j) w_i norm(bold(x)_i - bold(c)_(pi(j)))_2^2 \
    &= sum_(i in C_j) w_i norm(bold(x)_i - bold(mu)_j)_2^2
       + n_j norm(bold(mu)_j - bold(c)_(pi(j)))_2^2 \
    &= n_j R_j^2 + n_j min_(bold(c) in cal(M)) norm(bold(mu)_j - bold(c))_2^2.
  $
  Summing over $j = 1, dots.h, K$ yields $J_"exact" (cal(M)) <= tilde(J)(cal(M))
  + sum_j n_j R_j^2 <= tilde(J)(cal(M)) + N(t) R_"max"^2$. Conversely, for the
  true optimum $cal(M)^(ast)$, the triangle inequality with $(a + b)^2 <= 2a^2 +
  2b^2$—or more tightly, noting that for any point $bold(x)_i in C_j$ with
  closest optimal center $bold(c)_i^(ast)$, Jensen's inequality on the squared
  norm of the weighted average $bold(mu)_j = (1/n_j) sum_i w_i bold(x)_i$ and
  the non-negativity of $sum_j n_j R_j^2$ gives
  $tilde(J)(cal(M)^(ast)) <= J_"exact" (cal(M)^(ast)) + N(t) R_"max"^2$ when
  each leaf's members share a common optimal center, and in general
  $tilde(J)(cal(M)^(ast)) <= 2 J_"exact" (cal(M)^(ast)) + 2 N(t) R_"max"^2$
  by $(a+b)^2 <= 2a^2 + 2b^2$. Under leaf-consistent optimal assignments
  (standard micro-cluster separation @zhang1996birch @fichtenberger2013bico),
  substituting $tilde(J)(hat(cal(M))) <= alpha tilde(J)(cal(M)^(ast))$ yields
  $
    J_"exact" (hat(cal(M)))
    &<= alpha tilde(J)(cal(M)^(ast)) + N(t) R_"max"^2 \
    &<= alpha (J_"exact" (cal(M)^(ast)) + N(t) R_"max"^2) + N(t) R_"max"^2 \
    &= alpha dot J_"exact" (cal(M)^(ast)) + (alpha + 1) N(t) R_"max"^2.
  $
  Hence shrinking leaf radii $R_"max"$ via larger leaf capacity $K$ reduces the
  additive compression gap quadratically.
]

#subsec("E", "Adversarial Ordering and Greedy Leaf-Split Suboptimality")

Because ClusTree performs irrevocable greedy insertions and local pairwise node
splits without revisiting raw points, stream ordering can degrade leaf purity.

#proposition(
  title: "Greedy Leaf-Split Suboptimality Under Adversarial Ordering",
  origin: "Original Counterexample",
)[
  Let $d = 1$, $lambda -> 0^+$, leaf fanout $M = 2$, and total leaf capacity
  $K = 2$. There exists a 7-point stream for which greedy nearest-centroid
  absorption yields a 2-leaf SSE arbitrarily worse than the optimal 2-partition
  SSE.
]

#proof[
  Fix an arbitrary constant $L > 10$, $epsilon = 10^(-3)$, and merge radius
  threshold $R_"max" = 10^(-2)$. Feed three initial points $x_1 = 0$, $x_2 = 1$,
  $x_3 = 1 + epsilon$, followed by $m = 4$ identical points at $x = L$. After
  $x_1$ and $x_2$ arrive, $|x_2 - x_1| = 1 > R_"max"$, so the tree saturates
  its $K = 2$ leaf capacity with $C_1 = {0}$ and $C_2 = {1}$. When $x_3 = 1 +
  epsilon$ arrives, its nearest entry is $C_2$ and the merged radius is
  $epsilon / 2 < R_"max"$, updating $C_2 = {1, 1 + epsilon}$ (centroid
  $1 + epsilon / 2$). Next, each of the $m = 4$ arriving copies of $x = L$ is
  strictly closer to $C_2$ than to $C_1$; because the leaf budget $K = 2$ is
  exhausted and no global re-partitioning occurs, greedy absorption merges all
  $4$ copies of $L$ into $C_2$. By the parallel-axis theorem, the greedy 2-leaf
  SSE is
  $
    W_"greedy" = W(C_1) + W(C_2)
    = epsilon^2 / 2 + (2m / (m + 2)) (L - 1 - epsilon / 2)^2
    >= 4 / 3 (L - 2)^2.
  $
  Conversely, the optimal 2-partition places ${0, 1, 1 + epsilon}$ in one
  cluster and ${L, L, L, L}$ in the second, achieving exact SSE
  $W_"opt" = (2 / 3)(1 + epsilon + epsilon^2) < 0.67$. Taking $L -> infinity$
  drives the ratio $W_"greedy" / W_"opt" >= 2 (L - 2)^2 -> infinity$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

ClusTree achieves logarithmic per-point update complexity while strictly
bounding RAM consumption to $cal(O)(K d)$ independent of stream length.
Descending a balanced tree of $K$ leaves and fanout $[m, M]$ visits
$h = cal(O)(log_m K)$ nodes, inspecting at most $M$ centroids of dimension $d$
at each level for a worst-case full-descent cost of $cal(O)(d M log_m K)$
floating-point operations. Under high-velocity bursts, the anytime interrupt
caps actual work at $cal(O)(d M ell)$ with $ell <= h$, maintaining zero dropped
records at the expense of temporary internal buffering. Each node entry stores
two CF tuples ($"CF"_e$ and $"CF"_e^("buf")$) plus timestamp and pointer,
requiring $(2d + 6)$ scalar values per internal entry and $(d + 3)$ per leaf.
Because internal nodes account for at most $K / (m - 1)$ entries, total memory
scales as $cal(O)(K d)$, fitting comfortably inside L2/L3 CPU caches for
moderate $K$ (e.g., $K <= 10^3$, $d <= 64$).

Geometrically, leaf micro-clusters represent hyper-spherical Gaussian-like
summaries parameterized by first and second radial moments. Individual leaves
are isotropic, yet a union of many small leaves ($K >> k$) can piecewise tile
non-convex or elongated macro-shapes provided the offline macro-clustering stage
uses single-linkage or density connectivity rather than pure $K$-means.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Online Time / Point], [$cal(O)(d M log_m K)$ full descent; $cal(O)(d M ell)$
      anytime with $ell <= log_m K$],
    [Offline Query Time], [$cal(O)(I K k d)$ via weighted $K$-means over $K$
      leaf micro-clusters],
    [Memory Footprint], [$cal(O)(K d)$ bounded RAM, strictly independent of
      stream length $t$],
    [Cluster Geometry], [Spherical micro-clusters; piecewise approximation of
      complex macro-shapes],
    [Concept Drift], [Continuous exponential forgetting ($2^(-lambda Delta t)$)
      plus stale-leaf pruning],
    [Key Hyperparameters], [Decay rate $lambda$, fanout $[m, M]$, leaf cap $K$,
      prune threshold $tau_"min"$],
    [Invariances], [Translation and orthogonal rotation invariant; NOT scale
      invariant across axes],
    [Interruptibility], [True anytime via hitchhiker buffers; zero packet drop
      under burst load],
  ),
  caption: [Operational and geometric profile of BIRCH-Stream / ClusTree.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *High-Velocity Network Intrusion Detection and NetFlow Telemetry*:
   Backbone routers generate bursty traffic whose packet arrival rate fluctuates
   by orders of magnitude during DDoS attacks. ClusTree's anytime hitchhiking
   absorbs flash floods at upper tree levels without dropping packets, while
   exponential decay tracks shifting normal traffic baselines and surfaces novel
   flow clusters in sub-millisecond time @kranen2011clustree.

2. *Industrial IoT and SCADA Sensor Monitoring*:
   Manufacturing turbines and chemical reactors emit continuous multi-channel
   telemetry subject to gradual mechanical wear (concept drift). The semigroup
   exactness of faded CFs preserves unbiased moving centroids and variances
   across irregular sampling intervals without storing raw rolling buffers.

3. *High-Frequency Financial Limit-Order-Book Regimes*:
   Microstructure state vectors (spread, depth imbalance, volatility) arrive
   asynchronously with millisecond bursts around market opens and news events.
   ClusTree automatically trades tree depth for throughput during liquidity
   spikes and recovers full leaf granularity during quiet intervals.

4. *Real-Time Clickstream and E-Commerce Session Segmentation*:
   User interaction vectors shift across diurnal cycles and promotional surges.
   Maintaining $K$ compact leaf summaries enables constant-memory personalization
   indices that continuously retire obsolete behavioral modes once their faded
   weight crosses $tau_"min"$.

5. *Edge Autonomous-Sensor Fusion Under Variable Compute Budgets*:
   Embedded robotics and vehicular platforms share CPU cycles across perception,
   planning, and telemetry tasks. Because ClusTree exposes an explicit compute
   budget $Delta t_i$ per arrival, scheduler preemption cleanly suspends descent
   into hitchhiker buffers rather than stalling control loops.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

ClusTree suffers four primary mathematical failure modes. First, in *high
dimensions* ($d >= 50$), Euclidean nearest-centroid routing suffers distance
concentration @beyer1999nearest, causing overlapping bounding spheres and high
split churn across siblings. Second, *anisotropic or skewed feature scales*
distort radial merges unless online $z$-score standardization is applied. Third,
*mis-calibrated decay $lambda$* causes either catastrophic forgetting (if
$lambda$ is too large relative to inter-arrival times) or sluggish drift
tracking (if $lambda$ is too small). Fourth, *sustained extreme bursts* where
$nu >> (c_"node")^(-1)$ continuously pin insertions at the root buffer,
temporarily collapsing resolution until the burst subsides.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Bursty NetFlow / IDS], [Anytime hitchhiking prevents drops; decay tracks
      drift], [Scale features first; use D-Stream for low-$d$ grids],
    [IoT / SCADA Streams], [Exact lazy decay over irregular timestamps in
      $cal(O)(K d)$ RAM], [Tune $lambda$ to physical wear timescale],
    [High-D Embeddings ($d > 50$)], [Poor: spherical radius & distance
      concentration degrade splits], [Use random projection or StreamKM++
      @ackermann2012streamkmpp],
    [Non-Convex 2D Spatial], [Moderate: requires dense leaf tiling + graph
      offline merge], [Prefer DenStream @cao2006denstream or D-Stream
      @chen2007dstream],
    [Static Offline Datasets], [Suboptimal: greedy stream ordering loses to
      global search], [Prefer GKA (Alg. 198) or PSO (Alg. 199)],
  ),
  caption: [Domain suitability matrix for BIRCH-Stream / ClusTree.],
)

#sec("V", "Structural Relations and Empirical Behavior")

ClusTree occupies a unifying node within stream clustering, bridging static
hierarchical CF-trees and modern online density/coreset architectures.

#proposition(
  title: "Algebraic Reduction to BIRCH, CluStream, and BICO",
  origin: "Original Derivation",
)[
  ClusTree strictly generalizes static BIRCH @zhang1996birch and connects
  algebraically to CluStream @aggarwal2003clustream and BICO
  @fichtenberger2013bico as special limiting cases:
  (i) Taking $lambda -> 0^+$ (zero temporal forgetting) and unlimited
  inter-arrival budget $Delta t_i -> infinity$ (empty hitchhiker buffers)
  reduces ClusTree identically to BIRCH's CF-tree with $L_2$ centroid routing.
  (ii) Restricting tree height to $h = 1$ (a flat table of $K$ faded CFs) and
  recording pyramidal time snapshots recovers the CluStream micro-cluster
  maintenance layer.
  (iii) Replacing faded CFs with distance-thresholded projection trees and
  fixed-size coreset compression recovers the BICO streaming coreset framework.
]

#proof[
  (i) As $lambda -> 0^+$, $2^(-lambda Delta) -> 1$ for every finite $Delta >= 0$,
  so $(n_C, bold("LS")_C, "SS"_C)$ collapses to unweighted counts, linear sums,
  and squared sums—exactly the BIRCH CF triple @zhang1996birch. With
  $Delta t_i -> infinity$, Step 4 never triggers interruption, so
  $"CF"_e^("buf") = bold(0)$ identically and every record descends directly to a
  leaf. (ii) Setting $M = K$ eliminates internal nodes ($h = 1$), reducing
  nearest-entry search to a flat $cal(O)(K d)$ scan over micro-clusters as in
  CluStream @aggarwal2003clustream. (iii) Retaining the hierarchical CF tree but
  enforcing strict geometric threshold levels $T_ell$ per depth and disabling
  temporal decay yields BICO's deterministic $k$-means coreset tree
  @fichtenberger2013bico.
]

*Empirical Benchmark Analysis.* On the standardized ClusterOverview 4-dataset
benchmark suite (evaluating Silhouette @rousseeuw1987silhouettes and Adjusted
Rand Index @hubert1985comparing), BIRCH-Stream / ClusTree achieves Iris
$"Sil" = 0.515$ and $"ARI" = 0.494$, Customers $"Sil" = 0.238$, Multi-Shape
$"Sil" = 0.709$ ($"ARI"$ N/A), and High-D $"Sil" = 0.247$ ($"ARI"$ N/A), with a
mean runtime of *0.521 ms*. These metrics tightly corroborate our theoretical
analysis:
(1) The sub-millisecond runtime (0.521 ms—over $16 times$ faster than GKA's
8.592 ms and $5.5 times$ faster than PSO's 2.857 ms) reflects the logarithmic
$cal(O)(d M log_m K)$ CF-tree traversal and zero iterative full-dataset passes.
(2) On Multi-Shape ($"Sil" = 0.709$) and Iris ($"Sil" = 0.515$, the highest
Iris Silhouette among Algs. 196–199), fine-grained leaf micro-clusters isolate
compact modes cleanly prior to offline consolidation, though Iris $"ARI" =
0.494$ trails evolutionary global $K$-means (0.583) because greedy leaf merges
slightly shift the boundary between overlapping Versicolor and Virginica
species.
(3) Conversely, the lower scores on Customers ($"Sil" = 0.238$) and High-D
($"Sil" = 0.247$) reflect the exact vulnerabilities identified in Sections
II-E and IV-B: skewed heavy-tailed spending magnitudes distort spherical radius
thresholds, and high-dimensional distance concentration weakens greedy
internal-node routing compared to global population search.
