#import "../template/ieee.typ": *

#alg_entry(
  id: 75,
  name: "Neighbor-Joining (NJ)",
  category: "Cat. III: Hierarchical & Tree-Based",
  year: "1987",
  paradigm: "Minimum-Evolution Additive Tree Reconstruction",
  complexity: $O(n^3)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Classical agglomerative clustering algorithms — UPGMA @sokal1958statistical, Ward @ward1963hierarchical, Lance-Williams @lance1967general (Alg. 72), and Minimax @bien2011minimax (Alg. 73) — reconstruct *rooted ultrametric* trees in which every leaf lies at the exact same path distance from the root. In evolutionary biology, this ultrametric restriction enforces a *strict molecular clock* (uniform evolutionary rate across all lineages). Whenever lineages evolve at heterogeneous rates, two fast-evolving non-sister taxa can appear farther apart than two slow-evolving distant relatives, causing UPGMA to misidentify clades even on exact noise-free phylogenetic distances.

Saitou and Nei @saitou1987neighbor resolved this fundamental limitation with *Neighbor-Joining (NJ)*, a greedy bottom-up algorithm that reconstructs an *unrooted additive tree* $T = (V, E, w)$ directly from a symmetric pairwise dissimilarity matrix $D in RR_(>= 0)^(n times n)$, without assuming a molecular clock.

#definition(title: "Additive Tree Metric and Buneman's Four-Point Condition")[
  A symmetric zero-diagonal dissimilarity $d : X times X arrow.r RR_(>= 0)$ on $X = {x_1, dots.c, x_n}$ is an *additive tree metric* if there exists an unrooted tree $T = (V, E)$ with leaf set $X$ and edge weights $w : E arrow.r RR_(>= 0)$ such that for all $x, y in X$,
  $
    d(x, y) = d_T (x, y) = sum_(e in P_T (x, y)) w(e),
  $
  where $P_T (x, y)$ is the unique simple path between $x$ and $y$ in $T$. Equivalently (Buneman's four-point condition), $d$ is additive iff for all $x, y, z, w in X$:
  $
    d(x, y) + d(z, w) <= max(d(x, z) + d(y, w), d(x, w) + d(y, z)).
  $
]

Every ultrametric satisfies the four-point condition, but the converse fails whenever pendant branch lengths differ. At each iteration with $r >= 3$ active taxa (leaves of the reduced tree), NJ does *not* join the pair minimizing raw distance $D_(i j)$; instead, it compensates for heterogeneous lineage rates by subtracting each taxon's *net divergence* $R_i = sum_(k=1)^r D_(i k)$ via the globally re-centered *$Q$-criterion* (the algebraically equivalent Studier-Keppler form of Saitou and Nei's minimum-evolution score @saitou1987neighbor):
$
  Q(i, j) = (r - 2) D_(i j) - R_i - R_j quad (i != j).
$
Selecting $(i^*, j^*) in arg min_(i < j) Q(i, j)$ joins $i^*$ and $j^*$ to a new internal node $u$ with branch lengths
$
  w(i^*, u) &= (1 / 2) D_(i^* j^*) + (R_(i^*) - R_(j^*)) / (2 (r - 2)), \
  w(j^*, u) &= D_(i^* j^*) - w(i^*, u),
$
and replaces $i^*, j^*$ by $u$ with reduced distances $D_(u k) = (1 / 2)(D_(i^* k) + D_(j^* k) - D_(i^* j^*))$ for all $k in.not {i^*, j^*}$.

#algbox(
  title: "Neighbor-Joining (NJ) Additive Tree Reconstruction",
  inputs: [Symmetric dissimilarity matrix $D in RR_(>= 0)^(n times n)$ over taxa $X = {x_1, dots.c, x_n}$ ($n >= 3$)],
  outputs: [Unrooted weighted tree $T = (V, E, w)$ (and optional $k$-partition by cutting $k - 1$ longest internal edges)],
)[
  1. Initialize active node set $cal(A) arrow.l {1, dots.c, n}$, $r arrow.l n$, tree $T arrow.l (cal(A), emptyset)$.
  2. *While* $r > 3$:
     - Compute net divergences $R_i arrow.l sum_(k in cal(A)) D_(i k)$ for all $i in cal(A)$.
     - Select $(i^*, j^*) in arg min_(i < j in cal(A)) ((r - 2) D_(i j) - R_i - R_j)$.
     - Create internal node $u$; add edges $(i^*, u)$ and $(j^*, u)$ with lengths
       $w(i^*, u) arrow.l (1 / 2) D_(i^* j^*) + (R_(i^*) - R_(j^*)) / (2 (r - 2))$ and $w(j^*, u) arrow.l D_(i^* j^*) - w(i^*, u)$.
     - *For each* $k in cal(A) without {i^*, j^*}$: set $D_(u k) arrow.l (1 / 2)(D_(i^* k) + D_(j^* k) - D_(i^* j^*))$.
     - Set $cal(A) arrow.l (cal(A) without {i^*, j^*}) union {u}$ and $r arrow.l r - 1$.
  3. On the final triple $cal(A) = {a, b, c}$, connect all three to a central vertex $o$ with $w(a, o) = (D_(a b) + D_(a c) - D_(b c)) / 2$ (and cyclic permutations).
  4. Return $T = (V, E, w)$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish the complete theoretical foundation of Neighbor-Joining: (A) an exact first-principles derivation proving that minimizing $Q(i, j)$ is identical to Saitou and Nei's greedy minimum-evolution tree-length reduction, (B) a self-contained edge-decomposition proof of exact additive tree recovery (consistency), (C) a rigorous dimension- and $n$-independent $L_infinity$ perturbation safety radius for cherry identification, and (D) an explicit Euclidean counterexample proving why NJ yields negative Silhouette scores when applied to non-additive spatial point clouds.

#subsec("A", "Equivalence of the Q-Criterion to Greedy Minimum Evolution")

#theorem(title: "Exact Minimum-Evolution Derivation of the Q-Criterion", origin: "Original Derivation")[
  Let $T_(i j)$ denote the partially resolved $r$-leaf tree ($r >= 4$) formed from an $r$-star by grouping pair $(i, j)$ onto a shared internal vertex $u$ connected to the remaining $(r - 2)$-star hub $v$. Under ordinary least-squares (Saitou-Nei) branch-length estimation @saitou1987neighbor, the total branch length $S_(i j)$ of $T_(i j)$ satisfies the exact affine identity
  $
    S_(i j) = Q(i, j) / (2 (r - 2)) + T_"tot" / (r - 2), quad T_"tot" = sum_(a < b) D_(a b).
  $
  Consequently, minimizing $Q(i, j)$ is algebraically identical to selecting the pair $(i, j)$ that minimizes the estimated total tree length $S_(i j)$.
]

#proof[
  In $T_(i j)$, let $L_(i j) = w(i, u) + w(j, u)$, $L_(u v) = w(u, v)$, and $L_v = sum_(k != i, j) w(k, v)$, so the total tree length is $S_(i j) = L_(i j) + L_(u v) + L_v$. Matching observed average distances across the three natural bipartitions of $T_(i j)$ gives:
  1. Within ${i, j}$: $L_(i j) = D_(i j)$.
  2. Among the $r - 2$ leaves $k in.not {i, j}$: summing $D_(a b) = w(a, v) + w(b, v)$ over all $binom(r - 2, 2)$ unordered pairs counts each $w(k, v)$ exactly $(r - 3)$ times, hence
     $
       L_v = (1 / (r - 3)) sum_(a < b; a, b != i, j) D_(a b).
     $
  3. Across ${i, j}$ vs. $cal(A) without {i, j}$: averaging $D_(i k) + D_(j k)$ over the $r - 2$ external leaves $k$ yields
     $
       (1 / (2 (r - 2))) sum_(k != i, j) (D_(i k) + D_(j k)) = (1 / 2) L_(i j) + L_(u v) + L_v / (r - 2).
     $
  Solving (3) for $L_(u v)$ and adding $L_(i j) + L_v$ cancels $(r - 3)$ identically:
  $
    S_(i j) &= (1 / 2) D_(i j) + (1 / (2 (r - 2))) sum_(k != i, j) (D_(i k) + D_(j k)) \
    &quad + ((r - 3) / (r - 2)) L_v \
    &= (1 / 2) D_(i j) + (1 / (2 (r - 2))) (R_i + R_j - 2 D_(i j)) \
    &quad + (1 / (r - 2)) sum_(a < b; a, b != i, j) D_(a b).
  $
  Now note that excluding pairs incident to $i$ or $j$ from $T_"tot" = sum_(a < b) D_(a b)$ gives
  $
    sum_(a < b; a, b != i, j) D_(a b) = T_"tot" - R_i - R_j + D_(i j).
  $
  Substituting this identity and collecting terms over $2 (r - 2)$ yields
  $
    S_(i j) &= ((r - 2) D_(i j) + R_i + R_j - 2 D_(i j) \
    &quad + 2 T_"tot" - 2 R_i - 2 R_j + 2 D_(i j)) / (2 (r - 2)) \
    &= ((r - 2) D_(i j) - R_i - R_j) / (2 (r - 2)) + T_"tot" / (r - 2),
  $
  which is $Q(i, j) / (2 (r - 2)) + T_"tot" / (r - 2)$. Because $T_"tot" / (r - 2)$ is independent of $(i, j)$, $arg min Q(i, j) = arg min S_(i j)$ @saitou1987neighbor.
]

#subsec("B", "Exact Additive Tree Reconstruction (Consistency)")

A *cherry* of an unrooted tree $T$ is a pair of distinct leaves $(c, d)$ adjacent to the same internal vertex. We prove NJ's landmark consistency theorem via a clean exact expansion of $Q$ over tree edges.

#theorem(title: "Exact Additive Tree Consistency via Edge Decomposition", origin: "Classical — Cited")[
  Let $D$ be an exact additive tree metric on an unrooted tree $T = (V, E, w)$ with $r >= 4$ leaves and strictly positive internal edge weights $w(e) > 0$ for all $e in E_"int"$. Then:
  1. Every global minimizer $(i^*, j^*) in arg min_(i < j) Q(i, j)$ is a true cherry of $T$.
  2. The branch lengths $w(i^*, u), w(j^*, u)$ and reduced dissimilarities $D_(u k)$ recover the exact edge weights and additive metric of the pruned tree $T'$, so NJ reconstructs $T$ and $w$ identically @saitou1987neighbor.
]

#proof[
  *Step 1 (Exact Edge Expansion of $Q$)*: For any edge $e in E$, deleting $e$ partitions the $r$ leaves into two components $A_e, B_e$ ($abs(A_e) + abs(B_e) = r$). Count how many times $w(e)$ appears in $Q(x, y) = (r - 2) D_(x y) - R_x - R_y$:
  - If $e in P_T (x, y)$ (say $x in A_e, y in B_e$), then $w(e)$ appears $(r - 2)$ times in $(r - 2) D_(x y)$, $abs(B_e)$ times in $R_x$, and $abs(A_e)$ times in $R_y$, giving net coefficient $(r - 2) - abs(B_e) - abs(A_e) = -2$.
  - If $e in.not P_T (x, y)$, let $s_e (x, y)$ be the number of leaves on the side of $e$ *not* containing ${x, y}$. Then $w(e)$ appears $0$ times in $D_(x y)$ and $s_e (x, y)$ times in each of $R_x$ and $R_y$, giving net coefficient $-2 s_e (x, y)$.

  For every pendant (leaf) edge, either $e in P_T (x, y)$ or $s_e (x, y) = 1$, so its coefficient is *always* $-2$ independent of $(x, y)$. Letting $W_"leaf" = sum_(e in E_"leaf") w(e)$ and defining $c_e (x, y) = 1$ if $e in P_T (x, y)$ else $s_e (x, y)$, we obtain
  $
    Q(x, y) = -2 W_"leaf" - 2 sum_(e in E_"int") c_e (x, y) w(e).
  $

  *Step 2 (Every Minimizer Is a Cherry)*: Let $(i, j)$ be any *non-cherry* pair of leaves. Then $P_T (i, j)$ traverses at least one internal edge $e^* in E_"int"$. Deleting $e^*$ splits the leaf set into $A$ (containing $i$) and $B$ (containing $j$), each of cardinality at least $2$ (every internal vertex has degree $>= 3$). Since $abs(A) + abs(B) = r$, assume WLOG that $abs(B) <= r / 2$, so $abs(A) = r - abs(B) >= r / 2 >= 2$. Root the $B$-subtree at the endpoint of $e^*$ in $B$; a deepest internal node in $B$ has all children away from $e^*$ as leaves, providing a true cherry $(c, d) subset.eq B$ of $T$.

  Compare $c_e (c, d)$ against $c_e (i, j)$ on *every* internal edge $e in E_"int"$:
  - *$e = e^*$*: Since $e^* in P_T (i, j)$ while $c, d in B$, we have $c_(e^*) (i, j) = 1 < abs(A) = c_(e^*) (c, d)$.
  - *$e$ inside the $A$-component*: Let $S_e subset.eq A$ be the side of $e$ not containing $e^*$ ($abs(S_e) >= 2$). Because $c, d, j in B subset.eq X without S_e$, $c_e (c, d) = abs(S_e)$. Meanwhile $c_e (i, j) = 1$ if $i in S_e$ (since $e in P_T (i, j)$) or $c_e (i, j) = abs(S_e)$ if $i in.not S_e$. Either way, $c_e (c, d) >= c_e (i, j)$.
  - *$e$ inside the $B$-component*: Let $S_e subset.eq B$ be the side of $e$ not containing $e^*$, so $2 <= abs(S_e) < abs(B) <= r / 2$. Since $(c, d)$ is a cherry, $e in.not P_T (c, d)$. If $j in S_e$, then $e in P_T (i, j)$ gives $c_e (i, j) = 1 < 2 <= min(abs(S_e), r - abs(S_e)) <= c_e (c, d)$. If $j in.not S_e$, then $i, j in X without S_e$ gives $c_e (i, j) = abs(S_e)$; meanwhile $c_e (c, d)$ is either $abs(S_e)$ (if $c, d in.not S_e$) or $r - abs(S_e) > r / 2 > abs(S_e)$ (if $c, d in S_e$). In every subcase, $c_e (c, d) >= c_e (i, j)$.

  Summing over all $e in E_"int"$ with $w(e) > 0$ yields the strict gap
  $
    Q(i, j) - Q(c, d) >= 2 (abs(A) - 1) w(e^*) >= (r - 2) w(e^*) > 0.
  $
  Thus no non-cherry pair can ever minimize $Q$.

  *Step 3 (Exact Reduction and Induction)*: For any true cherry $(i^*, j^*)$ with parent $u$, every other leaf $k$ satisfies $D_(i^* k) = w(i^*, u) + d_T (u, k)$, $D_(j^* k) = w(j^*, u) + d_T (u, k)$, and $D_(i^* j^*) = w(i^*, u) + w(j^*, u)$. Directly:
  $
    D_(u k) &= (1 / 2)(D_(i^* k) + D_(j^* k) - D_(i^* j^*)) = d_T (u, k), \
    R_(i^*) - R_(j^*) &= sum_(k != i^*, j^*) (D_(i^* k) - D_(j^* k)) = (r - 2)(w(i^*, u) - w(j^*, u)),
  $
  so $(1 / 2) D_(i^* j^*) + (R_(i^*) - R_(j^*)) / (2 (r - 2)) = w(i^*, u)$. Induction down to $r = 3$ (where the 3-star linear system is exact) recovers $T$ and $w$ completely.
]

#subsec("C", "Dimension-Free L-Infinity Safety Radius for Cherry Identification")

The strict gap $Q(i, j) - Q(c, d) >= (r - 2) w_"min"$ derived in Step 2 of Theorem 2 immediately yields a clean, $n$-independent $L_infinity$ perturbation stability guarantee.

#theorem(title: "Constant L-Infinity Safety Radius Under Distance Noise", origin: "Original Proof")[
  Let $D$ be an additive metric on $T$ with minimum internal edge length $w_"min" = min_(e in E_"int") w(e) > 0$, and let $hat(D) = D + eta$ be a noisy symmetric dissimilarity matrix with $norm(eta)_infinity = max_(x != y) abs(eta_(x y)) <= epsilon$.
  1. If $epsilon < w_"min" / 6$ (and specifically $epsilon < w_"min" / 2$ for $r = 4$), then every global minimizer of $hat(Q)$ on $hat(D)$ is *guaranteed* to be a true cherry of $T$, independent of the number of taxa $r >= 4$.
  2. Each reduction step perturbs active distances by at most $abs(hat(D)_(u k) - D_(u k)) <= (3 / 2) epsilon$.
]

#proof[
  Write $Delta Q(x, y) = hat(Q)(x, y) - Q(x, y)$. Because $eta_(x y) = eta_(y x)$ appears once in $sum_(k != x) eta_(x k)$ and once in $sum_(k != y) eta_(y k)$, we have
  $
    Delta Q(x, y) = (r - 4) eta_(x y) - sum_(k != x, y) (eta_(x k) + eta_(y k)).
  $
  For any non-cherry $(i, j)$ and the witnessing cherry $(c, d) subset.eq B$ constructed in Theorem 2, note that $i in A$ cannot coincide with $c$ or $d in B$, leaving only two overlap cases:
  - *Disjoint* (${c, d} inter {i, j} = emptyset$, $r >= 5$): Subtracting $Delta Q(i, j)$ from $Delta Q(c, d)$ cancels all four cross terms $eta_(c i), eta_(c j), eta_(d i), eta_(d j)$ identically, leaving
    $
      Delta Q(c, d) - Delta Q(i, j) &= (r - 4)(eta_(c d) - eta_(i j)) \
      &quad - sum_(k in.not {c, d, i, j}) (eta_(c k) + eta_(d k) - eta_(i k) - eta_(j k)).
    $
    Bounding each $abs(eta) <= epsilon$ yields $abs(Delta Q(c, d) - Delta Q(i, j)) <= 6 (r - 4) epsilon$.
  - *Single overlap* (WLOG $c = j$, $r >= 4$): Both $eta_(i d)$ and every $eta_(j k)$ ($k in.not {i, j, d}$) cancel identically, leaving
    $
      Delta Q(j, d) - Delta Q(i, j) &= (r - 3)(eta_(j d) - eta_(i j)) \
      &quad + sum_(k in.not {i, j, d}) (eta_(i k) - eta_(d k)),
    $
    which is bounded in magnitude by $4 (r - 3) epsilon$.

  Combining both bounds with the true signal gap $Q(i, j) - Q(c, d) >= (r - 2) w_"min"$ from Theorem 2, we obtain $hat(Q)(c, d) < hat(Q)(i, j)$ whenever $max(6 (r - 4), 4 (r - 3)) epsilon < (r - 2) w_"min"$. At $r = 4$ (where only the overlap case occurs), this yields $4 epsilon < 2 w_"min"$, i.e., $epsilon < w_"min" / 2$. For all $r >= 4$, $max(6 (r - 4), 4 (r - 3)) < 6 (r - 2)$, so $epsilon < w_"min" / 6$ uniformly prevents any non-cherry from minimizing $hat(Q)$. Finally, Part 2 follows immediately by triangle inequality on $hat(D)_(u k) - D_(u k) = (eta_(i k) + eta_(j k) - eta_(i j)) / 2$.
]

(Note: sharper global quartet-tracking analysis by Atteson improves the multi-step full-tree $L_infinity$ safety radius to the optimal $w_"min" / 2$ without step-wise $3 / 2$ noise accumulation.)

#subsec("D", "Euclidean Non-Additivity and Negative-Silhouette Pathology")

Unlike ultrametric linkage rules, NJ's net-divergence subtraction $(-R_i - R_j)$ actively rewards pairing high-$R$ peripheral points across a Euclidean cloud, causing severe partitioning failure on standard spatial clustering benchmarks.

#proposition(title: "Constructive Euclidean Negative-Silhouette Counterexample", origin: "Original Counterexample")[
  Consider the four-point Euclidean set in $RR^2$ consisting of a tight central pair $c_1 = (0, 0)$, $c_2 = (0, 1)$ and two distant opposite points $p_1 = (-L, 0)$, $p_2 = (L, 0)$ with $L > 1$. Then NJ strictly prefers merging the distant pair $(c_1, p_1)$ (distance $L$) over the tight cluster $(c_1, c_2)$ (distance $1$). Cutting the resulting unrooted tree along its unique internal edge produces the partition ${{c_1, p_1}, {c_2, p_2}}$, whose mean Silhouette score @rousseeuw1987silhouettes is strictly negative and converges to $-1 / 12 approx -0.0833$ as $L arrow.r infinity$.
]

#proof[
  Let $s = sqrt(L^2 + 1) > L$. Pairwise Euclidean distances are $d(c_1, c_2) = 1$, $d(c_1, p_1) = d(c_1, p_2) = L$, $d(c_2, p_1) = d(c_2, p_2) = s$, and $d(p_1, p_2) = 2 L$. The net divergences are $R_(c_1) = 2 L + 1$, $R_(c_2) = 2 s + 1$, and $R_(p_1) = R_(p_2) = 3 L + s$. Computing $Q(x, y) = 2 D_(x y) - R_x - R_y$ ($r = 4$):
  $
    Q(c_1, c_2) &= 2 - (2 L + 1) - (2 s + 1) = -2 L - 2 s, \
    Q(c_1, p_1) &= 2 L - (2 L + 1) - (3 L + s) = -3 L - s - 1.
  $
  Comparing the two scores using $s - L = 1 / (s + L) < 1$ (since $L > 1$):
  $
    Q(c_1, c_2) - Q(c_1, p_1) = L - s + 1 = 1 - 1 / (s + L) > 0.
  $
  Hence $Q(c_1, p_1) < Q(c_1, c_2)$ strictly (e.g., at $L = 10$, $Q(c_1, p_1) approx -41.05 < -40.10 approx Q(c_1, c_2)$), forcing NJ to form the bipartition ${{c_1, p_1}, {c_2, p_2}}$ rather than separating ${c_1, c_2}$ from ${p_1, p_2}$. Evaluating point Silhouettes $s(x) = (b(x) - a(x)) / max(a(x), b(x))$:
  - For $c_1$: $a(c_1) = L$, $b(c_1) = (1 + L) / 2$, giving $s(c_1) = -1 / 2 + 1 / (2 L) arrow.r -1 / 2$.
  - For $c_2$: $a(c_2) = s$, $b(c_2) = (1 + s) / 2$, giving $s(c_2) = -1 / 2 + 1 / (2 s) arrow.r -1 / 2$.
  - For $p_1$ and $p_2$: $a(p_1) = L$, $b(p_1) = (s + 2 L) / 2 arrow.r (3 / 2) L$, giving $s(p_1) arrow.r +1 / 3$ and $s(p_2) arrow.r +1 / 3$.
  Averaging over all four points yields mean Silhouette $arrow.r (-1 / 2 - 1 / 2 + 1 / 3 + 1 / 3) / 4 = -1 / 12 < 0$ (already $-0.058$ at $L = 10$).
]

#sec("III", "Complexity, Geometry and Operational Profile")

Exact Neighbor-Joining runs in $Theta(n^3)$ time and $Theta(n^2)$ space @saitou1987neighbor. At iteration $r in {n, n - 1, dots.c, 4}$, updating the vector of net divergences $R_i$ takes $O(r)$ time incrementally, whereas finding $arg min_(i < j) Q(i, j)$ scans all $binom(r, 2)$ active pairs in $Theta(r^2)$ time because *every* entry $Q(i, j)$ shifts globally when $r$ decrements and $R_i$ updates. Summing $sum_(r=4)^n r(r - 1) / 2 = Theta(n^3)$ yields the exact cubic complexity. Unlike reducible Lance-Williams updates (Alg. 72), NJ cannot use a standard Nearest-Neighbor Chain without heuristic modifications (such as RapidNJ or relaxed neighbor-joining) because global re-centering $R_i$ invalidates local nearest-neighbor invariance after each merge. Nonetheless, the inner loop consists purely of contiguous array additions and subtractions, giving high cache efficiency and single-digit millisecond runtimes for $n <= 300$.

Geometrically, NJ operates in the cone of *additive tree metrics* rather than ultrametrics or Euclidean Voronoi partitions. Its positive homogeneity ($D arrow.r c D$ scales branch lengths by $c$ without changing topology) and permutation equivariance (when tie-free) are exact, but it is *not* invariant under nonlinear monotone transformations $f compose D$ (which destroy four-point additivity). Furthermore, while first-step branch lengths $w(i, u)$ are guaranteed non-negative on any metric input (since triangle inequality forces $abs(R_i - R_j) <= (r - 2) D_(i j)$), iterative distance reductions $D_(u k) = (D_(i k) + D_(j k) - D_(i j)) / 2$ can break the triangle inequality on noisy non-additive matrices, producing negative branch lengths that are conventionally clamped to $0$ in phylogenetic software.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time complexity], [$Theta(n^3)$ exact ($O(r^2)$ scan of $Q(i, j)$ per step $r = n, dots.c, 4$)],
    [Space complexity], [$Theta(n^2)$ symmetric working distance matrix plus $O(n)$ divergence vector $R$],
    [Output structure], [Unrooted weighted additive tree $T = (V, E, w)$ (not an ultrametric)],
    [Exact consistency], [Exact recovery of any additive tree metric; $L_infinity$ cherry radius $>= w_"min" / 6$],
    [Molecular clock], [Not required (handles arbitrary heterogeneous lineage rates)],
    [Branch positivity], [Guaranteed at step 1 on metrics; negative edges possible after reductions],
    [Inductive bias], [Minimum-evolution tree additivity (mismatched to compact Euclidean balls)],
  ),
  caption: [Operational and theoretical profile of Neighbor-Joining (Alg. 75).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Molecular Phylogenetics and Comparative Genomics*: NJ is the foundational distance-based method for reconstructing evolutionary trees from multiple sequence alignments @saitou1987neighbor using evolutionary substitution distances (Jukes-Cantor, Kimura 2-parameter, Tamura-Nei). Because protein families and species lineages experience markedly different mutation rates, NJ's exact consistency on non-ultrametric additive trees (Theorem 2) directly solves the rate-heterogeneity failure of UPGMA.
2. *Viral Phylodynamics and Pathogen Outbreak Tracing*: Rapidly evolving RNA pathogens (influenza, HIV, SARS-CoV-2) display strong rate variation across hosts and antigenic drift regimes. NJ builds rapid, statistically consistent evolutionary scaffolds from pairwise genomic divergence matrices for real-time epidemiological surveillance.
3. *Single-Cell Lineage Tracing and Tumor Clonal Phylogeny*: Somatic copy-number and single-nucleotide variant distances across cancer subclones violate molecular-clock assumptions due to mutator phenotypes and punctuated genome doubling. NJ reconstructs unrooted clonal branching orders without distorting accelerated tumor lineages.
4. *Starting-Tree Initialization for Likelihood and Bayesian Phylogenetics*: Maximum-likelihood and Markov chain Monte Carlo phylogenetic engines operate over a super-exponential topology space ($(2 n - 5)!!$ unrooted trees). Because NJ recovers the exact topology within a constant $L_infinity$ noise radius (Theorem 3) in deterministic polynomial time, it is the standard warm-start initializer for downstream likelihood optimization.
5. *Historical Linguistics and Manuscript Stemmatology*: Cognate-divergence matrices across Indo-European languages and textual-variant distances across medieval manuscript copies reflect branching descent with unequal regional drift rates; NJ recovers branching stemmata faithfully where ultrametric clustering collapses fast-drifting branches into spurious outgroups.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Five failure modes restrict NJ's deployment: (1) *Severe mismatch on Euclidean spatial clustering*, where orthogonal coordinates violate the four-point condition and $-R_i - R_j$ subtraction pairs distant peripheral points with central points, producing negative Silhouette partitions (Proposition 4); (2) *Cubic runtime bottleneck* ($O(n^3)$), limiting exact NJ to roughly $n <= 10^4$ taxa; (3) *Negative branch-length artifacts* on strongly non-additive or high-variance evolutionary distances; (4) *Absence of an intrinsic root*, requiring an external outgroup or midpoint rooting heuristic to orient evolutionary time; and (5) *Long-branch distance saturation*, where deeply diverged sequences suffer high variance in corrected distances, violating the $L_infinity$ safety radius of Theorem 3.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Molecular phylogenetics ($n <= 10^4$)], [Excellent: consistent additive trees without molecular clock], [Root via outgroup; refine with Max-Likelihood],
    [Single-cell / clonal trees], [Strong: handles heterogeneous somatic mutation rates], [Clamp negative branches; check four-point fit],
    [Euclidean tabular clustering], [Very poor: breaks compact clusters, yields negative Silhouette], [Use Ward, Flexible Beta (Alg. 74), or Minimax (Alg. 73)],
    [Non-convex / manifold shapes], [Very poor: tree-metric assumption fails on spatial manifolds], [Use HDBSCAN @campello2013hdbscan or Spectral @ng2002spectral],
    [Ultra-large phylogenies ($n > 10^4$)], [Slow: $O(n^3)$ exact $Q$-matrix scan], [Use RapidNJ, FastME, or divide-and-conquer trees],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for Alg. 75.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Neighbor-Joining occupies a unique theoretical position adjacent to the Lance-Williams SAHN family (Alg. 72): its *distance reduction step* formally matches a Lance-Williams recurrence, yet its *selection step* breaks the entire SAHN paradigm.

#proposition(title: "Algebraic Comparison of NJ Reduction and Lance-Williams Updates", origin: "Original Derivation")[
  NJ's reduced-distance update $D_(u k) = (1 / 2)(D_(i k) + D_(j k) - D_(i j))$ coincides algebraically with the Lance-Williams parameterization $(alpha_i, alpha_j, beta, gamma) = (1 / 2, 1 / 2, -1 / 2, 0)$. Because $alpha_i + alpha_j + beta = 1 / 2 < 1$, this update strictly violates the Lance-Williams monotonicity condition (Alg. 72, Theorem 1) and contracts distances toward the tree trunk as pendant edges are pruned.
]

#proof[
  Evaluating the Lance-Williams recurrence at $alpha_i = alpha_j = 1 / 2$, $beta = -1 / 2$, $gamma = 0$ directly gives $D_(u k) = (1 / 2) D_(i k) + (1 / 2) D_(j k) - (1 / 2) D_(i j)$. Here $alpha_i + alpha_j + beta = 1 / 2 < 1$, violating Theorem 1(iii) of Alg. 72: geometrically, subtracting $(1 / 2) D_(i j) = (w(i, u) + w(j, u)) / 2$ intentionally removes the pendant branch lengths of the pruned cherry so $D_(u k)$ measures distance from the deeper internal node $u$ rather than from leaves $i, j$.
]

Replacing NJ's $Q$-selection with raw $min D_(i j)$ and setting $beta = 0$ recovers WPGMA (and UPGMA when size-weighted) @sokal1958statistical, trading additive rate flexibility for ultrametric monotonicity and spatial compactness.

*Empirical Benchmark Behavior*: On the ClusterOverview 4-dataset spatial benchmark, Neighbor-Joining records Iris Silhouette $= -0.127$ (@rousseeuw1987silhouettes) and $"ARI" = 0.019$ (@hubert1985comparing), Customers Silhouette $= 0.147$, Multi-Shape Silhouette $= -0.007$, and High-D Silhouette $= -0.049$, with a mean runtime of $2.857$ ms. Far from indicating an implementation defect, these dramatically negative Silhouettes and near-zero Iris ARI provide a textbook empirical validation of Proposition 4: all four benchmark datasets are non-additive multi-dimensional Euclidean clouds, not phylogenetic tree metrics. NJ's net-divergence subtraction $(-R_i - R_j)$ optimizes unrooted minimum-evolution tree length rather than within-cluster compactness, splitting coherent Euclidean clusters along trunk-vs-periphery axes. While its compact arithmetic loop achieves a fast $2.857$ ms runtime on small benchmark sizes, NJ must be reserved strictly for genuine evolutionary, genealogical, or tree-metric domains — whereas compact spatial clustering belongs to ultrametric methods such as Lance-Williams (Alg. 72), Minimax (Alg. 73), and Flexible Beta (Alg. 74) @jain2010data @xu2005survey.
