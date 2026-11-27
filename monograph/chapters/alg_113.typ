#import "../template/ieee.typ": *

#alg_entry(
  id: 113,
  name: "Infomap",
  category: "Cat. V: Graph, Spectral & Network",
  year: "2008",
  paradigm: "Information-Theoretic Flow Compression (Map Equation)",
  complexity: $cal(O)(m log n)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Infomap, introduced by Rosvall and Bergstrom @rosvall2008infomap,
recasts community detection as a Minimum Description Length (MDL)
source-coding problem @wallace1968snob on a random walk. Rather than
maximizing static edge density against a configuration null model,
Infomap defines communities as regions of a graph that *trap*
probability flow for long intervals, enabling a two-level Huffman
codebook to compress an infinite random-walk trajectory into fewer
bits per step.

Let $G = (V, E, w)$ be a finite weighted graph (directed or
undirected) with $n = abs(V)$ vertices, $m = abs(E)$ edges, and
weights $w_(u v) >= 0$ with out-strengths
$d_u^"out" = sum_(v in V) w_(u v)$. Consider a Markov chain
$(X_t)_(t >= 0)$ on $V$ with transition matrix $P in [0, 1]^(n times n)$.
For a connected undirected graph, the unregularized walk uses
$P_(u v) = w_(u v) / d_u$ and possesses the closed-form stationary
distribution $pi_u = d_u / "vol"(V)$, where
$"vol"(V) = sum_(v in V) d_v$. For general directed or disconnected
graphs, teleportation with damping rate $tau in (0, 1)$ (typically
$tau = 0.15$) produces the ergodic PageRank kernel
$
  P_(u v) = (1 - tau) (w_(u v)) / (d_u^"out") + (tau) / (n)
$
(teleporting uniformly from dangling nodes with $d_u^"out" = 0$),
which admits a unique strictly positive stationary distribution
$pi P = pi$ ($pi_u > 0$, $sum_(u in V) pi_u = 1$).

Let $M = {M_1, dots, M_K}$ be a hard partition of $V$ into $K$
non-empty modules. For each module $i in {1, dots, K}$, define its
stationary occupancy $pi_i = sum_(u in M_i) pi_u$, its *exit flow*
$q_i^"out" = sum_(u in M_i) sum_(v in.not M_i) pi_u P_(u v)$, and its
*entry flow* $q_i^"in" = sum_(u in.not M_i) sum_(v in M_i) pi_u P_(u v)$.
By stationarity ($sum_u pi_u P_(u v) = pi_v$ and
$sum_v P_(u v) = 1$), exact flow conservation holds on both
undirected and directed graphs:
$
  pi_i - q_i^"in"
    = sum_(u, v in M_i) pi_u P_(u v)
    = pi_i - q_i^"out" quad => quad q_i^"in" = q_i^"out".
$
Let $q_"out" = sum_(i = 1)^K q_i^"out"$ be the total inter-module
transition rate and $p_i^"mod" = q_i^"out" + pi_i$ be the total
codebook usage rate of module $i$.

#definition(title: "Two-Level Map Equation")[
  The *two-level Map Equation* @rosvall2008infomap of a partition
  $M = {M_1, dots, M_K}$ under ergodic Markov flow $(P, pi)$ is the
  theoretical lower bound on the per-step description length (in
  bits) of a two-level prefix code consisting of one global *index
  codebook* $cal(Q)$ and $K$ local *module codebooks*
  ${cal(P)^i}_(i = 1)^K$:
  $
    L(M) = q_"out" H(cal(Q)) + sum_(i = 1)^K p_i^"mod" H(cal(P)^i),
  $
  where (with $0 log_2 0 := 0$)
  $
    H(cal(Q)) &= - sum_(i = 1)^K (q_i^"out") / (q_"out")
      log_2 (q_i^"out") / (q_"out"), \
    H(cal(P)^i) &= -(q_i^"out") / (p_i^"mod")
      log_2 (q_i^"out") / (p_i^"mod")
      - sum_(u in M_i) (pi_u) / (p_i^"mod")
      log_2 (pi_u) / (p_i^"mod").
  $
]

Intuitively, $cal(Q)$ names module entries whenever the walker jumps
between modules (at rate $q_"out"$), while $cal(P)^i$ reuses short
local codewords to encode both within-module node visits ($pi_u$) and
the special module-exit codeword ($q_i^"out"$).

#algbox(
  title: "Infomap Multilevel Map-Equation Optimizer",
  inputs: [Weighted graph $G = (V, E, w)$, teleportation $tau$
    ($tau = 0$ if connected undirected), tolerance $epsilon > 0$.],
  outputs: [Flow-compressing partition $M^* = {M_1, dots, M_K}$ and
    code length $L(M^*)$.],
)[
  1. *Stationary Flow:* Compute $pi$ analytically ($pi_u = d_u /
     "vol"(V)$ if undirected, $tau = 0$) or via power iteration; set
     edge flows $f_(u v) = pi_u P_(u v)$.
  2. *Initialize:* Assign each vertex to its own singleton module
     $M_u = {u}$ and maintain per-module scalars $(pi_i, q_i^"out")$
     and global $q_"out"$.
  3. *Local Greedy Sweeps:* Repeatedly visit vertices in random
     order; for each $u in M_a$, tally flow to neighboring modules
     $M_b$ and move $u -> M_b$ maximizing $-Delta L(u : a -> b) > 0$
     via Theorem 1 until no move decreases $L(M)$ by more than
     $epsilon$.
  4. *Multilevel Coarsening:* Contract each module $M_i$ into a
     super-node with aggregated inter-module flows
     $F_(i j) = sum_(u in M_i, v in M_j) f_(u v)$ and recurse Step 3
     on the coarse network until $L(M)$ stabilizes, followed by
     fine-level refinement sweeps @blondel2008louvain.
  5. *Return:* Best partition $M^*$ and code length $L(M^*)$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We present five formal mathematical results for Infomap: the exact
algebraic expansion enabling $cal(O)(1)$ incremental move scoring,
the Shannon source-coding optimality and singleton-boundary theorem,
an original KL-divergence decomposition proving strict zero-overfitting
on memoryless flow, an exact closed-form resolution threshold on
rings of cliques, and a constructive proof of directed-orientation
sensitivity.

#subsec("A", "Expanded Entropy Form and O(1) Local Move Scoring")

#theorem(
  title: "Exact Algebraic Expansion and Local Delta-L Decomposition",
  origin: "Original Derivation",
)[
  Define $h(x) = x log_2 x$ for $x > 0$ and $h(0) = 0$. For any
  ergodic Markov chain $(P, pi)$ and any partition
  $M = {M_1, dots, M_K}$, the Map Equation admits the exact expanded
  identity
  $
    L(M) = h(q_"out") - 2 sum_(i = 1)^K h(q_i^"out")
      - sum_(u in V) h(pi_u) + sum_(i = 1)^K h(p_i^"mod").
  $
  Consequently, $-sum_(u in V) h(pi_u) = H(pi)$ is invariant across
  all partitions, and moving a vertex $u$ from module $a$ to module
  $b$ changes $L(M)$ by
  $
    Delta L(u : a -> b) &= Delta h(q_"out")
      - 2 Delta h(q_a^"out") - 2 Delta h(q_b^"out") \
      &quad + Delta h(p_a^"mod") + Delta h(p_b^"mod"),
  $
  which is evaluable in $cal(O)(1)$ arithmetic operations once
  $u$'s flow into $M_a$ and $M_b$ is known.
]

#proof[
  Expand the two terms of Definition 1 using
  $log_2 (x / y) = log_2 x - log_2 y$. First, for the index codebook:
  $
    q_"out" H(cal(Q)) &= - sum_(i = 1)^K q_i^"out"
      (log_2 q_i^"out" - log_2 q_"out") \
      &= - sum_(i = 1)^K h(q_i^"out")
      + (sum_(i = 1)^K q_i^"out") log_2 q_"out" \
      &= h(q_"out") - sum_(i = 1)^K h(q_i^"out").
  $
  Second, for each module codebook $i in {1, dots, K}$, using
  $p_i^"mod" = q_i^"out" + sum_(u in M_i) pi_u$:
  $
    p_i^"mod" H(cal(P)^i) &= -q_i^"out"
      (log_2 q_i^"out" - log_2 p_i^"mod") \
      &quad - sum_(u in M_i) pi_u (log_2 pi_u - log_2 p_i^"mod") \
      &= -h(q_i^"out") - sum_(u in M_i) h(pi_u) + h(p_i^"mod").
  $
  Summing $p_i^"mod" H(cal(P)^i)$ over $i = 1, dots, K$ (noting that
  ${M_i}$ partitions $V$) and adding $q_"out" H(cal(Q))$ gives the
  claimed expanded identity. Because $pi_u$ depends solely on the
  Markov chain $(P, pi)$ and not on $M$, the term
  $-sum_(u in V) h(pi_u) = H(pi)$ cancels identically in
  $L(M') - L(M)$. Moving $u$ from $M_a$ to $M_b$ leaves
  $(q_i^"out", p_i^"mod")$ unchanged for every $i in.not {a, b}$,
  altering only $q_"out"$ and the four scalars
  $(q_a^"out", q_b^"out", p_a^"mod", p_b^"mod")$, yielding the five-term
  $cal(O)(1)$ formula for $Delta L(u : a -> b)$.
]

#subsec("B", "Shannon Coding Bounds and Boundary Evaluations")

#theorem(
  title: "Shannon Prefix-Code Optimality and Boundary Code Lengths",
  origin: "Classical — Cited",
)[
  Following @rosvall2008infomap:
  (i) For any two-level prefix-free code matching partition $M$, the
  expected per-step codeword length $bar(ell)$ satisfies
  $bar(ell) >= L(M)$, while symbol-by-symbol Huffman coding achieves
  $L(M) <= bar(ell)_"Huff" < L(M) + 1 + 2 q_"out"$.
  (ii) For the trivial 1-module partition $M_1 = {V}$,
  $L(M_1) = H(pi)$; for the $n$-singleton partition
  $M_n = {{v} : v in V}$ on a loopless walk ($P_(u u) = 0$),
  $L(M_n) = H(pi) + 2$.
]

#proof[
  For (i), each step of the stationary walk triggers an emission from
  codebook $cal(Q)$ with probability $q_"out"$ and from module
  codebook $cal(P)^i$ with probability $p_i^"mod"$ (one emission per
  node visit in $M_i$, rate $pi_i$, plus one exit emission per exit,
  rate $q_i^"out"$). By Shannon's noiseless coding theorem, any
  prefix-free code on distribution $cal(D)$ has expected length
  $bar(ell)_cal(D) >= H(cal(D))$, and Huffman coding satisfies
  $bar(ell)_cal(D) < H(cal(D)) + 1$. Taking the usage-weighted sum
  gives $L(M) <= bar(ell)_"Huff" < L(M) + q_"out" + sum_(i = 1)^K
  p_i^"mod"$. Since $sum_i p_i^"mod" = sum_i (pi_i + q_i^"out") =
  1 + q_"out"$, the upper slack is strictly bounded by
  $1 + 2 q_"out"$.

  For (ii), when $M_1 = {V}$, $q_1^"out" = q_"out" = 0$ and
  $p_1^"mod" = 1$, so Theorem 1 gives
  $L(M_1) = 0 - 0 + H(pi) + h(1) = H(pi)$. When $M_n = {{v}}$ and
  $P_(u u) = 0$, every visit to $u$ immediately exits ${u}$, so
  $q_u^"out" = pi_u$, $q_"out" = 1$, and $p_u^"mod" = 2 pi_u$. Using
  $h(2 pi_u) = 2 pi_u + 2 h(pi_u)$ in Theorem 1 yields
  $
    L(M_n) &= h(1) - 2 sum_(u in V) h(pi_u) - sum_(u in V) h(pi_u)
      + sum_(u in V) (2 pi_u + 2 h(pi_u)) \
      &= 2 sum_(u in V) pi_u - sum_(u in V) h(pi_u) = H(pi) + 2.
  $
]

#subsec("C", "Exact KL Decomposition and Zero Overfitting on Null Flow")

Unlike modularity, which finds high-$Q$ spurious partitions even on
complete or unstructured random graphs, Infomap admits an exact
information-theoretic anti-overfitting guarantee.

#theorem(
  title: "Exact Sub-Probability KL Decomposition and Null Immunity",
  origin: "Original Proof",
)[
  Let $(P, pi)$ be any irreducible Markov chain on $V$ with entropy
  rate $H(cal(X)) = -sum_(u in V) pi_u sum_(v in V) P_(u v) log_2
  P_(u v)$. For any partition $M = {M_1, dots, M_K}$, define for
  $u in M_i$ the two-level code-induced kernel
  $
    tilde(P)_M (v | u) = cases(
      pi_v / p_i^"mod" & "if" v in M_i,
      (q_i^"out" / p_i^"mod") (q_j^"out" / q_"out") (pi_v / p_j^"mod")
        & "if" v in M_j\, j != i,
    )
  $
  its row sum $Z_i = sum_(v in V) tilde(P)_M (v | u) in (0, 1]$, and
  normalized kernel $hat(P)_M (v | u) = tilde(P)_M (v | u) / Z_i$.
  Then:
  (i) Exact divergence identity:
  $
    L(M) = H(cal(X)) + sum_(i = 1)^K pi_i log_2 (1) / (Z_i)
      + sum_(u in V) pi_u D_"KL" (P(dot.op | u) || hat(P)_M (dot.op | u)).
  $
  (ii) Strict universal bound: $L(M) >= H(cal(X))$, with strict
  inequality $L(M) > H(cal(X))$ whenever $K >= 2$.
  (iii) Exact null-model immunity: if the walk is memoryless
  ($P_(u v) = pi_v$ for all $u, v in V$), then $M_1 = {V}$ is the
  *strict unique global minimizer* of $L(M)$: $L(M) > L(M_1)$ for
  every partition with $K >= 2$.
]

#proof[
  Compute the cross-entropy of $P$ against the sub-probability kernel
  $tilde(P)_M$:
  $
    E_M := - sum_(u in V) pi_u sum_(v in V) P_(u v)
      log_2 tilde(P)_M (v | u).
  $
  Splitting the inner sum into $v in M_i$ and $v in M_j$ ($j != i$),
  the term $-log_2 (pi_v / p_j^"mod")$ appears for every $(u, v)$
  with $v in M_j$, weighted by $sum_(u in V) pi_u P_(u v) = pi_v$;
  the term $-log_2 (q_i^"out" / p_i^"mod")$ appears for all
  $u in M_i, v in.not M_i$, weighted by $q_i^"out"$; and the term
  $-log_2 (q_j^"out" / q_"out")$ appears for all
  $u in.not M_j, v in M_j$, weighted by $q_j^"in" = q_j^"out"$.
  Collecting all three sums yields exactly
  $
    E_M = sum_(i = 1)^K p_i^"mod" H(cal(P)^i) + q_"out" H(cal(Q))
      = L(M).
  $
  Subtracting $H(cal(X)) = -sum_(u, v) pi_u P_(u v) log_2 P_(u v)$
  and substituting $tilde(P)_M (v | u) = Z_i hat(P)_M (v | u)$ gives
  $
    L(M) - H(cal(X)) &= sum_(i = 1)^K sum_(u in M_i) pi_u sum_(v in V)
      P_(u v) log_2 (P_(u v)) / (Z_i hat(P)_M (v | u)) \
      &= sum_(i = 1)^K pi_i log_2 (1) / (Z_i)
      + sum_(u in V) pi_u D_"KL" (P(dot.op | u) || hat(P)_M (dot.op | u)),
  $
  proving (i). For (ii), $D_"KL" >= 0$ by Gibbs' inequality, and for
  every $i$,
  $
    Z_i = (pi_i) / (p_i^"mod") + (q_i^"out") / (p_i^"mod")
      sum_(j != i) (q_j^"out") / (q_"out") (pi_j) / (p_j^"mod").
  $
  When $K >= 2$ and $P$ is irreducible, $q_j^"out" > 0$ for all $j$,
  so $pi_j / p_j^"mod" = pi_j / (pi_j + q_j^"out") < 1$ and
  $sum_(j != i) q_j^"out" / q_"out" < 1$, implying strictly
  $Z_i < (pi_i + q_i^"out") / p_i^"mod" = 1$, hence
  $log_2 (1 / Z_i) > 0$ and $L(M) > H(cal(X))$.
  For (iii), when $P_(u v) = pi_v$, $H(cal(X)) = H(pi) = L(M_1)$ by
  Theorem 2(ii); therefore (ii) immediately enforces
  $L(M) > H(pi) = L(M_1)$ for every partition with $K >= 2$.
]

#subsec("D", "Exact Ring-of-Cliques Resolution Threshold")

Although Infomap resists modularity's square-root resolution limit
@fortunato2007resolution, two-level compression still incurs a global
index-codebook cost that induces an *exponentially weaker*
resolution threshold on rings of cliques.

#theorem(
  title: "Exact Exponential Resolution Threshold of Two-Level Infomap",
  origin: "Original Derivation",
)[
  Let $G_(r, c)$ be the ring of $r$ identical cliques $K_c$
  ($c >= 3$, $r$ even) connected cyclically by single bridges, with
  $ell = binom(c, 2) + 1$ edges per clique ($m = r ell$) and
  unregularized undirected walk ($tau = 0$). Comparing the true
  $r$-clique partition $M_r$ against the pairwise-merged partition
  $M_(r / 2)$, Infomap strictly prefers merging pairs
  ($L(M_(r / 2)) < L(M_r)$) if and only if
  $
    r > r^* (ell) := ((2 ell + 1)^(2 ell + 1)) /
      (2 (ell + 1)^(2 ell + 2))
      tilde.op (2^(2 ell)) / (e dot.op ell).
  $
]

#proof[
  In $G_(r, c)$, total volume is $2 m = 2 r ell$, each clique has
  volume $2 ell$ ($pi_i = 1 / r$), and each clique has two outgoing
  bridges, so $q_i^"out" = 2 / (2 r ell) = 1 / (r ell)$,
  $q_"out" = 1 / ell$, and $p_i^"mod" = (ell + 1) / (r ell)$.
  Evaluating Theorem 1 on $M_r$ gives
  $
    L(M_r) - H(pi) &= -(1) / (ell) log_2 ell
      + (2) / (ell) log_2 (r ell) \
      &quad + (ell + 1) / (ell) log_2 ((ell + 1) / (r ell)) \
      &= (1 - ell) / (ell) log_2 r - log_2 ell
      + (ell + 1) / (ell) log_2 (ell + 1).
  $
  In the pairwise-merged partition $M_(r / 2)$ ($r / 2$ modules of
  two adjacent cliques each), each super-module still has exactly two
  outgoing bridges, giving $pi'_j = 2 / r$,
  $q'_j^"out" = 1 / (r ell)$, $q'_"out" = 1 / (2 ell)$, and
  $p'_j^"mod" = (2 ell + 1) / (r ell)$. Evaluating Theorem 1 again:
  $
    L(M_(r / 2)) - H(pi) &= (1 - 2 ell) / (2 ell) log_2 r
      - (1) / (2 ell) - log_2 ell \
      &quad + (2 ell + 1) / (2 ell) log_2 (2 ell + 1).
  $
  Subtracting the two expressions cancels $H(pi)$ and $-log_2 ell$:
  $
    2 ell (L(M_r) - L(M_(r / 2)))
      &= log_2 r + 1 + 2(ell + 1) log_2 (ell + 1) \
      &quad - (2 ell + 1) log_2 (2 ell + 1).
  $
  Thus $L(M_(r / 2)) < L(M_r)$ holds if and only if
  $
    log_2 r > (2 ell + 1) log_2 (2 ell + 1)
      - 2(ell + 1) log_2 (ell + 1) - 1,
  $
  which exponentiates to $r > (2 ell + 1)^(2 ell + 1) /
  (2 (ell + 1)^(2 ell + 2)) = r^*(ell)$. Rewriting
  $r^*(ell) = (2^(2 ell) / (ell + 1)) dot.op
  (1 + 1 / (2 ell))^(2 ell + 1) / (1 + 1 / ell)^(2 ell + 2)$ and
  taking $ell -> infinity$ yields $(e / e^2) dot.op 2^(2 ell) / ell
  = 2^(2 ell) / (e ell)$.
]

Theorem 4 quantifies Infomap's dramatic resolution advantage over
modularity (Alg. 112, Theorem 3): for a ring of $K_5$ cliques
($ell = 11$), modularity merges pairs at $r > 22$, whereas two-level
Infomap resolves all cliques intact up to $r^* (11) approx 1.31
times 10^5$—an exponential separation (and hierarchical multilevel
Infomap removes even this two-level threshold).

#subsec("E", "Directed Orientation Sensitivity and Scale Invariance")

#proposition(
  title: "Exact Orientation Sensitivity on Identical Undirected Skeletons",
  origin: "Original Counterexample",
)[
  (a) On connected undirected graphs ($tau = 0$), $L(M)$ is strictly
  invariant under positive weight scaling $w |-> alpha w$
  ($alpha > 0$) and vertex relabeling.
  (b) There exist two strongly connected unweighted directed graphs
  $G_A, G_B$ sharing the exact same 6-vertex, 8-edge undirected
  skeleton such that Infomap ($tau = 0$) strictly partitions $G_A$
  into $C_1 = {1, 2, 3}, C_2 = {4, 5, 6}$ ($L_A ({C_1, C_2}) <
  L_A (M_1)$) while strictly keeping $G_B$ undivided
  ($L_B ({C_1, C_2}) > L_B (M_1)$).
]

#proof[
  Part (a) follows immediately because $P_(u v) = w_(u v) / d_u$ and
  $pi_u = d_u / "vol"(V)$ are homogeneous of degree zero in $w$.

  For (b), take the 6-vertex skeleton consisting of the 6-cycle
  $1 - 2 - 3 - 4 - 5 - 6 - 1$ plus chords ${1, 3}$ and ${4, 6}$.
  In $G_A$, orient the outer ring cyclically
  $1 -> 2 -> 3 -> 4 -> 5 -> 6 -> 1$ and orient the chords *backward*
  as $3 -> 1$ and $6 -> 4$ (forming two trapping directed 3-cycles).
  Out-degrees are $d^+(1) = d^+(2) = d^+(4) = d^+(5) = 1$ and
  $d^+(3) = d^+(6) = 2$; direct substitution verifies the uniform
  stationary law $pi_u = 1 / 6$, so $L_A (M_1) = log_2 6 approx
  2.585$ bits. For $M = {C_1, C_2}$, $q_1^"out" = q_2^"out" =
  (1 / 6)(1 / 2) = 1 / 12$, $q_"out" = 1 / 6$, and
  $p_1^"mod" = p_2^"mod" = 7 / 12$, yielding
  $
    L_A ({C_1, C_2}) = (7) / (6) log_2 7 - (5) / (6)
      approx 2.442 < 2.585 = L_A (M_1).
  $
  In $G_B$, keep the exact same 8 undirected edges but flip the two
  chords *forward* to $1 -> 3$ and $4 -> 6$. Now $d^+(1) = d^+(4) = 2$
  and all other out-degrees are $1$. Solving $pi P = pi$ yields
  $pi_1 = pi_3 = pi_4 = pi_6 = 1 / 5$ and $pi_2 = pi_5 = 1 / 10$, so
  $L_B (M_1) = H(pi) = log_2 5 + 1 / 5 approx 2.522$ bits. For
  ${C_1, C_2}$, $q_1^"out" = q_2^"out" = pi_3 dot.op 1 = 1 / 5$,
  $q_"out" = 2 / 5$, and $p_1^"mod" = p_2^"mod" = 7 / 10$, giving
  $
    L_B ({C_1, C_2}) = (7) / (5) log_2 7 - (4) / (5)
      approx 3.130 > 2.522 = L_B (M_1).
  $
  Undirected modularity and spectral cuts assign identical partitions
  to $G_A$ and $G_B$, whereas Infomap cleanly separates cyclic flow
  traps ($G_A$) from feedforward flow chains ($G_B$).
]

#sec("III", "Complexity, Geometry and Operational Profile")

Exact global minimization of $L(M)$ over all partitions is NP-hard,
so Infomap employs a Louvain-style multilevel greedy heuristic
@blondel2008louvain. By Theorem 1, evaluating a candidate move of
vertex $u$ requires tallying its incident edge flows to neighboring
modules in $cal(O)("deg"(u))$ time and scoring five scalar $h(dot.op)$
differences in $cal(O)(1)$ time per candidate module. One sweep over
$V$ takes $cal(O)(m)$ time; across $cal(O)(log n)$ coarsening levels,
runtime scales as $cal(O)(m log n)$ with $cal(O)(n + m)$ space. Its
constant factor exceeds Chinese Whispers and Louvain due to per-move
`log2` evaluations and multilevel flow-table maintenance.

Geometrically, Infomap partitions graphs into *metastable Markov
basins* rather than Euclidean balls. While avoiding modularity's
square-root resolution limit, two-level Infomap exhibits a
*field-of-view limit* on high-diameter lattices and chains: because a
walker requires $Theta(L^2)$ steps to cross diameter $L$, Infomap
over-segments extended geometric ribbons into compact diffusive tiles.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time complexity],
    [$cal(O)(m log n)$ empirical multilevel greedy ($cal(O)(m)$ per
      sweep via Thm. 1); higher constant factor due to `log2` calls.],
    [Space complexity],
    [$cal(O)(n + m)$ for weighted edge flows, module scalars
      $(pi_i, q_i^"out")$, and coarse hierarchy levels.],
    [Hyperparameters],
    [None required on connected undirected graphs; teleportation
      $tau in (0, 1)$ (default $0.15$) for directed graphs.],
    [Determinism],
    [Stochastic vertex-visitation order in greedy sweeps; standard
      implementations retain best $L(M)$ across $R$ trials.],
    [Directed support],
    [Native first-class support; distinguishes cyclic traps from
      feedforward flow (Prop. 1b).],
    [Null-model behavior],
    [Strictly refuses to split memoryless uniform flow (Thm. 3iii);
      strong protection against noise overfitting.],
    [Resolution profile],
    [Exponentially weaker ring resolution threshold
      $r^*(ell) tilde.op 2^(2 ell) / (e ell)$ (Thm. 4); subject to
      field-of-view over-splitting on high-diameter grids.],
  ),
  caption: [Operational and theoretical profile of Infomap.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

*1. Scientific Citation Networks and Knowledge Mapping:*
In Rosvall and Bergstrom's flagship study @rosvall2008infomap,
weighted directed citation graphs trace how scholarly attention
circulates within disciplines; Infomap's orientation-sensitive flow
objective (Proposition 1b) isolates disciplinary traps cleanly.

*2. Ecological Food Webs and Biogeography:*
Species-occurrence and ocean-current larval transport graphs encode
asymmetric directed flows, allowing Infomap to extract
biogeographical provinces as flow-retaining basins directly.

*3. Financial Payment Flows and Trade Corridors:*
Interbank settlement and trade graphs span six orders of magnitude in
volume; scale invariance (Proposition 1a) and exponential resolution
(Theorem 4) uncover regional clearing blocs including small economies.

*4. Web Clickstream and Software Call Graphs:*
User navigation and control-flow traces follow ergodic PageRank walks
with teleportation $tau$; minimizing $L(M)$ directly compresses
clickstream or execution-trace logs into a two-level hierarchy.

*5. High-Dimensional kNN Graphs with Heterogeneous Scales:*
On high-dimensional tabular or embedding graphs with heterogeneous
densities, Infomap's exponential resolution (Theorem 4) and
zero-overfitting bias (Theorem 3iii) preserve small dense clusters
without fragmenting uniform noise clouds.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Four failure modes govern Infomap:
(i) *Field-of-view over-segmentation*: on high-diameter grids or
elongated non-convex manifolds, local mixing chops single ribbons
into diffusive tiles.
(ii) *Sparse detectability gap*: on ultra-sparse SBM graphs near the
Kesten–Stigum bound @abbe2018sbmreview, $H(cal(Q))$ overhead can
favor $M_1 = {V}$ over weak communities.
(iii) *Teleportation distortion*: on sink-heavy DAGs, $tau$ shifts
$pi_u$ and can obscure downstream modules.
(iv) *Greedy local minima*: non-convexity requires multi-trial
multilevel refinement, elevating latency on dense graphs.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Citation / web / financial flow graphs],
    [Excellent: native directed flow compression (Prop. 1b).],
    [Tune $tau$ carefully on sink-heavy DAGs.],
    [Multi-scale networks with tiny cliques],
    [Excellent: exponential resolution $2^(2 ell) / (e ell)$ (Thm. 4).],
    [Use hierarchical Infomap for nested levels.],
    [Noisy / low-structure null graphs],
    [High: provably refuses memoryless splits (Thm. 3iii).],
    [May under-detect ultra-sparse weak SBM blocks.],
    [Long ribbons / planar grids / meshes],
    [Poor: field-of-view limit over-splits long paths.],
    [Use Girvan–Newman or Spectral @shi2000ncut.],
    [Sub-millisecond real-time triage],
    [Moderate: `log2` and multi-trial overhead.],
    [Use Chinese Whispers or Leiden @traag2019leiden.],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for
    Infomap.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Infomap occupies a distinctive information-theoretic niche relative
to Normalized Cut @shi2000ncut, Walktrap @pons2006walktrap, Markov
Clustering (MCL) @vandongen2000mcl, and Louvain / Leiden
@blondel2008louvain @traag2019leiden. The following exact reformulation
links $L(M)$ explicitly to classical graph conductance
@chung1997spectral @vonluxburg2007tutorial.

#proposition(
  title: "Exact Conductance–Entropy Decomposition of the Map Equation",
  origin: "Original Derivation",
)[
  Let $G = (V, E, w)$ be a connected undirected weighted graph
  ($tau = 0$). For each module $M_i$ in partition
  $M = {M_1, dots, M_K}$, let
  $phi(M_i) = "cut"(M_i, bar(M)_i) / "vol"(M_i) in [0, 1]$ denote its
  graph conductance, and let $H_2 (p) = -p log_2 p - (1 - p) log_2 (1 - p)$
  be binary entropy. Then $q_i^"out" = pi_i phi(M_i)$ and the Map
  Equation decomposes exactly as
  $
    L(M) - H(pi) &= q_"out" H(cal(Q)) - H(pi_1, dots, pi_K) \
      &quad + sum_(i = 1)^K pi_i (1 + phi(M_i))
      H_2 ((phi(M_i)) / (1 + phi(M_i))).
  $
]

#proof[
  On a connected undirected graph ($tau = 0$),
  $pi_u = d_u / "vol"(V)$, $pi_i = "vol"(M_i) / "vol"(V)$, and
  $
    q_i^"out" = sum_(u in M_i, v in.not M_i) (w_(u v)) / ("vol"(V))
      = ("cut"(M_i, bar(M)_i)) / ("vol"(V)) = pi_i phi(M_i),
  $
  so $p_i^"mod" = pi_i (1 + phi(M_i))$. Let $X tilde.op pi$ be a
  stationary vertex and $I = M(X)$ its module index. By the chain
  rule of Shannon entropy,
  $
    H(pi) = H(I) + sum_(i = 1)^K pi_i H(X | I = i),
  $
  where $H(I) = H(pi_1, dots, pi_K)$. Inside module codebook
  $cal(P)^i$, the exit event occurs with probability
  $alpha_i = q_i^"out" / p_i^"mod" = phi(M_i) / (1 + phi(M_i))$ and
  a node emission occurs with probability $1 - alpha_i =
  1 / (1 + phi(M_i))$, conditioned on which $X$ follows
  $(pi_u / pi_i)_(u in M_i)$. By the grouping axiom of Shannon
  entropy,
  $
    H(cal(P)^i) = H_2 (alpha_i) + (1 - alpha_i) H(X | I = i).
  $
  Multiplying by $p_i^"mod" = pi_i (1 + phi(M_i))$ gives
  $p_i^"mod" H(cal(P)^i) = pi_i (1 + phi(M_i)) H_2 (alpha_i) +
  pi_i H(X | I = i)$. Summing over $i = 1, dots, K$, adding
  $q_"out" H(cal(Q))$, and subtracting $H(pi)$ cancels
  $sum_i pi_i H(X | I = i)$ identically, proving the formula.
]

Proposition 2 illuminates why Infomap generalizes Normalized Cut
@shi2000ncut ($"Ncut"(M) = sum_i phi(M_i)$): while Ncut sums
unweighted conductances and requires fixing $K$ externally, Infomap
weights each conductance nonlinearly through
$pi_i (1 + phi(M_i)) H_2 (phi(M_i) / (1 + phi(M_i)))$ and balances
the module-partition entropy reward $-H(pi_1, dots, pi_K)$ against
the inter-module switching cost $q_"out" H(cal(Q))$, automatically
determining $K$.

*Empirical Benchmark Behavior:* On the ClusterOverview 4-dataset
benchmark suite (Silhouette @rousseeuw1987silhouettes and Adjusted
Rand Index @hubert1985comparing), Infomap achieves the highest
*High-D* Silhouette (*0.318*, ARI = N/A) in the Algorithms
111–114 group and ties Chinese Whispers for the best *Customers*
score (*0.430*), corroborating Theorem 3 and Theorem 4: in noisy,
high-dimensional kNN graphs, Infomap's information-theoretic
null-model protection prevents fragmenting diffuse clouds into
spurious micro-clusters. On *Multi-Shape*, Infomap achieves Sil =
*0.478* (ARI = N/A), trailing Girvan–Newman's 0.782 because
field-of-view over-segmentation can subdivide elongated non-convex
manifolds into diffusive chunks. On *Iris*, it records Sil =
*0.339* and ARI = *0.400*, reflecting conservative merging across
the overlapping Versicolor–Virginica transition zone where
inter-class conductance $phi(M_i)$ is high. Its mean runtime of
*95.352 ms* is the slowest of the four chapters on these compact
benchmarks—reflecting flow-table initialization, logarithmic
evaluations, and multilevel trial restarts—even though its
asymptotic $cal(O)(m log n)$ scaling overtakes Girvan–Newman's
$cal(O)(m^2 n)$ and Walktrap's $cal(O)(m n H)$ on large sparse
graphs.
