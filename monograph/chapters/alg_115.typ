#import "../template/ieee.typ": *

#alg_entry(
  id: 115,
  name: "Fast Greedy Modularity (CNM)",
  category: "Cat. V: Graph, Spectral & Network",
  year: "2004",
  paradigm: "Agglomerative Modularity Maximization",
  complexity: $cal(O)(m d log n)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Let $G = (V, E, w)$ be an undirected weighted graph with vertex set $V = {1, dots, n}$, edge set $E$ of cardinality $m = |E|$, and symmetric nonnegative adjacency matrix $A in bb(R)^(n times n)$ with entries $A_(i j) = w(i, j) > 0$ for $(i, j) in E$ and $A_(i j) = 0$ otherwise (for unweighted simple graphs, $A_(i j) in {0, 1}$ and $A_(i i) = 0$). Let $k_i = sum_(j=1)^n A_(i j)$ denote the weighted degree of vertex $i$, and let $2 M = sum_(i=1)^n k_i$ be twice the total edge weight (reducing to $2 m$ in the unweighted setting; henceforth we write $m$ for total weight $M$ so unweighted and weighted formulas coincide). A *community partition* $cal(C) = {C_1, dots, C_K}$ of $V$ assigns each vertex $i$ a label $c_i in {1, dots, K}$.

Newman and Girvan @newman2004finding quantified community strength via the *modularity* functional
$ Q(cal(C)) = 1 / (2 m) sum_(i=1)^n sum_(j=1)^n (A_(i j) - (k_i k_j) / (2 m)) delta(c_i, c_j), $
which compares the fraction of edge weight falling inside communities against its expectation under the degree-preserving configuration null model $P_(i j) = k_i k_j / (2 m)$. Defining the symmetric inter-community weight fraction $e_(r s)$ and community volume fraction $a_r$ by
$ e_(r s) = 1 / (2 m) sum_(i in C_r) sum_(j in C_s) A_(i j), \
  a_r = sum_(s=1)^K e_(r s) = ("vol"(C_r)) / (2 m), $
where $"vol"(C_r) = sum_(i in C_r) k_i$, modularity compresses into the community-level quadratic form
$ Q(cal(C)) = sum_(r=1)^K (e_(r r) - a_r^2). $

Directly maximizing $Q(cal(C))$ over all Bell-number partitions of $V$ is strongly NP-hard @brandes2008modularity. Earlier greedy agglomeration @newman2004finding maintained a dense $n times n$ matrix of merge gains at $cal(O)(n (m + n))$ cost. Clauset, Newman, and Moore (CNM) @clauset2004cnm observed that on sparse networks ($m = cal(O)(n)$), most community pairs share no edges and cannot yield positive gain; storing only active adjacent pairs in balanced search trees and max-heaps accelerates greedy maximization by orders of magnitude.

#definition(title: "CNM Sparse Modularity Gain Structure")[
  Given a current partition $cal(C)$ with active community index set $cal(I)$, the *merge modularity gain* for any distinct pair $r, s in cal(I)$ is
  $ Delta Q_(r s) = Q(cal(C)_(r union s)) - Q(cal(C)) = 2 (e_(r s) - a_r a_s), $
  where $cal(C)_(r union s)$ merges $C_r$ and $C_s$ into a single cluster. The CNM data structure maintains: (i) for each $r in cal(I)$, a sparse row $cal(R)_r = {(s, Delta Q_(r s)) : s != r, e_(r s) > 0}$ stored as a balanced binary search tree (and row max-heap) keyed by neighbor id $s$; (ii) a global max-heap $cal(H)$ storing each community's best pair $(r, s^*(r), Delta Q_(r, s^*(r)))$ where $s^*(r) = arg max_(s in cal(R)_r) Delta Q_(r s)$; and (iii) the degree-fraction array $a_r = "vol"(C_r) / (2 m)$.
]

#algbox(
  title: "Clauset-Newman-Moore (CNM) Fast Greedy Modularity",
  inputs: [Undirected graph $G = (V, E, w)$ with positive degrees $k_i > 0$ and total weight $m = 1 / 2 sum_i k_i$.],
  outputs: [Dendrogram of merges and optimal greedy partition $hat(cal(C))$ maximizing $Q$ along the greedy trajectory.],
)[
  1. *Singleton initialization:* Set $C_i <- {i}$ and $a_i <- k_i / (2 m)$ for all $i in V$, with initial modularity $Q^((0)) <- -sum_(i=1)^n a_i^2$.
  2. *Sparse gain table:* For each edge $(i, j) in E$ ($i != j$), insert $(j, Delta Q_(i j))$ into $cal(R)_i$ and $(i, Delta Q_(i j))$ into $cal(R)_j$, where
     $ Delta Q_(i j) = A_(i j) / m - (k_i k_j) / (2 m^2). $
     For each $i in V$, insert $(i, s^*(i), Delta Q_(i, s^*(i)))$ into global max-heap $cal(H)$.
  3. *Greedy agglomeration loop:* For step $t = 1, dots, n - 1$ while $cal(H)$ is non-empty:
     - Pop maximum-gain pair $(i, j, Delta Q_(i j))$ from $cal(H)$.
     - Record $Q^((t)) <- Q^((t-1)) + Delta Q_(i j)$. If early stopping is enabled and $Delta Q_(i j) <= 0$, terminate and return current partition.
     - Merge the smaller-degree row into the larger-degree row: assume $|cal(R)_i| <= |cal(R)_j|$ (otherwise swap $i, j$), absorbing $C_i$ into $C_j$.
     - For every neighbor $l in cal(R)_i union cal(R)_j without {i, j}$, update $Delta Q'_(j l)$ (and symmetric entry in $cal(R)_l$):
       $ Delta Q'_(j l) = cases(
         Delta Q_(i l) + Delta Q_(j l) & "if" l in cal(R)_i inter cal(R)_j,
         Delta Q_(i l) - 2 a_j a_l & "if" l in cal(R)_i without cal(R)_j,
         Delta Q_(j l) - 2 a_i a_l & "if" l in cal(R)_j without cal(R)_i.
       ) $
     - Delete row $cal(R)_i$, set $a_j <- a_i + a_j$, $a_i <- 0$, and refresh affected maxima in $cal(H)$.
  4. *Return:* Partition $hat(cal(C)) = cal(C)^((t^*))$ attaining $t^* = arg max_t Q^((t))$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

#subsec("A", "Exactness of Sparse Updates and Connectedness")

We first prove that the three-case CNM recurrence is algebraically exact and that restricting attention solely to adjacent community pairs cannot miss any modularity-improving merge.

#lemma(title: "Exactness of the CNM Recurrence and Adjacent-Merge Sufficiency", origin: "Original Proof")[
  Let $cal(C)$ be any partition of $G$ with $k_i > 0$ for all $i in V$. Then:
  (i) Merging distinct communities $C_i$ and $C_j$ into $C'_j = C_i union C_j$ updates the exact modularity gain to every other community $C_l$ ($l in.not {i, j}$) according to the three-case formula in Step 3 of the CNM algorithm;
  (ii) For any non-adjacent pair $e_(r s) = 0$, $Delta Q_(r s) = -2 a_r a_s < 0$, so every positive-gain merge joins two communities connected by at least one edge;
  (iii) Early-stopped CNM (halting once $max_(r, s) Delta Q_(r s) <= 0$) increases $Q$ strictly monotonically and produces solely connected induced subgraphs.
]

#proof[
  For (i), let $C'_j = C_i union C_j$. By additivity of edge weights and volumes over disjoint sets,
  $ e'_(j l) = e_(i l) + e_(j l), \
    a'_j = a_i + a_j. $
  Substituting into the definition $Delta Q'_(j l) = 2(e'_(j l) - a'_j a_l)$ yields
  $ Delta Q'_(j l) = 2(e_(i l) - a_i a_l) + 2(e_(j l) - a_j a_l). $
  If $l in cal(R)_i inter cal(R)_j$, both $e_(i l) > 0$ and $e_(j l) > 0$, so both summands are stored in $cal(R)_i$ and $cal(R)_j$, giving $Delta Q'_(j l) = Delta Q_(i l) + Delta Q_(j l)$. If $l in cal(R)_i without cal(R)_j$, then $e_(j l) = 0$, so the second term equals $-2 a_j a_l$ and $Delta Q'_(j l) = Delta Q_(i l) - 2 a_j a_l$. The third case $l in cal(R)_j without cal(R)_i$ follows symmetrically with $e_(i l) = 0$. Furthermore, for any community $l$ adjacent to neither $i$ nor $j$, $e'_(j l) = 0$, so $l$ stays outside $cal(R)'_j$. For any pair $(u, v)$ disjoint from ${i, j}$, neither $e_(u v)$ nor $a_u, a_v$ changes, leaving $Delta Q_(u v)$ invariant.

  For (ii), since every vertex has $k_v > 0$, every non-empty community satisfies $a_r = "vol"(C_r) / (2 m) > 0$. If $e_(r s) = 0$, then $Delta Q_(r s) = 2(0 - a_r a_s) = -2 a_r a_s < 0$. Hence any pair achieving $Delta Q_(r s) > 0$ must have $e_(r s) > 0$ and is present in $cal(R)_r$ and $cal(H)$.

  For (iii), each accepted step satisfies $Delta Q_(i j) > 0$, so $Q^((t)) > Q^((t-1))$. By induction on $t$: at $t = 0$, every singleton ${v}$ is connected; at step $t$, since $Delta Q_(i j) > 0$ implies $e_(i j) > 0$ by (ii), at least one edge bridges the connected subgraphs $G[C_i]$ and $G[C_j]$, so $G[C_i union C_j]$ is connected.
]

#subsec("B", "Constructive Suboptimality via Degree-Product Bias")

At iteration $t = 1$ on any unweighted simple graph, $Delta Q_(u v) = 1 / m - (k_u k_v) / (2 m^2)$ is strictly decreasing in the degree product $k_u k_v$. Thus CNM always prioritizes merging the lowest-degree-product edge first—even when that edge crosses the true community boundary. Because CNM is irreversibly agglomerative, an early cross-boundary merge permanently traps the algorithm in a strictly suboptimal partition with spurious clusters.

#theorem(title: "Infinite-Family Constructive Suboptimality of CNM", origin: "Original Counterexample")[
  For every integer $q >= 6$, let $r = floor(q / 2)$ and construct the unweighted simple graph $G_q$ on $n = 2 q + 2$ vertices consisting of two disjoint cliques $A_0 tilde.op K_q$ and $B_0 tilde.op K_q$, two peripheral vertices $x_A, x_B$ joined by a single bridge edge $(x_A, x_B)$, edges joining $x_A$ to $r$ distinct vertices of $A_0$, and edges joining $x_B$ to $r$ distinct vertices of $B_0$. Then $G_q$ has $m = q(q - 1) + 2 r + 1$ edges, and:
  (i) Early-stopped CNM deterministically terminates at the 3-community partition $hat(cal(C)) = {A_0, B_0, {x_A, x_B}}$;
  (ii) Relative to the natural 2-community partition $cal(C)^* = {A_0 union {x_A}, B_0 union {x_B}}$, CNM suffers an exact strictly positive modularity deficit
  $ Q(cal(C)^*) - Q(hat(cal(C))) = (r - 2) / m + (3 (r + 1)^2) / (2 m^2) > 0. $
]

#proof[
  Compute all vertex degrees in $G_q$: $k_(x_A) = k_(x_B) = r + 1$; the $r$ attached vertices in $A_0$ (and $B_0$) have degree $q$; and the remaining $q - r = ceil(q / 2) >= 3$ unattached vertices in each clique have degree $q - 1$. The edge count is $m = 2 binom(q, 2) + 2 r + 1 = q(q - 1) + 2 r + 1$.

  *Step 1 (forced cross-cut merge):* Since $q >= 6$ and $r = floor(q / 2)$, we have $r + 1 <= q / 2 + 1 < q - 1$. Compare degree products across all edge types:
  $ k_(x_A) k_(x_B) &= (r + 1)^2, \
    k_(x_A) k_a &= (r + 1) q > (r + 1)^2, \
    k_a k_(a') &>= (q - 1)^2 > (r + 1)^2. $
  Hence $(x_A, x_B)$ uniquely minimizes $k_u k_v$ over $E$, uniquely maximizes $Delta Q_(u v) = 1 / m - k_u k_v / (2 m^2) > 0$, and is merged at $t = 1$ into $S = {x_A, x_B}$ with $"vol"(S) = 2(r + 1)$.

  *Step 2 (unattached clique seeding):* After step 1, $S$ has a single edge to each attached vertex $a_i in A_0$, giving
  $ Delta Q(S, {a_i}) = 1 / m - (2(r + 1) q) / (2 m^2). $
  Meanwhile, since $q - r >= 3$, unattached pairs $(u, v)$ exist inside $A_0$ and $B_0$ with
  $ Delta Q({u}, {v}) = 1 / m - ((q - 1)^2) / (2 m^2). $
  Since $2(r + 1) >= q + 1$, we have $2(r + 1) q >= q^2 + q > (q - 1)^2$, so every unattached intra-clique pair strictly beats every edge incident to $S$ (and beats any intra-clique edge touching an attached vertex of degree $q$). Thus step 2 merges an unattached pair inside $A_0$ (or $B_0$ by symmetry) into $U_2$.

  *Steps $3, dots, q$ (complete clique coalescence):* Suppose $U_s subset A_0$ has size $2 <= s <= q - 1$. Any vertex $w in A_0 without U_s$ shares exactly $s$ edges with $U_s$, while $"vol"(U_s) <= s q$ and $k_w <= q$. Therefore,
  $ Delta Q(U_s, {w}) >= s / m - (s q^2) / (2 m^2) = s / m (1 - q^2 / (2 m)). $
  Observe that $2 m = 2 q^2 - 2 q + 4 r + 2 >= 2 q^2$ because $4 r + 2 >= 2(q - 1) + 2 = 2 q$. Hence $1 - q^2 / (2 m) >= 1 / 2$, giving
  $ Delta Q(U_s, {w}) >= s / (2 m) >= 1 / m > Delta Q(S, {a_i}) $
  and also strictly exceeding $Delta Q({b}, {b'}) < 1 / m$ for any singleton pair in $B_0$. By induction, $U_s$ absorbs all remaining vertices of $A_0$ consecutively until $C = A_0$ is formed.

  *Step $q + 1$ onward (isolation of $S$ and coalescence of $B_0$):* Now $A_0$ has $r$ edges to $S$ and $"vol"(A_0) = q(q - 1) + r = m - r - 1$. Its merge gain with $S$ is
  $ Delta Q(A_0, S) &= r / m - ((m - r - 1) dot 2(r + 1)) / (2 m^2) \
                    &= (r m - (m - r - 1)(r + 1)) / m^2 \
                    &= ((r + 1)^2 - m) / m^2. $
  Since $(r + 1)^2 <= (q / 2 + 1)^2$ while $m >= q^2$ for all $q >= 6$ (e.g., at $q = 6$, $(r+1)^2 = 16 < 37 = m$), we obtain $Delta Q(A_0, S) < 0$. Thus $A_0$ never absorbs $S$. By identical reasoning to Steps 2--$q$, $B_0$ now coalesces completely before any edge to $S$ is chosen, reaching ${A_0, B_0, S}$. At this state, $Delta Q(A_0, S) = Delta Q(B_0, S) < 0$ and $Delta Q(A_0, B_0) = -2 a_(A_0) a_(B_0) < 0$, so early-stopped CNM halts at $hat(cal(C)) = {A_0, B_0, S}$.

  *Modularity gap:* For $cal(C)^*$, only $(x_A, x_B)$ is cut and both volumes equal $m$, so $a_1 = a_2 = 1 / 2$ and $Q(cal(C)^*) = (m - 1) / m - 1 / 2 = 1 / 2 - 1 / m$. For $hat(cal(C))$, $2 r$ edges are cut and the volume fractions are $(m - r - 1) / (2 m)$, $(m - r - 1) / (2 m)$, and $(r + 1) / m$, yielding
  $ Q(hat(cal(C))) &= 1 - (2 r) / m - 2 ((m - r - 1) / (2 m))^2 - ((r + 1) / m)^2 \
                   &= 1 / 2 - (r - 1) / m - (3 (r + 1)^2) / (2 m^2). $
  Subtracting yields $Q(cal(C)^*) - Q(hat(cal(C))) = (r - 2) / m + (3 (r + 1)^2) / (2 m^2) > 0$ since $r = floor(q / 2) >= 3$.
]

#subsec("C", "Dendrogram-Depth Complexity and Super-Hub Degradation")

#theorem(title: "Dendrogram-Depth Time Complexity of CNM", origin: "Classical — Cited")[
  Let $d$ denote the height (depth) of the binary merge dendrogram produced by CNM @clauset2004cnm on a graph with $n$ vertices and $m$ edges. Using balanced binary search trees for the sparse rows $cal(R)_r$ and a binary max-heap for $cal(H)$, the total runtime of CNM is $cal(O)(m d log n)$. Consequently, on sparse graphs ($m = cal(O)(n)$) with balanced dendrograms ($d = cal(O)(log n)$), CNM runs in $cal(O)(n log^2 n)$ time, whereas on skewed dendrograms ($d = Theta(n)$) it degrades to $cal(O)(n^2 log n)$ time.
]

#proof[
  Consider a single merge of community $i$ into $j$ with $|cal(R)_i| <= |cal(R)_j|$. Iterating over $l in cal(R)_i$ performs $|cal(R)_i|$ balanced-tree lookups/insertions in $cal(R)_j$ and $cal(R)_l$, each taking $cal(O)(log n)$ time, plus at most $|cal(R)_i| + 1$ max-heap updates in $cal(H)$ at $cal(O)(log n)$ each. For $l in cal(R)_j without cal(R)_i$, Step 3 subtracts $2 a_i a_l$ from $Delta Q_(j l)$, which in the standard CNM implementation scans $cal(R)_j$ in $cal(O)(|cal(R)_j| log n)$ time.

  Thus one merge costs $cal(O)((|cal(R)_i| + |cal(R)_j|) log n)$. Each original edge $(u, v) in E$ can contribute to at most two row entries at any level of the merge tree, and each time either endpoint's community participates in a merge, that community moves up one level toward the root of the dendrogram. Because the dendrogram has height $d$, each of the $m$ edges is processed at most $2 d$ times across all merges (counting both endpoints). Summing over all at most $n - 1$ merges yields total work
  $ T_("CNM") &= sum_(t=1)^(n-1) cal(O)((|cal(R)_(i_t)| + |cal(R)_(j_t)|) log n) \
              &= cal(O)(m d log n). $
  When merges are balanced ($d = cal(O)(log n)$) and $m = cal(O)(n)$, this gives $cal(O)(n log^2 n)$. Conversely, when a dominant mega-community absorbs singletons sequentially—a known pathology on power-law hubs and caterpillars where $a_j$ grows linearly—$d = Theta(n)$ and repeated scans of $cal(R)_j$ incur $Theta(n^2 log n)$ operations.
]

#subsec("D", "Resolution Limit and Global Weight-Scaling Stability")

#theorem(title: "Modularity Resolution Limit and NP-Hardness Barrier", origin: "Classical — Cited")[
  Exact maximization of $Q(cal(C))$ is strongly NP-hard @brandes2008modularity. Furthermore, as established by Fortunato and Barthélemy @fortunato2007resolution, for a ring of $R >= 3$ identical cliques $K_p$ ($p >= 3$, each with $l = binom(p, 2)$ internal edges) connected by single bridge edges ($m = R(l + 1)$ total edges), merging any adjacent pair of cliques strictly increases $Q$ whenever $l + 1 < sqrt(m / 2)$ (equivalently $R > 2(l + 1)$), so neither CNM nor an exact global optimizer can resolve the $R$ cliques.
]

#proof[
  In the ring of $R$ cliques, each clique $C_r tilde.op K_p$ has internal edge count $l$ and two incident inter-clique bridges, so $"vol"(C_r) = 2 l + 2$ and $a_r = (l + 1) / m$. Two adjacent cliques $C_1, C_2$ share one bridge edge, giving $e_(1 2) = 1 / (2 m)$. The exact modularity change from merging $C_1$ and $C_2$ is
  $ Delta Q_(1 2) = 2 (e_(1 2) - a_1 a_2) = 1 / m - (2 (l + 1)^2) / m^2. $
  Hence $Delta Q_(1 2) > 0$ if and only if $2 (l + 1)^2 < m$, i.e., $l + 1 < sqrt(m / 2)$. Because $m = R(l + 1)$, this reduces to $R > 2(l + 1)$. Regardless of the algorithmic search strategy, the true $R$-clique partition is strictly dominated in $Q$ by the $R / 2$-pair partition (for even $R$). Strong NP-hardness of maximizing $Q$ follows via reduction from 3-Partition @brandes2008modularity.
]

#proposition(title: "Positive Scale Invariance and Edge-Perturbation Sensitivity", origin: "Original Derivation")[
  Let $tilde(A) = alpha A$ for any $alpha > 0$. Then $Q(cal(C); tilde(A)) = Q(cal(C); A)$ and $Delta Q_(r s) (tilde(A)) = Delta Q_(r s) (A)$ for every partition $cal(C)$ and pair $(r, s)$, so CNM is strictly invariant to positive weight scaling and vertex relabeling (modulo tie-breaking). Moreover, if a single edge weight $A_(u v)$ is increased by $epsilon > 0$ ($tilde(m) = m + epsilon$), then for any partition $cal(C)$,
  $ |Q(cal(C); tilde(A)) - Q(cal(C); A)| <= (3 epsilon) / (m + epsilon). $
]

#proof[
  Scaling $A -> alpha A$ scales $k_i -> alpha k_i$ and $m -> alpha m$, leaving $e_(r s) = (2 alpha m)^(-1) sum alpha A_(i j)$ and $a_r = (2 alpha m)^(-1) sum alpha k_i$ identically unchanged. Vertex permutation $P A P^T$ relabels indices without altering $e_(r s)$ or $a_r$.

  For the perturbation bound, write $Q = E_("in") - sum_r a_r^2$ and set $theta = epsilon / (m + epsilon) in (0, 1)$. Increasing $A_(u v)$ by $epsilon$ shifts the internal edge fraction $E_("in") = m_("in") / m in [0, 1]$ to $(m_("in") + chi epsilon) / (m + epsilon)$ with $chi = delta(c_u, c_v) in {0, 1}$, so
  $ |tilde(E)_("in") - E_("in")| = |(chi m - m_("in")) epsilon| / (m (m + epsilon)) <= theta. $
  For the volume fractions, $tilde(a)_r = ("vol"(C_r) + nu_r epsilon) / (2(m + epsilon)) = (1 - theta) a_r + theta b_r$, where $nu_r in {0, 1, 2}$ counts the endpoints of $(u, v)$ in $C_r$ and $b_r = nu_r / 2$ forms a probability vector ($sum_r b_r = 1$). Consequently, $||tilde(a) - a||_1 = theta ||b - a||_1 <= 2 theta$. Since $sum_r (tilde(a)_r - a_r) = 0$ and $tilde(a)_r + a_r in [0, 2]$ implies $||tilde(a) + a - bold(1)||_infinity <= 1$, Hölder's inequality yields
  $ |sum_r tilde(a)_r^2 - sum_r a_r^2| &= |chevron.l tilde(a) - a, tilde(a) + a - bold(1) chevron.r| \
    &<= ||tilde(a) - a||_1 dot ||tilde(a) + a - bold(1)||_infinity <= 2 theta. $
  Combining both bounds gives $|Q(cal(C); tilde(A)) - Q(cal(C); A)| <= 3 theta = (3 epsilon) / (m + epsilon)$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

*Time and Space Complexity.* Initializing the sparse trees $cal(R)_i$ and global max-heap $cal(H)$ takes $cal(O)(m log n)$ time and $cal(O)(m + n)$ memory—one node per undirected edge in each endpoint's tree plus $n$ entries in $cal(H)$ and $a$. Because merges only collapse existing edges and never create new non-zeros beyond the union of neighbor sets, peak memory remains $cal(O)(m + n)$ throughout execution, enabling million-edge graphs to fit in RAM where dense $cal(O)(n^2)$ modularity matrices fail @clauset2004cnm. Runtime scales as $cal(O)(m d log n)$; however, cache locality is poor because pointer-heavy balanced trees (red-black or AVL trees) suffer random memory accesses during row unions.

*Dendrogram Imbalance and Super-Hub Bottleneck.* Empirical analyses of scale-free networks reveal that CNM frequently creates one early "super-community" that greedily absorbs low-degree nodes sequentially, driving $d -> Theta(n)$ and causing $cal(O)(n^2 log n)$ slowdowns. Consolidation-ratio balancing or multilevel heuristics such as Louvain @blondel2008louvain and Leiden @traag2019leiden avoid this pathology by alternating local vertex moves with graph coarsening.

*Geometric Inductive Bias.* When applied to metric datasets via a symmetric $k$NN graph, CNM makes no convexity or isotropy assumptions in ambient Euclidean space: any non-convex manifold whose $k$NN edges remain dense internally and sparse across gaps achieves high modularity. Its primary geometric weaknesses are the global resolution limit @fortunato2007resolution and irreversible early bridging across sparse noise edges (Theorem II.B).

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$cal(O)(m d log n)$; $cal(O)(n log^2 n)$ balanced sparse, $cal(O)(n^2 log n)$ skewed],
    [Space Complexity], [$cal(O)(m + n)$ via sparse BST rows $cal(R)_r$ and global heap $cal(H)$],
    [Cluster Count $K$], [Selected automatically at global peak $arg max_t Q^((t))$],
    [Determinism], [Deterministic modulo tie-breaking in $cal(H)$],
    [Connectedness], [Guaranteed connected communities under $Delta Q > 0$ stopping],
    [Invariances], [Positive weight-scale invariant; permutation equivariant],
    [Key Limitation], [Irreversible greedy merges (Thm. II.B); resolution limit $sqrt(m / 2)$],
  ),
  caption: [Operational and complexity profile of Fast Greedy Modularity (CNM).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *E-Commerce Co-Purchasing and Recommender Graphs:* In Clauset, Newman, and Moore's original study of the Amazon co-purchasing network ($n > 4 times 10^5$ products) @clauset2004cnm, product graphs are extremely sparse ($m / n < 10$) and exhibit hierarchical topical clusters. CNM's $cal(O)(m + n)$ memory footprint and automatic $K$ selection recover niche genre communities without pre-specifying taxonomy cardinality.
2. *Citation and Co-Authorship Bibliometrics:* Scientific citation networks possess heavy-tailed degree distributions where the configuration null $k_i k_j / (2 m)$ explicitly discounts incidental citations to ubiquitous survey papers or foundational textbooks. CNM produces a full merge dendrogram capturing subfield-to-discipline hierarchies.
3. *Biological Protein-Protein Interaction (PPI) Modules:* Sparse physical interaction networks encode functional protein complexes as densely interconnected subgraphs. Because early-stopped CNM guarantees connected induced subgraphs (Lemma II.A), recovered complexes correspond directly to physically contiguous molecular pathways.
4. *Non-Convex Metric Data via $k$NN Graph Reduction:* By converting point clouds into symmetric $k$NN graphs, CNM partitions arbitrarily shaped, crescent, or nested manifolds solely from local neighborhood connectivity, bypassing the Voronoi hyperplanes of $k$-means @lloyd1982least.
5. *Power-Grid and Infrastructure Zone Partitioning:* Transmission grids are planar-like sparse networks with bounded average degree and clear regional bottlenecks. CNM's dendrogram identifies nested islanding zones for cascading-failure containment and decentralized state estimation.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

CNM fails systematically in four regimes: (i) *extreme scale heterogeneity* where small dense communities of weight $l < sqrt(m / 2)$ are merged by the modularity resolution limit @fortunato2007resolution; (ii) *low-degree inter-cluster bridges* where degree-product minimization triggers irreversible cross-boundary merges (Theorem II.B); (iii) *hub-dominated scale-free graphs* where sequential absorption into a single mega-community degrades runtime to $Theta(n^2 log n)$ and skews cluster sizes; and (iv) *overlapping or bipartite networks* where hard agglomerative partitions cannot represent multi-membership nodes or anti-modular disassortative blocks.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Sparse Co-Purchase / Web], [High: $cal(O)(m + n)$ RAM; dendrogram], [Super-hub slowdown; use Louvain @blondel2008louvain],
    [Multi-Shape $k$NN Clouds], [High: non-convex topological fit], [Sensitive to $k$NN bridge noise],
    [Multi-Scale Micro-Cliques], [Low: resolution limit $l < sqrt(m / 2)$], [Use Spinglass @reichardt2006spinglass or Leiden @traag2019leiden],
    [Overlapping Social Ties], [Low: strict hard hierarchy], [Use CPM @palla2005cpm or BIGCLAM @yang2013bigclam],
    [Dense / Complete Graphs], [Moderate: $m = Theta(n^2)$ slows BSTs], [Use Leading Eigenvector @newman2006leading],
  ),
  caption: [Domain suitability matrix for Fast Greedy Modularity (CNM).],
)

#sec("V", "Structural Relations and Empirical Behavior")

CNM sits at the root of the modularity-optimization lineage alongside Louvain @blondel2008louvain, Leiden @traag2019leiden, Leading Eigenvector (Alg. 116) @newman2006leading, and Spinglass (Alg. 117) @reichardt2006spinglass. Whereas Girvan-Newman @girvan2002community is top-down divisive via $cal(O)(m^2 n)$ edge-betweenness recomputation, CNM is bottom-up agglomerative. Its precise algebraic relationship to Louvain is formalized below.

#proposition(title: "CNM as a Monotone Agglomeration Subspace of Multilevel Louvain", origin: "Original Proof")[
  Every merge step of early-stopped CNM corresponds to a valid super-node aggregation step in Phase 2 of Louvain @blondel2008louvain preceded by a single-pair community merge, with identical modularity objective $Q$. Conversely, Louvain strictly expands CNM's feasible transitions by allowing individual vertex reallocations ($c_i : r -> s$) that escape the cross-boundary trap of Theorem II.B and recover $cal(C)^*$ strictly.
]

#proof[
  In Louvain, contracting communities ${C_r}$ into super-nodes with edge weights $W_(r s) = sum_(i in C_r, j in C_s) A_(i j)$ and self-loops $W_(r r) = sum_(i, j in C_r) A_(i j)$ preserves total weight $2 m$, super-node degrees $K_r = "vol"(C_r)$, and modularity $Q$ identically. Moving super-node $i$ into adjacent super-node $j$ changes $Q$ by $2(W_(i j) / (2 m) - (K_i K_j) / (2 m)^2) = Delta Q_(i j)$, matching CNM Step 3 exactly. Unlike CNM, Louvain's Phase 1 moves individual vertices $u in C_r$ to any neighbor community $C_s$ whenever $Delta Q_(u: r -> s) > 0$. On the counterexample graph $G_q$ ($q >= 6, r = floor(q / 2) >= 3$) of Theorem II.B, starting from CNM's trapped state $hat(cal(C)) = {A_0, B_0, S}$ with $S = {x_A, x_B}$, reallocating $x_A$ from $S$ into $A_0$ removes $1$ internal edge $(x_A, x_B)$ and adds $r$ internal edges to $A_0$, while updating community volumes from $(m - r - 1, 2(r + 1))$ to $(m, r + 1)$. The exact modularity gain is
  $ Delta Q_(x_A -> A_0) = (r - 1) / m - ((r + 1)(m - 2 r - 2)) / (2 m^2). $
  Multiplying by $2 m^2$ gives $2(r - 1) m - (r + 1)(m - 2 r - 2) = (r - 3) m + 2(r + 1)^2 > 0$ for all $q >= 6$. Subsequently reallocating the remaining singleton $x_B$ into $B_0$ gains $r / m - (r + 1)(m - r - 1) / (2 m^2) > 0$, strictly ascending to the optimal partition $cal(C)^*$.
]

*Empirical Behavior on ClusterOverview Benchmarks.* On the standardized 4-dataset ClusterOverview suite evaluated via Silhouette @rousseeuw1987silhouettes and Adjusted Rand Index (ARI) @hubert1985comparing, CNM attains *Iris* $"Sil" = 0.323$, $"ARI" = 0.373$; *Customers* $"Sil" = 0.386$; *Multi-Shape* $"Sil" = 0.782$; and *High-D* $"Sil" = 0.303$, with a mean runtime of $3.417$ ms. These measurements mirror the theoretical analysis directly: CNM dominates its Category V peers on *Multi-Shape* ($0.782$ vs. $0.267$ for Alg. 116 and $0.088$ for Alg. 117) because sparse $k$NN manifolds have high intra-manifold modularity and few cross-manifold bridges, enabling greedy connected merges (Lemma II.A) to trace non-convex geometry faithfully. Conversely, on *Iris*, dense bridge edges between the overlapping Versicolor and Virginica clouds induce early irreversible cross-cluster agglomeration (Theorem II.B) and resolution-limit coalescence (Theorem II.D), depressing ARI to $0.373$ relative to spectral bisection ($0.495$) and stochastic Potts annealing ($0.636$). On *Customers* ($0.386$, highest of the four) and *High-D* ($0.303$), CNM provides a fast, parameter-free baseline for sparse non-convex $k$NN topologies when exact resolution tuning is unnecessary.
