#import "../template/ieee.typ": *

#alg_entry(
  id: 77,
  name: "GRINCH: Grafting and Rotation for Incremental Hierarchical Clustering",
  category: "Cat. III: Hierarchical & Tree-Based",
  year: "2019",
  paradigm: "Incremental Non-Local Tree Restructuring",
  complexity: $O(N d log N)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Online hierarchical clustering algorithms such as PERCH @kobren2017perch (Alg. 76) and BIRCH @zhang1996birch rely on local updates or fixed bounding boxes: when an early insertion splits a latent cluster across distant branches of the dendrogram, local sibling-aunt rotations cannot cross multi-hop structural barriers, and axis-aligned bounding boxes degrade in high-dimensional or non-Euclidean similarity spaces. *GRINCH* (*Grafting and Rotation for Incremental Hierarchical Clustering*) @monath2019grinch overcomes both limitations by operating over an arbitrary symmetric linkage similarity function $f: 2^X times 2^X -> bb(R)$ and combining *local rotations* with *non-local grafting* (re-attaching an entire subtree to a distant compatible branch) plus bounded-depth *restructuring*.

Let $X = \{bold(x)_1, dots.c, bold(x)_N\} subset.eq bb(R)^d$ arrive sequentially, and let $cal(T) = (V, E)$ be a rooted binary cluster tree with leaf set $"lvs"(v) subset.eq X$ under each node $v in V$. Denote the sibling of $v$ by $"sib"(v)$, its parent by $"par"(v)$, and its two children (for internal $u$) by $L(u), R(u)$. Each internal node $u$ carries a merge score $"score"(u) = f("lvs"(L(u)), "lvs"(R(u)))$. For additive linkages such as average inner-product (cosine on unit spheres) and squared Euclidean UPGMA @sokal1958statistical, each node stores compact sufficient statistics $(n_u, bold(s)_u, q_u)$ where
$ n_u = |"lvs"(u)|, quad bold(s)_u = sum_(bold(x) in "lvs"(u)) bold(x), \
  q_u = sum_(bold(x) in "lvs"(u)) || bold(x) ||_2^2. $

#definition(title: "Linkage Consistency, Rotation, and Grafting Primitives")[
  Let $f: 2^X times 2^X -> bb(R)$ be a symmetric linkage similarity and $cal(C)^* = \{C_1^*, dots.c, C_K^*\}$ a ground-truth partition.
  1. Linkage $f$ is *cluster-consistent* with $cal(C)^*$ @monath2019grinch if for every cluster $C_k^*$, every pair of nonempty disjoint subsets $A, B subset.eq C_k^*$, and every nonempty subset $S subset.eq X$ containing at least one foreign element ($S without C_k^* != emptyset$),
     $ f(A, B) > f(A, S). $
  2. *Local Rotation*: Given node $v$ with sibling $s = "sib"(v)$, parent $p = "par"(v)$, and aunt $a = "sib"(p)$, GRINCH rotates $v$ with $a$ (swapping $s$ and $a$ so $v, a$ share parent $p$) whenever
     $ f("lvs"(v), "lvs"(a)) > f("lvs"(v), "lvs"(s)). $
  3. *Non-Local Grafting*: Given node $v$ and a non-ancestor, non-descendant candidate node $u$ (with $u != "sib"(v)$), GRINCH detaches $v$ from $"par"(v)$ (replacing $"par"(v)$ with $"sib"(v)$) and inserts a new internal parent $p_"new" = (v, u)$ in place of $u$ whenever
     $ f("lvs"(v), "lvs"(u)) > max(&"score"("par"(v)), \
       &"score"("par"(u))). $
]

In addition, GRINCH optionally caps internal tree height via a bounded-size *restructure* step that runs exact batch HAC on a small sub-frontier around the newly formed node whenever a graft introduces a merge-score inversion or violates a maximum node-size threshold.

#algbox(
  title: "GRINCH (Incremental Rotation and Non-Local Grafting)",
  inputs: [Stream $X = (bold(x)_1, dots.c, bold(x)_N)$, linkage $f$, nearest-neighbor candidate budget $k_"nn"$, optional max-graft size $M_"max"$.],
  outputs: [Rooted binary dendrogram $cal(T)$ with sufficient statistics at all nodes.],
)[
  1. Initialize $cal(T)$ with root leaf $r$ holding $bold(x)_1$ and statistics $(1, bold(x)_1, || bold(x)_1 ||_2^2)$.
  2. *For* $t = 2, dots.c, N$:
     - *Nearest-Leaf Search*: Query exact or $k_"nn"$-approximate index for nearest leaf $ell^*$ maximizing $f(\{bold(x)_t\}, \{bold(y)_ell\})$.
     - *Insert Leaf*: Replace $ell^*$ with parent $p$ having children $ell^*$ and leaf $v_t = \{bold(x)_t\}$; update $(n, bold(s), q)$ along ancestors of $p$.
     - *Rotate Pass*: Set $v <- v_t$ (and check $ell^*$). *While* $"par"(v)$ and $"par"("par"(v))$ exist and $f("lvs"(v), "lvs"(a)) > f("lvs"(v), "lvs"(s))$ where $s = "sib"(v), a = "sib"("par"(v))$:
       - Swap $s$ and $a$; update $(n, bold(s), q)$ and $"score"$ at $p$ and $g = "par"(p)$.
     - *Graft Pass*: Starting from $v <- "par"(v_t)$ and walking upward (while $n_v <= M_"max"$):
       - Retrieve $k_"nn"$ nearest leaves to $v$ outside $"lvs"(v)$ and collect their ancestor candidates $cal(U)_v$ (excluding ancestors/descendants of $v$ and $"sib"(v)$).
       - Choose $u^* in arg max_(u in cal(U)_v) f("lvs"(v), "lvs"(u))$.
       - *If* $f("lvs"(v), "lvs"(u^*)) > max("score"("par"(v)), "score"("par"(u^*)))$: detach $v$ from $"par"(v)$, create $p_"new" = (v, u^*)$ at $u^*$, update statistics along paths from $"par"(v)$ and $u^*$ up to $w = "lca"(v, u^*)$, and continue upward from $v <- p_"new"$.
       - *Else* advance $v <- "par"(v)$.
  3. *Return* $cal(T)$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five formal results for GRINCH: exact $O(d)$ bilinear linkage evaluation and LCA-scoped updates, an exact Lyapunov equivalence between average-linkage rotations and Dasgupta's hierarchical cost, an explicit unit-sphere counterexample proving non-local barriers that trap PERCH while GRINCH escapes in one graft, exact dendrogram purity under cluster consistency, and capped-search complexity bounds.

#subsec("A", "Additive Sufficient Statistics and LCA-Scoped Updates")

#lemma(title: "Constant-Time Bilinear Linkage and Path-Scoped Invariance", origin: "Original Derivation")[
  For any two disjoint nodes $u, v in cal(T)$ storing $(n_u, bold(s)_u, q_u)$ and $(n_v, bold(s)_v, q_v)$, both the average inner-product similarity $f_"avg" (u, v)$ and average squared Euclidean distance $d_"avg"^2 (u, v)$ are computed exactly in $O(d)$ time via
  $ f_"avg" (u, v) &= (chevron.l bold(s)_u, bold(s)_v chevron.r) / (n_u n_v), \
    d_"avg"^2 (u, v) &= q_u / n_u + q_v / n_v - (2 chevron.l bold(s)_u, bold(s)_v chevron.r) / (n_u n_v). $
  Furthermore, when subtree $v$ is detached from $"par"(v)$ and grafted onto $u$ with $w = "lca"(v, u)$, the sufficient statistics of *every* ancestor strictly above $w$ (and of $w$ itself) remain invariant; updates are required solely along the two upward paths from $"par"(v)$ and $u$ to $w$.
]
#proof[
  By bilinearity of the Euclidean inner product,
  $ f_"avg" (u, v) &= 1 / (n_u n_v) sum_(bold(x) in "lvs"(u)) sum_(bold(y) in "lvs"(v)) chevron.l bold(x), bold(y) chevron.r \
                   &= (chevron.l bold(s)_u, bold(s)_v chevron.r) / (n_u n_v). $
  Expanding $|| bold(x) - bold(y) ||_2^2 = || bold(x) ||_2^2 + || bold(y) ||_2^2 - 2 chevron.l bold(x), bold(y) chevron.r$ and averaging over all $n_u n_v$ cross-pairs yields
  $ d_"avg"^2 (u, v) = (n_v q_u + n_u q_v - 2 chevron.l bold(s)_u, bold(s)_v chevron.r) / (n_u n_v), $
  which simplifies immediately to the stated $O(d)$ expression.

  For the graft update, let $P_v$ be the strict ancestors of $v$ strictly below $w$ and $P_u$ the strict ancestors of $u$ strictly below $w$. Detaching $v$ removes $"lvs"(v)$ from every node in $P_v$ (subtracting $(n_v, bold(s)_v, q_v)$), and attaching $v$ under $p_"new" = (v, u)$ adds $"lvs"(v)$ to every node in $P_u$ (adding $(n_v, bold(s)_v, q_v)$). Because $v$ and $u$ both reside strictly inside $w$, the net leaf set $"lvs"(w)$ is preserved identically, so $w$ and all strict ancestors of $w$ retain unchanged $(n, bold(s), q)$.
]

#subsec("B", "Exact Equivalence to Dasgupta Cost Descent Under Rotations")

A central theoretical question for incremental tree rotations is whether local swaps admit a global Lyapunov function guaranteeing cycle-free termination. We prove that under average linkage $f_"avg"$, GRINCH's local rotation rule is *exact coordinate descent* on Dasgupta's hierarchical clustering cost @dasgupta2016cost (equivalently strict ascent on Moseley--Wang revenue @moseley2017approximation).

#theorem(title: "Exact Dasgupta Cost Identity for Local Rotations", origin: "Original Proof")[
  Let $w(bold(x), bold(y)) >= 0$ be any symmetric pairwise similarity matrix on $X$, and let $f_"avg" (A, B) = (|A| |B|)^(-1) sum_(x in A, y in B) w(x, y)$. Define Dasgupta's cost @dasgupta2016cost and Moseley--Wang revenue @moseley2017approximation of a dendrogram $cal(T)$ by
  $ "Cost"(cal(T)) &= sum_(\{x, y\} subset.eq X) w(x, y) |"lvs"("lca"_cal(T) (x, y))|, \
    "Rev"(cal(T)) &= sum_(\{x, y\} subset.eq X) w(x, y) (N - |"lvs"("lca"_cal(T) (x, y))|). $
  If $cal(T)'$ is obtained from $cal(T)$ by rotating node $v$ with aunt $a$ (swapping sibling $s = "sib"(v)$ and aunt $a = "sib"("par"(v))$), then
  $ "Cost"(cal(T)') - "Cost"(cal(T)) &= -n_v n_s n_a (f_"avg" (v, a) \
                                     &quad - f_"avg" (v, s)), \
    "Rev"(cal(T)') - "Rev"(cal(T)) &= +n_v n_s n_a (f_"avg" (v, a) \
                                   &quad - f_"avg" (v, s)). $
  Consequently, every average-linkage GRINCH rotation strictly decreases $"Cost"(cal(T))$, guaranteeing finite cycle-free termination of any rotation sequence.
]
#proof[
  Let $p = "par"(v) = (v, s)$ and $g = "par"(p) = (p, a)$ in $cal(T)$, and let $p' = (v, a)$ and $g = (p', s)$ in $cal(T)'$. Partition all unordered pairs $\{x, y\} subset.eq X$ by their endpoint subtrees:
  - Pairs within $"lvs"(v)$, within $"lvs"(s)$, within $"lvs"(a)$, across $"lvs"(s) times "lvs"(a)$ (whose LCA is $g$ in both trees), and pairs with at least one point outside $"lvs"(g)$ all retain the exact same LCA leaf-set cardinality $|"lvs"("lca"(x, y))|$.
  - Every pair $\{x, y\} in "lvs"(v) times "lvs"(s)$ shifts its LCA from $p$ (size $n_v + n_s$) to $g$ (size $n_v + n_s + n_a$), increasing $|"lvs"("lca"(x, y))|$ by $+n_a$. Summing $w(x, y)$ over $"lvs"(v) times "lvs"(s)$ yields total weight $W(v, s) = n_v n_s f_"avg" (v, s)$, contributing $+n_a n_v n_s f_"avg" (v, s)$ to $"Cost"(cal(T)') - "Cost"(cal(T))$.
  - Every pair $\{x, y\} in "lvs"(v) times "lvs"(a)$ shifts its LCA from $g$ (size $n_v + n_s + n_a$) to $p'$ (size $n_v + n_a$), decreasing $|"lvs"("lca"(x, y))|$ by $-n_s$. Summing $w(x, y)$ over $"lvs"(v) times "lvs"(a)$ yields $W(v, a) = n_v n_a f_"avg" (v, a)$, contributing $-n_s n_v n_a f_"avg" (v, a)$.
  Adding both contributions gives
  $ "Cost"(cal(T)') - "Cost"(cal(T)) \
    = -n_v n_s n_a (f_"avg" (v, a) - f_"avg" (v, s)). $
  Since $"Rev"(cal(T)) = N sum_(\{x, y\}) w(x, y) - "Cost"(cal(T))$, the revenue increment has the opposite sign. Because GRINCH triggers a rotation iff $f_"avg" (v, a) > f_"avg" (v, s)$ and $n_v, n_s, n_a >= 1$, every rotation strictly reduces $"Cost"(cal(T))$. Since only finitely many ($(2 N - 3)!!$) rooted binary trees exist on $N$ leaves, no tree can repeat and any rotation pass terminates in finitely many steps.
]

#subsec("C", "Constructive Non-Local Barrier Separating Grafting From Rotations")

We now prove mathematically that local rotations alone—even if exhaustively tested at *every* node in the entire tree—cannot match the representational reach of GRINCH's grafting operator.

#theorem(title: "Strict Local-Rotation Trap Resolved by a Single Graft", origin: "Original Counterexample")[
  There exists a dataset $X subset.eq bb(S)^4 subset.eq bb(R)^5$ of $N = 10$ unit vectors partitioned into $K = 5$ clusters under cosine average linkage $f_"avg"$, together with an initial binary dendrogram $cal(T)_0$, such that:
  1. Every single node $v in cal(T)_0$ strictly fails the local rotation test ($f_"avg" (v, "aunt") < f_"avg" (v, "sib")$), making $cal(T)_0$ a strict local optimum of Dasgupta's cost under all single-edge rotations with $"DP"(cal(T)_0, cal(C)^*) = 31 / 35 < 1$.
  2. A single GRINCH graft step immediately produces a dendrogram $cal(T)_1$ with exact purity $"DP"(cal(T)_1, cal(C)^*) = 1$.
]
#proof[
  Define five unit prototype vectors $bold(u)_1, dots.c, bold(u)_5 in bb(S)^4 subset.eq bb(R)^5$ whose pairwise inner products form the $5 times 5$ Gram matrix
  $ S = mat(
    1, 0.6, 0.6, 0, 0;
    0.6, 1, 0, 0.5, 0;
    0.6, 0, 1, 0, 0.5;
    0, 0.5, 0, 1, 0;
    0, 0, 0.5, 0, 1
  ). $
  We verify $S$ is strictly positive definite so that $bold(u)_1, dots.c, bold(u)_5 in bb(S)^4$ exist via Cholesky factorization. Under the orthogonal symmetry split into antisymmetric vectors ($x_1 = 0, x_3 = -x_2, x_5 = -x_4$) and symmetric vectors ($x_3 = x_2, x_5 = x_4$), $S$ block-diagonalizes into the $2 times 2$ block $mat(1, 0.5; 0.5, 1)$ (eigenvalues $1.5, 0.5 > 0$) and the $3 times 3$ symmetric block
  $ M = mat(1, 0.6 sqrt(2), 0; 0.6 sqrt(2), 1, 0.5; 0, 0.5, 1), $
  whose leading principal minors are $M_1 = 1$, $M_2 = 1 - 0.72 = 0.28 > 0$, and $det(M) = 1 (1 - 0.25) - 0.72 = 0.03 > 0$. Hence $S succ 0$.

  Construct $X$ of size $N = 10$ with five clusters: $C_1^* = \{a_1, a_2\}$ (two copies of $bold(u)_1$), $C_2^* = \{b\}$ (copy of $bold(u)_2$), $C_3^* = \{c\}$ (copy of $bold(u)_3$), $C_4^* = D$ (three copies of $bold(u)_4$), and $C_5^* = E$ (three copies of $bold(u)_5$). Consider the dendrogram $cal(T)_0 = (((a_1, b), D), ((a_2, c), E))$, where $D$ and $E$ are internally pure binary subtrees. Check every node on the left branch for local rotations:
  - At $v = a_1$: sibling $b$ has $f_"avg" (a_1, b) = 0.6$, while aunt $D$ has $f_"avg" (a_1, D) = 0 < 0.6$.
  - At $v = b$: sibling $a_1$ has $f_"avg" (b, a_1) = 0.6$, while aunt $D$ has $f_"avg" (b, D) = 0.5 < 0.6$.
  - At $v = (a_1, b)$ (size $2$): sibling $D$ (size $3$) satisfies
    $ f_"avg" ((a_1, b), D) = (0 + 3 dot 0.5) / (2 dot 3) = 1.5 / 6 = 0.25, $
    whereas aunt $A = ((a_2, c), E)$ (size $5$) has total cross-weight $w(a_1, a_2) + w(a_1, c) + w(b, a_2) = 1.0 + 0.6 + 0.6 = 2.2$, yielding
    $ f_"avg" ((a_1, b), A) = 2.2 / (2 dot 5) = 0.22 < 0.25. $
  - At $v = D$: sibling $(a_1, b)$ has $f_"avg" = 0.25$, while aunt $A$ has $f_"avg" (D, A) = 0 < 0.25$.
  By exact left-right symmetry, every node on the right branch also strictly fails the rotation condition, and nodes inside $D$ or $E$ have constant similarity $1.0$. By Theorem II-B, $cal(T)_0$ is a strict local minimum of Dasgupta's cost under all local rotations. Its unordered cluster pairs number $|cal(P)^*| = 1 + 3 + 3 = 7$: all $6$ pairs in $D$ and $E$ have purity $1$, while $\{a_1, a_2\}$ has LCA at the root (size $10$) with purity $2 / 10$, giving $"DP"(cal(T)_0, cal(C)^*) = (6 + 0.2) / 7 = 31 / 35$.

  Now apply GRINCH's graft test to $v = a_2$ and $u = a_1$:
  $ f_"avg" (a_2, a_1) &= 1.0 \
                       &> max(f_"avg" (a_2, c), f_"avg" (a_1, b)) = 0.6. $
  GRINCH immediately detaches $a_2$ from $(a_2, c)$ and grafts it onto $a_1$, forming $cal(T)_1 = ((((a_1, a_2), b), D), (c, E))$. Every cluster $C_1^*, dots.c, C_5^*$ is now a pure maximal subtree, so $"DP"(cal(T)_1, cal(C)^*) = 1$.
]

#subsec("D", "Exact Dendrogram Purity Under Cluster Consistency")

#theorem(title: "Order-Independent Perfect Purity of Unconstrained GRINCH", origin: "Classical — Cited")[
  Suppose the linkage function $f$ is cluster-consistent with $cal(C)^*$ (Definition I.1) @monath2019grinch. Then for any arrival permutation $pi in S_N$, unconstrained GRINCH (testing all non-ancestor/descendant graft candidates without size cap $M_"max"$) constructs a dendrogram $cal(T)$ with $"DP"(cal(T), cal(C)^*) = 1$, even if the initial nearest-leaf placement is imperfect.
]
#proof[
  We induct on $t = 1, dots.c, N$ maintaining the invariant that every currently observed cluster $C_(k,t)^* = C_k^* inter \{bold(x)_1, dots.c, bold(x)_t\}$ appears as a pure maximal subtree $v_k$ in $cal(T)_t$. For $t = 1$ the claim is immediate.

  Assume the invariant holds at $t - 1$, and let $bold(x)_t in C_k^*$ be inserted next to any leaf $ell^* in C_j^*$ (whether $j = k$ or, due to approximate search, $j != k$).
  - If $j = k$, then $bold(x)_t$ already enters $v_k$, expanding it to $C_(k,t)^*$; consistency forbids any foreign subtree $u$ from exceeding the internal score of any edge inside $v_k$, preserving all pure maximal subtrees.
  - If $j != k$, splitting $ell^*$ places $v_t = \{bold(x)_t\}$ inside $v_j$. During the upward rotate pass from $ell^*$, every sibling-aunt step inside $v_j$ compares $f(A, B)$ for $A, B subset.eq C_j^*$ against $f(A, S)$ where $S$ contains $bold(x)_t in.not C_j^*$. By consistency, $f(A, B) > f(A, S)$, so local rotations eject $v_t$ completely outside $v_j$, restoring $v_j$ as a pure maximal subtree of $C_j^*$.
  - Next, if $C_(k,t-1)^* = emptyset$, then $v_t$ is already the pure maximal subtree of $C_(k,t)^*$. Otherwise, $C_(k,t-1)^*$ exists as a pure maximal subtree $v_k$ elsewhere in $cal(T)$. At node $v = v_t$, both $"par"(v_t)$ and $"par"(v_k)$ merge a pure subset of $C_k^*$ with a foreign subtree containing points outside $C_k^*$. By cluster consistency,
    $ f("lvs"(v_t), "lvs"(v_k)) > max(&"score"("par"(v_t)), \
                                      &"score"("par"(v_k))). $
    Thus the graft condition triggers for $u = v_k$ (or an even tighter subset of $v_k$), detaching $v_t$ and attaching it directly alongside $v_k$ to form the reunited pure maximal subtree $v'_k = (v_t, v_k)$ with $"lvs"(v'_k) = C_(k,t)^*$.
  At $t = N$, every $C_k^*$ forms a pure maximal subtree, hence $"DP"(cal(T), cal(C)^*) = 1$.
]

#subsec("E", "Capped-Search Complexity and Bounded-Depth Analysis")

#proposition(title: "Per-Insertion Complexity Under Capped Graft Search", origin: "Original Proof")[
  Suppose the tree height is $h$, the $k_"nn"$ nearest-leaf index answers queries in $T_"nn" (N, d)$ time, and the graft search restricts candidates to ancestors of the $k_"nn"$ nearest leaves with size cap $M_"max"$. If at most $G_"max"$ grafts are permitted per insertion, then each insertion requires at most
  $ O(T_"nn" (N, d) + (1 + G_"max") h k_"nn" d) $
  time under additive bilinear linkage, and $O(N d)$ total space across all $2 N - 1$ tree nodes.
]
#proof[
  The tree maintains $2 N - 1$ nodes, each storing $(n_v, bold(s)_v, q_v)$ in $d + 2$ scalars plus $O(1)$ pointers, giving $O(N d)$ total memory. Each insertion executes one nearest-leaf query ($T_"nn" (N, d)$), one leaf split, and at most $h$ ancestor rotations ($O(h d)$ time via Lemma II-A). At each of at most $h$ levels during the upward graft pass, collecting ancestors of $k_"nn"$ candidate leaves yields at most $h k_"nn"$ unique candidate nodes; evaluating $f_"avg"$ takes $O(d)$ per candidate via Lemma II-A, and executing a graft updates at most $2 h$ ancestors in $O(h d)$ time (Lemma II-A). Bounding restarts by $G_"max"$ yields $O(T_"nn" (N, d) + (1 + G_"max") h k_"nn" d)$ time per point; with $h = O(log N)$, $T_"nn" = O(d log N)$, and constant $k_"nn", G_"max"$, total time scales as $O(N d log N)$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

*Time, Space, and Parallelism.* Unlike PERCH's $2 d$-scalar bounding boxes, GRINCH stores only $(n_v, bold(s)_v, q_v)$—nearly half the per-node memory footprint ($O(N d)$)—and exploits exact $O(d)$ dot-product linkage updates without maintaining axis-aligned bounds. With balanced height $h = O(log N)$ and small $k_"nn"$ (typically $5$--$25$), GRINCH achieves $O(N d log N)$ practical scaling while repairing multi-hop mis-attachments. Its inherently sequential insertion loop limits intra-stream parallelism compared with round-based SCC @monath2021scc (Alg. 79), though $k_"nn"$ candidate scoring vectorizes cleanly via dense BLAS matrix-vector products.

*Geometric Invariances.* Under cosine average linkage $f_"avg" (u, v) = chevron.l bold(s)_u, bold(s)_v chevron.r / (n_u n_v)$ on unit-normalized vectors, GRINCH is *strictly invariant to orthogonal rotations* ($bold(x) mapsto bold(Q) bold(x)$) and positive radial rescaling—a decisive geometric advantage over PERCH's axis-aligned boxes. Under squared Euclidean average linkage $d_"avg"^2$, GRINCH is invariant to rigid motions (translations and orthogonal rotations) and equivariant to uniform scaling.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$O(N d log N)$ expected with $k_"nn"$ / $M_"max"$ caps; $O(N^2 d)$ unconstrained worst-case],
    [Space Complexity], [$O(N d)$ exact ($d + 2$ floats per node for additive sufficient statistics)],
    [Linkage Generality], [Arbitrary linkage $f(A, B)$: cosine UPGMA, squared Euclidean, Ward, learned neural scorers],
    [Invariances], [Orthogonal rotation invariance (unlike PERCH); rigid-motion invariance under $d_"avg"^2$],
    [Key Hyperparameters], [Candidate neighbor count $k_"nn"$, graft size cap $M_"max"$, restructure size cap],
    [Tree Operations], [Leaf split + local rotation + non-local grafting + optional bounded restructure],
  ),
  caption: [Operational and geometric profile of GRINCH (Alg. 77).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Streaming Cross-Document Entity Coreference.* Entity mentions embedded by transformer encoders lie on high-dimensional hyperspheres ($d = 256$--$1024$) where cosine average linkage matches semantic equivalence @monath2019grinch. GRINCH's orthogonal invariance and non-local grafting recover split entity mentions even under adversarial arrival order.
2. *Large-Scale Citation and Patent Disambiguation.* Multi-million-record scholarly streams arrive chronologically (inducing severe order bias as author topics drift over decades). Grafting merges late-career publications back into early-career branches separated by many hops in the dendrogram.
3. *Online Visual Instance and Product Deduplication.* E-commerce catalogs ingest noisy product photos continuously; cosine average linkage over normalized deep embeddings enables fast $O(d)$ graft checks that unify visually distinct angles of the same SKU once bridge images arrive.
4. *Incremental Speaker Diarization and Audio Clustering.* Real-time diarization streams short voice segments sequentially; GRINCH maintains an evolving speaker hierarchy with sub-millisecond incremental latency and corrects early short-utterance mis-assignments via grafting.
5. *Hierarchical Topic Discovery in Streaming Text.* News feeds and social streams exhibit both extreme cluster counts and multi-scale subtopic granularity; GRINCH supports exact $O(d)$ cosine UPGMA on TF-IDF or dense document vectors without storing an $N times N$ similarity matrix.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

GRINCH faces four operational failure regimes: (i) *Stale or low-recall $k_"nn"$ index*, where grafting fails to discover the true distant sibling $u^*$ if $u^*$ lies outside the top-$k_"nn"$ approximate neighbors; (ii) *Non-reducible or non-monotone linkages* (e.g., centroid linkage inversions), which can produce merge-score inversions along ancestor paths unless bounded restructuring is enabled; (iii) *Strongly non-convex manifold chains*, where average linkage prefers spherical clusters over elongated density filaments; and (iv) *Extreme distributed scale* ($N > 10^8$), where sequential point-by-point tree mutation bottlenecks distributed throughput relative to SCC (Alg. 79).

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [High-$d$ cosine entity coreference], [Rotation-invariant $O(d)$ UPGMA + non-local graft repair], [Tune $k_"nn"$ recall for very large $N$],
    [Order-biased temporal streams], [Grafting repairs multi-hop early insertion splits], [Enable bounded restructure if inversions appear],
    [Low-latency online clustering], [Sub-millisecond per-point update; compact $O(N d)$ state], [Use BICO (Alg. 78) under strict fixed-RAM cap],
    [Billion-scale batch distributed clustering], [Limited fit: inherently sequential insertion loop], [Use parallel round-based SCC (Alg. 79)],
    [Non-convex / manifold density chains], [Moderate fit: average linkage favors compact clusters], [Use HDBSCAN or single-linkage SLINK],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for GRINCH.],
)

#sec("V", "Structural Relations and Empirical Behavior")

GRINCH occupies the exact midpoint between PERCH @kobren2017perch (Alg. 76) and SCC @monath2021scc (Alg. 79). Setting $k_"nn" = 0$ (disabling grafting) and choosing AABB masking reduces GRINCH directly to PERCH; conversely, replacing sequential grafting with parallel batch reciprocal nearest-neighbor contractions yields SCC. Under average linkage, Theorem II-B also links GRINCH directly to Dasgupta's global hierarchical objective @dasgupta2016cost and Moseley--Wang's $1 / 3$-approximation analysis of batch UPGMA @moseley2017approximation.

#proposition(title: "Stationary Equivalence With Batch UPGMA on Ultrametric Data", origin: "Original Proof")[
  If the pairwise similarity matrix $w(bold(x), bold(y))$ satisfies the strong ultrametric inequality $w(bold(x), bold(z)) >= min(w(bold(x), bold(y)), w(bold(y), bold(z)))$ for all $bold(x), bold(y), bold(z) in X$ with distinct merge heights, then the dendrogram $cal(T)$ produced by unconstrained GRINCH under average linkage $f_"avg"$ is isomorphic to the unique exact batch UPGMA dendrogram @carlsson2010characterization @mullner2011modern.
]
#proof[
  Every ultrametric induces a unique laminar family of closed similarity balls $B(bold(x), tau) = \{bold(y) in X : w(bold(x), bold(y)) >= tau\}$ whose inclusion tree is the unique exact UPGMA dendrogram $cal(T)^*$ @carlsson2010characterization. At every threshold $tau$ separating an internal ball $B$ from its complement $X without B$, any $A, A' subset.eq B$ satisfy $f_"avg" (A, A') >= tau$, whereas any $S$ intersecting $X without B$ contains pairs with similarity strictly below $tau$, yielding $f_"avg" (A, S) < tau <= f_"avg" (A, A')$. Hence $f_"avg"$ is simultaneously cluster-consistent with *every* threshold partition of $cal(T)^*$. Applying Theorem II-D inductively to every level cut of the ultrametric hierarchy forces every ball $B$ of $cal(T)^*$ to appear as a pure maximal subtree in GRINCH's output $cal(T)$, proving $cal(T) equiv cal(T)^*$.
]

*Empirical Benchmark Analysis.* On the ClusterOverview 4-dataset benchmark @rousseeuw1987silhouettes @hubert1985comparing, GRINCH achieves: *Iris* Silhouette $= 0.489$, ARI $= 0.560$; *Customers* Silhouette $= 0.311$; *Multi-Shape* Silhouette $= 0.609$ (ARI N/A); *High-D* Silhouette $= 0.234$ (ARI N/A); with mean runtime $= 0.209$ ms. Two empirical signatures stand out and verify our analysis. First, GRINCH lifts Iris ARI dramatically from PERCH's $0.219$ to $0.560$ (matching batch-oriented SCC's exact $0.489$ / $0.560$ pair!) and raises both Customers ($0.273 -> 0.311$) and High-D ($0.220 -> 0.234$) Silhouette, confirming Theorem II-C's non-local repair power and orthogonal invariance over AABBs. Second, its $0.209$ ms mean runtime is the fastest of Algs. 76--79—$3.5 times$ faster than PERCH and $7.6 times$ faster than SCC—reflecting the low constant factor of Lemma II-A's $O(d)$ dot-product sufficient statistics.
