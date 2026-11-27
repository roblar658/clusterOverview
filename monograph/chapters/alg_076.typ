#import "../template/ieee.typ": *

#alg_entry(
  id: 76,
  name: "PERCH: Purity Enhancing Rotations for Cluster Hierarchies",
  category: "Cat. III: Hierarchical & Tree-Based",
  year: "2017",
  paradigm: "Incremental Extreme Hierarchical Clustering",
  complexity: $O(N d log N)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Classical agglomerative clustering algorithms such as SLINK @sibson1973slink, CLINK @defays1977clink, and Ward's method @ward1963hierarchical compute $Omega(N^2)$ pairwise affinities and scale poorly in the *extreme clustering* regime, where both the sample count $N$ and the number of latent clusters $K$ grow simultaneously ($K = Omega(N^alpha)$ for $alpha in (0, 1]$). Online conceptual clustering trees such as COBWEB @fisher1987cobweb and CLASSIT @gennari1989classit process observations sequentially, yet lack geometric search guarantees and accumulate irreversible insertion-order errors. *PERCH* (*Purity Enhancing Rotations for Cluster Hierarchies*) @kobren2017perch resolves both bottlenecks through an incremental binary cluster tree augmented with axis-aligned bounding boxes (AABBs), exact A\* nearest-leaf search, local tree rotations that repair early mis-attachments, and optional leaf-budget collapsing.

Let $X = \{bold(x)_1, dots.c, bold(x)_N\} subset.eq bb(R)^d$ arrive sequentially as a stream, partitioned by a ground-truth clustering $cal(C)^* = \{C_1^*, dots.c, C_K^*\}$. PERCH maintains a rooted binary tree $cal(T) = (V, E)$ whose leaves store points (or collapsed mini-clusters) with leaf set $"lvs"(v) subset.eq X$ under each node $v in V$. Each node $v$ maintains cardinality $n_v = |"lvs"(v)|$ and its tight axis-aligned bounding box
$ B_v = product_(j=1)^d [ell_(v,j), u_(v,j)], $
where $ell_(v,j) = min_(bold(x) in "lvs"(v)) x_j$ and $u_(v,j) = max_(bold(x) in "lvs"(v)) x_j$.

#definition(title: "Dendrogram Purity, AABB Bounds, and Masking")[
  Let $cal(P)_k^* = \{\{bold(x), bold(y)\} subset.eq C_k^* : bold(x) != bold(y)\}$ be the unordered intra-cluster pairs of $C_k^*$ and $cal(P)^* = union.big_(k=1)^K cal(P)_k^*$.
  1. The *dendrogram purity* of $cal(T)$ with respect to $cal(C)^*$ @kobren2017perch is
     $ "DP"(cal(T), cal(C)^*) &= 1 / (|cal(P)^*|) sum_(k=1)^K sum_(\{x, y\} in cal(P)_k^*) \
       &quad (|"lvs"("lca"_cal(T) (x, y)) inter C_k^*|) / (|"lvs"("lca"_cal(T) (x, y))|). $
  2. For a point $bold(x) in bb(R)^d$ and bounding box $B_v$, the exact minimum and maximum squared Euclidean distances from $bold(x)$ to $B_v$ are
     $ d_"min" (bold(x), B_v)^2 &= sum_(j=1)^d max(0, ell_(v,j) - x_j, x_j - u_(v,j))^2, \
       d_"max" (bold(x), B_v)^2 &= sum_(j=1)^d max(|x_j - ell_(v,j)|, |x_j - u_(v,j)|)^2. $
     Between two boxes $B_u$ and $B_v$, the minimum and maximum inter-box distances are defined analogously via coordinate intervals $[ell_(u,j), u_(u,j)]$ and $[ell_(v,j), u_(v,j)]$:
     $ d_"min" (B_u, B_v)^2 &= sum_(j=1)^d max(0, ell_(u,j) - u_(v,j), \
       &quad ell_(v,j) - u_(u,j))^2, \
       d_"max" (B_u, B_v)^2 &= sum_(j=1)^d max(|u_(u,j) - ell_(v,j)|, \
       &quad |u_(v,j) - ell_(u,j)|)^2. $
  3. Let node $v$ have sibling $s = "sib"(v)$, parent $p = "par"(v)$, grandparent $g = "par"(p)$, and aunt $a = "sib"(p)$. Node $s$ *masks* $a$ from $v$ whenever
     $ d_"max" (B_v, B_a) < d_"min" (B_v, B_s). $
     A *local rotation* resolves masking by swapping subtrees $s$ and $a$ so that $v$ and $a$ become siblings under $p$, while $s$ becomes the sibling of $p$ under $g$.
]

To prevent degenerate $O(N)$-depth chains on skewed streams, PERCH additionally applies a *balance rotation* when a node's sibling and aunt are equidistant within box bounds (or non-masking) and swapping reduces subtree depth imbalance without introducing masking. Under strict memory constraints with a maximum leaf budget $L$, *collapse mode* merges the internal node $u^*$ with minimum bounding-box diagonal $"diag"(B_u) = || bold(u)_u - bold(ell)_u ||_2$ into a single composite leaf whenever $|"lvs"(cal(T))| > L$.

#algbox(
  title: "PERCH (Incremental Hierarchical Insertion with Rotations)",
  inputs: [Stream $X = (bold(x)_1, dots.c, bold(x)_N) subset.eq bb(R)^d$, optional leaf budget $L in bb(N) union \{infinity\}$, balance flag.],
  outputs: [Binary cluster tree $cal(T)$ with AABB summaries $B_v$ at every node $v$.],
)[
  1. Initialize $cal(T)$ with root leaf $r$ holding $bold(x)_1$ and $B_r = product_j [x_(1,j), x_(1,j)]$.
  2. *For* $t = 2, dots.c, N$:
     - *Exact A\* Nearest-Leaf Search*: Maintain a min-heap $Q$ ordered by $d_"min" (bold(x)_t, B_u)$, initialized with $(d_"min" (bold(x)_t, B_"root"), "root")$. Pop node $u$ with smallest key; if $u$ is a leaf, set $ell^* <- u$ and stop; otherwise push both children $c_1, c_2$ of $u$ with keys $d_"min" (bold(x)_t, B_(c_i))$.
     - *Split Leaf*: Replace $ell^*$ with a new internal node $p$ whose children are $ell^*$ and a new leaf $v_t$ storing $bold(x)_t$. Set $B_p <- "hull"(B_(ell^*), B_(v_t))$.
     - *Upward Box Update and Rotation Pass*: Set $v <- ell^*$. *While* $"par"(v)$ and $"par"("par"(v))$ exist:
       - Let $p <- "par"(v)$, $s <- "sib"(v)$, $g <- "par"(p)$, $a <- "sib"(p)$.
       - Update $B_p <- "hull"(B_v, B_s)$ and $B_g <- "hull"(B_p, B_a)$.
       - *If* $d_"max" (B_v, B_a) < d_"min" (B_v, B_s)$ (or balance condition holds without violating non-masking): swap $s$ and $a$ (re-parenting $a$ under $p$ and $s$ under $g$), recompute $B_p <- "hull"(B_v, B_a)$, and continue from $v <- p$ (or re-check $v$ against its new aunt).
       - *Else* set $v <- p$.
     - *Collapse Mode*: If $|"lvs"(cal(T))| > L$, find internal node $u^* in arg min_u "diag"(B_u)$ and collapse $u^*$ into a composite leaf retaining $B_(u^*)$ and $n_(u^*)$.
  3. *Return* $cal(T)$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five formal properties governing PERCH: admissibility of its A\* bounding-box search, an exact algebraic formula for dendrogram purity gain under local rotations, perfect purity under box separability, a constructive high-dimensional counterexample revealing AABB inflation failure, and collapse-mode complexity and distortion bounds.

#subsec("A", "Admissibility and Consistency of AABB Best-First Search")

#lemma(title: "Admissibility and Monotonicity of Bounding-Box Search", origin: "Original Proof")[
  For any query $bold(x) in bb(R)^d$ and any node $u in cal(T)$, the AABB bounds satisfy
  $ d_"min" (bold(x), B_u) &<= min_(bold(y) in "lvs"(u)) || bold(x) - bold(y) ||_2 \
    &<= max_(bold(y) in "lvs"(u)) || bold(x) - bold(y) ||_2 <= d_"max" (bold(x), B_u). $
  Moreover, for any child $c$ of $u$, $d_"min" (bold(x), B_c) >= d_"min" (bold(x), B_u)$. Consequently, the A\* search in PERCH expands each node at most once and returns an exact Euclidean nearest leaf $ell^* in arg min_(ell in "lvs"(cal(T))) || bold(x) - bold(y)_ell ||_2$ upon popping the first leaf.
]
#proof[
  Fix $bold(y) in "lvs"(u)$. By definition of the tight AABB $B_u$, every coordinate satisfies $ell_(u,j) <= y_j <= u_(u,j)$ for $j in [d]$. Let $pi_j (x_j) = min(u_(u,j), max(ell_(u,j), x_j))$ be the metric projection of $x_j$ onto $[ell_(u,j), u_(u,j)]$. Then $|x_j - pi_j (x_j)| = max(0, ell_(u,j) - x_j, x_j - u_(u,j)) <= |x_j - y_j|$. Meanwhile, by convexity of $z mapsto |x_j - z|$ on $[ell_(u,j), u_(u,j)]$, its maximum is attained at an endpoint: $|x_j - y_j| <= max(|x_j - ell_(u,j)|, |x_j - u_(u,j)|)$. Squaring and summing across $j = 1, dots.c, d$ yields $d_"min" (bold(x), B_u)^2 <= || bold(x) - bold(y) ||_2^2 <= d_"max" (bold(x), B_u)^2$. Taking minimum and maximum over $bold(y) in "lvs"(u)$ proves the chain of inequalities.

  For monotonicity, $"lvs"(c) subset.eq "lvs"(u)$ implies $B_c subset.eq B_u$, so $d_"min" (bold(x), B_c) = inf_(bold(z) in B_c) || bold(x) - bold(z) ||_2 >= inf_(bold(z) in B_u) || bold(x) - bold(z) ||_2 = d_"min" (bold(x), B_u)$. Thus heap keys never decrease along root-to-leaf paths, preventing re-expansions. Suppose leaf $ell^*$ is the first leaf popped from $Q$. Because $B_(ell^*) = \{bold(y)_(ell^*)\}$, its key is exact: $d_"min" (bold(x), B_(ell^*)) = || bold(x) - bold(y)_(ell^*) ||_2$. For any other leaf $ell' in "lvs"(cal(T))$, either $ell'$ itself or a unique unexpanded ancestor $w$ of $ell'$ currently resides in $Q$. By heap optimality and admissibility,
  $ || bold(x) - bold(y)_(ell^*) ||_2 <= d_"min" (bold(x), B_w) <= || bold(x) - bold(y)_(ell') ||_2. $
  Hence $ell^*$ is an exact Euclidean nearest leaf.
]

#subsec("B", "Exact Algebraic Purity Change Under Local Rotation")

We derive the exact finite-sample change in dendrogram purity produced by a single PERCH rotation for arbitrary cluster compositions.

#theorem(title: "Exact Purity Increment of a Single Local Rotation", origin: "Original Derivation")[
  Let $v, s, a$ be pairwise disjoint subtrees in $cal(T)$ with $p = "par"(v) = "par"(s)$ and $g = "par"(p) = "par"(a)$, and let $cal(T)'$ be obtained by swapping $s$ and $a$ so that $p'$ has children $v, a$ and $g$ has children $p', s$. Denote $n_u = |"lvs"(u)|$ and $c_(u,k) = |"lvs"(u) inter C_k^*|$ for $u in \{v, s, a\}$ and $k in [K]$. Then
  $ "DP"(cal(T)', cal(C)^*) - "DP"(cal(T), cal(C)^*) = 1 / (|cal(P)^*|) sum_(k=1)^K Delta_k, $
  where each cluster's exact numerator contribution is
  $ Delta_k &= c_(v,k) c_(a,k) (c_(v,k) + c_(a,k)) / (n_v + n_a) \
            &quad - c_(v,k) c_(s,k) (c_(v,k) + c_(s,k)) / (n_v + n_s) \
            &quad + c_(v,k) (c_(s,k) - c_(a,k)) \
            &quad dot (c_(v,k) + c_(s,k) + c_(a,k)) / (n_v + n_s + n_a). $
  In particular, if $"lvs"(v), "lvs"(a) subset.eq C_k^*$ and $"lvs"(s) inter C_k^* = emptyset$, then
  $ "DP"(cal(T)', cal(C)^*) - "DP"(cal(T), cal(C)^*) = (n_v n_a n_s) / (|cal(P)^*| (n_v + n_s + n_a)) > 0. $
]
#proof[
  Consider every unordered pair $\{bold(x), bold(y)\} in cal(P)_k^*$ and determine whether its lowest common ancestor changes between $cal(T)$ and $cal(T)'$.
  - If both $bold(x), bold(y)$ lie strictly within one of $"lvs"(v)$, $"lvs"(s)$, or $"lvs"(a)$, their LCA lies inside that unchanged subtree.
  - If at least one point lies outside $"lvs"(g) = "lvs"(v) union "lvs"(s) union "lvs"(a)$, their LCA is a strict ancestor of $g$ whose leaf set is unchanged because $"lvs"(g)$ is preserved.
  - If one point lies in $"lvs"(s)$ and the other in $"lvs"(a)$, their LCA is $g$ in both $cal(T)$ and $cal(T)'$, with identical leaf set $"lvs"(g)$.
  Consequently, only two cross-sets change their LCA:
  1. The $c_(v,k) c_(s,k)$ pairs in $("lvs"(v) inter C_k^*) times ("lvs"(s) inter C_k^*)$ move their LCA from $p$ (purity $(c_(v,k) + c_(s,k)) / (n_v + n_s)$) to $g$ (purity $(c_(v,k) + c_(s,k) + c_(a,k)) / (n_v + n_s + n_a)$).
  2. The $c_(v,k) c_(a,k)$ pairs in $("lvs"(v) inter C_k^*) times ("lvs"(a) inter C_k^*)$ move their LCA from $g$ to $p'$ (purity $(c_(v,k) + c_(a,k)) / (n_v + n_a)$).
  Subtracting old purity contributions from new ones and collecting terms with denominator $n_v + n_s + n_a$ gives $Delta_k$ identically.

  For the pure-subtree specialization, $c_(v,k) = n_v$, $c_(a,k) = n_a$, $c_(s,k) = 0$, while $c_(v,j) = 0$ for all $j != k$. Hence $Delta_j = 0$ for $j != k$, and
  $ Delta_k &= n_v n_a (1 - (n_v + n_a) / (n_v + n_s + n_a)) \
            &= (n_v n_a n_s) / (n_v + n_s + n_a) > 0, $
  confirming strict purity improvement.
]

#subsec("C", "Order-Independent Perfect Purity Under Box Separability")

#theorem(title: "Perfect Dendrogram Purity Under AABB Separability", origin: "Classical — Cited")[
  Suppose the ground-truth partition $cal(C)^*$ satisfies *AABB cluster separability* @kobren2017perch: for every cluster $C_k^*$ and every pair of subsets $A, A' subset.eq C_k^*$ and $S subset.eq X without C_k^*$, the tight bounding boxes satisfy
  $ d_"max" (B_A, B_(A')) < d_"min" (B_A, B_S). $
  Then for any arrival permutation $pi in S_N$, non-collapsed PERCH ($L = infinity$) produces a tree $cal(T)$ with $"DP"(cal(T), cal(C)^*) = 1$.
]
#proof[
  We prove by induction on $t = 1, dots.c, N$ the *pure-maximal-subtree invariant*: after inserting $bold(x)_t$, every cluster $C_k^*$ having non-empty intersection $C_(k,t)^* = C_k^* inter \{bold(x)_1, dots.c, bold(x)_t\}$ appears as the exact leaf set $"lvs"(v_k) = C_(k,t)^*$ of some node $v_k in cal(T)_t$. For $t = 1$, the singleton root satisfies the invariant trivially.

  Assume the invariant holds at $t - 1$ and let $bold(x)_t in C_k^*$ arrive.
  *Case 1: $C_(k,t-1)^* != emptyset$.* For any $bold(y) in C_(k,t-1)^*$ and $bold(z) in.not C_k^*$, setting $A = \{bold(x)_t\}$, $A' = \{bold(y)\}$, and $S = \{bold(z)\}$ gives $|| bold(x)_t - bold(y) ||_2 = d_"max" (B_A, B_(A')) < d_"min" (B_A, B_S) = || bold(x)_t - bold(z) ||_2$. By Lemma II-A, exact A\* search returns a leaf $ell^* in C_(k,t-1)^* = "lvs"(v_k)$. Splitting $ell^*$ inserts $bold(x)_t$ inside $v_k$, expanding $"lvs"(v_k)$ to $C_(k,t)^*$. Every upward step within $v_k$ preserves $"lvs"(v_k)$ as a set, and once the upward pass reaches $v_k$, for its sibling $s$ and aunt $a$ (both disjoint from $C_k^*$), no foreign points enter $v_k$. Thus every $v_j$ remains pure and maximal.

  *Case 2: $C_(k,t-1)^* = emptyset$.* Let $ell^* in C_j^*$ ($j != k$) be the nearest leaf returned by A\* search, belonging to the pure maximal subtree $v_j$ of $C_j^*$. Splitting $ell^*$ creates parent $p$ with children $ell^*$ and $v_t = \{bold(x)_t\}$. If $ell^* = v_j$ (a singleton cluster), both $v_j$ and $v_t$ are immediately pure maximal subtrees. Otherwise, $ell^*$ is a strict descendant of $v_j$, so its aunt $a = "sib"(p)$ satisfies $"lvs"(a) subset.eq C_j^*$. By AABB separability with $A = "lvs"(ell^*) subset.eq C_j^*$, $A' = "lvs"(a) subset.eq C_j^*$, and $S = \{bold(x)_t\} subset.eq X without C_j^*$,
  $ d_"max" (B_(ell^*), B_a) < d_"min" (B_(ell^*), B_(v_t)). $
  Hence the masking condition fires, swapping $v_t$ and $a$: the reunited node $p' = (ell^*, a)$ is pure in $C_j^*$ while $v_t$ rises one level closer to $v_j$. Repeating this argument inductively at each ancestor strictly below $v_j$ bubbles $v_t$ completely out of $v_j$ until $v_t$ becomes the sibling of $v_j$. All clusters $C_(m,t)^*$ again occupy pure maximal subtrees. By induction, at $t = N$ every pair in $cal(P)_k^*$ has its LCA inside the pure subtree $v_k$, giving $"DP"(cal(T), cal(C)^*) = 1$.
]

#subsec("D", "High-Dimensional AABB Inflation Counterexample")

Theorem II-C relies on *AABB* separability rather than plain Euclidean separability. Because the diagonal of a $d$-dimensional bounding box can inflate by $Theta(sqrt(d))$ over pairwise point distances, we prove that Euclidean separability alone is strictly insufficient for PERCH's AABB masking rule whenever $d >= 3$.

#proposition(title: "Dimension-Induced AABB Masking Failure", origin: "Original Counterexample")[
  For every dimension $d >= 3$, there exists a finite dataset $X subset.eq bb(R)^d$ partitioned into two clusters $C_1^*, C_2^*$ satisfying strict Euclidean cluster separability,
  $ &max_(bold(x), bold(x)' in C_k^*) || bold(x) - bold(x)' ||_2 \
    &quad < min_(bold(x) in C_k^*, bold(y) in.not C_k^*) || bold(x) - bold(y) ||_2, $
  together with an arrival order under which PERCH's AABB masking test fails to trigger and yields
  $ "DP"(cal(T), cal(C)^*) = 1 - 2 / (d (d + 1)) < 1. $
]
#proof[
  Let $C_1^* = \{bold(e)_1, dots.c, bold(e)_d\} subset.eq bb(R)^d$ be the standard basis vectors and $C_2^* = \{bold(z)\}$ with $bold(z) = (1 + delta) bold(e)_1$, choosing $delta = 3 / 2$. For all $i != j$, intra-cluster distance is $|| bold(e)_i - bold(e)_j ||_2 = sqrt(2) approx 1.414$, while the minimum inter-cluster distance is
  $ min_(i in [d]) || bold(e)_i - bold(z) ||_2 = || bold(e)_1 - bold(z) ||_2 = 3 / 2 > sqrt(2), $
  since $|| bold(e)_j - bold(z) ||_2 = sqrt((5 / 2)^2 + 1) = sqrt(29) / 2$ for $j >= 2$. Thus strict Euclidean separability holds.

  Present $bold(e)_1, bold(e)_2, dots.c, bold(e)_d$ first, breaking equal-distance ties ($sqrt(2)$) toward the most recently inserted leaf. After $bold(e)_1$ and $bold(e)_2$ form the root, each $bold(e)_m$ ($m = 3, dots.c, d$) splits $bold(e)_(m-1)$ inside the right branch. Throughout these insertions, no masking inequality triggers because all pairwise distances equal $sqrt(2)$, yielding a root $g$ whose left child is $v = \{bold(e)_1\}$ and whose right child is subtree $a$ with $"lvs"(a) = \{bold(e)_2, dots.c, bold(e)_d\}$ and bounding box $B_a = \{0\} times [0, 1]^(d-1)$.

  Now $bold(z)$ arrives. Its unique nearest leaf is $bold(e)_1$ at distance $3 / 2$, so $bold(z)$ splits $bold(e)_1$ into $p = (v, s)$ with $s = \{bold(z)\}$ and aunt $a$. The farthest corner of $B_a$ from $bold(e)_1$ is $(0, 1, dots.c, 1)$, giving
  $ d_"max" (B_v, B_a) &= sqrt(1 + (d - 1)) = sqrt(d) \
                       &>= sqrt(3) > 3 / 2 = d_"min" (B_v, B_s) $
  for all $d >= 3$. Symmetrically, $d_"max" (B_s, B_a) = sqrt((5 / 2)^2 + d - 1) > 3 / 2 = d_"min" (B_s, B_v)$. Hence no masking rotation triggers, and $g$ has no grandparent. In the resulting tree, $C_1^*$ has $binom(d, 2)$ unordered pairs: the $binom(d - 1, 2)$ pairs inside $a$ have purity $1$, whereas the $d - 1$ pairs $\{bold(e)_1, bold(e)_j\}$ ($j >= 2$) have LCA $g$ of size $d + 1$ with purity $d / (d + 1)$. Since $C_2^*$ has zero pairs,
  $ "DP"(cal(T), cal(C)^*) &= (binom(d - 1, 2) + (d - 1) d / (d + 1)) / (binom(d, 2)) \
                           &= 1 - 2 / (d (d + 1)) < 1. $
]

#subsec("E", "Collapse-Mode Complexity and Approximation Distortion")

#proposition(title: "Leaf-Budget Complexity and Collapse Distortion Bound", origin: "Original Derivation")[
  In collapse mode with leaf budget $L >= 2$, PERCH maintains at most $2 L - 1$ nodes and uses $O(L d)$ memory. If tree height stays $h = O(log L)$ and collapse replaces internal node $u$ with a single representative $bold(c)_u in B_u$, then for any query $bold(x) in bb(R)^d$, the distance error introduced across all points hidden inside $u$ satisfies
  $ &max_(bold(y) in "lvs"(u)) | || bold(x) - bold(y) ||_2 - || bold(x) - bold(c)_u ||_2 | \
    &quad <= "diag"(B_u) = || bold(u)_u - bold(ell)_u ||_2. $
]
#proof[
  A full binary tree with at most $L$ leaves has at most $L - 1$ internal nodes, hence at most $2 L - 1$ nodes in total. Each node stores $2 d$ floats for $B_u$ plus $O(1)$ pointers and counts, giving strict $O(L d)$ space. Every insertion updates at most $h$ ancestors ($O(h d)$ time) and expands at most $2 L - 1$ nodes during A\* search (reduced to $O(d log L)$ when boxes are well-separated). For any $bold(y) in "lvs"(u) subset.eq B_u$ and $bold(c)_u in B_u$, coordinate containment $|y_j - c_(u,j)| <= u_(u,j) - ell_(u,j)$ implies $|| bold(y) - bold(c)_u ||_2 <= || bold(u)_u - bold(ell)_u ||_2 = "diag"(B_u)$. Applying the reverse triangle inequality $| || bold(x) - bold(y) ||_2 - || bold(x) - bold(c)_u ||_2 | <= || bold(y) - bold(c)_u ||_2$ completes the proof.
]

#sec("III", "Complexity, Geometry and Operational Profile")

*Time and Space Complexity.* In well-separated, low-intrinsic-dimension streams with balance rotations keeping height $h = O(log N)$, A\* search visits $O(log N)$ nodes per insertion and the upward update/rotation pass costs $O(d log N)$, yielding expected total runtime $O(N d log N)$ and $O(N d)$ space (or $O(L d)$ in collapse mode). Under adversarial overlapping high-dimensional data, A\* pruning weakens and search degenerates toward $O(N d)$ per insertion ($O(N^2 d)$ total).

*Geometric Invariances.* Because AABBs compute coordinate-wise minima and maxima, PERCH is invariant to translations ($bold(x) mapsto bold(x) + bold(t)$), uniform positive scaling ($bold(x) mapsto alpha bold(x)$), and coordinate permutations. However, AABBs are *not* rotation-invariant: an orthogonal rotation $bold(x) mapsto bold(Q) bold(x)$ alters bounding-box volume and diagonal lengths by up to $sqrt(d)$, directly affecting both A\* pruning efficiency and masking triggers.

*Cache and Memory Profile.* Tree traversal induces irregular pointer chasing across $2 N - 1$ nodes; contiguous arena allocation of node records $(bold(ell)_v, bold(u)_v, n_v)$ substantially improves L2/L3 cache locality and vectorizes $d_"min" / d_"max"$ loops via SIMD instructions.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$O(N d log N)$ expected balanced; $O(N^2 d)$ worst-case uncollapsed; $O(N L d)$ worst-case collapsed],
    [Space Complexity], [$O(N d)$ exact mode; $O(L d)$ bounded-memory collapse mode],
    [Streaming Support], [Native single-pass incremental insertion with online queryable hierarchy],
    [Invariances], [Translation, uniform scaling, coordinate permutation; *not* orthogonal rotation invariant],
    [Primary Hyperparameters], [Leaf budget $L$, balance-rotation toggle, collapse criterion],
    [Output Representation], [Rooted binary dendrogram supporting arbitrary $K$-cut extraction],
  ),
  caption: [Operational and geometric profile of PERCH (Alg. 76).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Extreme Author Disambiguation and Entity Resolution.* Bibliographic archives contain millions of mentions spanning $K > 10^5$ distinct authors @kobren2017perch. Batch HAC's $Omega(N^2)$ distance matrix is infeasible, whereas PERCH inserts new publications incrementally in $O(d log N)$ expected time and repairs early entity collisions via local rotations.
2. *Large-Scale Knowledge Base Coreference.* Streaming news and web extraction pipelines resolve entity mentions into fine-grained Wikipedia/Wikidata concepts where the cluster count grows nearly linearly with corpus size ($K = Omega(N^alpha)$) and dendrogram purity directly mirrors entity-linking F1.
3. *Online Face and Biometric Identity Indexing.* Metric-learned identity embeddings ($L_2$-normalized low-to-moderate dimensions after PCA) exhibit tight intra-person clusters; PERCH's A\* bounding-box search enables fast online enrollment without rebuilding global indexes.
4. *Fine-Grained Taxonomic and Metagenomic Binning.* Amplicon and k-mer abundance profiles span thousands of operational taxonomic units (OTUs) organized naturally as a multi-resolution phylogenetic hierarchy, matching PERCH's variable-granularity tree representation.
5. *Streaming Event Threading and Log Deduplication.* Continuous telemetry and breaking-news streams spawn new topics dynamically; PERCH's collapse mode bounds memory to $O(L d)$ while keeping fine-grained recent events uncollapsed until bounding boxes tighten.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

PERCH degrades predictably under four regimes: (i) *High ambient dimensionality* ($d >> 20$), where Proposition II-D proves AABB diagonals inflate by $sqrt(d)$, blunting both A\* pruning and masking triggers @beyer1999nearest; (ii) *Non-local tree barriers*, where a mis-attached subtree sits $>= 2$ hops away from its true cluster and local sibling-aunt rotations cannot jump across branches; (iii) *Strongly correlated oblique clusters*, where axis-aligned boxes overlap heavily even when clusters are linearly separable; and (iv) *Aggressive collapsing* ($L << K$), which irreversibly merges distinct minority clusters.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Extreme entity resolution ($K = Omega(N^alpha)$)], [Sub-quadratic streaming insertion; purity-preserving rotations], [Prefer GRINCH (Alg. 77) for non-local errors],
    [Low-$d$ spatially separated clusters], [Tight AABBs give fast A\* pruning and exact purity], [Pre-whiten / PCA if axes are strongly rotated],
    [Streaming bounded-memory indexing], [Collapse mode guarantees $O(L d)$ peak memory], [Use BICO (Alg. 78) if $k$-Means cost is primary],
    [High-$d$ dense embeddings ($d >= 50$)], [Weak fit: AABB diagonal inflates by $O(sqrt(d))$], [Use GRINCH (Alg. 77) or SCC (Alg. 79)],
    [Overlapping non-separable clusters], [Weak fit: local rotations freeze on early splits], [Use batch average-linkage HAC or SCC (Alg. 79)],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for PERCH.],
)

#sec("V", "Structural Relations and Empirical Behavior")

PERCH bridges classical online tree builders (BIRCH @zhang1996birch, ClusTree @kranen2011clustree, COBWEB @fisher1987cobweb) and modern scalable hierarchical clusterers (GRINCH @monath2019grinch, Alg. 77; SCC @monath2021scc, Alg. 79). Whereas BIRCH freezes micro-cluster assignments permanently without tree rotations, PERCH introduces local sibling-aunt rotations; GRINCH subsequently lifts PERCH's locality restriction via arbitrary-depth subtree grafting and arbitrary linkage functions.

#proposition(title: "Equivalence to Linkage Cuts Under AABB Separability", origin: "Original Proof")[
  Under AABB cluster separability (Theorem II-C), the $K$-cluster partition induced by the $K$ maximal pure subtrees of PERCH coincides identically with the $K$-cluster flat cuts of exact single-linkage (SLINK @sibson1973slink), complete-linkage (CLINK @defays1977clink), and average-linkage (UPGMA @mullner2011modern) hierarchical clustering.
]
#proof[
  AABB separability implies strict Euclidean separability: for all $k != j$, $max_(x, x' in C_k^*) || x - x' ||_2 < min_(x in C_k^*, y in C_j^*) || x - y ||_2$. In single-, complete-, and average-linkage HAC, every intra-cluster merge distance is bounded above by $max_k "diam"(C_k^*)$, whereas every cross-cluster distance is strictly larger. By induction on merges, all three batch algorithms complete every cluster $C_k^*$ before performing any cross-cluster merge, producing identical $K$-cluster cuts $cal(C)^*$. By Theorem II-C, PERCH also represents each $C_k^*$ as a pure maximal subtree, so extracting the $K$ pure maximal nodes recovers $cal(C)^*$ identically.
]

*Empirical Benchmark Analysis.* On the ClusterOverview 4-dataset benchmark @rousseeuw1987silhouettes @hubert1985comparing, PERCH records: *Iris* Silhouette $= 0.383$, ARI $= 0.219$; *Customers* Silhouette $= 0.273$; *Multi-Shape* Silhouette $= 0.709$ (ARI N/A); *High-D* Silhouette $= 0.220$ (ARI N/A); with mean runtime $= 0.722$ ms. These measurements align tightly with our theory: PERCH achieves the strongest Multi-Shape Silhouette ($0.709$) among Algs. 76--79 because spatially separated components yield tight, non-overlapping AABBs that enable accurate A\* attachment and purity-restoring rotations. Conversely, on Iris (where Versicolor and Virginica overlap, violating separability), local rotations cannot undo early cross-species merges, depressing ARI to $0.219$ compared with GRINCH's $0.560$ and BICO's $0.638$. Finally, High-D records the lowest Silhouette ($0.220$), directly corroborating Proposition II-D's $sqrt(d)$ AABB diagonal inflation mechanism.
