#import "../template/ieee.typ": *

#alg_entry(
  id: 60,
  name: "WPGMA (Weighted Pair Group Method with Arithmetic Mean)",
  category: "Cat. III: Hierarchical & Tree-Based",
  year: "1958",
  paradigm: "Agglomerative Ultrametric Linkage",
  complexity: $cal(O)(n^2)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Let $X = \{x_1, dots, x_n\}$ be a finite set of $n >= 2$ objects equipped with a symmetric, zero-diagonal dissimilarity matrix $D = (d_(i j)) in bb(R)_(>= 0)^(n times n)$ satisfying $d_(i i) = 0$ and $d_(i j) = d_(j i) >= 0$. Agglomerative hierarchical clustering constructs a sequence of nested partitions $cal(C)^((0)), cal(C)^((1)), dots, cal(C)^((n-1))$ starting from the singleton partition $cal(C)^((0)) = \{\{x_1\}, dots, \{x_n\}\}$ and ending at the trivial partition $cal(C)^((n-1)) = \{X\}$, together with merge heights $h_1, dots, h_(n-1) in bb(R)_(>= 0)$ encoding a rooted binary dendrogram $T$ on leaf set $X$.

Introduced by Sokal and Michener @sokal1958statistical alongside UPGMA (and later rediscovered as McQuitty's similarity analysis and unified in the Lance--Williams framework @lance1967general), the *Weighted Pair Group Method with Arithmetic Mean* (*WPGMA*) defines the dissimilarity between a newly merged cluster $A union B$ and any remaining cluster $C$ as the *unweighted arithmetic mean* of the two parent-branch dissimilarities $d_W (A, C)$ and $d_W (B, C)$, irrespective of the cardinalities $|A|$ and $|B|$. The term *weighted* reflects the resulting effect on individual leaves: because equal weight $1/2$ is assigned to each merging subtree rather than each leaf, elements belonging to smaller subtrees receive strictly higher effective weights in subsequent merges than elements buried deeper inside larger subtrees.

#definition(title: "WPGMA Linkage Recurrence and Cophenetic Ultrametric")[
  Let $cal(C)^((t-1))$ be the active cluster partition prior to merge $t in \{1, dots, n - 1\}$, initialized with $d_W (\{x_i\}, \{x_j\}) = d_(i j)$. At step $t$, WPGMA merges a minimizing pair
  $ (A_t, B_t) in arg min_(U != V in cal(C)^((t-1))) d_W (U, V), $
  records merge height $h_t = d_W (A_t, B_t)$, forms $cal(C)^((t)) = (cal(C)^((t-1)) without \{A_t, B_t\}) union \{A_t union B_t\}$, and updates inter-cluster dissimilarities to every $C in cal(C)^((t-1)) without \{A_t, B_t\}$ via the Lance--Williams rule with $alpha_A = alpha_B = 1/2$, $beta = 0$, $gamma = 0$:
  $ d_W (A_t union B_t, C) = 1 / 2 d_W (A_t, C) + 1 / 2 d_W (B_t, C). $
  In phylogenetic reconstruction under a molecular clock, the corresponding internal node height is scaled to $h_t / 2$ so that root-to-leaf path lengths equal $h_(n-1) / 2$. The induced *cophenetic dissimilarity* is $u(x_i, x_j) = h_(t(i, j))$ for $i != j$ (and $u(x_i, x_i) = 0$), where $t(i, j)$ is the unique step uniting $x_i$ and $x_j$.
]

Because WPGMA satisfies the Bruynooghe reducibility property (proved in Section II-B), the exact hierarchy can be computed in $cal(O)(n^2)$ worst-case time via the Nearest-Neighbor Chain (*NN-chain*) algorithm @mullner2011modern, avoiding the $cal(O)(n^3)$ cost of naive all-pairs scans.

#algbox(
  title: "WPGMA via Nearest-Neighbor Chain (Exact)",
  inputs: [Symmetric dissimilarity matrix $D = (d_(i j)) in bb(R)_(>= 0)^(n times n)$.],
  outputs: [Merge list $\{(A_t, B_t, h_t)\}_(t=1)^(n-1)$ defining dendrogram $T$.],
)[
  1. Initialize active set $cal(C) <- \{\{x_1\}, dots, \{x_n\}\}$, matrix $d_W (\{x_i\}, \{x_j\}) <- d_(i j)$, empty chain stack $S <- []$, and merge list $M <- []$.
  2. While $|cal(C)| > 1$:
     a. If $S$ is empty, push an arbitrary cluster $A in cal(C)$ onto $S$.
     b. Let $U$ be the top of $S$, and find its nearest active neighbor $V in arg min_(C in cal(C) without \{U\}) d_W (U, C)$ (breaking ties consistently, preferring the predecessor of $U$ in $S$ if tied).
     c. If $V$ is the predecessor of $U$ on $S$ (so $(U, V)$ is a reciprocal nearest-neighbor pair):
        - Pop $U$ and $V$ from $S$, record $(U, V, d_W (U, V))$ in $M$, and form $W <- U union V$.
        - For every $C in cal(C) without \{U, V\}$, set
          $ d_W (W, C) <- 1 / 2 d_W (U, C) + 1 / 2 d_W (V, C). $
        - Update $cal(C) <- (cal(C) without \{U, V\}) union \{W\}$.
     d. Else push $V$ onto $S$.
  3. Sort the $n - 1$ merge triples in $M$ by non-decreasing height $h_t$ and return $M$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

#subsec("A", "Dyadic Tree-Depth Representation")

Unlike UPGMA, whose inter-cluster dissimilarity is purely a set function of $(A, B)$, WPGMA depends intrinsically on the internal rooted binary merge subtrees $T_A$ and $T_B$ constructed prior to the comparison of $A$ and $B$. Our first lemma gives the exact closed-form expansion of $d_W (A, B)$ in terms of leaf depths.

#lemma(title: "Dyadic Tree-Depth Weight Expansion", origin: "Original Derivation")[
  Let $T_A$ and $T_B$ denote the rooted binary merge subtrees of two disjoint active clusters $A, B in cal(C)^((t))$. For each leaf $x in A$, let $ell(x, A) in bb(Z)_(>= 0)$ denote the path length (number of edges) from leaf $x$ to the root of $T_A$, with $ell(x, \{x\}) = 0$. Define leaf weights $mu_A (x) = 2^(-ell(x, A))$. Then:
  1. $mu_A$ is a probability distribution on $A$: $sum_(x in A) mu_A (x) = 1$.
  2. The WPGMA dissimilarity admits the exact bilinear expectation
     $ d_W (A, B) = sum_(x in A) sum_(y in B) 2^(-ell(x, A) - ell(y, B)) d(x, y). $
]

#proof[
  We prove both statements by induction on $|A| + |B| >= 2$. For the base case $|A| = |B| = 1$ with $A = \{x\}$ and $B = \{y\}$, we have $ell(x, A) = ell(y, B) = 0$, so $mu_A (x) = 2^0 = 1$ and the double sum reduces to $2^0 d(x, y) = d_W (\{x\}, \{y\})$.

  For the inductive step, assume without loss of generality that $|A| >= 2$, so $A$ was formed by merging two disjoint subtrees $A = A_1 union A_2$ with $|A_1| + |B| < |A| + |B|$ and $|A_2| + |B| < |A| + |B|$. Every leaf $x in A_1$ increases its depth by exactly one edge upon attaching $T_(A_1)$ beneath the root of $T_A$, so $ell(x, A) = ell(x, A_1) + 1$, and symmetrically $ell(x, A) = ell(x, A_2) + 1$ for $x in A_2$. Therefore,
  $ sum_(x in A) 2^(-ell(x, A)) &= 1 / 2 sum_(x in A_1) 2^(-ell(x, A_1)) + 1 / 2 sum_(x in A_2) 2^(-ell(x, A_2)) \ &= 1 / 2 (1) + 1 / 2 (1) = 1, $
  establishing Kraft's equality (1). Substituting the inductive hypothesis for $(A_1, B)$ and $(A_2, B)$ into the WPGMA recurrence gives
  $ d_W (A, B) &= 1 / 2 d_W (A_1, B) + 1 / 2 d_W (A_2, B) \ &= sum_(x in A_1) sum_(y in B) 2^(-(ell(x, A_1)+1) - ell(y, B)) d(x, y) \ &quad + sum_(x in A_2) sum_(y in B) 2^(-(ell(x, A_2)+1) - ell(y, B)) d(x, y) \ &= sum_(x in A) sum_(y in B) 2^(-ell(x, A) - ell(y, B)) d(x, y). $
  By symmetry in $A$ and $B$, the identity holds for all disjoint active clusters.
]

#subsec("B", "Reducibility, Monotonicity, and Exact Ultrametricity")

A central theoretical advantage of WPGMA over centroid methods (Algs. 61--62) is that its dendrogram is guaranteed to be free of height inversions and computable in quadratic time.

#theorem(title: "Bruynooghe Reducibility, Inversion-Free Ultrametricity, and NN-Chain Exactness", origin: "Classical — Cited")[
  WPGMA satisfies the Bruynooghe reducibility condition @lance1967general @mullner2011modern: for any three distinct clusters $A, B, C$ with $d_W (A, B) <= min(d_W (A, C), d_W (B, C))$,
  $ d_W (A union B, C) >= min(d_W (A, C), d_W (B, C)) >= d_W (A, B). $
  Consequently: (i) merge heights are monotonically non-decreasing, $h_1 <= h_2 <= dots <= h_(n-1)$; (ii) the cophenetic dissimilarity $u$ is an ultrametric on $X$, satisfying $u(x, z) <= max(u(x, y), u(y, z))$ for all $x, y, z in X$; and (iii) the NN-chain algorithm constructs an exact WPGMA dendrogram in $cal(O)(n^2)$ worst-case time.
]

#proof[
  Since $d_W (A union B, C) = (1/2) d_W (A, C) + (1/2) d_W (B, C)$ is a convex combination of $d_W (A, C)$ and $d_W (B, C)$, it satisfies
  $ min(d_W (A, C), d_W (B, C)) <= d_W (A union B, C) <= max(d_W (A, C), d_W (B, C)), $
  proving reducibility directly.

  For (i), at step $t$ the merged pair $(A_t, B_t)$ minimizes $d_W$ globally over $cal(C)^((t-1))$, so every untouched pair $(U, V)$ has $d_W (U, V) >= h_t$ and every newly formed pair $(A_t union B_t, C)$ has $d_W (A_t union B_t, C) >= min(d_W (A_t, C), d_W (B_t, C)) >= h_t$. Taking the minimum over $cal(C)^((t))$ yields $h_(t+1) >= h_t$.

  For (ii), fix any $x, y, z in X$ and let $m = max(t(x, y), t(y, z))$. Immediately after step $m$, $x$ and $y$ lie in the same cluster of $cal(C)^((m))$, as do $y$ and $z$. By transitivity of equivalence classes, $x$ and $z$ already lie in the same cluster, so $t(x, z) <= m$. Monotonicity (i) then gives $u(x, z) = h_(t(x, z)) <= h_m = max(u(x, y), u(y, z))$.

  For (iii), reducibility implies that merging any reciprocal nearest-neighbor pair $(U, V)$ cannot decrease the distance from any cluster $C$ on the chain stack $S$ below its existing distance to its stack successor, preserving the stack invariant that each element points to its nearest remaining neighbor @mullner2011modern. Each cluster is pushed onto $S$ at most once and popped at most once; each of the at most $3n - 3$ push/pop events performs one $cal(O)(n)$ scan or row update, giving $cal(O)(n^2)$ total time.
]

#subsec("C", "Exact Ultrametric Recovery and Perturbation Stability")

In phylogenetic systematics, where evolutionary distances arise from a molecular clock corrupted by finite-sequence sampling noise, the stability of tree recovery under $L_infinity$ perturbations is decisive.

#theorem(title: "Exact Ultrametric Recovery and Tight L-Infinity Stability", origin: "Original Proof")[
  Let $u^*$ be an ultrametric on $X$ with positive level-separation gap
  $ Delta = min \{ u^*(x, z) - u^*(x, y) : x, y, z in X, \, u^*(x, y) < u^*(x, z) \} > 0. $
  Suppose the observed dissimilarity matrix satisfies $d = u^* + epsilon$ with $norm(epsilon)_infinity = max_(i, j) |d_(i j) - u^*_(i j)| < Delta / 2$. Then WPGMA recovers a rooted binary tree $T$ compatible with every clade of $u^*$ (uniquely identical to the tree of $u^*$ if $u^*$ is binary) and achieves the tight cophenetic error bound
  $ norm(u - u^*)_infinity <= norm(epsilon)_infinity. $
  In particular, when $epsilon = 0$, WPGMA recovers $u^*$ identically: $u = u^*$.
]

#proof[
  Let $eta = norm(epsilon)_infinity < Delta / 2$. Say an active cluster $A subset.eq X$ is *$u^*$-consistent* if for every $x, x' in A$ and every $z in X without A$, $u^*(x, x') <= u^*(x, z)$, with the consequence that $u^*(x, z)$ depends only on $A$ and $z$ (by the ultrametric inequality, $u^*(x, z) = u^*(x', z)$ whenever $u^*(x, x') < u^*(x, z)$, and equality also holds when $u^*(x, x') = u^*(x, z) = h$ inside an unresolved multifurcation at level $h$). More precisely, call a partition $cal(C)^((t))$ *level-consistent* with $u^*$ if for every distinct $A, B in cal(C)^((t))$, $u^*(x, y)$ takes a single constant value $u^*(A, B)$ for all $(x, y) in A times B$, and $u^*(x, x') <= min_(C != A) u^*(A, C)$ for all $x, x' in A$.

  We prove by induction on $t = 0, dots, n - 1$ that $cal(C)^((t))$ is level-consistent with $u^*$ and that $|h_t - u^*(A_t, B_t)| <= eta$. The base case $t = 0$ is immediate since singletons satisfy both conditions vacuously.

  Assume $cal(C)^((t-1))$ is level-consistent with $u^*$. For any distinct $A, B in cal(C)^((t-1))$, Lemma 1 yields $sum_(x in A, y in B) mu_A (x) mu_B (y) = 1$ with non-negative weights, so
  $ |d_W (A, B) - u^*(A, B)| &= | sum_(x in A, y in B) mu_A (x) mu_B (y) (d(x, y) - u^*(A, B)) | \ &<= sum_(x in A, y in B) mu_A (x) mu_B (y) |epsilon(x, y)| <= eta. $
  Let $h^*_t = min_(U != V in cal(C)^((t-1))) u^*(U, V)$. Any pair $(A, B)$ with $u^*(A, B) = h^*_t$ satisfies $d_W (A, B) <= h^*_t + eta$. Conversely, any pair $(U, V)$ with $u^*(U, V) > h^*_t$ must satisfy $u^*(U, V) >= h^*_t + Delta$ by definition of $Delta$, and hence
  $ d_W (U, V) >= h^*_t + Delta - eta > h^*_t + eta, $
  using $eta < Delta / 2$. Consequently, WPGMA is forced to select a pair $(A_t, B_t)$ attaining $u^*(A_t, B_t) = h^*_t$.

  To verify level-consistency of $cal(C)^((t))$, consider any $C in cal(C)^((t-1)) without \{A_t, B_t\}$. For all $x in A_t$, $y in B_t$, and $z in C$, we have $u^*(x, y) = h^*_t <= min(u^*(A_t, C), u^*(B_t, C))$. By the ultrametric inequality on $(x, y, z)$, the two larger sides among $u^*(x, y), u^*(x, z), u^*(y, z)$ must coincide, so $u^*(A_t, C) = u^*(B_t, C)$; hence $u^*$ is constant across $(A_t union B_t) times C$. Moreover, every internal pair inside $A_t union B_t$ has $u^*$-distance at most $h^*_t <= u^*(A_t union B_t, C)$, preserving level-consistency. Finally, every pair $(x, y) in A_t times B_t$ satisfies $|u(x, y) - u^*(x, y)| = |h_t - u^*(A_t, B_t)| = |d_W (A_t, B_t) - u^*(A_t, B_t)| <= eta$. Since every unordered pair of distinct points is united at exactly one step $t$, $norm(u - u^*)_infinity <= eta$, and $eta = 0$ gives $u = u^*$.
]

#subsec("D", "Imbalanced Weight Concentration and Divergence from UPGMA")

Because leaf weights decay exponentially as $2^(-ell(x, A))$, WPGMA can deviate markedly from UPGMA whenever unbalanced ("caterpillar") subtrees form. We determine the exact minimal sample size at which their topologies diverge.

#proposition(title: "Caterpillar Weight Concentration and Minimal n = 5 Divergence from UPGMA", origin: "Original Counterexample")[
  1. For any dataset of cardinality $n <= 4$, WPGMA and UPGMA (Alg. 59) generate identical merge topologies (up to tie-breaking), differing at most in the final merge height $h_3$.
  2. The bound $n = 5$ is strictly minimal for topological divergence: on the 1D Euclidean dataset $X = \{0, 1, 8, 16, 28\} subset bb(R)$ with $d(x, y) = |x - y|$, WPGMA and UPGMA produce distinct dendrogram topologies and distinct 2-cluster partitions.
  3. On a sequential caterpillar cluster $A = (dots((\{x_1\} union \{x_2\}) union \{x_3\}) dots union \{x_m\})$ compared against a singleton $\{y\}$, the WPGMA--UPGMA discrepancy obeys
     $ d_W (A, \{y\}) - d_U (A, \{y\}) &= (1 / 2 - 1 / m) d(x_m, y) \ &quad + sum_(k=2)^(m-1) (2^(-(m-k+1)) - 1 / m) d(x_k, y) \ &quad + (2^(-(m-1)) - 1 / m) d(x_1, y), $
     placing mass $1/2$ on the single most recently joined element $x_m$ and vanishing mass $2^(-(m-1))$ on $x_1$.
]

#proof[
  For (1), when $|A| = |B| = 1$, $d_W (A union B, C) = (1/2) d(A, C) + (1/2) d(B, C) = d_U (A union B, C)$ identically. Thus WPGMA and UPGMA cannot disagree on any active distance until a cluster of size at least $3$ is formed. When $n <= 4$, after the first merge creates sizes $(2, 1, 1)$, either the two singletons merge into $(2, 2)$ or the size-$2$ cluster merges with a singleton into $(3, 1)$; in both cases only two active clusters remain, so the final merge is forced identically for both methods.

  For (2), trace both algorithms on $X = \{0, 1, 8, 16, 28\} subset bb(R)$:
  - *Step 1 (identical)*: The unique minimum distance is $d(0, 1) = 1$, merging $C_(1 2) = \{0, 1\}$ at $h_1 = 1$. Distances from $C_(1 2)$ are $d(C_(1 2), 8) = (8 + 7)/2 = 7.5$, $d(C_(1 2), 16) = 15.5$, and $d(C_(1 2), 28) = 27.5$.
  - *Step 2 (identical)*: Comparing $7.5$ against $d(8, 16) = 8$, $d(16, 28) = 12$, and $d(8, 28) = 20$, both methods merge $S = C_(1 2) union \{8\}$ at $h_2 = 7.5$.
  - *Step 3 (divergence)*: Remaining clusters are $S$, $\{16\}$, and $\{28\}$ with $d(16, 28) = 12$. Under WPGMA,
    $ d_W (S, \{16\}) = 1 / 2 (15.5) + 1 / 2 (8) = 11.75 < 12, $
    so WPGMA merges $S union \{16\}$ at $h_3 = 11.75$, yielding the 2-cluster cut $\{\{0, 1, 8, 16\}, \{28\}\}$. Under UPGMA,
    $ d_U (S, \{16\}) = 2 / 3 (15.5) + 1 / 3 (8) = 13 > 12, $
    so UPGMA merges $\{16, 28\}$ at $h_3 = 12$, yielding the distinct 2-cluster cut $\{\{0, 1, 8\}, \{16, 28\}\}$.

  For (3), in the caterpillar tree $T_A$, leaf $x_m$ has depth $1$, leaf $x_k$ ($2 <= k <= m - 1$) has depth $m - k + 1$, and leaf $x_1$ has depth $m - 1$. Subtracting the uniform UPGMA weight $1/m$ from the dyadic weight $2^(-ell(x_k, A))$ in Lemma 1 immediately gives the exact formula.
]

#subsec("E", "Metric Preservation and Transformation Equivariance")

We next establish which algebraic invariances WPGMA preserves across dissimilarity spaces, complementing the axiomatic characterization of hierarchical clustering in @carlsson2010characterization.

#proposition(title: "Positive Affine Equivariance, Metric Preservation, and Nonlinear Non-Invariance", origin: "Original Proof")[
  1. *Positive affine equivariance*: Under $tilde(d)(x, y) = a d(x, y) + b$ ($a > 0$, $b >= 0$, $x != y$), WPGMA yields the identical dendrogram topology with merge heights $tilde(h)_t = a h_t + b$.
  2. *Metric preservation*: If $d$ is a metric on $X$, then at every step $t$, $d_W$ (with $d_W (A, A) := 0$) remains a valid metric on $cal(C)^((t))$.
  3. *Non-invariance to monotone warping*: Unlike single and complete linkage, WPGMA is not invariant under strictly increasing nonlinear transforms $phi: bb(R)_(>= 0) -> bb(R)_(>= 0)$.
]

#proof[
  For (1), by Lemma 1, since $sum_(x in A, y in B) mu_A (x) mu_B (y) = 1$,
  $ tilde(d)_W (A, B) = sum_(x in A, y in B) mu_A (x) mu_B (y) (a d(x, y) + b) = a d_W (A, B) + b. $
  Because $r |-> a r + b$ is strictly increasing, the minimizing pair at every step is preserved identically and $tilde(h)_t = a h_t + b$.

  For (2), symmetry and positivity ($d_W (A, B) > 0$ for disjoint $A, B$) follow immediately from Lemma 1. For any three distinct clusters $A, B, C in cal(C)^((t))$, let $mu_A, mu_B, mu_C$ be their dyadic leaf probability distributions. Integrating the point-level triangle inequality $d(x, z) <= d(x, y) + d(y, z)$ against the product distribution $mu_A times.o mu_B times.o mu_C$ yields
  $ d_W (A, C) &= sum_(x, y, z) mu_A (x) mu_B (y) mu_C (z) d(x, z) \ &<= sum_(x, y, z) mu_A (x) mu_B (y) mu_C (z) (d(x, y) + d(y, z)) \ &= d_W (A, B) + d_W (B, C), $
  since each marginal sums to $1$.

  For (3), evaluate WPGMA on $X = \{0, 1, 8, 16, 28\}$ under $phi(r) = r^2$. Step 1 merges $C_(1 2) = \{0, 1\}$ at $1^2 = 1$. Step 2 compares $d_W^phi (C_(1 2), 8) = (64 + 49)/2 = 56.5$ against $8^2 = 64$, still merging $S = C_(1 2) union \{8\}$. At Step 3, $d_W^phi (C_(1 2), 16) = (256 + 225)/2 = 240.5$ and $d_W^phi (8, 16) = 64$, so
  $ d_W^phi (S, \{16\}) = 1 / 2 (240.5) + 1 / 2 (64) = 152.25 > 144 = phi(12) = d^phi (16, 28). $
  Hence WPGMA under $phi(r) = r^2$ merges $\{16, 28\}$ rather than $S union \{16\}$, strictly altering the topology relative to $phi(r) = r$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

WPGMA operates purely in dissimilarity space and requires no ambient vector coordinates. Because the update $d_W (A union B, C) = (d_W (A, C) + d_W (B, C)) / 2$ needs neither cluster sizes nor multiplications (reducing in binary arithmetic to an addition and a single exponent decrement), it has the lowest arithmetic constant factor among all average-based Lance--Williams linkages @lance1967general @mullner2011modern.

In the NN-chain implementation (Algorithm box, Section I), the upper-triangular dissimilarity matrix consumes $n(n - 1)/2$ floating-point words ($cal(O)(n^2)$ space), while the chain stack and active-cluster index occupy $cal(O)(n)$ space. Each merge overwrites row/column $min(U, V)$ in place and deactivates $max(U, V)$, scanning only contiguous active rows. Across all $n - 1$ merges, at most $3n - 3$ nearest-neighbor searches occur, yielding exact $cal(O)(n^2)$ worst-case time (versus $cal(O)(n^2 log n)$ for global min-heap priority queues and $cal(O)(n^3)$ for naive matrix scans). For large $n > 50,000$, the $cal(O)(n^2)$ memory footprint ($approx 10$ GB in double precision for $n = 50,000$) is the primary operational bottleneck, with cache performance governed by row-major streaming during the $cal(O)(n)$ neighbor scans.

Geometrically, WPGMA acts as a hybrid between UPGMA and local-branch linkage: by assigning weight $1/2$ to each coalescing lineage, it prevents massive clusters from completely drowning out small sister lineages during hierarchical aggregation, while still averaging pairwise distances to suppress single-link chaining bridges. Unlike centroid linkages, WPGMA never synthesizes off-manifold geometric representatives, allowing it to preserve moderately curved or non-spherical shapes when inter-manifold gaps exceed intra-manifold branch widths.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Input Space], [Arbitrary symmetric dissimilarity matrix $D in bb(R)_(>= 0)^(n times n)$],
    [Lance--Williams], [$alpha_A = alpha_B = 1/2$, $beta = 0$, $gamma = 0$ (combinatorial)],
    [Time Complexity], [$cal(O)(n^2)$ exact via NN-chain; $cal(O)(n^3)$ naive scan],
    [Space Complexity], [$cal(O)(n^2)$ in-place condensed dissimilarity array],
    [Reducibility], [Satisfied; guaranteed inversion-free ultrametric output],
    [Invariance], [Permutation and positive affine ($a d + b$); not nonlinear monotone],
    [Inductive Bias], [Equal lineage/branch weighting; resists large-clade masking],
    [Ultrametric Fit], [Exact recovery; $norm(u - u^*)_infinity <= norm(epsilon)_infinity$ for $norm(epsilon)_infinity < Delta/2$],
  ),
  caption: [Operational and mathematical profile of WPGMA (Alg. 60).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Phylogenetic Systematics with Unequal Taxon Sampling*: In classical numerical taxonomy @sokal1958statistical and microbial phylogenetics, different evolutionary clades are frequently sampled at wildly unequal densities (e.g., dozens of strains from a clinical outbreak versus a single isolate from a wild sister lineage). By assigning equal weight $1/2$ to each coalescing clade regardless of sample size, WPGMA prevents over-sampled clades from biasing divergence-time estimates toward their dense subclades, directly matching an equal-lineage evolutionary prior under an approximate molecular clock.
2. *Comparative Linguistics and Lexical Glottochronology*: Word-list cognate distance matrices across language families exhibit both strong ultrametric structure (divergent branching from proto-languages) and highly uneven documentation (many dialects recorded for dominant branches, single records for endangered sister languages). WPGMA's metric-preserving dyadic averaging recovers rooted family trees without allowing densely documented dialect continua to distort deep splits.
3. *Bioinformatics Microarray and Gene-Coexpression Clustering*: When grouping gene expression profiles from dissimilarity matrices (e.g., $1 - |r|$ correlation dissimilarities), housekeeping pathways contain hundreds of co-regulated transcripts whereas specialized regulatory modules contain only a few genes. WPGMA's branch-equalized linkage preserves small functional modules as distinct branches rather than allowing large housekeeping clusters to dominate inter-module distances @aggarwal2013data.
4. *Ecological Community and Biogeographical Classification*: In phytosociology and marine ecology, site-by-site Bray--Curtis or Jaccard dissimilarity matrices are non-Euclidean and often sample rare habitats sparsely alongside common biomes. Because WPGMA is purely combinatorial, preserves the triangle inequality (Proposition 2), and treats rare and abundant habitat branches symmetrically, it produces interpretable hierarchical biotope classifications without requiring Euclidean coordinates.
5. *Software Modularization and Dependency Tree Recovery*: In legacy software re-engineering, pairwise module coupling dissimilarities combine structural call graphs and semantic token overlap. Because utility libraries contain far more files than specialized domain controllers, equal branch weighting via WPGMA recovers balanced architectural hierarchies that avoid absorbing small domain subsystems into sprawling utility monoliths.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

WPGMA suffers from four well-characterized mathematical failure modes:
- *Non-ultrametric evolutionary rates (long-branch attraction)*: Whenever molecular evolution violates the constant-rate clock (additive trees that are not ultrametric), WPGMA systematically misgroups fast-evolving lineages with distant outgroups; Neighbor-Joining @saitou1987neighbor eliminates this bias by subtracting total-divergence terms.
- *Late-outlier sensitivity on caterpillar trees*: By Proposition 1(3), if a noise point or outlier attaches late to an imbalanced cluster of size $m$, it instantly receives weight $1/2$ in all subsequent merges, distorting higher-level branch distances far more strongly than in UPGMA (where its weight is $1/(m+1)$).
- *High-dimensional distance concentration*: In high-dimensional Euclidean spaces without low-dimensional intrinsic structure, pairwise distances concentrate around a common mean @beyer1999nearest, eroding the gap $Delta$ required by Theorem 2 and collapsing dendrogram contrast.
- *Quadratic memory barrier*: For $n > 10^5$ observations, storing the $cal(O)(n^2)$ dissimilarity matrix becomes infeasible, necessitating streaming pre-clustering such as BIRCH @zhang1996birch or approximate hierarchical trees @kobren2017perch.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Clock-like Phylogenies (Uneven Sampling)], [High: exact ultrametric recovery with equal clade weights], [Non-clock additive trees: use Neighbor-Joining @saitou1987neighbor],
    [Non-Euclidean Ecological Matrices], [High: metric-preserving combinatorial updates], [Noisy outliers: use UPGMA or HDBSCAN @campello2013hdbscan],
    [Gene Co-expression Hierarchies], [Moderate--High: protects small pathways from large blocks], [Late outlier distortion: use UPGMA or Ward @ward1963hierarchical],
    [Non-Convex Manifolds & Bridges], [Moderate: outperforms centroids, lags pure density methods], [Interlocking spirals: use DBSCAN @ester1996dbscan / HDBSCAN],
    [Large-Scale Data ($n > 10^5$)], [Low: $cal(O)(n^2)$ memory barrier for distance matrix], [Use BIRCH @zhang1996birch or PERCH @kobren2017perch],
  ),
  caption: [Domain suitability and failure-mode matrix for WPGMA.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Within the Lance--Williams hierarchical family @lance1967general, WPGMA occupies the exact midpoint between unweighted average linkage (UPGMA, Alg. 59) and weighted centroid linkage (WPGMC, Alg. 62). The following proposition formalizes its sandwiching between single and complete linkage and identifies the exact condition under which WPGMA coincides with UPGMA.

#proposition(title: "Linkage Sandwich Inequality and Exact Balanced-Tree Coincidence with UPGMA", origin: "Original Derivation")[
  Let $d_"SL" (A, B) = min_(x in A, y in B) d(x, y)$ and $d_"CL" (A, B) = max_(x in A, y in B) d(x, y)$ denote single and complete linkage dissimilarities, and let $d_U (A, B) = (|A| |B|)^(-1) sum_(x in A, y in B) d(x, y)$ denote UPGMA. For any disjoint active clusters $A, B$ with merge subtrees $T_A, T_B$:
  1. $d_"SL" (A, B) <= d_W (A, B) <= d_"CL" (A, B)$.
  2. If both $T_A$ and $T_B$ are perfect (all-equal-depth) binary trees of sizes $|A| = 2^a$ and $|B| = 2^b$, then $d_W (A, B) = d_U (A, B)$ identically.
  3. On squared Euclidean dissimilarities $D_(i j) = norm(bold(x)_i - bold(x)_j)_2^2$, WPGMA strictly dominates WPGMC (Alg. 62) by the sum of dyadic within-cluster variances: $D_W^((2)) (A, B) = D_M (A, B) + V_"dyad" (A) + V_"dyad" (B) >= D_M (A, B)$.
]

#proof[
  (1) By Lemma 1, $d_W (A, B)$ is a convex combination of $\{d(x, y) : (x, y) in A times B\}$, so it lies between their minimum $d_"SL" (A, B)$ and maximum $d_"CL" (A, B)$. (2) When every leaf $x in A$ has identical depth $ell(x, A) = a$, $2^(-ell(x, A)) = 2^(-a) = 1/|A|$, and similarly $2^(-ell(y, B)) = 1/|B|$, reducing Lemma 1 directly to $d_U (A, B)$. (3) Letting $bold(X) ~ mu_A$ and $bold(Y) ~ mu_B$ be independent with means $bold(m)_A, bold(m)_B$ (the WPGMC pseudo-centroids), expanding $norm(bold(X) - bold(Y))_2^2$ around $bold(m)_A - bold(m)_B$ eliminates cross terms and yields $D_W^((2)) (A, B) = norm(bold(m)_A - bold(m)_B)_2^2 + V_"dyad" (A) + V_"dyad" (B)$.
]

*Empirical Behavior on the ClusterOverview Benchmark*: Across the standardized 4-dataset ClusterOverview suite (evaluated with Silhouette @rousseeuw1987silhouettes and Adjusted Rand Index @hubert1985comparing), WPGMA achieves *Iris* $"Sil" = 0.447$, $"ARI" = 0.565$; *Customers* $"Sil" = 0.353$; *Multi-Shape* $"Sil" = 0.609$ ($"ARI" =$ N/A); and *High-D* $"Sil" = 0.236$ ($"ARI" =$ N/A), with a mean runtime of $2.209$ ms. Two structural signatures stand out in direct alignment with our theory:
- First, on *Multi-Shape*, WPGMA achieves $"Sil" = 0.609$ — substantially outperforming both Centroid Linkage (`0.452`, Alg. 61), Median Linkage (`0.452`, Alg. 62), and DIANA (`0.420`, Alg. 63). Because non-convex shapes place geometric centroids into empty off-manifold regions, centroid linkages conflate distinct manifolds; WPGMA avoids synthetic centroids altogether (Lemma 1) while its equal-branch weighting prevents large curved components from absorbing small separated components prematurely.
- Second, on *Iris*, WPGMA's $"ARI" = 0.565$ sits close to UPGMC (`0.560`) and DIANA (`0.574`) while trailing WPGMC (`0.653`), reflecting moderate sensitivity to the asymmetric boundary between *Versicolor* and *Virginica* where imbalanced caterpillar merges concentrate $1/2$ weight on late boundary joins (Proposition 1). On *High-D*, $"Sil" = 0.236$ reflects high-dimensional distance concentration @beyer1999nearest, while the $2.209$ ms mean runtime matches the $cal(O)(n^2)$ NN-chain cost of bottom-up agglomerative schemes.
