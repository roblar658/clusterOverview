#import "../template/ieee.typ": *

#cmp_entry(
  roman: "V",
  title: "Intra-Category Comparative Analysis: Hierarchical and Tree-Based (Algs. 56-80)",
  scope_label: "Head-to-Head Theoretical & Practical Comparison of Algorithms 56-80",
  provenance: "Original Unifying Proofs + Cited Classical Theorems",
)

The twenty-five hierarchical and tree-based algorithms of Category III (Algs. 56--80) replace flat $k$-partitioning with multi-scale nested hierarchies: they map a finite dissimilarity space $(cal(X), d)$ or Euclidean point cloud $cal(X) subset.eq bb(R)^d$ into a rooted cluster tree $cal(T)$ whose antichain cuts reveal structure across all resolutions @jain2010data @xu2005survey @mullner2011modern. Formally, a rooted hierarchy with monotone merge or split heights $h : V(cal(T)) -> bb(R)_(>= 0)$ is equivalent via Johnson's bijection to an *ultrametric* $u : cal(X) times cal(X) -> bb(R)_(>= 0)$ obeying $u(x, z) <= max{u(x, y), u(y, z)}$, whereas unrooted phylogenetic trees generalize this geometry to *additive tree metrics* governed by Buneman's four-point condition @carlsson2010characterization @saitou1987neighbor.

Algorithms 56--80 divide along four structural axes: (1) *Lance--Williams recurrence algebra*, where parameter inequalities separate inversion-free ultrametric linkages and $O(n^2)$ Bruynooghe-reducible Nearest-Neighbor Chain execution (Algs. 56--60, 72, 74 for $beta <= 0$) from non-reducible or inversion-prone methods (Algs. 61--62, 74 for $beta > 0$); (2) *chaining, diameter, and global tree cost*, contrasting ultrametric extrema (SLINK \#56, CLINK \#57) against prototype-anchored Minimax Linkage (\#73), Gini-gated Genie (\#69), and Dasgupta/Moseley--Wang optimal UPGMA (\#58) @dasgupta2016cost @moseley2017approximation; (3) *streaming and hybrid compression*, spanning Euclidean CF-trees (BIRCH \#64, coreset-calibrated BICO \#78), metric CR-trees (BUBBLE \#71), shrunk representatives (CURE \#65), dynamic $k$-NN graphs (CHAMELEON \#66), categorical link kernels (ROCK \#67), and density-core pre-compression (HCA-DBSCAN \#70); and (4) *divisive, additive, and web-scale hierarchies*, uniting DIANA (\#63), binary MONA (\#68), Neighbor-Joining (\#75), online restructuring (PERCH \#76, GRINCH \#77), and parallel reciprocal-neighbor contractions (SCC \#79, FINCH \#80). Sections I--IV establish unifying theorems across these four axes, culminating in Section V's complete 25-algorithm reference matrix and decision guide.

// =============================================================================
// SECTION I
// =============================================================================
#sec("I", "The Lance--Williams Unification, Monotonicity, and Space Distortion")

Nine algorithms in Category III—Single Linkage (SLINK, Alg. 56 @sibson1973slink), Complete Linkage (CLINK, Alg. 57 @defays1977clink), Average Linkage (UPGMA, Alg. 58 @sokal1958statistical), Ward's Method (Alg. 59 @ward1963hierarchical), WPGMA (Alg. 60 @sokal1958statistical), Centroid Linkage (UPGMC, Alg. 61 @sokal1958statistical), Median Linkage (WPGMC, Alg. 62 @lance1967general), and Flexible Beta Linkage (Alg. 74 @lance1967general)—are unified algebraically by the *Lance--Williams recurrence* (Alg. 72, @lance1967general). Let $cal(C)^((t))$ be the active cluster partition at step $t in {0, dots, n - 2}$ with symmetric non-negative inter-cluster dissimilarities $D(A, B) >= 0$. Upon merging a minimizing pair $(A, B) in arg min_(U != V) D(U, V)$ into $A union B$, the updated dissimilarity to every remaining cluster $C in cal(C)^((t)) without {A, B}$ obeys
$
  D(A union B, C) &= alpha_A D(A, C) + alpha_B D(B, C) \
    &quad + beta D(A, B) + gamma |D(A, C) - D(B, C)|,
$
where $(alpha_A, alpha_B, beta, gamma)$ depend at most on cluster cardinalities $(|A|, |B|, |C|)$.

#subsec("A", "Exact Milligan Monotonicity, Bruynooghe Reducibility, and the NN-Chain Dichotomy")

Two distinct structural properties govern every Lance--Williams method: (i) *monotonicity* (zero dendrogram inversions), ensuring merge heights $h_1 <= h_2 <= dots <= h_(n-1)$ induce a valid ultrametric $u(x, y) = h(x or y)$ without intersecting dendrogram branches; and (ii) *Bruynooghe reducibility*, ensuring that merging two reciprocal nearest neighbors never brings a third cluster strictly closer than before, which is the exact algebraic precondition for the amortized $O(n^2)$ Nearest-Neighbor Chain (NN-Chain) algorithm @mullner2011modern. Although frequently conflated in applied literature, monotonicity and reducibility are *strictly distinct*: Flexible Beta (Alg. 74) separates them sharply at $beta = 0$.

#definition(title: "Universal Monotonicity and Bruynooghe Reducibility")[
  A Lance--Williams update rule $(alpha_A, alpha_B, beta, gamma)$ is:
  1. *Universally Monotone* if for every triple of non-negative dissimilarities $(D(A, B), D(A, C), D(B, C))$ satisfying $D(A, B) <= min{D(A, C), D(B, C)}$,
     $ D(A union B, C) >= D(A, B); $
  2. *Universally Reducible* (in the sense of Bruynooghe) if under the same hypothesis $D(A, B) <= min{D(A, C), D(B, C)}$,
     $ D(A union B, C) >= min{D(A, C), D(B, C)}. $
]

#theorem(title: "Exact Lance--Williams Monotonicity, Reducibility, and NN-Chain Complexity Dichotomy", origin: "Original Unified Proof")[
  Fix arbitrary cluster cardinalities $(|A|, |B|, |C|)$ and real Lance--Williams parameters $(alpha_A, alpha_B, beta, gamma)$.
  1. *(Necessary and Sufficient Monotonicity Conditions):* The recurrence is universally monotone on all non-negative dissimilarities if and only if
     $
       alpha_A + alpha_B + beta >= 1, quad alpha_A + alpha_B >= 0, \
       min{alpha_A, alpha_B} + gamma >= 0.
     $
     (The classical Milligan sufficient conditions $alpha_A >= 0$, $alpha_B >= 0$, $alpha_A + alpha_B + beta >= 1$, $min{alpha_A, alpha_B} + gamma >= 0$ imply these immediately.) Whenever these hold at every merge step, the merge heights satisfy $h_1 <= h_2 <= dots <= h_(n-1)$ and $u(x, y) := h(x or y)$ defines an exact ultrametric on $cal(X)$.
  2. *(Necessary and Sufficient Reducibility Conditions):* The recurrence is universally Bruynooghe-reducible if and only if
     $
       alpha_A + alpha_B + beta >= 1, quad alpha_A + alpha_B >= 1, \
       min{alpha_A, alpha_B} + gamma >= 0.
     $
     Every universally reducible rule is universally monotone, and executes identically to greedy global-minimum agglomeration via the NN-Chain algorithm in $O(n^2)$ worst-case time (and $O(n)$ memory when distances are recomputed on the fly, or $O(n^2)$ memory with an explicit matrix).
  3. *(Exact Classification Across Algs. 56--62 and 74):*
     - *Reducible & Strictly Ultrametric ($O(n^2)$ NN-Chain):* Single Linkage (\#56), Complete Linkage (\#57), UPGMA (\#58), Ward (\#59), WPGMA (\#60), and Flexible Beta (\#74) for all $beta <= 0$.
     - *Monotone & Ultrametric, Yet Non-Reducible ($O(n^2 log n)$ heap required):* Flexible Beta (\#74) for every $beta in (0, 1)$, where $alpha_A + alpha_B = 1 - beta in (0, 1)$ violates reducibility while $alpha_A + alpha_B + beta = 1$ preserves strict monotonicity.
     - *Non-Monotone & Non-Reducible (Inversions, $O(n^2 log n)$ heap required):* Centroid Linkage / UPGMC (\#61) with $alpha_A + alpha_B + beta = 1 - (|A| |B|) / (|A| + |B|)^2 < 1$, and Median Linkage / WPGMC (\#62) with $alpha_A + alpha_B + beta = 3 / 4 < 1$.
]

#proof[
  *Part 1:* Let $d_0 := D(A, B) >= 0$. By hypothesis $D(A, C) >= d_0$ and $D(B, C) >= d_0$, so we may parameterize $D(A, C) = d_0 + u$ and $D(B, C) = d_0 + v$ with arbitrary independent $u, v >= 0$. Consider the branch $u <= v$, writing $delta := v - u >= 0$ so $|D(A, C) - D(B, C)| = delta$. Substituting into the Lance--Williams formula yields the exact identity
  $
    D(A union B, C) - d_0 &= alpha_A (d_0 + u) + alpha_B (d_0 + u + delta) \
      &quad + beta d_0 + gamma delta - d_0 \
      &= (alpha_A + alpha_B + beta - 1) d_0 \
      &quad + (alpha_A + alpha_B) u + (alpha_B + gamma) delta.
  $
  Symmetrically, on the branch $v <= u$ with $delta' := u - v >= 0$,
  $
    D(A union B, C) - d_0 &= (alpha_A + alpha_B + beta - 1) d_0 \
      &quad + (alpha_A + alpha_B) v + (alpha_A + gamma) delta'.
  $
  Because $(d_0, u, delta) in bb(R)_(>= 0)^3$ (resp. $(d_0, v, delta') in bb(R)_(>= 0)^3$) range freely and independently over the non-negative orthant, $D(A union B, C) - d_0 >= 0$ holds universally iff all three linear coefficients are non-negative on both branches: $alpha_A + alpha_B + beta - 1 >= 0$, $alpha_A + alpha_B >= 0$, $alpha_B + gamma >= 0$, and $alpha_A + gamma >= 0$. When satisfied at step $t$ (merging $(A, B)$ at height $h_(t+1) = d_0$), every updated distance satisfies $D(A union B, C) >= d_0$, while every untouched pair $(C, C')$ already satisfied $D(C, C') >= d_0$ by minimality of $(A, B)$. Hence $h_(t+2) = min D^((t+1)) >= h_(t+1)$. Inductively $h_1 <= dots <= h_(n-1)$, and for any three leaves $x, y, z$, if $h(x or z)$ exceeds both $h(x or y)$ and $h(y or z)$, then $x, y$ and $y, z$ would already lie in the same connected cluster prior to height $h(x or z)$, contradicting $x or z$ forming later; thus $u(x, z) <= max{u(x, y), u(y, z)}$.

  *Part 2:* Let $m := min{D(A, C), D(B, C)} >= d_0$ and write $eta := m - d_0 >= 0$ (so $m = d_0 + eta$). On the branch $D(A, C) = m <= D(B, C) = m + delta$ ($delta >= 0$), substituting $D(A, B) = d_0$, $D(A, C) = d_0 + eta$, and $D(B, C) = d_0 + eta + delta$ gives
  $
    D(A union B, C) - m &= (alpha_A + alpha_B + beta - 1) d_0 \
      &quad + (alpha_A + alpha_B - 1) eta + (alpha_B + gamma) delta,
  $
  and symmetrically for $D(B, C) = m <= D(A, C) = m + delta'$:
  $
    D(A union B, C) - m &= (alpha_A + alpha_B + beta - 1) d_0 \
      &quad + (alpha_A + alpha_B - 1) eta + (alpha_A + gamma) delta'.
  $
  Since $(d_0, eta, delta) in bb(R)_(>= 0)^3$ vary freely and independently, universal reducibility $D(A union B, C) - m >= 0$ holds iff $alpha_A + alpha_B + beta >= 1$, $alpha_A + alpha_B >= 1$, and $min{alpha_A, alpha_B} + gamma >= 0$. Because $alpha_A + alpha_B >= 1 > 0$, reducibility strictly implies monotonicity.

  Under reducibility, let the NN-Chain maintain a stack $S = (C_1, dots, C_r)$ where $C_(j+1)$ is a nearest neighbor of $C_j$. When a reciprocal nearest-neighbor (RNN) pair $(C_(r-1), C_r)$ merges into $M = C_(r-1) union C_r$, reducibility guarantees for every predecessor $C_j$ ($j <= r - 2$) that $D(C_j, M) >= min{D(C_j, C_(r-1)), D(C_j, C_r)} >= D(C_j, C_(j+1))$ (as distances strictly decrease along the stack), so the surviving prefix $(C_1, dots, C_(r-2))$ remains a valid NN-chain without rebuilding. Moreover, by reducibility no future merge can produce a dissimilarity strictly smaller than $D(C_(r-1), C_r)$ involving $C_(r-1)$ or $C_r$, so merging $(C_(r-1), C_r)$ immediately yields an exact dendrogram (identical up to ties to global-minimum agglomeration). Amortizing stack operations via the potential $Phi = |S| + 2 |cal(C)|$ shows at most $3 n - 2$ nearest-neighbor scans are executed in total, each costing $O(n)$, establishing $O(n^2)$ total time @mullner2011modern.

  *Part 3:* Direct substitution of Table I verifies all three regimes:
  - *Reducible ($O(n^2)$ NN-Chain):* SLINK ($alpha_A = alpha_B = 1 / 2, beta = 0, gamma = -1 / 2$), CLINK ($gamma = +1 / 2$), UPGMA ($alpha_A = |A| / (|A| + |B|), beta = gamma = 0$), WPGMA ($alpha_A = alpha_B = 1 / 2, beta = gamma = 0$), Ward ($alpha_A = (|A| + |C|) / N, beta = -|C| / N, gamma = 0 => alpha_A + alpha_B = 1 + |C| / N > 1$), and Flexible Beta ($alpha_A = alpha_B = (1 - beta) / 2, gamma = 0$) for $beta <= 0$ ($1 - beta >= 1$).
  - *Monotone, Non-Reducible ($O(n^2 log n)$):* Flexible Beta for $beta in (0, 1)$, where $alpha_A + alpha_B + beta = 1$ yet $alpha_A + alpha_B = 1 - beta < 1$.
  - *Non-Monotone ($O(n^2 log n)$):* UPGMC ($beta = -|A| |B| / (|A| + |B|)^2$) and WPGMC ($beta = -1 / 4$), both violating $alpha_A + alpha_B + beta >= 1$.
]

#subsec("B", "Geometry of Centroid and Median Inversions vs. 1D Dimensional Immunity")

Why do Centroid Linkage (UPGMC, Alg. 61) and Median Linkage (WPGMC, Alg. 62) exhibit dendrogram inversions in $bb(R)^2$, yet remain strictly ultrametric on the real line $bb(R)^1$? And how do their inversion depths behave as cluster sizes diverge? The next theorem provides exact sharp answers.

#theorem(title: "Sharp Inversion Bounds, Size-Masking Asymmetry, and 1D Monotonicity of UPGMC and WPGMC", origin: "Original Unified Proof")[
  Let $cal(X) subset.eq bb(R)^d$ be endowed with squared Euclidean dissimilarity $D(A, B) = norm(tilde(mu)_A - tilde(mu)_B)_2^2$, where $tilde(mu)_C$ is the exact arithmetic centroid $mu_C$ for UPGMC (Alg. 61) and the recursively midpoint-averaged pseudo-centroid for WPGMC (Alg. 62).
  1. *(Sharp Single-Step Inversion Floor):* At any merge $(A, B)$ minimizing $D$ across active clusters, every remaining cluster $C$ obeys
     $
       D_"UPGMC" (A union B, C) &>= (1 - (|A| |B|) / (|A| + |B|)^2) D(A, B) >= 3 / 4 D(A, B), \
       D_"WPGMC" (A union B, C) &>= 3 / 4 D(A, B),
     $
     corresponding to a worst-case unsquared Euclidean height drop factor of $sqrt(3) / 2 approx 0.8660$. Both bounds are sharp on an isosceles triangle $x_1 = (0, 0)$, $x_2 = (1, 0)$, $x_3 = (1 / 2, h)$ with $h in (sqrt(3) / 2, 1)$ (giving $h_1 = 1 > h_2 = h^2 -> (3 / 4)^+$ as $h -> (sqrt(3) / 2)^+$).
  2. *(Asymptotic Size-Imbalance Asymmetry):* Under highly unbalanced merges ($|A| >> |B|$, let $rho = |B| / |A| -> 0^+$), UPGMC's inversion margin vanishes quadratically in $sqrt(rho)$ as $1 - (|A| |B|) / (|A| + |B|)^2 = 1 - rho + O(rho^2) -> 1$ (recovering near-reducibility via *centroid masking* of the smaller cluster $B$), whereas WPGMC's inversion factor remains *persistently* $3 / 4$ regardless of $|A| / |B|$ (eliminating size-domination masking at the cost of persistent inversions at every scale).
  3. *(Exact 1D Reducibility & Non-Inversion):* When $d = 1$ ($cal(X) subset.eq bb(R)$), *both* UPGMC (Alg. 61) and WPGMC (Alg. 62) satisfy Bruynooghe reducibility $D(A union B, C) >= min{D(A, C), D(B, C)}$ at every step, admit zero dendrogram inversions, and coincide with $O(n^2)$ NN-Chain execution.
]

#proof[
  *Part 1:* By Stewart's identity in Hilbert space (setting $lambda = |A| / (|A| + |B|)$ for UPGMC, or $lambda = 1 / 2$ for WPGMC, so $tilde(mu)_(A union B) = lambda tilde(mu)_A + (1 - lambda) tilde(mu)_B$), expanding $norm(lambda (tilde(mu)_A - tilde(mu)_C) + (1 - lambda) (tilde(mu)_B - tilde(mu)_C))_2^2$ via $2 chevron.l u, v chevron.r = norm(u)_2^2 + norm(v)_2^2 - norm(u - v)_2^2$ yields
  $
    D(A union B, C) &= lambda D(A, C) + (1 - lambda) D(B, C) \
      &quad - lambda (1 - lambda) D(A, B).
  $
  Because $D(A, C) >= D(A, B)$ and $D(B, C) >= D(A, B)$,
  $
    D(A union B, C) &>= (lambda + (1 - lambda) - lambda (1 - lambda)) D(A, B) \
      &= (1 - lambda (1 - lambda)) D(A, B).
  $
  Since $lambda (1 - lambda) <= 1 / 4$ with equality at $lambda = 1 / 2$, both linkages satisfy $D(A union B, C) >= (3 / 4) D(A, B)$ ($sqrt(3) / 2$ in unsquared norm). On $x_1 = (0, 0)$, $x_2 = (1, 0)$, $x_3 = (1 / 2, h)$ with $h in (sqrt(3) / 2, 1)$, we have $D(x_1, x_2) = 1 < 1 / 4 + h^2 = D(x_1, x_3) = D(x_2, x_3)$, so both merge $\{x_1, x_2\}$ at squared height $h_1 = 1$ with $tilde(mu)_(1 2) = (1 / 2, 0)$ and then absorb $x_3$ at $h_2 = h^2 < 1$, approaching $3 / 4$ as $h -> (sqrt(3) / 2)^+$.

  *Part 2:* For UPGMC with $|B| = rho |A|$, $lambda (1 - lambda) = rho / (1 + rho)^2 = rho + O(rho^2) -> 0$ and $tilde(mu)_(A union B) - tilde(mu)_A = (rho / (1 + rho)) (tilde(mu)_B - tilde(mu)_A) -> 0$: large cluster $A$ masks $B$ and suppresses inversions. In WPGMC, $lambda = 1 / 2$ identically, so placing $m$ copies at $x_1 = (0, 0)$, $1$ point at $x_2 = (1, 0)$, and $1$ point at $x_3 = (1 / 2, h)$ still sets $tilde(mu)_(1 2) = (1 / 2, 0)$ and reproduces the $h_2 / h_1 = h^2 -> 3 / 4$ inversion for arbitrarily large $m$.

  *Part 3:* Suppose $cal(X) subset.eq bb(R)$ and $(A, B)$ minimizes $D(A, B) = (tilde(mu)_A - tilde(mu)_B)^2$ with $tilde(mu)_A <= tilde(mu)_B$. No active representative $tilde(mu)_C$ can lie in $(tilde(mu)_A, tilde(mu)_B)$ (otherwise $|tilde(mu)_C - tilde(mu)_A| < |tilde(mu)_B - tilde(mu)_A|$, contradicting minimality). Hence either $tilde(mu)_C >= tilde(mu)_B$ or $tilde(mu)_C <= tilde(mu)_A$. Since $tilde(mu)_(A union B) in [tilde(mu)_A, tilde(mu)_B]$, when $tilde(mu)_C >= tilde(mu)_B$ we obtain $tilde(mu)_C - tilde(mu)_(A union B) >= tilde(mu)_C - tilde(mu)_B >= 0$, whence
  $
    D(A union B, C) &= (tilde(mu)_C - tilde(mu)_(A union B))^2 \
      &>= (tilde(mu)_C - tilde(mu)_B)^2 = min{D(A, C), D(B, C)}.
  $
  The case $tilde(mu)_C <= tilde(mu)_A$ is symmetric, establishing exact 1D reducibility and showing that inversions are strictly multidimensional ($d >= 2$).
]

#subsec("C", "Variance Conservation Laws and Space-Dilation Ordering")

Beyond monotonicity, Euclidean linkages differ dramatically in *space distortion*—whether merging $(A, B)$ makes the composite cluster appear artificially *closer* to third clusters (*space-contracting*, inducing chaining) or *farther* (*space-dilating*, favoring compact equal-size spheres). The next unifying theorem links UPGMA, UPGMC, Ward, WPGMA, WPGMC, and Flexible Beta through exact algebraic variance identities.

#theorem(title: "Huygens--Steiner Linkage Identities and Exact Space-Dilation Ordering", origin: "Original Unified Proof")[
  Let $cal(X) subset.eq bb(R)^d$, and for any finite cluster $C subset.eq cal(X)$ let $mu_C = |C|^(-1) sum_(x in C) x$ and $sigma^2 (C) = |C|^(-1) sum_(x in C) norm(x - mu_C)_2^2$ denote its mean and internal variance.
  1. *(Exact UPGMA--UPGMC--Ward Variance Triangle):* Under squared Euclidean base distances $d^2 (x, y) = norm(x - y)_2^2$, every disjoint pair $(A, B)$ obeys
     $
       D_"UPGMA"^((2)) (A, B) &= D_"UPGMC"^((2)) (A, B) + sigma^2 (A) + sigma^2 (B), \
       D_"Ward" (A, B) &= (|A| |B|) / (|A| + |B|) D_"UPGMC"^((2)) (A, B).
     $
  2. *(Linkage Sandwich & Balanced-Tree Coincidence):* Whenever $D(A, B) <= min{D(A, C), D(B, C)}$, both $D_"UPGMA"$ and $D_"WPGMA"$ lie inside $[D_"SL" (A union B, C), D_"CL" (A union B, C)]$. On perfectly balanced binary subtrees, WPGMA (\#60) coincides identically with UPGMA (\#58) and WPGMC (\#62) with UPGMC (\#61); on a right-caterpillar of depth $m$, WPGMA assigns leaf $x_1$ weight $2^(-(m-1))$ vs. $1 / 2$ for $x_m$, amplifying late-arriving outliers by $2^(m-2)$.
  3. *(Flexible Beta Monotone Dilation & Local Ward Matching):* For Flexible Beta (Alg. 74, $alpha_A = alpha_B = (1 - beta) / 2$, $gamma = 0$, $beta < 1$),
     $
       &(partial D_"Flex"^((beta)) (A union B, C)) / (partial beta) \
       &quad = D(A, B) - 1 / 2 (D(A, C) + D(B, C)) <= 0.
     $
     Hence decreasing $beta < 0$ monotonically dilates space (matching Ward locally at $beta = -|C| / (2 |A| + |C|)$ when $|A| = |B|$, giving $-1 / 3$ at equal sizes), $beta = 0$ recovers WPGMA as the *unique* value preserving the metric triangle inequality universally, and $beta in (0, 1)$ contracts space while forfeiting reducibility.
]

#proof[
  *Part 1:* For $x in A$ and $y in B$, decompose $x - y = (x - mu_A) - (y - mu_B) + (mu_A - mu_B)$. Expanding the squared norm and averaging over $(x, y) in A times B$, all three cross inner products vanish identically because $sum_(x in A) (x - mu_A) = bold(0)$ and $sum_(y in B) (y - mu_B) = bold(0)$, leaving
  $
    D_"UPGMA"^((2)) (A, B) &= 1 / (|A| |B|) sum_(x in A, y in B) norm(x - y)_2^2 \
      &= norm(mu_A - mu_B)_2^2 + sigma^2 (A) + sigma^2 (B).
  $
  For Ward's ESS increment, set $N = |A| + |B|$ and note $mu_(A union B) = (|A| mu_A + |B| mu_B) / N$, so $mu_A - mu_(A union B) = (|B| / N)(mu_A - mu_B)$ and $mu_B - mu_(A union B) = -(|A| / N)(mu_A - mu_B)$. Applying the parallel-axis theorem to $A$ and $B$ around $mu_(A union B)$ gives
  $
    Delta "ESS"(A, B) &= |A| norm(mu_A - mu_(A union B))_2^2 \
      &quad + |B| norm(mu_B - mu_(A union B))_2^2 \
      &= (|A| |B|) / (|A| + |B|) norm(mu_A - mu_B)_2^2.
  $

  *Part 2:* Since $D_"SL" (A union B, C) = min{D(A, C), D(B, C)}$ and $D_"CL" (A union B, C) = max{D(A, C), D(B, C)}$, any convex combination $lambda D(A, C) + (1 - lambda) D(B, C)$ with $lambda in [0, 1]$ (including UPGMA and WPGMA) lies inside $[D_"SL", D_"CL"]$. By induction on merge depth, WPGMA assigns leaf $x in C$ dyadic weight $2^(-"depth"_C (x))$ (summing to $1$ by Kraft's equality for full binary trees), whereas UPGMA assigns uniform weight $|C|^(-1)$. When all leaves of $C$ sit at equal depth $h$, $2^(-h) = 2^(-h) = |C|^(-1)$, making both the pairwise averages (WPGMA vs. UPGMA) and the weighted averages of leaf vectors (WPGMC vs. UPGMC) coincide identically. On a caterpillar $(((x_1, x_2), x_3), dots, x_m)$, $"depth"(x_1) = m - 1$ while $"depth"(x_m) = 1$, giving weight ratio $2^(-1) / 2^(-(m-1)) = 2^(m-2)$.

  *Part 3:* Differentiating $D_"Flex"^((beta)) (A union B, C) = (1 / 2)(D(A, C) + D(B, C)) + beta [D(A, B) - (1 / 2)(D(A, C) + D(B, C))]$ with respect to $beta$ gives $D(A, B) - (1 / 2)(D(A, C) + D(B, C)) <= 0$ whenever $D(A, B) <= min{D(A, C), D(B, C)}$. At $beta = 0$, $alpha_A = alpha_B = 1 / 2, beta = gamma = 0$ is exact WPGMA, which preserves the triangle inequality because averaging $d(a, c) <= d(a, z) + d(z, c)$ and $d(b, c) <= d(b, z) + d(z, c)$ preserves $<=$. Conversely, for every $beta < 0$, take four collinear points $A = {0}, B = {1}, C = {2}, E = {3}$ on $bb(R)$: merging $(A, B)$ at distance $1$ yields $D(A union B, C) = (3 - beta) / 2$, $D(A union B, E) = (5 - 3 beta) / 2$, and $D(C, E) = 1$, so $D(A union B, E) - D(A union B, C) - D(C, E) = -beta > 0$, strictly violating the triangle inequality. For every $beta in (0, 1)$, take $C = {-2}, A = {0}, B = {1 / 2}, E = {2}$ on $bb(R)$: merging the unique closest pair $(A, B)$ ($d = 1 / 2$) gives $D(A union B, C) = (9 - 7 beta) / 4$ and $D(A union B, E) = (7 - 5 beta) / 4$, whose sum $4 - 3 beta < 4 = D(C, E)$ again violates the triangle inequality. Thus $beta = 0$ is the unique triangle-inequality-preserving parameter. Finally, when $|A| = |B| = m$ and $|C| = c$, Ward sets $alpha_A = alpha_B = (m + c) / (2 m + c) = (1 - beta_"Ward") / 2$ with $beta_"Ward" = -c / (2 m + c)$ (equaling $-1 / 3$ at $m = c$).
]

Theorem 2(1) explains *both* Ward's celebrated practical strength on noisy spherical data and its classic weakness on imbalanced clusters: multiplying squared centroid separation by $(|A| |B|) / (|A| + |B|)$ penalizes merging two large clusters of size $m$ by factor $m / 2$, creating a strong equal-size bias that resists outlier chaining yet splits large elongated clusters prematurely before absorbing tiny neighboring clusters. Conversely, UPGMA adds *additive* variance $sigma^2 (A) + sigma^2 (B)$ without size multiplication, balancing compactness and cluster-size imbalance cleanly.

#wide_figure(
  caption: [*Table I:* Exact Lance--Williams Parameterization, Monotonicity, Bruynooghe Reducibility, Space Distortion, and Complexity Across All 9 Combinatorial Agglomerative Algorithms (Algs. 56--62, 72, 74). Here $n_A = |A|$, $n_B = |B|$, $n_C = |C|$, $N = n_A + n_B + n_C$.],
  {
    set text(size: 7.4pt)
    set table(inset: (x: 3.8pt, y: 2.8pt))
    table(
      columns: (1.22fr, 0.85fr, 0.85fr, 0.72fr, 0.56fr, 0.60fr, 0.70fr, 0.85fr),
      align: (left, center, center, center, center, center, center, center),
      [Algorithm], [$alpha_A$], [$alpha_B$], [$beta$], [$gamma$], [Monotone?], [Reducible?], [Best Time],
      [\#56 SLINK @sibson1973slink], [$1 / 2$], [$1 / 2$], [$0$], [$-1 / 2$], [Yes], [Yes], [$O(n^2)$ (SLINK)],
      [\#57 CLINK @defays1977clink], [$1 / 2$], [$1 / 2$], [$0$], [$+1 / 2$], [Yes], [Yes], [$O(n^2)$ (NN-Ch.)],
      [\#58 UPGMA @sokal1958statistical], [$n_A / (n_A + n_B)$], [$n_B / (n_A + n_B)$], [$0$], [$0$], [Yes], [Yes], [$O(n^2)$ (NN-Ch.)],
      [\#59 Ward @ward1963hierarchical], [$(n_A + n_C) / N$], [$(n_B + n_C) / N$], [$-n_C / N$], [$0$], [Yes], [Yes], [$O(n^2)$ (NN-Ch.)],
      [\#60 WPGMA @sokal1958statistical], [$1 / 2$], [$1 / 2$], [$0$], [$0$], [Yes], [Yes], [$O(n^2)$ (NN-Ch.)],
      [\#61 UPGMC @sokal1958statistical], [$n_A / (n_A + n_B)$], [$n_B / (n_A + n_B)$], [$-n_A n_B / (n_A + n_B)^2$], [$0$], [No ($d >= 2$)], [No ($d >= 2$)], [$O(n^2 log n)$],
      [\#62 WPGMC @lance1967general], [$1 / 2$], [$1 / 2$], [$-1 / 4$], [$0$], [No ($d >= 2$)], [No ($d >= 2$)], [$O(n^2 log n)$],
      [\#74 Flex. $beta <= 0$ @lance1967general], [$(1 - beta) / 2$], [$(1 - beta) / 2$], [$beta <= 0$], [$0$], [Yes], [Yes], [$O(n^2)$ (NN-Ch.)],
      [\#74 Flex. $beta in (0, 1)$], [$(1 - beta) / 2$], [$(1 - beta) / 2$], [$beta > 0$], [$0$], [Yes], [No], [$O(n^2 log n)$],
      [\#72 Lance--Williams], [General], [General], [General], [General], [Thm. 1(1)], [Thm. 1(2)], [$O(n^2)$ / $O(n^2 log n)$],
    )
  },
)

// =============================================================================
// SECTION II
// =============================================================================
#sec("II", "Chaining, Diameter, Anti-Chaining, and Dasgupta Cost Guarantees")

While Section I separated linkages by local recurrence algebra, practical hierarchical performance is dictated by *global geometric and combinatorial pathologies*: Single Linkage's zero-breakdown chaining across thin noise filaments, Complete Linkage's greedy diameter suboptimality, Minimax Linkage's escape from the Lance--Williams class via prototype-anchored metric balls @bien2011minimax, Genie's economic inequality control over MST contractions @gagolewski2016genie, and global hierarchy optimality measured by Dasgupta's cost objective @dasgupta2016cost and Moseley--Wang's dual revenue @moseley2017approximation.

#subsec("A", "Ultrametric Extremality, Greedy Diameter Failure, and Minimax Non-Combinatoriality")

For any finite metric space $(cal(X), d)$, let $cal(U)_-(d) = {u " ultrametric" : u <= d}$ and $cal(U)_+(d) = {u " ultrametric" : u >= d}$ denote the lower and upper ultrametric cones. Single Linkage (Alg. 56) outputs the *unique maximum element* of $cal(U)_-(d)$—the *subdominant ultrametric* $u_"SL" (x, y) = min_(P : x -> y) max_(e in P) d(e)$, realized equivalently by the bottleneck path along any Minimum Spanning Tree (MST) @gower1969minimum @carlsson2010characterization. Complete Linkage (Alg. 57) outputs a *minimal element* of $cal(U)_+(d)$ (a superdominant ultrametric where every height-$h$ cluster is a maximal clique in the $h$-threshold graph), which is *not* unique in general. Sitting strictly between diameter-based Complete Linkage and prototype-based medoids is *Minimax Linkage* (Alg. 73, @bien2011minimax), which scores every candidate merge by the minimum covering-ball radius anchored at an actual member of the merged cluster:
$
  r_"MM" (C) := min_(p in C) max_(x in C) d(p, x), quad D_"MM" (A, B) := r_"MM" (A union B).
$

#theorem(title: "Minimax Reducibility, Non-Lance--Williams Impossibility, and CLINK Diameter Suboptimality", origin: "Original Unified Proof")[
  Let $(cal(X), d)$ be any finite metric space.
  1. *(Unconditional Reducibility & Tight Factor-2 Metric Sandwich of Minimax Linkage):* Minimax Linkage (Alg. 73) satisfies Bruynooghe reducibility $D_"MM" (A union B, C) >= min{D_"MM" (A, C), D_"MM" (B, C)}$ unconditionally on every symmetric dissimilarity space (without requiring the triangle inequality), hence is strictly ultrametric with zero inversions and monotone-transform invariant. On any metric space, every cluster $C$ satisfies
     $
       1 / 2 "diam"(C) <= r_"MM" (C) <= "diam"(C),
     $
     with both bounds tight, and every level-$h$ cut yields clusters covered by closed metric balls $B(p_r, h)$ centered at data prototypes $p_r in C_r$.
  2. *(Strict Non-Membership of Minimax Linkage in the Lance--Williams Family):* There exists *no* combinatorial update function $F(D(A, B), D(A, C), D(B, C), |A|, |B|, |C|)$—and in particular no Lance--Williams parameterization—that computes $D_"MM" (A union B, C)$.
  3. *(Constructive Chaining Collapse of SLINK vs. Greedy Diameter Suboptimality of CLINK):* Single Linkage (\#56) has finite-sample breakdown point $0$ (a single equispaced chain of $m$ noise points bridges two clusters of arbitrary separation $L$ at merge height $L / (m + 1)$). Conversely, Complete Linkage (\#57) avoids single-point chaining yet fails to minimize maximum cluster diameter even on $n = 4$ points in $bb(R)^1$: for every $epsilon in (0, 1 / 2)$, on $cal(X)_epsilon = {0, 2 - epsilon, 3 - 2 epsilon, 5 - 3 epsilon}$ with $k = 2$, CLINK produces maximum diameter $3 - 2 epsilon$ against optimal $2 - epsilon$, attaining ratio $(3 - 2 epsilon) / (2 - epsilon) -> (3 / 2)^-$.
]

#proof[
  *Part 1:* Let $p^* in A union B union C$ achieve $r_"MM" (A union B union C) = max_(x in A union B union C) d(p^*, x)$. Since $p^*$ lies in either $A union C$ or $B union C$:
  - If $p^* in A union C$, then $r_"MM" (A union C) <= max_(x in A union C) d(p^*, x) <= max_(x in A union B union C) d(p^*, x) = r_"MM" (A union B union C)$.
  - If $p^* in B union C$, symmetrically $r_"MM" (B union C) <= r_"MM" (A union B union C)$.
  Hence $D_"MM" (A union B, C) >= min{D_"MM" (A, C), D_"MM" (B, C)}$ holds *unconditionally*—even without assuming $D_"MM" (A, B) <= min{D_"MM" (A, C), D_"MM" (B, C)}$! By Theorem 1(2), Minimax Linkage is strictly inversion-free and ultrametric. For any $x, y in C$ and optimal prototype $p^* in C$, the triangle inequality gives $d(x, y) <= d(x, p^*) + d(p^*, y) <= 2 r_"MM" (C)$, proving $(1 / 2) "diam"(C) <= r_"MM" (C)$; conversely, $r_"MM" (C) <= max_(x, y in C) d(x, y) = "diam"(C)$. Equality $r_"MM" = (1 / 2) "diam"$ holds on $\{-1, 0, 1\} subset.eq bb(R)$ (prototype $0$), while $r_"MM" = "diam"$ holds on any regular simplex with equal pairwise distances.

  *Part 2:* Construct two 4-point finite metric spaces $(cal(X), d_1)$ and $(cal(X), d_2)$ partitioned identically into $A = {a}$, $B = {b}$, and $C = {c_1, c_2}$ (so $|A| = 1$, $|B| = 1$, $|C| = 2$):
  - In $d_1$: set $d_1 (a, b) = 1$ and all other five pairwise distances equal to $2$. Every triangle has side lengths $(1, 2, 2)$ or $(2, 2, 2)$, so $d_1$ is a valid metric. Direct evaluation gives $D_"MM"^((1)) (A, B) = 1$, $D_"MM"^((1)) (A, C) = 2$, $D_"MM"^((1)) (B, C) = 2$, and $D_"MM"^((1)) (A union B, C) = 2$ (witnessed by prototype $a$, since $max{d_1 (a, b), d_1 (a, c_1), d_1 (a, c_2)} = 2$).
  - In $d_2$: set $d_2 (a, b) = 1$, $d_2 (c_1, c_2) = 2$, $d_2 (c_1, a) = d_2 (c_2, b) = 2$, and $d_2 (c_1, b) = d_2 (c_2, a) = 3$. All four triangles have side triples $(1, 2, 3)$ or $(2, 2, 3)$, obeying the triangle inequality identically. Compute:
    - $D_"MM"^((2)) (A, B) = d_2 (a, b) = 1$;
    - $D_"MM"^((2)) (A, C) = r_"MM" ({a, c_1, c_2}) = 2$ (uniquely achieved at prototype $c_1$, since $max{d_2 (c_1, a), d_2 (c_1, c_2)} = 2$);
    - $D_"MM"^((2)) (B, C) = r_"MM" ({b, c_1, c_2}) = 2$ (uniquely achieved at prototype $c_2$, since $max{d_2 (c_2, b), d_2 (c_2, c_1)} = 2$);
    - Yet in $A union B union C = {a, b, c_1, c_2}$, *every* candidate prototype has distance $3$ to one other member ($d_2 (a, c_2) = d_2 (b, c_1) = 3$), yielding $D_"MM"^((2)) (A union B, C) = 3 != 2 = D_"MM"^((1)) (A union B, C)$!
  Because $(D(A, B), D(A, C), D(B, C), |A|, |B|, |C|) = (1, 2, 2, 1, 1, 2)$ is identical across both metrics while $D_"MM" (A union B, C)$ differs, no combinatorial update rule $F$—and a fortiori no Lance--Williams formula—can represent Minimax Linkage.

  *Part 3:* For SLINK, placing two dense clusters $K_1, K_2$ of internal spacing $delta_0$ at separation $L >> delta_0$ and inserting $m = ceil(L / delta_0)$ equally spaced bridge points at step $delta < delta_0$ connects $K_1$ and $K_2$ in the MST at bottleneck height $delta < delta_0$, merging $K_1$ and $K_2$ before either cluster forms internally. For CLINK on $cal(X)_epsilon = {0, 2 - epsilon, 3 - 2 epsilon, 5 - 3 epsilon}$ ($epsilon in (0, 1 / 2)$), the adjacent gaps are $(2 - epsilon, 1 - epsilon, 2 - epsilon)$. Step 1 merges the strictly closest central pair ${2 - epsilon, 3 - 2 epsilon}$ at distance $1 - epsilon$. At step 2, merging either extreme endpoint yields diameter $3 - 2 epsilon$, producing the $k = 2$ partition ${{0, 2 - epsilon, 3 - 2 epsilon}, {5 - 3 epsilon}}$ of max diameter $3 - 2 epsilon$, whereas the optimal $k = 2$ partition ${{0, 2 - epsilon}, {3 - 2 epsilon, 5 - 3 epsilon}}$ achieves max diameter $2 - epsilon$.
]

#remark(title: "Computational Price of Non-Combinatorial Prototypes")[
  Theorem 3(2) explains why Minimax Linkage (Alg. 73) runs in $O(n^3)$ worst-case time rather than $O(n^2)$: because $r_"MM" (A union B)$ cannot be updated in $O(1)$ from pairwise cluster distances alone, computing a trial merge requires scanning distances across all candidate prototypes in $A union B$ (though storing per-point maximum distances within each cluster accelerates practical runs substantially via NN-Chain, since Minimax *is* reducible!).
]

#subsec("B", "Gini-Controlled Anti-Chaining in Genie (Alg. 69)")

Genie (Alg. 69, @gagolewski2016genie) preserves Single Linkage's non-convex manifold tracking and MST speed while suppressing outlier isolation via the normalized cluster-size Gini index
$
  G(c) &= (sum_(i=1)^k sum_(j=1)^k |c_i - c_j|) / (2 (k - 1) n) \
    &= (sum_(r=1)^k (2 r - k - 1) c_((r))) / ((k - 1) n).
$
Whenever $G(c) <= g$, Genie performs an unconstrained Single-Linkage step; whenever $G(c) > g$, it restricts the minimum single-linkage merge to pairs containing a minimum-size cluster ($min{c_i, c_j} = c_"min"$).

#theorem(title: "Genie MST Sufficiency, Exact Outlier-Absorption Threshold, and Inevitable Inversions", origin: "Original Unified Proof")[
  Let $(cal(X), d)$ have distinct pairwise dissimilarities (or fixed tie-breaking) with unique MST $T^* = (cal(X), E^*)$, and fix $g in (0, 1]$.
  1. *(Exact MST Reduction & $O(n log n)$ Post-MST Complexity):* For every $g in (0, 1]$ and step $t in {0, dots, n - 2}$, Genie's selected merge $(C_(i^*), C_(j^*))$ is witnessed by an edge $e^* in E^*$ of $T^*$; hence once $T^*$ is computed, the entire Genie hierarchy across *any* $g$ runs in $O(n log n)$ time @gagolewski2016genie.
  2. *(Exact Anti-Chaining / Outlier-Absorption Threshold):* Suppose $cal(X)$ comprises two clean clusters $C_1, C_2$ of size $m >= 2$ (internal MST edges $< Delta_0$) separated by gap $Delta > Delta_0$, contaminated by $r >= 1$ isolated singleton outliers at distance $> Delta$ from $C_1 union C_2$ ($n = 2 m + r$). Whereas Single Linkage ($g = 1$) merges $C_1$ and $C_2$ before absorbing any outlier, for *every* threshold
     $
       g < min { (n - 3) / (2 n), space (2 r (m - 1)) / ((r + 1) n) }
     $
     (in particular, the default $g = 0.3$ succeeds for *every* $m >= 4$ and $1 <= r <= m$, where the bound exceeds $1 / 3$), Genie strictly blocks the merge $(C_1, C_2)$ until every outlier is absorbed, recovering the exact two-cluster macro-partition at $K = 2$.
  3. *(Inevitable Inversions & Balanced-Bridge Vulnerability):* For every $g in (0, 1)$, Genie admits strict merge-height inversions $h_(t+1) < h_t$ (destroying ultrametricity). Moreover, because $G(c)$ inspects only cardinalities, a sparse bridge growing at the same rate as two clusters ($G(c) <= g$) still induces Single-Linkage chaining.
]

#proof[
  *Part 1:* Inductively assume all contracted edges prior to step $t$ lie in $E^*$. Let $e^* = (u^*, v^*)$ realize the merge $(C_(i^*), C_(j^*))$. If $G(c) <= g$, $e^*$ is the globally lightest edge crossing $cal(C)^((t))$; if $G(c) > g$ with $c_(i^*) = c_"min"$, $e^*$ is the lightest edge leaving $C_(i^*)$ across $(C_(i^*), cal(X) without C_(i^*))$. In both cases the MST cut property forces $e^* in E^*$. Maintaining the Gini numerator $S = sum_(r=1)^k (2 r - k - 1) c_((r))$ in an order-statistic BST and incident $T^*$-edges in meldable heaps updates each contraction in $O(log n)$ amortized time.

  *Part 2:* Because all internal MST edges of $C_1, C_2$ have weight $< Delta_0 < Delta$, the cross-cluster bridge edge of weight $Delta$ between $C_1$ and $C_2$ can only become eligible once $C_1$ and $C_2$ have assembled into two components of sizes $m_1, m_2 >= m$ alongside $r' in {0, 1, dots, r}$ unabsorbed singleton outliers ($m_1 + m_2 = n - r'$). If $r' = 0$, all outliers are already absorbed and the $K = 2$ partition is exact. If $r' >= 1$, the active size vector $(1, dots, 1, m_1, m_2)$ ($k = r' + 2$) has exact Gini index
  $
    G(1, dots, 1, m_1, m_2) &= (2 r' (m_1 + m_2 - 2) + 2 |m_1 - m_2|) / (2 (r' + 1) n) \
      &>= (r' (n - r' - 2)) / ((r' + 1) n) =: phi(r').
  $
  Rewriting $phi(x) = n^(-1)[(n - 1) - (x + 1) - (n - 1) / (x + 1)]$ gives $phi''(x) = -2 (n - 1) / (n (x + 1)^3) < 0$, so $phi$ is strictly concave on $[1, r]$ and attains its minimum at an endpoint: $phi(1) = (n - 3) / (2 n)$ or $phi(r) = (2 r (m - 1)) / ((r + 1) n)$. When $g < min{phi(1), phi(r)}$, every state with $r' >= 1$ satisfies $G > g$, excluding $(C_1, C_2)$ (since $min{m_1, m_2} >= m > 1 = c_"min"$) and forcing singleton absorption until $r' = 0$. For $m >= 4$ and $1 <= r <= m$, $phi(r) - phi(1) = (2 r m - 2 m - r^2 - 2 r + 3) / (2(r + 1) n) >= 0$ and $phi(1) >= 1 / 3 > 0.3$.

  *Part 3:* Fix $g in (0, 1)$ and choose $r = 1$ (if $g < 1 / 2$) or large $r, m$ so $min{phi(1), phi(r)} > g$, placing all $r$ outliers at distance $D_"out" > Delta$ from $C_1, C_2$. By Part 2, Genie contracts an outlier edge at height $h_t = D_"out"$ before contracting $(C_1, C_2)$ at the final step of strictly smaller height $h_(n-1) = Delta < D_"out"$, producing a strict merge-height inversion. Conversely, linking two clusters of size $m$ by a chain of $m$ bridge points with spacing $delta < Delta_0$ lets the chain grow in lockstep with $C_1, C_2$ while keeping $G(c) <= g$, bypassing the Gini gate completely.
]

#subsec("C", "Dasgupta Cost & Moseley--Wang Revenue Separation")

How do classical linkages compare under a principled *global* objective over the entire tree topology, rather than cut-by-cut heuristics? For a symmetric non-negative similarity graph $G = (cal(X), w)$ ($w(i, j) >= 0$), Dasgupta @dasgupta2016cost defined the hierarchical cost and Moseley and Wang @moseley2017approximation defined its exact dual revenue:
$
  cal(L)_"Dasg" (cal(T); w) &= sum_(i < j) w(i, j) |cal(T)[i or j]|, \
  cal(R)_"MW" (cal(T); w) &= sum_(i < j) w(i, j) (n - |cal(T)[i or j]|),
$
where $cal(T)[i or j]$ is the subtree rooted at the lowest common ancestor $i or j$ of leaves $i$ and $j$. Because $cal(L)_"Dasg" (cal(T)) + cal(R)_"MW" (cal(T)) = n sum_(i < j) w(i, j)$ is constant across all trees, minimizing $cal(L)_"Dasg"$ and maximizing $cal(R)_"MW"$ share the exact same global optimum $cal(T)^*$, yet differ sharply in multiplicative approximability. The next theorem provides a complete, self-contained separation proving UPGMA's exact $1 / 3$-revenue guarantee while shattering *both* Single Linkage and Complete Linkage simultaneously on both objectives.

#theorem(title: "UPGMA 1/3-Revenue Optimality vs. Polynomial Collapse of Single and Complete Linkage", origin: "Original Unified Proof")[
  Let $G = (cal(X), w)$ be any $n$-vertex similarity graph with $w(i, j) >= 0$, and let $cal(T)^*$ be an optimal tree minimizing $cal(L)_"Dasg"$ (equivalently maximizing $cal(R)_"MW"$).
  1. *(Exact $1 / 3$-Approximation of UPGMA, @moseley2017approximation):* Average Linkage (UPGMA, Alg. 58, maximizing $w(A, B) / (|A| |B|)$ at each step) unconditionally satisfies
     $ cal(R)_"MW" (cal(T)_"UPGMA"; w) >= 1 / 3 cal(R)_"MW" (cal(T)^*; w). $
  2. *(Polynomial Collapse of Single Linkage on Both Objectives):* For every $m >= 4$, setting $M = ceil(m^(3 / 2))$ and $n = M + m = Theta(m^(3 / 2))$, there exists a deterministic similarity graph $G_"SL" (m)$ where UPGMA achieves $(1 - o(1)) cal(R)_"MW" (cal(T)^*)$, whereas Single Linkage (Alg. 56) simultaneously suffers
     $
       (cal(L)_"Dasg" (cal(T)_"SL")) / (cal(L)_"Dasg" (cal(T)^*)) &= Omega(n^(1 / 3)), \
       (cal(R)_"MW" (cal(T)_"SL")) / (cal(R)_"MW" (cal(T)^*)) &= O(n^(-1 / 3)).
     $
  3. *(Polynomial Collapse of Complete Linkage on Both Objectives):* For every $m >= 4$, setting $M = m^2$ and $n = M + 2 m = Theta(m^2)$, there exists a deterministic tie-free similarity graph $G_"CL" (m)$ where UPGMA achieves $(1 - o(1)) cal(R)_"MW" (cal(T)^*)$, whereas Complete Linkage (Alg. 57) simultaneously suffers
     $
       (cal(L)_"Dasg" (cal(T)_"CL")) / (cal(L)_"Dasg" (cal(T)^*)) &= Omega(sqrt(n)), \
       (cal(R)_"MW" (cal(T)_"CL")) / (cal(R)_"MW" (cal(T)^*)) &= O(n^(-1 / 2)).
     $
]

#proof[
  *Part 1:* For any tree $cal(T)$ and pair $\{i, j\}$, $n - |cal(T)[i or j]|$ counts third leaves $k in cal(X) without {i, j}$ outside $cal(T)[i or j]$, giving $cal(R)_"MW" (cal(T)) = sum_({i, j, k}) "rev"_cal(T) (i, j, k)$, where $"rev"_cal(T) (i, j, k) = w(i, j)$ if $k in.not cal(T)[i or j]$. Every unordered triple $\{i, j, k\}$ is resolved at the unique UPGMA merge $(A, B)$ joining two members ($i in A, j in B$) while the third lies in an active cluster $C in cal(C) without {A, B}$. Across $A times B times C$, UPGMA earns exact revenue $"rev"_"UPGMA" (A, B; C) = |C| w(A, B)$, while any tree $cal(T)^*$ earns at most $|C| w(A, B) + |B| w(A, C) + |A| w(B, C)$ (since $w >= 0$). UPGMA's greedy choice $(A, B)$ over $(A, C)$ and $(B, C)$ enforces
  $
    w(A, B) / (|A| |B|) >= max { w(A, C) / (|A| |C|), space w(B, C) / (|B| |C|) },
  $
  so $|C| w(A, B) >= |B| w(A, C)$ and $|C| w(A, B) >= |A| w(B, C)$. Summing with $|C| w(A, B) = |C| w(A, B)$ yields $3 "rev"_"UPGMA" (A, B; C) >= "rev"_(cal(T)^*) (A, B, C)$; summing over all merges $(A, B)$ and third clusters $C$ partitions $binom(cal(X), 3)$ disjointly and proves $cal(R)_"MW" (cal(T)_"UPGMA") >= (1 / 3) cal(R)_"MW" (cal(T)^*)$ @moseley2017approximation.

  *Part 2 (Single Linkage):* Define $G_"SL" (m)$ on $V = P union.sq Q$ ($|P| = M = ceil(m^(3 / 2))$, $|Q| = m$, $epsilon = 1 / 4$) where $P$ is a path with weights $w(p_r, p_(r+1)) = 1 + 2 epsilon$, $Q$ is a unit clique $K_m$, and $p_M$ connects to every $q_j in Q$ with weight $1 + epsilon$ (all other weights $0$). Separating $P$ and $Q$ at the root and bisecting $P$ evenly gives $cal(L)_"Dasg" (cal(T)^*) = (m^3 - m) / 3 + O(M log M) + m (M + m)(1 + epsilon) = (1 / 3 + o(1)) m^3$ and $cal(R)_"MW" (cal(T)^*) >= binom(m, 2) M = Omega(m^(7 / 2))$. Single Linkage first merges all of $P$ (weight $1 + 2 epsilon$) and then absorbs $q_1, dots, q_m$ one by one into $P$ via $(p_M, q_j)$ (weight $1 + epsilon > 1$), giving every $Q$-edge LCA size $>= M + 2$: thus $cal(L)_"Dasg" (cal(T)_"SL") >= binom(m, 2) M = Omega(m^(7 / 2))$, while revenue is bounded by $binom(m, 2) m + M (M + m)(1 + 2 epsilon) + m^2 (1 + epsilon) = O(m^3)$, yielding ratios $Omega(m^(1 / 2)) = Omega(n^(1 / 3))$ and $O(n^(-1 / 3))$. Conversely, once UPGMA pairs adjacent path nodes in $P$, their average similarity to neighboring path nodes drops to $(1 + 2 epsilon) / 2 = 3 / 4 < 1$, so UPGMA immediately merges $Q$ intact and attains $(1 - o(1)) cal(R)_"MW" (cal(T)^*)$.

  *Part 3 (Complete Linkage):* Define $G_"CL" (m)$ on $V = Q union.sq I$ where $Q = L union.sq R$ is $K_(m, m)$ ($|L| = |R| = m$, cross weights $1$, internal $0$) and $I$ has size $M = m^2$ with internal weights $2 delta$ and $I times Q$ weights $delta := n^(-4)$ (total $delta$-weight $o(1)$). Separating $Q$ from $I$ at the root gives $cal(L)_"Dasg" (cal(T)^*) <= 2 m^3 + o(1)$ and $cal(R)_"MW" (cal(T)^*) >= m^4$. Complete Linkage first merges $m$ disjoint pairs $P_r = {u_r, v_r}$ (sim. $1$); because $w(u_r, u_s) = 0$, every pair $(P_r, P_s)$ drops to minimum similarity $0 < delta < 2 delta$, forcing Complete Linkage to merge all of $I$ (sim. $2 delta$) and then absorb $P_1, dots, P_m$ one by one into $I$ (sim. $delta$). Every one of the $m^2 - m$ cross-pair edges of $K_(m, m)$ has LCA size $> m^2$, giving $cal(L)_"Dasg" (cal(T)_"CL") = Omega(m^4)$ and $cal(R)_"MW" (cal(T)_"CL") <= m n + 2 m^3 + o(1) = O(m^3)$, hence ratios $Omega(sqrt(n))$ and $O(n^(-1 / 2))$. Under UPGMA, pairs $(P_r, P_s)$ retain average similarity $1 / 2 >> 2 delta$, recovering $Q$ intact.
]

Theorem 5 illuminates a deep structural lesson: although Ward's method (Alg. 59) minimizes greedy SSE locally in Euclidean space (and acts as a 2-approximation counterexample in 1D when greedy merges split an optimal cluster), and Minimax/CLINK control metric diameters, *UPGMA (Alg. 58)* stands uniquely among classical linkages as the only rule enjoying a constant-factor global tree-revenue guarantee across arbitrary similarity graphs @moseley2017approximation.

#wide_figure(
  caption: [*Table II:* Head-to-Head Comparison of Chaining Resistance, Ultrametric Geometry, Prototype Anchoring, and Global Dasgupta / Moseley--Wang Guarantees (Algs. 56--59, 69, 73).],
  {
    set text(size: 7.4pt)
    set table(inset: (x: 3.8pt, y: 2.8pt))
    table(
      columns: (1.05fr, 0.90fr, 0.72fr, 0.95fr, 1.05fr, 1.15fr),
      align: (left, center, center, center, center, left),
      [Algorithm], [Ultrametric Lattice], [Chaining / Brk.], [Diameter Control], [Dasg. / MW Ratio], [Distinctive Geometric Role],
      [\#56 SLINK @sibson1973slink], [Unique subdom. $u_-$], [Severe ($0$)], [Unbounded], [$Omega(n^(1/3)) "/" O(n^(-1/3))$], [Exact MST / HDBSCAN\* core],
      [\#57 CLINK @defays1977clink], [Min. superdom. $u_+$], [Immune], [Greedy ($>= 3/2$)], [$Omega(sqrt(n)) "/" O(n^(-1/2))$], [Clique-cover compactness],
      [\#58 UPGMA @sokal1958statistical], [Interior ultrametric], [Moderate], [Avg. variance], [$>= 1/3$ MW rev. @moseley2017approximation], [Unique const.-factor tree rev.],
      [\#59 Ward @ward1963hierarchical], [Dilated ultrametric], [Strong], [Spherical ESS], [SSE-greedy (Euclid.)], [Gaussian MLE / balanced spheres],
      [\#69 Genie @gagolewski2016genie], [Non-ultrametric], [Gini-blocked], [MST + Gini], [Tracks SL on $G <= g$], [Non-convex + outlier absorption],
      [\#73 Minimax @bien2011minimax], [Prototype ultrametric], [Immune], [Factor-2 ball cover], [Non-LW ($O(n^3)$)], [Data-anchored $B(p_r, h)$ balls],
    )
  },
)

// =============================================================================
// SECTION III
// =============================================================================
#sec("III", "Scalable Multi-Stage and Metric-Space Hierarchical Architectures")

Every exact pairwise linkage in Sections I--II demands $Omega(n^2)$ distance evaluations, rendering direct agglomeration infeasible once $n > 10^5$. Seven algorithms in Category III break this quadratic wall via multi-stage compression and domain-adapted structural summaries: three height-balanced streaming summary trees—BIRCH (Alg. 64 @zhang1996birch), BICO (Alg. 78 @fichtenberger2013bico), and BUBBLE (Alg. 71 @ganti1999bubble)—and four non-spherical, graph-dynamic, categorical, or density-core hybrids—CURE (Alg. 65 @guha1998cure), CHAMELEON (Alg. 66 @karypis1999chameleon), ROCK (Alg. 67 @guha2000rock), and HCA-DBSCAN (Alg. 70).

#subsec("A", "Streaming Sufficient-Statistic Trees: BIRCH vs. BICO vs. BUBBLE")

BIRCH, BICO, and BUBBLE all scan massive datasets incrementally into a bounded-memory tree of $s << n$ micro-clusters before global Phase 2 clustering, yet diverge sharply in *coreset guarantees* and *ambient geometry*: *BIRCH (\#64, @zhang1996birch)* stores exact additive Euclidean *Clustering Features* $"CF"(C) = (N, bold("LS"), "SS")$ ($bold("LS") = sum x$, $"SS" = sum norm(x)_2^2$) gated by a heuristic uniform radius threshold $r <= T$; *BICO (\#78, @fichtenberger2013bico)* replaces the uniform radius test with a coreset-calibrated multilevel projection tree enforcing a global 1-means residual budget $R <= (epsilon^2 / 9) "OPT"_k (cal(X))$ @harpeled2004coresets; and *BUBBLE (\#71, @ganti1999bubble)* generalizes to *arbitrary metric spaces* $(cal(X), d)$ via *Cluster Representations* $"CR"(C) = (N, "Rep"(C), hat(o)_C, r(C), {"RowSum"_C (o)})$ anchored at the discrete clustroid $hat(o)_C in arg min_(o in C) sum_(x in C) d(o, x)^2$, with cross-row-sum surrogate $hat(S)_2 (o) = |C_2| (d(o, hat(o)_2)^2 + r(C_2)^2)$ and optional FastMap routing (BUBBLE-FM).

#theorem(title: "Strong Coreset Separation (BICO vs. BIRCH) and Metric Clustroid Envelopes (BUBBLE)", origin: "Original Unified Proof")[
  Let $cal(X) = union.sq.big_(j=1)^s C_j$ be compressed into $s$ leaf summaries.
  1. *(Exact Two-Sided Strong Coreset Bound for Weighted CFs, @fichtenberger2013bico):* In $bb(R)^d$, let $cal(S) = {(n_j, mu_j)}_(j=1)^s$ be the weighted leaf centroids with total 1-means compression residual $R := sum_(j=1)^s sum_(x in C_j) norm(x - mu_j)_2^2$. For *every* set of $k$ centers $Q subset.eq bb(R)^d$, writing $"cost"(cal(X), Q) := sum_(x in cal(X)) min_(q in Q) norm(x - q)_2^2$ and $"cost"(cal(S), Q) := sum_(j=1)^s n_j min_(q in Q) norm(mu_j - q)_2^2$,
     $
       "cost"(cal(X), Q) &<= "cost"(cal(S), Q) + R, \
       "cost"(cal(S), Q) &<= "cost"(cal(X), Q) + R + 2 sqrt(R dot "cost"(cal(X), Q)).
     $
     Hence, whenever BICO calibrates its rebuild thresholds so that $R <= (epsilon^2 / 9) "OPT"_k (cal(X))$ for $epsilon in (0, 1]$, $cal(S)$ is an exact *strong $(k, epsilon)$-coreset*: $(1 - epsilon) "cost"(cal(X), Q) <= "cost"(cal(S), Q) <= (1 + epsilon) "cost"(cal(X), Q)$ uniformly for all $Q subset.eq bb(R)^d$ of size $k$. Running weighted $k$-means++ @arthur2007kmeanspp on $cal(S)$ therefore yields an end-to-end $O((1 + epsilon) log k)$-approximation in $O(n d log s + s^2 d)$ time.
  2. *(Unbounded Distortion Failure of Uniform-Threshold BIRCH):* For every distortion target $M > 0$ and any uniform radius threshold $T > 0$ (even with $s >= k$ leaves available), there exists a 1D Euclidean dataset on which BIRCH's uniform-radius absorption merges points into leaf centroids yielding $k$-means cost ratio $("cost"(cal(X), Q_"BIRCH")) / "OPT"_k (cal(X)) > M$, whereas BICO achieves $(1 + epsilon)$ approximation.
  3. *(Universal Metric Clustroid Bound & Surrogate Envelope of BUBBLE):* In any metric space $(cal(X), d)$, the clustroid $hat(o)_C$ and metric radius $r(C) = (N^(-1) sum_(x in C) d(hat(o)_C, x)^2)^(1 / 2)$ satisfy
     $
       r(C) <= 2 inf_(z in cal(X)) (1 / N sum_(x in C) d(z, x)^2)^(1 / 2)
     $
     (improving to $sqrt(2)$ in Euclidean space), and BUBBLE's cross row-sum surrogate $hat(S)_2 (o) = |C_2| (d(o, hat(o)_2)^2 + r(C_2)^2)$ satisfies the unconditional two-sided metric envelope
     $
       1 / 2 S_2 (o) <= hat(S)_2 (o) <= 2 S_2 (o) + 3 |C_2| r(C_2)^2
     $
     for every $o in cal(X)$, becoming *exact* ($hat(S)_2 (o) = S_2 (o)$) in Euclidean space when $hat(o)_2 = mu_(C_2)$.
]

#proof[
  *Part 1:* For any $Q subset.eq bb(R)^d$ with $|Q| = k$, denote $q(x) := arg min_(q in Q) norm(x - q)_2$ and $q_j := arg min_(q in Q) norm(mu_j - q)_2$. Assigning every $x in C_j$ to $q_j$ (which cannot beat $q(x)$) and applying the parallel-axis theorem inside each $C_j$ yields
  $
    "cost"(cal(X), Q) &<= sum_(j=1)^s sum_(x in C_j) norm(x - q_j)_2^2 \
      &= sum_(j=1)^s ( n_j norm(mu_j - q_j)_2^2 + sum_(x in C_j) norm(x - mu_j)_2^2 ) \
      &= "cost"(cal(S), Q) + R.
  $
  Conversely, for every $x in C_j$, optimality of $q_j$ and the triangle inequality give $norm(mu_j - q_j)_2 <= norm(x - q(x))_2 + norm(x - mu_j)_2$. Squaring, summing over $j, x$, and applying Cauchy--Schwarz yields
  $
    "cost"(cal(S), Q) &<= "cost"(cal(X), Q) + R \
      &quad + 2 sum_(j=1)^s sum_(x in C_j) norm(x - q(x))_2 norm(x - mu_j)_2 \
      &<= "cost"(cal(X), Q) + R + 2 sqrt("cost"(cal(X), Q) dot R).
  $
  When $R <= (epsilon^2 / 9) "OPT"_k (cal(X)) <= (epsilon^2 / 9) "cost"(cal(X), Q)$ and $epsilon in (0, 1]$, these two inequalities yield $(1 - epsilon) "cost"(cal(X), Q) <= "cost"(cal(S), Q) <= (1 + epsilon) "cost"(cal(X), Q)$.

  *Part 2:* Fix $k = 2$, $T > 0$, and $M > 1$. Choose an integer $m > 81 M$ and construct the 1D dataset consisting of $m$ points at $-T$, $m$ points at $+T$, and a single isolated point at $10 T$. Placing $k = 2$ centers at $\{-T, +T\}$ incurs optimal cost $"OPT"_2 (cal(X)) <= (9 T)^2 = 81 T^2$. Now stream the $2 m$ points alternating $(-T, +T, -T, +T, dots)$ before $10 T$: at every absorption step, the running leaf centered at $0$ maintains radius $<= T$ (with final centroid $0$ and exact radius $T$), so BIRCH absorbs all $2 m$ points into a *single* leaf $(2 m, 0)$ and places $10 T$ in a second leaf $(1, 10 T)$. Phase 2 clustering on the two surviving leaves places centers at $\{0, 10 T\}$, incurring true cost $2 m T^2$ and distortion ratio $>= (2 m T^2) / (81 T^2) > 2 M$. Because BIRCH gates leaves by *unweighted* radius $r <= T$ instead of *weighted* residual cost $n_j r^2 <= T_"res"$ (as BICO enforces), a single high-mass leaf hides $O(m T^2)$ distortion.

  *Part 3:* Fix any $z in cal(X)$. By the triangle inequality and $(a + b)^2 <= 2(a^2 + b^2)$, $d(x, y)^2 <= 2 d(x, z)^2 + 2 d(z, y)^2$ for all $x, y in C$. Averaging over $C times C$ yields
  $
    r(C)^2 &<= 1 / N^2 sum_(x, y in C) d(x, y)^2 \
      &<= 4 / N sum_(x in C) d(x, z)^2;
  $
  taking square roots and infimum over $z$ gives the factor-$2$ bound ($sqrt(2)$ in Euclidean space, since $N^(-2) sum norm(x - y)_2^2 = 2 sigma^2 (C)$). For BUBBLE's surrogate, summing $d(o, y)^2 <= 2 d(o, hat(o)_2)^2 + 2 d(hat(o)_2, y)^2$ over $y in C_2$ gives
  $
    S_2 (o) &= sum_(y in C_2) d(o, y)^2 \
      &<= 2 |C_2| (d(o, hat(o)_2)^2 + r(C_2)^2) = 2 hat(S)_2 (o),
  $
  proving $(1 / 2) S_2 (o) <= hat(S)_2 (o)$. Symmetrically, summing $d(o, hat(o)_2)^2 <= 2 d(o, y)^2 + 2 d(hat(o)_2, y)^2$ gives $|C_2| d(o, hat(o)_2)^2 <= 2 S_2 (o) + 2 |C_2| r(C_2)^2$, whence $hat(S)_2 (o) <= 2 S_2 (o) + 3 |C_2| r(C_2)^2$. In Euclidean space with $hat(o)_2 = mu_(C_2)$, the parallel-axis theorem makes $hat(S)_2 (o) = S_2 (o)$ exact.
]

#subsec("B", "Non-Spherical, Dynamic-Graph, Categorical, and Density-Core Hybrids")

BIRCH, BICO, and BUBBLE all summarize clusters by single spherical centers (centroids or clustroids), inheriting spherical bias in Phase 2 unless combined with a shape-aware linkage. Four specialized Category III architectures—CURE (Alg. 65 @guha1998cure), CHAMELEON (Alg. 66 @karypis1999chameleon), ROCK (Alg. 67 @guha2000rock), and HCA-DBSCAN (Alg. 70)—tackle non-spherical manifolds, multi-density bridges, discrete categorical transactions, and noisy backgrounds directly.

#theorem(title: "Geometry, Scale Invariance, Link-Kernel Equivalence, and Bridge Immunity Across CURE, CHAMELEON, ROCK, and HCA-DBSCAN", origin: "Original Unified Proof")[
  1. *(CURE Shrinkage Interpolation & Concentric-Manifold Collapse, Alg. 65):* Representing each cluster $C$ by $c$ well-scattered points shrunk toward $mu_C$ via $r'_j = (1 - alpha) r_j + alpha mu_C$ ($alpha in [0, 1]$) recovers exact Single Linkage at $(c = |C|, alpha = 0)$ and Centroid Linkage at $alpha = 1$ (or $c = 1$). Shrinking suppresses radial outlier chaining by shifting an outlier at distance $R$ from $mu_C$ inward by $alpha R$; *however*, on two concentric spheres $S_1, S_2 subset.eq bb(R)^d$ of radii $R_1 < R_2$ centered at the origin, shrinking $S_2$'s representatives by any $alpha >= 1 - R_1 / R_2$ moves them strictly inside or onto $S_1$, causing CURE to falsely collapse $S_1$ and $S_2$ regardless of how large $c$ is.
  2. *(CHAMELEON Local Homothety Invariance & Bridge Rejection, Alg. 66):* In CHAMELEON's Phase 2, let $cal(G)_k$ be the weighted $k$-NN graph ($w_(i j) = 1 / d(x_i, x_j)$) partitioned by Phase 1 METIS min-cut bisections @karypis1998metis into subclusters, merged by maximizing $"RI"(C_i, C_j) dot "RC"(C_i, C_j)^alpha$ with $w_i := |C_i| / (|C_i| + |C_j|)$:
     $
       "RI"(C_i, C_j) &= (2 |"EC"(C_i, C_j)|) / (|"EC"(C_i)| + |"EC"(C_j)|), \
       "RC"(C_i, C_j) &= (overline("EC")(C_i, C_j)) / (w_i overline("EC")(C_i) + w_j overline("EC")(C_j)).
     $
     On separated clusters scaled independently by arbitrary positive factors $lambda_r > 0$ preserving intra-cluster $k$-NN lists, *both* $"RI"$ and $"RC"$ are identically invariant to $(lambda_1, dots, lambda_K)$. Moreover, a thin bridge of $b << min{|"EC"(C_i)|, |"EC"(C_j)|}$ edges between two dense clusters yields $"RI" = O(b / |"EC"|) -> 0$, rejecting chaining without centroid shrinkage (at the cost of non-monotonicity and non-ultrametric merge scores).
  3. *(ROCK Positive-Semidefinite Link Kernel & Exact Phase Transition at $f(theta) = 1 / 2$, Alg. 67):* Let $A in {0, 1}^(n times n)$ be the symmetric neighbor adjacency matrix ($A_(i j) = bold(1)["sim"(x_i, x_j) >= theta]$), so the shared-neighbor link matrix is $L = A^2 succ.eq 0$. ROCK's goodness measure $g(C_i, C_j) = "link"(C_i, C_j) / ((n_i + n_j)^(1 + 2 f(theta)) - n_i^(1 + 2 f(theta)) - n_j^(1 + 2 f(theta)))$ is a monotone greedy ascent step on $E_ell = sum_r (bold(1)_(C_r)^top A^2 bold(1)_(C_r)) / (n_r^(2 f(theta)))$ @guha2000rock. At $f(theta) = 1 / 2$, maximizing $E_ell$ is *algebraically identical* to minimizing Kernel $k$-Means @dhillon2004kernel with PSD link kernel $K = A^2$:
     $
       &sum_(r=1)^k sum_(x in C_r) norm(phi(x) - mu_r^phi)_2^2 \
       &quad = "Tr"(A^2) - sum_(r=1)^k (bold(1)_(C_r)^top A^2 bold(1)_(C_r)) / (|C_r|).
     $
     For $f(theta) > 1 / 2$, the denominator grows super-quadratically, penalizing giant clusters; for $f(theta) < 1 / 2$, it grows sub-quadratically, inducing snowball cluster absorption.
  4. *(HCA-DBSCAN Truncated Core Equivalence & Unbounded Noise-Bridge Immunity, Alg. 70):* Extracting DBSCAN cores $cal(M) = {M_1, dots, M_m}$ at radius $epsilon$ and $"MinPts" = mu$ before Stage 2 agglomeration annihilates every sparse noise bridge whose points satisfy $|N_epsilon (x)| < mu$ (achieving unbounded bridge immunity vs. SLINK's $0$-breakdown chaining). When Stage 2 executes Single Linkage over $cal(M)$, the resulting dendrogram is isomorphic to $epsilon$-truncated Single Linkage on $cal(X)_"core"$; when Stage 2 executes weighted Ward linkage over core summaries $(|M_a|, overline(x)_(M_a))$, it matches exact point-level Huygens--Steiner ESS merges across cores in $O(m^2)$ time ($m << n$) while resolving DBSCAN's border-point order ambiguity via deterministic nearest-core projection.
]

#proof[
  *Part 1:* At $(c = |C|, alpha = 0)$, $r'_j = r_j$ ranges over all points of $C$, so $D_"CURE" (A, B) = min_(a in A, b in B) norm(a - b)_2 = D_"SL" (A, B)$. At $alpha = 1$, every representative collapses to $r'_j = mu_C$, yielding $D_"CURE" (A, B) = norm(mu_A - mu_B)_2 = D_"UPGMC" (A, B)$. For an outlier $z in C$ at displacement $z - mu_C$ with norm $R$, shrinking yields $norm(z' - mu_C)_2 = (1 - alpha) R$, increasing its distance to any point $y$ outside the ball $B(mu_C, R)$ along the outward ray by exactly $alpha R$. On two symmetric concentric spheres $S_1 = R_1 bb(S)^(d-1)$ and $S_2 = R_2 bb(S)^(d-1)$ ($R_1 < R_2$) densely sampled so $mu_(S_2) = bold(0)$, every representative $r_j in S_2$ has shrunk radius $norm(r'_j)_2 = (1 - alpha) R_2$. Whenever $alpha >= 1 - R_1 / R_2$, we obtain $(1 - alpha) R_2 <= R_1$, pushing $S_2$'s representatives directly onto or inside the inner sphere $S_1$ so that cross-cluster representative distance collapses below intra-sphere spacing.

  *Part 2:* Under positive cluster-wise scaling $tilde(x) = c_r + lambda_r (x - c_r)$ preserving each cluster's $k$-NN topology, every intra-cluster edge weight $w_(u v) = 1 / d(u, v)$ inside cluster $r$ scales uniformly by $lambda_r^(-1)$. Consequently, for any two subclusters $C_i, C_j$ belonging to true cluster $r$, the edge-cut sums $|"EC"(C_i, C_j)|$, $|"EC"(C_i)|$, and $|"EC"(C_j)|$ all scale by the identical factor $lambda_r^(-1)$, as do their mean edge weights $overline("EC")$. Every factor $lambda_r^(-1)$ cancels identically between numerator and denominator of both $"RI"(C_i, C_j)$ and $"RC"(C_i, C_j)$, proving exact multi-scale homothety invariance. For two dense cliques of size $m$ bridged by $b$ edges of equal weight, $|"EC"(C_i)| = Theta(m^2)$ while $|"EC"(C_i, C_j)| = b$, yielding $"RI" = O(b / m^2) -> 0$ regardless of bridge proximity, whereas Single Linkage merges along the single shortest bridge edge. Non-monotonicity follows because merging two tightly knit halves increases internal bisection weight $|"EC"(C_i union C_j)|$, which can either raise or lower $"RI"$ relative to subsequent merges.

  *Part 3:* Since $A$ is real symmetric, $L = A^2 = A^top A$ satisfies $v^top L v = norm(A v)_2^2 >= 0$ for all $v in bb(R)^n$, proving $L succ.eq 0$ with exact Gram feature map $phi(x_i) = A_(i, :)$ (the $i$-th neighbor indicator vector) and $L_(i j) = chevron.l A_(i, :), A_(j, :) chevron.r = |N_theta (x_i) inter N_theta (x_j)|$. At $f(theta) = 1 / 2$, $2 f(theta) = 1$, so ROCK's objective becomes $E_ell = sum_(r=1)^k (bold(1)_(C_r)^top A^2 bold(1)_(C_r)) / |C_r|$. Expanding Kernel $k$-Means distortion in feature space $phi$ with centroid $mu_r^phi = |C_r|^(-1) sum_(x in C_r) phi(x)$ gives
  $
    &sum_(r=1)^k sum_(i in C_r) norm(phi(x_i) - mu_r^phi)_2^2 \
      &quad = sum_(i=1)^n norm(phi(x_i))_2^2 - sum_(r=1)^k |C_r| norm(mu_r^phi)_2^2 \
      &quad = "Tr"(A^2) - sum_(r=1)^k (bold(1)_(C_r)^top A^2 bold(1)_(C_r)) / (|C_r|) = "Tr"(A^2) - E_ell,
  $
  establishing exact algebraic equivalence! For general $f(theta)$, observe that merging $C_i$ and $C_j$ replaces $n_i^(1 + 2 f) + n_j^(1 + 2 f)$ by $(n_i + n_j)^(1 + 2 f)$ in the expected intra-cluster link denominator: because $n |-> n^(1 + 2 f)$ is strictly convex iff $1 + 2 f > 1$ (and grows faster than $n^2$, the max possible pair count, iff $f > 1 / 2$), $f(theta) > 1 / 2$ increasingly penalizes large merges while $f(theta) < 1 / 2$ rewards absorbing nodes into already large clusters.

  *Part 4:* Every bridge point $z$ with $|N_epsilon (z)| < mu$ is excluded from $cal(X)_"core"$, so removing arbitrarily many sparse bridge chains leaves the core graph $G_epsilon = (cal(X)_"core", E_epsilon)$ and all core components $cal(M)$ completely untouched. Because $M_1, dots, M_m$ are exactly the connected components of the $epsilon$-threshold graph on $cal(X)_"core"$, contracting every edge of weight $<= epsilon$ in the MST of $cal(X)_"core"$ produces the exact MST over $(cal(M), D_"SL")$, proving dendrogram isomorphism to $epsilon$-truncated Single Linkage on $cal(X)_"core"$. For Stage 2 Ward linkage, Theorem 2(1) shows that $Delta "ESS"(M_a, M_b)$ depends exclusively on $(|M_a|, overline(x)_(M_a))$ and $(|M_b|, overline(x)_(M_b))$, so running Ward NN-Chain on $cal(M)$ computes exact point-level ESS increases across all macro-merges in $O(m^2)$ time; assigning each border point $y in cal(X)_"border"$ to $arg min_(x in cal(X)_"core") d(y, x)$ (breaking ties by deterministic index) removes DBSCAN's BFS visitation-order dependence completely.
]

#wide_figure(
  caption: [*Table III:* Head-to-Head Comparison of Scalable Multi-Stage and Metric-Space Hierarchical Architectures (Algs. 64--67, 70--71, 78).],
  {
    set text(size: 7.3pt)
    set table(inset: (x: 3.5pt, y: 2.7pt))
    table(
      columns: (1.05fr, 1.15fr, 0.92fr, 0.72fr, 0.78fr, 1.15fr),
      align: (left, left, center, center, center, left),
      [Algorithm], [Compression / Summary], [Time Complexity], [Space / Geometry], [Formal Guarantee], [Primary Failure Mode],
      [\#64 BIRCH @zhang1996birch], [Euclidean CF $(N, bold("LS"), "SS")$], [$O(n d B log_B m)$], [Euclidean], [Exact CF additivity], [Unbounded $k$-means distortion; spheres],
      [\#78 BICO @fichtenberger2013bico], [Coreset-calibrated CF tree], [$O(n d log s + s^2 d)$], [Euclidean], [Strong $(k, epsilon)$-coreset], [Spherical $k$-means geometry only],
      [\#71 BUBBLE @ganti1999bubble], [Metric CR + clustroid / FastMap], [$O(n(B log_B n_L + s) + n_L^2)$], [Any metric], [Factor-2 clustroid bound], [Expensive metric evals; FastMap dist.],
      [\#65 CURE @guha1998cure], [$c$ shrunk reps + random sample], [$O(n_s^2 log n_s + n k c d)$], [Euclidean], [Chernoff sample + margin], [Concentric manifold collapse (Thm. 7.1)],
      [\#66 CHAMELEON @karypis1999chameleon], [$k$-NN METIS bisections + RI/RC], [$O(n^2 d + n m + m^2 log m)$], [Metric / $k$-NN], [Multi-scale homothety inv.], [Non-monotone; high $d$ $k$-NN hubness],
      [\#67 ROCK @guha2000rock], [Shared-neighbor link matrix $A^2$], [$O(n^2 + n m_m m_a + n^2 log n)$], [Categorical / set], [Kernel $k$-Means at $f = 1 / 2$], [Sensitive to $theta$; $O(n^2)$ memory],
      [\#70 HCA-DBSCAN], [DBSCAN cores $cal(M)$ + LW tree], [$O(n log n + m^2)$ / $O(n^2)$], [Metric / Euclid.], [Unbounded bridge immunity], [Inherits Stage 1 single-scale $epsilon$],
    )
  },
)

// =============================================================================
// SECTION IV
// =============================================================================
#sec("IV", "Divisive, Phylogenetic, and Extreme-Scale Incremental Hierarchies")

The remaining seven algorithms of Category III depart from classical bottom-up batch agglomeration along three distinct frontiers: (A) *top-down divisive splitting*—polythetic DIANA (Alg. 63 @kaufman1990finding) versus monothetic binary MONA (Alg. 68 @kaufman1990finding); (B) *non-ultrametric additive phylogeny*—Neighbor-Joining (Alg. 75 @saitou1987neighbor) contrasted against molecular-clock UPGMA/WPGMA; and (C) *web-scale incremental and parallel round-based hierarchies*—PERCH (Alg. 76 @kobren2017perch), GRINCH (Alg. 77 @monath2019grinch), SCC (Alg. 79 @monath2021scc), and FINCH (Alg. 80 @sarfraz2019finch).

#subsec("A", "Top-Down Divisive Trees: Polythetic DIANA vs. Monothetic MONA")

Whereas agglomerative methods make irreversible greedy decisions at the *finest* scale (pairwise distances) where noise dominates, divisive methods examine *global* cluster structure at the root first. DIANA (Alg. 63) splits the active cluster $C$ of largest diameter $"diam"(C)$ by seeding a splinter group $S = {x^*}$ at the point of largest mean dissimilarity in $C$ and iteratively relocating any $i in R = C without S$ maximizing $D_i = overline(d)(i, R without {i}) - overline(d)(i, S) > 0$. MONA (Alg. 68) restricts to binary feature matrices $X in {0, 1}^(n times d)$ and splits $C$ along the single coordinate $f^* in arg max_f A_C (f) = arg max_f sum_(g != f) |a_(f g) d_(f g) - b_(f g) c_(f g)|$, producing exact rule-interpretable subcubes.

#theorem(title: "DIANA Ultrametric Failure vs. MONA Covariance Equivalence and Parity Blindness", origin: "Original Unified Proof")[
  1. *(DIANA Dominating Ultrametric, Silhouette Sign, and $n = 6$ Exact-Ultrametric Failure, Alg. 63):* Because $C' subset.eq C => "diam"(C') <= "diam"(C)$, DIANA's split heights are monotonically non-increasing along every root-to-leaf path, inducing a *superdominant ultrametric* $u_"DIANA" >= d$ (in the same upper cone $cal(U)_+(d)$ as Complete Linkage, at worst-case time $O(n^3)$ vs. CLINK's $O(n^2)$). Moreover, for every candidate $i in R$, DIANA's transfer gain $D_i = overline(d)(i, R without {i}) - overline(d)(i, S)$ has the exact opposite sign of Rousseeuw's 2-cluster Silhouette @rousseeuw1987silhouettes ($D_i > 0 <=> s(i) < 0$). *Nevertheless*, there exists an exact $n = 6$ ultrametric $(cal(X), u)$ whose true root bipartition is recovered identically by SLINK, CLINK, UPGMA, WPGMA, and Minimax Linkage, yet is *shattered on the very first split* by DIANA.
  2. *(MONA Exact Covariance Identity, Outlier Damping, and Second-Order Parity Blindness, Alg. 68):* For any cluster $C$ of cardinality $m = |C|$ and binary attributes $f, g in {1, dots, d}$ with prevalences $p_f, p_g$, MONA's pairwise determinant association satisfies the exact identity
     $
       alpha_C (f, g) = |a_(f g) d_(f g) - b_(f g) c_(f g)| = m^2 |"Cov"_C (X_f, X_g)|,
     $
     which automatically damps rare-attribute outlier splits by $sqrt(p_f (1 - p_f))$ relative to Pearson $chi^2$ association analysis (and eliminates $0 / 0$ singularities). *However*, MONA is provably blind to higher-order interactions: on the $d = 3$ XOR parity dataset $cal(X)_"XOR"$ consisting of $r >= 1$ copies of each even-parity vertex ${(0, 0, 0), (0, 1, 1), (1, 0, 1), (1, 1, 0)}$ ($n = 4 r$), every pairwise association vanishes identically ($A_(cal(X)_"XOR") (f) = 0$ for all $f$), forcing MONA to terminate at the root with zero splits despite all four clusters being separated by Hamming distance $2$, whereas polythetic DIANA and Lance--Williams methods recover all four clusters cleanly.
]

#proof[
  *Part 1:* Monotonicity $C' subset.eq C => "diam"(C') <= "diam"(C)$ guarantees that assigning split height $h(C) = "diam"(C)$ yields a valid ultrametric $u_"DIANA" (x, y) = "diam"(C_(x or y)) >= d(x, y)$. For any $i in R$, $a(i) = overline(d)(i, R without {i})$ and $b(i) = overline(d)(i, S)$, so Rousseeuw's Silhouette $s(i) = (b(i) - a(i)) / max{a(i), b(i)} = -D_i / max{a(i), b(i)}$ is strictly negative iff $D_i > 0$.

  Now construct the exact $n = 6$ ultrametric $u$ on $cal(X) = {a_1} union {a_2, a_3, a_4} union {b_1, b_2}$ generated by the rooted tree: root splits at height $1.0$ into $A = {a_1, a_2, a_3, a_4}$ and $B = {b_1, b_2}$; $A$ splits at height $0.95$ into ${a_1}$ and $A' = {a_2, a_3, a_4}$; both $A'$ and $B$ merge internally at height $0.1$. Every pairwise distance satisfies the strong triangle inequality identically (so SLINK, CLINK, UPGMA, WPGMA, and Minimax all recover the exact root split ${A, B}$). Running DIANA's root split on $cal(X)$:
  - Mean dissimilarity of each element to the other $5$:
    $
      overline(d)(a_1) &= (3 dot 0.95 + 2 dot 1.0) / 5 = 0.97, \
      overline(d)(b_r) &= (1 dot 0.1 + 4 dot 1.0) / 5 = 0.82, \
      overline(d)(a_j) &= (2 dot 0.1 + 0.95 + 2 dot 1.0) / 5 = 0.63 quad (j >= 2).
    $
  - DIANA uniquely selects $a_1$ as splinter seed: $S = {a_1}$, $R = {a_2, a_3, a_4, b_1, b_2}$.
  - Relocating test $D_x = overline(d)(x, R without {x}) - overline(d)(x, S)$ for every $x in R$:
    - For $b_r in B$: $D_(b_r) = (3 dot 1.0 + 0.1) / 4 - 1.0 = -0.225 < 0$.
    - For $a_j in A'$: $D_(a_j) = (2 dot 0.1 + 2 dot 1.0) / 4 - 0.95 = -0.400 < 0$.
  Because every $D_x < 0$, DIANA halts immediately with the erroneous root split ${{a_1}, {a_2, a_3, a_4, b_1, b_2}}$, breaking the exact ultrametric clade $A$ permanently.

  *Part 2:* Substitute $b_(f g) = s_f - a_(f g)$, $c_(f g) = s_g - a_(f g)$, and $d_(f g) = m - s_f - s_g + a_(f g)$ into $a_(f g) d_(f g) - b_(f g) c_(f g)$:
  $
    a_(f g) d_(f g) - b_(f g) c_(f g) &= a_(f g) (m - s_f - s_g + a_(f g)) \
      &quad - (s_f - a_(f g))(s_g - a_(f g)) \
      &= m a_(f g) - s_f s_g \
      &= m^2 (a_(f g) / m - p_f p_g) = m^2 "Cov"_C (X_f, X_g).
  $
  By Cauchy--Schwarz, $|"Cov"_C (X_f, X_g)| <= sqrt(p_f (1 - p_f) p_g (1 - p_g))$, which scales as $O(sqrt(p_f))$ as $p_f -> 0^+$, damping single-object outlier splits; by contrast, $chi^2_(f g) = m "Cov"_C (X_f, X_g)^2 / (p_f (1 - p_f) p_g (1 - p_g))$ divides out $p_f (1 - p_f)$, allowing rare splits to dominate and blowing up when $p_f in {0, 1}$. Finally, on $cal(X)_"XOR"$, for every distinct pair $(f, g) in binom({1, 2, 3}, 2)$, the four bit pairs $(0, 0), (0, 1), (1, 0), (1, 1)$ each occur exactly $r$ times, so $a_(f g) = b_(f g) = c_(f g) = d_(f g) = r$, giving $alpha(f, g) = |r^2 - r^2| = 0$ and $A(f) = 0$ for every coordinate $f$, stalling MONA at the root. Under Hamming distance, every distinct pair of XOR vertices lies at distance $2$ while duplicates lie at distance $0$, so DIANA and every agglomerative linkage isolate all four clusters with zero error.
]

#subsec("B", "Ultrametric vs. Additive Phylogeny: Neighbor-Joining vs. UPGMA / WPGMA")

In computational phylogenetics and evolutionary taxonomy, molecular evolutionary rates frequently vary across lineages, violating the *molecular clock* (ultrametricity) while still generating an exact *additive tree metric* $d(i, j) = sum_(e in P_(i j)) ell_e$ along an unrooted tree $T = (V, E, ell)$ with positive edge lengths $ell_e > 0$. UPGMA (Alg. 58) and WPGMA (Alg. 60) greedily merge the pair minimizing raw $d(i, j)$, implicitly assuming equal root-to-leaf path lengths; *Neighbor-Joining* (NJ, Alg. 75 @saitou1987neighbor) corrects for unequal branch lengths by minimizing the rate-adjusted $Q$-criterion over $r >= 4$ active taxa:
$
  Q(i, j) := (r - 2) d(i, j) - S_i - S_j, quad S_i := sum_(k=1)^r d(i, k),
$
contracting $(i, j)$ to a parent node $u$ with distances $d(u, k) = (1 / 2)(d(i, k) + d(j, k) - d(i, j))$ (identical to a Lance--Williams update with $alpha_i = alpha_j = 1 / 2, beta = -1 / 2, gamma = 0$, yet coupled with a non-local global selection criterion $Q$).

#theorem(title: "Exact Additive Consistency of Neighbor-Joining vs. Molecular-Clock Collapse of UPGMA and Euclidean Distortion", origin: "Original Unified Proof")[
  Let $T = (V, E, ell)$ be an unrooted binary phylogenetic tree on $r >= 4$ leaves with strictly positive edge lengths $ell_e > 0$ and induced additive metric $d$.
  1. *(Exact Edge-Weight Decomposition of $Q$ & Global Consistency of NJ, @saitou1987neighbor):* Let $L(T) = sum_(e in E) ell_e$ be the total tree length and $E_"int" (T)$ the set of internal edges. For every distinct leaf pair $(i, j)$,
     $
       Q(i, j) &= -2 L(T) \
         &quad - 2 sum_(e in E_"int" (T) without P_(i j)) (r - m_e (i, j) - 1) ell_e,
     $
     where $m_e (i, j) in {2, dots, r - 2}$ is the leaf cardinality of the component of $T without {e}$ containing both $i$ and $j$. Consequently, *every* pendant branch length cancels out identically across all pairs, every minimizer $(x, y) in arg min_(i != j) Q(i, j)$ is strictly a *cherry* (adjacent pair of leaves) of $T$, and NJ reconstructs the exact topology and every branch length of $T$ in $O(n^3)$ time (retaining exact topology whenever additive noise satisfies $norm(hat(d) - d)_infinity < (1 / 2) min_(e in E_"int") ell_e$).
  2. *(Long-Branch Attraction Failure of UPGMA / WPGMA / SLINK / CLINK):* On the 4-taxon additive quartet $((1, 2) | (3, 4))$ with pendant lengths $(ell_1, ell_2, ell_3, ell_4) = (5, 1, 5, 1)$ and internal edge $ell_e = 1$, *all four* ultrametric linkages (SLINK, CLINK, UPGMA, WPGMA) falsely merge the short-branch non-cherry $(2, 4)$ ($d(2, 4) = 3 < 6 = d(1, 2) = d(3, 4)$), whereas NJ yields $Q(1, 2) = Q(3, 4) = -28 < -26 = Q(2, 4)$ and recovers $((1, 2) | (3, 4))$ identically.
  3. *(Euclidean Clustering Misfit of NJ):* Conversely, on non-tree Euclidean point clouds, NJ's subtraction of $-S_i - S_j$ systematically favors merging two *distant extreme outliers* (whose large row sums $S_i, S_j$ depress $Q(i, j)$) ahead of dense interior clusters, producing negative Silhouette partitions where Ward and UPGMA succeed cleanly.
]

#proof[
  *Part 1:* Fix any distinct leaves $i, j$ and examine the coefficient $c_e (i, j)$ of edge length $ell_e$ inside $Q(i, j) = (r - 2) d(i, j) - sum_k d(i, k) - sum_k d(j, k)$. Removing $e$ splits the $r$ leaves into two sides $A_e, B_e$.
  - If $e in P_(i j)$ (say $i in A_e$ of size $m$, $j in B_e$ of size $r - m$): $e$ contributes $+(r - 2)$ via $(r - 2) d(i, j)$, $-(r - m)$ via $-S_i$ (paths from $i$ to $B_e$), and $-m$ via $-S_j$ (paths from $j$ to $A_e$), giving exact coefficient $(r - 2) - (r - m) - m = -2$.
  - If $e in.not P_(i j)$ (so $i, j$ both lie in $A_e$ of size $m_e (i, j)$): $e$ contributes $0$ to $d(i, j)$, $-(r - m_e (i, j))$ via $-S_i$, and $-(r - m_e (i, j))$ via $-S_j$, giving exact coefficient $-2(r - m_e (i, j)) = -2 - 2(r - m_e (i, j) - 1)$.
  Notice that for every pendant edge $e_k$ ($k in.not {i, j}$), $i$ and $j$ lie on the side of size $m_(e_k) (i, j) = r - 1$, so $r - m_(e_k) - 1 = 0$ and *every* pendant edge has constant coefficient $-2$ regardless of $(i, j)$! Summing over all $e in E(T)$ yields the exact identity.

  Now suppose $(i, j)$ is *not* a cherry of $T$. Then the unique path $P_(i j)$ traverses at least one internal edge; let $e_1$ be the first internal edge leaving $i$'s neighbor $v_i$ toward $j$, and $e_"last"$ the last internal edge entering $j$'s neighbor $v_j$. Let $C_i$ be the component of $T without {e_1}$ containing $i$, and $C_j$ the component of $T without {e_"last"}$ containing $j$. Because $v_i, v_j$ have degree $>= 3$, $|C_i| >= 2$ and $|C_j| >= 2$; moreover $C_i inter C_j = emptyset$, so $min{|C_i|, |C_j|} <= r / 2$. Without loss of generality assume $|C_i| <= r / 2$. Because $C_i$ is a finite tree with $>= 2$ leaves, rooting $C_i$ at $v_i$ and picking two leaves descending from a deepest internal node of $C_i$ produces a cherry $(x, y)$ of $T$ strictly inside $C_i$. Compare $Q(x, y)$ with $Q(i, j)$ edge by edge using our exact decomposition (noting $P_(x y)$ contains *zero* internal edges):
  - For every internal edge $e in.not C_i$ with $e in.not P_(i j)$, all of $C_i union {j}$ (hence $x, y, i, j$) sit on the same side of $e$, so $m_e (x, y) = m_e (i, j)$.
  - For any internal edge $e in C_i$ (necessarily $e in.not P_(i j)$): if $e$ does not separate $\{x, y\}$ from $i$, then $m_e (x, y) = m_e (i, j)$; if $e$ separates $\{x, y\}$ from $i$, then the side containing $\{x, y\}$ is strictly contained inside $C_i$, so $m_e (x, y) <= |C_i| - 1 < r / 2 < r - m_e (x, y) = m_e (i, j)$, whence $r - m_e (x, y) - 1 > r - m_e (i, j) - 1$.
  - Finally, $e_1 in P_(i j)$ contributes $0$ to the bonus of $(i, j)$ yet is absent from $P_(x y)$ with $m_(e_1) (x, y) = |C_i| <= r / 2 <= r - 2$ (since $r >= 4$), contributing at least $2(r - |C_i| - 1) ell_(e_1) >= 2 ell_(e_1) > 0$ to $-Q(x, y)$!
  Hence $Q(x, y) <= Q(i, j) - 2 ell_(e_1) < Q(i, j)$ strictly, proving every global minimizer of $Q$ is a cherry. Contracting cherry $(x, y)$ to its parent $u$ sets $d(u, k) = (1 / 2)((ell_x + d(v, k)) + (ell_y + d(v, k)) - (ell_x + ell_y)) = d(v, k)$, preserving exact additivity inductively down to $r = 3$.

  *Part 2:* On $((1, 2) | (3, 4))$ with $(ell_1, ell_2, ell_3, ell_4, ell_e) = (5, 1, 5, 1, 1)$, pairwise distances are $d(1, 2) = d(3, 4) = 6$, $d(2, 4) = 3$, $d(1, 4) = d(2, 3) = 7$, $d(1, 3) = 11$. Because $d(2, 4) = 3$ is the unique global minimum distance, SLINK, CLINK, UPGMA, and WPGMA all merge $(2, 4)$ on step 1, corrupting the quartet split. By Part 1, $L(T) = 5 + 1 + 5 + 1 + 1 = 13$, giving $Q(1, 2) = Q(3, 4) = -2(13) - 2(1) = -28 < -26 = -2(13) = Q(2, 4)$, recovering the exact split.

  *Part 3:* Place six planar points in $bb(R)^2$: four distant boundary vertices $x_1 = (L, 0)$, $x_2 = (0, L)$, $x_3 = (-L, 0)$, $x_4 = (0, -L)$ plus two tight interior points $u = (0, 0)$, $v = (epsilon, 0)$ with $L = 10$ and $epsilon = 0.1$ ($r = 6$, so $r - 2 = 4$). Direct Euclidean distances give $d(u, v) = epsilon$, $d(x_2, x_3) = sqrt(2) L$, and row sums:
  $
    S_u &= 4 L + epsilon, \
    S_v &= 2 L + 2 sqrt(L^2 + epsilon^2) + epsilon \
      &= 4 L + epsilon + epsilon^2 / L + O(epsilon^4 / L^3), \
    S_2 &= (4 + 2 sqrt(2)) L + epsilon^2 / (2 L) + O(epsilon^4 / L^3), \
    S_3 &= (4 + 2 sqrt(2)) L + epsilon.
  $
  Evaluating $Q(i, j) = 4 d(i, j) - S_i - S_j$ yields
  $
    Q(u, v) &= -8 L + 2 epsilon - epsilon^2 / L + O(epsilon^4 / L^3) > -8 L + 1.9 epsilon, \
    Q(x_2, x_3) &= -8 L - epsilon - epsilon^2 / (2 L) + O(epsilon^4 / L^3) < -8 L - epsilon.
  $
  Hence $Q(x_2, x_3) < Q(u, v)$ strictly: NJ merges the two distant boundary points $x_2 = (0, 10)$ and $x_3 = (-10, 0)$ at Euclidean distance $10 sqrt(2) approx 14.14$ (which lie closer to $\{u, v\}$ at distance $10$, yielding negative Silhouette $(10 - 14.14) / 14.14 approx -0.293 < 0$) *before* merging the interior pair $\{u, v\}$ separated by $0.1$.
]

#subsec("C", "Extreme-Scale Incremental & Parallel Hierarchies: PERCH vs. GRINCH vs. SCC vs. FINCH")

Modern extreme-clustering pipelines ($N >= 10^6$, $K >= 10^4$) shatter $O(n^2)$ linkage barriers via two paradigms: (1) *online incremental restructuring*—*PERCH (\#76, @kobren2017perch)* maintains axis-aligned bounding boxes (AABBs) and repairs early misplacements via local sibling--aunt rotations ($d_"max" (v, a) < d_"max" (v, s)$), whereas *GRINCH (\#77, @monath2019grinch)* tracks $O(d)$ bilinear average-linkage summaries and adds *non-local grafting* (detaching subtree $v$ to its globally nearest non-ancestor $v'$ when $"sim"(v, v') > max{"sim"(v, "sib"(v)), "sim"(v', "sib"(v'))}$); and (2) *parallel round-based contraction*—parameter-free *FINCH (\#80, @sarfraz2019finch)* contracts all connected components of the first-neighbor graph $A(i, j) = bold(1)[j = kappa_i^1 or i = kappa_j^1 or kappa_i^1 = kappa_j^1]$ per round, whereas *SCC (\#79, @monath2021scc)* contracts reciprocal nearest-neighbor (RNN) 2-cycles under a monotone similarity schedule $tau_1 > dots > tau_R$, bridging FINCH and batch UPGMA.

#theorem(title: "Local-Rotation Trap Separation (GRINCH vs. PERCH) and Parallel Contraction Bounds (FINCH vs. SCC)", origin: "Original Unified Proof")[
  1. *(Strict Local-Rotation Trap Resolved by a Single GRINCH Graft, Algs. 76--77):* There exists a symmetric 2D Euclidean dataset and binary tree $cal(T)$ separating two true twin points $a_1 = (0, 0), a_2 = (0.1, 0)$ (cluster $C_1$) across opposite root subtrees alongside four pure background clusters $B_1, D_1, B_2, D_2$ of size $m = 10$ each such that:
     - *Every* local sibling rotation at *every* depth is strictly rejected under both PERCH's AABB max-linkage and GRINCH's average-linkage criteria, trapping PERCH permanently at $C_1$-purity $approx 0.524$;
     - A *single non-local graft* in GRINCH detaches $a_1$ and pairs it with $a_2$, immediately restoring exact dendrogram purity $1.0$ across all five clusters.
     Moreover, in $bb(R)^d$, PERCH's AABB bounding boxes suffer axis-rotation masking (expanding volume by up to $d^(d / 2)$ for diagonally oriented manifolds), whereas GRINCH's bilinear average-linkage summary is strictly rotation-invariant ($O(d)$ orthogonal equivariance).
  2. *(Deterministic Logarithmic Depth & Asynchronous Chaining of FINCH, Alg. 80):* In every round $t$ of FINCH on $K_t >= 2$ clusters (with unique 1-NNs), every connected component of $A$ has cardinality $>= 2$, forcing deterministic geometric compression $K_(t+1) <= floor(K_t / 2)$ and exact round ceiling $R <= floor(log_2 n)$ (tight on dyadic nested pairs). However, because FINCH contracts *one-way* 1-NN chains $i -> j -> ell$ into a single component, an asymmetric-density tail pointing sequentially into a dense core triggers *asynchronous chaining* that merges distinct clusters in a single coarse round.
  3. *(Unique 2-Cycle Pseudotrees, Order-Invariance, and Exact Reductions of SCC, Alg. 79):* Under deterministic tie-breaking, every weakly connected component of the directed 1-NN graph on $K$ active clusters is a functional pseudotree containing *exactly one* reciprocal nearest-neighbor 2-cycle $\{u, v\}$ ($kappa_u^1 = v$ and $kappa_v^1 = u$). Contracting strictly reciprocal pairs (or $tau_r$-gated components) in parallel is *identically order-invariant* under data permutations, provably recovers separated clusters when intra-cluster similarities strictly exceed inter-cluster similarities across thresholds, reduces to exact sequential UPGMA as the threshold schedule refines ($Delta tau -> 0^+$ with average linkage), and recovers FINCH's unrestricted 1-NN contraction when $tau = -infinity$ with one-way arcs admitted.
]

#proof[
  *Part 1:* Place $a_1 = (0, 0)$, $a_2 = (0.1, 0)$ (cluster $C_1$, distance $0.1$) and four clusters of $m = 10$ coincident points each: $B_1 = (10, 0)$, $D_1 = (10, 12)$, $B_2 = (-10, 0)$, $D_2 = (-10, 12)$, arranged in the 42-leaf tree $cal(T) = (p_1, p_2)$ with $p_1 = ((a_1, B_1), D_1)$ and $p_2 = ((a_2, B_2), D_2)$ (writing $u_1 = (a_1, B_1)$). Along $a_1$'s ancestry:
  - At $v = a_1$ (sibling $B_1$, aunt $D_1$): $d(a_1, B_1) = 10 < sqrt(244) approx 15.62 = d(a_1, D_1)$, blocking rotation under both max and average linkage.
  - At $v = u_1$ (sibling $D_1$, aunt $p_2$): $d_"avg" (u_1, D_1) = (10 dot 12 + sqrt(244)) / 11 approx 12.33$ ($d_"max" approx 15.62$), whereas the $200$ pairs in $B_1 times (B_2 union D_2)$ out of $231$ total pairs between $u_1$ and $p_2$ (distances $20$ and $sqrt(544) approx 23.32$) enforce $d_"avg" (u_1, p_2) > (2000 + 2332) / 231 approx 18.75 > 12.33$ and $d_"max" = sqrt(544) approx 23.32 > 15.62$, strictly blocking rotation under *both* criteria.
  By mirror symmetry $x |-> -x$ (up to the $0.1$ shift of $a_2$), both rotations on $p_2$ are equally blocked, trapping PERCH at $a_1 or a_2 = "root"$ ($C_1$-purity $(1 / 2)(1 + 2 / 42) approx 0.524$). A single GRINCH graft from $a_1$ finds nearest node $a_2$ ($d = 0.1 << 10$), detaches $a_1$, and forms $(a_1, a_2)$, restoring purity $1.0$ across all five clusters. Moreover, rotating a unit-thickness slab diagonally along $(1, dots, 1) / sqrt(d) in bb(R)^d$ inflates AABB corner distance by $sqrt(d)$, whereas GRINCH's bilinear summary $chevron.l mu_A, mu_B chevron.r$ is strictly $O(d)$-invariant.

  *Part 2:* Every active cluster $i in {1, dots, K_t}$ has an incident 1-NN edge in $A$, so every connected component has size $>= 2$, giving $2 K_(t+1) <= K_t => K_(t+1) <= floor(K_t / 2)$ and exact round ceiling $R <= floor(log_2 n)$ (tight on separated dyadic pairs). Conversely, on $x_j = 2^j - 1$ ($j = 0, dots, m - 1$), $kappa_(x_j)^1 = x_(j-1)$ for all $j >= 1$, collapsing all $m$ exponentially spaced points into one component in round 1.

  *Part 3:* In the directed 1-NN graph on $K$ active clusters (out-degree $1$ under canonical tie-breaking), following $i -> kappa_i^1 -> dots$ yields non-increasing symmetric weights; canonical tie-breaking forbids directed cycles of length $>= 3$, so every weakly connected component drains into a *unique* reciprocal 2-cycle $u <-> v$. Contracting RNN 2-cycles (or $tau_r$-gated components) synchronously depends only on the symmetric similarity matrix, ensuring strict permutation invariance. Refining $tau_1 > dots > tau_R$ so one RNN pair crosses per step recovers sequential NN-Chain UPGMA (Theorem 1.2), while $tau = -infinity$ with one-way arcs recovers FINCH identically.
]

#wide_figure(
  caption: [*Table IV:* Head-to-Head Comparison of Divisive, Phylogenetic, and Extreme-Scale Incremental / Parallel Hierarchies (Algs. 63, 68, 75--77, 79--80).],
  {
    set text(size: 7.3pt)
    set table(inset: (x: 3.4pt, y: 2.7pt))
    table(
      columns: (1.02fr, 1.05fr, 0.90fr, 0.72fr, 0.78fr, 1.25fr),
      align: (left, left, center, center, center, left),
      [Algorithm], [Tree Mechanism], [Time Complexity], [Tree Geometry], [Online / Par.?], [Key Theoretical Tradeoff],
      [\#63 DIANA @kaufman1990finding], [Top-down splinter bisection], [$O(n^3)$ (worst)], [Superdom. ultra.], [Batch], [Fails exact ultrametric at $n = 6$ (Thm. 8.1)],
      [\#68 MONA @kaufman1990finding], [Monothetic 1-var. binary split], [$O(n d^3)$], [Rule subcubes], [Batch], [Exact rule interpretability; XOR-blind (Thm. 8.2)],
      [\#75 NJ @saitou1987neighbor], [Min-evolution $Q$-cherry merge], [$O(n^3)$], [Additive unrooted], [Batch], [Exact variable-rate phylogeny; poor on Euclid.],
      [\#76 PERCH @kobren2017perch], [AABB search + local rotations], [$O(N d log N)$], [Purity binary tree], [Online], [Local-rotation traps & high-$d$ AABB masking],
      [\#77 GRINCH @monath2019grinch], [Bilinear avg. + rotate + graft], [$O(N d log N)$], [Purity binary tree], [Online], [Escapes rotation traps; $O(d)$ orthogonal inv.],
      [\#79 SCC @monath2021scc], [Scheduled parallel RNN rounds], [$O(N d log N + R |E|)$], [Multi-level graph], [Parallel], [Order-invariant; bridges FINCH & UPGMA],
      [\#80 FINCH @sarfraz2019finch], [Param.-free 1-NN contraction], [$O(n log n)$ (ANN)], [Coarse $log_2 n$ lvls], [Parallel], [Zero hyperparams; asymmetric chain collapse],
    )
  },
)

// =============================================================================
// SECTION V
// =============================================================================
#sec("V", "Complete 25-Algorithm Head-to-Head Comparison Matrix and Decision Guide")

Synthesizing Sections I--IV, Table V consolidates the mathematical properties, algorithmic complexity, ultrametricity, reducibility, and structural regimes of *all twenty-five* Category III hierarchical and tree-based algorithms (Algs. 56--80), followed by a theorem-driven practitioner's decision guide and open theoretical frontiers.

#wide_figure(
  caption: [*Table V:* Master Head-to-Head Theoretical and Practical Comparison Matrix of All 25 Hierarchical and Tree-Based Algorithms (Category III, Algorithms 56--80). Abbreviations: *LW* = Lance--Williams Combinatorial, *AC/PR* = Anti-Chaining / Prototype Linkage, *SC* = Scalable Multi-Stage / Metric / Hybrid, *DV/PH* = Divisive / Phylogenetic, *EX* = Extreme-Scale Incremental / Parallel; *Red.* = Bruynooghe Reducible; *Ult.* = Inversion-Free Ultrametric.],
  {
    set text(size: 6.8pt)
    set table(inset: (x: 3.0pt, y: 2.2pt))
    table(
      columns: (1.20fr, 0.44fr, 1.04fr, 0.50fr, 0.44fr, 0.64fr, 1.48fr),
      align: (left, center, center, center, center, center, left),
      [Algorithm], [Grp.], [Time Complexity], [Ult.?], [Red.?], [Domain], [Defining Mathematical Property / Failure Mode],
      [\#56 SLINK @sibson1973slink], [LW], [$O(n^2)$], [Yes], [Yes], [Dissim.], [Unique subdom. ultrametric / MST; $0$-brk. chaining],
      [\#57 CLINK @defays1977clink], [LW], [$O(n^2)$], [Yes], [Yes], [Dissim.], [Min. superdom. clique hierarchy; $Omega(sqrt(n))$ Dasg. cost],
      [\#58 UPGMA @sokal1958statistical], [LW], [$O(n^2)$], [Yes], [Yes], [Dissim.], [Exact $1 / 3$ MW revenue; additive variance; clock req.],
      [\#59 Ward @ward1963hierarchical], [LW], [$O(n^2)$], [Yes], [Yes], [Euclid.], [ESS-greedy space dilation; equal-size spherical bias],
      [\#60 WPGMA @sokal1958statistical], [LW], [$O(n^2)$], [Yes], [Yes], [Dissim.], [Dyadic $2^(-"depth")$ weights; caterpillar outlier amp.],
      [\#61 UPGMC @sokal1958statistical], [LW], [$O(n^2 log n)$], [No], [No], [Euclid.], [$>= 3 / 4$ inversions ($d >= 2$); size-domination masking],
      [\#62 WPGMC @lance1967general], [LW], [$O(n^2 log n)$], [No], [No], [Euclid.], [Anti-masking midpoint; persistent $3 / 4$ inversions],
      [\#63 DIANA @kaufman1990finding], [DV/PH], [$O(n^3)$], [Yes], [N/A], [Dissim.], [Top-down superdom. tree; fails ultrametric at $n = 6$],
      [\#64 BIRCH @zhang1996birch], [SC], [$O(n d B log_B m)$], [Ph. 2], [N/A], [Euclid.], [Single-pass additive CF-tree; unbounded $k$-means dist.],
      [\#65 CURE @guha1998cure], [SC], [$O(n_s^2 log n_s + n k c d)$], [No], [No], [Euclid.], [$c$ shrunk reps; collapses concentric manifolds],
      [\#66 CHAMELEON @karypis1999chameleon], [SC], [$O(n^2 d + m^2 log m)$], [No], [No], [Metric], [METIS + RI/RC scale invariance; rejects thin bridges],
      [\#67 ROCK @guha2000rock], [SC], [$O(n^2 + n m_m m_a)$], [No], [No], [Categ.], [PSD link kernel $A^2$; Kernel $k$-Means at $f = 1 / 2$],
      [\#68 MONA @kaufman1990finding], [DV/PH], [$O(n d^3)$], [Step], [N/A], [Binary], [Exact covariance 1-var. rules; blind to XOR parity],
      [\#69 Genie @gagolewski2016genie], [AC/PR], [$O(n^2)$ / $O(n log n)$], [No], [No], [Dissim.], [MST edges + Gini gate absorbs outliers; inversions],
      [\#70 HCA-DBSCAN @ester1996dbscan], [SC], [$O(n log n + m^2)$], [Ph. 2], [Ph. 2], [Metric], [Core extraction gives unbounded bridge immunity],
      [\#71 BUBBLE @ganti1999bubble], [SC], [$O(n(B log_B n_L + s))$], [Ph. 2], [N/A], [Metric], [Metric CR-tree + FastMap; factor-2 clustroid bound],
      [\#72 Lance--Williams @lance1967general], [LW], [$O(n^2)$ / $O(n^2 log n)$], [Thm. 1], [Thm. 1], [Dissim.], [Universal 4-param. combinatorial linkage umbrella],
      [\#73 Minimax @bien2011minimax], [AC/PR], [$O(n^3)$], [Yes], [Yes], [Dissim.], [Non-LW prototype $B(p_r, h)$ balls; factor-2 diam.],
      [\#74 Flexible Beta @lance1967general], [LW], [$O(n^2)$ / $O(n^2 log n)$], [$beta < 1$], [$beta <= 0$], [Dissim.], [Continuous dilation dial; metric iff $beta = 0$],
      [\#75 Neighbor-Joining @saitou1987neighbor], [DV/PH], [$O(n^3)$], [Addit.], [N/A], [Addit.], [Exact variable-rate additive tree; poor on Euclid.],
      [\#76 PERCH @kobren2017perch], [EX], [$O(N d log N)$], [No], [N/A], [Euclid.], [Online AABB + local rotations; rotation-trap stuck],
      [\#77 GRINCH @monath2019grinch], [EX], [$O(N d log N)$], [No], [N/A], [Euclid.], [Online bilinear avg. + graft escapes local traps],
      [\#78 BICO @fichtenberger2013bico], [SC], [$O(n d log s + s^2 d)$], [Ph. 2], [N/A], [Euclid.], [Streaming CF-tree with strong $(k, epsilon)$-coreset],
      [\#79 SCC @monath2021scc], [EX], [$O(N d log N + R |E|)$], [No], [N/A], [Metric], [Order-invariant parallel RNN rounds; bridges UPGMA],
      [\#80 FINCH @sarfraz2019finch], [EX], [$O(n log n)$ (ANN)], [No], [N/A], [Metric], [Parameter-free $<= floor(log_2 n)$ rounds; chain risk],
    )
  },
)

#subsec("A", "Mathematically Grounded Practitioner's Decision Guide")

Choosing effectively among Algorithms 56--80 follows directly from the invariants established in Theorems 1--10:

1. *Regime 1 — Moderate-Scale General Dissimilarities ($n <= 5 times 10^4$):* Default to *UPGMA (\#58)* ($O(n^2)$ NN-Chain) for an inversion-free ultrametric with provable $1 / 3$ Moseley--Wang global revenue (Theorem 5.1). Use *Minimax Linkage (\#73)* when $n <= 5 times 10^3$ and every cluster needs a data-anchored covering ball $B(p_r, h)$ (Theorem 3.1). Choose *Genie (\#69)* ($g = 0.3$) on non-convex manifolds with isolated outliers for fast $O(n log n)$ post-MST sweeps (Theorem 4.2), or *Flexible Beta (\#74)* with $beta in [-0.35, -0.20]$ (strictly reducible, Theorem 1.3) for Ward-like dilation on non-Euclidean dissimilarities.
2. *Regime 2 — Moderate-Scale Euclidean Spherical Data ($n <= 5 times 10^4$):* Deploy *Ward's Method (\#59)* ($O(n^2)$ NN-Chain) for noisy balanced Gaussian spheres (Theorem 2.1), upgrading to *HCA-DBSCAN (\#70)* under dense noise bridges (Theorem 7.4). Avoid *UPGMC (\#61)* and *WPGMC (\#62)* in $d >= 2$ due to dendrogram inversions and $O(n^2 log n)$ heap overhead (Theorems 1.3 and 2).
3. *Regime 3 — Streaming & Coreset Compression ($n >= 10^6$):* Prefer *BICO (\#78)* over *BIRCH (\#64)* in Euclidean pipelines to retain single-pass $O(n d log s)$ speed while guaranteeing an exact strong $(k, epsilon)$-coreset (Theorem 6.1--6.2). In non-Euclidean metric spaces (edit distance, DTW @sakoe1978dtw, graphs), deploy *BUBBLE / BUBBLE-FM (\#71)* (Theorem 6.3).
4. *Regime 4 — Non-Spherical Shapes, Multi-Density Manifolds, & Categorical Tables:* Use *CURE (\#65)* for non-nested elongated Euclidean shapes; switch to *CHAMELEON (\#66)* (homothety invariance and bridge rejection, Theorem 7.2) or *HCA-DBSCAN (\#70)* whenever concentric shells or multi-density contrasts occur. For sparse categorical/transactional sets, deploy *ROCK (\#67)* near $f(theta) approx 1 / 2$ (Theorem 7.3); for rule-interpretable binary trees without parity structure, deploy *MONA (\#68)* (Theorem 8.2).
5. *Regime 5 — Non-Clock Phylogenetics:* Deploy *Neighbor-Joining (\#75)* whenever evolutionary rates vary across lineages (exact additive recovery vs. UPGMA/WPGMA long-branch attraction, Theorem 9.1--9.2), and never apply NJ to Euclidean geometric clustering (Theorem 9.3).
6. *Regime 6 — Web-Scale Extreme Hierarchies ($N >= 10^6$, $K >= 10^4$):* For online insertion, prefer *GRINCH (\#77)* over *PERCH (\#76)* for orthogonal invariance and graft escape from local rotation traps (Theorem 10.1). For parallel batch ANN graphs, use *FINCH (\#80)* for zero-hyperparameter coarse hierarchies ($<= floor(log_2 n)$ levels, Theorem 10.2) and *SCC (\#79)* for order-invariant multi-resolution control bridging FINCH and UPGMA (Theorem 10.3).
7. *Regime 7 — Explicit Elimination Rules:* Avoid *SLINK (\#56)* on unfiltered noisy point clouds ($0$-breakdown chaining, Theorem 3.3), *CLINK (\#57)* when global Dasgupta cost matters (Theorem 5.3), *DIANA (\#63)* above $n > 3000$ or on unbalanced clades (Theorem 8.1), *Flexible Beta (\#74)* with $beta > 0$ under NN-Chain (Theorem 1.3), and *CURE (\#65)* with $alpha >= 1 - R_1 / R_2$ on concentric shells (Theorem 7.1).

#subsec("B", "Synthesis of Open Theoretical Frontiers Across Category III")

Unifying Algorithms 56--80 highlights four open mathematical frontiers:
1. *Agglomerative Dasgupta Cost Approximation:* While UPGMA (\#58) guarantees $1 / 3$ Moseley--Wang revenue @moseley2017approximation and SLINK/CLINK fail polynomially (Theorem 5), designing an $O(n^2)$ or subquadratic bottom-up linkage achieving $O(1)$ or $O(sqrt(log n))$ approximation to *Dasgupta cost* $cal(L)_"Dasg"$ without top-down sparsest-cut solvers @dasgupta2016cost remains open.
2. *Finite-Budget Online Restructuring Guarantees:* Establishing worst-case $cal(L)_"Dasg"$ or $cal(R)_"MW"$ bounds for GRINCH (\#77) under capped graft-search budgets and adversarial arrival orders remains unresolved.
3. *Strong Coreset Hierarchies Beyond Euclidean $k$-Means:* Extending BICO's (\#78) streaming coreset theory from spherical $k$-means to preserve full Ward, Minimax, or Dasgupta hierarchical objectives in arbitrary metric spaces (unifying BICO and BUBBLE \#71) is an active frontier.
4. *Cross-Paradigm Density--Kernel--Spectral Bridges:* SLINK (\#56) anchors HDBSCAN\*'s mutual-reachability MST @campello2015hierarchical (Comparison Chapter IV), ROCK (\#67) at $f(theta) = 1 / 2$ reduces exactly to Kernel $k$-Means @dhillon2004kernel (Theorem 7.3), and CHAMELEON's METIS cuts connect agglomeration to Normalized Cut spectral relaxations @shi2000ncut; these cross-paradigm equivalences and Kleinberg's axiomatic trade-offs @kleinberg2002impossibility are synthesized in Comparison Chapters I, II, and XII--XIII of Part IX.
