#import "../template/ieee.typ": *

#alg_entry(
  id: 62,
  name: "Median Linkage (WPGMC)",
  category: "Cat. III: Hierarchical & Tree-Based",
  year: "1967",
  paradigm: "Agglomerative Weighted-Centroid Linkage",
  complexity: $cal(O)(n^2 log n)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Let $X = \{bold(x)_1, dots, bold(x)_n\} subset bb(R)^p$ be a finite Euclidean point set with initial squared pairwise dissimilarity matrix $D^((0))_(i j) = norm(bold(x)_i - bold(x)_j)_2^2$. In Centroid Linkage (UPGMC, Alg. 61), merging a large cluster $A$ with a small cluster $B$ places the new barycenter close to $bold(mu)_A$, causing the geometry of the large cluster to mask that of the small cluster. To eliminate this cardinality-domination pathology, Gower (1967) and Lance and Williams @lance1967general introduced *Median Linkage*, formally the *Weighted Pair Group Method with Centroid* (*WPGMC*).

Rather than weighting each merging cluster by its sample size, WPGMC defines the representative of $A union B$ as the *exact Euclidean midpoint* of the parent representatives $bold(m)_A$ and $bold(m)_B$. Crucially, the historical designation *median* does *not* refer to the coordinate-wise $L_1$ statistical median; instead, it refers to the *geometric median of a triangle* in Euclidean geometry (the line segment connecting a third cluster's representative $bold(m)_C$ to the midpoint $bold(m)_(A union B)$ of segment $overline(bold(m)_A bold(m)_B)$), whose squared length is governed exactly by Apollonius's median theorem. Because equal weight $1/2$ is assigned to each merging branch regardless of cluster cardinality, individual observations receive unequal, tree-depth-dependent dyadic weights — hence the name *weighted* centroid linkage.

#definition(title: "WPGMC Pseudo-Centroid and Apollonius Linkage Recurrence")[
  Define the *dyadic pseudo-centroid* $bold(m)_A in bb(R)^p$ recursively along the merge tree $T_A$ of cluster $A subset.eq X$ by
  $ bold(m)_(\{bold(x)_i\}) = bold(x)_i, quad bold(m)_(A union B) = 1 / 2 bold(m)_A + 1 / 2 bold(m)_B. $
  The Euclidean and squared-Euclidean WPGMC dissimilarities between disjoint clusters $A, B$ are
  $ d_M (A, B) = norm(bold(m)_A - bold(m)_B)_2, quad D_M (A, B) = norm(bold(m)_A - bold(m)_B)_2^2. $
  Starting from $cal(C)^((0)) = \{\{bold(x)_1\}, dots, \{bold(x)_n\}\}$, at step $t in \{1, dots, n - 1\}$ WPGMC merges a minimizing pair
  $ (A_t, B_t) in arg min_(U != V in cal(C)^((t-1))) D_M (U, V), $
  records merge height $h_t = sqrt(D_M (A_t, B_t))$ (squared height $H_t = h_t^2$), and updates squared dissimilarities to every remaining cluster $C$ via the size-independent Lance--Williams recurrence ($alpha_A = alpha_B = 1/2$, $beta = -1/4$, $gamma = 0$):
  $ D_M (A_t union B_t, C) = 1 / 2 D_M (A_t, C) + 1 / 2 D_M (B_t, C) - 1 / 4 D_M (A_t, B_t). $
]

Because $beta = -1/4 < 0$ violates Bruynooghe reducibility unconditionally, WPGMC cannot use the $cal(O)(n^2)$ Nearest-Neighbor Chain algorithm @mullner2011modern and instead executes deterministically in $cal(O)(n^2 log n)$ worst-case time via indexed min-heaps (or near-$cal(O)(n^2)$ average time via Müllner's generic priority-queue algorithm).

#algbox(
  title: "Exact Median Linkage (WPGMC) via Indexed Heaps",
  inputs: [Points $X = \{bold(x)_1, dots, bold(x)_n\} subset bb(R)^p$ (or squared Euclidean matrix $D^((0))$).],
  outputs: [Merge list $\{(A_t, B_t, h_t)\}_(t=1)^(n-1)$ defining dendrogram $T$.],
)[
  1. Set $D_M (\{bold(x)_i\}, \{bold(x)_j\}) <- norm(bold(x)_i - bold(x)_j)_2^2$ for all $1 <= i < j <= n$; initialize active set $cal(C) <- \{\{bold(x)_1\}, dots, \{bold(x)_n\}\}$.
  2. Build per-cluster indexed min-heaps $Q_U$ over $\{(D_M (U, V), V) : V in cal(C) without \{U\}\}$ and global min-heap $G$ over $\{(min(Q_U), U) : U in cal(C)\}$.
  3. For $t = 1, dots, n - 1$:
     a. Extract global minimizing pair $(A, B)$ with squared distance $H_t = D_M (A, B)$ from $G$; set $h_t <- sqrt(H_t)$ and append $(A, B, h_t)$.
     b. Form $W <- A union B$ and remove $A, B$ from $cal(C)$ and $G$.
     c. For every remaining $C in cal(C)$:
        - Delete $A$ and $B$ from $Q_C$.
        - Update $D_M (W, C) <- (1/2) D_M (A, C) + (1/2) D_M (B, C) - (1/4) H_t$.
        - Insert $(D_M (W, C), W)$ into $Q_C$, update $C$'s key in $G$, and append $(D_M (W, C), C)$ to buffer $L_W$.
     d. Heapify $L_W$ into $Q_W$ in $cal(O)(|cal(C)|)$ time, insert $(min(Q_W), W)$ into $G$, and add $W$ to $cal(C)$.
  4. Return $\{(A_t, B_t, h_t)\}_(t=1)^(n-1)$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

#subsec("A", "Apollonius Parallelogram Identity and Dyadic Convexity")

Our first lemma proves the exact geometric foundation of WPGMC and characterizes its pseudo-centroid as a convex dyadic average over cluster leaves.

#lemma(title: "Apollonius Median Identity and Dyadic Pseudo-Centroid Expansion", origin: "Original Derivation")[
  Let $T_A$ be the rooted binary merge tree of cluster $A subset.eq X$, and let $ell(bold(x), A)$ be the edge depth of leaf $bold(x) in A$ below the root of $T_A$. Then:
  1. The WPGMC pseudo-centroid admits the exact convex-hull representation
     $ bold(m)_A = sum_(bold(x) in A) 2^(-ell(bold(x), A)) bold(x) in "conv"(A), quad sum_(bold(x) in A) 2^(-ell(bold(x), A)) = 1. $
  2. For any disjoint clusters $A, B, C subset bb(R)^p$, the squared pseudo-centroid distance obeys the exact size-independent Lance--Williams recurrence
     $ D_M (A union B, C) = 1 / 2 D_M (A, C) + 1 / 2 D_M (B, C) - 1 / 4 D_M (A, B). $
]

#proof[
  For (1), proceed by induction on $|A|$. For $|A| = 1$, $ell(bold(x), \{bold(x)\}) = 0$ and $bold(m)_(\{bold(x)\}) = 2^0 bold(x) = bold(x)$. If $A = A_1 union A_2$, every $bold(x) in A_1$ satisfies $ell(bold(x), A) = ell(bold(x), A_1) + 1$ and every $bold(x) in A_2$ satisfies $ell(bold(x), A) = ell(bold(x), A_2) + 1$, so
  $ bold(m)_A &= 1 / 2 bold(m)_(A_1) + 1 / 2 bold(m)_(A_2) \ &= sum_(bold(x) in A_1) 2^(-(ell(bold(x), A_1) + 1)) bold(x) + sum_(bold(x) in A_2) 2^(-(ell(bold(x), A_2) + 1)) bold(x) \ &= sum_(bold(x) in A) 2^(-ell(bold(x), A)) bold(x). $
  Weights are strictly positive and sum to $(1/2)(1) + (1/2)(1) = 1$ by induction, guaranteeing $bold(m)_A in "conv"(A)$.

  For (2), define $bold(u) = bold(m)_A - bold(m)_C$ and $bold(v) = bold(m)_B - bold(m)_C$. Then $bold(m)_(A union B) - bold(m)_C = (bold(u) + bold(v)) / 2$ and $bold(m)_A - bold(m)_B = bold(u) - bold(v)$. Applying the Hilbert-space parallelogram law $norm(bold(u) + bold(v))_2^2 + norm(bold(u) - bold(v))_2^2 = 2 norm(bold(u))_2^2 + 2 norm(bold(v))_2^2$ and dividing by $4$ yields
  $ norm((bold(u) + bold(v)) / 2)_2^2 = 1 / 2 norm(bold(u))_2^2 + 1 / 2 norm(bold(v))_2^2 - 1 / 4 norm(bold(u) - bold(v))_2^2, $
  which is precisely $D_M (A union B, C) = (1/2) D_M (A, C) + (1/2) D_M (B, C) - (1/4) D_M (A, B)$.
]

#subsec("B", "Exact Dyadic-Variance Gap Linking WPGMC and WPGMA")

Just as UPGMC differs from squared-Euclidean UPGMA by the sum of within-cluster variances, WPGMC satisfies an exact dyadic-variance conservation law linking it to squared-Euclidean WPGMA (Alg. 60).

#theorem(title: "Dyadic Variance Decomposition and Conservation Law with WPGMA", origin: "Original Derivation")[
  For any cluster $A$ with merge subtree $T_A$ and dyadic leaf measure $mu_A (bold(x)) = 2^(-ell(bold(x), A))$, define the *dyadic trace variance*
  $ V_"dyad" (A) = sum_(bold(x) in A) mu_A (bold(x)) norm(bold(x) - bold(m)_A)_2^2. $
  Let $D_W^((2)) (A, B) = sum_(bold(x) in A, bold(y) in B) mu_A (bold(x)) mu_B (bold(y)) norm(bold(x) - bold(y))_2^2$ be the squared-Euclidean WPGMA dissimilarity evaluated on the same subtrees $T_A, T_B$. Then:
  1. $D_M (A, B) = D_W^((2)) (A, B) - V_"dyad" (A) - V_"dyad" (B)$.
  2. $V_"dyad"$ satisfies the exact recursive merge law
     $ V_"dyad" (A union B) = 1 / 2 V_"dyad" (A) + 1 / 2 V_"dyad" (B) + 1 / 4 D_M (A, B). $
  3. Adding (2) and Lemma 1(2) cancels the $-1/4$ and $+1/4$ terms identically, recovering the squared-WPGMA recurrence $D_W^((2)) (A union B, C) = (1/2) D_W^((2)) (A, C) + (1/2) D_W^((2)) (B, C)$.
]

#proof[
  Let $bold(X) ~ mu_A$ and $bold(Y) ~ mu_B$ be independent random vectors, so $bb(E)[bold(X)] = bold(m)_A$, $bb(E)[bold(Y)] = bold(m)_B$, and $bb(E)[norm(bold(X) - bold(m)_A)_2^2] = V_"dyad" (A)$. Writing $bold(X) - bold(Y) = (bold(X) - bold(m)_A) + (bold(m)_A - bold(m)_B) - (bold(Y) - bold(m)_B)$ and taking expectations of the squared norm, independence and zero-mean centering annihilate all cross inner products, yielding
  $ D_W^((2)) (A, B) = V_"dyad" (A) + norm(bold(m)_A - bold(m)_B)_2^2 + V_"dyad" (B), $
  proving (1).

  For (2), note that $mu_(A union B) = (1/2) mu_A + (1/2) mu_B$ is an equal mixture of $mu_A$ and $mu_B$. By the law of total variance,
  $ V_"dyad" (A union B) &= 1 / 2 (V_"dyad" (A) + norm(bold(m)_A - bold(m)_(A union B))_2^2) \ &quad + 1 / 2 (V_"dyad" (B) + norm(bold(m)_B - bold(m)_(A union B))_2^2). $
  Since $bold(m)_A - bold(m)_(A union B) = (bold(m)_A - bold(m)_B)/2$ and $bold(m)_B - bold(m)_(A union B) = -(bold(m)_A - bold(m)_B)/2$, both squared shifts equal $(1/4) D_M (A, B)$, giving $V_"dyad" (A union B) = (1/2) V_"dyad" (A) + (1/2) V_"dyad" (B) + (1/4) D_M (A, B)$.

  For (3), summing Lemma 1(2) with (2) and $V_"dyad" (C) = (1/2) V_"dyad" (C) + (1/2) V_"dyad" (C)$ cancels $-(1/4) D_M (A, B)$ with $+(1/4) D_M (A, B)$ to give $D_W^((2)) (A union B, C) = (1/2) D_W^((2)) (A, C) + (1/2) D_W^((2)) (B, C)$.
]

#subsec("C", "Size-Invariant Reversals, Non-Negativity, and Simplex Inversion Profile")

Unlike UPGMC, whose negative Lance--Williams coefficient $beta = -n_A n_B / (n_A + n_B)^2$ vanishes on imbalanced merges, WPGMC locks $beta = -1/4$ across all cluster cardinalities. We characterize both its persistence of reversals under arbitrarily skewed cluster sizes and its surprising non-negativity preservation on general dissimilarity matrices.

#theorem(title: "Size-Invariant Inversion Persistence, Non-Negativity, and Simplex Height Decay", origin: "Original Proof")[
  1. *Universal non-negativity & sharp single-step bound*: Executed greedily on *any* symmetric non-negative dissimilarity matrix $D^((0)) in bb(R)_(>= 0)^(n times n)$ (even non-Euclidean), WPGMC preserves $D_M >= 0$ at every step and obeys the sharp bound
     $ H_(t+1) >= 3 / 4 H_t quad (h_(t+1) >= sqrt(3) / 2 h_t). $
     In $p = 1$, WPGMC is strictly inversion-free ($h_(t+1) >= h_t$); for all $p >= 2$, inversions persist even when $|A| / |B| -> infinity$.
  2. *Unbalanced inversion separation from UPGMC*: On $A = \{(-1, 0)^((m))\}$ ($m$ coinciding points), $B = \{(1, 0)\}$, $C = \{(0, 1.8)\}$ in $bb(R)^2$, WPGMC suffers a fixed $19\%$ squared-height inversion ($H_2 = 3.24 < 4.00 = H_1$) for every $m >= 1$, whereas UPGMC is completely inversion-free for all $m >= 15$.
  3. *Regular-simplex inversion asymptotics*: On the regular $(n - 1)$-simplex of unit edge length in $bb(R)^n$ (with infinitesimal tie-breaking), sequential caterpillar merges in WPGMC generate $n - 2$ consecutive strict inversions obeying the exact closed form
     $ H_m^"WPGMC" = 2 / 3 + 1 / 3 dot 4^(-(m-1)) arrow.b 2 / 3 quad (h_m arrow.b sqrt(2 / 3)), $
     strictly exceeding UPGMC's simplex floor $H_m^"UPGMC" = (m + 1) / (2 m) arrow.b 1 / 2$.
]

#proof[
  For (1), prove $D_M >= 0$ and $H_(t+1) >= (3/4) H_t$ simultaneously by induction on $t$. At step $t$, $(A_t, B_t)$ minimizes $D_M$ over $cal(C)^((t-1))$, so $D_M (A_t, C) >= H_t$ and $D_M (B_t, C) >= H_t$ for all $C in cal(C)^((t-1)) without \{A_t, B_t\}$. Thus
  $ D_M (A_t union B_t, C) &= 1 / 2 D_M (A_t, C) + 1 / 2 D_M (B_t, C) - 1 / 4 H_t \ &>= (1 / 2 + 1 / 2 - 1 / 4) H_t = 3 / 4 H_t >= 0, $
  while every untouched pair already satisfies $D_M (U, V) >= H_t >= (3/4) H_t >= 0$. Hence $H_(t+1) >= (3/4) H_t >= 0$. In $p = 1$, since $m_(A_t union B_t) = (m_(A_t) + m_(B_t)) / 2 in [m_(A_t), m_(B_t)]$ and no cluster representative $m_C$ can lie strictly inside $(m_(A_t), m_(B_t))$ without contradicting minimality of $|m_(B_t) - m_(A_t)|$, we have $|m_C - m_(A_t union B_t)| >= h_t$ for all $C$.

  For (2), after merging the $m$ copies of $(-1, 0)$ at height $0$, active centers are $A$ at $(-1, 0)$ (size $m$), $B$ at $(1, 0)$ (size $1$), and $C$ at $(0, 1.8)$ (size $1$), with $D(A, B) = 4.00$ and $D(A, C) = D(B, C) = 4.24$. Step $m$ merges $A union B$ at $H = 4.00$. Under WPGMC, $D_M (A union B, C) = 4.24 - (1/4)(4.00) = 3.24 < 4.00$ for *every* $m >= 1$ ($19\%$ drop). Under UPGMC, $D_C (A union B, C) = 4.24 - 4 m / (m + 1)^2$; for every $m >= 15$, $4 m / (m + 1)^2 <= 60 / 256 = 0.2344 < 0.24$, so $D_C (A union B, C) > 4.00$ and no inversion occurs.

  For (3), place $bold(x)_i = (1 / sqrt(2)) bold(e)_i in bb(R)^n$ ($norm(bold(x)_i - bold(x)_j)_2^2 = 1$). Let $a_m = norm(bold(m)_m)_2^2$ for the caterpillar cluster $S_m = (dots(\{bold(x)_1\} union \{bold(x)_2\}) dots union \{bold(x)_m\})$. Since $bold(m)_(m+1) = (1/2) bold(m)_m + (1/2) bold(x)_(m+1)$ with $bold(m)_m perp bold(x)_(m+1)$ and $norm(bold(x)_(m+1))_2^2 = 1/2$, we obtain the linear recurrence
  $ a_(m+1) = 1 / 4 a_m + 1 / 8, quad a_1 = 1 / 2. $
  Its unique solution is $a_m = 1/6 + (1/3) dot 4^(-(m-1))$ (verified at $m = 1$: $1/6 + 1/3 = 1/2$, and $a_(m+1) = 1/24 + (1/3) dot 4^(-m) + 1/8 = 1/6 + (1/3) dot 4^(-m)$). Adding $norm(bold(x)_(m+1))_2^2 = 1/2$ yields $H_m^"WPGMC" = a_m + 1/2 = 2/3 + (1/3) dot 4^(-(m-1)) arrow.b 2/3$.
]

#subsec("D", "Anti-Masking Step Bound and Minimal Divergence from UPGMC")

We now quantify the exact geometric tradeoff between WPGMC's anti-masking property and its heightened sensitivity to late-joining outliers, and establish the minimal cardinality at which WPGMC and UPGMC separate topologically.

#proposition(title: "Anti-Masking Midpoint Shift, Late-Outlier Amplification, and n = 5 Divergence", origin: "Original Counterexample")[
  1. *Anti-masking and outlier amplification*: Regardless of cluster cardinalities $|A|$ and $|B|$, merging $A$ and $B$ shifts each representative by exactly half their separation:
     $ norm(bold(m)_(A union B) - bold(m)_A)_2 = norm(bold(m)_(A union B) - bold(m)_B)_2 = 1 / 2 d_M (A, B). $
     Consequently, if a single late-merging outlier $bold(y)$ joins a cluster $A$ of size $m$, WPGMC amplifies the representative shift by a factor of $(m + 1) / 2$ relative to UPGMC.
  2. *Minimal $n = 5$ divergence from UPGMC*: For all $n <= 4$, WPGMC and UPGMC produce identical merge topologies (up to ties). For $n = 5$, on the integer line dataset $X = \{0, 2, 9, 18, 32\} subset bb(R)$, WPGMC and UPGMC yield strictly distinct dendrograms and 2-cluster cuts.
]

#proof[
  (1) Immediate from $bold(m)_(A union B) - bold(m)_A = (bold(m)_B - bold(m)_A) / 2$, giving shift $(1/2) d_M (A, B)$ versus UPGMC's $(1 / (m + 1)) d_C (A, \{bold(y)\})$, whose ratio is $(m + 1) / 2$.

  (2) Whenever $|A| = |B| = 1$, $bold(m)_(A union B) = bold(mu)_(A union B)$, so WPGMC and UPGMC coincide identically until a cluster of size $>= 3$ competes against at least two remaining clusters, which requires $n >= 5$. On $X = \{0, 2, 9, 18, 32\} subset bb(R)$:
  - *Step 1 (identical)*: Unique minimum $d(0, 2) = 2$ merges $C_(1 2) = \{0, 2\}$ with center $m_(1 2) = mu_(1 2) = 1$.
  - *Step 2 (identical)*: Active distances are $d(1, 9) = 8 < d(9, 18) = 9 < d(18, 32) = 14$, so both merge $S = C_(1 2) union \{9\}$. Their representatives now diverge: $m_S = (1 + 9)/2 = 5$, whereas $mu_S = (0 + 2 + 9)/3 = 11/3 approx 3.667$.
  - *Step 3 (divergence)*: Comparing distance from $S$ to $\{18\}$ against $d(18, 32) = 14$:
    $ d_M (S, \{18\}) = |18 - 5| = 13 < 14, quad d_C (S, \{18\}) = |18 - 11 / 3| = 43 / 3 approx 14.333 > 14. $
    Thus WPGMC merges $S union \{18\}$ to output the 2-cluster cut $\{\{0, 2, 9, 18\}, \{32\}\}$, whereas UPGMC merges $\{18, 32\}$ to output $\{\{0, 2, 9\}, \{18, 32\}\}$.
]

#subsec("E", "NN-Chain Failure and Exact Complexity Barrier")

#proposition(title: "Invalidity of NN-Chain for WPGMC and Deterministic Complexity", origin: "Classical — Cited")[
  Because WPGMC has $beta = -1/4 < 0$, it fails Bruynooghe reducibility @mullner2011modern, and the NN-chain algorithm fails to recover the exact WPGMC dendrogram on $X = \{(-2, 0), (2, 0), (0, 7.2), (0, 3.5)\} subset bb(R)^2$. Exact WPGMC executes in $cal(O)(n^2 log n)$ worst-case time and $cal(O)(n^2)$ space via indexed min-heaps.
]

#proof[
  Compute squared distances on $X$: $D(bold(x)_3, bold(x)_4) = 3.7^2 = 13.69$, $D(bold(x)_1, bold(x)_2) = 16.00$, $D(bold(x)_1, bold(x)_4) = D(bold(x)_2, bold(x)_4) = 4 + 12.25 = 16.25$, and $D(bold(x)_1, bold(x)_3) = D(bold(x)_2, bold(x)_3) = 4 + 51.84 = 55.84$. Exact WPGMC merges $C_(3 4) = \{bold(x)_3, bold(x)_4\}$ first ($H_1 = 13.69$, $bold(m)_(3 4) = (0, 5.35)$), then $C_(1 2) = \{bold(x)_1, bold(x)_2\}$ ($H_2 = 16.00 < 32.6225$), giving $\{\{bold(x)_1, bold(x)_2\}, \{bold(x)_3, bold(x)_4\}\}$. Conversely, starting NN-chain at $bold(x)_1$ finds reciprocal pair $(bold(x)_1, bold(x)_2)$ ($16.00 < 16.25$), merges $C_(1 2)$ with midpoint $(0, 0)$, and yields $D_M (C_(1 2), bold(x)_4) = 3.5^2 = 12.25 < 13.69$, causing NN-chain to merge $C_(1 2) union \{bold(x)_4\}$ and return the wrong 2-cluster cut. Indexed-heap accounting mirrors Algorithm 61, taking $cal(O)((n - t) log n)$ per merge and $cal(O)(n^2 log n)$ total time.
]

#sec("III", "Complexity, Geometry and Operational Profile")

WPGMC shares UPGMC's $cal(O)(n^2 log n)$ worst-case complexity and $cal(O)(n^2)$ memory footprint, yet its inner update $D_M (A union B, C) = (D_M (A, C) + D_M (B, C)) / 2 - D_M (A, B) / 4$ requires no cluster-size tracking or division (only additions and exponent shifts). Standard libraries such as `fastcluster` @mullner2011modern execute Müllner's generic priority-queue schedule, achieving near-$cal(O)(n^2)$ typical runtime while degrading toward $cal(O)(n^3)$ under frequent reversals; indexed per-cluster heaps eliminate that cubic worst case unconditionally.

Geometrically, each cluster is represented by $bold(m)_A in "conv"(A)$ that remains equally responsive to newly merged substructures regardless of size imbalance (Proposition 1(1)), eliminating UPGMC's masking on clean unequal-size boundaries. Conversely, because $beta = -1/4$ never decays with cluster-size imbalance (Theorem 2(2)), WPGMC induces more frequent inversions than UPGMC, requiring monotonic post-processing (e.g., $tilde(h)_t = max_(s <= t) h_s$) before threshold cuts.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Input Space], [Euclidean $bb(R)^p$ (combinatorial on squared $L_2$; $D_M >= 0$ for all $D^((0)) >= 0$)],
    [Lance--Williams], [$alpha_A = alpha_B = 1/2$, $beta = -1/4$, $gamma = 0$ (size-independent)],
    [Time Complexity], [$cal(O)(n^2 log n)$ worst-case (indexed heaps); $cal(O)(n^3)$ naive],
    [Space Complexity], [$cal(O)(n^2)$ condensed distance array + heap pointers],
    [Reducibility], [Always violated ($beta = -1/4$); inversions persist as $|A|/|B| -> infinity$],
    [Reversal Floor], [$h_(t+1) >= (sqrt(3)/2) h_t$ single-step; $h_m arrow.b sqrt(2/3)$ on unit simplex],
    [Invariance], [Rigid motions and isotropic scaling; non-invariant to affine shear],
    [Inductive Bias], [Equal-branch midpoint geometry; anti-masking, outlier-sensitive],
  ),
  caption: [Operational and mathematical profile of Median Linkage / WPGMC (Alg. 62).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Imbalanced Morphometric and Biological Taxonomy*: When classifying species or anatomical subspecies from continuous Euclidean morphometric measurements (such as flower sepal/petal dimensions or cranial landmarks), rare subspecies often have far fewer specimens than abundant sibling species. WPGMC's anti-masking midpoint rule (Proposition 1(1)) ensures that merging a densely sampled morphotype with an intermediate subpopulation preserves equal geometric representation for the sparser lineage, preventing the abundant group from absorbing the boundary.
2. *Multi-Scale Sensor Array Calibration and Spatial Interpolation*: In environmental monitoring networks where urban zones have dense sensor clusters and rural zones have sparse stations, fusing adjacent zones via unweighted centroids (UPGMC) pulls regional centers overwhelmingly into urban cores. WPGMC places hierarchical regional centers at the geometric midpoint of coalescing subregions, producing spatially balanced multi-resolution sensor hierarchies.
3. *Color Palette Quantization for Rare Accent Tones*: In perceptual CIELAB Euclidean color space, dominant background pixels vastly outnumber salient foreground accent hues. Standard centroid or Ward aggregation shifts combined prototypes toward the dominant background, washing out rare accents. WPGMC's equal-branch midpoint rule retains perceptual influence for low-frequency color clusters across intermediate merge levels.
4. *Balanced Hierarchical Topic and Document Prototype Trees*: When documents are embedded into normalized semantic vector spaces and sub-topics exhibit heavy-tailed document counts, UPGMC collapses niche sub-topics into mega-topics without shifting the mega-topic centroid. WPGMC keeps merged semantic prototypes centered between coalescing sub-topics, improving retrieval diversity across taxonomy branches @aggarwal2013data.
5. *Hierarchical Anomaly-Cleaned Prototype Synthesis*: In industrial quality control where known outlier batches are pre-screened via density filtering, WPGMC constructs hierarchical Euclidean prototypes across manufacturing regimes of unequal production volume, assigning equal geometric weight to low-volume specialty runs and high-volume baseline runs.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

WPGMC is constrained by four primary failure modes:
- *Severe late-outlier trajectory hijacking*: By Proposition 1(1), if an unremoved outlier joins a clean cluster of $m = 1,000$ points, WPGMC drags the cluster's pseudo-centroid halfway to the outlier — a $500.5 times$ greater distortion than UPGMC — catastrophically corrupting all subsequent merges.
- *Persistent dendrogram inversions under imbalance*: Unlike UPGMC, whose reversals dampen as $|A| / |B| -> infinity$, WPGMC maintains $beta = -1/4$ at every merge (Theorem 2(2)), generating frequent height reversals that break threshold-based dendrogram cutting.
- *Off-manifold pseudo-centroids on non-convex shapes*: Although $bold(m)_A in "conv"(A)$ (Lemma 1(1)), the convex hull of a curved crescent or ring includes empty interior space, causing WPGMC to mismerge interlocking shapes just like UPGMC.
- *Dependence on arbitrary merge history*: Because leaf weights $2^(-ell(bold(x), A))$ depend exponentially on internal tree balance, minor tie-breaking changes early in the hierarchy can alter downstream pseudo-centroids substantially.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Clean Imbalanced Morphometrics], [High: anti-masking tracks unequal-size class boundaries], [Outliers present: pre-filter or use UPGMA / Ward @ward1963hierarchical],
    [Unevenly Sampled Spatial Regions], [Moderate--High: balances sparse vs. dense zones], [Monotone dendrogram needed: use WPGMA (Alg. 60)],
    [Noisy / Heavy-Tailed Measurements], [Low: $(m+1)/2$ late-outlier shift amplification], [Use PAM / medoids @kaufman1990finding or HDBSCAN @campello2013hdbscan],
    [Non-Convex Manifolds & Rings], [Low: convex-hull midpoints land off-manifold], [Use WPGMA (Alg. 60), DBSCAN @ester1996dbscan, or Spectral @shi2000ncut],
    [Arbitrary Non-Euclidean Metrics], [Low: Apollonius identity requires Hilbert geometry], [Use WPGMA (Alg. 60) or DIANA (Alg. 63)],
  ),
  caption: [Domain suitability and failure-mode matrix for Median Linkage (WPGMC).],
)

#sec("V", "Structural Relations and Empirical Behavior")

WPGMC completes the $2 times 2$ combinatorial matrix of Sokal--Michener--Gower linkages indexed by (*unweighted* vs. *weighted* leaf schemes) $times$ (*average pairwise* vs. *centroid* geometry): UPGMA, WPGMA (Alg. 60), UPGMC (Alg. 61), and WPGMC (Alg. 62). The following proposition establishes two exact equivalence regimes connecting WPGMC to WPGMA and UPGMC.

#proposition(title: "Exact 1D Ordered Equivalence with WPGMA and Balanced-Tree Coincidence with UPGMC", origin: "Original Proof")[
  1. *1D ordered coincidence with WPGMA*: On $X subset bb(R)^1$ with $d(x, y) = |x - y|$, whenever two active clusters $A, B$ are linearly separated ($max(A) <= min(B)$) and share identical internal merge subtrees $T_A, T_B$, Euclidean WPGMC coincides identically with Euclidean WPGMA (Alg. 60): $d_M (A, B) = d_W (A, B)$.
  2. *Balanced-tree coincidence with UPGMC*: Whenever $T_A$ and $T_B$ are perfect binary trees (all leaves at equal depth inside each subtree), $bold(m)_A = bold(mu)_A$ and $bold(m)_B = bold(mu)_B$, so $D_M (A, B) = D_C (A, B)$ identically.
]

#proof[
  (1) Since $y - x = |y - x| >= 0$ for all $x in A, y in B$, Lemma 1 of WPGMA and Lemma 1 of WPGMC yield
  $ d_W (A, B) &= sum_(x in A) sum_(y in B) mu_A (x) mu_B (y) (y - x) \ &= sum_(y in B) mu_B (y) y - sum_(x in A) mu_A (x) x = m_B - m_A = d_M (A, B), $
  using $sum_(x in A) mu_A (x) = sum_(y in B) mu_B (y) = 1$.
  (2) If every leaf $x in A$ lies at identical depth $a = log_2 |A|$, then $mu_A (x) = 2^(-a) = 1 / |A|$, so $bold(m)_A = |A|^(-1) sum_(bold(x) in A) bold(x) = bold(mu)_A$, and symmetrically $bold(m)_B = bold(mu)_B$, giving $D_M (A, B) = D_C (A, B)$.
]

*Empirical Behavior on the ClusterOverview Benchmark*: Across the ClusterOverview 4-dataset benchmark (Silhouette @rousseeuw1987silhouettes and ARI @hubert1985comparing), Median Linkage (WPGMC) records *Iris* $"Sil" = 0.426$, $"ARI" = 0.653$; *Customers* $"Sil" = 0.353$; *Multi-Shape* $"Sil" = 0.452$ ($"ARI" =$ N/A); and *High-D* $"Sil" = 0.234$ ($"ARI" =$ N/A), with a mean runtime of $2.225$ ms. Comparing WPGMC directly against UPGMC (Alg. 61) reveals a textbook demonstration of our theoretical results:
- *Highest Iris ARI (`0.653`) via Anti-Masking*: Although WPGMC's *Iris* Silhouette (`0.426`) is lower than UPGMC's (`0.489`), WPGMC achieves the *highest* external label accuracy (*Iris* $"ARI" = 0.653$) across all four hierarchical methods in this block (surpassing WPGMA `0.565`, UPGMC `0.560`, and DIANA `0.574`). Because *Iris* is outlier-free and exhibits an asymmetric overlap between *Versicolor* and *Virginica*, WPGMC's equal-branch midpoint rule (Proposition 1(1)) avoids UPGMC's large-cluster masking and tracks the true inter-species boundary significantly better.
- *Coincidence on Customers, Multi-Shape, and High-D*: WPGMC reproduces UPGMC's exact Silhouette scores on *Customers* (`0.353`), *Multi-Shape* (`0.452`), and *High-D* (`0.234`). On well-balanced top-level merges, Proposition 3(2) shows $bold(m)_A approx bold(mu)_A$; on *Multi-Shape*, both centroid methods suffer identically from off-manifold convex-hull representatives compared to WPGMA (`0.609`); and on *High-D*, distance concentration @beyer1999nearest and simplex inversions (Theorem 2(3)) compress cluster separability uniformly. The $2.225$ ms runtime slightly exceeds UPGMC (`2.104` ms), consistent with WPGMC's higher inversion frequency (Theorem 2(2)) triggering extra priority-queue refreshes.
